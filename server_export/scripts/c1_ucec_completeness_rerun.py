#!/usr/bin/env python3
"""
c1_ucec_completeness_rerun.py -- POST HOC completeness-stratified rerun of the
already-published C1-UCEC confirmatory test.

Authority: review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 0 (pre-signature
decisions recorded 2026-09-29), decision D4 ("Already-published C1-UCEC: two post
hoc items, decided before any LUAD result exists") and the completeness-stratified
AUC wording fixed under D1's "Co-reported" bullet:

    Completeness-stratified AUC. Recompute with its permutation null by a
    deterministic rerun: same seed crc32(b"c1_ucec|perm"), same B = 1000, same
    cases, genes and K. The observed value from the published gene table is 0.561,
    against the pooled 0.597. Label it post hoc and report it in Results whatever
    the result.

This script does NOT re-derive the C1-UCEC test. It reproduces c1_run_test.py's
run_test() EXACTLY up to the genes x (1+perms) matrix M -- same prepare(), same
PCA call, same permutation seed, same SP.run_pool call, same usable filter, same
genes_u order, same sel_mask -- by importing c1_run_test and split_replication_power
rather than reimplementing them, so "identical machinery" is enforced by
construction (mirrors c1_run_test.py's own docstring convention). It then:

  1. VALIDATES the reproduced auc_obs / auc_null_mean / auc_null_sd / p_auc /
     rho_obs / p_rho / n_usable / n_selected / n_background against the already-
     published c1_ucec_result.json (abs tol 1e-9 for floats, exact for counts).
     FATAL on any mismatch, printing both sets of values, UNLESS --no-validate
     (which exists only for the synthetic dry run -- see the wrapper / dry-run
     script for why real validation cannot pass on synthetic data).
  2. Splits the usable genes into two strata by their DISCOVERY per-gene n (from
     --discovery-results' own "n" column): "complete" = n equals that file's modal
     n (printed), "incomplete" = the rest.
  3. Computes AUC_strat = sum_s U_s / sum_s (n1_s * n0_s), ranks computed WITHIN
     each stratum, for every column of M (observed + all permutation nulls), plus
     each stratum's own AUC with its own null (same B=1000 columns), stratum sizes
     (selected/background) and the fraction of confirmatory increments > 0 per
     stratum.
  4. Saves the genes x (1+B) matrix M, sel_mask and strata to an .npz, and writes
     c1_ucec_completeness_result.json (post_hoc: true, decision_date, npz sha256,
     this script's sha256, and the same code/discovery hashes c1_run_test.py
     records).

Nothing here can change the frozen §5 reading of the original test; this is a
descriptive, pre-declared-wording, post hoc addendum only (D1 "Co-reported", not a
gate).

Run (smp; env comes from lsf_c1_ucec_strat.sh so no MORPHO_* leaks in):
  bsub -q smp -n 8 -R "span[hosts=1]" -o c1_ucec_strat.out \
    lsf_c1_ucec_strat.sh --perms 1000
Smoke first (minutes; validates auc_obs only, refuses to print any stratified p):
  lsf_c1_ucec_strat.sh --perms 5 --tag smoke
"""
import argparse
import hashlib
import json
import os
import sys
import time
import zlib

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import residual_analysis as RA  # noqa: E402
import split_replication_power as SP  # noqa: E402
import c1_run_test as C1  # noqa: E402  (prepare, assemble, MIN_N -- not reimplemented)

FLOAT_TOL = 1e-9


def sha256(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return "unreadable"


def auc_strat(values, sel_mask, strata):
    """Stratified rank-sum AUC: sum_s U_s / sum_s (n1_s * n0_s), with ranks computed
    WITHIN each stratum (not globally). Pure function -- no I/O, no globals -- so it
    can be unit-tested directly.

    values : (n_genes,) or (n_genes, n_cols) array to rank (e.g. the M matrix: column
             0 observed, columns 1.. permutation nulls).
    sel_mask : (n_genes,) bool array, True = selected gene.
    strata : (n_genes,) array of stratum labels (any hashable, e.g. "complete" /
             "incomplete").

    Returns (auc, per_stratum):
      auc         : float (if values was 1-D) or (n_cols,) array -- the pooled
                    stratified statistic, sum_s U_s / sum_s (n1_s * n0_s).
      per_stratum : {label: {"n1": int, "n0": int, "auc": float | (n_cols,) array |
                    None}} -- None when a stratum has zero selected or zero
                    background genes (den contribution 0, excluded from the pooled
                    sum, exactly as the sum-of-U-over-sum-of-n1n0 formula implies).
    """
    values = np.asarray(values, dtype=float)
    sel_mask = np.asarray(sel_mask, dtype=bool)
    strata = np.asarray(strata)
    if len(values) != len(sel_mask) or len(values) != len(strata):
        raise ValueError(
            f"length mismatch: values={len(values)} sel_mask={len(sel_mask)} "
            f"strata={len(strata)}")
    squeeze = values.ndim == 1
    if squeeze:
        values = values[:, None]
    n_cols = values.shape[1]
    num = np.zeros(n_cols, dtype=float)
    den = 0.0
    per_stratum = {}
    for lab in sorted(set(strata.tolist())):
        idx = strata == lab
        n1 = int(sel_mask[idx].sum())
        n0 = int((~sel_mask[idx]).sum())
        if n1 == 0 or n0 == 0:
            per_stratum[lab] = dict(n1=n1, n0=n0, auc=None)
            continue
        ranks = rankdata(values[idx], axis=0)
        U = ranks[sel_mask[idx]].sum(axis=0) - n1 * (n1 + 1) / 2.0
        num += U
        den += n1 * n0
        auc_s = U / (n1 * n0)
        per_stratum[lab] = dict(n1=n1, n0=n0,
                                 auc=float(auc_s[0]) if squeeze else auc_s)
    auc = num / den if den > 0 else np.full(n_cols, np.nan)
    if squeeze:
        auc = float(auc[0])
    return auc, per_stratum


def emp_p(obs, null):
    """One-sided upper-tail empirical p-value, same convention as c1_run_test.py /
    split_replication_power.py: (1 + #{null >= obs}) / (1 + B)."""
    null = np.asarray(null, dtype=float)
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))


# The published c1_ucec_result.json (2026-09-16 run). Pinned literally so a later,
# possibly-different regeneration of that file cannot silently change what this
# rerun is checked against -- the whole point of "validate against the published
# result" is that the reference is fixed, not whatever currently sits on disk.
PUBLISHED_RESULT_SHA256 = "53fcb86138465e90ef2c47e3519786fc0fb4a534497f4a83850ee584ee6e653c"


def validate(log, computed, published_path, tol=FLOAT_TOL):
    """FATAL exit (printing both sets of values) unless every field of `computed`
    matches the corresponding field of the published JSON at `published_path`,
    AND the published file itself is byte-identical to the 2026-09-16 run this
    script is meant to reproduce (checked by sha256, not just by its field values)."""
    if not os.path.exists(published_path):
        sys.exit(f"FATAL: --validate-against file does not exist: {published_path}")
    got_file_hash = sha256(published_path)
    if got_file_hash != PUBLISHED_RESULT_SHA256:
        sys.exit(f"FATAL: {published_path} has sha256={got_file_hash}, expected the "
                 f"pinned 2026-09-16 result {PUBLISHED_RESULT_SHA256}. Refusing to "
                 f"validate against a file that is not byte-identical to the "
                 f"published C1-UCEC result -- if it was legitimately regenerated, "
                 f"update this literal deliberately, do not just let it drift.")
    with open(published_path, "r", encoding="utf-8") as f:
        published = json.load(f)

    float_fields = ["auc_obs", "auc_null_mean", "auc_null_sd", "p_auc",
                    "rho_obs", "p_rho"]
    count_fields = [("n_usable", "n_usable"), ("n_selected", "n_selected"),
                    ("n_background", "n_background"),
                    ("n_cases", "n_cases"), ("K", "K"),
                    ("n_genes_tested", "n_genes_tested")]
    hash_fields = ["sha256_discovery_results", "sha256_residual_analysis",
                   "sha256_split_replication_power"]

    mismatches = []
    for f in float_fields:
        got, want = computed[f], published[f]
        if not (np.isfinite(got) and np.isfinite(want) and abs(got - want) <= tol):
            mismatches.append((f, got, want))
    for got_key, want_key in count_fields:
        got, want = computed[got_key], published[want_key]
        if got != want:
            mismatches.append((got_key, got, want))
    for f in hash_fields:
        got, want = computed[f], published[f]
        if got != want:
            mismatches.append((f, got, want))

    if mismatches:
        log("\nFATAL: reproduced statistics do NOT match the published result -- "
            f"stopping before any stratified analysis is run. ({published_path})")
        log(f"{'field':<28}{'this run':>24}{'published':>24}")
        for f, got, want in mismatches:
            log(f"{f:<28}{got!r:>24}{want!r:>24}")
        sys.exit(1)
    log(f"VALIDATED: all {len(float_fields) + len(count_fields) + len(hash_fields)} "
        f"statistics and code/discovery hashes match {published_path} "
        f"(abs tol {tol:g} for floats, exact for counts and hashes; file itself "
        f"sha256-pinned to {PUBLISHED_RESULT_SHA256[:12]}...).")


def main():
    ap = argparse.ArgumentParser()
    # Identical defaults to c1_run_test.py's ArgumentParser (mode "test" branch) --
    # copied, not derived, so a change to c1_run_test.py's defaults does not
    # silently change this script's without a matching edit here.
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=RA.N_WORKERS)
    ap.add_argument("--tag", default="")
    ap.add_argument("--protein", default="/public/home/fjhui/ZW/ucec_c1/protein_ucec_c1.csv")
    ap.add_argument("--discovery-results",
                    default="/public/home/fjhui/ZW/scripts/pinned/ucec/residual_results_tumoronly.csv")
    ap.add_argument("--discovery-slide-map",
                    default="/public/home/fjhui/ZW/scripts/pinned/slide_type_map_ucec.tsv")
    # New in this script only.
    ap.add_argument("--validate-against", default="",
                     help="published c1_ucec_result.json to check the reproduced "
                          "statistics against before running the stratified "
                          "analysis. Default: <MORPHO_OUT>/c1_ucec_result.json -- "
                          "under an identical env to lsf_c1_test.sh (which this "
                          "wrapper's env block matches) that IS the original run's "
                          "output directory.")
    ap.add_argument("--no-validate", action="store_true",
                     help="skip the validate-against-published-JSON step. ONLY for "
                          "the synthetic dry run (dryrun_heldout.py), where the "
                          "monkeypatched loaders make the real published numbers "
                          "unreachable by construction. Never pass this against "
                          "real data.")
    args = ap.parse_args()

    def log(msg):
        print(msg, flush=True)

    validate_against = args.validate_against or str(RA.OUT_DIR / "c1_ucec_result.json")

    log("=== resolved paths (verify before trusting any number below) ===")
    log(f"  MORPHO_ROOT        : {RA.ROOT}")
    log(f"  MORPHO_RNA_MANIFEST: {RA.RNA_MANIFEST}")
    log(f"  MORPHO_SLIDE_MAP   : {RA.SLIDE_TYPE_MAP}")
    log(f"  MORPHO_WSI_EMB_DIR : {RA.WSI_EMB_DIR}")
    log(f"  MORPHO_OUT         : {RA.OUT_DIR}")
    log(f"  protein            : {args.protein}")
    log(f"  discovery results  : {args.discovery_results}")
    log(f"  validate against   : {validate_against}"
        + ("  [SKIPPED: --no-validate]" if args.no_validate else ""))
    log(f"  perms / workers    : {args.perms} / {args.workers}")
    log("=" * 74)

    t0 = time.time()

    # ---- reproduce c1_run_test.run_test() exactly up to M -------------------
    common, tested, selected, disc_incr, disc_fdr, R, P, W, K = C1.prepare(args, log)
    n = len(common)

    pcs = PCA(n_components=K, svd_solver="full").fit_transform(W)
    rng = np.random.RandomState(zlib.crc32(b"c1_ucec|perm"))
    perms = [rng.permutation(n) for _ in range(args.perms)]
    log(f"per-gene pass: {len(tested)} genes x (1 observed + {args.perms} permutations) "
        f"on {args.workers} workers ...")
    res = SP.run_pool(SP._gene_test, len(tested), args.workers, R, P, pcs, perms, C1.MIN_N)

    usable = [j for j in range(len(tested)) if res[j] is not None
              and np.isfinite(disc_incr.get(tested[j], np.nan))]
    dropped = len(tested) - len(usable)
    if dropped:
        log(f"{dropped} genes dropped (fewer than {C1.MIN_N} confirmatory protein values)")
    if len(usable) < 100:
        sys.exit(f"FATAL: only {len(usable)} usable genes")

    genes_u = [tested[j] for j in usable]
    M = np.vstack([res[j] for j in usable])
    xs = np.array([disc_incr[g] for g in genes_u])
    sel_mask = np.array([g in selected for g in genes_u])
    if sel_mask.sum() == 0 or (~sel_mask).sum() == 0:
        sys.exit(f"FATAL: degenerate gene sets after intersection "
                 f"({int(sel_mask.sum())} selected, {int((~sel_mask).sum())} background)")
    # ---- end exact reproduction ----------------------------------------------

    auc, rho, n1, n0 = C1.assemble(M, xs, sel_mask)
    auc_obs, auc_null = auc[0], auc[1:]
    rho_obs, rho_null = rho[0], rho[1:]
    # Computed here (not at the end) so validate() can check n_cases/K/n_genes_tested
    # and the three code/discovery hashes against the published JSON's own fields,
    # not just the 6 floats + 3 counts the original version checked.
    this_script_hash = sha256(os.path.abspath(__file__))
    discovery_hash = sha256(args.discovery_results)
    residual_analysis_hash = sha256(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "residual_analysis.py"))
    split_power_hash = sha256(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "split_replication_power.py"))
    computed = dict(
        auc_obs=float(auc_obs), auc_null_mean=float(auc_null.mean()),
        auc_null_sd=float(auc_null.std()),
        p_auc=emp_p(auc_obs, auc_null),
        rho_obs=float(rho_obs), p_rho=emp_p(rho_obs, rho_null),
        n_usable=len(usable), n_selected=int(n1), n_background=int(n0),
        n_cases=n, K=K, n_genes_tested=len(tested),
        sha256_discovery_results=discovery_hash,
        sha256_residual_analysis=residual_analysis_hash,
        sha256_split_replication_power=split_power_hash,
    )
    log(f"\n  reproduced AUC observed  : {computed['auc_obs']:.6f}")
    log(f"  reproduced p (AUC)       : {computed['p_auc']:.6f}")

    p_floor = 1.0 / (args.perms + 1)
    smoke_short_run = args.perms < 1000
    if args.no_validate:
        log("\nVALIDATION SKIPPED (--no-validate): this run is not checked against "
            "any published result. Only ever pass this for the synthetic dry run.")
    elif smoke_short_run:
        # A short (smoke) run cannot match the B=1000 published null -- same
        # resolution-limit reasoning c1_run_test.py's own verdict block applies to
        # p_auc. Validate only auc_obs (deterministic given the same data/PCA/K,
        # independent of B) and refuse to print any stratified p, exactly as the
        # wrapper script's header promises.
        if not os.path.exists(validate_against):
            sys.exit(f"FATAL: --validate-against file does not exist: {validate_against}")
        with open(validate_against, "r", encoding="utf-8") as f:
            published = json.load(f)
        if abs(computed["auc_obs"] - published["auc_obs"]) > FLOAT_TOL:
            log(f"\nFATAL (smoke): auc_obs mismatch. this run={computed['auc_obs']!r} "
                f"published={published['auc_obs']!r}")
            sys.exit(1)
        log(f"SMOKE VALIDATED: auc_obs matches {validate_against} "
            f"(abs tol {FLOAT_TOL:g}). B={args.perms} smallest attainable p is "
            f"{p_floor:.4f} -- refusing to print any stratified p for a smoke run.")
    else:
        validate(log, computed, validate_against)

    # ---- completeness strata --------------------------------------------------
    disc = pd.read_csv(args.discovery_results)
    modal_n = int(disc["n"].mode().iloc[0])
    log(f"\ndiscovery file modal per-gene n = {modal_n} "
        f"(--discovery-results: {args.discovery_results})")
    # The stratum boundary is a data-derived mode, but the boundary this script must
    # actually use is the UCEC discovery cohort size (100), not whatever a different
    # file vintage happens to report as its mode. Assert them equal so a wrong-vintage
    # file cannot silently move where "complete" is drawn.
    if modal_n != 100:
        sys.exit(f"FATAL: discovery file's modal per-gene n is {modal_n}, expected "
                 f"100 (the UCEC discovery cohort size) -- wrong file vintage, "
                 f"stopping rather than stratifying at the wrong boundary.")
    n_by_gene = disc.set_index("gene")["n"]
    missing = [g for g in genes_u if g not in n_by_gene.index]
    if missing:
        sys.exit(f"FATAL: {len(missing)} usable genes have no row in "
                 f"{args.discovery_results} (e.g. {missing[:5]}) -- cannot stratify.")
    strata = np.array(["complete" if n_by_gene[g] == modal_n else "incomplete"
                       for g in genes_u])
    # Per-gene CONFIRMATORY n (not the discovery n above) -- D1's "Co-reported" bullet
    # asks for this alongside M so a reader can check completeness on either side.
    confirm_n = (~np.isnan(P[:, usable])).sum(axis=0).astype(int)

    auc_strat_all, strat_info = auc_strat(M, sel_mask, strata)
    auc_strat_obs, auc_strat_null = auc_strat_all[0], auc_strat_all[1:]

    per_stratum_out = {}
    for lab, info in sorted(strat_info.items()):
        frac_pos = float((M[strata == lab, 0] > 0).mean())
        entry = dict(n_selected=info["n1"], n_background=info["n0"],
                     frac_positive_confirmatory=frac_pos,
                     median_confirm_n=int(np.median(confirm_n[strata == lab])))
        if info["auc"] is not None:
            arr = info["auc"]
            entry["auc_obs"] = float(arr[0])
            entry["auc_null_mean"] = float(arr[1:].mean())
            entry["auc_null_sd"] = float(arr[1:].std())
        else:
            entry["auc_obs"] = None
        per_stratum_out[lab] = entry
        log(f"  stratum {lab:<11} n_sel={info['n1']:<6} n_bg={info['n0']:<6} "
            f"frac_positive={frac_pos:.4f} "
            + (f"AUC={entry['auc_obs']:.4f} (null {entry['auc_null_mean']:.4f} "
               f"+/- {entry['auc_null_sd']:.4f})" if info["auc"] is not None
               else "AUC=n/a (degenerate stratum)"))

    result = dict(
        post_hoc=True,
        decision_date="2026-09-29",
        authority="review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 0, D1 "
                  "Co-reported + D4",
        cohort="ucec_c1", n_cases=n, K=K, perms=args.perms,
        modal_discovery_n=modal_n,
        n_usable=computed["n_usable"], n_selected=computed["n_selected"],
        n_background=computed["n_background"],
        pooled_auc_obs=computed["auc_obs"],
        auc_strat_obs=float(auc_strat_obs),
    )

    if args.no_validate or smoke_short_run:
        log("\nNOT computing a stratified p-value: "
            + ("--no-validate run" if args.no_validate else
               f"B={args.perms} smoke run (see refusal above)."))
        result["p_auc_strat"] = None
        result["validated"] = (False if args.no_validate else "smoke_auc_obs_only")
    else:
        p_strat = emp_p(auc_strat_obs, auc_strat_null)
        result["auc_strat_null_mean"] = float(auc_strat_null.mean())
        result["auc_strat_null_sd"] = float(auc_strat_null.std())
        result["p_auc_strat"] = p_strat
        result["validated"] = True
        log(f"\n  AUC_strat observed       : {auc_strat_obs:.6f}")
        log(f"  AUC_strat null (B={args.perms})".ljust(29)
            + f": {auc_strat_null.mean():.6f} +/- {auc_strat_null.std():.6f}")
        log(f"  p (AUC_strat, one-sided) : {p_strat:.6f}")
        if computed["p_auc"] < 0.05 and p_strat >= 0.05:
            log("\n  FIXED WORDING (D1 Co-reported): \"replicated at the pooled "
                "level; not robust to protein-completeness stratification\".")

    result["per_stratum"] = per_stratum_out
    result["minutes"] = round((time.time() - t0) / 60.0, 1)

    suffix = f"_{args.tag}" if args.tag else ""
    npz_fp = RA.OUT_DIR / f"c1_ucec_completeness{suffix}.npz"
    np.savez_compressed(npz_fp,
                        genes=np.array(genes_u),
                        sel_mask=sel_mask,
                        strata=strata,
                        M=M,
                        confirm_n=confirm_n)
    result["sha256_npz"] = sha256(str(npz_fp))
    result["sha256_this_script"] = this_script_hash
    result["sha256_discovery_results"] = computed["sha256_discovery_results"]
    result["sha256_residual_analysis"] = computed["sha256_residual_analysis"]
    result["sha256_split_replication_power"] = computed["sha256_split_replication_power"]
    result["sha256_c1_run_test"] = sha256(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "c1_run_test.py"))

    res_fp = RA.OUT_DIR / f"c1_ucec_completeness_result{suffix}.json"
    with open(res_fp, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    log(f"\n  -> {npz_fp}\n  -> {res_fp}\n  elapsed {result['minutes']} min")


if __name__ == "__main__":
    main()
