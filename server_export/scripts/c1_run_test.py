#!/usr/bin/env python3
"""
c1_run_test.py -- the realised C1-UCEC confirmatory test, exactly as frozen in
review/C1_UCEC_FROZEN_RULE_2026-09-15.md (with Addenda 1-2). This is the one real
analysis; it is not a pre-check and it has no tunable reading.

Structure, and how it differs from split_replication_power.py (which this imports):
  selection  = FIXED, already computed on the discovery cohort (n=100). The selected
               set is `fdr < 0.05 & incremental_r2 > 0` in the pinned discovery result
               file: 2,566 of 10,512 tested genes (both counts asserted below -- a
               different vintage of that file stops the run). Nothing is re-derived
               from discovery patients here (frozen rule section 2).
  test       = the WHOLE confirmatory cohort (no split). PCA(K) is fit on the
               confirmatory cohort's OWN WSI embeddings, K = round(20 * n / 100)
               via SP.k_for -- the formula is frozen, the resulting K is not
               (section 3).
  statistic  = AUC (primary) and Spearman rho vs the discovery increments
               (secondary), one-sided upper tail, computed with the same rank code
               as the pre-checks (section 3).
  null       = B = 1000 permutations of the patient <-> WSI-PC-row correspondence,
               every gene recomputed together inside each permutation so gene-gene
               correlation is preserved (section 4).
  reading    = p(AUC) < 0.05 replicated / 0.05-0.15 suggestive / >= 0.15 not
               replicated; < 45 usable cases -> abandoned, not run (section 5).

Per-gene increments, the CV, the ridge and the permutation machinery are imported
from split_replication_power.py / residual_analysis.py rather than reimplemented, so
"identical machinery" is enforced by construction. What is NOT imported is
RA.load_protein_matrix: the confirmatory proteome comes from PDC000439's
quantDataMatrix as a case-indexed CSV (frozen rule section 1 specifies that pull
method), not as the discovery cohort's raw "<aliquot> Log Ratio" TSV.

Modes:
  test       observed statistic + permutation null + gene-level table   (default)
  bootstrap  patient-level bootstrap 95% CI for the AUC (section 5). Separate mode
             because each resample costs one full per-gene pass, i.e. about what the
             1000 permutations cost; it is run after the headline test, on the same
             data, and cannot change the section-5 reading.

Run (smp, 16 cores; env comes from lsf_c1_test.sh so no MORPHO_* leaks in):
  bsub -q smp -n 16 -R "span[hosts=1]" -o c1_test.out lsf_c1_test.sh test --perms 1000
  bsub -q smp -n 16 -R "span[hosts=1]" -o c1_boot.out lsf_c1_test.sh bootstrap --boots 1000
Smoke first (minutes): lsf_c1_test.sh test --perms 5 --tag smoke
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
import split_replication_power as SP  # noqa: E402  (same _gene_test / run_pool / k_for)

N_TESTED_EXPECTED = 10512   # frozen rule section 2
N_SELECTED_EXPECTED = 2566  # frozen rule section 2
MIN_USABLE_CASES = 45       # frozen rule section 5: below this the test is abandoned
PCS_PER_100 = 20.0          # frozen rule section 3
MIN_N = 30                  # same floor as residual_analysis._process_gene


def sha256(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return "unreadable"


def load_confirmatory_protein(path, log):
    df = pd.read_csv(path, index_col=0).apply(pd.to_numeric, errors="coerce")
    log(f"protein (PDC000439 quantDataMatrix, case-indexed): "
        f"{df.shape[0]} cases x {df.shape[1]} genes  <- {path}")
    return df


def build_gene_sets(discovery_csv, confirm_genes, log):
    d = pd.read_csv(discovery_csv)
    if not d.gene.is_unique:
        sys.exit(f"FATAL: {discovery_csv} has duplicate gene rows -- the selection set would "
                 f"be ambiguous. Stopping.")
    if len(d) != N_TESTED_EXPECTED:
        sys.exit(f"FATAL: discovery result file has {len(d)} tested genes, frozen rule says "
                 f"{N_TESTED_EXPECTED}. Wrong file vintage -- stopping rather than testing a "
                 f"different selection set than the one pre-registered.")
    sel = d[(d.fdr < 0.05) & (d.incremental_r2 > 0)]
    if len(sel) != N_SELECTED_EXPECTED:
        sys.exit(f"FATAL: selection set is {len(sel)} genes, frozen rule says "
                 f"{N_SELECTED_EXPECTED}. Stopping.")
    log(f"discovery selection set: {len(sel)}/{len(d)} genes "
        f"({len(sel) / len(d):.1%}) <- {discovery_csv}")
    disc_incr = d.set_index("gene")["incremental_r2"]
    disc_fdr = d.set_index("gene")["fdr"]
    tested = [g for g in d.gene if g in confirm_genes]
    selected = set(sel.gene) & set(tested)
    log(f"usable universe (discovery-tested AND present in confirmatory RNA+protein): "
        f"{len(tested)} genes, of which {len(selected)} selected / "
        f"{len(tested) - len(selected)} background")
    return tested, selected, disc_incr, disc_fdr


def assemble(M, xs, sel_mask):
    """AUC and Spearman rho for column 0 (observed) and columns 1.. (null draws).

    Copied from split_replication_power.py so the statistic is computed by the same
    code as every pre-check: rank-based Mann-Whitney AUC, and Spearman as the mean
    product of standardised ranks.
    """
    n1, n0 = int(sel_mask.sum()), int((~sel_mask).sum())
    ranks = rankdata(M, axis=0)
    auc = (ranks[sel_mask].sum(axis=0) - n1 * (n1 + 1) / 2.0) / (n1 * n0)
    rs = rankdata(xs)
    rs = (rs - rs.mean()) / rs.std()
    rc = (ranks - ranks.mean(axis=0)) / ranks.std(axis=0)
    rho = (rs[:, None] * rc).mean(axis=0)
    return auc, rho, n1, n0


def prepare(args, log):
    rna = RA.load_rna_matrix()
    wsi = RA.load_wsi_embeddings()
    protein = load_confirmatory_protein(args.protein, log)

    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    log(f"cases with all three modalities: {len(common)}")
    if len(common) < MIN_USABLE_CASES:
        sys.exit(f"ABANDONED per frozen rule section 5: only {len(common)} usable cases "
                 f"(< {MIN_USABLE_CASES}). Report the availability finding; do not run the test.")

    disc_cases = set()
    if os.path.exists(args.discovery_slide_map):
        dm = pd.read_csv(args.discovery_slide_map, sep="\t")
        disc_cases = set(dm["case_submitter_id"])
    overlap = sorted(set(common) & disc_cases)
    if overlap:
        sys.exit(f"FATAL: {len(overlap)} confirmatory cases also appear in the discovery "
                 f"slide map ({overlap[:5]}) -- the independence claim in section 1 fails.")
    log(f"case-set independence vs discovery: 0 overlapping cases (checked against "
        f"{len(disc_cases)} discovery cases)")

    genes_present = set(rna.columns) & set(protein.columns)
    tested, selected, disc_incr, disc_fdr = build_gene_sets(args.discovery_results,
                                                            genes_present, log)

    R = rna.loc[common, tested].values.astype(float)
    P = protein.loc[common, tested].values.astype(float)
    W = wsi.loc[common].values.astype(float)
    K = SP.k_for(len(common), PCS_PER_100)
    log(f"capacity: K = round({PCS_PER_100:g} * {len(common)} / 100) = {K} PCs, "
        f"PCA fit on the confirmatory cohort's own embeddings ({W.shape[1]} dims)")
    return common, tested, selected, disc_incr, disc_fdr, R, P, W, K


def run_test(args, log):
    t0 = time.time()
    common, tested, selected, disc_incr, disc_fdr, R, P, W, K = prepare(args, log)
    n = len(common)

    pcs = PCA(n_components=K, svd_solver="full").fit_transform(W)
    rng = np.random.RandomState(zlib.crc32(b"c1_ucec|perm"))
    perms = [rng.permutation(n) for _ in range(args.perms)]
    log(f"per-gene pass: {len(tested)} genes x (1 observed + {args.perms} permutations) "
        f"on {args.workers} workers ...")
    res = SP.run_pool(SP._gene_test, len(tested), args.workers, R, P, pcs, perms, MIN_N)

    usable = [j for j in range(len(tested)) if res[j] is not None
              and np.isfinite(disc_incr.get(tested[j], np.nan))]
    dropped = len(tested) - len(usable)
    if dropped:
        log(f"{dropped} genes dropped (fewer than {MIN_N} confirmatory protein values)")
    if len(usable) < 100:
        sys.exit(f"FATAL: only {len(usable)} usable genes")

    genes_u = [tested[j] for j in usable]
    M = np.vstack([res[j] for j in usable])
    xs = np.array([disc_incr[g] for g in genes_u])
    sel_mask = np.array([g in selected for g in genes_u])
    if sel_mask.sum() == 0 or (~sel_mask).sum() == 0:
        sys.exit(f"FATAL: degenerate gene sets after intersection "
                 f"({int(sel_mask.sum())} selected, {int((~sel_mask).sum())} background)")

    auc, rho, n1, n0 = assemble(M, xs, sel_mask)
    auc_obs, auc_null = auc[0], auc[1:]
    rho_obs, rho_null = rho[0], rho[1:]
    p_auc = float((np.sum(auc_null >= auc_obs) + 1) / (len(auc_null) + 1))
    p_rho = float((np.sum(rho_null >= rho_obs) + 1) / (len(rho_null) + 1))
    sd = float(auc_null.std())

    out = dict(
        cohort="ucec_c1", n_cases=n, K=K, n_genes_tested=len(tested), n_usable=len(usable),
        n_selected=n1, n_background=n0, perms=args.perms,
        auc_obs=float(auc_obs), auc_null_mean=float(auc_null.mean()), auc_null_sd=sd,
        z_auc=float((auc_obs - auc_null.mean()) / sd) if sd > 0 else float("nan"),
        p_auc=p_auc,
        rho_obs=float(rho_obs), rho_null_mean=float(rho_null.mean()),
        rho_null_sd=float(rho_null.std()), p_rho=p_rho,
        confirm_negative_frac=float((M[:, 0] < 0).mean()),
        sha256_discovery_results=sha256(args.discovery_results),
        sha256_residual_analysis=sha256(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     "residual_analysis.py")),
        sha256_split_replication_power=sha256(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "split_replication_power.py")),
        minutes=round((time.time() - t0) / 60.0, 1),
    )

    suffix = f"_{args.tag}" if args.tag else ""
    gene_tbl = pd.DataFrame({
        "gene": genes_u,
        "selected": sel_mask,
        "discovery_incremental_r2": xs,
        "discovery_fdr": [disc_fdr[g] for g in genes_u],
        "confirmatory_incremental_r2": M[:, 0],
    }).sort_values("confirmatory_incremental_r2", ascending=False)
    gene_fp = RA.OUT_DIR / f"c1_ucec_gene_level{suffix}.csv"
    gene_tbl.to_csv(gene_fp, index=False)
    res_fp = RA.OUT_DIR / f"c1_ucec_result{suffix}.json"
    with open(res_fp, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    log("\n" + "=" * 74)
    log("C1-UCEC CONFIRMATORY TEST (frozen rule 2026-09-15, Addenda 1-2)")
    log("=" * 74)
    log(f"  cases / K                : {n} / {K}")
    log(f"  genes usable             : {len(usable)}  ({n1} selected, {n0} background)")
    log(f"  AUC observed             : {auc_obs:.4f}")
    log(f"  AUC null (B={args.perms})".ljust(29) + f": {auc_null.mean():.4f} +/- {sd:.4f}"
        f"   [the null is NOT 0.5 -- do not use 0.5 as a baseline]")
    log(f"  z (AUC)                  : {out['z_auc']:.2f}")
    log(f"  p (AUC, one-sided)       : {p_auc:.4f}   <- PRIMARY")
    log(f"  rho observed / null      : {rho_obs:.4f} / {rho_null.mean():.4f} "
        f"+/- {rho_null.std():.4f}")
    log(f"  p (rho, one-sided)       : {p_rho:.4f}   (secondary, not a second bar)")
    log(f"  confirmatory genes with negative increment: {out['confirm_negative_frac']:.1%}")
    p_floor = 1.0 / (args.perms + 1)
    if args.perms < 1000:
        # A short run cannot produce the frozen reading: its smallest attainable p is
        # 1/(B+1), so a p at that floor means "no permutation beat the observation",
        # NOT "not replicated". This project already published-then-caught that exact
        # censoring error once (transcriptome permutation FDR, PDAC's 0/1572).
        at_floor = p_auc <= p_floor + 1e-12
        verdict = (f"NO READING -- B={args.perms} is not the frozen analysis (section 4 fixes "
                   f"B=1000). The smallest attainable p here is {p_floor:.4f}; the observed "
                   f"p={p_auc:.4f} is "
                   + ("AT that floor, i.e. no permutation reached the observed statistic "
                      "(censored, and emphatically not a negative result)."
                      if at_floor else "resolution-limited.")
                   + " Section 5 applies only to the B=1000 run.")
    elif p_auc < 0.05:
        verdict = ("REPLICATED (p < 0.05). Report AUC + 95% bootstrap CI (run --mode "
                   "bootstrap), rho, and the gene-level table as the pre-declared external "
                   "replication; it becomes the headline external-validation result.")
    elif p_auc < 0.15:
        verdict = ("SUGGESTIVE (0.05 <= p < 0.15). NOT declared as replication. Report point "
                   "estimate, CI and p as-is; the TCGA-RPPA line stays the primary external-"
                   "validation narrative.")
    else:
        verdict = ("NOT REPLICATED (p >= 0.15). Report plainly alongside the RPPA/BRCA nulls. "
                   "The pre-checks' batch-blocked power stays in the record; at n=138 this "
                   "failure is itself a finding, and is not to be spun.")
    log(f"\n  PRE-DECLARED READING     : {verdict}")
    log(f"\n  -> {res_fp}\n  -> {gene_fp}\n  elapsed {out['minutes']} min")


def run_bootstrap(args, log):
    t0 = time.time()
    common, tested, selected, disc_incr, _fdr, R, P, W, K = prepare(args, log)
    n = len(common)
    log(f"patient-level bootstrap: {args.boots} resamples of {n} cases (with replacement), "
        f"PCA refit per resample, K held at {K} (the frozen formula's value for the full "
        f"cohort). Caveat recorded with the result: duplicated patients can fall in both "
        f"the train and the test fold of the 5-fold CV, which inflates absolute R2; the AUC "
        f"is a contrast between two gene sets computed on the same resample, so the effect "
        f"largely cancels, but the CI is not bias-free.")

    base_usable = None
    rng = np.random.RandomState(zlib.crc32(b"c1_ucec|boot"))
    aucs = []
    for b in range(args.boots):
        idx = rng.randint(0, n, n)
        pcs_b = PCA(n_components=K, svd_solver="full").fit_transform(W[idx])
        res = SP.run_pool(SP._gene_test, len(tested), args.workers, R[idx], P[idx], pcs_b,
                          [], MIN_N)
        usable = [j for j in range(len(tested)) if res[j] is not None
                  and np.isfinite(disc_incr.get(tested[j], np.nan))]
        if base_usable is None:
            base_usable = len(usable)
        if len(usable) < 100:
            continue
        genes_u = [tested[j] for j in usable]
        M = np.vstack([res[j] for j in usable])
        xs = np.array([disc_incr[g] for g in genes_u])
        sel_mask = np.array([g in selected for g in genes_u])
        auc, _rho, _n1, _n0 = assemble(M, xs, sel_mask)
        aucs.append(float(auc[0]))
        if (b + 1) % 25 == 0:
            log(f"  {b + 1}/{args.boots} resamples, running AUC mean "
                f"{np.mean(aucs):.4f}, {(time.time() - t0) / 60:.0f} min")

    a = np.array(aucs)
    lo, hi = np.percentile(a, [2.5, 97.5])
    out = dict(cohort="ucec_c1", n_cases=n, K=K, boots=len(a),
               auc_boot_mean=float(a.mean()), auc_boot_sd=float(a.std()),
               auc_ci95_lo=float(lo), auc_ci95_hi=float(hi),
               minutes=round((time.time() - t0) / 60.0, 1))
    suffix = f"_{args.tag}" if args.tag else ""
    fp = RA.OUT_DIR / f"c1_ucec_bootstrap{suffix}.json"
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    log(f"\n  AUC bootstrap mean       : {a.mean():.4f} (sd {a.std():.4f}, {len(a)} resamples)")
    log(f"  95% CI (percentile)      : [{lo:.4f}, {hi:.4f}]")
    log(f"  -> {fp}\n  elapsed {out['minutes']} min")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="test", choices=["test", "bootstrap"])
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--boots", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=RA.N_WORKERS)
    ap.add_argument("--tag", default="")
    ap.add_argument("--protein", default="/public/home/fjhui/ZW/ucec_c1/protein_ucec_c1.csv")
    ap.add_argument("--discovery-results",
                    default="/public/home/fjhui/ZW/scripts/pinned/ucec/residual_results_tumoronly.csv")
    ap.add_argument("--discovery-slide-map",
                    default="/public/home/fjhui/ZW/scripts/pinned/slide_type_map_ucec.tsv")
    args = ap.parse_args()

    def log(msg):
        print(msg, flush=True)

    log("=== resolved paths (verify before trusting any number below) ===")
    log(f"  MORPHO_ROOT        : {RA.ROOT}")
    log(f"  MORPHO_RNA_MANIFEST: {RA.RNA_MANIFEST}")
    log(f"  MORPHO_SLIDE_MAP   : {RA.SLIDE_TYPE_MAP}")
    log(f"  MORPHO_WSI_EMB_DIR : {RA.WSI_EMB_DIR}")
    log(f"  MORPHO_OUT         : {RA.OUT_DIR}")
    log(f"  protein            : {args.protein}")
    log(f"  discovery results  : {args.discovery_results}")
    log(f"  mode / workers     : {args.mode} / {args.workers}")
    log("=" * 74)

    if args.mode == "test":
        run_test(args, log)
    else:
        run_bootstrap(args, log)


if __name__ == "__main__":
    main()
