#!/usr/bin/env python3
"""
c1_ucec_stratified_reread.py -- POST HOC C1-UCEC batch-stratified re-read, D4 item 3
step 3-4 (D6b).

Authority: review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 0, decision D4,
item 3, steps 3-4:

    3. Re-read. Recompute the published 2,566-set AUC under a within-operator
       null: B = 1,000, seed crc32(b"c1_ucec|perm_op"), with singleton operators
       and unlabelled cases held fixed.
    4. Reporting. It is labelled post hoc and reported whatever the result.
    - If the attribute does not validate, or any D1 criterion fails on the UCEC
      analysed set, Methods state that C1-UCEC was tested under the unrestricted
      null only and carries no batch qualifier.
    In every case, UCEC's published section 5 reading is unchanged.

This does NOT re-derive the C1-UCEC test. Exactly like c1_ucec_completeness_rerun.py
(the D1 "Co-reported" post hoc script this one is modelled on), it reproduces
c1_run_test.py's run_test() EXACTLY up to the genes x (1+perms) matrix M -- same
prepare(), same PCA call, same permutation seed, same SP.run_pool call, same usable
filter, same genes_u order, same sel_mask -- by importing c1_run_test and
split_replication_power rather than reimplementing them. It then:

  1. VALIDATES the reproduced auc_obs / auc_null_mean / auc_null_sd / p_auc /
     rho_obs / p_rho / n_usable / n_selected / n_background against the already-
     published c1_ucec_result.json (abs tol 1e-9 for floats, exact for counts and
     hashes). FATAL on any mismatch, UNLESS --no-validate (synthetic dry run only).
  2. Reads idc_ucec_derive_batch_labels.py's summary JSON (--d4-summary) for its
     "qualifies" verdict (validated_operator AND all four D1 criteria hold on the
     published 138-case analysed set). This script never re-derives that verdict;
     it only consumes it.
  3a. If qualifies: builds a within-operator permutation null (B, seed
      crc32(b"c1_ucec|perm_op")) from --batch-labels, with singleton operators AND
      unlabelled cases held fixed (never pooled into a permutable "_missing"
      level -- this is departure (i) from residual_analysis_sitepack.py's main(),
      reimplemented here exactly as c1_run_test_luad.py's
      build_within_operator_perms does for LUAD, just with the UCEC seed tag).
      Reruns SP.run_pool with this null on the SAME PCs, SAME 2,566-selection
      genes, SAME discovery increments as the validated reproduction -- so the
      recomputed AUC's OBSERVED value is asserted bit-identical to auc_obs above
      (only the null differs; same assertion pattern as c1_run_test_luad.py's
      (a)(i)-vs-H1 check). This IS D4 item 3 step 3, read literally: the rule's
      own text does not residualise any PC for this item (contrast LUAD's H2,
      which residualises PCs for a DIFFERENT, batch-corrected 726-gene selection
      set -- UCEC never had a batch-corrected selection set to test, so D4 item 3
      has no H2 analogue). p_auc_strat is this step's p-value.
  3b. ALSO computes a second, clearly-separated, NON-GATING descriptive variant
      that additionally residualises the confirmatory K PCs on operator dummies
      before recomputing under the same within-operator null -- an LUAD-(a)(ii)-
      style extra. This was requested by the task text that commissioned this
      script (which asked for `build_operator_dummy_design`-style residualisation
      "and recompute the published 2,566-gene AUC under a within-operator null"),
      but is NOT asked for by the rule's own D4 item 3 text (read in full above,
      and again in review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md before relying
      on this docstring's paraphrase). Because there is no rule-frozen dummy list
      for UCEC (D4 was decided 2026-09-29, for a re-read of an ALREADY-PUBLISHED
      result -- no frozen-12-operator table like c1_run_test_luad.py's
      FROZEN_DUMMY_OPERATORS was ever written for UCEC), this variant reimplements
      residual_analysis_sitepack.py's ORIGINAL DYNAMIC dummy-selection rule
      exactly -- Counter.most_common(MAX_LEVELS=12), kept if size >= MIN_GROUP=3,
      full-rank guard -- computed from the realised (labelled-only; unlabelled
      cases are still never a dummy level, departure (i) again) UCEC operator
      labels, not the frozen LUAD list. This is the ONE place in this project
      where reusing sitepack's dynamic selection (rather than a frozen list) is
      correct, precisely because no rule text ever froze a UCEC dummy set. Every
      output field for this variant is prefixed `extra_residualized_` and carries
      `"descriptive_not_required_by_d4_item3_text": true` in the result JSON, so
      it can never be mistaken for the primary (gating) re-read above. It CANNOT
      change D4 item 3's reading, does not enter the fallback-clause decision, and
      is reported only because it was asked for -- not because the rule requires
      it. Flag this explicitly to the corresponding author before this script's
      output is cited anywhere: confirm whether the extra variant belongs in the
      paper at all, or should be dropped to keep D4 item 3 to its literal text.
  4. If NOT qualifies: prints the section-0-mandated fallback sentence verbatim
     and exits cleanly (rc=0) WITHOUT computing any stratified AUC, per D4 item
     3's own fallback clause -- this is an anticipated outcome, not an error.
  5. Writes the post hoc result JSON (post_hoc: true, decision_date: 2026-09-29,
     authority: this rule section/D6b, plus every hash
     c1_ucec_completeness_rerun.py's result JSON already records) and, when
     qualifies, the M / M_op / M_op_resid matrices (.npz).

Nothing here can change the frozen section 5 reading of the original C1-UCEC test
("In every case, UCEC's published section 5 reading is unchanged" -- D4's closing
sentence). This is a descriptive, pre-declared-wording, post hoc addendum only.

Run (smp; env from lsf_c1_ucec_d4_strat.sh so no MORPHO_* leaks in -- mirrors
lsf_c1_ucec_strat.sh's env block exactly):
  scp upload list (batch labels + D4 summary from idc_ucec_derive_batch_labels.py,
  run locally beforehand, plus this script and its wrapper -- everything else is
  already on the cluster from the original C1-UCEC run). Two destinations, not one:
    scp server_export/scripts/pinned/c1_ucec_batch_labels.tsv \\
        review/C1_UCEC_D4_2026-09-29/c1_ucec_d4_summary.json \\
        user@cluster:/public/home/fjhui/ZW/scripts/pinned/
    scp server_export/scripts/c1_ucec_stratified_reread.py \\
        server_export/scripts/lsf_c1_ucec_d4_strat.sh \\
        user@cluster:/public/home/fjhui/ZW/scripts/
  Smoke (minutes; validates auc_obs only, refuses to print any stratified p):
    bsub -q smp -n 8 -R "span[hosts=1]" -o c1_ucec_d4_strat_smoke.out \\
      lsf_c1_ucec_d4_strat.sh --perms 5 --tag smoke
  Full run (B=1000, the real post hoc analysis):
    bsub -q smp -n 8 -R "span[hosts=1]" -o c1_ucec_d4_strat.out \\
      lsf_c1_ucec_d4_strat.sh --perms 1000
"""
import argparse
import hashlib
import json
import os
import sys
import time
import zlib
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import residual_analysis as RA  # noqa: E402
import split_replication_power as SP  # noqa: E402
import c1_run_test as C1  # noqa: E402  (prepare, assemble, MIN_N -- not reimplemented)

FLOAT_TOL = 1e-9

# sitepack literals (residual_analysis_sitepack.py L89/98), copied for the
# dynamic-selection extra variant (step 3b docstring above). MIN_PERMUTABLE is not
# duplicated here: the D1 permutable-share check already ran, on this same axis, in
# idc_ucec_derive_batch_labels.py's evaluate_d1_criteria -- this script only acts on
# its qualifies verdict (read from --d4-summary below) and does not re-gate on it.
MAX_LEVELS = 12
MIN_GROUP = 3

# The published c1_ucec_result.json (2026-09-16 run). Pinned literally, same value
# and same rationale as c1_ucec_completeness_rerun.py's own copy of this constant:
# so a later, possibly-different regeneration of that file cannot silently change
# what this re-read is checked against.
PUBLISHED_RESULT_SHA256 = "53fcb86138465e90ef2c47e3519786fc0fb4a534497f4a83850ee584ee6e653c"


def sha256(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return "unreadable"


def emp_p(obs, null):
    """One-sided upper-tail empirical p-value: (1 + #{null >= obs}) / (1 + B).
    Same convention as c1_run_test.py / split_replication_power.py /
    c1_ucec_completeness_rerun.py."""
    null = np.asarray(null, dtype=float)
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))


def validate(log, computed, published_path, tol=FLOAT_TOL):
    """FATAL exit (printing both sets of values) unless every field of `computed`
    matches the corresponding field of the published JSON at `published_path`,
    AND that file is byte-identical to the 2026-09-16 run this script reproduces
    (checked by sha256). Identical logic to
    c1_ucec_completeness_rerun.py's validate() -- duplicated, not imported, so
    this script's own hash covers its complete validation logic (this project's
    convention: c1_run_test_luad.py similarly duplicates rather than imports
    c1_run_test.py's assemble/sha256)."""
    if not os.path.exists(published_path):
        sys.exit(f"FATAL: --validate-against file does not exist: {published_path}")
    got_file_hash = sha256(published_path)
    if got_file_hash != PUBLISHED_RESULT_SHA256:
        sys.exit(f"FATAL: {published_path} has sha256={got_file_hash}, expected the "
                 f"pinned 2026-09-16 result {PUBLISHED_RESULT_SHA256}. Refusing to "
                 f"validate against a file that is not byte-identical to the "
                 f"published C1-UCEC result.")
    with open(published_path, "r", encoding="utf-8") as f:
        published = json.load(f)

    float_fields = ["auc_obs", "auc_null_mean", "auc_null_sd", "p_auc", "rho_obs", "p_rho"]
    count_fields = [("n_usable", "n_usable"), ("n_selected", "n_selected"),
                    ("n_background", "n_background"), ("n_cases", "n_cases"),
                    ("K", "K"), ("n_genes_tested", "n_genes_tested")]
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


def load_operator_labels(common, batch_labels_path, log):
    """case_submitter_id -> operator_uuid, "" for any case with no row / an empty
    label. Same shape as c1_run_test_luad.py's load_operator_labels, ported to
    UCEC's batch-labels file (idc_ucec_derive_batch_labels.py's output, identical
    schema to pinned/c1_luad_batch_labels.tsv)."""
    bl = pd.read_csv(batch_labels_path, sep="\t", dtype=str).fillna("")
    lab = bl.set_index("case_submitter_id")["operator_uuid"].to_dict()
    raw = np.array([lab.get(c, "") for c in common])
    n_unlabelled = int((raw == "").sum())
    log(f"operator labels: {len(common) - n_unlabelled}/{len(common)} cases labelled "
        f"({n_unlabelled} unlabelled, held fixed in the within-operator null) "
        f"<- {batch_labels_path}")
    return raw


def build_within_operator_perms(raw, B, seed_bytes, log, tag):
    """Departure (i) from sitepack (same departure as c1_run_test_luad.py's
    build_within_operator_perms, reimplemented here with UCEC's own seed tag):
    unlabelled cases ("") are NOT a level -- excluded from every stratum, never
    permuted. Singleton strata are likewise held fixed (a permutation within a
    stratum of size 1 is the identity). Strata are processed in sorted label
    order for a deterministic RNG draw sequence."""
    n = len(raw)
    n_labelled = int((raw != "").sum())
    levels = sorted(set(v for v in raw if v))
    members = {lv: np.where(raw == lv)[0] for lv in levels}
    n_strata_ge3 = sum(1 for m in members.values() if len(m) >= 3)
    largest = max((len(m) for m in members.values()), default=0)
    fixed = int((raw == "").sum()) + sum(len(m) for m in members.values() if len(m) <= 1)
    permutable = 1.0 - fixed / n if n else 0.0
    coverage = n_labelled / n if n else 0.0
    largest_frac = (largest / n_labelled) if n_labelled else 1.0
    log(f"  [{tag}] operator stratum diagnostics: coverage {coverage:.1%}, "
        f"strata>=3cases {n_strata_ge3}, largest {largest}/{n_labelled}="
        f"{largest_frac:.1%}, permutable {permutable:.1%}, {len(levels)} levels "
        f"(cross-checked against idc_ucec_derive_batch_labels.py's D1 report, "
        f"which is the one this script's qualification GATE actually uses)")
    rng = np.random.RandomState(zlib.crc32(seed_bytes))
    perms = []
    for _ in range(B):
        p = np.arange(n)
        for lv in levels:  # sorted already
            m = members[lv]
            if len(m) > 1:
                p[m] = rng.permutation(m)
        perms.append(p)
    diag = dict(coverage=coverage, n_strata_ge3=n_strata_ge3, largest_frac=largest_frac,
                permutable=permutable, n_levels=len(levels))
    return perms, diag


def build_operator_dummy_design_dynamic(raw, pcs, log):
    """EXTRA variant only (step 3b of the module docstring; NOT part of D4 item
    3's literal text). Reimplements residual_analysis_sitepack.py L279-322's
    ORIGINAL dynamic dummy-selection rule exactly:
      dummy_levels = [lv for lv, k in counts.most_common(MAX_LEVELS) if k >= MIN_GROUP]
      if len(dummy_levels) == len(levels) and dummy_levels: dummy_levels = dummy_levels[:-1]
    computed from Counter.most_common on the REALISED UCEC operator labels -- not
    a frozen table (none exists for UCEC; c1_run_test_luad.py's
    FROZEN_DUMMY_OPERATORS has no UCEC counterpart, because D4 was decided
    2026-09-29 for a re-read of an already-published result, not pre-registered
    ahead of any data).

    One departure from sitepack's literal main(), kept consistent with (i) above:
    unlabelled cases ("") are EXCLUDED from `counts`/`levels` entirely (never a
    candidate dummy level), rather than pooled into sitepack's default "_missing"
    level via `lab.get(c, "") or "_missing"`. This matches build_within_operator_perms
    above -- unlabelled cases get no dummy and are simply absorbed by the
    intercept, exactly as c1_run_test_luad.py's build_operator_dummy_design
    departure (i) documents for LUAD.
    """
    counts = Counter(v for v in raw if v)
    levels = sorted(counts)
    dummy_levels = [lv for lv, k in counts.most_common(MAX_LEVELS) if k >= MIN_GROUP]
    if len(dummy_levels) == len(levels) and dummy_levels:
        dummy_levels = dummy_levels[:-1]   # keep the design full-rank (sitepack L305-306)
    n = len(raw)
    D = np.zeros((n, len(dummy_levels)))
    for j, lv in enumerate(dummy_levels):
        D[:, j] = (raw == lv).astype(float)
    Bdes = np.hstack([np.ones((n, 1)), D])
    coef, *_ = np.linalg.lstsq(Bdes, pcs, rcond=None)
    pcs_resid = pcs - Bdes @ coef
    kept = pcs_resid.var(axis=0).sum() / max(pcs.var(axis=0).sum(), 1e-12)
    log(f"  [extra_residualized] dynamic dummy selection (sitepack original rule, "
        f"MAX_LEVELS={MAX_LEVELS}, MIN_GROUP={MIN_GROUP}): {len(dummy_levels)} of "
        f"{len(levels)} realised operator levels kept as dummies "
        f"{dummy_levels}, PC variance retained {kept:.1%} (n={n})")
    return pcs_resid, dummy_levels


def main():
    ap = argparse.ArgumentParser()
    # Identical defaults to c1_run_test.py's ArgumentParser (mode "test" branch)
    # and to c1_ucec_completeness_rerun.py's own copy of them -- copied, not
    # derived, per that script's own stated rationale.
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=RA.N_WORKERS)
    ap.add_argument("--tag", default="")
    ap.add_argument("--protein", default="/public/home/fjhui/ZW/ucec_c1/protein_ucec_c1.csv")
    ap.add_argument("--discovery-results",
                    default="/public/home/fjhui/ZW/scripts/pinned/ucec/residual_results_tumoronly.csv")
    ap.add_argument("--discovery-slide-map",
                    default="/public/home/fjhui/ZW/scripts/pinned/slide_type_map_ucec.tsv")
    # New in this script.
    ap.add_argument("--validate-against", default="",
                     help="published c1_ucec_result.json to validate the "
                          "reproduced statistics against. Default: "
                          "<MORPHO_OUT>/c1_ucec_result.json, matching "
                          "c1_ucec_completeness_rerun.py's own default (same env, "
                          "same original-run output directory).")
    ap.add_argument("--no-validate", action="store_true",
                     help="skip the validate-against-published-JSON step. ONLY for "
                          "a synthetic dry run, where monkeypatched loaders make "
                          "the real published numbers unreachable by construction. "
                          "Never pass this against real data.")
    ap.add_argument("--batch-labels",
                    default="/public/home/fjhui/ZW/scripts/pinned/c1_ucec_batch_labels.tsv",
                    help="idc_ucec_derive_batch_labels.py's output (uploaded from "
                         "server_export/scripts/pinned/c1_ucec_batch_labels.tsv)")
    ap.add_argument("--d4-summary",
                    default="/public/home/fjhui/ZW/scripts/pinned/c1_ucec_d4_summary.json",
                    help="idc_ucec_derive_batch_labels.py's qualification summary "
                         "(uploaded from "
                         "review/C1_UCEC_D4_2026-09-29/c1_ucec_d4_summary.json). "
                         "This script trusts its \"qualifies\" field; it does not "
                         "re-derive it.")
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
    log(f"  batch labels       : {args.batch_labels}")
    log(f"  D4 summary         : {args.d4_summary}")
    log(f"  perms / workers    : {args.perms} / {args.workers}")
    log("=" * 78)

    t0 = time.time()

    # ---- reproduce c1_run_test.run_test() exactly up to M ---------------------
    common, tested, selected, disc_incr, disc_fdr, R, P, W, K = C1.prepare(args, log)
    n = len(common)

    pcs = PCA(n_components=K, svd_solver="full").fit_transform(W)
    rng = np.random.RandomState(zlib.crc32(b"c1_ucec|perm"))
    perms = [rng.permutation(n) for _ in range(args.perms)]
    log(f"per-gene pass: {len(tested)} genes x (1 observed + {args.perms} permutations) "
        f"on {args.workers} workers (unrestricted null, reproducing the published run) ...")
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
    # ---- end exact reproduction -----------------------------------------------

    auc, rho, n1, n0 = C1.assemble(M, xs, sel_mask)
    auc_obs, auc_null = auc[0], auc[1:]
    rho_obs, rho_null = rho[0], rho[1:]
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
            "any published result. Only ever pass this for a synthetic dry run.")
    elif smoke_short_run:
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
            f"{p_floor:.4f} -- refusing to compute any stratified null for a smoke run.")
    else:
        validate(log, computed, validate_against)

    # ---- D4 item 3: read the qualification verdict -----------------------------
    if not os.path.exists(args.d4_summary):
        sys.exit(f"FATAL: --d4-summary file does not exist: {args.d4_summary}. Run "
                 f"idc_ucec_derive_batch_labels.py first (locally; it needs no HPC "
                 f"data) and upload its output, or point --d4-summary at a summary "
                 f"file for a synthetic dry run.")
    with open(args.d4_summary, "r", encoding="utf-8") as f:
        d4_summary = json.load(f)
    qualifies = bool(d4_summary.get("qualifies"))
    validated_operator = bool(d4_summary.get("validated_operator"))
    log(f"\nD4 item 3 qualification (from {args.d4_summary}): "
        f"validated_operator={validated_operator} qualifies={qualifies}")

    result = dict(
        post_hoc=True,
        decision_date="2026-09-29",
        authority="review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 0, D4 "
                  "item 3 (D6b)",
        cohort="ucec_c1", n_cases=n, K=K, perms=args.perms,
        n_usable=computed["n_usable"], n_selected=computed["n_selected"],
        n_background=computed["n_background"],
        pooled_auc_obs=computed["auc_obs"], pooled_p_auc=computed["p_auc"],
        d4_summary_path=args.d4_summary,
        validated_operator=validated_operator, qualifies=qualifies,
    )

    FALLBACK_SENTENCE = ("Methods state that C1-UCEC was tested under the "
                          "unrestricted null only and carries no batch qualifier.")

    if not qualifies:
        log("\n" + "=" * 78)
        log("D4 item 3 fallback (operator axis did not validate, or a D1 "
            "criterion failed, on the published UCEC analysed set):")
        log(f"  \"{FALLBACK_SENTENCE}\"")
        log("=" * 78)
        log("Not computing a stratified AUC (D4 item 3's fallback clause). This is "
            "an anticipated, pre-declared outcome, not an error.")
        result["p_auc_strat"] = None
        result["fallback_wording"] = FALLBACK_SENTENCE
        result["validated"] = (False if args.no_validate else True)
        result["minutes"] = round((time.time() - t0) / 60.0, 1)
        result["sha256_this_script"] = this_script_hash
        result["sha256_discovery_results"] = computed["sha256_discovery_results"]
        result["sha256_residual_analysis"] = computed["sha256_residual_analysis"]
        result["sha256_split_replication_power"] = computed["sha256_split_replication_power"]
        result["sha256_c1_run_test"] = sha256(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "c1_run_test.py"))
        result["sha256_d4_summary"] = sha256(args.d4_summary)
        suffix = f"_{args.tag}" if args.tag else ""
        res_fp = RA.OUT_DIR / f"c1_ucec_d4_stratified_result{suffix}.json"
        with open(res_fp, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        log(f"\n  -> {res_fp}\n  elapsed {result['minutes']} min")
        return

    if args.no_validate or smoke_short_run:
        log("\nNOT computing a stratified AUC: "
            + ("--no-validate run" if args.no_validate else
               f"B={args.perms} smoke run (see refusal above)."))
        result["p_auc_strat"] = None
        result["validated"] = (False if args.no_validate else "smoke_auc_obs_only")
        result["minutes"] = round((time.time() - t0) / 60.0, 1)
        result["sha256_this_script"] = this_script_hash
        result["sha256_discovery_results"] = computed["sha256_discovery_results"]
        result["sha256_residual_analysis"] = computed["sha256_residual_analysis"]
        result["sha256_split_replication_power"] = computed["sha256_split_replication_power"]
        result["sha256_c1_run_test"] = sha256(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "c1_run_test.py"))
        result["sha256_d4_summary"] = sha256(args.d4_summary)
        suffix = f"_{args.tag}" if args.tag else ""
        res_fp = RA.OUT_DIR / f"c1_ucec_d4_stratified_result{suffix}.json"
        with open(res_fp, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        log(f"\n  -> {res_fp}\n  elapsed {result['minutes']} min")
        return

    # ---- qualifies AND B=1000 (or a real, non-smoke B): compute the re-read ----
    # Fail-closed provenance: the labels file this run is about to use must be the
    # SAME file idc_ucec_derive_batch_labels.py actually produced and recorded in
    # this same --d4-summary JSON, not a stale or hand-edited copy at the same path.
    expected_labels_hash = d4_summary.get("sha256_output_labels")
    if not expected_labels_hash:
        sys.exit(f"FATAL: {args.d4_summary} has no sha256_output_labels field (an "
                 f"older idc_ucec_derive_batch_labels.py run, before that field "
                 f"existed) -- regenerate it before trusting --batch-labels.")
    got_labels_hash = sha256(args.batch_labels)
    if got_labels_hash != expected_labels_hash:
        sys.exit(f"FATAL: --batch-labels sha256 mismatch: got={got_labels_hash} "
                 f"want={expected_labels_hash} (from {args.d4_summary}). The labels "
                 f"file does not match the one the qualification verdict was "
                 f"computed from -- refusing to residualise against it.")
    log(f"  hash OK: --batch-labels matches sha256_output_labels in "
        f"{args.d4_summary} ({got_labels_hash[:12]}...)")
    result["sha256_batch_labels"] = got_labels_hash

    raw_op = load_operator_labels(common, args.batch_labels, log)

    log("\n=== D4 item 3 step 3 (PRIMARY, rule-literal): same 2,566-set PCs, "
        "within-operator null ===")
    perms_op, diag_op = build_within_operator_perms(raw_op, args.perms, b"c1_ucec|perm_op",
                                                     log, tag="primary")
    res_op = SP.run_pool(SP._gene_test, len(tested), args.workers, R, P, pcs, perms_op, C1.MIN_N)
    M_op = np.vstack([res_op[j] for j in usable])
    auc_op, rho_op, *_ = C1.assemble(M_op, xs, sel_mask)
    auc_obs_op = float(auc_op[0])
    if abs(auc_obs_op - computed["auc_obs"]) > 1e-9:
        sys.exit(f"FATAL: within-operator-null observed AUC {auc_obs_op!r} != "
                 f"unrestricted-null observed AUC {computed['auc_obs']!r} -- these "
                 f"must be bit-identical (same PCs, same gene sets; only the null "
                 f"differs). This is a bug.")
    p_auc_strat = emp_p(auc_obs_op, auc_op[1:])
    log(f"  AUC_strat observed (== pooled AUC, by construction): {auc_obs_op:.6f}")
    log(f"  AUC_strat null (B={args.perms})".ljust(29)
        + f": {auc_op[1:].mean():.6f} +/- {auc_op[1:].std():.6f}")
    log(f"  p (AUC_strat, one-sided)                : {p_auc_strat:.6f}")
    if computed["p_auc"] < 0.05 and p_auc_strat >= 0.05:
        log("  Note: replicated under the unrestricted null but not under the "
            "operator-stratified null (LUAD clause [B], form 2 wording pattern).")

    result["auc_strat_obs"] = auc_obs_op
    result["auc_strat_null_mean"] = float(auc_op[1:].mean())
    result["auc_strat_null_sd"] = float(auc_op[1:].std())
    result["p_auc_strat"] = p_auc_strat
    result["operator_stratum_diagnostics"] = diag_op
    result["validated"] = True

    log("\n=== EXTRA (descriptive, NOT required by D4 item 3's literal text; see "
        "module docstring step 3b): residualised PCs, dynamic sitepack-original "
        "dummy selection, same within-operator null ===")
    pcs_resid, dummy_levels = build_operator_dummy_design_dynamic(raw_op, pcs, log)
    res_resid = SP.run_pool(SP._gene_test, len(tested), args.workers, R, P, pcs_resid,
                            perms_op, C1.MIN_N)
    M_resid = np.vstack([res_resid[j] for j in usable])
    auc_resid, rho_resid, *_ = C1.assemble(M_resid, xs, sel_mask)
    auc_obs_resid = float(auc_resid[0])
    p_auc_strat_resid = emp_p(auc_obs_resid, auc_resid[1:])
    log(f"  AUC observed (residualised PCs; NOT expected to equal the pooled AUC): "
        f"{auc_obs_resid:.6f}")
    log(f"  null (B={args.perms})".ljust(29)
        + f": {auc_resid[1:].mean():.6f} +/- {auc_resid[1:].std():.6f}")
    log(f"  p (one-sided)            : {p_auc_strat_resid:.6f}")

    result["extra_residualized_variant"] = dict(
        descriptive_not_required_by_d4_item3_text=True,
        rationale="Requested by the task that commissioned this script; the rule's "
                  "own D4 item 3 text does not residualise any PC for this item "
                  "(see module docstring). Does not gate or replace the primary "
                  "re-read above.",
        dummy_selection_rule="residual_analysis_sitepack.py original dynamic "
                             f"selection: Counter.most_common(MAX_LEVELS={MAX_LEVELS}), "
                             f"kept if size>=MIN_GROUP={MIN_GROUP}, full-rank guard.",
        dummy_levels=dummy_levels, n_dummy_levels=len(dummy_levels),
        auc_obs=auc_obs_resid,
        auc_null_mean=float(auc_resid[1:].mean()), auc_null_sd=float(auc_resid[1:].std()),
        p_auc=p_auc_strat_resid,
    )

    minutes = round((time.time() - t0) / 60.0, 1)
    result["minutes"] = minutes
    result["sha256_this_script"] = this_script_hash
    result["sha256_discovery_results"] = computed["sha256_discovery_results"]
    result["sha256_residual_analysis"] = computed["sha256_residual_analysis"]
    result["sha256_split_replication_power"] = computed["sha256_split_replication_power"]
    result["sha256_c1_run_test"] = sha256(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "c1_run_test.py"))
    result["sha256_batch_labels"] = sha256(args.batch_labels)
    result["sha256_d4_summary"] = sha256(args.d4_summary)

    suffix = f"_{args.tag}" if args.tag else ""
    npz_fp = RA.OUT_DIR / f"c1_ucec_d4_stratified{suffix}.npz"
    np.savez_compressed(npz_fp, genes=np.array(genes_u), sel_mask=sel_mask,
                        M_unrestricted=M, M_within_operator=M_op,
                        M_within_operator_residualized=M_resid,
                        operator_raw=raw_op)
    result["sha256_npz"] = sha256(str(npz_fp))

    res_fp = RA.OUT_DIR / f"c1_ucec_d4_stratified_result{suffix}.json"
    with open(res_fp, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    log(f"\n  -> {res_fp}\n  -> {npz_fp}\n  elapsed {minutes} min")


if __name__ == "__main__":
    main()
