#!/usr/bin/env python3
"""
c1_run_test_luad.py -- the C1-LUAD confirmatory test, per
review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md (the signed rule; section references
below are to that file: D1/D4/D5/D6 are in section 0). This is NOT a pre-check; it has
no tunable reading once B=1000 is used.

Structure, and how it differs from c1_run_test.py (the executed C1-UCEC test, which
this script mirrors in overall shape):
  primary set = FIXED, already computed on the discovery cohort (n=105). H1 is the
                uncorrected fdr<0.05 & incremental_r2>0 set: 2,310 of 10,724 tested
                genes (luad__residual_results_tumoronly.csv). H2 is the batch-corrected
                fdr_batchresid<0.05 & incremental_r2_batchresid>0 set: 726 genes
                (luad__sitepack_operator.csv). Both counts are asserted below.
  family     = a THREE-STEP fixed sequence at one-sided alpha 0.05, stopping at the
               first non-rejection (D1, D6a): H1 -> (a)(i) batch qualifier -> H2.
               (a)(ii) and the completeness-stratified AUC are co-reported, descriptive,
               and cannot move any reading.
  capacity   = K = round(20 * n / 100), PCA(K) fit on the CONFIRMATORY cohort's own
               embeddings (§3).
  null       = TWO null families, both B=1,000, both drawing gene-gene-correlation-
               preserving joint permutations (every gene recomputed together per draw):
                 - unrestricted (H1): RandomState(crc32(b"c1_luad|perm"))
                 - within-operator (a(i)/H2/(a)(ii)): RandomState(crc32(b"c1_luad|perm_op")),
                   built from the FROZEN operator labels in --batch-labels. Unlabelled
                   cases are held fixed (never pooled into a permutable level) -- this
                   is departure (i) from residual_analysis_sitepack.py's main(), whose
                   design cannot otherwise be imported (its dummy/residualisation/
                   permutation logic is inline in main(), bound to discovery-only data
                   paths via batch_leak_check.cohort_paths/find_proteome). Departure
                   (ii): the 12 residualisation dummy operators are FROZEN (below), not
                   recomputed at runtime from Counter.most_common.
  residual.  = confirmatory K PCs residualised on the frozen 12-operator dummy design
               by least squares with an intercept (sitepack L307-322, reimplemented,
               not imported -- see departures above). If a frozen operator has fewer
               than MIN_GROUP=3 realised cases, it loses its column; no replacement is
               added.

Modes (mirrors c1_run_test.py's test/bootstrap split, widened because the LUAD
confirmatory family is larger; running each as its own LSF job keeps every job inside
a sane walltime and lets a smoke run cover only the piece being changed):
  primary       H1, the (a)(i) batch qualifier, H2, (a)(ii), and the completeness-
                stratified AUC (D1, section 3, section 4, section 5). Writes the M
                matrix (genes x (1+B)) and per-gene confirmatory n as .npz.
  sensitivity   (s) and (s2), H1 only, on case subsets (section 3).
  bootstrap     Patient-level bootstrap intervals: H1 always; H2 only with
                --h2-claimed (pass this only after a `primary` run gave p1<0.05 AND
                p_strat<0.05) (section 5.4).
  descriptive   Prints the invocation for c1_luad_pergene.py (b1)/(b2) -- that is a
                separate script/module per the rule's Appendix item 3, not part of
                this file, because it does not reuse split_replication_power.py's
                per-gene machinery.

Every path is a REQUIRED argument; none has a default. That is deliberate: a default
is exactly how a stale environment silently redirects one cohort's run into another
cohort's output (this project already hit that once, UCEC into pdac/results,
2026-09-15). --rule and --rule-sha256 pin this run to one exact, signed copy of the
frozen rule; every other --*-sha256 argument pins one input file the same way. All are
FATAL-asserted before anything is computed. The result JSON records the sha256 of the
rule file, of this script itself, and of every imported module, next to the input
hashes (section 7.2).

--blind (section 5.6, D5g) stops immediately after `prepare()` computes and logs n, K
and the tested/selected gene counts. Nothing past that point runs: no AUC, rho, p,
gene-level table, matrix or result JSON is computed, written or printed. Use this for
every confirmatory smoke test; never smoke-test on confirmatory data without it.

Run (smp; env from lsf_c1_luad_test.sh so no MORPHO_* leaks in):
  bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_smoke.out \
    lsf_c1_luad_test.sh primary --perms 5 --tag smoke --blind
  bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_primary.out \
    lsf_c1_luad_test.sh primary --perms 1000
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

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
# split_replication_power.py's canonical HPC location is beside this script
# (${ZW}/scripts/, per section 7.2: "HPC scripts/"), matching c1_run_test.py's own
# convention and its documented sha256 lookup. Locally it sits at the repository top
# level instead (one directory up from server_export/scripts/), so the repo root is
# also added to sys.path -- import succeeds in both places, but the hash recorded in
# every result JSON always looks beside this script, matching what ships to the HPC;
# a local run legitimately reads "unreadable" for that lookup, exactly as
# c1_run_test.py already does (documented in c1_ucec_completeness_rerun.py's own
# dry-run notes).
sys.path.insert(1, os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
import residual_analysis as RA  # noqa: E402
import split_replication_power as SP  # noqa: E402  (k_for, _gene_test, run_pool)

# ---------------------------------------------------------------------------
# Frozen constants (section 2, section 3, section 4). Nothing here is derived
# from data seen after 2026-09-29; every count is asserted against these.
# ---------------------------------------------------------------------------
N_TESTED_EXPECTED = 10724     # discovery-tested genes, luad__residual_results_tumoronly.csv
N_SELECTED_EXPECTED = 2310    # H1: fdr<0.05 & incremental_r2>0
N_SELECTED_H2_EXPECTED = 726  # H2: fdr_batchresid<0.05 & incremental_r2_batchresid>0
MIN_USABLE_CASES = 45         # section 5.1/5.5: below this the test is abandoned
PCS_PER_100 = 20.0            # section 3
MIN_N = 30                    # same floor as residual_analysis._process_gene
N_CASE_IDS_EXPECTED = 113
MAX_LEVELS = 12               # sitepack literal (section 4)
MIN_GROUP = 3                 # sitepack literal (section 4)
MIN_PERMUTABLE = 0.50         # sitepack literal (section 4)
DISCOVERY_MODAL_N_EXPECTED = 105  # LUAD discovery cohort size (completeness strata)

# Frozen residualisation dummy set (section 4): the 12 operators ranked first by case
# count on the 112 expected analysed cases, ties broken by first appearance in sorted
# case order. 925cc759-f2ba-47ea-b2dc-fb1e39fa77da (3 cases) is deliberately NOT in
# this list -- it is left to the intercept, not a 13th dummy. If, at runtime, a listed
# operator has fewer than MIN_GROUP realised cases, it loses its column; no other
# operator is substituted in its place.
FROZEN_DUMMY_OPERATORS = [
    "9c62eb1b-888b-4319-a031-5ee965150a3e",   # 16 cases (of 112 expected)
    "888d279c-9467-4931-8ead-381499dea990",   # 11
    "975375e8-786f-45d0-bfdd-9e70146a3eb8",   # 10
    "fc28642d-0894-46ef-9386-1ba3b5bdd99c",   # 8
    "6ed758c2-29d4-4ccb-945c-a7454ee330f5",   # 6
    "f3d46998-d392-4ce1-bd72-cf82622a2333",   # 5
    "00000000-0000-0000-0000-000000000000",   # 5 (all-zero UUID; D5a: real stratum, literal)
    "c8152987-eedd-4cf7-9b33-fee82a532dbc",   # 4
    "56d2978b-1a1d-4c5a-a722-3daa9816a1f5",   # 4
    "c1a7d89d-1759-40c1-af7a-8ba9c7901cec",   # 4
    "052aeb46-e02a-4dae-ae21-c92500e0b9d4",   # 3
    "b698ee8d-5921-4ead-bf67-ee1d6ed62b52",   # 3
]

# The 5 non-standard-grid cases removed by sensitivity (s): 2 at 40x/0.25um (Leica),
# 3 at 0.32um with no operator (PixelMed). Section 1.3/3.
SENSITIVITY_S_DROP_CASES = frozenset({
    "C3L-02513", "C3L-02515", "C3L-03642", "C3L-03717", "C3L-03721",
})
# The 12 cases with no SS1553 slide (9 on a device absent from LUAD discovery, 3 with
# no recorded scanner). Section 1.3/3/6.
SENSITIVITY_S2_DROP_CASES = frozenset({
    "C3L-02513", "C3L-02515", "C3L-03642",                      # no scanner recorded
    "C3L-03717", "C3L-03721",                                   # SS7320, 40x
    "C3L-03985", "C3L-04757", "C3N-03765", "C3N-04157",
    "C3N-04168", "C3N-04176", "C3N-04180",                      # SS7559
})

# section 7.2: frozen sha256 of every local-mirror input this script reads. The HPC
# copy at each --*-path argument is asserted identical to these before anything else
# runs. A mismatch is FATAL -- it means the wrong file vintage, not a difference this
# script may silently absorb.
EXPECTED_SHA256 = {
    "case_ids": "3629b243a3aa21057dae24e7f5c1b47e8a2571e6ffbcff6f3170333e1ca9fe44",
    "discovery_results": "63aec1992a77e952b96c4e25e27372e2bee5c52294560c56ef2b229196e87be7",
    "discovery_sitepack": "aeef7232682c24d209024ba73efea8033cb7486779831bda1a1ba3c8bda45d9f",
    "batch_labels": "0a840533a3601490c360cae8e5fcfc5e6f426e5314e5ac02d1e1b00730f1b3d4",
    "discovery_slide_map": "3de16eb77f615908c1999806388d9871b48ed8bfc645a6366a7f64b7df9a843e",
    "discovery_xwalk": "967cabf18cb6b6c3646aaf05005c46fcdccd7b55eb8dc9882348e84ca4baa2bf",
    "residual_analysis": "a54a0a56494d315e15ef075425fe2133ecd9f05193c22656c1145bffdaac4a0a",
    "split_replication_power": "79530c075c0f7e75bf06d3c2f1edfe8792cbf4d6c71ae1d32e60bf44f31a0397",
}


def sha256(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return "unreadable"


def assert_hash(path, key, log, allow_skip_env=None):
    """FATAL unless sha256(path) == EXPECTED_SHA256[key]. allow_skip_env names an env
    var that, if set to a nonempty value, downgrades the mismatch to a loud warning --
    used only for discovery_xwalk/discovery_slide_map on the HPC where the pinned copy
    may legitimately carry an rna-manifest suffix the local mirror does not; unset by
    default, so the check is fail-closed unless someone deliberately opts out."""
    got = sha256(path)
    want = EXPECTED_SHA256[key]
    if got == want:
        log(f"  hash OK: {key} = {got[:12]}... ({path})")
        return got
    msg = (f"sha256 mismatch for {key}: got={got} want={want} ({path}). "
           f"This is the wrong file vintage -- stopping rather than testing against "
           f"an unverified input.")
    if allow_skip_env and os.environ.get(allow_skip_env):
        log(f"  WARNING (skip forced via {allow_skip_env}): {msg}")
        return got
    sys.exit(f"FATAL: {msg}")


class CheckFatal(Exception):
    pass


def load_case_ids(path, log):
    assert_hash(path, "case_ids", log)
    with open(path, encoding="utf-8") as f:
        ids = sorted({line.strip() for line in f if line.strip()})
    if len(ids) != N_CASE_IDS_EXPECTED:
        sys.exit(f"FATAL: {path} has {len(ids)} case ids, expected {N_CASE_IDS_EXPECTED}.")
    return set(ids)


def load_discovery_case_union(slide_map_path, xwalk_path, rna_manifest_path, log):
    """Case-level union of the three HPC discovery pinned files (section 7.2). Each
    must carry exactly 111 discovery cases; the union is what the independence check
    below must not overlap with the frozen 113 confirmatory cases."""
    assert_hash(slide_map_path, "discovery_slide_map", log, allow_skip_env="MORPHO_ALLOW_HPC_SLIDE_MAP_HASH")
    assert_hash(xwalk_path, "discovery_xwalk", log, allow_skip_env="MORPHO_ALLOW_HPC_XWALK_HASH")
    sm = pd.read_csv(slide_map_path, sep="\t", dtype=str)
    xw = pd.read_csv(xwalk_path, sep="\t", dtype=str)
    rm = pd.read_csv(rna_manifest_path, sep="\t", dtype=str)
    cases_sm = set(sm["case_submitter_id"]) if "case_submitter_id" in sm.columns else set()
    xwalk_case_col = "case_id" if "case_id" in xw.columns else "case_submitter_id"
    cases_xw = set(xw[xwalk_case_col]) if xwalk_case_col in xw.columns else set()
    cases_rm = set(rm["case_submitter_id"]) if "case_submitter_id" in rm.columns else set()
    for name, s in [("slide map", cases_sm), ("aliquot crosswalk", cases_xw), ("rna manifest", cases_rm)]:
        if len(s) != 111:
            log(f"  WARNING: discovery {name} has {len(s)} distinct cases, expected 111 "
                f"({slide_map_path if name == 'slide map' else xwalk_path if name == 'aliquot crosswalk' else rna_manifest_path})")
    union = cases_sm | cases_xw | cases_rm
    log(f"  discovery case union (slide map + xwalk + rna manifest): {len(union)} cases")
    return union


def independence_check(common, frozen_113, discovery_union, log, selftest=False):
    """Fail-closed. `selftest=True` re-adds 7 known bridge cases in memory only and
    must raise CheckFatal; used by --selftest, never on the real confirmatory case
    set."""
    cases = set(common)
    if selftest:
        cases = cases | {"C3L-02348", "C3L-02350", "C3N-00545", "C3N-00738",
                         "C3N-01023", "C3N-01024", "C3N-02003"}
    off_frozen = sorted(cases - frozen_113)
    if off_frozen:
        raise CheckFatal(f"{len(off_frozen)} analysed cases are not on the frozen "
                         f"113-case list: {off_frozen[:5]}")
    overlap = sorted(cases & discovery_union)
    if overlap:
        raise CheckFatal(f"{len(overlap)} analysed cases overlap the discovery case "
                         f"union: {overlap[:5]}")
    log(f"  independence check{' (selftest, +7 bridge cases)' if selftest else ''}: "
        f"0 overlap with discovery, 0 off the frozen list ({len(cases)} cases checked)")


def load_confirmatory_protein(path, log):
    df = pd.read_csv(path, index_col=0).apply(pd.to_numeric, errors="coerce")
    log(f"confirmatory protein (case-indexed, from c1_pull_omics_luad.py): "
        f"{df.shape[0]} cases x {df.shape[1]} genes  <- {path}")
    return df


def load_discovery_n(discovery_results_path):
    d = pd.read_csv(discovery_results_path)
    return d.set_index("gene")["n"]


def build_gene_sets(discovery_results_path, discovery_sitepack_path, confirm_genes, log):
    d = pd.read_csv(discovery_results_path)
    if not d.gene.is_unique:
        sys.exit(f"FATAL: {discovery_results_path} has duplicate gene rows -- the "
                 f"selection set would be ambiguous.")
    if len(d) != N_TESTED_EXPECTED:
        sys.exit(f"FATAL: discovery result file has {len(d)} tested genes, expected "
                 f"{N_TESTED_EXPECTED}. Wrong file vintage -- stopping.")
    sel_h1_all = d[(d.fdr < 0.05) & (d.incremental_r2 > 0)]
    if len(sel_h1_all) != N_SELECTED_EXPECTED:
        sys.exit(f"FATAL: H1 selection set is {len(sel_h1_all)} genes, expected "
                 f"{N_SELECTED_EXPECTED}. Stopping.")

    s = pd.read_csv(discovery_sitepack_path).drop_duplicates("gene")
    if not {"fdr_batchresid", "incremental_r2_batchresid"}.issubset(s.columns):
        sys.exit(f"FATAL: {discovery_sitepack_path} is missing fdr_batchresid / "
                 f"incremental_r2_batchresid.")
    sel_h2_all = s[(s.fdr_batchresid < 0.05) & (s.incremental_r2_batchresid > 0)]
    if len(sel_h2_all) != N_SELECTED_H2_EXPECTED:
        sys.exit(f"FATAL: H2 selection set is {len(sel_h2_all)} genes (from "
                 f"fdr_batchresid/incremental_r2_batchresid), expected "
                 f"{N_SELECTED_H2_EXPECTED}. Do not silently fall back to "
                 f"fdr_both/incremental_r2_both, a different, non-frozen estimator "
                 f"(that pair gives 802, not 726). Stopping.")

    disc_incr = d.set_index("gene")["incremental_r2"]
    disc_fdr = d.set_index("gene")["fdr"]
    disc_incr_bres = s.set_index("gene")["incremental_r2_batchresid"]

    tested = [g for g in d.gene if g in confirm_genes]
    selected_h1 = set(sel_h1_all.gene) & set(tested)
    selected_h2 = set(sel_h2_all.gene) & set(tested)
    log(f"discovery H1 selection: {len(sel_h1_all)}/{len(d)} genes; usable universe "
        f"(discovery-tested AND present in confirmatory RNA+protein): {len(tested)} "
        f"genes, of which {len(selected_h1)} H1-selected / "
        f"{len(tested) - len(selected_h1)} H1-background")
    log(f"discovery H2 (batch-corrected) selection: {len(sel_h2_all)} genes; "
        f"{len(selected_h2)} present in the usable universe, "
        f"{len(set(sel_h1_all.gene) & set(sel_h2_all.gene))} overlap with H1's "
        f"discovery selection")
    return tested, selected_h1, selected_h2, disc_incr, disc_fdr, disc_incr_bres


def prepare(args, log, drop_cases=frozenset()):
    """Shared load + intersect + independence-check path for every mode. drop_cases,
    if given, removes those cases from `common` (used by (s)/(s2) sensitivity) --
    applied as a drop set, not a probe-then-filter, so sensitivity modes only load the
    data once."""
    rna = RA.load_rna_matrix()
    wsi = RA.load_wsi_embeddings()
    protein = load_confirmatory_protein(args.protein, log)

    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    log(f"cases with all three modalities: {len(common)}")

    frozen_113 = load_case_ids(args.case_ids, log)
    discovery_union = load_discovery_case_union(
        args.discovery_slide_map, args.discovery_xwalk, args.discovery_rna_manifest, log)

    # Self-test FIRST, on the untouched `common` set, before it is trusted for anything:
    # re-adding the 7 known bridge cases in memory must raise CheckFatal. If it does
    # not, the overlap logic itself is broken and nothing downstream can be trusted.
    try:
        independence_check(common, frozen_113, discovery_union, log, selftest=True)
    except CheckFatal:
        pass
    else:
        sys.exit("FATAL: independence-check selftest did not raise CheckFatal -- the "
                 "overlap logic is broken. Refusing to proceed on unverified guard code.")

    # Real check, on the real analysed set. This must NOT raise.
    try:
        independence_check(common, frozen_113, discovery_union, log)
    except CheckFatal as e:
        sys.exit(f"FATAL: independence check failed on the real analysed set: {e}")

    if drop_cases:
        before = len(common)
        present_drop = set(drop_cases) & set(common)
        common = sorted(set(common) - present_drop)
        log(f"  sensitivity drop applied: {before} -> {len(common)} cases "
            f"({len(present_drop)} of {len(drop_cases)} listed cases were present)")

    if len(common) < MIN_USABLE_CASES:
        sys.exit(f"ABANDONED per section 5.1/5.5: only {len(common)} usable cases "
                 f"(< {MIN_USABLE_CASES}). Report the coverage table; do not run the "
                 f"test.")

    genes_present = set(rna.columns) & set(protein.columns)
    tested, sel_h1, sel_h2, disc_incr, disc_fdr, disc_incr_bres = build_gene_sets(
        args.discovery_results, args.discovery_sitepack, genes_present, log)

    R = rna.loc[common, tested].values.astype(float)
    P = protein.loc[common, tested].values.astype(float)
    W = wsi.loc[common].values.astype(float)
    K = SP.k_for(len(common), PCS_PER_100)
    log(f"capacity: K = round({PCS_PER_100:g} * {len(common)} / 100) = {K} PCs, "
        f"PCA fit on the confirmatory cohort's own embeddings ({W.shape[1]} dims)")
    return common, tested, sel_h1, sel_h2, disc_incr, disc_fdr, disc_incr_bres, R, P, W, K


def load_operator_labels(common, batch_labels_path, log):
    assert_hash(batch_labels_path, "batch_labels", log)
    bl = pd.read_csv(batch_labels_path, sep="\t", dtype=str).fillna("")
    lab = bl.set_index("case_submitter_id")["operator_uuid"].to_dict()
    raw = np.array([lab.get(c, "") for c in common])
    n_unlabelled = int((raw == "").sum())
    log(f"operator labels: {len(common) - n_unlabelled}/{len(common)} cases labelled "
        f"({n_unlabelled} unlabelled, held fixed in the within-operator null)")
    return raw


def build_within_operator_perms(raw, B, log):
    """Departure (i) from sitepack: unlabelled cases ("") are NOT a level; they are
    excluded from every stratum and never permuted. Strata are processed in sorted
    label order for a deterministic RNG draw sequence."""
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
    qualifies = (coverage >= 0.80 and n_strata_ge3 >= 2
                and largest_frac <= 0.90 and permutable >= MIN_PERMUTABLE)
    log(f"operator stratum diagnostics (D1 criteria): coverage {coverage:.1%} "
        f"(need >=80%), strata>=3cases {n_strata_ge3} (need >=2), largest "
        f"{largest}/{n_labelled}={largest_frac:.1%} (need <=90%), permutable "
        f"{permutable:.1%} (need >=50%) -> {'QUALIFIES' if qualifies else 'NOT TESTABLE'}")
    rng = np.random.RandomState(zlib.crc32(b"c1_luad|perm_op"))
    perms = []
    for _ in range(B):
        p = np.arange(n)
        for lv in levels:  # sorted already
            m = members[lv]
            if len(m) > 1:
                p[m] = rng.permutation(m)
        perms.append(p)
    return perms, qualifies, permutable, levels, members


def build_operator_dummy_design(raw, pcs, log):
    """Reimplements residual_analysis_sitepack.py L307-322 exactly (design matrix,
    intercept, least-squares residualisation), with departure (ii): dummy_levels is
    the FROZEN list filtered to operators with >= MIN_GROUP realised cases, not
    recomputed from Counter.most_common(MAX_LEVELS)."""
    dummy_levels = [u for u in FROZEN_DUMMY_OPERATORS if int((raw == u).sum()) >= MIN_GROUP]
    dropped = [u for u in FROZEN_DUMMY_OPERATORS if u not in dummy_levels]
    if dropped:
        log(f"  frozen dummy operators dropped (below {MIN_GROUP} realised cases): "
            f"{dropped}")
    n = len(raw)
    D = np.zeros((n, len(dummy_levels)))
    for j, lv in enumerate(dummy_levels):
        D[:, j] = (raw == lv).astype(float)
    Bdes = np.hstack([np.ones((n, 1)), D])
    coef, *_ = np.linalg.lstsq(Bdes, pcs, rcond=None)
    pcs_resid = pcs - Bdes @ coef
    kept = pcs_resid.var(axis=0).sum() / max(pcs.var(axis=0).sum(), 1e-12)
    log(f"  operator residualisation: {len(dummy_levels)} dummy columns (of 12 frozen "
        f"candidates), PC variance retained {kept:.1%} (n={n})")
    return pcs_resid, dummy_levels


def assemble(M, xs, sel_mask):
    n1, n0 = int(sel_mask.sum()), int((~sel_mask).sum())
    ranks = rankdata(M, axis=0)
    auc = (ranks[sel_mask].sum(axis=0) - n1 * (n1 + 1) / 2.0) / (n1 * n0)
    rs = rankdata(xs)
    rs = (rs - rs.mean()) / rs.std()
    rc = (ranks - ranks.mean(axis=0)) / ranks.std(axis=0)
    rho = (rs[:, None] * rc).mean(axis=0)
    return auc, rho, n1, n0


def emp_p(obs, null):
    null = np.asarray(null, dtype=float)
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))


def auc_strat(M, sel_mask, strata):
    """Sum-of-U over sum-of-(n1*n0), ranks computed WITHIN each stratum. Column 0 of
    the returned array is the observed statistic; the rest is the null. Matches the
    independently-verified formula in c1_ucec_completeness_rerun.py."""
    labels = sorted(set(strata))
    per_stratum = {}
    U_sum = np.zeros(M.shape[1])
    denom_sum = 0.0
    for lab in labels:
        idx = strata == lab
        sub = M[idx]
        sub_sel = sel_mask[idx]
        n1, n0 = int(sub_sel.sum()), int((~sub_sel).sum())
        if n1 == 0 or n0 == 0:
            per_stratum[lab] = dict(n1=n1, n0=n0, auc=None)
            continue
        ranks = rankdata(sub, axis=0)
        U = ranks[sub_sel].sum(axis=0) - n1 * (n1 + 1) / 2.0
        per_stratum[lab] = dict(n1=n1, n0=n0, auc=(U / (n1 * n0)))
        U_sum += U
        denom_sum += n1 * n0
    combined = U_sum / denom_sum if denom_sum > 0 else np.full(M.shape[1], np.nan)
    return combined, per_stratum


def run_family(tag, tested, selected, R, P, pcs, perms, disc_incr_map, workers, log):
    res = SP.run_pool(SP._gene_test, len(tested), workers, R, P, pcs, perms, MIN_N)
    usable = [j for j in range(len(tested))
              if res[j] is not None and tested[j] in disc_incr_map.index
              and np.isfinite(disc_incr_map.get(tested[j], np.nan))]
    dropped = len(tested) - len(usable)
    if dropped:
        log(f"  [{tag}] {dropped} genes dropped (fewer than {MIN_N} confirmatory "
            f"protein values, or no discovery comparator)")
    if len(usable) < 100:
        sys.exit(f"FATAL [{tag}]: only {len(usable)} usable genes.")
    genes_u = [tested[j] for j in usable]
    M = np.vstack([res[j] for j in usable])
    xs = np.array([disc_incr_map[g] for g in genes_u])
    sel_mask = np.array([g in selected for g in genes_u])
    n1, n0 = int(sel_mask.sum()), int((~sel_mask).sum())
    if n1 == 0 or n0 == 0:
        sys.exit(f"FATAL [{tag}]: degenerate gene sets after intersection "
                 f"({n1} selected, {n0} background).")
    auc, rho, n1, n0 = assemble(M, xs, sel_mask)
    auc_obs, auc_null = float(auc[0]), auc[1:]
    rho_obs, rho_null = float(rho[0]), rho[1:]
    p_auc = emp_p(auc_obs, auc_null)
    p_rho = emp_p(rho_obs, rho_null)
    log(f"  [{tag}] n_usable={len(usable)} n_sel={n1} n_bg={n0}  "
        f"AUC={auc_obs:.4f} (null {auc_null.mean():.4f}+/-{auc_null.std():.4f}) "
        f"p={p_auc:.4f}  |  rho={rho_obs:.4f} p_rho={p_rho:.4f}")
    return dict(tag=tag, genes=genes_u, M=M, xs=xs, sel_mask=sel_mask, n_usable=len(usable),
               n_selected=n1, n_background=n0, auc_obs=auc_obs,
               auc_null_mean=float(auc_null.mean()), auc_null_sd=float(auc_null.std()),
               p_auc=p_auc, rho_obs=rho_obs, rho_null_mean=float(rho_null.mean()),
               rho_null_sd=float(rho_null.std()), p_rho=p_rho,
               negative_frac=float((M[:, 0] < 0).mean()))


def strip_arrays(d):
    """JSON-safe copy of a run_family() result (drops M/xs/sel_mask/genes)."""
    if d is None:
        return None
    return {k: v for k, v in d.items() if k not in ("M", "xs", "sel_mask", "genes")}


def clause_B(p1, p_strat, op_qualifies):
    if p1 is None or p1 >= 0.05:
        return None
    if not op_qualifies:
        return ("; acquisition-batch (operator)-stratified robustness not testable "
               "in this cohort")
    if p_strat < 0.05:
        return (f"; including under an acquisition-batch (operator)-stratified null "
               f"(p = {p_strat:.4f})")
    return (f" under an unrestricted permutation null but not under an "
           f"acquisition-batch (operator)-stratified null (p = {p_strat:.4f})")


def clause_D(step3_reached, p2, auc2):
    if not step3_reached:
        return None
    if p2 < 0.05:
        return (f"; the batch-corrected (operator) discovery set (726) also "
               f"outranked background under the operator-stratified null "
               f"(AUC {auc2:.4f}; p = {p2:.4f})")
    return (f"; the batch-corrected (operator) discovery set (726) was not shown "
           f"to outrank background under the operator-stratified null "
           f"(AUC {auc2:.4f}; p = {p2:.4f})")


def bucket_wording(p1, d2_reading=None):
    if p1 < 0.05:
        return "REPLICATED (set-level ranking, see clauses)"
    if p1 < 0.15:
        return "SUGGESTIVE (0.05<=p1<0.15); NOT declared a replication"
    tail = {
        "GO": "not replicated",
        "MARGINAL": "not replicated at marginal pre-checked power; not evidence of absence",
        "NO-GO": "uninformative at pre-checked power",
    }.get(d2_reading, "not replicated (D2 reading not supplied -- fill in before citing)")
    return f"NOT REPLICATED (p1>=0.15): {tail}"


def run_primary(args, log):
    t0 = time.time()
    common, tested, sel_h1, sel_h2, disc_incr, disc_fdr, disc_incr_bres, R, P, W, K = \
        prepare(args, log)
    n = len(common)

    if args.blind:
        log(f"\nBLIND SMOKE: n={n} K={K} n_tested={len(tested)} "
            f"n_h1_selected(usable)={len(sel_h1)} n_h2_selected(usable)={len(sel_h2)}")
        log("Nothing further is computed or printed in --blind mode (section 5.6).")
        log(f"elapsed {round((time.time() - t0) / 60.0, 1)} min")
        return

    raw_op = load_operator_labels(common, args.batch_labels, log)
    perms_op, op_qualifies, permutable, op_levels, op_members = \
        build_within_operator_perms(raw_op, args.perms, log)

    pcs = PCA(n_components=K, svd_solver="full").fit_transform(W)
    rng = np.random.RandomState(zlib.crc32(b"c1_luad|perm"))
    perms_unrestricted = [rng.permutation(n) for _ in range(args.perms)]

    p_floor = 1.0 / (args.perms + 1)
    exact_b = (args.perms == 1000)

    log("\n=== STEP 1: H1 (2,310-gene uncorrected set, unrestricted null) ===")
    h1 = run_family("H1", tested, sel_h1, R, P, pcs, perms_unrestricted, disc_incr,
                    args.workers, log)
    p1 = h1["p_auc"]

    log("\n=== STEP 2: (a)(i) batch qualifier (same 2,310 set, within-operator null) ===")
    if op_qualifies:
        a_i = run_family("a(i)", tested, sel_h1, R, P, pcs, perms_op, disc_incr,
                         args.workers, log)
        if abs(a_i["auc_obs"] - h1["auc_obs"]) > 1e-9:
            sys.exit(f"FATAL: (a)(i) observed AUC {a_i['auc_obs']!r} != H1 observed "
                     f"AUC {h1['auc_obs']!r} -- these must be bit-identical (same "
                     f"pcs, same gene sets; only the null differs). This is a bug.")
        p_strat = a_i["p_auc"]
    else:
        log("  operator stratum does NOT qualify on the realised analysed set -- "
            "(a)(i) is 'not testable' (D1/section 4); no other axis is substituted.")
        a_i = None
        p_strat = None

    step3_reached = bool(op_qualifies and p1 < 0.05 and p_strat is not None and p_strat < 0.05)
    log(f"\n=== STEP 3: H2 (726-gene batch-corrected set) "
        f"{'[REACHED: p1<0.05 and p_strat<0.05]' if step3_reached else '[not reached -- descriptive only]'} ===")
    pcs_resid, dummy_levels = build_operator_dummy_design(raw_op, pcs, log)
    if op_qualifies:
        h2 = run_family("H2", tested, sel_h2, R, P, pcs_resid, perms_op, disc_incr_bres,
                        args.workers, log)
        a_ii = run_family("a(ii)", tested, sel_h1, R, P, pcs_resid, perms_op,
                          disc_incr_bres, args.workers, log)
    else:
        log("  operator stratum not testable -- H2 and (a)(ii) computed descriptively "
            "with the H1 (unresidualised, unrestricted-null) estimator (D1 Case 2).")
        h2 = run_family("H2(Case2,descriptive)", tested, sel_h2, R, P, pcs,
                        perms_unrestricted, disc_incr_bres, args.workers, log)
        a_ii = None
    p2 = h2["p_auc"]

    log("\n=== Completeness-stratified AUC (co-reported, reuses H1's M) ===")
    disc_n = load_discovery_n(args.discovery_results)
    modal_n = int(disc_n.mode().iloc[0])
    if modal_n != DISCOVERY_MODAL_N_EXPECTED:
        sys.exit(f"FATAL: discovery file's modal per-gene n is {modal_n}, expected "
                 f"{DISCOVERY_MODAL_N_EXPECTED} -- wrong file vintage, stopping rather "
                 f"than stratifying at the wrong boundary.")
    # h1["genes"] is stripped only when building the JSON-safe `result` dict below;
    # the live h1 dict returned by run_family still carries it.
    genes_u = h1["genes"]
    missing = [g for g in genes_u if g not in disc_n.index]
    if missing:
        sys.exit(f"FATAL: {len(missing)} H1-usable genes have no row in "
                 f"{args.discovery_results} (e.g. {missing[:5]}) -- cannot stratify.")
    strata = np.array(["complete" if disc_n[g] == modal_n else "incomplete" for g in genes_u])
    confirm_n = (~np.isnan(P[:, [tested.index(g) for g in genes_u]])).sum(axis=0)
    auc_compl_col, compl_info = auc_strat(h1["M"], h1["sel_mask"], strata)
    auc_compl_obs, auc_compl_null = float(auc_compl_col[0]), auc_compl_col[1:]
    p_compl = emp_p(auc_compl_obs, auc_compl_null)
    # pooled rho within-stratum (informational): rank-based rho computed within each
    # stratum's own xs/M[:,0], then case-weighted by stratum size
    within_rho = {}
    for lab in sorted(set(strata)):
        idx = strata == lab
        if idx.sum() < 3:
            continue
        rr = rankdata(h1["xs"][idx])
        cc = rankdata(h1["M"][idx, 0])
        within_rho[lab] = float(np.corrcoef(rr, cc)[0, 1])
    log(f"  AUC_compl observed={auc_compl_obs:.4f} null={auc_compl_null.mean():.4f}"
        f"+/-{auc_compl_null.std():.4f} p_compl={p_compl:.4f}")
    for lab, info in sorted(compl_info.items()):
        fp = float((h1["M"][strata == lab, 0] > 0).mean())
        med_n = float(np.median(confirm_n[strata == lab]))
        log(f"    stratum {lab:<11} n_sel={info['n1']:<5} n_bg={info['n0']:<5} "
            f"frac_positive={fp:.4f} median_confirm_n={med_n:g} within_rho={within_rho.get(lab, float('nan')):.4f} "
            + (f"AUC={info['auc'][0]:.4f}" if info["auc"] is not None else "AUC=n/a"))

    # ---- wording ----
    log("\n" + "=" * 78)
    log("C1-LUAD CONFIRMATORY TEST -- SECTION 5 READING")
    log("=" * 78)
    if not exact_b:
        log(f"NO READING: B={args.perms} is not the frozen analysis (section 5.6 "
            f"requires B=1000 exactly). Smallest attainable p is {p_floor:.4f}; "
            f"p1={p1:.4f} is "
            + ("AT that floor (censored, not a negative result)."
               if p1 <= p_floor + 1e-12 else "resolution-limited."))
    else:
        log(f"p1 = {p1:.4f}  -> {bucket_wording(p1, args.d2_reading)}")
        cB = clause_B(p1, p_strat, op_qualifies)
        if cB is not None:
            log(f"clause [B]: {cB}")
        if p1 < 0.05 and p_compl >= 0.05:
            log(f"clause [C]: ; not robust to protein-completeness stratification "
                f"(p = {p_compl:.4f})")
        cD = clause_D(step3_reached, p2, h2["auc_obs"])
        if cD is not None:
            log(f"clause [D]: {cD}")
        elif p1 < 0.05:
            log("H2 not tested in the confirmatory sequence (step 3 not reached).")

    minutes = round((time.time() - t0) / 60.0, 1)
    result = dict(
        rule_path=args.rule, rule_sha256=sha256(args.rule), rule_sha256_expected=args.rule_sha256,
        cohort="luad_c1", n_cases=n, K=K, perms=args.perms, exact_b=exact_b,
        n_genes_tested=len(tested),
        h1=strip_arrays(h1), p1=p1,
        a_i=strip_arrays(a_i), p_strat=p_strat,
        operator_qualifies=bool(op_qualifies), operator_permutable_frac=float(permutable),
        operator_dummy_levels=dummy_levels,
        step3_reached=step3_reached,
        h2=strip_arrays(h2), p2=p2,
        a_ii=strip_arrays(a_ii),
        completeness=dict(auc_compl_obs=auc_compl_obs,
                          auc_compl_null_mean=float(auc_compl_null.mean()),
                          auc_compl_null_sd=float(auc_compl_null.std()),
                          p_compl=p_compl,
                          per_stratum={lab: dict(n_selected=info["n1"], n_background=info["n0"],
                                                 auc_obs=(float(info["auc"][0]) if info["auc"] is not None else None),
                                                 frac_positive=float((h1["M"][strata == lab, 0] > 0).mean()),
                                                 median_confirm_n=float(np.median(confirm_n[strata == lab])),
                                                 within_stratum_rho=within_rho.get(lab))
                                      for lab, info in compl_info.items()}),
        d2_reading=args.d2_reading, blind=False,
        minutes=minutes,
        sha256_this_script=sha256(os.path.abspath(__file__)),
        sha256_discovery_results=sha256(args.discovery_results),
        sha256_discovery_sitepack=sha256(args.discovery_sitepack),
        sha256_residual_analysis=sha256(os.path.join(_SCRIPT_DIR, "residual_analysis.py")),
        sha256_split_replication_power=sha256(os.path.join(_SCRIPT_DIR, "split_replication_power.py")),
        sha256_batch_labels=sha256(args.batch_labels),
        sha256_case_ids=sha256(args.case_ids),
    )

    suffix = f"_{args.tag}" if args.tag else ""
    npz_fp = RA.OUT_DIR / f"c1_luad_primary{suffix}.npz"
    np.savez_compressed(npz_fp, genes=np.array(genes_u), sel_mask=h1["sel_mask"],
                        strata=strata, M=h1["M"], confirm_n=confirm_n)
    result["sha256_npz"] = sha256(str(npz_fp))

    gene_tbl = pd.DataFrame({
        "gene": genes_u,
        "h1_selected": h1["sel_mask"],
        "h2_selected": np.array([g in sel_h2 for g in genes_u]),
        "discovery_incremental_r2": h1["xs"],
        "discovery_fdr": [disc_fdr.get(g, float("nan")) for g in genes_u],
        "confirmatory_incremental_r2": h1["M"][:, 0],
        "confirmatory_n": confirm_n,
        "completeness_stratum": strata,
    }).sort_values("confirmatory_incremental_r2", ascending=False)
    gene_fp = RA.OUT_DIR / f"c1_luad_primary_gene_level{suffix}.csv"
    gene_tbl.to_csv(gene_fp, index=False)

    res_fp = RA.OUT_DIR / f"c1_luad_primary_result{suffix}.json"
    with open(res_fp, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=lambda o: (o.tolist() if hasattr(o, "tolist") else str(o)))
    log(f"\n  -> {res_fp}\n  -> {gene_fp}\n  -> {npz_fp}\n  elapsed {minutes} min")


def run_sensitivity(args, log):
    t0 = time.time()
    if args.variant not in ("s", "s2"):
        sys.exit("FATAL: --variant must be 's' or 's2'")
    drop = SENSITIVITY_S_DROP_CASES if args.variant == "s" else SENSITIVITY_S2_DROP_CASES
    seed_tag = b"c1_luad|perm_s" if args.variant == "s" else b"c1_luad|perm_s2"

    common, tested, sel_h1, sel_h2, disc_incr, disc_fdr, disc_incr_bres, R, P, W, K = \
        prepare(args, log, drop_cases=drop)
    n = len(common)
    if args.blind:
        log(f"\nBLIND SMOKE ({args.variant}): n={n} K={K} n_tested={len(tested)}")
        return

    pcs = PCA(n_components=K, svd_solver="full").fit_transform(W)
    rng = np.random.RandomState(zlib.crc32(seed_tag))
    perms = [rng.permutation(n) for _ in range(args.perms)]
    h1s = run_family(f"H1({args.variant})", tested, sel_h1, R, P, pcs, perms, disc_incr,
                     args.workers, log)

    minutes = round((time.time() - t0) / 60.0, 1)
    result = dict(cohort="luad_c1", variant=args.variant, n_cases=n, K=K, perms=args.perms,
                 dropped_cases=sorted(drop), h1=strip_arrays(h1s), p_auc=h1s["p_auc"],
                 minutes=minutes, sha256_this_script=sha256(os.path.abspath(__file__)))
    suffix = f"_{args.tag}" if args.tag else ""
    res_fp = RA.OUT_DIR / f"c1_luad_sensitivity_{args.variant}_result{suffix}.json"
    with open(res_fp, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=lambda o: (o.tolist() if hasattr(o, "tolist") else str(o)))
    log(f"\n  -> {res_fp}\n  elapsed {minutes} min")


def run_bootstrap(args, log):
    t0 = time.time()
    common, tested, sel_h1, sel_h2, disc_incr, disc_fdr, disc_incr_bres, R, P, W, K = \
        prepare(args, log)
    n = len(common)
    if args.blind:
        log(f"\nBLIND SMOKE (bootstrap): n={n} K={K} n_tested={len(tested)}")
        return

    raw_op = load_operator_labels(common, args.batch_labels, log)
    log(f"patient-level bootstrap: {args.boots} resamples of {n} cases (with "
        f"replacement), PCA refit per resample, K held at {K}. H2 resample: each "
        f"resample's K PCs are residualised on the FROZEN operator dummies (section 4) "
        f"using that resample's (duplicated) case labels. Caveat: a duplicated case can "
        f"fall in both the train and test fold of the 5-fold CV, which inflates "
        f"absolute R2; the AUC is a contrast between two gene sets on the same "
        f"resample, so the effect largely cancels, but the CI is not bias-free.")

    rng = np.random.RandomState(zlib.crc32(b"c1_luad|boot"))
    aucs_h1, aucs_h2 = [], []
    for b in range(args.boots):
        idx = rng.randint(0, n, n)
        pcs_b = PCA(n_components=K, svd_solver="full").fit_transform(W[idx])
        res = SP.run_pool(SP._gene_test, len(tested), args.workers, R[idx], P[idx],
                          pcs_b, [], MIN_N)
        usable = [j for j in range(len(tested)) if res[j] is not None
                  and np.isfinite(disc_incr.get(tested[j], np.nan))]
        if len(usable) >= 100:
            genes_u = [tested[j] for j in usable]
            M = np.vstack([res[j] for j in usable])
            xs = np.array([disc_incr[g] for g in genes_u])
            sel_mask = np.array([g in sel_h1 for g in genes_u])
            if sel_mask.any() and (~sel_mask).any():
                auc, *_ = assemble(M, xs, sel_mask)
                aucs_h1.append(float(auc[0]))

        if args.h2_claimed:
            raw_b = raw_op[idx]
            pcs_b_resid, _ = build_operator_dummy_design(raw_b, pcs_b, lambda *_a: None)
            res2 = SP.run_pool(SP._gene_test, len(tested), args.workers, R[idx], P[idx],
                               pcs_b_resid, [], MIN_N)
            usable2 = [j for j in range(len(tested)) if res2[j] is not None
                      and np.isfinite(disc_incr_bres.get(tested[j], np.nan))]
            if len(usable2) >= 100:
                genes_u2 = [tested[j] for j in usable2]
                M2 = np.vstack([res2[j] for j in usable2])
                xs2 = np.array([disc_incr_bres[g] for g in genes_u2])
                sel_mask2 = np.array([g in sel_h2 for g in genes_u2])
                if sel_mask2.any() and (~sel_mask2).any():
                    auc2, *_ = assemble(M2, xs2, sel_mask2)
                    aucs_h2.append(float(auc2[0]))

        if (b + 1) % 25 == 0:
            msg = f"  {b + 1}/{args.boots} resamples, H1 AUC mean {np.mean(aucs_h1):.4f}"
            if args.h2_claimed and aucs_h2:
                msg += f", H2 AUC mean {np.mean(aucs_h2):.4f}"
            log(msg + f", {(time.time() - t0) / 60:.0f} min")

    def ci(a):
        a = np.array(a)
        if len(a) == 0:
            return dict(boots=0, mean=None, sd=None, ci95_lo=None, ci95_hi=None)
        lo, hi = np.percentile(a, [2.5, 97.5])
        return dict(boots=len(a), mean=float(a.mean()), sd=float(a.std()),
                   ci95_lo=float(lo), ci95_hi=float(hi))

    result = dict(cohort="luad_c1", n_cases=n, K=K, requested_boots=args.boots,
                 h2_claimed=bool(args.h2_claimed),
                 h1=ci(aucs_h1), h2=(ci(aucs_h2) if args.h2_claimed else None),
                 minutes=round((time.time() - t0) / 60.0, 1),
                 sha256_this_script=sha256(os.path.abspath(__file__)))
    suffix = f"_{args.tag}" if args.tag else ""
    fp = RA.OUT_DIR / f"c1_luad_bootstrap_result{suffix}.json"
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    log(f"\n  H1 bootstrap: mean={result['h1']['mean']} CI=[{result['h1']['ci95_lo']}, "
        f"{result['h1']['ci95_hi']}] ({result['h1']['boots']} resamples)")
    if args.h2_claimed:
        log(f"  H2 bootstrap: mean={result['h2']['mean']} CI=[{result['h2']['ci95_lo']}, "
            f"{result['h2']['ci95_hi']}] ({result['h2']['boots']} resamples)")
    log(f"  -> {fp}\n  elapsed {result['minutes']} min")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["primary", "sensitivity", "bootstrap", "descriptive"])
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--boots", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=RA.N_WORKERS)
    ap.add_argument("--tag", default="")
    ap.add_argument("--blind", action="store_true",
                    help="section 5.6/D5g: stop after n/K/gene-count logging; compute "
                         "nothing else. Required for any confirmatory smoke test.")
    ap.add_argument("--variant", choices=["s", "s2"], default="",
                    help="sensitivity mode only")
    ap.add_argument("--h2-claimed", action="store_true",
                    help="bootstrap mode only: also bootstrap H2. Pass this ONLY after "
                         "a `primary` run gave p1<0.05 AND p_strat<0.05 (step 3 reached).")
    ap.add_argument("--d2-reading", choices=["GO", "MARGINAL", "NO-GO"], default=None,
                    help="D2 pre-check outcome, for the >=0.15 bucket wording (section "
                         "5.1). Not required to run; only affects the printed wording.")
    # --- required paths, no defaults (Appendix item 1) ---
    ap.add_argument("--rule", required=True, help="path to the SIGNED frozen rule file")
    ap.add_argument("--rule-sha256", required=True,
                    help="expected sha256 of --rule (from section 7.3)")
    ap.add_argument("--case-ids", required=True, help="pinned/c1_case_ids_luad.txt")
    ap.add_argument("--protein", required=True, help="protein_luad_c1.csv (c1_pull_omics_luad.py output)")
    ap.add_argument("--discovery-results", required=True, help="luad__residual_results_tumoronly.csv")
    ap.add_argument("--discovery-sitepack", required=True, help="luad__sitepack_operator.csv")
    ap.add_argument("--discovery-slide-map", required=True, help="HPC scripts/pinned/slide_type_map_luad.tsv")
    ap.add_argument("--discovery-xwalk", required=True, help="HPC scripts/pinned/aliquot_to_case_tumor_luad.tsv")
    ap.add_argument("--discovery-rna-manifest", required=True, help="HPC scripts/pinned/manifest_rna_tumor_luad.tsv")
    ap.add_argument("--batch-labels", required=True, help="pinned/c1_luad_batch_labels.tsv")
    args = ap.parse_args()

    def log(msg):
        print(msg, flush=True)

    log("=== resolved paths (verify before trusting any number below) ===")
    log(f"  MORPHO_ROOT        : {RA.ROOT}")
    log(f"  MORPHO_RNA_MANIFEST: {RA.RNA_MANIFEST}")
    log(f"  MORPHO_SLIDE_MAP   : {RA.SLIDE_TYPE_MAP}")
    log(f"  MORPHO_WSI_EMB_DIR : {RA.WSI_EMB_DIR}")
    log(f"  MORPHO_OUT         : {RA.OUT_DIR}")
    log(f"  rule               : {args.rule}")
    log(f"  protein            : {args.protein}")
    log(f"  discovery results  : {args.discovery_results}")
    log(f"  discovery sitepack : {args.discovery_sitepack}")
    log(f"  batch labels       : {args.batch_labels}")
    log(f"  mode / workers     : {args.mode} / {args.workers}  blind={args.blind}")
    log("=" * 78)

    got_rule_hash = sha256(args.rule)
    if got_rule_hash != args.rule_sha256:
        sys.exit(f"FATAL: --rule sha256 mismatch: got={got_rule_hash} "
                 f"want={args.rule_sha256}. Refusing to run against an unverified "
                 f"copy of the frozen rule.")
    log(f"  hash OK: signed rule = {got_rule_hash[:12]}...")

    if args.mode == "primary":
        run_primary(args, log)
    elif args.mode == "sensitivity":
        run_sensitivity(args, log)
    elif args.mode == "bootstrap":
        run_bootstrap(args, log)
    elif args.mode == "descriptive":
        log("\n'descriptive' is not implemented in this script -- run "
            "c1_luad_pergene.py directly for (b1)/(b2) (Appendix item 3). This mode "
            "exists only so `c1_run_test_luad.py descriptive` gives a clear pointer "
            "instead of an argparse error.")
        sys.exit(2)


if __name__ == "__main__":
    main()
