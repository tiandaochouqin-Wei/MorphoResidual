#!/usr/bin/env python3
"""
c1_ucec_plex_reread.py -- POST HOC TMT-plex-stratified re-read of the published
C1-UCEC confirmatory test (PDC000439; 138 cases).

STATUS AND AUTHORITY
  POST HOC. Governed by review/POSTHOC_UCEC_C1_PLEX_PRESPEC_2026-10-08.md (sidecar
  .sha256), released on GitHub before the --perms 1000 run; the --blind and --perms 5
  smoke steps may run earlier and count only if this script's sha256 in their LSF .out
  equals the sidecar value. Nothing here changes the published C1-UCEC reading (AUC
  0.597443, unrestricted null 0.498195 +/- 0.036296, p = 0.003996) or the D4
  within-operator re-read. No confirmatory null of the C1-UCEC test was stratified by TMT
  plex before this script (POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md, section 3).

WHAT IT DOES (adapted from c1_luad_plex_reread.py and c1_ucec_stratified_reread.py)
  The frozen machinery of c1_run_test.py (prepare(), assemble(), MIN_N), of
  split_replication_power.py (run_pool, _gene_test, k_for) and of residual_analysis.py is
  IMPORTED, not reimplemented. New code: the within-plex permutation null and the plex
  dummy design (both copied from c1_luad_plex_reread.py), and the driver.

  0. Pre-flight, before any data is read: MORPHO_* environment (must be the ucec_c1 tree),
     sha256 of every pinned input and of the three frozen modules. Any mismatch is FATAL.
  1. Reproduction gate (FATAL on any mismatch). Recompute the published C1-UCEC test
     exactly (c1_run_test.run_test up to the genes x (1+B) matrix: same prepare(), PCA
     call, seed crc32(b"c1_ucec|perm"), run_pool, usable filter, gene order, selection
     mask) and validate auc_obs, auc_null_mean, auc_null_sd, p_auc, rho_obs, p_rho (abs
     tol 1e-9) and n_usable, n_selected, n_background, n_cases, K, n_genes_tested and
     the code/discovery hashes against the sha256-pinned c1_ucec_result.json. If the D4
     matrix file c1_ucec_d4_stratified.npz is present on the host, its gene order and
     observed column are compared too (same pinned sha256 as in the repository).
  2. P1 (primary). The published 2,566-gene set (2,536 usable), the SAME PCs, under a
     within-plex null: B = 1000 draws from RandomState(crc32(b"c1_ucec|perm_plex")),
     strata in sorted label order, singleton plexes and unlabelled cases held fixed (here
     0 and 0). The observed AUC and rho are asserted identical to step 1 (only the null
     changes). Reported: AUC/rho, null mean/SD, z, one-sided p, null-SD ratio against
     the unrestricted null.
  3. P2 (descriptive). Same set, the K confirmatory PCs residualised on plex dummies
     (every plex with >= 3 cases in Counter.most_common order; the last is dropped for
     full rank when the dummies cover every case; least squares with intercept), SAME
     within-plex draws. This is the analogue of the LUAD P2 (and of D4's descriptive
     "extra" arm). It uses the PUBLISHED set because UCEC has no registered
     batch-corrected set; LUAD P2 used its registered 726-gene H2 set.
  4. P3 (sensitivity). P1 without C3L-05571 and C3N-02978: each has two Primary Tumor
     aliquots, one Disqualified (plex 04 / 12) and a Qualified re-run in plex 16, and the
     published protein matrix averaged every non-Normal aliquot of a case
     (c1_pull_omics.py), so these two protein values may mix two plexes. Rows are
     dropped after prepare(), which is identical to dropping them inside it; K is
     recomputed (round(20 * 136 / 100) = 27) and the PCA refit; the within-plex null
     comes from a fresh RandomState(crc32(b"c1_ucec|perm_plex")), i.e. P1's procedure
     repeated literally.

MODES
  --perms 1000 (default)  the analysis; writes result.json, .npz, .log.
  --perms N, 1 <= N < 1000  smoke: reproduces step 1 at B = N, checks auc_obs and the
                          counts against the published result, then stops. No plex null
                          is computed, nothing is written.
  --blind                 stops after prepare() and the plex-label load; prints counts
                          only (n, K, genes, plex levels and sizes, dummies). Nothing
                          is written.

OUTPUTS (in $MORPHO_OUT = /public/home/fjhui/ZW/ucec_c1/results/posthoc; never
  overwritten; JSON and npz are opened with mode "x", the log is append-only)
  c1_ucec_plex_reread_result[_<tag>].json
  c1_ucec_plex_reread[_<tag>].log
  c1_ucec_plex_reread[_<tag>].npz   genes x (1+B) matrices for step 1, P1, P2, P3
No name collides with c1_ucec_result.json, c1_ucec_gene_level.csv or the D4 files.

RUN: through lsf_c1_ucec_plex_reread.sh (sets MORPHO_*). See that wrapper's header.
Expected wall time at B = 1000 on 8 cores: ~3.4 h (four passes of ~51 min; D4's three
passes took 153.5 min).
"""
import argparse
import hashlib
import json
import os
import platform
import re
import socket
import sys
import time
import zlib
from collections import Counter

import numpy as np
import pandas as pd

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

POST_HOC_GOVERNED_BY = "review/POSTHOC_UCEC_C1_PLEX_PRESPEC_2026-10-08.md"

# ---------------------------------------------------------------------------
# Pinned constants. Changing any of these changes the analysis.
# ---------------------------------------------------------------------------
FLOAT_TOL = 1e-9
FROZEN_B = 1000
MIN_GROUP = 3                       # residual_analysis_sitepack.py / c1_luad MIN_GROUP
SEED_UNRESTRICTED = b"c1_ucec|perm"       # c1_run_test.run_test, reproduced
SEED_PLEX = b"c1_ucec|perm_plex"          # P1, P2 (shared draws) and P3 (fresh RandomState)
P3_DROP_CASES = frozenset({"C3L-05571", "C3N-02978"})

PUBLISHED_RESULT_SHA256 = "53fcb86138465e90ef2c47e3519786fc0fb4a534497f4a83850ee584ee6e653c"
D4_NPZ_SHA256 = "ef75f1c06d4fc33843a79e9a8d1bb4f4224acd201d30784badb97cf3be13384a"
PLEX_MAP_SHA256 = "4fe057a9ac57277db886365330fd3ae254f9dbcaff4581f28f68e2a57589dd26"
FROZEN_INPUT_SHA256 = {
    "case_ids": "d777fc88c1dd7b9eb817b0a6a1366b1b909386bd34dea857e4bc53cffc5807dc",
    "discovery_results": "100d8f8581c78905c2a98caf4bb2713e376d9264c6f00276446a9893b987c4a6",
}
FROZEN_CODE_SHA256 = {
    "residual_analysis": "a54a0a56494d315e15ef075425fe2133ecd9f05193c22656c1145bffdaac4a0a",
    "split_replication_power": "79530c075c0f7e75bf06d3c2f1edfe8792cbf4d6c71ae1d32e60bf44f31a0397",
    # local copy of the repository mirror of the HPC scripts; verify with sha256sum on the
    # HPC before the smoke run (the smoke run stops at pre-flight if it differs)
    "c1_run_test": "fc578c705e803ddb118ff1f65a61079723167985a66adf7285834bf4709391b5",
}
# What the reproduction must give (c1_ucec_result.json, 2026-09-16; also in the sha-pinned file)
GATE_EXPECTED = {"n_cases": 138, "K": 28, "n_genes_tested": 10153, "n_usable": 10146,
                 "n_selected": 2536, "n_background": 7610}

REQUIRED_ENV = ("MORPHO_ROOT", "MORPHO_RNA_MANIFEST", "MORPHO_SLIDE_MAP",
                "MORPHO_WSI_EMB_DIR", "MORPHO_OUT")
_TAG = re.compile(r"^[A-Za-z0-9_-]{0,40}$")
_CASE = re.compile(r"^C3[LN]-\d{5}$")


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------
def sha256(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return "unreadable"


def fatal(msg):
    sys.exit(f"FATAL: {msg}")


def assert_sha(path, want, label, log):
    got = sha256(path)
    if got != want:
        fatal(f"sha256 mismatch for {label}: got={got} want={want} ({path}). Wrong file "
              f"vintage; stopping before anything is computed.")
    log(f"  hash OK: {label} = {got[:12]}... ({path})")
    return got


class Tee:
    """print() plus, once open() is called, a line-buffered APPEND-mode copy to a log
    file, so a failed attempt stays on record."""

    def __init__(self):
        self.fh = None

    def open(self, path):
        self.fh = open(path, "a", encoding="utf-8", newline="\n")

    def __call__(self, msg=""):
        print(msg, flush=True)
        if self.fh is not None:
            self.fh.write(f"{msg}\n")
            self.fh.flush()

    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None


def _z(obs, mean, sd):
    return float((obs - mean) / sd) if sd > 0 else float("nan")


def _utc():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def emp_p(obs, null):
    """One-sided upper-tail empirical p: (1 + #{null >= obs}) / (1 + B), the convention
    of c1_run_test.py and c1_ucec_stratified_reread.py."""
    null = np.asarray(null, dtype=float)
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))


# ---------------------------------------------------------------------------
# Plex labels, within-plex null, plex dummy design (copied from c1_luad_plex_reread.py)
# ---------------------------------------------------------------------------
def read_plex_map(path):
    """case -> plex label. First two columns must be case_id (or case_submitter_id) and
    plex. Duplicate case rows and malformed case IDs are FATAL. An empty plex cell means
    unlabelled."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    cols = list(df.columns)
    if len(cols) < 2 or cols[0] not in ("case_id", "case_submitter_id") or cols[1] != "plex":
        fatal(f"{path}: expected the first two columns to be case_id (or "
              f"case_submitter_id) and plex; got {cols[:2]}")
    cases = [c.strip() for c in df[cols[0]].tolist()]
    plex = [p.strip() for p in df[cols[1]].tolist()]
    bad = [c for c in cases if not _CASE.match(c)]
    if bad:
        fatal(f"{path}: {len(bad)} malformed case IDs, e.g. {bad[:5]}")
    dup = sorted({c for c, k in Counter(cases).items() if k > 1})
    if dup:
        fatal(f"{path}: {len(dup)} case IDs appear more than once, e.g. {dup[:5]}")
    return dict(zip(cases, plex))


def labels_for(cases, mapping):
    """np.ndarray of plex labels aligned to `cases`; '' for a case with no row or an empty
    label (held fixed in the within-plex null, never a level)."""
    return np.array([mapping.get(c, "") for c in cases], dtype=object).astype(str)


def strata_diagnostics(raw):
    raw = np.asarray(raw)
    n = len(raw)
    labelled = [v for v in raw.tolist() if v]
    counts = Counter(labelled)
    n_unlabelled = n - len(labelled)
    singletons = sorted(lv for lv, k in counts.items() if k == 1)
    n_fixed = n_unlabelled + len(singletons)
    return dict(
        n_cases=n, n_levels=len(counts), n_unlabelled=n_unlabelled,
        n_singleton_strata=len(singletons), singleton_strata=singletons,
        n_cases_fixed=n_fixed,
        frac_cases_fixed=(n_fixed / n) if n else float("nan"),
        frac_permutable=(1.0 - n_fixed / n) if n else float("nan"),
        n_strata_ge3=sum(1 for k in counts.values() if k >= 3),
        largest_stratum=max(counts.values()) if counts else 0,
        smallest_stratum=min(counts.values()) if counts else 0,
        stratum_sizes={lv: int(counts[lv]) for lv in sorted(counts)},
    )


def build_within_strata_perms(raw, B, seed_bytes):
    """B within-stratum permutations of range(len(raw)). Unlabelled cases ('') are not a
    level and never move; a stratum of size 1 never moves. Strata are processed in sorted
    label order inside every draw, from ONE RandomState(zlib.crc32(seed_bytes)). This is
    the loop of c1_ucec_stratified_reread.build_within_operator_perms."""
    raw = np.asarray(raw)
    n = len(raw)
    levels = sorted(set(v for v in raw if v))
    members = {lv: np.where(raw == lv)[0] for lv in levels}
    rng = np.random.RandomState(zlib.crc32(seed_bytes))
    perms = []
    for _ in range(B):
        p = np.arange(n)
        for lv in levels:
            m = members[lv]
            if len(m) > 1:
                p[m] = rng.permutation(m)
        perms.append(p)
    return perms


def plex_dummy_levels(raw, min_group):
    """Levels with >= min_group analysed cases, in Counter.most_common order (count
    descending, ties by first appearance in `raw`, i.e. sorted case order). If those
    dummies cover every analysed case (no unlabelled case and no small level left to the
    intercept) the LAST one is dropped so the design is full rank. Unlabelled cases are
    never a level."""
    raw = np.asarray(raw)
    labelled = [v for v in raw.tolist() if v]
    counts = Counter(labelled)
    dummy = [lv for lv, k in counts.most_common() if k >= min_group]
    left_to_intercept = sorted(lv for lv in counts if lv not in set(dummy))
    dropped = None
    if dummy and len(labelled) == len(raw) and not left_to_intercept:
        dropped = dummy[-1]
        dummy = dummy[:-1]
    return dummy, dropped, left_to_intercept


def build_plex_dummy_design(raw, pcs, min_group):
    """Residualise pcs on [1, plex dummies] by least squares (np.linalg.lstsq). FATAL
    unless the design has full column rank."""
    raw = np.asarray(raw)
    dummy, dropped, left = plex_dummy_levels(raw, min_group)
    n = len(raw)
    D = np.zeros((n, len(dummy)))
    for j, lv in enumerate(dummy):
        D[:, j] = (raw == lv).astype(float)
    Bdes = np.hstack([np.ones((n, 1)), D])
    rank = int(np.linalg.matrix_rank(Bdes))
    if rank != Bdes.shape[1]:
        fatal(f"plex dummy design is rank deficient ({rank} < {Bdes.shape[1]} columns).")
    coef, *_ = np.linalg.lstsq(Bdes, pcs, rcond=None)
    pcs_resid = pcs - Bdes @ coef
    kept = float(pcs_resid.var(axis=0).sum() / max(pcs.var(axis=0).sum(), 1e-12))
    info = dict(n_dummies=len(dummy), dummy_levels=dummy, dropped_for_full_rank=dropped,
                levels_left_to_intercept=left, design_rank=rank,
                pc_variance_retained=kept, min_group=int(min_group))
    return pcs_resid, info


# ---------------------------------------------------------------------------
# Pre-flight
# ---------------------------------------------------------------------------
def env_guard():
    env = {}
    for k in REQUIRED_ENV:
        v = os.environ.get(k, "")
        if not v:
            fatal(f"{k} is not set. Run through lsf_c1_ucec_plex_reread.sh, which exports "
                  f"every MORPHO_* before python starts (residual_analysis.py reads them "
                  f"at import time).")
        env[k] = v
    root = os.path.normpath(env["MORPHO_ROOT"])
    out = os.path.normpath(env["MORPHO_OUT"])
    if os.path.basename(root) != "ucec_c1":
        fatal(f"basename(MORPHO_ROOT)={os.path.basename(root)!r}, expected 'ucec_c1'.")
    if not (out == root or out.startswith(root + os.sep)):
        fatal(f"MORPHO_OUT={out!r} does not lie inside the ucec_c1 tree "
              f"(MORPHO_ROOT={root!r}).")
    return env


def import_frozen(env, log):
    """Hash the three frozen modules from the files the import will actually load, then
    import. Pre-flight comes before the import because residual_analysis reads the
    environment at import time."""
    import importlib.util
    code_hashes = {}
    origins = {}
    for name, want in FROZEN_CODE_SHA256.items():
        spec = importlib.util.find_spec(name)
        if spec is None or not spec.origin:
            fatal(f"cannot locate module {name} on sys.path {sys.path[:3]}")
        origins[name] = os.path.realpath(spec.origin)
        code_hashes[name] = assert_sha(origins[name], want, f"frozen module {name}", log)
    import c1_run_test as C1  # noqa: E402  (imports residual_analysis, split_replication_power)
    RA, SP = C1.RA, C1.SP
    for name, mod in (("c1_run_test", C1), ("residual_analysis", RA),
                      ("split_replication_power", SP)):
        if os.path.realpath(mod.__file__) != origins[name]:
            fatal(f"imported {name} from {mod.__file__}, but hashed {origins[name]}.")
    for attr, key in (("ROOT", "MORPHO_ROOT"), ("RNA_MANIFEST", "MORPHO_RNA_MANIFEST"),
                      ("SLIDE_TYPE_MAP", "MORPHO_SLIDE_MAP"),
                      ("WSI_EMB_DIR", "MORPHO_WSI_EMB_DIR"), ("OUT_DIR", "MORPHO_OUT")):
        if os.path.normpath(str(getattr(RA, attr))) != os.path.normpath(env[key]):
            fatal(f"residual_analysis.{attr}={getattr(RA, attr)} != ${key}={env[key]}.")
    return C1, RA, SP, code_hashes


def family(C1, SP, tested, selected, disc_incr, R, P, pcs, perms, workers, log, name):
    """One genes x (1+B) pass exactly as c1_run_test.run_test builds M, then AUC and rho
    through C1.assemble. Returns a dict with the arrays."""
    log(f"[{name}] per-gene pass: {len(tested)} genes x (1 observed + {len(perms)} "
        f"permutations) on {workers} workers ...")
    t = time.time()
    res = SP.run_pool(SP._gene_test, len(tested), workers, R, P, pcs, perms, C1.MIN_N)
    usable = [j for j in range(len(tested)) if res[j] is not None
              and np.isfinite(disc_incr.get(tested[j], np.nan))]
    if len(usable) < 100:
        fatal(f"only {len(usable)} usable genes")
    genes_u = [tested[j] for j in usable]
    M = np.vstack([res[j] for j in usable])
    xs = np.array([disc_incr[g] for g in genes_u])
    sel_mask = np.array([g in selected for g in genes_u])
    if sel_mask.sum() == 0 or (~sel_mask).sum() == 0:
        fatal(f"degenerate gene sets ({int(sel_mask.sum())} selected, "
              f"{int((~sel_mask).sum())} background)")
    auc, rho, n1, n0 = C1.assemble(M, xs, sel_mask)
    log(f"[{name}] done in {(time.time() - t) / 60.0:.1f} min: n_usable={len(usable)} "
        f"n_selected={n1} n_background={n0} AUC_obs={auc[0]:.6f}")
    return dict(genes=genes_u, M=M, xs=xs, sel_mask=sel_mask, auc=auc, rho=rho,
                n_usable=len(usable), n_selected=int(n1), n_background=int(n0))


def stats_of(fam):
    auc, rho = fam["auc"], fam["rho"]
    a_obs, a_null = float(auc[0]), auc[1:]
    r_obs, r_null = float(rho[0]), rho[1:]
    return dict(
        n_usable=fam["n_usable"], n_selected=fam["n_selected"], n_background=fam["n_background"],
        perms=int(len(a_null)),
        auc_obs=a_obs, auc_null_mean=float(a_null.mean()), auc_null_sd=float(a_null.std()),
        z_auc=_z(a_obs, a_null.mean(), a_null.std()), p_auc=emp_p(a_obs, a_null),
        exceedances_auc=int(np.sum(a_null >= a_obs)),
        margin_auc=float(a_obs - a_null.mean()),
        rho_obs=r_obs, rho_null_mean=float(r_null.mean()), rho_null_sd=float(r_null.std()),
        z_rho=_z(r_obs, r_null.mean(), r_null.std()), p_rho=emp_p(r_obs, r_null),
        exceedances_rho=int(np.sum(r_null >= r_obs)),
        n_neg_obs_frac=float((fam["M"][:, 0] < 0).mean()),
    )


def validate(log, computed, published_path, tol=FLOAT_TOL):
    """FATAL unless every field of `computed` matches the sha-pinned published
    c1_ucec_result.json (identical to c1_ucec_stratified_reread.validate)."""
    if not os.path.exists(published_path):
        fatal(f"--validate-against file does not exist: {published_path}")
    got_hash = sha256(published_path)
    if got_hash != PUBLISHED_RESULT_SHA256:
        fatal(f"{published_path} has sha256={got_hash}, expected the pinned 2026-09-16 "
              f"result {PUBLISHED_RESULT_SHA256}.")
    with open(published_path, "r", encoding="utf-8") as f:
        published = json.load(f)
    float_fields = ["auc_obs", "auc_null_mean", "auc_null_sd", "p_auc", "rho_obs", "p_rho"]
    count_fields = ["n_usable", "n_selected", "n_background", "n_cases", "K", "n_genes_tested"]
    hash_fields = ["sha256_discovery_results", "sha256_residual_analysis",
                   "sha256_split_replication_power"]
    mism = []
    for f in float_fields:
        got, want = computed[f], published[f]
        if not (np.isfinite(got) and np.isfinite(want) and abs(got - want) <= tol):
            mism.append((f, got, want))
    for f in count_fields + hash_fields:
        if computed[f] != published[f]:
            mism.append((f, computed[f], published[f]))
    if mism:
        log("\nFATAL: reproduced statistics do NOT match the published result -- stopping "
            f"before any plex analysis. ({published_path})")
        log(f"{'field':<28}{'this run':>24}{'published':>24}")
        for f, got, want in mism:
            log(f"{f:<28}{got!r:>24}{want!r:>24}")
        sys.exit(1)
    log(f"VALIDATED: all {len(float_fields) + len(count_fields) + len(hash_fields)} "
        f"statistics and code/discovery hashes match {published_path} (abs tol {tol:g} for "
        f"floats, exact for counts and hashes; file sha256-pinned to "
        f"{PUBLISHED_RESULT_SHA256[:12]}...).")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1], allow_abbrev=False)
    ap.add_argument("--perms", type=int, default=FROZEN_B)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--tag", default="")
    ap.add_argument("--blind", action="store_true",
                    help="stop after prepare() and the plex-label load; counts only")
    # defaults = the paths c1_ucec_stratified_reread.py / lsf_c1_ucec_d4_strat.sh used
    ap.add_argument("--protein", default="/public/home/fjhui/ZW/ucec_c1/protein_ucec_c1.csv")
    ap.add_argument("--discovery-results",
                    default="/public/home/fjhui/ZW/scripts/pinned/ucec/residual_results_tumoronly.csv")
    ap.add_argument("--discovery-slide-map",
                    default="/public/home/fjhui/ZW/scripts/pinned/slide_type_map_ucec.tsv")
    ap.add_argument("--validate-against",
                    default="/public/home/fjhui/ZW/ucec_c1/results/c1_ucec_result.json")
    ap.add_argument("--d4-npz",
                    default="/public/home/fjhui/ZW/ucec_c1/results/c1_ucec_d4_stratified.npz",
                    help="optional cross-check against the D4 matrices (used if present)")
    ap.add_argument("--case-ids",
                    default="/public/home/fjhui/ZW/scripts/pinned/c1_case_ids_ucec.txt")
    ap.add_argument("--plex-map",
                    default="/public/home/fjhui/ZW/scripts/pinned/case_to_plex_ucec_c1.tsv")
    args = ap.parse_args(argv)
    if not (1 <= args.perms <= FROZEN_B):
        ap.error(f"--perms must be in 1..{FROZEN_B} ({FROZEN_B} = the analysis; fewer = smoke)")
    if not _TAG.match(args.tag):
        ap.error("--tag may only contain letters, digits, '_' and '-' (max 40)")
    return args


def output_paths(out_dir, tag):
    suffix = f"_{tag}" if tag else ""
    return dict(
        result=os.path.join(out_dir, f"c1_ucec_plex_reread_result{suffix}.json"),
        log=os.path.join(out_dir, f"c1_ucec_plex_reread{suffix}.log"),
        npz=os.path.join(out_dir, f"c1_ucec_plex_reread{suffix}.npz"),
    )


def run(args, log):
    t0 = time.time()
    real = (args.perms == FROZEN_B) and not args.blind
    env = env_guard()
    out = output_paths(env["MORPHO_OUT"], args.tag)
    if real:
        existing = [p for p in (out["result"], out["npz"]) if os.path.exists(p)]
        if existing:
            fatal(f"refusing to overwrite {existing}; use a new --tag.")
        os.makedirs(env["MORPHO_OUT"], exist_ok=True)
        log.open(out["log"])

    mode = "BLIND" if args.blind else (f"ANALYSIS (B={FROZEN_B})" if real
                                       else f"SMOKE (B={args.perms})")
    log("=" * 78)
    log("C1-UCEC PLEX-STRATIFIED RE-READ -- POST HOC (new attempt)")
    log(f"  governed by {POST_HOC_GOVERNED_BY}; does not change the published C1-UCEC "
        f"reading or the D4 re-read")
    log(f"  mode={mode} host={socket.gethostname()} lsf_job={os.environ.get('LSB_JOBID', '')} "
        f"start={_utc()} python={platform.python_version()}")
    log("=" * 78)

    log("pre-flight hashes:")
    hashes = {}
    hashes["plex_map"] = assert_sha(args.plex_map, PLEX_MAP_SHA256, "plex map", log)
    hashes["case_ids"] = assert_sha(args.case_ids, FROZEN_INPUT_SHA256["case_ids"],
                                    "case ids", log)
    hashes["discovery_results"] = assert_sha(
        args.discovery_results, FROZEN_INPUT_SHA256["discovery_results"],
        "discovery results", log)
    hashes["published_result"] = assert_sha(args.validate_against, PUBLISHED_RESULT_SHA256,
                                            "published c1_ucec_result.json", log)
    hashes["protein_csv_recorded_not_asserted"] = sha256(args.protein)
    hashes["rna_manifest_recorded_not_asserted"] = sha256(env["MORPHO_RNA_MANIFEST"])
    hashes["confirmatory_slide_map_recorded_not_asserted"] = sha256(env["MORPHO_SLIDE_MAP"])
    hashes["this_script"] = sha256(os.path.abspath(__file__))
    log(f"  recorded (not asserted; the reproduction gate is their integrity check):"
        f" protein {hashes['protein_csv_recorded_not_asserted'][:12]}..., rna manifest "
        f"{hashes['rna_manifest_recorded_not_asserted'][:12]}..., slide map "
        f"{hashes['confirmatory_slide_map_recorded_not_asserted'][:12]}...")
    log(f"  this script sha256 = {hashes['this_script']}")

    C1, RA, SP, code_hashes = import_frozen(env, log)
    from sklearn.decomposition import PCA  # noqa: E402
    workers = args.workers if args.workers else RA.N_WORKERS

    log("=== resolved paths ===")
    for k in REQUIRED_ENV:
        log(f"  {k:<20}: {env[k]}")
    log(f"  plex map            : {args.plex_map}")
    log(f"  protein             : {args.protein}")
    log(f"  discovery results   : {args.discovery_results}")
    log(f"  validate against    : {args.validate_against}")
    log(f"  perms / workers     : {args.perms} / {workers}  blind={args.blind}")
    log("=" * 78)

    # --- frozen prepare() ---------------------------------------------------------
    common, tested, selected, disc_incr, disc_fdr, R, P, W, K = C1.prepare(args, log)
    n = len(common)
    for key, got in (("n_cases", n), ("K", K), ("n_genes_tested", len(tested))):
        if got != GATE_EXPECTED[key]:
            fatal(f"reproduction gate: {key}={got}, expected {GATE_EXPECTED[key]}.")
    with open(args.case_ids, encoding="utf-8") as f:   # hash-asserted above
        frozen_cases = sorted({line.strip() for line in f if line.strip()})
    if common != frozen_cases:
        fatal("the analysed case list differs from the frozen case-id file.")
    log(f"  gate n / K / genes tested: {n} / {K} / {len(tested)} -> OK; analysed cases == "
        f"frozen 138-case file")

    plex_map = read_plex_map(args.plex_map)
    raw = labels_for(common, plex_map)
    diag = strata_diagnostics(raw)
    dummy_prev, dropped_prev, left_prev = plex_dummy_levels(raw, MIN_GROUP)
    p3_present = sorted(P3_DROP_CASES & set(common))
    log(f"plex labels: {len(plex_map)} rows in map; analysed {n - diag['n_unlabelled']}/{n} "
        f"labelled ({diag['n_unlabelled']} unlabelled)")
    log(f"  {diag['n_levels']} plex levels, sizes {diag['smallest_stratum']}-"
        f"{diag['largest_stratum']}, {diag['n_strata_ge3']} strata >= 3 cases, "
        f"{diag['n_singleton_strata']} singleton strata; cases held fixed "
        f"{diag['n_cases_fixed']}/{n} = {diag['frac_cases_fixed']:.1%}")
    log(f"  P2 design: {len(dummy_prev)} dummies (levels >= {MIN_GROUP} cases; one dropped "
        f"for full rank: {dropped_prev}), {len(left_prev)} levels left to the intercept")
    log(f"  P3: {len(p3_present)}/{len(P3_DROP_CASES)} listed cases analysed {p3_present}; "
        + ", ".join(f"{c} -> {plex_map.get(c, '(no row)')}" for c in sorted(P3_DROP_CASES)))
    if diag["n_unlabelled"]:
        log(f"  WARNING: {diag['n_unlabelled']} analysed cases have no plex label; they are "
            f"held fixed and get no dummy: {[c for c, v in zip(common, raw) if not v]}")

    if args.blind:
        log(f"\nBLIND: n={n} K={K} n_tested={len(tested)} n_selected(usable by universe)="
            f"{len(selected)} plex_levels={diag['n_levels']} dummies={len(dummy_prev)} "
            f"p3_cases_present={len(p3_present)}")
        log("Nothing further is computed, printed or written in --blind mode.")
        return 0

    # --- step 1: reproduce the published test exactly (c1_run_test.run_test up to M) -----
    log("\n=== STEP 1 REPRODUCTION (published 2,566 set, unrestricted null, published seed) ===")
    pcs = PCA(n_components=K, svd_solver="full").fit_transform(W)
    rng = np.random.RandomState(zlib.crc32(SEED_UNRESTRICTED))
    perms_u = [rng.permutation(n) for _ in range(args.perms)]
    step1 = family(C1, SP, tested, selected, disc_incr, R, P, pcs, perms_u, workers, log,
                   "step1")
    s1 = stats_of(step1)
    computed = dict(
        auc_obs=s1["auc_obs"], auc_null_mean=s1["auc_null_mean"], auc_null_sd=s1["auc_null_sd"],
        p_auc=s1["p_auc"], rho_obs=s1["rho_obs"], p_rho=s1["p_rho"],
        n_usable=s1["n_usable"], n_selected=s1["n_selected"], n_background=s1["n_background"],
        n_cases=n, K=K, n_genes_tested=len(tested),
        sha256_discovery_results=hashes["discovery_results"],
        sha256_residual_analysis=code_hashes["residual_analysis"],
        sha256_split_replication_power=code_hashes["split_replication_power"])
    if args.perms == FROZEN_B:
        validate(log, computed, args.validate_against)
    else:
        with open(args.validate_against, "r", encoding="utf-8") as f:
            published = json.load(f)
        bad = [k for k in ("n_usable", "n_selected", "n_background") if computed[k] != published[k]]
        if abs(computed["auc_obs"] - published["auc_obs"]) > FLOAT_TOL or bad:
            fatal(f"smoke: auc_obs {computed['auc_obs']!r} vs published "
                  f"{published['auc_obs']!r}, count mismatches {bad}.")
        log(f"SMOKE VALIDATED: auc_obs and counts match the published result (abs tol "
            f"{FLOAT_TOL:g}). B={args.perms}: smallest attainable p is "
            f"{1.0 / (args.perms + 1):.4f}; the null statistics are not compared.")

    d4 = dict(checked=False)
    if os.path.exists(args.d4_npz):
        got = sha256(args.d4_npz)
        if got != D4_NPZ_SHA256:
            log(f"  NOTE: {args.d4_npz} has sha256 {got[:12]}..., not the pinned "
                f"{D4_NPZ_SHA256[:12]}...; D4 cross-check skipped.")
        else:
            z = np.load(args.d4_npz, allow_pickle=True)
            if list(z["genes"]) != step1["genes"]:
                fatal("gene order differs from the D4 matrix file.")
            dobs = float(np.max(np.abs(z["M_unrestricted"][:, 0] - step1["M"][:, 0])))
            if dobs > FLOAT_TOL:
                fatal(f"observed column differs from D4's M_unrestricted by {dobs!r}.")
            ncol = min(z["M_unrestricted"].shape[1], step1["M"].shape[1])
            dnull = float(np.max(np.abs(z["M_unrestricted"][:, :ncol] - step1["M"][:, :ncol])))
            d4 = dict(checked=True, max_abs_diff_observed=dobs,
                      max_abs_diff_shared_columns=dnull, shared_columns=int(ncol))
            log(f"  D4 cross-check: gene order identical; observed column max |diff| = "
                f"{dobs:.2e}; shared columns ({ncol}) max |diff| = {dnull:.2e}")
    else:
        log(f"  D4 cross-check skipped ({args.d4_npz} not present).")

    if not real:
        log(f"\nSMOKE: B={args.perms}. The step-1 reproduction was checked on this host; no "
            f"within-plex null is computed and nothing is written.")
        log(f"elapsed {round((time.time() - t0) / 60.0, 1)} min")
        return 0

    # --- P1 --------------------------------------------------------------------------------
    log("\n=== P1: published set, same PCs, within-plex null (post hoc) ===")
    perms_plex = build_within_strata_perms(raw, args.perms, SEED_PLEX)
    p1 = family(C1, SP, tested, selected, disc_incr, R, P, pcs, perms_plex, workers, log, "P1")
    sp1 = stats_of(p1)
    for key in ("auc_obs", "rho_obs"):
        if abs(sp1[key] - s1[key]) > FLOAT_TOL:
            fatal(f"P1 observed {key} {sp1[key]!r} != step 1 {s1[key]!r}: they must be "
                  f"identical (same PCs, same genes; only the null differs).")
    if p1["genes"] != step1["genes"]:
        fatal("P1 usable gene list differs from step 1.")

    # --- P2 --------------------------------------------------------------------------------
    log("\n=== P2 (descriptive): published set, PCs residualised on plex dummies, same "
        "within-plex draws ===")
    pcs_resid, design = build_plex_dummy_design(raw, pcs, MIN_GROUP)
    log(f"  plex residualisation: {design['n_dummies']} dummies (dropped for full rank: "
        f"{design['dropped_for_full_rank']}; left to intercept: "
        f"{design['levels_left_to_intercept']}), PC variance retained "
        f"{design['pc_variance_retained']:.1%} (n={n})")
    p2 = family(C1, SP, tested, selected, disc_incr, R, P, pcs_resid, perms_plex, workers,
                log, "P2")
    sp2 = stats_of(p2)
    sp2.update(design)
    if p2["genes"] != step1["genes"]:
        fatal("P2 usable gene list differs from step 1.")

    # --- P3 --------------------------------------------------------------------------------
    log(f"\n=== P3 (sensitivity): P1 without {sorted(P3_DROP_CASES)} (two aliquots, two "
        f"plexes) ===")
    if p3_present:
        keep = np.array([c not in P3_DROP_CASES for c in common])
        common3 = [c for c, k in zip(common, keep) if k]
        n3 = len(common3)
        K3 = SP.k_for(n3, C1.PCS_PER_100)
        R3, P3, W3 = R[keep], P[keep], W[keep]
        raw3 = labels_for(common3, plex_map)
        diag3 = strata_diagnostics(raw3)
        log(f"  n = {n3}, K = {K3}, plex sizes {diag3['smallest_stratum']}-"
            f"{diag3['largest_stratum']}, cases held fixed {diag3['n_cases_fixed']}")
        pcs3 = PCA(n_components=K3, svd_solver="full").fit_transform(W3)
        perms3 = build_within_strata_perms(raw3, args.perms, SEED_PLEX)
        p3 = family(C1, SP, tested, selected, disc_incr, R3, P3, pcs3, perms3, workers, log,
                    "P3")
        sp3 = stats_of(p3)
        sp3["n_cases"], sp3["K"], sp3["dropped_cases"] = n3, K3, p3_present
        sp3["strata"] = diag3
        if p3["genes"] != step1["genes"]:
            log("  NOTE: P3 usable gene list differs from step 1 "
                f"({len(p3['genes'])} vs {len(step1['genes'])}).")
    else:
        log("  neither case is in the analysed set; P3 would equal P1 and is not run.")
        p3, sp3 = None, None

    # --- summaries -------------------------------------------------------------------------
    sp1["frac_cases_fixed"] = diag["frac_cases_fixed"]
    sp1["null_sd_ratio_vs_unrestricted"] = (
        sp1["auc_null_sd"] / s1["auc_null_sd"] if s1["auc_null_sd"] > 0 else float("nan"))
    sp1["null_mean_shift_vs_unrestricted"] = sp1["auc_null_mean"] - s1["auc_null_mean"]

    log("\n" + "=" * 78)
    log("POST HOC SUMMARY (does not change the published C1-UCEC reading)")
    log("=" * 78)
    log(f"step1 (reproduced): AUC={s1['auc_obs']:.6f} null {s1['auc_null_mean']:.4f}"
        f"+/-{s1['auc_null_sd']:.4f} p={s1['p_auc']:.6f}")
    for name, d in (("P1", sp1), ("P2", sp2), ("P3", sp3)):
        if d is None:
            log(f"{name}: not run")
            continue
        log(f"{name}: AUC={d['auc_obs']:.6f} null {d['auc_null_mean']:.4f}+/-"
            f"{d['auc_null_sd']:.4f} z={d['z_auc']:.2f} p={d['p_auc']:.6f} "
            f"({d['exceedances_auc']}/{d['perms']}) | rho={d['rho_obs']:.4f} null "
            f"{d['rho_null_mean']:.4f}+/-{d['rho_null_sd']:.4f} z={d['z_rho']:.2f} "
            f"p_rho={d['p_rho']:.6f} | n_sel={d['n_selected']} n_bg={d['n_background']}")
    log(f"P1 null SD ratio vs unrestricted: {sp1['null_sd_ratio_vs_unrestricted']:.3f}; "
        f"null mean shift {sp1['null_mean_shift_vs_unrestricted']:+.4f}; P2 dummies "
        f"{sp2['n_dummies']}, PC variance retained {sp2['pc_variance_retained']:.1%}")

    minutes = round((time.time() - t0) / 60.0, 1)
    with open(out["npz"], "xb") as fh:
        np.savez_compressed(
            fh, common=np.array(common), plex_raw=raw, genes=np.array(step1["genes"]),
            sel_mask=step1["sel_mask"], M_unrestricted=step1["M"], M_within_plex=p1["M"],
            M_within_plex_residualized=p2["M"],
            genes_p3=np.array(p3["genes"] if p3 else []),
            sel_mask_p3=(p3["sel_mask"] if p3 else np.array([], dtype=bool)),
            M_p3=(p3["M"] if p3 else np.empty((0, 0))))
    result = dict(
        post_hoc=True, governed_by=POST_HOC_GOVERNED_BY,
        published_reading_unchanged=True,
        note=("Post hoc plex-stratified re-read of the published C1-UCEC test. Nothing here "
              "changes the published section 5 reading or the D4 within-operator re-read."),
        cohort="ucec_c1", host=socket.gethostname(), lsf_job=os.environ.get("LSB_JOBID", ""),
        finished_utc=_utc(), minutes=minutes, perms=args.perms, workers=workers,
        seeds=dict(step1=SEED_UNRESTRICTED.decode(), p1_p2_p3=SEED_PLEX.decode()),
        n_cases=n, K=K, n_genes_tested=len(tested),
        reproduction_gate=dict(validated_against=args.validate_against,
                               published_sha256=PUBLISHED_RESULT_SHA256, d4_crosscheck=d4),
        step1=s1, plex_strata=diag, P1=sp1, P2=sp2, P3=sp3,
        paths=dict(plex_map=args.plex_map, protein=args.protein, case_ids=args.case_ids,
                   discovery_results=args.discovery_results,
                   discovery_slide_map=args.discovery_slide_map,
                   confirmatory_rna_manifest=env["MORPHO_RNA_MANIFEST"],
                   confirmatory_slide_map=env["MORPHO_SLIDE_MAP"],
                   wsi_emb_dir=env["MORPHO_WSI_EMB_DIR"], out=env["MORPHO_OUT"],
                   log=out["log"], npz=out["npz"]),
        sha256_inputs=hashes,
        sha256_code=dict(this_script=hashes["this_script"], **code_hashes),
        sha256_npz=sha256(out["npz"]),
        versions=dict(python=platform.python_version(), numpy=np.__version__,
                      pandas=pd.__version__, sklearn=__import__("sklearn").__version__),
    )
    with open(out["result"], "x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=2,
                  default=lambda o: (o.tolist() if hasattr(o, "tolist") else str(o)))
    log(f"\n  -> {out['result']}\n  -> {out['npz']}\n  -> {out['log']}\n  elapsed {minutes} min")
    return 0


def main(argv=None):
    args = parse_args(argv)
    log = Tee()
    try:
        return run(args, log)
    except SystemExit as e:
        if isinstance(e.code, str):
            log(e.code)
        raise
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
