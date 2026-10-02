#!/usr/bin/env python3
"""
c1_luad_plex_reread.py -- POST HOC TMT-plex-stratified re-read of the registered
C1-LUAD test.

STATUS AND AUTHORITY
  POST HOC. Governed by review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md, which every
  run pins by --prespec / --prespec-sha256. This is NOT part of the signed C1-LUAD rule
  (review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md, sha256 a68dfbec...). That rule has
  no plex stratum: its section 6 said none was possible. review/C1_LUAD_ADDENDA.md
  Entry 7 records that the claim was wrong, and that any plex-stratified re-read is
  post hoc, is specified in a dated file before it runs, and is reported whatever its
  result. Nothing here changes the registered section 5 reading (p1 and clauses
  [B], [C], [D]). Every statistic below is descriptive and post hoc.

WHAT IT DOES
  The frozen machinery of c1_run_test_luad.py (sha256 1d5ee05d...) is IMPORTED, not
  reimplemented: prepare(), run_family() (which calls assemble() and SP.run_pool),
  strip_arrays(), load_operator_labels(), EXPECTED_SHA256 and MIN_GROUP.
  residual_analysis (RA) and split_replication_power (SP) come in through it. The new
  code is only (a) the within-plex permutation null, the operator loop of
  build_within_operator_perms with a new seed, and (b) the plex dummy design.

  NAMING. P1, P2 and P3 below are the prespec's L1, L2 and L3 (its section L). They are
  not its section P, which is the discovery plex axis. The log and the result JSON
  carry this mapping.

  0. Pre-flight, before any data is read. The MORPHO_* environment must be set before
     python starts; residual_analysis.py reads it at import. MORPHO_ROOT must be the
     luad_c1 tree, with MORPHO_OUT inside it. The sha256 of every input file and of
     the three frozen modules is asserted; any mismatch is FATAL.
  1. Reproduction gate (FATAL on any mismatch). Step 1 of the primary run is recomputed
     exactly: prepare(); PCA(K, svd_solver="full") on the confirmatory embeddings; B
     unrestricted permutations from RandomState(zlib.crc32(b"c1_luad|perm")); then
     run_family on the 2,310-gene H1 set. Required: n == 112, K == 22,
     n_usable == 10406, n_selected == 2287, |auc_obs - 0.6439668501223573| <= 1e-9,
     and p1 == 1/1001 (checked only at B = 1000).
  2. P1. The same 2,310 set and the same PCs, under a within-plex null: B draws from
     RandomState(zlib.crc32(b"c1_luad|perm_plex")), strata in sorted label order,
     singleton plexes and unlabelled cases held fixed. The observed AUC and rho are
     asserted identical to step 1, since only the null changes. Reported for AUC and
     rho: null mean, SD, z and one-sided p; also the fraction of cases held fixed.
  3. P2. The 726-gene H2 set, exactly as in step 3 of the primary run (selection,
     background and incremental_r2_batchresid discovery increments). Its PCs are the
     confirmatory PCs residualised on plex dummies, and its null is the SAME
     within-plex draws as P1. Dummies: every plex with >= MIN_GROUP (3) analysed cases,
     in Counter.most_common order. If the dummies cover every analysed case, the last
     is dropped for full rank. This is residual_analysis_sitepack.py's rule with no
     MAX_LEVELS cap. The residuals do not depend on which level is dropped, because the
     column space is unchanged. The fit is least squares with an intercept, as
     build_operator_dummy_design does for operators. Reported: AUC, null mean/SD, z, p,
     rho, p_rho, PC variance retained, number of dummies.
  4. P3. Sensitivity: P1 repeated without C3N-02230 and C3N-02240, each of which has one
     aliquot measured in two plexes. It uses prepare(drop_cases=...), the frozen path
     that (s)/(s2) used. K is recomputed and the PCA refit. The within-plex null comes
     from a fresh RandomState(zlib.crc32(b"c1_luad|perm_plex")), so P1's procedure is
     repeated literally, seed included.

INPUTS (all required CLI arguments, no defaults, every file sha256-asserted)
  --rule / --rule-sha256       signed C1-LUAD rule; the hash must also equal the signed
                               value a68dfbec... below
  --prespec / --prespec-sha256 the signed post hoc prespec governing this script
  --case-ids                   pinned/c1_case_ids_luad.txt                (constant)
  --protein / --protein-sha256 luad_c1/protein_luad_c1.csv (no hash was frozen in 7.2,
                               so it is passed in; the gate cross-checks the content)
  --discovery-results          pinned/luad/residual_results_tumoronly.csv (constant)
  --discovery-sitepack         results/luad__sitepack_operator.csv        (constant)
  --discovery-slide-map        pinned/slide_type_map_luad.tsv             (constant)
  --discovery-xwalk            pinned/aliquot_to_case_tumor_luad.tsv      (constant)
  --discovery-rna-manifest     pinned/manifest_rna_tumor_luad.tsv         (constant)
  --batch-labels               pinned/c1_luad_batch_labels.tsv            (constant;
                               hash-asserted and coverage logged, not otherwise used)
  --plex-map / --plex-map-sha256  pinned/case_to_plex_luad_c1.tsv (module A); the
                               first two columns are case_id (or case_submitter_id)
                               and plex
  --rna-manifest-sha256        hash of $MORPHO_RNA_MANIFEST (the luad_c1 RNA manifest)
  $MORPHO_SLIDE_MAP            pinned/slide_type_map_luad_c1.tsv          (constant)
  Per-aliquot RNA files under $MORPHO_ROOT/omics/rna and the .pt embeddings under
  $MORPHO_WSI_EMB_DIR are not hashed one by one. The reproduction gate (auc_obs to
  1e-9 plus four exact counts) is their integrity check.

OUTPUTS (in $MORPHO_OUT = /public/home/fjhui/ZW/luad_c1/results; never overwritten.
The JSON and npz are opened with mode "x", and the run refuses to start if either
exists. The log is append-only)
  c1_luad_plex_reread_result[_<tag>].json  every statistic above, the gate table, the
                                           sha256 of every input, of this script and
                                           of each imported frozen module
  c1_luad_plex_reread[_<tag>].log          a copy of everything printed. Each attempt
                                           appends under its own header, so a failed
                                           attempt stays on record
  c1_luad_plex_reread[_<tag>].npz          genes x (1+B) matrices for P1, P2 and P3,
                                           for follow-ups that need no new permutation
  None of these names collides with the primary run's c1_luad_primary_* files.

MODES
  --perms 1000 (default)  the post hoc analysis; writes the three outputs.
  --perms N, 1 <= N < 1000  smoke run. Checks the gate fields except p1 on this host,
                          then stops. No within-plex null is computed and nothing is
                          written.
  --blind                 stops after prepare() and the plex-label load. Prints counts
                          only (n, K, genes, plex levels and sizes, dummies, P3 cases
                          present). Nothing is written.

RUN: through lsf_c1_luad_plex_reread.sh, which sets the MORPHO_* environment. Pin the
host and core count of the primary run (job 75682618: s002, 4 cores); see the wrapper
header for why. Smoke first, then the full run:
  bsub -q smp -n 4 -m s002 -R "span[hosts=1]" -o c1_luad_plex_smoke.out \
    lsf_c1_luad_plex_reread.sh --perms 5 [... required paths/hashes]
  bsub -q smp -n 4 -m s002 -R "span[hosts=1]" -o c1_luad_plex_reread.out \
    lsf_c1_luad_plex_reread.sh --perms 1000 [... required paths/hashes]
Expected wall time at B = 1000 on 4 cores: about 5 h, i.e. four run_family passes
(step 1, P1, P2, P3). The primary run's four passes took 4.7 h.

No openslide/torch/timm import happens here directly. residual_analysis imports torch;
openslide is never imported in this process.
"""
import argparse
import hashlib
import importlib.util
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
# Same search path as c1_run_test_luad.py: this directory first (HPC layout, where
# split_replication_power.py sits beside the scripts), then the repository top level
# (local layout: server_export/scripts -> ../..). Each module's hash is taken from the
# file that the import will actually load (importlib.util.find_spec), not from an
# assumed location.
sys.path.insert(0, _SCRIPT_DIR)
sys.path.insert(1, os.path.dirname(os.path.dirname(_SCRIPT_DIR)))

POST_HOC_GOVERNED_BY = "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md"

# ---------------------------------------------------------------------------
# Pinned constants. Changing any of these changes the analysis.
# ---------------------------------------------------------------------------
SIGNED_RULE_SHA256 = "a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022"

# Frozen code imported by this script (rule section 7.2).
FROZEN_CODE_SHA256 = {
    "c1_run_test_luad": "1d5ee05d57f80af4d90870aecb10e574e4c910a8ce1eab743a067bd9a9278742",
    "residual_analysis": "a54a0a56494d315e15ef075425fe2133ecd9f05193c22656c1145bffdaac4a0a",
    "split_replication_power": "79530c075c0f7e75bf06d3c2f1edfe8792cbf4d6c71ae1d32e60bf44f31a0397",
}

# Frozen inputs (rule section 7.2). The first five duplicate c1_run_test_luad.py's
# EXPECTED_SHA256. They are checked against it after import, and asserted here BEFORE
# import with no skip switch (c1_run_test_luad.assert_hash allows MORPHO_ALLOW_HPC_*
# overrides for two of them; this script does not).
FROZEN_INPUT_SHA256 = {
    "case_ids": "3629b243a3aa21057dae24e7f5c1b47e8a2571e6ffbcff6f3170333e1ca9fe44",
    "discovery_results": "63aec1992a77e952b96c4e25e27372e2bee5c52294560c56ef2b229196e87be7",
    "discovery_sitepack": "aeef7232682c24d209024ba73efea8033cb7486779831bda1a1ba3c8bda45d9f",
    "batch_labels": "0a840533a3601490c360cae8e5fcfc5e6f426e5314e5ac02d1e1b00730f1b3d4",
    "discovery_slide_map": "3de16eb77f615908c1999806388d9871b48ed8bfc645a6366a7f64b7df9a843e",
    "discovery_xwalk": "967cabf18cb6b6c3646aaf05005c46fcdccd7b55eb8dc9882348e84ca4baa2bf",
    "discovery_rna_manifest": "cacb2858814584a555991829b522807211f0cfe3204c3fc2f37b72e60815e82b",
    "confirmatory_slide_map": "d564b6bc5682bc6332aa3ad82387ed300a87d8142bba2d794650281350347d60",
}

FROZEN_B = 1000
SEED_UNRESTRICTED = b"c1_luad|perm"     # primary step 1 (reproduced)
SEED_PLEX = b"c1_luad|perm_plex"        # P1 and P2 (shared draws)
SEED_P3 = b"c1_luad|perm_plex"          # P3: P1 repeated literally, fresh RandomState

P3_DROP_CASES = frozenset({"C3N-02230", "C3N-02240"})

# P1/P2/P3 here are the prespec's section-L analyses, not its section P (discovery plex).
PRESPEC_LABELS = {"P1": "L1", "P2": "L2", "P3": "L3"}

# Reproduction gate: the published primary run (job 75682618, 2026-10-02).
GATE_EXPECTED = {
    "n": 112,
    "K": 22,
    "n_usable": 10406,
    "n_selected": 2287,
    "auc_obs": 0.6439668501223573,
    "p1": 1.0 / 1001,
}
GATE_FLOAT_TOL = 1e-9
OBS_IDENTITY_TOL = 1e-9     # same tolerance c1_run_test_luad.py uses for (a)(i) vs H1

REQUIRED_ENV = ("MORPHO_ROOT", "MORPHO_RNA_MANIFEST", "MORPHO_SLIDE_MAP",
                "MORPHO_WSI_EMB_DIR", "MORPHO_OUT")

_HEX64 = re.compile(r"^[0-9a-f]{64}$")
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
    """print() plus, once open() is called, a line-buffered copy to a log file. The log
    is opened in APPEND mode, so it is never truncated: a failed attempt (pre-flight or
    gate FATAL) stays on record, and a later attempt under the same tag appends below
    it, after its own header line."""

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


# ---------------------------------------------------------------------------
# Plex labels, within-plex null, plex dummy design (the only new statistical code)
# ---------------------------------------------------------------------------
def read_plex_map(path):
    """case -> plex label. First two columns must be case_id (or case_submitter_id) and
    plex. Duplicate case rows and malformed case IDs are FATAL. An empty plex cell
    means unlabelled."""
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
        fatal(f"{path}: {len(dup)} case IDs appear more than once (ambiguous plex), "
              f"e.g. {dup[:5]}")
    return dict(zip(cases, plex))


def labels_for(cases, mapping):
    """np.ndarray of plex labels aligned to `cases`; '' for a case with no row or an
    empty label (held fixed in the within-plex null, never a level)."""
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
    level and never move; a stratum of size 1 never moves. Strata are processed in
    sorted label order within every draw, from ONE RandomState(zlib.crc32(seed_bytes)),
    which fixes the RNG draw sequence. This is the loop of
    c1_run_test_luad.build_within_operator_perms with the seed as an argument; it
    reproduces that function's draws exactly when given b"c1_luad|perm_op"."""
    raw = np.asarray(raw)
    n = len(raw)
    levels = sorted(set(v for v in raw if v))
    members = {lv: np.where(raw == lv)[0] for lv in levels}
    rng = np.random.RandomState(zlib.crc32(seed_bytes))
    perms = []
    for _ in range(B):
        p = np.arange(n)
        for lv in levels:  # sorted already
            m = members[lv]
            if len(m) > 1:
                p[m] = rng.permutation(m)
        perms.append(p)
    return perms


def plex_dummy_levels(raw, min_group):
    """residual_analysis_sitepack.py's dummy rule without the MAX_LEVELS cap.
    Levels with >= min_group analysed cases, in Counter.most_common order: count
    descending, ties by first appearance in `raw`, which is sorted case order. If
    those dummies cover every analysed case (no unlabelled case and no small level
    left to the intercept), the LAST one is dropped so the design is full rank.
    Unlabelled cases are never a level (departure (i), as for operators)."""
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
    """Residualise pcs on [1, plex dummies] by least squares (np.linalg.lstsq), as
    c1_run_test_luad.build_operator_dummy_design does for operators. FATAL unless the
    design has full column rank."""
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
# Pre-flight: environment, frozen code, inputs
# ---------------------------------------------------------------------------
def env_guard():
    env = {}
    for k in REQUIRED_ENV:
        v = os.environ.get(k, "")
        if not v:
            fatal(f"{k} is not set. Run through lsf_c1_luad_plex_reread.sh, which exports "
                  f"every MORPHO_* before python starts (residual_analysis.py reads them "
                  f"at import time).")
        env[k] = v
    root = os.path.normpath(env["MORPHO_ROOT"])
    out = os.path.normpath(env["MORPHO_OUT"])
    if os.path.basename(root) != "luad_c1":
        fatal(f"basename(MORPHO_ROOT)={os.path.basename(root)!r}, expected 'luad_c1'.")
    if not (out == root or out.startswith(root + os.sep)):
        fatal(f"MORPHO_OUT={out!r} does not lie inside the luad_c1 tree "
              f"(MORPHO_ROOT={root!r}).")
    return env


def locate_frozen_modules():
    """Origin file of each frozen module that the import below would load, found
    without executing it."""
    found = {}
    for name in FROZEN_CODE_SHA256:
        spec = importlib.util.find_spec(name)
        if spec is None or not spec.origin:
            fatal(f"cannot locate module {name} on sys.path {sys.path[:3]}")
        found[name] = os.path.realpath(spec.origin)
    return found


def import_frozen(env, log):
    origins = locate_frozen_modules()
    code_hashes = {}
    for name, origin in origins.items():
        code_hashes[name] = assert_sha(origin, FROZEN_CODE_SHA256[name],
                                       f"frozen module {name}", log)
    import c1_run_test_luad as C1L  # noqa: E402  (imports residual_analysis, split_replication_power)

    for name, mod in (("c1_run_test_luad", C1L), ("residual_analysis", C1L.RA),
                      ("split_replication_power", C1L.SP)):
        if os.path.realpath(mod.__file__) != origins[name]:
            fatal(f"imported {name} from {mod.__file__}, but hashed {origins[name]}.")
    for key in ("case_ids", "discovery_results", "discovery_sitepack", "batch_labels",
                "discovery_slide_map", "discovery_xwalk"):
        if C1L.EXPECTED_SHA256[key] != FROZEN_INPUT_SHA256[key]:
            fatal(f"pinned hash for {key} disagrees with c1_run_test_luad.EXPECTED_SHA256.")
    # residual_analysis read the environment at import. Confirm that it saw the same
    # environment that env_guard() checked (also true if it was imported earlier).
    RA = C1L.RA
    for attr, key in (("ROOT", "MORPHO_ROOT"), ("RNA_MANIFEST", "MORPHO_RNA_MANIFEST"),
                      ("SLIDE_TYPE_MAP", "MORPHO_SLIDE_MAP"),
                      ("WSI_EMB_DIR", "MORPHO_WSI_EMB_DIR"), ("OUT_DIR", "MORPHO_OUT")):
        if os.path.normpath(str(getattr(RA, attr))) != os.path.normpath(env[key]):
            fatal(f"residual_analysis.{attr}={getattr(RA, attr)} != ${key}={env[key]}: "
                  f"the environment changed after residual_analysis was imported.")
    return C1L, code_hashes


def check_gate(repro, perms, log):
    """Compare the reproduced step-1 statistics with GATE_EXPECTED. Returns the table;
    exits FATAL if any checked field mismatches."""
    table = {}
    bad = []
    for key in ("n", "K", "n_usable", "n_selected", "auc_obs", "p1"):
        got, want = repro[key], GATE_EXPECTED[key]
        if key == "p1" and perms != FROZEN_B:
            table[key] = dict(reproduced=got, expected=want, checked=False, ok=None)
            log(f"  gate {key:<10}: reproduced={got!r} expected={want!r} -> not checked "
                f"(B={perms} != {FROZEN_B})")
            continue
        if key in ("auc_obs", "p1"):
            ok = bool(np.isfinite(got) and abs(got - want) <= GATE_FLOAT_TOL)
        else:
            ok = (int(got) == int(want))
        table[key] = dict(reproduced=got, expected=want, checked=True, ok=ok)
        log(f"  gate {key:<10}: reproduced={got!r} expected={want!r} -> "
            f"{'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(key)
    if bad:
        fatal(f"reproduction gate failed on {bad}. This run does not reproduce the "
              f"published C1-LUAD step 1, so no post hoc null is computed. If only "
              f"auc_obs differs, check the host first: the primary ran on s002 with 4 "
              f"cores, and LAPACK SVD differs across hosts.")
    log("  reproduction gate PASSED" + ("" if perms == FROZEN_B else
                                         " (smoke: p1 not checked)"))
    return table


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    # allow_abbrev=False: a flag must be spelled in full (no prefix matching).
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1], allow_abbrev=False)
    ap.add_argument("--perms", type=int, default=FROZEN_B)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--tag", default="")
    ap.add_argument("--blind", action="store_true",
                    help="stop after prepare() and the plex-label load; counts only; "
                         "nothing written")
    ap.add_argument("--rule", required=True)
    ap.add_argument("--rule-sha256", required=True)
    ap.add_argument("--prespec", required=True,
                    help="signed review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md")
    ap.add_argument("--prespec-sha256", required=True)
    ap.add_argument("--case-ids", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--protein-sha256", required=True)
    ap.add_argument("--rna-manifest-sha256", required=True,
                    help="sha256 of $MORPHO_RNA_MANIFEST (luad_c1 RNA manifest)")
    ap.add_argument("--discovery-results", required=True)
    ap.add_argument("--discovery-sitepack", required=True)
    ap.add_argument("--discovery-slide-map", required=True)
    ap.add_argument("--discovery-xwalk", required=True)
    ap.add_argument("--discovery-rna-manifest", required=True)
    ap.add_argument("--batch-labels", required=True)
    ap.add_argument("--plex-map", required=True, help="pinned/case_to_plex_luad_c1.tsv")
    ap.add_argument("--plex-map-sha256", required=True)
    args = ap.parse_args(argv)
    if not (1 <= args.perms <= FROZEN_B):
        ap.error(f"--perms must be in 1..{FROZEN_B} ({FROZEN_B} = the analysis; fewer = smoke)")
    if not _TAG.match(args.tag):
        ap.error("--tag may only contain letters, digits, '_' and '-' (max 40)")
    for k in ("rule_sha256", "prespec_sha256", "protein_sha256", "rna_manifest_sha256",
              "plex_map_sha256"):
        if not _HEX64.match(getattr(args, k)):
            ap.error(f"--{k.replace('_', '-')} must be 64 lowercase hex characters")
    return args


def output_paths(out_dir, tag):
    suffix = f"_{tag}" if tag else ""
    paths = dict(
        result=os.path.join(out_dir, f"c1_luad_plex_reread_result{suffix}.json"),
        log=os.path.join(out_dir, f"c1_luad_plex_reread{suffix}.log"),
        npz=os.path.join(out_dir, f"c1_luad_plex_reread{suffix}.npz"),
    )
    for p in paths.values():
        if os.path.basename(p).startswith("c1_luad_primary"):
            fatal(f"output name {p} would collide with the primary run's files.")
    return paths


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
    log("C1-LUAD PLEX-STRATIFIED RE-READ -- POST HOC (new attempt)")
    log(f"  governed by {POST_HOC_GOVERNED_BY}; does not change the registered "
        f"section 5 reading")
    log(f"  mode={mode} host={socket.gethostname()} lsf_job={os.environ.get('LSB_JOBID', '')} "
        f"start={_utc()} python={platform.python_version()}")
    log(f"  P1/P2/P3 in this output = prespec section L analyses "
        f"{PRESPEC_LABELS['P1']}/{PRESPEC_LABELS['P2']}/{PRESPEC_LABELS['P3']} "
        f"(not section P, the discovery plex axis)")
    log("=" * 78)

    # --- pre-flight: every input hash, before any data is read ---------------------
    log("pre-flight hashes:")
    hashes = {}
    if args.rule_sha256 != SIGNED_RULE_SHA256:
        fatal(f"--rule-sha256 {args.rule_sha256} is not the signed rule's "
              f"{SIGNED_RULE_SHA256}.")
    hashes["rule"] = assert_sha(args.rule, args.rule_sha256, "signed C1-LUAD rule", log)
    hashes["prespec"] = assert_sha(args.prespec, args.prespec_sha256, "post hoc prespec", log)
    hashes["plex_map"] = assert_sha(args.plex_map, args.plex_map_sha256, "plex map", log)
    hashes["protein"] = assert_sha(args.protein, args.protein_sha256,
                                   "confirmatory protein", log)
    hashes["confirmatory_rna_manifest"] = assert_sha(
        env["MORPHO_RNA_MANIFEST"], args.rna_manifest_sha256,
        "confirmatory RNA manifest ($MORPHO_RNA_MANIFEST)", log)
    hashes["confirmatory_slide_map"] = assert_sha(
        env["MORPHO_SLIDE_MAP"], FROZEN_INPUT_SHA256["confirmatory_slide_map"],
        "confirmatory slide map ($MORPHO_SLIDE_MAP)", log)
    for key, path in (("case_ids", args.case_ids),
                      ("discovery_results", args.discovery_results),
                      ("discovery_sitepack", args.discovery_sitepack),
                      ("discovery_slide_map", args.discovery_slide_map),
                      ("discovery_xwalk", args.discovery_xwalk),
                      ("discovery_rna_manifest", args.discovery_rna_manifest),
                      ("batch_labels", args.batch_labels)):
        hashes[key] = assert_sha(path, FROZEN_INPUT_SHA256[key], key, log)
    hashes["this_script"] = sha256(os.path.abspath(__file__))
    log(f"  this script sha256 = {hashes['this_script']}")

    C1L, code_hashes = import_frozen(env, log)
    from sklearn.decomposition import PCA  # noqa: E402
    workers = args.workers if args.workers else C1L.RA.N_WORKERS

    log("=== resolved paths ===")
    for k in REQUIRED_ENV:
        log(f"  {k:<20}: {env[k]}")
    log(f"  plex map            : {args.plex_map}")
    log(f"  protein             : {args.protein}")
    log(f"  perms / workers     : {args.perms} / {workers}  blind={args.blind}")
    log("=" * 78)

    # --- frozen prepare(): load, intersect, independence check, gene sets ----------
    common, tested, sel_h1, sel_h2, disc_incr, disc_fdr, disc_incr_bres, R, P, W, K = \
        C1L.prepare(args, log)
    n = len(common)
    for key, got in (("n", n), ("K", K)):
        if got != GATE_EXPECTED[key]:
            fatal(f"reproduction gate: {key}={got}, expected {GATE_EXPECTED[key]}.")
    log(f"  gate n / K: {n} / {K} -> OK")

    C1L.load_operator_labels(common, args.batch_labels, log)  # hash + coverage only
    plex_map = read_plex_map(args.plex_map)
    raw = labels_for(common, plex_map)
    diag = strata_diagnostics(raw)
    dummy_preview, dropped_preview, left_preview = plex_dummy_levels(raw, C1L.MIN_GROUP)
    with open(args.case_ids, encoding="utf-8") as f:   # hash-asserted above
        frozen_cases = {line.strip() for line in f if line.strip()}
    p3_present = sorted(P3_DROP_CASES & set(common))
    log(f"plex labels: {len(plex_map)} rows in map; {len(frozen_cases & set(plex_map))}/"
        f"{len(frozen_cases)} frozen cases covered; analysed {n - diag['n_unlabelled']}/{n} "
        f"labelled ({diag['n_unlabelled']} unlabelled)")
    log(f"  {diag['n_levels']} plex levels, sizes {diag['smallest_stratum']}-"
        f"{diag['largest_stratum']}, {diag['n_singleton_strata']} singleton strata; "
        f"cases held fixed {diag['n_cases_fixed']}/{n} = {diag['frac_cases_fixed']:.1%}")
    log(f"  P2 design: {len(dummy_preview)} dummies (levels >= {C1L.MIN_GROUP} cases, "
        f"one dropped for full rank: {dropped_preview}), "
        f"{len(left_preview)} levels left to the intercept")
    log(f"  P3: {len(p3_present)}/{len(P3_DROP_CASES)} listed cases analysed {p3_present}; "
        + ", ".join(f"{c} -> {plex_map.get(c, '(no row)')}" for c in sorted(P3_DROP_CASES)))
    if diag["n_unlabelled"]:
        log(f"  WARNING: {diag['n_unlabelled']} analysed cases have no plex label; they are "
            f"held fixed and get no dummy: "
            f"{[c for c, v in zip(common, raw) if not v]}")

    if args.blind:
        log(f"\nBLIND: n={n} K={K} n_tested={len(tested)} n_h1_selected(usable)="
            f"{len(sel_h1)} n_h2_selected(usable)={len(sel_h2)} plex_levels="
            f"{diag['n_levels']} dummies={len(dummy_preview)} p3_cases_present="
            f"{len(p3_present)}")
        log("Nothing further is computed, printed or written in --blind mode.")
        return 0

    # --- step 1 reproduction -------------------------------------------------------
    log("\n=== STEP 1 REPRODUCTION (2,310 set, unrestricted null, primary seed) ===")
    pcs = PCA(n_components=K, svd_solver="full").fit_transform(W)
    rng = np.random.RandomState(zlib.crc32(SEED_UNRESTRICTED))
    perms_u = [rng.permutation(n) for _ in range(args.perms)]
    step1 = C1L.run_family("step1", tested, sel_h1, R, P, pcs, perms_u, disc_incr,
                           workers, log)
    repro = dict(n=n, K=K, n_usable=step1["n_usable"], n_selected=step1["n_selected"],
                 auc_obs=step1["auc_obs"], p1=step1["p_auc"])
    gate = check_gate(repro, args.perms, log)

    if not real:
        log(f"\nSMOKE: B={args.perms}. The gate fields above were checked on this host; "
            f"no within-plex null is computed and nothing is written.")
        log(f"elapsed {round((time.time() - t0) / 60.0, 1)} min")
        return 0

    # --- P1 ------------------------------------------------------------------------
    log("\n=== P1 (prespec L1): 2,310 set, same PCs, within-plex null (post hoc) ===")
    perms_plex = build_within_strata_perms(raw, args.perms, SEED_PLEX)
    p1 = C1L.run_family("P1", tested, sel_h1, R, P, pcs, perms_plex, disc_incr,
                        workers, log)
    for key in ("auc_obs", "rho_obs"):
        if abs(p1[key] - step1[key]) > OBS_IDENTITY_TOL:
            fatal(f"P1 observed {key} {p1[key]!r} != step 1 {step1[key]!r}: these must be "
                  f"identical (same PCs, same genes; only the null differs). This is a bug.")
    if p1["genes"] != step1["genes"]:
        fatal("P1 usable gene list differs from step 1. This is a bug.")

    # --- P2 ------------------------------------------------------------------------
    log("\n=== P2 (prespec L2): 726 set, PCs residualised on plex dummies, same "
        "within-plex draws ===")
    pcs_resid, design = build_plex_dummy_design(raw, pcs, C1L.MIN_GROUP)
    log(f"  plex residualisation: {design['n_dummies']} dummies (dropped for full rank: "
        f"{design['dropped_for_full_rank']}; left to intercept: "
        f"{design['levels_left_to_intercept']}), PC variance retained "
        f"{design['pc_variance_retained']:.1%} (n={n})")
    p2 = C1L.run_family("P2", tested, sel_h2, R, P, pcs_resid, perms_plex, disc_incr_bres,
                        workers, log)

    # --- P3 ------------------------------------------------------------------------
    log("\n=== P3 (prespec L3): P1 without C3N-02230 and C3N-02240 (aliquot in two "
        "plexes) ===")
    if p3_present:
        (common3, tested3, sel_h1_3, _sel_h2_3, disc_incr3, _f3, _b3,
         R3, P3, W3, K3) = C1L.prepare(args, log, drop_cases=P3_DROP_CASES)
        if tested3 != tested or sel_h1_3 != sel_h1:
            fatal("P3 gene universe differs from step 1. This is a bug.")
        raw3 = labels_for(common3, plex_map)
        diag3 = strata_diagnostics(raw3)
        pcs3 = PCA(n_components=K3, svd_solver="full").fit_transform(W3)
        perms3 = build_within_strata_perms(raw3, args.perms, SEED_P3)
        p3 = C1L.run_family("P3", tested3, sel_h1_3, R3, P3, pcs3, perms3, disc_incr3,
                            workers, log)
        p3_meta = dict(n_cases=len(common3), K=K3, dropped_cases=p3_present,
                       strata=diag3)
    else:
        log("  neither case is in the analysed set; P3 would equal P1 and is not run.")
        p3, p3_meta = None, dict(n_cases=n, K=K, dropped_cases=[], strata=None)

    # --- summaries -----------------------------------------------------------------
    def summ(fam):
        if fam is None:
            return None
        d = C1L.strip_arrays(fam)
        d["z_auc"] = _z(d["auc_obs"], d["auc_null_mean"], d["auc_null_sd"])
        d["z_rho"] = _z(d["rho_obs"], d["rho_null_mean"], d["rho_null_sd"])
        return d

    s1, sp1, sp2, sp3 = summ(step1), summ(p1), summ(p2), summ(p3)
    sp1["frac_cases_fixed"] = diag["frac_cases_fixed"]
    sp1["null_sd_ratio_vs_unrestricted"] = (
        sp1["auc_null_sd"] / s1["auc_null_sd"] if s1["auc_null_sd"] > 0 else float("nan"))
    sp2.update(design)
    if sp3 is not None:
        sp3["frac_cases_fixed"] = p3_meta["strata"]["frac_cases_fixed"]

    log("\n" + "=" * 78)
    log("POST HOC SUMMARY (does not change the registered C1-LUAD reading)")
    log("=" * 78)
    log(f"step1 (reproduced): AUC={s1['auc_obs']:.6f} null {s1['auc_null_mean']:.4f}"
        f"+/-{s1['auc_null_sd']:.4f} p1={s1['p_auc']:.6f}")
    for name, d in (("P1", sp1), ("P2", sp2), ("P3", sp3)):
        if d is None:
            log(f"{name} (prespec {PRESPEC_LABELS[name]}): not run")
            continue
        log(f"{name} (prespec {PRESPEC_LABELS[name]}): AUC={d['auc_obs']:.6f} "
            f"null {d['auc_null_mean']:.4f}+/-"
            f"{d['auc_null_sd']:.4f} z={d['z_auc']:.2f} p={d['p_auc']:.6f} | "
            f"rho={d['rho_obs']:.4f} null {d['rho_null_mean']:.4f}+/-{d['rho_null_sd']:.4f} "
            f"z={d['z_rho']:.2f} p_rho={d['p_rho']:.6f} | n_sel={d['n_selected']} "
            f"n_bg={d['n_background']}")
    log(f"P1 cases held fixed: {sp1['frac_cases_fixed']:.1%}; P2 dummies "
        f"{sp2['n_dummies']}, PC variance retained {sp2['pc_variance_retained']:.1%}")

    minutes = round((time.time() - t0) / 60.0, 1)
    with open(out["npz"], "xb") as fh:
        np.savez_compressed(
            fh, common=np.array(common), plex_raw=raw,
            genes_p1=np.array(p1["genes"]), sel_mask_p1=p1["sel_mask"], M_p1=p1["M"],
            genes_p2=np.array(p2["genes"]), sel_mask_p2=p2["sel_mask"], M_p2=p2["M"],
            genes_p3=np.array(p3["genes"] if p3 else []),
            sel_mask_p3=(p3["sel_mask"] if p3 else np.array([], dtype=bool)),
            M_p3=(p3["M"] if p3 else np.empty((0, 0))))
    result = dict(
        post_hoc=True,
        governed_by=POST_HOC_GOVERNED_BY,
        registered_reading_unchanged=True,
        note=("Post hoc plex-stratified re-read of the registered C1-LUAD test. The "
              "signed rule has no plex stratum (C1_LUAD_ADDENDA.md Entry 7); nothing "
              "here changes p1 or clauses [B], [C], [D]."),
        prespec_section="L", prespec_labels=PRESPEC_LABELS,
        cohort="luad_c1", host=socket.gethostname(),
        lsf_job=os.environ.get("LSB_JOBID", ""),
        finished_utc=_utc(), minutes=minutes,
        perms=args.perms, workers=workers,
        seeds=dict(step1=SEED_UNRESTRICTED.decode(), p1_p2=SEED_PLEX.decode(),
                   p3=SEED_P3.decode()),
        n_cases=n, K=K, n_genes_tested=len(tested),
        reproduction_gate=gate,
        step1=s1,
        plex_strata=diag,
        P1=sp1, P2=sp2, P3=sp3, P3_meta=p3_meta,
        paths=dict(rule=args.rule, prespec=args.prespec, plex_map=args.plex_map,
                   protein=args.protein, case_ids=args.case_ids,
                   discovery_results=args.discovery_results,
                   discovery_sitepack=args.discovery_sitepack,
                   discovery_slide_map=args.discovery_slide_map,
                   discovery_xwalk=args.discovery_xwalk,
                   discovery_rna_manifest=args.discovery_rna_manifest,
                   batch_labels=args.batch_labels,
                   confirmatory_rna_manifest=env["MORPHO_RNA_MANIFEST"],
                   confirmatory_slide_map=env["MORPHO_SLIDE_MAP"],
                   wsi_emb_dir=env["MORPHO_WSI_EMB_DIR"], out=env["MORPHO_OUT"],
                   log=out["log"], npz=out["npz"]),
        sha256_inputs=hashes,
        sha256_code=dict(this_script=hashes["this_script"], **code_hashes),
        sha256_npz=sha256(out["npz"]),
        versions=dict(python=platform.python_version(), numpy=np.__version__,
                      pandas=pd.__version__,
                      sklearn=__import__("sklearn").__version__),
    )
    with open(out["result"], "x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=2,
                  default=lambda o: (o.tolist() if hasattr(o, "tolist") else str(o)))
    log(f"\n  -> {out['result']}\n  -> {out['npz']}\n  -> {out['log']}\n"
        f"  elapsed {minutes} min")
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
