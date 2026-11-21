#!/usr/bin/env python
"""Post hoc re-reads of C1-UCEC and C1-LUAD with transcriptome-derived covariates: the local stages.

Binding text: review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md, VERSION 3.1 (sections cited as
"PS n").  This script implements PS 3 (discovery analyses and sets), PS 5 (cells, the change D,
the capacity randomisation, the diagnostic), PS 6 (reading rules) and the local parts of PS 9
(self-tests, V1, V1b, V2, V3, the sidecar gate, Stage L1, Stage L2).

Modes (exactly one)
  --selftest                synthetic data only; its first line is "sha256 <hex> c1_covariates_local.py" (the
                            hash of the running script file).  Covers (kernel, BH, AUC, set ranking and set summaries,
                            partitions, canonical correlations, file round trips, parsers,
                            CLOSED_PARTS.txt, manifest, PS 6.1 types with closures, D and p_abs,
                            unlimited reruns, the sidecar check, and two full end-to-end runs of
                            Stage L1 and Stage L2 on tiny fake cohorts in a temporary folder that
                            never touches the real project files, including pre-flight refusals
                            that must create nothing).
  --validate-only           published quantities only (PS 9, before freezing; refused once the sidecar
                            exists, so that frozen json files are never rewritten):
      OP   operator_dummies on the real labels with the binding counts of PS 2 (12 indicators
           covering 73 of 100 / 75 of 105 cases; label counts only);
      S_BC the batch-corrected sets of the operator-axis tables (PS 3: 201 / 726 genes; the LUAD
           file must have the sha256 recorded in the registered result json);
      V1   the discovery code path, fed the published block, reproduces every fit gene's
           published increment (1e-9), own-mRNA R2 (1e-10), p and FDR (1e-12), n, and the
           published selection sets exactly;
      V1b  with the 20 expression components in the baseline the same code reproduces, for the
           published-set genes, r2_rna_pcs20 (1e-9 in both organs) and incr_over_pcs20 (1e-9 in
           UCEC, 1e-5 in LUAD) of scripts__transcriptome_baseline_<cohort>_perm.csv; the largest
           differences are recorded;
      V2   `assemble` on the saved confirmatory matrices reproduces the published readings
           (AUC, null mean, null SD, p, rho, p_rho to 1e-12; n selected / background).
  --gate-v3                 PS 9 "H0 and V3": the transferred discovery score tables reproduce the
                            four gain columns of results__composition_controls_<cohort>.csv (1e-9).
                            No morphology is involved.  Refused once the sidecar exists.
  --check-stage-h           format check only of whatever Stage H files are present (keys, dtypes,
                            shapes, gene order, finiteness, cases, score-table header and patients,
                            manifest syntax and hashes).  It computes no statistic.
  --check-sidecar           verifies the sidecar against the complete list of files the stages
                            require (listed + matching hash + pinned hashes) and prints the result.
  --stage L1 --manifest-sha <hex>
                            Stage L1 (PS 9; after the Stage H record).  Needs the frozen sidecar
                            and the Stage H manifest whose sha256 is given on the command line.
  --stage L2 --manifest-sha <hex>
                            Stage L2 (PS 9).  Needs the frozen sidecar, the same manifest and a
                            completed Stage L1 whose outputs have the hashes in result_L1.json.

Pre-flight (PS 9): every check that can fail on an interface or transfer slip runs before
<stage>.partial is created -- the sidecar, CLOSED_PARTS.txt, the manifest's sha256; for L1 the
score tables, V3 and V2; for L2 the Stage L1 hashes and the existence, manifest hash and format
of every Stage H file that is not closed.  A refusal creates nothing.

Run safety (PS 9 "Interrupted runs"): a stage writes into <stage>.partial and renames it to
<stage> only after result_<stage>.json is written; "completed" means that file exists; there is
no limit on attempts: a leftover <stage>.partial is moved to <stage>_failed<k> (k = 1, 2, ...)
and the stage runs again; a completed stage is never rerun.  Everything printed during a stage
is teed to <stage>.log.

The estimator is the published one (residual_analysis.py): per gene, on the patients with a
protein value, 5-fold KFold(shuffle, random_state 0), ridge alpha 1 on columns standardised on
the training fold, increment = R2(baseline + block) - R2(baseline), null = the block's rows
permuted over all patients with the published permutations (RandomState(0), 1,000 sequential
draws) and then restricted to the fit patients, one-sided p = (1 + #{null >= observed}) / 1001,
Benjamini-Hochberg over all fit genes.  `bh_fdr` is copied unchanged from residual_analysis.py,
`assemble` from c1_run_test.py and `cancorr` from transcriptome_baseline.py (_cancorr).  The
numerical kernel (ref_increment, batched_increments, _job, run_block) is that of the version
validated on 2026-11-21 (VERSION 2 of the text), unchanged.
"""
import os

# One BLAS thread per worker: the per-gene loop is parallelised over processes, and an
# oversubscribed BLAS would change nothing numerically but would slow every worker down.
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
for _v in THREAD_VARS:
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import collections  # noqa: E402
import csv  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import platform  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))                # MorphoResidual_paper/
PRESPEC_REL = "review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md"
OUT_REL = "server_export/results/c1_covariates"
STAGE_H_REL = OUT_REL + "/stage_h"
ANALYSIS_DIR_REL = "review/C1_COVARIATES_2026-11-21"
LOCAL_SCRIPT_REL = ANALYSIS_DIR_REL + "/c1_covariates_local.py"
CLUSTER_SCRIPT_REL = ANALYSIS_DIR_REL + "/c1_covariates_cluster.py"
WRAPPER_RELS = tuple(ANALYSIS_DIR_REL + "/lsf_c1_cov_%s.sh" % k for k in ("ucec", "luad", "discovery"))
CLOSED_REL = ANALYSIS_DIR_REL + "/CLOSED_PARTS.txt"
PINNED_REL = ANALYSIS_DIR_REL + "/PINNED_HOSTS.txt"
MAKE_ABSTRACT_REL = ANALYSIS_DIR_REL + "/make_abstract_variants.py"
ABSTRACT_VARIANTS_REL = ANALYSIS_DIR_REL + "/ABSTRACT_VARIANTS.md"
ABSTRACT_FINAL_REL = "review/QUALITY_REVIEW_2026-11-20/abstract_v12/ABSTRACT_FINAL.tex"
SIDECAR_LIST_REL = ANALYSIS_DIR_REL + "/SIDECAR_FILES.txt"
KERNEL_REL = "review/MEAS_ERROR_2026-10-08/theta_kernel.py"
KERNEL_DIR = os.path.join(PAPER, "review", "MEAS_ERROR_2026-10-08")

# PS 1 "Code the definitions come from" (PS 1 prints 8-hex prefixes; these are the full values of
# the repository copies, and the sidecar must give exactly these).
DEFINITION_MODULES = {
    "server_export/scripts/residual_analysis.py":
        "a54a0a56494d315e15ef075425fe2133ecd9f05193c22656c1145bffdaac4a0a",
    "split_replication_power.py":
        "79530c075c0f7e75bf06d3c2f1edfe8792cbf4d6c71ae1d32e60bf44f31a0397",
    "server_export/scripts/c1_run_test.py":
        "fc578c705e803ddb118ff1f65a61079723167985a66adf7285834bf4709391b5",
    "server_export/scripts/c1_run_test_luad.py":
        "1d5ee05d57f80af4d90870aecb10e574e4c910a8ce1eab743a067bd9a9278742",
    "server_export/scripts/composition_controls.py":
        "46b4a1d8a6973bb71b162697a04754e0e1c555eeaf3ff83dfc5c143087abdfa6",
    "server_export/scripts/transcriptome_baseline.py":
        "761df7d714285aa4941adcd358370c0836739578fbfe5874efd2aaee376359c0",
    "server_export/scripts/residual_analysis_sitepack.py":
        "d0cfeb01ef3d249455d4d26749bbec01b2c6bb3ada6509790f349073b0edbf52",
    KERNEL_REL:
        "6c95cd119d28e811e09967946e9f191c69f07f54f653991917e817ef99ec2734",
}
KERNEL_SHA = DEFINITION_MODULES[KERNEL_REL]
ABSTRACT_FINAL_SHA = "6a1bfede2f9a56dd811cf46b345ecd2989f7ffb5697b1578e574fa8fec298ddc"   # PS 1 (6a1bfede...)

sys.path.insert(0, KERNEL_DIR)
import theta_kernel as tk  # noqa: E402  (reference ridge = residual_analysis.cv_r2, folds, permutations)

B_PERM = 1000           # PS 3, PS 4: published number of permutations
FDR_LEVEL = 0.05        # PS 3: selection threshold
ALPHA = 0.05            # PS 6: "outranked" means p < 0.05 (p unrounded)
MIN_SET = 50            # PS 5: a set or background with fewer genes in the universe is untestable
MIN_FIT = 30            # published MIN_N: a gene needs 30 protein values
MAX_OP_LEVELS = 12      # PS 2: operator indicators
MIN_OP_GROUP = 3
N_CAP = 119             # PS 2: capacity placebos j = 1..119
CAP_COLS = 25           # PS 4: placebo files hold columns 0 and 1-25
TOL_REF = 1e-9          # batched kernel against the reference kernel
TOL_PUB = 1e-10         # own-mRNA R2 against the published value (PS 9 V1)
TOL_INC = 1e-9          # increments against published values (PS 9 V1, V1b R2 and UCEC increment, V3)
TOL_PFDR = 1e-12        # p and FDR against published values (PS 9 V1)
TOL_AUC = 1e-12         # `assemble` against the published confirmatory readings (PS 9 V2)
TOL_CANCORR = 1e-9      # self-test only: cancorr against a direct computation

SCORE_COLS = ("stroma", "immune", "proliferation")                      # PS 2, in this order
ANALYSES = ("C3", "R20", "OP", "M|C3", "M|R20")                         # PS 3 (M is published)
FORWARD_X = ("C3", "R20", "OP")                                         # PS 5 forward family
MAIN_ARMS = ("A_C3", "A_R20", "V_C3", "V_R20")                          # PS 4, 1,001 columns
ARM_TOKENS = ("A_C3", "A_R20", "A_R20cap", "V_C3", "V_R20")             # CLOSED_PARTS.txt arms
ANALYSIS_TOKENS = ("C3", "M|C3")                                        # CLOSED_PARTS.txt analyses: only these can be closed (PS 9)
MANIFEST_NAME = "c1_cov_stage_h_manifest.sha256"
GAUSS_SEED = 42         # PS 5 diagnostic: reference block RandomState(42).randn(n, K)
HEX64 = re.compile(r"[0-9a-fA-F]{64}")

REAL_COHORTS = {
    "ucec": dict(
        label="uterus",
        inputs="review/MEAS_ERROR_2026-10-08/export/ucec/inputs_ucec.npz",
        inputs_sha="e28b86fadd41ec86363577fffad1db1f2f8a5162e732c7acb0e4133c0982044e",
        conf="server_export/scripts/pinned/c1_ucec_d4_stratified.npz",
        conf_sha="ef75f1c06d4fc33843a79e9a8d1bb4f4224acd201d30784badb97cf3be13384a",
        arms=[("unrestricted", "M_unrestricted"), ("within-operator", "M_within_operator"),
              ("residualised", "M_within_operator_residualized")],
        published="server_export/scripts/pinned/c1_ucec_d4_followup.json", published_kind="ucec",
        compctl="server_export/results/results__composition_controls_ucec.csv",
        txperm="server_export/results/scripts__transcriptome_baseline_ucec_perm.csv",
        n_disc=100, n_fit=10512, n_sel=2566, n_univ=10146, n1=2536, n0=7610,
        n_op_dummies=12, n_op_covered=73, n_conf_cases=138, K=28,
        tol_v1b_inc=1e-9,
        sitepack="server_export/results/ucec__sitepack_operator.csv", n_bc=201, n_bc_universe=198,
        sitepack_prefix=None, sitepack_pub_key=None, sitepack_sha=None),
    "luad": dict(
        label="lung",
        inputs="review/MEAS_ERROR_2026-10-08/export/luad/inputs_luad.npz",
        inputs_sha="b64cdf3b31717a32e0d42c42f3f1dde74786d5598008f0d214297a71adaba095",
        conf="server_export/scripts/pinned/luad_c1/c1_luad_primary_primary.npz",
        conf_sha="2034e71143c2e065675d3159832b1a7be17aa38c7feeaacfc82862c411ea2aab",
        arms=[("unrestricted", "M")],
        published="server_export/scripts/pinned/luad_c1/c1_luad_primary_result_primary.json",
        published_kind="luad",
        compctl="server_export/results/results__composition_controls_luad.csv",
        txperm="server_export/results/scripts__transcriptome_baseline_luad_perm.csv",
        n_disc=105, n_fit=10724, n_sel=2310, n_univ=10406, n1=2287, n0=8119,
        n_op_dummies=12, n_op_covered=75, n_conf_cases=112, K=22,
        tol_v1b_inc=1e-5,           # PS 9 V1b: wider in LUAD (reference run's own float noise)
        # PS 3 S_BC: 726 genes; the file's sha256 is recorded in the registered result json (read at run time)
        sitepack="server_export/results/luad__sitepack_operator.csv", n_bc=726, n_bc_universe=721,
        sitepack_prefix="aeef7232", sitepack_pub_key="sha256_discovery_sitepack",
        sitepack_sha="aeef7232682c24d209024ba73efea8033cb7486779831bda1a1ba3c8bda45d9f"),
}


class Fatal(Exception):
    """A refusal or a fatal error.  Raised (not sys.exit) so the self-test can check that bad
    inputs are refused, and so a stage records the failure in its log before it stops."""


def fatal(msg):
    raise Fatal(msg)


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stage_h_name(key, arm):
    return "c1_cov_%s_%s.npz" % (key, arm)


def cap_name(key, j):
    """PS 4 / shared interface: c1_cov_<cohort>_A_R20cap_<jjj>.npz, jjj zero-padded to 3 digits."""
    return "c1_cov_%s_A_R20cap_%03d.npz" % (key, j)


def scores_name(key):
    return "c1_cov_discovery_scores_%s.tsv" % key


def covariates_record_name(key):
    return "c1_cov_confirmatory_covariates_%s.npz" % key


class Ctx:
    """Every path and constant a stage uses.  The real run uses the project root; the self-test
    builds a Ctx on a temporary root with the same relative layout, so the very same stage code
    runs end to end on synthetic files and can never write into the real result folders."""

    def __init__(self, root, cohorts, B, synthetic, pins=None):
        self.root = os.path.abspath(root)
        self.cohorts = cohorts
        self.B = B
        self.synthetic = synthetic
        self.pins = dict(pins or {})                 # rel path -> full sha256 the sidecar must give
        self.script_path = self.P(LOCAL_SCRIPT_REL)
        self.kernel_path = self.P(KERNEL_REL)
        self.closed_path = self.P(CLOSED_REL)
        self.pinned_path = self.P(PINNED_REL)
        self.prespec = self.P(PRESPEC_REL)
        self.sidecar = self.prespec + ".sha256"
        self.out_dir = self.P(OUT_REL)
        self.stage_h_dir = self.P(STAGE_H_REL)

    def P(self, rel):
        return os.path.join(self.root, *rel.split("/"))

    def H(self, name):
        return os.path.join(self.stage_h_dir, name)

    def required_rels(self, closed):
        """Every file the stages need frozen in the sidecar (task item 10; PS 1, PS 9): the
        pre-specification, both scripts, the three wrappers, CLOSED_PARTS.txt, the eight
        definition modules, the Abstract input, generator and output, and per cohort the inputs,
        the confirmatory matrix, the published results, the operator-axis (sitepack) table and, unless
        both discovery analyses that use them are closed, the discovery score table; PINNED_HOSTS.txt."""
        rels = [PRESPEC_REL, LOCAL_SCRIPT_REL, CLUSTER_SCRIPT_REL] + list(WRAPPER_RELS) + \
            [CLOSED_REL, PINNED_REL, MAKE_ABSTRACT_REL, ABSTRACT_FINAL_REL, ABSTRACT_VARIANTS_REL] + \
            list(DEFINITION_MODULES)
        for key, c in self.cohorts.items():
            rels += [c["inputs"], c["conf"], c["published"], c["compctl"], c["txperm"], c["sitepack"]]
            if closed is None or needs_scores(closed, key):
                rels.append(STAGE_H_REL + "/" + scores_name(key))
        out = []
        for r in rels:
            if r not in out:
                out.append(r)
        return out


def real_ctx():
    pins = dict(DEFINITION_MODULES)
    pins[ABSTRACT_FINAL_REL] = ABSTRACT_FINAL_SHA
    for c in REAL_COHORTS.values():
        pins[c["inputs"]] = c["inputs_sha"]
        pins[c["conf"]] = c["conf_sha"]
        if c.get("sitepack_sha"):
            pins[c["sitepack"]] = c["sitepack_sha"]
    return Ctx(PAPER, REAL_COHORTS, B_PERM, False, pins)


def env_record():
    import scipy
    import sklearn
    return dict(python=sys.version.split()[0], platform=platform.platform(), numpy=np.__version__,
                scipy=scipy.__version__, sklearn=sklearn.__version__, cpu_count=os.cpu_count(),
                threads={v: os.environ.get(v) for v in THREAD_VARS})


# ------------------------------------------------------------------------------------------
# CLOSED_PARTS.txt (PS 6.1 "Closed readings", PS 9 "Closed parts")
# ------------------------------------------------------------------------------------------
def _control_lines(path):
    """The lines of a frozen control file (CLOSED_PARTS.txt, PINNED_HOSTS.txt).  Both local
    parsers and the cluster's parse these files identically (shared rules): UTF-8, no byte-order
    mark, LF line ends only (a CR anywhere is fatal)."""
    if not os.path.isfile(path):
        fatal("REFUSED: %s not found" % path)
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw.startswith(b"\xef\xbb\xbf"):
        fatal("REFUSED: %s starts with a byte-order mark" % path)
    if b"\r" in raw:
        fatal("REFUSED: %s contains a carriage return (LF line endings only)" % path)
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        fatal("REFUSED: %s is not UTF-8" % path)
    return text.split("\n")


def _control_tokens(line, path):
    """None for a comment ('#' in column 0) or an empty line; else exactly three tokens separated by
    single spaces, with no leading or trailing whitespace and no whitespace inside a token."""
    if line == "" or line.startswith("#"):
        return None
    parts = line.split(" ")
    if len(parts) != 3 or any(p == "" or re.search(r"\s", p) for p in parts):
        fatal("REFUSED: unknown token in %s: %r (exactly three tokens separated by single spaces)" % (path, line))
    return parts


def parse_closed(path, cohorts):
    """Shared strict rules (PS 9 "Closed parts"): '#' comment lines (the '#' in column 0) and empty
    lines are skipped; every other line is exactly three tokens separated by single spaces, with no
    leading or trailing whitespace: 'arm <cohort> <arm>' (arm in A_C3, A_R20, A_R20cap, V_C3, V_R20)
    or 'analysis <cohort> <X>' (X is C3 or the literal token M|C3; no other discovery analysis can
    be closed).  Anything else is fatal.  Consistency, fatal if broken: for a cohort A_R20cap is
    listed if and only if A_R20 is; 'analysis <cohort> C3' if and only if 'analysis <cohort> M|C3'.
    A file with comment lines only means that nothing is closed."""
    arms, analyses, tokens = set(), set(), []
    for ln in _control_lines(path):
        f = _control_tokens(ln, path)
        if f is None:
            continue
        if f[0] not in ("arm", "analysis") or f[1] not in cohorts:
            fatal("REFUSED: unknown token in %s: %r" % (path, ln))
        if f[0] == "arm" and f[2] not in ARM_TOKENS:
            fatal("REFUSED: unknown arm in %s: %r" % (path, ln))
        if f[0] == "analysis" and f[2] not in ANALYSIS_TOKENS:
            fatal("REFUSED: unknown analysis in %s: %r (only C3 and M|C3 can be closed)" % (path, ln))
        target = arms if f[0] == "arm" else analyses
        if (f[1], f[2]) in target:
            fatal("REFUSED: %r is listed twice in %s" % (ln, path))
        target.add((f[1], f[2]))
        tokens.append(ln)
    for key in cohorts:
        if ((key, "A_R20cap") in arms) != ((key, "A_R20") in arms):
            fatal("REFUSED: %s: for %s, A_R20cap is listed if and only if A_R20 is" % (path, key))
        if ((key, "C3") in analyses) != ((key, "M|C3") in analyses):
            fatal("REFUSED: %s: for %s, 'analysis C3' is listed if and only if 'analysis M|C3' is" % (path, key))
    return dict(arms=arms, analyses=analyses, tokens=tokens)


def parse_pinned(path, cohorts):
    """PINNED_HOSTS.txt (PS 9 H-V1): '#' comment lines (column 0), empty lines, and lines of exactly
    three single-space-separated tokens 'pinned <cohort> <host>' (cohort in `cohorts`); at most one
    line per cohort.  Returns {cohort: host}."""
    out = {}
    for ln in _control_lines(path):
        f = _control_tokens(ln, path)
        if f is None:
            continue
        if f[0] != "pinned" or f[1] not in cohorts:
            fatal("REFUSED: unknown token in %s: %r" % (path, ln))
        if f[1] in out:
            fatal("REFUSED: %s has two pinned lines for %s" % (path, f[1]))
        out[f[1]] = f[2]
    return out


def fully_closed(closed, key):
    """Every Stage H arm of the cohort is closed (then it needs no pinned host)."""
    return all((key, a) in closed["arms"] for a in ARM_TOKENS)


def check_pinned(ctx, closed):
    """Stage time (PS 9): PINNED_HOSTS.txt parses and every cohort that is not fully closed has a
    pinned line."""
    pinned = parse_pinned(ctx.pinned_path, ctx.cohorts)
    for key in ctx.cohorts:
        if not fully_closed(closed, key) and key not in pinned:
            fatal("REFUSED: %s has no pinned line for %s, which is not fully closed" % (PINNED_REL, key))
    return pinned


def refuse_if_frozen(ctx, mode):
    """--validate-only and --gate-v3 write json files that the sidecar freezes: after freezing they
    must not rewrite them."""
    if os.path.exists(ctx.sidecar):
        fatal("REFUSED: %s: the sidecar %s exists; after freezing this mode would rewrite frozen files"
              % (mode, ctx.sidecar))


def no_closed():
    return dict(arms=set(), analyses=set(), tokens=[])


def arm_closed(closed, key, arm):
    return (key, arm) in closed["arms"]


def an_closed(closed, key, X):
    return (key, X) in closed["analyses"]


def needs_scores(closed, key):
    """The discovery score table is read unless both analyses that use it are closed (PS 9 H0/V3)."""
    return not (an_closed(closed, key, "C3") and an_closed(closed, key, "M|C3"))


# ------------------------------------------------------------------------------------------
# fail-closed sidecar gate (PS 9: "every stage refuses to start unless every sidecar hash matches")
# ------------------------------------------------------------------------------------------
def sidecar_problems(ctx):
    """Returns (problems, listed {rel: sha256}, closed or None).  Every line '<sha256> *<path
    relative to MorphoResidual_paper/>' is verified; '#' lines are comments; every required file
    must be listed; a pinned file must be listed with its pinned hash; CLOSED_PARTS.txt must be
    listed and parse."""
    probs, listed = [], {}
    if not os.path.exists(ctx.sidecar):
        return ["sidecar %s not found -- the pre-specification is not frozen" % ctx.sidecar], listed, None
    with open(ctx.sidecar, "rb") as fh:
        text = fh.read().decode("utf-8", errors="replace")
    for ln in text.split("\n"):
        line = ln.rstrip("\r").strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2 or not HEX64.fullmatch(parts[0]):
            probs.append("malformed sidecar line %r" % line)
            continue
        h = parts[0].lower()
        rel = parts[1].strip()
        rel = rel[1:] if rel.startswith("*") else rel
        if rel.startswith("/") or "\\" in rel or ".." in rel.split("/") or not rel:
            probs.append("sidecar path %r is not relative to MorphoResidual_paper/ with forward slashes" % rel)
            continue
        if rel in listed:
            probs.append("%s is listed twice in the sidecar" % rel)
        path = ctx.P(rel)
        if not os.path.isfile(path):
            probs.append("%s is listed in the sidecar but missing" % rel)
        elif sha256(path) != h:
            probs.append("%s changed since freezing (sidecar hash does not match)" % rel)
        if rel in ctx.pins and h != ctx.pins[rel]:
            probs.append("%s is listed with %s..., not the pinned %s..." % (rel, h[:8], ctx.pins[rel][:8]))
        listed[rel] = h
    closed = None
    if CLOSED_REL not in listed:
        probs.append("%s is not listed in the sidecar (not frozen)" % CLOSED_REL)
    else:
        try:
            closed = parse_closed(ctx.closed_path, ctx.cohorts)
        except Fatal as e:
            probs.append(str(e))
    if PINNED_REL in listed:
        try:
            parse_pinned(ctx.pinned_path, ctx.cohorts)
        except Fatal as e:
            probs.append(str(e))
    for rel in ctx.required_rels(closed):
        if rel not in listed:
            probs.append("%s is not listed in the sidecar (not frozen)" % rel)
    return probs, listed, closed


def check_sidecar(ctx):
    """Stage gate: refuses on any problem.  Returns (closed, number of listed files)."""
    probs, listed, closed = sidecar_problems(ctx)
    if probs:
        fatal("REFUSED: sidecar check failed (%d problems): %s" % (len(probs), " | ".join(probs[:8])))
    log("sidecar verified: %d files; %d required files listed; CLOSED_PARTS.txt: %s"
        % (len(listed), len(ctx.required_rels(closed)), closed["tokens"] or "nothing closed"))
    return closed, len(listed)


def check_sidecar_mode(ctx):
    """--check-sidecar: prints every problem and the status of every required file."""
    probs, listed, closed = sidecar_problems(ctx)
    req = ctx.required_rels(closed)
    for rel in req:
        st = "listed" if rel in listed else "NOT LISTED"
        if rel in listed and os.path.isfile(ctx.P(rel)) and sha256(ctx.P(rel)) != listed[rel]:
            st = "HASH MISMATCH"
        print("  %-13s %s" % (st, rel))
    extra = sorted(set(listed) - set(req))
    print("  %d further files listed in the sidecar (not required by the stages)" % len(extra))
    if closed is not None:
        print("  CLOSED_PARTS.txt: %s" % (closed["tokens"] or "nothing closed"))
    for p in probs:
        print("  PROBLEM: %s" % p)
    print("check-sidecar: %s (%d required files, %d listed, %d problems)"
          % ("SIDECAR OK" if not probs else "SIDECAR FAILED", len(req), len(listed), len(probs)))
    return not probs


# ------------------------------------------------------------------------------------------
# copied unchanged
# ------------------------------------------------------------------------------------------
def bh_fdr(pvals):
    """Identical to residual_analysis.bh_fdr."""
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    fdr = np.empty(n)
    fdr[order] = np.clip(ranked, 0, 1)
    return fdr


def assemble(M, xs, sel_mask):
    """Identical to c1_run_test.assemble."""
    n1, n0 = int(sel_mask.sum()), int((~sel_mask).sum())
    ranks = rankdata(M, axis=0)
    auc = (ranks[sel_mask].sum(axis=0) - n1 * (n1 + 1) / 2.0) / (n1 * n0)
    rs = rankdata(xs)
    rs = (rs - rs.mean()) / rs.std()
    rc = (ranks - ranks.mean(axis=0)) / ranks.std(axis=0)
    rho = (rs[:, None] * rc).mean(axis=0)
    return auc, rho, n1, n0


def cancorr(A, B):
    """Copied from transcriptome_baseline.py (_cancorr, inside main(); sha256 761df7d7...),
    unchanged: canonical correlations of the column-standardised blocks via QR and SVD."""
    Az = (A - A.mean(0)) / np.where(A.std(0) == 0, 1.0, A.std(0))
    Bz = (B - B.mean(0)) / np.where(B.std(0) == 0, 1.0, B.std(0))
    Qa, _ = np.linalg.qr(Az); Qb, _ = np.linalg.qr(Bz)  # noqa: E702
    return np.clip(np.linalg.svd(Qa.T @ Qb, compute_uv=False), 0, 1)


def emp_p(obs, null):
    return float((1 + int(np.sum(null >= obs))) / (1 + len(null)))


# ------------------------------------------------------------------------------------------
# the estimator: reference (unbatched) and batched over permutations (validated kernel, unchanged)
# ------------------------------------------------------------------------------------------
def ref_increment(y, base, block_v, folds):
    """Reference: R2([base, block]) - R2(base) with theta_kernel.oof_ridge (= RA.cv_r2)."""
    r2_base = tk.r2_of(y, tk.oof_ridge(base, y, folds))
    r2_full = tk.r2_of(y, tk.oof_ridge(np.hstack([base, block_v]), y, folds))
    return r2_base, r2_full - r2_base


def batched_increments(y, base, block_full, v, idx, alpha=tk.RIDGE_ALPHA):
    """Increment of the block over the baseline for every row order in `idx`.

    y           (nv,)   protein of the gene's fit patients
    base        (nv,q)  baseline columns on the fit patients (never permuted)
    block_full  (n,k)   block rows of all cohort patients
    v           (n,)    boolean fit mask
    idx         (R,n)   row orders; row 0 must be the identity.  As in residual_analysis.py the
                        block is permuted over all n patients first and then restricted to the
                        fit patients: block_full[perm][v].
    Returns r2_base (reference kernel) and the R increments."""
    nv = int(v.sum())
    folds = tk.folds_for(nv)
    r2_base = tk.r2_of(y, tk.oof_ridge(base, y, folds))
    blk = block_full[idx][:, v, :]                               # (R, nv, k)
    R = blk.shape[0]
    q = base.shape[1]
    p = q + blk.shape[2]
    eye = alpha * np.eye(p)
    preds = np.empty((R, nv))
    for tr, te in folds:
        Xtr = np.empty((R, len(tr), p))
        Xte = np.empty((R, len(te), p))
        Xtr[:, :, :q] = base[tr]
        Xte[:, :, :q] = base[te]
        Xtr[:, :, q:] = blk[:, tr, :]
        Xte[:, :, q:] = blk[:, te, :]
        mu = Xtr.mean(axis=1, keepdims=True)
        sd = Xtr.std(axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        Xtr = (Xtr - mu) / sd
        Xte = (Xte - mu) / sd
        ym = y[tr].mean()
        Xt = Xtr.transpose(0, 2, 1)
        G = np.matmul(Xt, Xtr) + eye
        rhs = np.matmul(Xt, (y[tr] - ym)[None, :, None])
        w = np.linalg.solve(G, rhs)
        preds[:, te] = np.matmul(Xte, w)[:, :, 0] + ym
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot > 0:
        r2 = 1.0 - np.sum((y[None, :] - preds) ** 2, axis=1) / ss_tot
    else:
        r2 = np.zeros(R)
    return r2_base, r2 - r2_base


_S = {}


def _init(mrna, prot, base_extra, block, idx, n_refcheck):
    _S.update(mrna=mrna, prot=prot, base_extra=base_extra, block=block, idx=idx, n_refcheck=n_refcheck)


def _job(g):
    y_all = _S["prot"][:, g]
    v = ~np.isnan(y_all)
    nv = int(v.sum())
    y = y_all[v]
    base = _S["mrna"][v, g].reshape(-1, 1)
    if _S["base_extra"] is not None:
        base = np.hstack([base, _S["base_extra"][v]])
    r2_base, inc = batched_increments(y, base, _S["block"], v, _S["idx"])
    folds = tk.folds_for(nv)
    # the batched kernel is checked against the reference kernel on every gene: the observed
    # increment, and the first n_refcheck permutations
    worst = 0.0
    for r in range(min(1 + _S["n_refcheck"], len(inc))):
        ref = ref_increment(y, base, _S["block"][_S["idx"][r]][v], folds)[1]
        worst = max(worst, abs(ref - inc[r]))
    null = inc[1:]
    exceed = int(np.sum(null >= inc[0]))
    return g, nv, r2_base, float(inc[0]), exceed, float(null.mean()) if len(null) else np.nan, worst


def run_block(d, block, base_extra, perms, workers, tag, n_refcheck=2):
    """Per-gene increment of `block` over [own mRNA, base_extra], its permutation p and BH FDR
    over all fit genes of the cohort.  workers <= 1 runs in-process (same _job, same numbers):
    the self-test uses it to stay fast, and it avoids a process pool on tiny inputs."""
    n, G = d["protein"].shape
    if block.shape[0] != n or perms.shape[1] != n:
        fatal("FATAL: %s: block or permutations do not have the cohort's %d rows" % (tag, n))
    idx = np.vstack([np.arange(n, dtype=np.int64)[None, :], perms])
    t0 = time.time()
    out = np.full((G, 6), np.nan)
    initargs = (d["mrna_log2tpm1"], d["protein"], base_extra, block, idx, n_refcheck)
    if workers <= 1:
        _init(*initargs)
        for k, r in enumerate(map(_job, range(G))):
            out[r[0]] = r[1:]
    else:
        with mp.Pool(workers, initializer=_init, initargs=initargs) as pool:
            for k, r in enumerate(pool.imap_unordered(_job, range(G), chunksize=16)):
                out[r[0]] = r[1:]
                if (k + 1) % 2000 == 0:
                    log("  %s: %d/%d genes, %.0f s" % (tag, k + 1, G, time.time() - t0))
    if np.isnan(out[:, :5]).any():
        fatal("FATAL: %s: a gene returned no result" % tag)
    worst = float(np.nanmax(out[:, 5]))
    if not worst <= TOL_REF:
        fatal("FATAL: %s: batched kernel differs from the reference kernel by %.3e" % (tag, worst))
    B = perms.shape[0]
    pval = (1 + out[:, 3]) / (1 + B)
    fdr = bh_fdr(pval)
    res = dict(n=out[:, 0].astype(int), r2_base=out[:, 1], inc=out[:, 2], exceed=out[:, 3].astype(int),
               null_mean=out[:, 4], pval=pval, fdr=fdr, selected=(fdr < FDR_LEVEL) & (out[:, 2] > 0),
               worst_ref=worst, k_block=int(block.shape[1]),
               q_base=1 + (0 if base_extra is None else int(base_extra.shape[1])))
    floor = int(np.sum(out[:, 3] == 0))
    kstar = int(np.ceil(G / (FDR_LEVEL * (B + 1))))
    n_bh = int((fdr < FDR_LEVEL).sum())
    log("  %s: block %d cols over a %d-col baseline; selected %d/%d; BH rejections before the sign "
        "filter %d; genes at the permutation floor %d (k* = %d); batched-vs-reference worst %.2e; %.0f s"
        % (tag, res["k_block"], res["q_base"], int(res["selected"].sum()), G, n_bh, floor, kstar, worst,
           time.time() - t0))
    res.update(n_floor=floor, kstar=kstar, n_bh=n_bh)
    return res


# ------------------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------------------
def load_inputs(ctx, key):
    c = ctx.cohorts[key]
    path = ctx.P(c["inputs"])
    if sha256(path) != c["inputs_sha"]:
        fatal("FATAL: %s is not the frozen file" % c["inputs"])
    with np.load(path, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    # a synthetic file can never enter a real run, and the self-test can never read a real file
    if bool(d["is_synthetic"]) != ctx.synthetic:
        fatal("FATAL: %s: synthetic flag %s in a %s run" % (c["inputs"], bool(d["is_synthetic"]),
                                                            "synthetic" if ctx.synthetic else "real"))
    for k in ("patients", "genes", "label_operator"):
        d[k] = d[k].astype(str)
    d["pcs"] = d["pcs"].astype(np.float64)        # exact upcast, as the published hstack does
    d["rna_pcs20"] = d["rna_pcs20"].astype(np.float64)
    n, G = d["protein"].shape
    if (n, G) != (c["n_disc"], c["n_fit"]) or d["pcs"].shape != (n, 20) or d["rna_pcs20"].shape != (n, 20):
        fatal("FATAL: unexpected shapes in %s" % c["inputs"])
    if int(d["ra_n_perm"]) != ctx.B or int(d["ra_seed"]) != 0 or int(d["cv_folds"]) != tk.CV_FOLDS \
            or float(d["ridge_alpha"]) != tk.RIDGE_ALPHA:
        fatal("FATAL: estimator constants in the inputs file differ from the published ones")
    if len(set(d["genes"])) != G or np.sum(~np.isnan(d["protein"]), axis=0).min() < MIN_FIT:
        fatal("FATAL: duplicate genes or a gene with fewer than 30 protein values")
    # PS 2: covariates are built in the published case order (sorted case identifiers)
    if list(d["patients"]) != sorted(set(d["patients"])):
        fatal("FATAL: discovery patients are not sorted and unique")
    pub_sel = (d["pub_fdr"] < FDR_LEVEL) & (d["pub_incremental_r2"] > 0)
    if int(pub_sel.sum()) != c["n_sel"]:
        fatal("FATAL: the published selection set is not the frozen size")
    d["pub_sel"] = pub_sel
    return d


def load_conf(ctx, key, genes_only=False):
    c = ctx.cohorts[key]
    path = ctx.P(c["conf"])
    if sha256(path) != c["conf_sha"]:
        fatal("FATAL: %s is not the frozen file" % c["conf"])
    with np.load(path, allow_pickle=False) as z:
        genes = z["genes"].astype(str)
        out = dict(genes=genes, sel_mask=z["sel_mask"].astype(bool), M={})
        if not genes_only:
            for name, k in c["arms"]:
                M = z[k]
                if M.shape != (c["n_univ"], ctx.B + 1) or not np.isfinite(M).all():
                    fatal("FATAL: unexpected confirmatory matrix %s" % k)
                out["M"][name] = M
    if len(genes) != c["n_univ"] or len(set(genes)) != len(genes) or int(out["sel_mask"].sum()) != c["n1"]:
        fatal("FATAL: unexpected confirmatory universe in %s" % c["conf"])
    return out


def check_universe(d, conf):
    """Every confirmatory gene must be a discovery fit gene, so every selector has a discovery
    increment for it (rho is taken over the universe)."""
    fit = set(d["genes"])
    missing = [g for g in conf["genes"] if g not in fit]
    if missing:
        fatal("FATAL: %d confirmatory genes are not discovery fit genes (%s ...)" % (len(missing), missing[:3]))


def published_perms(ctx, d):
    return tk.build_perms(d["protein"].shape[0], ctx.B, "ra", 0)


def operator_dummies(labels, n_cols, n_covered):
    """PS 2: one indicator for each of the 12 largest operator levels with at least 3 cases;
    every other case, including an unlabelled one (''), has an all-zero row.  The counts of PS 2
    (12 / 73 UCEC, 12 / 75 LUAD) are binding and asserted; a tie at the cut-off is fatal.
    Unlike residual_analysis_sitepack.py no level is dropped for full rank: PS 2 does not drop
    one, and with an all-zero reference group there is nothing to drop."""
    labels = np.asarray(labels).astype(str)
    cnt = collections.Counter(x for x in labels if x != "")
    eligible = sorted(((k, m) for k, m in cnt.items() if m >= MIN_OP_GROUP), key=lambda kv: (-kv[1], kv[0]))
    if len(eligible) > MAX_OP_LEVELS and eligible[MAX_OP_LEVELS][1] == eligible[MAX_OP_LEVELS - 1][1]:
        fatal("FATAL: tie at the operator-indicator cut-off")
    levels = [k for k, _ in eligible[:MAX_OP_LEVELS]]
    D = np.zeros((len(labels), len(levels)))
    for j, k in enumerate(levels):
        D[:, j] = (labels == k).astype(float)
    if D.shape[1] != n_cols or int(D.sum()) != n_covered:
        fatal("FATAL: operator indicators %d cols / %d cases differ from PS 2 (%d / %d)"
              % (D.shape[1], int(D.sum()), n_cols, n_covered))
    return D, levels


def check_operator(ctx, key, d):
    """PS 2 binding counts on the real labels (12 indicators covering 73 of 100 / 75 of 105 cases;
    a tie at the cut-off is fatal).  Label counts only.  Returns (indicators, levels, info)."""
    c = ctx.cohorts[key]
    D, levels = operator_dummies(d["label_operator"], c["n_op_dummies"], c["n_op_covered"])
    labels = np.asarray(d["label_operator"]).astype(str)
    cnt = collections.Counter(labels.tolist())
    info = dict(n_indicators=int(D.shape[1]), n_covered=int(D.sum()), n_cases=int(len(labels)),
                n_unlabelled=int(cnt.get("", 0)), levels=[[str(k), int(cnt[k])] for k in levels])
    log("OP %s: %d indicators cover %d of %d cases (%d unlabelled); levels and counts %s"
        % (key, info["n_indicators"], info["n_covered"], info["n_cases"], info["n_unlabelled"], info["levels"]))
    return D, levels, info


def _fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_sbc(ctx, key, universe):
    """PS 3 S_BC (reported, not read): genes with fdr_batchresid < 0.05 and incremental_r2_batchresid > 0
    in server_export/results/<cohort>__sitepack_operator.csv.  Fatal unless the set has exactly the
    frozen size (201 UCEC, 726 LUAD) and, for LUAD, unless the file has the sha256 recorded in the
    registered result json (key sitepack_pub_key; read in full from that json and compared; its first
    8 hex digits must be sitepack_prefix).  The number of S_BC genes in the confirmatory universe is
    logged (expected 198 / 721; not fatal).  rho (PS 5) is taken from incremental_r2_batchresid over the
    universe if every universe gene has a finite value in the table, otherwise there is no rho (and
    which genes are missing is logged).  Only the published table is read."""
    c = ctx.cohorts[key]
    path = ctx.P(c["sitepack"])
    if not os.path.isfile(path):
        fatal("REFUSED: operator-axis table %s not found" % c["sitepack"])
    h = sha256(path)
    prefix = c.get("sitepack_prefix")
    if prefix:
        if not h.startswith(prefix):
            fatal("REFUSED: %s has sha256 %s..., not the recorded %s..." % (c["sitepack"], h[:8], prefix))
        pk = c.get("sitepack_pub_key")
        if pk:
            with open(ctx.P(c["published"]), encoding="utf-8") as fh:
                want = json.load(fh).get(pk)
            if not (isinstance(want, str) and HEX64.fullmatch(want) and want.lower().startswith(prefix)):
                fatal("REFUSED: %s has no usable %s beginning %s" % (c["published"], pk, prefix))
            if h != want.lower():
                fatal("REFUSED: %s does not have the sha256 recorded under %s in %s"
                      % (c["sitepack"], pk, c["published"]))
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    for col in ("gene", "fdr_batchresid", "incremental_r2_batchresid"):
        if not rows or col not in rows[0]:
            fatal("REFUSED: %s has no column %s" % (c["sitepack"], col))
    fdr, inc = {}, {}
    for r in rows:
        if r["gene"] in inc:
            fatal("REFUSED: %s lists gene %s twice" % (c["sitepack"], r["gene"]))
        fdr[r["gene"]] = _fnum(r["fdr_batchresid"])
        inc[r["gene"]] = _fnum(r["incremental_r2_batchresid"])
    S = set(g for g in inc if fdr[g] < FDR_LEVEL and inc[g] > 0)     # NaN compares False: not selected
    if len(S) != c["n_bc"]:
        fatal("REFUSED: S_BC of %s has %d genes, not the frozen %d" % (key, len(S), c["n_bc"]))
    uni = list(universe)
    n_uni = len(S & set(uni))
    missing = [g for g in uni if g not in inc]
    nonfin = [g for g in uni if g in inc and not np.isfinite(inc[g])]
    if missing or nonfin:
        xs, note = None, ("no rho: %d universe genes absent from the table, %d without a finite "
                          "incremental_r2_batchresid" % (len(missing), len(nonfin)))
    else:
        xs, note = {g: inc[g] for g in uni}, ""
    log("S_BC %s: %d genes (frozen %d); %d in the confirmatory universe (expected %s; not fatal); %s"
        % (key, len(S), c["n_bc"], n_uni, c.get("n_bc_universe") or "n/a", note or "rho from incremental_r2_batchresid"))
    return dict(S=S, xs=xs, n=len(S), n_in_universe=n_uni, expected_in_universe=c.get("n_bc_universe"),
                rho_available=xs is not None, rho_note=note, sitepack_sha256=h, n_rows=len(rows))


def sbc_info(sbc):
    return {k: sbc[k] for k in ("n", "n_in_universe", "expected_in_universe", "rho_available", "rho_note",
                                "sitepack_sha256", "n_rows")}


def gaussian_block(n, K):
    """PS 5 diagnostic reference: a standard Gaussian block of the size of the morphology components."""
    return np.random.RandomState(GAUSS_SEED).randn(n, K)


def load_scores(path, patients):
    """Discovery composition scores from Stage H0 (PS 2, PS 4): tab-separated, header
    `patient stroma immune proliferation`, optional leading '#' lines (kept as comments), floats
    at full precision, exactly one row per discovery patient.  CRLF line ends are accepted
    because the file crosses from the cluster through Windows."""
    if not os.path.exists(path):
        fatal("REFUSED: discovery score table %s not found" % path)
    with open(path, "rb") as fh:
        text = fh.read().decode("utf-8-sig")
    lines = [ln[:-1] if ln.endswith("\r") else ln for ln in text.split("\n")]
    i, comments = 0, []
    while i < len(lines) and lines[i].startswith("#"):
        comments.append(lines[i])
        i += 1
    if i >= len(lines) or lines[i].split("\t") != ["patient"] + list(SCORE_COLS):
        fatal("REFUSED: unexpected header in %s" % path)
    tab = {}
    for ln in lines[i + 1:]:
        if ln == "":
            continue
        f = ln.split("\t")
        if len(f) != 4 or f[0] in tab:
            fatal("REFUSED: malformed or duplicate row in %s: %r" % (path, ln[:60]))
        try:
            tab[f[0]] = [float(x) for x in f[1:]]
        except ValueError:
            fatal("REFUSED: non-numeric score in %s: %r" % (path, ln[:60]))
    if set(tab) != set(patients):
        fatal("REFUSED: %s does not hold exactly the cohort's discovery patients (%d rows, %d missing, %d extra)"
              % (path, len(tab), len(set(patients) - set(tab)), len(set(tab) - set(patients))))
    C = np.asarray([tab[p] for p in patients], dtype=np.float64)
    if not np.isfinite(C).all():
        fatal("REFUSED: non-finite composition score in %s" % path)
    return C, comments


# ------------------------------------------------------------------------------------------
# sets (PS 3)
# ------------------------------------------------------------------------------------------
def rank_order(genes, pval, inc):
    """PS 3: p ascending, then increment descending, then gene symbol ascending."""
    return sorted(range(len(genes)), key=lambda i: (float(pval[i]), -float(inc[i]), str(genes[i])))


def lowest_r2_order(genes, r2):
    """PS 3 T_G: published own-mRNA R2 ascending, then gene symbol."""
    return sorted(range(len(genes)), key=lambda i: (float(r2[i]), str(genes[i])))


def _top(order, G, n, S):
    T = np.zeros(G, bool)
    T[order[:n]] = True
    rank = np.empty(G, dtype=np.int64)
    rank[order] = np.arange(1, G + 1)
    return dict(T=T, rank=rank, S=S)


def build_sets(d, res):
    """T_X (fixed size n = |S_M| among fit genes, no sign filter) and S_X for every analysis that
    is not closed (res[X] None: closed, no sets); the rule reference T_M (PS 3: the same rule on
    the published analysis M: published p ascending, published increment descending, symbol);
    and T_G."""
    genes = d["genes"]
    G = len(genes)
    n = int(d["pub_sel"].sum())
    out = {"M": _top(rank_order(genes, d["pub_pval"], d["pub_incremental_r2"]), G, n, d["pub_sel"].copy())}
    for X in ANALYSES:
        r = res[X]
        out[X] = None if r is None else _top(rank_order(genes, r["pval"], r["inc"]), G, n, r["selected"].copy())
    TG = np.zeros(G, bool)
    TG[lowest_r2_order(genes, d["pub_r2_rna"])[:n]] = True
    out["T_G"] = TG
    return out


def published_exceed(pub_pval, B):
    """#{null >= observed} of the published analysis M, recovered from its p = (1 + e) / (1 + B)."""
    e = np.rint(np.asarray(pub_pval, float) * (B + 1)).astype(int) - 1
    if np.max(np.abs((1 + e) / float(B + 1) - pub_pval)) > 1e-12 or e.min() < 0:
        fatal("FATAL: the published p-values are not of the form (1 + e) / %d" % (B + 1))
    return e


def t_summary(genes, T, rank, pval, inc, exceed, fdr, sm, tm, unis):
    """PS 7.3 tab:c1covsets, one T set: overlap with S_M and with T_M, the largest p in the set,
    the Benjamini-Hochberg q at that rank (BH over all fit genes of the analysis: the q of the
    last gene of the set), the number of its genes at the permutation floor, whether the floor
    alone filled the set (floor count of the analysis >= n), and the number of its genes with a
    non-positive discovery increment."""
    idx = np.flatnonzero(T)
    n = len(idx)
    last = int(np.flatnonzero(rank == n)[0])
    if not T[last]:
        fatal("FATAL: the last-ranked gene is not in its set")
    tset = set(genes[idx])
    largest_p = float(np.max(pval[idx]))
    if largest_p != float(pval[last]):
        fatal("FATAL: the set is not the n smallest p-values")
    n_floor_analysis = int(np.sum(exceed == 0))
    return dict(status="ok", n=n, n_in_universe=len(tset & unis), overlap_S_M=len(tset & sm),
                overlap_T_M=len(tset & tm), largest_p=largest_p, last_gene=str(genes[last]),
                bh_q_at_rank=float(fdr[last]), n_floor_in_set=int(np.sum(exceed[idx] == 0)),
                n_floor_analysis=n_floor_analysis, floor_alone_filled_set=bool(n_floor_analysis >= n),
                n_nonpositive_increment=int(np.sum(inc[idx] <= 0)))


def s_summary(S, fdr, exceed, B, unis, sm, genes):
    """PS 7.3 tab:c1covsets, one S set: its size, rejections before the sign filter, the floor
    count of the analysis and k*."""
    G = len(S)
    n_bh = int(np.sum(fdr < FDR_LEVEL))
    sset = set(genes[S])
    return dict(status="ok", n=int(S.sum()), n_in_universe=len(sset & unis), overlap_S_M=len(sset & sm),
                n_bh_before_sign_filter=n_bh, n_floor=int(np.sum(exceed == 0)),
                kstar=int(np.ceil(G / (FDR_LEVEL * (B + 1)))), note=s_note(n_bh, int(S.sum())))


SET_FIELDS = ("increment", "r2_base", "exceed", "pval", "fdr")


def write_sets(path, d, res, sets):
    """One row per discovery fit gene; floats written with repr so they read back exactly.  The
    columns of a closed analysis are empty strings."""
    cols = ["gene", "n", "pub_r2_rna", "pub_incremental_r2", "pub_pval", "pub_fdr", "S_M", "T_M", "T_M_rank", "T_G"]
    for X in ANALYSES:
        cols += ["%s_%s" % (X, f) for f in SET_FIELDS] + ["S_" + X, "T_" + X, "T_%s_rank" % X]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for i, g in enumerate(d["genes"]):
            row = [g, str(int(d["pub_n"][i])), repr(float(d["pub_r2_rna"][i])),
                   repr(float(d["pub_incremental_r2"][i])), repr(float(d["pub_pval"][i])),
                   repr(float(d["pub_fdr"][i])), str(int(d["pub_sel"][i])), str(int(sets["M"]["T"][i])),
                   str(int(sets["M"]["rank"][i])), str(int(sets["T_G"][i]))]
            for X in ANALYSES:
                r, s = res[X], sets[X]
                if r is None:
                    row += [""] * (len(SET_FIELDS) + 3)
                    continue
                row += [repr(float(r["inc"][i])), repr(float(r["r2_base"][i])), str(int(r["exceed"][i])),
                        repr(float(r["pval"][i])), repr(float(r["fdr"][i])),
                        str(int(s["S"][i])), str(int(s["T"][i])), str(int(s["rank"][i]))]
            w.writerow(row)
    return cols


def read_sets(path):
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], rows[1:]
    ix = {h: j for j, h in enumerate(head)}
    genes = [r[0] for r in body]

    def col_set(name):
        return set(r[0] for r in body if r[ix[name]] == "1")

    def col_float(name):
        return {r[0]: float(r[ix[name]]) for r in body}

    out = dict(genes=genes, S_M=col_set("S_M"), T_G=col_set("T_G"), pub_inc=col_float("pub_incremental_r2"),
               pub_r2_rna=col_float("pub_r2_rna"))
    out["M"] = dict(T=col_set("T_M"), S=out["S_M"], xs=out["pub_inc"])
    for X in ANALYSES:
        empty = [r[ix["S_" + X]] == "" for r in body]
        if all(empty):
            out[X] = None                                 # closed analysis: no sets
            continue
        if any(empty):
            fatal("FATAL: %s: analysis %s is partly empty" % (path, X))
        fdr = col_float(X + "_fdr")
        out[X] = dict(S=col_set("S_" + X), T=col_set("T_" + X), xs=col_float(X + "_increment"),
                      pval=col_float(X + "_pval"), fdr=fdr,
                      n_bh=int(sum(1 for v in fdr.values() if v < FDR_LEVEL)))
    return out


def s_note(n_bh, n_sel):
    """PS 3: how an FDR-defined set is reported when it is empty."""
    if n_bh == 0:
        return "empty at B = 1,000, limited by the permutation floor"
    if n_sel == 0:
        return "rejections but no positive increment"
    return ""


# ------------------------------------------------------------------------------------------
# cells (PS 5)
# ------------------------------------------------------------------------------------------
def cell(M, universe, set_genes, xs_by_gene=None, background=None, return_vec=False):
    """AUC of the set against its background on one confirmatory matrix.  background=None: every
    other gene of the universe (the C1-LUAD H2 rule).  The null is every column after column 0."""
    sel = np.fromiter((g in set_genes for g in universe), bool, len(universe))
    if background is None:
        keep = np.ones(len(universe), bool)
    else:
        keep = sel | np.fromiter((g in background for g in universe), bool, len(universe))
    n1, n0 = int(sel[keep].sum()), int((~sel[keep]).sum())
    if n1 < MIN_SET or n0 < MIN_SET:
        return dict(status="untestable", n1=n1, n0=n0)
    Mk = M[keep]
    if xs_by_gene is None:
        xs = np.arange(Mk.shape[0], dtype=float)
    else:
        xs = np.asarray([xs_by_gene[g] for g in universe[keep]], dtype=float)
    auc, rho, n1, n0 = assemble(Mk, xs, sel[keep])
    sd = float(auc[1:].std())
    out = dict(status="ok", n1=n1, n0=n0, n_null=int(len(auc) - 1), auc=float(auc[0]),
               null_mean=float(auc[1:].mean()), null_sd=sd,
               margin=float(auc[0] - auc[1:].mean()), z=float((auc[0] - auc[1:].mean()) / sd) if sd > 0 else None,
               p=emp_p(auc[0], auc[1:]), exceed=int(np.sum(auc[1:] >= auc[0])))
    if xs_by_gene is not None:
        out.update(rho=float(rho[0]), rho_null_mean=float(rho[1:].mean()), rho_null_sd=float(rho[1:].std()),
                   p_rho=emp_p(rho[0], rho[1:]))
    if return_vec:
        out["auc_vec"] = auc
    return out


PARTS = ("both", "S_M only", "T_X only")


def partition_cells(M, universe, sm, tx):
    """PS 5, reported and not read: genes in both S_M and T_X, in S_M only and in T_X only, each
    against the genes in neither.  If T_X has fewer than 50 genes in the universe the three parts
    are 'partition undefined'.  No rho: the parts have no single selector."""
    uni = set(universe)
    sm, tx = sm & uni, tx & uni
    parts = {"both": sm & tx, "S_M only": sm - tx, "T_X only": tx - sm}
    neither = uni - sm - tx
    if len(tx) < MIN_SET:
        return {k: dict(status="partition undefined", n1=len(s), n0=len(neither),
                        note="T_X has %d genes in the universe" % len(tx)) for k, s in parts.items()}
    return {k: cell(M, universe, s, background=neither) for k, s in parts.items()}


def d_change(auc_pub, auc_x):
    """PS 5: D_b = AUC_pub(b) - AUC_X(b) with the shared permutations, b = 0..B."""
    D = np.asarray(auc_pub, float) - np.asarray(auc_x, float)
    null = D[1:]
    sd = float(null.std())
    return dict(D_0=float(D[0]), D_null_mean=float(null.mean()), D_null_sd=sd,
                z_D=float((D[0] - null.mean()) / sd) if sd > 0 else None, p_D=emp_p(D[0], null))


def capacity_p(m_j, m_R):
    """PS 5: e = #{j : m_j <= m_R} over the 119 placebos; p_abs = (1 + e) / 120."""
    if len(m_j) != N_CAP:
        fatal("FATAL: capacity randomisation needs %d placebo margins, got %d" % (N_CAP, len(m_j)))
    e = int(sum(1 for m in m_j if m <= m_R))
    return dict(e=e, p_abs=(1 + e) / float(N_CAP + 1))


def tier(p):
    return "p<0.05" if p < 0.05 else ("0.05<=p<0.15" if p < 0.15 else "p>=0.15")


CELL_COLS = ["cohort", "family", "set", "matrix", "null", "background", "read", "status", "note",
             "n1", "n0", "n_null", "auc", "null_mean", "null_sd", "margin", "z", "p", "exceed", "tier",
             "rho", "rho_null_mean", "rho_null_sd", "p_rho"]
L2_EXTRA = ["class", "D_0", "D_null_mean", "D_null_sd", "z_D", "p_D", "f", "m_R", "e", "p_abs", "m_j"]


def make_row(cohort, family, set_name, matrix, null, read, r, background="rest of universe", note=""):
    row = dict(cohort=cohort, family=family, set=set_name, matrix=matrix, null=null, background=background,
               read=bool(read), note=note or r.get("note", ""))
    for k, v in r.items():
        if k not in ("auc_vec", "note"):
            row[k] = v
    if r.get("status") == "ok":
        row["tier"] = tier(r["p"])
    return row


def closed_cell(why):
    return dict(status="closed", note="closed (CLOSED_PARTS.txt: %s)" % why)


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, (bool, np.bool_)):
        return "True" if v else "False"
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    if isinstance(v, (list, tuple)):
        return ";".join(_fmt(x) for x in v)
    return str(v)


def write_cells(path, rows, cols):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([_fmt(r.get(k)) for k in cols])


def read_cells(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


# ------------------------------------------------------------------------------------------
# reading rules (PS 6)
# ------------------------------------------------------------------------------------------
def classify_cohort(c3_status, p_c3, r20_status, p_r, caps_status, p_abs):
    """PS 6.1 for one cohort; p values unrounded.  Statuses: 'ok' or 'closed' (any other status
    of a primary cell is fatal; after freezing nothing is 'not completed').

    C3 class: H if p < 0.05, N otherwise; 'closed' if A_C3 is closed.
    R20 class: 'closed' if A_R20 is closed; else H if p_R < 0.05; else 'closed' if the placebos
      are closed; else A if p_abs < 0.05; else U.
    Type: C if C3 class N; R if C3 class H or closed, and R20 class A; NC:C if C3 closed and R20
      class not A; NC:R if C3 class H and R20 closed; none otherwise.  Flagged = has a type (PS 6.4)."""
    for s in (c3_status, r20_status, caps_status):
        if s not in ("ok", "closed"):
            fatal("FATAL: unknown primary status %r" % s)
    c3 = "closed" if c3_status == "closed" else ("H" if p_c3 < ALPHA else "N")
    if r20_status == "closed":
        r20 = "closed"
    elif p_r < ALPHA:
        r20 = "H"
    elif caps_status == "closed":
        r20 = "closed"
    elif p_abs < ALPHA:
        r20 = "A"
    else:
        r20 = "U"
    if c3 == "N":
        typ = "C"
    elif r20 == "A":                               # c3 is H or closed here
        typ = "R"
    elif c3 == "closed":
        typ = "NC:C"
    elif r20 == "closed":                          # c3 is H here
        typ = "NC:R"
    else:
        typ = "none"
    return dict(C3_class=c3, R20_class=r20, type=typ, flagged=typ != "none", U=(r20 == "U"))


def case_of(flags):
    """PS 6.4: Case 0, 1 or 2 = number of flagged cohorts."""
    return int(sum(1 for f in flags if f))


def m4_count(entries):
    """PS 7.1 M4: k outranked of t that could be tested; closed and untestable cells are excluded
    from t.  entries are (name, status, p)."""
    tested = [(n, p) for n, s, p in entries if s == "ok"]
    out = [n for n, p in tested if p < ALPHA]
    return dict(k=len(out), t=len(tested), n_cells=len(entries), outranked=out,
                not_tested=[[n, s] for n, s, _ in entries if s != "ok"],
                cells=[dict(name=n, status=s, p=p) for n, s, p in entries])


# ------------------------------------------------------------------------------------------
# Stage H files (PS 4) and the manifest (PS 9)
# ------------------------------------------------------------------------------------------
def parse_manifest(path):
    """Only lines '<sha256> *<bare file name>' (blank lines are skipped; anything else, status
    lines included, is fatal).  Returns {name: sha256}."""
    hashed = {}
    with open(path, "rb") as fh:
        text = fh.read().decode("utf-8")
    for ln in text.split("\n"):
        line = ln.rstrip("\r")
        if not line.strip():
            continue
        m = re.fullmatch(r"([0-9a-fA-F]{64}) \*(\S+)", line)
        if not m:
            fatal("REFUSED: malformed manifest line %r (only '<sha256> *<file name>' lines are allowed)" % line)
        name = m.group(2)
        if "/" in name or "\\" in name:
            fatal("REFUSED: manifest names must be bare file names: %r" % name)
        if name in hashed:
            fatal("REFUSED: %s is listed twice in the manifest" % name)
        hashed[name] = m.group(1).lower()
    return hashed


def read_h_npz(path, genes_ref, ncols, n_cases):
    """Format check of one Stage H matrix (PS 4 "Files").  Returns (data or None, problems)."""
    probs = []
    try:
        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
    except Exception as e:  # an object array or a corrupt file: reported, never unpickled
        return None, ["cannot load without pickle: %s" % e]
    for k in ("genes", "M", "cases"):
        if k not in data:
            probs.append("key %s missing" % k)
    if probs:
        return None, probs
    g, M, cs = data["genes"], data["M"], data["cases"]
    if g.dtype.kind != "U" or g.ndim != 1:
        probs.append("genes is not a 1-d unicode array (%s)" % g.dtype)
    elif list(g) != list(genes_ref):
        probs.append("genes differ in value or order from the saved published matrix")
    if M.dtype != np.float64:
        probs.append("M is %s, not float64" % M.dtype)
    if M.shape != (len(genes_ref), ncols):
        probs.append("M has shape %s, expected %s" % (M.shape, (len(genes_ref), ncols)))
    elif not np.isfinite(M).all():
        probs.append("M has non-finite values")
    if cs.dtype.kind != "U" or cs.ndim != 1:
        probs.append("cases is not a 1-d unicode array (%s)" % cs.dtype)
    elif len(cs) != n_cases or list(cs) != sorted(set(cs)):
        probs.append("cases: %d entries, expected %d sorted unique" % (len(cs), n_cases))
    extra = sorted(set(data) - {"genes", "M", "cases"})
    if extra:
        probs.append("unexpected keys %s" % extra)
    return data, probs


def read_cov_record(path, n_cases, K, need_c3, need_r20):
    """Format check of c1_cov_confirmatory_covariates_<cohort>.npz (PS 4; shared interface): keys
    cases (sorted unique unicode), C3 (n x 3; absent or n x 0 allowed only if both C3 arms are
    closed), R20 (n x 20; may be absent only if every R20 arm is closed), pcs (n x K, the
    confirmatory morphology components), markers_json (JSON text).  Further keys are allowed.
    Returns (data or None, problems); data gains 'markers' (the parsed JSON)."""
    probs = []
    try:
        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
    except Exception as e:
        return None, ["cannot load without pickle: %s" % e]
    for k in ("cases", "pcs", "markers_json") + (("C3",) if need_c3 else ()) + (("R20",) if need_r20 else ()):
        if k not in data:
            probs.append("key %s missing" % k)
    if probs:
        return None, probs
    cs = data["cases"]
    if cs.dtype.kind != "U" or cs.ndim != 1 or len(cs) != n_cases or list(cs) != sorted(set(cs)):
        probs.append("cases: not %d sorted unique unicode entries" % n_cases)

    def fmat(name, shapes):
        a = data[name]
        if a.dtype != np.float64:
            probs.append("%s is %s, not float64" % (name, a.dtype))
        elif a.shape not in shapes:
            probs.append("%s has shape %s, expected %s" % (name, a.shape, " or ".join(str(s) for s in shapes)))
        elif not np.isfinite(a).all():
            probs.append("%s has non-finite values" % name)
    fmat("pcs", [(n_cases, K)])
    if "C3" in data:
        fmat("C3", [(n_cases, 3)] if need_c3 else [(n_cases, 3), (n_cases, 0)])
    if "R20" in data:
        fmat("R20", [(n_cases, 20)])
    mj = data["markers_json"]
    try:
        if mj.dtype.kind != "U" or mj.size != 1:
            raise ValueError("not a single unicode string")
        data["markers"] = json.loads(str(mj.reshape(-1)[0]))
    except Exception as e:
        probs.append("markers_json does not hold JSON text (%s)" % e)
    return data, probs


def stage_h_expected(ctx, key, closed):
    """(file name, number of columns, arm token) of every Stage H matrix of a cohort that is not
    closed."""
    out = [(stage_h_name(key, a), ctx.B + 1, a) for a in MAIN_ARMS if not arm_closed(closed, key, a)]
    if not arm_closed(closed, key, "A_R20cap"):
        out += [(cap_name(key, j), CAP_COLS + 1, "A_R20cap") for j in range(1, N_CAP + 1)]
    return out


def preflight_stage_h(ctx, closed, man):
    """PS 9 Stage L2: 'refuses to start unless every Stage H file it reads has the sha256 the
    manifest gives'.  For every file that is not closed: listed, present, hash, format (keys,
    dtypes, shapes, gene order), and one shared case list per cohort, the covariates record
    included.  Returns {cohort: {"hashes": {name: sha256}, "record": name or None}}."""
    plan = {}
    for key, c in ctx.cohorts.items():
        genes = load_conf(ctx, key, genes_only=True)["genes"]
        expected = stage_h_expected(ctx, key, closed)
        hashes, cases = {}, {}

        def verify(name):
            if name not in man:
                fatal("REFUSED: %s (not closed) is not listed in the Stage H manifest" % name)
            p = ctx.H(name)
            if not os.path.isfile(p):
                fatal("REFUSED: %s is listed in the manifest but missing from %s" % (name, ctx.stage_h_dir))
            h = sha256(p)
            if h != man[name]:
                fatal("REFUSED: %s does not have the sha256 the manifest gives" % name)
            return p, h
        for name, ncols, _ in expected:
            p, h = verify(name)
            data, probs = read_h_npz(p, genes, ncols, c["n_conf_cases"])
            if probs:
                fatal("REFUSED: %s: %s" % (name, "; ".join(probs)))
            hashes[name] = h
            cases[name] = tuple(data["cases"])
            del data
        rec = None
        rn = covariates_record_name(key)
        if expected or rn in man:
            p, h = verify(rn)
            need_c3 = not (arm_closed(closed, key, "A_C3") and arm_closed(closed, key, "V_C3"))
            need_r20 = not all(arm_closed(closed, key, a) for a in ("A_R20", "V_R20", "A_R20cap"))
            data, probs = read_cov_record(p, c["n_conf_cases"], c["K"], need_c3, need_r20)
            if probs:
                fatal("REFUSED: %s: %s" % (rn, "; ".join(probs)))
            hashes[rn] = h
            cases[rn] = tuple(data["cases"])
            rec = rn
        if len(set(cases.values())) > 1:
            fatal("REFUSED: the Stage H files of %s do not share one case list" % key)
        ignored = sorted(n for n in man if n.startswith("c1_cov_%s_" % key) and n not in hashes)
        log("pre-flight %s: %d Stage H files verified (manifest hash, format, one case list); closed arms %s; "
            "manifest entries not read: %d" % (key, len(hashes), sorted(a for k, a in closed["arms"] if k == key)
                                                or "none", len(ignored)))
        plan[key] = dict(hashes=hashes, record=rec)
    return plan


def load_verified_matrix(ctx, name, h, genes, ncols, n_cases):
    p = ctx.H(name)
    if sha256(p) != h:
        fatal("FATAL: %s changed after the pre-flight" % name)
    data, probs = read_h_npz(p, genes, ncols, n_cases)
    if probs:
        fatal("FATAL: %s: %s" % (name, "; ".join(probs)))
    return data["M"]


def check_stage_h(ctx, expect_cols, expect_cap_cols):
    """--check-stage-h: format only, no statistic.  Returns True if every present file passes."""
    ok = True
    closed = no_closed()
    if os.path.exists(ctx.closed_path):
        try:
            closed = parse_closed(ctx.closed_path, ctx.cohorts)
            log("CLOSED_PARTS.txt: %s" % (closed["tokens"] or "nothing closed"))
        except Fatal as e:
            ok = False
            log("CLOSED_PARTS.txt FAILED: %s" % e)
    for key, c in ctx.cohorts.items():
        d, conf = load_inputs(ctx, key), load_conf(ctx, key, genes_only=True)
        sp = ctx.H(scores_name(key))
        if os.path.exists(sp):
            try:
                _, comments = load_scores(sp, d["patients"])
                log("%s: %s OK (%d comment lines)" % (key, scores_name(key), len(comments)))
            except Fatal as e:
                ok = False
                log("%s: %s FAILED: %s" % (key, scores_name(key), e))
        else:
            log("%s: %s absent" % (key, scores_name(key)))
        cases_seen = {}
        names = [(n, expect_cols if a != "A_R20cap" else expect_cap_cols)
                 for n, _, a in stage_h_expected(ctx, key, closed)]
        absent = []
        for name, ncols in names:
            p = ctx.H(name)
            if not os.path.exists(p):
                absent.append(name)
                continue
            data, probs = read_h_npz(p, conf["genes"], ncols, c["n_conf_cases"])
            if data is not None and "cases" in data and data["cases"].dtype.kind == "U":
                cases_seen[name] = tuple(data["cases"])
            if probs:
                ok = False
                log("%s: %s FAILED: %s" % (key, name, "; ".join(probs)))
            else:
                log("%s: %s OK (%d x %d)" % (key, name, len(conf["genes"]), ncols))
        if absent:
            log("%s: %d of %d expected matrices absent (%s%s)" % (key, len(absent), len(names), ", ".join(absent[:4]),
                                                               " ..." if len(absent) > 4 else ""))
        rp = ctx.H(covariates_record_name(key))
        if os.path.exists(rp):
            need_c3 = not (arm_closed(closed, key, "A_C3") and arm_closed(closed, key, "V_C3"))
            need_r20 = not all(arm_closed(closed, key, a) for a in ("A_R20", "V_R20", "A_R20cap"))
            data, probs = read_cov_record(rp, c["n_conf_cases"], c["K"], need_c3, need_r20)
            if data is not None and data["cases"].dtype.kind == "U":
                cases_seen[covariates_record_name(key)] = tuple(data["cases"])
            if probs:
                ok = False
                log("%s: %s FAILED: %s" % (key, covariates_record_name(key), "; ".join(probs)))
            else:
                log("%s: %s OK (keys %s)" % (key, covariates_record_name(key), sorted(data)))
        else:
            log("%s: %s absent" % (key, covariates_record_name(key)))
        if len(set(cases_seen.values())) > 1:
            ok = False
            log("%s: FAILED: the case lists differ between Stage H files" % key)
    mp_ = ctx.H(MANIFEST_NAME)
    if os.path.exists(mp_):
        try:
            man = parse_manifest(mp_)
            for name, h in sorted(man.items()):
                p = ctx.H(name)
                if not os.path.exists(p):
                    ok = False
                    log("manifest: %s listed but absent" % name)
                elif sha256(p) != h:
                    ok = False
                    log("manifest: %s sha256 MISMATCH" % name)
            log("manifest: %d files; own sha256 %s" % (len(man), sha256(mp_)))
        except Fatal as e:
            ok = False
            log("manifest FAILED: %s" % e)
    else:
        log("manifest absent")
    log("check-stage-h: %s" % ("ALL PRESENT FILES PASS" if ok else "FAILED"))
    return ok


# ------------------------------------------------------------------------------------------
# validation (published quantities only)
# ------------------------------------------------------------------------------------------
def validate_v1(ctx, key, d, workers):
    res = run_block(d, d["pcs"], None, published_perms(ctx, d), workers, "V1 %s" % key)
    dif = dict(increment=float(np.max(np.abs(res["inc"] - d["pub_incremental_r2"]))),
               r2_rna=float(np.max(np.abs(res["r2_base"] - d["pub_r2_rna"]))),
               pval=float(np.max(np.abs(res["pval"] - d["pub_pval"]))),
               fdr=float(np.max(np.abs(res["fdr"] - d["pub_fdr"]))),
               n=int(np.max(np.abs(res["n"] - d["pub_n"]))))
    same_set = bool(np.array_equal(res["selected"], d["pub_sel"]))
    ok = (dif["r2_rna"] <= TOL_PUB and dif["increment"] <= TOL_INC and dif["pval"] <= TOL_PFDR
          and dif["fdr"] <= TOL_PFDR and dif["n"] == 0 and same_set)
    log("V1 %s: max abs differences %s; selection set identical: %s; selected %d -> %s"
        % (key, dif, same_set, int(res["selected"].sum()), "VALIDATED" if ok else "FAILED"))
    return dict(ok=ok, max_abs_diff=dif, selection_identical=same_set, n_selected=int(res["selected"].sum()),
                batched_vs_reference_worst=res["worst_ref"])


def validate_v1b(ctx, key, d):
    """PS 9 V1b: R2(own mRNA + rna_pcs20) by the reference kernel against r2_rna_pcs20 (1e-9 in both
    organs), and the observed increment of the 20 morphology components over [own mRNA, rna_pcs20]
    by the batched kernel (the code path of M|R20) against incr_over_pcs20 (1e-9 in UCEC, 1e-5 in
    LUAD), for the published-set genes.  The largest differences are recorded."""
    c = ctx.cohorts[key]
    tol_inc = c["tol_v1b_inc"]
    path = ctx.P(c["txperm"])
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    gi = {g: i for i, g in enumerate(d["genes"])}
    n = d["protein"].shape[0]
    ident = np.arange(n, dtype=np.int64)[None, :]
    dr2, dinc, dref, nmis, missing = [], [], 0.0, 0, []
    for r in rows:
        g = r["gene"]
        if g not in gi:
            missing.append(g)
            continue
        j = gi[g]
        v = ~np.isnan(d["protein"][:, j])
        y = d["protein"][v, j]
        folds = tk.folds_for(int(v.sum()))
        base = np.hstack([d["mrna_log2tpm1"][v, j].reshape(-1, 1), d["rna_pcs20"][v]])
        r2_ref = tk.r2_of(y, tk.oof_ridge(base, y, folds))
        _, inc = batched_increments(y, base, d["pcs"], v, ident)
        dref = max(dref, abs(ref_increment(y, base, d["pcs"][v], folds)[1] - inc[0]))
        dr2.append(abs(r2_ref - float(r["r2_rna_pcs20"])))
        dinc.append(abs(inc[0] - float(r["incr_over_pcs20"])))
        nmis += int(int(v.sum()) != int(float(r["n"])))
    pub_set = set(d["genes"][d["pub_sel"]])
    csv_set = set(r["gene"] for r in rows)
    dr2, dinc = np.asarray(dr2), np.asarray(dinc)
    res = dict(n_rows=len(rows), n_compared=int(len(dr2)), genes_not_in_inputs=missing[:20],
               n_genes_not_in_inputs=len(missing), published_set_genes_not_in_csv=len(pub_set - csv_set),
               csv_genes_not_in_published_set=len(csv_set - pub_set), n_mismatch=nmis,
               tol_r2_rna_pcs20=TOL_INC, tol_incr_over_pcs20=tol_inc,
               max_abs_diff_r2_rna_pcs20=float(dr2.max()) if len(dr2) else None,
               p999_abs_diff_r2_rna_pcs20=float(np.percentile(dr2, 99.9)) if len(dr2) else None,
               max_abs_diff_incr_over_pcs20=float(dinc.max()) if len(dinc) else None,
               p999_abs_diff_incr_over_pcs20=float(np.percentile(dinc, 99.9)) if len(dinc) else None,
               n_genes_over_1e9_r2=int(np.sum(dr2 > 1e-9)), n_genes_over_1e9_inc=int(np.sum(dinc > 1e-9)),
               n_genes_over_tol_inc=int(np.sum(dinc > tol_inc)),
               batched_vs_reference_worst=float(dref))
    res["ok"] = bool(len(dr2) > 0 and not missing and nmis == 0 and res["published_set_genes_not_in_csv"] == 0
                     and res["max_abs_diff_r2_rna_pcs20"] <= TOL_INC
                     and res["max_abs_diff_incr_over_pcs20"] <= tol_inc and dref <= TOL_REF)
    log("V1b %s: %d genes; max abs diff r2_rna_pcs20 %.3e (tol %.0e; 99.9th pct %.3e), incr_over_pcs20 %.3e "
        "(tol %.0e; 99.9th pct %.3e); genes over 1e-9 (R2 / increment): %d / %d; n mismatches %d; "
        "published-set genes missing from the csv %d; batched-vs-reference %.2e -> %s"
        % (key, res["n_compared"], res["max_abs_diff_r2_rna_pcs20"] or 0.0, TOL_INC,
           res["p999_abs_diff_r2_rna_pcs20"] or 0.0, res["max_abs_diff_incr_over_pcs20"] or 0.0, tol_inc,
           res["p999_abs_diff_incr_over_pcs20"] or 0.0, res["n_genes_over_1e9_r2"], res["n_genes_over_1e9_inc"],
           nmis, res["published_set_genes_not_in_csv"], dref, "VALIDATED" if res["ok"] else "FAILED"))
    return res


def validate_v2(ctx, key, d, conf):
    c = ctx.cohorts[key]
    with open(ctx.P(c["published"]), encoding="utf-8") as fh:
        pub = json.load(fh)
    xs = dict(zip(d["genes"], d["pub_incremental_r2"]))
    sm = set(d["genes"][d["pub_sel"]])
    mask = np.fromiter((g in sm for g in conf["genes"]), bool, len(conf["genes"]))
    ok = bool(np.array_equal(mask, conf["sel_mask"]))
    rows = {}
    for name, _ in c["arms"]:
        r = cell(conf["M"][name], conf["genes"], sm, xs_by_gene=xs)
        if c["published_kind"] == "ucec":
            w = pub["arms"][name]
            want = dict(auc=w["auc"]["obs"], null_mean=w["auc"]["null_mean"], null_sd=w["auc"]["null_sd"],
                        p=w["auc"]["p_one_sided"], rho=w["rho"]["obs"], p_rho=w["rho"]["p_one_sided"])
        else:
            w = pub["h1"]
            want = dict(auc=w["auc_obs"], null_mean=w["auc_null_mean"], null_sd=w["auc_null_sd"],
                        p=w["p_auc"], rho=w["rho_obs"], p_rho=w["p_rho"])
        worst = max(abs(r[k] - want[k]) for k in want)
        good = worst <= TOL_AUC and (r["n1"], r["n0"]) == (c["n1"], c["n0"])
        ok = ok and good
        rows[name] = dict(got={k: r[k] for k in want}, n1=r["n1"], n0=r["n0"], worst_abs_diff=worst, ok=good)
        log("V2 %s %s: AUC %.10f null %.10f SD %.10f p %.6f rho %.10f; n %d/%d; worst diff %.1e -> %s"
            % (key, name, r["auc"], r["null_mean"], r["null_sd"], r["p"], r["rho"], r["n1"], r["n0"], worst,
               "VALIDATED" if good else "FAILED"))
    return dict(ok=ok, mask_identical=bool(np.array_equal(mask, conf["sel_mask"])), arms=rows)


def validate_v3(ctx, key, d, C):
    """PS 9 V3: the transferred discovery scores reproduce the published composition controls
    (results__composition_controls_<cohort>.csv) for every published-set gene: the gain of each
    score alone and of the three together over own mRNA.  These involve no morphology."""
    c = ctx.cohorts[key]
    with open(ctx.P(c["compctl"]), encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    cols = [("%s_gain" % f, [j]) for j, f in enumerate(SCORE_COLS)] + [("ALL_gain", [0, 1, 2])]
    gi = {g: i for i, g in enumerate(d["genes"])}
    worst, m, skipped = 0.0, 0, 0
    for r in rows:
        if r["gene"] not in gi or any(r.get(k, "") == "" for k, _ in cols):
            skipped += 1
            continue
        g = gi[r["gene"]]
        v = ~np.isnan(d["protein"][:, g])
        y = d["protein"][v, g]
        folds = tk.folds_for(int(v.sum()))
        xr = d["mrna_log2tpm1"][v, g].reshape(-1, 1)
        r2_r = tk.r2_of(y, tk.oof_ridge(xr, y, folds))
        for k, which in cols:
            r2_c = tk.r2_of(y, tk.oof_ridge(np.hstack([xr, C[v][:, which]]), y, folds))
            worst = max(worst, abs((r2_c - r2_r) - float(r[k])))
        m += 1
    pub_set = set(d["genes"][d["pub_sel"]])
    covered = len(pub_set & set(r["gene"] for r in rows))
    ok = m == len(rows) and m > 0 and covered == len(pub_set) and worst <= TOL_INC
    log("V3 %s: %d of %d published composition-control genes compared (%d published-set genes covered of %d), "
        "worst abs difference %.2e -> %s" % (key, m, len(rows), covered, len(pub_set), worst,
                                             "VALIDATED" if ok else "FAILED"))
    return dict(ok=bool(ok), n_genes=m, n_rows=len(rows), n_skipped=skipped, published_set_covered=covered,
                worst_abs_diff=worst)


def gate_v3(ctx):
    out = {}
    for key in ctx.cohorts:
        d = load_inputs(ctx, key)
        C, comments = load_scores(ctx.H(scores_name(key)), d["patients"])
        out[key] = validate_v3(ctx, key, d, C)
        out[key]["score_file_sha256"] = sha256(ctx.H(scores_name(key)))
        out[key]["score_file_comments"] = comments
    return out


# ------------------------------------------------------------------------------------------
# stages
# ------------------------------------------------------------------------------------------
class _Tee(object):
    """Copies everything printed during a stage into <stage>.log (PS 9: runs are kept with
    their log)."""

    def __init__(self, a, b):
        self.a, self.b = a, b

    def write(self, s):
        self.a.write(s)
        self.b.write(s)

    def flush(self):
        self.a.flush()
        self.b.flush()


def refuse_if_not_runnable(ctx, stage):
    """Refusals about the stage folder itself (they create nothing)."""
    final = os.path.join(ctx.out_dir, stage)
    partial = final + ".partial"
    if os.path.exists(os.path.join(final, "result_%s.json" % stage)):
        fatal("REFUSED: Stage %s has completed (%s); a completed stage is never rerun" % (stage, final))
    if os.path.exists(final):
        fatal("REFUSED: %s exists without result_%s.json; resolve by hand" % (final, stage))
    if os.path.exists(os.path.join(partial, "result_%s.json" % stage)):
        fatal("REFUSED: %s holds a result file but was not renamed; resolve by hand" % partial)


def run_stage(ctx, stage, body, workers, record=None):
    """Run safety (PS 9 "Interrupted runs").  `body(out_dir)` writes the stage's outputs into
    out_dir and returns a dict for result_<stage>.json.  No limit on attempts: a leftover
    <stage>.partial is kept as <stage>_failed<k> and the stage runs again."""
    refuse_if_not_runnable(ctx, stage)
    final = os.path.join(ctx.out_dir, stage)
    partial = final + ".partial"
    if os.path.exists(partial):
        k = 1
        while os.path.exists("%s_failed%d" % (final, k)):
            k += 1
        failed = "%s_failed%d" % (final, k)
        os.rename(partial, failed)
        log("previous interrupted attempt kept as %s; Stage %s runs again" % (failed, stage))
    os.makedirs(partial)
    logf = open(os.path.join(partial, "%s.log" % stage), "w", encoding="utf-8", newline="\n")
    old = sys.stdout
    sys.stdout = _Tee(old, logf)
    t0 = time.time()
    try:
        log("Stage %s start; script %s; output %s" % (stage, ctx.script_path, partial))
        res = body(partial)
        outputs = {}
        for fn in sorted(os.listdir(partial)):
            if fn != "%s.log" % stage:      # the log is still being written; it is hashed in the results record
                outputs[fn] = sha256(os.path.join(partial, fn))
        res.update(record or {})
        res.update(stage=stage, outputs=outputs, script_sha256=sha256(ctx.script_path),
                   sidecar_sha256=sha256(ctx.sidecar), prespec_sha256=sha256(ctx.prespec),
                   kernel_sha256=sha256(ctx.kernel_path), closed_parts_sha256=sha256(ctx.closed_path),
                   env=env_record(), workers=workers,
                   attempts_before_this=sorted(fn for fn in os.listdir(ctx.out_dir)
                                               if re.fullmatch(r"%s_failed\d+" % re.escape(stage), fn)),
                   minutes=round((time.time() - t0) / 60.0, 2),
                   finished_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        with open(os.path.join(partial, "result_%s.json" % stage), "w", encoding="utf-8", newline="\n") as fh:
            json.dump(res, fh, indent=1, sort_keys=True)
        log("Stage %s complete in %.1f min; outputs %s" % (stage, res["minutes"], sorted(outputs)))
    except BaseException:
        log("Stage %s FAILED:\n%s" % (stage, traceback.format_exc()))
        raise
    finally:
        sys.stdout = old
        logf.close()
    os.rename(partial, final)
    log("Stage %s folder renamed to %s" % (stage, final))
    return final


def preflight_l1(ctx, closed):
    """L1 pre-flight: inputs, confirmatory matrices, the operator-indicator counts of PS 2, the
    batch-corrected sets, V2 on the saved matrices and, for every cohort whose score table is read,
    the score table and V3.  Any failure is a refusal."""
    pre = {}
    for key in ctx.cohorts:
        d = load_inputs(ctx, key)
        conf = load_conf(ctx, key)
        check_universe(d, conf)
        OP, levels, op_info = check_operator(ctx, key, d)             # PS 2 binding counts (task item 3)
        sbc = load_sbc(ctx, key, conf["genes"])
        v2 = validate_v2(ctx, key, d, conf)
        if not v2["ok"]:
            fatal("REFUSED: V2: the reference cells of %s do not reproduce the published readings" % key)
        C, comments, v3, ssha = None, None, None, None
        if needs_scores(closed, key):
            C, comments = load_scores(ctx.H(scores_name(key)), d["patients"])
            v3 = validate_v3(ctx, key, d, C)
            if not v3["ok"]:
                fatal("REFUSED: V3 failed for %s: the discovery scores are not the published ones" % key)
            ssha = sha256(ctx.H(scores_name(key)))
        else:
            log("%s: discovery analyses C3 and M|C3 are closed; no score table is read" % key)
        pre[key] = dict(d=d, conf=conf, C=C, comments=comments, v2=v2, v3=v3, score_sha=ssha,
                        OP=OP, levels=levels, op_info=op_info, sbc=sbc)
    return pre


def stage_l1(ctx, out_dir, workers, closed, pre):
    rows, summary = [], {}
    for key in ctx.cohorts:
        P = pre[key]
        d, conf, C = P["d"], P["conf"], P["C"]
        perms = published_perms(ctx, d)
        OP, levels, sbc = P["OP"], P["levels"], P["sbc"]          # checked in the pre-flight
        spec = {"C3": (C, None), "R20": (d["rna_pcs20"], None), "OP": (OP, None),
                "M|C3": (d["pcs"], C), "M|R20": (d["pcs"], d["rna_pcs20"])}
        res = {}
        for X in ANALYSES:
            if an_closed(closed, key, X):
                res[X] = None
                log("  %s %s: closed (CLOSED_PARTS.txt); no sets" % (key, X))
                continue
            blk, base = spec[X]
            res[X] = run_block(d, blk, base, perms, workers, "%s %s" % (key, X))
        sets = build_sets(d, res)
        write_sets(os.path.join(out_dir, "sets_%s.csv" % key), d, res, sets)

        genes, uni = d["genes"], conf["genes"]
        unis = set(uni)
        n = int(d["pub_sel"].sum())
        sm = set(genes[d["pub_sel"]])
        tm = set(genes[sets["M"]["T"]])
        T = {X: (set(genes[sets[X]["T"]]) if sets[X] is not None else None) for X in ANALYSES}
        S = {X: (set(genes[res[X]["selected"]]) if res[X] is not None else None) for X in ANALYSES}
        TG = set(genes[sets["T_G"]])
        xs = {X: (dict(zip(genes, res[X]["inc"])) if res[X] is not None else None) for X in ANALYSES}
        xs_pub = dict(zip(genes, d["pub_incremental_r2"]))

        # PS 7.3 tab:c1covsets
        exc_m = published_exceed(d["pub_pval"], ctx.B)
        ssum = {"T_M": t_summary(genes, sets["M"]["T"], sets["M"]["rank"], d["pub_pval"], d["pub_incremental_r2"],
                                 exc_m, d["pub_fdr"], sm, tm, unis),
                "S_M": s_summary(d["pub_sel"], d["pub_fdr"], exc_m, ctx.B, unis, sm, genes)}
        for X in ANALYSES:
            r = res[X]
            if r is None:
                ssum["T_" + X] = ssum["S_" + X] = dict(status="closed")
                continue
            ssum["T_" + X] = t_summary(genes, sets[X]["T"], sets[X]["rank"], r["pval"], r["inc"], r["exceed"],
                                       r["fdr"], sm, tm, unis)
            ssum["T_" + X].update(k_block=r["k_block"], q_base=r["q_base"], batched_vs_reference_worst=r["worst_ref"])
            ssum["S_" + X] = s_summary(r["selected"], r["fdr"], r["exceed"], ctx.B, unis, sm, genes)
        ssum["T_G"] = dict(status="ok", n=len(TG), n_in_universe=len(TG & unis), overlap_S_M=len(TG & sm),
                           overlap_T_M=len(TG & tm),
                           note="T_G is chosen by published own-mRNA R2, not by a p-value: the p-based "
                                "columns do not apply")
        summary[key] = dict(n_cases=int(d["protein"].shape[0]), n_fit=int(len(genes)), n_set=n,
                            operator_levels=levels, operator_indicators=P["op_info"], S_BC=sbc_info(sbc),
                            score_file_comments=P["comments"],
                            score_file_sha256=P["score_sha"], v3=P["v3"], v2=P["v2"], sets=ssum,
                            closed_analyses=sorted(X for X in ANALYSES if an_closed(closed, key, X)),
                            note_M_R20=("analysis M|R20 uses the permutation sequence of the main analysis "
                                        "(RandomState(0), 1,000 sequential draws) for every fit gene, so its "
                                        "p-values and floor counts are not those of the published "
                                        "transcriptome-baseline table (PS 3)"))

        for arm, _ in ctx.cohorts[key]["arms"]:
            M = conf["M"][arm]
            unres = arm == "unrestricted"
            rows.append(make_row(key, "reference", "S_M", "published", arm, False,
                                 cell(M, uni, sm, xs_by_gene=xs_pub)))
            rows.append(make_row(key, "rule reference", "T_M", "published", arm, False,
                                 cell(M, uni, tm, xs_by_gene=xs_pub)))
            for X in FORWARD_X:
                r = closed_cell("analysis %s %s" % (key, X)) if T[X] is None else cell(M, uni, T[X], xs_by_gene=xs[X])
                rows.append(make_row(key, "forward" if unres else "forward (descriptive)", "T_" + X, "published",
                                     arm, unres, r))
            # T_G has no discovery increment of its own: rho is not computed for it
            rows.append(make_row(key, "gene property", "T_G", "published", arm, False, cell(M, uni, TG)))
            for X in FORWARD_X:
                names = {"both": "S_M&T_%s" % X, "S_M only": "S_M-T_%s" % X, "T_X only": "T_%s-S_M" % X}
                if T[X] is None:
                    parts = {p: closed_cell("analysis %s %s" % (key, X)) for p in PARTS}
                else:
                    parts = partition_cells(M, uni, sm, T[X])
                for part in PARTS:
                    rows.append(make_row(key, "partition", names[part], "published", arm, False, parts[part],
                                         background="neither S_M nor T_%s" % X))
            if unres:
                # PS 5: S_C3, S_R20, S_OP on the published matrix, unrestricted arm only
                for X in FORWARD_X:
                    if S[X] is None:
                        r = closed_cell("analysis %s %s" % (key, X))
                    else:
                        r = cell(M, uni, S[X], xs_by_gene=xs[X])
                        r["note"] = ssum["S_" + X]["note"]
                    rows.append(make_row(key, "FDR-defined", "S_" + X, "published", arm, False, r))
                # PS 5: S_BC (reported, not read) on the published unrestricted arm; rho only if the
                # table has a finite value for every universe gene
                rows.append(make_row(key, "batch-corrected", "S_BC", "published", arm, False,
                                     cell(M, uni, sbc["S"], xs_by_gene=sbc["xs"]), note=sbc["rho_note"]))
    write_cells(os.path.join(out_dir, "cells_L1.csv"), rows, CELL_COLS)
    return dict(cohorts=summary, n_cells=len(rows))


def verify_l1(ctx, manifest_sha, closed):
    """Stage L2 refuses unless Stage L1 completed and its outputs have the recorded hashes, and
    unless Stage L1 ran on the same manifest and the same closed parts."""
    l1 = os.path.join(ctx.out_dir, "L1")
    rp = os.path.join(l1, "result_L1.json")
    if not os.path.exists(rp):
        fatal("REFUSED: Stage L1 has not completed (%s missing)" % rp)
    with open(rp, encoding="utf-8") as fh:
        r1 = json.load(fh)
    need = ["cells_L1.csv"] + ["sets_%s.csv" % k for k in ctx.cohorts]
    for fn in need:
        if fn not in r1.get("outputs", {}):
            fatal("REFUSED: %s is not recorded in result_L1.json" % fn)
    for fn, h in r1["outputs"].items():
        p = os.path.join(l1, fn)
        if not os.path.exists(p) or sha256(p) != h:
            fatal("REFUSED: Stage L1 output %s missing or changed since Stage L1" % fn)
    if r1.get("manifest_sha256") != manifest_sha.strip().lower():
        fatal("REFUSED: Stage L1 ran on another Stage H manifest")
    if sorted(r1.get("closed_parts", [])) != sorted(closed["tokens"]):
        fatal("REFUSED: Stage L1 ran with other closed parts")
    log("Stage L1 outputs verified against result_L1.json (%d files)" % len(r1["outputs"]))
    return r1, sha256(rp)


def verify_manifest(ctx, manifest_sha):
    mp_ = ctx.H(MANIFEST_NAME)
    if not manifest_sha:
        fatal("REFUSED: --manifest-sha is required (the value in the Stage H record)")
    if not os.path.exists(mp_):
        fatal("REFUSED: Stage H manifest %s not found" % mp_)
    if sha256(mp_) != manifest_sha.strip().lower():
        fatal("REFUSED: the manifest's sha256 is not the value given with --manifest-sha")
    man = parse_manifest(mp_)
    log("Stage H manifest verified: %d files" % len(man))
    return man


def stage_l2(ctx, out_dir, workers, closed, plan, r1_sha):
    l1 = os.path.join(ctx.out_dir, "L1")
    cells_l1 = read_cells(os.path.join(l1, "cells_L1.csv"))
    rows, per, diag = [], {}, {}
    matched_entries, forward_entries = [], []
    for key, c in ctx.cohorts.items():
        conf = load_conf(ctx, key)
        S = read_sets(os.path.join(l1, "sets_%s.csv" % key))
        uni = conf["genes"]
        sbc = load_sbc(ctx, key, uni)
        fit = set(S["genes"])
        if any(g not in fit for g in uni):
            fatal("FATAL: the sets file of %s does not cover the confirmatory universe" % key)
        for X in ANALYSES:
            if (S[X] is None) != an_closed(closed, key, X):
                fatal("FATAL: sets file of %s and CLOSED_PARTS.txt disagree on analysis %s" % (key, X))
        hashes = plan[key]["hashes"]
        n_cases = c["n_conf_cases"]
        H, st = {}, {}
        for arm in MAIN_ARMS:
            if arm_closed(closed, key, arm):
                H[arm], st[arm] = None, "closed"
            else:
                name = stage_h_name(key, arm)
                H[arm], st[arm] = load_verified_matrix(ctx, name, hashes[name], uni, ctx.B + 1, n_cases), "ok"
        rec = None
        if plan[key]["record"]:
            rn = plan[key]["record"]
            if sha256(ctx.H(rn)) != hashes[rn]:
                fatal("FATAL: %s changed after the pre-flight" % rn)
            need_c3 = not (arm_closed(closed, key, "A_C3") and arm_closed(closed, key, "V_C3"))
            need_r20 = not all(arm_closed(closed, key, a) for a in ("A_R20", "V_R20", "A_R20cap"))
            rec, probs = read_cov_record(ctx.H(rn), n_cases, c["K"], need_c3, need_r20)
            if probs:
                fatal("FATAL: %s: %s" % (rn, "; ".join(probs)))

        sm, xs_pub, tm = S["S_M"], S["pub_inc"], S["M"]["T"]
        ref = cell(conf["M"]["unrestricted"], uni, sm, xs_by_gene=xs_pub, return_vec=True)

        def rd(family, set_name, sel, arm, read, use_rho=True, note=""):
            """sel: (gene set, xs) or None for a closed analysis (with its name in `note`)."""
            if st[arm] == "closed":
                r = closed_cell("arm %s %s" % (key, arm))
            elif sel is None:
                r = closed_cell(note)
            else:
                r = cell(H[arm], uni, sel[0], xs_by_gene=sel[1] if use_rho else None, return_vec=True)
                if note and r.get("status") == "ok":
                    r["note"] = note
            row = make_row(key, family, set_name, arm, "unrestricted", read, r)
            rows.append(row)
            return row, r

        def tsel(X):
            return None if S[X] is None else (S[X]["T"], S[X]["xs"])

        def ssel(X):
            return None if S[X] is None else (S[X]["S"], S[X]["xs"])

        # primary cells, with the change D from the published matrix (PS 5)
        prim = {}
        for arm in ("A_C3", "A_R20"):
            row, r = rd("primary", "S_M", (sm, xs_pub), arm, True)
            if r["status"] == "ok":
                dd = d_change(ref["auc_vec"], r["auc_vec"])
                dd["f"] = r["margin"] / ref["margin"]
                row.update(dd)
            elif r["status"] != "closed":
                fatal("FATAL: primary cell %s %s is %s" % (key, arm, r["status"]))
            prim[arm] = (row, r)
        # capacity randomisation (PS 5): m_j from each placebo (columns 0 and 1-25), m_R from columns
        # 0 and 1-25 of A_R20
        caps_status = "closed" if arm_closed(closed, key, "A_R20cap") else "ok"
        m_j = []
        for j in range(1, N_CAP + 1):
            mname = "A_R20cap_%03d" % j
            if caps_status == "closed":
                rows.append(make_row(key, "capacity", "S_M", mname, "unrestricted (25 draws)", False,
                                     closed_cell("arm %s A_R20cap" % key)))
                continue
            name = cap_name(key, j)
            Mc = load_verified_matrix(ctx, name, hashes[name], uni, CAP_COLS + 1, n_cases)
            r = cell(Mc, uni, sm, xs_by_gene=xs_pub)
            del Mc
            if r["status"] != "ok":
                fatal("FATAL: capacity cell %s %s is %s" % (key, name, r["status"]))
            rows.append(make_row(key, "capacity", "S_M", mname, "unrestricted (25 draws)", False, r))
            m_j.append(r["margin"])
        m_R, cap = None, None
        if st["A_R20"] == "ok":
            m_R = cell(H["A_R20"][:, :CAP_COLS + 1], uni, sm)["margin"]
            if caps_status == "ok":
                cap = capacity_p(m_j, m_R)
        # reverse (read), matched (read), covariate-to-covariate, rule reference, FDR-defined and
        # gene-property cells (PS 5)
        rev = {X: rd("reverse", "S_M", (sm, xs_pub), "V_" + X, True)[1] for X in ("C3", "R20")}
        mat = {X: rd("matched", "T_M|" + X, tsel("M|" + X), "A_" + X, True,
                     note="" if S["M|" + X] is not None else "analysis %s M|%s" % (key, X))[1]
               for X in ("C3", "R20")}
        for X in ("C3", "R20"):
            rd("covariate-to-covariate", "T_" + X, tsel(X), "V_" + X, False,
               note="" if S[X] is not None else "analysis %s %s" % (key, X))
        for X in ("C3", "R20"):
            rd("rule reference", "T_M", (tm, xs_pub), "A_" + X, False)
        for X in ("C3", "R20"):
            Y = "M|" + X
            note = s_note(S[Y]["n_bh"], len(S[Y]["S"])) if S[Y] is not None else "analysis %s %s" % (key, Y)
            rd("FDR-defined", "S_" + Y, ssel(Y), "A_" + X, False, note=note)
        for X in ("C3", "R20"):
            rd("gene property", "T_G", (S["T_G"], None), "V_" + X, False, use_rho=False)
        for X in ("C3", "R20"):
            rd("batch-corrected", "S_BC", (sbc["S"], sbc["xs"]), "A_" + X, False, note=sbc["rho_note"])

        pc3, pr = prim["A_C3"][1], prim["A_R20"][1]
        cls = classify_cohort(pc3["status"], pc3.get("p"), pr["status"], pr.get("p"), caps_status,
                              cap["p_abs"] if cap else None)
        prim["A_C3"][0]["class"] = cls["C3_class"]
        prim["A_R20"][0]["class"] = cls["R20_class"]
        if m_R is not None:
            prim["A_R20"][0].update(m_R=m_R, m_j=list(m_j))
        if cap:
            prim["A_R20"][0].update(cap)

        # diagnostic (PS 5, not read): canonical correlations of the confirmatory morphology
        # components with C3 and with R20, from the covariates record
        dg = dict(record=plan[key]["record"], cancorr_pcs_C3=None, cancorr_pcs_R20=None,
                  cancorr_gauss_C3=None, cancorr_gauss_R20=None, gauss_seed=GAUSS_SEED, markers=None)
        if rec is not None:
            dg["markers"] = rec["markers"]
            dg["record_keys"] = sorted(k for k in rec if k != "markers")
            if "C3" in rec and rec["C3"].shape[1] == 3:
                dg["cancorr_pcs_C3"] = [float(x) for x in cancorr(rec["pcs"], rec["C3"])]
            if "R20" in rec:
                dg["cancorr_pcs_R20"] = [float(x) for x in cancorr(rec["pcs"], rec["R20"])]
            dg["K"] = int(rec["pcs"].shape[1])
            # reference: a standard Gaussian block of the size of the morphology components (PS 5)
            Gref = gaussian_block(int(rec["pcs"].shape[0]), dg["K"])
            if "C3" in rec and rec["C3"].shape[1] == 3:
                dg["cancorr_gauss_C3"] = [float(x) for x in cancorr(Gref, rec["C3"])]
            if "R20" in rec:
                dg["cancorr_gauss_R20"] = [float(x) for x in cancorr(Gref, rec["R20"])]
        diag[key] = dg

        def dstats(arm):
            return {k: prim[arm][0].get(k) for k in ("D_0", "D_null_mean", "D_null_sd", "z_D", "p_D", "f")}

        def pstats(r):
            return {k: r.get(k) for k in ("status", "n1", "n0", "auc", "null_mean", "null_sd", "margin", "z", "p",
                                          "rho", "p_rho")}
        per[key] = dict(
            label=c["label"], **cls,
            e=cap["e"] if cap else None, p_abs=cap["p_abs"] if cap else None,
            C3=dict(pstats(pc3), D=dstats("A_C3")),
            R20=dict(pstats(pr), D=dstats("A_R20"), margin_all_1000_columns=pr.get("margin"),
                     m_R_25_columns=m_R, m_j=list(m_j), placebos=caps_status,
                     e=cap["e"] if cap else None, p_abs=cap["p_abs"] if cap else None),
            published=dict(auc=ref["auc"], null_mean=ref["null_mean"], null_sd=ref["null_sd"], margin=ref["margin"],
                           p=ref["p"]),
            reverse={X: dict(status=r["status"], p=r.get("p"), margin=r.get("margin"),
                             outranked=(r["p"] < ALPHA) if r["status"] == "ok" else None) for X, r in rev.items()},
            markers_confirmatory=dg["markers"], diagnostic=dg, S_BC=sbc_info(sbc),
            closed=dict(arms=sorted(a for k, a in closed["arms"] if k == key),
                        analyses=sorted(a for k, a in closed["analyses"] if k == key)))
        for X in ("C3", "R20"):
            r = mat[X]
            matched_entries.append(("%s T_M|%s on A_%s" % (key, X, X), r["status"], r.get("p")))
    for r in cells_l1:
        if r["family"] == "forward" and r["read"] == "True":
            forward_entries.append(("%s %s on the published matrix" % (r["cohort"], r["set"]), r["status"],
                                    float(r["p"]) if r["p"] else None))
    if len(forward_entries) != 6 or len(matched_entries) != 4:
        fatal("FATAL: expected six forward and four matched read cells, got %d and %d"
              % (len(forward_entries), len(matched_entries)))
    classes = dict(rule="PS 6.1-6.4; p unrounded; outranked means p < 0.05; a cohort is flagged if it has a type",
                   cohorts=per, case=case_of(v["flagged"] for v in per.values()),
                   flagged=[k for k, v in per.items() if v["flagged"]],
                   M4=dict(forward=m4_count(forward_entries), matched=m4_count(matched_entries)),
                   closed_parts=closed["tokens"])
    write_cells(os.path.join(out_dir, "cells_L2.csv"), rows, CELL_COLS + L2_EXTRA)
    with open(os.path.join(out_dir, "classes.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(classes, fh, indent=1, sort_keys=True)
    for k, v in per.items():
        log("%s: C3 class %s (p %r); R20 class %s (p_R %r, e %r, p_abs %r); type %s; flagged %s"
            % (k, v["C3_class"], v["C3"]["p"], v["R20_class"], v["R20"]["p"], v["e"], v["p_abs"], v["type"],
               v["flagged"]))
    log("case %d; M4 forward %d of %d, matched %d of %d" % (classes["case"], classes["M4"]["forward"]["k"],
                                                            classes["M4"]["forward"]["t"],
                                                            classes["M4"]["matched"]["k"],
                                                            classes["M4"]["matched"]["t"]))
    return dict(n_cells=len(rows), stage_h_files={n: h for key in plan for n, h in plan[key]["hashes"].items()},
                result_L1_sha256=r1_sha, case=classes["case"], diagnostic=diag)


def do_stage(ctx, stage, workers, manifest_sha):
    """Pre-flight (every check that only refuses; it creates nothing), then the stage itself."""
    refuse_if_not_runnable(ctx, stage)
    closed, _ = check_sidecar(ctx)
    pinned = check_pinned(ctx, closed)
    log("PINNED_HOSTS.txt: %s" % (pinned or "no host pinned (every cohort fully closed)"))
    man = verify_manifest(ctx, manifest_sha)
    msha = manifest_sha.strip().lower()
    record = dict(manifest_sha256=msha, closed_parts=closed["tokens"], pinned_hosts=pinned)
    if stage == "L1":
        pre = preflight_l1(ctx, closed)
        log("pre-flight L1 passed")
        return run_stage(ctx, "L1", lambda od: stage_l1(ctx, od, workers, closed, pre), workers, record)
    if stage != "L2":
        fatal("FATAL: unknown stage %r" % stage)
    _, r1_sha = verify_l1(ctx, msha, closed)
    plan = preflight_stage_h(ctx, closed, man)
    log("pre-flight L2 passed")
    return run_stage(ctx, "L2", lambda od: stage_l2(ctx, od, workers, closed, plan, r1_sha), workers, record)


# ------------------------------------------------------------------------------------------
# self-test (synthetic only)
# ------------------------------------------------------------------------------------------
def _expect_fatal(fn, *a, **k):
    try:
        fn(*a, **k)
    except Fatal:
        return True
    raise AssertionError("expected a refusal from %s" % getattr(fn, "__name__", fn))


def _tree(path):
    """Every file and folder under `path` (to show that a refusal created nothing)."""
    out = []
    for dp, dn, fn in os.walk(path):
        out += [os.path.join(dp, x) for x in dn + fn]
    return sorted(out)


def _selftest_kernel():
    from scipy.stats import mannwhitneyu
    rng = np.random.RandomState(1)
    n, G, k = 60, 30, 7
    mrna = rng.randn(n, G)
    block = rng.randn(n, k)
    block[:, 0] = (rng.rand(n) < 0.1).astype(float)            # a sparse dummy-like column
    labels = np.array(["L%d" % (i % 9) for i in range(n - 2)] + ["", "X"])
    dummies = np.zeros((n, 8))
    for j in range(8):                                          # a pure indicator block (OP-like)
        dummies[:, j] = (labels == "L%d" % j).astype(float)
    extra3 = rng.randn(n, 3)
    extra20 = rng.randn(n, 20)
    prot = 0.5 * mrna + 0.4 * block[:, [1]] + rng.randn(n, G)
    prot[rng.rand(n, G) < 0.1] = np.nan
    perms = tk.build_perms(n, 25, "ra", 0)
    idx = np.vstack([np.arange(n)[None, :], perms])
    worst = 0.0
    for g in range(G):
        v = ~np.isnan(prot[:, g])
        y = prot[v, g]
        folds = tk.folds_for(int(v.sum()))
        xr = mrna[v, g].reshape(-1, 1)
        for base in (xr, np.hstack([xr, extra3[v]]), np.hstack([xr, extra20[v]])):
            for blk in (block, dummies):
                r2b, inc = batched_increments(y, base, blk, v, idx)
                for r in range(idx.shape[0]):
                    rb, ri = ref_increment(y, base, blk[idx[r]][v], folds)
                    worst = max(worst, abs(ri - inc[r]), abs(rb - r2b))
        nv, r2r, r2both, pinc = tk.published_increment(mrna[:, g], prot[:, g], block, v)
        worst = max(worst, abs(pinc - batched_increments(y, xr, block, v, idx[:1])[1][0]))
        null = tk.published_null_increments(mrna[:, g], prot[:, g], block, perms, v)
        worst = max(worst, float(np.max(np.abs(null - batched_increments(y, xr, block, v, idx)[1][1:]))))
    assert worst <= 1e-10, worst
    # run_block: the process pool and the in-process path give identical numbers
    d = dict(protein=prot, mrna_log2tpm1=mrna)
    a = run_block(d, block, extra3, perms, 1, "selftest serial")
    b = run_block(d, block, extra3, perms, 2, "selftest pool")
    for kk in ("inc", "r2_base", "pval", "fdr"):
        assert np.array_equal(a[kk], b[kk]), kk
    # bh_fdr against a direct implementation, with ties
    p = rng.rand(500)
    p[:40] *= 1e-4
    p[40:60] = p[40]
    m = len(p)
    direct = np.array([min(1.0, min(p[p >= pi] * m / np.array([np.sum(p <= pj) for pj in p[p >= pi]]))) for pi in p])
    assert np.max(np.abs(bh_fdr(p) - direct)) <= 1e-12
    # assemble against scipy's Mann-Whitney statistic, and cell() on top of it
    M = rng.randn(300, 11)
    sel = rng.rand(300) < 0.3
    auc, rho, n1, n0 = assemble(M, rng.randn(300), sel)
    for j in range(M.shape[1]):
        u = mannwhitneyu(M[sel, j], M[~sel, j]).statistic / (n1 * n0)
        assert abs(u - auc[j]) <= 1e-12
    uni = np.array(["g%03d" % i for i in range(300)])
    c = cell(M, uni, set(uni[sel]))
    assert abs(c["auc"] - auc[0]) <= 1e-12 and c["n1"] == n1 and c["p"] == emp_p(auc[0], auc[1:])
    assert cell(M, uni, set(uni[:10]))["status"] == "untestable"
    assert cell(M, uni, set(uni[:260]))["status"] == "untestable"      # background below 50
    print("  kernel, BH, assemble: OK (batched vs reference worst %.1e)" % worst)


def _selftest_cancorr():
    """cancorr against a direct computation: sqrt of the eigenvalues of
    Saa^-1 Sab Sbb^-1 Sba, on blocks of the confirmatory shapes (K = 28 / 22 with 3 and 20)."""
    rng = np.random.RandomState(11)
    for n, K, q in ((138, 28, 3), (112, 22, 20), (60, 5, 8)):
        A = rng.randn(n, K)
        Bm = 0.6 * A[:, :q] @ rng.randn(q, q) + rng.randn(n, q) if q <= K else rng.randn(n, q)
        Bm[:, 0] += 0.8 * A[:, 0]
        got = cancorr(A, Bm)
        Ac, Bc = A - A.mean(0), Bm - Bm.mean(0)
        Saa, Sbb, Sab = Ac.T @ Ac, Bc.T @ Bc, Ac.T @ Bc
        Mx = np.linalg.solve(Saa, Sab) @ np.linalg.solve(Sbb, Sab.T)
        ev = np.sort(np.clip(np.real(np.linalg.eigvals(Mx)), 0, None))[::-1][:min(K, q)]
        direct = np.sqrt(ev)
        assert got.shape == (min(K, q),) and np.max(np.abs(got - direct)) <= TOL_CANCORR, (n, K, q)
    # a zero-SD column is handled as in the published recipe (no division by zero)
    A = rng.randn(40, 4)
    A[:, 2] = 1.0
    assert np.isfinite(cancorr(A, rng.randn(40, 3))).all()
    # PS 5 diagnostic reference: a standard Gaussian block RandomState(42).randn(n, K) with C3 and with R20
    for n, K in ((138, 28), (112, 22)):
        G = gaussian_block(n, K)
        assert G.shape == (n, K) and np.array_equal(G, np.random.RandomState(42).randn(n, K))
        assert np.array_equal(G, gaussian_block(n, K))                 # deterministic
        for q in (3, 20):
            cg = cancorr(G, rng.randn(n, q))
            assert cg.shape == (q,) and np.all(cg >= 0) and np.all(cg <= 1 + 1e-9)
    print("  cancorr (copied from transcriptome_baseline.py) equals the direct computation; "
          "Gaussian reference block RandomState(42).randn(n, K): OK")


def _selftest_sets(tmp):
    # ranking rule with ties: p asc, increment desc, gene symbol asc
    genes = np.array(["b", "a", "c", "d", "e", "f"])
    p = np.array([0.01, 0.01, 0.01, 0.002, 0.5, 0.01])
    inc = np.array([0.2, 0.2, 0.3, 0.1, 0.9, -0.4])
    assert [genes[i] for i in rank_order(genes, p, inc)] == ["d", "c", "a", "b", "f", "e"]
    r2 = np.array([0.3, 0.1, 0.1, -0.2, 0.5, 0.1])
    assert [genes[i] for i in lowest_r2_order(genes, r2)] == ["d", "a", "c", "f", "b", "e"]
    # operator indicators: unlabelled case gets no indicator; tie at the cut-off is fatal
    labels = np.array(["A"] * 5 + ["B"] * 4 + ["C"] * 3 + ["D"] * 2 + [""] * 3)
    D, lv = operator_dummies(labels, 3, 12)
    assert lv == ["A", "B", "C"] and D[-3:].sum() == 0
    lab13 = np.array(sum((["L%02d" % i] * 3 for i in range(13)), []))
    _expect_fatal(operator_dummies, lab13, 12, 36)
    # partition backgrounds on a known split
    uni = np.array(["g%03d" % i for i in range(300)])
    rng = np.random.RandomState(3)
    M = rng.randn(300, 21)
    sm = set(uni[:120])
    tx = set(uni[60:180])
    parts = partition_cells(M, uni, sm, tx)
    keep_bg = np.arange(180, 300)
    for part, idx in (("both", np.arange(60, 120)), ("S_M only", np.arange(0, 60)), ("T_X only", np.arange(120, 180))):
        r = parts[part]
        rows = np.concatenate([idx, keep_bg])
        mask = np.zeros(len(rows), bool)
        mask[:len(idx)] = True
        auc = assemble(M[rows], np.arange(len(rows), dtype=float), mask)[0]
        assert r["status"] == "ok" and (r["n1"], r["n0"]) == (60, 120) and abs(r["auc"] - auc[0]) <= 1e-15, part
    und = partition_cells(M, uni, sm, set(uni[290:]) | {"not_in_universe"})
    assert all(v["status"] == "partition undefined" for v in und.values())
    # write / read round trip of the sets file, with one closed analysis; T_M by the PS 3 rule
    G = 80
    g = np.array(["G%03d" % i for i in range(G)])
    B = 999
    pub_exc = rng.randint(0, 6, G)
    d = dict(genes=g, pub_n=np.full(G, 40.0), pub_r2_rna=rng.randn(G), pub_incremental_r2=rng.randn(G) * 1e-3,
             pub_pval=(1 + pub_exc) / float(B + 1))
    d["pub_fdr"] = bh_fdr(d["pub_pval"])
    d["pub_sel"] = (d["pub_fdr"] < 0.3) & (d["pub_incremental_r2"] > 0)
    n = int(d["pub_sel"].sum())
    res = {}
    for X in ANALYSES:
        exc = rng.randint(0, 5, G)
        res[X] = dict(inc=rng.randn(G) * 10.0 ** rng.randint(-12, 3, G), r2_base=rng.randn(G), exceed=exc,
                      pval=(1 + exc) / float(B + 1))
        res[X]["fdr"] = bh_fdr(res[X]["pval"])
        res[X]["selected"] = (res[X]["fdr"] < 0.05) & (res[X]["inc"] > 0)
    res["OP"] = None                                                     # a closed analysis
    sets = build_sets(d, res)
    want_tm = set(g[rank_order(g, d["pub_pval"], d["pub_incremental_r2"])[:n]])
    assert set(g[sets["M"]["T"]]) == want_tm and sets["OP"] is None
    path = os.path.join(tmp, "sets_roundtrip.csv")
    write_sets(path, d, res, sets)
    back = read_sets(path)
    assert back["genes"] == list(g) and back["S_M"] == set(g[d["pub_sel"]]) and back["T_G"] == set(g[sets["T_G"]])
    assert back["M"]["T"] == want_tm and back["OP"] is None
    assert all(back["pub_inc"][x] == float(y) for x, y in zip(g, d["pub_incremental_r2"]))
    for X in ANALYSES:
        if X == "OP":
            continue
        assert back[X]["S"] == set(g[res[X]["selected"]]) and back[X]["T"] == set(g[sets[X]["T"]])
        assert len(back[X]["T"]) == n
        assert all(back[X]["xs"][x] == float(y) for x, y in zip(g, res[X]["inc"]))
        assert back[X]["n_bh"] == int((res[X]["fdr"] < 0.05).sum())
    # set summaries (PS 7.3) on a constructed analysis with known values
    gg = np.array(["k%02d" % i for i in range(12)])
    exc = np.array([0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8])
    pv = (1 + exc) / 1001.0
    ic = np.array([0.3, -0.1, 0.2, 0.0, 0.5, -0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    fd = bh_fdr(pv)
    top = _top(rank_order(gg, pv, ic), 12, 5, np.zeros(12, bool))
    assert set(gg[top["T"]]) == {"k00", "k01", "k02", "k03", "k04"}
    s5 = t_summary(gg, top["T"], top["rank"], pv, ic, exc, fd, {"k00", "k07"}, {"k01", "k02", "k09"}, set(gg[:10]))
    # order within the set: k00 (0.3), k02 (0.2), k01 (-0.1) at the floor, then k04 (0.5), k03 (0.0)
    assert s5["last_gene"] == "k03" and s5["largest_p"] == 2 / 1001.0 and s5["bh_q_at_rank"] == fd[3]
    assert (s5["overlap_S_M"], s5["overlap_T_M"], s5["n_floor_in_set"], s5["n_floor_analysis"]) == (1, 2, 3, 3)
    assert s5["floor_alone_filled_set"] is False and s5["n_nonpositive_increment"] == 2 and s5["n_in_universe"] == 5
    top2 = _top(rank_order(gg, pv, ic), 12, 2, np.zeros(12, bool))
    s2 = t_summary(gg, top2["T"], top2["rank"], pv, ic, exc, fd, set(), set(), set(gg))
    assert s2["floor_alone_filled_set"] is True and s2["last_gene"] == "k02" and s2["n_nonpositive_increment"] == 0
    ss = s_summary(fd < 0.05, fd, exc, 1000, set(gg), set(), gg)
    assert (ss["n_bh_before_sign_filter"], ss["n_floor"], ss["kstar"]) == (int((fd < 0.05).sum()), 3, 1)
    assert s_summary(np.zeros(12, bool), np.ones(12), exc, 1000, set(), set(), gg)["note"].startswith("empty")
    assert np.array_equal(published_exceed(pv, 1000), exc)
    _expect_fatal(published_exceed, pv + 1e-4, 1000)
    print("  ranking, T_M, operator indicators, partitions, sets round trip (closed analysis), set summaries: OK")


def _selftest_parsers(tmp):
    pats = ["P1", "P2", "P3"]
    vals = [[0.1, -1.0 / 3, 2.5e-17], [1.0, 2.0, 3.0], [-0.5, 0.25, 1e300]]
    body = "#exported on the cluster\r\n#markers stroma 24\r\npatient\tstroma\timmune\tproliferation\r\n"
    for p, v in zip(pats, vals):
        body += p + "\t" + "\t".join(repr(x) for x in v) + "\r\n"
    good = os.path.join(tmp, "scores_good.tsv")
    with open(good, "wb") as fh:
        fh.write(body.encode("utf-8"))
    C, comments = load_scores(good, pats)
    assert np.array_equal(C, np.array(vals)) and len(comments) == 2

    def bad(text, patients=pats):
        p = os.path.join(tmp, "scores_bad.tsv")
        with open(p, "wb") as fh:
            fh.write(text.encode("utf-8"))
        _expect_fatal(load_scores, p, patients)
    bad(body, pats + ["P4"])                                         # a missing patient fails
    bad(body, pats[:2])                                              # an extra patient fails
    bad(body.replace("immune", "immun"))                             # header
    bad(body + "P1\t1\t2\t3\r\n")                                    # duplicate
    bad(body.replace("2.0", "nan"))                                  # non-finite
    bad(body.replace("\t3.0", ""))                                   # short row
    # CLOSED_PARTS.txt and PINNED_HOSTS.txt: the shared strict rules
    cohorts = {"ucec": {}, "luad": {}}
    cp = os.path.join(tmp, "closed.txt")

    def wr(text, path=cp):
        with open(path, "wb") as fh:
            fh.write(text.encode("utf-8"))
    wr("# header only\n# nothing closed\n")
    c0 = parse_closed(cp, cohorts)
    assert c0["arms"] == set() and c0["analyses"] == set() and c0["tokens"] == []
    wr("# h\narm ucec A_C3\narm luad A_R20\narm luad A_R20cap\nanalysis luad C3\n\nanalysis luad M|C3\n")
    c1 = parse_closed(cp, cohorts)
    assert c1["arms"] == {("ucec", "A_C3"), ("luad", "A_R20"), ("luad", "A_R20cap")}
    assert c1["analyses"] == {("luad", "C3"), ("luad", "M|C3")} and len(c1["tokens"]) == 5
    assert arm_closed(c1, "ucec", "A_C3") and not arm_closed(c1, "luad", "A_C3") and an_closed(c1, "luad", "M|C3")
    for oktxt in ("arm ucec A_R20\narm ucec A_R20cap\n", "analysis ucec M|C3\nanalysis ucec C3\n",
                  "arm ucec A_C3\narm ucec V_C3\nanalysis ucec C3\nanalysis ucec M|C3"):     # no final newline
        wr(oktxt)
        assert parse_closed(cp, cohorts)["tokens"]
    for badtxt in ("arm ucec A_R20cap_001\n", "arm brca A_C3\n", "analysis ucec M\n", "arm ucec\n",
                   "closed ucec A_C3\n", "arm ucec A_C3\narm ucec A_C3\n", "arm ucec A_C3 extra\n",
                   "arm ucec A_C3 \n",                      # a trailing space
                   " arm ucec A_C3\n",                      # a leading space
                   "arm  ucec A_C3\n",                      # a double space
                   "arm\tucec\tA_C3\n",                     # tabs
                   "arm ucec A_C3\r\n",                     # a carriage return
                   "# ok\n  # an indented comment\n",       # an indented comment is not a comment
                   "# ok\n \n",                             # a whitespace-only line is not an empty line
                   "arm ucec A_R20cap\n",                   # A_R20cap without A_R20
                   "arm ucec A_R20\n",                      # A_R20 without A_R20cap
                   "analysis ucec C3\n",                    # C3 without M|C3
                   "analysis ucec M|C3\n",                  # M|C3 without C3
                   "analysis ucec R20\n", "analysis ucec M|R20\n", "analysis ucec OP\n"):   # not closable
        wr(badtxt)
        _expect_fatal(parse_closed, cp, cohorts)
    with open(cp, "wb") as fh:
        fh.write(b"\xef\xbb\xbf# bom\n")
    _expect_fatal(parse_closed, cp, cohorts)
    with open(cp, "wb") as fh:
        fh.write(b"# \xff\xfe not utf-8\n")
    _expect_fatal(parse_closed, cp, cohorts)
    _expect_fatal(parse_closed, os.path.join(tmp, "absent.txt"), cohorts)
    assert needs_scores(c1, "ucec") and not needs_scores(c1, "luad")
    assert fully_closed(dict(arms={("luad", a) for a in ARM_TOKENS}, analyses=set(), tokens=[]), "luad")
    assert not fully_closed(c1, "luad") and not fully_closed(no_closed(), "ucec")
    # PINNED_HOSTS.txt
    pp = os.path.join(tmp, "pinned.txt")
    wr("# header\n# comments only\n\n", pp)
    assert parse_pinned(pp, cohorts) == {}
    wr("# h\npinned ucec gpu02\n\npinned luad s006\n", pp)
    assert parse_pinned(pp, cohorts) == {"ucec": "gpu02", "luad": "s006"}
    for badtxt in ("pinned ucec gpu02\npinned ucec s006\n", "pinned brca gpu02\n", "pinned ucec\n",
                   "pinned ucec gpu02 x\n", "pinned ucec gpu02 \n", " pinned ucec gpu02\n", "pinned  ucec gpu02\n",
                   "arm ucec A_C3\n", "pinned ucec gpu02\r\n", "# a\n  # b\n", "# a\n \n"):
        wr(badtxt, pp)
        _expect_fatal(parse_pinned, pp, cohorts)
    _expect_fatal(parse_pinned, os.path.join(tmp, "absent.txt"), cohorts)
    # at stage time every cohort that is not fully closed needs a pinned line
    pctx = Ctx(os.path.join(tmp, "pc"), cohorts, 5, True)
    os.makedirs(os.path.dirname(pctx.pinned_path))
    all_luad = dict(arms={("luad", a) for a in ARM_TOKENS}, analyses=set(), tokens=[])
    all_both = dict(arms={(k, a) for k in cohorts for a in ARM_TOKENS}, analyses=set(), tokens=[])
    wr("# nothing pinned\n", pctx.pinned_path)
    _expect_fatal(check_pinned, pctx, no_closed())
    _expect_fatal(check_pinned, pctx, all_luad)                          # ucec is not fully closed
    assert check_pinned(pctx, all_both) == {}
    wr("pinned ucec gpu02\n", pctx.pinned_path)
    assert check_pinned(pctx, all_luad) == {"ucec": "gpu02"}
    _expect_fatal(check_pinned, pctx, no_closed())                       # luad has no line
    wr("pinned ucec gpu02\npinned luad s006\n", pctx.pinned_path)
    assert check_pinned(pctx, no_closed()) == {"ucec": "gpu02", "luad": "s006"}
    # manifest: only '<sha256> *<name>' lines; status lines, bad names and duplicates are refused
    genes = np.array(["a", "b", "c"])
    hdir = os.path.join(tmp, "mroot")
    ctx = Ctx(hdir, {}, 5, True)
    os.makedirs(ctx.stage_h_dir)
    for nm in ("x1.npz", "x2.npz"):
        np.savez_compressed(ctx.H(nm), genes=genes, M=np.zeros((3, 6)), cases=np.array(["c1", "c2"]))
    h1 = sha256(ctx.H("x1.npz"))
    with open(ctx.H(MANIFEST_NAME), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("%s *x1.npz\n%s *x2.npz\n" % (h1, sha256(ctx.H("x2.npz"))))
    man = parse_manifest(ctx.H(MANIFEST_NAME))
    assert man == {"x1.npz": h1, "x2.npz": sha256(ctx.H("x2.npz"))}
    for badtxt in ("%s *x1.npz\nnot_completed x3.npz\n" % h1, "%s *x1.npz\nuntestable x3.npz\n" % h1,
                   "%s *x1.npz\n%s *x1.npz\n" % (h1, h1), "%s *sub/x1.npz\n" % h1, "%s x1.npz\n" % h1,
                   "# comment\n%s *x1.npz\n" % h1, "%s  *x1.npz\n" % h1):
        with open(ctx.H("bad.sha256"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(badtxt)
        _expect_fatal(parse_manifest, ctx.H("bad.sha256"))
    _expect_fatal(verify_manifest, ctx, "f" * 64)
    _expect_fatal(verify_manifest, ctx, None)
    assert verify_manifest(ctx, sha256(ctx.H(MANIFEST_NAME))) == man
    data, probs = read_h_npz(ctx.H("x1.npz"), genes, 6, 2)
    assert not probs and data["M"].shape == (3, 6)
    assert read_h_npz(ctx.H("x1.npz"), genes[::-1], 6, 2)[1]          # gene order
    assert read_h_npz(ctx.H("x1.npz"), genes, 7, 2)[1]                # columns
    np.savez_compressed(ctx.H("f32.npz"), genes=genes, M=np.zeros((3, 6), np.float32), cases=np.array(["c1", "c2"]))
    assert read_h_npz(ctx.H("f32.npz"), genes, 6, 2)[1]               # dtype
    np.savez_compressed(ctx.H("obj.npz"), genes=np.array(["a", None, "c"], dtype=object), M=np.zeros((3, 6)),
                        cases=np.array(["c1", "c2"]))
    assert read_h_npz(ctx.H("obj.npz"), genes, 6, 2)[0] is None       # object arrays are never unpickled
    np.savez_compressed(ctx.H("xk.npz"), genes=genes, M=np.zeros((3, 6)), cases=np.array(["c1", "c2"]), extra=1)
    assert read_h_npz(ctx.H("xk.npz"), genes, 6, 2)[1]                # extra key in a matrix file
    # covariates record
    cs = np.array(["c1", "c2", "c3"])
    mk = np.array(json.dumps({"stroma": {"used": ["A", "B"]}}))
    np.savez_compressed(ctx.H("rec.npz"), cases=cs, C3=np.zeros((3, 3)), R20=np.ones((3, 20)), pcs=np.ones((3, 4)),
                        markers_json=mk, c3_testable=np.bool_(True))
    data, probs = read_cov_record(ctx.H("rec.npz"), 3, 4, True, True)
    assert not probs and data["markers"]["stroma"]["used"] == ["A", "B"]
    assert read_cov_record(ctx.H("rec.npz"), 3, 5, True, True)[1]      # wrong K
    np.savez_compressed(ctx.H("rec0.npz"), cases=cs, C3=np.zeros((3, 0)), R20=np.ones((3, 20)), pcs=np.ones((3, 4)),
                        markers_json=mk)
    assert not read_cov_record(ctx.H("rec0.npz"), 3, 4, False, True)[1]   # n x 0 only if C3 is closed
    assert read_cov_record(ctx.H("rec0.npz"), 3, 4, True, True)[1]
    np.savez_compressed(ctx.H("recn.npz"), cases=cs, R20=np.ones((3, 20)), pcs=np.ones((3, 4)), markers_json=mk)
    assert not read_cov_record(ctx.H("recn.npz"), 3, 4, False, True)[1]   # C3 absent: allowed only if closed
    assert read_cov_record(ctx.H("recn.npz"), 3, 4, True, True)[1]
    np.savez_compressed(ctx.H("recb.npz"), cases=cs, C3=np.zeros((3, 3)), R20=np.ones((3, 20)), pcs=np.ones((3, 4)),
                        markers_json=np.array("not json {"))
    assert read_cov_record(ctx.H("recb.npz"), 3, 4, True, True)[1]
    np.savez_compressed(ctx.H("recp.npz"), cases=cs, C3=np.zeros((3, 3)), R20=np.ones((3, 20)), markers_json=mk)
    assert read_cov_record(ctx.H("recp.npz"), 3, 4, True, True)[1]       # pcs missing
    print("  score table, CLOSED_PARTS.txt and PINNED_HOSTS.txt (strict, shared rules), manifest, Stage H matrix "
          "and covariates-record parsers: OK")


def _selftest_classes():
    """Every combination of the PS 6.1 inputs, with closures, against a table written out by hand."""
    c3_inputs = {"H": ("ok", 0.01), "N(0.05)": ("ok", 0.05), "N": ("ok", 0.30), "closed": ("closed", None)}
    # (r20 status, p_R, placebo status, p_abs) -> R20 class
    r20_inputs = [
        (("ok", 0.01, "ok", 0.50), "H"),
        (("ok", 0.30, "ok", 1 / 120.), "A"),
        (("ok", 0.30, "ok", 5 / 120.), "A"),                     # e = 4: p_abs 0.0417 < 0.05
        (("ok", 0.30, "ok", 6 / 120.), "U"),                     # e = 5: p_abs = 0.05 is not < 0.05
        (("ok", 0.05, "ok", 0.50), "U"),                         # p_R = 0.05 is not < 0.05
        (("closed", None, "ok", None), "closed"),                # A_R20 closed
        (("ok", 0.01, "closed", None), "H"),                     # placebos closed but p_R < 0.05: H
        (("ok", 0.30, "closed", None), "closed"),                # placebos closed and p_R >= 0.05
        (("closed", None, "closed", None), "closed"),
    ]
    want_type = {("H", "H"): "none", ("H", "A"): "R", ("H", "U"): "none", ("H", "closed"): "NC:R",
                 ("N", "H"): "C", ("N", "A"): "C", ("N", "U"): "C", ("N", "closed"): "C",
                 ("closed", "A"): "R", ("closed", "H"): "NC:C", ("closed", "U"): "NC:C",
                 ("closed", "closed"): "NC:C"}
    n = 0
    seen = set()
    for c3name, (c3s, pc3) in c3_inputs.items():
        c3cls = c3name.split("(")[0]
        for (rs, pr, cs, pa), r20cls in r20_inputs:
            r = classify_cohort(c3s, pc3, rs, pr, cs, pa)
            typ = want_type[(c3cls, r20cls)]
            got = (r["C3_class"], r["R20_class"], r["type"], r["flagged"], r["U"])
            want = (c3cls, r20cls, typ, typ != "none", r20cls == "U")
            assert got == want, (c3name, (rs, pr, cs, pa), got, want)
            seen.add(typ)
            n += 1
    assert seen == {"C", "R", "NC:C", "NC:R", "none"}
    _expect_fatal(classify_cohort, "untestable", None, "ok", 0.01, "ok", 0.5)
    _expect_fatal(classify_cohort, "ok", 0.01, "not_completed", None, "ok", 0.5)
    assert case_of([False, False]) == 0 and case_of([True, False]) == 1 and case_of([True, True]) == 2
    m = m4_count([("a", "ok", 0.01), ("b", "ok", 0.05), ("c", "untestable", None), ("d", "ok", 0.049),
                  ("e", "closed", None)])
    assert (m["k"], m["t"], m["n_cells"], m["outranked"]) == (2, 3, 5, ["a", "d"])
    assert m["not_tested"] == [["c", "untestable"], ["e", "closed"]]
    print("  PS 6.1 classes and types (%d combinations, closures included), case and M4 counts: OK" % n)


def _selftest_d_cap():
    # D on constructed AUC vectors
    a_pub = np.array([0.70, 0.50, 0.52, 0.48, 0.51, 0.49])
    a_x = np.array([0.60, 0.50, 0.50, 0.50, 0.50, 0.50])
    r = d_change(a_pub, a_x)
    Dn = np.array([0.0, 0.02, -0.02, 0.01, -0.01])
    assert abs(r["D_0"] - 0.10) < 1e-15 and abs(r["D_null_mean"] - Dn.mean()) < 1e-15
    assert abs(r["D_null_sd"] - Dn.std()) < 1e-15 and r["p_D"] == 1 / 6.
    r2 = d_change(a_x, a_pub)
    assert r2["p_D"] == 1.0                                          # D_0 below every D_b
    # e and p_abs over 119 placebos: ties count; class A needs e <= 4
    mj = [0.1] * 115 + [-0.2, -0.2, -0.3, -0.3]
    c = capacity_p(mj, -0.2)
    assert c["e"] == 4 and abs(c["p_abs"] - 5 / 120.) < 1e-15 and c["p_abs"] < ALPHA
    c5 = capacity_p(mj[:-1] + [-0.2], -0.2 + 1e-12)
    assert c5["e"] == 4
    c6 = capacity_p([0.1] * 114 + [-0.2] * 5, -0.2)
    assert c6["e"] == 5 and c6["p_abs"] == 6 / 120. and not c6["p_abs"] < ALPHA
    assert capacity_p([0.1] * 119, -0.5) == dict(e=0, p_abs=1 / 120.)
    _expect_fatal(capacity_p, [0.1] * 118, 0.0)
    _expect_fatal(capacity_p, [0.1] * 59, 0.0)
    # placebo seeds and draws of PS 2 (the cluster draws them; recorded here for the record)
    import zlib
    q = np.random.RandomState(zlib.crc32(b"c1_cov|%s|cap|%03d" % (b"ucec", 7))).permutation(138)
    assert (b"c1_cov|%s|cap|%03d" % (b"ucec", 7)) == b"c1_cov|ucec|cap|007" and sorted(q) == list(range(138))
    assert cap_name("luad", 7) == "c1_cov_luad_A_R20cap_007.npz" and cap_name("ucec", 119).endswith("_119.npz")
    # m_j and m_R on small constructed matrices: set genes rank first in column 0 only
    uni = np.array(["g%03d" % i for i in range(100)])
    rng = np.random.RandomState(5)
    Mc = rng.rand(100, 1 + CAP_COLS)
    Mc[:50, 0] = 2.0                                                # AUC 1 in column 0
    rc = cell(Mc, uni, set(uni[:50]))
    assert rc["auc"] == 1.0 and rc["n_null"] == CAP_COLS
    assert abs(rc["margin"] - (1.0 - assemble(Mc, np.arange(100.0), np.arange(100) < 50)[0][1:].mean())) < 1e-15
    # m_R from the first 25 permutation columns of a 1,001-column arm
    Mf = rng.rand(100, 1001)
    Mf[:50, 0] = 2.0
    mR = cell(Mf[:, :CAP_COLS + 1], uni, set(uni[:50]))["margin"]
    a = assemble(Mf, np.arange(100.0), np.arange(100) < 50)[0]
    assert abs(mR - (a[0] - a[1:26].mean())) < 1e-15
    print("  change D, e and p_abs over 119 placebos, m_R from 25 columns: OK")


def _selftest_run_stage(tmp):
    ctx = Ctx(os.path.join(tmp, "rs"), {}, 5, True)
    os.makedirs(os.path.dirname(ctx.sidecar))
    os.makedirs(os.path.dirname(ctx.closed_path))
    for f in (ctx.sidecar, ctx.prespec, ctx.closed_path):            # run_stage records these hashes
        with open(f, "w") as fh:
            fh.write("# test\n")
    for src, dst in ((os.path.abspath(__file__), ctx.script_path), (os.path.abspath(tk.__file__), ctx.kernel_path)):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src, dst)

    def boom(od):
        with open(os.path.join(od, "x.txt"), "w") as fh:
            fh.write("x")
        raise RuntimeError("simulated crash")

    def fine(od):
        with open(os.path.join(od, "x.txt"), "w") as fh:
            fh.write("y")
        return dict(n=1)

    for k in range(4):                                       # four interrupted attempts in a row
        try:
            run_stage(ctx, "LX", boom, 1)
        except RuntimeError:
            pass
    assert os.path.isdir(os.path.join(ctx.out_dir, "LX.partial"))
    assert all(os.path.isdir(os.path.join(ctx.out_dir, "LX_failed%d" % k)) for k in (1, 2, 3))
    run_stage(ctx, "LX", fine, 1)                             # no limit on attempts: the fifth completes
    assert os.path.isdir(os.path.join(ctx.out_dir, "LX_failed4")) and not os.path.exists(
        os.path.join(ctx.out_dir, "LX.partial"))
    with open(os.path.join(ctx.out_dir, "LX", "result_LX.json"), encoding="utf-8") as fh:
        r = json.load(fh)
    assert r["attempts_before_this"] == ["LX_failed1", "LX_failed2", "LX_failed3", "LX_failed4"]
    assert r["outputs"] == {"x.txt": sha256(os.path.join(ctx.out_dir, "LX", "x.txt"))} and "numpy" in r["env"]
    with open(os.path.join(ctx.out_dir, "LX", "LX.log"), encoding="utf-8") as fh:
        assert "complete" in fh.read()
    before = _tree(ctx.out_dir)
    _expect_fatal(run_stage, ctx, "LX", fine, 1)             # completed: never rerun
    assert _tree(ctx.out_dir) == before
    print("  run safety (.partial, _failed1.._failed4, unlimited reruns, completed never rerun, tee): OK")


def _selftest_sidecar_list():
    """SIDECAR_FILES.txt (beside this script) names every file the real stages require."""
    p = os.path.join(PAPER, *SIDECAR_LIST_REL.split("/"))
    if not os.path.exists(p):
        raise AssertionError("SIDECAR_FILES.txt not found beside the script")
    with open(p, encoding="utf-8") as fh:
        rels = [ln.strip() for ln in fh if ln.strip() and not ln.strip().startswith("#")]
    missing = [r for r in real_ctx().required_rels(None) if r not in rels]
    assert not missing, missing
    assert len(rels) == len(set(rels))
    print("  SIDECAR_FILES.txt names all %d files the real stages require: OK" % len(real_ctx().required_rels(None)))


def _selftest_make_sidecar():
    """make_sidecar.py compares the full required list (for the parsed CLOSED_PARTS.txt) with
    SIDECAR_FILES.txt and refuses if a required path is missing from the list."""
    sys.dont_write_bytecode = True
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import make_sidecar as ms
    this = sys.modules[__name__]
    with open(os.path.join(PAPER, *SIDECAR_LIST_REL.split("/")), encoding="utf-8") as fh:
        rels = [ln.strip() for ln in fh if ln.strip() and not ln.strip().startswith("#")]
    assert ms.required_problems(rels, this) == []
    cut = [r for r in rels if r != PINNED_REL and not r.endswith("luad__sitepack_operator.csv")]
    miss = ms.required_problems(cut, this)
    assert len(miss) == 2 and any(PINNED_REL in m for m in miss) and any("luad__sitepack" in m for m in miss), miss
    print("  make_sidecar.py: SIDECAR_FILES.txt covers the required list; a missing PINNED_HOSTS.txt or sitepack "
          "table is refused: OK")


DUMMY_RELS = [CLUSTER_SCRIPT_REL] + list(WRAPPER_RELS) + [MAKE_ABSTRACT_REL, ABSTRACT_FINAL_REL, ABSTRACT_VARIANTS_REL] + \
    [r for r in DEFINITION_MODULES if r != KERNEL_REL]


def write_sidecar(ctx, rels):
    with open(ctx.sidecar, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# written synthetic\n")
        for rel in rels:
            fh.write("%s *%s\n" % (sha256(ctx.P(rel)), rel))


def build_synthetic(root, closed_tokens, shifts, cap_shift, B=99, seed=7):
    """Tiny fake cohorts 'ucec' and 'luad' under `root` with the real relative layout: every file a
    stage reads, CLOSED_PARTS.txt, a frozen sidecar and a Stage H manifest (closed arms are not
    written to the manifest).  shifts: {cohort: {arm: shift of the set genes in column 0}};
    cap_shift: {cohort: shift in the placebo files}.  Returns (ctx, manifest sha256)."""
    rng = np.random.RandomState(seed)
    os.makedirs(root)
    for rel in [LOCAL_SCRIPT_REL, KERNEL_REL, PRESPEC_REL, CLOSED_REL, PINNED_REL] + DUMMY_RELS:
        os.makedirs(os.path.dirname(os.path.join(root, *rel.split("/"))), exist_ok=True)
    shutil.copyfile(os.path.abspath(__file__), os.path.join(root, *LOCAL_SCRIPT_REL.split("/")))
    shutil.copyfile(os.path.abspath(tk.__file__), os.path.join(root, *KERNEL_REL.split("/")))
    for rel in [PRESPEC_REL] + DUMMY_RELS:
        with open(os.path.join(root, *rel.split("/")), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("synthetic %s\n" % rel)
    with open(os.path.join(root, *CLOSED_REL.split("/")), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# synthetic CLOSED_PARTS.txt\n" + "".join(t + "\n" for t in closed_tokens))
    spec = {"ucec": dict(nconf=36, K=6, arms=[("unrestricted", "M_unrestricted"), ("within-operator", "M_within_operator"),
                                               ("residualised", "M_within_operator_residualized")], kind="ucec",
                         label="uterus"),
            "luad": dict(nconf=34, K=5, arms=[("unrestricted", "M")], kind="luad", label="lung")}
    n, G, nsel = 40, 160, 70
    cohorts = {}
    tmpctx = Ctx(root, {}, B, True)
    os.makedirs(tmpctx.stage_h_dir)
    for key, s in spec.items():
        cohorts[key] = {}
    closed = parse_closed(os.path.join(root, *CLOSED_REL.split("/")), cohorts)
    with open(os.path.join(root, *PINNED_REL.split("/")), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# synthetic PINNED_HOSTS.txt\n" + "".join("pinned %s host_%s\n" % (k, k) for k in spec
                                                           if not fully_closed(closed, k)))
    for key, s in spec.items():
        pats = np.array(["P%03d" % i for i in range(n)])
        genes = np.array(["G%03d" % i for i in range(G)])
        lab = sum((["L%02d" % i] * 3 for i in range(12)), []) + ["L12"] * 2 + [""] + ["L13"]
        labels = np.array(lab)[rng.permutation(n)]
        mrna = rng.randn(n, G)
        C3 = rng.randn(n, 3)
        pcs = rng.randn(n, 20).astype(np.float32)
        r20 = rng.randn(n, 20)
        prot = 0.5 * mrna + 0.5 * rng.randn(n, G)
        prot[:, :80] += 1.5 * C3[:, [0]]                             # strong C3 genes: S_C3 is not empty
        for j in range(G):
            prot[rng.choice(n, 3, replace=False), j] = np.nan
        pub_exc = rng.randint(5, 60, G)
        pub_inc = rng.randn(G) * 0.1
        selidx = np.arange(40, 40 + nsel)
        pub_exc[selidx] = rng.randint(0, 2, nsel)
        pub_inc[selidx] = np.abs(pub_inc[selidx]) + 1e-3
        pub_pval = (1 + pub_exc) / float(B + 1)
        pub_fdr = bh_fdr(pub_pval)
        pub_fdr[selidx] = np.minimum(pub_fdr[selidx], 0.01)
        pub_fdr[np.setdiff1d(np.arange(G), selidx)] = np.maximum(pub_fdr[np.setdiff1d(np.arange(G), selidx)], 0.06)
        inp = os.path.join(root, "inputs_%s.npz" % key)
        np.savez_compressed(inp, is_synthetic=np.array(True), patients=pats, genes=genes, mrna_log2tpm1=mrna,
                            protein=prot, pcs=pcs, rna_pcs20=r20, label_operator=labels, ra_n_perm=np.array(B),
                            ra_seed=np.array(0), cv_folds=np.array(5), ridge_alpha=np.array(1.0),
                            pub_n=(~np.isnan(prot)).sum(0).astype(float), pub_r2_rna=rng.randn(G),
                            pub_incremental_r2=pub_inc, pub_pval=pub_pval, pub_fdr=pub_fdr)
        drop = set(rng.choice(G, 10, replace=False).tolist())
        uni = np.array([g for i, g in enumerate(genes) if i not in drop])
        sm = set(genes[selidx])
        smask = np.array([g in sm for g in uni])
        Gu = len(uni)
        # the operator-axis (sitepack) table: S_BC = genes 20..79 (60 genes); decoys 100..104 (fdr < 0.05 but a
        # negative increment) and 105..109 (a positive increment but fdr 0.2) are not in it; in LUAD one
        # universe gene has no finite increment, so LUAD has no rho
        rs = np.random.RandomState(1000 + len(key))
        sp_rows = ["cohort,batch_axis,gene,incremental_r2_batchresid,fdr_batchresid"]
        for i, g in enumerate(genes):
            if 20 <= i < 80:
                f_, v_ = 0.01, abs(rs.randn()) * 0.05 + 1e-3
            elif 100 <= i < 105:
                f_, v_ = 0.01, -0.05
            elif 105 <= i < 110:
                f_, v_ = 0.2, 0.05
            else:
                f_, v_ = 0.5, rs.randn() * 0.05
            if key == "luad" and g == uni[-1]:
                v_ = float("nan")
            sp_rows.append("%s,operator,%s,%s,%s" % (key, g, repr(v_), repr(f_)))
        spp = os.path.join(root, "sitepack_%s.csv" % key)
        with open(spp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(sp_rows) + "\n")
        sp_sha = sha256(spp)

        def mat(shift, cols=B + 1):
            M = rng.randn(Gu, cols)
            M[smask, 0] += shift
            return M
        arms = {k: mat(1.5) for _, k in s["arms"]}
        conf = os.path.join(root, "conf_%s.npz" % key)
        np.savez_compressed(conf, genes=uni, sel_mask=smask, **arms)
        xs = dict(zip(genes, pub_inc))
        if s["kind"] == "ucec":
            pub = dict(arms={})
            for nm, k in s["arms"]:
                r = cell(arms[k], uni, sm, xs_by_gene=xs)
                pub["arms"][nm] = dict(auc=dict(obs=r["auc"], null_mean=r["null_mean"], null_sd=r["null_sd"],
                                                p_one_sided=r["p"]), rho=dict(obs=r["rho"], p_one_sided=r["p_rho"]))
        else:
            r = cell(arms["M"], uni, sm, xs_by_gene=xs)
            pub = dict(h1=dict(auc_obs=r["auc"], auc_null_mean=r["null_mean"], auc_null_sd=r["null_sd"],
                               p_auc=r["p"], rho_obs=r["rho"], p_rho=r["p_rho"]))
        pub["sha256_discovery_sitepack"] = sp_sha
        pubp = os.path.join(root, "pub_%s.json" % key)
        with open(pubp, "w", encoding="utf-8") as fh:
            json.dump(pub, fh)
        # composition controls and transcriptome-baseline tables for the published-set genes
        cc = ["gene,n,incr_orig,stroma_gain,immune_gain,proliferation_gain,ALL_gain"]
        tb = ["gene,n,r2_rna_pcs20,incr_over_pcs20"]
        for j in selidx:
            v = ~np.isnan(prot[:, j])
            y = prot[v, j]
            folds = tk.folds_for(int(v.sum()))
            xr = mrna[v, j].reshape(-1, 1)
            r2r = tk.r2_of(y, tk.oof_ridge(xr, y, folds))
            gains = [tk.r2_of(y, tk.oof_ridge(np.hstack([xr, C3[v][:, w]]), y, folds)) - r2r
                     for w in ([0], [1], [2], [0, 1, 2])]
            cc.append(",".join([genes[j], str(int(v.sum())), "0.0"] + [repr(float(x)) for x in gains]))
            base = np.hstack([xr, r20[v]])
            r2b = tk.r2_of(y, tk.oof_ridge(base, y, folds))
            r2w = tk.r2_of(y, tk.oof_ridge(np.hstack([base, pcs.astype(np.float64)[v]]), y, folds))
            tb.append(",".join([genes[j], str(int(v.sum())), repr(float(r2b)), repr(float(r2w - r2b))]))
        ccp, tbp = os.path.join(root, "compctl_%s.csv" % key), os.path.join(root, "txperm_%s.csv" % key)
        for p_, lines in ((ccp, cc), (tbp, tb)):
            with open(p_, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("\n".join(lines) + "\n")
        if needs_scores(closed, key):
            body = "# synthetic H0 export\r\npatient\tstroma\timmune\tproliferation\r\n" + "".join(
                "%s\t%s\r\n" % (p, "\t".join(repr(float(x)) for x in C3[i])) for i, p in enumerate(pats))
            with open(tmpctx.H(scores_name(key)), "wb") as fh:
                fh.write(body.encode("utf-8"))
        cases = np.array(["C%03d" % i for i in range(s["nconf"])])
        for arm in MAIN_ARMS:                     # a closed arm's file may exist; it is never read
            np.savez_compressed(tmpctx.H(stage_h_name(key, arm)), genes=uni, M=mat(shifts[key][arm]), cases=cases)
        if not arm_closed(closed, key, "A_R20cap"):
            for j in range(1, N_CAP + 1):
                np.savez_compressed(tmpctx.H(cap_name(key, j)), genes=uni, M=mat(cap_shift[key], CAP_COLS + 1),
                                    cases=cases)
        nconf, K = s["nconf"], s["K"]
        pcs_c = rng.randn(nconf, K)
        c3_closed = arm_closed(closed, key, "A_C3") and arm_closed(closed, key, "V_C3")
        C3c = np.empty((nconf, 0)) if c3_closed else pcs_c[:, :3] * 0.5 + rng.randn(nconf, 3)
        np.savez_compressed(tmpctx.H(covariates_record_name(key)), cases=cases, C3=C3c, R20=rng.randn(nconf, 20),
                            pcs=pcs_c, c3_testable=np.bool_(not c3_closed),
                            markers_json=np.array(json.dumps({"stroma": {"used": ["COL1A1"]}}, sort_keys=True)))
        cohorts[key] = dict(label=s["label"], inputs=os.path.basename(inp), inputs_sha=sha256(inp),
                            conf=os.path.basename(conf), conf_sha=sha256(conf), arms=s["arms"],
                            published=os.path.basename(pubp), published_kind=s["kind"],
                            compctl=os.path.basename(ccp), txperm=os.path.basename(tbp),
                            n_disc=n, n_fit=G, n_sel=nsel, n_univ=Gu, n1=int(smask.sum()), n0=int((~smask).sum()),
                            n_op_dummies=12, n_op_covered=36, n_conf_cases=nconf, K=K, tol_v1b_inc=1e-9,
                            sitepack=os.path.basename(spp), n_bc=60, n_bc_universe=None,
                            sitepack_prefix=sp_sha[:8] if key == "luad" else None,
                            sitepack_pub_key="sha256_discovery_sitepack" if key == "luad" else None,
                            sitepack_sha=None)
    pins = {KERNEL_REL: sha256(os.path.join(root, *KERNEL_REL.split("/"))),
            DUMMY_RELS[-1]: sha256(os.path.join(root, *DUMMY_RELS[-1].split("/")))}
    ctx = Ctx(root, cohorts, B, True, pins)
    man = []
    for key in cohorts:
        for name, _, _ in stage_h_expected(ctx, key, closed):
            man.append("%s *%s" % (sha256(ctx.H(name)), name))
        man.append("%s *%s" % (sha256(ctx.H(covariates_record_name(key))), covariates_record_name(key)))
    with open(ctx.H(MANIFEST_NAME), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(man) + "\n")
    write_sidecar(ctx, ctx.required_rels(closed))
    return ctx, sha256(ctx.H(MANIFEST_NAME))


def _selftest_sidecar_check(ctx):
    """--check-sidecar / the stage gate on a temporary tree."""
    assert check_sidecar_mode(ctx) and not sidecar_problems(ctx)[0]
    closed = sidecar_problems(ctx)[2]
    rels = ctx.required_rels(closed)
    with open(ctx.sidecar, "rb") as fh:
        good = fh.read()

    def expect_problem(substr):
        probs = sidecar_problems(ctx)[0]
        assert any(substr in p for p in probs), (substr, probs)
        assert not check_sidecar_mode(ctx)
        _expect_fatal(check_sidecar, ctx)
    # a required file not listed
    write_sidecar(ctx, [r for r in rels if r != WRAPPER_RELS[2]])
    expect_problem("lsf_c1_cov_discovery.sh is not listed")
    write_sidecar(ctx, [r for r in rels if r != CLOSED_REL])
    expect_problem("CLOSED_PARTS.txt is not listed")
    write_sidecar(ctx, [r for r in rels if r != "server_export/scripts/transcriptome_baseline.py"])
    expect_problem("transcriptome_baseline.py is not listed")
    # a listed file changed after freezing
    with open(ctx.sidecar, "wb") as fh:
        fh.write(good)
    tgt = ctx.P(ABSTRACT_VARIANTS_REL)
    with open(tgt, "rb") as fh:
        orig = fh.read()
    with open(tgt, "ab") as fh:
        fh.write(b"x")
    expect_problem("changed since freezing")
    with open(tgt, "wb") as fh:
        fh.write(orig)
    # a pinned file listed with another hash (file and sidecar agree, the pin does not)
    pinned = DUMMY_RELS[-1]
    with open(ctx.P(pinned), "rb") as fh:
        porig = fh.read()
    with open(ctx.P(pinned), "ab") as fh:
        fh.write(b"y")
    write_sidecar(ctx, rels)
    expect_problem("not the pinned")
    with open(ctx.P(pinned), "wb") as fh:
        fh.write(porig)
    # a malformed CLOSED_PARTS.txt (listed with its hash) is a problem
    with open(ctx.closed_path, "rb") as fh:
        corig = fh.read()
    with open(ctx.closed_path, "ab") as fh:
        fh.write(b"arm ucec A_X\n")
    write_sidecar(ctx, rels)
    expect_problem("unknown arm")
    with open(ctx.closed_path, "wb") as fh:
        fh.write(corig)
    # a malformed PINNED_HOSTS.txt (listed with its hash) is a problem
    with open(ctx.pinned_path, "rb") as fh:
        porig2 = fh.read()
    with open(ctx.pinned_path, "ab") as fh:
        fh.write(b"pinned ucec another_host\n")
    write_sidecar(ctx, rels)
    expect_problem("two pinned lines")
    with open(ctx.pinned_path, "wb") as fh:
        fh.write(porig2)
    # a malformed line and a missing sidecar
    with open(ctx.sidecar, "wb") as fh:
        fh.write(good + b"nothex *x\n")
    expect_problem("malformed sidecar line")
    os.remove(ctx.sidecar)
    expect_problem("not found")
    with open(ctx.sidecar, "wb") as fh:
        fh.write(good)
    assert not sidecar_problems(ctx)[0]
    print("  --check-sidecar on a temporary tree (unlisted, changed, pinned, CLOSED_PARTS, PINNED_HOSTS, malformed, "
          "missing): OK")
    # --validate-only and --gate-v3 refuse when the sidecar exists (through main(), on this temporary tree) and
    # write nothing; without a sidecar the refusal does not apply
    _expect_fatal(refuse_if_frozen, ctx, "--validate-only")
    global real_ctx
    keep_rc, keep_argv = real_ctx, sys.argv
    real_ctx = lambda: ctx
    try:
        for mode in ("--validate-only", "--gate-v3"):
            before = _tree(ctx.root)
            sys.argv = ["c1_covariates_local.py", mode]
            try:
                main()
            except SystemExit as e:
                assert e.code == 2, (mode, e.code)
            else:
                raise AssertionError("%s was not refused although the sidecar exists" % mode)
            assert _tree(ctx.root) == before, "a refused %s wrote something" % mode
    finally:
        real_ctx, sys.argv = keep_rc, keep_argv
    os.rename(ctx.sidecar, ctx.sidecar + ".away")
    refuse_if_frozen(ctx, "--validate-only")
    refuse_if_frozen(ctx, "--gate-v3")
    os.rename(ctx.sidecar + ".away", ctx.sidecar)
    print("  --validate-only and --gate-v3 refuse when the sidecar exists and write nothing: OK")


def _selftest_end_to_end_1(tmp):
    """Scenario 1, with closures: ucec C3 H, R20 closed (A_R20, the placebos and V_R20 together: type NC:R);
    luad C3 closed (A_C3 and V_C3 with the analyses C3 and M|C3, as the marker rule would), R20 H
    (type NC:C).  Pre-flight refusals create nothing; leftover partials are kept and the stage runs again."""
    root = os.path.join(tmp, "synthetic_root_1")
    tokens = ["arm ucec A_R20", "arm ucec A_R20cap", "arm ucec V_R20", "arm luad A_C3", "arm luad V_C3",
              "analysis luad C3", "analysis luad M|C3"]
    shifts = {"ucec": dict(A_C3=1.5, A_R20=0.0, V_C3=1.5, V_R20=0.0),
              "luad": dict(A_C3=0.0, A_R20=1.5, V_C3=0.0, V_R20=1.5)}
    ctx, msha = build_synthetic(root, tokens, shifts, cap_shift={"ucec": 1.5, "luad": 1.5})
    real_out = os.path.normcase(os.path.abspath(os.path.join(PAPER, *OUT_REL.split("/"))))
    assert not os.path.normcase(ctx.out_dir).startswith(real_out) and ctx.synthetic
    assert not os.path.normcase(ctx.root).startswith(os.path.normcase(PAPER) + os.sep)
    _selftest_sidecar_check(ctx)
    # PINNED_HOSTS.txt: every cohort that is not fully closed needs a pinned line (checked at stage time)
    closed0 = parse_closed(ctx.closed_path, ctx.cohorts)
    assert check_pinned(ctx, closed0) == {"ucec": "host_ucec", "luad": "host_luad"}
    with open(ctx.pinned_path, "rb") as fh:
        pk = fh.read()
    with open(ctx.pinned_path, "wb") as fh:
        fh.write(pk.replace(b"pinned luad host_luad\n", b""))
    _expect_fatal(check_pinned, ctx, closed0)
    with open(ctx.pinned_path, "wb") as fh:
        fh.write(pk)
    # pre-freeze gates on the synthetic files
    d_u0 = load_inputs(ctx, "ucec")          # luad has no score table here: its analyses C3 and M|C3 are closed
    assert validate_v3(ctx, "ucec", d_u0, load_scores(ctx.H(scores_name("ucec")), d_u0["patients"])[0])["ok"]
    for key in ctx.cohorts:
        d = load_inputs(ctx, key)
        assert validate_v1b(ctx, key, d)["ok"]
        assert validate_v2(ctx, key, d, load_conf(ctx, key))["ok"]
    assert check_stage_h(ctx, ctx.B + 1, CAP_COLS + 1)
    assert not check_stage_h(ctx, 6, CAP_COLS + 1)                    # wrong --expect-cols is caught

    def refused_creates_nothing(*a):
        before = _tree(ctx.out_dir)
        _expect_fatal(do_stage, ctx, *a)
        assert _tree(ctx.out_dir) == before, "a refusal created something"
    # --- L1 pre-flight refusals ---
    refused_creates_nothing("L1", 1, None)                           # no --manifest-sha
    refused_creates_nothing("L1", 1, "0" * 64)                        # wrong manifest sha
    refused_creates_nothing("L2", 1, msha)                            # L2 before L1
    with open(ctx.sidecar, "rb") as fh:
        sc = fh.read()
    os.rename(ctx.sidecar, ctx.sidecar + ".away")
    refused_creates_nothing("L1", 1, msha)                            # no sidecar
    os.rename(ctx.sidecar + ".away", ctx.sidecar)
    sp = ctx.H(scores_name("ucec"))
    with open(sp, "rb") as fh:
        sorig = fh.read()
    txt = sorig.decode("utf-8").split("\r\n")                         # a changed score table ...
    f0 = txt[2].split("\t")
    f0[1] = repr(float(f0[1]) + 0.5)
    txt[2] = "\t".join(f0)
    with open(sp, "wb") as fh:
        fh.write("\r\n".join(txt).encode("utf-8"))
    refused_creates_nothing("L1", 1, msha)                            # ... is refused by the sidecar
    lines = sc.decode("utf-8").split("\n")                            # ... and by V3 if re-frozen
    with open(ctx.sidecar, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(("%s *%s" % (sha256(sp), STAGE_H_REL + "/" + scores_name("ucec")))
                           if ln.endswith(scores_name("ucec")) else ln for ln in lines))
    assert not sidecar_problems(ctx)[0]
    d0 = load_inputs(ctx, "ucec")
    Cbad = load_scores(sp, d0["patients"])[0]
    assert not validate_v3(ctx, "ucec", d0, Cbad)["ok"]
    refused_creates_nothing("L1", 1, msha)                            # V3 fails in the pre-flight
    with open(sp, "wb") as fh:
        fh.write(sorig)
    with open(ctx.sidecar, "wb") as fh:
        fh.write(sc)
    # --- L1 runs (after a leftover partial from an interrupted attempt), and is never rerun ---
    os.makedirs(os.path.join(ctx.out_dir, "L1.partial"))
    do_stage(ctx, "L1", 1, msha)
    assert os.path.isdir(os.path.join(ctx.out_dir, "L1_failed1"))
    refused_creates_nothing("L1", 1, msha)
    l1 = os.path.join(ctx.out_dir, "L1")
    with open(os.path.join(l1, "result_L1.json"), encoding="utf-8") as fh:
        r1 = json.load(fh)
    assert sorted(r1["outputs"]) == ["cells_L1.csv", "sets_luad.csv", "sets_ucec.csv"]
    assert r1["manifest_sha256"] == msha and sorted(r1["closed_parts"]) == sorted(tokens)
    su = r1["cohorts"]["ucec"]["sets"]
    for nm in ("T_M", "T_C3", "T_R20", "T_M|C3", "T_M|R20"):
        for k in ("overlap_S_M", "overlap_T_M", "largest_p", "bh_q_at_rank", "n_floor_in_set",
                  "floor_alone_filled_set", "n_nonpositive_increment"):
            assert k in su[nm], (nm, k)
    assert su["T_OP"]["status"] == "ok" and su["S_OP"]["status"] == "ok" and su["T_C3"]["status"] == "ok"
    sl = r1["cohorts"]["luad"]["sets"]
    for nm in ("T_C3", "S_C3", "T_M|C3", "S_M|C3"):
        assert sl[nm] == {"status": "closed"}, nm
    assert sl["T_OP"]["status"] == "ok" and sl["T_M|R20"]["status"] == "ok"
    assert su["T_M"]["overlap_S_M"] <= 70 and su["T_M"]["overlap_T_M"] == 70
    for nm in ("S_M", "S_C3", "S_R20", "S_M|C3", "S_M|R20"):
        assert all(k in su[nm] for k in ("n", "n_bh_before_sign_filter", "n_floor", "kstar")), nm
    bu_, bl_ = r1["cohorts"]["ucec"]["S_BC"], r1["cohorts"]["luad"]["S_BC"]
    assert bu_["n"] == 60 and bu_["rho_available"] is True and bu_["rho_note"] == ""
    assert bl_["n"] == 60 and bl_["rho_available"] is False and "no rho" in bl_["rho_note"]
    assert r1["cohorts"]["ucec"]["operator_indicators"]["n_indicators"] == 12
    assert r1["cohorts"]["ucec"]["operator_indicators"]["n_covered"] == 36
    c1 = read_cells(os.path.join(l1, "cells_L1.csv"))
    # per arm: reference 1, rule reference 1, forward 3, gene property 1, partition 9; on the unrestricted arm
    # also FDR-defined 3 and S_BC 1; ucec has three saved arms, luad one
    assert len(c1) == 15 * 4 + 4 * 2
    fw = [r for r in c1 if r["family"] == "forward"]
    assert len(fw) == 6 and all(r["read"] == "True" and r["null"] == "unrestricted" for r in fw)
    assert [r["status"] for r in fw if r["cohort"] == "luad" and r["set"] == "T_C3"] == ["closed"]
    assert all(r["status"] == "ok" for r in fw if not (r["cohort"] == "luad" and r["set"] == "T_C3"))
    assert sum(1 for r in c1 if r["read"] == "True") == 6
    assert sum(1 for r in c1 if r["family"] == "forward (descriptive)") == 6     # ucec's two other arms
    assert sum(1 for r in c1 if r["family"] == "rule reference") == 4             # 3 ucec arms + 1 luad
    assert all(r["rho"] != "" for r in c1 if r["family"] == "rule reference")
    assert sum(1 for r in c1 if r["family"] == "FDR-defined") == 6
    assert sum(1 for r in c1 if r["family"] == "partition") == 3 * 3 * 4         # 3 X x 3 parts x 4 arms
    assert sum(1 for r in c1 if r["family"] == "partition" and r["status"] == "closed") == 3   # luad C3
    assert all(r["rho"] == "" for r in c1 if r["family"] in ("gene property", "partition"))
    bc = [r for r in c1 if r["family"] == "batch-corrected"]
    assert len(bc) == 2 and all(r["set"] == "S_BC" and r["matrix"] == "published" and r["null"] == "unrestricted"
                                and r["read"] == "False" and r["status"] == "ok" for r in bc)
    ub = [r for r in bc if r["cohort"] == "ucec"][0]
    lb = [r for r in bc if r["cohort"] == "luad"][0]
    assert ub["rho"] != "" and lb["rho"] == "" and "no rho" in lb["note"] and int(ub["n1"]) >= 50
    S = read_sets(os.path.join(l1, "sets_ucec.csv"))
    assert all(len(S[X]["T"]) == 70 for X in ANALYSES) and S["OP"] is not None
    assert len(S["T_G"]) == 70 and len(S["S_M"]) == 70 and len(S["M"]["T"]) == 70
    assert len(S["C3"]["S"]) > 0                                                  # the strong C3 genes
    Sl_ = read_sets(os.path.join(l1, "sets_luad.csv"))
    assert Sl_["C3"] is None and Sl_["M|C3"] is None and Sl_["OP"] is not None and Sl_["M|R20"] is not None
    # an independent re-computation of T_M on one arm
    conf = load_conf(ctx, "ucec")
    tmr = [r for r in c1 if r["family"] == "rule reference" and r["cohort"] == "ucec" and r["null"] == "unrestricted"][0]
    mask = np.array([g in S["M"]["T"] for g in conf["genes"]])
    xsv = np.array([S["pub_inc"][g] for g in conf["genes"]])
    a, rho_, _, _ = assemble(conf["M"]["unrestricted"], xsv, mask)
    assert float(tmr["auc"]) == a[0] and float(tmr["rho"]) == rho_[0]
    # an independent re-computation of the S_BC cell (ucec: with rho) from the sitepack table
    with open(os.path.join(ctx.root, "sitepack_ucec.csv"), encoding="utf-8", newline="") as fh:
        sp = list(csv.DictReader(fh))
    sbc_set = set(r["gene"] for r in sp if float(r["fdr_batchresid"]) < 0.05 and float(r["incremental_r2_batchresid"]) > 0)
    sbc_inc = {r["gene"]: float(r["incremental_r2_batchresid"]) for r in sp}
    assert len(sbc_set) == 60 and "G100" not in sbc_set and "G105" not in sbc_set
    a_bc, rho_bc, _, _ = assemble(conf["M"]["unrestricted"], np.array([sbc_inc[g] for g in conf["genes"]]),
                                  np.array([g in sbc_set for g in conf["genes"]]))
    assert float(ub["auc"]) == a_bc[0] and float(ub["rho"]) == rho_bc[0]
    # --- L2 pre-flight refusals (each creates nothing) ---
    cp = os.path.join(l1, "cells_L1.csv")
    with open(cp, "rb") as fh:
        orig = fh.read()
    with open(cp, "ab") as fh:
        fh.write(b"\n")
    refused_creates_nothing("L2", 1, msha)                            # an L1 output changed
    with open(cp, "wb") as fh:
        fh.write(orig)
    refused_creates_nothing("L2", 1, "0" * 64)                        # wrong manifest sha

    def tamper(name, how):
        p = ctx.H(name)
        with open(p, "rb") as fh:
            keep = fh.read()
        if how == "append":
            with open(p, "ab") as fh:
                fh.write(b"x")
        else:
            os.rename(p, p + ".away")
        refused_creates_nothing("L2", 1, msha)
        if how != "append":
            os.rename(p + ".away", p)
        with open(p, "wb") as fh:
            fh.write(keep)
    tamper(stage_h_name("ucec", "V_C3"), "append")                    # wrong hash
    tamper(cap_name("luad", 77), "append")                            # wrong hash, a placebo
    tamper(cap_name("luad", 119), "remove")                           # listed but missing
    tamper(covariates_record_name("luad"), "remove")                  # the covariates record missing
    # format slips, unlisted files and case lists: the pre-flight itself, with a manifest dict
    closed = parse_closed(ctx.closed_path, ctx.cohorts)
    man = parse_manifest(ctx.H(MANIFEST_NAME))
    assert set(preflight_stage_h(ctx, closed, man)) == {"ucec", "luad"}
    for drop in (cap_name("luad", 1), stage_h_name("luad", "V_R20"), covariates_record_name("ucec")):
        _expect_fatal(preflight_stage_h, ctx, closed, {k: v for k, v in man.items() if k != drop})   # unlisted
    uni = load_conf(ctx, "ucec", genes_only=True)["genes"]
    cases = np.array(["C%03d" % i for i in range(36)])
    variants = {
        "f32": dict(genes=uni, M=np.zeros((len(uni), ctx.B + 1), np.float32), cases=cases),
        "order": dict(genes=uni[::-1], M=np.zeros((len(uni), ctx.B + 1)), cases=cases),
        "cols": dict(genes=uni, M=np.zeros((len(uni), ctx.B)), cases=cases),
        "cases": dict(genes=uni, M=np.zeros((len(uni), ctx.B + 1)), cases=np.array(["D%03d" % i for i in range(36)])),
        "keys": dict(genes=uni, M=np.zeros((len(uni), ctx.B + 1)), cases=cases, extra=np.zeros(1)),
    }
    target = stage_h_name("ucec", "A_C3")
    keep_path = ctx.H(target) + ".keep"
    shutil.copyfile(ctx.H(target), keep_path)
    for nm, payload in variants.items():
        np.savez_compressed(ctx.H(target), **payload)
        _expect_fatal(preflight_stage_h, ctx, closed, dict(man, **{target: sha256(ctx.H(target))}))
    shutil.move(keep_path, ctx.H(target))
    rn = covariates_record_name("ucec")
    shutil.copyfile(ctx.H(rn), ctx.H(rn) + ".keep")
    with np.load(ctx.H(rn)) as z:
        recd = {k: z[k] for k in z.files}
    for nm, payload in (("noC3", {k: v for k, v in recd.items() if k != "C3"}),
                        ("K", dict(recd, pcs=np.zeros((36, 7)))),
                        ("cases", dict(recd, cases=np.array(["D%03d" % i for i in range(36)])))):
        np.savez_compressed(ctx.H(rn), **payload)
        _expect_fatal(preflight_stage_h, ctx, closed, dict(man, **{rn: sha256(ctx.H(rn))}))
    shutil.move(ctx.H(rn) + ".keep", ctx.H(rn))
    assert not os.path.exists(os.path.join(ctx.out_dir, "L2.partial"))
    assert not [f for f in os.listdir(ctx.out_dir) if f.startswith("L2")]
    # --- L2 runs after two interrupted attempts (unlimited reruns) ---
    os.makedirs(os.path.join(ctx.out_dir, "L2_failed1"))
    os.makedirs(os.path.join(ctx.out_dir, "L2.partial"))
    do_stage(ctx, "L2", 1, msha)
    assert os.path.isdir(os.path.join(ctx.out_dir, "L2_failed2"))
    l2 = os.path.join(ctx.out_dir, "L2")
    with open(os.path.join(l2, "classes.json"), encoding="utf-8") as fh:
        cl = json.load(fh)
    with open(os.path.join(l2, "result_L2.json"), encoding="utf-8") as fh:
        r2 = json.load(fh)
    u, lu = cl["cohorts"]["ucec"], cl["cohorts"]["luad"]
    dkeys = ("D_0", "D_null_mean", "D_null_sd", "z_D", "p_D", "f")
    assert (u["C3_class"], u["R20_class"], u["type"], u["flagged"], u["U"]) == ("H", "closed", "NC:R", True, False), u
    assert u["e"] is None and u["p_abs"] is None and u["R20"]["placebos"] == "closed"
    assert u["R20"]["m_R_25_columns"] is None and u["R20"]["status"] == "closed"
    assert all(u["C3"]["D"][k] is not None for k in dkeys) and all(u["R20"]["D"][k] is None for k in dkeys)
    assert (lu["C3_class"], lu["R20_class"], lu["type"], lu["flagged"]) == ("closed", "H", "NC:C", True), lu
    assert lu["R20"]["placebos"] == "ok" and lu["e"] is not None and 0 < lu["p_abs"] <= 1
    assert len(lu["R20"]["m_j"]) == 119 and lu["R20"]["m_R_25_columns"] is not None
    assert lu["R20"]["margin_all_1000_columns"] is not None
    assert all(lu["R20"]["D"][k] is not None for k in dkeys) and all(lu["C3"]["D"][k] is None for k in dkeys)
    assert lu["reverse"]["C3"]["status"] == "closed" and lu["reverse"]["C3"]["outranked"] is None
    assert u["reverse"]["C3"]["status"] == "ok" and isinstance(u["reverse"]["C3"]["outranked"], bool)
    assert u["reverse"]["R20"]["status"] == "closed" and isinstance(lu["reverse"]["R20"]["outranked"], bool)
    assert cl["case"] == 2 and sorted(cl["flagged"]) == ["luad", "ucec"]
    assert cl["M4"]["forward"]["n_cells"] == 6 and cl["M4"]["forward"]["t"] == 5        # luad T_C3 closed
    assert cl["M4"]["matched"]["n_cells"] == 4 and cl["M4"]["matched"]["t"] == 2        # ucec A_R20, luad A_C3
    assert sorted(cl["closed_parts"]) == sorted(tokens)
    # diagnostic: canonical correlations from the covariates record, against a direct call, and the Gaussian
    # reference block RandomState(42).randn(n, K) with C3 and with R20
    with np.load(ctx.H(covariates_record_name("ucec"))) as z:
        want = cancorr(z["pcs"], z["C3"])
        want20 = cancorr(z["pcs"], z["R20"])
        gref = np.random.RandomState(42).randn(z["pcs"].shape[0], z["pcs"].shape[1])
        gwant, gwant20 = cancorr(gref, z["C3"]), cancorr(gref, z["R20"])
    assert np.array_equal(np.array(u["diagnostic"]["cancorr_pcs_C3"]), want)
    assert np.array_equal(np.array(r2["diagnostic"]["ucec"]["cancorr_pcs_R20"]), want20)
    assert np.array_equal(np.array(u["diagnostic"]["cancorr_gauss_C3"]), gwant)
    assert np.array_equal(np.array(r2["diagnostic"]["ucec"]["cancorr_gauss_R20"]), gwant20)
    assert u["diagnostic"]["gauss_seed"] == 42
    assert lu["diagnostic"]["cancorr_pcs_C3"] is None and lu["diagnostic"]["cancorr_pcs_R20"] is not None
    assert lu["diagnostic"]["cancorr_gauss_C3"] is None and lu["diagnostic"]["cancorr_gauss_R20"] is not None
    assert u["markers_confirmatory"] == {"stroma": {"used": ["COL1A1"]}}
    assert u["S_BC"]["n"] == 60 and lu["S_BC"]["rho_available"] is False
    c2 = read_cells(os.path.join(l2, "cells_L2.csv"))
    assert len(c2) == 2 * (2 + 119 + 2 + 2 + 2 + 2 + 2 + 2 + 2)
    reads = [r for r in c2 if r["read"] == "True"]
    assert sorted(set(r["family"] for r in reads)) == ["matched", "primary", "reverse"] and len(reads) == 12
    fam = collections.Counter(r["family"] for r in c2)
    assert fam["capacity"] == 2 * N_CAP and fam["rule reference"] == 4 and fam["FDR-defined"] == 4
    assert fam["gene property"] == 4 and fam["covariate-to-covariate"] == 4 and fam["batch-corrected"] == 4
    assert sorted(set(r["set"] for r in c2 if r["family"] == "FDR-defined")) == ["S_M|C3", "S_M|R20"]
    assert sorted(set(r["matrix"] for r in c2 if r["family"] == "gene property")) == ["V_C3", "V_R20"]
    assert all(r["rho"] == "" for r in c2 if r["family"] == "gene property")
    # S_BC on A_C3 and A_R20 (reported, not read); rho only where the table has a finite value for every universe gene
    bcl = [r for r in c2 if r["family"] == "batch-corrected"]
    assert all(r["set"] == "S_BC" and r["read"] == "False" and r["null"] == "unrestricted" for r in bcl)
    assert sorted((r["cohort"], r["matrix"], r["status"]) for r in bcl) == [
        ("luad", "A_C3", "closed"), ("luad", "A_R20", "ok"), ("ucec", "A_C3", "ok"), ("ucec", "A_R20", "closed")]
    assert all(r["rho"] != "" for r in bcl if r["cohort"] == "ucec" and r["status"] == "ok")
    assert all(r["rho"] == "" for r in bcl if r["cohort"] == "luad")
    # closed: per cohort the A_R20 / A_C3 arm cells (primary, matched, rule reference, FDR-defined, S_BC = 5), the
    # V arm cells (reverse, covariate-to-covariate, gene property = 3) and, in ucec, the 119 placebos
    assert collections.Counter(r["cohort"] for r in c2 if r["status"] == "closed") == {"ucec": 5 + 3 + 119, "luad": 5 + 3}
    prim = [r for r in c2 if r["family"] == "primary" and r["status"] == "ok"]
    assert len(prim) == 2 and all(r["D_0"] != "" and r["class"] != "" for r in prim)
    # an independent re-computation of one primary cell, of one S_BC cell and of m_j for one placebo
    with np.load(ctx.H(stage_h_name("ucec", "A_C3"))) as z:
        Mx = z["M"]
    S = read_sets(os.path.join(l1, "sets_ucec.csv"))
    mask = np.array([g in S["S_M"] for g in conf["genes"]])
    auc = assemble(Mx, np.arange(len(mask), dtype=float), mask)[0]
    assert abs(u["C3"]["auc"] - auc[0]) < 1e-15 and u["C3"]["p"] == emp_p(auc[0], auc[1:])
    a_bc, rho_bc, _, _ = assemble(Mx, np.array([sbc_inc[g] for g in conf["genes"]]),
                                  np.array([g in sbc_set for g in conf["genes"]]))
    row = [r for r in bcl if r["cohort"] == "ucec" and r["matrix"] == "A_C3"][0]
    assert float(row["auc"]) == a_bc[0] and float(row["rho"]) == rho_bc[0]
    confl = load_conf(ctx, "luad")
    Sl = read_sets(os.path.join(l1, "sets_luad.csv"))
    maskl = np.array([g in Sl["S_M"] for g in confl["genes"]])
    with np.load(ctx.H(cap_name("luad", 42))) as z:
        a42 = assemble(z["M"], np.arange(len(maskl), dtype=float), maskl)[0]
    assert abs(lu["R20"]["m_j"][41] - (a42[0] - a42[1:].mean())) < 1e-15
    refused_creates_nothing("L2", 1, msha)                            # completed: never rerun
    print("  end-to-end 1 (closures; types NC:R and NC:C; case 2; S_BC cells; refusals create nothing; reruns): OK")


def _selftest_end_to_end_2(tmp):
    """Scenario 2, nothing closed: ucec C3 N (type C), R20 U; luad C3 H, R20 A (type R); case 2.  Also the
    S_BC loader (frozen size, sign rule, LUAD sha256 against the registered json, rho availability) and the
    operator-indicator counts on this tree."""
    import copy
    root = os.path.join(tmp, "synthetic_root_2")
    shifts = {"ucec": dict(A_C3=-1.5, A_R20=-1.5, V_C3=1.5, V_R20=0.0),
              "luad": dict(A_C3=1.5, A_R20=-1.5, V_C3=0.0, V_R20=0.0)}
    ctx, msha = build_synthetic(root, [], shifts, cap_shift={"ucec": -3.0, "luad": 1.5}, seed=8)
    assert parse_closed(ctx.closed_path, ctx.cohorts)["tokens"] == []
    assert check_pinned(ctx, no_closed()) == {"ucec": "host_ucec", "luad": "host_luad"}
    assert all(v["ok"] for v in gate_v3(ctx).values())                  # --gate-v3 body, both score tables present
    do_stage(ctx, "L1", 1, msha)
    do_stage(ctx, "L2", 1, msha)
    with open(os.path.join(ctx.out_dir, "L2", "classes.json"), encoding="utf-8") as fh:
        cl = json.load(fh)
    u, lu = cl["cohorts"]["ucec"], cl["cohorts"]["luad"]
    assert (u["C3_class"], u["R20_class"], u["type"], u["flagged"], u["U"]) == ("N", "U", "C", True, True), u
    assert (lu["C3_class"], lu["R20_class"], lu["type"], lu["flagged"], lu["U"]) == ("H", "A", "R", True, False), lu
    assert u["e"] == 119 and u["p_abs"] == 1.0 and lu["e"] == 0 and lu["p_abs"] == 1 / 120.
    assert u["reverse"]["C3"]["outranked"] is True
    assert cl["case"] == 2 and sorted(cl["flagged"]) == ["luad", "ucec"]
    assert cl["M4"]["forward"]["t"] == 6 and cl["M4"]["matched"]["t"] == 4
    c2 = read_cells(os.path.join(ctx.out_dir, "L2", "cells_L2.csv"))
    assert not any(r["status"] == "closed" for r in c2)
    # --- the S_BC loader ---
    def alt(key, **ov):
        co = copy.deepcopy(ctx.cohorts)
        co[key].update(ov)
        return Ctx(ctx.root, co, ctx.B, True, ctx.pins)
    uni_u = load_conf(ctx, "ucec", genes_only=True)["genes"]
    uni_l = load_conf(ctx, "luad", genes_only=True)["genes"]
    bu, bl = load_sbc(ctx, "ucec", uni_u), load_sbc(ctx, "luad", uni_l)
    assert bu["n"] == 60 and bu["rho_available"] and bu["xs"] is not None and bu["rho_note"] == ""
    assert bl["n"] == 60 and not bl["rho_available"] and bl["xs"] is None and "1 without a finite" in bl["rho_note"]
    assert "G020" in bu["S"] and "G079" in bu["S"] and "G019" not in bu["S"] and "G080" not in bu["S"]
    assert "G100" not in bu["S"] and "G105" not in bu["S"]                       # the sign rule and the FDR rule
    assert bu["n_in_universe"] == len(bu["S"] & set(uni_u)) and bu["n_in_universe"] >= 50
    assert load_sbc(ctx, "ucec", list(uni_u) + ["NOT_IN_TABLE"])["xs"] is None   # an absent universe gene: no rho
    _expect_fatal(load_sbc, alt("ucec", n_bc=61), "ucec", uni_u)                 # not the frozen size
    _expect_fatal(load_sbc, alt("luad", n_bc=59), "luad", uni_l)
    _expect_fatal(load_sbc, alt("luad", sitepack_prefix="00000000"), "luad", uni_l)   # not the recorded prefix
    with open(ctx.P(ctx.cohorts["luad"]["published"]), encoding="utf-8") as fh:
        pubj = json.load(fh)
    h_ok = pubj["sha256_discovery_sitepack"]
    assert h_ok == sha256(ctx.P(ctx.cohorts["luad"]["sitepack"]))
    for nm, val in (("badsha", h_ok[:-1] + ("0" if h_ok[-1] != "0" else "1")), ("nokey", None)):
        pj = dict(pubj)
        if val is None:
            del pj["sha256_discovery_sitepack"]
        else:
            pj["sha256_discovery_sitepack"] = val
        with open(os.path.join(ctx.root, "pub_luad_%s.json" % nm), "w", encoding="utf-8") as fh:
            json.dump(pj, fh)
        _expect_fatal(load_sbc, alt("luad", published="pub_luad_%s.json" % nm), "luad", uni_l)
    # --- the operator-indicator counts of PS 2 on this tree ---
    d_u = load_inputs(ctx, "ucec")
    Dop, lv, info = check_operator(ctx, "ucec", d_u)
    assert info["n_indicators"] == 12 and info["n_covered"] == 36 and info["n_unlabelled"] == 1 and len(lv) == 12
    _expect_fatal(check_operator, alt("ucec", n_op_covered=37), "ucec", d_u)
    _expect_fatal(check_operator, alt("ucec", n_op_dummies=11), "ucec", d_u)
    print("  end-to-end 2 (nothing closed; types C (class U) and R (class A); case 2), S_BC loader, "
          "operator-indicator counts: OK")


def selftest(keep_root=None):
    print("sha256 %s c1_covariates_local.py" % sha256(os.path.abspath(__file__)))
    t0 = time.time()
    tmp = keep_root or tempfile.mkdtemp(prefix="c1cov_selftest_")
    if keep_root:
        if os.path.exists(keep_root):
            fatal("REFUSED: --synthetic-root %s exists; give a new folder" % keep_root)
        os.makedirs(keep_root)
    if os.path.normcase(os.path.abspath(tmp)).startswith(os.path.normcase(PAPER) + os.sep):
        fatal("REFUSED: the synthetic root may not lie inside MorphoResidual_paper")
    print("selftest (synthetic only) in %s" % tmp)
    _selftest_kernel()
    _selftest_cancorr()
    _selftest_sets(tmp)
    _selftest_parsers(tmp)
    _selftest_classes()
    _selftest_d_cap()
    _selftest_run_stage(tmp)
    _selftest_sidecar_list()
    _selftest_make_sidecar()
    _selftest_end_to_end_1(tmp)
    _selftest_end_to_end_2(tmp)
    if not keep_root:
        shutil.rmtree(tmp, ignore_errors=True)
    print("SELFTEST PASSED in %.0f s" % (time.time() - t0))


# ------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--selftest", action="store_true")
    mode.add_argument("--validate-only", action="store_true")
    mode.add_argument("--gate-v3", action="store_true")
    mode.add_argument("--check-stage-h", action="store_true")
    mode.add_argument("--check-sidecar", action="store_true")
    mode.add_argument("--stage", choices=["L1", "L2"])
    ap.add_argument("--synthetic-root", help="--selftest only: keep the synthetic files in this new folder")
    ap.add_argument("--expect-cols", type=int, default=B_PERM + 1,
                    help="--check-stage-h: columns of the four main arms (6 for a --perms 5 smoke export)")
    ap.add_argument("--expect-cap-cols", type=int, default=CAP_COLS + 1,
                    help="--check-stage-h: columns of the placebo files")
    ap.add_argument("--manifest-sha", help="--stage L1 / L2: sha256 of the Stage H manifest (from the Stage H record)")
    ap.add_argument("--workers", type=int, default=max(1, min(24, (os.cpu_count() or 2) - 2)))
    args = ap.parse_args()
    try:
        if args.selftest:
            return selftest(args.synthetic_root)
        if sha256(tk.__file__) != KERNEL_SHA:
            fatal("FATAL: theta_kernel.py is not the frozen kernel")
        ctx = real_ctx()
        if args.check_sidecar:
            sys.exit(0 if check_sidecar_mode(ctx) else 1)
        if args.validate_only:
            refuse_if_frozen(ctx, "--validate-only")
            out, t0 = {}, time.time()
            for key in ctx.cohorts:
                d, conf = load_inputs(ctx, key), load_conf(ctx, key)
                check_universe(d, conf)
                _, _, op_info = check_operator(ctx, key, d)                # PS 2 counts: fatal if wrong
                out[key] = dict(OP=op_info, S_BC=sbc_info(load_sbc(ctx, key, conf["genes"])))
                out[key].update(V2=validate_v2(ctx, key, d, conf), V1b=validate_v1b(ctx, key, d))
                out[key]["V1"] = validate_v1(ctx, key, d, args.workers)
            ok = all(v["V1"]["ok"] and v["V1b"]["ok"] and v["V2"]["ok"] for v in out.values())
            os.makedirs(ctx.out_dir, exist_ok=True)
            with open(os.path.join(ctx.out_dir, "validate_only.json"), "w", encoding="utf-8", newline="\n") as fh:
                json.dump(dict(ok=ok, script_sha256=sha256(ctx.script_path), kernel_sha256=sha256(ctx.kernel_path),
                               env=env_record(), workers=args.workers, minutes=round((time.time() - t0) / 60.0, 1),
                               tolerances=dict(ref=TOL_REF, r2_rna=TOL_PUB, increment=TOL_INC, p_fdr=TOL_PFDR,
                                               auc=TOL_AUC, v1b_r2_rna_pcs20=TOL_INC,
                                               v1b_incr_over_pcs20={k: c["tol_v1b_inc"]
                                                                    for k, c in ctx.cohorts.items()}),
                               cohorts=out), fh, indent=1, sort_keys=True)
            log("V1 %s; V1b %s; V2 %s" % tuple("ok" if all(v[k]["ok"] for v in out.values()) else "FAILED"
                                             for k in ("V1", "V1b", "V2")))
            log("VALIDATED" if ok else "VALIDATION FAILED")
            sys.exit(0 if ok else 1)
        if args.gate_v3:
            refuse_if_frozen(ctx, "--gate-v3")
            res = gate_v3(ctx)
            ok = all(v["ok"] for v in res.values())
            os.makedirs(ctx.out_dir, exist_ok=True)
            with open(os.path.join(ctx.out_dir, "gate_v3.json"), "w", encoding="utf-8", newline="\n") as fh:
                json.dump(dict(ok=ok, script_sha256=sha256(ctx.script_path), env=env_record(), cohorts=res), fh,
                          indent=1, sort_keys=True)
            log("V3 VALIDATED" if ok else "V3 FAILED")
            sys.exit(0 if ok else 1)
        if args.check_stage_h:
            sys.exit(0 if check_stage_h(ctx, args.expect_cols, args.expect_cap_cols) else 1)
        if not args.manifest_sha:
            ap.error("--stage %s needs --manifest-sha" % args.stage)
        do_stage(ctx, args.stage, args.workers, args.manifest_sha)
    except Fatal as e:
        print(str(e), flush=True)
        sys.exit(2)


if __name__ == "__main__":
    main()
