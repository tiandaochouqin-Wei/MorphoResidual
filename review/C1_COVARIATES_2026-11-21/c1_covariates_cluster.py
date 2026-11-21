#!/usr/bin/env python3
"""
c1_covariates_cluster.py -- cluster side (Stage H and its gates) of the POST HOC re-reads
of the two confirmatory tests with transcriptome-derived covariates.

GOVERNED BY review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md, VERSION 3.1 (sections cited
"PS n"). Nothing here changes the published C1-UCEC reading or the registered C1-LUAD reading.

WHAT IS IMPORTED UNCHANGED (sha256 asserted before import, PS 1 / PS 9)
  residual_analysis.py       load_rna_matrix, load_protein_matrix, load_wsi_embeddings,
                             cv_r2, CV_FOLDS (read MORPHO_* at import time)
  split_replication_power.py k_for (via the published prepare functions)
  c1_run_test.py             prepare() and assemble() for C1-UCEC
  c1_run_test_luad.py        prepare() (incl. the frozen 113-case list, the independence
                             check and its hash assertions), assert_hash, assemble() for C1-LUAD
  composition_controls.py    comp_score, MARKERS (the three composition scores, PS 2)
  transcriptome_baseline.py is NOT importable as a function (its R20 lines sit inside
  main()); r20_recipe re-implements exactly the lines PS 2 quotes, and gate-r20 (H-V2)
  checks it against the published rna_pcs20 arrays.

THE NEW PER-GENE FUNCTION (PS 4) is _gene_test_cov: a copy of
  split_replication_power._gene_test with two additions, baseline covariates E (stay with
  their patients) and a general block B (rows permuted under the null):
      r2b      = cv_r2(hstack([xr[v], E[v]]), y[v])
      column 0 = cv_r2(hstack([xr[v], E[v], B[v]]), y[v]) - r2b
      column b = cv_r2(hstack([xr[v], E[v], B[perms[b-1]][v]]), y[v]) - r2b
  run through a multiprocessing pool built exactly like split_replication_power.run_pool.
  With E empty and B = the published morphology components it is the published function
  (gate H-V1 checks that on the saved matrix).

MODES (--mode; --cohort ucec|luad)
  selftest        synthetic data only, no project import (run locally; must pass).
  export-patients LOCAL helper: writes the patient list of an inputs_<cohort>.npz (one id per
                  line) for upload; reads patient ids only.
  blind           confirmatory tree, counts only (PS 9 "Blind run"): no increment of any kind;
                  prints the markers found per family and their SDs (PS 2 C3 rule).
  gate            H-V1: E empty, B = published components, columns 0..5 vs the saved matrix.
                  Strict verdict (identical genes, |diff| <= 1e-9) and fallback quantities
                  (largest deviation; AUC of the published set from the recomputed column 0
                  against the published AUC), both written into the gate json.
  gate-r20        H-V2, DISCOVERY tree: R20 rebuilt on the inputs file's patients vs rna_pcs20.
  export-scores   H0, DISCOVERY tree: the three composition scores of the discovery patients.
  arm             production (after the release): one of A_C3, A_R20, V_C3, V_R20.
  placebo         production: A_R20cap_j for j in --draws a-b (1 <= a <= b <= 119), one file
                  per draw; a draw whose file exists is skipped.
  manifest        production: ONE c1_cov_stage_h_manifest.sha256 over the output folders of
                  both confirmatory trees (run through the UCEC wrapper), lines
                  "<sha256> *<bare file name>" only.

PRODUCTION GUARDS (PS 6.1, PS 9 "Scripts", "Interrupted runs")
  - --sidecar with the sha256 given by --sidecar-sha256; every sidecar entry present on the
    cluster must match: the pre-specification (--prespec), this script, the three wrappers,
    the five imported modules, the cohort's saved matrix, CLOSED_PARTS.txt (--closed-parts)
    and, in the modes arm and placebo, PINNED_HOSTS.txt (--pinned-hosts).
  - CLOSED_PARTS.txt is read strictly (PS 9 "Closed parts"); a closed arm is refused (exit 4).
  - PINNED_HOSTS.txt (--pinned-hosts; a sidecar-listed file) must give this host for the
    cohort ("pinned <cohort> <host>"), and a gate json of this host (H-V1) written by the SAME
    script sha256 on the SAME inputs (recorded input hashes) must have passed_strict or
    passed_fallback. There is no free host argument.
  - No limit on attempts. Each attempt writes <stem>.attempt<k>.json before it computes;
    partial files left by an interrupted attempt are kept as
    <file>.partial.interrupted_before_attempt<k>; a completed output is never recomputed or
    overwritten (exclusive create + hard link). A per-run lock (flock, same host) stops two
    jobs of the same run from overlapping; placebo ranges of one cohort may not overlap.
  - A stop by the scheduler (SIGTERM/SIGUSR2/SIGINT) is logged as such: resubmit unchanged.
    Any exception is logged with the PS 9 "Defects" rule: report to the drafting session
    before any resubmission; an environment failure (disk, memory, file system) is then
    resubmitted unchanged and recorded; anything else is handled as a script defect.

FILES (in MORPHO_OUT = <tree>/results/posthoc/c1_covariates)
  c1_cov_<cohort>_<arm>.npz                 genes (unicode), M (float64, G x 1001), cases
  c1_cov_<cohort>_A_R20cap_<jjj>.npz        genes, M (G x 26: col 0 + perms 1-25), cases
  c1_cov_confirmatory_covariates_<cohort>.npz  cases, C3 (n x 3; n x 0 if C3 is closed),
                                            R20 (n x 20), pcs (n x K), markers_json
  <stem>.attempt<k>.json / .attempt<k>.log / .run<k>.json (/ .attempt<k>.exception.json)
  c1_cov_<cohort>_gate_hv1_<host>_<stamp>.json, c1_cov_<cohort>_blind_<host>_<stamp>.json ...
  c1_cov_stage_h_manifest.sha256            (UCEC tree; covers both trees)
  discovery tree: c1_cov_discovery_scores_<cohort>.tsv, gate_r20 / export_scores jsons.

SELF-TEST  --mode selftest prints "sha256 <hex> c1_covariates_cluster.py" as its first line.

EXIT CODES  0 ok / gate passed strictly; 3 gate not passed (H-V1: strict test failed, the
  fallback verdict is in the json); 4 refused (closed arm, completed output, overlapping
  range, a leftover file, an incomplete output set for the manifest); 143 stopped by the
  scheduler; anything else = a crash (FATAL text on stderr, or a Python traceback).

RUN only through lsf_c1_cov_ucec.sh / lsf_c1_cov_luad.sh / lsf_c1_cov_discovery.sh, which
set MORPHO_* before python starts. See SERVER_RUNBOOK.md.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
import platform
import re
import signal
import socket
import sys
import tempfile
import time
import traceback
import zlib

import numpy as np

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------------------
# Pinned constants (PS 1, 2, 4, 9). Changing any of these changes the analysis.
# ---------------------------------------------------------------------------------------
PRESPEC_REL = "review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md"
SCRIPT_REL = "review/C1_COVARIATES_2026-11-21/c1_covariates_cluster.py"
CLOSED_PARTS_REL = "review/C1_COVARIATES_2026-11-21/CLOSED_PARTS.txt"
PINNED_HOSTS_REL = "review/C1_COVARIATES_2026-11-21/PINNED_HOSTS.txt"
WRAPPER_REL = {
    "ucec": "review/C1_COVARIATES_2026-11-21/lsf_c1_cov_ucec.sh",
    "luad": "review/C1_COVARIATES_2026-11-21/lsf_c1_cov_luad.sh",
    "discovery": "review/C1_COVARIATES_2026-11-21/lsf_c1_cov_discovery.sh",
}

# PS 1 prints 8-hex prefixes; the full values are those of the repository copies (the
# prefixes are asserted below at import time). The cluster copies must be byte-identical.
FROZEN_CODE_SHA256 = {
    "residual_analysis": "a54a0a56494d315e15ef075425fe2133ecd9f05193c22656c1145bffdaac4a0a",
    "split_replication_power": "79530c075c0f7e75bf06d3c2f1edfe8792cbf4d6c71ae1d32e60bf44f31a0397",
    "c1_run_test": "fc578c705e803ddb118ff1f65a61079723167985a66adf7285834bf4709391b5",
    "c1_run_test_luad": "1d5ee05d57f80af4d90870aecb10e574e4c910a8ce1eab743a067bd9a9278742",
    "composition_controls": "46b4a1d8a6973bb71b162697a04754e0e1c555eeaf3ff83dfc5c143087abdfa6",
}
PS1_PREFIX = {
    "residual_analysis": "a54a0a56", "split_replication_power": "79530c07",
    "c1_run_test": "fc578c70", "c1_run_test_luad": "1d5ee05d",
    "composition_controls": "46b4a1d8",
}
FROZEN_CODE_SIDECAR_REL = {
    "residual_analysis": "server_export/scripts/residual_analysis.py",
    "split_replication_power": "split_replication_power.py",
    "c1_run_test": "server_export/scripts/c1_run_test.py",
    "c1_run_test_luad": "server_export/scripts/c1_run_test_luad.py",
    "composition_controls": "server_export/scripts/composition_controls.py",
}
THETA_KERNEL_SHA256 = "6c95cd119d28e811e09967946e9f191c69f07f54f653991917e817ef99ec2734"

ZW = "/public/home/fjhui/ZW"
COHORTS = {
    "ucec": dict(
        root_base="ucec_c1", disc_root_base="ucec", n=138, K=28, n_tested=10153,
        n_usable=10146, perm_seed=b"c1_ucec|perm", saved_key="M_unrestricted",
        saved_sha="ef75f1c06d4fc33843a79e9a8d1bb4f4224acd201d30784badb97cf3be13384a",
        saved_rel="server_export/scripts/pinned/c1_ucec_d4_stratified.npz",
        inputs_sha="e28b86fadd41ec86363577fffad1db1f2f8a5162e732c7acb0e4133c0982044e",
        n_disc=100, n_sel_usable=2536, n_bg_usable=7610,
        # published quantity: arms.unrestricted.auc.obs of pinned/c1_ucec_d4_followup.json
        published_auc=0.5974425564900907),
    "luad": dict(
        root_base="luad_c1", disc_root_base="luad", n=112, K=22, n_tested=10417,
        n_usable=10406, perm_seed=b"c1_luad|perm", saved_key="M",
        saved_sha="2034e71143c2e065675d3159832b1a7be17aa38c7feeaacfc82862c411ea2aab",
        saved_rel="server_export/scripts/pinned/luad_c1/c1_luad_primary_primary.npz",
        inputs_sha="b64cdf3b31717a32e0d42c42f3f1dde74786d5598008f0d214297a71adaba095",
        n_disc=105, n_sel_usable=2287, n_bg_usable=8119,
        # published quantity: h1.auc_obs of pinned/luad_c1/c1_luad_primary_result_primary.json
        published_auc=0.6439668501223573),
}
# C1-UCEC: path defaults of c1_run_test.py / c1_ucec_plex_reread.py (the published runs'
# paths); the two hashes are the plex re-read's pinned values.
UCEC_DEFAULTS = dict(
    protein=ZW + "/ucec_c1/protein_ucec_c1.csv",
    discovery_results=ZW + "/scripts/pinned/ucec/residual_results_tumoronly.csv",
    discovery_slide_map=ZW + "/scripts/pinned/slide_type_map_ucec.tsv",
    case_ids=ZW + "/scripts/pinned/c1_case_ids_ucec.txt",
    saved_matrix=ZW + "/ucec_c1/results/c1_ucec_d4_stratified.npz",
)
UCEC_INPUT_SHA256 = {
    "case_ids": "d777fc88c1dd7b9eb817b0a6a1366b1b909386bd34dea857e4bc53cffc5807dc",
    "discovery_results": "100d8f8581c78905c2a98caf4bb2713e376d9264c6f00276446a9893b987c4a6",
}
LUAD_RULE_SHA256 = "a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022"
LUAD_REQUIRED = ("rule", "rule_sha256", "case_ids", "protein", "protein_sha256",
                 "rna_manifest_sha256", "discovery_results", "discovery_sitepack",
                 "discovery_slide_map", "discovery_xwalk", "discovery_rna_manifest",
                 "batch_labels", "saved_matrix")

N_PERM = 1000             # PS 4: the published unrestricted permutations
N_GATE_PERMS = 5          # PS 9 H-V1: columns 0..5
N_CAP = 119               # PS 2: R20cap_j, j = 1..119
N_CAP_PERMS = 25          # PS 4: placebo arms use columns 0 and 1-25
N_R20 = 20                # PS 2
MIN_MARKERS = 6           # PS 2
MIN_N = 30                # PS 4 (= c1_run_test.MIN_N)
HV1_TOL = 1e-9            # PS 9 H-V1 strict
HV1_FALLBACK_DEV = 1e-6   # PS 9 H-V1 fallback: largest deviation below 1e-6 ...
HV1_FALLBACK_AUC = 1e-6   # ... and |AUC(S_M, recomputed column 0) - published AUC| <= 1e-6
HV2_TOL = 1e-6            # PS 9
SD_TINY = 1e-8            # PS 9 blind run
FAMILIES = ("stroma", "immune", "proliferation")   # PS 2, in this order
ARMS = {                  # PS 4: (E, B); None = no columns
    "A_C3": ("C3", "pcs"),
    "A_R20": ("R20", "pcs"),
    "V_C3": (None, "C3"),
    "V_R20": (None, "R20"),
}
CLOSABLE_ARMS = ("A_C3", "A_R20", "A_R20cap", "V_C3", "V_R20")
CLOSABLE_ANALYSES = ("C3", "M|C3")      # PS 9: no other discovery analysis can be closed
COVARIATES_KEYS = ("cases", "C3", "R20", "pcs", "markers_json")
MANIFEST_NAME = "c1_cov_stage_h_manifest.sha256"
REQUIRED_ENV = ("MORPHO_ROOT", "MORPHO_RNA_MANIFEST", "MORPHO_SLIDE_MAP",
                "MORPHO_WSI_EMB_DIR", "MORPHO_OUT")
DEFAULT_LIST_ROOTS = [ZW + "/ucec_c1/results", ZW + "/luad_c1/results", ZW + "/results",
                      ZW + "/ucec/results", ZW + "/luad/results", ZW + "/meas_error",
                      ZW + "/c1_covariates"]
MANIFEST_DIRS = {"ucec": ZW + "/ucec_c1/results/posthoc/c1_covariates",
                 "luad": ZW + "/luad_c1/results/posthoc/c1_covariates"}
LISTING_FLAG = re.compile(r"(comp|txome|transcriptome|rna_?pc|pcs20|c1_cov|c1_spec|operator|"
                          r"batchonly|stroma|immune|prolif|purity|estimate|xcell)", re.I)

for _k, _v in FROZEN_CODE_SHA256.items():
    assert _v.startswith(PS1_PREFIX[_k]), _k


class Refused(Exception):
    """A refusal that is not a crash: exit code 4."""


class SchedulerStop(BaseException):
    """Raised by the signal handler when the scheduler stops a production run."""


# ---------------------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------------------
def sha256(path):
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return "unreadable"
    return h.hexdigest()


def fatal(msg):
    sys.exit("FATAL: " + msg)


def utc(t=None):
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t))


def stamp():
    """UTC time plus process id, so two runs started in the same second never collide."""
    t = time.time()
    return time.strftime("%Y%m%dT%H%M%SZ", time.gmtime(t)) + "-%d-%06d" % (
        os.getpid(), int((t % 1) * 1e6))


def host():
    return socket.gethostname().split(".")[0]


def jdefault(o):
    if hasattr(o, "tolist"):
        return o.tolist()
    if isinstance(o, bytes):
        return o.decode("latin-1")
    return str(o)


class Tee:
    """print() plus, once open() is called, an append-mode copy in a log file."""

    def __init__(self):
        self.fh = None

    def open(self, path):
        self.fh = open(path, "a", encoding="utf-8", newline="\n")

    def __call__(self, msg=""):
        print(msg, flush=True)
        if self.fh is not None:
            self.fh.write(str(msg) + "\n")
            self.fh.flush()

    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None


# ---------------------------------------------------------------------------------------
# Exclusive writes, locks, attempts (PS 4 Files, PS 6.1, PS 9 "Interrupted runs")
# ---------------------------------------------------------------------------------------
def _publish(partial, final):
    """Make `partial` appear as `final` without ever replacing an existing `final`."""
    try:
        os.link(partial, final)
    except FileExistsError:
        raise Refused("refusing to overwrite " + final)
    except (AttributeError, OSError):
        if os.path.exists(final):
            raise Refused("refusing to overwrite " + final)
        os.rename(partial, final)
        return
    os.remove(partial)


def write_exclusive(final, writer, partial=None):
    """writer(fh) writes bytes to a fresh partial file (exclusive create; default
    <final>.partial); on success the partial is published as <final>. An existing <final>
    or partial is refused."""
    if os.path.exists(final):
        raise Refused("refusing to overwrite " + final)
    partial = partial or final + ".partial"
    try:
        fd = os.open(partial, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0),
                     0o644)
    except FileExistsError:
        raise Refused("a leftover " + partial + " exists; see SERVER_RUNBOOK")
    with os.fdopen(fd, "wb") as fh:
        writer(fh)
        fh.flush()
        os.fsync(fh.fileno())
    _publish(partial, final)
    return sha256(final)


def write_text_exclusive(final, text):
    return write_exclusive(final, lambda fh: fh.write(text.encode("utf-8")))


def write_json_exclusive(final, obj):
    return write_text_exclusive(final, json.dumps(obj, indent=2, default=jdefault) + "\n")


def save_matrix_npz(final, genes, M, cases):
    g = np.array(list(genes))
    c = np.array(list(cases))
    if g.dtype.kind != "U" or c.dtype.kind != "U":
        raise ValueError("genes/cases must be numpy unicode arrays")
    M = np.asarray(M, dtype=np.float64)
    if M.shape[0] != len(g) or not np.isfinite(M).all():
        raise ValueError("matrix shape/finite check failed before writing " + final)
    return write_exclusive(final, lambda fh: np.savez_compressed(fh, genes=g, M=M, cases=c))


class RunLock:
    """An exclusive flock on <out_dir>/.c1_cov_lock_<stem> held for the whole run. Production
    runs only on the pinned host, so a host-local lock is enough to stop a resubmitted job
    from touching the files of a job of the same run that is still alive. Not available on
    Windows (self-test only): then it is a no-op."""

    def __init__(self, out_dir, stem):
        self.path = os.path.join(out_dir, ".c1_cov_lock_" + stem)
        self.fh = None

    def acquire(self):
        try:
            import fcntl
        except ImportError:
            return False
        self.fh = open(self.path, "a")
        try:
            fcntl.flock(self.fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self.fh.close()
            self.fh = None
            raise Refused("another job of this run holds %s (is it still running? check bjobs; "
                          "resubmit only after it has ended)" % self.path)
        return True

    def release(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None


ATTEMPT_RE = r"^%s\.attempt(\d+)\.json$"


def attempts_recorded(out_dir, stem):
    pat = re.compile(ATTEMPT_RE % re.escape(stem))
    return sorted(int(m.group(1)) for m in (pat.match(f) for f in os.listdir(out_dir)) if m)


def begin_attempt(out_dir, stem, finals, extra_partials, record):
    """PS 6.1 / PS 9 "Interrupted runs". Refuses if every final output exists (a completed
    run is never repeated). Otherwise this is attempt k = (highest recorded attempt) + 1, with
    no limit. Partial files left by an interrupted attempt (of the finals and the run's own
    extra partials, e.g. its covariates-record partial) are kept under a distinct name,
    <partial>.interrupted_before_attempt<k>. Writes <stem>.attempt<k>.json (exclusive)."""
    if finals and all(os.path.exists(f) for f in finals):
        raise Refused("every output of " + stem + " exists: a completed run is never repeated")
    done = attempts_recorded(out_dir, stem)
    k = (max(done) + 1) if done else 1
    partials = [f + ".partial" for f in finals] + list(extra_partials)
    kept = []
    for p in partials:
        if os.path.exists(p):
            if k == 1:
                raise Refused("leftover " + p + " without any attempt record of " + stem +
                              "; stop and report")
            dst = p + ".interrupted_before_attempt%d" % k
            if os.path.exists(dst):
                raise Refused(dst + " exists already; stop and report")
            os.rename(p, dst)
            kept.append(os.path.basename(dst))
    rec = dict(record)
    rec.update(attempt=k, previous_attempts=done, interrupted_partials_kept=kept,
               attempt_start_utc=utc())
    write_json_exclusive(os.path.join(out_dir, stem + ".attempt%d.json" % k), rec)
    return k, kept


# ---------------------------------------------------------------------------------------
# CLOSED_PARTS.txt (PS 6.1, PS 9)
# ---------------------------------------------------------------------------------------
def parse_closed_parts_text(text):
    """CLOSED_PARTS.txt, parsed strictly (PS 9 "Closed parts"; the local script parses it
    identically). A line starting with '#' in column 0 is a comment; an empty line is
    ignored; every other line is exactly three tokens separated by single spaces, with no
    leading or trailing whitespace: 'arm <cohort> <arm>' (cohort ucec|luad; arm A_C3 | A_R20 |
    A_R20cap | V_C3 | V_R20) or 'analysis <cohort> <X>' (X: C3 or M|C3). A line of
    whitespace only is not an empty line and is refused. Consistency (per cohort): A_R20cap is
    listed if and only if A_R20 is; 'analysis C3' if and only if 'analysis M|C3' is. Raises
    ValueError on anything else."""
    arms, analyses = set(), set()
    for i, line in enumerate(text.split("\n"), 1):
        if line == "" or line.startswith("#"):
            continue
        parts = line.split(" ")
        if len(parts) != 3 or parts[1] not in COHORTS:
            raise ValueError("CLOSED_PARTS line %d not understood: %r" % (i, line))
        kind, cohort, what = parts
        if kind == "arm" and what in CLOSABLE_ARMS:
            item, bucket = (cohort, what), arms
        elif kind == "analysis" and what in CLOSABLE_ANALYSES:
            item, bucket = (cohort, what), analyses
        else:
            raise ValueError("CLOSED_PARTS line %d not understood: %r" % (i, line))
        if item in bucket:
            raise ValueError("CLOSED_PARTS line %d repeats %r" % (i, line))
        bucket.add(item)
    for cohort in sorted(COHORTS):
        if ((cohort, "A_R20cap") in arms) != ((cohort, "A_R20") in arms):
            raise ValueError("CLOSED_PARTS: for %s, 'arm %s A_R20cap' is listed if and only if "
                             "'arm %s A_R20' is (PS 9 'Closed parts')" % (cohort, cohort, cohort))
        if ((cohort, "C3") in analyses) != ((cohort, "M|C3") in analyses):
            raise ValueError("CLOSED_PARTS: for %s, 'analysis %s C3' is listed if and only if "
                             "'analysis %s M|C3' is (PS 9 'Closed parts')"
                             % (cohort, cohort, cohort))
    return dict(arms=arms, analyses=analyses)


def _read_control_file(path, what, flag):
    if not path or not os.path.isfile(path):
        fatal("%s must name %s (got %r)" % (flag, what, path))
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw.startswith(b"\xef\xbb\xbf") or b"\r" in raw:
        fatal(path + " must be UTF-8 without BOM, LF line endings")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError as e:
        fatal("%s is not UTF-8 (%s)" % (path, e))


def read_closed_parts(path):
    text = _read_control_file(path, "CLOSED_PARTS.txt", "--closed-parts")
    try:
        return parse_closed_parts_text(text)
    except ValueError as e:
        fatal(str(e))


def parse_pinned_hosts_text(text):
    """PINNED_HOSTS.txt, parsed strictly (PS 9 "Scripts", H-V1). A line starting with '#' in
    column 0 is a comment; an empty line is ignored; every other line is exactly three
    tokens separated by single spaces, with no leading or trailing whitespace:
    'pinned <cohort> <host>' (cohort ucec|luad), at most one line per cohort. Returns
    {cohort: host}. Raises ValueError on anything else."""
    pinned = {}
    for i, line in enumerate(text.split("\n"), 1):
        if line == "" or line.startswith("#"):
            continue
        parts = line.split(" ")
        if (len(parts) != 3 or parts[0] != "pinned" or parts[1] not in COHORTS
                or not re.fullmatch(r"\S+", parts[2])):
            raise ValueError("PINNED_HOSTS line %d not understood: %r" % (i, line))
        if parts[1] in pinned:
            raise ValueError("PINNED_HOSTS line %d: a second line for %s" % (i, parts[1]))
        pinned[parts[1]] = parts[2]
    return pinned


def read_pinned_hosts(path):
    text = _read_control_file(path, "PINNED_HOSTS.txt", "--pinned-hosts")
    try:
        return parse_pinned_hosts_text(text)
    except ValueError as e:
        fatal(str(e))


def check_pinned_host(pinned, cohort, this_host):
    """Production runs only on the host PINNED_HOSTS.txt gives for the cohort."""
    if cohort not in pinned:
        fatal("PINNED_HOSTS.txt has no 'pinned %s <host>' line: no host is pinned for %s, so "
              "no production run is possible (PS 9)" % (cohort, cohort))
    if pinned[cohort] != this_host:
        fatal("this host is %s, but PINNED_HOSTS.txt pins %s for %s; production runs only on "
              "the pinned host (PS 9)" % (this_host, pinned[cohort], cohort))


def arm_closed(closed, cohort, arm):
    return (cohort, arm) in closed["arms"]


def closed_summary(closed):
    return dict(arms=sorted("%s %s" % a for a in closed["arms"]),
                analyses=sorted("%s %s" % a for a in closed["analyses"]))


# ---------------------------------------------------------------------------------------
# Manifest (PS 4 Files, PS 9 Stage H)
# ---------------------------------------------------------------------------------------
def cap_name(cohort, j):
    return "c1_cov_%s_A_R20cap_%03d.npz" % (cohort, j)


def covariates_name(cohort):
    return "c1_cov_confirmatory_covariates_%s.npz" % cohort


def expected_outputs(cohort, closed):
    """The final outputs a cohort must have before the manifest is written: every arm not
    closed, every placebo draw unless A_R20cap is closed, and the covariates record unless
    every arm of the cohort is closed."""
    names = ["c1_cov_%s_%s.npz" % (cohort, a) for a in ARMS if not arm_closed(closed, cohort, a)]
    if not arm_closed(closed, cohort, "A_R20cap"):
        names += [cap_name(cohort, j) for j in range(1, N_CAP + 1)]
    if names:
        names.append(covariates_name(cohort))
    return names


def closed_outputs(cohort, closed):
    names = ["c1_cov_%s_%s.npz" % (cohort, a) for a in ARMS if arm_closed(closed, cohort, a)]
    if arm_closed(closed, cohort, "A_R20cap"):
        names += [cap_name(cohort, j) for j in range(1, N_CAP + 1)]
    return names


def build_manifest_text(dirs_by_cohort, closed):
    """Lines '<sha256> *<file name>' for every regular file named c1_cov_* in the output
    folders of both cohorts (manifest files excluded), sorted by name. Refuses while any
    *.partial exists, on a file name that occurs in two folders, when an expected output is
    missing, or when an output of a closed arm exists. No status lines (PS 4, PS 9)."""
    entries, where = {}, {}
    for cohort, d in sorted(dirs_by_cohort.items()):
        if not os.path.isdir(d):
            raise Refused("manifest directory does not exist: " + d)
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            if not f.startswith("c1_cov_") or not os.path.isfile(p) or f.startswith(MANIFEST_NAME):
                continue
            if f.endswith(".partial"):
                raise Refused(p + " is in progress (or an interrupted run left it: resubmit that "
                              "run unchanged, PS 6.1); no manifest now")
            if f in entries:
                raise Refused("file name " + f + " occurs in two manifest directories")
            entries[f] = sha256(p)
            where[f] = cohort
    missing, forbidden = [], []
    for cohort in sorted(dirs_by_cohort):
        missing += [n for n in expected_outputs(cohort, closed)
                    if n not in entries or where[n] != cohort]
        forbidden += [n for n in closed_outputs(cohort, closed) if n in entries]
    if forbidden:
        raise Refused("outputs of closed arms exist: %s; stop and report" % forbidden[:5])
    if missing:
        raise Refused("%d expected output(s) missing, e.g. %s: every arm not closed must "
                      "complete first (resubmit unchanged, PS 6.1)" % (len(missing), missing[:5]))
    lines = ["%s *%s" % (entries[f], f) for f in sorted(entries)]
    return "\n".join(lines) + "\n", entries


def write_manifest(final, text, force_rewrite):
    """Exclusive creation. A leftover <manifest>.partial (an interrupted manifest run) blocks
    the write; --force-rewrite moves it aside, and is accepted only while no manifest exists."""
    if os.path.exists(final):
        raise Refused(final + " exists; a manifest is written once and never rewritten")
    partial = final + ".partial"
    if os.path.exists(partial):
        if not force_rewrite:
            raise Refused("leftover %s; rerun with --force-rewrite (allowed only while no "
                          "manifest exists)" % partial)
        os.rename(partial, partial + ".interrupted_" + stamp())
    return write_text_exclusive(final, text)


# ---------------------------------------------------------------------------------------
# Definitions (PS 2, PS 4)
# ---------------------------------------------------------------------------------------
def cap_seed(cohort, j):
    """PS 2: crc32(b"c1_cov|<cohort>|cap|%03d" % j)."""
    return zlib.crc32(b"c1_cov|%s|cap|%03d" % (cohort.encode("ascii"), j))


def cap_perm(cohort, j, n):
    return np.random.RandomState(cap_seed(cohort, j)).permutation(n)


def placebo_from(R20, q):
    """PS 2: row i of the placebo is row q[i] of R20."""
    return R20[q]


def published_perms(seed_bytes, n, B=N_PERM):
    """PS 4: rng = RandomState(crc32(seed)); perms = [rng.permutation(n) for _ in range(B)],
    with no other draw from that generator."""
    rng = np.random.RandomState(zlib.crc32(seed_bytes))
    return [rng.permutation(n) for _ in range(B)]


def r20_recipe(rna_df):
    """PS 2, the lines of transcriptome_baseline.py (re-implemented; not importable):
    rv = rna.values.astype(float); sd = rv.std(0);
    rz = (rv - rv.mean(0)) / np.where(sd == 0, 1.0, sd);
    R20 = PCA(n_components=20, svd_solver="full").fit_transform(rz)."""
    from sklearn.decomposition import PCA
    rv = rna_df.values.astype(float)
    sd = rv.std(0)
    rz = (rv - rv.mean(0)) / np.where(sd == 0, 1.0, sd)
    pca = PCA(n_components=N_R20, svd_solver="full")
    R20 = pca.fit_transform(rz)
    return R20, pca.explained_variance_ratio_.astype(np.float64)


def build_c3(rna_df, comp_score, markers):
    """PS 2: composition_controls.comp_score per family, columns stroma, immune,
    proliferation. C3 is untestable if a family has fewer than MIN_MARKERS markers present
    (as columns of rna) OR no present marker with a non-zero standard deviation (pandas,
    ddof 1, as comp_score). Returns (C3 or None, info); info holds, per family, the markers
    listed, present, missing, the SD of every present marker, the markers that enter the
    score (non-zero SD) and whether the family passes."""
    info = {}
    cols = []
    testable = True
    for fam in FAMILIES:
        present = [m for m in markers[fam] if m in rna_df.columns]
        sds = {m: float(rna_df[m].std()) for m in present}
        nonzero = [m for m in present if sds[m] != 0]
        fam_ok = len(present) >= MIN_MARKERS and len(nonzero) >= 1
        info[fam] = dict(n_listed=len(markers[fam]), n_present=len(present), present=present,
                         missing=[m for m in markers[fam] if m not in rna_df.columns],
                         sd=sds, n_zero_sd=len(present) - len(nonzero), used=nonzero,
                         family_testable=bool(fam_ok))
        if not fam_ok:
            testable = False
            continue
        s, mk = comp_score(rna_df, markers[fam])
        if s is None or list(mk) != present:
            raise ValueError("comp_score disagrees with the marker count of family " + fam)
        cols.append(np.asarray(s, dtype=np.float64).reshape(-1, 1))
    if not testable:
        return None, info
    C3 = np.hstack(cols)
    if not np.isfinite(C3).all():
        raise ValueError("non-finite composition score although every family has a marker "
                         "with non-zero SD")
    return C3, info


def rna_column_counts(rna_df):
    rv = rna_df.values.astype(float)
    return dict(rna_columns=int(rv.shape[1]), all_zero_columns=int((rv == 0).all(0).sum()),
                **{"sd_below_1e-8_columns": int((rv.std(0) < SD_TINY).sum())})


def max_diff_sign_aligned(A, ref):
    """max |s_k A[:,k] - ref[:,k]| with s_k = sign of <A[:,k], ref[:,k]> (PCA signs)."""
    A = np.asarray(A, dtype=np.float64)
    ref = np.asarray(ref, dtype=np.float64)
    if A.shape != ref.shape:
        return float("inf"), []
    signs = []
    worst = 0.0
    for k in range(A.shape[1]):
        s = 1.0 if float(A[:, k] @ ref[:, k]) >= 0 else -1.0
        signs.append(int(s))
        worst = max(worst, float(np.max(np.abs(s * A[:, k] - ref[:, k]))))
    return worst, signs


def log_markers(log, minfo):
    for fam in FAMILIES:
        i = minfo[fam]
        log("  markers %-13s: %d/%d present, %d zero-SD, family %s"
            % (fam, i["n_present"], i["n_listed"], i["n_zero_sd"],
               "passes" if i["family_testable"] else "FAILS the PS 2 rule"))
        log("      SD (ddof 1) : " + (" ".join("%s=%.6g" % (m, i["sd"][m]) for m in i["present"])
                                       or "-"))
        log("      missing     : " + (",".join(i["missing"]) or "-"))


# ---------------------------------------------------------------------------------------
# The per-gene function (PS 4) and its pool (as split_replication_power.run_pool)
# ---------------------------------------------------------------------------------------
_W = {}


def _init_cov(rna, prot, E, B, perms, min_n, cv_fn, folds):
    if mp.current_process().name != "MainProcess":
        # pool workers: default signal behaviour (the stop handlers belong to the main process)
        for name in ("SIGTERM", "SIGUSR2", "SIGINT", "SIGXCPU"):
            sig = getattr(signal, name, None)
            if sig is not None:
                try:
                    signal.signal(sig, signal.SIG_DFL)
                except (OSError, ValueError):
                    pass
    _W.update(rna=rna, prot=prot, E=E, B=B, perms=perms, min_n=min_n, cv=cv_fn, folds=folds)


def _gene_test_cov(j):
    """split_replication_power._gene_test with baseline covariates E and a block B."""
    y = _W["prot"][:, j]
    v = ~np.isnan(y)
    if v.sum() < _W["min_n"]:
        return j, None
    yv = y[v]
    xr = _W["rna"][v, j].reshape(-1, 1)
    E = _W["E"]
    B = _W["B"]
    cv = _W["cv"]
    folds = _W["folds"]
    Ev = E[v]
    r2b = cv(np.hstack([xr, Ev]), yv, folds)
    out = np.empty(len(_W["perms"]) + 1)
    out[0] = cv(np.hstack([xr, Ev, B[v]]), yv, folds) - r2b
    for b, p in enumerate(_W["perms"], 1):
        out[b] = cv(np.hstack([xr, Ev, B[p][v]]), yv, folds) - r2b
    return j, out


def run_pool_cov(n_genes, workers, rna, prot, E, B, perms, min_n, cv_fn, folds):
    """workers >= 1: mp.Pool exactly as split_replication_power.run_pool. workers == 0:
    serial (self-test only)."""
    if workers <= 0:
        _init_cov(rna, prot, E, B, perms, min_n, cv_fn, folds)
        return dict(_gene_test_cov(j) for j in range(n_genes))
    with mp.Pool(workers, initializer=_init_cov,
                 initargs=(rna, prot, E, B, perms, min_n, cv_fn, folds)) as pool:
        return dict(pool.imap_unordered(_gene_test_cov, range(n_genes), chunksize=64))


def usable_indices(tested, res, disc_incr):
    """The published usable filter (c1_run_test.run_test / c1_run_test_luad.run_family):
    >= MIN_N protein values (res not None) and a finite discovery increment."""
    return [j for j in range(len(tested)) if res[j] is not None
            and np.isfinite(disc_incr.get(tested[j], np.nan))]


def usable_by_count(tested, P, disc_incr, min_n=MIN_N):
    """Same filter from counts only (blind run): no increment is computed."""
    nn = (~np.isnan(P)).sum(axis=0)
    return [j for j in range(len(tested)) if nn[j] >= min_n
            and np.isfinite(disc_incr.get(tested[j], np.nan))]


def matrix_from(tested, res, usable):
    genes_u = [tested[j] for j in usable]
    M = np.vstack([res[j] for j in usable]) if usable else np.empty((0, 0))
    return genes_u, M


def check_matrix(genes_u, M, saved_genes, ncols):
    if list(genes_u) != list(saved_genes):
        raise ValueError("usable gene list/order differs from the saved matrix's genes "
                         "(%d vs %d genes)" % (len(genes_u), len(saved_genes)))
    if M.shape != (len(saved_genes), ncols):
        raise ValueError("matrix shape %r, expected %r" % (M.shape, (len(saved_genes), ncols)))
    if not np.isfinite(M).all():
        raise ValueError("non-finite value in the matrix (fatal, PS 4)")


def compute_matrix(tested, R, P, E, B, perms, disc_incr, saved_genes, workers, cv_fn, folds):
    """One Stage H matrix (PS 4): every usable gene of the published test, columns 0 and
    1..len(perms). Raises ValueError if the gene list, shape or finiteness check fails."""
    res = run_pool_cov(len(tested), workers, R, P, E, B, perms, MIN_N, cv_fn, folds)
    usable = usable_indices(tested, res, disc_incr)
    genes_u, M = matrix_from(tested, res, usable)
    check_matrix(genes_u, M, saved_genes, len(perms) + 1)
    return genes_u, M


# ---------------------------------------------------------------------------------------
# Gate H-V1 verdicts (PS 9)
# ---------------------------------------------------------------------------------------
def hv1_auc(assemble_fn, M, genes, selected, disc_incr):
    """AUC of the published selection set from column 0 of M (the published `assemble`)."""
    sel = np.array([g in selected for g in genes])
    xs = np.array([float(disc_incr[g]) for g in genes])
    if sel.sum() == 0 or (~sel).sum() == 0:
        return float("nan"), int(sel.sum()), int((~sel).sum())
    auc, _rho, n1, n0 = assemble_fn(np.asarray(M)[:, :1], xs, sel)
    return float(auc[0]), int(n1), int(n0)


def hv1_verdicts(same_genes, M, saved_M, auc_new, auc_published):
    """passed_strict: identical gene list and columns 0..5 within HV1_TOL.
    passed_fallback: identical gene list, largest deviation < HV1_FALLBACK_DEV and
    |AUC - published| <= HV1_FALLBACK_AUC."""
    M = np.asarray(M, dtype=np.float64)
    saved_M = np.asarray(saved_M, dtype=np.float64)
    if same_genes and M.shape == saved_M.shape and np.isfinite(M).all():
        d = np.abs(M - saved_M)
        dmax = float(d.max()) if d.size else 0.0
        d999 = float(np.percentile(d, 99.9)) if d.size else 0.0
        dcol = [float(x) for x in d.max(axis=0)] if d.size else []
    else:
        dmax = d999 = float("inf")
        dcol = []
    auc_dev = abs(float(auc_new) - float(auc_published)) if np.isfinite(auc_new) else float("inf")
    strict = bool(same_genes and dmax <= HV1_TOL)
    fallback = bool(same_genes and dmax < HV1_FALLBACK_DEV and auc_dev <= HV1_FALLBACK_AUC)
    return dict(passed_strict=strict, passed_fallback=fallback, max_abs_diff=dmax,
                p999_abs_diff=d999, max_abs_diff_by_column=dcol, auc_abs_dev=auc_dev)


HV1_MATCH_KEYS = ("sha256_this_script", "sha256_saved_matrix", "sha256_inputs",
                  "sha256_frozen_modules", "cohort", "host")


def find_accepted_hv1(out_dir, cohort, this_host, expect):
    """PS 9: production accepts the pinned host only if a gate json of THIS host, written by
    the SAME script sha256 on the SAME inputs, has passed_strict or passed_fallback. Returns
    (file name, record) or (None, list of reasons)."""
    pat = re.compile(r"^c1_cov_%s_gate_hv1_%s_\d{8}T\d{6}Z-\d+-\d{6}\.json$"
                     % (re.escape(cohort), re.escape(this_host)))
    reasons = []
    for f in sorted(os.listdir(out_dir)):
        if not pat.match(f):
            continue
        try:
            with open(os.path.join(out_dir, f), encoding="utf-8") as fh:
                rec = json.load(fh)
        except (OSError, ValueError) as e:
            reasons.append("%s unreadable (%s)" % (f, e))
            continue
        bad = [k for k in HV1_MATCH_KEYS if rec.get(k) != expect.get(k)]
        if rec.get("gate") != "H-V1":
            bad.append("gate")
        if bad:
            reasons.append("%s: differs in %s" % (f, ",".join(bad)))
            continue
        if rec.get("passed_strict") is True or rec.get("passed_fallback") is True:
            return f, rec
        reasons.append("%s: neither passed_strict nor passed_fallback" % f)
    return None, reasons


# ---------------------------------------------------------------------------------------
# Covariates record (PS 4 Files)
# ---------------------------------------------------------------------------------------
def covariates_payload(cases, C3, R20, pcs, minfo, c3_closed):
    n = len(cases)
    C3a = np.empty((n, 0), dtype=np.float64) if (c3_closed or C3 is None) else \
        np.asarray(C3, dtype=np.float64)
    meta = dict(families=list(FAMILIES), c3_stored="closed" if C3a.shape[1] == 0 else "computed",
                markers=minfo)
    return dict(cases=np.array(list(cases)), C3=C3a, R20=np.asarray(R20, dtype=np.float64),
                pcs=np.asarray(pcs, dtype=np.float64),
                markers_json=np.array(json.dumps(meta, sort_keys=True, default=jdefault)))


def same_payload(path, payload):
    with np.load(path, allow_pickle=False) as z:
        if sorted(z.files) != sorted(COVARIATES_KEYS):
            return False
        for k in COVARIATES_KEYS:
            a, b = z[k], payload[k]
            if a.shape != b.shape or a.dtype.kind != b.dtype.kind or not np.array_equal(a, b):
                return False
    return True


def covariates_partial(out_dir, stem):
    return os.path.join(out_dir, stem + ".covariates_record.partial")


def ensure_covariates_file(out_dir, cohort, stem, payload, log):
    """c1_cov_confirmatory_covariates_<cohort>.npz, written by whichever production job gets
    there first, through a partial file owned by this run (<stem>.covariates_record.partial,
    so an interrupted write is kept by the run's next attempt); every job checks that its
    covariates are identical to the published record."""
    final = os.path.join(out_dir, covariates_name(cohort))
    if not os.path.exists(final):
        partial = covariates_partial(out_dir, stem)
        try:
            write_exclusive(final, lambda fh: np.savez_compressed(fh, **payload), partial=partial)
            log("  wrote " + os.path.basename(final))
        except Refused:
            if os.path.exists(partial):     # another job published first: drop our copy
                os.remove(partial)
            if not os.path.exists(final):
                raise
    if not same_payload(final, payload):
        fatal(final + " exists and differs from this job's covariates (keys " +
              ",".join(COVARIATES_KEYS) + "; different host or inputs?)")
    log("  covariates identical to " + os.path.basename(final))
    return sha256(final)


# ---------------------------------------------------------------------------------------
# Discovery score table (H0)
# ---------------------------------------------------------------------------------------
def scores_tsv_text(patients, C3, comments):
    lines = ["# " + c for c in comments]
    lines.append("patient\t" + "\t".join(FAMILIES))
    for p, row in zip(patients, C3):
        lines.append(p + "\t" + "\t".join(repr(float(x)) for x in row))
    return "\n".join(lines) + "\n"


def read_scores_tsv(path):
    rows = [ln.rstrip("\r\n").split("\t") for ln in open(path, encoding="utf-8")
            if ln.strip() and not ln.startswith("#")]
    if rows[0] != ["patient"] + list(FAMILIES):
        raise ValueError("bad header " + repr(rows[0]))
    return [r[0] for r in rows[1:]], np.array([[float(x) for x in r[1:]] for r in rows[1:]])


def read_patient_file(path):
    with open(path, encoding="utf-8") as fh:
        ids = [ln.strip() for ln in fh if ln.strip()]
    if ids != sorted(ids) or len(set(ids)) != len(ids):
        fatal(path + " must hold unique patient ids in sorted order")
    return ids


def load_inputs_patients(npz_path, cohort):
    got = sha256(npz_path)
    if got != COHORTS[cohort]["inputs_sha"]:
        fatal("%s sha256=%s, PS 1 expects %s" % (npz_path, got, COHORTS[cohort]["inputs_sha"]))
    z = np.load(npz_path, allow_pickle=False)
    return z, [str(p) for p in z["patients"].tolist()]


# ---------------------------------------------------------------------------------------
# Listing of result folders (PS 0 item 7) and embedding counts (blind run)
# ---------------------------------------------------------------------------------------
def list_tree(roots, max_depth=3, cap=50000):
    rows, truncated = [], False
    for root in roots:
        if not os.path.isdir(root):
            rows.append((root, "missing", "", ""))
            continue
        base_depth = root.rstrip("/").count("/")
        for dirpath, dirnames, filenames in os.walk(root):
            depth = dirpath.rstrip("/").count("/") - base_depth
            dirnames.sort()
            if depth >= max_depth:
                dirnames[:] = []
            for name in sorted(filenames) + [d + "/" for d in dirnames]:
                p = os.path.join(dirpath, name.rstrip("/"))
                try:
                    st = os.stat(p)
                    rows.append((p, str(st.st_size), utc(st.st_mtime),
                                 "flag" if LISTING_FLAG.search(name) else ""))
                except OSError:
                    rows.append((p, "unreadable", "", ""))
                if len(rows) >= cap:
                    return rows, True
    return rows, truncated


def embedding_counts(RA, cases):
    import pandas as pd
    slide_df = pd.read_csv(RA.SLIDE_TYPE_MAP, sep="\t")
    stype = slide_df.set_index("slide_submitter_id")["sample_type"].to_dict()
    pt = sorted(s for s, t in stype.items() if t == "Primary Tumor")
    pt_case = sorted(s for s in pt if s[:9] in set(cases))
    slide_re = re.compile(r"^(C3[A-Z]-\d{5}-\d+)")
    files = sorted(os.path.basename(str(f)) for f in RA.WSI_EMB_DIR.glob("*.pt"))
    found_pt, n_normal, n_unmatched = set(), 0, 0
    for f in files:
        m = slide_re.match(f[:-3])
        t = stype.get(m.group(1)) if m else None
        if t == "Primary Tumor":
            found_pt.add(m.group(1))
        elif t is None:
            n_unmatched += 1
        else:
            n_normal += 1
    miss = [s for s in pt if s not in found_pt]
    miss_case = [s for s in pt_case if s not in found_pt]
    return dict(pt_files_found=len(files), primary_tumor_slides_in_map=len(pt),
                primary_tumor_slides_with_embedding=len(found_pt),
                primary_tumor_slides_missing=len(miss), missing_examples=miss[:20],
                analysed_case_pt_slides_expected=len(pt_case),
                analysed_case_pt_slides_missing=len(miss_case),
                files_non_tumor=n_normal, files_unmatched=n_unmatched)


# ---------------------------------------------------------------------------------------
# Pre-flight: environment, frozen modules, sidecar
# ---------------------------------------------------------------------------------------
def env_guard(cohort, tree):
    env = {}
    for k in REQUIRED_ENV:
        v = os.environ.get(k, "")
        if not v:
            fatal(k + " is not set. Run through the lsf_c1_cov_*.sh wrappers (residual_analysis "
                  "reads MORPHO_* at import time).")
        env[k] = v
    root = os.path.normpath(env["MORPHO_ROOT"])
    out = os.path.normpath(env["MORPHO_OUT"])
    if tree == "confirmatory":
        want = COHORTS[cohort]["root_base"]
        if os.path.basename(root) != want:
            fatal("basename(MORPHO_ROOT)=%r, expected %r" % (os.path.basename(root), want))
        if not out.startswith(root + os.sep) or os.path.basename(out) != "c1_covariates":
            fatal("MORPHO_OUT=%r must be <MORPHO_ROOT>/.../c1_covariates" % out)
    else:
        want = COHORTS[cohort]["disc_root_base"]
        if os.path.basename(root) != want:
            fatal("basename(MORPHO_ROOT)=%r, expected the discovery tree %r"
                  % (os.path.basename(root), want))
        if os.path.basename(out) != cohort or \
                os.path.basename(os.path.dirname(out)) != "discovery":
            fatal("MORPHO_OUT=%r must be .../c1_covariates/discovery/%s" % (out, cohort))
    return env


def import_frozen(env, log):
    """Hash the five frozen modules from the files the imports will load, then import them
    (the precedent c1_ucec_plex_reread.import_frozen)."""
    import importlib
    import importlib.util
    if _SCRIPT_DIR not in sys.path:
        sys.path.insert(0, _SCRIPT_DIR)
    origins, hashes = {}, {}
    for name, want in FROZEN_CODE_SHA256.items():
        spec = importlib.util.find_spec(name)
        if spec is None or not spec.origin:
            fatal("cannot locate module %s on sys.path %r" % (name, sys.path[:3]))
        origins[name] = os.path.realpath(spec.origin)
        got = sha256(origins[name])
        if got != want:
            fatal("sha256 mismatch for frozen module %s: got=%s want=%s (%s)"
                  % (name, got, want, origins[name]))
        hashes[name] = got
        log("  hash OK: frozen module %s = %s... (%s)" % (name, got[:12], origins[name]))
    mods = {}
    for name in ("residual_analysis", "split_replication_power", "c1_run_test",
                 "c1_run_test_luad", "composition_controls"):
        mods[name] = importlib.import_module(name)
        if os.path.realpath(mods[name].__file__) != origins[name]:
            fatal("imported %s from %s but hashed %s" % (name, mods[name].__file__, origins[name]))
    RA = mods["residual_analysis"]
    for attr, key in (("ROOT", "MORPHO_ROOT"), ("RNA_MANIFEST", "MORPHO_RNA_MANIFEST"),
                      ("SLIDE_TYPE_MAP", "MORPHO_SLIDE_MAP"),
                      ("WSI_EMB_DIR", "MORPHO_WSI_EMB_DIR"), ("OUT_DIR", "MORPHO_OUT")):
        if os.path.normpath(str(getattr(RA, attr))) != os.path.normpath(env[key]):
            fatal("residual_analysis.%s=%s != $%s=%s" % (attr, getattr(RA, attr), key, env[key]))
    for m in ("c1_run_test", "c1_run_test_luad", "split_replication_power",
              "composition_controls"):
        if getattr(mods[m], "RA", RA) is not RA:
            fatal(m + " imported a different residual_analysis module object")
    if mods["c1_run_test"].MIN_N != MIN_N or mods["c1_run_test_luad"].MIN_N != MIN_N:
        fatal("published MIN_N differs from PS 4 (30)")
    return mods, hashes, origins


def read_sidecar(path):
    entries = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2 or not re.match(r"^[0-9a-fA-F]{64}$", parts[0]):
                fatal("malformed sidecar line: " + repr(line))
            rel = parts[1].strip()
            if rel.startswith("*"):
                rel = rel[1:]
            rel = rel.replace("\\", "/")
            if rel in entries:
                fatal("the sidecar lists %s twice" % rel)
            entries[rel] = parts[0].lower()
    return entries


def check_sidecar_entries(entries, cluster_files):
    """Every file of `cluster_files` ({sidecar path: cluster path}) must be listed in the
    sidecar and have the listed sha256. Returns (checked, entries not on the cluster)."""
    checked = {}
    for rel, path in sorted(cluster_files.items()):
        if rel not in entries:
            fatal("the sidecar does not list " + rel)
        if not path or not os.path.isfile(path):
            fatal("%s (sidecar entry %s) is not on the cluster" % (path, rel))
        g = sha256(path)
        if g != entries[rel]:
            fatal("sidecar mismatch for %s: %s has %s, the sidecar says %s"
                  % (rel, path, g, entries[rel]))
        checked[rel] = g
    return checked, sorted(set(entries) - set(checked))


def cluster_sidecar_files(args, cohort, origins):
    """{sidecar path: cluster path} of every sidecar entry a cluster run uses (PS 9
    "Scripts"): this file and the pre-specification, the three wrappers, the five imported
    modules, the cohort's saved matrix, CLOSED_PARTS.txt and (when given) PINNED_HOSTS.txt."""
    cluster = {PRESPEC_REL: args.prespec, SCRIPT_REL: os.path.abspath(__file__),
               COHORTS[cohort]["saved_rel"]: args.saved_matrix,
               CLOSED_PARTS_REL: args.closed_parts}
    if args.pinned_hosts:
        cluster[PINNED_HOSTS_REL] = args.pinned_hosts
    for rel in WRAPPER_REL.values():
        cluster[rel] = os.path.join(_SCRIPT_DIR, os.path.basename(rel))
    for name, rel in FROZEN_CODE_SIDECAR_REL.items():
        cluster[rel] = origins[name]
    return cluster


def verify_sidecar(args, cohort, origins, log):
    """PS 9 "Scripts": the sidecar has the sha256 given on the command line, and every
    sidecar entry present on the cluster matches: the pre-specification, this script, the
    three wrappers, the five imported modules, the cohort's saved matrix, CLOSED_PARTS.txt
    and, in the production modes arm and placebo (and in the manifest mode if given),
    PINNED_HOSTS.txt. Every other entry (local inputs, local scripts, gate logs, Abstract
    variants) is listed as 'not on the cluster'."""
    for a in ("sidecar", "sidecar_sha256", "prespec", "closed_parts"):
        if not getattr(args, a):
            fatal("production needs --sidecar, --sidecar-sha256, --prespec and --closed-parts")
    if args.mode in ("arm", "placebo") and not args.pinned_hosts:
        fatal("production needs --pinned-hosts (PINNED_HOSTS.txt, listed in the sidecar)")
    if not re.match(r"^[0-9a-fA-F]{64}$", args.sidecar_sha256):
        fatal("--sidecar-sha256 must be 64 hex digits")
    got = sha256(args.sidecar)
    if got != args.sidecar_sha256.lower():
        fatal("sidecar sha256 %s != --sidecar-sha256 %s" % (got, args.sidecar_sha256))
    entries = read_sidecar(args.sidecar)
    wrapper = os.environ.get("C1COV_WRAPPER", "")
    own = os.path.join(_SCRIPT_DIR, os.path.basename(WRAPPER_REL[cohort]))
    if not wrapper or os.path.realpath(wrapper) != os.path.realpath(own):
        fatal("C1COV_WRAPPER=%r is not %s; run through the wrapper beside this script"
              % (wrapper, own))
    cluster = cluster_sidecar_files(args, cohort, origins)
    checked, not_here = check_sidecar_entries(entries, cluster)
    log("  sidecar %s... verified: %d entries match on this cluster, %d entries are not on "
        "the cluster (not checked here)" % (got[:12], len(checked), len(not_here)))
    return dict(sidecar_sha256=got, checked=checked, not_on_cluster=not_here)


# ---------------------------------------------------------------------------------------
# Confirmatory preparation (published prepare functions, unchanged)
# ---------------------------------------------------------------------------------------
def resolve_paths(args):
    if args.cohort == "ucec":
        for k, v in UCEC_DEFAULTS.items():
            if not getattr(args, k):
                setattr(args, k, v)
    else:
        miss = [k for k in LUAD_REQUIRED if not getattr(args, k)]
        if miss:
            fatal("C1-LUAD needs every published path/hash flag; missing: " +
                  ", ".join("--" + m.replace("_", "-") for m in miss))


def confirmatory_prepare(args, mods, env, log):
    """Calls the published prepare() of the cohort's script unchanged; returns a dict."""
    cohort = args.cohort
    c = COHORTS[cohort]
    hashes = {}
    if cohort == "ucec":
        C1 = mods["c1_run_test"]
        for key in ("case_ids", "discovery_results"):
            g = sha256(getattr(args, key))
            if g != UCEC_INPUT_SHA256[key]:
                fatal("sha256 mismatch for %s: got=%s want=%s" % (key, g, UCEC_INPUT_SHA256[key]))
            hashes[key] = g
            log("  hash OK: %s = %s..." % (key, g[:12]))
        # recorded here; production compares them with the H-V1 record of this host
        hashes["protein"] = sha256(args.protein)
        hashes["rna_manifest"] = sha256(env["MORPHO_RNA_MANIFEST"])
        hashes["slide_map"] = sha256(env["MORPHO_SLIDE_MAP"])
        hashes["discovery_slide_map"] = sha256(args.discovery_slide_map)
        common, tested, selected, disc_incr, disc_fdr, R, P, W, K = C1.prepare(args, log)
        with open(args.case_ids, encoding="utf-8") as fh:
            frozen_cases = sorted({ln.strip() for ln in fh if ln.strip()})
        if list(common) != frozen_cases:
            fatal("the analysed case list differs from the frozen case-id file")
        assemble_fn = C1.assemble
    else:
        C1L = mods["c1_run_test_luad"]
        g = sha256(args.rule)
        if g != args.rule_sha256.lower() or g != LUAD_RULE_SHA256:
            fatal("--rule sha256 %s does not equal --rule-sha256 and the signed value %s"
                  % (g, LUAD_RULE_SHA256))
        hashes["rule"] = g
        for key, path, want in (("protein", args.protein, args.protein_sha256),
                                ("rna_manifest", env["MORPHO_RNA_MANIFEST"],
                                 args.rna_manifest_sha256)):
            g = sha256(path)
            if g != want.lower():
                fatal("sha256 mismatch for %s: got=%s want=%s (%s)" % (key, g, want, path))
            hashes[key] = g
            log("  hash OK: %s = %s... (%s)" % (key, g[:12], path))
        for key, path in (("discovery_results", args.discovery_results),
                          ("discovery_sitepack", args.discovery_sitepack),
                          ("batch_labels", args.batch_labels)):
            hashes[key] = C1L.assert_hash(path, key, log)
        (common, tested, selected, _sel_h2, disc_incr, disc_fdr, _bres, R, P, W,
         K) = C1L.prepare(args, log)
        for key, path in (("case_ids", args.case_ids),
                          ("discovery_slide_map", args.discovery_slide_map),
                          ("discovery_xwalk", args.discovery_xwalk),
                          ("discovery_rna_manifest", args.discovery_rna_manifest),
                          ("slide_map", env["MORPHO_SLIDE_MAP"])):
            hashes[key] = sha256(path)
        assemble_fn = C1L.assemble
    n = len(common)
    for key, got in (("n", n), ("K", K), ("n_tested", len(tested))):
        if got != c[key]:
            fatal("published preparation gave %s=%s, expected %s" % (key, got, c[key]))
    if list(common) != sorted(common):
        fatal("cases are not in sorted order")
    log("  prepare(): n=%d K=%d genes tested=%d (as published)" % (n, K, len(tested)))
    return dict(common=list(common), tested=list(tested), selected=set(selected),
                disc_incr=disc_incr, R=R, P=P, W=W, K=K, input_hashes=hashes,
                assemble=assemble_fn)


def load_saved(args, cohort, log, columns=None):
    c = COHORTS[cohort]
    g = sha256(args.saved_matrix)
    if g != c["saved_sha"]:
        fatal("saved matrix %s sha256=%s, PS 1 expects %s" % (args.saved_matrix, g, c["saved_sha"]))
    z = np.load(args.saved_matrix, allow_pickle=False)
    genes = [str(x) for x in z["genes"].tolist()]
    if len(genes) != c["n_usable"]:
        fatal("saved matrix has %d genes, expected %d" % (len(genes), c["n_usable"]))
    M = None
    if columns is not None:
        M = np.asarray(z[c["saved_key"]][:, :columns], dtype=np.float64)
    log("  hash OK: saved matrix = %s... (%s; key %s; %d genes)"
        % (g[:12], args.saved_matrix, c["saved_key"], len(genes)))
    return genes, M, g


def morphology_pcs(W, K):
    from sklearn.decomposition import PCA
    return PCA(n_components=K, svd_solver="full").fit_transform(W)


def confirmatory_covariates(RA, CC, prep, log):
    """PS 2 on the confirmatory cases: rna = RA.load_rna_matrix() restricted to the
    analysed cases (every gene; the target gene is not removed)."""
    rna = RA.load_rna_matrix()
    rc = rna.loc[prep["common"]]
    chk = rna.loc[prep["common"], prep["tested"]].values.astype(float)
    if chk.shape != prep["R"].shape or not np.array_equal(chk, prep["R"]):
        fatal("the second RNA load differs from the matrix prepare() used")
    try:
        C3, minfo = build_c3(rc, CC.comp_score, CC.MARKERS)
    except ValueError as e:
        fatal(str(e))
    R20, evr = r20_recipe(rc)
    if not np.isfinite(R20).all():
        fatal("non-finite R20")
    counts = rna_column_counts(rc)
    log_markers(log, minfo)
    log("  C3 %s; R20 %s, variance explained by 20 components %.4f; RNA columns %d (all-zero "
        "%d, SD<1e-8 %d)" % ("testable" if C3 is not None else "UNTESTABLE (marker rule, PS 2)",
                              R20.shape, float(evr.sum()), counts["rna_columns"],
                              counts["all_zero_columns"], counts["sd_below_1e-8_columns"]))
    return dict(C3=C3, R20=R20, evr=evr, markers=minfo, counts=counts)


# ---------------------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------------------
def base_record(args, mode):
    return dict(script="c1_covariates_cluster.py", governed_by=PRESPEC_REL, post_hoc=True,
                mode=mode, cohort=args.cohort, host=host(), fqdn=socket.gethostname(),
                lsf_job=os.environ.get("LSB_JOBID", ""), start_utc=utc(),
                argv=sys.argv[1:], args=vars(args),
                sha256_this_script=sha256(os.path.abspath(__file__)),
                wrapper=os.environ.get("C1COV_WRAPPER", ""),
                sha256_wrapper=sha256(os.environ["C1COV_WRAPPER"])
                if os.environ.get("C1COV_WRAPPER") else None,
                versions=versions())


def versions():
    v = dict(python=platform.python_version(), numpy=np.__version__)
    for m in ("scipy", "sklearn", "pandas"):
        try:
            v[m] = __import__(m).__version__
        except Exception:  # noqa: BLE001
            v[m] = None
    v["threads"] = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                   "MKL_NUM_THREADS")}
    return v


def finish_record(rec, t0):
    rec["end_utc"] = utc()
    rec["minutes"] = round((time.time() - t0) / 60.0, 2)
    return rec


def parse_draws(s):
    m = re.match(r"^(\d{1,3})-(\d{1,3})$", s or "")
    if not m:
        fatal("--draws must look like 1-30")
    a, b = int(m.group(1)), int(m.group(2))
    if not (1 <= a <= b <= N_CAP):
        fatal("--draws must lie within 1-%d" % N_CAP)
    return a, b


def placebo_stem(cohort, a, b):
    return "c1_cov_%s_A_R20cap_%03d-%03d" % (cohort, a, b)


def check_draw_overlap(out_dir, cohort, a, b):
    """Placebo ranges of one cohort may not overlap (each draw belongs to one run)."""
    pat = re.compile(r"^c1_cov_%s_A_R20cap_(\d{3})-(\d{3})\.attempt\d+\.json$" % re.escape(cohort))
    for f in sorted(os.listdir(out_dir)):
        m = pat.match(f)
        if m:
            a2, b2 = int(m.group(1)), int(m.group(2))
            if (a2, b2) != (a, b) and a2 <= b and a <= b2:
                raise Refused("draws %d-%d overlap the earlier placebo run %03d-%03d (%s); "
                              "resubmit that run unchanged instead" % (a, b, a2, b2, f))


def scheduler_stop_text(signum, k, run_stem):
    return ("STOPPED BY THE SCHEDULER (signal %s) during attempt %d of %s: not a script "
            "defect. Resubmit the same command unchanged (only the wall-clock limit may be "
            "raised; PS 6.1)." % (signum, k, run_stem))


def exception_text(k, run_stem):
    return ("EXCEPTION in attempt %d of %s (PS 9 'Defects'): REPORT TO THE DRAFTING SESSION "
            "BEFORE ANY RESUBMISSION. Do not resubmit yet. An environment failure (disk, "
            "memory, file system) is resubmitted unchanged and the exception is recorded; "
            "anything else is handled as a script defect." % (k, run_stem))


def _raise_stop(signum, _frame):
    raise SchedulerStop(signum)


def install_stop_handlers():
    for name in ("SIGTERM", "SIGUSR2", "SIGINT", "SIGXCPU"):
        sig = getattr(signal, name, None)
        if sig is not None:
            try:
                signal.signal(sig, _raise_stop)
            except (OSError, ValueError):
                pass


def mode_confirmatory(args, log):
    t0 = time.time()
    cohort = args.cohort
    c = COHORTS[cohort]
    env = env_guard(cohort, "confirmatory")
    out_dir = env["MORPHO_OUT"]
    os.makedirs(out_dir, exist_ok=True)
    resolve_paths(args)
    tag = "%s_%s_%s" % (args.mode.replace("-", ""), host(), stamp())
    production = args.mode in ("arm", "placebo")
    closed = None
    if production:
        pinned = read_pinned_hosts(args.pinned_hosts)
        check_pinned_host(pinned, cohort, host())
        closed = read_closed_parts(args.closed_parts)
        if args.mode == "arm":
            if args.arm not in ARMS:
                fatal("--arm must be one of " + ", ".join(ARMS))
            if arm_closed(closed, cohort, args.arm):
                raise Refused("arm %s %s is closed in CLOSED_PARTS.txt (PS 6.1, PS 9); it is "
                              "not run" % (cohort, args.arm))
            run_stem = "c1_cov_%s_%s" % (cohort, args.arm)
            finals = [os.path.join(out_dir, run_stem + ".npz")]
        else:
            if arm_closed(closed, cohort, "A_R20cap"):
                raise Refused("arm %s A_R20cap is closed in CLOSED_PARTS.txt (PS 6.1, PS 9); no "
                              "placebo is run" % cohort)
            a, b = parse_draws(args.draws)
            run_stem = placebo_stem(cohort, a, b)
            finals = [os.path.join(out_dir, cap_name(cohort, j)) for j in range(a, b + 1)]
            check_draw_overlap(out_dir, cohort, a, b)
        if all(os.path.exists(f) for f in finals):
            raise Refused("every output of %s exists: a completed run is never repeated" % run_stem)
    rec = base_record(args, args.mode)
    rec.update(env=env)
    if not production:
        log.open(os.path.join(out_dir, "c1_cov_%s_%s.log" % (cohort, tag)))
    log("=" * 78)
    log("C1 COVARIATES (POST HOC, %s) -- cohort=%s mode=%s host=%s lsf_job=%s start=%s"
        % (PRESPEC_REL, cohort, args.mode, host(), rec["lsf_job"], rec["start_utc"]))
    log("  this script sha256 = " + rec["sha256_this_script"])
    log("=" * 78)
    mods, code_hashes, origins = import_frozen(env, log)
    RA, CC = mods["residual_analysis"], mods["composition_controls"]
    rec["sha256_frozen_modules"] = code_hashes
    if production:
        rec["sidecar"] = verify_sidecar(args, cohort, origins, log)
        rec["closed_parts"] = closed_summary(closed)
        rec["pinned_hosts"] = pinned
        log("  CLOSED_PARTS: %s" % json.dumps(rec["closed_parts"]))
        log("  PINNED_HOSTS: %s (this host: %s)" % (json.dumps(pinned), host()))
    saved_genes, saved_M, saved_sha = load_saved(
        args, cohort, log, columns=(1 + N_GATE_PERMS) if args.mode == "gate" else None)
    rec["sha256_saved_matrix"] = saved_sha
    for k in REQUIRED_ENV:
        log("  %-20s: %s" % (k, env[k]))
    prep = confirmatory_prepare(args, mods, env, log)
    rec["sha256_inputs"] = prep["input_hashes"]
    n, K, tested = len(prep["common"]), prep["K"], prep["tested"]
    rec["counts"] = dict(n_cases=n, K=K, n_tested=len(tested))

    if production:
        expect = {k: rec.get(k) for k in HV1_MATCH_KEYS}
        hv1, why = find_accepted_hv1(out_dir, cohort, host(), expect)
        if hv1 is None:
            fatal("no H-V1 gate json of host %s, written by this script (%s...) on these inputs, "
                  "with passed_strict or passed_fallback, in %s: %s"
                  % (host(), rec["sha256_this_script"][:12], out_dir, why or "none found"))
        rec["hv1_record"] = dict(file=hv1, passed_strict=why["passed_strict"],
                                 passed_fallback=why["passed_fallback"],
                                 max_abs_diff=why.get("max_abs_diff"))
        log("  H-V1 record accepted for this host: %s (strict %s, fallback %s, max |diff| %s)"
            % (hv1, why["passed_strict"], why["passed_fallback"], why.get("max_abs_diff")))

    if args.mode == "blind":
        return mode_blind(args, log, RA, CC, prep, saved_genes, rec, out_dir, tag, t0)

    W = prep["W"]
    pcs = morphology_pcs(W, K)
    perms = published_perms(c["perm_seed"], n, N_PERM)
    workers = args.workers if args.workers else RA.N_WORKERS
    cv_fn, folds = RA.cv_r2, RA.CV_FOLDS

    if args.mode == "gate":
        return mode_gate(args, log, prep, pcs, perms, saved_genes, saved_M, workers, cv_fn,
                         folds, rec, out_dir, tag, t0)

    # ---------------- production ----------------
    cov = confirmatory_covariates(RA, CC, prep, log)
    c3_closed = arm_closed(closed, cohort, "A_C3") and arm_closed(closed, cohort, "V_C3")
    if cov["C3"] is None and not c3_closed:
        fatal("C3 is untestable in %s by the marker rule of PS 2, but CLOSED_PARTS.txt does not "
              "close both A_C3 and V_C3 (the blind run should have shown this before freezing): "
              "stop and report" % cohort)
    payload = covariates_payload(prep["common"], cov["C3"], cov["R20"], pcs, cov["markers"],
                                 c3_closed)
    lock = RunLock(out_dir, run_stem)
    lock.acquire()
    try:
        k, kept = begin_attempt(out_dir, run_stem, finals, [covariates_partial(out_dir, run_stem)],
                                dict(rec, note="attempt record (PS 9 'Interrupted runs')"))
        log.open(os.path.join(out_dir, run_stem + ".attempt%d.log" % k))
        log("attempt %d of %s (no limit on attempts; previous attempts' partial files kept: %s)"
            % (k, run_stem, ", ".join(kept) or "none"))
        rec["attempt"] = k
        install_stop_handlers()
        try:
            rec["sha256_covariates_file"] = ensure_covariates_file(out_dir, cohort, run_stem,
                                                                   payload, log)
            rec["covariates"] = dict(markers=cov["markers"], rna_counts=cov["counts"],
                                     r20_explained_variance_ratio=cov["evr"],
                                     c3_testable=cov["C3"] is not None, c3_closed=c3_closed)
            outputs = run_production(args, log, cohort, prep, cov, pcs, perms, saved_genes,
                                     finals, workers, cv_fn, folds)
        except SchedulerStop as e:
            log(scheduler_stop_text(e.args[0] if e.args else "?", k, run_stem))
            return 143
        except SystemExit:
            # fatal() inside an attempt is an exception of the run too (PS 9 "Defects")
            log(exception_text(k, run_stem))
            raise
        except Exception:  # noqa: BLE001
            tb = traceback.format_exc()
            log(tb)
            log(exception_text(k, run_stem))
            try:
                write_json_exclusive(os.path.join(out_dir, run_stem + ".attempt%d.exception.json"
                                                  % k), finish_record(dict(rec, traceback=tb), t0))
            except Exception:  # noqa: BLE001
                pass
            raise
        rec.update(outputs=outputs, n_usable=len(saved_genes), workers=workers)
        write_json_exclusive(os.path.join(out_dir, run_stem + ".run%d.json" % k),
                             finish_record(rec, t0))
        log("DONE: %s completed (attempt %d), %d file(s) written, %.1f min"
            % (run_stem, k, len(outputs), rec["minutes"]))
        return 0
    finally:
        lock.release()


def run_production(args, log, cohort, prep, cov, pcs, perms, saved_genes, finals, workers,
                   cv_fn, folds):
    n = len(prep["common"])
    E0 = np.empty((n, 0), dtype=np.float64)
    blocks = dict(pcs=pcs, C3=cov["C3"], R20=cov["R20"])
    outputs = {}
    if args.mode == "arm":
        e_name, b_name = ARMS[args.arm]
        E = blocks[e_name] if e_name else E0
        B = blocks[b_name]
        log("\n=== arm %s: E=%s (%d cols), B=%s (%d cols), %d genes x (1 + %d) on %d workers ==="
            % (args.arm, e_name or "none", E.shape[1], b_name, B.shape[1], len(prep["tested"]),
               N_PERM, workers))
        t1 = time.time()
        try:
            genes_u, M = compute_matrix(prep["tested"], prep["R"], prep["P"], E, B, perms,
                                        prep["disc_incr"], saved_genes, workers, cv_fn, folds)
        except ValueError as e:
            fatal(str(e))
        outputs[os.path.basename(finals[0])] = save_matrix_npz(finals[0], genes_u, M,
                                                               prep["common"])
        log("  wrote %s (%d genes x %d) in %.1f min" % (os.path.basename(finals[0]), M.shape[0],
                                                         M.shape[1], (time.time() - t1) / 60))
        return outputs
    a, b = parse_draws(args.draws)
    for j, final in zip(range(a, b + 1), finals):
        if os.path.exists(final):
            log("  draw %03d: %s exists, skipped" % (j, os.path.basename(final)))
            continue
        q = cap_perm(cohort, j, n)
        E = placebo_from(cov["R20"], q)
        t1 = time.time()
        try:
            genes_u, M = compute_matrix(prep["tested"], prep["R"], prep["P"], E, pcs,
                                        perms[:N_CAP_PERMS], prep["disc_incr"], saved_genes,
                                        workers, cv_fn, folds)
        except ValueError as e:
            fatal("draw %03d: %s" % (j, e))
        outputs[os.path.basename(final)] = save_matrix_npz(final, genes_u, M, prep["common"])
        log("  draw %03d: wrote %s in %.1f min (seed crc32 %d)"
            % (j, os.path.basename(final), (time.time() - t1) / 60, cap_seed(cohort, j)))
    return outputs


def mode_gate(args, log, prep, pcs, perms, saved_genes, saved_M, workers, cv_fn, folds, rec,
              out_dir, tag, t0):
    cohort = args.cohort
    c = COHORTS[cohort]
    n = len(prep["common"])
    tested = prep["tested"]
    log("\n=== H-V1: new per-gene function, E empty, B = published components, columns 0..%d ==="
        % N_GATE_PERMS)
    E0 = np.empty((n, 0), dtype=np.float64)
    res = run_pool_cov(len(tested), workers, prep["R"], prep["P"], E0, pcs,
                       perms[:N_GATE_PERMS], MIN_N, cv_fn, folds)
    usable = usable_indices(tested, res, prep["disc_incr"])
    genes_u, M = matrix_from(tested, res, usable)
    same_genes = list(genes_u) == list(saved_genes)
    # published quantities only: the AUC of the published selection set (S_M) from column 0
    auc_new, n1, n0 = hv1_auc(prep["assemble"], M, genes_u, prep["selected"], prep["disc_incr"])
    auc_saved, _s1, _s0 = hv1_auc(prep["assemble"], saved_M, saved_genes, prep["selected"],
                                  prep["disc_incr"])
    v = hv1_verdicts(same_genes, M, saved_M, auc_new, c["published_auc"])
    rec.update(gate="H-V1", gate_host=host(), tolerance_strict=HV1_TOL,
               tolerance_fallback_deviation=HV1_FALLBACK_DEV,
               tolerance_fallback_auc=HV1_FALLBACK_AUC, identical_gene_list=same_genes,
               n_usable=len(genes_u), n_selected_usable=n1, n_background_usable=n0,
               auc_recomputed_col0=auc_new, auc_published=c["published_auc"],
               auc_saved_col0=auc_saved, workers=workers, **v)
    rec["production_eligible"] = bool(v["passed_strict"] or v["passed_fallback"])
    log("  host                    : " + host())
    log("  usable genes            : %d (saved %d); identical list and order: %s"
        % (len(genes_u), len(saved_genes), same_genes))
    log("  largest |diff| cols 0..5: %.3e   99.9th percentile: %.3e   by column: %s"
        % (v["max_abs_diff"], v["p999_abs_diff"],
           " ".join("%.2e" % x for x in v["max_abs_diff_by_column"])))
    log("  AUC of S_M (%d / %d)    : recomputed col 0 %r; published %r; |diff| %.3e; saved col 0 %r"
        % (n1, n0, auc_new, c["published_auc"], v["auc_abs_dev"], auc_saved))
    log("  strict  (genes identical, |diff| <= 1e-9)                  : %s"
        % ("PASSED" if v["passed_strict"] else "FAILED"))
    log("  fallback (|diff| < 1e-6 and |AUC - published| <= 1e-6)      : %s"
        % ("PASSED" if v["passed_fallback"] else "FAILED"))
    log("H-V1 PASSED (strict)" if v["passed_strict"] else
        "H-V1 STRICT FAILED (fallback %s; try the next host of the PS 9 order)"
        % ("passed" if v["passed_fallback"] else "failed"))
    write_json_exclusive(os.path.join(out_dir, "c1_cov_%s_gate_hv1_%s_%s.json"
                                      % (cohort, host(), tag.split("_")[-1])),
                         finish_record(rec, t0))
    return 0 if v["passed_strict"] else 3


def mode_blind(args, log, RA, CC, prep, saved_genes, rec, out_dir, tag, t0):
    cohort = args.cohort
    tested = prep["tested"]
    usable = usable_by_count(tested, prep["P"], prep["disc_incr"])
    genes_u = [tested[j] for j in usable]
    rna = RA.load_rna_matrix()
    rc = rna.loc[prep["common"]]
    try:
        C3, minfo = build_c3(rc, CC.comp_score, CC.MARKERS)
    except ValueError as e:
        fatal(str(e))
    _R20, evr = r20_recipe(rc)
    emb = embedding_counts(RA, prep["common"])
    roots = args.list_roots.split(",") if args.list_roots else DEFAULT_LIST_ROOTS
    rows, truncated = list_tree(roots)
    lst = os.path.join(out_dir, "c1_cov_%s_blind_listing_%s_%s.tsv"
                       % (cohort, host(), tag.split("_")[-1]))
    write_text_exclusive(lst, "path\tsize\tmtime_utc\tname_flag\n" +
                         "".join("\t".join(r) + "\n" for r in rows))
    counts = dict(n_cases=len(prep["common"]), K=prep["K"], genes_tested=len(tested),
                  genes_usable_by_count=len(genes_u),
                  usable_identical_to_saved_matrix=list(genes_u) == list(saved_genes),
                  n_selected_usable=len(set(genes_u) & prep["selected"]),
                  rna=rna_column_counts(rc),
                  r20_explained_variance_ratio=[float(x) for x in evr],
                  r20_explained_variance_total=float(evr.sum()),
                  markers=minfo, c3_testable=C3 is not None, embeddings=emb,
                  listing=dict(roots=roots, rows=len(rows), truncated=truncated,
                               flagged=sum(1 for r in rows if r[3] == "flag"),
                               file=os.path.basename(lst)))
    rec["counts"] = counts
    log("\nBLIND (counts only; no increment of any kind):")
    log("  cases %d, K %d, genes tested %d, usable %d (identical to saved matrix: %s), "
        "published set genes among them %d"
        % (counts["n_cases"], counts["K"], counts["genes_tested"],
           counts["genes_usable_by_count"], counts["usable_identical_to_saved_matrix"],
           counts["n_selected_usable"]))
    log("  RNA columns %d, all-zero %d, SD < 1e-8 %d" % (counts["rna"]["rna_columns"],
                                                       counts["rna"]["all_zero_columns"],
                                                       counts["rna"]["sd_below_1e-8_columns"]))
    log("  R20 explained variance: total %.4f; by component %s"
        % (counts["r20_explained_variance_total"],
           " ".join("%.4f" % x for x in counts["r20_explained_variance_ratio"])))
    log_markers(log, minfo)
    log("  C3 testable (PS 2: >= %d markers present and >= 1 with non-zero SD, every family): %s"
        % (MIN_MARKERS, counts["c3_testable"]))
    log("  embeddings: %s" % json.dumps(emb))
    log("  result-folder listing: %d rows (%d name-flagged), truncated=%s -> %s"
        % (len(rows), counts["listing"]["flagged"], truncated, os.path.basename(lst)))
    write_json_exclusive(os.path.join(out_dir, "c1_cov_%s_blind_%s_%s.json"
                                      % (cohort, host(), tag.split("_")[-1])),
                         finish_record(rec, t0))
    log("BLIND DONE")
    return 0


def mode_discovery(args, log):
    t0 = time.time()
    cohort = args.cohort
    env = env_guard(cohort, "discovery")
    out_dir = env["MORPHO_OUT"]
    os.makedirs(out_dir, exist_ok=True)
    tag = "%s_%s_%s" % (args.mode.replace("-", ""), host(), stamp())
    log.open(os.path.join(out_dir, "c1_cov_%s_%s.log" % (cohort, tag)))
    rec = base_record(args, args.mode)
    rec.update(env=env, gate_host=host())
    log("=" * 78)
    log("C1 COVARIATES (POST HOC) DISCOVERY TREE -- cohort=%s mode=%s host=%s start=%s"
        % (cohort, args.mode, host(), rec["start_utc"]))
    log("  this script sha256 = " + rec["sha256_this_script"])
    log("=" * 78)
    if not args.patients or not args.inputs_npz:
        fatal("--patients and --inputs-npz are required")
    z, npz_patients = load_inputs_patients(args.inputs_npz, cohort)
    patients = read_patient_file(args.patients)
    if patients != npz_patients:
        fatal("--patients differs from the 'patients' array of the inputs file")
    if len(patients) != COHORTS[cohort]["n_disc"]:
        fatal("%d patients, PS 2 expects %d" % (len(patients), COHORTS[cohort]["n_disc"]))
    rec["sha256_inputs"] = dict(inputs_npz=sha256(args.inputs_npz), patients=sha256(args.patients))
    mods, code_hashes, _origins = import_frozen(env, log)
    rec["sha256_frozen_modules"] = code_hashes
    RA, CC = mods["residual_analysis"], mods["composition_controls"]
    rec["sha256_inputs"]["rna_manifest"] = sha256(env["MORPHO_RNA_MANIFEST"])
    rna = RA.load_rna_matrix()
    miss = [p for p in patients if p not in rna.index]
    if miss:
        fatal("%d patients have no RNA in this tree: %s" % (len(miss), miss[:5]))
    rc = rna.loc[patients]

    if args.mode == "gate-r20":
        R20, evr = r20_recipe(rc)
        ref = np.asarray(z["rna_pcs20"], dtype=np.float64)
        dmax, signs = max_diff_sign_aligned(R20, ref)
        passed = bool(dmax <= HV2_TOL)
        rec.update(gate="H-V2", passed=passed, tolerance=HV2_TOL, max_abs_diff_sign_aligned=dmax,
                   column_signs=signs, rna=rna_column_counts(rc),
                   r20_explained_variance_ratio=evr, n_patients=len(patients))
        log("  host %s; R20 %s vs rna_pcs20 %s; max |diff| after sign alignment %.3e; signs %s"
            % (host(), R20.shape, ref.shape, dmax, signs))
        log("H-V2 PASSED" if passed else "H-V2 FAILED")
        write_json_exclusive(os.path.join(out_dir, "c1_cov_%s_gate_r20_%s_%s.json"
                                          % (cohort, host(), tag.split("_")[-1])),
                             finish_record(rec, t0))
        return 0 if passed else 3

    # export-scores (H0)
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    if common != patients:
        fatal("sorted(common) of the discovery tree (%d) differs from the inputs file's "
              "patients (%d)" % (len(common), len(patients)))
    rec["sha256_inputs"].update(protein=sha256(os.environ.get("MORPHO_PROTEIN_TSV", "") or os.devnull),
                                slide_map=sha256(env["MORPHO_SLIDE_MAP"]))
    try:
        C3, minfo = build_c3(rc, CC.comp_score, CC.MARKERS)
    except ValueError as e:
        fatal(str(e))
    log_markers(log, minfo)
    rec.update(markers=minfo, n_patients=len(patients), rna=rna_column_counts(rc))
    final = os.path.join(out_dir, "c1_cov_discovery_scores_%s.tsv" % cohort)
    if C3 is None:
        rec.update(untestable=True)
        write_json_exclusive(os.path.join(out_dir, "c1_cov_%s_export_scores_untestable_%s.json"
                                          % (cohort, tag.split("_")[-1])), finish_record(rec, t0))
        raise Refused("C3 untestable in discovery %s (marker rule of PS 2); no score table "
                      "written" % cohort)
    comments = ["c1_covariates_cluster.py export-scores (H0), cohort %s, %d patients" % (cohort, len(patients)),
                "host %s; lsf_job %s; written %s" % (host(), rec["lsf_job"], utc()),
                "script_sha256 %s" % rec["sha256_this_script"],
                "composition_controls_sha256 %s" % code_hashes["composition_controls"]]
    for fam in FAMILIES:
        comments.append("markers %s %d/%d present (zero-SD %d): %s"
                        % (fam, minfo[fam]["n_present"], minfo[fam]["n_listed"],
                           minfo[fam]["n_zero_sd"], ",".join(minfo[fam]["used"])))
    h = write_text_exclusive(final, scores_tsv_text(patients, C3, comments))
    rec.update(outputs={os.path.basename(final): h})
    write_json_exclusive(os.path.join(out_dir, "c1_cov_%s_export_scores_%s.json"
                                      % (cohort, tag.split("_")[-1])), finish_record(rec, t0))
    log("  wrote %s (sha256 %s)" % (final, h))
    log("EXPORT DONE")
    return 0


def mode_manifest(args, log):
    if args.cohort != "ucec":
        fatal("the manifest is written once, through the UCEC wrapper (it covers both trees)")
    env = env_guard(args.cohort, "confirmatory")
    if os.path.normpath(env["MORPHO_OUT"]) != os.path.normpath(MANIFEST_DIRS["ucec"]):
        fatal("MORPHO_OUT is not %s" % MANIFEST_DIRS["ucec"])
    resolve_paths(args)
    closed = read_closed_parts(args.closed_parts)
    log("C1 COVARIATES manifest -- host=%s start=%s" % (host(), utc()))
    _mods, _hashes, origins = import_frozen(env, log)
    verify_sidecar(args, args.cohort, origins, log)
    log("  CLOSED_PARTS: %s" % json.dumps(closed_summary(closed)))
    text, entries = build_manifest_text(MANIFEST_DIRS, closed)
    final = os.path.join(env["MORPHO_OUT"], MANIFEST_NAME)
    h = write_manifest(final, text, args.force_rewrite)
    log("manifest: %d files -> %s" % (len(entries), final))
    log("manifest sha256 = " + h)
    log("(record this value: Stage L1 and Stage L2 are given it on the command line, PS 9)")
    return 0


def mode_export_patients(args, log):
    if not args.inputs_npz or not args.out:
        fatal("export-patients needs --inputs-npz and --out")
    _z, patients = load_inputs_patients(args.inputs_npz, args.cohort)
    if patients != sorted(patients) or len(set(patients)) != len(patients):
        fatal("inputs patients are not unique and sorted")
    h = write_text_exclusive(args.out, "".join(p + "\n" for p in patients))
    log("wrote %d patient ids to %s (sha256 %s)" % (len(patients), args.out, h))
    return 0


# ---------------------------------------------------------------------------------------
# Self-test (synthetic data only; no project import)
# ---------------------------------------------------------------------------------------
def _cv_r2_copy(X, y, folds):
    """Verbatim copy of residual_analysis.cv_r2 (self-test only; production uses RA.cv_r2)."""
    from sklearn.model_selection import KFold
    kf = KFold(n_splits=folds, shuffle=True, random_state=0)
    preds = np.zeros_like(y, dtype=float)
    ridge_I = 1.0 * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        X_tr, X_te, y_tr = X[tr], X[te], y[tr]
        x_mean = X_tr.mean(axis=0)
        x_std = X_tr.std(axis=0)
        x_std[x_std == 0] = 1.0
        Xtr_s = (X_tr - x_mean) / x_std
        Xte_s = (X_te - x_mean) / x_std
        y_mean = y_tr.mean()
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + ridge_I, Xtr_s.T @ (y_tr - y_mean))
        preds[te] = Xte_s @ w + y_mean
    ss_res = np.sum((y - preds) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def _comp_score_copy(rna, markers):
    """Verbatim copy of composition_controls.comp_score (self-test only)."""
    mk = [m for m in markers if m in rna.columns]
    if len(mk) < 6:
        return None, mk
    z = rna[mk].sub(rna[mk].mean()).div(rna[mk].std().replace(0, np.nan))
    return z.mean(axis=1).values.reshape(-1, 1), mk


def _assemble_copy(M, xs, sel_mask):
    """Verbatim copy of c1_run_test.assemble (self-test only)."""
    from scipy.stats import rankdata
    n1, n0 = int(sel_mask.sum()), int((~sel_mask).sum())
    ranks = rankdata(M, axis=0)
    auc = (ranks[sel_mask].sum(axis=0) - n1 * (n1 + 1) / 2.0) / (n1 * n0)
    rs = rankdata(xs)
    rs = (rs - rs.mean()) / rs.std()
    rc = (ranks - ranks.mean(axis=0)) / ranks.std(axis=0)
    rho = (rs[:, None] * rc).mean(axis=0)
    return auc, rho, n1, n0


def _ref_r2_sklearn(X, y):
    """Independent implementation: sklearn KFold + StandardScaler + Ridge(alpha=1)."""
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold
    from sklearn.preprocessing import StandardScaler
    pred = np.zeros(len(y))
    for tr, te in KFold(n_splits=5, shuffle=True, random_state=0).split(X):
        sc = StandardScaler().fit(X[tr])
        m = Ridge(alpha=1.0, fit_intercept=True, solver="cholesky").fit(sc.transform(X[tr]), y[tr])
        pred[te] = m.predict(sc.transform(X[te]))
    return 1.0 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)


def _ref_gene(R, P, E, B, perms, j, r2fn):
    """Straightforward per-gene implementation of PS 4 (no shared code with the worker)."""
    y_all = P[:, j]
    keep = [i for i in range(len(y_all)) if not np.isnan(y_all[i])]
    if len(keep) < MIN_N:
        return None
    keep = np.array(keep)
    y = y_all[keep]
    base = np.column_stack([R[keep, j]] + [E[keep, c] for c in range(E.shape[1])])
    r2b = r2fn(base, y)
    vals = []
    for p in [np.arange(len(y_all))] + list(perms):
        Bp = np.array([B[p[i]] for i in range(len(y_all))])
        X = np.column_stack([base] + [Bp[keep, c] for c in range(B.shape[1])])
        vals.append(r2fn(X, y) - r2b)
    return np.array(vals)


def _published_gene_test_copy(R, P, pcs, perms, j):
    """Structure of split_replication_power._gene_test (self-test only)."""
    y = P[:, j]
    v = ~np.isnan(y)
    if v.sum() < MIN_N:
        return None
    yv = y[v]
    xr = R[v, j].reshape(-1, 1)
    r2r = _cv_r2_copy(xr, yv, 5)
    out = np.empty(len(perms) + 1)
    out[0] = _cv_r2_copy(np.hstack([xr, pcs[v]]), yv, 5) - r2r
    for b, p in enumerate(perms, 1):
        out[b] = _cv_r2_copy(np.hstack([xr, pcs[p][v]]), yv, 5) - r2r
    return out


def _try_theta_kernel():
    for d in (os.path.join(_SCRIPT_DIR, "..", "MEAS_ERROR_2026-10-08"), _SCRIPT_DIR):
        p = os.path.join(d, "theta_kernel.py")
        if os.path.isfile(p):
            if sha256(p) != THETA_KERNEL_SHA256:
                return None, "theta_kernel.py found but its sha256 differs from PS 1"
            import importlib.util
            spec = importlib.util.spec_from_file_location("theta_kernel_c1cov", p)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod, p
    return None, "theta_kernel.py not found"


def _raises(exc, fn, *a, **kw):
    import contextlib
    import io
    try:
        with contextlib.redirect_stderr(io.StringIO()):    # argparse usage text stays out of the log
            fn(*a, **kw)
    except exc:
        return True
    return False


def selftest(args, log):
    log("sha256 %s c1_covariates_cluster.py" % sha256(os.path.abspath(__file__)))
    import pandas as pd
    import shutil
    from sklearn.decomposition import PCA
    t0 = time.time()
    rs = np.random.RandomState(20261121)
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond)))
        log("  [%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))

    # ---- synthetic cohort ----
    n, G, K = 60, 14, 6
    R = rs.randn(n, G)
    P = R * 0.5 + rs.randn(n, G)
    P[rs.rand(n, G) < 0.12] = np.nan
    P[:35, 3] = np.nan                       # gene 3: < 30 protein values -> skipped
    P[:, 7] = np.where(np.arange(n) < 5, np.nan, P[:, 7])
    pcs = rs.randn(n, K)
    C3 = rs.randn(n, 3)
    R20s = rs.randn(n, 20)
    perms = published_perms(b"c1_selftest|perm", n, 4)
    E0 = np.empty((n, 0))

    log("-- 1. per-gene function (PS 4)")
    for label, E, B in (("A-type (E=C3, B=pcs)", C3, pcs), ("A-type (E=R20, B=pcs)", R20s, pcs),
                        ("V-type (E=none, B=C3)", E0, C3), ("V-type (E=none, B=R20)", E0, R20s)):
        res = run_pool_cov(G, 0, R, P, E, B, perms, MIN_N, _cv_r2_copy, 5)
        worst, none_ok = 0.0, True
        for j in range(G):
            ref = _ref_gene(R, P, E, B, perms, j, _ref_r2_sklearn)
            if (ref is None) != (res[j] is None):
                none_ok = False
                continue
            if ref is not None:
                worst = max(worst, float(np.max(np.abs(ref - res[j]))))
        ok("per-gene function == independent sklearn implementation, " + label,
           none_ok and worst <= 1e-9, "(max |diff| %.2e)" % worst)
    ok("gene with < 30 protein values is skipped", res[3] is None)
    res = run_pool_cov(G, 0, R, P, E0, pcs, perms, MIN_N, _cv_r2_copy, 5)
    exact = all((res[j] is None and _published_gene_test_copy(R, P, pcs, perms, j) is None) or
                np.array_equal(res[j], _published_gene_test_copy(R, P, pcs, perms, j))
                for j in range(G))
    ok("E empty, B = pcs: bit-identical to the published _gene_test structure", exact)
    res_pool = run_pool_cov(G, 2, R, P, C3, pcs, perms, MIN_N, _cv_r2_copy, 5)
    res_ser = run_pool_cov(G, 0, R, P, C3, pcs, perms, MIN_N, _cv_r2_copy, 5)
    ok("multiprocessing pool == serial (bit-identical)",
       all((res_pool[j] is None and res_ser[j] is None) or np.array_equal(res_pool[j], res_ser[j])
           for j in range(G)))
    tk, where = _try_theta_kernel()
    if tk is None:
        log("  [SKIP] theta_kernel comparison: " + where)
    else:
        worst = 0.0
        resA = run_pool_cov(G, 0, R, P, C3, pcs, perms, MIN_N, _cv_r2_copy, 5)
        for j in range(G):
            if resA[j] is None:
                continue
            v = ~np.isnan(P[:, j])
            y = P[v, j]
            f = tk.folds_for(int(v.sum()))
            base = np.hstack([R[v, j].reshape(-1, 1), C3[v]])
            r2b = tk.r2_of(y, tk.oof_ridge(base, y, f))
            for b, p in enumerate([np.arange(n)] + perms):
                inc = tk.r2_of(y, tk.oof_ridge(np.hstack([base, pcs[p][v]]), y, f)) - r2b
                worst = max(worst, abs(inc - resA[j][b]))
        ok("per-gene function == tk.oof_ridge increments (" + where + ")", worst <= 1e-12,
           "(max |diff| %.2e)" % worst)

    log("-- 2. usable filter, matrix assembly, compute_matrix")
    tested = ["G%02d" % j for j in range(G)]
    disc = {g: (np.nan if g == "G05" else 0.1 * (j + 1)) for j, g in enumerate(tested)}
    use = usable_indices(tested, res, disc)
    ok("usable filter = (>= 30 protein values) and finite discovery increment",
       use == [j for j in range(G) if j not in (3, 5)])
    ok("count-based usable filter (blind) == increment-based filter",
       usable_by_count(tested, P, disc) == use)
    genes_u, M = matrix_from(tested, res, use)
    saved_genes = [tested[j] for j in use]
    good = not _raises(ValueError, check_matrix, genes_u, M, saved_genes, len(perms) + 1)
    bad = _raises(ValueError, check_matrix, genes_u[::-1], M, saved_genes, len(perms) + 1)
    ok("matrix checks accept the right gene order and reject a different one", good and bad)
    big_perms = published_perms(b"c1_selftest|perm", n, 40)
    gA, MA = compute_matrix(tested, R, P, C3, pcs, big_perms, disc, saved_genes, 0, _cv_r2_copy, 5)
    gC, MC = compute_matrix(tested, R, P, R20s[cap_perm("ucec", 7, n)], pcs,
                            big_perms[:N_CAP_PERMS], disc, saved_genes, 0, _cv_r2_copy, 5)
    ok("arm matrix: G x (1 + B); placebo matrix: G x 26 = columns 0 and 1-25 of the same perms",
       MA.shape == (len(saved_genes), 41) and MC.shape == (len(saved_genes), N_CAP_PERMS + 1)
       and N_CAP_PERMS + 1 == 26)
    ok("compute_matrix refuses a gene list that differs from the saved matrix's",
       _raises(ValueError, compute_matrix, tested, R, P, C3, pcs, perms, disc,
               saved_genes[:-1], 0, _cv_r2_copy, 5))

    log("-- 3. R20 recipe and H-V2 comparison")
    rna = pd.DataFrame(np.abs(rs.randn(40, 300)) * 3, index=["C3L-%05d" % i for i in range(40)],
                       columns=["g%d" % i for i in range(300)])
    rna.iloc[:, 5] = 0.0
    rna.iloc[:, 9] = 2.5                                   # constant non-zero column
    R20, evr = r20_recipe(rna)
    rv = rna.values.astype(float)
    sd = rv.std(0)
    rz = (rv - rv.mean(0)) / np.where(sd == 0, 1.0, sd)
    U, S, Vt = np.linalg.svd(rz - rz.mean(0), full_matrices=False)
    ref = U[:, :20] * S[:20]
    d, _ = max_diff_sign_aligned(R20, ref)
    ok("R20 recipe == SVD of the standardised matrix (up to column sign)", d <= 1e-8,
       "(max |diff| %.2e)" % d)
    ok("R20 explained variance == SVD", np.allclose(evr, (S ** 2 / np.sum(S ** 2))[:20], atol=1e-12))
    ok("zero-SD columns enter as zeros", np.all(rz[:, 5] == 0) and np.all(rz[:, 9] == 0))
    flipped = ref * np.array([(-1) ** k for k in range(20)])
    d2, signs = max_diff_sign_aligned(flipped, ref)
    ok("sign-aligned comparison handles flipped columns", d2 == 0.0 and signs[1] == -1)
    cc = rna_column_counts(rna)
    ok("RNA column counts (all-zero, SD < 1e-8)",
       cc["all_zero_columns"] == 1 and cc["sd_below_1e-8_columns"] == 2)

    log("-- 4. C3 rule (PS 2)")
    mk = {f: ["%s%d" % (f[:3].upper(), i) for i in range(8)] for f in FAMILIES}
    cols = sum(mk.values(), [])
    rr = pd.DataFrame(rs.randn(30, len(cols)), columns=cols)
    rr["STR0"] = 1.0                                        # zero-SD marker -> skipped
    C, info = build_c3(rr, _comp_score_copy, mk)
    z = rr[mk["stroma"][1:]]
    man = ((z - z.mean()) / z.std(ddof=1)).mean(axis=1).values
    ok("C3: three columns in order; zero-SD marker skipped; ddof 1",
       C is not None and C.shape == (30, 3) and np.allclose(C[:, 0], man, atol=1e-14)
       and info["stroma"]["n_zero_sd"] == 1 and "STR0" not in info["stroma"]["used"])
    ok("C3: the SD of every present marker is reported (pandas, ddof 1)",
       info["stroma"]["sd"]["STR0"] == 0.0 and
       abs(info["immune"]["sd"]["IMM2"] - float(rr["IMM2"].std(ddof=1))) < 1e-15
       and len(info["proliferation"]["sd"]) == 8)
    rr2 = rr.drop(columns=mk["immune"][:3])                 # 5 immune markers left
    C2, info2 = build_c3(rr2, _comp_score_copy, mk)
    ok("C3 untestable with fewer than 6 markers present in a family",
       C2 is None and info2["immune"]["n_present"] == 5 and not info2["immune"]["family_testable"])
    rr3 = rr.copy()
    for m in mk["proliferation"]:
        rr3[m] = 3.0
    C3b, info3 = build_c3(rr3, _comp_score_copy, mk)
    ok("C3 untestable when no present marker of a family has a non-zero SD (not a crash)",
       C3b is None and info3["proliferation"]["n_present"] == 8
       and info3["proliferation"]["used"] == [] and info3["stroma"]["family_testable"])
    rr4 = rr.copy()
    for m in mk["immune"][:7]:
        rr4[m] = 0.5
    C4, info4 = build_c3(rr4, _comp_score_copy, mk)
    ok("one marker with non-zero SD is enough (6+ present)",
       C4 is not None and info4["immune"]["used"] == ["IMM7"]
       and np.allclose(C4[:, 1], ((rr["IMM7"] - rr["IMM7"].mean()) / rr["IMM7"].std()).values))

    log("-- 5. placebo convention (PS 2, PS 4)")
    seeds_ok = all(cap_seed("ucec", j) == zlib.crc32(b"c1_cov|ucec|cap|%03d" % j) and
                   cap_seed("luad", j) == zlib.crc32(b"c1_cov|luad|cap|%03d" % j)
                   for j in range(1, N_CAP + 1))
    ok("placebo seeds = crc32(b'c1_cov|<cohort>|cap|%03d' % j), j = 1..119",
       seeds_ok and N_CAP == 119)
    ok("seed bytes for j = 7 and j = 119",
       (b"c1_cov|%s|cap|%03d" % (b"ucec", 7)) == b"c1_cov|ucec|cap|007"
       and cap_seed("luad", 119) == zlib.crc32(b"c1_cov|luad|cap|119"))
    for coh, nn in (("ucec", 138), ("luad", 112)):
        qs = [cap_perm(coh, j, nn) for j in range(1, N_CAP + 1)]
        ok("119 distinct permutations of 0..%d (%s)" % (nn - 1, coh),
           all(sorted(q.tolist()) == list(range(nn)) for q in qs) and
           len({tuple(q.tolist()) for q in qs}) == N_CAP)
    ok("placebo file names: 3-digit draw index",
       cap_name("ucec", 1) == "c1_cov_ucec_A_R20cap_001.npz"
       and cap_name("luad", 119) == "c1_cov_luad_A_R20cap_119.npz")
    qs = [cap_perm("ucec", j, 138) for j in (1, 2)]
    Rbig = rs.randn(138, 20)
    pl = placebo_from(Rbig, qs[0])
    ok("row i of the placebo is row q[i] of R20",
       all(np.array_equal(pl[i], Rbig[qs[0][i]]) for i in range(138)))
    ok("placebo keeps the column covariance of R20",
       np.allclose(np.cov(pl.T), np.cov(Rbig.T), atol=1e-12))
    p1 = published_perms(b"c1_ucec|perm", 138, N_PERM)
    p2 = published_perms(b"c1_ucec|perm", 138, N_PERM)[:N_CAP_PERMS]
    ok("published permutations: first 25 = columns 1-25 of the 1,000",
       N_CAP_PERMS == 25 and all(np.array_equal(a, b) for a, b in zip(p1[:N_CAP_PERMS], p2)))
    rng = np.random.RandomState(zlib.crc32(b"c1_luad|perm"))
    ok("published permutation recipe (luad)",
       np.array_equal(published_perms(b"c1_luad|perm", 112, 3)[2],
                      [rng.permutation(112) for _ in range(3)][2]))
    ok("--draws parsing: 1-30, 91-119, 7-7 accepted",
       parse_draws("1-30") == (1, 30) and parse_draws("91-119") == (91, 119)
       and parse_draws("7-7") == (7, 7))
    ok("--draws parsing: 0-5, 100-120, 30-1, 1-1000, 'a-b' refused",
       all(_raises(SystemExit, parse_draws, s) for s in ("0-5", "100-120", "30-1", "1-1000", "a-b")))
    ok("placebo run stem: c1_cov_<cohort>_A_R20cap_<aaa>-<bbb>",
       placebo_stem("ucec", 1, 30) == "c1_cov_ucec_A_R20cap_001-030")

    log("-- 6. gate H-V1 verdicts (PS 9)")
    gsyn = ["S%05d" % i for i in range(400)]
    Msyn = rs.randn(400, 6) * 0.05
    sel = set(gsyn[:120])
    dsyn = {g: float(x) for g, x in zip(gsyn, rs.randn(400))}
    a_pub, _n1, _n0 = hv1_auc(_assemble_copy, Msyn, gsyn, sel, dsyn)
    sel_mask = np.array([g in sel for g in gsyn])
    from scipy.stats import mannwhitneyu
    u = mannwhitneyu(Msyn[sel_mask, 0], Msyn[~sel_mask, 0], alternative="two-sided").statistic
    ok("hv1_auc == Mann-Whitney U / (n1 n0) on column 0", abs(a_pub - u / (120 * 280)) < 1e-15
       and (_n1, _n0) == (120, 280))
    v = hv1_verdicts(True, Msyn.copy(), Msyn, a_pub, a_pub)
    ok("identical matrix: strict and fallback pass", v["passed_strict"] and v["passed_fallback"]
       and v["max_abs_diff"] == 0.0)
    Mn = Msyn + 3e-10 * rs.choice([-1, 1], size=Msyn.shape)
    v = hv1_verdicts(True, Mn, Msyn, hv1_auc(_assemble_copy, Mn, gsyn, sel, dsyn)[0], a_pub)
    ok("deviation 3e-10: strict passes", v["passed_strict"] and v["passed_fallback"])
    Mf = Msyn.copy()
    Mf[5, 2] += 4e-7
    Mf[11, 0] += 2e-8
    a_f = hv1_auc(_assemble_copy, Mf, gsyn, sel, dsyn)[0]
    v = hv1_verdicts(True, Mf, Msyn, a_f, a_pub)
    ok("deviation 4e-7, AUC unchanged: strict fails, fallback passes",
       not v["passed_strict"] and v["passed_fallback"] and abs(v["max_abs_diff"] - 4e-7) < 1e-12)
    v = hv1_verdicts(True, Mf, Msyn, a_pub + 2e-6, a_pub)
    ok("deviation 4e-7 but |AUC - published| = 2e-6: neither passes",
       not v["passed_strict"] and not v["passed_fallback"])
    v = hv1_verdicts(True, Mf, Msyn, a_pub + 9e-7, a_pub)
    v2 = hv1_verdicts(True, Mf, Msyn, a_pub - 1.1e-6, a_pub)
    ok("|AUC - published| = 9e-7 is within the fallback, 1.1e-6 is not",
       v["passed_fallback"] and not v2["passed_fallback"])
    Mg = Msyn.copy()
    Mg[3, 4] += 2e-6
    v = hv1_verdicts(True, Mg, Msyn, a_pub, a_pub)
    ok("deviation 2e-6 (not below 1e-6), AUC equal: neither passes",
       not v["passed_strict"] and not v["passed_fallback"])
    v = hv1_verdicts(False, Msyn.copy(), Msyn, a_pub, a_pub)
    ok("different gene list: neither passes, deviation infinite",
       not v["passed_strict"] and not v["passed_fallback"] and v["max_abs_diff"] == float("inf"))
    Mx = Msyn.copy()
    Mx[0, 0] = np.nan
    v = hv1_verdicts(True, Mx, Msyn, a_pub, a_pub)
    ok("a non-finite recomputed value: neither passes", not v["passed_strict"] and not v["passed_fallback"])

    tmp = tempfile.mkdtemp(prefix="c1cov_selftest_")
    try:
        log("-- 7. production acceptance of H-V1 records (same host, script, inputs)")
        gdir = os.path.join(tmp, "gates")
        os.makedirs(gdir)
        expect = dict(sha256_this_script="a" * 64, sha256_saved_matrix="b" * 64,
                      sha256_inputs=dict(protein="c" * 64, case_ids="d" * 64),
                      sha256_frozen_modules=dict(FROZEN_CODE_SHA256), cohort="ucec", host="s006")

        def gate_json(name, **over):
            r = dict(expect, gate="H-V1", passed_strict=False, passed_fallback=False)
            r.update(over)
            with open(os.path.join(gdir, name), "w", encoding="utf-8") as fh:
                json.dump(r, fh)

        st = "20261121T000000Z-1-000001"
        gate_json("c1_cov_ucec_gate_hv1_s006_%s.json" % st)
        ok("a failed gate json is not accepted", find_accepted_hv1(gdir, "ucec", "s006", expect)[0] is None)
        gate_json("c1_cov_ucec_gate_hv1_s006_%s.json" % st.replace("-1-", "-2-"),
                  passed_fallback=True, sha256_this_script="e" * 64)
        ok("a passed gate json of another script sha256 is not accepted",
           find_accepted_hv1(gdir, "ucec", "s006", expect)[0] is None)
        gate_json("c1_cov_ucec_gate_hv1_s006_%s.json" % st.replace("-1-", "-3-"),
                  passed_strict=True, sha256_inputs=dict(protein="f" * 64, case_ids="d" * 64))
        ok("a passed gate json on other inputs is not accepted",
           find_accepted_hv1(gdir, "ucec", "s006", expect)[0] is None)
        gate_json("c1_cov_ucec_gate_hv1_s002_%s.json" % st.replace("-1-", "-4-"),
                  passed_strict=True, host="s002")
        ok("a passed gate json of another host is not accepted",
           find_accepted_hv1(gdir, "ucec", "s006", expect)[0] is None)
        gate_json("c1_cov_luad_gate_hv1_s006_%s.json" % st.replace("-1-", "-5-"),
                  passed_strict=True, cohort="luad")
        ok("a passed gate json of the other cohort is not accepted",
           find_accepted_hv1(gdir, "ucec", "s006", expect)[0] is None)
        gate_json("c1_cov_ucec_gate_hv1_s006_%s.json" % st.replace("-1-", "-6-"),
                  passed_fallback=True)
        f_ok, r_ok = find_accepted_hv1(gdir, "ucec", "s006", expect)
        ok("a fallback-passed gate json of this host, script and inputs is accepted",
           f_ok is not None and f_ok.endswith("-6-000001.json") and r_ok["passed_fallback"])

        log("-- 8. CLOSED_PARTS.txt and PINNED_HOSTS.txt (PS 6.1, PS 9)")
        for nm_, rd_ in (("CLOSED_PARTS.txt", read_closed_parts),
                         ("PINNED_HOSTS.txt", read_pinned_hosts)):
            here = os.path.join(_SCRIPT_DIR, nm_)
            if os.path.isfile(here):
                val_ = rd_(here)
                ok("the project's %s parses" % nm_, isinstance(val_, dict),
                   "(%s)" % json.dumps(closed_summary(val_) if nm_.startswith("CLOSED")
                                       else val_))
        empty = parse_closed_parts_text("# header\n# only comments\n\n")
        ok("comments only: nothing closed", empty == dict(arms=set(), analyses=set()))
        cpx = parse_closed_parts_text(
            "# h\narm ucec A_C3\narm ucec V_C3\narm luad A_R20cap\narm luad A_R20\n"
            "analysis luad M|C3\nanalysis luad C3\n\n# end\n")
        ok("tokens parse: arms and analyses",
           cpx["arms"] == {("ucec", "A_C3"), ("ucec", "V_C3"), ("luad", "A_R20cap"),
                           ("luad", "A_R20")}
           and cpx["analyses"] == {("luad", "M|C3"), ("luad", "C3")})
        bads = ["arm ucec A_X\n", "arm brca A_C3\n", "arms ucec A_C3\n", "arm ucec  A_C3\n",
                "arm ucec A_C3 extra\n", "analysis luad M\n", "arm ucec A_C3\narm ucec A_C3\n",
                " arm ucec A_C3\n", "arm ucec A_C3 \n", "arm\tucec\tA_C3\n",
                "  # indented comment\n", " \n", "\t\n", "arm ucec A_C3\r\n",
                "analysis luad R20\n", "analysis luad M|R20\n", "analysis ucec OP\n",
                "analysis ucec M|C3 \n", "arm ucec\n", "ARM ucec A_C3\n"]
        ok("malformed, repeated, whitespace-padded or non-closable lines are refused (%d cases)"
           % len(bads), all(_raises(ValueError, parse_closed_parts_text, b) for b in bads))
        inc = ["arm ucec A_R20cap\n", "arm luad A_R20\n", "analysis ucec C3\n",
               "analysis luad M|C3\n", "arm ucec A_R20cap\narm luad A_R20\n"]
        ok("consistency: A_R20cap iff A_R20, and analysis C3 iff analysis M|C3, per cohort",
           all(_raises(ValueError, parse_closed_parts_text, b) for b in inc))
        ok("closing A_R20 and A_R20cap (or C3 and M|C3) together is allowed, per cohort",
           parse_closed_parts_text("arm ucec A_R20\narm ucec A_R20cap\n")["arms"]
           == {("ucec", "A_R20"), ("ucec", "A_R20cap")}
           and parse_closed_parts_text("analysis luad C3\nanalysis luad M|C3\n")["analyses"]
           == {("luad", "C3"), ("luad", "M|C3")})
        fcp = os.path.join(tmp, "CLOSED_PARTS_crlf.txt")
        with open(fcp, "wb") as fh:
            fh.write(b"# h\r\narm ucec A_C3\r\n")
        ok("a CRLF CLOSED_PARTS file is fatal", _raises(SystemExit, read_closed_parts, fcp))
        fbom = os.path.join(tmp, "CLOSED_PARTS_bom.txt")
        with open(fbom, "wb") as fh:
            fh.write(b"\xef\xbb\xbf# h\n")
        ok("a BOM or a non-UTF-8 control file is fatal",
           _raises(SystemExit, read_closed_parts, fbom)
           and _raises(SystemExit, read_pinned_hosts, fbom))
        fbad = os.path.join(tmp, "CLOSED_PARTS_latin1.txt")
        with open(fbad, "wb") as fh:
            fh.write(b"# h \xe9\n")
        ok("a non-UTF-8 control file is fatal, not a traceback",
           _raises(SystemExit, read_closed_parts, fbad))
        ok("arm_closed", arm_closed(cpx, "luad", "A_R20") and not arm_closed(cpx, "luad", "A_C3"))

        # PINNED_HOSTS.txt (PS 9 "Scripts", H-V1)
        ok("PINNED_HOSTS: comments only pins nothing",
           parse_pinned_hosts_text("# h\n\n# only comments\n") == {})
        ok("PINNED_HOSTS: one line per cohort",
           parse_pinned_hosts_text("# h\npinned ucec s006\npinned luad s002\n")
           == {"ucec": "s006", "luad": "s002"})
        pbad = ["pinned ucec\n", "pinned ucec s006 x\n", "pinned ucec  s006\n",
                "pinned ucec s006 \n", " pinned ucec s006\n", "pinned brca s006\n",
                "pinned ucec \n", "pin ucec s006\n", "pinned ucec s006\npinned ucec s002\n",
                "pinned\tucec\ts006\n", "  # indented comment\n", " \n",
                "pinned ucec s006\r\n", "pinned ucec s0\t6\n"]
        ok("PINNED_HOSTS: malformed, padded, repeated or tabbed lines are refused (%d cases)"
           % len(pbad), all(_raises(ValueError, parse_pinned_hosts_text, b) for b in pbad))
        ok("PINNED_HOSTS: the running host must equal the host the file gives for the cohort",
           check_pinned_host({"ucec": host()}, "ucec", host()) is None
           and _raises(SystemExit, check_pinned_host, {"ucec": "s006"}, "ucec", "s002")
           and _raises(SystemExit, check_pinned_host, {"luad": host()}, "ucec", host())
           and _raises(SystemExit, check_pinned_host, {}, "ucec", host()))

        # sidecar: production verifies PINNED_HOSTS.txt too
        class _A:
            pass
        a_sc = _A()
        a_sc.prespec, a_sc.saved_matrix = "p", "m"
        a_sc.closed_parts, a_sc.pinned_hosts = "c", "h"
        csf = cluster_sidecar_files(a_sc, "ucec", {k_: "o_" + k_ for k_ in FROZEN_CODE_SIDECAR_REL})
        ok("the sidecar check covers PINNED_HOSTS.txt and CLOSED_PARTS.txt (production)",
           csf.get(PINNED_HOSTS_REL) == "h" and csf.get(CLOSED_PARTS_REL) == "c"
           and PRESPEC_REL in csf and SCRIPT_REL in csf and len(csf) == 2 + 3 + 5 + 2 + 1)
        pfile = os.path.join(tmp, "PINNED_HOSTS.txt")
        with open(pfile, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# test\npinned ucec %s\n" % host())
        pfile_other = os.path.join(tmp, "other", "PINNED_HOSTS.txt")
        os.makedirs(os.path.dirname(pfile_other))
        with open(pfile_other, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# test\npinned ucec not-this-host\n")
        ent_p = {PINNED_HOSTS_REL: sha256(pfile)}
        ok("sidecar check: a PINNED_HOSTS.txt with another hash is fatal",
           check_sidecar_entries(ent_p, {PINNED_HOSTS_REL: pfile})[0] == {PINNED_HOSTS_REL: sha256(pfile)}
           and _raises(SystemExit, check_sidecar_entries, ent_p, {PINNED_HOSTS_REL: pfile_other})
           and _raises(SystemExit, check_sidecar_entries, {}, {PINNED_HOSTS_REL: pfile}))

        # PS 9 "Defects": the texts logged on a stop and on an exception
        et = exception_text(3, "c1_cov_ucec_A_R20")
        st = scheduler_stop_text(15, 3, "c1_cov_ucec_A_R20")
        ok("an exception is reported before any resubmission; environment failures are "
           "resubmitted unchanged and recorded; anything else is a script defect",
           "REPORT TO THE DRAFTING SESSION BEFORE ANY RESUBMISSION" in et
           and "environment failure (disk, memory, file system)" in et
           and "resubmitted unchanged and the exception is recorded" in et
           and "handled as a script defect" in et and "Resubmit" not in et)
        ok("a scheduler stop keeps the 'resubmit unchanged' message",
           "Resubmit the same command unchanged" in st and "not a script defect" in st)

        # production refusal of a closed arm, driven through mode_confirmatory with a fake env
        croot = os.path.join(tmp, "fakezw", "ucec_c1")
        cout = os.path.join(croot, "results", "posthoc", "c1_covariates")
        os.makedirs(cout)
        cpf = os.path.join(tmp, "CLOSED_PARTS.txt")
        with open(cpf, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# test\narm ucec A_C3\narm ucec A_R20cap\narm ucec A_R20\n")
        saved_env = {k: os.environ.get(k) for k in REQUIRED_ENV}
        os.environ.update(MORPHO_ROOT=croot, MORPHO_RNA_MANIFEST=os.path.join(tmp, "m"),
                          MORPHO_SLIDE_MAP=os.path.join(tmp, "s"),
                          MORPHO_WSI_EMB_DIR=os.path.join(tmp, "w"), MORPHO_OUT=cout)
        try:
            refused = []
            for argv in (["--mode", "arm", "--cohort", "ucec", "--arm", "A_C3"],
                         ["--mode", "placebo", "--cohort", "ucec", "--draws", "1-30"]):
                a_ = parse_args(argv + ["--closed-parts", cpf, "--pinned-hosts", pfile])
                refused.append(_raises(Refused, mode_confirmatory, a_, lambda *_x: None))
            ok("production refuses a closed arm and closed placebos (before any import)",
               all(refused) and os.listdir(cout) == [])
            a_ = parse_args(["--mode", "arm", "--cohort", "ucec", "--arm", "A_C3",
                             "--closed-parts", cpf, "--pinned-hosts", pfile_other])
            ok("production refuses on a host other than the one PINNED_HOSTS.txt gives",
               _raises(SystemExit, mode_confirmatory, a_, lambda *_x: None) and os.listdir(cout) == [])
            ok("production refuses a cohort that PINNED_HOSTS.txt does not pin",
               _raises(SystemExit, check_pinned_host, read_pinned_hosts(pfile), "luad", host()))
            a_ = parse_args(["--mode", "arm", "--cohort", "ucec", "--arm", "A_R20",
                             "--pinned-hosts", pfile])
            ok("production refuses without --closed-parts",
               _raises(SystemExit, mode_confirmatory, a_, lambda *_x: None))
            ok("production modes need --pinned-hosts; there is no free --pinned-host",
               _raises(SystemExit, parse_args, ["--mode", "arm", "--cohort", "ucec", "--arm",
                                                "A_C3", "--closed-parts", cpf])
               and _raises(SystemExit, parse_args, ["--mode", "arm", "--cohort", "ucec", "--arm",
                                                    "A_C3", "--closed-parts", cpf,
                                                    "--pinned-host", host()])
               and _raises(SystemExit, parse_args, ["--mode", "gate", "--cohort", "ucec",
                                                    "--pinned-hosts", pfile]))
        finally:
            for k_, v_ in saved_env.items():
                if v_ is None:
                    os.environ.pop(k_, None)
                else:
                    os.environ[k_] = v_

        log("-- 9. files: exclusive writes, attempts without limit, covariates record")
        f1 = os.path.join(tmp, "c1_cov_ucec_A_C3.npz")
        h = save_matrix_npz(f1, genes_u, np.zeros((len(genes_u), 3)), ["a", "b"])
        with np.load(f1, allow_pickle=False) as zz:
            ok("npz round trip: unicode genes, float64 M, cases",
               zz["genes"].dtype.kind == "U" and list(zz["genes"]) == list(genes_u)
               and zz["M"].dtype == np.float64 and list(zz["cases"]) == ["a", "b"]
               and sorted(zz.files) == ["M", "cases", "genes"]
               and not os.path.exists(f1 + ".partial") and h == sha256(f1))
        ok("an existing output is never overwritten",
           _raises(Refused, save_matrix_npz, f1, genes_u, np.zeros((len(genes_u), 3)), ["a"]))
        ok("a non-finite matrix is never written",
           _raises(ValueError, save_matrix_npz, os.path.join(tmp, "c1_cov_bad.npz"), genes_u,
                   np.full((len(genes_u), 3), np.nan), ["a"])
           and not os.path.exists(os.path.join(tmp, "c1_cov_bad.npz")))
        adir = os.path.join(tmp, "attempts")
        os.makedirs(adir)
        stem = "c1_cov_luad_V_R20"
        fin = os.path.join(adir, stem + ".npz")
        cpart = covariates_partial(adir, stem)
        ks, kept_all = [], []
        for i in range(5):                                   # five interrupted attempts
            k_, kept = begin_attempt(adir, stem, [fin], [cpart], dict(x=1))
            ks.append(k_)
            kept_all += kept
            with open(fin + ".partial", "wb") as fh:
                fh.write(b"half %d" % i)
            if i % 2 == 0:
                with open(cpart, "wb") as fh:
                    fh.write(b"cov %d" % i)
        k6, kept6 = begin_attempt(adir, stem, [fin], [cpart], dict(x=1))
        kept_all += kept6
        ok("attempts are not limited: 1..6, each recorded",
           ks + [k6] == [1, 2, 3, 4, 5, 6] and attempts_recorded(adir, stem) == [1, 2, 3, 4, 5, 6])
        ok("partial files of every interrupted attempt are kept under distinct names",
           len(kept_all) == 5 + 3 and len(set(kept_all)) == len(kept_all)
           and all(os.path.exists(os.path.join(adir, x)) for x in kept_all)
           and open(os.path.join(adir, stem + ".npz.partial.interrupted_before_attempt3"),
                    "rb").read() == b"half 1"
           and not os.path.exists(fin + ".partial") and not os.path.exists(cpart))
        with open(os.path.join(adir, stem + ".attempt6.json"), encoding="utf-8") as fh:
            a6 = json.load(fh)
        ok("the attempt record lists the previous attempts and the kept partials",
           a6["attempt"] == 6 and a6["previous_attempts"] == [1, 2, 3, 4, 5]
           and a6["interrupted_partials_kept"] == kept6)
        save_matrix_npz(fin, genes_u, np.zeros((len(genes_u), 2)), ["a"])
        ok("a completed run is never repeated",
           _raises(Refused, begin_attempt, adir, stem, [fin], [cpart], {}))
        lone = os.path.join(adir, "c1_cov_luad_V_C3.npz")
        with open(lone + ".partial", "wb") as fh:
            fh.write(b"?")
        ok("a leftover partial without any attempt record is refused (stop and report)",
           _raises(Refused, begin_attempt, adir, "c1_cov_luad_V_C3", [lone], [], {}))
        os.remove(lone + ".partial")
        with open(os.path.join(tmp, "c1_cov_x.partial"), "wb") as fh:
            fh.write(b"x")
        ok("a leftover .partial blocks a write",
           _raises(Refused, write_text_exclusive, os.path.join(tmp, "c1_cov_x"), "y"))
        os.remove(os.path.join(tmp, "c1_cov_x.partial"))
        # placebo range overlap
        pdir = os.path.join(tmp, "placebo")
        os.makedirs(pdir)
        pfin = [os.path.join(pdir, cap_name("ucec", j)) for j in range(1, 31)]
        begin_attempt(pdir, placebo_stem("ucec", 1, 30), pfin, [], {})
        ok("placebo ranges may not overlap; disjoint ranges and the same range are allowed",
           _raises(Refused, check_draw_overlap, pdir, "ucec", 25, 40)
           and _raises(Refused, check_draw_overlap, pdir, "ucec", 1, 119)
           and not _raises(Refused, check_draw_overlap, pdir, "ucec", 31, 60)
           and not _raises(Refused, check_draw_overlap, pdir, "ucec", 1, 30)
           and not _raises(Refused, check_draw_overlap, pdir, "luad", 1, 119))
        # run lock (POSIX only)
        try:
            import fcntl  # noqa: F401
            l1 = RunLock(pdir, "c1_cov_ucec_A_C3")
            l1.acquire()
            l2 = RunLock(pdir, "c1_cov_ucec_A_C3")
            blocked = _raises(Refused, l2.acquire)
            l1.release()
            again = l2.acquire()
            l2.release()
            ok("a run lock stops a second job of the same run", blocked and again)
        except ImportError:
            log("  [SKIP] run lock (no fcntl on this platform; the cluster is Linux)")
        # covariates record
        cdir = os.path.join(tmp, "cov")
        os.makedirs(cdir)
        cases = ["C3L-%05d" % i for i in range(n)]
        minfo = dict(stroma=dict(used=["A"]), immune=dict(used=["B"]),
                     proliferation=dict(used=["C"]))
        pay = covariates_payload(cases, C3, R20s, pcs, minfo, False)
        hcov = ensure_covariates_file(cdir, "ucec", "c1_cov_ucec_A_C3", pay, lambda *_x: None)
        cfile = os.path.join(cdir, covariates_name("ucec"))
        with np.load(cfile, allow_pickle=False) as zc:
            meta = json.loads(str(zc["markers_json"]))
            ok("covariates record: keys cases, C3 (n x 3), R20 (n x 20), pcs (n x K), markers_json",
               sorted(zc.files) == sorted(COVARIATES_KEYS) and zc["C3"].shape == (n, 3)
               and zc["R20"].shape == (n, 20) and zc["pcs"].shape == (n, K)
               and np.array_equal(zc["pcs"], pcs) and zc["cases"].dtype.kind == "U"
               and meta["markers"]["immune"]["used"] == ["B"] and meta["c3_stored"] == "computed"
               and hcov == sha256(cfile))
        ok("covariates record: a second job with identical covariates passes",
           ensure_covariates_file(cdir, "ucec", "c1_cov_ucec_V_R20", pay, lambda *_x: None) == hcov)
        pay2 = covariates_payload(cases, C3, R20s, pcs + 1e-12, minfo, False)
        ok("covariates record: different pcs are fatal",
           _raises(SystemExit, ensure_covariates_file, cdir, "ucec", "c1_cov_ucec_V_C3", pay2,
                   lambda *_x: None))
        ok("covariates record: no partial left behind",
           [f for f in os.listdir(cdir) if f.endswith(".partial")] == [])
        pay3 = covariates_payload(cases, None, R20s, pcs, minfo, True)
        ok("covariates record: C3 stored as n x 0 when C3 is closed", pay3["C3"].shape == (n, 0)
           and json.loads(str(pay3["markers_json"]))["c3_stored"] == "closed")
        # a crashed covariates write is kept by the next attempt of its own run (minor 8)
        c2dir = os.path.join(tmp, "cov2")
        os.makedirs(c2dir)
        st2 = "c1_cov_luad_A_R20"
        f2 = os.path.join(c2dir, st2 + ".npz")
        begin_attempt(c2dir, st2, [f2], [covariates_partial(c2dir, st2)], {})
        with open(covariates_partial(c2dir, st2), "wb") as fh:
            fh.write(b"crashed mid-write")
        _k, kept = begin_attempt(c2dir, st2, [f2], [covariates_partial(c2dir, st2)], {})
        h2 = ensure_covariates_file(c2dir, "luad", st2, pay, lambda *_x: None)
        ok("an interrupted covariates write is kept aside by the run's next attempt and the "
           "record is then written", kept == [st2 + ".covariates_record.partial"
                                              ".interrupted_before_attempt2"]
           and h2 == sha256(os.path.join(c2dir, covariates_name("luad"))))

        log("-- 10. manifest (one file for both trees; hashes only)")
        mu, ml = os.path.join(tmp, "mu"), os.path.join(tmp, "ml")
        os.makedirs(mu)
        os.makedirs(ml)
        none_closed = dict(arms=set(), analyses=set())
        for coh, d_ in (("ucec", mu), ("luad", ml)):
            for nm in expected_outputs(coh, none_closed):
                write_text_exclusive(os.path.join(d_, nm), nm)
            write_text_exclusive(os.path.join(d_, "c1_cov_%s_A_C3.attempt1.json" % coh), "{}")
        ok("expected outputs per cohort: 4 arms + 119 placebos + covariates record",
           len(expected_outputs("ucec", none_closed)) == 4 + 119 + 1)
        with open(os.path.join(mu, "unrelated.txt"), "w") as fh:
            fh.write("not listed")
        dirs = dict(ucec=mu, luad=ml)
        text, entries = build_manifest_text(dirs, none_closed)
        lines = text.rstrip("\n").split("\n")
        ok("manifest: lines '<sha256> *<bare file name>' only, both trees, sorted",
           all(re.match(r"^[0-9a-f]{64} \*c1_cov_[^/\\ ]+$", ln) for ln in lines)
           and len(lines) == 2 * (124 + 1) and "unrelated.txt" not in text
           and lines == sorted(lines, key=lambda s: s.split("*", 1)[1])
           and "c1_cov_luad_A_R20cap_119.npz" in entries
           and not any(ln.startswith(("not_completed", "untestable")) for ln in lines))
        mf = os.path.join(mu, MANIFEST_NAME)
        with open(mf + ".partial", "wb") as fh:
            fh.write(b"interrupted manifest")
        ok("a leftover manifest partial blocks the write without --force-rewrite",
           _raises(Refused, write_manifest, mf, text, False))
        hm = write_manifest(mf, text, True)
        ok("--force-rewrite moves the leftover aside while no manifest exists",
           hm == sha256(mf) and not os.path.exists(mf + ".partial"))
        text2, _ = build_manifest_text(dirs, none_closed)
        ok("the manifest does not list itself or its leftovers", MANIFEST_NAME not in text2)
        ok("the manifest is never rewritten, with or without --force-rewrite",
           _raises(Refused, write_manifest, mf, text, False)
           and _raises(Refused, write_manifest, mf, text, True))
        pp = os.path.join(ml, "c1_cov_luad_V_C3.npz.partial")
        with open(pp, "wb") as fh:
            fh.write(b"running")
        ok("no manifest while a .partial exists", _raises(Refused, build_manifest_text, dirs, none_closed))
        os.remove(pp)
        write_text_exclusive(os.path.join(ml, "c1_cov_ucec_A_R20.attempt1.json"), "{}")
        write_text_exclusive(os.path.join(mu, "c1_cov_ucec_A_R20.attempt1.json"), "{}")
        ok("a file name in both trees is refused", _raises(Refused, build_manifest_text, dirs, none_closed))
        os.remove(os.path.join(ml, "c1_cov_ucec_A_R20.attempt1.json"))
        os.remove(os.path.join(ml, cap_name("luad", 77)))
        ok("a missing output (not closed) is refused",
           _raises(Refused, build_manifest_text, dirs, none_closed))
        cl = dict(arms={("luad", "A_R20cap")}, analyses=set())
        ok("an output of a closed arm that exists is refused",
           _raises(Refused, build_manifest_text, dirs, cl))
        for j in range(1, N_CAP + 1):
            if os.path.exists(os.path.join(ml, cap_name("luad", j))):
                os.remove(os.path.join(ml, cap_name("luad", j)))
        t3, e3 = build_manifest_text(dirs, cl)
        ok("closed placebos are neither expected nor listed",
           not any("luad_A_R20cap" in k for k in e3) and "c1_cov_ucec_A_R20cap_119.npz" in e3)
        cl_all = dict(arms={("luad", a) for a in CLOSABLE_ARMS}, analyses=set())
        ok("a cohort with every arm closed expects nothing (not even the covariates record)",
           expected_outputs("luad", cl_all) == [])

        log("-- 11. sidecar")
        sc = os.path.join(tmp, "side.sha256")
        fa, fb = os.path.join(tmp, "fa"), os.path.join(tmp, "fb")
        for p_, t_ in ((fa, "A"), (fb, "B")):
            with open(p_, "w") as fh:
                fh.write(t_)
        with open(sc, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# written 2026-11-21T00:00:00Z\n%s *%s\n%s  %s\n%s *%s\n"
                     % (sha256(fa), PRESPEC_REL, sha256(fb).upper(), SCRIPT_REL, "0" * 64,
                        "server_export/results/local_only.csv"))
        ent = read_sidecar(sc)
        ok("sidecar parsing ('*' and two-space forms, comments, case)",
           ent == {PRESPEC_REL: sha256(fa), SCRIPT_REL: sha256(fb),
                   "server_export/results/local_only.csv": "0" * 64})
        chk, nh = check_sidecar_entries(ent, {PRESPEC_REL: fa, SCRIPT_REL: fb})
        ok("sidecar check: present entries match, others reported as not on the cluster",
           sorted(chk) == sorted([PRESPEC_REL, SCRIPT_REL]) and nh == ["server_export/results/local_only.csv"])
        ok("sidecar check: a mismatch is fatal",
           _raises(SystemExit, check_sidecar_entries, ent, {PRESPEC_REL: fb}))
        ok("sidecar check: an unlisted required file is fatal",
           _raises(SystemExit, check_sidecar_entries, ent, {CLOSED_PARTS_REL: fa}))
        ok("sidecar check: a listed required file missing on the cluster is fatal",
           _raises(SystemExit, check_sidecar_entries, ent, {PRESPEC_REL: os.path.join(tmp, "nope")}))
        with open(sc, "a", encoding="utf-8", newline="\n") as fh:
            fh.write("%s *%s\n" % (sha256(fa), PRESPEC_REL))
        ok("a sidecar that lists a path twice is fatal", _raises(SystemExit, read_sidecar, sc))

        log("-- 12. score table, patients, listing")
        pts = ["C3L-%05d" % i for i in range(10)]
        S3 = rs.randn(10, 3) / 7.0
        tsv = os.path.join(tmp, "c1_cov_discovery_scores_ucec.tsv")
        write_text_exclusive(tsv, scores_tsv_text(pts, S3, ["markers stroma 24/24", "host x"]))
        rp, rS = read_scores_tsv(tsv)
        raw = open(tsv, "rb").read()
        ok("score TSV: '#' lines, header, repr floats round-trip exactly, LF, no BOM",
           rp == pts and np.array_equal(rS, S3) and b"\r" not in raw and not raw.startswith(b"\xef\xbb\xbf"))
        npz = os.path.join(tmp, "inputs_syn.npz")
        np.savez_compressed(npz, patients=np.array(pts), rna_pcs20=rs.randn(10, 20))
        saved = dict(COHORTS["ucec"])
        COHORTS["ucec"]["inputs_sha"] = sha256(npz)
        try:
            pf = os.path.join(tmp, "c1_cov_discovery_patients_ucec.txt")
            ns = argparse.Namespace(cohort="ucec", inputs_npz=npz, out=pf)
            mode_export_patients(ns, lambda *_a: None)
            ok("export-patients writes the inputs file's patients, one per line",
               read_patient_file(pf) == pts and open(pf, "rb").read() == ("\n".join(pts) + "\n").encode())
        finally:
            COHORTS["ucec"] = saved
        rows, trunc = list_tree([tmp, os.path.join(tmp, "missing_dir")], max_depth=2)
        ok("result-folder listing (names, sizes, mtimes; missing roots recorded)",
           any(r[0].endswith("c1_cov_ucec_A_C3.npz") and r[3] == "flag" for r in rows)
           and any(r[1] == "missing" for r in rows) and not trunc)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    log("-- 13. constants and this file")
    ok("published AUCs (published quantities) and usable set sizes pinned",
       COHORTS["ucec"]["published_auc"] == 0.5974425564900907
       and COHORTS["luad"]["published_auc"] == 0.6439668501223573
       and (COHORTS["ucec"]["n_sel_usable"], COHORTS["ucec"]["n_bg_usable"]) == (2536, 7610)
       and (COHORTS["luad"]["n_sel_usable"], COHORTS["luad"]["n_bg_usable"]) == (2287, 8119))
    raw = open(os.path.abspath(__file__), "rb").read()
    ok("this script: UTF-8 without BOM, LF line endings",
       not raw.startswith(b"\xef\xbb\xbf") and b"\r" not in raw)
    _ = PCA  # imported above to fail early if sklearn is missing

    n_fail = sum(1 for _n, c in checks if not c)
    log("\nselftest: %d checks, %d failed, %.1f s" % (len(checks), n_fail, time.time() - t0))
    log("SELFTEST PASSED" if n_fail == 0 else "SELFTEST FAILED")
    return 0 if n_fail == 0 else 1


# ---------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------
MODES = ("selftest", "export-patients", "blind", "gate", "gate-r20", "export-scores", "arm",
         "placebo", "manifest")


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="post hoc C1 covariate re-reads, cluster side",
                                 allow_abbrev=False)
    ap.add_argument("--mode", required=True, choices=MODES)
    ap.add_argument("--cohort", choices=sorted(COHORTS), default=None)
    ap.add_argument("--arm", choices=sorted(ARMS), default=None)
    ap.add_argument("--draws", default=None, help="placebo draws a-b within 1-119, e.g. 1-30")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--saved-matrix", default=None)
    ap.add_argument("--sidecar", default=None)
    ap.add_argument("--sidecar-sha256", default=None)
    ap.add_argument("--prespec", default=None)
    ap.add_argument("--closed-parts", default=None, help="CLOSED_PARTS.txt (listed in the sidecar)")
    ap.add_argument("--pinned-hosts", default=None,
                    help="PINNED_HOSTS.txt (listed in the sidecar; modes arm, placebo)")
    ap.add_argument("--force-rewrite", action="store_true",
                    help="manifest: move a leftover manifest .partial aside (only while no "
                         "manifest exists)")
    ap.add_argument("--list-roots", default=None, help="comma-separated; blind listing roots")
    ap.add_argument("--patients", default=None)
    ap.add_argument("--inputs-npz", default=None)
    ap.add_argument("--out", default=None, help="export-patients output file")
    # published data arguments (c1_run_test.py / c1_run_test_luad.py / plex re-read wrappers)
    ap.add_argument("--protein", default=None)
    ap.add_argument("--discovery-results", default=None)
    ap.add_argument("--discovery-slide-map", default=None)
    ap.add_argument("--case-ids", default=None)
    ap.add_argument("--rule", default=None)
    ap.add_argument("--rule-sha256", default=None)
    ap.add_argument("--protein-sha256", default=None)
    ap.add_argument("--rna-manifest-sha256", default=None)
    ap.add_argument("--discovery-sitepack", default=None)
    ap.add_argument("--discovery-xwalk", default=None)
    ap.add_argument("--discovery-rna-manifest", default=None)
    ap.add_argument("--batch-labels", default=None)
    args = ap.parse_args(argv)
    if args.mode != "selftest" and not args.cohort:
        ap.error("--cohort is required for mode " + args.mode)
    if args.pinned_hosts and args.mode not in ("arm", "placebo", "manifest"):
        ap.error("--pinned-hosts belongs to the modes arm, placebo and manifest")
    if args.mode in ("arm", "placebo") and not args.pinned_hosts:
        ap.error("--mode %s needs --pinned-hosts" % args.mode)
    if args.mode == "arm" and not args.arm:
        ap.error("--mode arm needs --arm")
    if args.mode == "placebo" and not args.draws:
        ap.error("--mode placebo needs --draws a-b")
    if args.force_rewrite and args.mode != "manifest":
        ap.error("--force-rewrite belongs to --mode manifest")
    return args


def main(argv=None):
    args = parse_args(argv)
    log = Tee()
    try:
        if args.mode == "selftest":
            return selftest(args, log)
        if args.mode == "export-patients":
            return mode_export_patients(args, log)
        if args.mode == "manifest":
            return mode_manifest(args, log)
        if args.mode in ("gate-r20", "export-scores"):
            return mode_discovery(args, log)
        return mode_confirmatory(args, log)
    except Refused as e:
        log("REFUSED: %s" % e)
        return 4
    except SystemExit as e:
        if isinstance(e.code, str):
            log(e.code)
        raise
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
