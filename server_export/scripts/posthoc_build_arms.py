#!/usr/bin/env python3
"""
posthoc_build_arms.py -- build the variant-input ARM files for the M2 (embedding
medium) and Q (disqualified-aliquot) post hoc re-runs of the discovery estimator
(residual_analysis.py), per review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md
(being written now -- the frozen doc number is a fixed reference, not a claim that
file already exists). This script DOES NOT MODIFY residual_analysis.py, does not
call its main(), and does not compute any residual/incremental-R2/permutation
statistic; it only writes per-arm copies of the two files residual_analysis.py
reads through MORPHO_SLIDE_MAP / MORPHO_ALIQUOT_XWALK, so that a later, separate
job (lsf_posthoc_residual.sh + residual_analysis.py, unmodified) can be pointed at
one arm at a time.

POST HOC. This whole family of runs is post hoc relative to the discovery
analysis already reported in the paper: every arm is reported whatever the
result (section 0 of the PRESPEC doc), no arm is added after a result is seen,
and nothing here can change a C1-UCEC or C1-LUAD confirmatory reading (those are
governed separately by review/C1_*_FROZEN_RULE_*.md and never read these files).

Arms written, under <--out-root>/<cohort>/ (default <--zw-root>/scripts/pinned/posthoc/<cohort>/):
  R           control re-run: unmodified copies of the CURRENT pinned slide map
              and aliquot crosswalk for this cohort (resolved via morpho_env.sh,
              never reimplemented here -- see resolve_morpho_env()). Built for
              EVERY cohort.
  F           (built where PRESPEC M2 qualifies the cohort: k >= 3) the slide map with every OCT Primary
              Tumor slide's sample_type changed to 'Excluded_OCT'. Verified
              against residual_analysis.load_wsi_embeddings(): that function
              keeps a (case, slide) row only when sample_type == 'Primary Tumor'
              and the slide has a matching embedding file; anything else
              (including 'Excluded_OCT') falls into its "else: normal/other ->
              exclude" branch (residual_analysis.py, the per-slide loop inside
              load_wsi_embeddings). A mixed patient (>=1 FFPE + >=1 OCT analysed
              slide) keeps its FFPE slide(s) and a non-empty pooled embedding; an
              OCT-only patient loses every analysed slide and drops out of
              `common` in residual_analysis.py's own intersection logic -- this
              script does not simulate that drop, it only relabels the source
              row, exactly as the PRESPEC's M2 calls for.
  D01..D19    (built with F, on the same M2 condition) all 19 arms relabel the Primary Tumor
              slides of k RANDOMLY drawn analysed patients to 'Excluded_D' (same
              mechanism as F, different patients), k = the number of patients
              classified oct_only under the CURRENT map (see classify below).
              Draw for arm s (s = 0..18): np.random.RandomState(1000 + s)
              .choice(sorted_analysed_case_ids, k, replace=False), exactly as
              specified (no other seed, no reseeding, no replacement). A cohort with
              k < 3, or with a medium join below 95%, gets no F or D arm and a
              DISCLOSURES.txt naming the criterion it failed (PRESPEC M3).
  Q           (built where >=1 analysed aliquot is Disqualified) the aliquot crosswalk with every row
              whose aliquot_id has PDC aliquot_status 'Disqualified' (per
              --pdc-aliquot-status) removed. A patient loses its row only if ALL
              of its xwalk rows were Disqualified; a patient with >=1 remaining
              Qualified aliquot keeps that row.

'Analysed patients', for sizing k and for the D-arm draw pool, is defined EXACTLY
as the task spec: the Primary Tumor slides (from the CURRENT, unmodified slide
map) that have a matching embedding file under MORPHO_WSI_EMB_DIR, rolled up to
case_submitter_id. This script computes that set TWO independent ways and FATALs
if they disagree:
  (a) a local, read-only scan of MORPHO_WSI_EMB_DIR using the identical filename
      regex residual_analysis.load_wsi_embeddings() uses (copied verbatim as
      SLIDE_RE below; residual_analysis.py itself is never imported for this
      -- no torch load, no mean-pool, just a filename match), rolled up via the
      slide map's case_submitter_id column;
  (b) importing residual_analysis.py (frozen, read-only -- see the hard rule in
      the governing task; never edited) and calling its own
      load_wsi_embeddings() with every MORPHO_* env var resolved by
      morpho_env.sh, taking the index of the returned per-case matrix.
This is NOT the final analysed-n for any arm's run (that also intersects RNA and
protein, which this script never loads) -- the docstring note in the task is
honoured literally: "the final n is reported by the run itself" (residual_analysis.py's
own `common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))`).

Inputs (every one is sha256-asserted before use; FATAL on mismatch, no silent
fallback):
  - residual_analysis.py, morpho_env.sh       (frozen code; hardcoded constants)
  - the cohort's CURRENT MORPHO_SLIDE_MAP / MORPHO_ALIQUOT_XWALK, resolved by
    sourcing morpho_env.sh (never reimplemented)        (--slide-map-sha256 / --xwalk-sha256)
  - pinned/slide_embedding_medium.tsv                   (--slide-embedding-medium + its --*-sha256;
                                                           required only for cohort in {ucec, luad})
  - pinned/pdc_aliquot_status.tsv (module D's output)    (--pdc-aliquot-status + its --*-sha256;
                                                           required only for cohort in {ucec, luad})

Outputs: <out-root>/<cohort>/arm_<ARM>_{slide_map,xwalk}.tsv (LF, UTF-8, same
columns as the source file, nothing reordered) plus one
<out-root>/<cohort>/MANIFEST.tsv with columns
  cohort  arm  field  file  sha256  n_patients_dropped  dropped_case_ids  note
(field is the MORPHO_* env var that arm's file overrides -- MORPHO_SLIDE_MAP or
MORPHO_ALIQUOT_XWALK -- which is exactly what lsf_posthoc_residual.sh needs to
know which variable to export). The manifest is also printed to stdout.

No permutation, no residualisation, no PCA happens in this file -- the only
randomness anywhere is the D-arm patient draw, and every draw uses a seeded
numpy RandomState as specified above; nothing here is unseeded.

Import-order note: this script never imports openslide, torch or timm directly.
It DOES import residual_analysis.py (which imports torch at module scope) for
the cross-check in (b) above; residual_analysis.py's own env-vars-at-import-time
requirement is honoured (every MORPHO_* var is set via os.environ BEFORE that
import, inside main(), not at module top level -- see compute_analysed_cases_via_RA()).

Run (HPC; this script takes no --rule/--mode -- it is not the confirmatory test,
just an input builder):
  /public/home/fjhui/miniconda3/bin/python -u posthoc_build_arms.py \\
    --cohort ucec \\
    --slide-map-sha256     <sha256 of the CURRENT pinned/slide_type_map_ucec.tsv> \\
    --xwalk-sha256          <sha256 of the CURRENT pinned/aliquot_to_case_tumor_ucec.tsv> \\
    --slide-embedding-medium        /public/home/fjhui/ZW/scripts/pinned/slide_embedding_medium.tsv \\
    --slide-embedding-medium-sha256 <sha256 of that file> \\
    --pdc-aliquot-status            /public/home/fjhui/ZW/scripts/pinned/pdc_aliquot_status.tsv \\
    --pdc-aliquot-status-sha256     <sha256 of that file>

For a cohort outside {ucec, luad} (ccrcc, gbm, pdac), --slide-embedding-medium /
--pdc-aliquot-status and their --*-sha256 may be omitted; only arm R is built,
and the script says so.

--dry-run computes and prints everything (including the manifest) without
writing any file under --out-root; use it to sanity-check k and the case
breakdown before committing arm files other jobs will read.
"""
import argparse
import hashlib
import os
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Frozen constants
# ---------------------------------------------------------------------------
ZW_DEFAULT = "/public/home/fjhui/ZW"
COHORTS = ("ccrcc", "luad", "ucec", "gbm", "pdac")
M2_MIN_OCT_ONLY = 3    # PRESPEC M2: a cohort qualifies for the medium arms at k >= 3
M0_MIN_JOIN = 0.95     # PRESPEC M0: below this the cohort gets a disclosure, not arms
N_D_ARMS = 19          # frozen: s = 0..18, D01..D19
D_SEED_BASE = 1000     # frozen: RandomState(1000 + s)
EXCLUDED_OCT_LABEL = "Excluded_OCT"
EXCLUDED_D_LABEL = "Excluded_D"
DISQUALIFIED_STATUS = "Disqualified"
RESULT_SLIDE_MAP_FIELD = "MORPHO_SLIDE_MAP"
RESULT_XWALK_FIELD = "MORPHO_ALIQUOT_XWALK"

# sha256 of the local-mirror copies as read when this script was written
# (2026-10-02). Hardcoded, not a CLI flag -- these two are frozen CODE files
# (residual_analysis.py is explicitly in the "import, never edit" list; this
# script never writes to either), so a legitimate future change to either must
# update this constant deliberately, the same discipline c1_run_test_luad.py's
# own EXPECTED_SHA256["residual_analysis"] already uses (identical value).
EXPECTED_SHA256_RESIDUAL_ANALYSIS = "a54a0a56494d315e15ef075425fe2133ecd9f05193c22656c1145bffdaac4a0a"
EXPECTED_SHA256_MORPHO_ENV = "1cf031e05b088bea2dd850f6376c058d134660e4d62fdc780fcc29fd87b8bff7"

# Mirrors residual_analysis.py's load_wsi_embeddings() slide_re EXACTLY (copied,
# not imported, because that function has no accessor at slide granularity --
# it only returns the final per-case mean-pooled matrix). Do not change this
# without re-reading that function; the cross-check in main() will catch drift.
SLIDE_RE = re.compile(r"^(C3[A-Z]-\d{5}-\d+)")


# ---------------------------------------------------------------------------
# sha256 / IO helpers
# ---------------------------------------------------------------------------
def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_sha256(path, expected, label):
    if not os.path.isfile(path):
        sys.exit(f"FATAL: {label} not found at {path!r}.")
    got = sha256_of(path)
    if got != expected:
        sys.exit(f"FATAL: sha256 mismatch for {label} ({path}): got={got} want={expected}. "
                  f"Wrong file vintage -- refusing to build arms from an unverified input.")
    print(f"  hash OK: {label} = {got[:12]}... ({path})")
    return got


def write_tsv(df, path):
    """LF, UTF-8, columns unchanged -- no file this script writes should ever
    differ from its source except in the one column (sample_type) one arm's
    rows were relabelled in."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        df.to_csv(f, sep="\t", index=False, lineterminator="\n")


def read_pinned_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, comment="#").fillna("")


# ---------------------------------------------------------------------------
# Pure logic (unit-tested on synthetic data in the scratchpad -- see report)
# ---------------------------------------------------------------------------
def relabel_oct_slides(slide_map_df, medium_lookup, label=EXCLUDED_OCT_LABEL):
    """Arm F. medium_lookup: slide_submitter_id -> 'FFPE' | 'OCT' (already
    normalised; unknown slides are simply absent from the dict and left alone)."""
    df = slide_map_df.copy()
    is_pt = df["sample_type"] == "Primary Tumor"
    med = df["slide_submitter_id"].map(medium_lookup.get)
    is_oct = is_pt & (med == "OCT")
    relabelled_slides = sorted(df.loc[is_oct, "slide_submitter_id"].tolist())
    df.loc[is_oct, "sample_type"] = label
    return df, relabelled_slides


def draw_d_arm_cases(analysed_cases_sorted, k, seed):
    """Literal translation of the spec: RandomState(seed).choice(sorted_ids, k,
    replace=False). analysed_cases_sorted must already be sorted by the caller
    (not re-sorted here) -- that sortedness is itself part of the frozen draw."""
    rng = np.random.RandomState(seed)
    drawn = rng.choice(analysed_cases_sorted, k, replace=False)
    return sorted(drawn.tolist())


def relabel_cases(slide_map_df, case_ids, label=EXCLUDED_D_LABEL):
    """Arms D01..D19. Relabels ALL Primary Tumor rows of the given cases (not
    just OCT ones) -- these patients are meant to drop out entirely, mirroring
    a same-sized random exclusion rather than a medium-selective one."""
    df = slide_map_df.copy()
    case_set = set(case_ids)
    is_target = df["case_submitter_id"].isin(case_set) & (df["sample_type"] == "Primary Tumor")
    df.loc[is_target, "sample_type"] = label
    return df, int(is_target.sum())


def filter_disqualified(xwalk_df, status_lookup, disqualified_value=DISQUALIFIED_STATUS):
    """Arm Q. status_lookup: xwalk's aliquot_id value -> PDC aliquot_status
    string (join key chosen by the caller -- see load_status_lookup)."""
    df = xwalk_df.copy()
    status = df["aliquot_id"].map(status_lookup.get)
    is_disq = status == disqualified_value
    kept = df.loc[~is_disq].reset_index(drop=True)
    dropped_aliquots = sorted(df.loc[is_disq, "aliquot_id"].tolist())
    before_cases = set(df["case_id"])
    after_cases = set(kept["case_id"])
    dropped_cases = sorted(before_cases - after_cases)
    return kept, dropped_cases, dropped_aliquots


def classify_patients_by_medium(case_to_present_slides, medium_lookup):
    """For each analysed case, the set of known media among its PRESENT
    (embedding-bearing) Primary Tumor slides decides oct_only / ffpe_only /
    mixed / unknown. A case with some labelled and some unlabelled slides is
    classified from the labelled ones only (the unlabelled slide neither proves
    nor disproves OCT-only-ness, so it is not allowed to turn a true oct_only
    case into 'unknown' nor vice versa)."""
    out = {}
    for case, slides in case_to_present_slides.items():
        labels = {medium_lookup[s] for s in slides if s in medium_lookup}
        if labels == {"OCT"}:
            out[case] = "oct_only"
        elif labels == {"FFPE"}:
            out[case] = "ffpe_only"
        elif labels == {"OCT", "FFPE"}:
            out[case] = "mixed"
        elif not labels:
            out[case] = "unknown"
        else:
            out[case] = "unexpected:" + ",".join(sorted(labels))
    return out


def present_slide_ids_from_emb_dir(emb_dir):
    """Filesystem scan mirroring residual_analysis.load_wsi_embeddings()'s own
    filename match (SLIDE_RE above) -- which slide_submitter_ids have an
    embedding file at all, independent of sample_type. Read-only; no .pt file
    is opened, only its name is matched."""
    out = set()
    for f in sorted(Path(emb_dir).glob("*.pt")):
        m = SLIDE_RE.match(f.stem)
        if m:
            out.add(m.group(1))
    return out


# ---------------------------------------------------------------------------
# HPC glue (not exercised by local syntax-check / unit tests: needs the real
# pinned tree, the real embeddings, morpho_env.sh's hardcoded /public/home/...
# paths, and torch)
# ---------------------------------------------------------------------------
def resolve_morpho_env(morpho_env_sh, cohort):
    """Source morpho_env.sh (never reimplemented) and harvest every MORPHO_*
    var it exports for this cohort, by running `source <file> <cohort> && env`
    in a fresh bash -- exactly the sourcing contract morpho_env.sh's own header
    documents ("Source this; do not execute it.")."""
    cmd = f"source {shlex.quote(morpho_env_sh)} {shlex.quote(cohort)} && env"
    proc = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"FATAL: `source {morpho_env_sh} {cohort}` failed (rc={proc.returncode}). "
                  f"stderr:\n{proc.stderr}")
    env = {}
    for line in proc.stdout.splitlines():
        if line.startswith("MORPHO_") and "=" in line:
            k, _, v = line.partition("=")
            env[k] = v
    required = ["MORPHO_ROOT", "MORPHO_RNA_MANIFEST", "MORPHO_ALIQUOT_XWALK",
                "MORPHO_SLIDE_MAP", "MORPHO_PROTEIN_TSV", "MORPHO_WSI_EMB_DIR"]
    missing = [k for k in required if k not in env]
    if missing:
        sys.exit(f"FATAL: morpho_env.sh {cohort} did not export {missing}; got keys {sorted(env)}.")
    return env


def compute_analysed_cases_via_ra(env, scripts_dir, scratch_out):
    """Cross-check (b): import the frozen residual_analysis.py (read-only;
    never edited) with every MORPHO_* var it reads at import time set from
    morpho_env.sh's resolution, MORPHO_OUT overridden to a private scratch
    subdirectory so merely counting cases never touches the published
    results/ tree, and return sorted(load_wsi_embeddings().index)."""
    os.environ["MORPHO_ROOT"] = env["MORPHO_ROOT"]
    os.environ["MORPHO_RNA_MANIFEST"] = env["MORPHO_RNA_MANIFEST"]
    os.environ["MORPHO_ALIQUOT_XWALK"] = env["MORPHO_ALIQUOT_XWALK"]
    os.environ["MORPHO_SLIDE_MAP"] = env["MORPHO_SLIDE_MAP"]
    os.environ["MORPHO_PROTEIN_TSV"] = env["MORPHO_PROTEIN_TSV"]
    os.environ["MORPHO_WSI_EMB_DIR"] = env["MORPHO_WSI_EMB_DIR"]
    os.environ["MORPHO_OUT"] = str(scratch_out)  # never the published results dir
    os.environ["MORPHO_N_WORKERS"] = env.get("MORPHO_N_WORKERS", "1")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    import residual_analysis as RA  # noqa: E402  (env must be set first; see above)
    wsi = RA.load_wsi_embeddings()
    return sorted(wsi.index.tolist())


# ---------------------------------------------------------------------------
# Loaders for the two new pinned inputs (schema not frozen anywhere yet --
# validated defensively, FATAL on anything unexpected rather than a silent
# wrong join)
# ---------------------------------------------------------------------------
def load_medium_lookup(path, cohort):
    df = read_pinned_tsv(path)
    if "cohort" in df.columns:
        df = df[df["cohort"].str.lower() == cohort.lower()]
    medium_col = next((c for c in ("embedding_medium", "Embedding_Medium", "medium")
                       if c in df.columns), None)
    if "slide_submitter_id" not in df.columns or medium_col is None:
        sys.exit(f"FATAL: {path} is missing required columns; have {list(df.columns)}, "
                  f"need 'slide_submitter_id' and one of "
                  f"embedding_medium/Embedding_Medium/medium.")
    lookup, bad = {}, []
    for sid, m in zip(df["slide_submitter_id"], df[medium_col]):
        mu = (m or "").strip().upper()
        if mu in ("", "UNKNOWN", "UNKNOWN_SLIDE_LEVEL", "NA", "N/A"):
            continue
        if mu not in ("FFPE", "OCT"):
            bad.append((sid, m))
            continue
        if sid in lookup and lookup[sid] != mu:
            sys.exit(f"FATAL: {path} has conflicting medium labels for slide {sid!r}: "
                      f"{lookup[sid]} vs {mu}.")
        lookup[sid] = mu
    if bad:
        sys.exit(f"FATAL: {path} has {len(bad)} row(s) with an embedding-medium value "
                  f"that is neither FFPE/OCT nor a recognised unknown marker, e.g. "
                  f"{bad[:5]}. Fix the file or extend this loader deliberately -- "
                  f"do not guess.")
    return lookup


def load_status_lookup(path, cohort):
    df = read_pinned_tsv(path)
    if "cohort" in df.columns:
        df = df[df["cohort"].str.lower() == cohort.lower()]
    # Join-key risk: channel_map_<cohort>.tsv (the likely source of this file)
    # carries BOTH an 'aliquot_id' column (a PDC UUID) and an 'aliquot_submitter_id'
    # column (the CPT-prefixed human id) -- and the aliquot crosswalk this script
    # joins against (build_aliquot_xwalk.py's 'aliquot_id' column) holds the
    # CPT-prefixed id, NOT the UUID. Prefer aliquot_submitter_id; the coverage
    # check in main() catches a wrong choice (or a genuinely bad file) either way.
    id_col = next((c for c in ("aliquot_submitter_id", "aliquot_id") if c in df.columns), None)
    if id_col is None or "aliquot_status" not in df.columns:
        sys.exit(f"FATAL: {path} is missing required columns; have {list(df.columns)}, "
                  f"need one of aliquot_submitter_id/aliquot_id plus aliquot_status.")
    lookup = {}
    for a, st in zip(df[id_col], df["aliquot_status"]):
        if a in lookup and lookup[a] != st:
            sys.exit(f"FATAL: {path} has conflicting aliquot_status for {id_col}={a!r}: "
                      f"{lookup[a]!r} vs {st!r}.")
        lookup[a] = st
    return lookup, id_col


# ---------------------------------------------------------------------------
# Arm builders -> manifest rows
# ---------------------------------------------------------------------------
def manifest_row(cohort, arm, field, path, n_dropped, dropped_ids, note, write):
    sha = sha256_of(path) if write else ""
    return dict(cohort=cohort, arm=arm, field=field, file=str(path),
                sha256=sha, n_patients_dropped=n_dropped,
                dropped_case_ids=";".join(dropped_ids), note=note)


def build_arm_r(cohort, slide_map_df, xwalk_df, out_dir, write):
    sm_path, xw_path = out_dir / "arm_R_slide_map.tsv", out_dir / "arm_R_xwalk.tsv"
    if write:
        write_tsv(slide_map_df, sm_path)
        write_tsv(xwalk_df, xw_path)
    note = "control re-run: unmodified copy of the current pinned file"
    return [manifest_row(cohort, "R", RESULT_SLIDE_MAP_FIELD, sm_path, 0, [], note, write),
            manifest_row(cohort, "R", RESULT_XWALK_FIELD, xw_path, 0, [], note, write)]


def build_arm_f(cohort, slide_map_df, medium_lookup, oct_only_cases, out_dir, write):
    df, relabelled_slides = relabel_oct_slides(slide_map_df, medium_lookup)
    path = out_dir / "arm_F_slide_map.tsv"
    if write:
        write_tsv(df, path)
    note = (f"{len(relabelled_slides)} OCT Primary Tumor slide row(s) relabelled "
            f"{EXCLUDED_OCT_LABEL}; {len(oct_only_cases)} oct_only patient(s) lose "
            f"every analysed slide")
    return [manifest_row(cohort, "F", RESULT_SLIDE_MAP_FIELD, path,
                         len(oct_only_cases), oct_only_cases, note, write)]


def build_arm_d(cohort, slide_map_df, analysed_cases_sorted, k, s, out_dir, write):
    seed = D_SEED_BASE + s
    drawn = draw_d_arm_cases(analysed_cases_sorted, k, seed)
    df, n_relab = relabel_cases(slide_map_df, drawn)
    arm = f"D{s + 1:02d}"
    path = out_dir / f"arm_{arm}_slide_map.tsv"
    if write:
        write_tsv(df, path)
    note = f"RandomState({seed}).choice(analysed, k={k}, replace=False); {n_relab} PT slide row(s) relabelled"
    return [manifest_row(cohort, arm, RESULT_SLIDE_MAP_FIELD, path, k, drawn, note, write)]


def build_arm_q(cohort, xwalk_df, status_lookup, out_dir, write):
    kept, dropped_cases, dropped_aliquots = filter_disqualified(xwalk_df, status_lookup)
    path = out_dir / "arm_Q_xwalk.tsv"
    if write:
        write_tsv(kept, path)
    shown = ",".join(dropped_aliquots[:10]) + ("..." if len(dropped_aliquots) > 10 else "")
    note = f"{len(dropped_aliquots)} Disqualified aliquot row(s) removed ({shown})"
    return [manifest_row(cohort, "Q", RESULT_XWALK_FIELD, path,
                         len(dropped_cases), dropped_cases, note, write)]


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("--cohort", required=True, choices=COHORTS)
    ap.add_argument("--zw-root", default=ZW_DEFAULT)
    ap.add_argument("--scripts-dir", default=None,
                    help="defaults to this script's own directory")
    ap.add_argument("--morpho-env-sh", default=None,
                    help="defaults to <scripts-dir>/morpho_env.sh")
    ap.add_argument("--out-root", default=None,
                    help="defaults to <zw-root>/scripts/pinned/posthoc")
    ap.add_argument("--slide-map-sha256", required=True,
                    help="expected sha256 of the CURRENT MORPHO_SLIDE_MAP morpho_env.sh resolves")
    ap.add_argument("--xwalk-sha256", required=True,
                    help="expected sha256 of the CURRENT MORPHO_ALIQUOT_XWALK morpho_env.sh resolves")
    ap.add_argument("--slide-embedding-medium", default=None,
                    help="pinned/slide_embedding_medium.tsv; required for cohort in {ucec, luad}")
    ap.add_argument("--slide-embedding-medium-sha256", default=None)
    ap.add_argument("--pdc-aliquot-status", default=None,
                    help="pinned/pdc_aliquot_status.tsv (module D's output); required for cohort in {ucec, luad}")
    ap.add_argument("--pdc-aliquot-status-sha256", default=None)
    ap.add_argument("--dry-run", action="store_true",
                    help="compute and print everything; write nothing under --out-root")
    args = ap.parse_args()

    scripts_dir = args.scripts_dir or os.path.dirname(os.path.abspath(__file__))
    morpho_env_sh = args.morpho_env_sh or os.path.join(scripts_dir, "morpho_env.sh")
    out_root = args.out_root or os.path.join(args.zw_root, "scripts", "pinned", "posthoc")
    if (args.slide_embedding_medium is None or args.slide_embedding_medium_sha256 is None
            or args.pdc_aliquot_status is None or args.pdc_aliquot_status_sha256 is None):
        sys.exit("FATAL: arms F/D/Q need --slide-embedding-medium(+sha256) and "
                  "--pdc-aliquot-status(+sha256). Whether a cohort qualifies for them is "
                  "decided from its own data (PRESPEC M2 and Q), so these are required for "
                  "every cohort -- a cohort that does not qualify gets the M3 disclosure, "
                  "which still has to be computed.")

    print(f"=== posthoc_build_arms: cohort={args.cohort} dry_run={args.dry_run} ===")
    print("--- hash checks ---")
    assert_sha256(os.path.join(scripts_dir, "residual_analysis.py"),
                  EXPECTED_SHA256_RESIDUAL_ANALYSIS, "residual_analysis.py")
    assert_sha256(morpho_env_sh, EXPECTED_SHA256_MORPHO_ENV, "morpho_env.sh")

    print("--- resolving cohort paths via morpho_env.sh (not reimplemented) ---")
    env = resolve_morpho_env(morpho_env_sh, args.cohort)
    for k in ("MORPHO_ROOT", "MORPHO_SLIDE_MAP", "MORPHO_ALIQUOT_XWALK", "MORPHO_WSI_EMB_DIR"):
        print(f"  {k} = {env[k]}")

    slide_map_path, xwalk_path = env["MORPHO_SLIDE_MAP"], env["MORPHO_ALIQUOT_XWALK"]
    assert_sha256(slide_map_path, args.slide_map_sha256, "MORPHO_SLIDE_MAP (current)")
    assert_sha256(xwalk_path, args.xwalk_sha256, "MORPHO_ALIQUOT_XWALK (current)")

    slide_map_df = read_pinned_tsv(slide_map_path)
    xwalk_df = read_pinned_tsv(xwalk_path)
    for col in ("slide_submitter_id", "sample_type", "case_submitter_id"):
        if col not in slide_map_df.columns:
            sys.exit(f"FATAL: {slide_map_path} is missing column {col!r}; has {list(slide_map_df.columns)}.")
    for col in ("aliquot_id", "case_id"):
        if col not in xwalk_df.columns:
            sys.exit(f"FATAL: {xwalk_path} is missing column {col!r}; has {list(xwalk_df.columns)}.")

    print("--- analysed patients (Primary Tumor slides with an embedding file) ---")
    present_ids = present_slide_ids_from_emb_dir(env["MORPHO_WSI_EMB_DIR"])
    pt = slide_map_df[slide_map_df["sample_type"] == "Primary Tumor"]
    case_to_present = defaultdict(list)
    for sid, cid in zip(pt["slide_submitter_id"], pt["case_submitter_id"]):
        if sid in present_ids:
            case_to_present[cid].append(sid)
    analysed_from_slides = sorted(c for c, lst in case_to_present.items() if lst)
    print(f"  (a) filename-scan rollup: {len(analysed_from_slides)} analysed case(s)")

    scratch_out = Path(env["MORPHO_ROOT"]) / "results" / "posthoc" / "_build_arms_scratch"
    analysed_via_ra = compute_analysed_cases_via_ra(env, scripts_dir, scratch_out)
    print(f"  (b) residual_analysis.load_wsi_embeddings() rollup: {len(analysed_via_ra)} analysed case(s)")
    # (b) is the definition, because (b) IS the function every arm's run calls.
    # The two disagree asymmetrically and the two directions mean different things:
    #   (b)-only -- load_wsi_embeddings() analyses a case this script's scan missed.
    #       That is a defect in SLIDE_RE or the sample_type filter here, and no arm
    #       built from this run could be trusted, so it stays FATAL.
    #   (a)-only -- this script's scan sees an embedding file for a case that
    #       load_wsi_embeddings() then drops (it is unmatched against the pinned
    #       slide map, or fails the sample_type filter). Such a case is simply not
    #       analysed; keeping it would put non-analysed patients into the pool the
    #       D arms draw k from, and into the denominators of k and the medium join.
    #       So intersect down to (b) and report what was dropped.
    only_a = sorted(set(analysed_from_slides) - set(analysed_via_ra))
    only_b = sorted(set(analysed_via_ra) - set(analysed_from_slides))
    if only_b:
        sys.exit(f"FATAL: residual_analysis.load_wsi_embeddings() analyses {len(only_b)} "
                 f"case(s) this script's filename scan did not find: {only_b[:10]}. "
                 f"SLIDE_RE or the sample_type filter here no longer mirrors that "
                 f"function; re-read it before trusting any arm built from this run.")
    analysed = sorted(analysed_via_ra)
    if only_a:
        for c in only_a:
            case_to_present.pop(c, None)
        print(f"  note: {len(only_a)} case(s) have a Primary Tumor embedding file but are "
              f"not analysed by load_wsi_embeddings() (unmatched against the pinned slide "
              f"map, or filtered on sample_type), and are excluded here too: "
              f"{only_a[:10]}{' ...' if len(only_a) > 10 else ''}")
    print(f"  analysed patients taken from (b): {len(analysed)} case(s)")

    cohort_out = Path(out_root) / args.cohort
    write = not args.dry_run
    if write:
        cohort_out.mkdir(parents=True, exist_ok=True)

    rows = build_arm_r(args.cohort, slide_map_df, xwalk_df, cohort_out, write)

    # Which arms a cohort gets is decided from its own data, exactly as sections M2
    # and Q of the pre-specification state ("in every cohort with k >= 3 OCT-only
    # patients", "in every cohort with at least one analysed Disqualified aliquot").
    # A cohort that fails a criterion gets the section M3 disclosure naming the
    # criterion, not silence.
    disclosures = []
    print("--- medium / aliquot-status inputs ---")
    assert_sha256(args.slide_embedding_medium, args.slide_embedding_medium_sha256,
                  "slide_embedding_medium.tsv")
    assert_sha256(args.pdc_aliquot_status, args.pdc_aliquot_status_sha256,
                  "pdc_aliquot_status.tsv")
    medium_lookup = load_medium_lookup(args.slide_embedding_medium, args.cohort)
    status_lookup, status_id_col = load_status_lookup(args.pdc_aliquot_status, args.cohort)

    n_xwalk = len(xwalk_df)
    matched = sum(1 for a in xwalk_df["aliquot_id"] if a in status_lookup)
    frac = matched / n_xwalk if n_xwalk else 0.0
    print(f"  aliquot_status join (xwalk.aliquot_id <-> status.{status_id_col}): "
          f"{matched}/{n_xwalk} = {frac:.1%} matched")
    if frac < 0.95:
        unmatched = sorted(set(xwalk_df["aliquot_id"]) - set(status_lookup))[:5]
        sys.exit(f"FATAL: only {frac:.1%} of xwalk aliquots matched {args.pdc_aliquot_status} "
                 f"on column {status_id_col!r} (<95% threshold) -- likely a UUID-vs-submitter-id "
                 f"join-key mismatch (see load_status_lookup's docstring note). "
                 f"Unmatched examples: {unmatched}. Refusing to build arm Q on an unreliable join.")

    # Section M0's own gate, on the same footing as the arm-Q join check above: a
    # wholesale medium join failure would silently drive k to 0 and make the whole
    # M2 comparison vacuous. M0 answers a failed join with a disclosure, not a
    # FATAL, so this skips the medium arms rather than killing arm R and arm Q.
    present_slides = sorted({s for slides in case_to_present.values() for s in slides})
    med_matched = sum(1 for s in present_slides if s in medium_lookup)
    med_frac = med_matched / len(present_slides) if present_slides else 0.0
    print(f"  medium join (present Primary Tumor slides <-> slide_embedding_medium): "
          f"{med_matched}/{len(present_slides)} = {med_frac:.1%} matched")

    classify = classify_patients_by_medium(case_to_present, medium_lookup)
    from collections import Counter
    breakdown = Counter(classify.get(c, "not_analysed") for c in analysed)
    print(f"  patient medium breakdown (of {len(analysed)} analysed): {dict(breakdown)}")
    oct_only_cases = sorted(c for c in analysed if classify.get(c) == "oct_only")
    k = len(oct_only_cases)
    print(f"  k (oct_only analysed patients) = {k}")

    if med_frac < 0.95:
        disclosures.append(
            f"M0/M2: medium arms not built for {args.cohort}: only {med_frac:.1%} of its "
            f"{len(present_slides)} analysed Primary Tumor slides joined to "
            f"slide_embedding_medium.tsv (M0 requires more than 95%). k was computed as "
            f"{k} but is not trustworthy under this join rate.")
    elif k < M2_MIN_OCT_ONLY:
        disclosures.append(
            f"M2: medium arms not built for {args.cohort}: it has k = {k} OCT-only analysed "
            f"patients, below the pre-specified bar of {M2_MIN_OCT_ONLY}. Medium breakdown of "
            f"its {len(analysed)} analysed patients: {dict(breakdown)}.")
    else:
        rows += build_arm_f(args.cohort, slide_map_df, medium_lookup, oct_only_cases,
                            cohort_out, write)
        for s in range(N_D_ARMS):
            rows += build_arm_d(args.cohort, slide_map_df, analysed, k, s, cohort_out, write)

    dq_analysed = sorted({a for a in xwalk_df["aliquot_id"]
                          if str(status_lookup.get(a, "")).strip().lower() == "disqualified"})
    print(f"  analysed aliquots flagged Disqualified at PDC: {len(dq_analysed)}")
    if dq_analysed:
        rows += build_arm_q(args.cohort, xwalk_df, status_lookup, cohort_out, write)
    else:
        disclosures.append(
            f"Q: arm Q not built for {args.cohort}: none of its {n_xwalk} analysed aliquots is "
            f"flagged Disqualified at PDC, so the arm would be identical to arm R.")

    if disclosures:
        print("\n=== DISCLOSURES (pre-specification section M3) ===")
        for d in disclosures:
            print("  " + d)
        if write:
            disc_path = cohort_out / "DISCLOSURES.txt"
            disc_path.write_text("\n".join(disclosures) + "\n", encoding="utf-8", newline="\n")
            print(f"wrote {disc_path}")

    manifest_df = pd.DataFrame(rows)
    print("\n=== MANIFEST ===")
    print(manifest_df.to_string(index=False))
    if write:
        manifest_path = cohort_out / "MANIFEST.tsv"
        write_tsv(manifest_df, manifest_path)
        print(f"\nwrote {manifest_path}")
    else:
        print("\n(dry run: no files or manifest written)")


if __name__ == "__main__":
    main()
