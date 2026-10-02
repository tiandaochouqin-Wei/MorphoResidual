#!/usr/bin/env python3
"""
c1_pull_omics_luad.py -- C1-LUAD confirmatory cohort: pull PDC000489 (TMT11) protein +
GDC STAR-Counts RNA for the 113 confirmatory cases, per the (not-yet-signed) frozen rule
review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 1.5 and Appendix item 2.

THIS SCRIPT DOES NOT RUN AGAINST REAL DATA YET. It is written now, against the draft
rule, so it is ready the moment the rule is signed; it must not be run until then
(section 1.6: RNA/protein are the LAST modality pulled, after slide download, reader QC
and full extraction). Writing it does not download, open, or embed any confirmatory
value -- every number in this file is either a frozen count from already-public
metadata (pinned/*.tsv, pinned/*.txt) or a structural constant (a study id, a column
name, a file name).

Modelled on c1_pull_omics.py (the EXECUTED C1-UCEC pull), but NOT a copy-paste with
cohort strings swapped -- LUAD's crosswalk has to be stricter, for a reason UCEC's
script never had to handle. c1_pull_omics.py's build_aliquot_xwalk() keeps a
biospecimenPerStudy row unless "Normal" is a *substring* of its sample_type. That
filter does NOT exclude PDC's own QC/bridge/pooled pseudo-case rows, whose sample_type
is "Not Reported", not "Normal": re-running that exact filter against PDC000153
(discovery) on 2026-09-29 kept 115 "cases" (111 real + 4 pseudo-case labels --
'Internal Reference - Pooled Sample', 'Normal Only IR', 'Taiwanese IR', 'Tumor Only
IR'), and 3 of those 4 pseudo-case labels turned out to have a real column in that
study's quantDataMatrix -- they would have leaked into a case-indexed protein CSV as
fake "cases" (review/C1_LUAD_STEP0B_2026-09-29/STEP0B_WORKFLOW_RESULT.json, branch
"protein", discrepancy (1); independently re-verified there by a second reviewer pass,
same finding). PDC000489 (the confirmatory study this script actually queries) has its
own 11 such pseudo-case labels over 12 rows (rule section 1.2 item 1, reproduced by
c1_luad_derive_cases.py, EXPECTED_PSEUDO_LABELS below): CR.C1, CR.C1.newly.combined,
CR.C2, CR.D, CR.D.MAM, CR.D.newly.combined, CR.D.pool.2, CR.LSCC, Taiwanese IR, Tumor
Only CONF, Tumor Only DISC. build_aliquot_xwalk() below closes this by requiring
sample_type == "Primary Tumor" EXACTLY (not "no 'Normal' substring"), together with
case_submitter_id in the frozen 113-case set -- both conditions independently exclude
every pseudo-case row, since a pseudo row's sample_type is never "Primary Tumor" and
its case label (e.g. "CR.C1") is never one of the 113 real case ids. KNOWN_PSEUDO_LABELS
and the C3L-/C3N- prefix check are a third, structural belt-and-suspenders check on top
of both (mirrors c1_luad_derive_cases.py's own non_c3 check), not the primary guard.

Also stricter than c1_pull_omics.py in a second way the UCEC pull never needed: FATAL
if any kept confirmatory aliquot_submitter_id also appears in the DISCOVERY crosswalk
(--discovery-xwalk) -- the section 1.5 independence guard. UCEC had no prior-wave bridge
cases to worry about; LUAD does (7 patients excluded at the case level in
c1_luad_derive_cases.py, section 1.2), and this is the aliquot-level backstop on top of
that case-level exclusion, checked fresh against whatever PDC actually returns at pull
time rather than trusted from the Step 0b snapshot.

Like UCEC, the case set is NOT rediscovered here: it is the 113 cases already fixed by
--case-ids (pinned/c1_case_ids_luad.txt). This script only pulls the other two
modalities for that fixed set -- it must not add or drop cases. Unlike c1_pull_omics.py,
--case-ids and --discovery-xwalk have NO default (Appendix item 1's "every path is
required" discipline, applied here because c1_pull_omics.py's own defaulted, UCEC-
specific --case-ids/--out-dir were flagged by the Step 0b review as exactly the kind of
stale-default risk this project already got burned by once, 2026-09-15, UCEC output
silently landing in pdac/results).

quantDataMatrix's column-header aliquot format is NOT assumed to match
cptac2_replicate.py's TCGA-barcode parsing (see c1_pull_omics.py's own header comment);
it is printed raw and joined against the aliquot->case crosswalk built here. If fewer
than half the crosswalk's aliquots match a matrix column, the script stops rather than
silently emitting a near-empty matrix -- this guards against a wrong header-format
assumption, a different failure mode than low case coverage (which the coverage table
records, not FATALs on).

Usage:
  python c1_pull_omics_luad.py protein [--pdc-study PDC000489] [--data-type log2_ratio] \
      --case-ids <path> --discovery-xwalk <path> --out-dir <path>
  python c1_pull_omics_luad.py rna     --case-ids <path> --discovery-xwalk <path> --out-dir <path>
  python c1_pull_omics_luad.py both    --case-ids <path> --discovery-xwalk <path> --out-dir <path>

  --case-ids         pinned/c1_case_ids_luad.txt (113 ids). REQUIRED, no default.
  --discovery-xwalk  aliquot_to_case_tumor_luad.tsv (the DISCOVERY crosswalk, columns
                      aliquot_id/case_id/sample_type -- the same file c1_run_test_luad.py
                      hash-checks as "discovery_xwalk"). REQUIRED, no default. Used only
                      for the independence guard; never merged into the confirmatory
                      output.
  --out-dir          REQUIRED, no default. Recommended HPC path (unchecked here, this
                      script does not create cohort directories on its own):
                      /public/home/fjhui/ZW/luad_c1

Outputs (relative to --out-dir):
  aliquot_to_case_tumor_luad_c1.tsv   confirmatory aliquot->case crosswalk actually used
                                       (Primary Tumor aliquots of the 113 frozen cases
                                       only -- for audit, not for reuse as a discovery
                                       crosswalk)
  protein_luad_c1.csv                 case-indexed protein matrix (index_col=0 on read),
                                       matching exactly what
                                       c1_run_test_luad.load_confirmatory_protein expects
  protein_luad_c1_gene_nan_counts.tsv per-gene NaN counts (section 1.5: "The pull writes
                                       a per-case coverage table and per-gene NaN
                                       counts")
  manifest_rna_tumor_luad_c1.tsv      RNA manifest in residual_analysis.py's expected
                                       tumor-manifest format (file_id/file_name/
                                       case_submitter_id/sample_type)
  gdc_rna_download_manifest_luad_c1.tsv   gdc_download.py-compatible download manifest
  coverage_luad_c1.tsv                per-case coverage: has_protein / has_rna / both,
                                       for all 113 frozen cases (section 1.6 step 5:
                                       "the coverage table"). Running `protein` then
                                       `rna` (or vice versa) in two separate invocations
                                       merges into the same file rather than clobbering
                                       the other modality's column.

This script does not import residual_analysis, torch, numpy or any HPC/GPU-only
module -- only argparse/hashlib/json/os/sys/time/urllib + pandas -- so it can be syntax-
checked and its pure filtering logic unit-tested (see build_aliquot_xwalk) on any
machine, including one with no torch and no network access.
"""
import argparse
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from collections import Counter

import pandas as pd

PDC_URL = "https://pdc.cancer.gov/graphql"
GDC_FILES = "https://api.gdc.cancer.gov/files"

CONFIRMATORY_STUDY_DEFAULT = "PDC000489"
# The DISCOVERY study id. This script pulls the CONFIRMATORY study only; --pdc-study
# is checked against this at startup so a copy-paste mistake (or a bad default someone
# adds later) cannot point this script at the study c1_pull_discovery scripts already
# cover. Real protection, not decorative: checked against the actual --pdc-study value
# used, before any query is issued (contrast the vacuous "_source_study" assert flagged
# in protein_route_check.py's review, STEP0B_WORKFLOW_RESULT.json branch "protein").
FORBIDDEN_DISCOVERY_STUDY = "PDC000153"

N_CASE_IDS_EXPECTED = 113

# The 11 QC/bridge/pooled pseudo-case labels in PDC000489's biospecimenPerStudy
# (sample_type = "Not Reported"; rule section 1.2 item 1, reproduced exactly by
# c1_luad_derive_cases.py's EXPECTED_PSEUDO_LABELS). The PRIMARY guard against these is
# structural: sample_type == "Primary Tumor" and case_submitter_id in the frozen 113
# (build_aliquot_xwalk, below) -- a pseudo row satisfies neither. This literal set is a
# second, independent check: if one of these exact labels is ever seen with
# sample_type == "Primary Tumor" (i.e. PDC's own labelling has changed underneath this
# rule), pull_protein() FATALs rather than silently trusting the structural filter
# alone.
KNOWN_PSEUDO_CASE_LABELS = frozenset({
    "CR.C1", "CR.C1.newly.combined", "CR.C2", "CR.D", "CR.D.MAM",
    "CR.D.newly.combined", "CR.D.pool.2", "CR.LSCC", "Taiwanese IR",
    "Tumor Only CONF", "Tumor Only DISC",
})

# section 7.2 / Appendix item 1 discipline: hash-assert the frozen inputs this script
# reads before trusting anything in them. --discovery-xwalk's HPC copy may legitimately
# live at a path the local mirror does not (c1_run_test_luad.py's own note on
# discovery_xwalk/discovery_slide_map), so it honours the same skip env var name for a
# consistent operator experience across both scripts.
EXPECTED_SHA256 = {
    "case_ids": "3629b243a3aa21057dae24e7f5c1b47e8a2571e6ffbcff6f3170333e1ca9fe44",
    "discovery_xwalk": "967cabf18cb6b6c3646aaf05005c46fcdccd7b55eb8dc9882348e84ca4baa2bf",
}

COVERAGE_FILENAME = "coverage_luad_c1.tsv"


def sha256(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return "unreadable"


def assert_hash(path, key, log, allow_skip_env=None):
    """FATAL unless sha256(path) == EXPECTED_SHA256[key], mirroring
    c1_run_test_luad.assert_hash exactly (same message shape, same skip-env
    mechanism) so an operator who knows one script's hash-mismatch behaviour already
    knows this one's."""
    got = sha256(path)
    want = EXPECTED_SHA256[key]
    if got == want:
        log(f"  hash OK: {key} = {got[:12]}... ({path})")
        return got
    msg = (f"sha256 mismatch for {key}: got={got} want={want} ({path}). This is the "
           f"wrong file vintage -- stopping rather than pulling omics against an "
           f"unverified case list or crosswalk.")
    if allow_skip_env and os.environ.get(allow_skip_env):
        log(f"  WARNING (skip forced via {allow_skip_env}): {msg}")
        return got
    sys.exit(f"FATAL: {msg}")


def pdc_query(q, tries=6):
    req = urllib.request.Request(
        PDC_URL, data=json.dumps({"query": q}).encode(),
        headers={"Content-Type": "application/json"})
    for a in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                d = json.load(r)
            if d.get("errors"):
                sys.exit(f"PDC graphql error: {str(d['errors'])[:300]}")
            return d["data"]
        except urllib.error.URLError as e:
            if a == tries - 1:
                raise
            print(f"  PDC attempt {a + 1}/{tries} failed: {e}", flush=True)
            time.sleep(4 * (a + 1))


def gdc_post(url, payload, tries=6):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
    for a in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except urllib.error.URLError as e:
            if a == tries - 1:
                raise
            print(f"  GDC attempt {a + 1}/{tries} failed: {e}", flush=True)
            time.sleep(4 * (a + 1))


def load_case_ids(path, log):
    """The frozen 113. NOT rediscovered -- see module docstring. Hash-checked and
    count-checked; either mismatch is FATAL (section 1.2: "The case set is not
    re-filtered or expanded after signature")."""
    assert_hash(path, "case_ids", log)
    with open(path, encoding="utf-8") as f:
        ids = sorted({line.strip() for line in f if line.strip()})
    if len(ids) != N_CASE_IDS_EXPECTED:
        sys.exit(f"FATAL: {path} has {len(ids)} case ids, expected "
                 f"{N_CASE_IDS_EXPECTED}. The case set is fixed by file, not "
                 f"rederived here -- stopping rather than pulling omics for the wrong "
                 f"set.")
    return ids


def load_discovery_aliquots(path, log):
    """Loads the DISCOVERY aliquot->case crosswalk's aliquot id column, for the
    independence guard (section 1.5: "FATAL if any kept aliquot_submitter_id appears
    in the discovery crosswalk"). Read-only; never merged into any confirmatory
    output."""
    assert_hash(path, "discovery_xwalk", log, allow_skip_env="MORPHO_ALLOW_HPC_XWALK_HASH")
    df = pd.read_csv(path, sep="\t", dtype=str)
    col = "aliquot_id" if "aliquot_id" in df.columns else "aliquot_submitter_id"
    aliquots = set(df[col].dropna())
    log(f"  discovery crosswalk: {len(aliquots)} aliquots ({path})")
    return aliquots


def fetch_confirmatory_biospecimen(pdc_study_id):
    """I/O only: queries PDC for biospecimenPerStudy(pdc_study_id). Kept separate from
    build_aliquot_xwalk() below so the filtering/leak-detection logic can be unit-
    tested with a synthetic, PDC-shaped `rows` list and zero network access."""
    q = ('{ biospecimenPerStudy(pdc_study_id: "%s" acceptDUA: true) '
         '{ aliquot_submitter_id case_submitter_id sample_type } }' % pdc_study_id)
    rows = pdc_query(q)["biospecimenPerStudy"] or []
    return rows


def build_aliquot_xwalk(rows, case_ids, discovery_aliquots, log=print):
    """PURE (no I/O, no sys.exit, no dependence on any file on disk): filters a PDC
    biospecimenPerStudy-shaped `rows` list (each item a dict with
    aliquot_submitter_id/case_submitter_id/sample_type, exactly what
    pdc_query(...)["biospecimenPerStudy"] returns) to Primary Tumor aliquots of the
    frozen `case_ids`, per rule section 1.5:

      "Crosswalk: restricted, before any averaging, to Primary Tumor aliquots of the
      113 frozen cases. Pseudo-case rows and the seven excluded patients never enter."

    Unlike c1_pull_omics.build_aliquot_xwalk (UCEC), which keeps a row unless
    "Normal" is a *substring* of sample_type -- a filter that does NOT exclude PDC's
    "Not Reported" QC/bridge/pooled pseudo-case rows, 3 of which were found to carry a
    real quantDataMatrix column on PDC000153 (module docstring; STEP0B_WORKFLOW_
    RESULT.json, branch "protein", discrepancy (1)) -- this function requires
    sample_type == "Primary Tumor" EXACTLY, and separately requires
    case_submitter_id to be one of the frozen `case_ids`. Both conditions
    independently exclude every pseudo-case row: a pseudo row's sample_type is never
    "Primary Tumor", and its case label (e.g. "CR.C1") is never in `case_ids`, which
    holds only real C3L-/C3N- ids. KNOWN_PSEUDO_CASE_LABELS and the C3L-/C3N- prefix
    check below are a third, structural belt-and-suspenders check on top of both.

    Also checks the independence guard here (structurally, not with sys.exit -- the
    caller decides what to do with a non-empty `leaked_aliquots`): a kept aliquot id
    that also appears in `discovery_aliquots` is recorded, never silently dropped or
    silently kept.

    Returns (xwalk, diag).
      xwalk: {aliquot_submitter_id: case_submitter_id}, one entry per kept Primary
             Tumor aliquot. A case with >1 such aliquot gets >1 entry, all pointing to
             that same case (diag['multi_aliquot_cases'] records which; not an error,
             matches discovery/UCEC "average >1 tumour aliquot per case").
      diag:  dict for logging and FATAL decisions in the caller:
               sample_type_tally        {sample_type: row count}, all rows
               n_rows                   len(rows)
               n_kept_rows              rows passing both filters
               n_kept_cases             distinct cases in xwalk
               multi_aliquot_cases      {case: [aliquot, ...]} where len > 1
               conflicts                [(aliquot, case_a, case_b), ...] -- SAME
                                        aliquot id claimed by two DIFFERENT kept
                                        cases; a real data-integrity break, always
                                        FATAL in the caller if non-empty
               pseudo_labels_seen       sorted list of KNOWN_PSEUDO_CASE_LABELS that
                                        appeared with sample_type == "Primary Tumor"
                                        in the raw rows (should always be empty)
               non_c3_case_ids          sorted list of kept case ids that are not
                                        C3L-/C3N- prefixed (should always be empty)
               leaked_aliquots          sorted list of kept aliquot ids also present
                                        in `discovery_aliquots` (should always be
                                        empty; non-empty is the section 1.5 FATAL)
    """
    case_set = set(case_ids)
    tally = Counter(r.get("sample_type") for r in rows)

    kept_rows = [r for r in rows
                 if (r.get("sample_type") or "").strip() == "Primary Tumor"
                 and r.get("case_submitter_id") in case_set]

    pseudo_labels_seen = sorted({
        r.get("case_submitter_id") for r in rows
        if r.get("case_submitter_id") in KNOWN_PSEUDO_CASE_LABELS
        and (r.get("sample_type") or "").strip() == "Primary Tumor"
    })

    xwalk = {}
    conflicts = []
    case_aliquots = {}
    for r in kept_rows:
        a = r.get("aliquot_submitter_id")
        c = r.get("case_submitter_id")
        if not a:
            continue
        if a in xwalk and xwalk[a] != c:
            conflicts.append((a, xwalk[a], c))
            continue
        xwalk[a] = c
        case_aliquots.setdefault(c, []).append(a)
    multi_aliquot_cases = {c: sorted(set(al)) for c, al in case_aliquots.items()
                            if len(set(al)) > 1}

    non_c3_case_ids = sorted({c for c in xwalk.values()
                              if not (c.startswith("C3L-") or c.startswith("C3N-"))})

    leaked_aliquots = sorted(set(xwalk) & set(discovery_aliquots))

    diag = dict(
        sample_type_tally=dict(tally),
        n_rows=len(rows),
        n_kept_rows=len(kept_rows),
        n_kept_cases=len(set(xwalk.values())),
        multi_aliquot_cases=multi_aliquot_cases,
        conflicts=conflicts,
        pseudo_labels_seen=pseudo_labels_seen,
        non_c3_case_ids=non_c3_case_ids,
        leaked_aliquots=leaked_aliquots,
    )
    log(f"[xwalk] {len(rows)} biospecimen rows; sample_type tally: {dict(tally)}")
    log(f"[xwalk] {len(kept_rows)} Primary Tumor rows restricted to the frozen "
        f"{len(case_set)}-case set -> {len(xwalk)} aliquots covering "
        f"{diag['n_kept_cases']}/{len(case_set)} frozen cases")
    if multi_aliquot_cases:
        log(f"  {len(multi_aliquot_cases)} case(s) carry >1 Primary Tumor aliquot "
            f"(averaged downstream, matches discovery convention): "
            f"{multi_aliquot_cases}")
    return xwalk, diag


def pull_protein(args, log):
    case_ids = load_case_ids(args.case_ids, log)
    discovery_aliquots = load_discovery_aliquots(args.discovery_xwalk, log)
    log(f"{len(case_ids)} frozen confirmatory case ids loaded from {args.case_ids}")

    if args.pdc_study == FORBIDDEN_DISCOVERY_STUDY:
        sys.exit(f"FATAL: --pdc-study is the DISCOVERY study ({FORBIDDEN_DISCOVERY_STUDY}). "
                 f"This script pulls the CONFIRMATORY study only; refusing to run.")

    log(f"\n[protein] fetching biospecimenPerStudy({args.pdc_study}) ...")
    rows = fetch_confirmatory_biospecimen(args.pdc_study)

    xwalk, diag = build_aliquot_xwalk(rows, case_ids, discovery_aliquots, log)

    if diag["conflicts"]:
        sys.exit(f"FATAL: {len(diag['conflicts'])} aliquot_submitter_id(s) claimed by "
                 f"two different kept cases -- a real data-integrity break, cannot "
                 f"proceed: {diag['conflicts'][:5]}")
    if diag["pseudo_labels_seen"]:
        sys.exit(f"FATAL: {len(diag['pseudo_labels_seen'])} KNOWN pseudo-case "
                 f"label(s) carry sample_type == 'Primary Tumor' in the raw PDC "
                 f"response -- {diag['pseudo_labels_seen']}. PDC's own labelling has "
                 f"apparently changed since the rule was drafted; this needs "
                 f"corresponding-author review before proceeding, not a silent "
                 f"pass-through.")
    if diag["non_c3_case_ids"]:
        sys.exit(f"FATAL: {len(diag['non_c3_case_ids'])} kept case id(s) are not "
                 f"C3L-/C3N- CPTAC ids -- {diag['non_c3_case_ids']}. Refusing: this "
                 f"would mean a pseudo-case or malformed label passed every other "
                 f"filter.")
    if diag["leaked_aliquots"]:
        sys.exit(f"FATAL (independence guard, section 1.5): "
                 f"{len(diag['leaked_aliquots'])} kept confirmatory "
                 f"aliquot_submitter_id(s) also appear in the discovery crosswalk "
                 f"({args.discovery_xwalk}) -- {diag['leaked_aliquots'][:5]}. This "
                 f"would mean a discovery aliquot re-entered the confirmatory pull; "
                 f"refusing to build a protein matrix from it.")
    log(f"[protein] independence guard OK: 0 of {len(xwalk)} kept aliquots overlap "
        f"the {len(discovery_aliquots)}-aliquot discovery crosswalk")

    xwalk_path = os.path.join(args.out_dir, "aliquot_to_case_tumor_luad_c1.tsv")
    with open(xwalk_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("aliquot_id\tcase_id\tsample_type\n")
        for a, c in sorted(xwalk.items()):
            f.write(f"{a}\t{c}\tPrimary Tumor\n")
    log(f"[protein] wrote {len(xwalk)} rows -> {xwalk_path}")

    q = '{ quantDataMatrix(pdc_study_id: "%s", data_type: "%s") }' % (
        args.pdc_study, args.data_type)
    mat = pdc_query(q)["quantDataMatrix"]
    header = mat[0]
    log(f"[pdc] quantDataMatrix {args.pdc_study}/{args.data_type}: "
        f"{len(mat) - 1} genes x {len(header) - 1} columns")
    log(f"  raw header sample (first 5): {header[1:6]}")

    # section 1.5: "FATAL on duplicate gene symbols." Checked on the RAW response,
    # before any dict/DataFrame construction could silently let a later duplicate
    # overwrite an earlier one.
    gene_counts = Counter(r[0] for r in mat[1:])
    dup_genes = sorted(g for g, n in gene_counts.items() if n > 1)
    if dup_genes:
        sys.exit(f"FATAL: {len(dup_genes)} duplicate gene symbol(s) in the "
                 f"quantDataMatrix response -- {dup_genes[:10]}. A case-indexed "
                 f"matrix keyed on gene symbol would silently overwrite one row's "
                 f"values with another's; refusing to build it.")
    log(f"  0 duplicate gene symbols among {len(gene_counts)} genes")

    # observed header format (verified 2026-09-16 against PDC000439, and again against
    # PDC000153 in Step 0b, branch "protein"): "<aliquot_uuid>:<aliquot_submitter_id>"
    # -- take the part after the colon.
    matched = sum(1 for c in header[1:] if c.split(":")[-1] in xwalk)
    n_aliquots = len(xwalk)
    log(f"  {matched}/{len(header) - 1} matrix columns match a kept confirmatory "
        f"aliquot ({n_aliquots} aliquots in the crosswalk)")
    if n_aliquots and matched < n_aliquots / 2:
        sys.exit(
            "FATAL: fewer than half the crosswalk's aliquots matched a matrix column "
            "-- the header format assumption (raw aliquot_submitter_id after ':') is "
            "probably wrong for this study. Inspect the raw header above before "
            "adapting the parsing; do not guess a TCGA-barcode transform (that is for "
            "the CPTAC-2 legacy cohorts only).")

    cols, keep_idx = [], []
    for i, c in enumerate(header[1:]):
        case = xwalk.get(c.split(":")[-1])
        if case:
            cols.append(case)
            keep_idx.append(i)

    rows_by_gene = {}
    for r in mat[1:]:
        gene = r[0]
        vals = r[1:]
        rows_by_gene[gene] = [vals[i] for i in keep_idx]
    df = pd.DataFrame(rows_by_gene, index=cols).apply(pd.to_numeric, errors="coerce")
    df = df.groupby(level=0).mean()  # >1 kept aliquot per case -> average, matches discovery

    case_set = set(case_ids)
    covered = set(df.index) & case_set
    log(f"[protein] {df.shape[0]} cases x {df.shape[1]} genes; "
        f"{len(covered)}/{len(case_set)} of the frozen 113 present")

    out_path = os.path.join(args.out_dir, "protein_luad_c1.csv")
    df.to_csv(out_path)
    log(f"[protein] wrote -> {out_path}")

    # section 1.5: "The pull writes a per-case coverage table and per-gene NaN counts."
    nan_path = os.path.join(args.out_dir, "protein_luad_c1_gene_nan_counts.tsv")
    n_cases = df.shape[0]
    with open(nan_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("gene\tn_present\tn_nan\tfrac_nan\n")
        for gene in df.columns:
            n_nan = int(df[gene].isna().sum())
            n_present = n_cases - n_nan
            frac = (n_nan / n_cases) if n_cases else float("nan")
            f.write(f"{gene}\t{n_present}\t{n_nan}\t{frac:.6f}\n")
    log(f"[protein] wrote per-gene NaN counts ({df.shape[1]} genes) -> {nan_path}")

    cov = load_coverage(args.out_dir, case_ids)
    for c in case_ids:
        cov[c]["has_protein"] = "True" if c in covered else "False"
    write_coverage(args.out_dir, case_ids, cov)
    return covered


def pull_rna(args, log):
    case_ids = load_case_ids(args.case_ids, log)
    log(f"{len(case_ids)} frozen confirmatory case ids loaded from {args.case_ids}")

    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.submitter_id", "value": case_ids}},
            {"op": "in", "content": {"field": "data_type",
                                      "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type",
                                      "value": ["STAR - Counts"]}},
        ],
    }
    fields = "file_id,file_name,md5sum,file_size,cases.submitter_id,cases.samples.sample_type"
    payload = {"filters": filters, "size": 5000, "format": "json", "fields": fields}
    hits = gdc_post(GDC_FILES, payload)["data"]["hits"]
    log(f"[gdc] {len(hits)} STAR-Counts files across the {len(case_ids)} frozen "
        f"confirmatory cases (before sample_type filtering)")

    kept, dropped_other, dropped_nocase = [], 0, 0
    manifest_path = os.path.join(args.out_dir, "manifest_rna_tumor_luad_c1.tsv")
    dl_manifest_path = os.path.join(args.out_dir, "gdc_rna_download_manifest_luad_c1.tsv")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fm, \
         open(dl_manifest_path, "w", encoding="utf-8", newline="\n") as fd:
        fm.write("file_id\tfile_name\tcase_submitter_id\tsample_type\n")
        fd.write("id\tfilename\tmd5\tsize\tstate\n")
        for h in hits:
            cs = (h.get("cases") or [{}])[0]
            case = cs.get("submitter_id")
            stype = (cs.get("samples", [{}])[0] or {}).get("sample_type", "")
            if not case:
                dropped_nocase += 1
                continue
            # exact match, not just "no 'Normal' substring" -- matches
            # residual_analysis.load_rna_matrix()'s own defensive filter
            # (manifest["sample_type"] == "Primary Tumor") and this rule's discipline
            # throughout (section 1.5: "RNA. GDC STAR-Counts, tumour only").
            if stype != "Primary Tumor":
                dropped_other += 1
                continue
            fm.write(f"{h['file_id']}\t{h['file_name']}\t{case}\t{stype}\n")
            fd.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum', '')}\t"
                     f"{h.get('file_size', 0)}\treleased\n")
            kept.append(case)

    covered = set(kept) & set(case_ids)
    log(f"[rna] kept {len(kept)} Primary Tumor files (dropped {dropped_other} "
        f"non-Primary-Tumor, {dropped_nocase} no-case), covering "
        f"{len(covered)}/{len(case_ids)} frozen cases")
    dup = len(kept) - len(set(kept))
    if dup:
        log(f"  note: {dup} case(s) have >1 Primary Tumor STAR file (multiple RNA "
            f"aliquots); load_rna_matrix() averages these per gene, matching discovery")
    missing = set(case_ids) - covered
    if missing:
        log(f"  {len(missing)} frozen cases have no Primary Tumor STAR-Counts file: "
            f"{sorted(missing)[:10]}{' ...' if len(missing) > 10 else ''}")
    log(f"[rna] wrote {manifest_path}\n[rna] wrote {dl_manifest_path}")
    log(f"\nNEXT: nohup python gdc_download.py {dl_manifest_path} "
        f"{args.out_dir}/omics/rna 6 > {args.out_dir}/rna_dl.log 2>&1 &")

    cov = load_coverage(args.out_dir, case_ids)
    for c in case_ids:
        cov[c]["has_rna"] = "True" if c in covered else "False"
    write_coverage(args.out_dir, case_ids, cov)
    return covered


def load_coverage(out_dir, case_ids):
    """Merge-safe read: if protein and rna are pulled in two separate invocations
    (the common case -- rna via GDC is usually run separately from the protein pull),
    this preserves whichever column the previous run already wrote instead of
    clobbering it with 'unknown'."""
    path = os.path.join(out_dir, COVERAGE_FILENAME)
    cov = {c: {"has_protein": "unknown", "has_rna": "unknown"} for c in case_ids}
    if os.path.exists(path):
        try:
            existing = pd.read_csv(path, sep="\t", dtype=str)
        except Exception as e:
            print(f"  WARNING: could not read existing {path} ({e}); starting a "
                  f"fresh coverage table", flush=True)
            return cov
        for _, r in existing.iterrows():
            c = r.get("case_submitter_id")
            if c in cov:
                if pd.notna(r.get("has_protein")):
                    cov[c]["has_protein"] = r["has_protein"]
                if pd.notna(r.get("has_rna")):
                    cov[c]["has_rna"] = r["has_rna"]
    return cov


def write_coverage(out_dir, case_ids, cov):
    """section 1.6 step 5: 'the coverage table.' One row per frozen case, always, in
    frozen sorted order, whatever the coverage."""
    path = os.path.join(out_dir, COVERAGE_FILENAME)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("case_submitter_id\thas_protein\thas_rna\tboth\n")
        for c in case_ids:
            hp, hr = cov[c]["has_protein"], cov[c]["has_rna"]
            both = "True" if (hp == "True" and hr == "True") else "False"
            f.write(f"{c}\t{hp}\t{hr}\t{both}\n")
    n_both = sum(1 for c in case_ids if cov[c]["has_protein"] == "True"
                 and cov[c]["has_rna"] == "True")
    print(f"  coverage: {n_both}/{len(case_ids)} cases have both modalities so far "
          f"-> {path}", flush=True)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["protein", "rna", "both"])
    ap.add_argument("--pdc-study", default=CONFIRMATORY_STUDY_DEFAULT,
                    help=f"frozen confirmatory study id (rule section 1.1); the "
                         f"default ({CONFIRMATORY_STUDY_DEFAULT}) is the frozen value "
                         f"and should not normally be overridden")
    ap.add_argument("--data-type", default="log2_ratio")
    ap.add_argument("--case-ids", required=True,
                    help="pinned/c1_case_ids_luad.txt, the frozen 113. REQUIRED, no "
                         "default (Appendix item 1 discipline).")
    ap.add_argument("--discovery-xwalk", required=True,
                    help="the DISCOVERY aliquot->case crosswalk "
                         "(aliquot_to_case_tumor_luad.tsv), for the section 1.5 "
                         "independence guard. REQUIRED, no default.")
    ap.add_argument("--out-dir", required=True,
                    help="output directory (created if missing). REQUIRED, no "
                         "default.")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    def log(msg):
        print(msg, flush=True)

    log("=== c1_pull_omics_luad.py ===")
    log(f"  mode             : {args.mode}")
    log(f"  pdc-study        : {args.pdc_study}")
    log(f"  data-type        : {args.data_type}")
    log(f"  case-ids         : {args.case_ids}")
    log(f"  discovery-xwalk  : {args.discovery_xwalk}")
    log(f"  out-dir          : {args.out_dir}")
    log("=" * 78)

    if args.mode in ("protein", "both"):
        pull_protein(args, log)
    if args.mode in ("rna", "both"):
        pull_rna(args, log)

    log(f"\ncoverage table -> {os.path.join(args.out_dir, COVERAGE_FILENAME)}")


if __name__ == "__main__":
    main()
