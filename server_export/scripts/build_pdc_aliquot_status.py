#!/usr/bin/env python3
"""
build_pdc_aliquot_status.py -- PDC aliquot status table for the Disqualified-aliquot
sensitivity (review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md, item "Q. PDC-Disqualified
aliquots (sub-item c)"). POST HOC and reported whatever the result (prespec section 0):
nothing here changes a published count, it only states which analysed aliquots PDC flags
as Disqualified so the Q/R arm comparison in the prespec can be run.

PURPOSE
  1. Build one row per Primary Tumor aliquot, across the five discovery studies
     (ccrcc/luad/ucec/gbm/pdac), from the PDC biospecimenPerStudy responses already
     saved to review/BATCH_AXES_SCOPING_2026-10-02/raw/ by fetch_biospecimen.py on
     2026-10-02 (see that script and build_plex_maps.py for how they were fetched and
     parsed). This script makes NO network requests; it only reads those saved files.
  2. Cross-check the local tumour-aliquot crosswalks
     (server_export/scripts/aliquot_to_case_tumor*.tsv -- "which aliquot(s) each
     analysed case actually uses", per build_aliquot_xwalk_v2.py / morpho_env.sh
     MORPHO_ALIQUOT_XWALK) against that status table, and report, per cohort:
       - analysed cases (figures/figdata/scores_<cohort>.csv index) whose crosswalk
         aliquot is Disqualified;
       - analysed cases with more than one crosswalk row (more than one analysed tumour
         aliquot).
     This confirms or corrects the scoping finding quoted in the prespec: "UCEC:
     C3L-00157, C3L-00356, C3L-00938, C3L-01247, C3L-01253 Disqualified, and C3N-01825
     with a second, Disqualified aliquot; LUAD: C3N-00545", and extends the same check
     to CCRCC, GBM and PDAC, which the scoping had not covered.

INPUTS (all local, read-only; every one is sha256-asserted below -- FATAL on mismatch)
  review/BATCH_AXES_SCOPING_2026-10-02/raw/r3_biospecimen_<PDC_ID>.json
      One per study (PDC000127 ccrcc, PDC000153 luad, PDC000125 ucec, PDC000204 gbm,
      PDC000270 pdac). Each is the raw GraphQL response of
      `{ biospecimenPerStudy(pdc_study_id: "<PDC_ID>" acceptDUA: true)
         { aliquot_id aliquot_submitter_id case_submitter_id sample_type
           aliquot_status } }`, saved verbatim by fetch_biospecimen.py. Not refetched.
  server_export/scripts/aliquot_to_case_tumor[_<cohort>].tsv
      The "local (shared-dir) crosswalk copies" build_plex_maps.py also reads (its
      docstring), columns aliquot_id, case_id, sample_type, built by
      build_aliquot_xwalk_v2.py with KEEP={"Primary Tumor"}. This is the join key for
      "which aliquot(s) an analysed case uses" -- NOT the verified pinned/ copy that
      morpho_env.sh's MORPHO_ALIQUOT_XWALK points HPC jobs at (scripts/pinned/
      aliquot_to_case_tumor<suffix>.tsv), which is HPC-only and not mirrored locally;
      see "WHAT I COULD NOT FULLY SETTLE" below.
  figures/figdata/scores_<cohort>.csv
      Index column = the analysed discovery cases for that cohort (same file
      build_plex_maps.py's analysed_cases() reads).

OUTPUT
  server_export/scripts/pinned/pdc_aliquot_status.tsv
      One row per Primary Tumor aliquot found in the five raw biospecimen snapshots,
      columns: study, case_submitter_id, aliquot_submitter_id, aliquot_status,
      sample_type ("Primary Tumor" on every row by construction). `study` is the short
      cohort key (ccrcc/luad/ucec/gbm/pdac) used throughout this project's --cohort
      arguments (residual_analysis_sitepack.py, lsf_sitepack.sh, etc.), not the PDC
      accession; the cohort -> PDC_ID map is STUDIES below. Rows are sorted by
      (study, case_submitter_id, aliquot_submitter_id) for a deterministic diff.
  Per-cohort cross-check summary (stdout only; not written to a file): disqualified
  analysed-case list, multi-aliquot analysed-case list, and the table's own sha256.

Randomness: none. This script only joins and reformats already-fixed tables; there is
no permutation, bootstrap or sampling step, so no RandomState is used anywhere in it.

WHAT I COULD NOT FULLY SETTLE
  - The local server_export/scripts/aliquot_to_case_tumor_pdac.tsv this script is told
    to treat as "which aliquots each analysed case actually uses" carries 14 rows that
    are not real patient Primary-Tumor aliquots (6 Cell Lines, 5 "Tumor", 3 Not
    Reported -- QC2/QC3/QC6, KoreanReference1/2/3, two pooled-sample rows, and two Not
    Reported rows whose case_id IS a real-looking case id, C3L-03350 and C3N-01900).
    This matches exactly the Sep-2026 contamination event build_aliquot_xwalk_v2.py's
    own docstring and morpho_env.sh describe (the shared scripts/ dir vs. the verified
    scripts/pinned/ copy the HPC job actually reads). Neither C3L-03350 nor C3N-01900
    is an analysed PDAC case (figures/figdata/scores_pdac.csv), so this contamination
    does not change the Disqualified/multi-aliquot verdict for any analysed PDAC case
    -- but it does mean this local file is not provably identical to the verified
    scripts/pinned/aliquot_to_case_tumor_pdac.tsv the HPC job reads, which is not
    mirrored to this machine and so cannot be checked directly. Likewise UCEC's local
    crosswalk carries one harmless non-patient row ("Ref"/"Ref", Not Reported), also
    not an analysed case. Both are reported but not corrected here (not this module's
    file to edit).
  - Whether aliquot_status can change over time at PDC (i.e. whether an aliquot flagged
    Disqualified today was already Disqualified when the discovery run that produced
    scores_<cohort>.csv executed). This script states PDC's CURRENT flag (as fetched
    2026-10-02) against a crosswalk file with no fetch timestamp of its own; it is the
    same status data the scoping already relied on for the UCEC/LUAD finding, so this
    is a pre-existing limitation of that finding, not one this script introduces.
"""
import csv
import hashlib
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
RAW_DIR = os.path.join(PAPER, "review", "BATCH_AXES_SCOPING_2026-10-02", "raw")
FIGDATA_DIR = os.path.join(PAPER, "figures", "figdata")
OUT_PATH = os.path.join(HERE, "pinned", "pdc_aliquot_status.tsv")

# cohort -> (PDC study id, raw biospecimen JSON filename, local tumour-crosswalk filename)
STUDIES = {
    "ccrcc": ("PDC000127", "r3_biospecimen_PDC000127.json", "aliquot_to_case_tumor.tsv"),
    "luad":  ("PDC000153", "r3_biospecimen_PDC000153.json", "aliquot_to_case_tumor_luad.tsv"),
    "ucec":  ("PDC000125", "r3_biospecimen_PDC000125.json", "aliquot_to_case_tumor_ucec.tsv"),
    "gbm":   ("PDC000204", "r3_biospecimen_PDC000204.json", "aliquot_to_case_tumor_gbm.tsv"),
    "pdac":  ("PDC000270", "r3_biospecimen_PDC000270.json", "aliquot_to_case_tumor_pdac.tsv"),
}
COHORT_ORDER = ["ccrcc", "luad", "ucec", "gbm", "pdac"]

# Frozen sha256 of every input this script reads (2026-10-02 snapshot). A mismatch is
# FATAL: it means the raw PDC response, the local crosswalk, or the analysed-case list
# has changed vintage since this script and its report were written, and the cohort's
# cross-check below would silently be run against a different set of facts.
EXPECTED_SHA256 = {
    "biospecimen_ccrcc": "aadddf369941e06561b02566cf97133ce405cdd904683904a6d60bcc52751a9b",
    "biospecimen_luad":  "acc8b5fa18f179d873c1025304048575cfe7659f617da86420e3e731986428f7",
    "biospecimen_ucec":  "7ba43253498c5c611fac795a9f6762c3a10589dbc61d4553dcf0ef2aadbbb2f7",
    "biospecimen_gbm":   "9c9c5d5cff2b961c84dc55b08fb21c77486665ec6de5a953fcecf56a89125cf2",
    "biospecimen_pdac":  "87da63ba120934b94a46dc62b91e253c7238191147074719caee3a4e9dc530e1",
    "xwalk_ccrcc": "e823ce38eb24ffade8cb31098cf5378b84c59a897d5b5238e92288ad0f755830",
    "xwalk_luad":  "967cabf18cb6b6c3646aaf05005c46fcdccd7b55eb8dc9882348e84ca4baa2bf",
    "xwalk_ucec":  "8e1d6f412c60b2f720b9875f9e813136e51773457f21d651912311e9d2d59fc4",
    "xwalk_gbm":   "6a9458de508a60471877daf5bcc60f3bb3f595ec9d97b7fe660af254d76f1907",
    "xwalk_pdac":  "9c17fa2fc3d5279d1c4c30d38f5c58439c4c467af2d77bcaed1289d9d340bb56",
    "scores_ccrcc": "f76f7b82222ff55807add12ac5e837471f5c266c88216c515200ae457b196455",
    "scores_luad":  "f06ff9a3a1f354140c7c9e7269b81045cf16e578adc77804e4bbe8b09998682a",
    "scores_ucec":  "04035b4efd5886d65d2a701dfbdb7cf7ee6d91d5dd4c6f1ac5f6451af7d13545",
    "scores_gbm":   "73294e8571f983398ef800cb1f057b90dab9b7cf2529028a0db1b5ec1b584a83",
    "scores_pdac":  "93885f0f87308be6312511af4975066ee9fdb5ed89f1ae87e61c210b0fff6acf",
}


def sha256(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return "unreadable"


def assert_hash(path, key):
    """FATAL unless sha256(path) == EXPECTED_SHA256[key]. No skip env: every input here
    is a fixed historical artifact committed in this repo, not a path that legitimately
    differs between the HPC and a local mirror."""
    got = sha256(path)
    want = EXPECTED_SHA256[key]
    if got != want:
        sys.exit(f"FATAL: sha256 mismatch for {key}: got={got} want={want} ({path}). "
                  f"This is the wrong file vintage -- stopping rather than reporting "
                  f"against an unverified input.")
    return got


def load_biospecimen(coh, pid):
    """Return (status_by_aliquot, pt_rows): status_by_aliquot maps every
    aliquot_submitter_id in the raw snapshot to its aliquot_status; pt_rows is the list
    of raw rows with sample_type == 'Primary Tumor' (one row per distinct aliquot --
    verified no aliquot_submitter_id repeats within a study's Primary Tumor rows)."""
    path = os.path.join(RAW_DIR, f"r3_biospecimen_{pid}.json")
    assert_hash(path, f"biospecimen_{coh}")
    with open(path, encoding="utf-8") as f:
        payload = json.load(f)
    rows = payload["data"]["biospecimenPerStudy"]
    status_by_aliquot = {}
    pt_rows = []
    seen_pt_aliquots = set()
    for r in rows:
        status_by_aliquot[r["aliquot_submitter_id"]] = r.get("aliquot_status")
        if (r.get("sample_type") or "").strip() == "Primary Tumor":
            a = r["aliquot_submitter_id"]
            if a in seen_pt_aliquots:
                sys.exit(f"FATAL: {coh} ({pid}): aliquot {a!r} appears more than once "
                          f"among Primary Tumor rows -- the 'one row per aliquot' "
                          f"assumption this script relies on does not hold; needs "
                          f"inspection before this table can be trusted.")
            seen_pt_aliquots.add(a)
            pt_rows.append(dict(case_submitter_id=r["case_submitter_id"],
                                 aliquot_submitter_id=a,
                                 aliquot_status=r.get("aliquot_status"),
                                 sample_type="Primary Tumor"))
    return status_by_aliquot, pt_rows


def load_xwalk(coh, xw_name):
    """Return case_id -> sorted list of aliquot_id, from the local tumour crosswalk."""
    path = os.path.join(HERE, xw_name)
    assert_hash(path, f"xwalk_{coh}")
    case_to_aliquots = defaultdict(list)
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            case_to_aliquots[row["case_id"]].append(row["aliquot_id"])
    return {c: sorted(a) for c, a in case_to_aliquots.items()}


def load_analysed_cases(coh):
    """Return the sorted analysed-case list for a cohort, from the first column of
    figures/figdata/scores_<cohort>.csv (no pandas dependency needed for one column)."""
    path = os.path.join(FIGDATA_DIR, f"scores_{coh}.csv")
    assert_hash(path, f"scores_{coh}")
    cases = []
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)  # header
        for row in reader:
            if row and row[0]:
                cases.append(row[0])
    return sorted(set(cases))


def main():
    all_rows = []
    summary = {}

    for coh in COHORT_ORDER:
        pid, _json_name, xw_name = STUDIES[coh]
        status_by_aliquot, pt_rows = load_biospecimen(coh, pid)
        for r in pt_rows:
            all_rows.append(dict(study=coh, **r))

        case_to_aliquots = load_xwalk(coh, xw_name)
        cases = load_analysed_cases(coh)

        unmapped = []
        disqualified = []
        multi_aliquot = {}
        for c in cases:
            aliquots = case_to_aliquots.get(c, [])
            if not aliquots:
                unmapped.append(c)
                continue
            if len(aliquots) > 1:
                multi_aliquot[c] = aliquots
            statuses = {a: status_by_aliquot.get(a, "NOT_IN_SNAPSHOT") for a in aliquots}
            if any(v == "Disqualified" for v in statuses.values()):
                disqualified.append((c, statuses))

        summary[coh] = dict(
            pdc_study_id=pid, n_analysed_cases=len(cases),
            n_unmapped_analysed_cases=len(unmapped), unmapped_cases=unmapped,
            disqualified_analysed_cases=[c for c, _ in disqualified],
            disqualified_detail={c: st for c, st in disqualified},
            multi_aliquot_analysed_cases=multi_aliquot,
        )

    all_rows.sort(key=lambda r: (r["study"], r["case_submitter_id"], r["aliquot_submitter_id"]))

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    fieldnames = ["study", "case_submitter_id", "aliquot_submitter_id", "aliquot_status", "sample_type"]
    with open(OUT_PATH, "w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r[k] for k in fieldnames})

    out_hash = sha256(OUT_PATH)
    print(f"wrote {OUT_PATH}")
    print(f"  rows: {len(all_rows)}")
    print(f"  sha256: {out_hash}")
    print()
    print("Per-cohort cross-check (analysed cases vs. local tumour crosswalk vs. PDC status):")
    for coh in COHORT_ORDER:
        s = summary[coh]
        print(f"\n=== {coh} ({s['pdc_study_id']}) ===")
        print(f"  analysed cases: {s['n_analysed_cases']}")
        if s["unmapped_cases"]:
            print(f"  ANALYSED CASES WITH NO CROSSWALK ROW AT ALL: {s['unmapped_cases']}")
        else:
            print("  every analysed case has >=1 crosswalk row: yes")
        print(f"  disqualified analysed-aliquot cases ({len(s['disqualified_analysed_cases'])}): "
              f"{s['disqualified_analysed_cases']}")
        for c, st in s["disqualified_detail"].items():
            print(f"      {c}: {st}")
        print(f"  analysed cases with >1 crosswalk aliquot ({len(s['multi_aliquot_analysed_cases'])}): "
              f"{s['multi_aliquot_analysed_cases']}")

    return summary


if __name__ == "__main__":
    main()
