#!/usr/bin/env python3
"""
C1-LUAD confirmatory/discovery independence check (Step 0b, branch pdc,
2026-09-29). METADATA ONLY. Fail-closed: any missing file, any unexpected
discovery case count, or any overlap between the 113 confirmatory-independent
cases and the discovery cohort is FATAL (non-zero exit). This is the general
risk check Step 0 §0 flagged ("a general risk to check for every future C1
cohort"), applied here as its own script per the task, mirroring
`c1_run_test.py`'s `prepare()` FATAL guard for UCEC.

Compares the 113-case list (server_export/scripts/pinned/c1_case_ids_luad.txt,
written by c1_luad_derive_cases.py) against:
  (a) the discovery slide map, server_export/scripts/slide_type_map_luad.tsv
      (column case_submitter_id; must hold exactly 111 cases) -- case-level
      overlap check.
  (b) the discovery aliquot crosswalk,
      server_export/scripts/aliquot_to_case_tumor_luad.tsv (columns
      aliquot_id, case_id, sample_type) -- overlap check at BOTH the aliquot
      level (PDC000489 Primary Tumor aliquot_submitter_ids of the 113 cases,
      read from the same raw JSON c1_luad_derive_cases.py reads, vs. this
      file's aliquot_id column) and the case level (case_id column vs. the
      113 list).

Every invocation runs BOTH checks and reports both: (1) the normal check on the
pinned 113-case list, which must PASS, and (2) a self-test that deliberately
re-adds the 7 known bridge cases to the list (in memory only -- it never
touches the pinned file) and asserts the check then exits FATAL, as a check
that this script's own overlap logic actually fires rather than passing
vacuously. There is no separate opt-in mode: `--selftest` on the command line
is accepted but has no effect, since the self-test always runs alongside the
normal check.

Usage:
  python c1_luad_independence_check.py   # runs both checks; exits non-zero if
                                          # either does not behave as expected
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PINNED_DIR = os.path.join(HERE, "pinned")
RAW_DIR = os.path.join(
    HERE, "..", "..", "review", "C1_LUAD_STEP0B_2026-09-29", "raw")

CASE_IDS_PATH = os.path.join(PINNED_DIR, "c1_case_ids_luad.txt")
SLIDE_MAP_PATH = os.path.join(HERE, "slide_type_map_luad.tsv")
ALIQUOT_XWALK_PATH = os.path.join(HERE, "aliquot_to_case_tumor_luad.tsv")
CONF_BIOSPECIMEN_PATH = os.path.join(RAW_DIR, "biospecimenPerStudy_PDC000489.json")

EXPECTED_N_FINAL = 113
EXPECTED_DISCOVERY_CASES = 111
BRIDGE_CASES_FOR_SELFTEST = {
    "C3L-02348", "C3L-02350", "C3N-00545", "C3N-00738",
    "C3N-01023", "C3N-01024", "C3N-02003",
}


class CheckFatal(SystemExit):
    pass


def fatal(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    raise CheckFatal(1)


def require_file(path, what):
    if not os.path.exists(path):
        fatal(f"missing {what}: {path}")


def load_case_ids(path):
    with open(path, encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def load_slide_map_cases(path):
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if "case_submitter_id" not in rows[0]:
        fatal(f"{path} has no case_submitter_id column (got {list(rows[0].keys())})")
    return {r["case_submitter_id"] for r in rows}


def load_aliquot_xwalk(path):
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for col in ("aliquot_id", "case_id"):
        if col not in rows[0]:
            fatal(f"{path} has no {col} column (got {list(rows[0].keys())})")
    aliquots = {r["aliquot_id"] for r in rows}
    cases = {r["case_id"] for r in rows}
    return aliquots, cases


def load_confirmatory_tumor_aliquots(case_ids):
    """PDC000489 Primary Tumor aliquot_submitter_id for exactly the given cases,
    read from the same raw JSON c1_luad_derive_cases.py reads. No network call."""
    require_file(CONF_BIOSPECIMEN_PATH, "confirmatory raw biospecimen JSON")
    with open(CONF_BIOSPECIMEN_PATH, encoding="utf-8") as f:
        d = json.load(f)
    rows = d.get("data", {}).get("biospecimenPerStudy")
    if rows is None:
        fatal(f"{CONF_BIOSPECIMEN_PATH} has no data.biospecimenPerStudy")
    case_set = set(case_ids)
    aliquots = {r["aliquot_submitter_id"] for r in rows
                if r.get("sample_type") == "Primary Tumor"
                and r.get("case_submitter_id") in case_set}
    return aliquots


def run_check(case_ids, label):
    print(f"\n=== independence check [{label}]: {len(case_ids)} confirmatory case(s) ===")
    case_set = set(case_ids)

    # --- (a) discovery slide map ---
    require_file(SLIDE_MAP_PATH, "discovery slide map")
    disc_slide_cases = load_slide_map_cases(SLIDE_MAP_PATH)
    print(f"  discovery slide map: {len(disc_slide_cases)} distinct case_submitter_id")
    if len(disc_slide_cases) != EXPECTED_DISCOVERY_CASES:
        fatal(f"discovery slide map case count = {len(disc_slide_cases)}, expected "
              f"{EXPECTED_DISCOVERY_CASES}")
    overlap_a = case_set & disc_slide_cases
    print(f"  case-level overlap vs. slide map: {len(overlap_a)} {sorted(overlap_a) if overlap_a else ''}")
    if overlap_a:
        fatal(f"[{label}] {len(overlap_a)} confirmatory case(s) appear in the "
              f"discovery slide map (case-level independence violated): {sorted(overlap_a)}")

    # --- (b) discovery aliquot crosswalk: case level and aliquot level ---
    require_file(ALIQUOT_XWALK_PATH, "discovery aliquot crosswalk")
    disc_xwalk_aliquots, disc_xwalk_cases = load_aliquot_xwalk(ALIQUOT_XWALK_PATH)
    print(f"  discovery aliquot crosswalk: {len(disc_xwalk_cases)} distinct case_id, "
          f"{len(disc_xwalk_aliquots)} distinct aliquot_id")
    if len(disc_xwalk_cases) != EXPECTED_DISCOVERY_CASES:
        fatal(f"discovery aliquot crosswalk case count = {len(disc_xwalk_cases)}, "
              f"expected {EXPECTED_DISCOVERY_CASES}")

    overlap_b_case = case_set & disc_xwalk_cases
    print(f"  case-level overlap vs. aliquot crosswalk: {len(overlap_b_case)} "
          f"{sorted(overlap_b_case) if overlap_b_case else ''}")
    if overlap_b_case:
        fatal(f"[{label}] {len(overlap_b_case)} confirmatory case(s) appear in the "
              f"discovery aliquot crosswalk (case-level independence violated): "
              f"{sorted(overlap_b_case)}")

    conf_tumor_aliquots = load_confirmatory_tumor_aliquots(case_ids)
    print(f"  confirmatory PDC000489 Primary Tumor aliquots for these "
          f"{len(case_ids)} case(s): {len(conf_tumor_aliquots)}")
    overlap_b_aliquot = conf_tumor_aliquots & disc_xwalk_aliquots
    print(f"  aliquot-level overlap vs. aliquot crosswalk: {len(overlap_b_aliquot)} "
          f"{sorted(overlap_b_aliquot) if overlap_b_aliquot else ''}")
    if overlap_b_aliquot:
        fatal(f"[{label}] {len(overlap_b_aliquot)} confirmatory Primary Tumor "
              f"aliquot(s) are identical to a discovery Primary Tumor aliquot "
              f"(aliquot-level independence violated): {sorted(overlap_b_aliquot)}")

    print(f"  [{label}] PASS: zero overlap at case level (slide map + aliquot "
          f"crosswalk) and at aliquot level.")
    return True


def main():
    require_file(CASE_IDS_PATH, "113-case pinned list")
    case_ids = load_case_ids(CASE_IDS_PATH)
    print(f"loaded {len(case_ids)} case ids from {CASE_IDS_PATH}")
    if len(case_ids) != EXPECTED_N_FINAL:
        fatal(f"pinned case list has {len(case_ids)} ids, expected {EXPECTED_N_FINAL}")

    normal_ok = None
    try:
        normal_ok = run_check(case_ids, "normal")
    except CheckFatal:
        normal_ok = False
        print("\n[normal] result: FATAL (see above) -- this should NOT happen; "
              "the pinned 113-case list is supposed to already exclude all "
              "discovery/bridge overlap.")

    selftest_ok = None
    if True:
        # Always run the self-test in the same invocation so both results are
        # reported in one place, per the task's "run both ... and report both".
        contaminated = sorted(set(case_ids) | BRIDGE_CASES_FOR_SELFTEST)
        print(f"\n--- --selftest: re-adding the 7 known bridge cases in memory "
              f"only ({len(contaminated)} cases total) ---")
        try:
            run_check(contaminated, "selftest (deliberately contaminated)")
            selftest_ok = False
            print("\n[selftest] result: DID NOT RAISE FATAL -- this is a bug in "
                  "the overlap logic; the selftest is supposed to fail closed.")
        except CheckFatal:
            selftest_ok = True
            print("\n[selftest] result: FATAL as expected (overlap logic correctly "
                  "fires on the deliberately contaminated list).")

    print("\n=== SUMMARY ===")
    print(f"  normal mode (113 pinned cases):        {'PASS' if normal_ok else 'FATAL'}")
    print(f"  --selftest mode (+7 bridge cases):     "
          f"{'FATAL as expected (correct)' if selftest_ok else 'DID NOT FATAL (BUG)'}")

    if not normal_ok or not selftest_ok:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
