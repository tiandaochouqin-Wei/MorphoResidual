#!/usr/bin/env python3
"""
C1-LUAD TMT design metadata, run/plex level (Step 0b, branch pdc, task item 4).
METADATA ONLY. Reads ONLY the saved PDC `studyExperimentalDesign` raw JSON
under review/C1_LUAD_STEP0B_2026-09-29/raw/ -- no network call here.

IMPORTANT SCOPE NOTE (read before using this file): this produces a RUN-level
summary (one row per TMT plex/run), not a per-case or per-aliquot plex
assignment. After exhaustive query/field discovery (introspection is
disabled on PDC's endpoint; five rounds of deliberately-wrong-name "did you
mean" probing logged in query_log_pdc.tsv, covering studyExperimentalDesign,
uiStudyExperimentalDesign -- same underlying type, no additional fields --,
a nonexistent aliquotRunMetadata query, and filesPerStudy), no queryable
field links an aliquot_submitter_id/case_submitter_id to a
study_run_metadata_submitter_id (plex/run), and `aliquot_is_ref` (the
candidate reference-channel flag) is null on every one of the 27 confirmatory
and 25 discovery rows. So:
  - aliquot -> plex/run assignment for the 113 confirmatory cases: NOT
    derivable from PDC's public GraphQL metadata without a quantitative
    query (quantDataMatrix's column headers), which this branch's task does
    not authorise for PDC000489 (the confirmatory study) and which was
    therefore not attempted. `c1_luad_plex_map.tsv` (case -> plex) is NOT
    written by this script for that reason -- writing it would require
    fabricating the linkage.
  - reference-channel identity per plex: NOT derivable either (aliquot_is_ref
    is uniformly null in the metadata PDC returns for both studies).
  - What IS derivable and is written here: the run/plex-level design table
    (one row per TMT run: run name/date, TMT chemistry, analyte,
    experiment_number, number_of_fractions) for both studies, and the
    plex/run COUNT per study.

Usage:
  python c1_luad_plex_runs_summary.py
Outputs:
  server_export/scripts/pinned/c1_luad_plex_runs_summary.tsv
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(
    HERE, "..", "..", "review", "C1_LUAD_STEP0B_2026-09-29", "raw")
PINNED_DIR = os.path.join(HERE, "pinned")

EXPECTED_N_RUNS = {"PDC000489": 27, "PDC000153": 25}
EXPECTED_EXPERIMENT_TYPE = {"PDC000489": "TMT11", "PDC000153": "TMT10"}


def fatal(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def load_design(pdc_study_id):
    path = os.path.join(RAW_DIR, f"studyExperimentalDesign_{pdc_study_id}.json")
    if not os.path.exists(path):
        fatal(f"missing raw file {path} -- run fetch_run_level_design.py first")
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    if d.get("errors"):
        fatal(f"{path} is a saved GraphQL ERROR response: {d['errors']}")
    rows = d.get("data", {}).get("studyExperimentalDesign")
    if rows is None:
        fatal(f"{path} has no data.studyExperimentalDesign")
    return rows


def main():
    os.makedirs(PINNED_DIR, exist_ok=True)
    out_path = os.path.join(PINNED_DIR, "c1_luad_plex_runs_summary.tsv")

    all_rows = []
    for pid, tag in [("PDC000489", "confirmatory"), ("PDC000153", "discovery")]:
        rows = load_design(pid)
        print(f"{pid} [{tag}]: {len(rows)} run rows")
        if len(rows) != EXPECTED_N_RUNS[pid]:
            fatal(f"{pid} run count = {len(rows)}, expected {EXPECTED_N_RUNS[pid]}")
        exp_types = {r.get("experiment_type") for r in rows}
        if exp_types != {EXPECTED_EXPERIMENT_TYPE[pid]}:
            fatal(f"{pid} experiment_type values = {exp_types}, expected "
                  f"{{'{EXPECTED_EXPERIMENT_TYPE[pid]}'}}")
        ref_vals = {r.get("aliquot_is_ref") for r in rows}
        if ref_vals != {None}:
            print(f"  NOTE: aliquot_is_ref is NOT uniformly null for {pid} "
                  f"({ref_vals}) -- re-examine; earlier probing found it always "
                  f"null and this summary's scope note assumed that.")
        for r in sorted(rows, key=lambda r: r["study_run_metadata_submitter_id"]):
            all_rows.append({
                "pdc_study_id": pid,
                "cohort": tag,
                "study_run_metadata_submitter_id": r["study_run_metadata_submitter_id"],
                "experiment_number": r.get("experiment_number"),
                "experiment_type": r.get("experiment_type"),
                "analyte": r.get("analyte"),
                "number_of_fractions": r.get("number_of_fractions"),
                "aliquot_is_ref": r.get("aliquot_is_ref"),
            })

    cols = ["pdc_study_id", "cohort", "study_run_metadata_submitter_id",
            "experiment_number", "experiment_type", "analyte",
            "number_of_fractions", "aliquot_is_ref"]
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(cols) + "\n")
        for r in all_rows:
            f.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
    print(f"\nwrote {len(all_rows)} run rows -> {out_path}")

    print("\n=== SUMMARY ===")
    print(f"  confirmatory (PDC000489): {EXPECTED_N_RUNS['PDC000489']} TMT11 plexes/runs")
    print(f"  discovery   (PDC000153): {EXPECTED_N_RUNS['PDC000153']} TMT10 plexes/runs")
    print("  reference-channel identity: NOT populated in metadata (aliquot_is_ref "
          "null on every row, both studies)")
    print("  case/aliquot -> plex assignment: NOT derivable from any metadata-only "
          "PDC GraphQL field found (see module docstring); c1_luad_plex_map.tsv is "
          "deliberately NOT produced.")


if __name__ == "__main__":
    main()
