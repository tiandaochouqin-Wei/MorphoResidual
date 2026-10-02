#!/usr/bin/env python3
"""
Build pinned per-cohort case -> TMT-plex maps for the batch_leak_check 'plex' axis.

Purpose
  Writes columns 1-2 (case_id, plex) of each
  review/BATCH_AXES_SCOPING_2026-10-02/plex_map_<cohort>.tsv as a minimal
  2-column, LF, sorted-by-case_id TSV under server_export/scripts/pinned/, for
  upload to the HPC pinned/ directory that the patched
  batch_leak_check.case_batch_labels reads. LOCAL-ONLY tool: its inputs
  (review/BATCH_AXES_SCOPING_2026-10-02/plex_map_*.tsv) exist only in this local
  checkout, built by that folder's own build_plex_maps.py from cached PDC JSON;
  nothing here reads or writes any HPC-only RNA/protein/WSI-embedding data, and
  this script is not meant to run on the HPC.

  Post hoc, and governed by review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md
  (section P, "Inputs"). Does not compute, read or alter any discovery or
  confirmatory result.

Inputs (sha256-asserted against EXPECTED_SHA256 below; FATAL on mismatch)
  review/BATCH_AXES_SCOPING_2026-10-02/plex_map_ccrcc.tsv
  review/BATCH_AXES_SCOPING_2026-10-02/plex_map_luad.tsv
  review/BATCH_AXES_SCOPING_2026-10-02/plex_map_ucec.tsv
  review/BATCH_AXES_SCOPING_2026-10-02/plex_map_gbm.tsv
  review/BATCH_AXES_SCOPING_2026-10-02/plex_map_pdac.tsv
  review/BATCH_AXES_SCOPING_2026-10-02/plex_map_luad_conf.tsv

Outputs (2 columns, header "case_id<TAB>plex", case_id sorted, LF line endings)
  server_export/scripts/pinned/case_to_plex_ccrcc.tsv
  server_export/scripts/pinned/case_to_plex_luad.tsv
  server_export/scripts/pinned/case_to_plex_ucec.tsv
  server_export/scripts/pinned/case_to_plex_gbm.tsv
  server_export/scripts/pinned/case_to_plex_pdac.tsv
  server_export/scripts/pinned/case_to_plex_luad_c1.tsv   (from plex_map_luad_conf.tsv)

Usage
  python build_pinned_plex_maps.py
  No arguments: both the source and destination directories are fixed relative
  to this file's own location, since this checkout is the only place the source
  scoping maps exist.

Each output row's plex value is exactly the single resolved plex the source
plex_map_<cohort>.tsv already picked for a multi-plex case (that folder's
build_plex_maps.py: "majority by number of measured tumor channels, ties ->
lowest plex"), e.g. UCEC C3N-01825 -> plex 10 (its non-Disqualified aliquot),
per review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md section P.
"""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
SCOPING = os.path.join(PAPER, "review", "BATCH_AXES_SCOPING_2026-10-02")
PINNED = os.path.join(HERE, "pinned")

# cohort -> source scoping file (review/BATCH_AXES_SCOPING_2026-10-02/)
COHORTS = {
    "ccrcc": "plex_map_ccrcc.tsv",
    "luad": "plex_map_luad.tsv",
    "ucec": "plex_map_ucec.tsv",
    "gbm": "plex_map_gbm.tsv",
    "pdac": "plex_map_pdac.tsv",
    "luad_c1": "plex_map_luad_conf.tsv",
}

# Frozen sha256 of every scoping input this script reads (computed 2026-10-02,
# the day build_plex_maps.py wrote them). A mismatch means the scoping maps
# changed since this pinning was written -- FATAL rather than silently pinning a
# different plex assignment than the one that was reviewed.
EXPECTED_SHA256 = {
    "plex_map_ccrcc.tsv": "f37e5e6e4030b7ede5c6c1b10fa7420720f3964dd40e448d66d2fbe5a5aec2b2",
    "plex_map_luad.tsv": "c3037a8a48d36620b8b0b6386d023999c5abb7bf6e5c47a500790e8e7b7e3087",
    "plex_map_ucec.tsv": "e17cbafb5bf231ea3e2ea219fa977cef46390f1a5cdfb208c60911a4aad01a03",
    "plex_map_gbm.tsv": "351dbae3a106fa4f37b667ca686265b9fb9860f5c4f4d89aec4a2558c8ec8d07",
    "plex_map_pdac.tsv": "b4fa6d75043b5da50dd7912000822b857d16fb01f6db600737801794409e3a75",
    "plex_map_luad_conf.tsv": "13397cef6054e284df5e9ff1f4e56b78ca23605d263671373be1c9ee7e0fab3c",
}


def sha256_of(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def assert_hash(path, name):
    got = sha256_of(path)
    want = EXPECTED_SHA256[name]
    if got != want:
        sys.exit(f"FATAL: sha256 mismatch for {name}: got={got} want={want} "
                  f"({path}). The scoping plex map changed since this pinning was "
                  f"written -- re-derive EXPECTED_SHA256 and re-check the plex "
                  f"assignment before trusting this output.")
    return got


def build_one(src_name):
    """Read columns 1-2 of a scoping plex_map_<cohort>.tsv; return sorted
    [(case_id, plex), ...], FATAL on a missing plex or a case with conflicting
    plex values across duplicate rows."""
    src_path = os.path.join(SCOPING, src_name)
    assert_hash(src_path, src_name)
    with open(src_path, encoding="utf-8-sig", newline="") as f:
        lines = f.read().splitlines()
    if not lines:
        sys.exit(f"FATAL: {src_path} is empty.")
    header = lines[0].split("\t")
    if header[:2] != ["case_id", "plex"]:
        sys.exit(f"FATAL: unexpected header in {src_path}: {header[:2]!r} "
                  f"(expected ['case_id', 'plex']).")
    seen = {}
    for ln in lines[1:]:
        if not ln.strip():
            continue
        cols = ln.split("\t")
        case_id, plex = cols[0].strip(), cols[1].strip()
        if not case_id or not plex:
            sys.exit(f"FATAL: {src_path} has a case ({case_id!r}) with no "
                      f"resolved plex ({plex!r}) -- investigate before pinning.")
        if case_id in seen and seen[case_id] != plex:
            sys.exit(f"FATAL: {src_path} case {case_id} has conflicting plex "
                      f"values {seen[case_id]!r} vs {plex!r}.")
        seen[case_id] = plex
    return sorted(seen.items())


def write_pinned(rows, out_path):
    os.makedirs(PINNED, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("case_id\tplex\n")
        for case_id, plex in rows:
            f.write(f"{case_id}\t{plex}\n")


def main():
    for cohort, src_name in COHORTS.items():
        rows = build_one(src_name)
        out_path = os.path.join(PINNED, f"case_to_plex_{cohort}.tsv")
        write_pinned(rows, out_path)
        print(f"{cohort}: {len(rows)} cases -> {out_path}  "
              f"sha256={sha256_of(out_path)}")


if __name__ == "__main__":
    main()
