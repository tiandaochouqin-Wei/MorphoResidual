#!/usr/bin/env python3
"""
Build the pinned per-slide embedding-medium (FFPE/OCT) table for the
batch_leak_check 'medium' axis.

Purpose
  Extracts slide_submitter_id, case_submitter_id, disease_type and
  embedding_medium from the public PathDB per-slide metadata file for every row
  whose Slide_ID matches a GDC slide_submitter_id and whose Disease_Type is one
  of the five analysed cancers, and writes it as a sorted, LF, UTF-8 TSV under
  server_export/scripts/pinned/, for upload to the HPC pinned/ directory that
  the patched batch_leak_check.case_batch_labels reads. LOCAL-ONLY tool: its
  input is a public PathDB download cached in this checkout; nothing here reads
  or writes any HPC-only RNA/protein/WSI-embedding data, and this script is not
  meant to run on the HPC.

  Post hoc, and governed by review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md
  (section M, "M0. Labels"). Does not compute, read or alter any discovery or
  confirmatory result. PathDB Specimen_Type is never used (per M0).

Input (sha256-asserted against EXPECTED_SHA256 below; FATAL on mismatch)
  review/BATCH_AXES_SCOPING_2026-10-02/pathdb/cptac_metadata_07-09-2024.csv
  (downloaded 2026-10-02T12:06:15Z; this is the same file and hash cited in
  review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md section M0 and in this task's
  computed instructions)

Output (4 columns, header
  "slide_submitter_id<TAB>case_submitter_id<TAB>disease_type<TAB>embedding_medium",
  sorted by slide_submitter_id, LF line endings)
  server_export/scripts/pinned/slide_embedding_medium.tsv

Row selection
  Slide_ID matches ^C3[LN]-\\d{5}-\\d+$ (a GDC slide_submitter_id) AND
  Disease_Type in {Clear Cell Renal Cell Carcinoma, Lung Adenocarcinoma,
  Uterine Corpus Endometrial Carcinoma, Glioblastoma,
  Pancreatic Ductal Adenocarcinoma}. embedding_medium is PathDB's
  Embedding_Medium value verbatim (including the empty string, where PathDB
  itself leaves it blank -- the downstream join in batch_leak_check.py treats
  blank the same as missing, it is not fabricated here).

FATAL conditions
  - a selected Slide_ID appears twice with two different non-empty
    Embedding_Medium values (a genuine PathDB conflict -- not observed when
    this script was written: see the report), or
  - EXPECTED_SHA256 mismatch on the input file.

Usage
  python build_pinned_slide_medium.py
  No arguments: both the source and destination paths are fixed relative to
  this file's own location.
"""
import csv
import hashlib
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(PAPER, "review", "BATCH_AXES_SCOPING_2026-10-02", "pathdb",
                    "cptac_metadata_07-09-2024.csv")
PINNED = os.path.join(HERE, "pinned")
OUT = os.path.join(PINNED, "slide_embedding_medium.tsv")

# sha256 of review/BATCH_AXES_SCOPING_2026-10-02/pathdb/cptac_metadata_07-09-2024.csv,
# the public PathDB download named in review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md
# section M0 (downloaded 2026-10-02T12:06:15Z). FATAL on mismatch rather than
# silently joining against a different vintage of this file.
EXPECTED_SHA256_PATHDB = (
    "6f999cec4381ea34f245c2faabdbc828837b0125f6a9950bb5d42ec4febe6da0"
)

SLIDE_ID_RE = re.compile(r"^C3[LN]-\d{5}-\d+$")
DISEASES = {
    "Clear Cell Renal Cell Carcinoma",
    "Lung Adenocarcinoma",
    "Uterine Corpus Endometrial Carcinoma",
    "Glioblastoma",
    "Pancreatic Ductal Adenocarcinoma",
}


def sha256_of(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def assert_hash(path, expected, label):
    got = sha256_of(path)
    if got != expected:
        sys.exit(f"FATAL: sha256 mismatch for {label}: got={got} want={expected} "
                  f"({path}). Refusing to join against an unverified PathDB "
                  f"vintage.")
    return got


def main():
    assert_hash(SRC, EXPECTED_SHA256_PATHDB, "cptac_metadata_07-09-2024.csv")

    with open(SRC, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for col in ("Slide_ID", "Case_ID", "Disease_Type", "Embedding_Medium"):
            if col not in reader.fieldnames:
                sys.exit(f"FATAL: {SRC} is missing expected column {col!r}.")
        rows = list(reader)

    seen = {}  # slide_submitter_id -> (case_submitter_id, disease_type, medium)
    for r in rows:
        slide = (r["Slide_ID"] or "").strip()
        if not SLIDE_ID_RE.match(slide):
            continue
        disease = (r["Disease_Type"] or "").strip()
        if disease not in DISEASES:
            continue
        case = (r["Case_ID"] or "").strip()
        medium = (r["Embedding_Medium"] or "").strip()
        if slide in seen:
            prev_case, prev_disease, prev_medium = seen[slide]
            if prev_medium and medium and prev_medium != medium:
                sys.exit(f"FATAL: {SRC} slide {slide} has conflicting "
                          f"Embedding_Medium values {prev_medium!r} vs "
                          f"{medium!r}.")
            if prev_case != case or prev_disease != disease:
                sys.exit(f"FATAL: {SRC} slide {slide} has conflicting "
                          f"Case_ID/Disease_Type across duplicate rows: "
                          f"({prev_case!r}, {prev_disease!r}) vs "
                          f"({case!r}, {disease!r}).")
            # keep whichever row carries a non-empty medium, if either does
            seen[slide] = (case, disease, prev_medium or medium)
        else:
            seen[slide] = (case, disease, medium)

    os.makedirs(PINNED, exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("slide_submitter_id\tcase_submitter_id\tdisease_type\t"
                "embedding_medium\n")
        for slide in sorted(seen):
            case, disease, medium = seen[slide]
            f.write(f"{slide}\t{case}\t{disease}\t{medium}\n")

    counts = Counter((d, m) for _, d, m in seen.values())
    print(f"wrote {OUT}  ({len(seen)} slides)")
    print("counts per disease_type x embedding_medium:")
    for (disease, medium), n in sorted(counts.items()):
        print(f"  {disease:<40} {medium or '(blank)':<8} {n}")
    print(f"sha256={sha256_of(OUT)}")


if __name__ == "__main__":
    main()
