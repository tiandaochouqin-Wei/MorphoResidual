#!/usr/bin/env python3
"""
Go/no-go check BEFORE spending GPU time: do RNA-seq, WES, protein, and WSI
actually share overlapping CPTAC case IDs? Run this after download_gdc.py,
download_wsi.sh, and the manual PDC step -- before submit_extract.lsf.

CPTAC case IDs look like C3L-01978 / C3N-00495 (project-prefix + 5 digits).
This scans filenames/manifests across all four sources with that pattern
rather than assuming a fixed column layout, since we haven't verified PDC's
sample.txt schema live (it requires the manual portal export first).
"""
import os
import re
import sys
from pathlib import Path

CASE_RE = re.compile(r"C3[A-Z]-\d{5}")
OMICS_DIR = Path(os.environ.get("MORPHO_OMICS_DIR", "./omics"))
WSI_DIR = Path(os.environ.get("MORPHO_WSI_DIR", "./WSI"))


def cases_from_manifest(path):
    if not path.exists():
        return set()
    return set(CASE_RE.findall(path.read_text(errors="replace")))


def cases_from_dir_names(d):
    if not d.exists():
        return set()
    out = set()
    for p in d.rglob("*"):
        out |= set(CASE_RE.findall(p.name))
    return out


ALIQUOT_RE = re.compile(r"CPT\d{7,}")


def load_aliquot_crosswalk():
    # PDC's tmt10.tsv header uses TMT-plex aliquot IDs (e.g. CPT0079430001),
    # not case barcodes -- confirmed live 2026-07-15 via head -1 on the real
    # file. The crosswalk (aliquot_submitter_id -> case_submitter_id) comes
    # from PDC's biospecimenPerStudy GraphQL query for study PDC000127
    # (dbe94609-1fb3-11e9-b7f8-0a80fada099c); 194 patient rows, 14 QC/NCI7
    # reference channels excluded.
    path = Path(__file__).parent / "aliquot_to_case.tsv"
    mapping = {}
    if not path.exists():
        return mapping
    with path.open(encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            aliquot, case = line.rstrip("\n").split("\t")
            mapping[aliquot] = case
    return mapping


def cases_from_protein_file(d):
    crosswalk = load_aliquot_crosswalk()
    for p in d.glob("*.tsv"):
        if "tmt10" not in p.name:
            continue
        header = p.open(errors="replace").readline()
        aliquots = set(ALIQUOT_RE.findall(header))
        if not aliquots:
            continue
        return {crosswalk[a] for a in aliquots if a in crosswalk}
    return set()


def main():
    rna = cases_from_manifest(OMICS_DIR / "manifest_rna.tsv")
    wes = cases_from_manifest(OMICS_DIR / "manifest_wes.tsv")
    protein = cases_from_protein_file(OMICS_DIR / "protein")
    # ccrcc_real = the actual CPTAC-CCRCC Discovery pathology cohort (390 SVS,
    # verified live 2026-07-15 via Faspex). The old "ccrcc" dir was a wrong
    # download (radiology CT/MR pulled from the wrong .tcia manifest) -- do
    # not point back at it.
    wsi = cases_from_dir_names(WSI_DIR / "raw" / "ccrcc_real")
    if not wsi:
        wsi = cases_from_dir_names(WSI_DIR / "emb" / "ccrcc_real")

    print(f"RNA-seq cases:  {len(rna)}")
    print(f"WES cases:      {len(wes)}")
    print(f"Protein cases:  {len(protein)}  "
          f"{'(WARNING: 0 found -- check tmt10.tsv header format by hand)' if not protein else ''}")
    print(f"WSI cases:      {len(wsi)}  "
          f"{'(WARNING: 0 found -- slide filenames may not embed case ID; check manually)' if not wsi else ''}")

    all_four = rna & wes & wsi & (protein if protein else rna & wes & wsi)
    three_no_protein = rna & wes & wsi
    print(f"\nOverlap RNA∩WES∩WSI (excl. protein, which needs manual-step data): {len(three_no_protein)}")
    if protein:
        print(f"Overlap across all 4 modalities: {len(all_four)}")

    n = len(all_four) if protein else len(three_no_protein)
    print(f"\n{'PASS' if n >= 50 else 'FAIL'}: {n} patients with usable multi-modal overlap "
          f"({'>=50, proceed to feature extraction' if n >= 50 else '<50, do not spend GPU time yet -- diagnose the ID mismatch first'})")
    sys.exit(0 if n >= 50 else 1)


if __name__ == "__main__":
    main()
