# MorphoResidual — CCRCC pilot data pipeline

All endpoints below were tested live against the real APIs on 2026-07-13, not
guessed from docs. Two things turned out differently than the obvious
approach and materially shape these scripts:

1. **CPTAC-CCRCC pathology slides are not in the standard NBIA REST API.**
   They ship as a static `.tcia` manifest + require the official NBIA Data
   Retriever CLI (its `DownloadServlet` 500s on a hand-built curl request —
   it needs the client's internal auth handshake). `download_wsi.sh` extracts
   the retriever from its RPM without root and drives it with the manifest.
2. **PDC (protein) has no stable download URL via the API.** GraphQL gives
   full metadata (confirmed: study PDC000127 → `Protein Assembly` category →
   the 79 MB `tmt10.tsv` gene-level abundance matrix we actually need, not
   the 400 GB of raw spectra also in that study). But byte access requires a
   7-day signed URL only obtainable by exporting a manifest from the PDC web
   portal. `fetch_pdc_protein.md` is a 5-minute manual step for this reason —
   not a shortcut I skipped, a real constraint on PDC's side.

RNA-seq and WES download cleanly end-to-end via `download_gdc.py` (pure
Python, no `gdc-client` binary needed — GDC's `/data/<uuid>` direct endpoint
works for open-access files). **Filter precisely as written**: the sibling
"Splice Junction Quantification" files under the same RNA-seq workflow are
`access=controlled` and will 403 if the filter is loosened.

## Order of operations

Run all of this from `mn02` (the login node — confirmed to have outbound
internet; GPU nodes may not).

```bash
export MORPHO_OMICS_DIR=/public/home/fjhui/ZW/omics
export MORPHO_WSI_DIR=/public/home/fjhui/ZW/WSI

# 1. RNA-seq + WES (fully automated, ~1 GB, minutes)
python download_gdc.py

# 2. WSI manifest + slides (semi-automated, ~150-400 GB, hours)
#    Read the script's step-1 note first: check for a newer .tcia manifest
#    version than the one hardcoded (v11_20230818) before a big run.
bash download_wsi.sh

# 3. Protein abundance (manual, ~85 MB, 5 minutes) -- see fetch_pdc_protein.md

# 4. Go/no-go: do the four sources actually share patients?
python verify_pilot_data.py
# If this reports FAIL, stop and diagnose ID mismatches before spending GPU time.

# 5. One-time env setup, then submit feature extraction
bash setup_env.sh
# requires HF_TOKEN with accepted MahmoodLab/UNI license -- see extract_features.py docstring
bsub < submit_extract.lsf
# if the interactive queue rejects non -Is submission, fall back to running
# extract_features.py by hand inside `bsub -q interactive -gpu "num=1" -Is bash`
```

## What you get at the end

- `omics/rna/*.tsv` — per-sample GENCODE-v36 STAR gene counts (open access)
- `omics/mutation/*.maf` — per-case masked somatic mutations
- `omics/protein/*.tsv` — gene-level TMT10 protein abundance matrix
- `WSI/emb/ccrcc/*.pt` — per-slide `{embeddings: [N_tiles, D], coords}` bags
  from the frozen UNI encoder, ready for the ABMIL aggregation + residual
  regression in `plan.md` §2.

`WSI/raw/ccrcc` holds the source slides until you pass `--delete-raw` to
`extract_features.py` (off by default — kept until you've spot-checked a few
embeddings look sane).

## Known gaps to watch

- WSI file format after NBIA download is unconfirmed (DICOM-wrapped vs native
  SVS) — `extract_features.py` auto-detects via `openslide` or `wsidicom`,
  but this hasn't been exercised against a real downloaded file yet. If
  `load_slide()` throws on the first real file, that's the thing to debug.
- `verify_pilot_data.py`'s protein/WSI case-ID extraction is regex-based
  (`C3[A-Z]-\d{5}`) because PDC's `sample.txt` schema and the WSI filename
  convention weren't directly inspectable without the manual PDC step. If it
  reports 0 cases found for either, open the file and check the ID format by
  eye rather than trusting the regex blindly.
