# Fetching CCRCC protein abundance from PDC (manual, one-time, ~85 MB)

**Why this step can't be scripted end-to-end:** PDC's GraphQL API (verified
live) gives full file *metadata* for a study, but not a stable download URL —
byte-level access only exists behind a **time-limited signed URL** minted when
you export a manifest from the web portal. Per PDC's own docs, that exported
manifest **expires after 7 days**. This is a portal-side design choice, not a
gap in our scripting — there's no API-only path around it.

The good news: we only need 3 small files (~85 MB total), not the raw mass
spec data, so this is a five-minute manual step.

## What we verified we need (via live GraphQL query against PDC000127)

Study **PDC000127** = "CPTAC CCRCC Discovery Study - Proteome" (110-tumor TMT10
discovery cohort). Under `data_category = "Protein Assembly"`:

| file | size | what it is |
|---|---|---|
| `CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv` | 79.2 MB | **gene-level TMT10 protein abundance matrix — this is the file the residual analysis needs** |
| `CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.summary.tsv` | 4.0 MB | per-protein/sample summary stats |
| `CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.sample.txt` | 6 KB | sample/aliquot ID mapping (needed to join to GDC case IDs) |

**Do NOT select "Raw Mass Spectra" or "Peptide Spectral Matches"** — those are
575 + 1150 files, ~400 GB, and are not what we need (confirmed live; the study
has them but they're instrument-level intermediate data, not the harmonized
abundance table).

## Steps

1. Go to https://proteomic.datacommons.cancer.gov/pdc/browse (or the current
   PDC portal) and search/filter for study **PDC000127**.
2. On the study's file list, filter `Data Category = Protein Assembly`.
3. Select the 3 files above (or just select all 4 in that category if the
   `.peptides.tsv` filter is inconvenient — it's only 174 MB, still cheap).
4. Use **"Export File Manifest"** (CSV or TSV) — download it locally.
5. Download the **PDC Data Download Client** from
   https://proteomic.datacommons.cancer.gov/pdc/data-download-documentation
   (same family of tool as `gdc-client`, standalone, no install needed).
6. Within 7 days of exporting the manifest, run it (on the HPC login node, or
   locally then `scp`/WinSCP the 3 small files over):
   ```bash
   ./pdc-client download -m /path/to/exported_manifest.tsv -d ./omics/protein
   ```
7. Confirm you got `CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv`
   (~79 MB) in `omics/protein/`.

If you'd rather not install the PDC client at all: since it's only ~85 MB,
downloading the 3 files straight from the portal's browser "Download" button
on each file page and transferring via WinSCP is equally fine for a pilot.
