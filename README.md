# MorphoResidual

Analysis code for a study of whether tumour H&E morphology predicts protein
abundance beyond what is explained by a gene's own mRNA, evaluated across
five CPTAC cancer cohorts (CCRCC, LUAD, UCEC, GBM, PDAC) with replication
attempts in independent TCGA/CPTAC-3 cohorts.

Manuscript: submitted to *Nature Communications*; a citation and DOI will be
added here on publication.

## What is (and isn't) in this repository

This repository holds the residual-estimation pipeline, the robustness and
replication analyses, and the figure-generation scripts referenced in the
manuscript's Methods and Code availability sections.

It does **not** hold: the manuscript source; raw whole-slide images, protein
spectra or sequence data (available from PDC/GDC/TCIA under their own
terms, see the manuscript's Data availability statement); or the large
derived tables and slide/patient manifests referenced as Supplementary Data
(deposited separately, e.g. via Zenodo, alongside the accepted manuscript).

## Layout

- `figures/*.py` — one script per main/supplementary figure or table,
  reading precomputed result tables and writing the published figure PDFs.
  `figures/mrstyle.py` is the shared plotting style module.
- `server_export/scripts/` — the residual-estimation, acquisition-robustness
  and cross-cohort replication pipeline, run on the source institution's
  HPC cluster against CPTAC/TCGA/GDC data.
- Top-level `*.py` — cohort-level analysis scripts (clinical association,
  composition controls, subsampling/ablation diagnostics, encoder
  comparisons, external-cohort replication, etc.) referenced by file name
  throughout the Methods.

## Requirements and reproducibility

Software versions, random seeds and foundation-model encoder weights used
are listed in the manuscript's Methods ("Software, versions, seeds and
encoder weights"). Third-party model weights (e.g. UNI, Phikon) are not
redistributed here and remain subject to their own providers' licences.

## License

MIT — see [LICENSE](LICENSE). This licence covers the code in this
repository only; it does not extend to third-party model weights it loads.
