# MorphoResidual

Analysis code for a study of whether tumour H&E morphology predicts protein
abundance beyond what is explained by a gene's own mRNA, evaluated across
five CPTAC cancer cohorts (CCRCC, LUAD, UCEC, GBM, PDAC) with replication
attempts in independent TCGA/CPTAC-3 cohorts.

Manuscript: submitted to *Nature Communications*; a citation and DOI will be
added here on publication.

## Reproducibility status (read this first)

This is the actual pipeline code, not a cleaned-up reference reimplementation,
and it is **not click-to-run outside the source HPC environment**. Concretely:

- The `server_export/scripts/` and top-level `*.py` scripts were run against
  restricted-access CPTAC/TCGA/GDC data on the source institution's HPC
  cluster. Getting a given number to reproduce requires (a) obtaining that
  data yourself under the terms in the manuscript's Data availability
  statement, (b) installing `requirements.txt`, and (c) pointing each script
  at your own data layout via the environment variables documented in its
  header (most commonly `MORPHO_SIG_CSV`, `MORPHO_OUT`, `MORPHO_WSI_EMB_DIR`,
  `MORPHO_COHORT` and `MORPHO_ROOT`; grep a script's docstring for the full
  list it uses). Several
  scripts still default to the source institution's absolute path
  (`/public/home/fjhui/ZW`) when a variable is unset — that default will not
  resolve on another machine and is not meant to.
- The `figures/*.py` scripts are the easiest entry point for inspection:
  each one is a self-contained matplotlib script that reads a small,
  already-computed table (from `figures/figdata/` or `server_export/`,
  deposited separately as Supplementary Data, see below) and writes one
  published figure. No GPU or restricted data is needed to *read* them, only
  to *run* them once the corresponding data table is in hand.
- Exact package version pins used for the reported analyses are not yet
  captured in this repository (see the manuscript's Limitations); versions in
  `requirements.txt` are unpinned.

## What is (and isn't) in this repository

This repository holds the residual-estimation pipeline, the robustness and
replication analyses, and the figure-generation scripts referenced in the
manuscript's Methods and Code availability sections.

It does **not** hold: the manuscript source; raw whole-slide images, protein
spectra or sequence data (available from PDC/GDC/TCIA under their own terms,
see below); or the large derived tables and slide/patient manifests
referenced as Supplementary Data (deposited separately, e.g. via Zenodo,
alongside the accepted manuscript). It also omits the manuscript's own
LaTeX-editing scripts and third-party literature-metadata dumps, which are
not analysis code.

## Which script makes which published display item

The 13 figures actually embedded in the current manuscript, and the script
that produces each one:

| Manuscript figure | Script | Output |
|---|---|---|
| Fig. 1 (flagship discovery) | `figures/make_master_composite.py` | `Fig_master.pdf` |
| Fig. 2 (capacity/composition controls) | `figures/make_fig_controls.py` | `Fig_controls.pdf` |
| Fig. 3 (biology vignettes) | `figures/make_fig_biology.py` | `Fig_biology.pdf` |
| Fig. 4 (robustness composite) | `figures/make_composite2.py` | `Fig_master2.pdf` |
| Fig. 5 (clinical) | `figures/make_fig_clinical.py` | `Fig_clinical.pdf` |
| Fig. 6 (mechanism) | `figures/make_fig_mechanism.py` | `Fig_mechanism.pdf` |
| Methods pipeline diagram | `figures/make_fig_methods_pipeline.py` | `Fig_methods_pipeline.pdf` |
| Supp. Fig. (reproducibility) | `figures/make_fig_supp1.py` | `Fig_supp1_reproducibility.pdf` |
| Supp. Fig. (scores) | `figures/make_fig_supp2.py` | `Fig_supp2_scores.pdf` |
| Supp. Fig. (mTOR) | `figures/mtor_mechanism.py` | `Fig_mtor.pdf` |
| Supp. Fig. (pan-organ core) | `figures/make_fig_pan_organ_core.py` | `Fig_pan_organ_core.pdf` |
| Supp. Fig. (replication) | `figures/make_fig_replication.py` | `Fig_replication.pdf` |
| Supp. Fig. (WSI spatial) | `figures/make_fig_wsi_spatial_supp.py` | `Fig_wsi_spatial_supp.pdf` |

Some supplementary tables are produced by a matching `make_supptable_*.py`
script in the repository root (e.g. `make_supptable_nuisance.py` →
`SuppTable_nuisance.tex`); others, including the two main-text tables (atlas
and enrichment survival), are compiled by hand from the underlying result
files rather than by a dedicated script.

`figures/` contains roughly twenty additional `make_*`/plotting scripts
(e.g. `make_fig1.py`, `make_figs.py`, `make_realdata_figs*.py`,
`make_txome_fig.py`) that are **earlier drafts of these same figures, kept
for provenance and not referenced by the current manuscript** — several
carry an explicit `# superseded by ...` note. If a script's output filename
does not appear in the table above, it is not currently live; check the
table, not the script list, to find the source of a specific published
number.

## Layout

- `figures/*.py` — one script per main/supplementary figure or table (see
  the mapping above), reading precomputed result tables and writing the
  published figure PDFs. `figures/mrstyle.py` is the shared plotting style
  module every other figure script imports.
- `server_export/scripts/` — the residual-estimation, acquisition-robustness
  and cross-cohort replication pipeline, run on the source institution's
  HPC cluster against CPTAC/TCGA/GDC data.
- Top-level `*.py` — cohort-level analysis scripts (clinical association,
  composition controls, subsampling/ablation diagnostics, encoder
  comparisons, external-cohort replication, etc.) referenced by file name
  throughout the Methods; `make_supptable_*.py` scripts build the
  Supplementary Tables.

## Data

No data is redistributed here. The five discovery cohorts and their
Proteomic Data Commons study identifiers: CCRCC (PDC000127), LUAD
(PDC000153), UCEC (PDC000125), GBM (PDC000204), PDAC (PDC000270); matched
RNA-seq is from the Genomic Data Commons and H&E whole-slide images from the
corresponding TCIA *CPTAC Pathology Discovery Cohort* collections. Full
accessions, release dates and replication-cohort sources (TCGA-RPPA,
CPTAC-2 TCGA-BRCA, CPTAC-3 confirmatory cohorts) are listed in the
manuscript's Data availability statement.

## Requirements and reproducibility

`requirements.txt` lists the third-party packages imported across this
repository (collected by scanning every committed script; versions
unpinned, see Reproducibility status above). Random seeds, cross-validation
splits and foundation-model encoder weights used are documented in the
manuscript's Methods ("Software, versions, seeds and encoder weights").
Third-party model weights (UNI, Phikon, Phikon-v2) are not redistributed
here and remain subject to their own providers' licences.

## License

MIT — see [LICENSE](LICENSE). This licence covers the code in this
repository only; it does not extend to third-party model weights it loads.
