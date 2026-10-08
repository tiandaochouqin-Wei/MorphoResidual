# Post hoc BRCA external replication on frozen tissue sections, with a set-level readout: pre-specification (2026-10-08)

**Status.** Written 2026-10-08 (times UTC), before any frozen Tissue Slide (TS) image was
downloaded and before any set-level statistic was computed for either BRCA arm. The
published BRCA arm (CPTAC-2 BRCA iTRAQ4 proteome, PDC000173, with GDC STAR-Counts RNA and
TCGA Diagnostic Slides; Phikon, stride 0; 102 patients) reported 0 of 10,031 genes and
99.1% negative increments (main text, external replication). The 2026-09-11 plan (item C3
of `REMAINING_WORK_2026-09-11.md`) proposed a TS arm with a note, dated before the DX
result, that both arms would be reported; the DX result is now published, so that note can
bind only the TS arm and the set-level readouts defined here. On 2026-10-08 the
corresponding author asked for this analysis. Frozen by the sha256 of this file, the
scripts and the gene-set files (sidecar `POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256`) and
released publicly on GitHub before any residual-analysis job (`RUNBOOK.sh` step 6 onward) is
submitted. To save time, steps 0-5 (case list, tested-gene list, TS manifest, download, md5
check, feature extraction, slide QC) were started on the server on 2026-10-08 before freezing;
they read file lists, slide headers and images only and compute no protein- or mRNA-related
statistic.

## 0. What it asks, and scope

The DX images are FFPE diagnostic sections, not adjacent to the proteome aliquot (all DX
slides carry vial code 01Z). TCGA's frozen Tissue Slides are top/bottom sections of the
frozen block from which the analytes were taken (the proteome aliquot's vial carries a
primary-tumour TS slide for 100 of 102 patients). The TS arm asks whether the external null
reflects the diagnostic slides' distance from the analysed tissue. A set-level readout is
added to both arms because a per-gene FDR count at n = 102 is a weak test.

Post hoc; reported whatever the result; no change to any published number. Prior exposure:
none. No TS image exists locally or on the server; no set-level statistic has been computed
for the DX arm.

## 1. Data and arms

- **Cases.** The 102 patients of the published arm (intersection of PDC000173 tumour
  aliquots, GDC Diagnostic Slides and GDC STAR-Counts RNA; re-derived on 2026-10-08, the
  scripts assert equality with the published case list). RNA and protein inputs are the
  published ones, unchanged (including the two patients whose metastatic-sample RNA is
  averaged in, as published).
- **TS slides (primary rule).** Every GDC 'Tissue Slide' file with a primary-tumour sample
  code (01x) for these 102 patients: 189 slides (95 top, 86 bottom, 8 middle sections), 102
  of 102 patients, 59.3 GB, open access; metastatic (06A) and normal (11A) slides excluded.
- **Extraction and pooling.** Phikon, 256-px tiles at level 0, stride 0, the published
  tissue filter (`extract_phikon.py` with `--manifest`). Per patient: mean over tiles
  within a slide, then the equal-weight mean over the patient's slides, as in the published
  DX arm (`--pool slide`).
- **Arms.**
  - **TS (primary).**
  - **DX (reproduction).** The published DX arm re-run with the patched script. Gate: it
    must reproduce the published per-gene table exactly (10,031 tested, 0 significant,
    negative fraction 0.991) or the job stops; its matrix supplies the DX set-level readout.
  - **TS-tile (sensitivity).** TS with tile-weighted pooling (`--pool tile`), the CPTAC-3
    discovery convention.
- **Estimator and null.** Unchanged from the published arm (`cptac2_residual.py`): 20 PCs
  (full SVD, all patients), 5-fold `KFold(shuffle, random_state=0)`, closed-form ridge
  lambda = 1, B = 1,000 permutations with `RandomState(i)`, BH-FDR; significant = FDR < 0.05
  and increment > 0. The patched script also saves the genes x 1,001 increment matrix
  (column 0 observed) so that the set-level readout uses the same permutations.
- **Slide QC.** `ts_slide_qc.py` records microns per pixel, apparent magnification,
  scanner and tile counts per slide; reported descriptively, no stratification.

## 2. Readouts

- **Per gene (as published):** n, genes tested, significant, percentage, median increment
  of significant genes, fraction of tested genes with a negative increment.
- **Set level (new).** AUC of the BRCA increment of a discovery gene set against all other
  BRCA-tested genes (rank Mann-Whitney, the CPTAC-3 confirmatory statistic), one-sided p
  against the same 1,000 permutations; rho between the discovery score and the BRCA
  increment, reported descriptively.
  - **Primary set:** `pub_k3`, genes significant (FDR < 0.05, positive increment) in at
    least three of the five published CPTAC-3 discovery cohorts (777 genes). Breast is a new
    organ, so a cross-organ set is used; k >= 3 has the largest excess over the
    independence expectation (777 observed vs about 485 expected) among k = 1..5, while the
    170-gene core (k >= 4) is small.
  - **Secondary sets** (Holm within each arm, descriptive): `pub_k4` (170-gene core),
    `pub_k2`, `op_k2` (operator batch-corrected, >= 2 cohorts), `phi_k3` (Phikon discovery,
    >= 3 cohorts).
  - **Symbols:** exact matching between the CPTAC-3 discovery tables and the BRCA tables;
    if fewer than 50% of a set's genes are among the BRCA-tested genes the computation stops,
    and any remedy is a dated amendment made before any set-level statistic is read.

## 3. Readings and wording (tiers: p < 0.05 shown; 0.05 <= p < 0.15 suggestive; p >= 0.15 not shown)

Results, after the sentence "On iTRAQ, BRCA showed no morphology-predictable proteins: 0 of
10,031 tested, 99.1% with a negative increment.", exactly one of:

- **TS primary p < 0.05:** "Post hoc, on frozen sections from the tumour blocks of the
  same 102 patients, [s] of [t] genes reached FDR < 0.05 ([f]% with a negative increment),
  and the 777 genes significant in at least three CPTAC-3 cohorts ([k] tested) ranked above
  the other BRCA genes (AUC [a]; null [m] +/- [sd], p = [p]); on the diagnostic slides the
  same set gave AUC [a_dx] (p = [p_dx])."
- **TS primary 0.05 <= p < 0.15:** the same sentence with "ranked above" replaced by "were
  suggestively ranked above".
- **TS primary p >= 0.15:** "Post hoc, on frozen sections from the tumour blocks of the
  same 102 patients, [s] of [t] genes reached FDR < 0.05 ([f]% with a negative increment),
  and the 777 genes significant in at least three CPTAC-3 cohorts ([k] tested) did not rank
  above the other BRCA genes (AUC [a], p = [p]; diagnostic slides AUC [a_dx], p = [p_dx]);
  the external null is therefore not explained by the diagnostic slides' distance from the
  analysed tissue."

Abstract: unchanged unless the TS primary p < 0.05, in which case "a whole-proteome cohort
none." becomes "a whole-proteome cohort none, although on frozen sections from the analysed
blocks a cross-organ gene set ranked above background (post hoc)." Nothing else enters the
Abstract or Discussion.

Supplement: a table with every arm (TS, DX, TS-tile) x every set (primary and secondary,
raw and Holm p), the per-gene readouts, the slide QC summary, and the DX reproduction check.

## 4. Freeze and run order

1. Hashed into the sidecar and released on GitHub: this file; in
   `review/BRCA_TS_ARM_2026-10-08/`: `scripts/cptac2_replicate.py`,
   `scripts/cptac2_residual.py`, `scripts/extract_phikon.py`, `scripts/ts_slide_qc.py`,
   `scripts/brca_tested_genes.py`, the LSF wrappers, `RUNBOOK.sh`, `local_brca_setlevel.py`,
   and the gene-set files `setlevel/candidate_sets/pub_k3.tsv`, `pub_k4.tsv`, `pub_k2.tsv`,
   `op_k2.tsv`, `phi_k3.tsv`.
2. On the server, following `RUNBOOK.sh`: deploy to a separate scripts folder; build the
   TS manifest; download and md5-check the 189 slides; smoke-extract three slides; extract
   all; run the DX reproduction (gate) and the TS and TS-tile arms; run `ts_slide_qc.py` and
   `brca_tested_genes.py`.
3. Brought back: the three increment matrices, the per-gene tables and summaries, the
   pooling manifests, the slide QC table, the tested-gene list and the job logs (sha256
   recorded on arrival); then `local_brca_setlevel.py` is run once.
