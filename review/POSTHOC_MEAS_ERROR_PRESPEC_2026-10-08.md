# Post hoc own-mRNA measurement-error (regression-dilution) test of the beyond-mRNA signal: pre-specification (file date 2026-10-08, local time UTC+8)

**Dates and times.** The date in this file's name (file date 2026-10-08, local time UTC+8) is
the local date of the drafting session. Every time and date stated inside this file is UTC and
carries its UTC date. The derivation, the local analysis code and this text were drafted on
2026-10-07 between about 21:00Z and 23:59Z; a first complete version was saved on 2026-10-08
at about 00:05Z, and it was revised from about 00:20Z on 2026-10-08 after a cross-check of the
three post hoc pre-specifications of this file date. The sidecar's `# written` line records
when it was frozen. Server-side times below are the cluster's clock as printed in the run
manifests (cluster time is UTC+8, consistent with the archive headers and the local transfer
time); each is converted to UTC.

**Status.** Drafted after the published discovery results, the post hoc batch-axis
analyses (`POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md`), the published C1-UCEC and registered
C1-LUAD confirmatory readings and the alternative-selection-set re-read
(`POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md`), and before any statistic that relates
morphology to protein or mRNA beyond the published incremental R^2 was computed on project
data for this question by anyone. The corresponding author asked for the own-mRNA
measurement-error question of Limitations (v) to be tested (2026-10-07 UTC).

- **Export runs on the cluster, all before freezing, as allowed.** `measurement_error_checks.py`
  v2 (sha256 `66b2a7fb...`, unchanged, with `theta_kernel.py` `6c95cd11...` and the published
  `residual_analysis.py` `a54a0a56...`) was run through `lsf_measurement_error.sh`
  (`7bc94395...`) with `--B 0` (default). Every submission made for this analysis is listed
  here (`--out` relative to `/public/home/fjhui/ZW/`; production jobs submitted 2026-10-07 at
  about 23:25Z; their LSF job IDs were not recorded locally):

  | LSF job | Cohort | Host | `--out` | Exit | Reproduction verdict (max abs difference) | Used |
  |---|---|---|---|---|---|---|
  | 75737658 (smoke, `--max-genes 40`, before the production jobs) | UCEC | s002 | `meas_error/smoke_ucec` | 3 | CLOSE (2.15e-5) | no |
  | production (ID not recorded) | CCRCC | s002 | `meas_error/ccrcc` (wrapper default) | 0 | EXACT (0.0) | yes |
  | production (ID not recorded) | LUAD | s006 | `meas_error/luad` (wrapper default) | 3 | CLOSE (5.14e-5) | no |
  | production (ID not recorded) | UCEC | s006 | `meas_error/ucec` (wrapper default) | 0 | EXACT (0.0) | yes |
  | production (ID not recorded) | GBM | s004 | `meas_error/gbm` (wrapper default) | 3 | CLOSE (2.85e-5) | no |
  | production (ID not recorded) | PDAC | s004 | `meas_error/pdac` (wrapper default) | 0 | EXACT (0.0) | yes |
  | 75737669 (host rerun) | LUAD | s002 | `meas_error/luad_s002` | 0 | EXACT (0.0) | yes |
  | 75737670 (host rerun) | LUAD | s004 | not recorded locally | 0 | EXACT | no |
  | 75737671 (host rerun) | GBM | s002 | not recorded locally | 3 | CLOSE | no |
  | 75737672 (host rerun) | GBM | s006 | `meas_error/gbm_s006` | 0 | EXACT (0.0) | yes |

  Exit 3 is the script's exit on any verdict other than EXACT; CLOSE means a maximum
  difference above 1e-10 but at most 1e-4 from the published per-gene values, and the script
  prints "do NOT freeze until explained". Explanation: the differences are host-dependent
  numerics of the morphology PCA (same inputs and code; the same cohort reproduces exactly on
  another host). The four host reruns were submitted unchanged except for `--out` to find,
  for LUAD and GBM, a host whose numerics match the host class of each cohort's 2026-09-01
  published run; that is why the two folders used carry a host suffix. The host choice was made
  on the reproduction of published quantities only. Of the two exact LUAD reruns (s002 and
  s004), s002 was packed; both reproduced the published quantities exactly, the choice was
  made before any quantity of this analysis existed, and s002 is also the CCRCC host; the
  s004 export stays on the cluster and is not used. The five exports used are the EXACT ones:
  CCRCC (s002), LUAD (`luad_s002`), UCEC (s006), GBM (`gbm_s006`), PDAC (s004); they finished at
  2026-10-07T23:30:45Z, 23:37:37Z, 23:28:34Z, 23:37:39Z and 23:28:14Z (run manifests). With
  `--B 0` the script computes no statistic relating morphology to protein or mRNA except the
  **published** per-gene incremental R^2, which it re-derives from the exported arrays: verdict
  `EXACT (<=1e-10)` in all five used exports, max |difference| 0.0 for n, R^2_rna, R^2_both and
  incremental R^2 over 9,635 / 10,724 / 10,512 / 10,785 / 9,897 fit genes (CCRCC / LUAD / UCEC
  / GBM / PDAC). It also wrote RNA-only and protein-completeness descriptors (raw-count depth,
  per-gene replicate ICC, Poisson reliability bound, missingness, transcriptome-wide RNA PCs,
  batch labels), none of which involves morphology. Only the five used exports were packed into
  `export/meas_export_5cohorts.tgz` (folders `ccrcc/`, `luad_s002/`, `ucec/`, `gbm_s006/`,
  `pdac/`), transferred at about 2026-10-07T23:41Z and extracted to `M/export/<c>/` (the two
  suffixed folders as `luad/` and `gbm/`; file names and bytes unchanged). On 2026-10-07 at
  about 23:45-23:50Z this session ran `sha256sum -c export_sha256_<c>.txt` in each cohort
  folder: every line OK (re-checked 2026-10-08 about 00:25Z, with the two renamed folders'
  `inputs_*.npz` also hashed straight from the archive: identical).
- **Apart from the runs listed above and F2 (below), nothing has run on project data for this analysis.**
  `local_theta.py` has not been run on any exported file, not even `--repro-only` (no
  `local_theta/` folder exists in any cohort folder at the time of writing); the published-
  quantity check F1 is step 1 of the run order after freezing (section 10). The full run
  refuses real data until the sidecar below exists and matches.
- Frozen by the sha256 of this file and of every file listed in
  `review/MEAS_ERROR_2026-10-08/SIDECAR_FILES.txt` (section 10), written to the sidecar
  `review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md.sha256`, and released publicly on GitHub
  under the release rule of section 10
  before the first real-data invocation of `local_theta.py` (including `--repro-only`). The only
  real-data step allowed before freezing is F2 (side files; section 10), which reads no
  morphology-omics relation.
- **The rules are those of `theory/DERIVATION.md` section 9 ("Final reading rules")**, revised
  after an independent theory verification (TAG verify-theory; its report `review_verify.md`
  stays in the drafting session's scratch folder and is not released; the changes it caused are
  DERIVATION section 8.1, and its scripts and logs are released in `M/theory/final-rules/`) and
  calibrated on synthetic data only (DERIVATION sections 7.6 and 7.8). Where this file fills a
  gap the derivation leaves open, the text says so and takes the side less favourable to the
  paper's reading. Every remaining difference from DERIVATION section 9, with its direction,
  is listed at the end of section 7 ("Departures from DERIVATION section 9").

**F2 (side files) was run on 2026-10-08 at 00:43:29Z** (`M/side/build_side_files.py`, console output in `M/side/build_side_files.log`; a first invocation at 00:42:55Z stopped on the CCRCC purity sheet, which holds two `TumorPurity` rows, and was repeated after the script was corrected to use the first, RNA-seq ESTIMATE, row). It reads no morphology, protein or RNA values (only the GMT, the keyword map, the covariate tables and the `patients` array of each inputs file). `families_reactome2022.tsv` has 6,303 rows (genes per family: translation 1,908, secretion 2,065, splicing 250, ptm 1,535, folding 233, ecm 312; terms matched: 39, 61, 9, 30, 15, 15 of 1,816 GMT terms). Purity coverage, patients with a value of the inputs patients: CCRCC 103 of 103, LUAD 103 of 105 (not the 104 expected: C3L-01890 has no value in the source and C3N-00545 is not in it), UCEC 96 of 100, PDAC 97 of 137 (the source holds `epithelial_cancer_deconv` for 100 of 140 cases); GBM has no purity file.

## 0. Question, scope and prior exposure

**The question (main.tex Limitations (v)).** A single noisy RNA-seq measurement M = X + e
(reliability lambda = var X / var M) leaves part of the true transcript X in the protein
residual r_p = P - beta M. Morphology could predict that residual by proxying X ("regression
dilution", H_noise) rather than post-transcriptional variation u (H_post). Limitations (v)
currently states that a denoised or replicate-averaged mRNA baseline "would be the decisive
test and has not been run", and the Abstract says "with each transcript measured once, part
may be true expression, not post-transcriptional regulation". This analysis does not run a
denoised baseline. It tests a signature: under H_noise the morphology-only prediction of the
residual is aligned with the measured transcript, with alignment theta = lambda / (1 -
lambda) > 0; under H_post theta = 0; if morphology reads the mRNA error itself, or a
transcript component that does not reach protein, theta = -1 (DERIVATION sections 1-2).

- **Post hoc and reported whatever the result.** Nothing below was pre-registered. After
  freezing, no set, stratum, arm, seed, threshold, reading or sentence is added, dropped or
  changed; a later change is a dated amendment released before the affected run, and a run
  made with different hashes is reported as such.
- **No effect on published or registered readings.** Table 1 counts, the published C1-UCEC
  and registered C1-LUAD readings and every post hoc re-read stand as reported. Nothing here
  is pooled with them.
- **What the test can and cannot separate (fixed now and repeated in every sentence that
  reports a favourable reading).** theta separates "enters protein through the measured
  transcript, linearly" from "does not". It cannot separate biological post-transcriptional
  variation from a protein-side technical artefact that morphology can read (for example
  interference or reference-pool effects shared within a TMT plex): such an artefact enters
  r_p exactly as u does (DERIVATION section 3 item 4, section 6). It does not exclude a
  non-linear transcript contribution: a non-linear mRNA-to-protein map with lambda = 1 gives
  theta = 0.00-0.03 (DERIVATION section 6). The spline sensitivity (section 6) is a partial
  check of the second point; nothing here checks the first beyond the descriptive within-plex
  arm (section 3). The fixed qualifier Q of section 8.5 carries this into every favourable
  sentence.

**Prior exposure (disclosed).** Known to the session that wrote this file when it was
written:
1. All published per-gene results (incremental R^2, R^2_rna, permutation p, FDR), Fig. 1f
   (R^2_rna against increment), the held-out seed-reliability increments, the
   transcriptome-wide baseline (median 40% of the noise-adjusted advantage retained), the
   CCRCC purity-adjusted increments (r = 0.85 with the unadjusted, 96% of the median
   retained), the composition controls, and the plex, medium and operator batch-axis counts.
2. **The 2026-09-27 dilution-bound analysis**, recovered byte-identical to
   `review/recovered_2026-09-27_dilution_bound/` (generated 2026-09-27 about 03:49-04:05Z by
   agent `a3b5cce95bf975a94`, workflow `wf_d5422fb4-b59`; never used in the manuscript). It
   used only published per-gene tables and, for its simulations, the real local UNI slide
   embeddings with **synthetic** X, M and P (protein depending only on true mRNA, morphology a
   perfect proxy for it). It computed no alignment, sign test or any statistic relating
   morphology to measured protein or mRNA. Its numbers (final version) are recorded only in
   that agent's session log (lines 191/195/199 of its jsonl), not in the released folder:
   lambda* = R^2_rna / (R^2_rna + held-out increment), the largest reliability
   at which perfect-proxy dilution could produce a gene's increment, had medians 0.312 /
   0.705 / 0.384 / 0.785 / 0.511 over the significant sets (CCRCC / LUAD / UCEC / GBM /
   PDAC; p25-p75 0.000-0.736, 0.281-0.895, 0.000-0.753, 0.451-0.925, 0.000-0.826), with
   59.5 / 34.5 / 57.1 / 28.5 / 49.3% of significant genes requiring lambda <= 0.5;
   restricted to the bottom seven deciles of R^2_rna the medians were 0.061 / 0.536 / 0.190
   / 0.636 / 0.036. 63 / 46 / 74 / 46 / 60% of significant genes had a held-out increment
   above their own perfect-proxy ceiling at lambda = 0.90. In the semi-synthetic simulation
   the held-out median that dilution alone manufactured crossed the observed one at lambda
   about 0.90 (LUAD, GBM), 0.80-0.85 (CCRCC, PDAC) and 0.75 (UCEC). Its own README records
   that these crossovers do not exclude dilution (single-library reliability of low-count
   genes can lie in that range) and that its decile-shape argument is confounded by depth.
3. `theory/DERIVATION.md` and all its simulations (`sim_*.py`, `section7_tables.md`,
   `theory/verify/`, and the verification and post-review checks of the drafting session's
   scratch folders `phaseC/verify-theory/` and `phaseC/final-rules/`, copied byte-identical to
   `M/theory/final-rules/`): **synthetic data only**. The one external input is public GDC file
   metadata (numbers of RNA files per case, queried 2026-10-07 UTC), which shows that
   multi-file CPTAC-3 cases are almost always distinct tumour samples of one case, not
   re-runs of one aliquot. Every threshold below was calibrated on these simulations before
   any exported file existed.
4. The export (Status and section 1): the reproduction verdicts, maximum differences, exit
   codes and hosts of the ten jobs listed in Status, as reported to this session from their
   LSF output; this session read the five `repro_check_<c>.json` and `run_manifest_<c>.json`
   files of the used exports (published-quantity reproduction, label coverage, shapes, hosts,
   hashes), listed the array names inside each `inputs_<c>.npz` (names only; no array was
   loaded at that point) and read the column headers of `gene_table_<c>.csv`, `patient_table_<c>.csv` and
   `replicate_icc_<c>.csv`. While reading the last header it printed, unintentionally,
   the first data row of `replicate_icc_ccrcc.csv` (gene A1BG: ICC 0.39, 40 cases with
   replicates, mean 3.3 files per case). No set-level ICC, depth, missingness or reliability
   summary has been computed. Any other inspection of the exported files by anyone before
   freezing must be added to this item before the sidecar is written. Added under that
   rule: F2 (Status) loaded the `patients` array (patient identifiers only) of each
   `inputs_<c>.npz` to count purity coverage, and a parse check of the F2 outputs loaded
   the same arrays; no other array of an exported file was loaded.
5. C1-LUAD confirmatory results (AUC 0.644 on the 2,310 set; 726 set under the
   operator-stratified null) and the C1-UCEC re-reads are known; no measurement-error
   statistic exists for either confirmatory cohort.

No alignment statistic (N, N_c, theta, S, z_N), W-only prediction of a protein residual, or
any other statistic relating morphology to protein or to mRNA beyond the published increment
has been computed on project data by anyone.

## 1. Data

All paths relative to `MorphoResidual_paper/`; folder `M = review/MEAS_ERROR_2026-10-08/`.
Per cohort `<c>` in {ccrcc, luad, ucec, gbm, pdac}, the used export is in `M/export/<c>/`
(Status):

| File | Content |
|---|---|
| `inputs_<c>.npz` | patients x genes arrays: `mrna_log2tpm1`, `protein` (NaN = missing), raw `counts`, `tpm`; the published 20 morphology PCs `pcs` (PCA fitted once on all common patients, as published) and `wsi_pooled`; `rna_pcs20` (transcriptome-wide RNA PCs as in `transcriptome_baseline.py`); `label_<axis>` for operator, plex, quarter, scanner, mpp, medium, cohort_stream; `n_rna_files`; `icc_log2tpm1`; the published per-gene results `pub_*`; `is_synthetic = False` |
| `gene_table_<c>.csv` | per gene: depth, completeness, ICC, lambda_Poisson, published values, the export's re-derivation of the published increment |
| `patient_table_<c>.csv` | per patient: library size, STAR counters, number of RNA files |
| `replicate_icc_<c>.csv`, `replicate_files_<c>.npz` | CCRCC and UCEC only: per-gene one-way ICC(1) of log2(TPM+1) across a case's RNA files, and the per-file values from which `local_theta.py --check-icc` recomputes it |
| `repro_check_<c>.json`, `run_manifest_<c>.json`, `export_sha256_<c>.txt` | reproduction verdict, run record, transfer hashes |

| Cohort | Patients | Fit genes | Published significant set | Operator labels | RNA replicate cases |
|---|---|---|---|---|---|
| CCRCC | 103 | 9,635 | 2,191 | 103/103 | yes |
| LUAD | 105 | 10,724 | 2,310 | 104/105 | no |
| UCEC | 100 | 10,512 | 2,566 | 99/100 | yes |
| GBM | 99 | 10,785 | 702 | 99/99 | no |
| PDAC | 137 | 9,897 | 1,572 | 137/137 (scanner 137/137) | no |

**RNA replicates are distinct tumour samples.** Where a case has several RNA-seq files
(CCRCC and UCEC; the export's ICC table counts the cases), public GDC metadata show these are
different tumour samples of the case, not re-sequencing of one library, so the ICC counts
intra-tumour sampling variation as well as library noise. The published pipeline averages a
case's TPM over its files, so for those cases the analysed M is more reliable than one file.
In both respects the single-file ICC is a **lower bound** on the reliability lambda of the
analysed M (DERIVATION section 3 item 7), which is the direction every rule below needs.
LUAD, GBM and PDAC have no replicate design; their lambda_lo is transferred (section 5) and is
an assumption, not a measurement.

## 2. Statistics (per gene, then per set)

Exactly as implemented in `local_theta.py` (sha256 in section 10) on top of `theta_kernel.py`
(byte-identical to the server copy). Per gene g, on its fit patients (non-missing protein),
with M and P standardised over the fit patients so that beta-hat = delta-hat = rho-hat:
r_p = P - rho-hat M; r_m = M - rho-hat P; r_p-hat = W-only out-of-fold prediction of r_p from
the 20 raw morphology PCs with the published estimator (5-fold `KFold(shuffle, random_state=0)`
on the gene's fit patients, in-fold standardisation, closed-form ridge alpha = 1), never the
nested increment d_p.

- N_g = rho-hat cov(r_p-hat, M) (numerator of theta); D_g = cov(r_p-hat, r_p);
  S_g = sign(rho-hat) corr(r_p-hat, r_m) (the sign test, reported only).
- **Tested statistic:** N_c,set = sum over the set of sign(rho-hat_g) corr(r_p-hat_g, M_g),
  computed as the kernel's `Nc` (= rho-hat corr) divided by |rho-hat|, for the observed value
  and every null draw alike (rho-hat depends only on (M, P), which the null preserves). The
  rho-hat-weighted sum has a positive bias proportional to G var(rho-hat) under H_post and is
  reported only as the diagnostic column `Nc_rhow_*`.
- **theta of record:** theta_pc = (sum N - E0 sum N) / (sum D - E0 sum D), E0 = mean over the
  primary arm's null draws; lambda-hat = theta_pc / (1 + theta_pc). Raw theta_set = sum N /
  sum D is a sensitivity.
- **Interval:** patient bootstrap, 200 resamples (`--boot-seed 20261008`), one resample of
  the n patients shared by all genes, each patient keeping its fold in the gene's published
  KFold (integer weights equal duplicated rows), of (sum N_b - E0 sum N) / (sum D_b - E0 sum
  D) with **E0 held fixed at the observed-sample null mean** (a convention; see Departures);
  95% percentile interval and bootstrap SE.
- **Set-level z:** z = (observed set sum - null mean) / null SD over the B draws, for D (z_D),
  N_c (z_N) and S (z_S); permutation p-values are reported beside them. Thresholds apply to
  z, not to p (p cannot go below 1/201). A gene with a non-finite observed or null value in
  an arm is dropped from that arm's set sums (the script's rule; the count is in the output).
- **Per gene:** only the proportions of genes whose N_c lies below / above the 2.5 / 97.5%
  band of their own null (2.5% each expected); no gene-level calls.

## 3. Nulls

Joint patient<->slide permutation of the PC rows: one draw per permutation shared by all
genes and all targets, folds, standardisation and alpha fixed; B = 200; **fresh draws**,
`RandomState(20261007)` (`--seed 20261007`, the script default), not the `RandomState(0)`
sequence that produced the published selection. Arms:

- **Primary: within operator** (`aperio.User` from the slide headers, the labels the
  sitepack analyses used): patients permuted only within operator, levels in sorted order,
  singletons fixed; an unlabelled patient forms its own `_missing` stratum, which with one
  such patient (LUAD, UCEC) is a fixed singleton. **PDAC: within operator x scanner**
  (`operator_scanner`).
- **Sensitivity: unrestricted** (`none`), the published null's construction.
- **Descriptive only: within TMT plex** (`plex`, every cohort). No reading is taken from it.
  It is run because a plex-shared protein artefact read by morphology through between-plex
  differences would be absorbed by it; an artefact acting within plexes is not excluded by it,
  so no outcome of this arm separates biology from protein-side batch.
- Nothing is residualised on operator for theta (residualising biases theta upward under
  H_post; DERIVATION section 4.3).

**phi_op** (reported for every read set): (E0_within[sum D] - E0_unres[sum D]) /
(sum D_obs - E0_unres[sum D]), the share of the permutation-centred W-only signal that the
within-operator null keeps in its centre, i.e. attributable to between-operator structure;
the same for the nested increment Delta_p.

## 4. Gene eligibility, sets and strata

- **Eligibility (missingness).** Fitting each gene on its non-missing patients selects on P,
  a collider, and pushes theta and z_N negative even under H_post (5% truncation gave z_N of
  -3.5 / -4.4 at n_sel 65-95; DERIVATION section 3 item 1). The primary analysis therefore
  uses **complete genes only** (n_fit = n, no missing protein value in the cohort).
- **Primary set per cohort:** `pinned_sig` = published FDR < 0.05 and published incremental
  R^2 > 0, restricted to complete genes. If it has fewer than 20 complete genes it is
  labelled **MNAR-exposed** and read on `pinned_sig|miss<=5` (complete genes plus genes with
  0-5% missing protein); a negative outcome there goes to R4w, never R4. The complete-gene
  `pinned_sig` row (if it has >= 5 genes) is then the "complete-gene set" of rule 3.
- **Strata of each base set** (base sets need >= 5 genes for z; strata are reported only with
  >= 15 genes; theta is read only with >= 20): depth tertiles `|depthT1..3` (tertiles of the
  median raw count over the fit patients, cut on the base set's complete genes, complete genes
  only); missingness strata `|miss0-5` (0 < f <= 5%), `|miss5-20` (5 < f <= 20%), `|miss20+`
  (f > 20%) and `|miss<=5` (complete plus 0-5%).
- **Families (`--sets-file`).** Six signature families, the keyword map of
  `figures/enrichment_c2_levels.py` (`SIGNATURE_FAMILIES` with folding and ECM, and the
  `\becm\b` rule of `families_of`) applied to the term names of
  `figures/figdata/Reactome_2022.gmt`: family f = union of the genes of every term whose
  lowercased name matches f (a gene may belong to several families). Written as
  `M/side/families_reactome2022.tsv` (columns `set`, `gene`; section 10, F2). The script then
  forms `pinned_sig&<family>` (read) and the unselected `<family>` (descriptive), each on
  complete genes, each with its own strata.
- **Read strata (used by rules 6, 9 and 10).** Within a cohort: the depth tertiles of
  `pinned_sig`, the `pinned_sig&<family>` sets, and the depth tertiles of each
  `pinned_sig&<family>` set (the family x depth crossing counts), each with >= 20 genes so that
  theta is read. Missingness strata are the MNAR diagnostic and are not read (the `|miss<=5`
  row of an MNAR-exposed `pinned_sig` is that cohort's primary set, not a stratum); the
  descriptive sets are not read. **k** of a cohort = the number of its read strata, counted
  from its `theta_<c>_reading.csv` (rows of these kinds with `n_genes` >= 20).
- **Selection strength** s = n_sel / G_tested of the set or stratum itself: complete pinned
  genes in it over complete published-tested genes satisfying the same definition.
- **Descriptive only:** `all` (all complete fit genes) and the unselected family sets: z_D,
  z_N, theta are reported, no rule is applied (null genes make theta unstable there).
- **Not analysed:** the batch-corrected (operator, plex, medium) discovery sets of Table 1.
  Their selection used residualised components, and the script's selection-strength
  definition applies to subsets of `pinned_sig` only. The reading concerns the uncorrected
  significant sets; every sentence says so.

## 5. Reliability band, selection strength and kappa

- **lambda_lo** (lower bound on lambda). CCRCC and UCEC: the **median** per-gene replicate
  ICC(1) of log2(TPM+1) over the set or stratum (`median_icc`). LUAD, GBM and PDAC
  ("transferred", flagged in every sentence): for the whole set and for every non-tertile
  row, lambda_T = the minimum of the available values of the median ICC of CCRCC
  `pinned_sig` and of UCEC `pinned_sig`, passed as `--lambda-lo-transfer`; for depth tertile
  j, the minimum of the available values of the two cohorts' `median_icc` of
  `pinned_sig|depthT<j>`, passed as `--lambda-lo-transfer-tertiles` (tertile j matched to
  tertile j by within-set rank, not by absolute count). These values are read from the
  CCRCC and UCEC `theta_<c>_reading.csv` (the `median_icc` column does not depend on the arm)
  and written with four decimals into `M/export/lambda_transfer.txt` before the LUAD, GBM and
  PDAC runs. **Availability.** A cohort contributes only if its run completed (F1 passed and
  the run wrote its reading file). If only one of CCRCC and UCEC contributes, its value is used
  alone and flagged "transferred from one cohort" in every sentence; a tertile value missing
  in both is replaced by lambda_T, flagged. If neither contributes, both flags are omitted and
  rho-bar^2 (the script's fallback; a rigorous but weak bound, since rho^2 <= lambda) is used
  and flagged. The minimum of two is the conservative choice: a lower lambda_lo lowers
  theta_noise,min and raises the f_X bound. Values are clipped to [0.001, 0.99].
- **lambda_hi** (upper bound): the median per-gene Poisson bound
  lambda_P = 1 - mean_i(1/(c_i + 0.5)) / var_i(ln(c_i + 0.5)) over the fit patients, clipped
  to [0, 1]. It ignores library overdispersion and is loose. If lambda_hi < lambda_lo for a
  set, the f_X range is reported as its upper end only and flagged "bounds inconsistent";
  no reading uses lambda_hi.
- **kappa(s, lambda_lo)**, the winner's-curse floor of theta under pure H_noise after
  selection (minimum over seeds and effect sizes of theta_set,selected / theta_population,
  less 0.05, rounded down to 0.05; synthetic; bands lower-inclusive):

| s = n_sel / G_tested | lambda_lo <= 0.55 | 0.55 < lambda_lo <= 0.65 | lambda_lo > 0.65 |
|---|---|---|---|
| >= 50% | 0.65 | 0.65 | 0.55 |
| 20-50% | 0.55 | 0.45 | 0.35 |
| 10-20% | 0.45 | 0.40 | 0.25 |
| 5-10% | 0.40 | 0.35 | 0.25 |
| 2-5% | 0.35 | 0.35 | 0.25 |
| < 2% | 0.30 | 0.30 | 0.25 |

  kappa = 1 for unselected sets and for the confirmatory cohort (`--confirmatory`).
- **theta_noise,min** = kappa lambda_lo / (1 - lambda_lo): the smallest selected-set theta that
  pure dilution at reliability lambda_lo produced in simulation.
- **X-share range** f_X in [(1 - lambda_hi) theta_pc / lambda_hi, (1 - lambda_lo) theta_pc /
  (kappa lambda_lo)] intersected with [0, 1]; post-transcriptional share lower bound
  1 - upper end. **Added here (more conservative than DERIVATION):** for R1 the sentence quotes
  f_X,max = (1 - lambda_lo) max(theta_pc,CI-hi, 0) / (kappa lambda_lo), clipped to [0, 1], the
  largest dilution share compatible with the upper end of the interval (DERIVATION's R1 blind
  spot is a mixture with theta_pop < 0.1).

## 6. Sensitivities and diagnostics

Computed for every set with the primary arm's permutations; they change a reading only where
section 7 says so.

1. **Unrestricted null** (section 3); the rules of section 7 are also applied to its z_N,
   theta_pc and interval (from `theta_<c>_sets.csv`, arm `none`) for the operator rule.
2. **Purity:** M and P residualised (in-sample OLS with intercept, per gene, on fit patients
   with a purity value; patients without one leave the fit set) on a per-patient purity
   before rho-hat, r_p and r_p-hat are formed. Purity files (`M/side/purity_<c>.tsv`,
   columns `patient`, `purity`; section 10, F2), the sources used by
   `composition_pathologist.py`: CCRCC ESTIMATE `TumorPurity` (the first, RNA-seq-based, of the sheet's two `TumorPurity`
   rows, Status, F2; `subtype_ccrcc_raw.xlsx`,
   sheet "ESTIMATE scores", tumour aliquots mapped to cases by
   `server_export/scripts/aliquot_to_case_tumor.tsv`, mean per case); UCEC `Purity_Cancer`
   (`subtype_ucec.txt`, tumour rows); LUAD `TUMOR_PURITY_BYESTIMATE_RNASEQ`
   (`composition_luad_cbioportal.csv`, Gillette et al.; 104 of 105 cases expected, to be
   confirmed by F2, whose recorded count governs); PDAC
   `epithelial_cancer_deconv` (`subtype_pdac_raw.xlsx`, sheet "Molecular_phenotype_data",
   the CPTAC deconvolution; no ESTIMATE purity is available there). **GBM has no purity
   estimate** (only xCell stroma/immune scores): the variant is not run and the supplement
   says so.
3. **Spline baseline:** P ~ M replaced by a natural cubic spline with 3 df fitted in-sample on
   the gene's fit patients, knots at their minimum, terciles and maximum (a convention);
   theta^spline = cov(r_p-hat, f-hat(M)) / cov(r_p-hat, r_p^spline), centred on its own null
   mean, its z_N, and Delta_p^spline (published ridge on [spline basis] against [spline basis,
   PCs]) beside the published Delta_p.
4. **20-RNA-PC baseline:** P ~ M + `rna_pcs20` (in-sample OLS); own-transcript component
   b_M M_perp (M residualised on the RNA PCs); N, D, theta, z_N as in the linear case. This is
   the "shared versus gene-specific" discrimination of R4: an artefact shared across genes is
   absorbed by the RNA PCs, a gene-specific one is not.
5. **Raw theta_set** (not permutation-centred).
6. **Missingness strata** (section 4), the MNAR diagnostic: a monotone fall of z_N across
   0% > 0-5% > 5-20% > >20% (strata present) is the signature of selection on completeness.
7. **Diagnostics, reported only:** Wm = corr(W-only prediction of M, M) within fit sets
   (set mean, z); the sign test z_S; the rho-hat-weighted z; per-gene band proportions; arm
   information (levels, fixed singletons). Wm was underpowered in simulation and is not a
   decision quantity.

## 7. Reading rules (numeric; applied by hand to `theta_<c>_reading.csv`, which holds numbers only)

Constants: z_D gate 3; z_N thresholds -3 (anti-alignment), 2 (alignment; one-sided nominal
2.5% under a normal null, realised H_post rate about 2-5% in simulation) and 3 (strong;
2 <= z_N < 3 is "weak"); R1 band |theta_pc| <= 0.05 for sets of >= 100 genes, <= 0.10 for
20-99; R4 cut theta_pc <= -0.3; heterogeneity 2 bootstrap SE. All z, theta and intervals are
from the primary arm unless stated. Every comparison with a threshold uses the unrounded value
(see Formats). Applied to each cohort's primary set, then to each read stratum (section 4), in
this order:

1. **R0 / R6.** z_D < 3 => R6, "no W-only signal on this set"; stop and check the
   reproduction record. Also R6: sum D - E0 sum D <= 0.
2. **R4.** z_N <= -3 under both the primary and the unrestricted arm and theta_pc <= -0.3 =>
   anti-aligned with the measured transcript. Sub-reading by the 20-RNA-PC variant:
   "shared" if its z_N > -3, "not shared" if its z_N <= -3 and its theta_pc <= -0.3,
   otherwise "unresolved". *Gap filled here:* z_N <= -3 and theta_pc <= -0.3 under the primary
   arm but z_N > -3 unrestricted is also read R4, with the qualifier "within-operator only".
3. **R4w.** z_N <= -3 and -0.3 < theta_pc < 0 => weak anti-alignment, undecided
   (missingness, an anti-aligned channel or operator structure); never H_err. Resolution, in
   order: (a) if the set is MNAR-exposed (read on `|miss<=5`), the complete-gene set's z_N
   (section 4; available if it has >= 5 genes): z_N > -3 there => "confined to proteins with
   missing values"; otherwise, or if it is not available, (b); (b) the 20-RNA-PC variant:
   z_N > -3 there => "absorbed by the transcriptome-wide baseline"; otherwise "unresolved".
   An unresolved R4w on a set that is itself the complete-gene set (not MNAR-exposed) adds
   "persisted in fully quantified proteins"; that phrase is used in no other case.
4. **R1.** -3 < z_N < 2 and |theta_pc| within the band => not aligned with the measured
   transcript.
5. **R7.** -3 < z_N < 2 and |theta_pc| outside the band => undecidable (D noise); interval
   reported, no reading. *Gaps filled here:* z_N <= -3 with theta_pc >= 0, and z_N >= 2
   with theta_pc <= 0, are also R7.
6. **z_N >= 2 and theta_pc > 0:** **R3** if theta_pc < theta_noise,min and the
   interval's upper end < theta_noise,min; **R2** otherwise (including theta_pc >=
   theta_noise,min, and an interval that contains theta_noise,min). A weak call
   (2 <= z_N < 3) is interpreted only if theta_pc exceeds the R1 band, or the call replicates
   (DERIVATION section 9: "in another stratum or cohort"): any other read stratum of the same
   cohort (section 4, including the `pinned_sig&<family>` sets and the primary set when the
   call is on a stratum), or any other cohort's primary set, has z_N >= 2 under the primary
   arm. *Gap filled here:* an uninterpreted weak call is **R1w**, "weak alignment of
   negligible size, not interpreted".
7. **Operator rule (DERIVATION section 9, items 2 and 6).** For every set or stratum that
   receives a reading, **whenever the unrestricted arm's z_N >= 2 or <= -3**, an operator
   clause is written. Rules 1-6 are applied to the unrestricted arm's numbers (z_D, z_N,
   theta_pc and its interval from `theta_<c>_sets.csv`, arm `none`; kappa, lambda_lo and
   theta_noise,min unchanged), and the clause reports that unrestricted reading, its z_N and
   phi_op (D and Delta_p). If the primary arm's z_N is beyond the same threshold in the same
   direction, the clause is "same direction" (the unrestricted null agrees; phi_op is the
   share between operators). Otherwise it is "operator-mediated transcript alignment"
   (unrestricted z_N >= 2) or "operator-mediated RNA artefact" (unrestricted z_N <= -3).
   Neither form is a gene-level H_noise or H_err reading, and neither changes the primary
   reading.
8. **Arbitration (C1-LUAD).** Active only if the confirmatory export of section 10, F4 is
   listed in the sidecar at freezing; it cannot be activated after the first real-data run.
   If active, the confirmatory reading (kappa = 1, `--confirmatory`, lambda_lo transferred as
   for discovery LUAD) of the discovery-pinned 2,310 and 726 sets decides R2 versus R3 for
   LUAD's discovery primary set; the discovery reading is reported beside it.
9. **R5 (heterogeneity), a reporting rule.** Within a cohort, over read strata (section 4):
   pairs of depth tertiles of one parent (`theta_<c>_hetero.csv`, kind `depth`) and pairs of
   `pinned_sig&<family>` sets (computed the same way from `reading.csv`): if two strata fall
   in different rules, or |theta_pc,a - theta_pc,b| / max(SE_a, SE_b) > 2, the strata are
   reported separately and not pooled. Missingness strata are not part of R5; they are the
   MNAR diagnostic (section 6, item 6). R5 makes no inferential claim: the sentence states the
   number of comparisons and that they are not corrected.
10. **Which readings reach the main text.** The cohort's primary-set reading. A stratum
   reaches the main text only if it reads R2 or R4 with |z_N| >= 3 while the cohort's
   primary reading is neither; every other stratum, the descriptive sets and the plex arm go
   to the supplement. Its sentence names k (section 4).
11. **Sensitivity riders** (they never change the rule, only add a clause):
   - R1: "spline-stable" if theta^spline_pc is within the R1 band and its z_N < 2,
     otherwise "spline-sensitive".
   - R2 / R3: "carried by purity" if the purity variant's theta_pc is within the R1 band and
     its z_N < 2; "purity not available" in GBM.
   - Any reading: the missingness-stratum pattern (section 6, item 6).

**Departures from DERIVATION section 9** (every remaining difference, with its direction
relative to the paper's reading; "less favourable" means more readings that qualify the
beyond-mRNA claim):
1. Operator rule and weak-call replication: **reverted** to DERIVATION's broader triggers in
   this version (rules 6 and 7). Remaining difference: the operator clause is also written
   when both arms agree ("same direction"), which DERIVATION does not ask for. Direction: less
   favourable or neutral (more reporting).
2. R4 under the primary arm only (rule 2, gap filled): DERIVATION requires both nulls.
   Direction: less favourable (more R4 readings; R4 changes the Abstract).
3. R4w resolution (rule 3): DERIVATION's "complete-gene stratum" is made operational as the
   complete-gene set of an MNAR-exposed cohort (>= 5 genes, z only), since such a cohort has
   by definition fewer than 20 complete genes. Direction: neutral (R4w is "not decisive"
   whatever the branch).
4. R7 gap fills and R1w (rules 5 and 6): DERIVATION is silent on z_N <= -3 with
   theta_pc >= 0, z_N >= 2 with theta_pc <= 0, and on how an uninterpreted weak call is
   worded. Direction: neutral (no reading is taken; R1w has its own clause and is never
   grouped with R1).
5. Arbitration (rule 8): DERIVATION has the confirmatory reading decide R2 versus R3 "for
   those sets" (the 2,310 / 726 sets); here it decides LUAD's discovery primary set.
   Direction: two-sided (it can move LUAD either way). Inactive unless F4 is frozen; F4 is not
   written.
6. lambda_lo: DERIVATION's "set-level ICC" is implemented as the median per-gene ICC. Direction
   not known before the run; no alternative summary is computed before or after.
7. Transfer: DERIVATION matches "the raw-count depth stratum"; here tertiles are matched by
   within-set rank, and whole sets and non-tertile rows take the minimum of the two cohorts'
   `pinned_sig` medians. Direction not known before the run (it depends on how the cohorts'
   depth distributions compare); the minimum of two cohorts is the less favourable choice
   within this construction.
8. lambda_hi: 1/(c + 0.5) replaces DERIVATION's 1/count (undefined at zero counts). This raises
   lambda_hi and lowers the reported lower end of the f_X range. Direction: favourable in that
   reported number only; no reading uses lambda_hi.
9. f_X,max for R1 (section 5): added. Direction: less favourable.
10. Bootstrap interval with E0 fixed (section 2): DERIVATION does not specify how E0 enters.
   Holding it fixed omits the null mean's Monte Carlo uncertainty (B = 200), so the interval
   is somewhat narrower. Direction: favourable for R3 against R2 (R3 needs the upper end below
   theta_noise,min). Declared, not corrected.
11. Selection strength on complete genes (section 4): DERIVATION's s is "of the cohort's
   pinned set"; here it is complete pinned over complete tested genes. Direction not known
   before the run.
12. Strata: no three-way family x depth x missingness crossing; depth tertiles on complete
   genes only; missingness strata per base set. Direction: fewer strata, so fewer chances for
   a rule-10 promotion (favourable) and for an R5 split (neutral).
13. R5 restricted to depth-tertile pairs within a parent and family-set pairs. Direction:
   fewer heterogeneity reports; no reading changes.
14. Rule 10 (main text) is filled here; DERIVATION is silent. Direction: neutral (it decides
   placement, not readings).
15. Sentences: the supplement sentences of section 8.5 replace DERIVATION's wording. They add
   the fixed qualifier Q to every favourable sentence and drop DERIVATION's "therefore not an
   artefact of the gene's own mRNA measurement error" and "independent of the linear
   contribution". Direction: less favourable. DERIVATION's sign-test sentences are not
   written; z_S is reported in the table only. Direction: two-sided (one favourable and one
   unfavourable sentence dropped).

**Multiplicity.** Five primary readings are made with no correction. Under H_post, about one
weak call in eight runs of five cohorts is expected; the weak-call rule above, the strong
threshold of 3 for the main-text stratum rule and the -3 threshold for anti-alignment are the
protection. No reading is combined across cohorts into a pooled test.

**Formats.** theta, interval, lambda, kappa, theta_noise,min, f_X and phi_op to 2 decimals;
z to 1 decimal; gene counts as integers; cohort order CCRCC, LUAD, UCEC, GBM, PDAC.
Readings use unrounded values; if rounding would print a value on the other side of a
threshold, decimals are added until it does not (e.g. p = 50/1,001 printed as 0.04995;
z = 1.96 not "2.0").
Here: a z of 1.96 or 2.96 is printed 1.96 / 2.96, never "2.0" / "3.0"; a z of -2.96 is printed
-2.96, never "-3.0"; theta_pc = 0.054 against the 0.05 band is printed 0.054, never "0.05";
the same for theta_pc against -0.3, against theta_noise,min and for an interval end against
theta_noise,min; decimals are added until the printed value is on the correct side.

## 8. Manuscript wording, fixed now (filled from the outputs after the run, never before)

`[names]` = the cohorts in that class in the fixed order; numbers per cohort in the same
order, separated by semicolons. Only clauses whose class occurs are written. `[date]` = the
UTC date of the sidecar's `# written` line, written as "D Month 2026".

### 8.1 Results: a new final paragraph of `\S\ref{sec:txome}` (after "...invisible to transcriptomics.")

Fixed opening, every outcome:
"Post hoc, in an analysis specified on [date] in a dated file released publicly before it ran
(Supplementary Note~\ref{snote:measerr}), we asked whether the signal could be the gene's own
mRNA measurement error. If morphology recovered true expression that a single RNA-seq
measurement misses, its prediction of the protein residual would be aligned with the measured
transcript, by an amount $\theta$ that grows with the transcript's reliability $\lambda$
($\theta=\lambda/(1-\lambda)$ when the whole signal passes through the transcript, $0$ when
none does). On each cohort's fully quantified significant proteins we estimated $\theta$
against a within-operator permutation null and compared it with the smallest value that
dilution produced in simulation at the reliability these genes show between RNA-seq libraries
of distinct tumour samples of the same patient (a lower bound; measured in CCRCC and UCEC and
assumed to transfer to LUAD, GBM and PDAC)."

[If a cohort is MNAR-exposed, "fully quantified significant proteins" is followed by
" (in [names], significant proteins with up to 5\% missing values, because fewer than 20 were
fully quantified)".] [If lambda_T came from one cohort: "measured in [CCRCC/UCEC] and assumed
to transfer"; if from neither: "bounded by the squared mRNA-protein correlation where no
replicate libraries were available".]

Then, in this order, the clauses that occur:
- R1: "In [names] the alignment was indistinguishable from zero ($\theta=[\theta]$, 95\%
  interval [lo] to [hi]; $z=[z]$), whereas dilution would have required
  $\theta\geq[\theta_{\min}]$; at most [f_X,max] of the signal can be dilution at the
  replicate-based reliability."
- R1w: "In [names] a weak alignment ($z=[z]$) was of negligible size ($\theta=[\theta]$),
  did not replicate, and is not interpreted."
- R3: "In [names] part of the signal was aligned with the measured transcript
  ($\theta=[\theta]$, 95\% interval [lo] to [hi], below the dilution minimum
  $[\theta_{\min}]$): at most [f_X,hi] of it can be dilution at any reliability at or above
  the replicate-based bound."
- R2: "In [names] the alignment ($\theta=[\theta]$, 95\% interval [lo] to [hi]; $z=[z]$)
  reached what dilution at the replicate-based reliability would produce
  ($[\theta_{\min}]$), so recovery of true expression lost to RNA-seq noise is not excluded
  there[; it was carried by tumour purity ($\theta=[\theta^{pur}]$ after adjustment)]."
- R4: "In [names] the prediction was anti-aligned with the measured transcript
  ($\theta=[\theta]$; $z=[z]$[, within-operator only]); [shared:] it disappeared under the
  transcriptome-wide baseline ($z=[z^{rnapc}]$), as expected if morphology reads an RNA
  measurement artefact shared across genes [not shared:] it persisted under the
  transcriptome-wide baseline ($\theta=[\theta^{rnapc}]$), consistent with a gene-specific
  RNA artefact or a transcript component that does not reach protein [unresolved:] the
  transcriptome-wide baseline did not resolve whether it is shared ($z=[z^{rnapc}]$)."
- R4w: "In [names] a weak anti-alignment ($\theta=[\theta]$; $z=[z]$) [confined:] was confined
  to proteins with missing values (fully quantified proteins: $z=[z^{compl}]$) [absorbed:] was
  absorbed by the transcriptome-wide baseline ($z=[z^{rnapc}]$) [unresolved:] was not resolved
  by [MNAR-exposed: the fully quantified proteins ($z=[z^{compl}]$) or] the transcriptome-wide
  baseline ($z=[z^{rnapc}]$)[ complete-gene set: and persisted in fully quantified proteins];
  at this size it is not the signature of an RNA artefact ($\theta\approx-1$)."
- R6 / R7: "In [names] the test was not decisive ([R6: no morphology-only signal on the
  set, $z=[z_D]$] [R7: $\theta=[\theta]$, 95\% interval [lo] to [hi]])."
- Operator rule (rule 7), operator-mediated: "Under an unrestricted null [names] read [as
  aligned / as anti-aligned] ($z=[z^{unres}]$) but not within operators; that part is carried
  by differences between acquisition operators ([phi_op] of the signal)."
- Operator rule, same direction: "In [names] an unrestricted null gave the same direction
  ($z=[z^{unres}]$); [phi_op] of the signal lies between acquisition operators."
- Main-text stratum (rule 10): "Within [cohort], [stratum] ([n] proteins) read [aligned
  ($\theta=[\theta]$, $z=[z]$) / anti-aligned ($\theta=[\theta]$, $z=[z]$)], one of [k]
  strata examined without correction."

Fixed closing, every outcome (it carries the qualifier Q of section 8.5 into the main text):
"These are the batch-uncorrected significant sets. The test concerns the linear
contribution of the gene's own measured transcript: it cannot separate post-transcriptional
regulation from protein-side technical variation that enters the residual in the same way,
such as interference shared within a TMT plex, and it does not exclude non-linear transcript
effects[ R1 all spline-stable: , although replacing the linear mRNA baseline by a spline
left $\theta$ within $\pm[band]$][ any spline-sensitive: ; with a spline baseline [names]
gave $\theta=[\theta^{spl}]$ ($z=[z^{spl}]$)]."

### 8.2 Methods: one paragraph appended to `\S\ref{sec:txomemeth}`

"\emph{Own-mRNA measurement error (post hoc).} This analysis was specified on [date] in a
dated file released publicly before it ran (Supplementary Note~\ref{snote:measerr}). For each
fully quantified protein, the residual of protein on own mRNA (both standardised) was predicted
from the 20 morphology components alone with the published estimator and folds. The alignment
$\theta=\sum_g(\hat\rho_g\,\mathrm{cov}(\hat r_{p,g},M_g)-E_0)/\sum_g(\mathrm{cov}(\hat
r_{p,g},r_{p,g})-E_0)$ over each significant set was centred on a within-operator
patient$\leftrightarrow$slide permutation null (within operator and scanner in PDAC; 200
fresh permutations shared by all genes; unrestricted null as sensitivity), with a
patient-bootstrap interval (200 resamples), and
$\sum_g\mathrm{sign}(\hat\rho_g)\,\mathrm{corr}(\hat r_{p,g},M_g)$ was tested against the same
null. Reliability was bounded below by the one-way intraclass correlation of
$\log_2(\mathrm{TPM}+1)$ between a case's RNA-seq files (distinct tumour samples; CCRCC and
UCEC, transferred to the other cohorts as the smaller of the two) and above by a Poisson
count bound; the deflation of $\theta$ by selection on significance was calibrated by
simulation. Proteins with any missing value were excluded because selection on protein
completeness biases $\theta$ (Supplementary Note~\ref{snote:measerr})[MNAR-exposed: , except
in [names], where fewer than 20 significant proteins were fully quantified and proteins with
up to 5\% missing values were included, a choice that biases $\theta$ downwards]."

[If lambda_T came from one cohort: "transferred to the other cohorts from [CCRCC/UCEC]"; if
from neither: "...UCEC) or, in [names], by the squared mRNA-protein correlation".]

### 8.3 Limitations (v): replace the last sentence

Replace "A denoised or replicate-averaged mRNA baseline would be the decisive test and has
not been run; until then, ``beyond the transcript'' means beyond its measured abundance." (one
sentence) by:

"A post hoc test of its signature (\S\ref{sec:txome}; Supplementary
Note~\ref{snote:measerr}) asks whether morphology's prediction of the protein residual is
aligned with the measured transcript, as dilution requires: [the class clauses below, joined
by semicolons]. The test addresses the linear contribution of the measured transcript only,
cannot distinguish post-transcriptional regulation from protein-side technical variation
shared within TMT plexes, and rests on replicate libraries in CCRCC and UCEC and on the
assumption that their reliability transfers to the other three cohorts. A denoised or
replicate-averaged mRNA baseline, which would also address non-linear and gene-specific
effects, has not been run; until then, ``beyond the transcript'' means beyond its measured
abundance."

Class clauses:
- R1: "in [names] the alignment was indistinguishable from zero, so linear dilution
  of the measured transcript accounts for at most [max f_X,max] of the signal at the
  replicate-based reliability"
- R1w: "in [names] a weak alignment of negligible size was not interpreted"
- R3: "in [names] part of the signal was aligned with the measured transcript, but at most
  [f_X,hi] of it can be dilution"
- R2: "in [names] the alignment reached what dilution would produce, so recovery of true
  expression lost to RNA-seq noise is not excluded there"
- R4: "in [names] the prediction was anti-aligned with the measured transcript, the signature
  of an RNA measurement artefact read by morphology or of a transcript component that does
  not reach protein[, shared across genes / not shared / unresolved]"
- R4w / R6 / R7: "in [names] the test was not decisive"
- Operator rule (added to that cohort's clause): operator-mediated: "with part of it carried
  by differences between acquisition operators"; same direction: "in the same direction under
  an unrestricted null"

### 8.4 Abstract and Discussion

**Abstract policy (the same in the three post hoc pre-specifications of file date
2026-10-08).** No favourable post hoc result enters the Abstract. The only permitted
Abstract edits are qualifiers, fixed in advance, that stop an existing sentence from
overstating in the light of a post hoc result:
C1-UCEC within-plex re-read,
P1 p_AUC >= 0.05: the uterus clause "a set-level replication;" becomes "a set-level
replication, not shown under a post hoc within-TMT-plex null;"; BRCA frozen-section arm,
any set-level reading in its supplement table with p < 0.05: "a whole-proteome cohort
none." becomes "a whole-proteome cohort none at the gene level."; measurement-error checks:
the R2 and R4 clauses of its section 8.4. Nothing else. If the Abstract is later shortened,
a qualifier added under this policy is kept in substance.

Here: the R2 and R4 clauses below are the only Abstract edits of this analysis. Each
reports an outcome unfavourable to the paper's reading and qualifies the existing clause on
single transcript measurement, which would otherwise understate how much of the signal may
be true expression (R2) or omit that the signal behaves as an RNA artefact would (R4) and
so overstate the closing sentence's "protein-regulatory tissue state". No other outcome
changes the Abstract.

- **No cohort's primary reading is R2 or R4 (every other outcome, including R1 in all five):
  the Abstract is unchanged**, and so is the Discussion sentence "because that transcript is
  measured once, part of what morphology recovers may be true expression lost to measurement
  noise rather than post-transcriptional regulation (Limitations~(v))". Reason, fixed now:
  no favourable post hoc result enters the Abstract, and the clause stays literally true,
  because the test does not exclude non-linear transcript effects.
- **At least one cohort reads R2, none R4:** the Abstract clause "with each transcript
  measured once, part may be true expression, not post-transcriptional regulation." becomes
  "with each transcript measured once, part may be true expression, not post-transcriptional
  regulation, and a post hoc test could not exclude this for a substantial part of the signal
  in [names]."; the Discussion sentence gains "; a post hoc test could not exclude this in
  [names]", inserted immediately before " (Limitations~(v))".
- **At least one cohort reads R4, none R2:** the Abstract clause becomes "with each
  transcript measured once, part may be true expression, not post-transcriptional regulation;
  in [names] the signal was anti-aligned with the measured transcript, as expected if
  morphology reads an RNA artefact (post hoc)."; the Discussion sentence gains "; in [names]
  the signal was anti-aligned with the measured transcript (post hoc)", inserted immediately
  before " (Limitations~(v))".
- **Both R2 and R4 occur:** R2's clause first, then R4's: Abstract "with each transcript
  measured once, part may be true expression, not post-transcriptional regulation, and a post
  hoc test could not exclude this for a substantial part of the signal in [names R2]; in
  [names R4] the signal was anti-aligned with the measured transcript, as expected if
  morphology reads an RNA artefact (post hoc)."; Discussion "; a post hoc test could not
  exclude this in [names R2]; in [names R4] the signal was anti-aligned with the measured
  transcript (post hoc)", inserted immediately before " (Limitations~(v))".
- If the C1-LUAD arbitration is active, LUAD's class for this section is the arbitrated one.
- Operator clauses, strata and sensitivities never change the Abstract or the Discussion.

### 8.5 Supplement: new Supplementary Note `snote:measerr` and table `tab:measerr`

The note is headed "Own-mRNA measurement error (specified [date], UTC)". It gives the model
(one paragraph), the statistic, nulls, eligibility, the kappa table and the simulation
calibration (pointing to the public DERIVATION.md and `theory/final-rules/`), this file's
prior exposure item 2 in one sentence, the departures of section 7 in one sentence each, and
per cohort, then per read stratum, the sentence of the cohort's rule below, verbatim with
all numbers. `[lambda-source]` = ", transferred from CCRCC and UCEC" / ", transferred from
[CCRCC/UCEC] only" / ", the squared mRNA-protein correlation (no replicate libraries)" / empty
(measured).

**Fixed qualifier Q** (appended verbatim to every sentence marked +Q):
"This reading concerns the linear contribution of the gene's own measured transcript; it
cannot separate post-transcriptional regulation from protein-side technical variation that
morphology can read, such as interference shared within a TMT plex, and it does not exclude a
non-linear transcript contribution."

- R1 (+Q): "In [cohort], the morphology-predicted protein residual of the [n] [fully quantified]
  significant proteins was not aligned with the measured transcript ($\theta=[\theta]$, 95\%
  interval [lo] to [hi]; within-operator permutation $z=[z_N]$), whereas regression dilution
  at the reliability bound $\lambda_{lo}=[\lambda_{lo}]$[lambda-source] produced
  $\theta\geq[\theta_{\min}]$ after selection in simulation; at most [f_X,max] of the signal
  can be linear dilution at that reliability. With a spline mRNA baseline the increment was
  [dp_spline] (published [dp_linear]) and $\theta$ was [theta_spline] ($z=[z_spline]$)."
- R1w (+Q): "In [cohort], a weak alignment ($z=[z_N]$) of negligible size ($\theta=[\theta]$,
  95\% interval [lo] to [hi]) did not replicate in another stratum, family set or cohort and
  is not interpreted."
- R3 (+Q): "In [cohort], the alignment ($\theta=[\theta]$, 95\% interval [lo] to [hi]; $z=[z_N]$)
  lay below the dilution minimum $[\theta_{\min}]$ at $\lambda_{lo}=[\lambda_{lo}]$
  [lambda-source]: at most [f_X,hi] of the morphology--residual covariance can be linear
  dilution of the measured transcript at any reliability at or above that bound, and the
  remaining [post_share_lb] is not explained by it.[ It was carried by tumour purity
  ($\theta=[\theta^{pur}]$, $z=[z^{pur}]$ after adjustment).][ Purity was not available.]"
- R2: "In [cohort], the alignment ($\theta=[\theta]$, 95\% interval [lo] to [hi]; $z=[z_N]$;
  a reliability of $\hat\lambda=[\hat\lambda]$ if the signal were entirely dilution) reached
  the range that dilution produces at $\lambda_{lo}=[\lambda_{lo}]$[lambda-source]
  ($\theta_{\min}=[\theta_{\min}]$). We therefore cannot exclude that morphology recovers true
  transcript abundance lost to RNA-seq noise in [cohort], and there the estimand is protein
  beyond the measured transcript.[ It was carried by tumour purity ($\theta=[\theta^{pur}]$,
  $z=[z^{pur}]$ after adjustment).][ Purity was not available.]"
- R4: "In [cohort], the morphology-predicted residual was anti-aligned with the measured
  transcript ($\theta=[\theta]$, 95\% interval [lo] to [hi]; within-operator $z=[z_N]$,
  unrestricted $z=[z^{unres}]$)[, within operators only]. [shared:] Under the
  transcriptome-wide (20-RNA-PC) baseline the anti-alignment disappeared ($z=[z^{rnapc}]$): a
  shared RNA measurement artefact read by morphology is the leading explanation. [not shared:]
  It persisted under the transcriptome-wide baseline ($\theta=[\theta^{rnapc}]$,
  $z=[z^{rnapc}]$): a gene-specific RNA artefact, or a transcript component that morphology
  reads but that does not reach protein; $\theta$ does not separate the two. [unresolved:] The
  transcriptome-wide baseline did not resolve whether it is shared ($\theta=[\theta^{rnapc}]$,
  $z=[z^{rnapc}]$); an RNA artefact and a transcript component that does not reach protein
  both remain possible."
- R4w (+Q): "In [cohort], a weak anti-alignment ($\theta=[\theta]$; $z=[z_N]$) [confined:] was
  confined to proteins with missing values (the [n_c] fully quantified proteins: $z=[z^{compl}]$),
  consistent with selection on protein completeness [absorbed:] was absorbed by the
  transcriptome-wide baseline ($z=[z^{rnapc}]$) [unresolved:] was resolved neither by
  [MNAR-exposed: the [n_c] fully quantified proteins ($z=[z^{compl}]$) nor by] the
  transcriptome-wide baseline ($z=[z^{rnapc}]$)[ complete-gene set: ; it persisted in fully
  quantified proteins]; with $\theta>-0.3$ it is not the signature of an RNA artefact
  ($\theta\approx-1$), and no reading is taken."
- R5: "Within [cohort], [strata] fell in different readings or differed by more than two
  bootstrap standard errors of $\theta$ ([values]; [m] comparisons, not corrected) and are
  reported separately."
- R6 / R7: "In [cohort], no reading: [R6: no morphology-only signal on the set
  ($z_D=[z_D]$)] [R7: $\theta=[\theta]$, 95\% interval [lo] to [hi], $z=[z_N]$]."
- Operator rule, same direction: "Under the unrestricted null the set read [reading]
  ($z=[z^{unres}]$, $\theta=[\theta^{unres}]$), in the same direction as within operators;
  [phi_op_D] of the permutation-centred signal ([phi_op_dp] of the nested increment) lies
  between operators."
- Operator rule, operator-mediated: "Under the unrestricted null the set read [reading]
  ($z=[z^{unres}]$, $\theta=[\theta^{unres}]$), within operators [primary reading]
  ($z=[z_N]$): the difference is between-operator structure ($\phi_{op}=[phi_op_D]$ of the
  signal, [phi_op_dp] of the nested increment), read as operator-mediated [transcript
  alignment / RNA artefact], not as gene-level dilution or artefact."
- Main-text stratum (rule 10) in the supplement: the same sentence as in section 8.1, with k.

The table has one row per cohort x set or stratum with: n genes, s, lambda_lo (source),
lambda_hi, kappa, theta_noise,min, theta_pc and interval, z_D, z_N (primary, unrestricted,
plex), z_S, phi_op (D and Delta_p), purity / spline / 20-RNA-PC / raw-theta values,
missingness strata and the reading. Depth tertiles are described with the fixed sentence
"Dilution predicts that $\theta$ rises with depth and that few selected proteins sit in the
top tertile; $\theta$ by tertile was [values] and the top tertile held [n] of [N] selected
proteins."; the missingness pattern with "[A monotone fall of $z$ with missing protein
([values]), the signature of selection on completeness, confined to proteins the primary
analysis excludes / No monotone trend with missingness ([values])]." A cohort not analysed is
one row: "not run: [F1 failed / not completed (section 10)]"; GBM's purity cells read "not
available".

## 9. What this analysis cannot establish

- That the residual is biological: protein-side technical artefacts read by morphology look
  like u (section 0).
- That no transcript-borne component exists: non-linear mRNA-to-protein maps, and mixtures
  with f_X below about 0.1, can read R1.
- The reliability of LUAD, GBM and PDAC transcripts: transferred, not measured.
- A gene-level attribution: per-gene theta is heavy-tailed and no gene-level call is made.
- That an anti-alignment is an artefact: a transcript component that morphology reads but
  that does not reach protein gives theta = -1 as well.
- Anything about the batch-corrected discovery sets of Table 1 or, unless arbitration is
  active, about the confirmatory cohorts.

## 10. Freeze list and run order

**Line endings.** Hashes are of the bytes as they are. `theta_kernel.py` has CRLF line
terminators and its hash equals the server copy recorded in all five manifests; it must not
be normalised. `theory/DERIVATION.md`, `README_commands.txt` and several files of
`theory/final-rules/` are CRLF as well and are hashed as they are. `local_theta.py`,
`test_local_theta.py` and this file are LF.

**Files the gate requires** (`local_theta.py` refuses real data unless each is listed and
unchanged):

| File | sha256 at writing |
|---|---|
| `review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md` | (this file; in the sidecar) |
| `M/local_theta.py` | `2d10449a1d30b95cbebb6a73e13aee3ea88d77c6b8d739523f9c4b71dcae7632` |
| `M/theta_kernel.py` | `6c95cd119d28e811e09967946e9f191c69f07f54f653991917e817ef99ec2734` |
| `M/export/ccrcc/inputs_ccrcc.npz` | `68933ea6fd5ba174a360d7ced323d02394d7e105ed8dfe9b8fa016150cd67131` |
| `M/export/luad/inputs_luad.npz` | `b64cdf3b31717a32e0d42c42f3f1dde74786d5598008f0d214297a71adaba095` |
| `M/export/ucec/inputs_ucec.npz` | `e28b86fadd41ec86363577fffad1db1f2f8a5162e732c7acb0e4133c0982044e` |
| `M/export/gbm/inputs_gbm.npz` | `c35cb9fa00b3010dec7c8a7b5f154126124a5698a6884700cf5c5d7d9f5b55fd` |
| `M/export/pdac/inputs_pdac.npz` | `cf8497b2c097bebe655668384d9320c8a0f3bf347d69ec88ce1848da354cf64d` |
| `M/side/families_reactome2022.tsv` | built in F2 |
| `M/side/purity_ccrcc.tsv`, `purity_luad.tsv`, `purity_ucec.tsv`, `purity_pdac.tsv` | built in F2 |

No `theta_kernel_local.py` exists or is needed: the sign-weighted statistic is `Nc`/|rho-hat|
on the unchanged kernel's output.

**The sidecar list `M/SIDECAR_FILES.txt`** (one path per line, relative to
`MorphoResidual_paper/`; `__pycache__/` never listed) holds, besides the files above: itself;
every other exported file of every used export as named in `export_sha256_<c>.txt`
(`gene_table`, `patient_table`, `repro_check`, `run_manifest`, `export_sha256`, and for CCRCC
and UCEC `replicate_icc` and `replicate_files`; all verified against the transfer hashes,
Status); `M/export/meas_export_5cohorts.tgz` (`2b22fedf...`); `M/measurement_error_checks.py`
(`66b2a7fb...`), `M/lsf_measurement_error.sh` (`7bc94395...`), `M/test_local_theta.py`
(`694d6f74...`), `M/test_synthetic.py` (`32207f57...`), `M/README_commands.txt` (after F3);
the superseded drafts kept for the record (`M/measurement_error_checks_v1.py`,
`M/lsf_measurement_error_v1.sh`, `M/README_commands_v1.txt`, `M/test_synthetic_v1.py`,
`M/local_theta_pre-review.py`, `M/test_local_theta_pre-review.py`,
`M/theory/DERIVATION_pre-review.md`); `M/theory/DERIVATION.md` (`3881fba2...`) and every
simulation file of `M/theory/` (`sim_*.py`, `sim_*.csv`, `sim_*.log`, `section7_tables.md`,
`verify/*`); the calibration evidence copied from the drafting session's scratch folders into
`M/theory/final-rules/` (from `phaseC/final-rules/`: `vt_final.py`, `reclass.py`,
`reclass.log`, `vt_F1.log`, `vt_F2_F4.log`, `vt_F3a.log`, `vt_F3b.log`, `vt_F5.log`,
`vt_F5b.log`; from `phaseC/verify-theory/`: `vt_finite.py`, `vt_pop.py`, `vt_pop.log`,
`vt_E2.log` to `vt_E7.log`, `vt_E6.py`, `vt_E7.py`, `vt_bias.py`, `vt_bias.log`; byte-identical
copies; `vt_final.py` and `vt_E*.py` import `vt_finite` from their own folder; the editing
script `apply_edits.py` of that folder is not evidence and is not copied);
`server_export/scripts/residual_analysis.py` (`a54a0a56...`); the published per-gene results
the exports were checked against, `server_export/pinned/<c>/residual_results_tumoronly.csv`
(their sha256 equal the `sig_csv_sha256` of each run manifest); `composition_pathologist.py`
(the source of the purity definitions); the F2 builder `M/side/build_side_files.py` and its
sources (`figures/figdata/Reactome_2022.gmt` `c3ea9df7...`, `figures/enrichment_c2_levels.py`
`9116a732...`, `figures/figdata/subtype_ccrcc_raw.xlsx` `73137225...`,
`figures/figdata/subtype_ucec.txt` `7c3d6c36...`,
`figures/figdata/composition_luad_cbioportal.csv` `489ff1b0...`,
`figures/figdata/subtype_pdac_raw.xlsx` `252049de...`,
`server_export/scripts/aliquot_to_case_tumor.tsv` `e823ce38...`); and, if F4 is done, the
C1-LUAD export files. The F2 outputs are already in the list; one that is not built is removed
from the list before the sidecar is written, and the analysis that needs it is "not run (not
frozen)" and cannot be added later. Nothing else is released for this analysis.

**What is released (the same rule in the three post hoc pre-specifications of file date
2026-10-08).** This file, its sidecar and every file the sidecar lists are released on
GitHub before the run, except derived data files (`.npz` and `.tgz` files and the per-gene
result tables named below) and third-party source tables under `figures/figdata/`. These
follow the repository's data policy: each is fixed by its sha256 in the sidecar, the
derived data files are deposited with the other derived data at acceptance, and the
third-party tables are public at their sources. The release commit message lists every
withheld path.

Here, withheld: the five `M/export/<c>/inputs_<c>.npz`, `M/export/ccrcc/replicate_files_ccrcc.npz`,
`M/export/ucec/replicate_files_ucec.npz`, `M/export/meas_export_5cohorts.tgz`, the five
`M/export/<c>/gene_table_<c>.csv`, `M/export/ccrcc/replicate_icc_ccrcc.csv`,
`M/export/ucec/replicate_icc_ucec.csv`, the five
`server_export/pinned/<c>/residual_results_tumoronly.csv`, and the five source tables under
`figures/figdata/` named above. Everything else in the list is released.

**Before freezing (none computes an alignment statistic or reads a morphology-omics relation):**
- **F2.** Write `M/side/build_side_files.py` (reads only the GMT, the keyword map, the
  covariate tables of section 6 item 2 and the `patients` array of each inputs file, for
  coverage) and run it to produce the families file and the four purity files; record per
  file the number of genes or patients covered.
- **F3.** Correct `README_commands.txt` steps 6-7 (stale: they still give `--strata
  none,operator,plex`, say the RNA-PC variant is not implemented and name the old seed). Where
  it differs from this file, this file governs.
- **F4 (optional; enables rule 8).** A confirmatory export reusing
  `c1_run_test_luad.prepare()` (K = 22 PCs, registered fit sets and folds) that writes the
  `inputs_<c>.npz` layout with operator labels, raw counts, and the discovery 2,310 / 726
  memberships as a sets file, and reproduces the published C1-LUAD per-gene confirmatory
  increments to 1e-10. It is not written, and the server inputs
  (`/public/home/fjhui/ZW/luad_c1/`) have not been checked from this session.
- **F5.** Add any further exposure to section 0, item 4; complete `M/SIDECAR_FILES.txt`; write
  the sidecar with the command below; release this file, the listed files and the sidecar on
  GitHub under the release rule above; record the release commit and time.

```
cd MorphoResidual_paper
{ echo "# sha256 of the frozen files for POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md (paths relative to MorphoResidual_paper/)";
  echo "# written $(date -u +%Y-%m-%dT%H:%M:%SZ)";
  grep -v '^#' review/MEAS_ERROR_2026-10-08/SIDECAR_FILES.txt | grep -v '^[[:space:]]*$' | xargs sha256sum -b; } \
  > review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md.sha256
```

The command must end without a "No such file" message; every line has the form
`<sha256> *<path>`, which is what `local_theta.py` reads.

**Run order after freezing and release (one completed run per cohort; no smoke run on real
data; all runs local, on the machine that wrote the sidecar, from `M/`).**
1. **F1, published quantities only.** Per cohort: `python local_theta.py --inputs
   export/<c>/inputs_<c>.npz --repro-only --repro-pvals 20 --check-pca --check-icc --out
   export/<c>/local_theta` (it computes no alignment statistic). Required: exit 0,
   "REPRODUCES the published per-gene results (<=1e-10)", p-value mismatches 0, ICC difference
   <= 1e-10 in CCRCC and UCEC. The `--check-pca` difference is reported only (the analysis uses
   the exported PCs, not a refit, and PCA numerics are host-dependent, Status). A cohort that
   fails F1 is not analysed and is reported as "not run: F1 failed" with its output; the other
   cohorts proceed, and section 5's availability rule governs the transfer. The F1 outputs are
   hashed into the results record.
2. **Common arguments:** `--B 200 --seed 20261007 --boot 200 --boot-seed 20261008 --sets-file
   side/families_reactome2022.tsv --save-null --workers 8` (the worker count does not change
   any number).
3. CCRCC: `--strata operator,none,plex --primary-arm operator --variants spline,rnapc,purity
   --purity-file side/purity_ccrcc.tsv`.
4. UCEC: as CCRCC with `side/purity_ucec.tsv`.
5. Write lambda_T and the three tertile values (section 5, including its availability rule)
   into `M/export/lambda_transfer.txt` with four decimals and their sources.
6. LUAD: as CCRCC with `side/purity_luad.tsv --lambda-lo-transfer <lambda_T>
   --lambda-lo-transfer-tertiles <a,b,c>` (both flags omitted if neither CCRCC nor UCEC
   contributes).
7. GBM: as LUAD but `--variants spline,rnapc` and no purity file.
8. PDAC: as LUAD with `--strata operator_scanner,none,plex --primary-arm operator_scanner`
   and `side/purity_pdac.tsv`.
9. If F4 was frozen: C1-LUAD with `--confirmatory --strata operator,none --primary-arm operator
   --lambda-lo-transfer <lambda_T>` and the discovery-membership sets file.
10. Hash every output (`theta_<c>_{sets,reading,hetero,variants,pergene,manifest}`, `_boot`,
   `_null_<arm>`, the F1 outputs, `lambda_transfer.txt`) and the console log of every
   invocation into a results record; apply section 7 by hand to the reading files; fill the
   sentences of section 8 exactly as written.

**Failed or interrupted runs (the same wording in the three post hoc pre-specifications of
file date 2026-10-08).** A run that stops before writing its result file is reported with its
log and repeated once, unchanged, on the same host with a rerun tag; if the pinned host is
unavailable, a new smoke/gate is run on another host first (recorded) and the run is pinned
there; if the rerun also fails, the analysis is reported as not completed and no further run
is made. A local analysis invocation that exits FATAL before writing its result is not a run,
is recorded, and may be repeated after the input is corrected.

Here: the result file of a run is `theta_<c>_reading.csv`; the host is the machine that wrote
the sidecar; the rerun tag is `--out export/<c>/local_theta_rerun1`; the gate on another host
is F1 for that cohort plus `test_local_theta.py` ("ALL CHECKS PASSED") on that host, recorded
with its hostname; "exits FATAL" includes any `REFUSED:` exit of `local_theta.py` before
`theta_<c>_reading.csv` is written (sidecar mismatch, non-reproduction, a wrong path or
argument), and "the input is corrected" means restoring the frozen bytes or correcting the
command line, never changing a frozen file. A cohort reported as not completed is one row of
the supplement table ("not run: not completed") and contributes to no transfer. A run that
completes is never repeated with different settings.
