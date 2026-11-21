# Post hoc re-reads of the two confirmatory tests with transcriptome-derived covariates (C1-UCEC and C1-LUAD): pre-specification (file date 2026-11-21, local time UTC+8)

**Dates and times.** The date in this file's name is the local date (UTC+8) on which
drafting began. Every time stated inside this file is UTC. The sidecar's `# written` line
records when the file was frozen; the manuscript cites that UTC date.

**Status.** Drafted from 2026-11-20 UTC, after the published discovery results, the
published C1-UCEC reading (2026-09-16), the registered C1-LUAD reading (2026-10-02), the
post hoc analyses pre-specified on 2026-10-02 and 2026-10-08, and the Abstract rewrite of
`C1_LUAD_ADDENDA.md` Entry 12 (2026-11-20). An internal review of 2026-11-20 asked whether
the confirmatory set-level ordering needs histology, or whether any reading of tissue state
taken from the same RNA-seq would give it; the corresponding author asked for that to be
tested (2026-11-20). A first draft was reviewed adversarially on 2026-11-20 by six
independent readers, a second on the same day by four, and a third by two before freezing;
this text reflects the three reviews.
The main changes they caused: one freeze instead of two; a capacity comparison over 119
random reassignments instead of one; sets of fixed size instead of FDR-defined sets;
wording fixed by the source of an unfavourable reading; the complete Abstract for every
outcome; and no reading closed after freezing. The first draft and the three sets of
findings are kept beside the scripts (`DRAFT1_*`, `AUDIT1_findings.txt`,
`AUDIT2_findings.txt`, `AUDIT3_findings.txt`) and are released with them.

This file is frozen **once**, by the sha256 of this file, every script (local and cluster),
every input and the gate logs of section 9, listed in the sidecar
`POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256`. It is released publicly on GitHub, as
an annotated tag and a Release, after the corresponding author's written go-ahead and
before any stage of section 9 (H production, L1, L2) is run. Before freezing, only the
checks listed in section 9 are run: they use synthetic data, published quantities, or
RNA-seq alone.

## 0. Question, scope and prior exposure

**Question.** Both confirmatory tests rank genes by how much histology adds to the gene's
own mRNA in new patients and ask whether the genes selected in discovery rank above the
rest. Their null breaks the pairing between a patient's slide and that patient's omics. It
cannot tell an ordering that needs histology from one that any measure of tissue state
would give, including one that the same RNA-seq already holds. This analysis asks three
things of each confirmatory cohort:

1. **Conditional (primary).** Does the published uncorrected selection set still outrank
   background when the confirmatory baseline already holds transcriptome-derived covariates?
2. **Reverse.** Does that set also outrank background when those covariates take the place
   of histology in the confirmatory cohort, with no slide used?
3. **Forward and matched.** Do sets chosen in discovery without histology outrank
   background under confirmatory histology, and do sets chosen by histology after the
   covariates entered the discovery baseline outrank background when the confirmatory
   baseline holds covariates built in the same way?

None of the three is a test of specificity. The first asks whether histology adds to the
covariates; the others describe how far the covariates give the same ordering.

**Scope.**

- **Post hoc, reported whatever the result.** Nothing below was pre-registered. After
  freezing, no covariate, set, arm, cell, background, threshold, reading or sentence is
  added, dropped or changed. The only later documents are: a correction under the Defects rule
  of section 9; the Stage H record and the results record of section 9; and the dated
  entries of `C1_LUAD_ADDENDA.md` that sections 7.4 and 7.5 require.
- **No effect on the registered or published readings.** The published C1-UCEC reading
  (AUC 0.597; null 0.498 +/- 0.036; p = 0.004) and the registered C1-LUAD reading (AUC
  0.644; null 0.524 +/- 0.019; p1 = 0.001; p_strat = 0.001; p2 = 0.001) stand as reported.
  Nothing here is pooled with either, and no wording of
  `C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md` section 5.2 changes.
- **No new cohort.** Section 5.8 of the C1-LUAD rule ("No further CPTAC-3 confirmatory
  cohort is attempted before submission") is not touched: these are re-reads of the two
  cohorts already read.
- **Only the uncorrected sets are read.** The batch-corrected sets (201 genes in UCEC, 726
  in LUAD) are reported on the new arms and not read (section 5): the registered reading of
  the 726-gene set rests on operator-residualised components and a within-operator null,
  which no new arm reproduces. The Abstract's clause on the 726-gene set therefore keeps no
  qualifier whatever the result. Every sentence of section 7 that reports a cell's result
  names the uncorrected set.
- **One null.** Every new confirmatory arm is read against the unrestricted
  patient<->slide permutation null of the published tests; no operator-stratified or
  plex-stratified null is run, to limit arms. A restricted null usually has a higher mean
  but can have a smaller spread (C1-UCEC published set: z = 2.73 unrestricted, 2.89 within
  operator), so nothing can be inferred about a stratified null in either direction.

**Departures from earlier rules and plans (each stated with its reason).**

1. C1-UCEC rule section 3 and C1-LUAD rule section 3 fixed the confirmatory estimator; the
   primary and reverse arms change its baseline or its block. The UCEC rule's section 5 also
   says "No other framework ... is to be substituted after this file exists". These re-reads
   are added beside the fixed tests and replace nothing.
2. Section 2 of both rules fixed the selection sets as not re-derived; the forward and
   matched cells use other discovery sets, as `POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md`
   did.
3. C1-LUAD rule 5.7: "The Abstract's LUAD wording is the section 5.2 clauses", item 0
   ("Room the realised branch does not use is not refilled after the LUAD result"; "No
   result clause is cut") and item 3 ("No other descriptive item enters the Abstract"). In
   an unfavourable outcome section 7.4 adds a sentence after those clauses and makes room
   by removing words fixed now. No protected content item (organ, n, set, AUC, interval, p,
   counts, retention, caveat, external nulls) is removed; words inside result clauses are
   (cuts c1, c3, c5, c6 of section 7.4). Cut c1 removes the cohort-size hedge from Entry
   12's clause, which makes an unfavourable statement stronger. The only recorded precedent
   is Entry 12, which added an unfavourable clause under the same reading.
4. The Abstract policy of the three post hoc pre-specifications of 2026-10-08 lists its
   permitted qualifiers and says "Nothing else." Entry 12 reads that policy as limiting the
   Abstract edits those three analyses may cause. This file follows that reading and adds a
   qualifier of the same kind for a new analysis.
5. C1-LUAD rule 5.7 item 6 ("Every claim names the set") is kept for every sentence that
   reports a cell's result. The set is named as "the uncorrected [uterine|lung] set" or
   "uncorrected discovery set", not always in the rule's own phrase ("genes selected by the
   uncorrected discovery analysis (2,310)"), as in the Abstracts of Entries 11 and 12. Item 2 is kept: "replicat\*" is not used for any result of this
   analysis.
6. `NC_ROADMAP.md` (2026-09-07; T1-1 and section 6, R4) rejected a whole-gene ridge
   permutation FDR as computationally infeasible and as a risk ("morphology is only a noisy
   proxy of the transcriptome"), and rejected a "symmetric" comparison of histology with RNA
   components on the published sets because a selection effect favours the block that made
   the selection; R11 rejected covariate scans as an empty table and a new multiplicity
   family. Sections 3 and 5 run related analyses. What changed: reviewers of the finished
   manuscript raised the question, so the risk is taken on purpose; a batched kernel makes
   the whole-gene permutation computable; the comparison is made on confirmatory patients,
   not on the selection data; the forward, matched and reverse cells are barred from any
   claim that one selector beats another (section 5); and the new readings are left
   uncorrected because every headline consequence follows from a failure (section 6.1).
7. The plan that prompted this analysis (`1_plan_from_five_reviews.md`, B1) had a third
   part: freezing the discovery PCA basis and ridge weights and applying them to the
   confirmatory patients. It is not done here. The discovery weights were not stored, the
   statistic would have no cross-validation and a different number of components (20
   against 28 and 22), and it bears on per-gene transfer, which the paper already reports as
   not holding, not on the question above.
8. The checked plan (`2_plan_checked_against_rules_FINAL.md`, D1) asked for the choice
   between two framings of the paper to be written into this file, and for a favourable
   result to enter Results only. Section 7.5 fixes the title and the opening sentences of
   the Abstract, the Introduction and the Discussion for every outcome. In the favourable
   case this file also adds one sentence to Limitations (iv) and rows to Fig. 6a (sections
   7.8, 7.9), and in every case an item in the Introduction's list of reported tests
   (section 7.6), because leaving the result out of the figure that shows every other re-read
   would be selective. Left to the authors after the results, and nothing else: the target
   journal; which free text is shortened if the main text exceeds the journal's length; the
   prose of the Supplementary Note beyond its fixed sentences; the layout and column order of
   the two supplementary tables, and the drawing and legend of the new figure rows; the
   wording of the dated entries of `C1_LUAD_ADDENDA.md` and of the provenance rows; a cut to
   the Abstract that a journal's limit forces, within the bounds of section 7.4; a further
   weakening of the framing under section 7.5; and the Methods wording that discloses a
   defect found after a run (section 9).
9. Entry 12 made the Abstract's opening sentence identical to a sentence of the
   Introduction, and Entries 11 and 12 record the journal's instruction that the Abstract
   contain "Here, we show". Cut c2 (Case 1) and the Case 2 openers of section 7.5 end both.
   The cut is made for length and the Case 2 change follows from the result.

**Prior exposure and expectations (disclosed).**

1. *What is expected.* The internal review that prompted this analysis put the chance of a
   favourable outcome slightly below one half and tied an unfavourable one to a change of
   target journal. The drafting session expects the lung set not to outrank background once
   the composition scores are in the baseline (item 2). For the uterine cohort it also
   regards class N as more likely than not: the published z is 2.73 and 36% of the median
   increment was retained in discovery, so Case 2, and with it the change of title, is the
   outcome the drafting session regards as most likely. Several choices were made under these expectations and could favour
   the paper: an unresolved reading under the expression components (class U) changes no headline
   wording (it is stated in Results, Discussion and Limitations); the capacity comparison is
   expected to tilt toward that class and has the resolution of 119 draws
   (section 5); only the unrestricted null is run; the batch-corrected sets are reported and
   not read, so the Abstract's clause on the 726-gene lung set stands unqualified in every
   outcome, and "the uncorrected lung set did not" can be read as if the batch-corrected set
   had been tested and had held; a reverse or forward cell that outranks changes no headline
   wording (it is stated in Results and, for reverse cells, in Limitations (iv)); the sets
   read in the forward family have
   a fixed size and are diluted when the selecting analysis finds few genes, which makes
   "outranked" less likely in a family where it is unfavourable (section 3); the title and
   the opening sentences change only if both cohorts are flagged, although class N requires
   a fall to about 0.26 of the published margin in the lung cohort and only to about 0.60 in
   the uterine cohort (section 6.3); a reading closed before freezing is worded "was not
   run", which is milder than "did not"; and the target journal is left open.
2. *Discovery, composition.* The published composition control (`sec:stroma`) added the
   three scores of section 2 to the discovery baseline: the median increment of the
   published sets retained 36% in UCEC and 6% in LUAD (68% CCRCC, 34% PDAC, 3% GBM). The
   per-gene values are in `server_export/results/results__composition_controls_<cohort>.csv`
   for every gene of the published sets: the gain of each score and of all three over own
   mRNA, and histology's increment over own mRNA plus the scores. In the other direction,
   published deconvolution estimates (ESTIMATE, xCell, Gillette et al.) left the
   morphology-residual link at 93-104% in the other four cohorts (`sec:stroma`).
3. *Discovery, expression components.* The published transcriptome-wide baseline
   (`sec:txome`) added the 20 components of section 2 to the discovery baseline: a median
   40% of the noise-adjusted advantage was retained (77% in UCEC, the highest; 40% in LUAD),
   and UCEC kept a positive median increment over that baseline (+0.075; Entry 12).
   Per-gene values and permutation p-values under that baseline exist for the genes of the
   published sets (`server_export/results/scripts__transcriptome_baseline_<cohort>_perm.csv`):
   227 (UCEC) and 55 (LUAD) genes pass FDR < 0.05 with a positive increment when
   Benjamini-Hochberg is taken over the published set, and within the published sets 99 and
   55 genes sit at the permutation floor. Taken over all fit genes, as section 3 does,
   Benjamini-Hochberg needs 211 and 215 genes at the floor, so the FDR-defined sets S_M|R20
   are expected (not computed) to be empty in both cohorts and S_M|C3 possibly in LUAD. That
   is why the sets that are read are of fixed size (section 3). The canonical correlations
   between the discovery morphology and expression components are published (`snote:txome`).
4. *Discovery, measurement-error analysis.* `POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md` ran,
   on the discovery data of both organs, sensitivity variants that adjusted the alignment
   statistic for tumour purity and for the 20 expression components (in-sample covariates,
   not a permuted block).
5. *Discovery, operator.* `incremental_r2_batchonly` in
   `server_export/results/<cohort>__sitepack_operator.csv` is the per-gene increment of
   operator dummies in place of histology, without a p-value (Supplementary Table posthoc,
   "B.o. vs n.m.").
6. *Confirmatory, other sets.* Morphology-selected sets other than the published ones have
   been read on the confirmatory matrices: the operator-corrected sets S1 and S2, the other
   sets of `POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md`, discovery top-k sets and the
   selected sets of the other four organs (UCEC; 2026-10-07 and 2026-10-08), and the
   operator-corrected 726 set (LUAD; registered step 3). These are sets chosen by histology
   after operator correction. The set T_OP of section 3 is different: it is chosen by
   operator indicators in place of histology. The per-gene confirmatory increments have been
   in the published gene-level tables since 2026-09-16 (UCEC) and 2026-10-02 (LUAD).
7. *Blind.* To the knowledge of the corresponding author and of the drafting session, and
   after a search of the local project (scripts, logs, result folders, session notes; the
   cluster's folders were not searched and are listed in the blind-run record of section 9):
   no statistic that uses the composition scores, the expression components, or a gene set
   chosen by the scores, by the components or by operator indicators in place of histology
   has been computed on either confirmatory cohort, other than the RNA-only quantities
   of the blind run of section 9 (markers found and their standard deviations, and the
   variance explained by the 20 components); and no discovery analysis has used the
   scores, the components or operator indicators as the permuted block.
8. *The drafting session.* It read array shapes, gene lists and label counts of the files of
   section 1, read the existing result files named above, and ran the checks of section 9.
   The validation run printed, for the published analysis M, the number of genes at the
   permutation floor (1,748; 1,354) and of rejections before the sign filter (3,263; 3,086);
   the fixed-size rule of section 3 was written with these known. A reviewer of the second
   draft ran a simulation on synthetic data only. Nothing else was computed on project data.

## 1. Inputs

Paths are relative to `MorphoResidual_paper/`. Hashes below are prefixes; the full sha256
of every file is in the sidecar, and the scripts assert the full values.

- **Discovery inputs.** `review/MEAS_ERROR_2026-10-08/export/ucec/inputs_ucec.npz`
  (`e28b86fa...`; 100 patients, 10,512 fit genes) and `.../luad/inputs_luad.npz`
  (`b64cdf3b...`; 105 patients, 10,724 fit genes): own-gene mRNA (log2(TPM+1)), protein,
  the 20 published morphology components, the 20 expression components (`rna_pcs20`),
  operator labels and the published per-gene results.
- **Discovery composition scores** (exported on the cluster before freezing, step H0 of
  section 9): `server_export/results/c1_covariates/stage_h/c1_cov_discovery_scores_ucec.tsv`
  and `..._luad.tsv`.
- **Confirmatory matrices.** `server_export/scripts/pinned/c1_ucec_d4_stratified.npz`
  (`ef75f1c0...`): 10,146 usable genes x 1,001 columns (column 0 the observed confirmatory
  increment, columns 1-1000 the permutations), arms `M_unrestricted`, `M_within_operator`,
  `M_within_operator_residualized`. `server_export/scripts/pinned/luad_c1/c1_luad_primary_primary.npz`
  (`2034e711...`): 10,406 usable genes x 1,001 columns, the registered H1 arm.
- **Published results, for validation.** `server_export/scripts/pinned/c1_ucec_d4_followup.json`;
  `server_export/scripts/pinned/luad_c1/c1_luad_primary_result_primary.json`;
  `server_export/results/results__composition_controls_ucec.csv` and `..._luad.csv`;
  `server_export/results/scripts__transcriptome_baseline_ucec_perm.csv` and `..._luad_perm.csv`.
- **The batch-corrected sets.** `server_export/results/ucec__sitepack_operator.csv` and
  `server_export/results/luad__sitepack_operator.csv` (columns `fdr_batchresid`,
  `incremental_r2_batchresid`).
- **Code the definitions come from:** `server_export/scripts/residual_analysis.py`
  (`a54a0a56...`), `split_replication_power.py` (`79530c07...`, at the paper root),
  `server_export/scripts/c1_run_test.py` (`fc578c70...`),
  `server_export/scripts/c1_run_test_luad.py` (`1d5ee05d...`),
  `server_export/scripts/composition_controls.py` (`46b4a1d8...`),
  `server_export/scripts/transcriptome_baseline.py` (`761df7d7...`),
  `server_export/scripts/residual_analysis_sitepack.py` (`d0cfeb01...`),
  `review/MEAS_ERROR_2026-10-08/theta_kernel.py` (`6c95cd11...`).
- **Files of this analysis**, in `review/C1_COVARIATES_2026-11-21/`:
  `c1_covariates_local.py`, `c1_covariates_cluster.py`, the LSF wrappers
  `lsf_c1_cov_*.sh`, `make_abstract_variants.py`, its input
  `review/QUALITY_REVIEW_2026-11-20/abstract_v12/ABSTRACT_FINAL.tex` (`6a1bfede...`, the
  Entry 12 Abstract) and its output `ABSTRACT_VARIANTS.md`, `CLOSED_PARTS.txt` and `PINNED_HOSTS.txt`
  (section 9), `make_sidecar.py` with its list `SIDECAR_FILES.txt`, and the gate logs in
  `gate_logs/`.
- **On the cluster only** (`/public/home/fjhui/ZW/`): the RNA-seq, protein and slide
  embeddings of both confirmatory cohorts and the full discovery RNA-seq matrices.

## 2. Covariates

Each covariate set is computed inside one cohort, on that cohort's analysed cases only and
in the published case order (sorted case identifiers): the 100 and 105 discovery patients
of the inputs files, and the 138 and 112 analysed confirmatory cases. `rna` below is
`residual_analysis.load_rna_matrix()` (every gene; log2(TPM+1); tumour files averaged per
case) restricted to those cases; no gene is removed, and the target gene is not removed from
any covariate, as in the published controls.

- **C3, three transcriptomic composition scores** (columns stroma, immune, proliferation,
  in this order): `composition_controls.comp_score(rna, composition_controls.MARKERS[f])`,
  imported unchanged. For each family the markers present as columns of `rna` are z-scored
  over the cases (pandas mean and standard deviation, ddof 1; a marker with zero standard
  deviation is skipped) and averaged. C3 is untestable in a cohort if a family has fewer
  than 6 markers present or no marker with a non-zero standard deviation (the published
  control would have dropped the family; here the three-score definition is kept whole).
  The markers used per family are reported.
- **R20, 20 expression components**, as `transcriptome_baseline.py` and the export that
  built `rna_pcs20`: `rv = rna.values.astype(float); sd = rv.std(0); rz = (rv - rv.mean(0)) /
  np.where(sd == 0, 1.0, sd); R20 = PCA(n_components=20, svd_solver="full").fit_transform(rz)`.
  In discovery R20 is the `rna_pcs20` array of the inputs file. The number stays 20 in the
  confirmatory cohorts.
- **OP, operator indicators (discovery only):** one indicator for each of the 12 largest
  scanning-operator levels with at least 3 cases; every other case, including the one
  unlabelled case in each cohort, has an all-zero row. Binding counts, from the labels
  alone: UCEC 12 indicators covering 73 of 100 cases; LUAD 12 covering 75 of 105. A tie at
  the cut-off would be fatal; there is none.
- **Capacity placebos R20cap_j (confirmatory only), j = 1..119:** `R20cap_j = R20[q_j]`
  (row i of the placebo is row `q_j[i]` of R20), with `q_j =
  RandomState(crc32(b"c1_cov|<cohort>|cap|%03d" % j)).permutation(n)`, `<cohort>` being
  `ucec` or `luad` and n 138 or 112. Each placebo has the size and covariance of R20 and
  carries no patient information.

Transcriptomic scores are the only composition measure in these four data sets; no
pathologist-read composition exists for them (Limitations (iv)). R20 is the published
baseline, with every gene standardised: it is not a tuned reader of tissue state, and it
also carries a de-noised copy of the gene's own transcript (Limitations (v)). It is called
"expression components", never a tissue covariate.

## 3. Discovery analyses and sets

Every discovery analysis uses the published estimator and changes exactly one thing, the
block or the baseline: per gene, on the patients with a protein value, 5-fold
`KFold(shuffle=True, random_state=0)`, ridge (alpha 1) on columns standardised within the
training fold; increment = R2(baseline + block) - R2(baseline); null = the block's rows
permuted over all patients with the published sequence (`RandomState(0)`, 1,000 sequential
draws), the baseline staying with its patients; one-sided p = (1 + #{null >= observed}) /
1001.

| Analysis | Baseline | Block (permuted) |
|---|---|---|
| M (published) | own mRNA | 20 morphology components |
| C3 | own mRNA | C3 |
| R20 | own mRNA | R20 |
| OP | own mRNA | OP |
| M\|C3 | own mRNA + C3 | 20 morphology components |
| M\|R20 | own mRNA + R20 | 20 morphology components |

The published permutation analysis under the R20 baseline (`transcriptome_baseline.py`)
drew its permutations as `RandomState(i).permutation(n)`, i = 0..999, and only for the genes
of the published sets. Analysis M|R20 uses the sequence of the main analysis for every fit
gene, so its p-values and floor counts are not those of the published table (99; 55), and
the supplementary table says so.

- **S_M**, the published uncorrected sets (2,566 UCEC; 2,310 LUAD): FDR < 0.05 and
  increment > 0 in analysis M. They are the only sets read in the primary and reverse cells.
- **Sets that are read, T_X** for X in C3, R20, OP, M|C3, M|R20: the n genes ranked first by
  analysis X, with n the size of S_M in that cohort (2,566; 2,310), ordered by p ascending,
  then by increment descending, then by gene symbol. No sign filter is applied. A fixed size
  is used because Benjamini-Hochberg at 1,000 permutations returns either no gene or at
  least k\* = ceil(m / (0.05 x 1001)) genes (211; 215), which is expected to empty the
  FDR-defined sets under the adjusted baselines (section 0, item 3), and because sets of one
  size can be set side by side. Two consequences are reported for every set and not
  corrected: when more than n genes of an analysis sit at the permutation floor, the set is
  the n floor genes with the largest increments, that is, it is ordered by increment alone;
  and when few genes carry signal the set runs into large p-values and is mostly noise,
  which lowers its AUC.
- **A rule reference, reported and not read, T_M:** the same rule applied to the published
  analysis M. It differs from S_M by the rule alone (S_M has a sign filter and an FDR
  threshold), and so separates what the rule does from what the block or the baseline does.
- **Sets that are reported and not read, S_X:** FDR < 0.05 (Benjamini-Hochberg over all fit
  genes) and increment > 0 in analysis X, with the number selected, the number of rejections
  before the sign filter, the number of genes at the permutation floor and k\*. An S_X with
  no rejection is "empty at B = 1,000, limited by the permutation floor"; one with
  rejections but no positive increment is reported as that.
- **The batch-corrected sets, reported and not read, S_BC:** `fdr_batchresid < 0.05` and
  `incremental_r2_batchresid > 0` in the operator-axis discovery table: 201 genes in UCEC
  (the Batch-corr. column of the published Table 1) and 726 in LUAD (the set of the
  registered step 3).
- **A gene-property set, reported and not read, T_G:** the n genes with the lowest
  published own-mRNA R2 in discovery (`pub_r2_rna`, ascending, then gene symbol). It asks
  whether genes poorly predicted by their own transcript rank above background whatever
  selects them, under histology and under the covariates alone.

## 4. Confirmatory arms (Stage H, on the cluster)

Each arm keeps everything of the published test: `cases = sorted(common)` from the published
script's own `prepare` (138; 112), its tested genes, MIN_N = 30, `residual_analysis.cv_r2`,
the morphology components `pcs = PCA(n_components=K, svd_solver="full").fit_transform(W)`
with `W = wsi.loc[cases].values.astype(float)` and K = 28 (UCEC) or 22 (LUAD), and the
published permutations: `rng = RandomState(crc32(b"c1_ucec|perm"))` or
`RandomState(crc32(b"c1_luad|perm"))`, `perms = [rng.permutation(n) for _ in range(1000)]`,
with no other draw from that generator. Permutation b is therefore the same reassignment of
rows in every arm and in the published matrix.

**The per-gene function** (a copy of `split_replication_power._gene_test` with two
additions, the baseline covariates E and a general block B). For gene j with fit mask `v`
(protein not missing; the gene is skipped if `v.sum() < 30`), own-mRNA column `xr`,
protein `y`, E of shape (n, q) and B of shape (n, k):

    r2b      = cv_r2(hstack([xr[v], E[v]]), y[v])
    column 0 = cv_r2(hstack([xr[v], E[v], B[v]]), y[v]) - r2b
    column b = cv_r2(hstack([xr[v], E[v], B[perms[b-1]][v]]), y[v]) - r2b

`B[p]` replaces row i of the block by row `p[i]`. E has no columns in the V arms and in the
gate that reproduces the published matrix.

| Arm | E (stays with its patients) | B (rows permuted under the null) | Columns |
|---|---|---|---|
| A_C3 | C3 | pcs | 0 and 1-1000 |
| A_R20 | R20 | pcs | 0 and 1-1000 |
| A_R20cap_j, j = 1..119 | R20cap_j | pcs | 0 and 1-25 |
| V_C3 | none | C3 | 0 and 1-1000 |
| V_R20 | none | R20 | 0 and 1-1000 |

Every arm is computed on exactly the published usable gene list (the `genes` array of the
saved matrix); no gene is dropped, imputed or re-selected, and a non-finite value is a fatal
error. The design is a ridge with alpha 1 on standardised columns, so no gene fails for
rank.

**Files.** For `<cohort>` in `ucec`, `luad`: `c1_cov_<cohort>_<arm>.npz` for the four
1,001-column arms and `c1_cov_<cohort>_A_R20cap_<jjj>.npz` for each placebo, each with
exactly the keys `genes` (unicode, equal in value and order to the saved matrix's), `M`
(float64; genes x 1,001, or genes x 26) and `cases`;
`c1_cov_confirmatory_covariates_<cohort>.npz` with `cases`, `C3`, `R20`, `pcs` and the
markers used; and one `c1_cov_stage_h_manifest.sha256` with the sha256 of every output of
both cohorts, whose file names are unique across the two.

## 5. Cells: statistics, background and null

A cell is one set read on one matrix. The statistic is the published one: the AUC of the
matrix's column-0 increments of the set against the background (rank Mann-Whitney) by
`assemble`, copied unchanged from `c1_run_test.py`; the null is the same statistic on the
permutation columns; one-sided p = (1 + #{null >= observed}) / (1 + B); margin = AUC minus
the null mean; z = margin / null SD. The null mean is never assumed to be 0.5. Spearman rho
over the universe between the selecting analysis's discovery increment and the confirmatory
increment is reported, as in the published tests, for every cell whose set is chosen by a
discovery increment (not for T_G or partition cells), and is not read.

- **Universe and background.** The usable confirmatory genes of the published test (10,146;
  10,406). The background of a set is every other gene of the universe (the C1-LUAD H2
  rule). A set or background with fewer than 50 genes in the universe is "untestable".

| Family | Cell (set on matrix), in each cohort | Role |
|---|---|---|
| Primary | S_M on A_C3; S_M on A_R20 | read (6.1) |
| Capacity | S_M on each A_R20cap_j | enters the R20 class (6.1) |
| Reverse | S_M on V_C3; S_M on V_R20 | read (6.2) |
| Forward | T_C3, T_R20, T_OP on the published matrix | read (6.2) |
| Matched | T_M\|C3 on A_C3; T_M\|R20 on A_R20 | read (6.2) |
| Covariate-to-covariate | T_C3 on V_C3; T_R20 on V_R20 | reported, not read |
| Rule reference | T_M on the published matrix, on A_C3 and on A_R20 | reported, not read |
| FDR-defined | S_C3, S_R20, S_OP on the published matrix; S_M\|C3 on A_C3; S_M\|R20 on A_R20 | reported, not read |
| Gene property | T_G on the published matrix, on V_C3 and on V_R20 | reported, not read |
| Batch-corrected sets | S_BC on the published matrix, on A_C3 and on A_R20 | reported, not read |
| Partition | for X in C3, R20, OP, on the published matrix: the genes in both S_M and T_X, in S_M only and in T_X only, each against the genes in neither | reported, not read |
| Reference | S_M on the published matrix | reproduces the published reading |

Cells on the published matrix use its unrestricted arm. The S_BC cells use the published
estimator with unresidualised components and the unrestricted null, which is not the
estimator of the 726-gene set's registered reading (operator-residualised components;
within-operator null). The 201-gene set's published reading required both the unrestricted
and the within-operator null on unresidualised components, so its cell on the published
matrix repeats the published unrestricted value. The cells show whether the adjustment
moves these sets as it moves S_M. For C1-UCEC the forward,
rule-reference, gene-property, partition and reference cells on the published matrix are
also computed on the two other saved arms (within-operator null; operator-residualised
components), reported and not read. Partition cells have a different background from every
other family and are not comparable with them; if T_X has fewer than 50 genes in the
universe the three parts are "partition undefined".

**The change from the published matrix (reported, not read).** For X in A_C3, A_R20 and the
set S_M, with the shared permutations: D_b = AUC_pub(b) - AUC_X(b) for b = 0..1000; D_0, the
mean and SD of D_1..D_1000, z_D, one-sided p_D = (1 + #{D_b >= D_0}) / 1001, and the retained
fraction f = margin_X / margin_pub. The reference distribution is exact only under the
hypothesis that slides are exchangeable with respect to protein, mRNA and the covariates
jointly, which the published tests already reject. p_D therefore describes D_0 against
permutation noise; it is not a test that the two arms carry equal signal, and for A_R20 it
does not separate what the components carry from what 20 extra predictors cost. It supports
no statement that a covariate explains an ordering.

**The capacity randomisation (enters the R20 class).** For each placebo j, m_j = AUC of S_M
on A_R20cap_j in column 0 minus its mean over columns 1-25; m_R is the same quantity from
columns 0 and 1-25 of A_R20. p_abs = (1 + #{j : m_j <= m_R}) / 120. It tests whether the
assignment of the 20 components to patients matters for the margin. m_R uses the first 25
permutation columns and can differ in the third decimal from the cell's margin, which uses
all 1,000; the supplementary table prints both. Real components also lower the residual
variance, which a placebo does not, so the margin under the real assignment is expected, if
anything, to be the larger one; the size of this tilt is not estimated. Real components are
also correlated with the morphology components in the observed column and not in the
permuted ones, which changes the cost of the morphology block in column 0 in a direction
that depends on the set's null mean. With 119 draws a real assignment whose randomisation p
is 0.02 is classed A about nine times in ten, and one whose p is 0.04 about half the time.
The tilt and the resolution both favour class U, which changes no headline wording.

**A diagnostic, not read.** The canonical correlations between the confirmatory morphology
components and C3, and R20, by the function `_cancorr` of `transcriptome_baseline.py`, from
the covariates file of section 4, and as a reference those of a standard Gaussian block of
the size of the morphology components (`RandomState(42).randn(n, K)`) with C3 and with R20.

**What the margins are not.** Sets differ in how they were chosen, and every discovery
selection carries its own selection effect. Margins are compared across cells descriptively
only. No cell is called an effect size, and no sentence says that one selector beats
another.

## 6. Reading rules

p is used unrounded. "Outranked" means p < 0.05; "not shown" means p >= 0.05. In the
Abstract only, "did not" stands for "was not shown to" (p >= 0.05); it claims neither that
the set ranked below background nor that the covariates explain the ordering. The tables
also print the tier of the registered rule (p < 0.05; 0.05 <= p < 0.15; p >= 0.15).

### 6.1 Primary readings (four: one per cohort and covariate set; no correction)

- **C3 class** of a cohort: **H** if S_M on A_C3 has p < 0.05; **N** otherwise.
- **R20 class** of a cohort, from p_R (S_M on A_R20) and p_abs: **H** if p_R < 0.05; **A**
  if p_R >= 0.05 and p_abs < 0.05 (with 119 placebos, at most four placebo margins at or
  below m_R); **U** otherwise. Three covariates cost almost no capacity, so C3 has no
  placebo; 20 do (43 predictors at about 90 training patients in LUAD, 49 at about 110 in
  UCEC).
- **Closed readings.** A primary reading can be closed without a result only before
  freezing, by the marker rule of section 2 or by a gate of section 9; `CLOSED_PARTS.txt`
  and the gate record name every closed part, and the file holds no entry if nothing is
  closed.
  After freezing no primary reading is closed: a run stopped by the scheduler (wall-clock
  limit, node failure, queue) is resubmitted unchanged until it completes, the wall-clock
  limit being the only setting that may be raised; an exception is handled under the Defects
  rule of section 9. The manuscript is not submitted until every arm not closed before
  freezing is complete, Stages L1 and L2 have run and the results record is released.
- **Type** of a cohort, from its C3 reading (rows) and its R20 reading (columns):

  | C3 \ R20 | H | A | U | closed |
  |---|---|---|---|---|
  | H | none | R | none | NC:R |
  | N | C | C | C | C |
  | closed | NC:C | R | NC:C | NC:C |
- **No correction, and why.** The four readings are four questions and are never summarised
  as one; each is printed with its p. Every consequence for the Abstract, the title and the
  opening sentences follows from a failure, for which leaving p uncorrected is the stricter
  choice. The favourable sentences say that the readings are uncorrected.

### 6.2 Secondary readings (each cell: "outranked" or "not shown"; counts are reported)

Reverse: four cells. Forward: six cells. Matched: four cells. No correction. An outcome
that outranks is unfavourable to the paper in the reverse and forward families and
favourable in the matched family; none of these cells changes the Abstract, the title or
the opening sentences.

### 6.3 What a class means, fixed now

- **C3 class H:** histology adds to these three transcriptomic scores for the published
  set's ordering. It does not show independence from composition measured in other ways.
- **C3 class N:** histology's contribution to the ordering could not be told apart from
  what the scores carry. Three covariates do not cost the power to see it: from the
  published margins and null SDs, N requires the margin to fall below about 0.6 of its
  published value in UCEC and about 0.26 in LUAD (0.61 and 0.25 from the empirical nulls).
  The two cohorts therefore do not have the same chance of N for the same proportional
  loss.
- **R20 class A:** the ordering is not shown with the 20 components in the baseline, and
  the real components leave less of the margin than random reassignments of them do.
- **R20 class U:** the set is not shown to outrank background with the 20 components in the
  baseline, and the comparison does not say whether the components or the 20 extra
  predictors are the reason. The question is unresolved in that cohort; it is not answered
  in the paper's favour.
- A reverse cell that outranks: the set's ordering can be obtained in the confirmatory
  cohort without a slide. A forward cell that outranks: genes chosen without histology also
  receive larger histology increments in new patients. Both are what one expects if
  histology reads tissue state that the covariates also record; neither shows that
  histology is redundant, and a cell that does not outrank does not show that it is needed.
- A matched cell that outranks: genes chosen by histology after the covariates entered the
  discovery baseline also rank above background when the confirmatory baseline holds
  covariates built in the same way.

### 6.4 Cases

A cohort is **flagged** if it has a type. **Case 0:** no cohort flagged. **Case 1:** one.
**Case 2:** both.

## 7. Manuscript wording, fixed now (filled from the outputs after Stage L2, never before)

[a], [m], [s], [p], [d] are the AUC, null mean, null SD, one-sided p and margin of the named
cell, printed to three decimals (p = 1/1,001 as 0.001); a p below 0.05 that would print as
0.050 is printed with as many decimals as show it to be below 0.05 (50/1,001 as 0.04995),
as in the post hoc pre-specifications of 2026-10-08; classes use unrounded values. [date] is
the UTC date of the sidecar's `# written` line, printed as "[day] [month] 2026 (UTC)". [X],
[Y] and [Z] are "uterus" or "lung". [P] is "transcriptomic
composition scores" for type C or NC:C and "20 expression components" for type R or NC:R.
No Stage L1 number is written into the manuscript before Stage L2 completes.

**The sentence QS**, used in 7.6 to 7.8, exactly one of:

- One cohort of type C or R (the other not flagged): "Post hoc, the uncorrected
  [uterine|lung] set was not shown to outrank background once [P] entered the confirmatory
  baseline (\S\ref{sec:c1ucecres})."
- Both of one type C or R: "Post hoc, neither uncorrected set was shown to outrank
  background once [P] entered the confirmatory baseline (\S\ref{sec:c1ucecres})."
- Types C and R: "Post hoc, the uncorrected uterine set was not shown to outrank
  background once [P of uterus] entered the confirmatory baseline, nor the uncorrected lung
  set once [P of lung] did (\S\ref{sec:c1ucecres})."
- One cohort of type NC (the other not flagged): "Post hoc, the test of the uncorrected
  [uterine|lung] set with [P] in the confirmatory baseline was not run
  (\S\ref{sec:c1ucecres})."
- Both of type NC: "Post hoc, the tests of the uncorrected sets with [P] in the
  confirmatory baseline were not run (\S\ref{sec:c1ucecres})."; with different [P]: "...
  with [P of uterus] (uterus) and [P of lung] (lung) in the confirmatory baseline were not
  run (\S\ref{sec:c1ucecres})."
- One cohort of type C or R and one of type NC: the first form for the one, followed by the
  fourth form for the other.

### 7.1 Results: one new paragraph at the end of `\S\ref{sec:c1ucecres}` (after "...so it bounds protein batch rather than isolating it.")

The paragraph is these sentences, in this order and in no other.

**M1 (C3), exactly one of.** Every form opens: "Post hoc, in an analysis specified on
[date] in a dated file released publicly before it ran (Supplementary
Note~\ref{snote:c1cov}), we added transcriptomic stromal, immune and proliferation scores
to the confirmatory baseline: " and continues:

- both cohorts H: "both uncorrected discovery sets (2{,}566; 2{,}310 genes) still outranked
  background (uterus: AUC [a] against [m] $\pm$ [s], one-sided $p=$[p], margin [d] against
  $0.099$ without the scores; lung: [a] against [m] $\pm$ [s], $p=$[p], margin [d] against
  $0.120$)."
- one H (X), one N (Y): "the uncorrected discovery set still outranked background in [X]
  ([2{,}566|2{,}310] genes; AUC [a] against [m] $\pm$ [s], one-sided $p=$[p], margin [d]
  against [$0.099$|$0.120$] without the scores) but was not shown to in [Y]
  ([2{,}566|2{,}310] genes; [a] against [m] $\pm$ [s], $p=$[p], margin [d] against
  [$0.099$|$0.120$])."
- both N: "neither uncorrected discovery set (2{,}566; 2{,}310 genes) was shown to outrank
  background (uterus: AUC [a] against [m] $\pm$ [s], one-sided $p=$[p], margin [d] against
  $0.099$ without the scores; lung: [a] against [m] $\pm$ [s], $p=$[p], margin [d] against
  $0.120$)."
- one cohort closed (Z), the other H or N (X): "in [X] the uncorrected discovery set
  ([2{,}566|2{,}310] genes) [still outranked background|was not shown to outrank background]
  (AUC [a] against [m] $\pm$ [s], one-sided $p=$[p], margin [d] against [$0.099$|$0.120$]
  without the scores); in [Z] the test was not run."
- both closed: the whole of M1 is "Post hoc, an analysis specified on [date] in a dated
  file released publicly before it ran (Supplementary Note~\ref{snote:c1cov}) was to add
  transcriptomic stromal, immune and proliferation scores to the confirmatory baseline in
  the tests of the uncorrected discovery sets (2{,}566; 2{,}310 genes); that test was run in
  neither cohort."

The published margins are 0.0992 (UCEC) and 0.1195 (LUAD), from the published AUC and null
mean.

**M2 (R20):** "With 20 expression components added instead, the two uncorrected discovery
sets gave one-sided $p=$[p] (uterus; margin [d]) and [p] (lung; margin [d]); [e] of 119
(uterus) and [e] of 119 (lung) placebo baselines of the same size left a margin as small or
smaller[tail]." [e] = #{j : m_j <= m_R}. The tail is, for each cohort in the order uterus,
lung: nothing if its R20 class is H; "; in [uterus|lung] the real components therefore took
more of the margin than 20 extra predictors cost" if A; "; in [uterus|lung] what the real
components took is therefore not separated from what 20 extra predictors cost" if U. If the
R20 reading of one cohort (Z) is closed, M2 is "With 20 expression components added
instead, the uncorrected [uterine|lung] set gave one-sided $p=$[p] (margin [d]), and [e] of
119 placebo baselines of the same size left a margin as small or smaller[tail]; in [Z] this
test was not run." If both are closed, M2 is "The test with 20 expression components added
instead was run in neither cohort."

**M3 (reverse):** "In place of histology, the three scores alone gave one-sided $p=$[p]
(uterus) and [p] (lung) for the same sets, and the 20 components alone [p] and [p]." A
covariate set whose reverse cell is closed in one cohort (Z) has the clause "[the three
scores|the 20 components] alone gave one-sided $p=$[p] in [X] and were not run in [Z]";
closed in both cohorts, "[the three scores|the 20 components] alone were not run".

**M4 (forward and matched):** "Of six discovery sets chosen without histology, each the
size of its cohort's published set (by the scores, the components or scanning operator),
[k] outranked background under
confirmatory histology; of four chosen by histology after the covariates had entered the
discovery baseline, [k] did so when the confirmatory baseline held covariates built in the
same way (the p-value at which each set was cut, and the same rule applied to histology
alone, are in Supplementary Table~\ref{tab:c1covsets})." A count is given "of [t] that
could be tested" if a cell is untestable or closed.

**M5 (the fixed limits; verbatim in every outcome):** "The covariates are transcriptomic
estimates, the null is unrestricted and the four readings are uncorrected. A set that still
outranks background shows that histology adds to these estimates, not that the ordering is
independent of composition measured in other ways or robust to operator; a set not shown to
means that histology's contribution could not be told apart from what the covariates carry,
as expected if histology reads the same tissue state, not that histology reads nothing
else. Covariates that order a set without a slide do not show that histology is redundant,
and covariates that fail to do not show that it is needed; nor does a set chosen without
histology that fails to outrank background, because sets of fixed size are diluted when
their analysis finds few genes; the batch-corrected sets were re-read descriptively only
(Supplementary Table~\ref{tab:c1cov})." If a primary reading is closed, "the four readings
are" becomes "the [three|two] readings that were run are" ("the one reading that was run
is"), here and in section 7.2, and the paragraph of section 7.2 gains the sentence "In [Z]
the test with [the scores|the components] in the baseline was not run[cause]." with the
cause clause of section 9.

If the main text then exceeds the journal's length, free text elsewhere is shortened; none
of M1 to M5 is.

### 7.2 Methods

One sentence appended to `\S\ref{sec:c1ucecmeth}`, after "...may combine two plexes
removed.": "A fifth post hoc analysis, specified on [date] in a dated file released
publicly before it ran, re-reads this test and the lung test with transcriptome-derived
covariates (\S\ref{sec:c1luadmeth}; Supplementary Note~\ref{snote:c1cov})."

One paragraph appended to `\S\ref{sec:c1luadmeth}`: "\emph{Transcriptome-derived
covariates in the two confirmatory tests (post hoc).} This analysis was specified on [date]
in a dated file released publicly before it ran (Supplementary Note~\ref{snote:c1cov}). In
each confirmatory cohort the published test was repeated for the uncorrected discovery set
with the baseline enlarged by three transcriptomic composition scores (stromal, immune,
proliferation; \S\ref{sec:stroma}) or by 20 expression components
(\S\ref{sec:txomemeth}), each computed on the confirmatory cases, and with the scores or
the components in place of the morphology components. To separate what the 20 components
carry from what 20 extra predictors cost, the components were also reassigned to patients
at random 119 times, and the margin of the real assignment was ranked among the 119. Cases,
genes, morphology components, folds and permutations were those of the published tests;
only the unrestricted null was run and the four readings are uncorrected. In discovery,
genes were ranked by the published estimator with the scores, the components or
scanning-operator indicators in place of morphology, and with morphology after the scores
or the components had entered the baseline; the sets read hold as many genes as the
published set. The batch-corrected discovery sets were also evaluated, descriptively, on the published
test and with the scores or the components in the baseline."

### 7.3 Supplement: Supplementary Note `snote:c1cov` and Supplementary Tables `tab:c1cov` and `tab:c1covsets`

`tab:c1cov` has one row per cell of section 5, in every outcome: cohort, family, set,
matrix, null, whether the cell is read, genes in the set and in the background, AUC, null
mean, null SD, margin, z, one-sided p, tier and rho; and, for the primary cells, D_0, its
null mean and SD, z_D, p_D, f, the class, and for R20 m_R, the 119 values m_j and p_abs.
`tab:c1covsets` gives, for every T set chosen by an analysis (T_X and T_M; for T_G only
the overlaps): its overlap with S_M and with T_M, the largest p in the set, the Benjamini-Hochberg q at that rank, the number of its genes at the permutation
floor, whether the floor alone filled the set, and the number of its genes with a
non-positive discovery increment; for every S set its size, rejections before the sign
filter, floor count and k\*; the note of section 3 on the permutation sequence of M|R20;
and the markers used per composition family in each of the four computations. The note
states sections 0, 5, 6.3, 8 and the order of section 9 in substance, gives the diagnostic
of section 5, and contains these sentences verbatim:

- "The capacity comparison ranks the margin under the real assignment of the 20 components
  among 119 random reassignments; it tests whether the assignment of the components to
  patients matters. Real components also lower the residual variance, which a reassignment
  does not, so the comparison is expected, if anything, to understate what the components
  take, and with 119 reassignments it can miss a moderate loss."
- "Margins are not compared between selectors: the sets were chosen in different ways and
  each selection carries its own selection effect."
- "In the Abstract, 'did not' stands for 'was not shown to' (one-sided $p\geq0.05$); it
  claims neither that the set ranked below background nor that the covariates explain the
  ordering."

### 7.4 Abstract

**Policy** (that of the post hoc pre-specifications of 2026-10-08, as Entry 12 reads it).
No favourable post hoc result enters the Abstract. The only Abstract edits of this
analysis are the qualifier sentence below, the cuts that make room for it and, in Case 2,
the opening sentence (7.5). The qualifier reports an outcome unfavourable to the paper's reading
and follows the sentence on the two later cohorts, which a reader would otherwise take to
show an ordering that histology adds to what the transcriptome records. If the Abstract is
later shortened, or the target journal sets another limit, the qualifier sentence is kept
word for word; if a limit makes that impossible, a dated entry of `C1_LUAD_ADDENDA.md`
records the change, and the sentence still says that the analysis was post hoc, which
covariates entered and which uncorrected set was not shown to outrank background.

**The texts.** `ABSTRACT_VARIANTS.md`, generated by `make_abstract_variants.py` and listed
in the sidecar, holds the complete Abstract for each of the 24 flagged outcomes and its
length under the four counting conventions of Entries 11 and 12 (141 / 143 / 149 / 151 for
the Entry 12 Abstract).

- **Case 0:** the Abstract is unchanged. Reason, fixed now: no favourable post hoc result
  enters the Abstract.
- **Case 1** (cohort of type C or R): "Post hoc, adjusted for [composition
  scores|expression components], the uncorrected [uterine|lung] set did not." is inserted
  before "Elsewhere,".
- **Case 2**, both of one type: "Post hoc, adjusted for [composition scores|expression
  components], neither uncorrected set did." Types C and R: "Post hoc, adjusted (uterus,
  [..]; lung, [..]), neither uncorrected set did."
- **A cohort of type NC:C or NC:R** has "the test of the uncorrected [uterine|lung] set
  adjusted for [..] was not run", in the forms given in `ABSTRACT_VARIANTS.md`.

**Room, fixed now.** Words are removed in this order until the Abstract is no longer than
the Entry 12 Abstract under each of the four conventions, and no further: (c1) "At these
cohort sizes,"; (c2, Case 1 only) "Here, we show that"; (c3) "Elsewhere,"; (c5) "genes"
after "2,310"; (c6) "patients" after "138". In Case 2 the opening sentence is replaced
first (7.5), and the cuts are then counted. Labels c4 and c7 of the second draft are
retired. No protected content item is removed (section 0, departure 3); "later", "from one
consortium" and Entry 12's clause "morphology adds no net accuracy over transcriptome-wide
baselines in four cohorts" are not touched. The resulting lengths: the two Case 2 variants
of one type equal the Entry 12 Abstract (141 / 143 / 149 / 151); the four Case 1 variants
and the two Case 2 variants of types C and R are one word longer (142 / 144 / 150 / 152)
after every cut. The sixteen variants with a closed reading are 3 to 14 words longer; they
can arise only if a reading is closed before freezing, in which case the gate record names
the variants that remain possible and their lengths. A variant is submitted at its stated
length; if the journal requires a cut, it is made by a dated entry of `C1_LUAD_ADDENDA.md`
that removes no protected content item and no word of the qualifier.

The Abstract edit, and any change under 7.5, is recorded in a dated entry of
`C1_LUAD_ADDENDA.md` that quotes the variant used and is released publicly, as an annotated
tag and a GitHub Release (it may be part of the results record), before the manuscript is
edited.

### 7.5 Title and opening sentences (the fork, fixed now)

- **Case 0 and Case 1:** the title and the opening sentences of the Introduction's summary
  ("Here, we show that histology carries a set-level signal of which proteins depart from
  their own transcript.") and of the Discussion are unchanged. In Case 1 the Abstract's
  opening sentence loses "Here, we show that" (cut c2); the cut is made for length and is
  not a change of claim.
- **Case 2:** the title becomes "Testing whether histology signals which proteins depart
  from their own transcript in five human cancers" (in the main text and on the title page
  of the Supplementary Information). The Abstract's opening sentence becomes
  "We test whether histology signals which proteins depart from their own transcript." The
  Introduction's sentence becomes "Here, we test whether histology signals which proteins
  depart from their own transcript." The Discussion's opening sentence becomes "Within one
  consortium, histology carries a set-level signal of which proteins depart from their own
  transcript; in two later cohorts the uncorrected sets were not shown to outrank
  background once [transcriptomic composition scores|20 expression
  components|transcriptome-derived covariates] entered the baseline (post hoc)." (both of
  type C; both of type R; types C and R), and
  its next sentence begins "The signal recurs across five organs, and" in place of "It
  recurs across five organs within one consortium, and". With both cohorts of type NC the
  second half reads "; in two later cohorts a test of the uncorrected sets with
  transcriptome-derived covariates in the baseline was not run (post hoc)."; with one of
  type NC, "; in one later cohort the uncorrected set was not shown to outrank background
  once [transcriptomic composition scores|20 expression components] entered the baseline,
  and in the other a test of that kind was not run (post hoc).".
- **The target journal** is the authors' decision after the results. Any further change to
  the framing may weaken the paper's claims and may not strengthen them; it is recorded in a
  dated entry of `C1_LUAD_ADDENDA.md`.

### 7.6 Introduction

- In every case the list in the last sentence gains an item: "...a post hoc test for
  transcript measurement error, a post hoc test against transcriptome-derived covariates
  and the nulls outside CPTAC-3 are all reported."
- Cases 1 and 2: QS is inserted after "...it also passed an acquisition-batch
  (operator)-stratified null."

### 7.7 Discussion

- Cases 1 and 2, after the sentence that ends "...the signal held unevenly
  (\S\ref{sec:c1res}).": QS; then, for a cohort of type C or R, "[In that cohort|In those
  cohorts] the ordering does not distinguish histology from [a transcriptome-derived estimate of tissue composition |
  what the transcriptome as a whole records | what the transcriptome records]." (type C;
  type R; types C and R). With one cohort of type C or R and one of type NC the order is:
  the result sentence of QS, the "In that cohort ..." sentence, the "not run" sentence of
  QS.
- Cases 1 and 2: the sentence that ends "...outranked background genes; the test ranked
  gene sets, not the protein of individual genes." is followed by QS.
- Any cohort with R20 class U, in every case: "With 20 expression components in the
  confirmatory baseline [the uncorrected uterine set was not|the uncorrected lung set was
  not|neither uncorrected set was] shown to outrank background, and the placebo
  comparison did not separate this from what 20 extra predictors cost; whether histology
  adds there to what the transcriptome as a whole records is therefore unresolved
  (\S\ref{sec:c1ucecres})." is inserted after the sentence that ends
  "...not the protein of individual genes.", after QS if QS is present there.
- Case 0 with no class U: the Discussion is unchanged.

### 7.8 Limitations (iv), in the main text and in the matching paragraph of the Supplementary limitations note: sentences appended, in this order

- Cases 1 and 2: QS. Case 0: "Post hoc, both uncorrected sets still outranked background
  with transcriptomic composition scores in the confirmatory baseline, which does not make the ordering
  independent of composition measured in other ways (\S\ref{sec:c1ucecres})."
- In every case, if a reverse cell outranks: "Without a slide, [the three scores alone
  ordered [both uncorrected sets|the uncorrected uterine set|the uncorrected lung set]
  above background][ and ][the 20 expression components alone ordered [both uncorrected
  sets|the uncorrected uterine set|the uncorrected lung set] above background]
  (\S\ref{sec:c1ucecres})." The first part is used if a reverse C3 cell outranks, the second
  if a reverse R20 cell outranks; each names the cohorts in which its own cell outranks.
- Any cohort with R20 class U: the class-U sentence of 7.7.

### 7.9 Figure and provenance table

- Fig. 6a gains, in every outcome and tagged post hoc, one row for each of the four primary
  cells and one row per cohort for the placebo margins.
- The supplementary provenance table gains a row for this file, for the Stage H record and
  for the results record.

## 8. What this analysis cannot establish

- **Independence from composition.** The scores and components come from the same RNA-seq
  as the baseline, and the scores track pathologist-read composition weakly in the one
  cohort where both exist (PDAC). Score error can leave a class H where a better measure of
  composition would not.
- **That a covariate explains an ordering.** Nothing here measures how much of a margin a
  covariate removes; the change D is printed with its own permutation spread and supports
  no such statement.
- **A conditional null.** The permutation null breaks the slide from everything else. In
  the observed column morphology and the covariates are dependent; in every permuted column
  they are independent. The test is not an exact test of "no association between histology
  and protein given the covariates".
- **Robustness to operator or plex** of any new arm: only the unrestricted null is run.
- **The 20 components in the uterine cohort.** With a published z of 2.73 the uterine test
  has little room: class U is possible there even if the components take much of the
  margin.
- **A sharp capacity comparison.** It is tilted toward class U and has the resolution of
  119 draws (section 5).
- **Undiluted forward and matched sets.** A fixed-size set is mostly noise when its
  analysis finds few genes (section 3).
- **Anything at the gene level.** Most confirmatory increments are negative; every reading
  is about the ordering of a set against its background.
- **Which way composition and histology run.** If histology reads composition and
  composition shapes the protein-to-transcript ratio, reverse and forward cells that outrank
  are what both the paper's reading and the critique predict.
- **Comparisons between selectors** (section 5), **anything beyond description about the
  batch-corrected sets**, and anything about CCRCC, GBM or PDAC.
- **Not done, and why:** 10, 22 or 28 expression components (20 is the published
  definition); a within-operator null for the new arms (each would be a further arm and a
  further fork); a placebo for C3 (three predictors).

## 9. Freeze, gates, stages and run order

**Scripts.** Local: `c1_covariates_local.py`. Cluster: `c1_covariates_cluster.py`, which
imports `residual_analysis`, `split_replication_power`, `c1_run_test`, `c1_run_test_luad`
and `composition_controls` unchanged and asserts their sha256, and is run through wrappers
cloned from the published confirmatory wrappers, with the same data arguments and hash
checks as the published runs. Locally every stage refuses to start unless every sidecar
hash matches. On the cluster a production run refuses unless the sidecar has the sha256
given on its command line; every sidecar entry the cluster run uses (this file, the cluster
script, the three wrappers, the five imported modules, the cohort's saved matrix,
`CLOSED_PARTS.txt` and `PINNED_HOSTS.txt`) matches; the host is the one `PINNED_HOSTS.txt`
names for the cohort; and that host has a gate record, written by the frozen script on the
frozen inputs, that passed strictly or by the fallback of H-V1.

**Before freezing** (synthetic data, published quantities, or RNA-seq alone; none computes
a new increment on a confirmatory cohort or a new discovery analysis). Logs go to
`review/C1_COVARIATES_2026-11-21/gate_logs/`.

- *Self-tests.* On synthetic data: the batched kernel equals the reference kernel
  (`theta_kernel.oof_ridge`, which is `residual_analysis.cv_r2`); `bh_fdr` and `assemble`
  equal direct implementations; set ranking, partition backgrounds, file round-trips and
  the rules of section 6 on constructed cases; the cluster per-gene function equals an
  independent implementation and the local kernel.
- *V1 (local).* The code path used for every discovery analysis, given the published block,
  reproduces the published per-gene increment (to 1e-9), own-mRNA R2 (1e-10), p-value and
  FDR (1e-12) of every fit gene and the published sets (2,566; 2,310) exactly.
- *V1b (local).* With R20 in the baseline the same code reproduces, for every gene of the
  published sets, the published R2 of own mRNA plus the 20 components (`r2_rna_pcs20`;
  1e-9 in both organs) and the published increment of morphology over that baseline
  (`incr_over_pcs20`; 1e-9 in UCEC, 1e-5 in LUAD). The LUAD tolerance is wider because the
  published LUAD transcriptome-baseline run used morphology components that differ from the
  published ones at single-precision level: its own column of increments over own mRNA
  differs from the published increments by up to 1.3e-6, while its own-mRNA R2 agrees to
  1e-15. The largest difference is recorded.
- *V2 (local).* `assemble` on the saved matrices reproduces the published readings (AUC,
  null mean, null SD, p, rho to 1e-12; 2,536 / 7,610 and 2,287 / 8,119 genes).
- *Blind run (cluster; counts only).* Per cohort: cases, K, genes tested and usable, RNA-seq
  columns, all-zero columns and columns with SD below 1e-8, variance explained by the 20
  components, markers found per family and their SDs, embedding files found; and a listing
  of the cluster's result folders for item 7 of section 0.
- *H-V1 (cluster).* The cluster per-gene function with no E and B = the published
  morphology components reproduces columns 0 to 5 of the saved matrix for every usable gene
  (1e-9) with an identical gene list. Hosts are tried in this order, read from the job logs
  on the cluster: UCEC, the host of the job that wrote the saved matrix, then the host of
  the published test, then s006; LUAD, the host of the primary run, then s006. The first
  host that passes is pinned for that cohort and written into `PINNED_HOSTS.txt`. If none passes at 1e-9, the host with the smallest
  largest deviation is pinned provided that deviation is below 1e-6 and the AUC of S_M from
  its recomputed column 0 equals the published AUC to 1e-6; the deviation is recorded.
  Otherwise every Stage H arm of that cohort is closed. The gate record names the host of
  every gate run.
- *H-V2 (cluster).* The R20 construction, applied to the discovery RNA-seq of each organ on
  the patients of the inputs file, on the host whose export produced that file, reproduces
  `rna_pcs20` up to the sign of each column (1e-6). Otherwise A_R20, the placebos and V_R20
  of that organ's confirmatory cohort are closed.
- *H0 and V3.* The cluster exports the discovery C3 scores of both organs, on exactly the
  patients of the inputs files, at full double precision. Locally they must reproduce the
  published composition controls for every gene of the published sets (`stroma_gain`,
  `immune_gain`, `proliferation_gain`, `ALL_gain`; 1e-9; these involve no morphology). If
  C3 is untestable in a discovery cohort, or V3 is not passed, no score table is used and
  the discovery analyses C3 and M|C3 of that organ are closed, with the cells that read
  their sets.
- *Closed parts.* C3 untestable in a confirmatory cohort closes its A_C3 and V_C3. Every
  closed arm, analysis and cell is listed in `CLOSED_PARTS.txt`, which the scripts read: a
  closed arm is not run, a closed cell has the status "closed", and a closed primary
  reading gives the types of section 6.1. No definition is changed to pass a gate. A reading
  is closed under H-V1 or H-V2 only after every host of the gate's list has run the gate and
  failed it; a host that cannot be reached is waited for and is not a failure. The gate
  record gives the log of every attempt, and the corresponding author's go-ahead names each
  closed reading. A closure never removes a flag: every cohort with a closed primary reading
  is flagged (section 6.1), so a closure can only replace "did not" by "was not run"
  (section 0, item 1). The placebo arms are closed only together with A_R20, and only the
  discovery analyses C3 and M|C3 can be closed. For each closed reading the gate record
  fixes a clause of at most twelve words giving the cause; it follows "was not run" or
  "was run in neither cohort" in M1 or M2 and in section 7.2, not in the Abstract. A closed
  form that section 7 does not give (for example for the sentence of section 7.2 on the 119
  reassignments) is fixed in the gate record before the sidecar is written. If `CLOSED_PARTS.txt` has an entry, the
  gate record states, before the sidecar is written, the form that each sentence of section
  7 then takes, using only the closed forms given there, and the Abstract variants that
  remain possible. If it has none, every closed form of section 7 and the sixteen Abstract
  variants with a closed reading are void. A production run that finds a condition for
  closing that `CLOSED_PARTS.txt` does not list stops and is handled as a defect.

The gate record at the end of this section is filled from the logs before the sidecar is
written. The self-test logs, `validate_only.json` and `gate_v3.json` that are frozen are
those written by the frozen versions of the scripts; the gate record gives the script
hashes they carry.

**Release.** After the corresponding author's written go-ahead, recorded with its UTC time:
this file, the sidecar, the scripts and wrappers, `make_abstract_variants.py` with its
input `ABSTRACT_FINAL.tex`, `ABSTRACT_VARIANTS.md`, `CLOSED_PARTS.txt`, `PINNED_HOSTS.txt`,
the gate logs, and the first draft with the review findings are committed, tagged `c1-covariates-posthoc-prespec-2026-11-21` and
published as a GitHub Release. Withheld, as derived data fixed by their sha256 in the
sidecar and deposited with the other derived data at acceptance: the inputs files, the
confirmatory matrices, the discovery score tables and the published per-gene result tables (including the
two operator-axis tables).
The commit message lists them.

**Stage H (cluster; production).** Per cohort, four 1,001-column arms and 119 placebo arms:
one LSF job per 1,001-column arm, and the placebos in one or more jobs, each writing one
file per draw and skipping a draw whose file exists; on the pinned host; each output
written once and never overwritten. Expected 1 to 3 hours per 1,001-column arm, and about
three such arms of computing for the placebos of a cohort. The manifest is written once every arm
not closed is complete and lists every file in the two output folders. The outputs and the
manifest are transferred unchanged to `server_export/results/c1_covariates/stage_h/`, and
the job logs to `server_export/results/c1_covariates/stage_h_joblogs/`.

**Stage H record.** After the corresponding author's written go-ahead, recorded with its
UTC time, the manifest and the job logs are committed,
tagged `c1-covariates-posthoc-stageh-<date>` and published as a GitHub Release, before any
local stage. The record holds file hashes and logs only.

**Stage L1 (local; one run, after the Stage H record).** It refuses to start unless the
Stage H manifest has the sha256 given on the command line (the value in the Stage H
record). Every discovery analysis and set of section 3, and every cell on the published
matrices. Outputs under `server_export/results/c1_covariates/L1/`: `sets_<cohort>.csv`,
`cells_L1.csv`, `result_L1.json` (with the sha256 of the other outputs), `L1.log`.

**Stage L2 (local; one run).** It refuses to start unless every Stage H file it reads has
the sha256 the manifest gives, the manifest has the sha256 given on the command line, and
the Stage L1 outputs have the hashes in `result_L1.json`; these checks are made before the
stage's folder is created. It computes every cell on the Stage H matrices, the change D,
the capacity randomisation and the diagnostic of section 5, and writes `cells_L2.csv`,
`classes.json` (the classes, types and case of section 6, the counts of M4 and the reverse
flags, with the unrounded values they rest on, and the markers used) and `result_L2.json`.
A result file records the sha256 of every output of its stage except the stage log, which
is hashed in the results record.

**Interrupted runs.** A run stopped before writing its result is kept with its log and
resubmitted unchanged (section 6.1); each attempt is recorded. An attempt counts from the
start of its computation; a refusal at the pre-flight checks is not an attempt. A completed
run is never repeated.

**Defects.** A defect found before a run is corrected by a dated amendment that contains
the diff of the script, gives the new script hash and a new sidecar, states which stages
had been run when the defect was found, and is released publicly, as an annotated tag and a
GitHub Release after the corresponding author's written go-ahead, before the run; it
changes no definition, reading or sentence. An exception during a run is reported before
any resubmission: if it is an environment failure (disk, memory, file system), the run is
resubmitted unchanged and the exception is recorded; otherwise it is treated as a defect
found before that run, and the amendment also reproduces what the stopped run had written
to its log and names every run already completed, whose outputs stand. If a correction
needs a choice that this file does not make, or a sentence of section 7 cannot be filled as
written, the amendment states the choice, is released before the affected run or before the
manuscript is edited, records which results were then known and, where two choices are
defensible, takes the one less favourable to the paper. If a pinned host is withdrawn from
the cluster after freezing, H-V1 is repeated on the next host of its list under such an
amendment; nothing else about the run changes. A defect found after a run is disclosed in Methods with
both results; section 7 is filled from the run made with the frozen code; the corrected run
is reported beside it if the defect made the frozen run depart from this file's text; and
if the two runs give different classes, the type, the case, the class-U sentences of 7.7
and 7.8, the reverse sentences of 7.8 and the counts of M4 are those of whichever run is
less favourable to the paper (type C before R before none; class A before U before H; a
reverse or forward cell that outranks before one that does not; a matched cell that does
not outrank before one that does), and this sentence is appended
to M5: "A script defect found after the run was corrected; the corrected run changed the
[C3|R20] class in [uterus|lung] from [..] to [..] (Methods)."

**Results record.** After the corresponding author's written go-ahead, recorded with its
UTC time, the outputs of both local stages, their logs
and their sha256 are committed, tagged `c1-covariates-posthoc-results-<date>` and published
as a GitHub Release. The sentences of section 7 are then filled in exactly as written, by
transcription from `classes.json` and the two cell tables.

**Gate record** (2026-11-20/21 UTC; every log named here is in `review/C1_COVARIATES_2026-11-21/gate_logs/`
and is listed in the sidecar; times are CST, UTC+8, as printed by the cluster).

*Scripts frozen.* `c1_covariates_local.py` sha256 `1ea40fa50759366d09f61360d6d2f80dce9703d6d95f6a8ce43ff0f0ae61e779`;
`c1_covariates_cluster.py` `da89be848e8dc6b64e8e9fb92fd6bc1a387b94aa3667374e6f1a8284331501dc`;
`lsf_c1_cov_ucec.sh` `23c9df61d31db78c9f6c6ef1039d056e5a2d8e715782605ba5dc0dcc85d320e0`;
`lsf_c1_cov_luad.sh` `bbb94eb7ff141533567a20c488e059704032b2219e64dab01398dfca0df06d39`;
`lsf_c1_cov_discovery.sh` `933a9c962cfbba73dee1bd8273fb00dd2b46412b1351bbad6fd9ce49333c103b`.
Every gate JSON of the cluster carries the cluster script's hash `da89be84...`.

*Self-tests* (`selftest_local.log`, `selftest_cluster.log`; the first line of each is the hash of the script
tested): local, SELFTEST PASSED in 27 s; cluster, 113 checks, 0 failed.

*V1, V1b, V2* (`validate_only.log`, `validate_only.json`, written by the local script above; exit 0, VALIDATED).
V1: the published per-gene increment is reproduced to 9.8e-15 (UCEC) and 2.1e-14 (LUAD), own-mRNA R2 to 4.4e-16,
p-values and FDR exactly, and the published sets (2,566; 2,310) are identical. V1b: R2 of own mRNA plus the 20
components to 1.6e-13 in both organs; increment of morphology over that baseline to 1.5e-13 (UCEC) and 4.0e-6
(LUAD; 2,304 of 2,310 genes above 1e-9, inside the 1e-5 tolerance of section 9 for the reason given there). V2: AUC,
null mean, null SD, p and rho of the four saved arms reproduced with worst difference 0.0. Operator indicators:
12 covering 73 of 100 (UCEC) and 75 of 105 (LUAD). Batch-corrected sets: 201 genes, 198 in the confirmatory
universe (UCEC); 726 genes, 721 in the universe, file hash equal to the registered one (LUAD).

*Blind runs* (counts only): UCEC job 76713695 on s006: 138 cases, K 28, genes tested 10,153, usable 10,146
(identical to the saved matrix), 2,536 set genes among them; 59,427 RNA columns, 2,396 all-zero; 20 components
explain 0.5844 of the variance; markers found 24/24, 25/25, 24/24, none with zero SD; C3 testable; 169 of 169
slide embeddings present. LUAD job 76712907 on s002: 112 cases, K 22, genes tested 10,417, usable 10,406
(identical), 2,287 set genes; 3,309 all-zero RNA columns; explained variance 0.4746; markers 24/24, 25/25, 24/24;
C3 testable; 135 of 136 slide embeddings present (the missing one, `C3N-02142-23`, is the case lost to IDC
coverage in the registered rule). The result-folder listings (`c1_cov_*_blind_listing_*.tsv`; 447 rows (UCEC)
and 440 (LUAD), not truncated; they cover the discovery trees, both confirmatory trees, the folders of the earlier post hoc and
measurement-error analyses and this analysis's own files) were read for section 0, item 7. In the two
confirmatory trees every file belongs to an analysis named in section 0: the C1 tests, their bootstraps and
sensitivities, the per-gene rediscoveries, the D4 stratified and completeness re-reads and the two plex re-reads.
The files whose names mention composition, stroma or operator are the discovery-tree composition and stroma
controls and operator-axis tables of 2026-08-31 to 2026-09-06, which section 0, items 2 and 5 disclose, and the
published UCEC completeness re-read. No file of any other kind of analysis of a confirmatory cohort was found.

*H-V1.* UCEC, host s001 (the host of the job that wrote the saved matrix; first of its list), job 76713696:
strict PASSED, largest difference 0 on columns 0 to 5, AUC of the published set from the recomputed column 0
0.5974425564900907 equal to the published value. LUAD, host s002 (host of the primary run; first of its list),
job 76712908: strict PASSED, largest difference 0, AUC 0.6439668501223573 equal to the published value. No other
host was tried or needed. **Pinned hosts: UCEC s001, LUAD s002** (`PINNED_HOSTS.txt`).

*H-V2.* UCEC on s006 (job 76713697) and LUAD on s002 (job 76712905), the hosts that produced the inputs files:
PASSED, largest difference after sign alignment 0 (all signs +1).

*H0 and V3.* Discovery score tables exported on the same hosts, UCEC job 76713698 (sha256 `865a8a35984c93ea6d6e3b8056b49ed45bff5432f60bb6052e01b1f3e62a7b2a`)
and LUAD job 76712906 (`95c986a0e302e16fd3321dcd39958710ffa1ce0fc76cbba18a8d0c568fdfea58`); the hashes were
re-computed locally after the transfer and agree. V3 (`gate_v3.log`): all 2,566 (UCEC) and all 2,310 (LUAD)
published composition-control genes reproduced, worst difference 6.7e-16 and 4.4e-16; VALIDATED.
`--check-stage-h` on the two tables (`check_stage_h_pre.log`): passed.

*Failed attempts, kept.* Eight first attempts (UCEC jobs 76712075, 76712076, 76712077, 76712079; LUAD jobs
76712078, 76712080, 76712082, 76712083; 2026-11-21 06:28-06:44 CST) stopped at the import of `sklearn` with
`ImportError: numpy.core.multiarray failed to import`, before any data were read: the conda base now holds numpy
2.0.2 (installed 2026-03-19) beside scipy 1.11.4 (installed 2026-09-27), which do not import together, whereas the
October 2026 runs used numpy 1.26.4. Their logs and the folder logs of those attempts are kept. The remedy,
introduced before any gate had run, is a private copy of numpy 1.26.4 (`numpy_1.26.4_tree.tar.gz`, sha256
`92f2aeb5603a18ed35b2242db5fa07f6ac43a96104e48a2b0ffa76f3663ef8ba`, unpacked by each wrapper into the node-local
temporary folder and checked to import as 1.26.4); the conda base and `~/.local` were not changed. A first attempt
to install it inside `/public` stopped because the account had reached its file-count quota (2,000,000 files); the
half-installed folder was removed and the conda package cache cleaned by the corresponding author, and about
7,400 file names were free afterwards. All gates above were run with the remedy in place (run JSONs: numpy 1.26.4).

*Closed parts.* None. `CLOSED_PARTS.txt` holds no entry, so every closed form of section 7 and the sixteen Abstract
variants with a closed reading are void. The Abstract variants that remain possible are the eight of
`ABSTRACT_VARIANTS.md` without a closed reading: 1-u[-]-l[C], 1-u[-]-l[R], 1-u[C]-l[-], 1-u[R]-l[-], 2-u[C]-l[C],
2-u[C]-l[R], 2-u[R]-l[C], 2-u[R]-l[R], besides the unchanged Abstract of Case 0.

*What these gates do not show.* No new increment of any kind was computed on a confirmatory cohort, and no
discovery analysis of section 3 was run, before this record was written.
