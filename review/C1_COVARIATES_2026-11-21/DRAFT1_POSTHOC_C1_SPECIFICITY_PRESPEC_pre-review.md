# Post hoc C1 re-reads for specificity to histology: tissue-covariate baselines and non-histology selectors (C1-UCEC and C1-LUAD): pre-specification (file date 2026-11-21, local time UTC+8)

**Dates and times.** The date in this file's name is the local date (UTC+8) of the drafting
session. Every time stated inside this file is UTC. The sidecar's `# written` line records
when the file was frozen; the manuscript cites that UTC date.

**Status.** Written on 2026-11-20 UTC, after the published discovery results, the published
C1-UCEC reading (2026-09-16), the registered C1-LUAD reading (2026-10-02), the post hoc
analyses pre-specified on 2026-10-02 and 2026-10-08, and the Abstract rewrite recorded in
`C1_LUAD_ADDENDA.md` Entry 12 (2026-11-20). An internal review of 2026-11-20 asked whether
the confirmatory set-level ordering needs histology, or whether any reading of tissue
state would give it. The corresponding author asked for that question to be tested
(2026-11-20). Frozen by the sha256 of this file, the local script and every local input,
listed in the sidecar `POSTHOC_C1_SPECIFICITY_PRESPEC_2026-11-21.md.sha256`; released
publicly on GitHub, as an annotated tag and a Release, before any run of section 9.

## 0. Question, scope and prior exposure

**Question.** Both confirmatory tests rank genes by how much histology adds to the gene's
own mRNA in new patients, and ask whether the genes selected in discovery rank above the
rest. Their null breaks the pairing between a patient's slide and that patient's omics. It
cannot tell an ordering that needs histology from one that any measure of tissue state
would give, including one that the transcriptome already holds. This analysis asks three
things of each confirmatory cohort:

1. **Conditional (primary).** Does the published selection set still outrank background
   when the confirmatory baseline already holds transcriptomic covariates of tissue state?
2. **Reverse.** Does the published selection set also outrank background when those
   covariates take the place of histology in the confirmatory cohort?
3. **Forward.** Do sets selected in discovery without histology outrank background under
   confirmatory histology?

**Scope.**

- **Post hoc, reported whatever the result.** Nothing below was pre-registered. After
  freezing, no covariate, set, arm, cell, background, threshold, reading or sentence is
  added, dropped or changed, except by a dated amendment released before the run it affects
  (section 9).
- **No effect on the registered or published readings.** The published C1-UCEC reading
  (AUC 0.597; null 0.498 +/- 0.036; p = 0.004) and the registered C1-LUAD reading (AUC
  0.644; null 0.524 +/- 0.019; p1 = 0.001; p_strat = 0.001; p2 = 0.001) stand as reported.
  Nothing here is pooled with either, and no wording of
  `C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md` section 5.2 changes.
- **No new cohort.** Section 5.8 of the C1-LUAD rule ("No further CPTAC-3 confirmatory
  cohort is attempted before submission") is not touched: these are re-reads of the two
  cohorts already read.
- **Departures from the two rules.** (i) The primary and reverse arms change the
  confirmatory estimator (the baseline, or the block), which section 3 of both rules fixed.
  (ii) The forward and matched cells use discovery selection sets other than the published
  ones, which section 2 of both rules fixed as not re-derived.
- **Kept from the C1-LUAD rule.** The word "replicat\*" is not used for any result of this
  analysis (rule 5.7 item 2), and every sentence that states a result names the set (rule
  5.7 item 6).
- **Unrestricted null only.** Every new confirmatory arm is read against the unrestricted
  patient<->slide permutation null of the published tests. No operator-stratified or
  plex-stratified null is run for a new arm. Because the unrestricted null is the easiest of
  the three to pass, an ordering not shown under it would not be shown under the others; an
  ordering that holds under it is not thereby robust to operator or plex (section 8).

**Prior exposure (disclosed).**

1. *Discovery, composition.* The published composition control (main text,
   `sec:stroma`) added the three scores of section 2 to the discovery baseline: the median
   increment of the published sets retained 36% in UCEC and 6% in LUAD (68% CCRCC, 34% PDAC,
   3% GBM). These numbers were known to everyone involved and are why LUAD is expected to be
   the harder cohort here.
2. *Discovery, expression components.* The published transcriptome-wide baseline
   (`sec:txome`) added the 20 components of section 2 to the discovery baseline: a median
   40% of the noise-adjusted advantage was retained (77% in UCEC, the highest), and 227
   (UCEC) and 55 (LUAD) of the published sets passed FDR < 0.05 under that baseline, with
   Benjamini-Hochberg taken over the published set only and a different permutation
   sequence. The sets S_M|R20 of section 3 are related to these counts but are not the same.
3. *Discovery, operator.* `incremental_r2_batchonly` in
   `server_export/results/<cohort>__sitepack_operator.csv` is the per-gene increment of
   operator dummies in place of histology, without a p-value; its summary is in the
   supplement (Supplementary Table posthoc, "B.o. vs n.m.").
4. *Confirmatory, other sets.* On the C1-UCEC matrix, unrestricted AUCs were computed on
   2026-10-07 for discovery top-k sets and for the selected sets of the other four discovery
   organs (0.564-0.585), as disclosed in `POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md`
   section 0, and that pre-specification's sets were then read. The per-gene confirmatory
   increments have been in the published gene-level tables since 2026-09-16 (UCEC) and
   2026-10-02 (LUAD).
5. *Blind.* To the knowledge of the corresponding author and of the session that wrote
   this file, and after a search of the project's scripts, logs and result folders: no
   statistic that uses the composition scores, the expression components or an
   operator-selected, score-selected or component-selected gene set has been computed on
   either confirmatory cohort by anyone; no discovery analysis has used the scores, the
   components or operator dummies as the block with a permutation p-value; and the sets
   S_M|C3 and S_M|R20 of section 3 have not been computed.
6. *The drafting session (2026-11-20).* It read array shapes, gene lists and label counts of
   the files of section 1, and ran `--selftest` (synthetic data) and `--validate-only`
   (published quantities; section 9). It computed nothing else on project data.

## 1. Inputs

Paths are relative to `MorphoResidual_paper/`; the sha256 of each is in the sidecar.

- **Discovery inputs.** `review/MEAS_ERROR_2026-10-08/export/ucec/inputs_ucec.npz`
  (`e28b86fa...`; 100 patients, 10,512 fit genes) and `.../luad/inputs_luad.npz`
  (`b64cdf3b...`; 105 patients, 10,724 fit genes), exported on the cluster on 2026-10-07 for
  the measurement-error analysis: own-gene mRNA (log2(TPM+1)), protein, the 20 published
  morphology components, the 20 expression components (`rna_pcs20`), operator labels and
  the published per-gene results. They reproduce the published per-gene increment exactly
  (section 9, V1).
- **Confirmatory matrices.** `server_export/scripts/pinned/c1_ucec_d4_stratified.npz`
  (`ef75f1c0...`): 10,146 usable genes x 1,001 columns (column 0 the observed confirmatory
  increment, columns 1-1000 the permutations), arms `M_unrestricted`, `M_within_operator`
  and `M_within_operator_residualized`. `server_export/scripts/pinned/luad_c1/c1_luad_primary_primary.npz`
  (`2034e711...`): 10,406 usable genes x 1,001 columns, the registered H1 arm (unrestricted
  null).
- **Published readings, for validation.** `server_export/scripts/pinned/c1_ucec_d4_followup.json`,
  `server_export/scripts/pinned/luad_c1/c1_luad_primary_result_primary.json`.
- **Published composition controls, for validation of the transferred scores.**
  `server_export/results/results__composition_controls_ucec.csv` and `..._luad.csv`.
- **Definitions, by reference.** `server_export/scripts/composition_controls.py`
  (`46b4a1d8...`: marker lists and score), `server_export/scripts/transcriptome_baseline.py`
  (`761df7d7...`: expression components), `server_export/scripts/residual_analysis.py`
  (`a54a0a56...`: estimator), `split_replication_power.py` (`79530c07...`),
  `server_export/scripts/c1_run_test.py` (`fc578c70...`) and
  `server_export/scripts/c1_run_test_luad.py` (`1d5ee05d...`: the confirmatory tests),
  `server_export/scripts/residual_analysis_sitepack.py` (operator dummies),
  `review/MEAS_ERROR_2026-10-08/theta_kernel.py` (`6c95cd11...`: the reference kernel the
  local script imports).
- **On the cluster only** (`/public/home/fjhui/ZW/`): the RNA-seq, protein and slide
  embeddings of both confirmatory cohorts and the full discovery RNA-seq matrices. Stage H
  (section 9) reads them there.

## 2. Covariates

Three covariate sets. Each is computed inside one cohort, on that cohort's analysed cases
only (discovery: the 100 and 105 patients of the inputs files; confirmatory: the 138 and 112
analysed cases of the published tests).

- **C3, three transcriptomic composition scores** (stroma, immune, proliferation, in this
  order), exactly as the published control: `composition_controls.comp_score` with
  `composition_controls.MARKERS`, unchanged. For each family, the markers present as columns
  of the cohort's RNA-seq matrix (`residual_analysis.load_rna_matrix`: all genes,
  log2(TPM+1), tumour files averaged per case) are z-scored over the analysed cases and
  averaged; a family with fewer than 6 markers present makes C3 untestable in that cohort.
  The target gene is not removed from a score it belongs to, as in the published control.
  The number of markers used per family is reported.
- **R20, 20 expression components**, exactly as the published transcriptome-wide
  baseline (`transcriptome_baseline.py`): every column of the cohort's RNA-seq matrix on the
  analysed cases is standardised (mean 0; SD with ddof 0; a zero SD is set to 1), and
  `PCA(n_components=20, svd_solver="full").fit_transform` gives the 20 columns. In discovery
  these are the `rna_pcs20` arrays of the inputs files. The number stays 20 in the
  confirmatory cohorts; it is not scaled with n as the morphology components are.
- **OP, operator dummies (discovery only)**, as `residual_analysis_sitepack.py`: one dummy
  for each of the 12 largest scanning-operator levels with at least 3 cases; every other
  case, including the one unlabelled case in each cohort, has an all-zero row. UCEC: 12
  dummies covering 73 of 100 cases. LUAD: 12 dummies covering 75 of 105 cases. (Counted from
  the labels alone.)
- **Capacity placebo R20cap (confirmatory only).** The R20 matrix with its rows reassigned
  to patients by one fixed permutation, `RandomState(crc32(b"c1_spec|ucec|cap")).permutation(138)`
  and `RandomState(crc32(b"c1_spec|luad|cap")).permutation(112)`. It keeps the size and
  covariance of R20 and carries no patient information. It is one draw.

Transcriptomic scores are the only composition measure available in these four cohorts; no
pathologist-read composition exists for them (Limitations (iv)).

## 3. Discovery sets

Every set uses the published discovery estimator and changes exactly one thing, the block
or the baseline: per gene, on the patients with a protein value, 5-fold
`KFold(shuffle=True, random_state=0)`, ridge (alpha 1) on columns standardised within the
training fold; increment = R2(baseline + block) - R2(baseline); null = the block's rows
permuted over all patients with the published sequence (`RandomState(0)`, 1,000 sequential
draws), the baseline staying with its patients; one-sided p = (1 + #{null >= observed}) /
1001; Benjamini-Hochberg over all fit genes of the cohort (10,512; 10,724); selected = FDR <
0.05 and increment > 0.

| Set | Baseline | Block | Stage |
|---|---|---|---|
| S_M (published; 2,566 UCEC, 2,310 LUAD) | own mRNA | 20 morphology components | published |
| S_C3 | own mRNA | C3 | L2 |
| S_R20 | own mRNA | R20 | L1 |
| S_OP | own mRNA | OP | L1 |
| S_M\|C3 | own mRNA + C3 | 20 morphology components | L2 |
| S_M\|R20 | own mRNA + R20 | 20 morphology components | L1 |

S_C3, S_R20 and S_OP are selected without histology. S_M|C3 and S_M|R20 are selected by
histology after the covariates have entered the discovery baseline. For each set the script
reports the number selected, the number of Benjamini-Hochberg rejections before the sign
filter, the number of genes at the permutation floor p = 1/1001, and k\* = ceil(m / (0.05 x
1001)): with 1,000 permutations Benjamini-Hochberg returns either no gene or at least k\*
genes, so an empty set is "empty at B = 1,000, limited by the permutation floor; not
evidence of absence".

## 4. Confirmatory arms (Stage H, on the cluster)

For each cohort the arm keeps everything of the published test: the analysed cases, the
gene universe and MIN_N = 30, the K morphology components (K = 28 for UCEC and 22 for LUAD,
fitted on the confirmatory cohort's own mean-pooled embeddings), `residual_analysis.cv_r2`,
and the published permutations, drawn in the same way from
`RandomState(crc32(b"c1_ucec|perm"))` and `RandomState(crc32(b"c1_luad|perm"))`, so that
permutation b is the same reassignment of rows in every arm and in the published matrix.
Each arm gives a genes x 1,001 matrix on the published usable universe (10,146; 10,406).

| Arm | Baseline | Block (rows permuted under the null) |
|---|---|---|
| A_C3 | own mRNA + C3 | K morphology components |
| A_R20 | own mRNA + R20 | K morphology components |
| A_R20cap | own mRNA + R20cap | K morphology components |
| V_C3 | own mRNA | C3 |
| V_R20 | own mRNA | R20 |

In the A arms the covariates stay with their patients in every permutation. In the V arms
no slide is used: the covariates are the block and their rows are permuted. Stage H also
exports the discovery C3 scores of both cohorts (three numbers per patient), which the
discovery sets S_C3 and S_M|C3 need.

## 5. Cells: statistic, background and null

A cell is one set read on one matrix. The statistic is the published one: the AUC of the
matrix's column-0 increments of the set against the background (rank Mann-Whitney) by
`assemble`, copied unchanged from `c1_run_test.py`; the null is the same statistic on
columns 1-1000; one-sided p = (1 + #{null >= observed}) / 1001; margin = AUC minus the null
mean; z = margin / null SD. The null mean is never assumed to be 0.5. Spearman rho between
the set's own discovery increment and the confirmatory increment is reported and not read.

- **Universe and background.** The usable confirmatory genes of the published test. The
  background of a set is every other gene of the universe (the C1-LUAD H2 rule).
- **Minimum size.** A set, or a background, with fewer than 50 genes in the universe is
  "untestable".

| Family | Cell (set on matrix), in each cohort | Role |
|---|---|---|
| Primary | S_M on A_C3; S_M on A_R20 | read, section 6.1 |
| Capacity | S_M on A_R20cap | enters the R20 reading only |
| Reverse | S_M on V_C3; S_M on V_R20 | read, section 6.2 |
| Forward | S_C3, S_R20, S_OP on the published matrix | read, section 6.2 |
| Matched | S_M\|C3 on A_C3; S_M\|R20 on A_R20 | read, section 6.2 |
| Covariate-to-covariate | S_C3 on V_C3; S_R20 on V_R20 | reported, not read |
| Partition | for Z in C3, R20, OP, on the published matrix: the genes in both S_M and S_Z, in S_M only, and in S_Z only, each against the genes in neither | reported, not read |
| Reference | S_M on the published matrix | reproduces the published reading |

For C1-UCEC the forward, partition and reference cells are also computed on the two other
saved arms (within-operator null; operator-residualised components); these are reported and
not read. For C1-LUAD only the unrestricted arm is saved.

**What the margins are not.** The sets differ in size and in how they were selected, and
every discovery selection carries its own selection effect. Margins are compared across
cells descriptively only. No cell is called an effect size, and no sentence says that one
selector "beats" another.

## 6. Reading rules

p is used unrounded. "Outranked" means p < 0.05; "not shown" means p >= 0.05. The table
also prints the tier of the registered rule (p < 0.05; 0.05 <= p < 0.15; p >= 0.15).

### 6.1 Primary readings (four, one per cohort and covariate set; no multiplicity correction)

- **C3 class** of a cohort: **H** if S_M on A_C3 has p < 0.05; **N** otherwise.
- **R20 class** of a cohort, from p_R (S_M on A_R20) and p_cap (S_M on A_R20cap):
  **H** if p_R < 0.05; **A** (absorbed) if p_R >= 0.05 and p_cap < 0.05; **U**
  (uninformative) if p_R >= 0.05 and p_cap >= 0.05. Three covariates cost almost no
  capacity, so C3 needs no placebo; 20 covariates do (43 predictors at about 90 training
  patients in LUAD, 49 at about 110 in UCEC), so a failure under R20 is read only against
  the placebo.
- **Flag F** of a cohort: true if its C3 class is N or its R20 class is A. F is what sections
  7.4 to 7.6 use.
- **Not completed, untestable.** A primary cell that is untestable, or not completed when
  the manuscript is submitted, is reported as such and counts as "not shown": C3 class N.
  If A_R20 is not completed the R20 class is "not completed" and F rests on C3 alone; if
  only A_R20cap is not completed, a p_R >= 0.05 is class A.
- **No correction, and why.** The four readings are four questions and are never summarised
  as one. Each is printed with its p. The only consequence that reaches the Abstract follows
  from a failure (section 7.4), for which leaving p uncorrected is the stricter choice.

### 6.2 Secondary readings (each cell: "outranked" or "not shown"; counts are reported)

- Reverse: four cells. Forward: six cells. Matched: four cells.
- No correction is applied and none of these cells changes the Abstract, the Discussion or
  the Limitations.

### 6.3 What a class means, fixed now

- C3 class H: the published set's ordering survives these three transcriptomic scores. It
  does not show independence from composition measured in other ways.
- C3 class N, R20 class A: the ordering is not shown beyond the covariates. It does not show
  that the covariates explain it.
- R20 class U: nothing is learned about the 20 components at this cohort size.
- A reverse cell that outranks: the set's ordering can be obtained in the confirmatory
  cohort without a slide.
- A forward cell that outranks: genes chosen without histology also receive larger
  histology increments in new patients.
- A matched cell that outranks: a selection made by histology after the covariates are in
  the baseline also ranks above background when the confirmatory baseline holds them.

## 7. Manuscript wording, fixed now (filled from the outputs after Stage L2, never before)

[a], [m], [s], [p] are the AUC, null mean, null SD and one-sided p of the named cell, printed
as in the confirmatory paragraphs (three decimals; p as computed, 0.001 at the floor). [date]
is the UTC date of the sidecar's `# written` line. A cell that is untestable or not completed
is named so in place of its numbers. Stage L1 numbers are not written into the manuscript
before Stage L2 completes or is closed as not completed.

### 7.1 Results: one new paragraph at the end of `\S\ref{sec:c1ucecres}` (after "...so it bounds protein batch rather than isolating it.")

The paragraph is these sentences, in this order and in no other.

**M1 (C3), exactly one of:**

- both cohorts H: "Post hoc, in an analysis specified on [date] in a dated file released
  publicly before it ran (Supplementary Note~\ref{snote:c1spec}), we added transcriptomic
  stromal, immune and proliferation scores to the confirmatory baseline: both uncorrected
  discovery sets (2{,}566; 2{,}310 genes) still outranked background (uterus: AUC [a]
  against [m] $\pm$ [s], one-sided $p=$[p]; lung: [a] against [m] $\pm$ [s], $p=$[p])."
- one H (X), one N (Y): "Post hoc, in an analysis specified on [date] in a dated file
  released publicly before it ran (Supplementary Note~\ref{snote:c1spec}), we added
  transcriptomic stromal, immune and proliferation scores to the confirmatory baseline: the
  uncorrected discovery set still outranked background in [X] ([n] genes; AUC [a] against
  [m] $\pm$ [s], one-sided $p=$[p]) but was not shown to in [Y] ([n] genes; [a] against [m]
  $\pm$ [s], $p=$[p])."
- both N: "Post hoc, in an analysis specified on [date] in a dated file released publicly
  before it ran (Supplementary Note~\ref{snote:c1spec}), we added transcriptomic stromal,
  immune and proliferation scores to the confirmatory baseline: neither uncorrected
  discovery set (2{,}566; 2{,}310 genes) was shown to outrank background (uterus: AUC [a]
  against [m] $\pm$ [s], one-sided $p=$[p]; lung: [a] against [m] $\pm$ [s], $p=$[p])."

**M2 (R20):** "With 20 expression components added instead, one-sided $p$ was [p] in
uterus and [p] in lung[tail]." The tail is, for each cohort in the order uterus, lung:
nothing if its R20 class is H; "; a placebo baseline of the same size left the
[uterine|lung] ordering in place ($p=$[p])" if A; "; a placebo baseline of the same size
also removed the [uterine|lung] ordering ($p=$[p]), so that result is uninformative" if U.

**M3 (reverse, C3), exactly one of:**

- both outrank: "In place of histology, the three scores alone also ordered both sets
  above background ($p=$[p]; [p])."
- one (X) outranks: "In place of histology, the three scores alone also ordered the
  [uterine|lung] set above background ($p=$[p]) but not the [lung|uterine] set ($p=$[p])."
- neither: "In place of histology, the three scores alone did not order either set above
  background ($p=$[p]; [p])."

**M4 (the fixed limits; verbatim in every outcome):** "These covariates are
transcriptomic estimates and the null is unrestricted: an ordering that holds is not
thereby independent of composition measured in other ways or robust to operator, and one
not shown is not thereby explained by the covariates (Supplementary
Table~\ref{tab:c1spec})."

If the main text then exceeds the journal's length, free text elsewhere is shortened; none
of M1 to M4 is.

### 7.2 Methods: one paragraph appended to `\S\ref{sec:c1luadmeth}`

"\emph{Specificity to histology (post hoc).} This analysis was specified on [date] in a
dated file released publicly before it ran (Supplementary Note~\ref{snote:c1spec}). In each
confirmatory cohort the published test was repeated with the baseline enlarged by three
transcriptomic composition scores (stromal, immune, proliferation; \S\ref{sec:confmeth}) or
by 20 expression components (\S\ref{sec:txomemeth}), each computed on the confirmatory
cases; with the 20 components reassigned to patients at random as a capacity placebo; and
with the scores or the components in place of the morphology components. Cases, genes,
morphology components, folds and the 1{,}000 permutations were those of the published
tests, and only the unrestricted null was run. In discovery, sets were also selected with
the scores, the components or scanning-operator indicators in place of morphology, and with
morphology after the scores or the components had entered the baseline, by the published
estimator with Benjamini--Hochberg control over all fitted genes."

### 7.3 Supplement: new Supplementary Note `snote:c1spec` and Supplementary Table `tab:c1spec`

The table has one row per cell of section 5, in every outcome: cohort, family, set, matrix,
null, genes in the set and in the background, AUC, null mean, null SD, margin, z, one-sided
p, tier, rho, and the class of section 6.1 where it applies. A second table gives every
discovery set of section 3 with the numbers listed there, the overlap of each with S_M, and
the markers used per composition family in each of the four cohort-by-wave computations.
The note states sections 0, 5 ("What the margins are not") and 8 of this file in substance,
and contains these sentences verbatim:

- Forward: "Of [t] discovery sets selected without histology that were large enough to
  test (by the three scores, by the 20 components or by scanning operator, in each organ),
  [k] outranked background under confirmatory histology[: names]."
- Matched: "Of [t] discovery sets selected by histology after the covariates had entered
  the discovery baseline that were large enough to test, [k] outranked background when the
  confirmatory baseline held the same covariates[: names]." If t = 0: "No discovery set
  selected by histology after the covariates had entered the discovery baseline was large
  enough to test."
- Reverse, R20: "In place of histology, the 20 components alone gave one-sided $p=$[p] in
  uterus and [p] in lung."
- Placebo: "The capacity placebo is a single random reassignment of the 20 components to
  patients."

### 7.4 Abstract

**Abstract policy** (that of the post hoc pre-specifications of file date 2026-10-08,
applied here as `C1_LUAD_ADDENDA.md` Entry 12 reads it). No favourable post hoc result
enters the Abstract. The only permitted Abstract edits are qualifiers, fixed in advance,
that stop an existing sentence from overstating in the light of a post hoc result. If the
Abstract is later shortened, a qualifier added under this policy is kept in substance.

Here: the clause below is the only Abstract edit of this analysis. It reports an outcome
unfavourable to the paper's reading and qualifies the sentence on the two later cohorts,
which a reader would otherwise take to show an ordering that histology adds to what the
transcriptome records of tissue state.

- **F is false in both cohorts** (every other outcome of the secondary and descriptive
  cells included): the Abstract is unchanged. Reason, fixed now: no favourable post hoc
  result enters the Abstract.
- **F is true in both cohorts:** the sentence on the two later cohorts gains, immediately
  before its final full stop, "; post hoc, neither ordering was shown beyond transcriptomic
  tissue covariates".
- **F is true in one cohort:** it gains "; post hoc, the [uterine|lung] ordering was not
  shown beyond transcriptomic tissue covariates".

The Abstract stands at 141 to 143 words against the journal's limit of 150, and the clause
has 10 or 12. Room is made as in Entries 11 and 12: by a dated entry of
`C1_LUAD_ADDENDA.md`, released before the edit, that cuts no item protected by C1-LUAD rule
5.7 or by Entry 11 and leaves the clause verbatim. Read literally, rule 5.7 item 3 ("No
other descriptive item enters the Abstract") is departed from by this clause, as it was by
the qualifiers of the 2026-10-08 pre-specifications; the departure weakens the paper's
claim and the same entry records it.

### 7.5 Discussion

- **F true in at least one cohort:** immediately after the sentence that ends "not the
  per-gene increment, which is negative for most genes there." this sentence is inserted:
  "Post hoc, [neither ordering was | the uterine ordering was not | the lung ordering was
  not] shown beyond transcriptomic tissue covariates (\S\ref{sec:c1ucecres})."
- **F false in both:** the Discussion is unchanged.

### 7.6 Limitations (iv): one sentence appended, exactly one of

- F true in at least one cohort: "Post hoc, [neither confirmatory ordering was | the
  uterine confirmatory ordering was not | the lung confirmatory ordering was not] shown with
  transcriptomic tissue covariates in the baseline (\S\ref{sec:c1ucecres})."
- F false in both: "Post hoc, both confirmatory orderings held with these scores in the
  baseline, which does not make them independent of composition measured in other ways
  (\S\ref{sec:c1ucecres})."

### 7.7 Figure, provenance table and what this rule does not fix

- Fig. 6a gains one row for each of the four primary cells in every outcome, tagged post
  hoc; its legend names them.
- The supplementary provenance table gains a row for this file, for each amendment and for
  the results record.
- **Minimum edits.** Sections 7.1 to 7.6 are made in every case and are not weakened.
  Whether the title, the opening sentence or the target journal change after an
  unfavourable reading is the authors' decision, taken after the results and recorded in a
  dated entry of `C1_LUAD_ADDENDA.md`; such a change may weaken the paper's claims and may
  not strengthen them.

## 8. What this analysis cannot establish

- **Independence from composition.** The scores and components come from the same RNA-seq
  as the baseline and track pathologist-read composition weakly in the one cohort where
  both exist (PDAC). A class H says the ordering survives these estimates.
- **That a covariate explains an ordering.** A class N or A can come from lost power. The
  test has no direct measure of how much of a margin a covariate removes; margins before
  and after are printed, not tested against each other.
- **Robustness to operator or plex** of any new arm: only the unrestricted null is run.
- **Anything at the gene level.** Most confirmatory increments are negative; every reading
  is about the ordering of a set against its background.
- **Which way composition and histology run.** If histology reads composition and
  composition shapes the protein-to-transcript ratio, a reverse or forward cell that
  outranks is what both the paper's reading and the critique predict.
- **The placebo's variability.** R20cap is one draw; class A against class U can turn on it.
- **Comparisons between selectors.** Set sizes, block sizes and selection effects differ
  (section 5).
- **Anything about CCRCC, GBM or PDAC**, which have no confirmatory cohort here.

## 9. Freeze, stages and run order

**Local script.** `review/C1_SPECIFICITY_2026-11-21/c1_specificity_local.py`. Stage runs
refuse to start unless every sidecar hash matches, and refuse to overwrite a completed
stage.

**Before freezing (published quantities and synthetic data only).**

- `--selftest`: on synthetic data the batched kernel equals the reference kernel
  (`theta_kernel.oof_ridge`, which is `residual_analysis.cv_r2`), `bh_fdr` equals a direct
  implementation and `assemble` equals scipy's Mann-Whitney statistic.
- `--validate-only`, V1: the code path used for every new discovery set, given the
  published block, reproduces the published per-gene increment, p-value and FDR of every fit
  gene and the published selection sets (2,566; 2,310). V2: `assemble` on the saved matrices
  reproduces the published readings. The log is `validate_only.log` beside the script; its
  result is recorded at the end of this section.

**Release.** This file, its sidecar, the script, `validate_only.log` and the sidecar list
are committed, tagged `c1-posthoc-specificity-prespec-2026-11-21` and published as a GitHub
Release before Stage L1. Withheld, as derived data fixed by their sha256 in the sidecar and
deposited with the other derived data at acceptance: the two inputs files, the two
confirmatory matrices and the two composition-control tables. The commit message lists
them.

**Stage L1 (local; one run).** `--stage L1`: sets S_R20, S_OP and S_M|R20 in both
cohorts; the reference, forward and partition cells for R20 and OP.

**Stage H (cluster).** The cluster scripts do not exist when this file is frozen, because
they must first pass a gate on the cluster, which needs the corresponding author's login.
They implement sections 2 and 4 literally by importing the published modules
(`residual_analysis`, `split_replication_power`, `c1_run_test`, `c1_run_test_luad`,
`composition_controls`), and are frozen by **Amendment 1**: a dated file that lists the
scripts and their sha256, the pinned host of each cohort, and the gate and blind outputs
below, and is released (tag and Release) before the first Stage H production run. Amendment
1 may add file names, hashes, paths, hosts and gate outputs. It may not change any
definition, arm, cell, threshold, reading or sentence of this file. Stage L1's result files
are hashed in Amendment 1, since Stage L1 will have been read by then.

- *Gate H-V1 (published quantities only).* On the pinned host, the Stage H script recomputes
  the published estimator's observed column and permutation columns 1 to 5 and must equal
  the saved confirmatory matrix (whose sha256 on the cluster must equal the value in section
  1) to 1e-9 for every usable gene, with an identical usable gene list.
- *Blind run (counts only).* Cases, K, genes tested, RNA-seq genes entering the components,
  and markers found per family; no increment of any new arm.
- *Production.* Five arms per cohort and the export of the discovery scores. Each arm is
  run once. Outputs: `c1_spec_<cohort>_<arm>.npz` (`genes`, `M`),
  `c1_spec_discovery_scores_<cohort>.tsv` and `c1_spec_confirmatory_scores_<cohort>.tsv`,
  transferred unchanged to `server_export/results/c1_specificity/stage_h/` with their
  sha256 from the cluster logs.

**Stage L2 (local; one run).** `--stage L2`. Gate V3 first: the transferred discovery
scores, written at full double precision, must reproduce the published composition
controls (for every published-set gene, the gain of each score alone and of the three
together over the gene's own mRNA, columns `stroma_gain`, `immune_gain`,
`proliferation_gain` and `ALL_gain`, to 1e-9; these involve no morphology component);
otherwise Stage L2 stops. Then sets S_C3 and S_M|C3, and every remaining cell of section 5.

**Failures.** A run that fails before writing its result (a crash, the walltime, a node
failure) is kept with its log and repeated once, unchanged; if the rerun also fails, that
stage or arm is reported as not completed. A completed run is never repeated. A defect in a
frozen script is corrected only by a dated amendment released before the affected run, and
both the defect and the fix are disclosed.

**Results record.** The result files, logs and their sha256 are committed and tagged
`c1-posthoc-specificity-results-<date>` after Stage L2. The sentences of section 7 are then
filled in exactly as written, by hand from `cells_L1.csv` and `cells_L2.csv`, which hold
numbers only.

**Validation record (2026-11-20 UTC, before freezing).** [filled from `validate_only.log`
before the sidecar is written]
