# DRAFT v3 — NOT FROZEN until signed (§7) and registered (D3); no confirmatory slide image, RNA or protein file may be downloaded before that

# C1-LUAD replication: draft analysis rule

- **Drafted 2026-09-27.** Written from the Step 0 metadata-only feasibility check
  (`review/C1_LUAD_STEP0_2026-09-27.md`). It mirrors
  `review/C1_UCEC_FROZEN_RULE_2026-09-15.md`.
- **Rewritten 2026-09-29 (v2).** Section 0 was recorded after the adversarial review
  (`review/C1_LUAD_RULE_REVIEW_2026-09-29.json`). Sections 1-7 were rewritten after
  the metadata-only Step 0b (`review/C1_LUAD_STEP0B_2026-09-29/`). The previous
  version is archived as `C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.v1_0704.md`, sha256
  `e8dcb240146b14ffb0ba6606923a58fe5a5fe5d2fe86c853b6bd3ab6c0ab9c0b`.
- **Revised 2026-09-29 (v3).** v3 applies an 8-agent verification of v2. It also records
  the author's decisions in D6. v2 is archived as
  `C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.v2.md`, sha256
  `769e451a288fce35309a601a4a0246ad8c7df4390344c413aa5e36ca9fc05fc6`. A closing
  4-agent check of v3 was then applied.
- **Data touched so far.**
  - No confirmatory RNA or protein value has been downloaded.
  - No confirmatory slide pixel has been decoded or viewed.
  - Step 0b read PDC and GDC metadata and the IDC index.
  - For batch labels and pixel spacing, Step 0b also range-read the start of the base
    VOLUME instance of each of the 135 confirmatory primary series (158 HTTP range
    requests, logged in `query_log_idc.tsv`). The final range was 256 KiB for 125
    series and 4 MiB for 10.
  - These fixed-size ranges ran past the end of the DICOM header (about 5 KB to
    1.7 MB) into the PixelData element: about 0.11-2.4 MB per series and about 40 MB
    in total. For at least 80 series this is more than any offset table could occupy,
    so compressed tile bytes were fetched.
  - Those bytes were parsed with pydicom `stop_before_pixels`, never decoded, and
    discarded. Only the tag dump `idc_step3b_confirmatory_headers.json` was kept.
  - The reads were authorised as header-only and in fact overran the header. This is
    disclosed here and in the deviations table.
- **Status.** This is a draft until the corresponding author signs it (§7) and the
  third-party timestamp exists (D3).

## 0. Pre-signature decisions (recorded 2026-09-29)

Recorded after the 11-agent adversarial review of this draft
(`review/C1_LUAD_RULE_REVIEW_2026-09-29.json`) and **before** any LUAD pre-check has been
run and before any confirmatory slide, RNA or protein file has been downloaded. The
corresponding author chose all four options below on 2026-09-29. D5 was added the same
day, after Step 0b. D6 was added after the v2 verification, and it amends D1, D2 and D4
where D6 says so. Sections 1-7 were rewritten to match. Wherever sections 1-7 conflict
with this section, this section governs. Nothing in this section may be changed after
any pre-check output (D2) exists. Any run of `lsf_precheck_luad.sh`, smoke mode
included, counts as D2 output.

Correction to a rationale the review relied on: the 2026-09-29 operator-width result
(`review/recalc/operator_width_summary.csv`) used an **equal-rank** Gaussian arm
(`operator_width_control.py` L166-168), not a variance-matched one. Per main.tex
Limitations (i) such an arm cannot price the operator residualisation. It is therefore
**not** cited anywhere in this rule as evidence of confounding. D1 rests on two facts:
(1) the paper's abstract headline is the batch-corrected set (726 for LUAD); and
(2) the paper rests its discovery batch verdict on the batch-stratified permutation null.

### D1. Confirmatory family: fixed sequence H1 -> (a)(i) -> H2, one-sided alpha 0.05 at each step (sequence amended by D6a)

- **Step 1, H1 (identical to C1-UCEC).**
  - Selection set: the 2,310-gene uncorrected set.
  - K PCs, formula as in §3, fitted on raw confirmatory embeddings.
  - Null: unrestricted patient<->slide permutation, B = 1,000.
  - Result: p1. The §5 buckets (0.05 / 0.15 / 45 cases) are unchanged and read on p1 alone.
- **Step 2, batch qualifier (a)(i); tested only if p1 < 0.05.**
  - The same 2,310-set AUC is recomputed under the within-operator permutation null
    (§4; B = 1,000), giving p_strat.
  - Its wording is clause [B] of §5.2.
  - If p_strat >= 0.05, Results add this sentence next to p_strat (Results only):
    "The within-stratum null keeps every case inside its stratum, so any association
    carried by differences between strata, including genuine biology aliased with
    accrual, is credited to chance, and the fixed cases lower its power. A
    non-significant result therefore means the data cannot separate the signal from
    stratum structure, not that it is an artefact."
  - If p1 >= 0.05, p_strat is reported descriptively.
  - Both nulls are reported with mean, SD and the fraction of cases held fixed.
- **Step 3, H2, the batch-corrected set; claimable only if p1 < 0.05 and
  p_strat < 0.05.**
  - Selection set: the 726 genes (`luad__sitepack_operator.csv`,
    `fdr_batchresid < 0.05 & incremental_r2_batchresid > 0`).
  - Background: every tested gene not in the 726. This keeps the 1,679 genes that are
    in the 2,310 set but not the 726. That makes the contrast harder to win but does
    not change the test's size.
  - ρ is computed against discovery `incremental_r2_batchresid`.
  - Increments come from `RA.cv_r2`, unclipped, as in H1.
  - Case 1, the operator stratum is testable (§4):
    - Confirmatory PCs are residualised on the frozen operator dummy set of §4. This is
      the sitepack design, reimplemented, because it is inline in sitepack's `main()`.
    - Null: the within-operator permutations of §4, giving p2.
  - Case 2, operator not testable: H2 is computed with the H1 estimator and is
    descriptive only.
  - Wording: clause [D] of §5.2, pass or fail, but only when step 3 is reached
    (p1 < 0.05, p_strat < 0.05 and Case 1). Otherwise no [D] clause is printed in the
    Abstract, and Results report H2 descriptively as "not tested in the confirmatory
    sequence".
  - If p1 >= 0.05, or p_strat >= 0.05, or the operator stratum is not testable, H2 is
    descriptive, with no claim of replication of the batch-corrected set.
- **Stratum.**
  - One axis only, never joint.
  - Operator is primary. Scan quarter is used only if operator fails at Step 0b.
    Operator passes (§4), so scan quarter is never used.
  - A case's label is its value in `pinned/c1_luad_batch_labels.tsv`: the mode over its
    header-read GDC Primary Tumor slides, with ties going to the first slide in sorted
    slide-ID order. The label is frozen before signature and is not recomputed from the
    slides that later yield embeddings.
  - The one tied case, C3N-02281 (slides -22 and -23 by different operators), is
    labelled `81c9cdc4-1329-4611-a782-6e2f9e83e35d`. That operator has one case, so
    C3N-02281 is held fixed.
  - An axis qualifies only if all of these hold:
    - its DICOM attribute is validated on discovery slides' IDC copies against the
      SVS-derived `aperio.User` / `aperio.Date`;
    - it is present for >= 80% of analysed cases;
    - there are >= 2 strata with >= 3 cases each;
    - the largest stratum holds <= 90% of cases;
    - >= 50% of cases are permutable (singletons and unlabelled cases stay fixed).
  - The attribute names and the qualification outcome of each axis are fixed in the
    signed rule from Step 0b metadata, never after data exposure.
  - The only later check is the §4 re-check on the realised analysed set, done from the
    coverage table before any statistic is computed. It can only downgrade (a)(i), H2
    and (a)(ii) to "not testable"; it never substitutes another axis.
- **Co-reported (descriptive, not a gate): completeness-stratified AUC.**
  - Strata: genes whose discovery per-gene n equals 105, versus the rest. The realised
    stratum sizes after intersection are reported.
  - AUC_compl = sum over strata of U_s, divided by the sum of n1_s x n0_s. Its one-sided
    p is p_compl. Its null uses the same B = 1,000 permutation columns as H1.
  - Also reported: the within-stratum pooled ρ, and the fraction of positive increments
    per stratum for selected and background genes.
  - Fixed wording: if p1 < 0.05 and p_compl >= 0.05, write "replicated at the
    pooled level; not robust to protein-completeness stratification". It goes verbatim
    in Results, and in the Abstract as clause [C] of §5.2.
  - `c1_run_test_luad.py` must save the genes x (1+B) matrix M (.npz, with its sha256 in
    the result JSON) and the per-gene confirmatory n.

### D2. LUAD split-half pre-check: discovery data only, run before signature

- **Script.** `split_replication_power.py`, unchanged (sha256 must equal the one recorded
  by the C1-UCEC run, 79530c07...), with its `--group-csv` mode.
- **P0.** Build `luad_batch_groups.csv` on the 105-case discovery analysis set
  (`figures/figdata/scores_luad.csv`). Its sha256 and build date go in §7.2 before any
  P1, P1b or P2 run. Columns:
  - `q_early_vs_late`: 2016Q3-2017Q2 = sel, 2017Q3-Q4 = test.
  - `q_2017Q3_vs_rest`: 2017Q3 = test.
  - `op_balanced_0..4`: greedy operator balance with seeds 0-4, operators kept whole.
  - `stream_c3n_vs_c3l` (descriptive, D6d): sel = C3N (71 cases), test = C3L (34).
  - Every column must meet the script's own floor, n_sel >= 40 and n_test >= 30.
  - The first seven columns follow the column schema of `ucec_batch_groups.csv`. No
    builder script for that UCEC file survives. The construction
    (`make_batch_groups.py`, `build_luad_batch_groups.py`) was therefore reconstructed
    by inference from the UCEC file's values and from the script's floors. Parity with
    the UCEC pre-check groups holds at the schema level, not the algorithm level; see §6.
  - Case labels use all local Primary Tumor slides. The script intersects the file with
    its own RNA ∩ protein ∩ WSI set on HPC.
- **P1.** 10 random splits: `--n-test 50 --splits 10 --perms 200 --top-frac 0.15 --seed 0`,
  run as one job.
- **P1b (descriptive, D6d).** 10 random splits at the stream split's size:
  `--n-test 34 --splits 10 --perms 200 --top-frac 0.15 --seed 1`, run as one job. This
  is the matched-n reference for `stream_c3n_vs_c3l`.
- **P2.** One run per group column of P0, 200 permutations each. The seven UCEC-schema
  columns gate the reading below. `stream_c3n_vs_c3l` is descriptive:
  - If its p_auc < 0.05, it is read as within-discovery transfer of the selection across
    the C3N/C3L accrual boundary at n_test = 34. It is never called a replication
    (§5.7 item 2).
  - If its p_auc >= 0.05 and P1b passes (>= 8/10 splits with p_auc < 0.05), it is read
    as evidence of accrual dependence.
  - If its p_auc >= 0.05 and P1b does not pass, it is read as uninformative at
    n_test = 34.
  - It never enters GO / MARGINAL / NO-GO.
- **Reading, fixed now** (P2 counts refer to the seven UCEC-schema columns).
  - GO: P1 >= 8/10 **and** P2 >= 6/7 with p_auc < 0.05.
  - NO-GO: P1 < 5/10 **or** P2 < 4/7.
  - MARGINAL: anything else.
- **Consequences.**
  - GO: proceed to signature.
  - MARGINAL: proceed to signature, with the §5 ">= 0.15" bullet reworded to "not
    replicated at marginal pre-checked power; not evidence of absence".
  - NO-GO: the corresponding author decides whether to download.
    - If C1-LUAD is not run, the paper reports that it was halted at the pre-check,
      with the pre-check numbers.
    - If it is run anyway, the ">= 0.15" bullet is pre-declared uninformative.
- **Scope.**
  - The pre-check gates only go / marginal / no-go. It cannot change D1, K, the null or
    any reading.
  - Every P1/P1b/P2 row is reported whatever the outcome.
  - No residualised-selection arm is run, because `split_replication_power.py` stays
    unchanged. The (a)(i) within-operator null is not rehearsed either. H2 and the
    (a)(i) qualifier therefore enter without any pre-check, and how often the qualifier
    fires at n of about 112 was unknown at signing.
  - Any `lsf_precheck_luad.sh` run, smoke mode included, is D2 output for the §0 lock.
  - Known asymmetries:
    - The selection half has ~55 patients with top-15% raw selection, versus FDR
      selection at n = 105 in the real test, so NO-GO is weaker evidence than GO.
    - Within-discovery splits share TMT plex.
    - UCEC's realised z was 2.73 against a pre-check median of 4.26.

### D3. Third-party timestamp: OSF registration plus GitHub Release, before any confirmatory download

- **Order of steps.**
  1. Sign (§7.3).
  2. Compute the signed file's SHA-256 and keep it outside the file.
  3. Register on OSF: "Secondary Data Preregistration" template, embargo allowed.
     Attach the signed rule, the Step-0/0b notes with their dated errata
     (`review/C1_LUAD_STEP0_ERRATA_2026-09-29.md`), the 113-case ID list and the hash
     manifest.
  4. Place the signed file and a `.sha256` sidecar at a tracked path in the public repo,
     push, and create a GitHub Release. Pushing requires the author's go-ahead.
  5. Write the OSF DOI, the release tag and the timestamps with timezone into a receipt
     file.
  6. Only then may confirmatory data be downloaded. Keep the PDC/GDC/IDC download logs.
- **Paper wording for LUAD.** "registered on OSF on <date> (<DOI>) after metadata checks
  that included partial byte-range reads of the confirmatory slide files to obtain DICOM
  header fields (pixel bytes within those ranges were discarded without decoding), and
  before any complete confirmatory slide, transcript or protein file was downloaded
  (download and range-read logs in the repository)".
- **UCEC.** Its existing disclosure (no third-party timestamp) stands and is not upgraded.

### D4. Already-published C1-UCEC: three post hoc items, decided before any LUAD result exists (item 3 added by D6b)

- **Completeness-stratified AUC.** Recompute with its permutation null by a
  deterministic rerun: same seed `crc32(b"c1_ucec|perm")`, same B = 1,000, same cases,
  genes and K. The observed value from the published gene table is 0.561, against the
  pooled 0.597. Label it post hoc and report it in Results whatever the result.
- **PathDB-323 slide sensitivity.** It was pre-declared in UCEC Addendum 2 item 1 but
  never executed. Methods states exactly that. It is not run now.
- **Batch-stratified re-read (RUN).**
  1. **Validate, then read.** On the IDC copies of UCEC discovery slides, validate that
     the DICOM ImageComments `User =` field equals `aperio.User` in
     `slide_acquisition_meta.tsv`, as §4 did for LUAD. Then read the same field from the
     headers of the UCEC confirmatory Primary Tumor slides. Those files are already on
     the HPC and were already analysed, so this is not a blinding question.
  2. **Labels.** A case's label is the mode over its slides in
     `scripts/pinned/slide_type_map_ucec_c1.tsv` (the 169 analysed Primary Tumor slides)
     whose header carries the field. Ties go to the first slide in sorted slide-ID
     order. Labels are frozen before the re-read. The D1 qualification criteria are
     evaluated on the published UCEC analysed set.
  3. **Re-read.** Recompute the published 2,566-set AUC under a within-operator null:
     B = 1,000, seed `crc32(b"c1_ucec|perm_op")`, with singleton operators and
     unlabelled cases held fixed.
  4. **Reporting.** It is labelled post hoc and reported whatever the result.
  - If the attribute does not validate, or any D1 criterion fails on the UCEC analysed
    set, Methods state that C1-UCEC was tested under the unrestricted null only and
    carries no batch qualifier.
- In every case, UCEC's published §5 reading is unchanged.

### D5. Decisions after Step 0b (recorded 2026-09-29, before any D2 output)

Step 0b ran on 2026-09-29 under the author's authorisation
(`review/C1_LUAD_STEP0B_2026-09-29/`). It read metadata, plus DICOM header range reads
that were authorised as header-only but ran past the header (see the preamble).

- **What it reproduced.** Every Step 0 count that it re-derived from source: case
  counts, GDC and IDC slide counts and coverage, and demographics. It also derived the
  operator and scan-quarter labels.
- **What it did not re-derive.** It did not re-query PathDB. So the §1.3
  embedding-medium, FFPE and OCT figures are carried over unchanged from Step 0
  (2026-09-27) and have not been independently re-verified. IDC's own
  `embeddingMedium_CodeMeaning` and `tissueFixative_CodeMeaning` fields are empty for
  every matched slide.

The corresponding author chose D5a-D5d on 2026-09-29. D5e-D5g are drafter defaults; they
stand unless the author overrides them before any D2 output exists.

- **D5a.** The all-zero operator UUID (5 confirmatory cases; §4) is treated as its own
  operator stratum, taken literally.
- **D5b.** The difference in recorded ancestry between the two waves is disclosed in
  §6. It is not adjusted for, and no ancestry-adjusted AUC is added.
- **D5c.** Independence is patient-level. The seven patients who are also discovery
  Primary Tumor cases are excluded. This includes C3N-00545, even though its
  confirmatory aliquot differs from its discovery aliquot. The case set stays 113.
- **D5d.** D2 is run only after this rewrite. P1 runs as one job (`--splits 10
  --seed 0`). The reason: split s uses `RandomState(seed*1000+s)` and its permutation
  seed is keyed on the split label, so per-seed parallel jobs would define different
  splits.
- **D5e.** Scale: no rescaling, which is the discovery and C1-UCEC policy. A
  descriptive sensitivity analysis (s) covers the five non-standard-grid cases (§1.3).
- **D5f.** The PathDB-pooled slide sensitivity analysis of the previous draft is
  dropped, for three reasons:
  - PathDB rows and IDC series are not the same slide sets (C3N-02282 has no PathDB
    row).
  - The analogous UCEC analysis was pre-declared and never run (D4).
  - The GDC criterion is the discovery criterion.
- **D5g.** Confirmatory smoke runs are allowed only with `--blind` (§5.6).

### D6. Decisions after the v2 verification (recorded 2026-09-29, before any D2 output)

The corresponding author chose these on 2026-09-29.

- **D6a.** The confirmatory family is a three-step fixed sequence: H1, then the (a)(i)
  qualifier, then H2, each at one-sided alpha 0.05, stopping at the first
  non-rejection. H2 is claimable only if p1 < 0.05 and p_strat < 0.05. This bounds the
  familywise error of every positive claim the §5.2 clauses can print.
- **D6b.** C1-UCEC gets a post hoc batch-stratified re-read (D4, item 3).
- **D6c.** The Abstract quotes the bootstrap interval for both cohorts (§5.4, §5.7).
- **D6d.** Three descriptive additions, none of which can change §5:
  - the stream-disjoint pre-check split with its matched-n reference (D2 P0, P1b, P2);
  - disclosure of the confirmatory pre-analytic field fill rates (§6), with no
    covariate analysis;
  - the single-scanner sensitivity analysis (s2) (§3).

## 1. Cohort, case set and data

### 1.1 Studies

- **Confirmatory**: CPTAC LUAD Confirmatory Study - Proteome, **PDC000489**. It is TMT11,
  with 27 TMT11 runs in PDC's run-level experimental design. It is the proteome
  fraction, not the phosphoproteome (PDC000490) or the acetylome (PDC000491).
- **Discovery**: PDC000153, TMT10, 25 runs.

### 1.2 Case set: 113 cases, frozen by file

- The case set is `server_export/scripts/pinned/c1_case_ids_luad.txt`: 113 IDs, sorted,
  sha256 `3629b243a3aa21057dae24e7f5c1b47e8a2571e6ffbcff6f3170333e1ca9fe44`.
- **Derivation.** `c1_luad_derive_cases.py` reads only the raw `biospecimenPerStudy`
  JSON saved on 2026-09-29 under `review/C1_LUAD_STEP0B_2026-09-29/raw/`, and every
  count below is FATAL-asserted.
  1. **Pseudo-cases.** PDC000489 lists 131 case labels over 245 aliquot rows. Eleven
     labels (12 rows, `sample_type` = "Not Reported") are QC, bridge or pooled
     pseudo-cases: `CR.C1`, `CR.C1.newly.combined`, `CR.C2`, `CR.D`, `CR.D.MAM`,
     `CR.D.newly.combined`, `CR.D.pool.2`, `CR.LSCC`, `Taiwanese IR`,
     `Tumor Only CONF`, `Tumor Only DISC`. That leaves 120 real cases with a Primary
     Tumor aliquot.
  2. **Discovery patients.** Seven of the 120 patients are also Primary Tumor cases of
     PDC000153 (intersection by `case_submitter_id`): C3L-02348, C3L-02350, C3N-00545,
     C3N-00738, C3N-01023, C3N-01024, C3N-02003.
     - Six of them carry the identical `aliquot_submitter_id` in both studies, i.e. a
       cross-batch bridge re-measurement.
     - **C3N-00545** carries a different confirmatory aliquot (CPT0066910003,
       Qualified) from its discovery aliquot (CPT0066890003, flagged Disqualified in
       PDC). It is among the 105 analysed discovery cases.
     - All seven are excluded. **Independence is patient-level, not aliquot-level
       (D5c).**
  3. **Result.** 113 cases remain, all C3L-/C3N-. The excluded labels and cases are
     in `pinned/c1_excluded_luad.tsv`.
- **Independence check, fail-closed** (`c1_luad_independence_check.py`).
  - Result: 0 of 113 cases overlap the discovery slide map (`slide_type_map_luad.tsv`,
    111 cases) or the discovery aliquot crosswalk (`aliquot_to_case_tumor_luad.tsv`,
    111 aliquots), at case level and at aliquot level.
  - Its built-in self-test re-adds the seven discovery patients in memory and must
    exit FATAL.
  - Run 2026-09-29: normal check PASS; self-test FATAL, as required.
  - The analysis script repeats this check before loading any confirmatory file, and
    asserts that the analysed cases are a subset of the frozen list.
- The case set is not re-filtered or expanded after signature.

### 1.3 Slides

- **Criterion.** A slide counts only if its GDC slide entity has
  `sample_type = Primary Tumor`. This is the discovery criterion (UCEC Addendum 2,
  item 1).
  - Confirmatory: 136 slides cover 113 of 113 cases (93 cases with one slide, 17 with
    two, 3 with three or more). Files: `pinned/c1_primary_slide_ids_luad.txt`, and the
    full entity table `pinned/c1_luad_gdc_slides.tsv` (136 Primary Tumor + 123 Solid
    Tissue Normal).
  - Discovery reproduced slide-for-slide: 141 Primary Tumor / 106 Solid Tissue Normal,
    identical to `slide_type_map_luad.tsv`.
- **IDC availability.** 135 of the 136 primary slides have an IDC series, joined by
  `ContainerIdentifier` across all `cptac_luad` SM series.
  - **C3N-02142-23**, the only primary slide of C3N-02142, has no IDC series. No
    alternative source (PathDB SVS, TCIA) is sought.
  - **Expected analysed n = 112**, so K = 22 by §3. C3N-02142 stays in the frozen
    list and is reported as lost to IDC coverage.
- **Pathways and pixel grid of the 135 slides.**

  | Pathway (Manufacturer / model) | Transfer syntax | Level-0 spacing | Slides |
  |---|---|---|---|
  | Leica Biosystems / "Aperio converted by com.pixelmed.convert.TIFFToDicom" | JPEG Baseline | 0.49 µm | 122 |
  | same | JPEG Baseline | 0.50 µm | 8 |
  | same | JPEG Baseline | 0.25 µm (40x): C3L-03717-21, C3L-03721-21 | 2 |
  | PixelMed / com.pixelmed.convert.TIFFToDicom | JPEG Baseline | 0.32 µm, no objective recorded: C3L-02513-23, C3L-02515-23, C3L-03642-21 | 3 |

  Each of the five slides not at 0.49-0.50 µm is its case's only primary slide.
- **Discovery slides in IDC.** All 141 discovery primary slides are Leica-pathway
  copies of the SVS files the discovery pipeline read: 139 JPEG Baseline at 0.49 µm,
  and 2 JPEG 2000 at 0.32 µm (C3L-01924). Their IDC spacing agrees with the SVS mpp
  to within 0.004 µm for 141 of 141. Their ImageComments tag carries the original
  Aperio description (§4). The Leica pathway is therefore a conversion of Aperio SVS
  scans, not a different format.
- **Scanners.** ScanScope IDs are taken from ImageComments.
  - Confirmatory primary slides: SS1553 122; SS7559 8 slides (7 cases, 0.50 µm,
    Nov-Dec 2018); SS7320 2 slides (2 cases, 0.25 µm, 40x). The 3 PixelMed slides carry
    no scanner ID.
  - Discovery primary slides: all SS1553 (139; 2 blank).
  - So 9 of the 112 expected analysed cases were scanned only on devices absent from
    LUAD discovery (SS7559: C3L-03985, C3L-04757, C3N-03765, C3N-04157, C3N-04168,
    C3N-04176, C3N-04180; SS7320: C3L-03717, C3L-03721). A further 3 cases (C3L-02513,
    C3L-02515, C3L-03642) have no recorded scanner.
  - The 9 cases fall in operator strata that contain no SS1553 slide:
    - the all-zero stratum (SS7320, 2 cases; SS7559, 3 cases);
    - `ad6974fd` and `e5e971cd` (SS7559, 2 cases each).
  - The within-operator null never exchanges these cases with SS1553 cases. The (a)(i)
    and H2 nulls therefore credit their difference from SS1553 cases to between-operator
    structure.
  - The 3 cases with no recorded scanner carry no operator label and are held fixed.
  - Sensitivity analysis (s2) removes all 12 cases.
- **Scale policy (D5e).**
  - Tiles are 256 px, at level 0, dense, with no rescaling. This is identical to
    discovery, which tiled its two 0.32 µm slides unrescaled, and to C1-UCEC.
  - The two 0.25 µm slides lie outside the discovery range. They are tiled unrescaled
    like the rest, and the descriptive sensitivity analysis (s) removes the five
    non-standard-grid cases.
  - `extract_features.py`'s `--level > 0` path is not used, because it passes
    level-k coordinates to openslide.
- **Slide preparation.** These figures are as reported by Step 0's 2026-09-27 PathDB
  join; Step 0b did not re-derive them (see §6).
  - Confirmatory: of the 136 primary slides, 135 match a PathDB row and all 135 are
    FFPE; C3N-02282's slide has no PathDB row.
  - Discovery primary slides: 126 FFPE and 15 OCT.
- **Extraction.**
  - Use `extract_c1_dicom.py` unchanged. It imports the tiling, tissue filter, UNI
    encoder and transform from `extract_features.py`.
  - Run with `--only-file pinned/c1_primary_slide_ids_luad_idc.txt`. That list is the
    136-ID list minus C3N-02142-23, which has no IDC series; the unchanged script exits
    on any listed ID it cannot find.
  - Before extraction, run a `--dry-run` over all downloaded series. The script
    resolves every series before its dry-run loop and stops at the first one it cannot
    resolve. Step 0b found ContainerIdentifier = GDC slide ID, with the PatientID
    prefix, for 135 of 135 series. An abort is therefore handled as a reader failure
    under §1.4.
  - The slide map is `pinned/slide_type_map_luad_c1.tsv`, built by
    `build_slide_map_c1.py` from the 136 IDs and pinned before signature. Its local
    run on 2026-09-29 gave 136 of 136 Primary Tumor.

### 1.4 Reader QC, before any extraction

- **Thresholds.** These are fixed inside `qc_reader_equivalence.py`: pooled mean |Δ|
  < 1 grey level, and every channel's mean signed Δ within ±0.5.
- **Leica pathway.**
  - (i) Confirmatory series C3L-00444-23 (0.49 µm) and C3L-03717-21 (0.25 µm, 40x):
    openslide-DICOM versus wsidicom on identical level-0 regions.
  - (ii) Discovery slide C3L-00001-21 (0.4942 µm, scanner SS1553): its original SVS
    versus its IDC Leica-pathway copy.
- **PixelMed pathway.** Arm (i) only, on C3L-02513-23. No discovery slide exists in
  this pathway, so arm (ii) is not available; this is disclosed (§6).
- **How the slides were chosen.** Each is the first in sorted order within its
  pathway x spacing category; for discovery, only C3-style IDs were considered.
- **Failure handling.**
  - A DIFFERENT result stops extraction until the cause is found.
  - A fix is a dated code correction. It is limited to the decoder or reader and must
    pass the same identical-grid QC.
  - A fix cannot change the slide list, the scale policy, the tiling, the tissue
    filter or the encoder.
  - No pathway is dropped because of a reader failure.
  - This is outcome-blind, because of the order of operations in §1.6.

### 1.5 Omics

- **Protein.**
  - Source: the PDC000489 `quantDataMatrix` with `data_type` log2_ratio, the same
    route as the published C1-UCEC test.
  - Crosswalk: restricted, before any averaging, to Primary Tumor aliquots of the 113
    frozen cases. Pseudo-case rows and the seven excluded patients never enter.
  - FATAL if any kept `aliquot_submitter_id` appears in the discovery crosswalk.
  - FATAL on duplicate gene symbols.
  - The pull writes a per-case coverage table and per-gene NaN counts.
  - Why a new script: `c1_pull_omics.py`'s own crosswalk filter keeps "Not Reported"
    pseudo-case rows (found 2026-09-29 on PDC000153; 3 of 4 matched real matrix
    columns). So the LUAD pull uses a new, hashed script, and the executed UCEC
    script is not re-parameterised.
- **Discovery route check, before signature and on discovery data only.**
  - `protein_route_check.py` compares the PDC000153 log2_ratio matrix pulled on
    2026-09-29 (`pinned/luad_discovery_protein_via_pdc.csv`, 110 cases x 11,029
    genes) with `RA.load_protein_matrix`.
  - C3N-00545 has no column in that matrix, although its CDAP TSV aliquot entered the
    discovery analysis.
  - Pre-declared consequences:
    - If missing values in the PDC route turn out to be encoded as 0 or another
      sentinel, they are recoded to NaN. This is a dated code correction, applied to
      C1-UCEC as well, and disclosed.
    - Otherwise the log2_ratio route is kept unchanged, whatever the concordance.
      Any systematic difference (shared vs unshared peptides, CDAP version) is
      disclosed for both C1 tests.
- **RNA.** GDC STAR-Counts, tumour only; 113 of 113 cases have a file (Step 0). It is
  loaded with `RA.load_rna_matrix`, unchanged.

### 1.6 Order of operations after signature (blinding)

1. D3 registration and GitHub release; the receipt is written.
2. Download the 135 primary DICOM series listed in
   `pinned/c1_primary_slide_ids_luad_idc.txt` and nothing else, keeping the download
   log.
3. Reader QC (§1.4).
4. Full extraction. The embedding manifest (per-slide `.pt` sha256, n_tiles, status)
   is written and dated.
5. Only then are RNA and protein pulled, and the coverage table written.
6. Blind smoke run (§5.6), then the B = 1,000 run, then the bootstrap.

No confirmatory RNA or protein file is opened before step 5.

## 2. Selection sets (discovery; not re-derived)

- **Files.**
  - `server_export/pinned/luad/residual_results_tumoronly.csv`, byte-identical to
    `server_export/results/luad__residual_results_tumoronly.csv`.
  - `server_export/results/luad__sitepack_operator.csv`.
  - The HPC copies must hash identically to these (§7.2).
- **H1 set.** `fdr < 0.05 & incremental_r2 > 0`: **2,310** of 10,724 tested genes
  (modal n = 105).
- **H2 set.** `fdr_batchresid < 0.05 & incremental_r2_batchresid > 0`: **726** genes.
  - The column pair is not `fdr_both`/`incremental_r2_both`, which gives 802 in the
    operator file.
  - The "845" in the previous draft came from the quarter file.
  - 631 of the 726 genes are in the 2,310 set; 95 are not.
- **Background.**
  - H1: every other discovery-tested gene in the usable universe.
  - H2: every usable gene not in the 726, including the 1,679 genes in the 2,310 set
    but not in the 726.
- **Usable universe.**
  - Discovery-tested genes with a confirmatory RNA column and at least 30 non-missing
    confirmatory protein values (MIN_N = 30, as in C1-UCEC).
  - Symbols are matched exactly, with no remapping.
  - The number of discovery-tested genes absent from the confirmatory data is
    reported.

## 3. Statistics (one-sided upper tail throughout)

- **Capacity.**
  - K = round(20 x n / 100), where n is the realised number of analysed cases
    (expected 112, giving K = 22).
  - This is the project convention K/n = 0.2 used by the pre-checks and by C1-UCEC.
    It is not a ratio match to LUAD discovery (20/105 = 0.19), and the Methods say so.
  - PCA(K, `svd_solver="full"`) is fitted on the confirmatory mean-pooled embeddings
    only.
- **Increments.** `RA.cv_r2` (5-fold KFold, shuffle, random_state 0; closed-form
  ridge, standardised), as in C1-UCEC.
- **H1.** Mann-Whitney AUC, selected versus background genes. The secondary statistic
  is Spearman ρ against discovery `incremental_r2`; it is reported, not a second bar.
- **(a)(i) qualifier (step 2).** H1's AUC under the within-operator null (§4), giving
  p_strat.
- **H2 (step 3).** The same statistic for the 726 set, with increments from `RA.cv_r2`
  on the operator-residualised K PCs (§4), under the within-operator null, giving p2.
  ρ is computed against discovery `incremental_r2_batchresid`.
- **Rho p-values.** Each ρ p-value uses the null of its own AUC.
- **Completeness-stratified AUC** (D1, co-reported).
  - Strata: genes whose discovery per-gene n = 105, versus n < 105.
  - Composition before intersection: selected 2,211 / 99, background 5,708 / 2,706.
    The realised sizes are reported.
  - AUC_compl = Σ_s U_s / Σ_s n1_s·n0_s, with ranks computed within each stratum, giving
    p_compl. Its null uses H1's permutation columns.
  - Also reported: the within-stratum pooled ρ, and the fraction of positive increments
    per stratum for selected and background genes.
- **Descriptive only.** These cannot change §5. Apart from the §5.4 sentence, they carry
  no wording trigger. They are reported in one table whatever the outcome.
  - **(a)(ii)** The 2,310 set on H2's estimator: residualised PCs, within-operator null,
    the same permutation draws as H2.
  - **(b) Per-gene rediscovery**, in two versions:
    - **(b1) Uncorrected discovery estimator.** K = 20 fixed (`residual_analysis.N_PCS`),
      with `residual_analysis._process_gene` and `bh_fdr` unchanged. B = 1,000
      permutations are shared across genes, as in discovery, drawn from
      `RandomState(crc32(b"c1_luad|pergene"))`. BH runs over the usable universe.
    - **(b2) Batch-corrected estimator.** K = 20 fixed (sitepack `N_PCS`). PCs are
      residualised on the frozen dummy set (§4). The estimator is sitepack's `_fit` with
      `plain_folds`, reimplemented rather than imported (Appendix item 1):
      - per gene, on that gene's valid cases;
      - 5-fold KFold (shuffle, random_state 0);
      - ridge with alpha 1.0 on features standardised within each training fold;
      - R² clipped to [-5, 1];
      - increment = R²(RNA + residualised PCs) - R²(RNA).

      This is the estimator that defined `incremental_r2_batchresid` and the 726.
      B = 1,000 within-operator permutations are drawn from
      `RandomState(crc32(b"c1_luad|pergene_op"))`, with unlabelled cases held fixed. BH
      runs over the usable universe. The count is fdr < 0.05 & increment > 0.

    For each version, report:
    - m (the number of genes tested);
    - the number of BH rejections before the increment > 0 filter, and the count after
      it;
    - the overlap of the count with the 2,310 set and with the 726 set, against the
      hypergeometric expectation;
    - whether the count reaches 50, the discovery positivity criterion (itself set post
      hoc);
    - a censoring diagnostic:
      - the number of genes at p = 1/(B+1), and the number at p <= 2/(B+1);
      - k* = ceil(m / (0.05·(B+1))). With B = 1,000, BH returns either no gene or at
        least k* genes before the sign filter, so reaching 50 mostly records whether BH
        rejects anything at all;
      - a zero count while any gene sits at the floor is reported as "zero at
        B = 1,000, limited by the permutation floor; not evidence of absence";
      - a non-zero count with genes at the floor is reported as "resolution-limited at
        B = 1,000".
  - **(s) Grid/pathway sensitivity.**
    - H1 only, under the unrestricted null.
    - Cases: the n_s cases other than the five non-standard-grid cases (C3L-02513,
      C3L-02515, C3L-03642, C3L-03717, C3L-03721).
    - K = round(20·n_s/100), which is 21 at 107.
    - B = 1,000 permutations of n_s rows, drawn afresh from
      `RandomState(crc32(b"c1_luad|perm_s"))`.
  - **(s2) Single-scanner sensitivity (D6d).**
    - H1 only, under the unrestricted null.
    - Cases: the n_s2 cases other than the 12 with no SS1553 slide (§1.3).
    - K = round(20·n_s2/100), which is 20 at 100.
    - B = 1,000 permutations of n_s2 rows, drawn afresh from
      `RandomState(crc32(b"c1_luad|perm_s2"))`.
  - **Bootstrap interval** (§5.4).

## 4. Nulls, strata and residualisation

- **Unrestricted null** (H1, completeness-stratified AUC).
  - B = 1,000 permutations of the case<->PC-row pairing, drawn once from
    `RandomState(crc32(b"c1_luad|perm"))` and shared by every gene.
  - p = (1 + #{null >= obs}) / (1 + B).
  - The null mean and SD are reported; 0.5 is never assumed.
- **Within-operator null** ((a)(i), H2, (a)(ii)).
  - Permutes the case<->PC-row pairing within operator strata, using the frozen labels.
    Singleton operators and unlabelled cases are held fixed.
  - B = 1,000 draws, taken once from `RandomState(crc32(b"c1_luad|perm_op"))`. Within
    each draw, strata are processed in sorted label order. The same draws are shared by
    (a)(i), H2 and (a)(ii).
  - What the null tests: it credits all between-operator association to chance and
    tests only within-operator association.
  - The fraction of cases held fixed is reported.
- **Operator labels, fixed before signature.**
  - Source: DICOM (0020,4000) ImageComments of each slide's base VOLUME instance
    carries the original Aperio description string. Operator = its `User =` field;
    scan date = its `Date =` field.
  - Validated on the discovery IDC copies against the SVS-derived `aperio.User` and
    `aperio.Date`: operator 139 of 139 (the two slides of C3L-01924 have no operator
    on either side); date 141 of 141.
  - Confirmatory: present on 132 of 135 slides; absent on exactly the three
    PixelMed-pathway slides.
  - Per-case label: the modal value over the case's header-read primary slides, with
    ties going to the first slide in sorted slide-ID order. The one tie is C3N-02281
    (D1).
  - Frozen in `pinned/c1_luad_batch_labels.tsv`. It was derived by
    `idc_step3c_derive_batch_labels.py` from the Step 0b range reads described in the
    preamble; the pixel bytes in those ranges were fetched but never decoded or stored
    (log: `query_log_idc.tsv`).
- **All-zero UUID** (`00000000-0000-0000-0000-000000000000`; C3L-03717, C3L-03721,
  C3N-03765, C3N-04176, C3N-04180). This value is present verbatim in the Aperio header
  of two scanners used in 2018. It is its own operator stratum (D5a).
- **Qualification**, by the D1 criteria, evaluated 2026-09-29 from Step 0b metadata.

  | D1 criterion | Required | 113 frozen cases | 112 expected analysed cases |
  |---|---|---|---|
  | Coverage | >= 80% | 109/113 = 96.5% | 109/112 = 97.3% |
  | Strata with >= 3 cases | >= 2 | 13 (of 32 operators) | 13 |
  | Largest stratum | <= 90% of cases | 16/109 = 14.7% | 16/109 = 14.7% |
  | Permutable cases | >= 50% | 98/113 = 86.7% | 98/112 = 87.5% |

  Unlabelled cases: C3L-02513, C3L-02515 and C3L-03642 (PixelMed only), and C3N-02142
  (no IDC slide; it is not among the 112).
  - **Operator qualifies; scan quarter is not used.**
  - Scan quarter, recorded for D1: it also qualifies (6 strata, 2017Q3-2018Q4; largest
    30/109 = 27.5%; 109/113 labelled).
  - The four operator criteria are re-checked on the realised analysed set, from the
    coverage table, before any statistic is computed. If any fails there, (a)(i), H2
    and (a)(ii) take the "not testable" wording (§5.2 clause [B], third form), and no
    other axis is substituted (D1).
- **Residualisation** (H2, (a)(ii), (b2)).
  - The sitepack design is inline in `residual_analysis_sitepack.py`'s `main()`, which
    runs only on discovery paths. It is therefore reimplemented in
    `c1_run_test_luad.py`, from sitepack L304-322, with two stated departures:
    1. **Unlabelled cases are not a level.** They get no dummy and are held fixed in
       every within-operator permutation. Sitepack would instead pool them into a
       dummy-coded, permutable `_missing` level.
    2. **The dummy set is frozen now,** from the frozen labels on the 112 expected
       cases. It is the 12 operators ranked first by case count, with ties broken by
       first appearance in sorted case order:

       | # | Operator | Cases |
       |---|---|---|
       | 1 | `9c62eb1b-888b-4319-a031-5ee965150a3e` | 16 |
       | 2 | `888d279c-9467-4931-8ead-381499dea990` | 11 |
       | 3 | `975375e8-786f-45d0-bfdd-9e70146a3eb8` | 10 |
       | 4 | `fc28642d-0894-46ef-9386-1ba3b5bdd99c` | 8 |
       | 5 | `6ed758c2-29d4-4ccb-945c-a7454ee330f5` | 6 |
       | 6 | `f3d46998-d392-4ce1-bd72-cf82622a2333` | 5 |
       | 7 | `00000000-0000-0000-0000-000000000000` | 5 |
       | 8 | `c8152987-eedd-4cf7-9b33-fee82a532dbc` | 4 |
       | 9 | `56d2978b-1a1d-4c5a-a722-3daa9816a1f5` | 4 |
       | 10 | `c1a7d89d-1759-40c1-af7a-8ba9c7901cec` | 4 |
       | 11 | `052aeb46-e02a-4dae-ae21-c92500e0b9d4` | 3 |
       | 12 | `b698ee8d-5921-4ead-bf67-ee1d6ed62b52` | 3 |

       `925cc759-f2ba-47ea-b2dc-fb1e39fa77da` (3 cases) is left to the intercept.
  - The K confirmatory PCs are residualised on these dummies by least squares, with an
    intercept.
  - If a listed operator falls below 3 realised cases, it loses its column and no
    replacement is added.
  - MAX_LEVELS = 12, MIN_GROUP = 3 and MIN_PERMUTABLE = 0.50 are sitepack's values,
    copied as literals. The MIN_PERMUTABLE warning is printed.
- **Cross-wave context.** Only 2 of the 32 confirmatory operators appear in discovery.
  Confirmatory scan quarters run 2017Q3-2018Q4; discovery quarters run 2016Q3-2017Q4.

## 5. Decision rule

### 5.1 Primary reading (thresholds identical to C1-UCEC §5), on p1 alone

- **Fewer than 45 usable cases (§5.5):** abandoned. The test is not run and the
  availability finding is reported.
- **p1 < 0.05:** set-level ranking replicated, worded only with the §5.2 clauses.
- **0.05 <= p1 < 0.15:** suggestive; not declared a replication.
- **p1 >= 0.15:** not replicated. The wording depends on the D2 outcome, which is
  recorded in §7.4 before signature:
  - GO: "not replicated";
  - MARGINAL: "not replicated at marginal pre-checked power; not evidence of absence";
  - NO-GO, but run anyway: "uninformative at pre-checked power".

Every D2 row is reported whatever the outcome, and the D2 reading is not revised after
the confirmatory result. This replaces the UCEC sentence about not questioning a
within-cohort pre-check.

### 5.2 Wording clauses (fixed now)

When p1 < 0.05, the LUAD reading is built from these clauses, in this order and in no
other. In the Abstract it sits in the same sentence as the UCEC reading (§5.7 item 1).

- **[A]** "lung (N patients): genes selected by the uncorrected discovery analysis
  (2,310) outranked background (AUC x, 95% patient-bootstrap interval a-b; one-sided
  p = p1), a set-level replication".
- **[B]** Exactly one of:
  1. ", including under an acquisition-batch (operator)-stratified null (p = p_strat)"
  2. " under an unrestricted permutation null but not under an acquisition-batch
     (operator)-stratified null (p = p_strat)"
  3. "; acquisition-batch (operator)-stratified robustness not testable in this cohort"
- **[C]** If p_compl >= 0.05: "; not robust to protein-completeness stratification
  (p = p_compl)".
- **[D]** Only when step 3 is reached (form 1 of [B]), pass or fail, exactly one of:
  - p2 < 0.05: "; the batch-corrected (operator) discovery set (726) also outranked
    background under the operator-stratified null (AUC y; p = p2)"
  - p2 >= 0.05: "; the batch-corrected (operator) discovery set (726) was not shown to
    outrank background under the operator-stratified null (AUC y; p = p2)"

  When step 3 is not reached, [D] is omitted from the Abstract, and Results report H2
  descriptively as "not tested in the confirmatory sequence".

When p1 >= 0.05, the Abstract carries only [A'], followed by the §5.1 bucket wording.
[A'] reads: "lung (N patients): genes selected by the uncorrected discovery analysis
(2,310) versus background, AUC x, 95% patient-bootstrap interval a-b, one-sided
p = p1".

### 5.3 H2 and the descriptive items

H2 is read as in D1. The descriptive items of §3 never enter the §5 reading.

### 5.4 Bootstrap interval (cannot alter any reading)

- Method: the `c1_run_test.py` bootstrap machinery, with LUAD parameters:
  - exactly 1,000 patient resamples;
  - seed `crc32(b"c1_luad|boot")`;
  - PCA refitted per resample, with K held at the full-cohort value;
  - plain KFold, as in C1-UCEC.
- Report the percentile 95% interval together with the bootstrap mean, SD, bias
  (bootstrap mean minus estimate) and the permutation null mean.
- If p < 0.05 while the interval contains the null mean, Results read: "established
  against the permutation null, not precisely estimated at this n" (the published UCEC
  wording).
- **H1.** Computed always. p = p1, and the null mean is that of the unrestricted null.
  The H1 interval is quoted in the Abstract, next to the UCEC interval (D6c).
- **H2.** Computed only when step 3 is reached. Each resample's K PCs are residualised
  on the frozen operator dummies (§4) of the resampled cases. p = p2, and the null mean
  is that of the within-operator null.

### 5.5 Usable cases and abandonment

- A case is usable if it is on the frozen list and has parsed RNA, a protein aliquot,
  and at least one primary-slide embedding with status ok.
- The run writes a 113-row coverage table.
- Abandonment may be read only if every missing case is attributed in that table either
  to a documented access or corrupt-file failure, or to a documented data property (next
  bullet).
- A case missing because of a documented, outcome-blind property of the data is
  recorded as such, is not "fixed", and counts as a documented loss. These properties
  are:
  - no IDC series for any of its primary slides (C3N-02142, documented here in advance);
  - no `quantDataMatrix` column for its aliquot;
  - no RNA file;
  - zero tissue tiles under the frozen filter.
- Only a failure to process an available, valid file is a pipeline error to fix, not
  an outcome.
- Every attribution is written in the coverage table before the first unblinded output
  (§1.6 step 5), and is not changed after it.

### 5.6 Run guards

- A §5 reading is printed only when B = 1,000 exactly; any other B prints "no
  reading".
- The bootstrap requires 1,000 completed resamples; any skipped resample is FATAL.
- Confirmatory smoke runs use `--blind` (D5g). Nothing is written to disk or printed
  except n, K, gene counts and timing: no AUC, ρ, p, gene-level table or result JSON.
- The time of the first unblinded output is recorded in the result JSON and in the
  post-signature addenda file (§7.4).
- All inputs and code are hash-asserted at start (§7.2); any mismatch is FATAL.

### 5.7 Manuscript consequences, fixed now

These are content rules. Wording that this rule gives in quotation marks for the
manuscript is used verbatim. That includes D1, D3, §3 (b), §5.1, §5.2, §5.4, and items 5
and 6 below. The Abstract's LUAD wording is the §5.2 clauses, with the §5.1 bucket
strings after [A']. The NC abstract limit is 200 words.

0. **Room in the Abstract, made now.**
   - Before signature, and without reference to any LUAD result, the Abstract is cut to
     leave room for the longest §5.2 branch ([A]+[B]+[C]+[D]).
   - The words removed are listed in §7.4 before signature.
   - Room the realised branch does not use is not refilled after the LUAD result.
   - No result clause is cut. Protected: the discovery counts; the transcriptome-baseline
     retention and its caveat; the UCEC reading (organ, n, AUC, interval, p); and the
     external protein-array and whole-proteome nulls.
1. **Abstract.** It reports the primary §5 reading of every C1 cohort that was run:
   organ, n, AUC, bootstrap interval and one-sided p, for UCEC and LUAD in the same
   sentence. A suggestive or not-replicated LUAD reading is never left out while UCEC
   is kept.
2. **"Replicat\*"** is used only for a cohort with p1 < 0.05, always as "set-level",
   and with the §5.2 clauses.
3. **Other items in the Abstract.** H2 enters the Abstract whenever step 3 is reached,
   pass or fail, as clause [D]. The completeness result enters only as clause [C]. The
   bootstrap interval is quoted for both cohorts (D6c). No other descriptive item
   enters the Abstract.
4. **Abandonment.** The Abstract is unchanged; Results and Methods state the
   availability finding.
5. **Results.** The heading is plural only if LUAD p1 < 0.05; otherwise it reads "in
   one of two confirmatory cohorts". LUAD uses UCEC's paragraph structure and table
   columns.
6. **Naming the set.** Every claim names the set: "genes selected by the uncorrected
   discovery analysis (2,310)". The UCEC abstract sentence is clarified the same way
   ("the uncorrected 2,566-gene set").
7. **Methods and Data availability** add PDC000489, IDC `cptac_luad` and the OSF
   registration.
8. **Pipeline wording.** For LUAD, Methods may say "same encoder and tiling". They may
   not say "same pipeline" or "same acquisition" unless they also state the 5 off-grid
   slides and the 12 cases without an SS1553 slide (§1.3).

### 5.8 Reporting commitments

- C1-LUAD is reported whatever the outcome, including abandonment or a D2 NO-GO halt.
- The C1 programme is C1-UCEC (done, 2026-09-16) and C1-LUAD (this rule). No further
  CPTAC-3 confirmatory cohort is attempted before submission. Any later attempt is
  reported whatever its outcome.
- No pooled or meta-analytic combination with C1-UCEC is presented as a replication.
- The UCEC post hoc items of D4 are reported whatever their result.

### 5.9 Code corrections

- **Before the first unblinded output:** a bug found and fixed is recorded in a dated
  entry of the addenda file (§7.4) with the bug, the fix and the new hashes, and is
  third-party-timestamped as in D3.
- **After unblinding:** any change is disclosed in Methods as a post-unblinding change,
  with both results.
  - The §5 reading, the §5.7 wording and the Abstract use the pre-change run, i.e. the
    code hashed in §7.2.
  - The corrected run is also reported, but only if the change corrects a defect that
    made the pre-change run depart from this rule's text.
  - Both readings appear in the Abstract if the two runs fall in different §5 buckets,
    or if they give different forms of §5.2 clause [B], [C] or [D].
  - In that case, "Replicat\*", the plural Results heading, form 1 of [B] and the pass
    form of [D] are used only if both runs support them. Where one is withheld, each
    run's statistic and p-value are printed in its place, labelled by run.

## 6. What this rule does not resolve (disclosed whatever the outcome)

- **Setting.** Same consortium, same collection design, same organ: this is an
  independent patient set, not an independent setting.
- **Accrual population (D5b).** Not adjusted for.

  | Recorded race (GDC) | Confirmatory (n = 113) | Discovery (n = 111) |
  |---|---|---|
  | White | 69 (61%) | 40 (36%) |
  | Asian | 40 | 64 |
  | Black | 4 | 1 |
  | American Indian | 0 | 1 |
  | Unknown | 0 | 5 |

  - C3N/C3L split: 72/41 in confirmatory, 73/34 in discovery. The four Taiwanese
    11LU cases appear in discovery only. All 40 Asian confirmatory patients are C3N.
  - In discovery, operator is aliased with ancestry (main.tex Limitations (ii)).
  - The unrestricted null cannot separate a morphology-protein association from
    accrual structure that recurs in both waves. The operator-stratified null credits
    all between-operator association to chance and tests only within-operator
    association.
- **Operators.** Only 2 of the 32 confirmatory operators appear in discovery, so the
  cross-wave comparison is largely operator-disjoint. Within-wave operator structure
  is handled only by (a)(i) and H2.
- **Scanners.** 9 of the 112 expected analysed cases were scanned only on devices
  absent from LUAD discovery (SS7559, SS7320), and 3 have no recorded scanner (§1.3).
  So the scanner natural experiment (main.tex L270-277) does not cover this wave. The 9
  fall in operator strata with no SS1553 slide, and the 3 unlabelled cases are held
  fixed (§1.3); all 12 are removed by (s2).
- **Protein batch.** PDC metadata do not give case-to-TMT-plex membership or the
  identity of the reference channel: `studyExperimentalDesign` is run-level only, and
  `aliquot_is_ref` is null on all 52 runs. No plex stratum is possible, and within-wave
  plex effects are not addressed.
- **Pre-analytic variables (D6d, disclosure only).** GDC sample-level fields, counted
  per case as a Primary Tumor sample with the field filled
  (`review/C1_LUAD_STEP0B_2026-09-29/step3b_preanalytic_fill_rates_both.json`):

  | Field | Confirmatory | Discovery |
  |---|---|---|
  | Excision-to-freezing time | 113/113 | 107/111 |
  | Clamping-to-freezing time | 62/113 | 71/111 |
  | days_to_collection | 0/113 | 0/111 |

  No pre-analytic covariate enters any analysis. The discovery analysis used none, and
  these times describe the frozen omics aliquot, not the slide.
- **Slide preparation.** Confirmatory is about 100% FFPE; discovery is 126 FFPE and
  15 OCT slides. These are Step 0's 2026-09-27 PathDB figures. Step 0b did not
  re-derive them, and IDC's own embedding-medium and fixative fields are empty for
  every matched slide, so no second metadata source corroborates them.
- **Acquisition grid.**
  - Two slides are at 40x (0.25 µm), outside the discovery range.
  - Three PixelMed-pathway slides are at 0.32 µm and carry no operator label.
  - These are handled only by sensitivity analyses (s) and (s2).
  - Reader-QC arm (ii) is not available for the PixelMed pathway.
- **Step 0b range reads.** Before signature, header range reads of all 135
  confirmatory primary series fetched compressed pixel bytes (preamble). The bytes were
  not decoded, and they carry no omics value.
- **C3N-00545.** It is excluded as a discovery patient even though its confirmatory
  aliquot is distinct. Its discovery aliquot is flagged Disqualified in PDC and has no
  column in the PDC000153 `quantDataMatrix`, yet it entered the discovery analysis
  through the CDAP TSV. The protein-route check characterises this.
- **C3N-02142.** Lost to IDC coverage.
- **Pre-checks.**
  - The within-discovery pre-checks (D2) share TMT batch and collection wave with the
    selection data.
  - UCEC's realised z (2.73) fell below its pre-check median (4.26).
  - The LUAD group construction was reconstructed by inference, so parity with the UCEC
    pre-check groups holds at the schema level only (D2 P0).
  - H2 and the (a)(i) qualifier enter without any pre-check.
- **Other external tests.** No LUAD RPPA or other cross-consortium test is attempted
  here.

## 7. Sign-off and freezing

### 7.1 Knowledge at signing

The signer states, in their own words, what they had seen before signing. At least:

- the C1-UCEC result and its bootstrap interval;
- the discovery LUAD n-fragility and operator-width analyses (the latter with an
  equal-rank noise arm);
- the D2 pre-check outcome;
- the discovery protein-route check;
- the D4 UCEC post hoc results, if they have been run;
- the Step 0 and Step 0b metadata: case lists, slide lists, header-derived batch
  labels and scanners, demographics, pre-analytic fill rates.

They had seen no confirmatory RNA or protein value, and no decoded confirmatory slide
pixel. The undecoded bytes fetched in Step 0b are disclosed in the preamble.

### 7.2 Frozen inputs and code (sha256)

- The code hashed in this table is part of the rule. Where the text leaves a detail open,
  the hashed code fixes it.
- A disagreement between text and hashed code is handled under §5.9: as a
  pre-unblinding addendum if found before the first unblinded output, and as a
  post-unblinding change if found after.
- Local mirrors are listed first. The HPC copies are asserted identical at run time.
- Rows marked "at signing" are filled in before signature.
- Rows marked "before the first P1, P1b or P2 run" or "before its first run on UCEC
  data" are filled in, with the recording time, before that run.
- The D2 builders already produced an earlier `luad_batch_groups.csv` (2026-09-29 07:27,
  without `stream_c3n_vs_c3l`). The hashes recorded are those of the builder version
  that produces the CSV used by P2.

| File | sha256 |
|---|---|
| `server_export/pinned/luad/residual_results_tumoronly.csv` | `63aec1992a77e952b96c4e25e27372e2bee5c52294560c56ef2b229196e87be7` |
| `server_export/results/luad__sitepack_operator.csv` | `aeef7232682c24d209024ba73efea8033cb7486779831bda1a1ba3c8bda45d9f` |
| `pinned/c1_case_ids_luad.txt` | `3629b243a3aa21057dae24e7f5c1b47e8a2571e6ffbcff6f3170333e1ca9fe44` |
| `pinned/c1_excluded_luad.tsv` | `5d3f089597fea2154cac82cf0401959baed10d2f30a40d0f01eda36932fb0a9d` |
| `pinned/c1_primary_slide_ids_luad.txt` (136) | `b7a102b590714a6c6c8739973b8696dd6f7aeb57d8a4860c469177f970a861f3` |
| `pinned/c1_primary_slide_ids_luad_idc.txt` (135) | `97c339d8c93a81b582ff40331ea88c90002c81f0a6ae3ece330eb608fdb62c9f` |
| `pinned/c1_luad_gdc_slides.tsv` | `a8e6a068eaf79b0a6b7aa67ab8debd5450833cdecaffb73629d7f9a055886134` |
| `pinned/c1_luad_batch_labels.tsv` | `0a840533a3601490c360cae8e5fcfc5e6f426e5314e5ac02d1e1b00730f1b3d4` |
| `pinned/luad_batch_groups.csv` (D2; rebuilt with `stream_c3n_vs_c3l`; built 2026-09-29, LF-only) | `57bb7441558cc9600196dd0cd97aaa6268a77a0b2b41ce5c230f7ae1091c1668` |
| `pinned/luad_batch_groups.csv.sha256` (sidecar) | `014d0257668f0b5df7794ffb95fa3dec5432131e5ddae206f9426ed9ea45c2d1` |
| `make_batch_groups.py` (D2 P0 builder) | `aa689821b0fe881e7b3dfb98028ed1f2799b977dc101caa77144df0d37a17125` |
| `build_luad_batch_groups.py` (D2 P0 builder) | `24eb6a061bacd203262883b97e2bef9a065dfaae4fc31763a2ecf0f720c9f981` |
| `pinned/luad_discovery_protein_via_pdc.csv` | `aef3270a117f6189678137e7527e9c945468ab42f6205596d8673bfac651f8ff` |
| `slide_type_map_luad.tsv` (discovery, local mirror) | `3de16eb77f615908c1999806388d9871b48ed8bfc645a6366a7f64b7df9a843e` |
| `aliquot_to_case_tumor_luad.tsv` (discovery, local mirror) | `967cabf18cb6b6c3646aaf05005c46fcdccd7b55eb8dc9882348e84ca4baa2bf` |
| `manifest_rna_tumor_luad.tsv` (discovery, local mirror) | `cacb2858814584a555991829b522807211f0cfe3204c3fc2f37b72e60815e82b` |
| HPC `scripts/pinned/slide_type_map_luad.tsv`, `aliquot_to_case_tumor_luad.tsv`, `manifest_rna_tumor_luad.tsv` (discovery) | confirmed identical to the three local-mirror rows above, verified 2026-09-30 (see 7.4) |
| `residual_analysis.py` | `a54a0a56494d315e15ef075425fe2133ecd9f05193c22656c1145bffdaac4a0a` |
| `split_replication_power.py` (repository top level, not `server_export/scripts/`; HPC `scripts/`) | `79530c075c0f7e75bf06d3c2f1edfe8792cbf4d6c71ae1d32e60bf44f31a0397` |
| `residual_analysis_sitepack.py` (read for reference; its residualisation design, `_fit` and `plain_folds` are reimplemented) | `d0cfeb01ef3d249455d4d26749bbec01b2c6bb3ada6509790f349073b0edbf52` |
| `extract_c1_dicom.py` | `a9a734fa33118c41158aed7a6ff6a7b924d5070e7517ae157afd1b510f2b73c2` |
| `extract_features.py` | `e3d3f95736de21dd612f89c90f36bc35e8c7223a22112d1f38cb3e9cba3018dd` |
| `qc_reader_equivalence.py` | `45aea4d84d768f95870cc937079afee947dac7a6b16ea594422ab12fd1e315dc` |
| `build_slide_map_c1.py` | `7b33174abf03a3d861e3a0664874dc1a811d5e7d7893e3de5ad37352b6239707` |
| `c1_luad_derive_cases.py` | `dcd09d8525debf2e918fd268f0a81a548422f0ed7fd74b41e1f348f5b0a4e34e` |
| `c1_luad_independence_check.py` | `9fdf4c75c6db8bebb334c8e8e8c5d95394481bb9eef65c3edf71ba4ee7d38495` |
| `idc_step3c_derive_batch_labels.py` (review/C1_LUAD_STEP0B_2026-09-29/) | `d0b03afed259451ae0ca8051bcb65761338252da00e9b5bf3c2ff51e4f0f2a5c` |
| `lsf_precheck_luad.sh` (D2; gains the P1b and stream modes) | `f74a1a6ccb59c434588fb821e4d44580f1cc8f16727fe87c42294833dfa9de39` |
| `protein_route_check.py` | `91996eab48331eb235d0407e5e03682b8a80279dababe02df056b6f63c939aa3` |
| `lsf_protein_route.sh` | `ae8f99ec43b829041806e6dca5ece2d2c784c6e12a789f79ed3b2ac4444c8adb` |
| `pinned/slide_type_map_luad_c1.tsv` (built 2026-09-30 14:57 UTC by `build_slide_map_c1.py` from the 136-ID list; GDC slide-entity metadata only) | `d564b6bc5682bc6332aa3ad82387ed300a87d8142bba2d794650281350347d60` |
| `build_slide_map_gdc.py` (imported by `build_slide_map_c1.py` for `fetch_slide_entities`) | `8a14f20b9458c6e15dee508224003d47f45c8455328dbf7e56cd34d72dee5ba8` |
| `c1_run_test_luad.py` | `1d5ee05d57f80af4d90870aecb10e574e4c910a8ce1eab743a067bd9a9278742` |
| `c1_luad_pergene.py` (the per-gene (b1)/(b2) module) | `505f7b1156ccb1dcb026b395d4c32d2aee3263bb1a6246e732300fa045b744f1` |
| `c1_pull_omics_luad.py` | `6c91f5f83bf4ce0541435b0942669c264051fa1f66ee696612cf528d4bc51384` |
| `lsf_extract_c1_luad.sh` | `c7ca2916286b756206a969e0f140a7be6964adfb40109753217c9f380fa05e54` |
| `lsf_c1_luad_test.sh` | `9ac035ecf2688ef10f6e83feb96fd67ae776621270692f9d8f1cd30cfaf8f73f` |
| `c1_ucec_completeness_rerun.py` (D4 item 1) | `c910efdcf080ddaa92b48aee35f33a0739a2cd61bd8ebcb2c72bcf2418960b16` |
| `lsf_c1_ucec_strat.sh` (D4 item 1 wrapper) | `f46eb2ea0b2f322fbce7fd7b4401b5a2c3e34f7d5ac105d184f5ebbc7f88a555` |
| `idc_ucec_derive_batch_labels.py` (D4 item 3, steps 1-2; local, no HPC) | `4ee868ca058c7040ac9ed43cf4c04a28171ea84111aa1347037aec5bd6576403` |
| `c1_ucec_stratified_reread.py` (D4 item 3, steps 3-4) | `b83a9322b91959f392fca0f9329ad582b2a17944c1842010c07cf0f37b0ab124` |
| `lsf_c1_ucec_d4_strat.sh` (D4 item 3 wrapper) | `10405b1ad8401a05c20cea5fd111bd22c9e1a70a31285817b2909a8baed9c169` |
| `c1_ucec_completeness_followup.py` (the two §7.4 follow-ups; imports `auc_strat` from `c1_ucec_completeness_rerun.py` unchanged) | `9c0edd5090810f4b3d842bfa8b2c3c2f7b638c617e135f5a7d0dfc9abd11b9b8` |
| `pinned/c1_ucec_completeness.npz` (D4 item 1 output, copied from the HPC) | `a85814f65a698da81c0e0bee6d294c4f7c2bac9ff9a7be4fb6c3bf3298142aa8` |
| `pinned/c1_ucec_completeness_followup.json` | `f6c1389ca105eaccbba2c508307a636b76a90b185fd1825416db14fdf06d88b3` |
| `c1_ucec_d4_followup.py` (the two §7.4 D4-item-3 follow-ups; imports `assemble` from `c1_run_test.py` unchanged) | `f010c67b2b535e5c2573fb0b1b6be3830ec8e5aab8698bf5528240ef19f98ce7` |
| `pinned/c1_ucec_d4_stratified.npz` (D4 item 3 output, copied from the HPC) | `ef75f1c06d4fc33843a79e9a8d1bb4f4224acd201d30784badb97cf3be13384a` |
| `pinned/c1_ucec_d4_followup.json` | `8f659f67f42b3383a06e539150dbc06674484ffa3545f3d4d34e8adfb6729a21` |
| `pinned/c1_ucec_batch_labels.tsv` (D4 item 3 labels; frozen 2026-09-30 before the re-read) | `7c8ac183151cf6f7a82b85047f2cdf39090ae82a62dad77104607797ae49767a` |
| `review/C1_UCEC_D4_2026-09-29/c1_ucec_d4_summary.json` | `d191f9a1e0c9a16a2f68e6c678a35e0e4e1fcfd2f3d881a7c5386b15ce6403d3` |
| `review/C1_UCEC_D4_2026-09-29/idc_ucec_confirmatory_per_slide_labels.tsv` | `872f2ed31076b95727dec7a4cfb9f06a96c43354d430d8360138517efdce0547` |
| `review/C1_UCEC_D4_2026-09-29/idc_ucec_discovery_validation.tsv` | `9c8b6196399b1066f28f1a655496ffb897a365945d8adac825260f1bd9199ad8` |
| `review/C1_UCEC_D4_2026-09-29/query_log_idc.tsv` | `4015252f6111f51d3889e9b0916852869d01536666c9bb622f68221d7c5259a8` |

### 7.3 Signature

- Signed by: Wei Zhang
- Date and time, with time zone: 2026-10-01 04:16 CST (UTC+8)
- sha256 of this file (recorded outside the file): see `review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md.sha256`, computed on this exact file state, immediately after this signature was saved and before any further edit
- OSF registration DOI:
- GitHub release tag:

### 7.4 Addenda log (pre-signature entries only)

- This file is never edited after its hash is recorded (D3 step 2).
- Post-signature entries go in a separate dated file, `review/C1_LUAD_ADDENDA.md`, each
  entry hashed and timestamped as in D3. These are the §5.6 unblinding time and any
  §5.9 addenda.

Pre-signature entries:

- **D2 input and code hashes, recorded 2026-09-30 (all before their first run):**
  `pinned/luad_batch_groups.csv` sha256 `57bb7441558cc9600196dd0cd97aaa6268a77a0b2b41ce5c230f7ae1091c1668`
  (built 2026-09-29, confirmed unchanged at run time via the sidecar check in
  `lsf_precheck_luad.sh`); `lsf_precheck_luad.sh` sha256
  `f74a1a6ccb59c434588fb821e4d44580f1cc8f16727fe87c42294833dfa9de39`;
  `split_replication_power.py` on the HPC confirmed identical to the local mirror,
  sha256 `79530c075c0f7e75bf06d3c2f1edfe8792cbf4d6c71ae1d32e60bf44f31a0397` (checked by
  `lsf_precheck_luad.sh`'s own unconditional hash guard on every run below, all
  passed). D4 scripts: none has yet been run on the HPC. Their hashes are recorded
  in §7.2 before any first run (item 1: `c1_ucec_completeness_rerun.py` /
  `lsf_c1_ucec_strat.sh`; item 3: `idc_ucec_derive_batch_labels.py` /
  `c1_ucec_stratified_reread.py` / `lsf_c1_ucec_d4_strat.sh`). Correction, same day:
  an earlier version of this sentence said item 1 "was already run ... on
  2026-09-16/2026-09-29". That was wrong. 2026-09-16 is the original C1-UCEC run, and
  item 1 had only been unit-tested locally; no completeness-stratified output exists.
- **D2 pre-check outcome, run 2026-09-30 on the HZAU HPC (`smp` queue, `-n 4
  -R 'span[hosts=1]'`; all nine jobs logged `hash OK` for both
  `split_replication_power.py` and `luad_batch_groups.csv` and `environment OK`
  before computing anything, per `lsf_precheck_luad.sh`'s own guards):**

  | Run | n_sel / n_test | splits with p_auc < 0.05 | median z (AUC) | output |
  |---|---|---|---|---|
  | P1 (`random`, seed 0) | 55 / 50 | **10/10** | 3.41 | `split_replication_power_luad_p1_random.csv` |
  | P1b (`random34`, seed 1) | 71 / 34 | **10/10** | 2.79 | `split_replication_power_luad_p1b_random34.csv` |
  | P2 `q_early_vs_late` | 55 / 49 | 1/1 | 4.35 | `split_replication_power_luad_p2_q_early_vs_late.csv` |
  | P2 `q_2017Q3_vs_rest` | 67 / 37 | 1/1 | 4.18 | `split_replication_power_luad_p2_q_2017Q3_vs_rest.csv` |
  | P2 `op_balanced_0` | 52 / 52 | 1/1 | 3.22 | `split_replication_power_luad_p2_op_balanced_0.csv` |
  | P2 `op_balanced_1` | 52 / 52 | 1/1 | 1.81 | `split_replication_power_luad_p2_op_balanced_1.csv` |
  | P2 `op_balanced_2` | 52 / 52 | 1/1 | 3.96 | `split_replication_power_luad_p2_op_balanced_2.csv` |
  | P2 `op_balanced_3` | 52 / 52 | 1/1 | 4.33 | `split_replication_power_luad_p2_op_balanced_3.csv` |
  | P2 `op_balanced_4` | 52 / 52 | 1/1 | 2.61 | `split_replication_power_luad_p2_op_balanced_4.csv` |
  | P2 `stream_c3n_vs_c3l` (descriptive; D6d) | 71 / 34 | 0/1 (z=0.99) | -- | `split_replication_power_luad_p2_stream_c3n_vs_c3l.csv` |

  Reading, applying section 0 D2's literals exactly: P1 = 10/10 >= 8/10; the seven
  UCEC-schema P2 columns = 7/7 >= 6/7. **GO.** Proceed to signature (D2
  Consequences). No P1/P1b/P2 setting was changed after seeing any of these numbers.
  Descriptive reading of `stream_c3n_vs_c3l` (D2 Scope, fixed before this run): its
  own split did not reach p_auc < 0.05 while its matched-n reference P1b passed
  (10/10 >= 8/10) -- per the pre-declared rule this is read as evidence of accrual
  dependence (not as under-powered noise), consistent with the recorded-ancestry
  difference between the two waves already disclosed in section 6.
- **Discovery protein-route check outcome, run 2026-09-30 on the HZAU HPC (`smp`
  queue, `-n 1`; discovery-only, PDC000153, never PDC000489 -- see
  `protein_route_check.py`'s own `FORBIDDEN_STUDY` assertion):** comparing route
  "pdc" (PDC GraphQL `quantDataMatrix(pdc_study_id="PDC000153")`, the same query
  method `c1_pull_omics_luad.py` uses for the confirmatory PDC000489 pull) against
  route "cdap" (`residual_analysis.load_protein_matrix()` under `source
  morpho_env.sh luad`, the CDAP TSV route the LUAD discovery pipeline actually used
  to produce `server_export/results/luad__residual_results_tumoronly.csv`), on the
  10,724 genes actually tested in the discovery run: n_pdc_cases/n_cdap_cases/
  n_shared_cases = 110/109/108; n_pdc_genes/n_cdap_genes = 11,029/11,032;
  gene_jaccard = 0.99973; tested genes missing from pdc / from cdap = 0 / 0 (all
  10,724); NaN-pattern agreement rate = 1.0 (0 nan-only-in-pdc, 0 nan-only-in-cdap,
  over 1,191,132 cells); median / p99.9 absolute delta over 1,104,387 shared valid
  cells = 2.499e-05 / 4.995e-05; per-gene Pearson r median / p5 (11,024 genes with
  r) = 0.999999997 / 0.999999869; pdc_exact_zeros / cdap_exact_zeros /
  cdap_non_numeric_sentinel_cells = 218 / 53 / 0. Reading: the two routes are
  equivalent for practical purposes -- near-total gene-set overlap, zero missing
  tested genes on either side, perfect NaN-position agreement, and per-gene
  correlations indistinguishable from 1 (deltas at floating-point/rounding scale
  relative to log2-ratio units). The only asymmetry (exact-zero count, 218 vs 53)
  is small relative to the 1.1M-cell comparison and does not move any of the
  agreement statistics above. This validates the PDC quantDataMatrix route -- the
  method `c1_pull_omics_luad.py` will use on the confirmatory cohort -- against the
  route the discovery pipeline actually ran, before that method is ever pointed at
  confirmatory data.
- **HPC hashes of the discovery pinned files (§7.2), verified 2026-09-30:** HPC
  `scripts/pinned/slide_type_map_luad.tsv` = `3de16eb77f615908c1999806388d9871b48
  ed8bfc645a6366a7f64b7df9a843e`, `aliquot_to_case_tumor_luad.tsv` =
  `967cabf18cb6b6c3646aaf05005c46fcdccd7b55eb8dc9882348e84ca4baa2bf`,
  `manifest_rna_tumor_luad.tsv` = `cacb2858814584a555991829b522807211f0cfe3204c3
  fc2f37b72e60815e82b` -- all three byte-identical to the local-mirror hashes
  recorded in 7.2.
- **Abstract words removed under §5.7 item 0, done 2026-10-01, with no LUAD result in
  existence and none consulted.**

  *Counting convention, fixed before any cutting and used for every figure below:*
  render the LaTeX (`\noindent` counted as one token, other commands dropped, `\&`->&,
  `\%`->%, each `$...$` collapsed to one space-free token), split on whitespace, and
  count an en-dash numeric range as two words. Under it the Abstract stood at **195**
  and the NC limit is 200.

  *Reservation:* the longest §5.2 branch [A]+[B] form 1+[C]+[D], typeset in the
  Abstract's own style and rendered with placeholder values, is **55** words including
  the joining semicolon. §5.7 item 6's mandated renaming of the UCEC set
  ("discovery-selected genes" -> "the uncorrected 2,566-gene set") costs a further
  **+2**, and was applied now rather than reserved. Target therefore 200 - 55 = 145,
  cut to **144** so that the worst case lands at 199/200 rather than 200/200. Net cut
  53 words. Room the realised branch does not use is not refilled (§5.7 item 0).

  *What was removed, in document order.* No result clause was cut; every protected item
  of §5.7 item 0 survives with its numbers unchanged, verified clause by clause.
  1. The whole opening framing sentence, 23 words: "A gene's own mRNA explains part of
     its protein abundance; we asked whether H\&E morphology predicts the residual that
     one RNA-seq measurement leaves." Its content survives compressed as "predicted
     protein abundance beyond the gene's own mRNA", which is the manuscript's title
     phrase.
  2. From the methods sentence: "we embedded whole-slide images with a pathology
     foundation model" and "measuring each gene's out-of-fold incremental $R^2$ over an
     mRNA-only baseline" (19 words) compress to "H\&E foundation-model embeddings" plus
     the surviving "incremental $R^2$" inside the parenthesis; "under permutation
     false-discovery control" (4) becomes "permutation FDR" (2); "99--137 patients each"
     becomes "$n=99$--137" (-1); "from one proteogenomic consortium" becomes "CPTAC"
     (-3).
  3. In the transcriptome caveat: "because each transcript is measured once" -> "with
     each transcript measured once" (-1); "part of what morphology recovers may be" ->
     "part may be" (-3); "rather than" -> ", not" (-1).
  4. In the UCEC sentence: "In a later same-consortium uterine cohort (138 patients),"
     -> "Later, within CPTAC, uterus (138 patients):" (-1), which also puts the sentence
     in the "organ (n):" frame §5.2 [A] requires, so the LUAD clause appends in parallel.
  5. The external nulls became their own sentence: "; beyond the consortium, four
     protein-array cohorts" -> ". Outside CPTAC, four protein-array cohorts" (-2), and
     "and a whole-proteome cohort none" -> "a whole-proteome cohort none" (-1). This is
     what frees the semicolon slot: the UCEC sentence now ends at "a set-level
     replication."
  6. Closing: "H\&E is therefore an accessible" -> "H\&E is an accessible" (-1); "a
     signal partly shared" -> "partly shared" (-1); "regulation that transcriptomics"
     -> "regulation transcriptomics" (-1).

  *How it was chosen.* Four independent compressions were drafted under different
  strategies (cut the methods hardest / cut the closing hardest / rewrite wholesale for
  density / minimal-edit surgery) and judged by two independent reviewers, one reading
  as an NC editor and one hunting for overclaim and for quietly altered protected
  numbers. Both recomputed every word count themselves, both ranked the minimal-edit
  version first, and both proposed the same two corrections, which are applied above:
  "pre-correction" -> "pre-batch-correction" (it sat beside "permutation FDR" and could
  be misread as pre-FDR) and "beyond a gene's own mRNA" -> "beyond the gene's own mRNA"
  (the title's wording). Drafts materially under 144 were penalised, since unused room
  is forfeited.

  *One decision left open for the corresponding author, flagged by both judges.* The
  chosen text uses "CPTAC" three times, and the manuscript never expands it anywhere
  (53 uses in main.tex, 24 in the supplement, zero expansions). Measured alternatives:
  glossing it once in the Abstract costs +4, and removing it entirely in favour of "one
  proteogenomic consortium" / "in the same consortium" / "Beyond the consortium" costs
  +6. At 144 there is no slack, so either would have to be paid by deleting
  "$n=99$--137" (2) and then 2 or 4 more from the organ list or "permutation FDR". The
  recommendation recorded here is to keep "CPTAC" in the Abstract and instead expand it
  at first use in the Introduction, which costs the Abstract nothing and fixes a defect
  the manuscript already has. If the author decides otherwise before signature, the
  change is made then, not after any LUAD result.
- **D4 item 3, steps 1-2 (labels and qualification), run 2026-09-30 locally
  (`idc_ucec_derive_batch_labels.py`, collection `cptac_ucec`, no `--limit`):**
  - Access. 582 HTTP range reads, all status 206; 20 escalated past 256 KiB. About
    200 MB was fetched in total, parsed in memory and discarded. None of it was
    decoded or written to disk. As in Step 0b, a fixed-length range can run past the
    header into PixelData. For UCEC this carries no blinding consequence, because its
    confirmatory data is published.
  - Step 1, validation. On 389 UCEC discovery Primary Tumor slides with both an
    ImageComments value and a `slide_acquisition_meta.tsv` row, `User =` matched
    `aperio.User` on 389/389 and `Date =` matched `aperio.Date` on 389/389.
    VALIDATED.
  - Step 2, labels. 167 of the 169 analysed confirmatory slides carry the field. The
    two without it are both slides of C3L-01064, the one unlabelled case. So 137 of
    138 cases are labelled (99.3%). An independent recomputation of every case's
    modal label, with ties to the first sorted slide ID, reproduced the labels file
    with 0 mismatches. The labels file is frozen at the §7.2 hash, before the
    re-read.
  - D1 criteria, operator axis, on the 138-case analysed set. 45 levels; 19 strata
    with >= 3 cases; largest stratum 14 cases (10.2% of labelled); 85.5% of cases
    permutable. All four criteria hold. **QUALIFIES.** The re-read (steps 3-4)
    therefore runs.
  - The all-zero UUID. It labels 8 UCEC cases and is used as one stratum, taken
    literally, as D5a does for LUAD. D4 item 3's text does not address this value.
    Treating it as unlabelled instead would not change the verdict: coverage 93.5%,
    18 strata >= 3, largest 10.9%, 79.7% permutable.
  - Scan-quarter axis, descriptive only (D4 item 3 never falls back to it): 10
    levels, 137/138 labelled, largest 2017Q3 (35.8%), all four criteria hold.
- **D4 item 1 (completeness-stratified rerun), 2026-09-30.** HPC smoke (`--perms 5`)
  passed. It reproduced n = 138, 169 slides, 10,153 genes and K = 28, and matched
  auc_obs 0.597443 to the published JSON within 1e-9. The observed completeness-
  stratified AUC combines the complete stratum (2,406 x 6,068, AUC 0.5598) and the
  incomplete stratum (130 x 1,542, AUC 0.6450), weighted by pair count, to 0.561, the
  value D4 states. The smoke run printed no p, by design.
- **D4 item 1, full run (`--perms 1000`), finished 2026-09-30 23:41 CST, 77.1 min,
  rc = 0.** Full validation passed: all 15 statistics and the code and discovery
  hashes match the published `c1_ucec_result.json`, floats to 1e-9 and counts
  exactly. It reproduced auc_obs 0.597443 and p_auc 0.003996, the published values.
  The completeness strata, split at the discovery file's modal per-gene n = 100:

  | Stratum | selected | background | frac. positive | AUC | permutation null |
  |---|---|---|---|---|---|
  | complete | 2,406 | 6,068 | 0.0785 | 0.5598 | 0.4654 +/- 0.0375 (z = 2.52) |
  | incomplete | 130 | 1,542 | 0.0173 | 0.6450 | 0.5557 +/- 0.0366 (z = 2.44) |

  Combined: **AUC_strat = 0.560924**, which is the 0.561 D4 pre-declared from the
  published gene table, against a stratified null of 0.466610 +/- 0.037404, one-sided
  **p = 0.006993** (7/1001, well clear of the 1/1001 floor).

  **Reading (rewritten 2026-10-01 after an adversarial check; the first version of
  this paragraph was wrong in three ways, each recorded below rather than silently
  replaced. Every figure here was independently recomputed from
  `server_export/pinned/ucec_c1/c1_ucec_gene_level.csv` joined to the discovery
  `n` column, which reproduces the published pooled AUC to 0.597442556490091 and both
  stratum AUCs exactly.)**

  *What the run establishes.* Restricted to the 8,474 genes quantified in all 100
  discovery patients -- 98.6% of the stratified statistic's gene pairs -- selected
  genes still outrank background: 0.5598 against a permutation null of 0.4654 +/-
  0.0375 (z = 2.52). The pair-weighted statistic is 0.560924 against 0.466610 +/-
  0.037404, one-sided p = 0.006993. The set-level ordering survives holding discovery
  completeness fixed.

  *What it does not establish.*
  - **(a) It is not two strata corroborating each other.** Pair weights are 0.98646
    complete and 0.01354 incomplete, so AUC_strat - AUC_complete = +0.0012: the
    statistic essentially *is* the complete stratum. The two nulls are about 0.82
    correlated (solving Var(AUC_strat) with the quoted SDs), because both come from
    the same 1,000 permutations. The withdrawn first version said "both strata point
    the same way at similar strength, so neither carries the result alone." That was
    wrong.
  - **(b) Completeness does explain the fall in the AUC's level.** The 0.597 -> 0.561
    fall is exactly the removal of the 23.31% of gene pairs that compare a complete
    with an incomplete gene; those pairs sit at an observed 0.7176. The four-block
    decomposition of the pooled AUC is 0.5598 / 0.7884 / 0.3843 / 0.6450 at weights
    0.75650 / 0.19224 / 0.04087 / 0.01039, summing to 0.597443. What survives
    stratification is the excess over the null, 0.0992 -> 0.0943 (-5%). The withdrawn
    version called it an error to read the drop as "completeness explained part of the
    signal"; that was too strong. It explains the level, not the excess.
  - **(c) The stated reason for the null displacement was wrong.** Inside the complete
    stratum the discovery n equals 100 for all 8,474 genes by the stratum's own
    definition, so per-gene n cannot displace that stratum's null to 0.4654. The
    driver is a selection-induced nuisance unrelated to completeness: selected genes
    rank *below* background on mRNA-only predictability (AUC(`r2_rna`) = 0.3998 within
    the complete stratum) and `r2_rna` predicts the confirmatory increment (Spearman
    +0.36 among complete-stratum background genes).

  *The structural point that matters most, and that the first reading missed.*
  Per-gene missingness is **invariant** under this permutation.
  `split_replication_power._gene_test` computes each gene's observed-patient mask
  `v = ~isnan(protein)` once, then permutes only the morphology rows (`pcs[p][v]`);
  the gene's own n, its mRNA-only baseline `r2r` and its capacity penalty are
  identical in the observed column and in all 1,000 null columns. Completeness was
  therefore structurally incapable of producing the pooled p. The null had already
  priced it, and the decomposition shows exactly where: the cross-completeness block's
  implied null is 0.602 against the complete block's 0.465, and they cancel to the
  published pooled null of 0.498. So D4 item 1 **confirms that the null was doing its
  job**; it does not remove a confound that was operating. (This also means the
  pooled null landing at 0.498 ~ 0.5 is cancellation, not design.)

  *Scope and caveats, all to be carried into any Results sentence.*
  - The strata are **discovery** completeness. The script computes and stores the
    confirmatory per-gene n but never stratifies on it, and the "incomplete" bin spans
    n = 30-99 in one level, where the gradient is undiminished: Spearman(discovery n,
    confirmatory increment) among background genes is +0.3564 inside that bin versus
    +0.3492 across all background genes. Nothing was adjusted within it.
  - Inside the incomplete stratum, discovery n alone separates selected from
    background at AUC 0.6965, above the statistic's own 0.6450; and that stratum's
    z falls to 1.89 if the top 7 of its 130 selected genes are dropped.
  - p = 0.0040 versus 0.0070 is 4 versus 7 exceedances in 1,001; the Monte Carlo SE at
    this level is about 0.002, so the difference carries no information. No interval
    was computed for the change in z (2.73 -> 2.52), so "nearly unchanged" should not
    be asserted without one.
  - No bootstrap interval exists for AUC_strat. The published pooled interval is
    0.48-0.65, wider than the entire above-null excess.
  - 93.2% of confirmatory increments are negative in both strata, so the AUC ranks
    degrees of harm from adding morphology, as the published text already says.
  - D1's "co-reported" bullet asked for the within-stratum pooled rho and the fraction
    of positive increments *per stratum for selected and background separately*. The
    script produced neither: its `frac_positive` is stratum-wide. Recomputed locally
    (observed values only, no nulls): frac positive selected/background = 0.1018
    (245/2,406) and 0.0692 (420/6,068) complete, 0.0615 (8/130) and 0.0136 (21/1,542)
    incomplete; within-stratum rho = 0.1876 complete and 0.3235 incomplete, against the
    published pooled 0.29357.
  - "Post hoc" applies to the null and p only. The point estimate 0.5609235 was already
    in `review/C1_LUAD_STEP0B_2026-09-29/local_checks/auc_strat_check_output.json`
    (2026-09-29) and D4 quotes 0.561 in advance, so there is no selective reporting.
  - The claim must not be allowed to widen into "not explained by technical artefact."
    A per-gene covariate is absorbed by this null for free; patient-level structure
    (site, plex, operator) is not, and that is D4 item 3's job.

  D4 requires this to be reported in Results whatever the result. **Done 2026-10-01**,
  together with D4 item 3: a paragraph in §`sec:c1ucecres`, a Methods paragraph in
  §`sec:c1ucecmeth`, and a new Supplementary Note `snote:c1posthoc` carrying the detail
  and every limit recorded in this section. Both documents recompile with 0 undefined
  references and 0 errors (main 39 pages, supplement 48). Every number in the new text
  was checked against `pinned/c1_ucec_completeness_followup.json` and
  `pinned/c1_ucec_d4_followup.json`.
- **Two follow-ups to D4 item 1, specified 2026-10-01 BEFORE running them, and
  reported whatever they show.** They are post hoc and beyond D4's text, which
  specified only the discovery-completeness rerun. They are recorded here in advance
  so that neither the analysis nor the reporting can be chosen after seeing the
  answer. They add no permutation: both are pure arithmetic on the already-saved
  `c1_ucec_completeness.npz` (`genes`, `sel_mask`, `strata`, `M` = 10,146 x 1,001,
  `confirm_n`), reusing the frozen `auc_strat` from `c1_ucec_completeness_rerun.py`
  unchanged. They cannot change any published UCEC reading (D4's closing sentence),
  and they create no new claim: their only job is to bound the two disclosure gaps
  named above.
  1. **Confirmatory-completeness stratification.** The rule fixed in advance, by exact
     analogy with the discovery rule (which splits at the discovery cohort size, 100):
     a gene is "complete" if `confirm_n` equals the confirmatory cohort size, 138, and
     "incomplete" otherwise. Report (a) the stratified AUC and its permutation p under
     that split, (b) the same under the joint 2x2 of discovery x confirmatory
     completeness, and (c) for reference, a split at every distinct value of the
     stratifier. Any cell with zero selected or zero background genes drops out, as
     `auc_strat`'s own formula already implies.
  2. **A paired null for the level drop.** `AUC_pooled` and `AUC_strat` are computed
     on the same 1,001 columns, so their difference D has its own permutation
     distribution. Report D observed, the mean and SD of D under the null, the
     standardised position of D observed within it, and the two-sided empirical p.
     The pre-stated interpretation: if D observed sits inside the null spread of D,
     the 0.597 -> 0.561 fall is the structural removal of cross-completeness pairs and
     nothing more; if D observed sits above it, the pooled statistic lost more than the
     null did, and the fall is not purely structural. Both outcomes are reported.
  **Result, run 2026-10-01 locally on a downloaded copy of the npz
  (`c1_ucec_completeness_followup.py`; npz sha256
  `a85814f65a698da81c0e0bee6d294c4f7c2bac9ff9a7be4fb6c3bf3298142aa8`, 10,146 genes x
  1,001 columns). The script first reproduced both reference statistics from that same
  M -- pooled 0.597443 / p 0.003996 and discovery-stratified 0.560924 / p 0.006993 --
  so it is operating on the same numbers the HPC run reported.**

  *1. Confirmatory completeness.* `confirm_n` spans 33-138 over 77 distinct values,
  with 8,761 of 10,146 genes at 138.

  | Stratification | AUC | null | z | one-sided p |
  |---|---|---|---|---|
  | none (pooled, published) | 0.597443 | 0.498195 +/- 0.036296 | 2.73 | 0.003996 |
  | discovery n (D4 item 1) | 0.560924 | 0.466610 +/- 0.037404 | 2.52 | 0.006993 |
  | **confirmatory n** | 0.559666 | 0.462130 +/- 0.037792 | **2.58** | 0.006993 |
  | **both, 2x2** | 0.554170 | 0.461248 +/- 0.037573 | **2.47** | 0.009990 |
  | confirmatory n, all 77 levels | 0.559210 | 0.461515 +/- 0.037940 | 2.57 | 0.006993 |

  The gap named above is closed: the ordering survives stratifying on completeness as
  measured in the cohort the AUC is actually computed in, on both completeness
  variables jointly, and on the finest available split of the confirmatory one. The
  same weighting caveat as before applies and must be stated rather than repeated as
  an error: the confirmatory-complete stratum carries 0.99314 of the pair weight, and
  the joint statistic is 0.99283 the complete|complete cell (2,381 selected vs 5,852
  background, AUC 0.553821, z 2.47). So the defensible sentence is "among genes
  quantified in every patient of both cohorts, selected genes still outrank
  background," not "four cells agree." That said, all four joint cells do point the
  same way -- complete|complete z 2.47, complete|incomplete z 1.62 (25 vs 216),
  incomplete|complete z 2.18 (73 vs 455), incomplete|incomplete z 2.18 (57 vs 1,087)
  -- which is descriptive support, not independent corroboration, since all four share
  the same 1,000 permutations.

  *2. Paired null for the level drop.* D = AUC_pooled - AUC_strat(discovery n):

  - D observed **+0.036519**
  - D under the null **+0.031585 +/- 0.003593**, 95% of the null in [+0.024642, +0.038293]
  - D observed sits at **z = +1.37** inside that null; two-sided empirical **p = 0.177**

  By the reading fixed in advance, this is the "inside the null spread" branch: the
  0.597 -> 0.561 fall is the structural removal of cross-completeness pairs, and the
  pooled statistic did not lose materially more than the null did. Equivalently, the
  excess over each statistic's own null falls +0.099248 -> +0.094314, i.e. by 0.004934
  or 5.0%, and 0.004934 is 1.37 SD of what the permutations themselves produce for
  that same difference. So "nearly unchanged" is now supported, with the honest
  qualifier that the point estimate is a 5% loss and that this comparison has limited
  power: D observed is inside the null but in its upper half, and the 95% null interval
  reaches +0.0383 against the observed +0.0365.

  Resolution caveat, unchanged: p values of 0.006993 and 0.009990 are 7/1001 and
  10/1001 against the pooled 4/1001, and the Monte Carlo SE at this level is about
  0.002, so the differences among them carry no information.
- **D4 item 3, steps 3-4 (the within-operator re-read), 2026-09-30.** Before
  submission the four uploaded files were confirmed on the HPC at their §7.2 hashes.
  The HPC smoke (`--perms 5`) passed. It reproduced the same n, K, gene counts and
  auc_obs 0.597443 as item 1, and read validated_operator = True and qualifies =
  True. By design it computed no stratified AUC. The labels-hash guard runs only on
  the full path, so the smoke did not exercise it.
- **D4 item 3, steps 3-4, full run (`--perms 1000`), finished 2026-10-01 01:27 CST,
  153.5 min, rc = 0.** Full validation passed against the published
  `c1_ucec_result.json`, and the labels-hash guard fired on a real run for the first
  time and passed (`hash OK`, `7c8ac183151c...`). 137/138 cases labelled; the one
  unlabelled case held fixed. Operator diagnostics reproduce the local D1 report
  exactly: 45 levels, 19 strata with >= 3 cases, largest 14/137 = 10.2%, 85.5%
  permutable.

  | Arm | AUC | null | one-sided p | exceedances |
  |---|---|---|---|---|
  | unrestricted (published) | 0.597443 | 0.498195 +/- 0.036296 | 0.003996 | 3/1000 |
  | **primary, rule-literal: same PCs, within-operator null** | 0.597443 | 0.505936 +/- 0.031649 | **0.003996** | 3/1000 |
  | descriptive extra: operator-residualised PCs, same null | 0.567199 | 0.488699 +/- 0.030618 | **0.005994** | 5/1000 |

  The primary arm's observed AUC is bit-identical to the pooled one, as the script's
  own assertion requires: only the null changed.

  **Reading. This entry is written after an adversarial check, and states only what
  survived it. Every quantitative claim below was recomputed locally from the frozen
  labels or by simulation.**

  *What the run establishes.* Under a null that confines every permutation inside an
  operator block, the set-level ordering still beats the null at p = 0.003996. An
  ordering that were itself a between-operator effect would drive the restricted null
  up towards the observed 0.597, because the restricted null preserves exactly the
  between-operator association; it sits at 0.506. **So the confirmatory ordering is
  not a between-operator effect.** Residualising the morphology PCs on operator
  dummies as well, which is not required by D4's text, leaves AUC 0.567 against
  0.489 at p = 0.005994.

  *What must not be claimed, with the reasons.*
  - **(a) The size of the null shift attributes nothing.** The shift is +0.007741,
    7.8% of the unrestricted excess, with a Monte Carlo SE of 0.0015 (so 5-11%). It is
    NOT a bound on how much of the ordering is operator-driven, because confining
    permutations to small blocks raises the null mechanically whatever the blocks
    mean: with 45 blocks the expected number of cases keeping their own morphology in
    a draw is one per block, 46/138 = 33.3%, so real signal leaks into the null by
    construction. In a matched simulation on the exact realised 45-level size profile
    with the project's own `cv_r2` and AUC, but with **operator-free random blocks**,
    the null rose by 18.9%, 19.4% and 31.5% of the excess across three seeds -- all
    larger than the observed 7.8%. A label-permuted, size-matched partition null on
    the real data was not run; until it is, the shift's magnitude is uninterpretable.
    An earlier draft of this entry read "raises the null only from 0.498 to 0.506 --
    7.8% of the excess -- so at most a small part of the ordering is attributable to
    operator." That inference is withdrawn.
  - **(b) No operator-adjusted effect size exists for the primary arm.** 0.597443 is
    unchanged by construction; only the reference moved. §4 states what this null
    tests: it credits all between-operator association to chance and tests only
    within-operator association. The arm licenses "survives an operator-stratified
    null", not "survives operator adjustment".
  - **(c) The three z values, 2.73 / 2.89 / 2.56, are not comparable and are not
    reported as evidence.** They are ratios against three different nulls; the
    exceedance counts are 3, 3 and 5 of 1000, and at this resolution the Monte Carlo
    SE is about 0.002, so neither the identical p nor the 0.004-versus-0.006
    difference carries information. This repeats the error already recorded and
    corrected for D4 item 1 in this same section.
  - **(d) "Residualised against operator" overstates the coverage.** Recomputed from
    the frozen labels under the script's own rule: the 12 dummies are the 12 largest
    levels, sizes 14, 14, 8, 7, 7, 6, 5, 5, 5, 4, 4, 4, covering 83 of 138 cases
    (60.1%). **55 cases (39.9%) in 33 levels receive no column**, including seven
    3-case levels that pass MIN_GROUP and are dropped only by the cap. The binding
    lever here is MAX_LEVELS, not MIN_GROUP, which is inert because every kept level
    has at least four cases -- the opposite of the discovery-side finding, which
    therefore does not transfer. PC variance retained is 66.4%, and the PCA is fit
    *before* residualisation, so the 28 retained directions were chosen to maximise
    variance that includes operator variance.
  - **(e) The cost of residualisation is 14.2%, not 21%.** Arms 2 and 3 share the same
    1,000 permutation draws, so the like-for-like contrast is arm 3 against arm 2:
    excess 0.078500 against 0.091507, i.e. 85.8% retained. Comparing arm 3 to the
    unrestricted excess charges the restriction penalty to residualisation as well.
    From the reported means alone, E[D_null] = 0.017237 against D_obs = 0.030244, so
    roughly 57% of the observed drop is what permuted pairings produce anyway; the
    proper paired figure is computable from the saved npz and is *pending*.
  - **(f) A rank statistic cannot answer the discovery-side operator finding.** That
    finding is a collapse in the *level and sign* of per-gene increments: on the real
    sweep (`review/recalc/operator_width_summary.csv`) UCEC discovery at the 12:3
    design goes from a median increment of +0.1183 to -0.1035, -87.5% of reference,
    while the same-rank noise arm keeps +0.0575, with 90.4% of genes below their own
    noise arm; UCEC is the worst of the five cohorts at that design. The confirmatory
    statistic is a Mann-Whitney AUC over ranks and 93.2% of confirmatory increments
    are already negative, so it is invariant to any monotone common shrinkage. This
    re-read bears on the ordering only, and must never be cited as answering the
    discovery-side sensitivity. (An earlier note in this file cited "increments at 33%
    of baseline against 71% for noise" as the discovery-side figure. That was the
    LUAD eight-gene smoke test, not the full sweep, and is withdrawn.)
  - **(g) Only one statistic of a two-statistic family was re-read.** The published
    result reports rho = 0.293568 at p = 0.004995 alongside the AUC.
    `c1_ucec_stratified_reread.py` computes `rho_op` (L535) and `rho_resid` (L565) and
    stores neither. Pending from the saved npz. Until then the correct wording is that
    the AUC arm survives, not that the replication survives.
  - **(h) Nothing external validates the new numbers.** The validation block checks the
    reproduction of the published unrestricted result. 0.505936, 0.031649, 0.567199,
    0.488699 and 0.030618 exist only in this run's own output.

  *A defect in the operator axis itself, found during this check and not previously
  recorded.* Of the 45 levels, **44 sit entirely inside one scan quarter; the sole
  exception is the all-zero UUID**, whose 8 cases span 2018Q4, 2019Q1 and 2019Q3. It
  also appears on **0 of the 393** discovery Primary Tumor slides, so the 389/389
  `User =` validation never covered it, and D5a's LUAD rationale for treating it
  literally was never checked for UCEC. Consequences: in the primary arm those 8 cases
  are permuted as one exchangeable block, so genuine batch structure spanning three
  quarters is destroyed rather than held fixed, which is anti-conservative on the axis
  under test; in the residualised arm it consumes one of the 12 columns to subtract a
  single shared mean from a three-quarter mixture, and displaces a real 3-case
  operator from the design. The qualification sensitivity is on record (treating it as
  unlabelled still qualifies), but the AUC was never recomputed that way, and D4 item
  3's text does not address the value. Not done; disclosed.

  *Two follow-ups to D4 item 3, specified 2026-10-01 BEFORE running them and reported
  whatever they show.* They close gaps (e) and (g) above and add no permutation: both
  are arithmetic on the saved `c1_ucec_d4_stratified.npz`, which holds all three arms'
  full 1,001-column matrices, reusing `assemble` from `c1_run_test.py` unchanged so the
  statistic cannot drift. As with the D4 item 1 follow-ups, they are post hoc and
  beyond D4's text, they cannot change any published UCEC reading, and they create no
  new claim.
  1. **The missing rho arm.** Report the observed rho, its null mean and SD and its
     one-sided p, for all three arms, using the discovery increments as `xs` exactly as
     the published run did. The published anchor to reproduce first is rho = 0.293568
     at p = 0.004995.
  2. **Paired nulls for both transitions.** Because all three arms share the same
     column index, D can be formed column by column. Report, for
     D_restrict = AUC(unrestricted) - AUC(within-operator) and
     D_resid = AUC(within-operator) - AUC(residualised), and for the same two
     differences in rho: D observed, the null mean and SD of D, D's standardised
     position in its own null, and the two-sided empirical p. Pre-stated reading: D
     observed inside the null spread means that transition moved the statistic no more
     than the permutations themselves move it; D observed above the null spread means
     the observed statistic lost more than the null did.

  **Result, run 2026-10-01 locally (`c1_ucec_d4_followup.py`; npz sha256
  `ef75f1c06d4fc33843a79e9a8d1bb4f4224acd201d30784badb97cf3be13384a`). The script
  refuses to report anything unless the unrestricted arm first reproduces BOTH published
  anchors; it reproduced AUC 0.5974425564900907 and rho 0.2935675639399924 exactly.**

  *1. The rho arm.* rho behaves like the AUC in every arm, so the published
  two-statistic family is now re-read in full:

  | Arm | rho | null | one-sided p | exceedances |
  |---|---|---|---|---|
  | unrestricted (published) | 0.293568 | 0.158599 +/- 0.053622 | 0.004995 | 4/1000 |
  | within-operator | 0.293568 | 0.167052 +/- 0.046313 | 0.003996 | 3/1000 |
  | residualised | 0.256896 | 0.143112 +/- 0.046151 | 0.007992 | 7/1000 |

  Limit (g) above is therefore lifted: both halves of the family survive both arms, and
  the wording may be "the replication survives an operator-stratified null" rather than
  "the AUC arm survives".

  *2. Paired nulls.*
  - **Restriction.** D observed is exactly 0 for both statistics, by construction. The
    paired null of that difference is -0.007742 +/- **0.049659** for the AUC and
    -0.008452 +/- 0.073114 for rho. So the draw-to-draw spread of the difference is
    6.4 times its own mean. This is a second, independent reason the +0.0077 null shift
    is not interpretable as an attribution, alongside the matched-block simulation in
    (a): the shift is small even relative to the noise of the comparison that produced
    it.
  - **Residualisation.** AUC: D observed +0.030243 against a null of +0.017237 +/-
    0.022753, z = +0.57, two-sided p = 0.545, with the null reproducing **57%** of the
    observed drop. rho: D observed +0.036671 against +0.023939 +/- 0.034207, z = +0.37,
    p = 0.711, null reproducing **65%**. By the reading fixed in advance, both sit
    inside the null spread: **residualising the morphology PCs on the 12 largest
    operator levels did not cost the statistic more than projecting out thirteen
    columns costs anyway.** That is the honest form of limit (e) -- the like-for-like
    retention is 85.8% of the within-operator excess, and most of even that 14.2% is
    generic projection cost rather than operator removal.

  *What is deliberately NOT done, and why.* Two checks the adversarial review asked for
  need a fresh per-gene pass and cannot be derived from any saved matrix: a
  label-permuted, size-matched partition null (which would calibrate the null shift),
  and an all-zero-UUID sensitivity (holding those 8 cases fixed instead of permuting
  them as a block). Each is roughly another 2.5 HPC hours per seed. The first is not
  required by anything claimed here, because the surviving claim -- that the ordering is
  not a between-operator effect -- rests on the restricted null sitting at 0.506 rather
  than near 0.597, not on the size of the shift, which is withdrawn as uninterpretable
  in (a). The second remains an open disclosure.

  *Two further scope limits, both verified locally.* In the confirmatory cohort the
  operator axis does not hold the slide-count confound fixed: eta-squared of slide
  count on operator is 0.312 against a label-permuted floor of 0.324 +/- 0.094, i.e.
  nothing beyond chance, so this re-read must not be cited as controlling the
  slide-count pathway that main.tex Limitations (ii) rests the UCEC operator concern
  on. And only 4 of the 45 confirmatory operators appear among the 27 discovery
  operators, covering 23 of 138 cases (16.7%), so the axis tested is largely
  wave-specific; a between-operator null is in any case blind to a scanning artefact
  shared across operators or drifting within one.
- **Confirmatory slide map (§1.3), built 2026-09-30 14:57 UTC locally, under the
  author's authorisation of the same day:** `build_slide_map_c1.py
  pinned/c1_primary_slide_ids_luad.txt pinned/slide_type_map_luad_c1.tsv`. This queried
  GDC case, sample and slide-entity metadata only; no pixel, RNA or protein value was
  touched. Input and builder hashes matched §7.2 before the run. Result: 136 of 136
  slide IDs have a GDC slide entity, all `Primary Tumor`, covering exactly the 113
  frozen cases (93 with one slide, 17 with two, 3 with three), as stated in §1.2. It
  agrees with Step 0's `pinned/c1_luad_gdc_slides.tsv` on every slide's sample type
  and case (0 disagreements). The output is LF-only; its hash is in §7.2.
- Code freeze. As of 2026-09-30 no §7.2 row is left "at signing".
- **Code-freeze hash re-verification, run 2026-10-01 locally, immediately before
  signature.** Every row of §7.2 (52 rows) was re-hashed from the file on disk and
  diffed against its recorded value. **52/52 match, 0 mismatches.** This is the final
  scripted check before §7.3 is signed; no file in §7.2 has drifted since its hash was
  recorded.

## Deviations from the C1-UCEC rule

| Item | C1-UCEC | C1-LUAD |
|---|---|---|
| Independence | no overlap existed | 7 discovery patients excluded, patient-level (§1.2) |
| Slide headers | read after the rule, before any statistic (UCEC Addendum 2) | range reads of all 135 confirmatory primary series before signature (Step 0b, authorised 2026-09-29 as header-only). They overran into undecoded pixel bytes. This went against the review verifiers' advice (batch-2, stat-5) and is disclosed in the preamble. |
| Reader QC | one pathway, arms (i)+(ii) | per pathway; PixelMed arm (ii) unavailable (§1.4) |
| Batch handling | unrestricted null only; a post hoc stratified re-read is now added (D4 item 3) | three-step fixed sequence H1 -> operator-stratified qualifier -> H2, plus a completeness co-report (D1, D6a) |
| Per-gene FDR | listed as a forbidden substitute | descriptive only (b); cannot substitute for §5 |
| Manuscript consequences | UCEC-only narrative | fixed clauses plus content rules covering both cohorts (§5.2, §5.7) |
| PathDB-pooled slide sensitivity | pre-declared, never run (D4) | dropped (D5f) |
| Guards | count asserts, B < 1,000 guard | hash asserts, fail-closed independence, coverage table, `--blind` smoke, exact-B guard, blinding order |
| Timestamp | none | OSF + GitHub release (D3) |
| Unchanged | | AUC/ρ statistics, 0.05/0.15/45 thresholds, K formula, B = 1,000, MIN_N = 30, one-sided |

## Appendix: implementation before signature (hash each into §7.2)

1. **`c1_run_test_luad.py`**, a new file.
   - **Imports** `residual_analysis` and `split_replication_power`. It does not import
     `residual_analysis_sitepack`: the §4 residualisation and within-operator
     permutation are reimplemented, with the frozen dummy list and MAX_LEVELS,
     MIN_GROUP and MIN_PERMUTABLE copied as literals.
   - **Paths.** All paths are required argparse arguments, with no defaults. The signed
     rule's path is one of them.
   - **Hashes.** It asserts the §7.2 hashes. The result JSON records the sha256 of the
     rule file, of the script itself and of every imported module, next to the input
     hashes.
   - **Checks.** It runs the fail-closed independence check against the union of the
     three HPC discovery pinned files, and writes the coverage table.
   - **Computes** H1, (a)(i), H2 (steps 1-3 in sequence), the completeness-stratified
     AUC, (a)(ii), (b), (s) and (s2).
   - **Saves** M (genes x (1+B)) and the per-gene confirmatory n as `.npz`, with the
     npz sha256 in the result JSON.
   - **Guards.** `--blind` and the exact-B guard.
   - **Bootstrap** mode follows §5.4.
   - **Output text.** The clause strings are copied verbatim from §5.2 and the bucket
     strings from §5.1. Seeds are as in §3 and §4.
2. **`c1_pull_omics_luad.py`**, a new file.
   - The crosswalk is restricted to the Primary Tumor aliquots of the frozen 113 cases.
   - FATALs as in §1.5.
   - Writes `protein_luad_c1.csv`, the RNA manifest and the coverage table.
3. **Per-gene (b1)/(b2) module.**
   - (b1) feeds the confirmatory CSV to `residual_analysis._process_gene` / `bh_fdr`
     unchanged.
   - (b2) reimplements sitepack's `_fit` and `plain_folds` exactly as stated in §3 (b2),
     without importing `residual_analysis_sitepack`, on PCs residualised with the frozen
     dummy set.
4. **`pinned/slide_type_map_luad_c1.tsv`**, built by `build_slide_map_c1.py` from
   `pinned/c1_primary_slide_ids_luad.txt`.
5. **Wrappers `lsf_extract_c1_luad.sh` and `lsf_c1_luad_test.sh`.**
   - Explicit paths. `unset` lists each variable by name, since bash `unset` does not
     accept a glob: `unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK
     MORPHO_SLIDE_MAP MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS
     MORPHO_SIG_CSV MORPHO_BATCH_AXIS MORPHO_MAX_LEVELS MORPHO_N_PERM`.
   - The extraction wrapper hard-codes `pinned/c1_primary_slide_ids_luad_idc.txt`.
   - No argument pass-through that could override frozen parameters.
   - `chmod +x` after upload.
6. **D2 files.**
   - `build_luad_batch_groups.py` adds `stream_c3n_vs_c3l`; the CSV and its sidecar are
     rebuilt and re-hashed.
   - `lsf_precheck_luad.sh` gains a `random34` mode (P1b:
     `--n-test 34 --splits 10 --perms 200 --top-frac 0.15 --seed 1`, one job) and
     accepts `stream_c3n_vs_c3l` in group mode.
7. **D4 item 3.** A UCEC operator-label extraction plus validation script, and a
   stratified re-read script built on `c1_ucec_completeness_rerun.py`'s exact
   reproduction of the published run.
8. **Abstract cut** under §5.7 item 0, listed in §7.4.
9. **Errata** for the Step 0 and Step 0b notes: `review/C1_LUAD_STEP0_ERRATA_2026-09-29.md`.
10. **Smoke tests.** Every new script is smoke-tested through its wrapper on discovery
    LUAD data or synthetic data. Confirmatory data are never used for a smoke test,
    except with `--blind`.
11. **Executed UCEC scripts.** `c1_run_test.py`, `c1_pull_omics.py` and `lsf_c1_test.sh`
    are not edited.
