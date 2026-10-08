# Post hoc BRCA external replication on frozen tissue sections, with a set-level readout: pre-specification (file date 2026-10-08, local time UTC+8)

**Dates and times.** The file name carries the local date (UTC+8), 2026-10-08. Every time
in this file is UTC and is given with its UTC date. Drafting began on 2026-10-07 UTC at
about 21:00Z; the first draft is dated below, and the revision described under Status was
made on 2026-10-08 UTC after 00:20Z. Wherever the paper text below says "specified on
[date]", [date] is the UTC date of the `# written` line of this file's sidecar, written
as for example "8 October 2026 (UTC)".

**Status.** This file was first drafted at 2026-10-07T22:41Z. That draft is kept unchanged as
`POSTHOC_BRCA_TS_PRESPEC_2026-10-08_pre-review.md`.
- *First review.* The same night (report received 2026-10-07T23:08Z) an adversarial review
  found nine must-fix items (M1-M9) and ten nice-to-have items (N1-N10). They are applied
  throughout; the labels M1-M9 and N1-N10 in parentheses refer to that review.
- *Second review.* A cross-file review of this file and of the two other post hoc
  pre-specifications of 2026-10-08 (the C1-UCEC within-plex re-read and the
  measurement-error checks) found nine further must-fix items for this file (B1-B9) and
  wording and file-list items. They are applied in this version, before freezing.
- *Third review.* On 2026-10-08 from about 00:50Z a final check of the three files
  made the Abstract policy paragraph byte-identical across them, unified the rounding
  sentence, the failed-run paragraph and the "specified on [date] in a dated file
  released publicly before it ran" wording, and added the minimum-group wording of
  section 3. No script changed.

The file is frozen by the sha256 of this file, every script and every input listed in
section 4, all in the sidecar `POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256`. The freeze
time is the UTC time that `make_setlevel_sidecar.py` writes into the sidecar header.
Everything is released publicly on GitHub before any residual-analysis job (`RUNBOOK.sh`
step 6 onward) is submitted, under the rule below. On 2026-10-07 UTC the corresponding
author asked for this analysis.

**What is released (the same rule in the three post hoc pre-specifications of file date
2026-10-08).** This file, its sidecar and every file the sidecar lists are released on
GitHub before the run, except derived data files (`.npz` and `.tgz` files and the per-gene
result tables named below) and third-party source tables under `figures/figdata/`. These
follow the repository's data policy: each is fixed by its sha256 in the sidecar, the
derived data files are deposited with the other derived data at acceptance, and the
third-party tables are public at their sources. The release commit message lists every
withheld path.

Here, nothing is withheld: the sidecar lists no `.npz` or `.tgz` file, no published
per-gene result table and no third-party table. The five candidate-set files and
`setlevel/universe_common.tsv` hold, per gene, the selection flags and the discovery score
(the mean of the published per-gene increments over the cohorts that tested the gene) that
the set-level readout uses, and `brca_tested_genes.{tsv,txt}` hold the BRCA-tested symbols
(the .tsv with each protein's non-missing count); all are released so that the readout can
be checked. The matrices and result files that the runs will
write, and the server's `dx_arm_inputs_outputs.sha256`, do not exist at freezing and are
not listed.

The published BRCA arm used the CPTAC-2 BRCA iTRAQ4 proteome (PDC000173), GDC STAR-Counts
RNA and TCGA Diagnostic Slides (DX), with Phikon at stride 0, in 102 patients. It reported
0 of 10,031 genes and 99.1% negative increments (main text, external replication).

The 2026-09-11 plan (item C3 of `REMAINING_WORK_2026-09-11.md`) asked for a TS arm, and for a
dated note, written before the DX result, saying that both arms would be reported. That
note was never written: the plan contains only the to-do item. The DX result is now
published. The commitment to report "whatever the result" therefore binds only the TS arm
and the set-level readouts defined here.

**What existed when this file was written, and what will exist when it is frozen (review M9).**

- *Local, before writing (all 2026-10-07 UTC):*
  - the GDC/PDC inventory and the TS dry-run selections (21:03-21:08Z);
  - the vial match (21:03Z);
  - the discovery selection-set inventory and the candidate set files (21:29Z);
  - the HGNC drift pre-check (21:31Z);
  - the scripts (21:06-21:37Z).
- *Server, before freezing (cross-review B1).*
  - Five of the seven deployed scripts (`cptac2_replicate.py`, `extract_phikon.py`,
    `ts_slide_qc.py`, `brca_tested_genes.py`, `lsf_cptac2_ts_extract.sh`) are deployed to
    `$ZW/scripts_ts_arm/` and are byte-identical to the sidecar copies.
  - The v2 `cptac2_residual.py` and `lsf_cptac2_ts_residual.sh` deployed with them have not
    run, are not kept locally, and are replaced by v3 at step 5b.
  - Steps 0, 0c and 1 of `RUNBOOK.sh` have run. Step 1 printed 189/232 TS files kept for
    102/102 cases, 59.3 GB. Step 0c printed `tested (n>=30, in RNA): 10031`, the published
    number.
  - Step 2 (download) was started before freezing. Nothing in this file depends on whether
    it has finished.
  - These steps ran on 2026-10-07 UTC after 21:37Z. Their exact times are in the server
    logs and shell history, which are brought back at step 8.
- *Step 0c output:* `brca_tested_genes.tsv` and `brca_tested_genes.txt` (BRCA-tested gene
  symbols with non-missing counts) arrived in the local folder at 2026-10-07T22:59Z.
  - Nobody has opened either file beyond the count line printed on the server.
  - They are hashed into the sidecar exactly as they arrived.
  - The overlap between the selection sets and these symbols has not been computed. It is
    computed only after freezing (section 2, Symbols), by a read-only command.
- *Did not exist at writing, and will not exist or will not have been read at freezing:*
  - any TS image locally;
  - any TS embedding (steps 3-4);
  - any step-5 slide-QC output (magnification, mpp, tile counts);
  - any residual matrix;
  - any set-level statistic for either arm.

  The freeze precedes reading any step-5 output. If freezing happens only after step-5
  output has been created, that output stays unopened until the sidecar exists.

**The blindness statement for steps 0-5, with its one exception.** Steps 0-5 compute no
statistic relating morphology to protein or mRNA. The one exception is step 0c, which reads
per-protein non-missing counts and RNA symbol labels only.
- Step 0 hashes the published inputs and results without opening them.
- Step 1 reads GDC/PDC metadata and the case column of `rna_case_file.tsv`.
- Step 2 downloads and md5-checks.
- Steps 3-4 process images only.
- Step 5 reads slide headers, and the shapes and coordinates of the embedding tensors.

## 0. What it asks, scope and prior exposure

- **Material.**
  - The DX images are FFPE diagnostic sections. All 102 DX slides carry vial code 01Z.
  - TCGA's frozen Tissue Slides (TS) are sections of frozen tumour tissue: top (TS), bottom
    (BS) or middle (MS) sections of the specimen.
  - The proteome aliquot's vial carries a primary TS-arm slide for 100 of 102 patients (the
    match is on the vial, `vial_match.py`).
  - The portion codes differ (slides 01-03, aliquot 11-61), so it is not established that
    the slides are adjacent to the analysed aliquot.
- **What it asks.**
  - Can a signal be recovered from frozen sections of the same tumours (the proteome
    aliquot's vial in 100 of 102 patients), where the FFPE diagnostic slides showed none?
  - A TS null cannot settle whether section type or distance contributes. Frozen sections
    carry their own artefacts (ice crystals, folds, thicker sections), their magnification
    may differ from DX, and n is again 102.
- **Why a set-level readout.** It is added to both arms because a per-gene FDR count at
  n = 102 is a weak test.
- **Scope.** The analysis is post hoc and is reported whatever the result. No published
  number changes.
- **Prior exposure (disclosed; times UTC).**
  1. The DX per-gene summary is public (0/10,031; 99.1% with a negative increment). It was
     known when the sets were chosen.
  2. The 170-gene cross-organ core and the UCEC and LUAD set-level readings are public (AUC
     0.597 and 0.644).
  3. On 2026-10-07 UTC (before 20:33Z) the same analysis lineage saw the selected sets of the
     other four discovery organs on the UCEC C1 confirmatory matrix: unrestricted AUC
     0.564-0.585 (`POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md`, section 0, item 2). Before
     the BRCA set was chosen, it was therefore known that cross-organ sets transfer within
     CPTAC.
  4. No BRCA morphology-protein statistic beyond the published summary has been seen.
     - The selection inventory (`brca_selection_set_inventory.py`) reads only the five
       discovery tables. Its BRCA-overlap columns are `PENDING` and `brca_gene_file` is
       null.
     - The choice of k = 3 used discovery counts only.
     - The published DX per-gene table has not been read gene by gene locally. On
       2026-10-07 at about 23:10Z, no BRCA CPTAC-2 per-gene result file was found in
       `MorphoResidual_paper/`: by file name, by content of its csv/tsv/json files, in its
       zip/tar archives, and in `source_data/*.xlsx`. None was found by file name elsewhere
       under `D:/claude/pajinsen/` either. The table reached the paper only as its summary.
       This check covers local files only.
     - No DX set-level statistic has ever been computed: the published job did not save its
       permutation draws.
  5. `brca_tested_genes.{tsv,txt}` are local but unopened (Status above).

## 1. Data and arms

- **Cases.** The 102 patients of the published arm: the intersection of PDC000173 tumour
  aliquots, GDC Diagnostic Slides and GDC STAR-Counts RNA. They were re-derived on
  2026-10-07 UTC, and the scripts assert equality with the published case list. The RNA and
  protein inputs are the published ones, unchanged. This includes the two patients whose
  metastatic-sample RNA is averaged in, as published.
- **TS slides: the frozen list (review M6).**
  - *The list.* The 189 files with `keep = 1` in
    `BRCA_TS_ARM_2026-10-08/dryrun/primary_all/ts_slide_selection.tsv`. They are identical
    to the files in `dryrun/primary_all/manifest_brca_ts_slides.txt`. The sha256 of the
    file names is `9576e733...e049`, where the hashed text is the 189 names, one per line,
    CR stripped, in `LC_ALL=C sort` order, each line newline-terminated (cross-review).
  - *Content.* 95 top, 86 bottom and 8 middle sections. Every GDC 'Tissue Slide' file with a
    primary-tumour sample code (01A/01B) for these patients, covering 102 of 102 patients.
    Open access. 59.3 GB, where GB = 10^9 bytes (59.35 x 10^9 bytes).
  - *Excluded.* Metastatic (06A, 2 files) and normal (11A, 41 files) slides.
  - *Server re-query.* The re-query at step 1 reproduced 189 files, 102 cases and 59.3 GB.
  - *GDC changes.* If GDC adds or removes files, the difference is reported and the
    analysis uses the intersection of the frozen list with the GDC manifest. Nothing is
    added.
  - *No substitution.* A slide whose download fails after retries, or that is skipped at
    extraction ("no tissue tiles", or still failing after resubmission), is not replaced.
    Every such slide is listed in the supplement with its reason (`ts_skipped_slides.txt`).
  - *Tile count.* No minimum tile count per slide: a slide with at least one tissue tile is
    used.
  - *Missing cases.* A patient without a usable TS slide leaves the TS arms only; the DX arm
    keeps 102. The placeholder [n] below is the actual number of TS patients.
    - n >= 95 (the coverage of the primary_ts_only rule): the TS arms are reported in the
      main text as in section 3.
    - n < 95: the TS arms go to the supplement only (section 3).
  - *Enforcement.* The residual wrapper (`lsf_cptac2_ts_residual.sh` v3) refuses to run the
    TS arm unless the embeddings equal exactly the frozen list (intersected with the GDC
    manifest) minus the recorded skips. It passes no argument to the dx job and only
    `--pool tile` to the ts job; any other argument stops it. `cptac2_residual.py`
    (patch v3) stops if any kept embedding is not a TS/BS/MS slide, or if the number kept
    differs from that number. `local_brca_setlevel.py` checks that ts and ts_tile hold the
    same cases, a subset of the 102 DX cases.
  - **Not done.** None of these is computed, now or later, for this paper, even though the
    embeddings make them possible after the fact:
    - the primary_ts_only and ms_vial selection rules;
    - any top-versus-bottom (TS/BS/MS) split;
    - any per-slide-type or per-magnification analysis.
- **Extraction and pooling.**
  - Phikon, 256-px tiles at level 0, stride 0, with the published tissue filter
    (`extract_phikon.py` with `--manifest`).
  - Per patient: the mean over tiles within a slide, then the equal-weight mean over the
    patient's slides, as in the published DX arm (`--pool slide`).
- **Magnification (N4).**
  - At level 0, a 256-px tile covers about 64 um at 40x (about 0.25 um/px) and about 128 um
    at 20x (about 0.5 um/px). That is a two-fold difference in field width, four-fold in
    area.
  - The TS and DX mpp distributions (`ts_slide_qc.py`, step 5) go to the supplement with
    this explanation.
  - Whatever they show, nothing is re-run, restricted or stratified by magnification.
- **Arms.**
  - **TS (primary).**
  - **DX (reproduction).**
    - *What it is.* The published DX arm re-run with the patched script, on the default
      code path, reading the DX embeddings in `$ROOT/emb_phikon`.
    - *Required server checks (review M7, RUNBOOK step 5b).* Before step 6, check that this
      folder is the one the published job read, and that its tiles are at stride 0 (the
      Methods), read from the saved tile coordinates (x-step 256 px; stride 512 gives
      512). The published copies of `cptac2_residual.py` and `extract_phikon.py` on the
      server must hash to the local published copies. If anything differs, a dated
      amendment to this file is written and released before step 6, for example to set
      `DX_EMB_DIR` for the dx job or to correct the Methods.
    - *Gate (review M3 and cross-review B2), applied locally by `local_brca_setlevel.py`.*
      It compares the re-run with the published per-gene table
      (`results/brca_protein_results.csv`) and the published summary
      (`results/summary_brca.csv`). Both published files must first pass the published-file
      check, which the script performs:
      - each hashes to its line (the path ending `results/brca_protein_results.csv` or
        `results/summary_brca.csv`) in `dx_arm_inputs_outputs.sha256`, written at step 0
        and brought back from the server;
      - the published summary reads n_cases 102, tested 10,031, significant 0, and
        negative_frac_all rounding to 0.991.

      If either published file fails this check, the gate is **fail**. Otherwise:
      - **exact**: the same gene set and max |published - re-run increment| <= 1e-12, and
        the summaries agree: n_cases, tested and significant identical, negative_frac_all
        identical to 4 decimals (0.01 percentage points, the precision of the Methods'
        99.11% vs 99.12%);
      - **host**: the same gene set, max |difference| <= 1e-4 (the full-SVD PCA agrees
        across hosts only to about 3e-5), and the same summary identity;
      - **fail**: anything else.

      Under exact or host, the DX arm supplies the DX segment of section 3. Under fail,
      the DX set-level rows are labelled "DX re-run (not identical to published)", go to
      the supplement only, and the DX segment reads "not reproduced" (section 3). The TS
      arms are reported as usual in every case. The gate tier is written into every DX
      output row and into the json.
  - **TS-tile (sensitivity).** TS with tile-weighted pooling (`--pool tile`), the CPTAC-3
    discovery convention.
- **Estimator and null.** Unchanged from the published arm (`cptac2_residual.py`):
  - 20 PCs (full SVD, all patients);
  - 5-fold `KFold(shuffle, random_state=0)`;
  - closed-form ridge with lambda = 1;
  - B = 1,000 permutations with `RandomState(i)`;
  - BH-FDR; significant = FDR < 0.05 and increment > 0.

  The patched script also saves the genes x 1,001 increment matrix (column 0 observed), so
  the set-level readout uses the same permutations. All three residual jobs run on one
  pinned host (`bsub -m`), which is recorded (N6) in `prerun_checks.txt` by
  `echo "HOST=$HOST"` at step 6.
- **Slide QC.** `ts_slide_qc.py` records microns per pixel, apparent magnification, scanner
  and tile counts per slide. They are reported descriptively, with no stratification.

## 2. Readouts

- **Per gene (as published):**
  - n;
  - genes tested;
  - significant;
  - percentage;
  - median increment of the significant genes;
  - fraction of tested genes with a negative increment.
- **Set level (new).** The universe is the BRCA-tested genes among the 8,309 genes tested
  in all five discovery cohorts (review M1); every frozen set file has exactly these 8,309
  rows.
  - *Statistic.* The AUC of the BRCA increment of a discovery gene set against the other
    BRCA-tested genes of this universe (rank Mann-Whitney, the CPTAC-3 confirmatory
    statistic), with a one-sided p against the same 1,000 permutations:
    p = (#null >= observed + 1)/1,001.
  - *Descriptive.* Spearman rho between the discovery score and the BRCA increment, and the
    null mean and SD.
  - *Minimum group.* A set with fewer than 10 selected genes, or fewer than 10 other
    (background) genes, among the BRCA-tested genes is untestable (`MIN_GROUP = 10` in the
    frozen script). This is not expected to occur. A primary set that were untestable for
    this reason is reported with the coverage-failure wording of section 3 (the frozen script
    records it as untestable and does not stop).
  - **Primary set:** `pub_k3`, the genes significant (FDR < 0.05, positive increment) in at
    least three of the five published CPTAC-3 discovery cohorts (777 genes). Breast is a new
    organ, so a cross-organ set is used.
    - Among k = 1..5, k >= 3 has the largest absolute excess over the independence
      expectation: 777 observed versus 484.6 expected.
    - The ratio is largest at k = 5 (5.4x; k = 4, 3.39x).
    - The 170-gene core (k >= 4) is small.
  - **Secondary sets** (Holm within each arm, descriptive):
    - `pub_k4` (the 170-gene core);
    - `pub_k2`;
    - `op_k2` (operator batch-corrected, >= 2 cohorts);
    - `phi_k3` (Phikon discovery, >= 3 cohorts).

    The sets are nested (`pub_k4` within `pub_k3` within `pub_k2`) and positively
    dependent, so Holm is conservative here (N9).
  - **Symbols:** exact matching between the CPTAC-3 discovery tables and the BRCA tables;
    no alias remapping. Each set's coverage among the BRCA-tested genes is reported.
    - *When.* It is computed after freezing and before step 6 (so that no residual job is
      wasted), by the read-only command in `RUNBOOK.sh` that reads the frozen set files and
      `brca_tested_genes.txt` and writes nothing. `brca_selection_set_inventory.py
      --brca-genes` is not used: it rewrites the frozen set and inventory files.
      `local_brca_setlevel.py` recomputes the same coverage from the matrices.
    - *Rules.* The computation changes nothing except through these rules:
      - If fewer than 50% of the **primary** set's genes are BRCA-tested, the computation
        stops. Any remedy is a dated amendment made before any set-level statistic is read.
      - A **secondary** set below 50% is reported as "untestable (coverage)", and Holm runs
        over the remaining secondary sets (N5).

## 3. Readings and wording (tiers: p < 0.05 shown; 0.05 <= p < 0.15 suggestive; p >= 0.15 not shown)

**Rounding.** Readings use unrounded values; if rounding would print a value on the other
side of a threshold, decimals are added until it does not (e.g. p = 50/1,001 printed as
0.04995; z = 1.96 not "2.0"). Here the thresholds are 0.05 and 0.15: p = 50/1,001 is printed
as 0.04995 (not 0.050 or 0.0500), and p = 150/1,001 as 0.1499 (not 0.150).

**Notation.**
- p_ts, a_ts: TS primary arm, primary set.
- p_tile: TS-tile, primary set.
- p_dx, a_dx: DX, primary set.
- [s] of [t], [f]: significant and tested genes and the percentage with a negative
  increment on TS.
- [k]: the number of the 777 genes that were tested.
- [n]: the number of TS patients.
- [m] +/- [sd], [m_dx] +/- [sd_dx]: the null mean and SD of the AUC, on TS and on DX.
- [date]: the UTC date of the sidecar's `# written` line.

**Results.** The text below follows the existing sentence "On iTRAQ, BRCA showed no
morphology-predictable proteins: 0 of 10,031 tested, 99.1% with a negative increment." It
is one sentence made of a TS segment and a DX segment (review M2). It is followed by this
sentence: "This frozen-section and set-level analysis was specified on [date] in a dated file
released publicly before it ran (Supplementary Note~\ref{snote:brcats})."

*TS segment, n >= 95, primary-set coverage at least 50%.* Exactly one of the following:
- p_ts < 0.05: "Post hoc, on frozen sections from the same tumours ([n] of the 102
  patients; the proteome aliquot's own vial in 100 of 102), [s] of [t] genes reached
  FDR < 0.05 ([f]% with a negative increment), and the 777 genes significant in at least
  three CPTAC-3 cohorts ([k] tested) ranked above the other BRCA-tested genes of the
  discovery universe (AUC [a_ts]; null [m] +/- [sd], p = [p_ts])"
- 0.05 <= p_ts < 0.15: the same, with "ranked above" replaced by "were suggestively ranked
  above".
- p_ts >= 0.15: the same opening, then "... did not rank above the other BRCA-tested genes
  of the discovery universe (AUC [a_ts]; null [m] +/- [sd], p = [p_ts])". The pre-review
  "therefore not explained by the diagnostic slides' distance" clause is deleted.

If s > 0, "[s] of [t] genes reached FDR < 0.05" is followed by "(versus 0 on the diagnostic
slides)".

*TS segment, n < 95 (cross-review B3).* "Post hoc, frozen sections were usable for only [n]
of the 102 patients, fewer than the 95 fixed in advance, so the frozen-section readouts are
given in Supplementary Table [x] only".

*TS segment, primary-set coverage below 50% with no amendment.* "Post hoc, on frozen
sections from the same tumours ([n] of the 102 patients; the proteome aliquot's own vial in
100 of 102), [s] of [t] genes reached FDR < 0.05 ([f]% with a negative increment); the
set-level comparison could not be made because only [c]% of the 777 genes were among the
BRCA-tested genes". In this case the DX segment is omitted. If the primary set is
untestable because of the minimum group (section 2) rather than coverage, "because only
[c]% of the 777 genes were among the BRCA-tested genes" is replaced by "because too few
genes were on one side of the comparison to test it".

*Both n < 95 and coverage below 50% (cross-review B3).* The two variants are combined:
"Post hoc, frozen sections were usable for only [n] of the 102 patients, fewer than the 95
fixed in advance, so the frozen-section readouts are given in Supplementary Table [x] only;
the set-level comparison could not be made on either slide type because only [c]% of the
777 genes were among the BRCA-tested genes". The DX segment is omitted.

*DX segment, appended after a semicolon (not in the coverage variants).* When the TS
segment is the n >= 95 variant, "the same set" refers to its set. When the TS segment is the
n < 95 variant, "the same set" is replaced by "the 777 genes significant in at least three
CPTAC-3 cohorts ([k] tested)" and "background" by "the other BRCA-tested genes of the
discovery universe". "Also" is used only when p_ts < 0.05 and n >= 95. Exactly one of:
- Gate exact or host, p_dx < 0.05: "on the diagnostic slides the same set [also] ranked
  above background (AUC [a_dx]; null [m_dx] +/- [sd_dx], p = [p_dx])."
- Gate exact or host, 0.05 <= p_dx < 0.15: "on the diagnostic slides the same set was
  suggestively ranked above background (AUC [a_dx]; null [m_dx] +/- [sd_dx], p = [p_dx])."
- Gate exact or host, p_dx >= 0.15: "on the diagnostic slides the same set gave AUC [a_dx]
  (null [m_dx] +/- [sd_dx], p = [p_dx])."
- Gate fail: "the diagnostic-slide arm was not reproduced exactly by the re-run (largest
  difference [d]), so its set-level reading is given in Supplementary Table [x] only."

*Optional closing sentence.* It is allowed only when every one of these holds:
- s = 0;
- n >= 95;
- the gate is exact or host;
- p_ts, p_tile and p_dx are all >= 0.15.

The sentence: "Frozen sections of the same tumours did not recover a signal either; because
frozen sections carry their own artefacts, this does not exclude a contribution of section
type."

*Supplement only.*
- The TS-tile arm and every secondary set never enter the main text.
- If TS-tile and TS primary disagree in tier (one p < 0.05 and the other >= 0.15), the
  supplement says so in one sentence. The main text does not change.

**Abstract (cross-review B4, B6; orchestrator decision D1).**

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

Here: the condition is that any set-level reading in the supplement table (any arm, any
set, raw one-sided p) has p < 0.05. Such a reading would leave the existing sentence "a
whole-proteome cohort none" overstating the BRCA null; the reading itself is not stated in
the Abstract, and no set-level result is. Otherwise the Abstract is unchanged.

**Consequences of the Abstract qualifier (B4).** If, and only if, the condition above holds,
the following five places in `main.tex` as of 2026-10-07 UTC gain "at the gene level". The
Abstract itself is line 63.
1. The Fig. 5d caption (`fig:master2`), near line 556: "... CPTAC-2 TCGA-BRCA cohort were
   tested the same way and are equally null at the gene level".
2. Results, near lines 868-869: "... and is not established at the gene level, at the
   sample sizes tested ($n=85$--451), across consortium and platform boundaries".
3. Discussion, near line 1122: "BRCA's whole-proteome null at the gene level is the
   sterner test".
4. Results, near line 861: "the whole-proteome BRCA null at the gene level is not subject
   to this limit".
5. Discussion, near lines 1010-1012: "one independent whole-proteome mass-spectrometry
   cohort found no morphology-predictable protein at the gene level" (replacing "at all").

Any further Discussion text is not pre-specified. It may not describe the post hoc
readings more strongly than the Results sentence.

**Methods, data availability and Supplement (cross-review B5).** These are unconditional:
the frozen-section arm is reported under any result.
- *Methods:* a new paragraph in `main.tex` after the sentence ending "with no panel
  base-rate comparison" (near line 1739), reading: "A post hoc arm on TCGA-BRCA frozen
  tissue slides, specified on [date] in a dated file released publicly before it
  ran (Supplementary Note~\ref{snote:brcats}), repeated the BRCA analysis on the 189
  primary-tumour frozen tissue slides (95 top, 86 bottom and 8 middle sections) of the same
  102 patients, a list fixed before any image was processed. Tiles were encoded with Phikon
  at stride 0, pooled by slide and then by patient as in the diagnostic-slide arm, and the
  estimand, cross-validation, ridge penalty, PCA, 1,000 permutations and BH-FDR were
  unchanged. The set-level readout was the AUC of the 777 genes significant in at least
  three of the five CPTAC-3 discovery cohorts against the other BRCA-tested genes among the
  8,309 genes tested in all five, with a one-sided p from the same permutations. The
  diagnostic-slide arm was re-run through the same code and had to reproduce the published
  result before its set-level reading was used."
- *Data availability:* near line 1906, "TCGA-BRCA diagnostic slides" becomes "TCGA-BRCA
  diagnostic and frozen tissue slides".
- *Supplement:* a new Supplementary Note (`snote:brcats`), placed after Supplementary
  Note~`snote:c1posthoc` in `supplement.tex`, holds one table with:
  - every arm (TS, DX, TS-tile) x every set (primary and secondary, raw and Holm p,
    coverage, untestable notes);
  - the per-gene readouts;
  - the DX gate (tier, maximum difference, published-file check, summary comparison);
  - the case counts and the TS patients missing relative to DX;
  - the skip list with reasons;
  - any GDC-list difference;
  - the slide-QC summary with the TS and DX mpp distributions and the field-of-view note.

  Wherever the paper counts its post hoc analyses (for example the sentence beginning "Four
  post hoc checks" near `supplement.tex` line 495), the count is updated by the editor, not
  by this file.

## 4. Freeze and run order

1. **Sidecar and release.**
   - Command: `python review/BRCA_TS_ARM_2026-10-08/make_setlevel_sidecar.py --prespec
     review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md`, run from `MorphoResidual_paper/`.
   - It writes the sha256 sidecar, with the UTC time in its header (review M4), of:
     - this file and its `_pre-review` copy;
     - `local_brca_setlevel.py`;
     - the five set files `setlevel/candidate_sets/{pub_k3,pub_k4,pub_k2,op_k2,phi_k3}.tsv`;
     - the inventory script and its json/tsv, `setlevel/{universe_common.tsv,
       symbol_qc.json,hgnc_drift_summary.tsv}`, `hgnc_drift_precheck.py`;
     - in `BRCA_TS_ARM_2026-10-08/scripts/`: `cptac2_replicate.py`, `cptac2_residual.py`,
       `extract_phikon.py`, `ts_slide_qc.py`, `brca_tested_genes.py`,
       `lsf_cptac2_ts_extract.sh`, `lsf_cptac2_ts_residual.sh`;
     - `RUNBOOK.sh`, `FREEZE_CHECKLIST.txt` (marked superseded by this file),
       `make_setlevel_sidecar.py`, `test_synthetic_pipeline.py`;
     - the frozen TS list `dryrun/primary_all/{ts_slide_selection.tsv,
       manifest_brca_ts_slides.txt, md5_brca_ts.txt}`;
     - `vial_match.py`, `vial_match.tsv` and its summary;
     - `gdc_inventory_ts.py` and the GDC inventory (`inventory_ts_files.tsv`,
       `inventory_ts_case_summary.tsv`, `inventory_summary.json`,
       `inventory_dx_cases.tsv`);
     - `brca_tested_genes.{tsv,txt}`, hashed as they arrived.
   - The GitHub release excludes `__pycache__/` and the superseded `*_v1.*` files (N8); the
     `_pre-review.md` copy of this file is released.
2. **Server steps 0-2.** Steps 0, 0c and 1 have run and step 2 was started (Status). The
   five deployed scripts the user is running are unchanged and are listed in the sidecar.
   `cptac2_residual.py` and `lsf_cptac2_ts_residual.sh` were revised to v3 after
   deployment and before any use. They are re-deployed at step 5b.
3. **Steps 3-5 (images, headers, shapes).** Step 4b writes the skip list. The step-5
   output is not read before the freeze.
4. **Step 5b, server checks (review M4 and M7).**
   - The sha256 of every deployed file in `$ZW/scripts_ts_arm/` is checked against the
     sidecar: the 7 scripts and the frozen list copy.
   - The DX embedding folder the published job read, and its stride from the tile
     coordinates.
   - The sha256 of `$S/cptac2_residual.py` (expected `1e12432f...`) and
     `$S/extract_phikon.py` (expected `4c616526...`).

   The output is kept in `prerun_checks.txt`. **If anything differs, a dated amendment to
   this file (UTC time, what differed, what changes) is written and released before step
   6.**
5. **Coverage (after freezing, before step 6).** The read-only command of section 2
   (Symbols) is run, locally, and its output is saved as `coverage_check.txt` (step 5c of
   `RUNBOOK.sh`) and brought back with the other files. The rules of section 2 apply.
6. **Steps 6-7.** The DX re-run and the TS and TS-tile arms, on one pinned host.
   - *Failed or interrupted runs (the same wording in the three post hoc
     pre-specifications of file date 2026-10-08).* A run that stops before writing its
     result file is reported with its log and repeated once, unchanged, on the same host
     with a rerun tag; if the pinned host is unavailable, a new smoke/gate is run on another
     host first (recorded) and the run is pinned there; if the rerun also fails, the
     analysis is reported as not completed and no further run is made. A local analysis
     invocation that exits FATAL before writing its result is not a run, is recorded, and
     may be repeated after the input is corrected.
   - *Here:* the result file of a residual job is its summary csv, the last file it writes
     (the per-gene csv and the increment-matrix `.npz` are written before it, so a job
     that stops in between leaves partial outputs and has not completed). The rerun tag is
     the log suffix `_rerun1` (`-o`/`-e` of the `bsub`); partial outputs of the failed
     attempt are renamed with the suffix `_failed1` before the rerun, kept, listed and
     brought back.
   - *Host unavailable.* The DX re-run is the gate for a host. If the pinned host is
     unavailable before the DX re-run has run, a new host is chosen and recorded. If it
     is unavailable after the DX re-run has run, the DX re-run (the new gate) is repeated
     first on another host, recorded, and all three jobs are rerun there; the earlier
     outputs are kept and listed but not analysed. This change of host is a dated amendment
     to this file.
   - At most one repeat of any job and at most one change of host are made. A DX-gate
     outcome of fail is not a failed run.
7. **Brought back:**
   - the three increment matrices, the per-gene tables and summaries;
   - the published `results/brca_protein_results.csv` and `results/summary_brca.csv`, and
     `dx_arm_inputs_outputs.sha256` (all three are used only by the gate);
   - the pooling manifests, the slide-QC tables, the skip list, `prerun_checks.txt`,
     `coverage_check.txt` and the job logs (including any `_rerun1` logs and `_failed1`
     partial outputs).

   Each file's sha256 is recorded on arrival.
8. **Local read-out.** `local_brca_setlevel.py` is run once, with
   `--dx-published-sha256 dx_arm_inputs_outputs.sha256`. It refuses to run without the
   verified sidecar, and it refuses to overwrite a result. An invocation that exits FATAL
   before writing `setlevel_result.*` is not a run, is recorded, and may be repeated after
   the input is corrected; nothing in this file or in the frozen script changes.

## 5. Notes

- Out of scope and unchanged: `main.tex` describes the sparse DX re-extraction as
  "every-other-tile", but stride 512 keeps about one tile in four (the existing F7 to-do,
  N10). The KIRC within-patient DX-versus-TS comparison of the C3 plan is not done.
