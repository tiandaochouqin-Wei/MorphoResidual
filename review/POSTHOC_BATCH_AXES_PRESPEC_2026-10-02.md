# Post hoc batch-axis analyses: pre-specification (2026-10-02)

**Status.** Written 2026-10-02 after the registered C1-LUAD test was unblinded and after
the scoping in `review/BATCH_AXES_SCOPING_2026-10-02/`, before any job below is
submitted. The corresponding author authorised these analyses in writing on 2026-10-02
(chat instruction to proceed with items 1-4, item 4 being this specification with its
sub-items a-c). Frozen by the sha256 of this file and the inputs and scripts in section F,
released publicly on GitHub before the first job runs. Answers the TMT-plex and
embedding-medium items of the NC ceiling review (`review/NC_CEILING_REVIEW_2026-09-27.md`)
and two findings of the scoping.

## 0. Scope and status of every analysis below

- **Post hoc and reported whatever the result.** Nothing below was pre-registered. No axis,
  arm, threshold or reading is added, dropped or changed after any of these jobs reports.
- **No effect on the confirmatory readings.** The registered C1-LUAD reading and the
  published C1-UCEC reading stand as reported. Analysis L re-reads C1-LUAD under a
  different null after its result is known: its observed statistic is known, its null is
  not, and it is labelled post hoc wherever it is reported.
- **No effect on the published discovery counts.** Table 1's Sig., Batch-corr. and
  Strictest columns stay as published; every new count is reported beside them.
- **Correction already made.** The signed C1-LUAD rule's section 6 statement that PDC gives
  no case-to-plex membership is wrong (C1-LUAD addenda Entry 7). The per-channel fields
  `tmt_126` to `tmt_131c` of PDC `studyExperimentalDesign` give it for all five discovery
  studies and PDC000489; a second PDC route agrees on all 1,343 channels.

## P. TMT plex as a batch axis (discovery, all five cohorts)

**Inputs.** `pinned/case_to_plex_<cohort>.tsv` (columns 1-2 of the scoping
`plex_map_<cohort>.tsv`), read by the patched `batch_leak_check.case_batch_labels`, which
aborts if an analysed case is unlabelled. UCEC C3N-01825 (aliquots in plexes 10 and 16) is
assigned to plex 10, the plex of its non-Disqualified aliquot.

**Computation.** `residual_analysis_sitepack.py <cohort>` unchanged, through
`lsf_sitepack.sh <cohort> plex`, with `MORPHO_N_PERM=1000` and every other setting as for the
operator axis, except **`MORPHO_MAX_LEVELS=30`**, so every plex with at least three analysed
cases gets a dummy (the default 12 would leave 43/51/22/58 cases in CCRCC/LUAD/UCEC/PDAC
uncorrected). Also `batch_leak_check.py <cohort> plex` for the eta-squared diagnostic.

**Reported per cohort.** Genes tested; Plex-corr. count (`fdr_batchresid < 0.05` and
`incremental_r2_batchresid > 0`, defined as Table 1's Batch-corr.); the Strictest-plex count
(`fdr_both`); PC variance retained; the batch-only control against its noise-matched
partner; overlap (count, Jaccard) with Sig. and with the operator Batch-corr. set;
*retained* = median residualised increment / median plain increment over Sig.; top-15
family enrichment as in Supplementary Table `tab:c2`.

**Readings.** P = Plex-corr. count; O = operator Batch-corr. count (815/726/201/680/682).
1. P >= 50 and P >= 0.5 O: the residual survives protein batch about as well as acquisition
   batch.
2. P >= 50 and P < 0.5 O: plex-sensitive. If *retained* >= 0.5 the drop is read as mostly
   power lost to strata of 3-8 cases; if < 0.5, plex-aligned morphology carried a
   substantial part of the estimate.
3. P < 50: the cohort is not positive under plex control, and the Introduction's
   "processing plex" sentence names it.
Plex is not aliased with operator (scoping, p = 0.12-0.64); in CCRCC it tracks accrual stream
(p = 0.002) and scan quarter (p = 0.03). A within-plex null credits all between-plex
association to chance, so no outcome separates protein batch from accrual case-mix.

**Cannot establish:** absence of within-plex confounding; channel-position, ratio
compression or reference-pool effects; a plex-corrected effect size.

## L. C1-LUAD confirmatory test under a within-plex null (post hoc, sub-item a)

**Script.** `c1_luad_plex_reread.py` via `lsf_c1_luad_plex_reread.sh`, importing the frozen
`c1_run_test_luad.py` machinery. Input `pinned/case_to_plex_luad_c1.tsv` (113/113 frozen cases;
27 plexes of 4-5 cases).

**Gate.** It first recomputes step 1 with the frozen code path and must reproduce
n = 112, K = 22, 10,406 usable genes (2,287 selected), AUC 0.6439668501223573 (to 1e-9) and
p1 = 1/1,001, or stop.

**Computation.** B = 1,000 within-plex permutations from
`RandomState(crc32(b"c1_luad|perm_plex"))`, strata in sorted label order, singletons and
unlabelled cases fixed. L1: the 2,310-set AUC and rho under that null. L2: the 726 set on
components residualised on plex dummies (every plex with >= 3 cases, one dropped for full
rank), same draws. L3: L1 without C3N-02230 and C3N-02240 (one aliquot measured in two
plexes).

**Readings.** L1 p < 0.05: "the set-level ordering also exceeds a within-plex null (post
hoc)". L1 p >= 0.05: "not shown under a within-plex null (post hoc); the registered reading
is unchanged", with the note that 27 strata of 4-5 cases make this null restrictive and
credit all between-plex association to chance. L2 and L3 are descriptive. None of them
enters the Abstract's registered clauses; a one-sentence report goes in Results.

## M. Embedding medium (OCT-frozen versus FFPE)

**M0. Labels.** `pinned/slide_embedding_medium.tsv`, built from the public PathDB file
`cptac_metadata_07-09-2024.csv` (downloaded 2026-10-02T12:06:15Z, sha256
`6f999cec4381ea34f245c2faabdbc828837b0125f6a9950bb5d42ec4febe6da0`), joined by exact
`slide_submitter_id` to each cohort's analysed Primary Tumor slides (the slides with an
embedding, as the batch suite defines them); PathDB `Specimen_Type` is never used. A
patient is FFPE, OCT or mixed. Reported per cohort: slide and patient counts per level,
join-failure rate, Cramer's V of medium against operator with a 500-shuffle null. If more
than 5% of a cohort's analysed slides fail to join, that cohort's medium analyses are
replaced by a disclosure.

**M1. Medium as a batch axis.** Run `lsf_sitepack.sh <cohort> medium` (default level cap) in
every cohort whose patient labels have >= 2 levels with >= 3 patients each and a largest
level <= 90%. Reported as in P. Readings: as in P with O the operator count, except in
UCEC, where operator is nested in OCT status at 0.99 and is the finer correction: there
M >= O reads "medium adds nothing beyond operator", M < O is reported without
interpretation, and M < 50 reads "UCEC is not positive under medium control".

**M2. FFPE-only re-run.** In every cohort with k >= 3 OCT-only patients, the discovery
estimator (`residual_analysis.py`, unchanged; B = 1,000) is re-run on arms that differ only in
the slide map or crosswalk, built by `posthoc_build_arms.py` and run by
`lsf_posthoc_residual.sh`:
- **R**: current inputs, re-run now (the control every arm is compared with);
- **F**: OCT slides excluded (mixed patients keep their FFPE slides; OCT-only patients drop);
- **D01-D19**: all slides of k random analysed patients excluded, drawn with
  `RandomState(1000+s).choice(sorted IDs, k, replace=False)`, s = 0..18.
Reported: counts for every arm, overlap of F with R, *retained* on R's Sig. set, and F's rank
among the D arms. Readings: F >= 50 and F >= median(D): OCT patients carry no more signal
than an equal-sized random set. F below all 19 D counts (one-sided p <= 0.05): OCT patients
carry disproportionate signal; in UCEC this is written as dependence on the
OCT/Ukraine/operator stratum (Supplementary Note `snote:ucec`), never on medium alone.
F < 50: the cohort is not positive on FFPE slides alone, stated in Results and Limitations.
Anything else: indeterminate (D mimics the patients lost, not the slides mixed patients
lose).

**M3. Cohorts that do not qualify** under M0-M2 get a disclosure in place of analysis,
stating which criterion failed (for example, GBM if every analysed slide is FFPE).

**Cannot establish:** whether any effect is freezing artefact rather than the accrual
factors nested with medium.

## Q. PDC-Disqualified aliquots (sub-item c)

**Finding.** PDC flags as Disqualified the analysed tumour aliquot of UCEC C3L-00157,
C3L-00356, C3L-00938, C3L-01247 and C3L-01253, a second aliquot of UCEC C3N-01825, and the
analysed aliquot of LUAD C3N-00545 (confirmed or corrected per cohort by
`build_pdc_aliquot_status.py`, which also checks CCRCC, GBM and PDAC).

**Computation.** In every cohort with at least one analysed Disqualified aliquot: arm **Q**
(the crosswalk without Disqualified aliquots) against arm **R**, through
`lsf_posthoc_residual.sh`.

**Reported.** Patients dropped, Sig. counts in R and Q, overlap (count, Jaccard), *retained*.
**Readings.** Q >= 50 and Jaccard(Q, R) >= the median Jaccard between R and the D arms of M2
in the same cohort (or, where M2 does not run, >= 0.8): "robust to excluding Disqualified
aliquots". Otherwise the drop is reported and the affected cohort's per-gene counts are
flagged in Results. The published counts stay as published either way; the Methods state
which aliquots are Disqualified.

## C. Corrections made whatever the results (done 2026-10-02)

1. main.tex Results and Supplementary Note: the false "no patient-to-plex assignment"
   sentences are corrected (C1-LUAD addenda Entry 7).
2. main.tex Results and Methods: the batch suite residualises the morphology components
   only, not protein; the text said both (code: `residual_analysis_sitepack.py`, PCs only).
3. main.tex Methods: CCRCC is described as the only cohort for which a plex label "had been
   assembled when this analysis was run", not the only one for which one exists.

## R. Where the results go

Plex: a Plex-corr. column in Table 1 (caption note on the level cap), one Results paragraph,
`sec:batchmeth` (source: PDC `tmt_*` fields, queried 2026-10-02; the cap), `sec:confmeth`
(replace "not run in the other four cohorts"), a Supplementary table with the P outputs and
the earlier CCRCC plex+purity run (1,340 of 9,635, labelled as run before this
specification). L: one Results sentence after the registered LUAD reading, and the
Supplementary Note `snote:c1luad`. M and Q: Supplementary Note `snote:ucec`, Supplementary
Methods `sm:census`, a Supplementary table, and one sentence each in `sec:batchmeth` and
`sec:c1res`.

## F. Freeze record

Every script and input below is fixed as of this file's release, before any job in
sections P, L, M or Q is submitted. This file's own sha256 is in its sidecar
`POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md.sha256`. Any later change to any of them is a
dated amendment released separately, naming what changed and why; a run made with a
different hash is reported as such.

None of these scripts is hashed in the signed C1-LUAD rule's section 7.2, and none of the
rule's section 7.2 scripts is modified: `c1_luad_plex_reread.py` imports the frozen
`c1_run_test_luad.py` rather than reimplementing it, and `batch_leak_check.py` (not a
section 7.2 file) keeps its scanner, operator, quarter, mpp and cohort_stream axes
unchanged. All 44 section 7.2 files were re-verified byte-identical at this freeze.

**Scripts (post hoc; none is hashed in the signed rule's section 7.2)**

| File | sha256 |
|---|---|
| `server_export/scripts/batch_leak_check.py` | `ac492a795a9997286eb23191d416cd47fc10e4955bc25cf25d9f9934213edc6b` |
| `server_export/scripts/build_pinned_plex_maps.py` | `885343a9ff8d520d8b852aff8c68f7d9cab3a9b671cf8a63d6fbbd9e7a8fa43c` |
| `server_export/scripts/build_pinned_slide_medium.py` | `32273d4f33ac1adf7c0337b216c48ed1083a23130e7980b76157bd2be17b2c57` |
| `server_export/scripts/build_pdc_aliquot_status.py` | `fb612b84ac25e73f6d7ea8a3a661cc64ec77974dfcf3b471f5dd71c78ec4cd3f` |
| `server_export/scripts/c1_luad_plex_reread.py` | `45a9317df29112f5efe7eb43f3300dcbde386800e71b0c3b43f27c163f707d2b` |
| `server_export/scripts/lsf_c1_luad_plex_reread.sh` | `aaae7cd54de6725e7bb638426467a27aab8724cf94b5f65069b95da60684cd66` |
| `server_export/scripts/posthoc_build_arms.py` | `b9b464f4c7fd70dd2fa273c8dbbe5cd305d3712c078072cc4a291a4aa99a01cc` |
| `server_export/scripts/posthoc_summarise_arms.py` | `335adf068364dce9b7544cbec53997a71e41e8443e1b6bacd32c027717a0eab5` |
| `server_export/scripts/lsf_posthoc_residual.sh` | `d46fac9e0d78961a30a21058ee47135c601c5db2c6049e7ae5253b6a8dc05d97` |

**Pinned inputs the scripts assert**

| File | sha256 |
|---|---|
| `server_export/scripts/pinned/case_to_plex_ccrcc.tsv` | `bbc400d33a5fea38380bc68226b2c21624ca7f31fbec07566264c2da2df401a6` |
| `server_export/scripts/pinned/case_to_plex_luad.tsv` | `014868ad61a5cfbe1cc112d603e5d2d669ddd82f03442e3f998cbd97170b616e` |
| `server_export/scripts/pinned/case_to_plex_ucec.tsv` | `e3f383d85f11a6e35a63676ffbdebbd7e7bbad71bed38a199a594772ae141f03` |
| `server_export/scripts/pinned/case_to_plex_gbm.tsv` | `e0aa9e78c7bf1edf8cc67fd7995037fc09327822e3163f882751beaf55589b3a` |
| `server_export/scripts/pinned/case_to_plex_pdac.tsv` | `c612c808a7a55134f8e393613c55eb56d9443a2c6863abb6b5c798a871710099` |
| `server_export/scripts/pinned/case_to_plex_luad_c1.tsv` | `491f7355b9ad5edfbf12d89cf7305713c38b54d6a8e91daa13970ed4c8829258` |
| `server_export/scripts/pinned/slide_embedding_medium.tsv` | `36783b539e14f3bad73d59d9769f7dd2d01ea21d08f0cd4d79dbfde48d098ef1` |
| `server_export/scripts/pinned/pdc_aliquot_status.tsv` | `61ff0547c9209b9ec6f28dfa8a26cf15ab0a99a0a9beacd6be194fdd2a382a02` |

**Public source of the medium labels**

| File | sha256 |
|---|---|
| `review/BATCH_AXES_SCOPING_2026-10-02/pathdb/cptac_metadata_07-09-2024.csv` | `6f999cec4381ea34f245c2faabdbc828837b0125f6a9950bb5d42ec4febe6da0` |

**Record of the pre-release code review.** The nine scripts were audited against this
file before the freeze, and four changes were made as a result, all before any job ran:
the plex axis's abort now applies to the analysed population the pinned maps were built
from and reports, rather than aborts, when a caller passes the wider WSI-embedding index;
the C1-LUAD plex re-read writes to `luad_c1/results/posthoc/` instead of the registered
test's own results directory; the medium and Disqualified-aliquot arms are now built where
each cohort's data meets the M2 and Q conditions, instead of for a hard-coded pair of
cohorts, with the M3 disclosure emitted where a condition fails; and the M2 arms are
gated on the specified k >= 3 and on M0's 95% join rate, both of which the code had not
enforced.
