# Results record: C1 covariate re-reads (post hoc), Stage L1 and Stage L2 (stages run 2026-11-22 UTC)

Record required by section 9 ("Results record") of `review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md`
(tag `c1-covariates-posthoc-prespec-2026-11-21`, commit fb15ed4). The Stage H record (commit bf15abb) was published
before either local stage ran: annotated tag `c1-covariates-posthoc-stageh-2026-11-22` created 04:40:56Z, pushed
04:41:00Z, GitHub Release published 04:41:18Z.

**Go-ahead.** The corresponding author's written go-ahead for releasing this record (the single word "push", given in
reply to a summary that said the record would hold the outputs of both stages, their logs and their sha256, with the
matrices withheld) was given in the drafting session at **2026-11-23T02:38:48Z** (the message's timestamp in the session
transcript). This text was written and then checked by three independent read-only reviewers after that go-ahead; it was
changed after their review only in the wording of its disclosures (the sha256 of the released `RESULTS_RECORD.md` is in
the commit message). Nothing in the manuscript has been edited; the sentences of section 7 are filled in after this record.

## What was run

Each local stage was run **once**, at its first attempt, with the frozen script (`c1_covariates_local.py`, sha256
`1ea40fa50759366d09f61360d6d2f80dce9703d6d95f6a8ce43ff0f0ae61e779`) and the Stage H manifest sha256 `ca39784874776c7fe0e8faf3c315b1ff8fcb642a9ea06992b9f53813db8c3e13` on the command line. No closed part
(`CLOSED_PARTS.txt` has no entry), pinned hosts UCEC s001 and LUAD s002, sidecar sha256 `e5391d15d70781b80d0e9aa1257cfa04a9bb725500c2924ae60e474a2ea9cf07`.
No attempt preceded either run (`attempts_before_this` is empty in both result files); neither run raised an exception;
no script, definition, set, arm, threshold, reading or sentence was changed after freezing. Both result files carry the
same script, pre-specification, sidecar, manifest and closed-parts hashes.

| Stage | Computation | Wall time | Finished (UTC, `finished_utc`) | Cells |
|---|---|---|---|---|
| L1 | every discovery analysis and set of section 3; every cell on the published matrices | 126.71 min | 2026-11-22T06:48:37Z | 68 |
| L2 | every cell on the Stage H matrices, the change D, the capacity randomisation, the diagnostic; classes, types, case | 1.65 min | 2026-11-22T06:52:20Z | 270 |

The console logs print the local clock of the machine, which is UTC+8 (its 14:48:37 is the 06:48:37Z of
`result_L1.json`). The first console line of Stage L1 (pre-flight) is 12:41:43 local, 04:41:43Z; its computation began
at 12:41:54 local, 04:41:54Z (the first line of `L1/L1.log`), both after the Stage H Release at 04:41:18Z. Local
environment as recorded in the result files: Python 3.13.13, numpy 2.3.5, scipy
1.17.1, scikit-learn 1.9.0, one thread per library, 24 worker processes. The cluster
jobs of Stage H used numpy 1.26.4 (Stage H record).

Pre-flight checks passed before each stage folder was created. Both stages verified the sidecar (100 files; 33 required
files listed) and the Stage H manifest (310 files); Stage L2 in addition verified the Stage L1 outputs against
`result_L1.json` (3 files) and, per cohort, 124 Stage H files by manifest hash, format and one case list. In Stage L1,
in each of the ten discovery analyses (five per cohort) the batched kernel and the reference kernel gave increments
with worst difference 0.00e+00 (`L1/L1.log`).

## Order of reading

Stage L2 was launched 35 seconds after Stage L1 finished (L1 finished at 14:48:37 local; the first line of the Stage L2
console is 14:49:12; the Stage L2 computation ran from 14:50:41 to 14:52:20). **Before Stage L2 had finished, the
assistant that ran the stages read the Stage L1 outputs** in the drafting session: the stage log `L1/L1.log`, the first
part of `L1/result_L1.json` (the summaries of the sets), and all 68 cells of `L1/cells_L1.csv`, which include the six
forward cells (T_C3, T_R20 and T_OP in each cohort) that section 6.2 reads. The corresponding author had seen no reading
before Stage L2 had finished. No cell on a Stage H matrix (primary, capacity, reverse or matched) and no class, type or
case existed before Stage L2 finished. Section 9 does not restrict reading the Stage L1 outputs after Stage L1; Stage L2
verified them against the hashes of `result_L1.json` before it started. The readings below are those the script wrote
(`L2/classes.json`); the manuscript was not edited between the stages or since.

## Readings as written by the frozen script (`L2/classes.json`)

Printed to three decimals as in section 7 (p = 1/1,001 as 0.001); classes use the unrounded values in `classes.json`.
All cells are one-sided and uncorrected, with the unrestricted null. Margin = AUC minus null mean. C3 class H: p < 0.05;
N: p >= 0.05. R20 class H: p_R < 0.05; A: p_R >= 0.05 and p_abs < 0.05; U: otherwise. m_R is the margin of the real
components over the first 25 permutation columns (section 5) and can differ in the third decimal from the margin of the
cell, which uses all 1,000.

| Cohort | Baseline of the confirmatory test | AUC | Null mean | Null SD | p | Margin | Reading |
|---|---|---|---|---|---|---|---|
| uterus (UCEC) | published (no covariates) | 0.597 | 0.498 | 0.036 | 0.004 | 0.099 | |
| uterus (UCEC) | A_C3: three composition scores | 0.551 | 0.515 | 0.035 | 0.161 | 0.036 | C3 class N |
| uterus (UCEC) | A_R20: 20 expression components | 0.546 | 0.508 | 0.041 | 0.199 | 0.038 | R20 class U; m_R 0.039; 46 of 119 placebo margins at or below m_R; p_abs 0.392 |
| lung (LUAD) | published (no covariates) | 0.644 | 0.524 | 0.019 | 0.001 | 0.120 | |
| lung (LUAD) | A_C3: three composition scores | 0.623 | 0.556 | 0.016 | 0.001 | 0.067 | C3 class H |
| lung (LUAD) | A_R20: 20 expression components | 0.598 | 0.564 | 0.017 | 0.029 | 0.033 | R20 class H; m_R 0.029; 0 of 119 placebo margins at or below m_R; p_abs 0.008 |

Set sizes. The published uncorrected sets S_M have 2,566 (uterus) and 2,310 (lung) genes among the
10,512 and 10,724 fit discovery genes. Of these, 2,536 and 2,287 are among the
10,146 and 10,406 genes usable in the confirmatory matrices, which is the n1 of every S_M cell in
`cells_L1.csv` and `cells_L2.csv`. Every fixed-size set (T_M, T_G, T_C3, T_R20, T_OP, T_M|C3, T_M|R20) has the size of
its cohort's published set in discovery; its size in the confirmatory universe is its n1.

Types (section 6.1):

- **uterus**: C3 class N, R20 class U, type **C**; flagged: yes.
- **lung**: C3 class H, R20 class H, type **none**; flagged: no.

**Case 1** (section 6.4): flagged cohorts: uterus (UCEC).

Reverse cells (the covariates in place of morphology; section 6.2):

| Cohort | Covariates | p | Margin | Reading |
|---|---|---|---|---|
| uterus | the three scores alone | 0.016 | 0.102 | outranked |
| uterus | the 20 components alone | 0.321 | 0.024 | not shown |
| lung | the three scores alone | 0.001 | 0.148 | outranked |
| lung | the 20 components alone | 0.003 | 0.054 | outranked |

Forward cells (sets chosen without histology, read on the published matrices): 4 of 6 outranked background.
Matched cells (sets chosen by histology after the covariates entered the discovery baseline): 1 of 4 did.
The cells and their p-values are in `L2/classes.json` and `L2/cells_L2.csv`; the sets and their cut points are in
`L1/result_L1.json` and `L1/cells_L1.csv`.

A class is read as section 6.3 fixes it; this record states no further interpretation. The limits of section 8 apply.

## Released here, with sha256

`RESULTS_SHA256.txt` lists the files below (`sha256sum -c RESULTS_SHA256.txt` from this folder). The result files
`L1/result_L1.json` and `L2/result_L2.json` record the sha256 of the other outputs of their stage; the stage logs
`L1/L1.log` and `L2/L2.log` are hashed here, as section 9 requires. The two console logs are the output of the
command lines of both stages as they ran (CRLF line endings kept, as written); every other file has LF endings. The
bytes are exact on GitHub (raw files) and in a clone with `core.autocrlf=false`; a clone with `core.autocrlf=true`
rewrites line endings and the check then fails.

| File | Bytes | sha256 |
|---|---|---|
| `L1/L1.log` | 4,822 | `eb510b5456fa9bb028bb278f044dacac988ced28e7247eb30f9636ce9161df2d` |
| `L1/result_L1.json` | 18,038 | `39be9aaa4053e5db0f5e606c810c78ebe547baed9311db8661d35322e9210f64` |
| `L1/cells_L1.csv` | 17,337 | `6d5512860da369f387b83ffdede1d023df39c60dc6cfb981acf837a57954fa08` |
| `L2/L2.log` | 939 | `cb9bc6c90899bf35e074d0d4a2172d8beb045fcca254b80d2fc8a33f755fc3fd` |
| `L2/result_L2.json` | 40,825 | `5aef7b1d9bdd53bb4db96ae8f2ebf743a814240b8473cf01a80e4b3f00d5e6b4` |
| `L2/cells_L2.csv` | 89,155 | `e26bbe610aa44c08b64a39ab1f28d00efadcf1e7f8e151e6c1c2222a9aaa5aca` |
| `L2/classes.json` | 38,711 | `ea6068039edf4a6759b488b0388d7850959688643c784e0d30cbf0d7c6e63149` |
| `console/stage_L1_console.log` | 7,771 | `5a7563c589d00613486c9dc277724a3151cbf16314b1372e3c3b25a8a5232dff` |
| `console/stage_L2_console.log` | 1,668 | `c5254a5e8cc2308018cd8dffa7a6eb12e294f90a62c7ea4e9b7a2861004320c5` |

`L1/result_L1.json` names, as written by the frozen script, the last-ranked gene (`last_gene`) of each of the twelve
re-ranked sets, and `L2/classes.json` lists the marker genes of the composition scores, which are those of the
pre-specification; no other gene-level content is released here.

## Two Stage L1 outputs not released here (a departure from section 9)

Section 9 (Results record) says that the outputs of both local stages are committed. Two outputs of Stage L1 are not:
`sets_ucec.csv` and `sets_luad.csv` (one row per fit discovery gene, 10,512 and 10,724 rows, 50 columns). This
departs from that sentence; section 9 (Release) names the published per-gene result tables, not Stage L1 outputs. Reason:
besides the columns Stage L1 wrote (the set flags S_M, T_M, T_M_rank and T_G and, for each of C3, R20, OP, M|C3 and
M|R20, the increment, baseline R2, exceedance count, p, FDR, set flags and rank), each file carries four columns copied
from the published per-gene discovery tables (`pub_r2_rna`, `pub_incremental_r2`, `pub_pval`, `pub_fdr`), which Release
withholds until acceptance, and the file cannot be split without changing its hash. Their sha256 are registered in
`outputs` of `L1/result_L1.json` and listed in `WITHHELD_SHA256.txt` (check them with `sha256sum -c WITHHELD_SHA256.txt`
in the folder that holds the two files, `server_export/results/c1_covariates/L1/` of the author's working copy); they
will be deposited with the other derived data at acceptance. Until then the size of every set, its overlap with S_M,
its floor count and k* are in `L1/result_L1.json` and every Stage L1 cell is in `L1/cells_L1.csv`; the genes of the sets
are not public.

| File | Bytes | sha256 |
|---|---|---|
| `sets_ucec.csv` | 5,896,441 | `73e4f6e92575611298bf1ce7c1cab010902eea507a6d5de16b537258e280971b` |
| `sets_luad.csv` | 6,021,639 | `90dcf0a4341e3a54c16b1d59ae5de963fecbe69fb93d632fb17562e65c1b161b` |

Also still withheld under section 9 (Release): the 248 Stage H files (246 matrices and 2 covariates records; fixed by
the manifest, as in the Stage H record), the two discovery score tables, the inputs files and the confirmatory matrices
(fixed by the sidecar), and the published per-gene result tables including the two operator-axis tables.

## Next

1. A dated entry of `C1_LUAD_ADDENDA.md` quotes the Abstract variant that section 7.4 assigns to the case and type
   above; it is released (tag and GitHub Release) before the manuscript is edited.
2. The sentences of section 7 are then filled in as written, by transcription from `classes.json` and the two cell
   tables.
