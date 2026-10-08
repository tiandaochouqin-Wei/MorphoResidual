# Post hoc C1-UCEC confirmatory test under a within-plex null: pre-specification

**Dates and times.** The date in this file's name (file date 2026-10-08, local time UTC+8)
is the local date of the drafting session. Every time and date stated inside this file is UTC and carries
its UTC date. The text was drafted on 2026-10-07 between about 21:00Z and 23:59Z. It was
revised on 2026-10-08 from about 00:20Z after a cross-check of the three post hoc
pre-specifications of this file date, and the blind and smoke record (section 7) was
completed at about 00:40Z, before the sidecar was written.

**Status.** Drafted after the published C1-UCEC reading (2026-09-16), the D4
within-operator re-read (2026-10-01) and the alternative-selection-set re-read
(`POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md`; results commit `ab96648`,
2026-10-07T20:55Z), and before any within-plex statistic for C1-UCEC was computed by
anyone. The corresponding author asked for this analysis on 2026-10-08 local time
(2026-10-07 UTC, before the drafting began). It is the UCEC counterpart of section L of
`POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md` (C1-LUAD under a within-plex null) and follows it
wherever the two cohorts allow; where it departs, the departure is named. The PDC queries
behind the labels ran at 2026-10-07T21:02:11Z-21:02:29Z, the plex map was built at
21:03:27Z, and the script and wrapper were finalised at 23:09Z (sha256 `339c2fd6...` and
`ce1bb071...`; full values in the sidecar). An adversarial review of the first draft of
this file (2026-10-07, before freezing) is reflected throughout; it did not lead to any
change in the script or the wrapper after 23:09Z. A cross-check of this file against the
two other post hoc pre-specifications of this file date (the BRCA frozen-section arm and
the measurement-error checks) followed on 2026-10-08 from about 00:20Z; it unified the
Abstract policy, the date convention, the rounding rule and the rule for failed runs across
the three files, and changed no script or wrapper.

Frozen by the sha256 of every file listed in
`review/C1_UCEC_PLEX_2026-10-08/SIDECAR_FILES.txt`, written to the sidecar
`POSTHOC_UCEC_C1_PLEX_PRESPEC_2026-10-08.md.sha256`, and released publicly on GitHub before
the `--perms 1000` run (section 7). The `--perms 5 --blind` and `--perms 5 --tag smoke`
steps may run before the sidecar is written; they count only if the `this script sha256`
line in their LSF .out equals the sidecar value; otherwise they do not count and are
repeated with the final script, before or after freezing.

## 0. Scope and prior exposure

- **Post hoc, reported whatever the result.** Nothing below was pre-registered. After
  freezing, no arm, case, seed, threshold, reading or wording is added, dropped or changed.
- **No effect on the published or registered readings.** The published C1-UCEC reading
  (AUC 0.597; unrestricted null 0.498 +/- 0.036, p = 0.004), the D4 within-operator
  re-read, the registered C1-LUAD reading and the alternative-selection-set re-read stand
  as reported. Nothing here is pooled with any of them, and no favourable result of this
  analysis enters the Abstract; section 4(f) fixes the one qualifier an unfavourable P1
  adds.
- **Prior exposure (disclosed).** The following were known to the session that wrote this
  file when it was written:
  1. The observed C1-UCEC statistics (AUC 0.597, rho 0.294) are published; the within-plex
     null is not known.
  2. D4 within-operator re-read: null 0.506 +/- 0.032, p = 0.004; with the morphology
     components residualised on operator dummies, AUC 0.567, p = 0.006.
  3. The alternative-selection-set re-read, including its two plex-corrected discovery
     sets (S4 Plex-corr., S5 Strictest-plex), both of which held under the unrestricted
     and within-operator nulls (Holm-adjusted). Those nulls were not stratified by plex.
  4. C1-LUAD section L: L1 AUC 0.644 against a within-plex null of 0.546 +/- 0.019,
     p = 0.001 (the null moved from 0.524 unrestricted to 0.546 within plex); L2 (726 set,
     plex-residualised components) 0.651 (0.6515) against 0.519 +/- 0.024, p = 0.001; L3
     also known.
  5. The UCEC discovery Plex-corr. count of Table 1 (2,124 genes, against 2,566 Sig.;
     section P of the 2026-10-02 specification).
  6. For C1-UCEC itself, only label counts have been computed from the plex labels: 16
     plexes, their sizes, the two multi-plex cases and the number of P2 dummies.

  No within-plex statistic and no plex-by-data association for C1-UCEC has been computed
  by anyone.

## 1. Plex labels

`case_to_plex_ucec_c1.tsv` (in `review/C1_UCEC_PLEX_2026-10-08/`, copied to the server as
`pinned/case_to_plex_ucec_c1.tsv`, sha256 `4fe057a9...`) is columns 1-2 of
`plex_map_ucec_c1.tsv`, built by `build_plex_map_ucec_c1.py` from the PDC000439
`studyExperimentalDesign` fields `tmt_126` to `tmt_131c` (16 TMT11 runs, 176 channel rows)
for the 138 C1-UCEC cases (`pinned/c1_case_ids_ucec.txt`, sha256 `d777fc88...`, the order
of the D4 matrices). The three raw PDC responses, the query log with their sha256, and
every script that produced the map are in the sidecar list. An independent route
(`paginatedCasesSamplesAliquots`, `crossval_route2.py`) agrees on all 176
run-channel-aliquot triples and on all 138 case labels. Result: 16 plexes of 8-9 cases
(six of 8, ten of 9); every case labelled; no singleton; no case held fixed.

**The two multi-plex cases.** C3L-05571 (CPT0275270003 in plex 04, PDC-Disqualified;
CPT0275270003.1 in plex 16) and C3N-02978 (CPT0244470003 in plex 12, PDC-Disqualified;
CPT0244470003.1 in plex 16) are assigned following section P's treatment of UCEC discovery
case C3N-01825 (plex of the non-Disqualified aliquot): 16; the builder used for C1-LUAD
would give 04/12 (`differs_from_original_rule = True` in `plex_map_ucec_c1.tsv`); chosen
from metadata before any statistic; P3 bounds this choice. (C1-LUAD's two multi-plex cases,
C3N-02230 and C3N-02240, were one aliquot measured in two plexes with no Disqualified
aliquot, and kept the original "majority of channels, ties to the smallest plex" rule; they
are not a precedent for this assignment.) The second route's assignment of these two cases
also used the aliquot status from `biospecimenPerStudy`, so for them the two routes are not
fully independent. The published protein matrix averaged every non-Normal aliquot of a
case (`c1_pull_omics.py`, listed in the sidecar; `groupby(level=0).mean()`), and both aliquots of each case are in
its pull rule, so these two protein values may combine two plexes; P3 handles this.

## 2. Gate

`c1_ucec_plex_reread.py`, via `lsf_c1_ucec_plex_reread.sh`, first checks the sha256 of the
plex map, the case list, the discovery results, the published result file and the three
frozen modules it imports (`c1_run_test.py`, `split_replication_power.py`,
`residual_analysis.py`); any mismatch stops the job before data are read. It then
recomputes the published test with the frozen code path, seed `crc32(b"c1_ucec|perm")` and
B = 1,000, and must reproduce n = 138, K = 28, 10,153 genes tested, 10,146 usable (2,536
selected, 7,610 background), AUC 0.5974425564900907 and rho 0.2935675639399924 (to 1e-9),
the published null mean, SD and p_AUC = 4/1,001 and p_rho, and the hashes recorded in the
sha256-pinned published result file, or stop. If the server copy of a frozen module
differs from the hash in the script, the job stops; the frozen code is not edited to pass
the gate.

## 3. Computation

B = 1,000 within-plex permutations from `RandomState(crc32(b"c1_ucec|perm_plex"))`, strata
in sorted label order, singletons and unlabelled cases held fixed (there are none).

The within-plex permutation, the plex-dummy design and the driver are adapted from
`c1_luad_plex_reread.py` and `c1_ucec_stratified_reread.py`; both are in the sidecar.

- **P1 (the reading).** The published 2,566-gene set (2,536 usable) on the same
  components: AUC and rho under the within-plex null. The observed AUC and rho equal the
  published ones; only the null changes.
- **P2 (descriptive).** The published set on components residualised on plex dummies
  (every plex with >= 3 cases; 15 dummies, plex 07 dropped for full rank; which plex is
  dropped does not change the residuals), same draws as P1. Its observed AUC differs from
  P1's. This departs from C1-LUAD L2, which used the registered 726-gene set; C1-UCEC has no
  registered batch-corrected set, so the published set is used, as in the D4 descriptive
  arm.
- **P3 (sensitivity).** P1 without C3L-05571 and C3N-02978 (n = 136, K = round(27.2) = 27,
  components refitted; plex 16 keeps 7 cases), from a fresh
  `RandomState(crc32(b"c1_ucec|perm_plex"))`.

All p values are one-sided, p = (1 + #{null >= observed})/(1 + B), floor 1/1,001.

## 4. Readings and wording, fixed now

**(a) What is read.** Only P1's p_AUC is read. rho and p_rho are reported for P1, P2 and P3
but are not read. P2 and P3 are not read except as stated in (c).

**(b) Two tiers, as in section L.** p_AUC < 0.05: "exceeds its null"; p_AUC >= 0.05: "not
shown to exceed its null". A value with 0.05 <= p_AUC < 0.15 takes the p >= 0.05 sentence
and is not called suggestive. This departs from the three tiers of the UCEC rule's section
5 and of the alternative-selection-set re-read (0.05-0.15 "suggestive"); the direct
precedent for a within-plex re-read of a confirmatory test is section L, which has two.

**Results wording** (main.tex, at the end of the alternative-selection-set paragraph,
after "...a note added after the result.", as a new sentence; not after the operator
sentences, whose opening sentence states a provenance that does not apply here). [date] is
the UTC date of the sidecar's `# written` line; [m], [s], [z], [p] are P1's null mean, null
SD, z and p_AUC (z = (AUC - [m])/[s], one-sided). Exactly one of:

- P1 p_AUC < 0.05: "A fourth post hoc analysis, specified on [date] in a dated file released
  publicly before it ran (Supplementary Note~\ref{snote:c1posthoc}), permutes patients only
  within TMT plex (16 plexes of 8--9 patients, none held fixed): the ordering also exceeds
  its null (AUC $0.597$ against [m] $\pm$ [s], $z=$[z], one-sided $p=$[p])."
- P1 p_AUC >= 0.05: "A fourth post hoc analysis, specified on [date] in a dated file
  released publicly before it ran (Supplementary Note~\ref{snote:c1posthoc}), permutes
  patients only within TMT plex (16 plexes of 8--9 patients, none held fixed): the ordering
  was not shown to exceed its null (AUC $0.597$ against [m] $\pm$ [s], $z=$[z], one-sided
  $p=$[p]); the published test against the unrestricted null stands as reported. With 16 strata of 8--9 cases this null credits all
  between-plex association to chance and lowers power, so a non-significant result means
  the data cannot separate the signal from plex structure, not that it is an artefact."

**(c) P1 and P3 disagree.** The P1 sentence stays as above. If P3's p_AUC falls on the
other side of 0.05 from P1's, the text "; without the two patients whose protein values may
combine two plexes, $p=$[p3]" is inserted immediately after "$p=$[p]" inside the
parenthesis of the P1 sentence. In every other case P3 is reported in the supplement only.

**(d) P2** is descriptive in every direction and changes no wording. Because its observed
AUC differs from P1's, it is reported with its margin (observed AUC minus its own null
mean).

**(e) Null shift.** P1's null mean shift and SD ratio against the reproduced unrestricted
null are reported in the supplement with the fixed note: "a shift is not an attribution;
block restriction raises a null mechanically" (as for the operator re-read).

**(f) Abstract and Discussion.** No sentence is added to the Discussion in any case.

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

Here: the qualifier is triggered by P1's p_AUC alone (section 4(a)); P2 and P3 never change
the Abstract. With P1 p_AUC < 0.05 the Abstract is unchanged. The alternative-selection-set
re-read's single abstract exception (its case (C) with n_below = 0) is governed by its own
file and is not changed here.

**(g) Format.** AUC, null mean, null SD, margin and p are given to 3 decimals in the text
and tables (p as 0.004, not 4/1,001), z to 2. Readings use unrounded values; if rounding
would print a value on the other side of a threshold, decimals are added until it does not
(e.g. p = 50/1,001 printed as 0.04995; z = 1.96 not "2.0"). Under this rule 0.050 or 0.0500
is never printed for a p below 0.05, and a p of exactly 0.05 or above takes the p >= 0.05
sentence.

## 5. What is reported for each arm

In a new last subsection of Supplementary Note `snote:c1posthoc` and a three-row
Supplementary Table (P1, P2, P3), beside the reproduced step-1 row:

- **Step 1 (reproduced):** AUC, null mean, SD, p_AUC; rho, p_rho; the gate verdict.
- **P1:** AUC; within-plex null mean and SD; z; p_AUC with the exceedance count; rho, its
  null mean and SD, p_rho; null mean shift and SD ratio against step 1 (with the note in
  4(e)); number of plexes, plex sizes, cases held fixed (0).
- **P2:** AUC; null mean and SD; margin; z; p_AUC; rho, p_rho; number of dummies (15),
  the plex dropped for full rank, PC variance retained.
- **P3:** n (136), K (27); AUC; null mean and SD; margin; z; p_AUC; rho, p_rho; plex sizes
  after removal; whether the usable gene list differs from step 1.

## 6. Manuscript changes that follow the result (applied after it, not now)

1. main.tex Results: the sentence of section 4, at the position stated there (with the
   P3 insertion of 4(c) if triggered).
2. main.tex Methods (the C1-UCEC post hoc paragraph, after "...the specification discloses
   both.", around lines 1571-1600): one sentence describing the fourth analysis: a null
   that permutes patients only within TMT plex (16 plexes of 8--9 patients from PDC's
   per-channel experimental-design fields), with the plex-residualised and two-patient-removed
   descriptive variants, specified on [date] in a dated file released publicly before it
   ran.
3. main.tex `sec:batchmeth`, the *Protein batch* sentence: add the provenance of the
   C1-UCEC labels (PDC000439 `tmt_126` to `tmt_131c`, retrieved 2026-10-07 UTC; a second
   PDC route agrees on all 176 channels and all 138 patients; the two multi-plex patients
   assigned to the plex of their non-Disqualified aliquot).
4. supplement.tex, opening of `snote:c1posthoc` (around lines 531-535): "the third was
   specified later, in a file of its own (last subsection)" becomes "the third and fourth
   were specified later, each in a file of its own (last two subsections)".
5. supplement.tex, alternative-selection-set subsection (around line 654): after the
   sentence "No confirmatory null is stratified by TMT plex or embedding medium, so a set
   that holds says nothing about the confirmatory ranking's robustness to those axes.",
   insert the cross-reference "(The published set alone is re-read under a within-plex
   null in the next subsection.)". That sentence is the wording fixed by the
   alternative-selection-set specification; only the cross-reference is added, the
   sentence itself is not changed.
6. supplement.tex: the new subsection and table of section 5.
7. Date convention of the third analysis (applied before this file, by the orchestrator;
   referenced here, not repeated). Every manuscript sentence says "specified on [UTC date of
   the sidecar's `# written` line] in a dated file released publicly before it ran". The
   alternative-selection-set sidecar was written 2026-10-07T20:33:56Z, so main.tex
   (around line 1592) and supplement.tex (the heading around line 613) now read "7 October 2026"
   (UTC) for the third analysis, the same convention as [date] above. The fourth analysis
   takes the UTC date of its own sidecar's `# written` line, which may differ from 7 October.
8. main.tex Abstract (the uterus clause, lines 54-56 as of 2026-10-08): only if P1 p_AUC
   >= 0.05, the qualifier of section 4(f), exactly as written there. Otherwise unchanged.

## 7. Freeze and run order

1. **Finalise.** This file, the script (`POST_HOC_GOVERNED_BY` = this file's name), the
   wrapper, the plex labels and their provenance are final. The wrapper header and the
   server checklist (`review/SERVER_RUN_2026-10-08.md`, step 4) add
   `pinned/c1_case_ids_ucec.txt` to the upload list if it is not on the cluster.
2. **Blind and smoke** on the server: `--perms 5 --blind` (counts only; computes no
   association statistic; it still runs `prepare()`, which includes the PCA) and
   `--perms 5 --tag smoke`, which computes no within-plex null; smoke
   recomputes the observed statistic and the first 5 draws of the published unrestricted
   null. From each LSF .out, the LSF job ID, the host and the `this script sha256` line are
   recorded in the table below. These are the only cells of this file filled after
   drafting, and they were filled before the sidecar was written. If a smoke run misses the
   gate, the smoke may be resubmitted on another host, as the wrapper header says. Every
   blind and smoke attempt, including one that failed before computing anything, gets its
   own row in the record table (job, host, script sha256, outcome). The analysis host is the
   host of the first smoke attempt that passed the gate with a `this script sha256` equal to
   the sidecar value. The wrapper header's `bsub` examples (without `bash`, the third without
   `-m`) and step 4c of the server checklist are superseded by item 4. Any blind, smoke or
   gate run made after the sidecar is written is recorded in the results record brought
   back under item 5 (with its LSF .out), not in this file, whose sha256 is frozen.
3. **Freeze and release.** The sidecar is written from the list in
   `review/C1_UCEC_PLEX_2026-10-08/SIDECAR_FILES.txt` (paths relative to
   `MorphoResidual_paper/`, `__pycache__/` excluded), with a `# written <UTC time>` line:

   ```
   cd MorphoResidual_paper
   { echo "# sha256 of the frozen files for POSTHOC_UCEC_C1_PLEX_PRESPEC_2026-10-08.md (paths relative to MorphoResidual_paper/)";
     echo "# written $(date -u +%Y-%m-%dT%H:%M:%SZ)";
     xargs sha256sum -b < review/C1_UCEC_PLEX_2026-10-08/SIDECAR_FILES.txt; } \
     > review/POSTHOC_UCEC_C1_PLEX_PRESPEC_2026-10-08.md.sha256
   ```

   The sidecar's script value must equal the `this script sha256` recorded from the counted
   blind and smoke .out files (a pre-final-script row of the table does not count; see the
   note under the table). This file, the listed files and the sidecar are then released on
   GitHub under the rule below.

   **What is released (the same rule in the three post hoc pre-specifications of file date
   2026-10-08).** This file, its sidecar and every file the sidecar lists are released on
   GitHub before the run, except derived data files (`.npz` and `.tgz` files and the per-gene
   result tables named below) and third-party source tables under `figures/figdata/`. These
   follow the repository's data policy: each is fixed by its sha256 in the sidecar, the
   derived data files are deposited with the other derived data at acceptance, and the
   third-party tables are public at their sources. The release commit message lists every
   withheld path.

   Here, withheld: `server_export/scripts/pinned/c1_ucec_d4_stratified.npz` and
   `server_export/pinned/ucec/residual_results_tumoronly.csv` (the published per-gene
   discovery results). No third-party table is listed.
4. **The analysis.** One `--perms 1000` run, pinned to the analysis host of item 2:
   `bsub -q smp -n 8 -R "span[hosts=1]" -m <analysis host> -o c1_ucec_plex_reread.out
   bash lsf_c1_ucec_plex_reread.sh --perms 1000` (C1-LUAD precedent: `-m s002`; `bash` is
   required, see the record table). Only one run that writes a result is reported.

   *Failed or interrupted runs (the same wording in the three post hoc pre-specifications
   of file date 2026-10-08).* A run that stops before writing its result file is reported
   with its log and repeated once, unchanged, on the same host with a rerun tag; if the
   pinned host is unavailable, a new smoke/gate is run on another host first (recorded) and
   the run is pinned there; if the rerun also fails, the analysis is reported as not
   completed and no further run is made. A local analysis invocation that exits FATAL
   before writing its result is not a run, is recorded, and may be repeated after the input
   is corrected.

   Here: the result file is `c1_ucec_plex_reread_result.json`; a run stops before writing
   it at the gate (for example through a LAPACK difference on the host) or for any other
   reason (a crash, the walltime, a node failure). The rerun tag is `--tag rerun1`. The
   smoke/gate on another host is `--perms 5 --tag smoke` there, which must print `SMOKE
   VALIDATED` with the sidecar's script sha256. "Recorded" means recorded in the results
   record of item 5, not in this file. This analysis has no local analysis invocation. No
   second run that passes the gate is ever made.
5. **Brought back:** `c1_ucec_plex_reread_result[_<tag>].json`,
   `c1_ucec_plex_reread[_<tag>].npz`, `c1_ucec_plex_reread[_<tag>].log` and the LSF `.out`
   logs of the blind, smoke and full runs; their sha256 are recorded when they arrive. The
   wording of section 4 and the changes of section 6 are then applied exactly as written.

**Blind and smoke record** (filled from the LSF .out before the sidecar is written)

| Step | LSF job ID | Host | `this script sha256` | Outcome |
|---|---|---|---|---|
| `--perms 5 --blind` and `--perms 5 --tag smoke`, attempt 1 | 75737659, 75737660 | not recorded | none (no computation) | exit 127: `bsub` was given the wrapper without `bash`; failed before any computation; no output; not a run |
| `--perms 5 --blind` (attempt 2) | 75737661 | s002 | `7800d3aa...` (pre-final script) | rc = 0; counts only; does not count (script sha differs from the sidecar value) |
| `--perms 5 --tag smoke`, smoke attempt 2 | 75737662 | s006 | `7800d3aa...` (pre-final script) | rc = 0, SMOKE VALIDATED; does not count (script sha differs from the sidecar value) |
| `--perms 5 --tag smoke`, smoke attempt 3 | 75737663 | s006 | `339c2fd6...` (final script) | rc = 0, SMOKE VALIDATED; counts; **analysis host: s006** |
| `--perms 5 --blind --tag blind2`, blind attempt 3 | 75737681 | s006 | `339c2fd6...` (final script) | rc = 0; counts; 2026-10-08T00:31:46Z-00:33:26Z; BLIND: n=138 K=28 n_tested=10153 n_selected(usable by universe)=2536 plex_levels=16 dummies=15 p3_cases_present=2 |

Note. Smoke attempt 1 failed at submission, not at the gate. Smoke attempt 3 repeats attempt 2
because the script was finalised (sha256 `339c2fd6...`) after attempt 2 had run, which the
rule in Status requires; blind attempt 3 repeats blind attempt 2 for the same reason, before
the sidecar was written, and printed the same counts. The analysis host is s006: the first
smoke that passed with the final script ran there, and so did the pre-final smoke and the
final blind. The pre-final rows are kept as run and are not used for any reading. The LSF
`-o` file of blind attempt 3 (`c1_ucec_plex_blind2.out`) also holds, above it, the output of
blind attempt 2, because LSF appends to an existing `-o` file; the two are separated by
their LSF headers. Full script sha256 values are in the LSF .out files, which are brought
back (item 5).
