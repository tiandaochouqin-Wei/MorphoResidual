# C1-LUAD post-signature addenda

Required by the signed rule (`review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md`, sha256
`a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022`), §7.4: "Post-signature
entries go in a separate dated file, `review/C1_LUAD_ADDENDA.md`, each entry hashed and
timestamped as in D3. These are the §5.6 unblinding time and any §5.9 addenda." The signed
file itself is not edited.

**Timestamp status:** entries 1-8 were written on 2026-10-02 and released in the public
GitHub release `c1-luad-addenda-2026-10-02` (2026-10-02 12:05:36 UTC; the corresponding
author authorised the push in writing on 2026-10-02). The file released there is the one
that hashes to sha256
`88041cd01b0a2aab5056d0825c6c482d75c4c3019ed23c8717192a3e86af6dc2`, and the GitHub
release time is the third-party timestamp for that version; Entry 4 fixed the handling of
the reader check before that check was run, which is what that timestamp establishes.

The present file is a later version and hashes differently; its sidecar carries the
current hash, not the released one. It differs from the released version in three places,
all of them additions made after the reader check: Entry 9, the sha256 table and
closing paragraph added to Entry 1, and the one "Result" line of Entry 4, which is marked
in place. Entries 9 and later carry the timestamp of the next release, not of
`c1-luad-addenda-2026-10-02`.

---

## Entry 1 (2026-10-02): run log of every confirmatory job

All on the HZAU HPC, `smp` queue, from `/public/home/fjhui/ZW/scripts`, through
`lsf_c1_luad_test.sh` (sha256 `9ac035ecf2688ef10f6e83feb96fd67ae776621270692f9d8f1cd30cfaf8f73f`),
every one with `--rule-sha256 a68dfbec...` and the rule hash check logged `hash OK`.
Times are CST (UTC+8) from the LSF logs.

| # | Job | Mode / flags | Host, cores | Start | End | rc | Printed unblinded values? |
|---|---|---|---|---|---|---|---|
| 0 | 75682543 | primary `--perms 5 --tag smoke --blind`, `-n 4` | -- | -- (PEND, killed before start) | -- | -- | no (never ran) |
| 1 | (first `-n 1` smoke) | primary `--perms 5 --tag smoke --blind` | s004, 1 | 2026-10-01 23:49:19 | 23:50:22 | 1 | no: stopped at `FileNotFoundError` for `pinned/slide_type_map_luad_c1.tsv` (Entry 3), after `RNA matrix ... 113 x 59427` and before any statistic |
| 2 | (second `-n 1` smoke) | same | s004, 1 | 2026-10-02 00:04:47 | 00:06:25 | 0 | no: `BLIND SMOKE: n=112 K=22 n_tested=10417 n_h1_selected(usable)=2287 n_h2_selected(usable)=721`; "Nothing further is computed or printed in --blind mode (section 5.6)." |
| 3 | 75682618 | primary `--perms 1000 --tag primary` | s002, 4 | 2026-10-02 00:13:47 | 04:55:32 | 0 | **yes: first unblinded run** (Entry 2) |
| 4 | 75683262 | bootstrap `--boots 1000 --h2-claimed --tag bootstrap` | s001, 4 | 2026-10-02 ~05:41 | 08:50:56 | 0 | yes (189.5 min) |
| 5 | 75686865 | sensitivity `--variant s --perms 1000 --tag sens_s` | s003, 4 | 09:50:37 | 10:58:23 | 0 | yes (67.3 min) |
| 6 | 75686867 | sensitivity `--variant s2 --perms 1000 --tag sens_s2` | s004, 4 | 09:50:40 | 10:56:41 | 0 | yes (65.6 min) |
| 7 | 75687748 | descriptive `--which b1 --perms 1000 --tag pergene_b1` | s003, 4 | 11:00:51 | 12:08:58 | 0 | yes (68.0 min) |
| 8 | 75687751 | descriptive `--which b2 --perms 1000 --tag pergene_b2` | s004, 4 | 11:00:53 | 11:43:03 | 0 | yes (42.0 min) |

Each of runs 3-8 was run exactly once. No run was repeated, and no parameter was changed
after any unblinded output. Outputs are under `/public/home/fjhui/ZW/luad_c1/results/`
(`c1_luad_primary_result_primary.json`, `..._gene_level_primary.csv`, `..._primary.npz`,
`c1_luad_bootstrap_result_bootstrap.json`, `c1_luad_sensitivity_{s,s2}_result_sens_{s,s2}.json`,
`c1_luad_pergene_result_pergene_{b1,b2}.json` and the per-gene CSVs). They were copied off the cluster on
2026-10-02, after the reader check of Entry 9, and are held at
`server_export/scripts/pinned/luad_c1/` in the working tree. They are **not** in this
public repository: the set is 78 MB, almost all of it the permutation array
`c1_luad_primary_primary.npz`, and large derived artefacts are deposited separately
(Zenodo at acceptance) rather than committed here. Their sha256 are recorded now so the
deposited files can be checked against this dated entry:

| File | sha256 |
|---|---|
| `c1_luad_primary_result_primary.json` | `b8e89e146de6343d94abd68c9befc0ba5a1bc8038a99b0a1e6dd5c461e34ca75` |
| `c1_luad_primary_gene_level_primary.csv` | `ae23a78219322570d5e3dab050db9495f6c83517273612ee75c786b67afde4f7` |
| `c1_luad_primary_primary.npz` | `2034e71143c2e065675d3159832b1a7be17aa38c7feeaacfc82862c411ea2aab` |
| `c1_luad_bootstrap_result_bootstrap.json` | `f02f11d36923473c9bfac2c9ac52e5133b701213eaf07d0c2ad791618ba69318` |
| `c1_luad_sensitivity_s_result_sens_s.json` | `00d156413a7a5ec794dfe861e9b8532ddd3f9545dcc3c4f42c44813bc4d29bca` |
| `c1_luad_sensitivity_s2_result_sens_s2.json` | `e22393ef65ddca60331818d3feb9527ea96d92efe55fc2289477dd69ecc44648` |
| `c1_luad_pergene_result_pergene_b1.json` | `a53efb879225f3ada9f0f4d32825c20b3777c67f510ff9c9d3e88796ea971216` |
| `c1_luad_pergene_result_pergene_b2.json` | `5f1f66f4af8d72e21c3f4c558ebddbc57dcfe7843bd0ea57f4c4518d4b68f8c7` |
| `c1_luad_pergene_b1_pergene_b1.csv` | `d8dadcf76eb14e9c6d6ef445fd03aa14bf2491ae8c33428ec458834381203c5b` |
| `c1_luad_pergene_b2_pergene_b2.csv` | `52631cab43472bc845f9f02c14d1d0206781b9e83d3a927d5ab91fb28c6c4492` |

The two quantities of Entry 6 were read from these files, with no new permutation: the
bootstrap SDs from the bootstrap JSON (0.024 at step 1, 0.033 at step 3), and the
per-stratum selected/background positive fractions from the gene table (0.270 against
0.139 complete, 0.129 against 0.038 incomplete). Reading the gene table reproduces three
values the frozen run had already written to the JSON — the 85.8% negative-increment
fraction, rho = 0.3978, and the 2,287/8,119 split — which is the check that the table was
read on the same terms the frozen script used.

## Entry 2 (2026-10-02): §5.6 unblinding time

§5.6 requires "The time of the first unblinded output is recorded in the result JSON and in
the post-signature addenda file (§7.4)." The frozen `c1_run_test_luad.py` (sha256
`1d5ee05d57f80af4d90870aecb10e574e4c910a8ce1eab743a067bd9a9278742`) does not write such a
field; this is an implementation gap in the frozen code, not corrected after unblinding
(§5.9: no post-unblinding code change was made).

Reconstruction from the LSF logs: the first unblinded output was produced by run 3 (job
75682618), which started 2026-10-02 00:13:47 CST and wrote its result files at 04:55:32 CST.
The step-1 to step-3 lines are printed during the run, without timestamps. The author first
saw unblinded values during that run, via `bpeek 75682618` (STEP 1-3 lines visible,
completeness block not yet printed); the exact time of that `bpeek` was not logged. The
first unblinded output is therefore bounded to 2026-10-02 00:13:47-04:55:32 CST
(2026-10-01 16:13:47-20:55:32 UTC), after the D3 timestamps (OSF 2026-09-30 22:00:24 UTC,
GitHub release 22:12:54 UTC).

## Entry 3 (2026-10-02): deployment events before the first unblinded output (no code or input change)

1. **Missing pinned file on the HPC.** Smoke run 1 stopped with `FileNotFoundError:
   /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad_c1.tsv`. The file had never been
   copied to the HPC. It was copied from the local mirror, and its HPC sha256
   `d564b6bc5682bc6332aa3ad82387ed300a87d8142bba2d794650281350347d60` was checked equal to the
   value fixed in §7.2 before smoke run 2. Nothing else changed.
2. **Files first copied to the HPC after signature, each sha256-verified against §7.2 before
   use:** `lsf_c1_luad_test.sh` (`9ac035ec...`), `c1_luad_pergene.py` (`505f7b11...`),
   `pinned/c1_luad_batch_labels.tsv` (`0a840533...`), `results/luad__sitepack_operator.csv`
   (`aeef7232...`), the signed rule (`a68dfbec...`, copied by WinSCP under its original
   filename), `pinned/slide_type_map_luad_c1.tsv` (`d564b6bc...`).
3. **Vestigial variable.** `lsf_c1_luad_test.sh` exports `MORPHO_ALIQUOT_XWALK` pointing at
   `${C1}/pinned_aliquot_to_case_tumor_luad_c1.tsv`, a file that does not exist (the pull
   wrote `aliquot_to_case_tumor_luad_c1.tsv`). Neither `c1_run_test_luad.py` nor
   `c1_luad_pergene.py` references the variable; `residual_analysis.py` only reads it into an
   unused module constant. No effect on any output. Not changed (the wrapper is hash-pinned).
4. **Smoke run 0** was submitted with `-n 4`, stayed PEND (no single `smp` host had 4 free
   slots) and was killed before starting; it produced no output.

## Entry 4 (2026-10-02): §1.4 reader QC was not run before extraction (deviation)

§1.6 fixes the order: download (2) -> Reader QC §1.4 (3) -> full extraction (4). Step 3 was
not carried out as specified. Before extraction only the `extract_c1_dicom.py --dry-run`
readability pass was run (135 readable, 0 unreadable; this is the §1.3 Extraction bullet's
dry-run, not the §1.4 pixel-equivalence check). The §1.4 check -- `qc_reader_equivalence.py`
(sha256 `45aea4d84d768f95870cc937079afee947dac7a6b16ea594422ab12fd1e315dc`), openslide-DICOM
vs `wsidicom` on identical level-0 regions of C3L-00444-23, C3L-03717-21 and C3L-02513-23
(arm i), and discovery C3L-00001-21 SVS vs its IDC Leica-pathway copy (arm ii) -- was
omitted by oversight (the extraction commands were prepared treating the dry-run as the
reader gate).

Consequence: the check is no longer outcome-blind. It is to be run now, with the frozen
script and its fixed thresholds, and its verdict reported in Methods and Supplementary Note
whatever it is. If any arm returns DIFFERENT, §1.4's failure handling and §5.9's
post-unblinding branch apply (the §5 reading and the Abstract stay on the pre-change run;
any corrected run is reported alongside, labelled). Concretely: if the two runs fall in
different §5 buckets or give different forms of clause [B], [C] or [D], both readings enter
the Abstract, and "Replicat\*", form 1 of [B] and the pass form of [D] are kept only if both
runs support them.

The check is run by `server_export/scripts/lsf_qc_reader_c1_luad.sh` (sha256
`c211ffc9cd6550c655960d497bbb0ed450bd8a379a75b2445b6a82046f5c6cb7`), which only
orchestrates: it hash-asserts the frozen `qc_reader_equivalence.py` and picks the base
VOLUME instance with the frozen extractor's own `resolve_series()`. Arm (ii) needs the IDC
copy of discovery slide C3L-00001-21 (series
`1.3.6.1.4.1.5962.99.1.250299211.778750893.1640927839051.2.0`), downloaded for this purpose.
This entry, the handling above and the wrapper are released publicly before the check is
submitted.

**Result** (this line alone was appended after the check, and is *not* part of the
2026-10-02 12:05:36 UTC release, where it read "**Result:** (to be appended in a later
entry)."; no other word of this entry is changed)**:** every arm EQUIVALENT; see Entry 9.

## Entry 5 (2026-10-02): §5.9 code corrections

None. No script, input or parameter was changed after the first unblinded output.

## Entry 7 (2026-10-02, after unblinding): §6 "Protein batch" is factually wrong

The signed rule's §6 says: "PDC metadata do not give case-to-TMT-plex membership or the
identity of the reference channel: `studyExperimentalDesign` is run-level only, and
`aliquot_is_ref` is null on all 52 runs. No plex stratum is possible". The second half of
the evidence is true (`aliquot_is_ref` is null), but the conclusion is wrong. The same
`studyExperimentalDesign(pdc_study_id, acceptDUA: true)` query exposes eleven per-channel
object fields, `tmt_126`, `tmt_127n`, `tmt_127c`, `tmt_128n`, `tmt_128c`, `tmt_129n`,
`tmt_129c`, `tmt_130n`, `tmt_130c`, `tmt_131`, `tmt_131c`, each with
`{aliquot_id aliquot_submitter_id aliquot_run_metadata_id}`. They are filled for every run
of PDC000127/153/125/204/270 and PDC000489. Step 0b (`query_log_pdc.tsv`) probed only scalar
field names and never requested them.

Found 2026-10-02 by a metadata-only scoping task (no confirmatory omics value used), with
every query logged in `review/BATCH_AXES_SCOPING_2026-10-02/pdc_query_log.tsv` and raw
responses saved. A second PDC route (`paginatedCasesSamplesAliquots` ->
`aliquot_run_metadata`) agrees on all 1,343 channel entries; the legacy CCRCC
`case_to_plex.tsv` agrees on 110/110 cases. For PDC000489: 113/113 frozen cases have a plex
(27 levels, 4-5 cases each); C3N-02230 and C3N-02240 have the same aliquot measured in two
runs (one in 2020).

Consequences: (i) the manuscript text that repeated the §6 claim is corrected (main
Results, Supplementary Note); (ii) the registered test had no plex stratum and none is
added to its reading; (iii) any plex-stratified re-read of C1-LUAD is post hoc, to be
specified in a dated file before it is run and reported whatever its result;
(iv) `c1_luad_plex_runs_summary.py`'s scope note carries the same error and is left as
written (it documents what Step 0b checked), with this entry as the correction.

## Entry 8 (2026-10-02): what this release adds to the public repository

D3's paper wording says the download and range-read logs are "in the repository". This
release adds:
- the D3 receipt `review/C1_LUAD_D3_RECEIPT_2026-10-01.md`;
- this file and its sha256 sidecar;
- `review/c1_luad_logs/step0b_2026-09-29/`: the Step 0b metadata-query and DICOM range-read
  logs (`query_log_idc.tsv` lists every range read), as copies in which only the Windows
  user name in a local cache path is replaced; its README gives the sha256 of each
  original and copy;
- the C1-LUAD, C1-UCEC D4 and pre-check scripts that were untracked. All 19 of them that
  are hashed in §7.2 match their recorded sha256 (checked before commit); three post-date
  §7.2 and are not hashed there: `download_c1_luad_dicom.py` (sha256 `5fdaa28c...`),
  `lsf_qc_reader_c1_luad.sh` (`c211ffc9...`) and `c1_luad_plex_runs_summary.py`
  (`60d3de7a...`, whose scope note is wrong; see Entry 7).
The post-signature HPC logs (DICOM download, extraction, omics pull, test runs) are on the
compute cluster and follow in a later release.

## Entry 9 (2026-10-02): §1.4 reader QC result, and the staging artefacts the first attempt found

Run as Entry 4 committed, after the public release `c1-luad-addenda-2026-10-02`
(2026-10-02 12:05 UTC) fixed its handling. Both attempts used
`lsf_qc_reader_c1_luad.sh` (`c211ffc9...`), which hash-asserts the frozen
`qc_reader_equivalence.py` (`45aea4d8...`); its thresholds, slide list and the two
readers were identical in both, and no code was changed between them.

**Attempt 1 — job 75698522, 2026-10-02 21:02:02–21:02:50 CST, rc=1.**
- Arm (ii), discovery C3L-00001-21 SVS vs its IDC copy: `VERDICT: EQUIVALENT`, 64/64
  regions bit-identical, pooled mean |Δ| 0.0000, max 0.
- Arm (i), C3L-00444-23, C3L-03717-21 and C3L-02513-23: **no verdict**. The script
  raised `wsidicom.errors.WsiDicomUidDuplicateError` inside `WsiDicom.open()`, before
  reading any region, so it never reached either of its two `VERDICT` statements
  (dimension mismatch, exit 2; threshold failure, exit 1). The exit status of 1 was the
  uncaught exception. §1.4's "a DIFFERENT result stops extraction" branch was therefore
  not entered, and arm (i) supplied no evidence either way.

**Cause.** `wsidicom` enumerates every file in the series directory. The confirmatory
staging tree `/public/home/fjhui/ZW/c1_luad_slides_dicom/` held 483 files named
`<uuid>.dcm<digits>` — partial objects left by the interrupted downloads that
`download_c1_luad_dicom.py` documents (its docstring records that a single
`download_from_selection()` over all 135 series "died twice in a row", which is why it
loops one series at a time). Each carries the SOPInstanceUID of the complete file
beside it, so `wsidicom` sees one identifier on two files and refuses.

**These files were not inputs to the extraction, and this is checkable.**
`extract_c1_dicom.resolve_series()` enumerates `series_dir.glob("*.dcm")` and the
series discovery uses `rglob("*.dcm")`; both require the name to end in `.dcm`, which
`<uuid>.dcm<digits>` does not. Independently, `resolve_series()` raises
`"two base-size VOLUME instances (concatenation?)"` when the two largest VOLUME
instances tie, so a duplicate that *did* match the glob would have stopped the run
rather than been chosen silently. No series raised it; all 135 extracted.

**What the 483 were.** 482 were strictly smaller than the complete file of the same
name; none was an orphan (every one had its complete `.dcm` present). The one
exception, in series C3N-02234, had exactly the same length (108,383,034 bytes) as
`1238c702-1576-4ece-ae02-ee7696e32ef8.dcm` but different content from byte 81,117,779
(74.8% in), and neither tail was zero-filled. Both parse and report the same
SOPInstanceUID, 2,989 frames and 108,045,934 bytes of JPEG pixel data, but only the
`.dcm` file decodes in full, to a (2989, 240, 240, 3) array; the suffixed copy yields
340,416,000 values, about 1,970 frames, and fails. The file the extractor read is
therefore the intact one, which the `*.dcm` glob had already guaranteed.

**Handling.** The 483 were moved, with their directory structure, to
`/public/home/fjhui/ZW/c1_luad_partial_downloads_2026-10-02/` (inventory in its
`_moved_list.txt`; nothing deleted, the move is reversible). This is not a §1.4 "fix":
no decoder, reader, threshold, slide list, scale policy, tiling, tissue filter or
encoder was touched, and no file that is an input to the extraction, or to any reported
analysis, was moved: the moved files entered nothing but `wsidicom`'s enumeration of the
series directory, which is exactly what attempt 1 tripped on and what the move cleared.
All 135 series were then re-resolved with `resolve_series()`: 135 of 135, 0 failures.

**Attempt 2 — job 75699834, 2026-10-02 21:50:50–21:51:45 CST, rc=0,
"every arm EQUIVALENT".** All four comparisons give `VERDICT: EQUIVALENT` with 64/64
regions bit-identical, pooled mean |Δ| 0.0000, 99.9th percentile 0.0, max 0, and mean
signed R/G/B difference +0.0000 on each channel. Arm (i): C3L-00444-23 (level-0
87,647 × 45,770), C3L-03717-21 (105,576 × 69,865), C3L-02513-23 (17,140 × 23,098);
both readers reported the same level-0 dimensions for each. Arm (ii): C3L-00001-21
(25,895 × 23,643). PNG strips are under `/public/home/fjhui/ZW/luad_c1/reader_qc/`.

**Reading.** §1.4's failure handling, and the Entry 4 pre-commitment released at
12:05 UTC, are not triggered: no arm failed the frozen thresholds, so the §5 reading
and the Abstract are unaffected and no clause form changes. The check remains not
outcome-blind, which is disclosed wherever it is reported.

## Entry 10 (2026-10-03): execution of the post hoc batch-axis programme, and a reporting correction its own audit forced

**What ran.** Everything `review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md` specifies, under
the released text plus its two dated amendments and the amendment-2 correction, all
timestamped before the corresponding jobs ran. Sections P and M1: seven sitepack jobs
(plex in five cohorts at the 30-level cap; medium in CCRCC and UCEC), all rc=0, 2026-10-03
00:45–00:55 CST. Sections M2 and Q: 88 arms through `lsf_posthoc_residual.sh` under a
6-concurrent job-group limit, all rc=0 by 07:31 CST; the unmodified control arm R
reproduced every published Sig. count exactly (2,191 / 2,310 / 2,566 / 702 / 1,572).
Section L: job 75702342 on s002 (3 cores), 385.4 min, rc=0; the reproduction gate matched
all six pinned values including p1 = 1/1,001, and P1/P2/P3 all returned p = 0.001
(within-plex null 0.5458±0.0194 against AUC 0.6440; details in the manuscript). One
earlier submission of the 21 CCRCC arms exited immediately with code 126 — the re-uploaded
wrapper had lost its execute bit — and produced no output; it was resubmitted unchanged.
Outputs: `results/posthoc/` under each cohort root and
`/public/home/fjhui/ZW/luad_c1/results/posthoc/c1_luad_plex_reread_result_plex.json`
(+ .npz, .log); sha256 to be recorded here when the files are mirrored into the
repository's pinned tree.

**Mirrored (2026-10-04), at `server_export/pinned/posthoc/` in the repository.** The L
result and log ship in the repository; the permutation array (.npz, 234 MB) is hashed
here and held for the Zenodo deposit with the other large derived artefacts, as are the
seven per-gene sitepack tables. The five arm summary tables ship alongside.

| File | sha256 |
|---|---|
| `c1_luad_plex_reread_result_plex.json` | `903c183bd9aa4672472e395e1541b81f9f7080bf64c7eaa2a9b473c6d734aa31` |
| `c1_luad_plex_reread_plex.log` | `b7ec0d2e022444e779cbbb2982410ecf9b140a270a484c4ae5ee4430bd96edbb` |
| `c1_luad_plex_reread_plex.npz` (not in the repository) | `3128195d4cf1e3d96307fa9b228f64e4f1c1d2553d72b508b14733823f33f208` |
| `posthoc_summary_ccrcc.csv` | `9560d6d73b3288581c5454990b7d5f1af53f53c7038e0e80fda39b26a9378f6f` |
| `posthoc_summary_luad.csv` | `5f983e8fa6e9f0f8592c417ca19f085bdcb994fce492089c62e54cf1725c481d` |
| `posthoc_summary_ucec.csv` | `323990565c78163a8384059d21d8d24a2b590144ada33c7c9f787707752a392c` |
| `posthoc_summary_gbm.csv` | `6805f724298f94f30bea39b140b366a2a8540c54e749a896b47b47790d5411ce` |
| `posthoc_summary_pdac.csv` | `e612922c0ba0bc63cdf4982b6eb3ad3c3de9e412d01150f04f917e0e052af412` |

**Reporting correction (disclosed, not silent).** After the first summary of these
results was drafted, an independent recomputation (three fresh reviewers over the raw
tables and the specification text) verified every reported number and every reading
assignment, and found that the draft reported only part of what the specification
mandates per cohort, with the omissions skewed favourable: the batch-only vs
noise-matched control (which in UCEC-medium shows the FFPE/OCT indicator alone
recovering a median increment of +0.111 over the published significant genes, nearly the
embedding's own +0.118); the Strictest counts (UCEC-medium = 0; CCRCC-medium a
by-construction copy of the corrected count, its 13-patient level being below the
15-case training-fold floor); the specified Jaccard overlaps (a one-sided percentage had
been substituted — e.g. UCEC-plex "50.7% of the operator set kept" where the specified
Jaccard is 0.046); M2's retained column (−0.43 and −0.06 in the two count-indeterminate
cohorts); Q's "patients dropped" (0 in GBM, whose arm Q is therefore identical to R and
its Jaccard of 1.000 an identity); the M3 disclosures owed to LUAD and PDAC
(medium-axis largest level 93.3% and 93.4%, above the 90% bar) and GBM's dual criterion
(one medium level, and k=0); and section L itself, which at that point had not been run
at all. Every mandated quantity but one is now computed, L has run, and the manuscript
carries the set (Table 1 caption, the Results paragraph in sec:c1res, Methods
sec:batchmeth, Supplementary Note snote:ucec, Supplementary Methods sm:census,
Supplementary Table tab:posthoc, and snote:c1luad for L). The one still outstanding is
section P's top-15 family enrichment, which runs through `enrichment_check.py` — the
same script, families and Enrichr libraries that produced the published per-axis table —
after the published operator output it would overwrite is set aside; its result is
appended below when it lands, whatever it shows, and the script's re-run of the
uncorrected baseline doubles as a drift control against the published column.

**Enrichment result (2026-10-03, appended as promised above).** Run on the login node
against Enrichr/Reactome with the script's four signature families. At least one family
sits in the top 15 of every corrected set on both axes: plex 8/13/7/8/4 of 15 terms
in-family (CCRCC/LUAD/UCEC/GBM/PDAC), medium 5 (CCRCC) and 9 (UCEC); UCEC's empty
strictest-medium set has nothing to enrich. The baseline re-run gave 7/12/11/6/3,
matching the published family pattern (the published table's six-family vocabulary adds
folding and ECM, so its counts sit one to three higher in CCRCC, GBM and PDAC). The
published operator output `scripts/enrichment_survival.tsv` (sha256 `6d2291125ba2...`)
was backed up before the runs and restored byte-identical afterwards; the new outputs
are `enrichment_survival_plex.tsv` and `enrichment_survival_medium.tsv` beside it. One M0 quantity, Cramér's V
between medium and operator, was computed locally over the pinned analysed-case lists
with a 500-shuffle null (V = 0.81/0.90/0.92/0.90 in CCRCC/LUAD/UCEC/PDAC, p ≤ 0.004;
GBM one level); those lists are exact for CCRCC and UCEC and two to three patients short
of the arm population elsewhere, which is stated where the number is used.

## Entry 6 (2026-10-02): reporting gaps in the frozen code, filled from saved outputs

The rule asks for two quantities that the frozen scripts do not print; they are computed
afterwards from saved outputs, with no new permutation, and labelled as such:
- D1 / §3: "the fraction of positive increments per stratum for selected and background
  genes" -- the script records only the stratum-wide fraction (as found for C1-UCEC, §7.4).
  Source: `c1_luad_primary_gene_level_primary.csv`.
- §5.4: bootstrap SD is in the bootstrap JSON but was not printed; read from
  `c1_luad_bootstrap_result_bootstrap.json`.
