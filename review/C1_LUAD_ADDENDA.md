# C1-LUAD post-signature addenda

Required by the signed rule (`review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md`, sha256
`a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022`), §7.4: "Post-signature
entries go in a separate dated file, `review/C1_LUAD_ADDENDA.md`, each entry hashed and
timestamped as in D3. These are the §5.6 unblinding time and any §5.9 addenda." The signed
file itself is not edited.

**Timestamp status:** entries 1-8 were written on 2026-10-02 and are released, with this
file's sha256 sidecar, in the public GitHub release `c1-luad-addenda-2026-10-02`
(the corresponding author authorised the push in writing on 2026-10-02). The release time
recorded by GitHub is the third-party timestamp for these entries. Later entries go in a
later release.

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
`c1_luad_pergene_result_pergene_{b1,b2}.json` and the per-gene CSVs). **To do:** copy them to
`server_export/pinned/luad_c1/` and record their sha256 here.

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

**Result:** (to be appended in a later entry).

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

## Entry 6 (2026-10-02): reporting gaps in the frozen code, filled from saved outputs

The rule asks for two quantities that the frozen scripts do not print; they are computed
afterwards from saved outputs, with no new permutation, and labelled as such:
- D1 / §3: "the fraction of positive increments per stratum for selected and background
  genes" -- the script records only the stratum-wide fraction (as found for C1-UCEC, §7.4).
  Source: `c1_luad_primary_gene_level_primary.csv`.
- §5.4: bootstrap SD is in the bootstrap JSON but was not printed; read from
  `c1_luad_bootstrap_result_bootstrap.json`.
