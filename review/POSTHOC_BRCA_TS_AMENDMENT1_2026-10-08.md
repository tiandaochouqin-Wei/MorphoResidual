# Amendment 1 to POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md: the diagnostic-slide embedding folder

**Written** after the server checks of `RUNBOOK.sh` steps 4-5c and before any residual-analysis
job (step 6 onward) was submitted; drafted 2026-10-08 from about 15:40Z (UTC), and frozen at
the UTC time of the `# written` line of its sidecar
(`POSTHOC_BRCA_TS_AMENDMENT1_2026-10-08.md.sha256`), which also lists the sha256 of every
evidence file named below. Amendment and sidecar are released publicly on GitHub before step 6.
The pre-specification and its sidecar (`POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256`, sha256
`abad8575...`) are unchanged.

**Why this amendment exists.** Section 1 of the pre-specification (DX arm, "Required server
checks"), section 4, item 4 and the step 5b DECISION of `RUNBOOK.sh` require that, before
step 6, the server confirm that the DX embedding folder the re-run reads (`$ROOT/emb_phikon`)
is the folder the published DX job read and that its tiles are at stride 0, and they require a
dated amendment, released before step 6, if anything differs ("for example to set
`DX_EMB_DIR` for the dx job"). Something differed.

## Evidence

Cluster time is UTC+8. The frozen step-5b block was run from `step5b_checks.sh` (a four-line
header, a shebang and three comment lines, followed by a byte-for-byte copy of lines 24-26
and 112-143 of the frozen `RUNBOOK.sh`; sha256 `d4b16b30...`), which appends to
`$ROOT/prerun_checks.txt`; that file was copied to `prerun_checks_step5b.txt` before step 6.
The read-only commands of section "Further checks" below were saved as
`$ROOT/step5b_extra_checks.txt`. Both files, `deployed_scripts_ts_arm.sha256`,
`frozen_scripts.sha256` and the two published LSF output files are kept in
`review/BRCA_TS_ARM_2026-10-08/step5b_evidence/` and hashed in this amendment's sidecar.

- **Deployed scripts (M4).** `sha256sum -c $ROOT/frozen_scripts.sha256` printed OK for all
  eight entries (the seven deployed scripts, including `cptac2_replicate.py`, and the frozen
  189-file list) and FAILED for none. The block ran twice, so `prerun_checks_step5b.txt`
  holds the eight lines twice, all OK.
- **Published script copies (M7 c).** `$S/cptac2_residual.py` hashes (CR-normalised) to
  `1e12432f...` and `$S/extract_phikon.py` to `4c616526...`, the local published copies.
- **Embedding folders and density (M7 a)**, from the saved tile coordinates of the same five
  slides in each folder:
  - `$ROOT/emb_phikon` (folder time 2026-09-11 14:24 +0800): x-step **512** on all five
    (8,041, 10,925, 10,974, 17,938 and 9,253 tiles), i.e. stride 512, the sparse extraction.
  - `$ROOT/emb_phikon_dense` (folder time 2026-09-14 08:18 +0800): x-step **256** on all
    five (32,146, 43,815, 44,021, 71,637 and 37,066 tiles), i.e. stride 0, full density. The
    dense/sparse tile ratio is 3.99 to 4.01 on the five slides.
  - `$ROOT/emb_phikon_ts`: x-step 256 on all five sampled slides (the frozen-section arm).
- **Which folder the published DX job read (M7 b).** The shell history holds two submission
  lines of `cptac2_residual.py --project TCGA-BRCA` from `$S` (history lines 132 and 184):
  one with the default folder `$ROOT/emb_phikon` (`-o cptac2_brca_residual.out`) and one with
  `MORPHO_WSI_EMB_DIR=$ROOT/emb_phikon_dense` (`-o cptac2_brca_residual_dense.out`).
  - The dense job ran from 2026-09-14 10:13:09 to 14:35:27 +0800 and reports writing
    `results/brca_protein_results.csv` and `results/summary_brca.csv`; both files carry the
    time 2026-09-14 14:35:18 +0800, inside that window, and nothing wrote them afterwards.
  - The sparse run's outputs are kept as `results/brca_protein_results_sparse.csv` and
    `results/summary_brca_sparse.csv` (file time 2026-09-14 10:13:46 +0800, 37 s after the
    dense job started; they were copied, not moved, before the dense job wrote).
  - `cptac2_brca_residual.out` (file time 2026-09-11 23:44:51 +0800) holds two LSF reports,
    because LSF appends to an `-o` file: a first attempt that ran 06:32:29-06:33:42 on
    2026-09-11 and "Exited with exit code 1" (73 s, before the last change to `emb_phikon` at
    14:24), and a run from 20:23:42 to 23:44:51 the same day, "Successfully completed". This
    bears only on the provenance of the sparse figure, not on which run wrote the published
    files. The dense log holds one report, "Successfully completed".
  - The published script copy `$S/cptac2_residual.py` carries the time 2026-09-11 05:56:19
    +0800, before both runs, so both ran the copy that hashes to `1e12432f...`.
  - The published summary reads n_cases 102, tested 10,031, significant 0,
    negative_frac_all 0.9911275047, which the paper's Methods report as the full-density
    result (99.11%); the sparse summary reads the same counts with negative_frac_all
    0.9912271957, the Methods' 99.12% for the sparse re-extraction.
  - Both folders hold 107 `.pt` files with identical names.

**Conclusion.** The published DX arm read `$ROOT/emb_phikon_dense` (stride 0), not
`$ROOT/emb_phikon` (stride 512) as section 1 of the pre-specification assumed. Nothing
points the other way. The gate is the falsification test of this inference: had the
published files come from the sparse embeddings, the dense re-run would differ from them per
gene by far more than the gate's host tolerance and the gate would read "fail".

The Methods' statement that the published BRCA arm used full tiling density (stride 0) is
correct and needs no correction. The tile counts above (dense/sparse ratio 3.99 to 4.01) also
show that the sparse re-extraction keeps about one tile in four, not "every other tile" as the
Methods say; that is the existing to-do F7 (pre-specification section 5), which this
amendment does not change.

## What changes

1. The dx job of step 6 is submitted with this line instead of `RUNBOOK.sh` line 172 (the
   variable is set inside the job command, as in the published submission, history line 184):

   ```
   bsub -q smp -n 8 -m $HOST -o $ROOT/dx_repro.log -e $ROOT/dx_repro.err "DX_EMB_DIR=/public/home/fjhui/ZW/cptac2_brca/emb_phikon_dense bash $ZW/scripts_ts_arm/lsf_cptac2_ts_residual.sh dx"
   ```

   Within its first minutes, `grep 'dx embeddings:' $ROOT/dx_repro.log` (or `bpeek`) must read
   `dx embeddings: /public/home/fjhui/ZW/cptac2_brca/emb_phikon_dense (107 .pt)`. If it names
   any other folder or count, the job is killed before it writes any result, its log is kept
   as `dx_repro_wrongdir.log`, and it is resubmitted with this line; that is a submission
   error, not a failed run under section 4, item 6. The TS and TS-tile lines (`RUNBOOK.sh`
   lines 187-188) are submitted unchanged except for a `bash` before the wrapper path, without
   `DX_EMB_DIR`.

   Superseded frozen text: pre-specification section 1, DX "What it is" ("`$ROOT/emb_phikon`"
   becomes "`$ROOT/emb_phikon_dense`"); `RUNBOOK.sh` step 5b DECISION ("the published DX job
   read $ROOT/emb_phikon" becomes "read $ROOT/emb_phikon_dense (Amendment 1)"); `RUNBOOK.sh`
   line 172 (becomes the line above). The wrapper is unchanged; its header hook and its
   `dx embeddings:` line are what this check uses. `RUNBOOK.sh` line 35 (step 0) is historical.
2. The DX header census of step 5 (`RUNBOOK.sh` line 101) took its tile counts from
   `$ROOT/emb_phikon` (sparse); its mpp and magnification columns come from the slide headers
   and stand. Before step 6 a second census, which reads headers and tensor shapes only, was
   run: `ts_slide_qc.py --raw-dir $ROOT/slides --emb-dir $ROOT/emb_phikon_dense --out
   $ROOT/results/dx_slide_qc_dense.tsv` (a new file; `dx_slide_qc.tsv` is kept). It found
   the same 107 slides and the same objective-power and mpp counts as the first census. The
   supplement's DX tile counts come from the new file; both files are brought back at step 8.
3. Nothing else changes: the TS and TS-tile jobs, the pinned host rule, the gate (its inputs,
   tiers and thresholds), the sets, the readings, the wording, the Abstract policy and the
   release rule stand as frozen.

## Notes

- `$ZW/scripts_ts_arm/` also holds `cptac2_residual_v1.py`, the superseded patch, against
  `RUNBOOK.sh` step 0's "do NOT copy the *_v1.* files". It is harmless: no job calls it (the
  wrapper checks that `cptac2_residual.py` is the v3 copy), and the M4 check excludes it
  (`grep -v '_v1\.'`). The folder also holds a copy of `step5b_checks.sh`, called by no job.
  Both appear in `deployed_scripts_ts_arm.sha256`.
- Step 5 (header census): the 189 frozen-section slides are all at 40x (mpp 0.248 for 150,
  0.250 for 39). The 107 diagnostic slides are mixed: apparent magnification 40x for 86 (mpp
  0.245-0.252) and 20x for 21, of which 8 report mpp 0.232 and 13 report 0.499; by pixel
  size, 94 slides are at 0.232-0.252 um/px and 13 at 0.499 um/px. As section 1 fixes, both
  distributions go to the supplement and nothing is restricted or stratified.
- Step 4: all 189 frozen-section slides were embedded (no skipped slide;
  `ts_skipped_slides.txt` is empty). Step 5c (local): 766 of the 777 primary-set genes
  (98.6%) and 98.1% to 99.1% of each secondary set are BRCA-tested; no coverage rule fires.
- Blindness: none of the steps above computes a statistic relating morphology to protein or
  mRNA. They read file names, folder and file times, tile coordinates, slide headers, hashes,
  the gene lists, and the one-row published summaries `results/summary_brca.csv` and
  `results/summary_brca_sparse.csv`. Reading those summaries departs from the "do not open"
  of `RUNBOOK.sh` steps 0 and 8 (lines 40 and 203) for the published DX files; they hold only n_cases, tested, significant,
  pct_sig, median_incr_sig and negative_frac_all, all already in the paper (99.11% and
  99.12%), so they add nothing about any readout. The published per-gene tables
  (`brca_protein_results.csv` and its `_sparse` copy) were not opened; only their file times
  were read.
- Paper: the Supplementary Note on this arm (`snote:brcats`, written after the read-out) will
  state that the diagnostic-slide re-run read the full-density (stride 0) embeddings that the
  published arm read, that the pre-specification had named the sparse folder, and that this
  was corrected by a dated amendment of [UTC date of this amendment's sidecar], released
  publicly before the re-run. The Methods are unchanged apart from the existing F7 item.
