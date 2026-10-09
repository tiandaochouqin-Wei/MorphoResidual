# Amendment 1 to POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md: genes without purity-matched patients in the purity variant

**Written** after the CCRCC run completed and the UCEC run failed, and before any further run
of `local_theta.py` on project data; drafted 2026-10-09 from about 18:50Z (UTC; the file date
2026-10-10 is the local date, UTC+8), and frozen at the UTC time of the `# written` line of its
sidecar (`POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md.sha256`), which lists the sha256 of every
file named below and of every file of the original sidecar. Amendment and sidecar are released
publicly on GitHub before any further run. The pre-specification and its sidecar
(`POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md.sha256`, `# written 2026-10-08T01:04:28Z`, released
in commit `bc66e6b` at 2026-10-08T01:31:27Z) are unchanged. `M/` below is
`review/MEAS_ERROR_2026-10-08/`.

**Why this amendment exists.** The UCEC run stopped on a deterministic defect of the frozen
`local_theta.py`, not on anything in the data. Section 10's failed-run rule would repeat it
once, unchanged; an unchanged repeat runs the same code on the same bytes and fails
identically, after which UCEC would be "not completed". The same code path is taken by the
LUAD and PDAC runs, whose purity files also lack patients (LUAD 103 of 105, PDAC 97 of 137;
F2), so they are expected to fail the same way (below). Following the frozen rules unchanged
would therefore remove three of the five cohorts for a reason unrelated to the data. This
amendment corrects the defect with two code lines and changes nothing else.

## Evidence

All runs ran locally on host `KIRC-CZC5497R6K`, the machine that wrote the sidecar, from `M/`.
Console logs are in `M/run_logs/` and are hashed in this amendment's sidecar.

- **F1 (published quantities only) passed for all five cohorts** (`F1_<c>.log`, 2026-10-08
  01:33:55Z to 01:35:48Z): each exit 0 with "REPRODUCES the published per-gene results
  (<=1e-10)", 0 of 20 p-value mismatches, and ICC difference 0.0 in CCRCC and UCEC.
- **CCRCC completed** (`RUN_ccrcc.log`): started 2026-10-08T01:37:06Z, ended 01:53:40Z,
  exit 0, with the frozen command line of section 10 steps 2-3, and wrote
  `theta_ccrcc_reading.csv` and the other outputs to `export/ccrcc/local_theta/` (file names
  listed, not opened).
- **UCEC failed** (`RUN_ucec.log`): started 2026-10-08T01:53:40Z, ended 02:09:17Z, exit 1.
  Command:

  ```
  python local_theta.py --inputs export/ucec/inputs_ucec.npz --B 200 --seed 20261007 --boot 200 --boot-seed 20261008 --sets-file side/families_reactome2022.tsv --save-null --workers 8 --strata operator,none,plex --primary-arm operator --variants spline,rnapc,purity --purity-file side/purity_ucec.tsv
  ```

  It ran the bootstrap, the three null arms, phi_op and the spline and RNA-PC variants, and
  stopped at its first step of the purity variant (log time 02:09:08Z), in the observed
  draws (`analyse`, `local_theta.py` line 978, `vo = run_draws(...)`), inside a pool worker:
  `_job` line 365 `fid = tk.fold_ids_for(nv)`, `theta_kernel.py` line 54 (`fold_ids_for`)
  and line 47 (`folds_for`), sklearn `_split.py` line 404: "ValueError: Cannot have number of
  splits n_splits=5 greater than the number of samples: n_samples=0." The only warnings are
  three numpy `RuntimeWarning: All-NaN slice encountered` (`_nanfunctions_impl.py` line 1617),
  logged after the third null arm and before phi_op, i.e. while the set rows of the three
  arms were assembled; they are unrelated to the failure.
  The run left three partial outputs in `export/ucec/local_theta/`, written at 02:04:22Z to
  02:04:31Z when the null arms finished: `theta_ucec_null_operator.npz`,
  `theta_ucec_null_none.npz` and `theta_ucec_null_plex.npz` (beside F1's
  `repro_only_ucec.json`). They were not opened; they are kept unchanged, are not used, and
  are hashed into the results record. No reading file was written.
- **Cause.** `purity_residualised()` sets to NaN, in both mRNA and protein, every gene with
  fewer than 30 patients having both protein and a purity value; its docstring says such
  patients "leave the gene's fit set" (section 6, item 2 of the pre-specification: "patients
  without one leave the fit set"), and for these genes the fit set is empty. UCEC has
  purity for 96 of its 100 patients (F2). `_job()` then computes `nv = 0` and calls
  `tk.fold_ids_for(nv)` before any branch on the statistic kind, and sklearn's `KFold`
  refuses zero samples. Downstream, `variant_rows()` already drops every gene whose variant
  statistics are not finite (`local_theta.py` line 800), so a gene returned without a
  statistic is handled; only the early `fold_ids_for` call fails. The exports keep only genes
  with at least 30 protein-quantified patients (`measurement_error_checks.py` v2, header), so
  on every other path `nv >= 30`; `nv = 0` arises only in the purity variant.
- **Which rows can be affected.** The base sets and their depth tertiles hold complete genes
  only (n_fit = n), whose purity-matched patient count equals the number of patients with
  purity (UCEC 96, LUAD 103, PDAC 97; CCRCC 103), at least 30, so none of them is dropped and
  their purity rows are unaffected. Only genes of the missingness strata (`|miss*`) can fall
  below 30 purity-matched patients, so only those strata can lose genes or their purity row.
- **Scope.** CCRCC has purity for all 103 patients, so no gene lost a patient and the run
  completed. GBM runs no purity variant. LUAD and PDAC would fail in the same way as soon as
  one analysed gene has fewer than 30 purity-matched patients; with 2 and 40 patients without
  purity, and the missingness strata holding genes near the 30-patient floor, this is
  expected. It was not checked, because checking would mean running the frozen script on
  their data.

## The change

`M/local_theta_a1.py` is a byte copy of the frozen `M/local_theta.py` (sha256 `2d10449a...`,
LF line endings, kept) with exactly two code changes and a header block at the top of its
docstring that names this amendment and the two changes:

```
@@ -94,7 +104,7 @@
 
 PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))          # MorphoResidual_paper/
 PRESPEC_REL = "review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md"
-SIDECAR = os.path.join(PAPER, "review", "POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md.sha256")
+SIDECAR = os.path.join(PAPER, "review", "POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md.sha256")
 REPRO_TOL = 1e-10
 BATCH_TOL = 1e-9
 MIN_SET = 5          # minimum set size for any z
@@ -361,6 +371,8 @@
     kind = _S["kind"]
     v = ~np.isnan(_S["prot"][:, gi])
     nv = int(v.sum())
+    if nv == 0:
+        return gi, None
     M, P = _S["mrna"][v, gi], _S["prot"][v, gi]
     fid = tk.fold_ids_for(nv)
     if kind == "stats":
```

1. A gene with no fit patients returns no statistic, as `purity_residualised()` intends. The
   gene is then NaN in the observed and null arrays of the purity variant and is dropped from
   every set sum of that variant by `variant_rows()`; a purity row with fewer than 5 remaining
   genes is not written (the existing `MIN_SET` rule).
2. The fail-closed gate of `local_theta_a1.py` reads this amendment's sidecar. `PRESPEC_REL`
   is unchanged, so the gate still refuses unless the original pre-specification is listed;
   the amendment sidecar lists it, the original sidecar, every file of the original sidecar
   and `local_theta_a1.py`.

`theta_kernel.py`, the frozen `local_theta.py` and every other frozen file are unchanged.

## Proof

**Synthetic tests** (`M/test_a1_guard.py`, console log `M/test_a1_guard.log`; inputs built
with `test_local_theta.py`'s own generator, scenario "purity", n = 120, G = 60, all flagged
`is_synthetic=True`; protein made missing so that 15 genes have 29 purity-matched patients
and the others 40 when purity is given for 40 patients; 33 checks, all PASS):

- T0: line diff of `local_theta_a1.py` against `local_theta.py`: the header block, the
  `SIDECAR` line and the two guard lines, nothing else; no CR byte.
- T1: with the partial purity file, the frozen script crashes with the same KFold
  `n_samples=0` error from `_job` at `tk.fold_ids_for(nv)`, both with `--variants purity` and
  with the production shape `--variants spline,rnapc,purity`; the latter leaves exactly the
  three `theta_syn_null_*.npz` files and no reading file, the UCEC pattern. `local_theta_a1.py`
  completes both; workers 2 and workers 1 give identical outputs. At gene level, the a1
  purity variant gives no statistic for exactly the 15 genes below 30 matched patients and
  finite statistics for all others, and `purity_residualised()` leaves every gene with either
  0 or at least 30 fit patients.
- T2: with a complete purity file, frozen and amended runs (`--variants
  spline,rnapc,purity`) write identical outputs: every `theta_*` csv byte-identical, every
  npz array-identical, the manifests identical apart from their timestamps, script hash and
  `--out`.
- T3: with the partial purity file, the amended run's purity rows use exactly the set genes
  with at least 30 matched patients (for example `all|miss5-20`: 15 of 30), a set left with
  none of them has no purity row, and every other output equals a frozen run with
  `--variants spline,rnapc`: sets, pergene, hetero csv byte-identical, null and bootstrap npz
  array-identical, `reading.csv` identical once its `purity_*` columns are dropped,
  `variants.csv` identical once its purity rows are dropped.
- T4: on a synthetic file not flagged synthetic, `local_theta_a1.py` refuses and names the
  amendment sidecar, which did not exist at test time.
- T5: `test_local_theta.py --quick` pointed at `local_theta_a1.py` (a scratch copy with four
  literal substitutions: its folder, the import, the script path in the gate test and the
  script-name filter of check G1g) ends "ALL CHECKS PASSED" (79 checks). The frozen suite on
  the frozen script fails only G1j ("production SIDECAR constant absent"), which fails by
  construction now that the original sidecar exists; every other check line and every printed
  number of the two runs is identical.

**CCRCC stands as run.** The guard changes behaviour only when `nv == 0`, where the frozen
code raises and cannot return. The frozen CCRCC run completed (exit 0, reading file
written), so no call on any path reached `nv = 0`, and every other line it executed is
byte-identical in `local_theta_a1.py`: the amended code gives CCRCC identical numbers. CCRCC
is not repeated.

## What changes in the run order (section 10)

1. Write the amendment sidecar from `MorphoResidual_paper/`, from the list
   `M/SIDECAR_FILES_A1.txt`:

   ```
   { echo "# sha256 of the frozen files for POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md (paths relative to MorphoResidual_paper/)";
     echo "# written $(date -u +%Y-%m-%dT%H:%M:%SZ)";
     grep -v '^#' review/MEAS_ERROR_2026-10-08/SIDECAR_FILES_A1.txt | grep -v '^[[:space:]]*$' | xargs sha256sum -b; } \
     > review/POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md.sha256
   ```

   It must end without a "No such file" message. Release the amendment, its sidecar and every
   listed file on GitHub under section 10's release rule, with these additions, which are
   hashed in the sidecar but withheld from this release: the two console logs
   `run_logs/RUN_ccrcc.log` and `run_logs/RUN_ucec.log`, and the outputs already written in
   `export/<c>/local_theta/` (CCRCC's complete outputs, UCEC's three partial null files and
   every cohort's F1 output `repro_only_<c>.json`). They hold results or partial results; the
   logs, csv and json files are released with the results record, the `.npz` files under section 10's
   rule for derived data (deposited at acceptance). Record the commit and time.
2. **UCEC** is rerun with `local_theta_a1.py`, the frozen UCEC command line above unchanged,
   and `--out export/ucec/local_theta_rerun1`. This is the one repeat allowed by the
   failed-run rule, made with the amended script instead of unchanged.
3. `M/export/lambda_transfer.txt` is written per section 5 (availability rule included) from
   `export/ccrcc/local_theta/theta_ccrcc_reading.csv` and
   `export/ucec/local_theta_rerun1/theta_ucec_reading.csv`.
4. **LUAD, GBM and PDAC** are run with `local_theta_a1.py` and their frozen command lines of
   section 10 steps 6-8 (default `--out export/<c>/local_theta`).
5. Console logs go to new files (`run_logs/RUN_ucec_rerun1.log`, `RUN_luad.log`,
   `RUN_gbm.log`, `RUN_pdac.log`); the logs listed in the amendment sidecar are never
   overwritten, since any change to them makes the gate refuse.

Every other argument, seed, rule, reading and wording stands as frozen. Superseded frozen
text: in section 10 steps 4 and 6-8 and in `README_commands.txt` step 7, the script for UCEC,
LUAD, GBM and PDAC is `local_theta_a1.py` (the gate of which reads this amendment's sidecar),
and UCEC's output folder is `export/ucec/local_theta_rerun1`; in the failed-run rule,
"repeated once, unchanged" becomes, for UCEC, "repeated once with `local_theta_a1.py`"; in
section 10's release rule, "Everything else in the list is released" does not apply to the
files withheld in step 1 above, and the release commit message lists them with the other
withheld paths.

## Limits

- If the UCEC rerun (`local_theta_rerun1`) fails, UCEC is reported as not completed ("not
  run: not completed", one supplement row) and contributes to no transfer. LUAD, GBM and PDAC
  keep section 10's failed-run rule unchanged: a failed or interrupted run is repeated once,
  unchanged, with `local_theta_a1.py` and `--out export/<c>/local_theta_rerun1`, and if that
  fails the cohort is not completed. No further amendment is made for any cohort. A
  `REFUSED:` exit before `theta_<c>_reading.csv` is written (section 10: sidecar mismatch,
  wrong path or argument) remains not a run, as frozen.
- If the UCEC rerun fails, section 5's availability rule applies and lambda_T comes from
  CCRCC only, flagged "transferred from CCRCC only".
- The amendment does not touch the readings: in a cohort where some genes lack purity, the
  purity sensitivity is computed on the genes with at least 30 purity-matched patients, as the
  frozen `purity_residualised()` already specified; the reported `purity_n_genes` gives the
  number. A missingness stratum with fewer than 5 genes having at least 30 purity-matched
  patients has no purity value, and its purity cell reads "not available (fewer than 30
  purity-matched patients)"; a purity theta with `purity_n_genes` below 20 is read as not
  readable (`MIN_THETA`).

## Blindness

No alignment statistic or reading number of UCEC, LUAD, GBM or PDAC has been seen. The five
F1 logs (published quantities only) were read in full. From the failed UCEC run were read:
its start and end lines, command line, warnings, traceback and exit line, the times and
stage tags (such as `[arm plex]`) of its progress lines, and, by the reviewer of this draft,
its lines 1-40 and 141-185 (set sizes, selection strengths, MNAR flags, permutation-strata
information, progress lines and the purity-coverage line); none holds an alignment
statistic or a reading number. From `RUN_ccrcc.log` were read its start line, command line,
end line with the exit code, its purity-coverage line ("purity available for 103/103
patients"); and when the end of that log was first checked, the last row of
the CCRCC reading table it prints (the `pinned_sig&translation` stratum) was displayed. CCRCC
had completed under the frozen rules, and this amendment changes no CCRCC number and no rule.
The partial UCEC outputs were not opened (only their names, sizes and file times were
listed). No CCRCC result file was opened. The synthetic tests read no project data.

## Disclosure

The supplement's measurement-error note (`snote:measerr`) will state that the first UCEC run
failed on a code defect (genes left without purity-matched patients in the purity
sensitivity), and that the defect was corrected by a dated amendment of [UTC date of this
amendment's sidecar], released publicly before UCEC was rerun and before the LUAD, GBM and
PDAC runs; that CCRCC ran with the original script, which gives identical numbers there; and
the purity sensitivity's gene counts where purity was incomplete.

## Notes

- `test_local_theta.py`'s G1j fails on the frozen script now that the original sidecar
  exists, and will fail on `local_theta_a1.py` once the amendment sidecar exists. Section 10's
  gate on another host is needed only if the pinned host is unavailable; it is then F1 plus
  the frozen `test_local_theta.py` (which tests `local_theta.py`, not the amended copy), and
  G1j is the one check expected to fail, by construction, and is recorded as such.
- After the amendment sidecar exists, checks T4a, T4b and T5a of `test_a1_guard.py` fail by
  construction (as G1j does), and the script is never again run with `--log` inside `M/`:
  `test_a1_guard.log` is listed in the sidecar, and rewriting it would make the gate refuse.
  Its default log path is the scratch folder it is given.
- The logged run of `test_a1_guard.py` used the scratch folder `D:/claude/pajinsen/_scratch_meas_a1`
  (outside the repository); its outputs are not released, and the script regenerates them.
