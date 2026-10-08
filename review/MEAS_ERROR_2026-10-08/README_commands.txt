Where this file differs from POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md, the pre-specification governs.
MEAS_ERROR_2026-10-08 -- exact commands, v2 (EXPORT job + local analysis).  The export jobs have run (pre-specification, Status); the local analysis has not.
(v1 files are kept as *_v1.* ; v1 stored nested-increment statistics from which the discriminator theta cannot be
 recovered -- DERIVATION.md 2.3 -- so the job is now an export job and every statistic is computed locally.)

FILES
  measurement_error_checks.py   EXPORT job (cluster).  Default = no discriminator is computed on the cluster.
  theta_kernel.py               shared numerics (published ridge, W-only OOF statistics, permutation/bootstrap draws)
  lsf_measurement_error.sh      LSF wrapper
  local_theta.py                LOCAL N, N_c, D, theta, lambda_hat_noise, S_cross (+ nulls, bootstrap CI); FAIL-CLOSED
  test_synthetic.py             plumbing test of the export job on a fake cohort (real RA loaders)
  test_local_theta.py           end-to-end test of local_theta.py on synthetic cohorts (population values, nulls, gate)
  theory/                       DERIVATION.md and the simulations

0. LOCAL PRE-FLIGHT (already green on 2026-10-08; rerun after ANY edit; each must end with "ALL CHECKS PASSED")
   python test_synthetic.py   <scratch dir>              # ~1 min
   python test_local_theta.py <scratch dir>              # ~1 min (add --quick for a 30 s version)

1. INSTALL on the login node (new filenames only, nothing overwritten; the two .py files must sit next to
   residual_analysis.py and batch_leak_check.py):
   scp measurement_error_checks.py theta_kernel.py lsf_measurement_error.sh fjhui@<host>:/public/home/fjhui/ZW/scripts/
   ssh fjhui@<host> "sed -i 's/\r$//' /public/home/fjhui/ZW/scripts/lsf_measurement_error.sh"

2. SMOKE (minutes; 40 genes).  Must print, in this order:
     "[depth] max |log2(TPM+1) re-read - RA matrix| = 0.00e+00"
     "[labels] {... 'operator': 'N/N labelled', 'plex': 'N/N labelled', ...}"   (no "_error" key; if there is one the
                                                                                 stratified local arms are impossible)
     "[selftest] theta_kernel.oof_ridge vs RA.cv_r2: max |diff| = 0.00e+00"
     "[repro] {... 'verdict': 'EXACT (<=1e-10)'}"
   bsub -q smp -n 8 -o meas_smoke_ucec.out -e meas_smoke_ucec.err bash lsf_measurement_error.sh ucec --max-genes 40 --out /public/home/fjhui/ZW/meas_error/smoke_ucec

3. PRODUCTION, one job per cohort (no --B: nothing is computed on the cluster except the published increments):
   for c in ccrcc luad ucec gbm pdac; do
     bsub -q smp -n 8 -o meas_$c.out -e meas_$c.err bash lsf_measurement_error.sh $c
   done
   exit code 0 = EXACT reproduction of the published per-gene results; 3 = does NOT reproduce (stop, do not use).
   Expected runtime per cohort: ESTIMATE 10-30 min, not measured -- the RNA files are read twice (RA loader + the
   raw-count pass), the .pt embeddings once, then about 20,000 ridge fits (~1 min) for the repro; the smoke job gives
   the real figure (extrapolate: the repro part scales with the number of fit genes).
   Expected sizes per cohort (estimate from the synthetic export, ~32 bytes per patient x gene entry before compression;
   check with ls -la): inputs_<c>.npz 25-45 MB, replicate_files_<c>.npz ~10-15 MB (CCRCC, UCEC only), gene_table
   <c>.csv 2-4 MB, the rest < 1 MB.  All five cohorts: roughly 250 MB.

4. BRING BACK (per cohort, from /public/home/fjhui/ZW/meas_error/<cohort>/ into review/MEAS_ERROR_2026-10-08/export/<cohort>/):
   inputs_<c>.npz  gene_table_<c>.csv  patient_table_<c>.csv  repro_check_<c>.json  run_manifest_<c>.json
   export_sha256_<c>.txt  replicate_icc_<c>.csv + replicate_files_<c>.npz (CCRCC, UCEC only; the job says
   "[icc] only N cases ..." and writes neither for the others)
   scp -r fjhui@<host>:/public/home/fjhui/ZW/meas_error/<c>/{inputs,gene_table,patient_table,repro_check,run_manifest,export_sha256,replicate_icc,replicate_files}_<c>.* export/<c>/
   then:  cd export/<c> && sha256sum -c export_sha256_<c>.txt      # every line must say OK (transfer integrity)
   Do NOT pull config_/done_ (internal checkpoints).  repro_check_<c>.json must say "EXACT (<=1e-10)" BEFORE anything else.

5. LOCAL REPRODUCTION CHECK (published quantities only; allowed without the sidecar; computes no discriminator):
   python local_theta.py --inputs export/<c>/inputs_<c>.npz --repro-only --repro-pvals 20 --check-pca --check-icc
   -> exit 0 and "REPRODUCES the published per-gene results (<=1e-10)"; the p-value line must show n_mismatch 0
      (re-derives 20 published permutation p-values with RA's own RandomState(0) sequence, ~1 min);
      PCA diff is informational (PCA is fitted once on ALL common patients in residual_analysis.py, so nothing
      is re-fitted in-fold); ICC diff must be <= 1e-10 for CCRCC/UCEC.

6. FREEZE (section 10 of the pre-specification governs; before the first statistic is computed):
   F2  write M/side/build_side_files.py and run it: side/families_reactome2022.tsv, side/purity_{ccrcc,luad,ucec,pdac}.tsv
       (reads only the GMT, covariate tables and the 'patients' array of each inputs file).
   F3  this file.   F4 (confirmatory C1-LUAD export, enables rule 8): optional, NOT done.
   F5  complete M/SIDECAR_FILES.txt, then from MorphoResidual_paper/ :
         { echo "# ..."; echo "# written $(date -u +%Y-%m-%dT%H:%M:%SZ)";
           grep -v '^#' review/MEAS_ERROR_2026-10-08/SIDECAR_FILES.txt | grep -v '^[[:space:]]*$' | xargs sha256sum -b; } \
           > review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md.sha256
       (every line "<sha256> *<path>"; must end without "No such file"); release on GitHub, record commit and time.
       Any later edit of any listed file makes local_theta.py REFUSE.  Do not normalise line endings before hashing.

7. LOCAL RUN after freezing, from M/, local only, one completed run per cohort, NO smoke run on real data:
   F1  per cohort <c> (ccrcc luad ucec gbm pdac):
       python local_theta.py --inputs export/<c>/inputs_<c>.npz --repro-only --repro-pvals 20 --check-pca --check-icc --out export/<c>/local_theta
       needs exit 0, "REPRODUCES ... (<=1e-10)", p-value mismatches 0, ICC diff <= 1e-10 (CCRCC, UCEC); --check-pca is informational.
       A cohort failing F1 is "not run: F1 failed"; the others proceed.
   COMMON = --B 200 --seed 20261007 --boot 200 --boot-seed 20261008 --sets-file side/families_reactome2022.tsv --save-null --workers 8
   CCRCC: python local_theta.py --inputs export/ccrcc/inputs_ccrcc.npz COMMON --strata operator,none,plex --primary-arm operator \
            --variants spline,rnapc,purity --purity-file side/purity_ccrcc.tsv
   UCEC:  as CCRCC with export/ucec/inputs_ucec.npz and side/purity_ucec.tsv
   then write lambda_T and the three depth-tertile values (pre-spec section 5, with its availability rule) into
        export/lambda_transfer.txt (four decimals, with sources)
   LUAD:  as CCRCC with luad files and side/purity_luad.tsv, plus --lambda-lo-transfer <lambda_T> --lambda-lo-transfer-tertiles <a,b,c>
          (both omitted if neither CCRCC nor UCEC contributes)
   GBM:   as LUAD but --variants spline,rnapc and no --purity-file
   PDAC:  as LUAD (same variants) but --strata operator_scanner,none,plex --primary-arm operator_scanner
          and --purity-file side/purity_pdac.tsv
   C1-LUAD confirmatory (--confirmatory --strata operator,none --primary-arm operator --lambda-lo-transfer <lambda_T>): only if F4 was frozen; it was not.
   Outputs in export/<c>/local_theta/: theta_<c>_{sets,reading,hetero,variants,pergene,manifest}, _boot, _null_<arm>.
   Hash every output and the console log of each invocation into a results record; a failed run follows pre-spec section 10
   ("Failed or interrupted runs"); the rerun tag is --out export/<c>/local_theta_rerun1.  A completed run is never repeated with other settings.
   Cost (laptop estimate): about 5 min per cohort on 8 workers for 3 arms + bootstrap; nulls ~30 MB per arm.  --secondary and --max-genes are not used.

8. OPTIONAL (speed fallback only; computes the discriminator ON THE CLUSTER, hence refused unless the frozen pre-spec
   hash is supplied and recorded):
   bash lsf_measurement_error.sh <c> --B 200 --prespec-sha256 $(sha256sum POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md | cut -d' ' -f1)
   writes perms_<c>.npy, wonly_obs_<c>.npy (G,5), wonly_perm_<c>.npy (G,B,4) [N,Nc,D,S_cross]; --perm-genes sig --perm-bg 2000
   restricts the null to the pinned set + random background genes; --resume after a kill.  The on-cluster path has the
   unstratified null only; the stratified arms exist only locally.  Timing on the node: python measurement_error_checks.py --bench --workers 8

WHAT THE TESTS ESTABLISH (synthetic data only; DERIVATION.md section numbers)
  - batched/weighted ridge == the published ridge (<=1e-15); bootstrap weights == duplicated rows; local statistics ==
    theory/sim_signtest.gene_stats on the same folds (<=1e-15); within-batch permutations == residual_analysis_sitepack.py's.
  - population values (n=4000, 40 genes, 20% missing protein): theta_set +1.50 / +4.10 (noise, lam .6/.8; population 1.50/4.00),
    -0.003 (post), -1.01 (err), +0.455 / +0.496 (mix f_X .3/.5; population .45/.50); S_cross sign flips where theta crosses R2_rna.
  - n=100, K=20, 120 genes, B=200 (30 seeds each, unselected genes): N_c,set z = +33 +- 1.7 (noise; 100 % of seeds p<0.05),
    -30 (err), and under H_post z = +0.54 +- 1.07 with two-sided p<0.05 in 10 % (unstratified arm) / 6.7 % (within-operator
    arm) of seeds -- the set-level test is MILDLY LIBERAL under H_post (small positive offset), so the "indistinguishable
    from null" threshold was calibrated on this; the frozen R1 is -3 < z_N < 2 (pre-specification, section 7).  D_set z ~ 11-29.
    Per-gene share of N_c outside the own-null 2.5/97.5 % band under H_post: 3.3 % / 3.4 % (null expectation 2.5 % each).
  - the exported file reproduces RA's published per-gene results exactly (0.0 difference) and the job exits 3 on a perturbed CSV.
