MEAS_ERROR_2026-10-08 -- exact commands (DRAFT; nothing here has been run on project data)

0. Install on the login node (new filenames only, nothing overwritten):
   scp measurement_error_checks.py lsf_measurement_error.sh fjhui@<host>:/public/home/fjhui/ZW/scripts/
   (then on the node: sed -i 's/\r$//' /public/home/fjhui/ZW/scripts/lsf_measurement_error.sh)

1. Calibrate on the target node (seconds, synthetic noise, no data):
   /public/home/fjhui/miniconda3/bin/python /public/home/fjhui/ZW/scripts/measurement_error_checks.py --bench --workers 8
   -> use its "ms per OOF ridge fit" t in:  wall_h = G * B * 2 * t / 3.6e6 / 8

2. Smoke (minutes) -- must print "[depth] max |log2(TPM+1) re-read - RA matrix| = 0.00e+00",
   "[selftest] ... max |diff| = 0.00e+00" and a [repro] MATCHES line on the pinned genes it covers:
   bsub -q smp -n 8 -o meas_smoke_ucec.out -e meas_smoke_ucec.err bash lsf_measurement_error.sh ucec --max-genes 40 --B 3 --out /public/home/fjhui/ZW/meas_error/smoke_ucec

3. Production, one job per cohort (B and gene scope to be fixed by the pre-spec):
   for c in ccrcc luad ucec gbm pdac; do
     bsub -q smp -n 8 -o meas_$c.out -e meas_$c.err bash lsf_measurement_error.sh $c --B 1000
   done
   cheaper variant (pinned-significant genes + 2000 random background genes get the null; all genes get the observed OOF):
     bash lsf_measurement_error.sh $c --B 1000 --perm-genes sig --perm-bg 2000
   observed-only (B=0): no null, minutes.
   resume after a kill: same command + --resume

4. Pull back (per cohort, from /public/home/fjhui/ZW/meas_error/<cohort>/): gene_table_<c>.csv,
   inputs_<c>.npz, oof_pred_<c>.npy, oof_fold_<c>.npy, perm_stats_<c>.npy, perms_<c>.npy,
   patient_table_<c>.csv, replicate_icc_<c>.csv (CCRCC, UCEC), repro_check_<c>.json, run_manifest_<c>.json.
   Do not pull obs_scalars_/config_/done_ (internal checkpoint files).
   Check repro_check_<c>.json says "MATCHES published pipeline" BEFORE looking at anything else.

Local plumbing test (synthetic cohort, ~1 min):  python test_synthetic.py <scratch dir>
