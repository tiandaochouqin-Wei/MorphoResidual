#!/bin/bash
# Measurement-error EXPORT job (v2) for ONE discovery cohort under LSF.  DRAFT 2026-10-08, not run.
# The job exports the inputs needed to compute every statistic LOCALLY (local_theta.py) after the
# pre-specification is frozen.  By default it computes NO morphology-vs-protein statistic other than
# the PUBLISHED incremental R^2 (used to prove the exported inputs reproduce the published results).
#
# Install (login node):  cp measurement_error_checks.py theta_kernel.py lsf_measurement_error.sh /public/home/fjhui/ZW/scripts/
#   (both python files must sit next to residual_analysis.py / batch_leak_check.py; the script adds its
#    own directory and MORPHO_SCRIPTS_DIR to sys.path.  Nothing existing is overwritten: new filenames only.)
#
# Usage (submitted, not run directly):
#   bsub -q smp -n 8 -o meas_<cohort>.out -e meas_<cohort>.err bash lsf_measurement_error.sh <cohort> [extra python args]
#   smoke  (minutes):  ... lsf_measurement_error.sh ucec --max-genes 40 --out $ZW/meas_error/smoke_ucec
#   resume: only the OPTIONAL --B path checkpoints; resubmit the SAME command with --resume appended.
# <cohort> in {ccrcc, luad, ucec, gbm, pdac}.
#
# Exit code: 0 = exported inputs reproduce the published per-gene results to 1e-10 (verdict EXACT);
# 3 = they do not (files are still written; do NOT use them); other = failure.
#
# Same conventions as lsf_residual.sh / lsf_c1_test.sh: a clean MORPHO_* environment, paths from
# morpho_env.sh (pinned manifests / crosswalks / slide maps), explicit MORPHO_OUT so that nothing can
# land in a cohort's results/ directory.
set -uo pipefail
CANCER=$1; shift
ZW=/public/home/fjhui/ZW

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
# bsub copies the submitting shell's environment: a stale MORPHO_OUT from another cohort once
# redirected a whole run (2026-09-15).  Clear before sourcing.
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP MORPHO_PROTEIN_TSV \
      MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV

# pinned significant-gene snapshot = the PUBLISHED per-gene results (Table 1 counts are asserted inside the script)
export MORPHO_SIG_CSV="${ZW}/scripts/pinned/results_snapshot_20260901/${CANCER}/residual_results_tumoronly.csv"

# fresh output directory, separate from every existing results dir
export MORPHO_OUT="${ZW}/meas_error/${CANCER}"
source "${ZW}/scripts/morpho_env.sh" "${CANCER}" || exit 1
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
export MORPHO_SCRIPTS_DIR="${ZW}/scripts"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "=== meas-error EXPORT v2: cohort=${CANCER} host=$(hostname) date=$(date) workers=${MORPHO_N_WORKERS} args=$* ==="
morpho_print_env
if ! morpho_check_env; then
  echo "FATAL: inputs missing, refusing to run"
  exit 1
fi
mkdir -p "${MORPHO_OUT}"
echo "  python        : $(/public/home/fjhui/miniconda3/bin/python -c 'import sys,numpy,sklearn,pandas;print(sys.version.split()[0],numpy.__version__,sklearn.__version__,pandas.__version__)')"
echo

cd "${ZW}/scripts" || exit 1
/public/home/fjhui/miniconda3/bin/python -u measurement_error_checks.py \
    --cohort "${CANCER}" --out "${MORPHO_OUT}" --workers "${MORPHO_N_WORKERS}" "$@"
rc=$?
echo
echo "=== done cohort=${CANCER} rc=${rc} date=$(date) ==="
ls -la "${MORPHO_OUT}"
cat "${MORPHO_OUT}/repro_check_${CANCER}.json" 2>/dev/null
exit ${rc}
