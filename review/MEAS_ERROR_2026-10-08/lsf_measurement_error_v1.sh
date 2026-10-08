#!/bin/bash
# Measurement-error raw material for ONE discovery cohort under LSF.  DRAFT 2026-10-08, not run.
#
# Install (login node):  cp measurement_error_checks.py lsf_measurement_error.sh /public/home/fjhui/ZW/scripts/
#   (the python file must sit next to residual_analysis.py; it adds its own directory and
#    MORPHO_SCRIPTS_DIR to sys.path.  Nothing existing is overwritten: new filenames only.)
#
# Usage (submitted, not run directly):
#   bsub -q smp -n 8 -o meas_<cohort>.out -e meas_<cohort>.err bash lsf_measurement_error.sh <cohort> [extra python args]
#   smoke  (minutes):  ... lsf_measurement_error.sh ucec --max-genes 40 --B 3 --out $ZW/meas_error/smoke_ucec
#   timing on the node (seconds, synthetic noise, no data):  python measurement_error_checks.py --bench --workers 8
#   resume after a wall-clock kill: resubmit the SAME command with --resume appended.
# <cohort> in {ccrcc, luad, ucec, gbm, pdac}.  Default permutation budget B=200, all genes;
# see README_commands.txt for the cheaper "pinned-significant + 2000 random background" option:
#   ... lsf_measurement_error.sh ccrcc --B 1000 --perm-genes sig --perm-bg 2000
#
# Same conventions as lsf_residual.sh / lsf_c1_test.sh: a clean MORPHO_* environment, paths
# from morpho_env.sh (pinned manifests / crosswalks / slide maps), explicit MORPHO_OUT so that
# nothing can land in a cohort's results/ directory, 8 workers (a 16-core span[hosts=1] request
# does not schedule on this cluster's smp queue).
set -uo pipefail
CANCER=$1; shift
ZW=/public/home/fjhui/ZW

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
# bsub copies the submitting shell's environment: a stale MORPHO_OUT from another cohort once
# redirected a whole run (2026-09-15).  Clear before sourcing.
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP MORPHO_PROTEIN_TSV \
      MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV

# pinned significant-gene snapshot (Table 1 counts are asserted inside the script)
export MORPHO_SIG_CSV="${ZW}/scripts/pinned/results_snapshot_20260901/${CANCER}/residual_results_tumoronly.csv"

# fresh output directory, separate from every existing results dir
export MORPHO_OUT="${ZW}/meas_error/${CANCER}"
source "${ZW}/scripts/morpho_env.sh" "${CANCER}" || exit 1
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
export MORPHO_SCRIPTS_DIR="${ZW}/scripts"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "=== meas-error: cohort=${CANCER} host=$(hostname) date=$(date) workers=${MORPHO_N_WORKERS} args=$* ==="
morpho_print_env
if ! morpho_check_env; then
  echo "FATAL: inputs missing, refusing to run"
  exit 1
fi
mkdir -p "${MORPHO_OUT}"
echo "  python        : $(/public/home/fjhui/miniconda3/bin/python -c 'import sys,numpy,sklearn;print(sys.version.split()[0],numpy.__version__,sklearn.__version__)')"
echo

cd "${ZW}/scripts" || exit 1
/public/home/fjhui/miniconda3/bin/python -u measurement_error_checks.py \
    --cohort "${CANCER}" --out "${MORPHO_OUT}" --workers "${MORPHO_N_WORKERS}" "$@"
rc=$?
echo
echo "=== done cohort=${CANCER} rc=${rc} date=$(date) ==="
ls -la "${MORPHO_OUT}"
exit ${rc}
