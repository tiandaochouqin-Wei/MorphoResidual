#!/bin/bash
# Operator-design width sweep for one cohort: does the increment survive wider operator
# correction, and is any loss confounding or the price of removing w dimensions?
#
# Uses batch_leak_check.cohort_paths for the data paths, exactly as
# residual_analysis_sitepack.py does, so morpho_env.sh is not needed -- but we still
# clear MORPHO_OUT so an inherited value from an interactive shell cannot redirect the
# result, which is how the covariate run ended up writing into luad/results.
set -uo pipefail
CANCER=$1
ZW=/public/home/fjhui/ZW

if [ "${CANCER}" = "ccrcc" ]; then
  export MORPHO_SIG_CSV="${ZW}/results/residual_results_tumoronly.csv"
else
  export MORPHO_SIG_CSV="${ZW}/${CANCER}/results/residual_results_tumoronly.csv"
fi
export MORPHO_COHORT="${CANCER}"
export MORPHO_WIDTHS="${MORPHO_WIDTHS:-0,12,20,full}"
unset MORPHO_OUT
export MORPHO_OPW_OUT="${ZW}/scripts/operator_width_${CANCER}.csv"

echo "=== operator width sweep: ${CANCER} host=$(hostname) date=$(date) ==="
echo "    widths=${MORPHO_WIDTHS}  sig=${MORPHO_SIG_CSV}"
if [ ! -f "${MORPHO_SIG_CSV}" ]; then
  echo "FATAL: ${MORPHO_SIG_CSV} not found"
  exit 1
fi
echo

cd "${ZW}/scripts"
/public/home/fjhui/miniconda3/bin/python -u operator_width_control.py
rc=$?
echo
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
