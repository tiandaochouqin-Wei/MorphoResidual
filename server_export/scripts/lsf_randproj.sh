#!/bin/bash
# Matched-capacity / random-projection control for one cohort.
#
# Answers the reviewer's "the gain is just 20 extra ridge columns" objection by
# comparing the real WSI PCs against 20 columns of pure noise and against 20
# random projections of the raw 1024-d embedding, over the cohort's own
# morphology-predictable gene set.
#
# MORPHO_SIG_CSV is set explicitly to the deterministic tumour-only baseline.
# The script's own default is residual_results_wsibatch.csv, which exists only
# for ccRCC -- leaving it unset would make four of five cohorts die on a missing
# file, and would compare cohorts on different gene sets even where it worked.
set -uo pipefail
CANCER=$1
ZW=/public/home/fjhui/ZW

if [ "${CANCER}" = "ccrcc" ]; then
  export MORPHO_SIG_CSV="${ZW}/results/residual_results_tumoronly.csv"
else
  export MORPHO_SIG_CSV="${ZW}/${CANCER}/results/residual_results_tumoronly.csv"
fi

source "${ZW}/scripts/morpho_env.sh" "${CANCER}" || exit 1

echo "=== randproj control: ${CANCER} host=$(hostname) date=$(date) ==="
morpho_print_env
if ! morpho_check_env; then
  echo "FATAL: inputs missing, refusing to run"
  exit 1
fi
echo

cd "${ZW}/scripts"
/public/home/fjhui/miniconda3/bin/python -u randproj_control.py
rc=$?
echo
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
