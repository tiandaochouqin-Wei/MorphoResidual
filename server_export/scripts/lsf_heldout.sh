#!/bin/bash
# Held-out effect size on the batch-corrected significant set (see held_out_corrected.py
# for the full rationale and the verified aggregation rule).
#
# Smoke-test through THIS wrapper: MORPHO_MAX_GENES=8 ./lsf_heldout.sh luad
set -uo pipefail
CANCER=$1
ZW=/public/home/fjhui/ZW

if [ "${CANCER}" = "ccrcc" ]; then
  export MORPHO_SIG_CSV="${ZW}/results/residual_results_tumoronly.csv"
  export MORPHO_SITEPACK_CSV="${ZW}/results/sitepack_operator.csv"
else
  export MORPHO_SIG_CSV="${ZW}/${CANCER}/results/residual_results_tumoronly.csv"
  export MORPHO_SITEPACK_CSV="${ZW}/${CANCER}/results/sitepack_operator.csv"
fi
export MORPHO_COHORT="${CANCER}"
unset MORPHO_OUT
export MORPHO_HELDOUT_OUT="${ZW}/scripts/held_out_corrected_${CANCER}.csv"

echo "=== held-out (batch-corrected set): ${CANCER} host=$(hostname) date=$(date) ==="
echo "    sig=${MORPHO_SIG_CSV}"
echo "    sitepack=${MORPHO_SITEPACK_CSV}"
if [ ! -f "${MORPHO_SITEPACK_CSV}" ]; then
  echo "FATAL: ${MORPHO_SITEPACK_CSV} not found"
  exit 1
fi
echo

cd "${ZW}/scripts"
/public/home/fjhui/miniconda3/bin/python -u held_out_corrected.py
rc=$?
echo
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
