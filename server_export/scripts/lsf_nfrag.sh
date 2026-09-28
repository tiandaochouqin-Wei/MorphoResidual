#!/bin/bash
# n-fragility: separate fold-assignment noise from genuine sensitivity to n.
#
# Smoke-test through THIS wrapper, not by hand-assembling a python command -- the
# operator-width run failed on all five cohorts because the wrapper's default disagreed
# with the script's and only the script had been exercised:
#     MORPHO_MAX_GENES=8 ./lsf_nfrag.sh luad
set -uo pipefail
CANCER=$1
ZW=/public/home/fjhui/ZW

if [ "${CANCER}" = "ccrcc" ]; then
  export MORPHO_SIG_CSV="${ZW}/results/residual_results_tumoronly.csv"
else
  export MORPHO_SIG_CSV="${ZW}/${CANCER}/results/residual_results_tumoronly.csv"
fi
export MORPHO_COHORT="${CANCER}"
unset MORPHO_OUT
export MORPHO_NFRAG_OUT="${ZW}/scripts/n_fragility_${CANCER}.csv"

echo "=== n-fragility: ${CANCER} host=$(hostname) date=$(date) ==="
echo "    seeds=${MORPHO_NFRAG_SEEDS:-10} draws=${MORPHO_NFRAG_DRAWS:-10} fracs=${MORPHO_NFRAG_FRACS:-0.95,0.90,0.85,0.80}"
echo "    sig=${MORPHO_SIG_CSV}"
if [ ! -f "${MORPHO_SIG_CSV}" ]; then
  echo "FATAL: ${MORPHO_SIG_CSV} not found"
  exit 1
fi
echo

cd "${ZW}/scripts"
/public/home/fjhui/miniconda3/bin/python -u n_fragility.py
rc=$?
echo
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
