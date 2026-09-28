#!/bin/bash
# Per-gene demographic controls (sex, age, race, and all three together) for one cohort,
# each with a dimension-matched Gaussian-noise counterpart so the "kept %" can be read.
#
# Supersedes lsf_sex.sh: same estimand, sex included, plus the capacity control.
# Which blocks run is decided by the Python script from the data, not here: UCEC has no
# sex variation, and CCRCC/UCEC have too few non-white patients for a race adjustment.
#
# MORPHO_SIG_CSV is set explicitly to the deterministic tumour-only baseline, the same
# convention as lsf_randproj.sh and lsf_sex.sh.
set -uo pipefail
CANCER=$1
ZW=/public/home/fjhui/ZW

if [ "${CANCER}" = "ccrcc" ]; then
  export MORPHO_SIG_CSV="${ZW}/results/residual_results_tumoronly.csv"
else
  export MORPHO_SIG_CSV="${ZW}/${CANCER}/results/residual_results_tumoronly.csv"
fi
export MORPHO_COHORT="${CANCER}"
export MORPHO_COVAR_TSV="${ZW}/scripts/covariates_by_case.tsv"

# morpho_env.sh sets MORPHO_OUT only if it is not already set, so a value inherited from
# an interactive shell that sourced a DIFFERENT cohort earlier would silently redirect the
# output. Clear it first; the data paths are set unconditionally and are not affected.
unset MORPHO_OUT

source "${ZW}/scripts/morpho_env.sh" "${CANCER}" || exit 1

# after the source: MORPHO_OUT does not exist before it, and set -u would abort
export MORPHO_COVAR_OUT="${MORPHO_OUT:-${ZW}/scripts}/covariate_control_${CANCER}.csv"

echo "=== covariate control: ${CANCER} host=$(hostname) date=$(date) ==="
morpho_print_env
if ! morpho_check_env; then
  echo "FATAL: inputs missing, refusing to run"
  exit 1
fi
if [ ! -f "${MORPHO_COVAR_TSV}" ]; then
  echo "FATAL: ${MORPHO_COVAR_TSV} not found -- upload covariates_by_case.tsv first"
  exit 1
fi
echo

cd "${ZW}/scripts"
/public/home/fjhui/miniconda3/bin/python -u covariate_control.py
rc=$?
echo
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
