#!/bin/bash
# Per-gene sex control for one cohort.
#
# Answers "is the morphology increment just sex?" by putting sex into the mRNA
# baseline and re-measuring, over the cohort's own morphology-predictable gene set,
# with the same estimator as the published run.
#
# UCEC is uniformly female and is refused by the Python script, not by this wrapper.
#
# MORPHO_SIG_CSV is set explicitly to the deterministic tumour-only baseline, the
# same convention as lsf_randproj.sh.
set -uo pipefail
CANCER=$1
ZW=/public/home/fjhui/ZW

if [ "${CANCER}" = "ccrcc" ]; then
  export MORPHO_SIG_CSV="${ZW}/results/residual_results_tumoronly.csv"
else
  export MORPHO_SIG_CSV="${ZW}/${CANCER}/results/residual_results_tumoronly.csv"
fi
export MORPHO_COHORT="${CANCER}"
export MORPHO_SEX_TSV="${ZW}/scripts/sex_by_case.tsv"

source "${ZW}/scripts/morpho_env.sh" "${CANCER}" || exit 1

# after the source: MORPHO_OUT does not exist before it, and set -u would abort.
# Results go to the cohort's own results directory, as randproj_control does, rather
# than into the shared scripts directory.
export MORPHO_SEX_OUT="${MORPHO_OUT:-${ZW}/scripts}/sex_control_${CANCER}.csv"

echo "=== sex control: ${CANCER} host=$(hostname) date=$(date) ==="
morpho_print_env
if ! morpho_check_env; then
  echo "FATAL: inputs missing, refusing to run"
  exit 1
fi
if [ ! -f "${MORPHO_SEX_TSV}" ]; then
  echo "FATAL: ${MORPHO_SEX_TSV} not found -- upload sex_by_case.tsv first"
  exit 1
fi
echo

cd "${ZW}/scripts"
/public/home/fjhui/miniconda3/bin/python -u sex_control.py
rc=$?
echo
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
