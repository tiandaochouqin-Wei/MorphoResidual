#!/bin/bash
# Run residual_analysis.py for one CPTAC cohort under LSF.
# Path resolution lives in morpho_env.sh so the ccRCC-vs-per-cancer layout split
# has exactly one definition.
#
# Usage (submitted, not run directly):
#   bsub -q smp -n 8 -o ... -e ... lsf_residual.sh <cancer>
# For bit-exact replication of a previous run, pin the node: bsub -m <host>.
# svd_solver="full" makes the pipeline deterministic on a given host; across
# hosts LAPACK's SVD takes different code paths and results agree to ~3e-5
# (identical significant set, but not byte-identical).
set -uo pipefail
CANCER=$1
ZW=/public/home/fjhui/ZW

source "${ZW}/scripts/morpho_env.sh" "${CANCER}" || exit 1

echo "=== job start: cancer=${CANCER} host=$(hostname) date=$(date) ==="
echo "  cores visible : $(nproc)"
morpho_print_env
if ! morpho_check_env; then
  echo "FATAL: inputs missing, refusing to run"
  exit 1
fi
echo

mkdir -p "${MORPHO_OUT}"
cd "${ZW}/scripts"
/public/home/fjhui/miniconda3/bin/python -u residual_analysis.py
rc=$?
echo
echo "=== job end: cancer=${CANCER} rc=${rc} date=$(date) ==="
# rc=1 only means the >=50-protein gate was not cleared; that is a scientific
# result, not a job failure. Report it plainly rather than masking it.
exit ${rc}
