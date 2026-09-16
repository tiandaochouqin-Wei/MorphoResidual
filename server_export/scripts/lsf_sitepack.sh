#!/bin/bash
set -uo pipefail
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
CANCER=$1
export MORPHO_BATCH_AXIS="${2:-operator}"
export MORPHO_N_PERM="${MORPHO_N_PERM:-1000}"
export MORPHO_N_WORKERS="${MORPHO_N_WORKERS:-8}"
echo "=== sitepack: ${CANCER} axis=${MORPHO_BATCH_AXIS} perms=${MORPHO_N_PERM} host=$(hostname) date=$(date) ==="
cd /public/home/fjhui/ZW/scripts
/public/home/fjhui/miniconda3/bin/python -u residual_analysis_sitepack.py "${CANCER}"
rc=$?
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
