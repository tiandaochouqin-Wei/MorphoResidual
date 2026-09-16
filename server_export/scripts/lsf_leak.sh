#!/bin/bash
set -uo pipefail
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
CANCER=$1
shift
echo "=== batch leak check: cancer=${CANCER} host=$(hostname) date=$(date) ==="
cd /public/home/fjhui/ZW/scripts
/public/home/fjhui/miniconda3/bin/python -u batch_leak_check.py "${CANCER}" "$@"
rc=$?
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
