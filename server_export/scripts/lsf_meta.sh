#!/bin/bash
# Acquisition-metadata census over every SVS header. Header reads only, so this
# is I/O bound rather than CPU bound; 8 workers is plenty.
set -uo pipefail
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
export META_WORKERS="${META_WORKERS:-8}"
echo "=== acquisition meta census: host=$(hostname) date=$(date) ==="
cd /public/home/fjhui/ZW/scripts
/public/home/fjhui/miniconda3/bin/python -u build_acquisition_meta.py
rc=$?
echo "=== done rc=${rc} date=$(date) ==="
exit ${rc}
