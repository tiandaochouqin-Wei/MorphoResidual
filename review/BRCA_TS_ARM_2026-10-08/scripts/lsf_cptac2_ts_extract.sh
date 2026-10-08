#!/bin/bash
# BRCA TS arm: Phikon, dense stride (--stride 0 = tile 256 px at level 0, same as the published DX arm).
# Submit from the login node (one job per shard; interactive queue = 24 h limit, 2 jobs/user, same as lsf_extract_c1.sh):
#   bsub -q interactive -gpu "num=1" -o ts_ext0.log -e ts_ext0.err \
#        /public/home/fjhui/ZW/scripts_ts_arm/lsf_cptac2_ts_extract.sh 0/2
#   bsub -q interactive -gpu "num=1" -o ts_ext1.log -e ts_ext1.err \
#        /public/home/fjhui/ZW/scripts_ts_arm/lsf_cptac2_ts_extract.sh 1/2
# Smoke test first (3 slides, prints tiles + wall time per slide -> calibrate the GPU-hour estimate):
#   bsub -q interactive -gpu "num=1" -o ts_smoke.log -e ts_smoke.err \
#        /public/home/fjhui/ZW/scripts_ts_arm/lsf_cptac2_ts_extract.sh 0/1 --limit 3
# Re-submitting the same shard resumes (finished slides are skipped by extract_phikon.py).
# Extra args after the shard are passed to extract_phikon.py.
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
# owkin/phikon must already be in ~/.cache/huggingface (it is, from the DX arm / RPPA arms); GPU nodes may have no internet.
export HF_HUB_OFFLINE=1
SHARD=${1:?usage: lsf_cptac2_ts_extract.sh <i/n> [extra extract_phikon.py args]}
shift
ROOT=/public/home/fjhui/ZW/cptac2_brca
echo "=== job start: host=$(hostname) date=$(date) shard=${SHARD} args=$* ==="
nvidia-smi -L
cd /public/home/fjhui/ZW/scripts_ts_arm || exit 1   # patched copies; published scripts untouched
for f in "${ROOT}/manifest_brca_ts_slides.txt" extract_phikon.py; do
  [ -s "$f" ] || { echo "FATAL: $f missing"; exit 1; }
done
grep -q -- "--manifest" extract_phikon.py || { echo "FATAL: scripts/extract_phikon.py is not the patched copy (no --manifest)"; exit 1; }
/public/home/fjhui/miniconda3/bin/python -u extract_phikon.py \
  --raw-dir "${ROOT}/slides_ts" --out-dir "${ROOT}/emb_phikon_ts" \
  --manifest "${ROOT}/manifest_brca_ts_slides.txt" --stride 0 --shard "${SHARD}" "$@"
rc=$?
echo "=== job end rc=${rc} date=$(date) ==="
exit ${rc}
