#!/bin/bash
# C1-UCEC confirmatory UNI extraction from the IDC DICOM copy. Submit from the login node:
#   bsub -q interactive -gpu "num=1" -o c1_extract_s0.log -e c1_extract_s0.err \
#        /public/home/fjhui/ZW/scripts/lsf_extract_c1.sh --shard 0 --nshards 2
# interactive queue: 24h limit, 2 jobs per user -> two shards. Re-submitting the same
# shard resumes (finished slides are skipped).
#
# HF_TOKEN is read from ~/.hf_token, not inherited from the submitting shell: the
# interactive queue was observed re-sourcing shell rc files inconsistently across
# submissions (HF_TOKEN_len 39/37/14 on three otherwise-identical bsub calls of the
# same --only smoke command on 2026-09-16), so an exported value cannot be trusted to
# survive to the job. `chmod 600 ~/.hf_token` beforehand; only this file is authoritative.
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
if [ ! -s "$HOME/.hf_token" ]; then
  echo "FATAL: $HOME/.hf_token missing or empty" >&2
  exit 1
fi
HF_TOKEN=$(cat "$HOME/.hf_token")
export HF_TOKEN
echo "=== job start: host=$(hostname) date=$(date) HF_TOKEN_len=${#HF_TOKEN} args=$* ==="
nvidia-smi -L
cd /public/home/fjhui/ZW/scripts || exit 1
/public/home/fjhui/miniconda3/bin/python -u extract_c1_dicom.py \
  --dicom-root /public/home/fjhui/ZW/c1_ucec_slides_dicom/cptac_ucec \
  --out-dir /public/home/fjhui/ZW/ucec_c1/WSI/emb "$@"
rc=$?
echo "=== job end rc=${rc} date=$(date) ==="
exit ${rc}
