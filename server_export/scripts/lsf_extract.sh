#!/bin/bash
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
CANCER=$1
echo "=== job start: cancer=$CANCER host=$(hostname) HF_TOKEN_len=${#HF_TOKEN} ==="
nvidia-smi -L
cd /public/home/fjhui/ZW/scripts
/public/home/fjhui/miniconda3/bin/python extract_features.py \
  --raw-dir /public/home/fjhui/ZW/${CANCER}/WSI/raw \
  --out-dir /public/home/fjhui/ZW/${CANCER}/WSI/emb
