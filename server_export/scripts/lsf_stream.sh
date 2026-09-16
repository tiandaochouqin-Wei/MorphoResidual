#!/bin/bash
# 流式抽取的 LSF 包装。用法(在 login02 上):
#   export HF_TOKEN=hf_xxx
#   bsub -q interactive -gpu "num=1" -o stream_gbm.log -e stream_gbm.err \
#        /public/home/fjhui/ZW/scripts/lsf_stream.sh gbm
#
# 为什么走 interactive 队列: gpu 队列 PEND 上万排不上，interactive 队列能收批处理
# 作业、24h 时限、每人限 2 个作业。CCRCC/LUAD/UCEC 都是这么跑通的。
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
CANCER=$1
if [ -z "$CANCER" ]; then echo "用法: lsf_stream.sh <癌种名>"; exit 1; fi

echo "=== job start: cancer=$CANCER host=$(hostname) HF_TOKEN_len=${#HF_TOKEN} ==="
nvidia-smi -L
echo "--- 开工前配额 ---"
diskquota

cd /public/home/fjhui/ZW/scripts || exit 1
/public/home/fjhui/miniconda3/bin/python stream_extract.py --cancer "$CANCER"
rc=$?

echo "--- 收工后配额 ---"
diskquota
echo "=== job end rc=$rc ==="
exit $rc
