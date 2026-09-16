#!/bin/bash
# 定位 GPFS 文件数(inode)配额的大户。只读，不删任何东西。
#
# 背景: 华农作重平台 /public 的配额是 空间 20T + 文件数 200 万 两条独立红线。
# 空间用了 9.2T 很宽松，文件数 196 万 /200 万 已经 98%，再涨就杀作业。
# 这个脚本回答"那 196 万个文件到底在哪"。
#
# 用法:  bash inode_audit.sh            # 审计 $HOME
#        bash inode_audit.sh /public/home/fjhui/ZW
#
# 注意: GPFS 上 find 遍历百万文件比较慢，预计几分钟到十几分钟。
#       建议丢后台: nohup bash inode_audit.sh > inode_audit.txt 2>&1 &

set -u
ROOT="${1:-$HOME}"
DEPTH="${2:-3}"

echo "=== inode audit: $ROOT (host=$(hostname), $(date)) ==="
echo

echo "--- 当前配额 ---"
diskquota 2>/dev/null || echo "(diskquota 不可用，跳过)"
echo

echo "--- 各子目录文件数 (深度 $DEPTH，只列 >1000 的，降序) ---"
# 一次遍历，按前 DEPTH 层路径聚合，避免对每个目录重复 find
find "$ROOT" -xdev -type f -printf '%p\n' 2>/dev/null \
  | awk -v root="$ROOT" -v d="$DEPTH" '
      BEGIN { rn = split(root, rp, "/") }
      {
        n = split($0, p, "/")
        key = ""
        lim = rn + d
        if (lim > n - 1) lim = n - 1
        for (i = 2; i <= lim; i++) key = key "/" p[i]
        if (key == "") key = "/"
        cnt[key]++
      }
      END { for (k in cnt) if (cnt[k] > 1000) printf "%10d  %s\n", cnt[k], k }
    ' | sort -rn
echo

echo "--- 常见 inode 黑洞快检 ---"
for d in \
  "$HOME/miniconda3/pkgs" \
  "$HOME/miniconda3/envs" \
  "$HOME/.cache/pip" \
  "$HOME/.cache/huggingface" \
  "$HOME/.cache/torch" \
  "$HOME/.conda" \
  "$HOME/.local/lib" \
  "$HOME/.singularity" \
  "$HOME/.nv"
do
  if [ -d "$d" ]; then
    printf "%10d  %s\n" "$(find "$d" -xdev -type f 2>/dev/null | wc -l)" "$d"
  fi
done
echo

echo "--- __pycache__ / .pyc 散件 ---"
printf "%10d  __pycache__ 目录下的文件\n" \
  "$(find "$ROOT" -xdev -type d -name __pycache__ -prune -exec find {} -type f \; 2>/dev/null | wc -l)"
echo

echo "--- 已抽完特征、raw 可回收的癌种 ---"
ZW="${ZW_ROOT:-/public/home/fjhui/ZW}"
if [ -d "$ZW" ]; then
  for c in "$ZW"/*/; do
    name=$(basename "$c")
    raw="$c/WSI/raw"; emb="$c/WSI/emb"
    [ -d "$raw" ] || continue
    nraw=$(find "$raw" -xdev -type f 2>/dev/null | wc -l)
    nemb=$(find "$emb" -maxdepth 1 -name '*.pt' 2>/dev/null | wc -l)
    printf "%-10s raw文件=%-9d emb切片=%-6d\n" "$name" "$nraw" "$nemb"
  done
else
  echo "(未找到 $ZW)"
fi
echo
echo "=== 审计结束 $(date) ==="
