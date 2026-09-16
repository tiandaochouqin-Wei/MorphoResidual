# 华农 /public 文件数配额告急 — 处置手册

## 问题
2026-09,平台管理员孟庆营提醒:账号 `fjhui` 的 GPFS 文件数配额已用 **1,964,635 / 2,000,000 (98.2%)**。
空间那条线(9.197T / 20T)很宽松,**卡住的是文件个数,不是容量**。
一旦运行中的作业把总数顶过 200 万,GPFS 会直接杀掉作业。

管理员另外说明:账号下有一点点数据是他们放的,他们之后自行清理;
要我们做的是——跑任务时随手清理自己产出的临时文件。

## 已排除
`extract_features.py` 一张切片只写一个 `.pt`。UCEC 393 张 = 393 个 inode。
抽取产物不是元凶。

## 最可疑的大户
1. **`WSI/raw/` 下的 DICOM 散件**。`download_wsi.py` 对每个 series 调 `zf.extractall()`,
   解压成一个目录装 N 个 `.dcm` + LICENSE + `.done`。五个癌种全量下载,轻松几十万 inode。
2. **`miniconda3/pkgs` 和 `~/.cache/pip`**。conda 解包后的 pkgs 目录是经典 inode 黑洞,
   十万级很常见,而且纯属可再生数据。
3. `__pycache__` 散件。

## 处置步骤

### 第一步 — 只看不动(先跑这个)
```bash
nohup bash inode_audit.sh > inode_audit.txt 2>&1 &
```
遍历百万文件在 GPFS 上要几分钟到十几分钟,所以丢后台。跑完 `cat inode_audit.txt`,
就知道那 196 万个文件具体压在哪个目录。

### 第二步 — 清可再生的垃圾(零风险)
```bash
bash inode_reclaim.sh caches          # 先空跑,看会清掉多少
bash inode_reclaim.sh caches --apply  # 确认后真清
bash inode_reclaim.sh pycache --apply
```
`caches` 刻意**没有**动 `~/.cache/huggingface/hub`,UNI 的权重在里面,删了要重下 1.3G。

### 第三步 — 打包已抽完癌种的 raw(数据不丢)
```bash
bash inode_reclaim.sh packraw ucec           # 空跑
bash inode_reclaim.sh packraw ucec --apply   # 真打包
```
把 `raw/` 打成单个 tar 再删散件,几万 inode 降到 1。原始数据留在 tar 里,
需要时 `tar -xf` 还原。脚本内置安全闸:`emb/` 里没有 `.pt` 就拒绝执行,
tar 校验不通过就不删原件。

**执行前请自行核对**该癌种抽取日志末尾是 `ok=<切片数> fail=0`。

### 第四步 — 以后别再堆(GBM / PDAC 用这个)
不要再「先全量下载,再全量抽取」。改用流式脚本:
```bash
export HF_TOKEN=hf_xxx
bsub -q interactive -gpu "num=1" -o stream_gbm.log -e stream_gbm.err \
     /public/home/fjhui/ZW/scripts/lsf_stream.sh gbm
```
`stream_extract.py` 逐张切片「下载 → 抽特征 → 立刻 rm -rf 原始文件」,
raw 峰值恒为一张片子。393 张片子净增 inode = 393。

特性:
- **断点续跑**:已有 `.pt` 的 series 自动跳过,作业被杀了直接重投。
- **不留半截文件**:先写 `.pt.partial` 再 rename。
- **启动先清残留**:上一轮被杀留下的临时目录会在开工时清掉。
- `--dry-run` 先看会处理哪些。

## 文件清单
| 文件 | 作用 |
|---|---|
| `inode_audit.sh` | 只读审计,定位 inode 大户 |
| `inode_reclaim.sh` | 回收,默认空跑,`--apply` 才动手 |
| `stream_extract.py` | 流式下载+抽取+删原件,inode 峰值 O(1) |
| `lsf_stream.sh` | 流式抽取的 LSF 提交包装 |
