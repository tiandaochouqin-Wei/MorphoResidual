# 服务器操作清单：C1 转录组协变量事后重读（Stage H 及其门控；对应预声明第 3.1 版）

依据：`review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md` **第 3.1 版**（下称 PS；"PS 9" 指其第 9 节）。
脚本：`c1_covariates_cluster.py`；包装脚本：`lsf_c1_cov_ucec.sh`、`lsf_c1_cov_luad.sh`、`lsf_c1_cov_discovery.sh`；关闭清单：`CLOSED_PARTS.txt`；固定主机清单：`PINNED_HOSTS.txt`（均在 `review/C1_COVARIATES_2026-11-21/`）。

**您需要书面同意的两次放行（缺一不可；同意书要记下 UTC 时间）**

1. **发布预声明**（第 3 部分第 3 步：tag `c1-covariates-posthoc-prespec-2026-11-21` 加 GitHub Release）。**这第一次放行的同意书必须点名每一个已关闭的读数**：`CLOSED_PARTS.txt` 里每一条对应的读数都要写出（例如"关闭 ucec 的 A_C3 与 V_C3"）；一条都没有就写"没有关闭任何读数"。每个关闭读数的原因（不超过十二个词）在门控记录里。
2. **发布 Stage H 记录**（第 5 部分第 1 步：tag `c1-covariates-posthoc-stageh-<日期>` 加 GitHub Release）。在此之前不跑本地 L1、L2。

（之后的结果记录是另一次放行，不在本清单范围内。）

**顺序（PS 9）**：第 1 部分上传 → 第 2 部分冻结前的门控 → **第 3 部分停**（我根据门控记录写 `PINNED_HOSTS.txt` 和 `CLOSED_PARTS.txt`、填门控记录、重新生成自测和验证日志、写 sidecar；您第一次书面同意后发布 tag + GitHub Release）→ 第 4 部分正式运行 → **第 5 部分再停**（您第二次书面同意后发布 Stage H 记录，然后才在本地跑 L1、L2）。

- 登录需要您本人操作（密码 + 动态口令）。下文 `USER@HOST` 是登录用户名和地址，端口 22333。
- 本地命令都在 Git Bash 里、从 `D:\claude\pajinsen\MorphoResidual_paper` 运行；服务器命令在登录后的 shell 里运行。
- 每个代码块只有一条命令，整块复制即可。尖括号 `<...>` 是要替换的占位符。
- 所有 `bsub` 都带 `bash`（10 月份不带 `bash` 时退出码是 127），日志名都带 `%J`（LSF 自动换成作业号），所以重交同一条命令时日志不会互相覆盖。
- 本地的 python 命令都带 `-B`（不在项目文件夹里留 `__pycache__`）。
- 冻结前的步骤只做计数、复现已发表的量（H-V1、H-V2）和导出发现队列的组成分数（H0），**不算任何新的增量**。
- **冻结前的每条命令都把日志写到能拷回 `review/C1_COVARIATES_2026-11-21/gate_logs/` 的地方**：本地命令用 `tee` 或 `>` 直接写进 `gate_logs/`；服务器上的作业用 `bsub -o` 写进 `/public/home/fjhui/ZW/c1_covariates/gate_logs/`，检查类命令用 `| tee` 写进同一个文件夹；2k 一次性取回。只做传输或建文件夹的命令（`scp`、`chmod`、`mkdir`）没有需要留存的输出。

---

## 0. 先核实的假设（任何一条不符：停下告诉我，不要改脚本）

1. **冻结模块**都在 `/public/home/fjhui/ZW/scripts/`，sha256 前缀依次为：`residual_analysis.py` a54a0a56、`split_replication_power.py` 79530c07、`c1_run_test.py` fc578c70、`c1_run_test_luad.py` 1d5ee05d、`composition_controls.py` 46b4a1d8。前四个由 10 月的 plex 重读核对过；**`composition_controls.py` 在服务器上的副本从未核对过**。
2. **本分析的脚本和三个包装脚本必须放在同一个文件夹** `/public/home/fjhui/ZW/scripts/`（与五个模块并列）：正式运行时脚本会在它自己旁边找三个包装脚本并逐个对照 sidecar。
3. **已发表的矩阵**：UCEC `/public/home/fjhui/ZW/ucec_c1/results/c1_ucec_d4_stratified.npz`（ef75f1c0…；s006 上的 plex 重读读过它）；LUAD `/public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz`（2034e711…；路径是从 `c1_run_test_luad.py` 的命名规则推出来的，**没有核对过**）。
4. **C1-UCEC 的输入**用 `c1_run_test.py` 和 plex 重读的默认路径：`ucec_c1/protein_ucec_c1.csv`、`scripts/pinned/ucec/residual_results_tumoronly.csv`（100d8f85…）、`scripts/pinned/slide_type_map_ucec.tsv`、`scripts/pinned/c1_case_ids_ucec.txt`（d777fc88…）。
5. **C1-LUAD 的输入**用 plex 重读（作业 75702342，s002）记录下来的路径和哈希：规则 `scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md`（a68dfbec…），蛋白 `luad_c1/protein_luad_c1.csv`（65c0daca…），RNA 清单 `luad_c1/manifest_rna_tumor_luad_c1.tsv`（b3f58f5d…），以及 `scripts/pinned/` 下的发现文件。
6. **发现队列**的环境用 `source morpho_env.sh ucec|luad`，与 10 月 7 日的测量误差导出相同。按 10 月的服务器记录，被采用的 `inputs_<cohort>.npz` 是 UCEC 在 **s006**、LUAD 在 **s002** 上导出的；PS 9 要求 H-V2 和 H0 在"导出这个 inputs 文件的主机"上跑，所以用这两台。如果您记得不是这两台，先告诉我。
7. **H-V1 的主机顺序**（PS 9）：UCEC，先是写出已保存矩阵的那个作业的主机，再是已发表检验的主机，最后是 s006；LUAD，先是主运行的主机，再是 s006。前两者要从旧作业日志里读出来（第 2b 步），我不猜。**一个读数只有在其列表里的每一台主机都跑过门控并且都没通过之后才算关闭；够不着的主机要等，不算失败。**
8. **队列** `smp` 的运行时限未知（第 2a 步查）。每个 1,001 列的臂预计 1–3 小时；119 个安慰剂合起来约等于 3 个臂的计算量，所以分成 4 个作业（1-30、31-60、61-90、91-119），每个约 1 个臂的量。
9. **Python** 为 `/public/home/fjhui/miniconda3/bin/python`，与已发表的运行相同；计算节点不需要联网；LSF 的 `-o` 支持 `%J`。
10. **新建的文件夹**：`/public/home/fjhui/ZW/c1_covariates/{inputs,frozen,discovery,gate_logs,stage_h_joblogs}`，以及两个确证树下的 `results/posthoc/c1_covariates`。脚本只在这些文件夹里写文件，从不覆盖已有文件。
11. **numpy 1.26.4 的私有压缩包**：`/public/home/fjhui/ZW/c1_covariates/numpy_1.26.4_tree.tar.gz`（sha256 `92f2aeb5603a18ed35b2242db5fa07f6ac43a96104e48a2b0ffa76f3663ef8ba`）。conda 里现在是 numpy 2.0.2，与 scipy 1.11.4 不能一起导入，所以三个包装脚本在每个作业开头把这个压缩包解压到计算节点本地的临时目录并检查版本为 1.26.4（不是就立刻退出）。**这个压缩包必须一直在**，正式运行前不要移动或删除它。你的账号文件数配额已满（2,000,000），所以我们只存这一个文件，而不是解开的目录。

---

## 1. 本地准备与上传（冻结前）

### 1a. 本地：从 inputs 文件导出两份患者名单（只读患者编号）
先建日志文件夹：
```bash
mkdir -p review/C1_COVARIATES_2026-11-21/gate_logs
```

```bash
python -B review/C1_COVARIATES_2026-11-21/c1_covariates_cluster.py --mode export-patients --cohort ucec --inputs-npz review/MEAS_ERROR_2026-10-08/export/ucec/inputs_ucec.npz --out review/C1_COVARIATES_2026-11-21/c1_cov_discovery_patients_ucec.txt 2>&1 | tee review/C1_COVARIATES_2026-11-21/gate_logs/export_patients_ucec.log
```

```bash
python -B review/C1_COVARIATES_2026-11-21/c1_covariates_cluster.py --mode export-patients --cohort luad --inputs-npz review/MEAS_ERROR_2026-10-08/export/luad/inputs_luad.npz --out review/C1_COVARIATES_2026-11-21/c1_cov_discovery_patients_luad.txt 2>&1 | tee review/C1_COVARIATES_2026-11-21/gate_logs/export_patients_luad.log
```

应分别显示 `wrote 100 patient ids` 和 `wrote 105 patient ids`。如果名单文件已经存在，脚本会拒绝覆盖（`REFUSED`），那就直接用已有的。

### 1b. 本地：记下要上传文件的哈希（上传后在服务器上对照）
```bash
sha256sum review/C1_COVARIATES_2026-11-21/c1_covariates_cluster.py review/C1_COVARIATES_2026-11-21/lsf_c1_cov_ucec.sh review/C1_COVARIATES_2026-11-21/lsf_c1_cov_luad.sh review/C1_COVARIATES_2026-11-21/lsf_c1_cov_discovery.sh review/C1_COVARIATES_2026-11-21/c1_cov_discovery_patients_ucec.txt review/C1_COVARIATES_2026-11-21/c1_cov_discovery_patients_luad.txt review/MEAS_ERROR_2026-10-08/export/ucec/inputs_ucec.npz review/MEAS_ERROR_2026-10-08/export/luad/inputs_luad.npz server_export/scripts/composition_controls.py 2>&1 | tee review/C1_COVARIATES_2026-11-21/gate_logs/upload_hashes_local.log
```

### 1c. 服务器上建文件夹
```bash
ssh -p 22333 USER@HOST "mkdir -p /public/home/fjhui/ZW/c1_covariates/inputs /public/home/fjhui/ZW/c1_covariates/frozen /public/home/fjhui/ZW/c1_covariates/gate_logs /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs"
```

### 1d. 上传脚本和三个包装脚本到 `scripts/`（WinSCP 请用二进制模式，不要转换换行）
```bash
scp -P 22333 review/C1_COVARIATES_2026-11-21/c1_covariates_cluster.py review/C1_COVARIATES_2026-11-21/lsf_c1_cov_ucec.sh review/C1_COVARIATES_2026-11-21/lsf_c1_cov_luad.sh review/C1_COVARIATES_2026-11-21/lsf_c1_cov_discovery.sh USER@HOST:/public/home/fjhui/ZW/scripts/
```

### 1e. 上传两份患者名单和两个 inputs 文件
```bash
scp -P 22333 review/C1_COVARIATES_2026-11-21/c1_cov_discovery_patients_ucec.txt review/C1_COVARIATES_2026-11-21/c1_cov_discovery_patients_luad.txt review/MEAS_ERROR_2026-10-08/export/ucec/inputs_ucec.npz review/MEAS_ERROR_2026-10-08/export/luad/inputs_luad.npz USER@HOST:/public/home/fjhui/ZW/c1_covariates/inputs/
```

### 1f. 服务器上（已登录）：核对上传的文件和冻结模块的哈希
```bash
cd /public/home/fjhui/ZW/scripts && sha256sum c1_covariates_cluster.py lsf_c1_cov_ucec.sh lsf_c1_cov_luad.sh lsf_c1_cov_discovery.sh residual_analysis.py split_replication_power.py c1_run_test.py c1_run_test_luad.py composition_controls.py morpho_env.sh ../c1_covariates/inputs/* ../ucec_c1/results/c1_ucec_d4_stratified.npz ../luad_c1/results/c1_luad_primary_primary.npz 2>&1 | tee /public/home/fjhui/ZW/c1_covariates/gate_logs/upload_check_server.log
```

- 脚本和三个包装脚本必须与 1b 的本地值一致。
- 五个模块的前缀必须依次是 a54a0a56、79530c07、fc578c70、1d5ee05d、46b4a1d8；`composition_controls.py` 还要与 1b 里本地 `server_export/scripts/composition_controls.py` 的完整哈希一致。
- 两份患者名单和两个 inputs 文件要与 1b 一致（inputs：e28b86fa…、b64cdf3b…）。
- 两个矩阵：ef75f1c0…、2034e711…。
- 如果 `composition_controls.py` 不存在：从本地 `server_export/scripts/composition_controls.py` 上传到 `scripts/`。如果存在但哈希不符：**停下告诉我，不要覆盖**。
- 如果某个矩阵不在这个路径：先 `ls` 找一下，把找到的路径告诉我，不要猜。
- 不要对这些文件跑 `sed -i 's/\r$//'`：它们已经是 LF 换行，冻结时记录的就是现在的哈希。

### 1g. 包装脚本加执行权限
```bash
chmod +x /public/home/fjhui/ZW/scripts/lsf_c1_cov_ucec.sh /public/home/fjhui/ZW/scripts/lsf_c1_cov_luad.sh /public/home/fjhui/ZW/scripts/lsf_c1_cov_discovery.sh
```

---

## 2. 冻结前的步骤（按顺序做）

**门控必须用最终版脚本跑。** 正式运行只接受"同一个脚本 sha256、同样的输入"写出的 H-V1 记录；冻结前脚本只要再改一个字节，2c–2j 就要重做。

### 2a. 查队列的运行时限
```bash
bqueues -l smp 2>&1 | tee /public/home/fjhui/ZW/c1_covariates/gate_logs/bqueues_smp.log
```

把 `RUNLIMIT`（默认和最大）告诉我（输出已存进 `gate_logs/bqueues_smp.log`，2k 取回）。如果默认时限短于 4 小时，正式运行时要在 `bsub` 后面加 `-W 时:分`（不超过最大时限）；具体数值我按截图告诉您。

### 2b. 先列出旧作业日志，再读出主机（H-V1 按这个顺序试）

先找出写了这三样东西的日志：已保存的 UCEC 矩阵（D4 作业）、已发表的 C1-UCEC 检验、C1-LUAD 主运行。
```bash
find /public/home/fjhui/ZW -maxdepth 3 -name '*.out' -newermt 2026-09-01 -print0 | xargs -0 -r grep -l -E 'c1_ucec_d4_stratified|C1-UCEC CONFIRMATORY TEST|STEP 1: H1' | xargs -r ls -la --time-style=long-iso 2>&1 | tee /public/home/fjhui/ZW/c1_covariates/gate_logs/old_job_logs_find.log
```

然后对上面找到的每个文件读出主机（把 `<日志文件>` 换成实际文件名，可以一次写多个）：
```bash
grep -H -m4 -E 'executed on host|job start: host=|Started at|c1_ucec_d4_stratified|C1-UCEC CONFIRMATORY TEST|STEP 1: H1' <日志文件> 2>&1 | tee /public/home/fjhui/ZW/c1_covariates/gate_logs/old_job_logs_hosts.log
```

- 两条命令的输出都已存进 `gate_logs/`，也发我；实际文件名和主机要写进门控记录。
- **UCEC 试主机的顺序（PS 9）**：(1) 写出 `c1_ucec_d4_stratified.npz` 的那个作业的主机；(2) 已发表的 C1-UCEC 检验（日志里有 `C1-UCEC CONFIRMATORY TEST`）的主机，如果和 (1) 不同；(3) s006，如果还没试过。
- **LUAD 试主机的顺序（PS 9）**：(1) C1-LUAD 主运行（日志里有 `STEP 1: H1`）的主机（应为 s002）；(2) s006，如果和 (1) 不同。
- 找不到日志：停下告诉我，不要猜主机。

### 2c. 空跑（只计数）：UCEC
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m s006 -o /public/home/fjhui/ZW/c1_covariates/gate_logs/c1_cov_ucec_blind_%J.out bash lsf_c1_cov_ucec.sh --mode blind
```

### 2d. 空跑（只计数）：LUAD
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m s002 -o /public/home/fjhui/ZW/c1_covariates/gate_logs/c1_cov_luad_blind_%J.out bash lsf_c1_cov_luad.sh --mode blind --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

空跑应显示：UCEC `cases 138, K 28, genes tested 10153, usable 10146 (identical to saved matrix: True), published set genes among them 2536`；LUAD `cases 112, K 22, genes tested 10417, usable 10406 (... True) ... 2287`；RNA 列数、全零列数、SD<1e-8 的列数；20 个成分的解释方差；**三个标记基因家族各自找到了多少个、每个标记的 SD、缺了哪些**；`C3 testable: True/False`；embedding 文件数；最后是 `BLIND DONE`。

- PS 2 的规则：某个家族找到的标记少于 6 个，或者没有一个标记的 SD 不为 0，这个队列的 C3 就"无法检验"，A_C3 和 V_C3 在冻结前关闭（写进 `CLOSED_PARTS.txt`，由我来写）。
- 结果文件夹清单写在 `c1_cov_<cohort>_blind_listing_*.tsv`，供 PS 0 第 7 条使用。

### 2e. H-V1：UCEC（按 2b 的顺序，一台一台试）
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <主机> -o /public/home/fjhui/ZW/c1_covariates/gate_logs/c1_cov_ucec_gate_<主机>_%J.out bash lsf_c1_cov_ucec.sh --mode gate
```

输出里有两个判定（都写进 `c1_cov_ucec_gate_hv1_<主机>_*.json`）：

- `strict`：基因列表完全相同，且第 0–5 列每个值的差 ≤ 1e-9。
- `fallback`：最大差 < 1e-6，且用重算的第 0 列算出的已发表基因集 AUC 与已发表的 AUC（UCEC 0.5974425564900907）相差 ≤ 1e-6。

怎么走（PS 9）：

- 末行是 `H-V1 PASSED (strict)`（退出码 0）：这台主机固定为 UCEC 的正式运行主机（我把它写进 `PINNED_HOSTS.txt`），列表里后面的主机不用再试。
- 末行是 `H-V1 STRICT FAILED ...`（退出码 3）：换列表里的下一台主机，重新提交同一条命令（`-m` 和 `-o` 里的主机名换掉）。
- 列表里的每台主机都没有 strict 通过：**后备规则**——在 fallback 通过的主机里选"最大差"最小的那台。这个选择由我根据各台的 json 来定。
- **只有列表里的每一台主机都跑过这个门控、并且都既没有 strict 也没有 fallback 通过，UCEC 的全部 Stage H 臂才在冻结前关闭**（我写进 `CLOSED_PARTS.txt`，您的第一次书面同意要点名）。
- **够不着的主机（节点下线、排不上队）要等，不算失败，也不能跳过**；等不到时告诉我，不要自己决定。
- 每台主机每次尝试的 `.out` 都已进 `gate_logs/`（2k 取回），门控记录要写出每一次尝试的主机和日志。每台大约几分钟到二十分钟。

### 2f. H-V1：LUAD（按 2b 的顺序：先主运行的主机，应为 s002；再 s006）
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <主机> -o /public/home/fjhui/ZW/c1_covariates/gate_logs/c1_cov_luad_gate_<主机>_%J.out bash lsf_c1_cov_luad.sh --mode gate --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

判定和走法同 2e（LUAD 已发表的 AUC 是 0.6439668501223573）：第一台主机 strict 通过就固定它，不必再试 s006；strict 不通过就试第二台；两台都没有 strict 通过时，按后备规则在 fallback 通过的主机里选最大差最小的那台。**只有两台主机（若主运行的主机就是 s006，则只有一台）都跑过并且都没通过，LUAD 的全部 Stage H 臂才在冻结前关闭**；够不着的主机要等。

### 2g. H-V2（gate-r20）：UCEC 发现队列，在 s006 上
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 1 -m s006 -o /public/home/fjhui/ZW/c1_covariates/gate_logs/c1_cov_disc_ucec_gate_r20_%J.out bash lsf_c1_cov_discovery.sh ucec --mode gate-r20 --patients /public/home/fjhui/ZW/c1_covariates/inputs/c1_cov_discovery_patients_ucec.txt --inputs-npz /public/home/fjhui/ZW/c1_covariates/inputs/inputs_ucec.npz
```

### 2h. H-V2（gate-r20）：LUAD 发现队列，在 s002 上
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 1 -m s002 -o /public/home/fjhui/ZW/c1_covariates/gate_logs/c1_cov_disc_luad_gate_r20_%J.out bash lsf_c1_cov_discovery.sh luad --mode gate-r20 --patients /public/home/fjhui/ZW/c1_covariates/inputs/c1_cov_discovery_patients_luad.txt --inputs-npz /public/home/fjhui/ZW/c1_covariates/inputs/inputs_luad.npz
```

两个都应以 `H-V2 PASSED` 结尾（按列对齐符号后，最大差 ≤ 1e-6）。`FAILED` 时不要换参数也不要换主机，告诉我：按 PS 9，这个门控的主机只有导出这个 inputs 文件的那一台，它跑过并且没有通过之后，这个器官确证队列的 A_R20、安慰剂和 V_R20 才在冻结前关闭（我写进 `CLOSED_PARTS.txt`）；够不着这台主机就等，不算失败。

### 2i. H0（export-scores）：UCEC
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 1 -m s006 -o /public/home/fjhui/ZW/c1_covariates/gate_logs/c1_cov_disc_ucec_export_scores_%J.out bash lsf_c1_cov_discovery.sh ucec --mode export-scores --patients /public/home/fjhui/ZW/c1_covariates/inputs/c1_cov_discovery_patients_ucec.txt --inputs-npz /public/home/fjhui/ZW/c1_covariates/inputs/inputs_ucec.npz
```

### 2j. H0（export-scores）：LUAD
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 1 -m s002 -o /public/home/fjhui/ZW/c1_covariates/gate_logs/c1_cov_disc_luad_export_scores_%J.out bash lsf_c1_cov_discovery.sh luad --mode export-scores --patients /public/home/fjhui/ZW/c1_covariates/inputs/c1_cov_discovery_patients_luad.txt --inputs-npz /public/home/fjhui/ZW/c1_covariates/inputs/inputs_luad.npz
```

两个都应以 `EXPORT DONE` 结尾，并打印每个家族的标记和 SD。出现 `sorted(common) ... differs`：发现树里的病例和 inputs 文件对不上，停下告诉我。出现 `REFUSED: C3 untestable in discovery`：不是故障，告诉我（这个器官的发现分析 C3 和 M|C3 会关闭）。

### 2k. 把门控记录和分数表传回本地

服务器上打包：①两个确证树的输出文件夹和发现树的输出文件夹；②`gate_logs/` 里的全部日志（作业 `.out` 和检查类命令的输出）。
```bash
cd /public/home/fjhui/ZW && tar czf c1_cov_gates.tgz ucec_c1/results/posthoc/c1_covariates luad_c1/results/posthoc/c1_covariates c1_covariates/discovery
```

```bash
tar czf /public/home/fjhui/ZW/c1_cov_gate_logs.tgz -C /public/home/fjhui/ZW/c1_covariates/gate_logs .
```

服务器上记下两个压缩包的哈希：
```bash
sha256sum /public/home/fjhui/ZW/c1_cov_gates.tgz /public/home/fjhui/ZW/c1_cov_gate_logs.tgz 2>&1 | tee /public/home/fjhui/ZW/c1_cov_gates_sha256.txt
```

本地取回（下面每块一条命令，按顺序）：
```bash
mkdir -p review/C1_COVARIATES_2026-11-21/gate_logs/server_logs server_export/results/c1_covariates/stage_h
```

```bash
scp -P 22333 USER@HOST:/public/home/fjhui/ZW/c1_cov_gates.tgz USER@HOST:/public/home/fjhui/ZW/c1_cov_gate_logs.tgz USER@HOST:/public/home/fjhui/ZW/c1_cov_gates_sha256.txt review/C1_COVARIATES_2026-11-21/gate_logs/
```

```bash
sha256sum review/C1_COVARIATES_2026-11-21/gate_logs/c1_cov_gates.tgz review/C1_COVARIATES_2026-11-21/gate_logs/c1_cov_gate_logs.tgz 2>&1 | tee review/C1_COVARIATES_2026-11-21/gate_logs/gates_tgz_sha256_local.log
```

本地的两个值必须与服务器上的（`gate_logs/c1_cov_gates_sha256.txt`）一致，再继续。
```bash
tar xzf review/C1_COVARIATES_2026-11-21/gate_logs/c1_cov_gates.tgz -C review/C1_COVARIATES_2026-11-21/gate_logs/
```

```bash
tar xzf review/C1_COVARIATES_2026-11-21/gate_logs/c1_cov_gate_logs.tgz -C review/C1_COVARIATES_2026-11-21/gate_logs/server_logs/
```

```bash
cp review/C1_COVARIATES_2026-11-21/gate_logs/c1_covariates/discovery/ucec/c1_cov_discovery_scores_ucec.tsv review/C1_COVARIATES_2026-11-21/gate_logs/c1_covariates/discovery/luad/c1_cov_discovery_scores_luad.tsv server_export/results/c1_covariates/stage_h/
```

### 2l. 本地：V3 和分数表格式检查（输出写进 `gate_logs/`）
```bash
python -B review/C1_COVARIATES_2026-11-21/c1_covariates_local.py --gate-v3 > review/C1_COVARIATES_2026-11-21/gate_logs/gate_v3.log 2>&1
```

```bash
python -B review/C1_COVARIATES_2026-11-21/c1_covariates_local.py --check-stage-h > review/C1_COVARIATES_2026-11-21/gate_logs/check_stage_h_pre.log 2>&1
```

两条命令把输出写进 `gate_v3.log` 和 `check_stage_h_pre.log`，屏幕上不显示；跑完后看末几行（`tail -n 5` 这两个文件）或发我。`gate_v3.log` 应以 `V3 VALIDATED` 结尾。`check_stage_h_pre.log` 这时只有两份分数表，其余 Stage H 文件显示缺失是正常的；分数表那两行必须是 OK。

---

## 3. 停：冻结与发布（服务器上什么都不跑）

1. 告诉我第 2 部分已经做完（日志已经在 `gate_logs/` 里；第 6 部分列的东西不用截图，告诉我文件名即可）。
2. 我按这个顺序做（全部在本地）：
   1. 根据 `gate_logs/` 里每台主机每次尝试的记录，填写预声明末尾的门控记录：self-tests、V1/V1b/V2、空跑、H-V1（试过的每台主机、每次的结果、固定了哪台）、H-V2 和 H0（各自的主机）、V3、关闭的部分及原因、仍可能出现的摘要变体。
   2. 写 `PINNED_HOSTS.txt`（每个队列一行 `pinned <队列> <主机>`，取自 H-V1 门控记录；某队列所有主机都没通过则没有这一行）和 `CLOSED_PARTS.txt`（什么都没关闭时只有注释行）。读数只有在其门控列表里的每一台主机都跑过并失败之后才能关闭。
   3. 用最终版脚本重新生成自测日志和验证日志，写进 `gate_logs/`；脚本从此不能再改一个字节。
   4. 写 sidecar `review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256`。
3. **您第一次书面同意**（记录 UTC 时间；同意书要点名每一个关闭的读数，没有就写"没有关闭任何读数"）后，提交并发布 tag `c1-covariates-posthoc-prespec-2026-11-21` 和 GitHub Release。
4. **发布之前**：不要跑 `--mode arm`、`--mode placebo`、`--mode manifest`，也不要做冒烟测试。脚本本身也有防护：这些模式需要冻结后的 sidecar、它的 sha256、预声明、`CLOSED_PARTS.txt` 和 `PINNED_HOSTS.txt`（二者都列在 sidecar 里，哈希必须一致），作业所在的主机必须等于 `PINNED_HOSTS.txt` 给这个队列的主机，并且本机有一份由同一脚本在同样输入上写出、`passed_strict` 或 `passed_fallback` 为真的 H-V1 记录，缺一项就拒绝运行。

---

## 4. 发布之后：正式运行（Stage H）

### 4a. 上传冻结的预声明、sidecar、`CLOSED_PARTS.txt` 和 `PINNED_HOSTS.txt`
```bash
scp -P 22333 review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 review/C1_COVARIATES_2026-11-21/CLOSED_PARTS.txt review/C1_COVARIATES_2026-11-21/PINNED_HOSTS.txt USER@HOST:/public/home/fjhui/ZW/c1_covariates/frozen/
```

```bash
sha256sum /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

- 第一行（sidecar 的 sha256）就是下面命令里的 `<SIDECAR_SHA>`，必须与发布记录里的值一致。另三行要与 sidecar 文件里列的值一致（脚本也会再查一遍）。
- `<UCEC主机>` 和 `<LUAD主机>` 是 `PINNED_HOSTS.txt` 里 `pinned ucec ...` 和 `pinned luad ...` 两行给出的主机，`-m` 用它；下面所有命令里 `--pinned-hosts` 指向同一个文件，脚本核对作业所在主机等于文件给的主机，不接受别的写法。
- 如果 `CLOSED_PARTS.txt` 关闭了某个臂（例如 `arm ucec A_C3`），对应的那条命令就不要提交（提交了脚本也会以 `REFUSED ... is closed` 拒绝，退出码 4）。`arm <cohort> A_R20cap` 表示这个队列的四个安慰剂作业都不提交（它和 `arm <cohort> A_R20` 总是一起出现）。
- 如果 `PINNED_HOSTS.txt` 没有某个队列的行，这个队列的全部臂都已在 `CLOSED_PARTS.txt` 里关闭，不提交它的任何作业。
- 如果 2a 显示需要 `-W`，在每条 `bsub` 后面加上我告诉您的 `-W 时:分`。
- 作业日志写进 `/public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/`（第 1c 步已建），4i 取回到 `server_export/results/c1_covariates/stage_h_joblogs/`。

### 4b. UCEC：四个 1,001 列的臂（每条命令一个作业，`-n 8`，每个约 1–3 小时）
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <UCEC主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_ucec_A_C3_%J.out bash lsf_c1_cov_ucec.sh --mode arm --arm A_C3 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <UCEC主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_ucec_A_R20_%J.out bash lsf_c1_cov_ucec.sh --mode arm --arm A_R20 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <UCEC主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_ucec_V_C3_%J.out bash lsf_c1_cov_ucec.sh --mode arm --arm V_C3 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <UCEC主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_ucec_V_R20_%J.out bash lsf_c1_cov_ucec.sh --mode arm --arm V_R20 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

### 4c. UCEC：119 个安慰剂基线，分 4 个作业（每个 draw 写一个文件；已存在的 draw 自动跳过）
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <UCEC主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_ucec_placebo_1-30_%J.out bash lsf_c1_cov_ucec.sh --mode placebo --draws 1-30 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <UCEC主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_ucec_placebo_31-60_%J.out bash lsf_c1_cov_ucec.sh --mode placebo --draws 31-60 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <UCEC主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_ucec_placebo_61-90_%J.out bash lsf_c1_cov_ucec.sh --mode placebo --draws 61-90 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 8 -R "span[hosts=1]" -m <UCEC主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_ucec_placebo_91-119_%J.out bash lsf_c1_cov_ucec.sh --mode placebo --draws 91-119 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

### 4d. LUAD：四个 1,001 列的臂（`-n 4`）
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <LUAD主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_luad_A_C3_%J.out bash lsf_c1_cov_luad.sh --mode arm --arm A_C3 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <LUAD主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_luad_A_R20_%J.out bash lsf_c1_cov_luad.sh --mode arm --arm A_R20 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <LUAD主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_luad_V_C3_%J.out bash lsf_c1_cov_luad.sh --mode arm --arm V_C3 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <LUAD主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_luad_V_R20_%J.out bash lsf_c1_cov_luad.sh --mode arm --arm V_R20 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

### 4e. LUAD：119 个安慰剂基线，分 4 个作业
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <LUAD主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_luad_placebo_1-30_%J.out bash lsf_c1_cov_luad.sh --mode placebo --draws 1-30 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <LUAD主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_luad_placebo_31-60_%J.out bash lsf_c1_cov_luad.sh --mode placebo --draws 31-60 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <LUAD主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_luad_placebo_61-90_%J.out bash lsf_c1_cov_luad.sh --mode placebo --draws 61-90 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 4 -R "span[hosts=1]" -m <LUAD主机> -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_luad_placebo_91-119_%J.out bash lsf_c1_cov_luad.sh --mode placebo --draws 91-119 --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt --rule /public/home/fjhui/ZW/scripts/pinned/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 65c0dacae5ffe22648cb2bf602bc0a2f02c029bd797a135d9dd284beb6145bf2 --rna-manifest-sha256 b3f58f5d2b5d319abb63a2177151694e46c59c253a685729223855b7ac856dfc --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv --saved-matrix /public/home/fjhui/ZW/luad_c1/results/c1_luad_primary_primary.npz
```

范围必须正好是这 4 个、互不重叠：脚本拒绝与已有范围重叠的新范围（例如已经交过 `1-30` 再交 `25-40`）。

### 4f. 查看进度和是否完成
```bash
bjobs -w
```

```bash
cd /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs && grep -H -E '^DONE:|REFUSED|STOPPED BY THE SCHEDULER|EXCEPTION|FATAL|Traceback|TERM_|job end' c1_cov_ucec_*.out c1_cov_luad_*.out
```

```bash
ls /public/home/fjhui/ZW/ucec_c1/results/posthoc/c1_covariates/c1_cov_ucec_A_R20cap_[0-9][0-9][0-9].npz /public/home/fjhui/ZW/luad_c1/results/posthoc/c1_covariates/c1_cov_luad_A_R20cap_[0-9][0-9][0-9].npz | wc -l
```

```bash
ls -la /public/home/fjhui/ZW/ucec_c1/results/posthoc/c1_covariates/c1_cov_*_A_C3.npz /public/home/fjhui/ZW/ucec_c1/results/posthoc/c1_covariates/c1_cov_*_A_R20.npz /public/home/fjhui/ZW/ucec_c1/results/posthoc/c1_covariates/c1_cov_*_V_*.npz /public/home/fjhui/ZW/ucec_c1/results/posthoc/c1_covariates/c1_cov_confirmatory_covariates_ucec.npz /public/home/fjhui/ZW/luad_c1/results/posthoc/c1_covariates/c1_cov_*_A_C3.npz /public/home/fjhui/ZW/luad_c1/results/posthoc/c1_covariates/c1_cov_*_A_R20.npz /public/home/fjhui/ZW/luad_c1/results/posthoc/c1_covariates/c1_cov_*_V_*.npz /public/home/fjhui/ZW/luad_c1/results/posthoc/c1_covariates/c1_cov_confirmatory_covariates_luad.npz
```

- 每个作业完成时，`.out` 里有 `DONE: <名字> completed (attempt k) ...`，输出文件夹里有 `<名字>.attempt<k>.json`、`.attempt<k>.log` 和 `.run<k>.json`。
- 什么都没关闭时，安慰剂计数应为 238（每个队列 119），每个队列有 4 个臂文件和 1 个 `c1_cov_confirmatory_covariates_<cohort>.npz`。

### 4g. 作业被停了怎么办（PS 6.1、PS 9 "Interrupted runs" 和 "Defects"）

- **调度器停的**（`.out` 里有 `STOPPED BY THE SCHEDULER`、`TERM_RUNLIMIT`、`TERM_OWNER`、节点故障、排队被删等）：等 `bjobs -w` 里这个作业已经消失，**把同一条命令原样再交一次**。只有运行时限可以调高（加或加大 `-W`），其余一字不改。次数**没有上限**，直到完成为止。
  - 脚本会自动编号为第 k 次尝试（写 `<名字>.attempt<k>.json`），把上一次留下的半成品改名为 `…partial.interrupted_before_attempt<k>` 保存下来，不删任何东西。
  - 安慰剂作业重交时，已写出的 draw 自动跳过，只补缺的。
  - 已经完成的输出永远不会重算或覆盖（再交会显示 `REFUSED: every output ... exists`，这表示已完成，不是故障）。
  - 如果显示 `another job of this run holds ...`：上一个同名作业还在跑，等它结束再交。
- **预检阶段的 FATAL**（哈希不符、主机不对、`PINNED_HOSTS.txt` 或 `CLOSED_PARTS.txt` 格式不对、sidecar 不对、找不到合格的 H-V1 记录）：这时还没开始计算，不算一次尝试。把 `.out` 发我，不要改任何文件。
- **`EXCEPTION` 或 Python `Traceback`**（PS 9 "Defects"）：**先报告，再决定能不能重交**。停下，不要自己重交，把 `.out` 和输出文件夹里的 `<名字>.attempt<k>.exception.json` 发我（`.out` 里会有 `REPORT TO THE DRAFTING SESSION BEFORE ANY RESUBMISSION`）。
  - **环境故障**（磁盘、内存、文件系统）：我确认后，原样重交同一条命令，并记录这次异常。
  - **其他任何情况**一律按脚本缺陷处理：需要一份带脚本 diff 的更正声明、新的 sidecar，您书面同意后发布，才能再跑；声明还要复现停掉的这次运行写到日志里的内容，并点名所有已完成的运行（它们的输出照旧有效）。
- **固定的主机从集群撤下了**：先告诉我。脚本不接受别的主机。按 PS 9，要在该队列列表里的下一台主机上重做 H-V1，这要走一份更正声明（新的 `PINNED_HOSTS.txt`、新的 sidecar，您书面同意后发布），其余一概不变。
- 在所有没有关闭的臂完成之前，稿件不投（PS 6.1）。

### 4h. 写 manifest（所有作业都完成后做一次；只通过 UCEC 包装脚本，覆盖两个队列的输出文件夹）
```bash
cd /public/home/fjhui/ZW/scripts && bsub -q smp -n 1 -o /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs/c1_cov_manifest_%J.out bash lsf_c1_cov_ucec.sh --mode manifest --sidecar /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 --sidecar-sha256 <SIDECAR_SHA> --prespec /public/home/fjhui/ZW/c1_covariates/frozen/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md --closed-parts /public/home/fjhui/ZW/c1_covariates/frozen/CLOSED_PARTS.txt --pinned-hosts /public/home/fjhui/ZW/c1_covariates/frozen/PINNED_HOSTS.txt
```

- 成功时末尾打印 `manifest sha256 = <64 位>`：截图发我（下文记作 `<MANIFEST_SHA>`，它进 Stage H 记录，L1 和 L2 都要用）。
- `REFUSED: ... expected output(s) missing`：还有没完成的臂或 draw，按 4g 重交，完成后再写 manifest。
- `REFUSED: ... is in progress`：还有作业在写文件，或者有被打断的作业还没重交：先按 4g 处理。
- `REFUSED: leftover ... manifest ... partial`（上一次写 manifest 时被打断）：在同一条命令末尾加 `--force-rewrite` 再交一次（只在 manifest 还不存在时允许）。
- manifest 一旦写出，就不能重写。

### 4i. 传回本地并核对
```bash
cd /public/home/fjhui/ZW && tar czf c1_cov_stage_h.tgz -C /public/home/fjhui/ZW/ucec_c1/results/posthoc/c1_covariates . -C /public/home/fjhui/ZW/luad_c1/results/posthoc/c1_covariates .
```

```bash
tar czf /public/home/fjhui/ZW/c1_cov_stage_h_joblogs.tgz -C /public/home/fjhui/ZW/c1_covariates/stage_h_joblogs .
```

```bash
sha256sum /public/home/fjhui/ZW/c1_cov_stage_h.tgz /public/home/fjhui/ZW/c1_cov_stage_h_joblogs.tgz /public/home/fjhui/ZW/ucec_c1/results/posthoc/c1_covariates/c1_cov_stage_h_manifest.sha256
```

本地：

```bash
scp -P 22333 USER@HOST:/public/home/fjhui/ZW/c1_cov_stage_h.tgz USER@HOST:/public/home/fjhui/ZW/c1_cov_stage_h_joblogs.tgz server_export/results/c1_covariates/
```

```bash
sha256sum server_export/results/c1_covariates/c1_cov_stage_h.tgz server_export/results/c1_covariates/c1_cov_stage_h_joblogs.tgz
```

两个压缩包的值必须与服务器上的一致。
```bash
tar xzf server_export/results/c1_covariates/c1_cov_stage_h.tgz -C server_export/results/c1_covariates/stage_h/
```

```bash
mkdir -p server_export/results/c1_covariates/stage_h_joblogs
```

```bash
tar xzf server_export/results/c1_covariates/c1_cov_stage_h_joblogs.tgz -C server_export/results/c1_covariates/stage_h_joblogs/
```

```bash
(cd server_export/results/c1_covariates/stage_h && sha256sum c1_cov_stage_h_manifest.sha256 && sha256sum -c --quiet c1_cov_stage_h_manifest.sha256 && echo ALL-OK)
```

- 第一行必须等于 `<MANIFEST_SHA>`；最后必须显示 `ALL-OK`（manifest 里每个文件的哈希都对得上）。
```bash
python -B review/C1_COVARIATES_2026-11-21/c1_covariates_local.py --check-stage-h
```

格式检查应全部 OK。之后告诉我。作业日志已在 `server_export/results/c1_covariates/stage_h_joblogs/`，Stage H 记录要用它们。

---

## 5. 再停：Stage H 记录，然后才是本地的 L1、L2

1. 我准备 Stage H 记录（manifest 和作业日志，只有文件哈希和日志）。**您第二次书面同意**（记录 UTC 时间）后，提交并发布 tag `c1-covariates-posthoc-stageh-<日期>` 和 GitHub Release。
2. **发布之后**才在本地跑（`<MANIFEST_SHA>` 用 Stage H 记录里的值）：
```bash
python -B review/C1_COVARIATES_2026-11-21/c1_covariates_local.py --stage L1 --manifest-sha <MANIFEST_SHA>
```

```bash
python -B review/C1_COVARIATES_2026-11-21/c1_covariates_local.py --stage L2 --manifest-sha <MANIFEST_SHA>
```

L1 完成后才跑 L2。两个都只跑一次；被打断时原样重跑（规则同 4g）。

---

## 6. 发我看的东西

冻结前的所有输出已经在服务器的 `gate_logs/` 和本地 `gate_logs/` 里；下面是要告诉我的结论和文件名，不需要截图。

1. 1f 的结果（`gate_logs/upload_check_server.log`）：哈希是否都对得上。
2. 2a 的 `RUNLIMIT`。
3. 2b 两条命令的结论（日志文件名和主机，列表里每台主机是哪台）。
4. 2c、2d 空跑输出中 `BLIND` 那一段（尤其是标记基因和 SD）。
5. 每台主机每次 H-V1：主机、largest |diff|、99.9 分位差、AUC 一行、strict 和 fallback 两行；有没有够不着的主机。
6. 2g–2j 的结尾几行，以及 2k 两个压缩包的哈希是否一致。
7. 2l 两个日志的末几行。
8. 正式运行后：4f 的 grep 输出和计数、4h 的 manifest sha256、4i 的 `ALL-OK`。
