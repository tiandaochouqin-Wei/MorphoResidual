#!/bin/bash
# C1-UCEC batch-stratified re-read (post hoc, D4 item 3, D6b --
# review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 0) under LSF.
#
# This wrapper is for c1_ucec_stratified_reread.py ONLY. Its sibling script,
# idc_ucec_derive_batch_labels.py, runs LOCALLY (it only probes the public IDC
# DICOM-header index over HTTP -- no HPC-resident RNA/protein/WSI data, no
# multiprocessing worker pool -- exactly like the LUAD Step 0b IDC scripts under
# review/C1_LUAD_STEP0B_2026-09-29/, none of which have an LSF wrapper either).
# Run it locally FIRST, then upload its two outputs here.
#
# scp upload list (everything else this script needs -- c1_run_test.py,
# split_replication_power.py, residual_analysis.py, the pinned discovery/protein
# data -- is already on the cluster from the original C1-UCEC run; only these four
# files need to go up):
#   scp server_export/scripts/pinned/c1_ucec_batch_labels.tsv \
#       review/C1_UCEC_D4_2026-09-29/c1_ucec_d4_summary.json \
#       user@cluster:/public/home/fjhui/ZW/scripts/pinned/
#   scp server_export/scripts/c1_ucec_stratified_reread.py \
#       server_export/scripts/lsf_c1_ucec_d4_strat.sh \
#       user@cluster:/public/home/fjhui/ZW/scripts/
#
# Submit from the login node:
#   Smoke (minutes; validates auc_obs ONLY against the published JSON and refuses
#   to compute any stratified AUC even if c1_ucec_d4_summary.json says
#   "qualifies": true -- B=5 permutations cannot support the frozen B=1000
#   within-operator null):
#     bsub -q smp -n 8 -R "span[hosts=1]" -o c1_ucec_d4_strat_smoke.out \
#       lsf_c1_ucec_d4_strat.sh --perms 5 --tag smoke
#   Full run (the real post hoc analysis, B=1000):
#     bsub -q smp -n 8 -R "span[hosts=1]" -o c1_ucec_d4_strat.out \
#       lsf_c1_ucec_d4_strat.sh --perms 1000
#
# morpho_env.sh is deliberately NOT sourced, for the same reason lsf_c1_test.sh and
# lsf_c1_ucec_strat.sh do not source it: the C1 proteome is PDC000439's
# quantDataMatrix pulled as a case-indexed CSV (frozen rule section 1), not the
# glob morpho_env.sh resolves. Every MORPHO_* is unset first -- bsub copies the
# submitting shell's environment, and a stale MORPHO_OUT from an earlier
# `source morpho_env.sh <other cohort>` in that shell silently redirected a whole
# UCEC run's output into pdac/results on 2026-09-15 (see lsf_c1_test.sh's own
# comment). Environment block below is byte-for-byte the same ucec_c1 exports as
# lsf_c1_test.sh / lsf_c1_ucec_strat.sh, on purpose: this script MUST read the
# exact same MORPHO_OUT as the original run, because that directory already holds
# c1_ucec_result.json -- c1_ucec_stratified_reread.py's --validate-against default
# finds it there without any extra path to keep in sync.
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV

C1=/public/home/fjhui/ZW/ucec_c1
export MORPHO_ROOT="${C1}"
export MORPHO_RNA_MANIFEST="${C1}/manifest_rna_tumor_ucec_c1.tsv"
export MORPHO_SLIDE_MAP=/public/home/fjhui/ZW/scripts/pinned/slide_type_map_ucec_c1.tsv
export MORPHO_WSI_EMB_DIR="${C1}/WSI/emb"
export MORPHO_OUT="${C1}/results"
# Worker count follows the slots LSF actually granted (see lsf_c1_test.sh's own
# comment on why a hardcoded count oversubscribes on this cluster's smp queue).
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
export MORPHO_PROTEIN_TSV="${C1}/protein_ucec_c1.csv"
export MORPHO_ALIQUOT_XWALK="${C1}/pinned_aliquot_to_case_tumor_ucec_c1.tsv"

echo "=== job start: host=$(hostname) date=$(date) cores=$(nproc) args=$* ==="
mkdir -p "${MORPHO_OUT}"
cd /public/home/fjhui/ZW/scripts || exit 1
/public/home/fjhui/miniconda3/bin/python -u c1_ucec_stratified_reread.py \
  --workers "${MORPHO_N_WORKERS}" \
  --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_ucec_batch_labels.tsv \
  --d4-summary /public/home/fjhui/ZW/scripts/pinned/c1_ucec_d4_summary.json \
  "$@"
rc=$?
echo "=== job end rc=${rc} date=$(date) ==="
exit ${rc}
