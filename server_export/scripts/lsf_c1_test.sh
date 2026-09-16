#!/bin/bash
# C1-UCEC confirmatory test under LSF. Submit from the login node:
#   bsub -q smp -n 16 -R "span[hosts=1]" -o c1_smoke.out  lsf_c1_test.sh test --perms 5 --tag smoke
#   bsub -q smp -n 16 -R "span[hosts=1]" -o c1_test.out   lsf_c1_test.sh test --perms 1000
#   bsub -q smp -n 16 -R "span[hosts=1]" -o c1_boot.out   lsf_c1_test.sh bootstrap --boots 1000
#
# morpho_env.sh is deliberately NOT sourced: it resolves the proteome by globbing
# omics/protein/*.tmt*.tsv and exits FATAL when that misses, but the C1 proteome is
# PDC000439's quantDataMatrix pulled as a case-indexed CSV (frozen rule section 1).
# Every MORPHO_* is unset first -- bsub copies the submitting shell's environment, and
# a stale MORPHO_OUT from an earlier `source morpho_env.sh <other cohort>` in that shell
# silently redirected a whole UCEC run's output into pdac/results on 2026-09-15.
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV

C1=/public/home/fjhui/ZW/ucec_c1
export MORPHO_ROOT="${C1}"
export MORPHO_RNA_MANIFEST="${C1}/manifest_rna_tumor_ucec_c1.tsv"
export MORPHO_SLIDE_MAP=/public/home/fjhui/ZW/scripts/pinned/slide_type_map_ucec_c1.tsv
export MORPHO_WSI_EMB_DIR="${C1}/WSI/emb"
export MORPHO_OUT="${C1}/results"
# Worker count follows the slots LSF actually granted, not a hardcoded 16: a 16-core
# span[hosts=1] request does not schedule on this cluster's smp queue (affinity: "not
# enough processor units ... 3 hosts"), so this job is normally submitted with fewer
# cores, and 16 workers on 8 slots would just oversubscribe.
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
# Read at import time by residual_analysis.py but never used by c1_run_test.py
# (load_protein_matrix is not called; the C1 proteome is read from protein_ucec_c1.csv).
export MORPHO_PROTEIN_TSV="${C1}/protein_ucec_c1.csv"
export MORPHO_ALIQUOT_XWALK="${C1}/pinned_aliquot_to_case_tumor_ucec_c1.tsv"

echo "=== job start: host=$(hostname) date=$(date) cores=$(nproc) args=$* ==="
mkdir -p "${MORPHO_OUT}"
cd /public/home/fjhui/ZW/scripts || exit 1
/public/home/fjhui/miniconda3/bin/python -u c1_run_test.py "$@"
rc=$?
echo "=== job end rc=${rc} date=$(date) ==="
exit ${rc}
