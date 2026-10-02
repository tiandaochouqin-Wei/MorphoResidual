#!/bin/bash
# Step 0b (2026-09-29) branch "protein": discovery-only PDC-vs-CDAP-TSV protein
# route equivalence check for LUAD. Uses PDC000153 (discovery) ONLY -- never
# PDC000489 (confirmatory). Metadata-only Step 0b authorisation, not the frozen
# rule; this script does not touch any confirmatory file.
#
# Before submitting, upload these two files (this run's local products) to the
# login node's scripts/ tree:
#   scp D:/claude/pajinsen/MorphoResidual_paper/server_export/scripts/protein_route_check.py \
#       <user>@<host>:/public/home/fjhui/ZW/scripts/
#   scp D:/claude/pajinsen/MorphoResidual_paper/server_export/scripts/pinned/luad_discovery_protein_via_pdc.csv \
#       <user>@<host>:/public/home/fjhui/ZW/scripts/pinned/
#
# Submit from the login node:
#   bsub -q smp -n 4 -R "span[hosts=1]" -o protein_route.out -e protein_route.err \
#        lsf_protein_route.sh
#
# Unlike lsf_c1_test.sh (which deliberately does NOT source morpho_env.sh, because
# the C1 confirmatory proteome is a pre-built case-indexed CSV, not a TMT TSV glob),
# THIS check's whole point is residual_analysis.load_protein_matrix()'s own
# discovery-CDAP-TSV route, which only morpho_env.sh resolves (it globs
# omics/protein/*.tmt*.tsv under the "luad" cohort root). So morpho_env.sh IS
# sourced here, with cohort "luad" (discovery), not any "*_c1" cohort string.
#
# Every MORPHO_* is unset first for the same reason lsf_c1_test.sh unsets them:
# bsub copies the submitting shell's environment, and a stale MORPHO_OUT left over
# from an earlier `source morpho_env.sh <other cohort>` in that same shell has
# silently redirected a whole run's output before (2026-09-15, UCEC output into
# pdac/results) -- this guard stays verbatim for the same reason lsf_c1_test.sh
# keeps it.
set -uo pipefail
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV

ZW=/public/home/fjhui/ZW
source "${ZW}/scripts/morpho_env.sh" luad || exit 1

echo "=== job start: host=$(hostname) date=$(date) cores=$(nproc) ==="
echo "resolved paths:"
morpho_print_env
if ! morpho_check_env; then
  echo "FATAL: inputs missing, refusing to run"
  exit 1
fi
echo

PDC_CSV="${ZW}/scripts/pinned/luad_discovery_protein_via_pdc.csv"
TESTED_CSV="${ZW}/scripts/pinned/luad/residual_results_tumoronly.csv"
if [ ! -e "${PDC_CSV}" ]; then
  echo "FATAL: ${PDC_CSV} missing -- scp it first (see header)"
  exit 1
fi
if [ ! -e "${TESTED_CSV}" ]; then
  echo "FATAL: ${TESTED_CSV} missing -- confirm the HPC-side pinned/luad/ copy "
  echo "  matches server_export/pinned/luad/residual_results_tumoronly.csv before relying on it"
  exit 1
fi

mkdir -p "${MORPHO_OUT}"
cd "${ZW}/scripts" || exit 1
/public/home/fjhui/miniconda3/bin/python -u protein_route_check.py \
    --pdc-csv "${PDC_CSV}" \
    --tested-genes-csv "${TESTED_CSV}" \
    --out-json "${MORPHO_OUT}/protein_route_check.json" \
    --out-csv "${MORPHO_OUT}/protein_route_check_pergene.csv"
rc=$?
echo
echo "=== job end rc=${rc} date=$(date) ==="
exit ${rc}
