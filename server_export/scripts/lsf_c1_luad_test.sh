#!/bin/bash
# C1-LUAD confirmatory test (c1_run_test_luad.py) / per-gene descriptive module
# (c1_luad_pergene.py) under LSF, per
# review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md, Appendix item 5.
#
# Unlike the looser C1-UCEC precedent lsf_c1_test.sh (whose c1_run_test.py hardcodes
# ucec_c1 defaults for --protein/--discovery-results/--discovery-slide-map), EVERY path
# argument of c1_run_test_luad.py and c1_luad_pergene.py is a REQUIRED argparse
# argument with NO default (Appendix item 1: "Paths. All paths are required argparse
# arguments, with no defaults."). This wrapper therefore does NOT hardcode
# --rule / --rule-sha256 / --case-ids / --protein / --discovery-results /
# --discovery-sitepack / --discovery-slide-map / --discovery-xwalk /
# --discovery-rna-manifest / --batch-labels (exact flag names read from
# c1_run_test_luad.py's argparse block, shared verbatim by c1_luad_pergene.py) -- it
# forwards "$@" untouched for all of them, and every bsub call below must supply them
# explicitly. Nothing here overrides a frozen statistical parameter (--perms, --boots,
# --variant, --h2-claimed, --d2-reading all come from "$@" too).
#
# Dispatch: the FIRST token of "$@" is the mode.
#   primary | sensitivity | bootstrap  -> c1_run_test_luad.py "$@" (mode token kept --
#                                          c1_run_test_luad.py's own argparse has a
#                                          required positional `mode` argument with
#                                          these same three choices, plus a
#                                          "descriptive" choice that only prints a
#                                          pointer to this file's `descriptive` branch).
#   descriptive                        -> c1_luad_pergene.py "$@" MINUS the leading
#                                          "descriptive" token (shifted off below):
#                                          c1_luad_pergene.py's argparse has NO
#                                          positional mode argument at all (confirmed
#                                          from its argparse block -- only --which
#                                          b1|b2|both selects which half runs), so
#                                          passing "descriptive" through to it would be
#                                          parsed as an unrecognized positional and
#                                          error out.
#
# scp upload (this file only -- c1_run_test_luad.py, c1_luad_pergene.py,
# residual_analysis.py and split_replication_power.py are uploaded separately and are
# NOT modified by this task):
#   scp server_export/scripts/lsf_c1_luad_test.sh user@cluster:/public/home/fjhui/ZW/scripts/
#   ssh cluster chmod +x /public/home/fjhui/ZW/scripts/lsf_c1_luad_test.sh
#
# Example bsub invocations. Paths below are the realistic pinned/results locations
# (HPC scripts/pinned/ mirrors this repo's server_export/scripts/pinned/); --rule and
# --rule-sha256 are placeholders -- fill in the actual uploaded path and the sha256
# recorded in section 7.3 once the rule is signed, before any of these run for real:
#
#   Blind smoke (section 5.6/D5g -- required for every confirmatory smoke test; prints
#   only n, K and gene counts, writes nothing):
#     bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_smoke.out \
#       lsf_c1_luad_test.sh primary --perms 5 --tag smoke --blind \
#         --rule /public/home/fjhui/ZW/scripts/review/C1_LUAD_FROZEN_RULE_SIGNED.md \
#         --rule-sha256 <sha256 recorded in section 7.3 at signature> \
#         --case-ids /public/home/fjhui/ZW/scripts/pinned/c1_case_ids_luad.txt \
#         --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv \
#         --discovery-results /public/home/fjhui/ZW/scripts/pinned/luad/residual_results_tumoronly.csv \
#         --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv \
#         --discovery-slide-map /public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad.tsv \
#         --discovery-xwalk /public/home/fjhui/ZW/scripts/pinned/aliquot_to_case_tumor_luad.tsv \
#         --discovery-rna-manifest /public/home/fjhui/ZW/scripts/pinned/manifest_rna_tumor_luad.tsv \
#         --batch-labels /public/home/fjhui/ZW/scripts/pinned/c1_luad_batch_labels.tsv
#
#   Primary run (B=1000, the frozen analysis -- section 5.6 prints a reading only at
#   exactly B=1000), same required paths as above:
#     bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_primary.out \
#       lsf_c1_luad_test.sh primary --perms 1000 --rule ... --rule-sha256 ... \
#         --case-ids ... --protein ... --discovery-results ... --discovery-sitepack ... \
#         --discovery-slide-map ... --discovery-xwalk ... --discovery-rna-manifest ... \
#         --batch-labels ...
#
#   Sensitivity ((s)/(s2), section 3), same required paths, plus --variant:
#     bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_sens_s.out \
#       lsf_c1_luad_test.sh sensitivity --variant s --perms 1000 [... same required paths]
#     bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_sens_s2.out \
#       lsf_c1_luad_test.sh sensitivity --variant s2 --perms 1000 [... same required paths]
#
#   Bootstrap (section 5.4), same required paths; add --h2-claimed ONLY after a
#   `primary` run gave p1<0.05 AND p_strat<0.05 (step 3 reached):
#     bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_boot.out \
#       lsf_c1_luad_test.sh bootstrap --boots 1000 [... same required paths]
#
#   Descriptive ((b1)/(b2), section 3 -- cannot enter the section 5 reading), same
#   required paths, no mode argument reaches c1_luad_pergene.py:
#     bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_pergene.out \
#       lsf_c1_luad_test.sh descriptive --perms 1000 [... same required paths]
#
# morpho_env.sh is deliberately NOT sourced (same reasoning as the UCEC precedent
# lsf_c1_test.sh): it resolves the proteome by globbing omics/protein/*.tmt*.tsv and
# exits FATAL when that misses, but the C1 proteome is a pulled CSV
# (c1_pull_omics_luad.py's protein_luad_c1.csv, forwarded via --protein above), not a
# TMT tsv morpho_env.sh would glob for.
#
# Every MORPHO_* is unset first, by LITERAL name -- bash `unset` does not accept a
# glob ('MORPHO_*' would not unset anything) -- and bsub copies the submitting shell's
# environment, so a stale MORPHO_OUT left over from an earlier
# `source morpho_env.sh <other cohort>` in that same shell would otherwise silently
# redirect this run's output (the exact 2026-09-15 incident lsf_c1_test.sh documents).
set -u

fatal() { echo "FATAL: $*" >&2; exit 1; }

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV \
      MORPHO_BATCH_AXIS MORPHO_MAX_LEVELS MORPHO_N_PERM

C1=/public/home/fjhui/ZW/luad_c1
export MORPHO_ROOT="${C1}"
export MORPHO_RNA_MANIFEST="${C1}/manifest_rna_tumor_luad_c1.tsv"
export MORPHO_SLIDE_MAP=/public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad_c1.tsv
export MORPHO_WSI_EMB_DIR="${C1}/WSI/emb"
export MORPHO_OUT="${C1}/results"
# Worker count follows the slots LSF actually granted, not a hardcoded 16 (same
# reasoning as lsf_c1_test.sh: a 16-core span[hosts=1] request does not schedule on
# this cluster's smp queue -- "not enough processor units ... 3 hosts" -- so this job
# is normally submitted with fewer cores).
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
# Read at import time by residual_analysis.py but not called directly by
# c1_run_test_luad.py / c1_luad_pergene.py (the C1 proteome comes from --protein,
# forwarded through "$@" above, per prepare()/load_confirmatory_protein()).
export MORPHO_PROTEIN_TSV="${C1}/protein_luad_c1.csv"
export MORPHO_ALIQUOT_XWALK="${C1}/pinned_aliquot_to_case_tumor_luad_c1.tsv"

# --- environment sanity: the exact bug this guards against ------------------------
# Same FATAL guard pattern lsf_precheck_luad.sh already uses for the D2 pre-check
# (basename(MORPHO_ROOT) check, MORPHO_OUT-inside-tree check), adapted to luad_c1: a
# stale MORPHO_OUT inherited from an earlier cohort's export in the same shell must not
# silently redirect this run's output into another cohort's results directory (the
# 2026-09-15 UCEC-into-pdac/results incident both precedents document).
root_base=$(basename "${MORPHO_ROOT}")
[ "${root_base}" = "luad_c1" ] || fatal "basename(MORPHO_ROOT)='${root_base}', expected 'luad_c1'."
case "${MORPHO_OUT}" in
  "${MORPHO_ROOT}"|"${MORPHO_ROOT}"/*) : ;;
  *) fatal "MORPHO_OUT='${MORPHO_OUT}' does not lie inside the luad_c1 tree (MORPHO_ROOT='${MORPHO_ROOT}')." ;;
esac
echo "environment OK: MORPHO_ROOT and MORPHO_OUT both confirmed under the luad_c1 tree."

MODE="${1:-}"
[ -n "${MODE}" ] || fatal "usage: lsf_c1_luad_test.sh primary|sensitivity|bootstrap|descriptive [args...] (see the header comment for the required --rule/--case-ids/--protein/... paths, none of which this wrapper defaults)"

echo "=== job start: host=$(hostname) date=$(date) cores=$(nproc) mode=${MODE} args=$* ==="
mkdir -p "${MORPHO_OUT}"
cd /public/home/fjhui/ZW/scripts || fatal "cannot cd to /public/home/fjhui/ZW/scripts"
PY=/public/home/fjhui/miniconda3/bin/python

case "${MODE}" in
  primary|sensitivity|bootstrap)
    "${PY}" -u c1_run_test_luad.py "$@"
    ;;
  descriptive)
    # c1_luad_pergene.py's argparse has no positional mode argument (--which
    # b1|b2|both is the only mode-like flag, and it defaults to "both") -- shift the
    # "descriptive" token off before forwarding, or argparse would reject it as an
    # unrecognized positional.
    shift
    "${PY}" -u c1_luad_pergene.py "$@"
    ;;
  *)
    fatal "unknown mode '${MODE}', expected primary|sensitivity|bootstrap|descriptive"
    ;;
esac
rc=$?
echo "=== job end rc=${rc} date=$(date) ==="
exit ${rc}
