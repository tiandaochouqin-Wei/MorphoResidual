#!/bin/bash
# lsf_c1_luad_plex_reread.sh -- LSF wrapper for c1_luad_plex_reread.py, the POST HOC
# TMT-plex-stratified re-read of the registered C1-LUAD test.
#
# POST HOC. Governed by review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md, which every run
# pins with --prespec/--prespec-sha256. Not part of the signed C1-LUAD rule
# (review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md); cannot change its section 5 reading
# (review/C1_LUAD_ADDENDA.md Entry 7).
#
# Purpose : run c1_luad_plex_reread.py in the same environment as the primary C1-LUAD
#           run (lsf_c1_luad_test.sh). The environment block and sanity guards below
#           copy that wrapper's, plus two changes:
#             - it also unsets MORPHO_ALLOW_HPC_SLIDE_MAP_HASH and
#               MORPHO_ALLOW_HPC_XWALK_HASH, so no hash check can be downgraded;
#             - it forwards only the documented flags listed below, once each and in
#               '--flag value' form. Anything else is FATAL: --workers (set here from
#               the LSF slot count), --flag=value, a repeated flag, or a positional.
#           The MORPHO_* paths are fixed in this file and cannot be changed by
#           arguments. Every input path is hash-asserted inside the script, either
#           against a constant or against a --*-sha256 argument.
# Inputs  : documented flags (all paths and hashes are required by the script; none
#           is defaulted here):
#             --perms N             1000 = the analysis; 1..999 = smoke (no null, no output)
#             --tag T               output-name suffix, [A-Za-z0-9_-]
#             --blind               counts only, nothing written
#             --rule P --rule-sha256 H          signed C1-LUAD rule (H must be a68dfbec...)
#             --prespec P --prespec-sha256 H    signed post hoc prespec
#             --case-ids P                      pinned/c1_case_ids_luad.txt
#             --protein P --protein-sha256 H    luad_c1/protein_luad_c1.csv
#             --rna-manifest-sha256 H           sha256 of ${C1}/manifest_rna_tumor_luad_c1.tsv
#             --discovery-results P --discovery-sitepack P --discovery-slide-map P
#             --discovery-xwalk P --discovery-rna-manifest P --batch-labels P
#                                               as for lsf_c1_luad_test.sh primary
#             --plex-map P --plex-map-sha256 H  pinned/case_to_plex_luad_c1.tsv
# Outputs : ${MORPHO_OUT} = /public/home/fjhui/ZW/luad_c1/results/posthoc/
#             c1_luad_plex_reread_result[_<tag>].json, c1_luad_plex_reread[_<tag>].log,
#             c1_luad_plex_reread[_<tag>].npz  (never overwritten; primary files untouched)
#
# Host pinning. The reproduction gate requires |auc_obs - 0.6439668501223573| <= 1e-9.
# PCA uses svd_solver="full", which is deterministic on one host, but LAPACK takes
# different code paths on different hosts (lsf_residual.sh: they agree to ~3e-5, not
# bit-for-bit). The primary run (job 75682618) ran on s002 with 4 cores, so submit with
# -m s002 -n 4. A run elsewhere that misses the gate stops FATAL; it never gives a silent
# difference. The smoke run (--perms 5, a few minutes) checks auc_obs on the chosen host
# before the ~5 h run is committed.
#
# Upload (the frozen scripts and inputs are already on the cluster from the primary run):
#   scp server_export/scripts/c1_luad_plex_reread.py \
#       server_export/scripts/lsf_c1_luad_plex_reread.sh \
#       user@cluster:/public/home/fjhui/ZW/scripts/
#   scp <module A output> case_to_plex_luad_c1.tsv user@cluster:/public/home/fjhui/ZW/scripts/pinned/
#   scp review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md (signed) to the same directory as the rule
#   ssh cluster chmod +x /public/home/fjhui/ZW/scripts/lsf_c1_luad_plex_reread.sh
#   ssh cluster sha256sum /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv \
#       /public/home/fjhui/ZW/luad_c1/manifest_rna_tumor_luad_c1.tsv   # record in prespec section F
#
# Submit (P = /public/home/fjhui/ZW/scripts/pinned). Blind first, then smoke, then the run:
#   bsub -q smp -n 4 -m s002 -R "span[hosts=1]" -o c1_luad_plex_blind.out \
#     lsf_c1_luad_plex_reread.sh --perms 5 --blind \
#       --rule <HPC path of the signed rule> \
#       --rule-sha256 a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022 \
#       --prespec <HPC path of the signed prespec> --prespec-sha256 <its sha256> \
#       --case-ids $P/c1_case_ids_luad.txt \
#       --protein /public/home/fjhui/ZW/luad_c1/protein_luad_c1.csv --protein-sha256 <...> \
#       --rna-manifest-sha256 <sha256 of luad_c1/manifest_rna_tumor_luad_c1.tsv> \
#       --discovery-results $P/luad/residual_results_tumoronly.csv \
#       --discovery-sitepack /public/home/fjhui/ZW/results/luad__sitepack_operator.csv \
#       --discovery-slide-map $P/slide_type_map_luad.tsv \
#       --discovery-xwalk $P/aliquot_to_case_tumor_luad.tsv \
#       --discovery-rna-manifest $P/manifest_rna_tumor_luad.tsv \
#       --batch-labels $P/c1_luad_batch_labels.tsv \
#       --plex-map $P/case_to_plex_luad_c1.tsv --plex-map-sha256 <module A's sha256>
#   bsub ... -o c1_luad_plex_smoke.out  lsf_c1_luad_plex_reread.sh --perms 5 [same flags, no --blind]
#   bsub ... -o c1_luad_plex_reread.out lsf_c1_luad_plex_reread.sh --perms 1000 [same flags]
# Use the same paths the primary run used (lsf_c1_luad_test.sh header). Every one is
# hash-asserted, so a wrong path stops FATAL.
#
# morpho_env.sh is deliberately NOT sourced, for the reason lsf_c1_luad_test.sh gives.
# Every MORPHO_* is unset by literal name before the luad_c1 values are exported.
set -u

fatal() { echo "FATAL: $*" >&2; exit 1; }

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV \
      MORPHO_BATCH_AXIS MORPHO_MAX_LEVELS MORPHO_N_PERM \
      MORPHO_ALLOW_HPC_SLIDE_MAP_HASH MORPHO_ALLOW_HPC_XWALK_HASH

C1=/public/home/fjhui/ZW/luad_c1
export MORPHO_ROOT="${C1}"
export MORPHO_RNA_MANIFEST="${C1}/manifest_rna_tumor_luad_c1.tsv"
export MORPHO_SLIDE_MAP=/public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad_c1.tsv
export MORPHO_WSI_EMB_DIR="${C1}/WSI/emb"
# A separate subtree, not ${C1}/results: that directory holds the registered test's
# own outputs, and a post hoc re-read must not be interleaved with them even though
# its file names differ and it writes with mode "x".
export MORPHO_OUT="${C1}/results/posthoc"
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
# Read at import time by residual_analysis.py but not used by this analysis (the
# protein comes from --protein); kept identical to lsf_c1_luad_test.sh.
export MORPHO_PROTEIN_TSV="${C1}/protein_luad_c1.csv"
export MORPHO_ALIQUOT_XWALK="${C1}/pinned_aliquot_to_case_tumor_luad_c1.tsv"

# --- environment sanity (same guards as lsf_c1_luad_test.sh) ------------------------
root_base=$(basename "${MORPHO_ROOT}")
[ "${root_base}" = "luad_c1" ] || fatal "basename(MORPHO_ROOT)='${root_base}', expected 'luad_c1'."
case "${MORPHO_OUT}" in
  "${MORPHO_ROOT}"|"${MORPHO_ROOT}"/*) : ;;
  *) fatal "MORPHO_OUT='${MORPHO_OUT}' does not lie inside the luad_c1 tree (MORPHO_ROOT='${MORPHO_ROOT}')." ;;
esac
echo "environment OK: MORPHO_ROOT and MORPHO_OUT both confirmed under the luad_c1 tree."

# --- argument whitelist ----------------------------------------------------------------
VALUE_FLAGS=" --perms --tag --rule --rule-sha256 --prespec --prespec-sha256 --case-ids \
--protein --protein-sha256 --rna-manifest-sha256 --discovery-results --discovery-sitepack \
--discovery-slide-map --discovery-xwalk --discovery-rna-manifest --batch-labels \
--plex-map --plex-map-sha256 "
BOOL_FLAGS=" --blind "
ARGS=()
SEEN=" "
while [ $# -gt 0 ]; do
  tok="$1"
  case "${tok}" in
    --*=*) fatal "use '--flag value', not '${tok}'." ;;
  esac
  case "${SEEN}" in
    *" ${tok} "*) fatal "flag ${tok} given more than once." ;;
  esac
  case "${VALUE_FLAGS}" in
    *" ${tok} "*)
      [ $# -ge 2 ] || fatal "${tok} needs a value."
      val="$2"
      case "${val}" in
        -*) fatal "${tok} value '${val}' looks like a flag." ;;
      esac
      if [ "${tok}" = "--perms" ]; then
        case "${val}" in
          ''|*[!0-9]*) fatal "--perms must be a positive integer, got '${val}'." ;;
        esac
      fi
      if [ "${tok}" = "--tag" ]; then
        case "${val}" in
          *[!A-Za-z0-9_-]*) fatal "--tag may only contain letters, digits, '_' and '-'." ;;
        esac
      fi
      ARGS+=("${tok}" "${val}")
      SEEN="${SEEN}${tok} "
      shift 2
      continue
      ;;
  esac
  case "${BOOL_FLAGS}" in
    *" ${tok} "*)
      ARGS+=("${tok}")
      SEEN="${SEEN}${tok} "
      shift
      continue
      ;;
  esac
  fatal "'${tok}' is not a documented flag of this wrapper (see the header). --workers is set here from the LSF slot count."
done
[ ${#ARGS[@]} -gt 0 ] || fatal "no arguments; see the header for the required flags."

echo "=== job start: host=$(hostname) date=$(date) cores=$(nproc) workers=${MORPHO_N_WORKERS} args=${ARGS[*]} ==="
mkdir -p "${MORPHO_OUT}"
cd /public/home/fjhui/ZW/scripts || fatal "cannot cd to /public/home/fjhui/ZW/scripts"
PY=/public/home/fjhui/miniconda3/bin/python

"${PY}" -u c1_luad_plex_reread.py --workers "${MORPHO_N_WORKERS}" "${ARGS[@]}"
rc=$?
echo "=== job end: host=$(hostname) rc=${rc} date=$(date) ==="
exit ${rc}
