#!/bin/bash
# lsf_c1_ucec_plex_reread.sh -- LSF wrapper for c1_ucec_plex_reread.py, the POST HOC
# TMT-plex-stratified re-read of the published C1-UCEC confirmatory test (PDC000439).
#
# POST HOC. Governed by review/POSTHOC_UCEC_C1_PLEX_PRESPEC_2026-10-08.md (sidecar .sha256),
# released on GitHub before the --perms 1000 run, which is submitted with bsub -m <the host
# that ran the smoke>. Cannot change the published C1-UCEC reading or the D4 re-read.
# Upload also pinned/c1_case_ids_ucec.txt if it is not already on the cluster.
#
# Purpose : run c1_ucec_plex_reread.py in the same environment as the published run and the
#           D4 re-read (lsf_c1_test.sh / lsf_c1_ucec_d4_strat.sh: byte-identical ucec_c1
#           exports), with two changes borrowed from lsf_c1_luad_plex_reread.sh:
#             - MORPHO_OUT is the subtree ${C1}/results/posthoc, so a post hoc re-read is not
#               interleaved with the published outputs in ${C1}/results (the published
#               c1_ucec_result.json is read from ${C1}/results via --validate-against);
#             - it forwards only the three documented flags below, once each, in
#               '--flag value' form. Anything else is FATAL (--workers is set here from
#               the LSF slot count; every input path and hash is fixed inside the python
#               script).
# Flags   : --perms N   1000 = the analysis; 1..999 = smoke (reproduction gate only, no null,
#                       nothing written)
#           --tag T     output-name suffix, [A-Za-z0-9_-]
#           --blind     counts only (n, K, plex levels, dummies), nothing written
# Outputs : /public/home/fjhui/ZW/ucec_c1/results/posthoc/
#             c1_ucec_plex_reread_result[_<tag>].json, c1_ucec_plex_reread[_<tag>].log,
#             c1_ucec_plex_reread[_<tag>].npz   (never overwritten)
#
# Upload (everything else is already on the cluster from the published run and D4):
#   scp review/C1_UCEC_PLEX_2026-10-08/case_to_plex_ucec_c1.tsv \
#       user@cluster:/public/home/fjhui/ZW/scripts/pinned/
#   scp review/C1_UCEC_PLEX_2026-10-08/c1_ucec_plex_reread.py \
#       review/C1_UCEC_PLEX_2026-10-08/lsf_c1_ucec_plex_reread.sh \
#       user@cluster:/public/home/fjhui/ZW/scripts/
#   ssh cluster 'chmod +x /public/home/fjhui/ZW/scripts/lsf_c1_ucec_plex_reread.sh'
#   ssh cluster 'cd /public/home/fjhui/ZW/scripts && sha256sum c1_run_test.py \
#       residual_analysis.py split_replication_power.py pinned/case_to_plex_ucec_c1.tsv \
#       pinned/c1_case_ids_ucec.txt pinned/ucec/residual_results_tumoronly.csv \
#       ../ucec_c1/results/c1_ucec_result.json ../ucec_c1/results/c1_ucec_d4_stratified.npz'
#   (compare with the constants at the top of c1_ucec_plex_reread.py; the script asserts them)
#
# Submit from the login node, in this order (blind and smoke take minutes):
#   bsub -q smp -n 8 -R "span[hosts=1]" -o c1_ucec_plex_blind.out \
#     lsf_c1_ucec_plex_reread.sh --perms 5 --blind
#   bsub -q smp -n 8 -R "span[hosts=1]" -o c1_ucec_plex_smoke.out \
#     lsf_c1_ucec_plex_reread.sh --perms 5 --tag smoke
#   bsub -q smp -n 8 -R "span[hosts=1]" -o c1_ucec_plex_reread.out \
#     lsf_c1_ucec_plex_reread.sh --perms 1000
# Core count / host: the gate needs |auc_obs - 0.5974425564900907| <= 1e-9. PCA uses
# svd_solver="full", deterministic on one host but LAPACK paths can differ across hosts.
# D4 passed the same gate with -n 8; use -n 8 again. If the smoke misses the gate on auc_obs,
# resubmit the smoke on another host before touching anything else.
#
# morpho_env.sh is deliberately NOT sourced (see lsf_c1_ucec_d4_strat.sh). Every MORPHO_* is
# unset first: bsub copies the submitting shell's environment, and a stale MORPHO_OUT once
# redirected a whole UCEC run into pdac/results (2026-09-15).
set -u

fatal() { echo "FATAL: $*" >&2; exit 1; }

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV \
      MORPHO_BATCH_AXIS MORPHO_MAX_LEVELS MORPHO_N_PERM \
      MORPHO_ALLOW_HPC_SLIDE_MAP_HASH MORPHO_ALLOW_HPC_XWALK_HASH

C1=/public/home/fjhui/ZW/ucec_c1
export MORPHO_ROOT="${C1}"
export MORPHO_RNA_MANIFEST="${C1}/manifest_rna_tumor_ucec_c1.tsv"
export MORPHO_SLIDE_MAP=/public/home/fjhui/ZW/scripts/pinned/slide_type_map_ucec_c1.tsv
export MORPHO_WSI_EMB_DIR="${C1}/WSI/emb"
export MORPHO_OUT="${C1}/results/posthoc"
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
export MORPHO_PROTEIN_TSV="${C1}/protein_ucec_c1.csv"
export MORPHO_ALIQUOT_XWALK="${C1}/pinned_aliquot_to_case_tumor_ucec_c1.tsv"

# --- environment sanity ------------------------------------------------------------------
root_base=$(basename "${MORPHO_ROOT}")
[ "${root_base}" = "ucec_c1" ] || fatal "basename(MORPHO_ROOT)='${root_base}', expected 'ucec_c1'."
case "${MORPHO_OUT}" in
  "${MORPHO_ROOT}"|"${MORPHO_ROOT}"/*) : ;;
  *) fatal "MORPHO_OUT='${MORPHO_OUT}' does not lie inside the ucec_c1 tree (MORPHO_ROOT='${MORPHO_ROOT}')." ;;
esac
echo "environment OK: MORPHO_ROOT and MORPHO_OUT both confirmed under the ucec_c1 tree."

# --- argument whitelist --------------------------------------------------------------------
VALUE_FLAGS=" --perms --tag "
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
[ ${#ARGS[@]} -gt 0 ] || fatal "no arguments; see the header for the flags."

echo "=== job start: host=$(hostname) date=$(date) cores=$(nproc) workers=${MORPHO_N_WORKERS} args=${ARGS[*]} ==="
mkdir -p "${MORPHO_OUT}"
cd /public/home/fjhui/ZW/scripts || fatal "cannot cd to /public/home/fjhui/ZW/scripts"
PY=/public/home/fjhui/miniconda3/bin/python

"${PY}" -u c1_ucec_plex_reread.py --workers "${MORPHO_N_WORKERS}" "${ARGS[@]}"
rc=$?
echo "=== job end: host=$(hostname) rc=${rc} date=$(date) ==="
exit ${rc}
