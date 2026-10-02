#!/bin/bash
# lsf_posthoc_residual.sh <cohort> <arm>
#
# Runs residual_analysis.py, COMPLETELY UNCHANGED (no edit, no monkeypatch, no
# argument it accepts -- it takes none), against one variant-input ARM built by
# posthoc_build_arms.py. Post hoc relative to the discovery analysis already
# reported in the paper, governed by review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md
# (being written now). Does not touch, does not read for its own CV/permutation
# logic, and cannot change any C1-UCEC / C1-LUAD confirmatory reading -- those
# are separate frozen rules this file never sources.
#
# Inputs: the MANIFEST.tsv that `posthoc_build_arms.py --cohort <cohort>`
# already wrote at ${ZW}/scripts/pinned/posthoc/<cohort>/MANIFEST.tsv (columns:
# cohort, arm, field, file, sha256, n_patients_dropped, dropped_case_ids, note).
# This script looks up every row for the requested <cohort>/<arm>, RE-HASHES
# each named file right now, and FATALs if it does not match the manifest's
# recorded sha256 -- an arm directory is meant to be built once and read many
# times by many bsub submissions, never hand-edited in between.
#
# Outputs: whatever residual_analysis.py's own main() writes (today, a single
# fixed-name results/residual_results_tumoronly.csv -- see that file for the
# authoritative name) under
#   ${MORPHO_ROOT}/results/posthoc/<arm>/
# which is NEVER ${MORPHO_ROOT}/results (the published discovery results dir);
# this wrapper FATALs before running if MORPHO_OUT would resolve to that path.
#
# Dispatch is positional and exactly 2 arguments: <cohort> <arm>. Nothing else
# is accepted -- residual_analysis.py's main() takes no argparse arguments at
# all, so there is nothing legitimate to forward, and forwarding anything could
# only ever be an attempt to change a frozen input through the back door. This
# mirrors lsf_c1_luad_test.sh's own "every MORPHO_* unset by literal name"
# discipline, extended here to "no extra argv at all".
#
# scp upload (this file and posthoc_build_arms.py only -- residual_analysis.py,
# morpho_env.sh are already on the HPC and are NOT modified by this task):
#   scp server_export/scripts/lsf_posthoc_residual.sh server_export/scripts/posthoc_build_arms.py \
#     user@cluster:/public/home/fjhui/ZW/scripts/
#   ssh cluster chmod +x /public/home/fjhui/ZW/scripts/lsf_posthoc_residual.sh
#
# Build the arms first (once per cohort; see posthoc_build_arms.py's own header
# for the full command, including the two new pinned inputs it needs):
#   /public/home/fjhui/miniconda3/bin/python -u posthoc_build_arms.py --cohort ucec [...]
#
# Example bsub invocations, one job per arm (R always exists; F/D01-D19/Q exist
# only for cohort in {ucec, luad} -- see posthoc_build_arms.py):
#   bsub -q smp -n 8 -R "span[hosts=1]" -o posthoc_ucec_R.out \
#     lsf_posthoc_residual.sh ucec R
#   bsub -q smp -n 8 -R "span[hosts=1]" -o posthoc_ucec_F.out \
#     lsf_posthoc_residual.sh ucec F
#   bsub -q smp -n 8 -R "span[hosts=1]" -o posthoc_ucec_D01.out \
#     lsf_posthoc_residual.sh ucec D01
#   ... D02 .. D19 ...
#   bsub -q smp -n 8 -R "span[hosts=1]" -o posthoc_ucec_Q.out \
#     lsf_posthoc_residual.sh ucec Q
#   (same four/one-arm pattern for luad; only R for ccrcc/gbm/pdac)
#
set -u

fatal() { echo "FATAL: $*" >&2; exit 1; }

[ $# -eq 2 ] || fatal "usage: lsf_posthoc_residual.sh <cohort> <arm>  (exactly 2 positional args -- residual_analysis.py's main() takes no arguments, so none are forwarded; nothing else is accepted here)"
COHORT="$1"
ARM="$2"

case "${COHORT}" in
  ccrcc|luad|ucec|gbm|pdac) : ;;
  *) fatal "unknown cohort '${COHORT}', expected ccrcc|luad|ucec|gbm|pdac" ;;
esac

case "${ARM}" in
  R|F|Q) : ;;
  D0[1-9]|D1[0-9]) : ;;  # D01..D19
  *) fatal "unknown arm '${ARM}'; expected R, F, Q, or D01..D19" ;;
esac

ZW=/public/home/fjhui/ZW
SCRIPTS="${ZW}/scripts"
MANIFEST="${SCRIPTS}/pinned/posthoc/${COHORT}/MANIFEST.tsv"
[ -f "${MANIFEST}" ] || fatal "no manifest at ${MANIFEST} -- run posthoc_build_arms.py --cohort ${COHORT} [...] first"

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}

# Every MORPHO_* unset first, by LITERAL name (same reasoning, and the same
# 2026-09-15 incident, lsf_c1_luad_test.sh's header documents: a stale value
# left over from an earlier `source morpho_env.sh <other cohort>` in the same
# submitting shell must not silently redirect this run).
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV \
      MORPHO_BATCH_AXIS MORPHO_MAX_LEVELS MORPHO_N_PERM

source "${SCRIPTS}/morpho_env.sh" "${COHORT}" || fatal "sourcing morpho_env.sh ${COHORT} failed"

# Pull this arm's (field, file, sha256) rows out of the manifest by column
# NAME (not position), so a future manifest column reorder cannot silently
# misalign field<->file<->sha256.
ROWS=$(awk -F'\t' -v cohort="${COHORT}" -v arm="${ARM}" '
  NR==1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
  $col["cohort"] == cohort && $col["arm"] == arm { print $col["field"] "\t" $col["file"] "\t" $col["sha256"] }
' "${MANIFEST}")
[ -n "${ROWS}" ] || fatal "manifest ${MANIFEST} has no row for cohort=${COHORT} arm=${ARM} (arm not built for this cohort -- F/D*/Q only exist for ucec and luad)"

echo "=== job start: host=$(hostname) date=$(date) cohort=${COHORT} arm=${ARM} ==="

while IFS=$'\t' read -r FIELD FILE SHA; do
  [ -n "${FIELD}" ] || continue
  [ -f "${FILE}" ] || fatal "arm file missing on disk: ${FILE} (field ${FIELD}) -- was the posthoc/ tree moved or partially rebuilt?"
  GOT=$(sha256sum "${FILE}" | awk '{print $1}')
  [ "${GOT}" = "${SHA}" ] || fatal "sha256 mismatch for ${FIELD} arm file ${FILE}: got=${GOT} want=${SHA}. The arm file was modified after posthoc_build_arms.py wrote it -- rebuild the arm, do not hand-edit it."
  case "${FIELD}" in
    MORPHO_SLIDE_MAP)     export MORPHO_SLIDE_MAP="${FILE}" ;;
    MORPHO_ALIQUOT_XWALK) export MORPHO_ALIQUOT_XWALK="${FILE}" ;;
    *) fatal "manifest names field '${FIELD}' for arm ${ARM}, which this wrapper does not know how to override (only MORPHO_SLIDE_MAP / MORPHO_ALIQUOT_XWALK are recognised)" ;;
  esac
  echo "  override OK: ${FIELD} = ${FILE} (sha256 ${SHA:0:12}...)"
done <<< "${ROWS}"

# MORPHO_OUT: a private, arm-specific subdirectory under this cohort's own
# results/ tree -- never ${MORPHO_ROOT}/results itself, which is the published
# discovery output residual_analysis.py's fixed filename would otherwise
# silently overwrite.
export MORPHO_OUT="${MORPHO_ROOT}/results/posthoc/${ARM}"
case "${MORPHO_OUT}" in
  "${MORPHO_ROOT}/results/posthoc/"*) : ;;
  *) fatal "MORPHO_OUT='${MORPHO_OUT}' does not lie under '${MORPHO_ROOT}/results/posthoc/' -- refusing" ;;
esac
[ "${MORPHO_OUT}" != "${MORPHO_ROOT}/results" ] || fatal "MORPHO_OUT resolved to the published results dir -- refusing to run"

# Worker count follows what LSF actually granted (same reasoning as
# lsf_c1_luad_test.sh / lsf_residual.sh: a fixed core count does not reliably
# get scheduled on this cluster's smp queue).
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-${MORPHO_N_WORKERS:-8}}"

echo "  cores visible : $(nproc)"
morpho_print_env
if ! morpho_check_env; then
  fatal "inputs missing after arm override, refusing to run"
fi

mkdir -p "${MORPHO_OUT}"
cd "${SCRIPTS}" || fatal "cannot cd to ${SCRIPTS}"
/public/home/fjhui/miniconda3/bin/python -u residual_analysis.py
rc=$?

echo "=== job end: host=$(hostname) date=$(date) cohort=${COHORT} arm=${ARM} rc=${rc} ==="
# rc=1 only means the pilot >=50-protein gate was not cleared (see
# residual_analysis.py's own exit code); that is a result, not a job failure,
# per lsf_residual.sh's identical convention. Report it plainly.
exit ${rc}
