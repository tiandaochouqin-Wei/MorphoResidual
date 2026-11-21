#!/bin/bash
# lsf_c1_cov_luad.sh -- LSF wrapper for c1_covariates_cluster.py on the C1-LUAD confirmatory
# tree (112 analysed cases). POST HOC, governed by
# review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md, VERSION 3.1 (sections "PS n"). Cannot
# change the registered C1-LUAD reading.
#
# Cloned from server_export/scripts/lsf_c1_luad_plex_reread.sh: the same luad_c1 exports,
# the same unset list (including the two hash-skip switches, so no hash check can be
# downgraded) and guards, and an argument whitelist. Changes:
#   - MORPHO_OUT is ${C1}/results/posthoc/c1_covariates;
#   - C1COV_WRAPPER is exported (the production modes hash this wrapper against the sidecar);
#   - the documented flags are those of c1_covariates_cluster.py. Every data path is a
#     required flag with no default, exactly as for the published run and the plex re-read,
#     and each is hash-asserted inside the script (against a constant of
#     c1_run_test_luad.py, PS 1, or a --*-sha256 flag).
#   - no manifest mode: the one manifest is written through lsf_c1_cov_ucec.sh and covers
#     this tree too.
# Flags (each at most once, '--flag value' form):
#   --mode M            blind | gate | arm | placebo
#   --arm A             A_C3 | A_R20 | V_C3 | V_R20                  (mode arm)
#   --draws a-b         1 <= a <= b <= 119, e.g. 1-30               (mode placebo)
#   --sidecar P --sidecar-sha256 H --prespec P --closed-parts P      (modes arm, placebo)
#   --pinned-hosts P    PINNED_HOSTS.txt (sidecar-listed; names the host)  (modes arm, placebo)
#   --list-roots L      blind-run listing roots
#   data (required for every mode; the published run's paths):
#   --rule P --rule-sha256 H --case-ids P --protein P --protein-sha256 H
#   --rna-manifest-sha256 H --discovery-results P --discovery-sitepack P
#   --discovery-slide-map P --discovery-xwalk P --discovery-rna-manifest P
#   --batch-labels P --saved-matrix P
# Submit with bash, -n 4 and the host chosen by H-V1 (the published primary run: s002, -n 4),
# for example (P=/public/home/fjhui/ZW/scripts/pinned; full commands in SERVER_RUNBOOK.md):
#   bsub -q smp -n 4 -R "span[hosts=1]" -m s002 -o c1_cov_luad_gate_s002_%J.out \
#     bash lsf_c1_cov_luad.sh --mode gate --rule $P/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md ...
#
# morpho_env.sh is deliberately NOT sourced, for the reason lsf_c1_luad_test.sh gives. Every
# MORPHO_* is unset by literal name before the luad_c1 values are exported.
set -u

fatal() { echo "FATAL: $*" >&2; exit 1; }

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV \
      MORPHO_BATCH_AXIS MORPHO_MAX_LEVELS MORPHO_N_PERM \
      MORPHO_ALLOW_HPC_SLIDE_MAP_HASH MORPHO_ALLOW_HPC_XWALK_HASH C1COV_WRAPPER

C1=/public/home/fjhui/ZW/luad_c1
export MORPHO_ROOT="${C1}"
export MORPHO_RNA_MANIFEST="${C1}/manifest_rna_tumor_luad_c1.tsv"
export MORPHO_SLIDE_MAP=/public/home/fjhui/ZW/scripts/pinned/slide_type_map_luad_c1.tsv
export MORPHO_WSI_EMB_DIR="${C1}/WSI/emb"
export MORPHO_OUT="${C1}/results/posthoc/c1_covariates"
export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
# Read at import time by residual_analysis.py but not used (the protein comes from
# --protein); kept identical to lsf_c1_luad_test.sh / lsf_c1_luad_plex_reread.sh.
export MORPHO_PROTEIN_TSV="${C1}/protein_luad_c1.csv"
export MORPHO_ALIQUOT_XWALK="${C1}/pinned_aliquot_to_case_tumor_luad_c1.tsv"
C1COV_WRAPPER="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
export C1COV_WRAPPER

# --- environment sanity (same guards as lsf_c1_luad_test.sh) ------------------------
root_base=$(basename "${MORPHO_ROOT}")
[ "${root_base}" = "luad_c1" ] || fatal "basename(MORPHO_ROOT)='${root_base}', expected 'luad_c1'."
case "${MORPHO_OUT}" in
  "${MORPHO_ROOT}"/*) : ;;
  *) fatal "MORPHO_OUT='${MORPHO_OUT}' does not lie inside the luad_c1 tree (MORPHO_ROOT='${MORPHO_ROOT}')." ;;
esac
echo "environment OK: MORPHO_ROOT and MORPHO_OUT both confirmed under the luad_c1 tree."

# --- argument whitelist ----------------------------------------------------------------
VALUE_FLAGS=" --mode --arm --draws --sidecar --sidecar-sha256 --prespec --closed-parts \
--pinned-hosts --list-roots --rule --rule-sha256 --case-ids --protein \
--protein-sha256 --rna-manifest-sha256 --discovery-results --discovery-sitepack \
--discovery-slide-map --discovery-xwalk --discovery-rna-manifest --batch-labels \
--saved-matrix "
ARGS=()
SEEN=" "
MODE=""
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
      case "${val}" in
        *[!A-Za-z0-9_.,/:-]*) fatal "${tok} value '${val}' contains a character outside [A-Za-z0-9_.,/:-]." ;;
      esac
      if [ "${tok}" = "--mode" ]; then
        case "${val}" in
          blind|gate|arm|placebo) MODE="${val}" ;;
          *) fatal "--mode must be blind|gate|arm|placebo on the luad_c1 tree (the manifest goes through lsf_c1_cov_ucec.sh), got '${val}'." ;;
        esac
      fi
      if [ "${tok}" = "--arm" ]; then
        case "${val}" in
          A_C3|A_R20|V_C3|V_R20) : ;;
          *) fatal "--arm must be A_C3|A_R20|V_C3|V_R20, got '${val}'." ;;
        esac
      fi
      if [ "${tok}" = "--draws" ]; then
        printf '%s' "${val}" | grep -Eq '^[0-9]{1,3}-[0-9]{1,3}$' || \
          fatal "--draws must look like a-b (e.g. 1-30), got '${val}'."
      fi
      if [ "${tok}" = "--sidecar-sha256" ]; then
        printf '%s' "${val}" | grep -Eq '^[0-9a-fA-F]{64}$' || \
          fatal "--sidecar-sha256 must be 64 hex digits, got '${val}'."
      fi
      ARGS+=("${tok}" "${val}")
      SEEN="${SEEN}${tok} "
      shift 2
      continue
      ;;
  esac
  fatal "'${tok}' is not a documented flag of this wrapper (see the header). --workers and --cohort are set here."
done
[ -n "${MODE}" ] || fatal "--mode is required; see the header."

echo "=== job start: host=$(hostname) date=$(date) cores=$(nproc) workers=${MORPHO_N_WORKERS} args=${ARGS[*]} ==="
mkdir -p "${MORPHO_OUT}"
cd /public/home/fjhui/ZW/scripts || fatal "cannot cd to /public/home/fjhui/ZW/scripts"
PY=/public/home/fjhui/miniconda3/bin/python

# --- numpy 1.26.4, the version of the published October 2026 runs. The conda base now holds numpy
# 2.0.2 next to scipy 1.11.4, which cannot import together (ImportError: numpy.core.multiarray failed
# to import). The tree is kept as ONE archive because the account is at its file-count quota; it is
# unpacked into this job's node-local temporary folder, put first on PYTHONPATH and removed at exit.
NP_TAR=/public/home/fjhui/ZW/c1_covariates/numpy_1.26.4_tree.tar.gz
[ -f "${NP_TAR}" ] || fatal "private numpy archive missing: ${NP_TAR}"
NP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/c1cov_np_XXXXXX") || fatal "cannot create a node-local temporary folder"
trap 'rm -rf "${NP_DIR}"' EXIT
tar xzf "${NP_TAR}" -C "${NP_DIR}" || fatal "cannot unpack ${NP_TAR}"
export PYTHONPATH="${NP_DIR}${PYTHONPATH:+:${PYTHONPATH}}"
echo "  numpy archive : $(sha256sum "${NP_TAR}" | cut -d' ' -f1)"
NPV=$("${PY}" -c 'import numpy;print(numpy.__version__)' 2>&1 | tail -1)
[ "${NPV}" = "1.26.4" ] || fatal "numpy is '${NPV}', expected 1.26.4 (private copy unpacked from ${NP_TAR})"
echo "  numpy         : ${NPV}"

"${PY}" -u c1_covariates_cluster.py --cohort luad --workers "${MORPHO_N_WORKERS}" "${ARGS[@]}"
rc=$?
echo "=== job end: host=$(hostname) rc=${rc} date=$(date) ==="
exit ${rc}
