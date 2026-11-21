#!/bin/bash
# lsf_c1_cov_discovery.sh -- LSF wrapper for c1_covariates_cluster.py on a DISCOVERY tree
# (UCEC n=100, LUAD n=105), for the two pre-freeze steps of PS 9 that read discovery RNA-seq:
#   gate-r20       H-V2: the R20 recipe on the inputs file's patients vs its rna_pcs20
#   export-scores  H0  : the three composition scores of the discovery patients (TSV)
# POST HOC, governed by review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md, VERSION 3.1.
# Both steps run before freezing and read RNA-seq only (no protein-linked statistic).
#
# Environment: as the measurement-error export that built inputs_<cohort>.npz
# (review/MEAS_ERROR_2026-10-08/lsf_measurement_error.sh): every MORPHO_* unset, MORPHO_OUT
# preset to a fresh folder, then `source morpho_env.sh <cohort>` (pinned manifests,
# crosswalks and slide maps), one BLAS thread. MORPHO_OUT is
# /public/home/fjhui/ZW/c1_covariates/discovery/<cohort>, outside every results/ folder.
#
# Usage: bash lsf_c1_cov_discovery.sh <ucec|luad> --mode <gate-r20|export-scores> \
#          --patients <file> --inputs-npz <file>
# Submit on the host whose export produced the inputs file (PS 9 H-V2; the 2026-10-07
# measurement-error exports adopted were UCEC s006, LUAD s002); the json records the host:
#   bsub -q smp -n 1 -m s006 -o c1_cov_disc_ucec_gate_r20_%J.out \
#     bash lsf_c1_cov_discovery.sh ucec --mode gate-r20 \
#       --patients /public/home/fjhui/ZW/c1_covariates/inputs/c1_cov_discovery_patients_ucec.txt \
#       --inputs-npz /public/home/fjhui/ZW/c1_covariates/inputs/inputs_ucec.npz
set -u

fatal() { echo "FATAL: $*" >&2; exit 1; }

CANCER="${1:-}"
case "${CANCER}" in
  ucec|luad) shift ;;
  *) fatal "first argument must be ucec or luad, got '${CANCER}'." ;;
esac
ZW=/public/home/fjhui/ZW

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV \
      MORPHO_BATCH_AXIS MORPHO_MAX_LEVELS MORPHO_N_PERM \
      MORPHO_ALLOW_HPC_SLIDE_MAP_HASH MORPHO_ALLOW_HPC_XWALK_HASH C1COV_WRAPPER

export MORPHO_OUT="${ZW}/c1_covariates/discovery/${CANCER}"
source "${ZW}/scripts/morpho_env.sh" "${CANCER}" || fatal "morpho_env.sh ${CANCER} failed"
export MORPHO_N_WORKERS=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
C1COV_WRAPPER="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
export C1COV_WRAPPER

root_base=$(basename "${MORPHO_ROOT}")
[ "${root_base}" = "${CANCER}" ] || fatal "basename(MORPHO_ROOT)='${root_base}', expected the discovery tree '${CANCER}'."
[ "${MORPHO_OUT}" = "${ZW}/c1_covariates/discovery/${CANCER}" ] || fatal "MORPHO_OUT was changed to '${MORPHO_OUT}'."

# --- argument whitelist ----------------------------------------------------------------
VALUE_FLAGS=" --mode --patients --inputs-npz "
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
        *[!A-Za-z0-9_.,/:-]*) fatal "${tok} value '${val}' contains a character outside [A-Za-z0-9_.,/:-]." ;;
      esac
      if [ "${tok}" = "--mode" ]; then
        case "${val}" in
          gate-r20|export-scores) MODE="${val}" ;;
          *) fatal "--mode must be gate-r20 or export-scores on a discovery tree, got '${val}'." ;;
        esac
      fi
      ARGS+=("${tok}" "${val}")
      SEEN="${SEEN}${tok} "
      shift 2
      continue
      ;;
  esac
  fatal "'${tok}' is not a documented flag of this wrapper (see the header)."
done
[ -n "${MODE}" ] || fatal "--mode is required; see the header."

echo "=== c1_cov discovery: cohort=${CANCER} host=$(hostname) date=$(date) args=${ARGS[*]} ==="
morpho_print_env
if ! morpho_check_env; then
  fatal "inputs missing, refusing to run"
fi
mkdir -p "${MORPHO_OUT}"
cd "${ZW}/scripts" || fatal "cannot cd to ${ZW}/scripts"
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
echo "  python        :$("${PY}" -c 'import sys,numpy,sklearn,pandas;print(sys.version.split()[0],numpy.__version__,sklearn.__version__,pandas.__version__)')"

"${PY}" -u c1_covariates_cluster.py --cohort "${CANCER}" "${ARGS[@]}"
rc=$?
echo "=== job end: host=$(hostname) rc=${rc} date=$(date) ==="
ls -la "${MORPHO_OUT}"
exit ${rc}
