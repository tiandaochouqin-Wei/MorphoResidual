#!/bin/bash
# BRCA residual analysis, one arm per call. Same estimator as the published arm (cptac2_residual.py:25,87-102:
# 5-fold KFold(shuffle, rs=0) OOF closed-form ridge lambda=1, 20 PCs (full SVD), 1000 patient<->slide permutations
# RandomState(i), BH-FDR). ~2 s/gene x 10,031 genes ~ 5-6 h on CPU (cptac2_residual.py:15-17), so submit, do not run inline.
#   bsub -q smp -n 8 -m <HOST> -o dx_repro.log -e dx_repro.err \
#        /public/home/fjhui/ZW/scripts_ts_arm/lsf_cptac2_ts_residual.sh dx                 # DX arm re-run = reproduction check
#   bsub -q smp -n 8 -m <HOST> -o ts_resid.log -e ts_resid.err \
#        /public/home/fjhui/ZW/scripts_ts_arm/lsf_cptac2_ts_residual.sh ts                 # TS arm, slide-mean pooling (as DX)
#   ... ts --pool tile                                                                      # pre-declared sensitivity
# Patch v2: each call also writes the genes x (1+B=1001) increment matrix brca_incr_matrix[_ts][_tilepool].npz (~80 MB) next to the csv.
# Pin ONE host (bsub -m <HOST>) for all three jobs: the full-SVD PCA agrees across hosts only to ~3e-5 (lsf_residual.sh header).
# v3 (2026-10-08, adversarial review M5/M6/M7):
#   * MORPHO_WSI_EMB_DIR and MORPHO_OUT inherited from the submitting shell are cleared (bsub copies the environment), so the
#     TS arm can never read the DX embeddings; MORPHO_OUT is set explicitly per arm.
#   * dx: embeddings = ${ROOT}/emb_phikon unless DX_EMB_DIR is set; DX_EMB_DIR may be set ONLY by a dated amendment to the
#     pre-specification (RUNBOOK step 5b, the published DX job's embedding folder/density).
#   * ts: the embedding folder must hold exactly the frozen 189-file list (deployed as frozen_manifest_brca_ts_slides.txt,
#     checked by the sha256 of its sorted file names) intersected with the GDC manifest of step 1, minus the slides recorded in
#     ${ROOT}/ts_skipped_slides.txt (stem<TAB>reason; reasons no_tissue / extract_failed / download_failed). Any other
#     difference stops the job; cptac2_residual.py is then called with --expect-kept <that number>.
#   * (cross-review) The only extra argument accepted is '--pool tile' on the ts arm; the dx arm takes none. Anything else (for
#     example --perms, --max-genes, --stride) stops the job, so the frozen estimator cannot be changed from the command line.
#   * A job that stops before writing its summary csv (the last file cptac2_residual.py writes) is repeated once, unchanged, on the
#     same host with log suffix _rerun1 (pre-specification section 4, item 6); partial outputs are renamed *_failed1 first.
set -uo pipefail
ARM=${1:?usage: lsf_cptac2_ts_residual.sh <dx|ts> [extra cptac2_residual.py args]}
shift
case "${ARM}:$*" in
  dx:|ts:|"ts:--pool tile") ;;
  *) echo "FATAL: unsupported arguments '$*' (dx takes none; ts takes none or exactly '--pool tile')"; exit 1 ;;
esac
unset MORPHO_WSI_EMB_DIR MORPHO_OUT
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-8} MKL_NUM_THREADS=${MKL_NUM_THREADS:-8}
ROOT=/public/home/fjhui/ZW/cptac2_brca
SD=/public/home/fjhui/ZW/scripts_ts_arm
FROZEN_LIST_SHA256=9576e7337466cb20960c281db252b1fa36a80dece76f8199ca1bed988364e049   # sorted file names of the 189 frozen TS files
N_FROZEN=189
echo "=== resid start: arm=${ARM} host=$(hostname) date=$(date -u +%Y-%m-%dT%H:%M:%SZ) args=$* ==="
cd "${SD}" || exit 1   # patched copies live here; the published scripts in ../scripts stay untouched
grep -q -- "--slide-type" cptac2_residual.py && grep -q -- "--expect-kept" cptac2_residual.py || { echo "FATAL: cptac2_residual.py is not the patch-v3 copy (needs --slide-type and --expect-kept)"; exit 1; }
[ -s "${ROOT}/protein_brca.csv" ] && [ -s "${ROOT}/rna_case_file.tsv" ] || { echo "FATAL: ${ROOT}/protein_brca.csv or rna_case_file.tsv missing"; exit 1; }
EXTRA=()
case "${ARM}" in
  dx)
    export MORPHO_OUT=${ROOT}/results_dx_repro        # never overwrite the published results/
    if [ -n "${DX_EMB_DIR:-}" ]; then
      export MORPHO_WSI_EMB_DIR=${DX_EMB_DIR}         # only by dated amendment (see header)
    fi
    echo "dx embeddings: ${MORPHO_WSI_EMB_DIR:-${ROOT}/emb_phikon} ($(ls "${MORPHO_WSI_EMB_DIR:-${ROOT}/emb_phikon}"/*.pt 2>/dev/null | wc -l) .pt)"
    ;;
  ts)
    export MORPHO_OUT=${ROOT}/results_ts
    # pin to the published DX-arm case set (cut from rna_case_file.tsv, RUNBOOK step 0)
    [ -s "${ROOT}/dx_case_list.tsv" ] || { echo "FATAL: ${ROOT}/dx_case_list.tsv missing"; exit 1; }
    FZ=${SD}/frozen_manifest_brca_ts_slides.txt
    GM=${ROOT}/manifest_brca_ts_slides.txt
    SK=${ROOT}/ts_skipped_slides.txt
    for f in "${FZ}" "${GM}"; do [ -s "$f" ] || { echo "FATAL: $f missing"; exit 1; }; done
    [ -f "${SK}" ] || { echo "FATAL: ${SK} missing (create it at RUNBOOK step 4b; an empty file means no slide was skipped)"; exit 1; }
    T=$(mktemp -d) || exit 1
    tail -n +2 "${FZ}" | cut -f2 | tr -d '\r' | sed 's/\.svs$//' | LC_ALL=C sort -u > "${T}/frozen"
    got=$(tail -n +2 "${FZ}" | cut -f2 | tr -d '\r' | LC_ALL=C sort | sha256sum | cut -c1-64)
    [ "${got}" = "${FROZEN_LIST_SHA256}" ] && [ "$(wc -l < "${T}/frozen")" -eq "${N_FROZEN}" ] || { echo "FATAL: ${FZ} is not the frozen 189-file list (sha256 ${got})"; exit 1; }
    tail -n +2 "${GM}" | cut -f2 | tr -d '\r' | sed 's/\.svs$//' | LC_ALL=C sort -u > "${T}/gdc"
    LC_ALL=C comm -12 "${T}/frozen" "${T}/gdc" > "${T}/analysis"
    if ! cmp -s "${T}/frozen" "${T}/gdc"; then
      echo "NOTE: GDC manifest differs from the frozen list (reported, not substituted; analysis = intersection):"
      LC_ALL=C comm -3 "${T}/frozen" "${T}/gdc" | sed 's/^/  /'
    fi
    ls "${ROOT}/emb_phikon_ts/" 2>/dev/null | grep '\.pt$' | sed 's/\.pt$//' | LC_ALL=C sort -u > "${T}/pt"
    cut -f1 "${SK}" | tr -d '\r' | grep -v '^$' | LC_ALL=C sort -u > "${T}/skip"
    extra_pt=$(LC_ALL=C comm -23 "${T}/pt" "${T}/analysis")
    [ -z "${extra_pt}" ] || { echo "FATAL: embeddings outside the analysed frozen list (move them out of emb_phikon_ts; reported, not analysed):"; echo "${extra_pt}"; exit 1; }
    LC_ALL=C comm -23 "${T}/analysis" "${T}/skip" > "${T}/expected"
    if ! cmp -s "${T}/expected" "${T}/pt"; then
      echo "FATAL: TS embeddings != frozen list minus recorded skips (incomplete shard, or stale ts_skipped_slides.txt):"
      echo "  missing without a recorded skip:"; LC_ALL=C comm -23 "${T}/expected" "${T}/pt" | sed 's/^/    /'
      echo "  recorded as skipped but present:"; LC_ALL=C comm -12 "${T}/skip" "${T}/pt" | sed 's/^/    /'
      exit 1
    fi
    NEXP=$(wc -l < "${T}/expected")
    echo "ts slide check OK: frozen ${N_FROZEN}, analysed list $(wc -l < "${T}/analysis"), recorded skips $(wc -l < "${T}/skip"), embeddings ${NEXP}"
    rm -rf "${T}"
    EXTRA=(--case-list "${ROOT}/dx_case_list.tsv" --expect-kept "${NEXP}")
    ;;
  *)
    echo "FATAL: arm must be dx or ts"; exit 1 ;;
esac
/public/home/fjhui/miniconda3/bin/python -u cptac2_residual.py --project TCGA-BRCA --slide-type "${ARM}" ${EXTRA[@]+"${EXTRA[@]}"} "$@"
rc=$?
echo "=== resid end: arm=${ARM} rc=${rc} date=$(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
exit ${rc}
