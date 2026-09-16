#!/bin/bash
# Resolve every MORPHO_* path for one cohort. Source this; do not execute it.
#
#   source morpho_env.sh <cancer>
#
# Exists because the layout is not uniform: ccRCC predates the per-cancer
# directory scheme and lives at the ZW root, with embeddings under
# WSI/emb/ccrcc_real and manifests carrying no cohort suffix, while the other
# four follow ZW/<cancer>/. That branch had already been copied into one wrapper
# and got it wrong (the ccRCC job died on a missing proteome). One definition,
# sourced everywhere, so the next wrapper cannot reintroduce it.
#
# A caller may pre-set MORPHO_OUT or MORPHO_SIG_CSV; both are honoured.

MORPHO_CANCER=$1
if [ -z "${MORPHO_CANCER}" ]; then
  echo "FATAL: morpho_env.sh needs a cohort name" >&2
  return 1 2>/dev/null || exit 1
fi

ZW=/public/home/fjhui/ZW

if [ "${MORPHO_CANCER}" = "ccrcc" ]; then
  _croot=${ZW}
  _suffix=""
  _emb=${ZW}/WSI/emb/ccrcc_real
else
  _croot=${ZW}/${MORPHO_CANCER}
  _suffix="_${MORPHO_CANCER}"
  _emb=${ZW}/${MORPHO_CANCER}/WSI/emb
fi

# Match the TMT proteome table specifically. A bare *.tsv glob is unsafe: the
# ccRCC protein dir also holds a .summary.tsv (spectral counts, no "Log Ratio"
# columns) and some stray manifests; picking one of those yields an empty
# protein matrix with no error at all.
_pdir=${_croot}/omics/protein
_protein=$(ls ${_pdir}/*Proteome*.tmt*.tsv 2>/dev/null | grep -v '\.summary\.' | head -1)
if [ -z "${_protein}" ]; then
  _protein=$(ls ${_pdir}/*.tmt*.tsv 2>/dev/null | grep -v '\.summary\.' | head -1)
fi
if [ -z "${_protein}" ]; then
  echo "FATAL: no TMT proteome tsv under ${_pdir}/" >&2
  return 1 2>/dev/null || exit 1
fi

# Manifests come from scripts/pinned/, not scripts/. The scripts/ directory is
# shared with other people on this account, and the legacy builders
# (build_slide_map.py, build_aliquot_xwalk.py) write the same filenames with
# looser rules. On 2026-09-01 a concurrent run of them replaced the UCEC and PDAC
# slide maps and crosswalks mid-session -- re-admitting 15 and 79 adjacent-normal
# slides and 6 PDAC cell-line aliquots -- which silently contaminated every job
# that started afterwards. pinned/ holds this pipeline's verified inputs so that
# cannot happen again, and nobody else's files are touched.
_pin=${ZW}/scripts/pinned
export MORPHO_ROOT="${_croot}"
export MORPHO_RNA_MANIFEST="${_pin}/manifest_rna_tumor${_suffix}.tsv"
export MORPHO_ALIQUOT_XWALK="${_pin}/aliquot_to_case_tumor${_suffix}.tsv"
export MORPHO_SLIDE_MAP="${_pin}/slide_type_map${_suffix}.tsv"
export MORPHO_PROTEIN_TSV="${_protein}"
export MORPHO_WSI_EMB_DIR="${_emb}"
# MORPHO_OUT must be explicit: every analysis script writes a FIXED filename, so
# an unset value would silently overwrite the ccRCC results at ZW/results/.
export MORPHO_OUT="${MORPHO_OUT:-${_croot}/results}"
export MORPHO_N_WORKERS="${MORPHO_N_WORKERS:-8}"

morpho_print_env() {
  echo "  cohort        : ${MORPHO_CANCER}"
  echo "  root          : ${MORPHO_ROOT}"
  echo "  rna manifest  : ${MORPHO_RNA_MANIFEST}"
  echo "  aliquot xwalk : ${MORPHO_ALIQUOT_XWALK}"
  echo "  slide map     : ${MORPHO_SLIDE_MAP}"
  echo "  protein       : ${MORPHO_PROTEIN_TSV}"
  echo "  wsi emb       : ${MORPHO_WSI_EMB_DIR}"
  echo "  out           : ${MORPHO_OUT}"
  echo "  workers       : ${MORPHO_N_WORKERS}"
  [ -n "${MORPHO_SIG_CSV:-}" ] && echo "  sig csv       : ${MORPHO_SIG_CSV}"
  return 0
}

morpho_check_env() {
  local missing=0 f
  for f in "${MORPHO_RNA_MANIFEST}" "${MORPHO_ALIQUOT_XWALK}" "${MORPHO_SLIDE_MAP}" \
           "${MORPHO_PROTEIN_TSV}" "${MORPHO_WSI_EMB_DIR}" ${MORPHO_SIG_CSV:-}; do
    if [ ! -e "${f}" ]; then echo "  MISSING: ${f}"; missing=1; fi
  done
  return ${missing}
}
