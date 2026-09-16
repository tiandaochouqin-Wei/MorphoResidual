#!/usr/bin/env bash
# Download the REAL CPTAC-CCRCC pathology WSIs via Aspera Faspex.
#
# Context: the earlier download_wsi.py pulled 727 series from a .tcia
# manifest that turned out to be the RADIOLOGY (CT/MR) sub-collection --
# verified live via getSeries (20/20 sampled series were CT/MR, 0 pathology).
# The real pathology data is distributed as a Faspex package, found via the
# "CPTAC-CCRCC Pathology Discovery Cohort" link on the TCIA collection page
# (https://www.cancerimagingarchive.net/collection/cptac-ccrcc/), not through
# the NBIA DICOM-series API at all.
#
# Downloads to a NEW directory (ccrcc_pathology_real) rather than overwriting
# the existing (wrong) raw/ccrcc, so nothing is lost if this needs retrying.
set -euo pipefail

ENV_NAME="aspera"
FASPEX_URL='https://faspex.cancerimagingarchive.net/aspera/faspex?context=eyJyZXNvdXJjZSI6InBhY2thZ2VzIiwidHlwZSI6ImV4dGVybmFsX2Rvd25sb2FkX3BhY2thZ2UiLCJpZCI6IjY4NiIsInBhc3Njb2RlIjoiNDMwZGIyZTM4NzFkMjFlNGNjYTEzMzJjNGY5NWZkMGZkNDc5NTNkYiIsInBhY2thZ2VfaWQiOiI2ODYiLCJlbWFpbCI6ImhlbHBAY2FuY2VyaW1hZ2luZ2FyY2hpdmUubmV0In0='
DEST_DIR="${MORPHO_WSI_DIR:-./WSI}/raw/ccrcc_pathology_real"

mkdir -p "${DEST_DIR}"
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME}"

echo "=== checking the exact receive command syntax this ascli version expects ==="
echo "(flag names have drifted across aspera-cli versions -- read this before"
echo " trusting the command below verbatim)"
ascli faspex5 packages receive --help || true
echo

echo "=== attempting download ==="
echo "destination: ${DEST_DIR}"
ascli faspex5 packages receive --url="${FASPEX_URL}" --to-folder="${DEST_DIR}"

echo
echo "=== done. quick sanity check ==="
find "${DEST_DIR}" -type f | wc -l
du -sh "${DEST_DIR}"
echo "spot-check one file's extension to confirm it's pathology (.svs/.dcm/.tif), not CT:"
find "${DEST_DIR}" -type f | head -5
