#!/usr/bin/env bash
# Download CPTAC-CCRCC pathology WSIs for the pilot cohort.
#
# Verified live on 2026-07-13:
#   - CPTAC-CCRCC pathology slides are NOT queryable through the standard
#     NBIA getSeries REST API (Collection=CPTAC-CCRCC only exposes
#     CT/MR/RTSTRUCT there). They are distributed as a static .tcia manifest
#     (a plain-text list of SeriesInstanceUIDs) hosted directly on
#     cancerimagingarchive.net.
#   - The manifest's actual byte-download endpoint (DownloadServlet) requires
#     the official client's internal auth handshake -- a hand-built curl GET
#     against it 500s with "password is null" in DownloadServlet.decrypt().
#     There is no working curl-only path; the official NBIA Data Retriever
#     CLI is required.
#   - The retriever ships only as RPM/DEB, but can be extracted without root
#     (no `yum install` needed) and run directly from the extracted tree.
set -euo pipefail

WSI_ROOT="${MORPHO_WSI_DIR:-./WSI}"
RAW_DIR="${WSI_ROOT}/raw/ccrcc"
TOOL_DIR="${WSI_ROOT}/tools/nbia-retriever"
MANIFEST_URL="https://www.cancerimagingarchive.net/wp-content/uploads/TCIA-CPTAC-CCRCC_v11_20230818.tcia"
MANIFEST_PATH="${WSI_ROOT}/manifests/ccrcc_pathology.tcia"
RPM_URL="https://github.com/CBIIT/NBIA-TCIA/releases/download/DR-4_4_3-TCIA-20240916-1/nbia-data-retriever-4.4.3-1.x86_64.rpm"
RPM_PATH="${WSI_ROOT}/tools/nbia-data-retriever.rpm"

mkdir -p "${RAW_DIR}" "${TOOL_DIR}" "$(dirname "${MANIFEST_PATH}")"

echo "=== 1. IMPORTANT: check for a newer manifest first ==="
echo "This script uses TCIA-CPTAC-CCRCC_v11_20230818.tcia (verified reachable"
echo "as of 2026-07-13). TCIA periodically posts newer versions. Before a"
echo "large run, check https://www.cancerimagingarchive.net/collection/cptac-ccrcc/"
echo "for a higher-numbered .tcia manifest and update MANIFEST_URL above if found."
echo

echo "=== 2. downloading manifest ==="
if [ -s "${MANIFEST_PATH}" ]; then
    echo "manifest already present at ${MANIFEST_PATH}, skipping download"
    echo "(delete it first if you want to force a fresh download)"
else
    curl -fSL --retry 3 --connect-timeout 15 -o "${MANIFEST_PATH}" "${MANIFEST_URL}"
fi
n_series=$(grep -c '^1\.' "${MANIFEST_PATH}" || true)
echo "manifest ready: ${MANIFEST_PATH} (${n_series} series/slides listed)"

echo "=== 3. setting up NBIA Data Retriever CLI (no root) ==="
# Confirmed live 2026-07-13: this RPM's internal layout uses the lowercase
# Debian-style path, not the CamelCase path the TCIA wiki documents for RPM.
RETRIEVER_BIN="${TOOL_DIR}/opt/nbia-data-retriever/bin/nbia-data-retriever"
if [ ! -x "${RETRIEVER_BIN}" ] && [ -d "${TOOL_DIR}/opt" ]; then
    # already extracted under a different layout -- auto-discover instead of failing
    found=$(find "${TOOL_DIR}" -iname "*retriever" -type f -perm -u+x 2>/dev/null | head -1)
    [ -n "${found}" ] && RETRIEVER_BIN="${found}"
fi
if [ ! -x "${RETRIEVER_BIN}" ]; then
    if [ -s "${RPM_PATH}" ]; then
        echo "RPM already present at ${RPM_PATH}, skipping download"
    else
        curl -fSL --retry 3 --connect-timeout 15 -o "${RPM_PATH}" "${RPM_URL}"
    fi
    if command -v rpm2cpio >/dev/null 2>&1; then
        ( cd "${TOOL_DIR}" && rpm2cpio "$(realpath "${RPM_PATH}")" | cpio -idmv )
    else
        echo "rpm2cpio not found on this system. Options:"
        echo "  a) 'module load rpm2cpio' or similar if your HPC provides it as a module"
        echo "  b) ask the platform admin to install it (it's a tiny, standard RPM utility)"
        echo "  c) extract the RPM on a machine you control and rsync the"
        echo "     ${TOOL_DIR}/opt/NBIADataRetriever tree over"
        exit 1
    fi
fi
if [ ! -x "${RETRIEVER_BIN}" ]; then
    echo "ERROR: expected binary not found at ${RETRIEVER_BIN} after extraction."
    echo "Run 'find ${TOOL_DIR} -iname \"*retriever*\"' to locate the real path"
    echo "and adjust RETRIEVER_BIN above -- RPM internal layout can vary by version."
    exit 1
fi
echo "retriever ready: ${RETRIEVER_BIN}"

echo "=== 4. downloading WSIs (this is the multi-hour, multi-GB step) ==="
echo "destination: ${RAW_DIR}"
"${RETRIEVER_BIN}" --cli "${MANIFEST_PATH}" -d "${RAW_DIR}" --agree-to-license -v

echo "=== done. verify counts ==="
find "${RAW_DIR}" -type f | wc -l
du -sh "${RAW_DIR}"
