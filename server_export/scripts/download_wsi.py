#!/usr/bin/env python3
"""
Download CPTAC-CCRCC pathology WSIs directly via the NBIA v4 REST API,
bypassing the official NBIA Data Retriever CLI entirely.

Why this replaces download_wsi.sh: the legacy DownloadServlet used by the
retriever CLI requires the client's internal encrypted auth handshake (a
bare curl GET 500s with "password is null" in DownloadServlet.decrypt()),
and the retriever CLI itself was verified to fail silently (exit 1, zero
output, no missing libs) on this HPC on 2026-07-13 -- likely a jpackage
runtime issue specific to this environment. The v4 API's getImage endpoint
was confirmed live the same day to work with a plain unauthenticated GET,
returning a ZIP of the series' DICOM instances (LICENSE + N .dcm files).
No Java, no RPM extraction, no root needed.

Manifest format (verified): key=value header lines, then
"ListOfSeriesToDownload=" followed by one SeriesInstanceUID per line.
"""
import io
import os
import sys
import time
import urllib.request
import urllib.error
import zipfile

API_URL = "https://services.cancerimagingarchive.net/nbia-api/services/v4/getImage"
MANIFEST_PATH = os.environ.get(
    "MORPHO_WSI_MANIFEST",
    os.path.join(os.environ.get("MORPHO_WSI_DIR", "./WSI"), "manifests", "ccrcc_pathology.tcia"),
)
RAW_DIR = os.path.join(os.environ.get("MORPHO_WSI_DIR", "./WSI"), "raw", "ccrcc")


def parse_manifest(path):
    uids = []
    in_list = False
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line == "ListOfSeriesToDownload=":
                in_list = True
                continue
            if in_list and line:
                uids.append(line)
    return uids


def download_series(uid, dest_dir, retries=4, timeout=180):
    out_dir = os.path.join(dest_dir, uid)
    marker = os.path.join(out_dir, ".done")
    if os.path.exists(marker):
        return "skip"

    url = f"{API_URL}?SeriesInstanceUID={uid}"
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            with zipfile.ZipFile(io.BytesIO(data)) as zf:
                os.makedirs(out_dir, exist_ok=True)
                zf.extractall(out_dir)
            open(marker, "w").close()
            n_dcm = sum(1 for n in os.listdir(out_dir) if n.lower().endswith(".dcm"))
            size_mb = len(data) / 1e6
            return f"ok ({n_dcm} dcm, {size_mb:.1f} MB)"
        except (urllib.error.URLError, zipfile.BadZipFile, TimeoutError) as e:
            print(f"    attempt {attempt}/{retries} failed for {uid}: {e}")
            time.sleep(3 * attempt)
    return "FAILED"


def main():
    if not os.path.exists(MANIFEST_PATH):
        sys.exit(f"manifest not found: {MANIFEST_PATH}")
    uids = parse_manifest(MANIFEST_PATH)
    print(f"{len(uids)} series in manifest")
    os.makedirs(RAW_DIR, exist_ok=True)

    ok = skip = fail = 0
    total_bytes_est = 0
    for i, uid in enumerate(uids, 1):
        result = download_series(uid, RAW_DIR)
        if result == "skip":
            skip += 1
            status = "skip (already done)"
        elif result == "FAILED":
            fail += 1
            status = "FAILED"
        else:
            ok += 1
            status = result
        print(f"[{i}/{len(uids)}] {uid}: {status}")

    print(f"\ndone. ok={ok} skip={skip} fail={fail} total={len(uids)}")
    print(f"raw slides at: {RAW_DIR}")
    if fail:
        print("re-run this script to retry only the failed ones (skip logic handles the rest)")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
