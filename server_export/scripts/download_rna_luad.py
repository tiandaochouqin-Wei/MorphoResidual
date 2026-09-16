#!/usr/bin/env python3
"""
LUAD RNA-seq downloader for the external-validation cohort. Reads the prebuilt
manifest_rna_tumor_luad.tsv (111 Primary-Tumor STAR-Counts file_ids, already
filtered to the PDC000153 LUAD case set) and pulls each open-access file directly
from GDC -- no gdc-client, no re-query needed. Run on a node with internet (mn02).

  files -> $MORPHO_LUAD_RNA (default /public/home/fjhui/ZW/luad/omics/rna)

Note: LUAD external validation of contribution A needs only RNA + protein + WSI.
WES/mutation is NOT downloaded (mediation arm C was dropped at pilot scale).
"""
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

GDC_DATA = "https://api.gdc.cancer.gov/data"
MANIFEST = Path(os.environ.get("MORPHO_LUAD_MANIFEST",
                               str(Path(__file__).parent / "manifest_rna_tumor_luad.tsv")))
OUT = Path(os.environ.get("MORPHO_LUAD_RNA", "/public/home/fjhui/ZW/luad/omics/rna"))


def download(file_id, file_name, retries=4):
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / file_name
    if dest.exists() and dest.stat().st_size > 0:
        return "skip"
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(urllib.request.Request(f"{GDC_DATA}/{file_id}"),
                                        timeout=120) as resp, open(dest, "wb") as f:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            return "ok"
        except urllib.error.HTTPError as e:
            print(f"  attempt {attempt} HTTP {e.code} for {file_name}")
            if e.code == 403:
                return "forbidden"
        except Exception as e:
            print(f"  attempt {attempt} failed for {file_name}: {e}")
        time.sleep(3 * attempt)
    return "gaveup"


def main():
    rows = [ln.rstrip("\n").split("\t") for ln in open(MANIFEST, encoding="utf-8")]
    header = {c: i for i, c in enumerate(rows[0])}
    data = rows[1:]
    print(f"LUAD RNA: {len(data)} files -> {OUT}")
    counts = {}
    for i, r in enumerate(data, 1):
        fid, fname = r[header["file_id"]], r[header["file_name"]]
        st = download(fid, fname)
        counts[st] = counts.get(st, 0) + 1
        if i % 20 == 0 or st not in ("ok", "skip"):
            print(f"  [{i}/{len(data)}] {fname} -> {st}")
    print(f"done: {counts}")
    sys.exit(0 if counts.get("gaveup", 0) + counts.get("forbidden", 0) == 0 else 1)


if __name__ == "__main__":
    main()
