#!/usr/bin/env python3
"""
Download matched RNA-seq (gene counts) + WES (masked somatic mutation) for the
CPTAC-CCRCC pilot cohort directly from GDC, no gdc-client binary needed.

Verified live against https://api.gdc.cancer.gov on 2026-07-13:
  - CPTAC-3 has no per-cancer disease_type label; isolate CCRCC via
    cases.primary_site = "Kidney" (261 cases under CPTAC-3).
  - RNA-seq: data_type="Gene Expression Quantification", workflow_type=
    "STAR - Counts" -> access="open" (521 files). The sibling
    "Splice Junction Quantification" files under the same workflow are
    access="controlled" and will 403 -- do NOT widen the data_type filter.
  - WES: data_category="Simple Nucleotide Variation", data_type=
    "Masked Somatic Mutation" (353 files, per GDC convention this masked
    tier is the open-access one; verify access="open" per-file below anyway).
  - Direct byte download works via GET https://api.gdc.cancer.gov/data/<uuid>
    for access="open" files -- no auth, no gdc-client required.
"""
import json
import os
import sys
import time
import urllib.request
import urllib.error

GDC_FILES_API = "https://api.gdc.cancer.gov/files"
GDC_DATA_API = "https://api.gdc.cancer.gov/data"

OUT_ROOT = os.environ.get("MORPHO_OMICS_DIR", "./omics")
RNA_DIR = os.path.join(OUT_ROOT, "rna")
WES_DIR = os.path.join(OUT_ROOT, "mutation")

RNA_FILTERS = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["CPTAC-3"]}},
        {"op": "in", "content": {"field": "cases.primary_site", "value": ["Kidney"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
    ],
}

WES_FILTERS = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["CPTAC-3"]}},
        {"op": "in", "content": {"field": "cases.primary_site", "value": ["Kidney"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
    ],
}

FIELDS = "file_id,file_name,file_size,access,cases.submitter_id,cases.samples.sample_type"


def gdc_query(filters, size=2000):
    payload = json.dumps({
        "filters": filters,
        "fields": FIELDS,
        "size": str(size),
    }).encode()
    req = urllib.request.Request(GDC_FILES_API, data=payload,
                                  headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = json.load(resp)
    hits = data["data"]["hits"]
    total = data["data"]["pagination"]["total"]
    if total > size:
        print(f"  WARNING: total={total} > page size={size}, results truncated; raise `size`.")
    return hits


def download_file(file_id, file_name, dest_dir, expected_size=None, retries=4):
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, file_name)
    if os.path.exists(dest) and (expected_size is None or os.path.getsize(dest) == int(expected_size)):
        print(f"  skip (already downloaded): {file_name}")
        return True
    url = f"{GDC_DATA_API}/{file_id}"
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=120) as resp, open(dest, "wb") as f:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            return True
        except urllib.error.HTTPError as e:
            body = e.read().decode(errors="replace")
            print(f"  attempt {attempt}/{retries} HTTP {e.code} for {file_name}: {body[:200]}")
            if e.code == 403:
                print(f"  -> file is not open-access, skipping permanently: {file_name}")
                return False
        except Exception as e:
            print(f"  attempt {attempt}/{retries} failed for {file_name}: {e}")
        time.sleep(3 * attempt)
    print(f"  GAVE UP on {file_name}")
    return False


def run(label, filters, dest_dir, manifest_path):
    print(f"\n=== {label}: querying GDC ===")
    hits = gdc_query(filters)
    print(f"  {len(hits)} files found")
    open_hits = [h for h in hits if h.get("access") == "open"]
    if len(open_hits) < len(hits):
        print(f"  skipping {len(hits) - len(open_hits)} controlled-access files")

    with open(manifest_path, "w", encoding="utf-8") as mf:
        mf.write("file_id\tfile_name\tfile_size\tcase_submitter_id\n")
        for h in open_hits:
            case_id = h["cases"][0]["submitter_id"] if h.get("cases") else "NA"
            mf.write(f"{h['file_id']}\t{h['file_name']}\t{h['file_size']}\t{case_id}\n")
    print(f"  manifest written: {manifest_path}")

    ok, fail = 0, 0
    for i, h in enumerate(open_hits, 1):
        print(f"[{i}/{len(open_hits)}] {h['file_name']} ({int(h['file_size'])/1e6:.1f} MB)")
        if download_file(h["file_id"], h["file_name"], dest_dir, h.get("file_size")):
            ok += 1
        else:
            fail += 1
    print(f"=== {label}: done. ok={ok} fail={fail} ===")
    return ok, fail


if __name__ == "__main__":
    os.makedirs(OUT_ROOT, exist_ok=True)
    rna_ok, rna_fail = run("RNA-seq gene counts (CPTAC-CCRCC pilot)", RNA_FILTERS,
                            RNA_DIR, os.path.join(OUT_ROOT, "manifest_rna.tsv"))
    wes_ok, wes_fail = run("WES masked somatic mutation (CPTAC-CCRCC pilot)", WES_FILTERS,
                            WES_DIR, os.path.join(OUT_ROOT, "manifest_wes.tsv"))
    total_fail = rna_fail + wes_fail
    print(f"\nTOTAL: rna_ok={rna_ok} wes_ok={wes_ok} failures={total_fail}")
    sys.exit(1 if total_fail else 0)
