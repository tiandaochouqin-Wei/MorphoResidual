#!/usr/bin/env python3
"""tcga_kirc_pull.py — B-11 step 2: (a) write a gdc-client manifest for the diagnostic
slides of the tri-modal TCGA-KIRC cases; (b) peek at the RPPA file format so we can
build the protein matrix. Runs on mn02 (GDC API).  python tcga_kirc_pull.py"""
import io, sys
try:
    import requests
except Exception as e:
    sys.exit(f"requests unavailable: {e}")

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data/"
OUT = "/public/home/fjhui/ZW/tcga_kirc"
import os; os.makedirs(OUT, exist_ok=True)


def f_and(*c): return {"op": "and", "content": list(c)}
def f_in(f, v): return {"op": "in", "content": {"field": f, "value": v}}
PROJ = f_in("cases.project.project_id", ["TCGA-KIRC"])


def query(filters, fields, size=20000):
    body = {"filters": filters, "size": size, "format": "json", "fields": fields}
    r = requests.post(GDC_FILES, json=body, timeout=180); r.raise_for_status()
    return r.json()["data"]["hits"]


def case_set(filters):
    hits = query(filters, "cases.submitter_id")
    return {c["submitter_id"] for h in hits for c in (h.get("cases") or []) if c.get("submitter_id")}

slides_cases = case_set(f_and(PROJ, f_in("data_type", ["Slide Image"]),
                              f_in("experimental_strategy", ["Diagnostic Slide"])))
rna_cases = case_set(f_and(PROJ, f_in("data_type", ["Gene Expression Quantification"]),
                           f_in("analysis.workflow_type", ["STAR - Counts"])))
rppa_cases = case_set(f_and(PROJ, f_in("data_type", ["Protein Expression Quantification"])))
tri = slides_cases & rna_cases & rppa_cases
print(f"tri-modal cases: {len(tri)}")

# ---- (a) slide manifest for tri-modal cases ----
sl = query(f_and(PROJ, f_in("data_type", ["Slide Image"]),
                 f_in("experimental_strategy", ["Diagnostic Slide"]),
                 f_in("cases.submitter_id", sorted(tri))),
           "file_id,file_name,md5sum,file_size,cases.submitter_id")
man = f"{OUT}/manifest_kirc_slides.txt"
tot = 0
with io.open(man, "w") as fh:
    fh.write("id\tfilename\tmd5\tsize\tstate\n")
    for h in sl:
        fh.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum','')}\t{h.get('file_size',0)}\treleased\n")
        tot += h.get("file_size", 0) or 0
print(f"[manifest] {len(sl)} slides, {tot/1e9:.1f} GB -> {man}")

# ---- (b) peek one RPPA file ----
rp = query(f_and(PROJ, f_in("data_type", ["Protein Expression Quantification"])),
           "file_id,file_name,cases.submitter_id", size=1)
if rp:
    fid = rp[0]["file_id"]
    txt = requests.get(GDC_DATA + fid, timeout=120).text
    print(f"\n[RPPA peek] {rp[0]['file_name']}  (case {rp[0]['cases'][0]['submitter_id']})")
    lines = txt.splitlines()
    print("  header:", lines[0][:200])
    for ln in lines[1:6]:
        print("   ", ln[:120])
    print(f"  total rows (antibodies): {len(lines)-1}")
print(f"\nNEXT: gdc-client download -m {man}  (476 GB, background)")
