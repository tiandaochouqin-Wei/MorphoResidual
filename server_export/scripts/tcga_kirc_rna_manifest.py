#!/usr/bin/env python3
"""tcga_kirc_rna_manifest.py — write a gdc-client-style manifest for TCGA-KIRC
tumour RNA-seq (STAR counts) + a case<->file crosswalk for the residual analysis.
Runs on mn02.  python tcga_kirc_rna_manifest.py
Then:  nohup python /public/home/fjhui/ZW/scripts/gdc_download.py \
             /public/home/fjhui/ZW/tcga_kirc/manifest_kirc_rna.txt \
             /public/home/fjhui/ZW/tcga_kirc/rna 6 > rna_dl.log 2>&1 &"""
import io, os, sys
try:
    import requests
except Exception as e:
    sys.exit(f"requests unavailable: {e}")

OUT = "/public/home/fjhui/ZW/tcga_kirc"
os.makedirs(OUT, exist_ok=True)
FILES = "https://api.gdc.cancer.gov/files"


def f_and(*c): return {"op": "and", "content": list(c)}
def f_in(f, v): return {"op": "in", "content": {"field": f, "value": v}}

hits = requests.post(FILES, json={
    "filters": f_and(
        f_in("cases.project.project_id", ["TCGA-KIRC"]),
        f_in("data_type", ["Gene Expression Quantification"]),
        f_in("analysis.workflow_type", ["STAR - Counts"])),
    "size": 5000, "format": "json",
    "fields": "file_id,file_name,md5sum,file_size,cases.submitter_id,cases.samples.sample_type"},
    timeout=180).json()["data"]["hits"]
print(f"{len(hits)} STAR files")

man = f"{OUT}/manifest_kirc_rna.txt"
xwalk = f"{OUT}/rna_case_file.tsv"
tot = ntum = 0
with io.open(man, "w") as fm, io.open(xwalk, "w") as fx:
    fm.write("id\tfilename\tmd5\tsize\tstate\n")
    fx.write("case\tsample_type\tfilename\n")
    for h in hits:
        cs = (h.get("cases") or [{}])[0]
        case = cs.get("submitter_id")
        stype = (cs.get("samples", [{}])[0] or {}).get("sample_type", "")
        if not case or "Normal" in stype:            # tumour only
            continue
        fm.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum','')}\t{h.get('file_size',0)}\treleased\n")
        fx.write(f"{case}\t{stype}\t{h['file_name']}\n")
        tot += h.get("file_size", 0) or 0
        ntum += 1
print(f"[rna] {ntum} tumour STAR files, {tot/1e9:.2f} GB")
print(f"  manifest -> {man}")
print(f"  crosswalk -> {xwalk}")
