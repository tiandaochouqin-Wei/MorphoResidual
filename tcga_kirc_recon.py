#!/usr/bin/env python3
"""tcga_kirc_recon.py — B-11 feasibility: how many TCGA-KIRC cases have RPPA + RNA +
diagnostic H&E slides, and how big is each modality. Runs on mn02 (GDC API).
    python tcga_kirc_recon.py"""
import json, sys
try:
    import requests
except Exception as e:
    sys.exit(f"requests unavailable: {e}")

GDC = "https://api.gdc.cancer.gov/files"


def cases_for(filters, label):
    body = {"filters": filters, "size": 20000, "format": "json",
            "fields": "cases.submitter_id,file_size,data_type,experimental_strategy"}
    r = requests.post(GDC, json=body, timeout=120)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    cases, total_bytes = set(), 0
    for h in hits:
        total_bytes += h.get("file_size", 0) or 0
        for c in h.get("cases", []) or []:
            if c.get("submitter_id"):
                cases.add(c["submitter_id"])
    print(f"[{label}] files={len(hits)}  cases={len(cases)}  total={total_bytes/1e9:.1f} GB")
    return cases


def f_and(*content):
    return {"op": "and", "content": list(content)}


def f_in(field, values):
    return {"op": "in", "content": {"field": field, "value": values}}


PROJ = f_in("cases.project.project_id", ["TCGA-KIRC"])

# 1) diagnostic H&E slides (SVS)
slides = cases_for(f_and(PROJ,
    f_in("data_type", ["Slide Image"]),
    f_in("experimental_strategy", ["Diagnostic Slide"])), "H&E diagnostic slides")

# 2) RNA-seq STAR counts
rna = cases_for(f_and(PROJ,
    f_in("data_type", ["Gene Expression Quantification"]),
    f_in("analysis.workflow_type", ["STAR - Counts"])), "RNA-seq (STAR)")

# 3) RPPA protein expression
rppa = cases_for(f_and(PROJ,
    f_in("data_type", ["Protein Expression Quantification"])), "RPPA protein")

tri = slides & rna & rppa
print("\n==================== B-11 FEASIBILITY (TCGA-KIRC) ====================")
print(f" slides {len(slides)} | RNA {len(rna)} | RPPA {len(rppa)}")
print(f" cases with ALL THREE (RPPA + RNA + H&E slide): {len(tri)}")
print(" -> if >~150 tri-modal cases, B-11 is well-powered; the slide set is the")
print("    download cost. RPPA antibody/gene count is checked once a file is pulled.")
