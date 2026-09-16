#!/usr/bin/env python3
"""Probe the GDC files API to work out the exact shape of Slide Image records
for the CPTAC-3 UCEC / PDAC cohorts before writing the real slide_type_map builder."""
import glob
import json
import os
import re
import sys
import urllib.request

GDC = "https://api.gdc.cancer.gov/files"


def q(payload, timeout=90):
    req = urllib.request.Request(
        GDC, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def local_cases(cancer):
    emb = f"/public/home/fjhui/ZW/{cancer}/WSI/emb"
    cases = set()
    slides = set()
    for p in glob.glob(emb + "/*.pt"):
        sid = os.path.basename(p)[:-3]
        slides.add(sid)
        m = re.match(r"^(C3[A-Z]-\d{5})", sid)
        if m:
            cases.add(m.group(1))
    return cases, slides


for cancer in ["ucec", "pdac"]:
    cases, slides = local_cases(cancer)
    print(f"\n########## {cancer}: {len(slides)} local slides, {len(cases)} cases")
    payload = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id",
                                         "value": ["CPTAC-3"]}},
                {"op": "in", "content": {"field": "cases.submitter_id",
                                         "value": sorted(cases)}},
                {"op": "in", "content": {"field": "data_type",
                                         "value": ["Slide Image"]}},
            ],
        },
        "fields": ",".join([
            "file_id", "file_name", "access",
            "cases.submitter_id",
            "cases.primary_site",
            "cases.samples.sample_type",
            "cases.samples.submitter_id",
            "cases.samples.portions.slides.submitter_id",
            "cases.samples.portions.slides.section_location",
        ]),
        "size": "5000",
    }
    d = q(payload)
    hits = d["data"]["hits"]
    total = d["data"]["pagination"]["total"]
    print(f"GDC returned {len(hits)} hits (total={total})")
    print("--- first 2 raw hits ---")
    print(json.dumps(hits[:2], indent=2)[:2500])

    # quick sample_type tally by file
    tally = {}
    for h in hits:
        for c in h.get("cases", []):
            for s in c.get("samples", []):
                tally[s.get("sample_type")] = tally.get(s.get("sample_type"), 0) + 1
    print("--- sample_type tally (per file x sample) ---", tally)
    sites = set()
    for h in hits:
        for c in h.get("cases", []):
            ps = c.get("primary_site")
            if isinstance(ps, list):
                sites.update(ps)
            elif ps:
                sites.add(ps)
    print("--- primary sites ---", sites)
