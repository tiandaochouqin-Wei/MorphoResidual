#!/usr/bin/env python3
"""Adjudicate the 3 LUAD mismatches: what does the authoritative GDC slide
entity (samples.portions.slides) say the parent sample of each slide is?"""
import json
import urllib.request

CASES = "https://api.gdc.cancer.gov/cases"


def q(case_ids):
    payload = {
        "filters": {"op": "in", "content": {"field": "submitter_id", "value": case_ids}},
        "expand": "samples,samples.portions,samples.portions.slides",
        "fields": "submitter_id,primary_site",
        "size": "50",
    }
    req = urllib.request.Request(
        CASES, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)["data"]["hits"]


for h in q(["C3L-02365", "C3N-02000", "C3L-00006", "C3L-00017"]):
    print("\n" + "=" * 72)
    print(f"case {h['submitter_id']}  site={h.get('primary_site')}")
    for s in h.get("samples", []):
        slides = []
        for p in s.get("portions", []):
            for sl in p.get("slides", []):
                slides.append(sl.get("submitter_id"))
        print(f"  sample {s.get('submitter_id'):<18} type={s.get('sample_type'):<22} "
              f"slides={slides}")
