#!/usr/bin/env python3
"""Second GDC probe: find out how CPTAC-3 slide images are actually indexed."""
import json
import urllib.request

GDC = "https://api.gdc.cancer.gov/files"


def q(payload, timeout=120):
    req = urllib.request.Request(
        GDC, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def head(label, payload):
    try:
        d = q(payload)
    except Exception as e:
        print(f"\n### {label}: ERROR {e}")
        return None
    hits = d["data"]["hits"]
    print(f"\n### {label}: total={d['data']['pagination']['total']} shown={len(hits)}")
    print(json.dumps(hits[:2], indent=2)[:1800])
    return hits


# A. What data_types exist under CPTAC-3 at all? use facets via aggregations
payload_a = {
    "filters": {"op": "in", "content": {"field": "cases.project.project_id",
                                        "value": ["CPTAC-3"]}},
    "facets": "data_type,data_category,experimental_strategy",
    "fields": "file_id",
    "size": "0",
}
try:
    d = q(payload_a)
    for f in ["data_type", "data_category", "experimental_strategy"]:
        buckets = d["data"]["aggregations"][f]["buckets"]
        print(f"\n=== CPTAC-3 {f} ===")
        for b in buckets[:30]:
            print(f"  {b['doc_count']:>8}  {b['key']}")
except Exception as e:
    print("facet query failed:", e)

# B. Slide Image files under CPTAC-3, no case filter
head("CPTAC-3 Slide Image (no case filter)", {
    "filters": {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["CPTAC-3"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Slide Image"]}},
    ]},
    "fields": "file_id,file_name,access,cases.submitter_id,cases.primary_site,"
              "cases.samples.sample_type,cases.samples.portions.slides.submitter_id,"
              "cases.samples.portions.slides.section_location",
    "size": "3",
})

# C. Does a known CCRCC slide id resolve? (that map was built successfully before)
head("file_name = C3N-00646-22.svs", {
    "filters": {"op": "in", "content": {"field": "file_name",
                                        "value": ["C3N-00646-22.svs"]}},
    "fields": "file_id,file_name,data_type,cases.submitter_id,cases.primary_site,"
              "cases.samples.sample_type",
    "size": "3",
})

# D. A known UCEC slide
head("file_name = C3L-00006-21.svs", {
    "filters": {"op": "in", "content": {"field": "file_name",
                                        "value": ["C3L-00006-21.svs"]}},
    "fields": "file_id,file_name,data_type,cases.submitter_id,cases.primary_site,"
              "cases.samples.sample_type",
    "size": "3",
})

# E. Is case C3L-00006 in GDC at all, and under which project/site?
head("case C3L-00006 any file", {
    "filters": {"op": "in", "content": {"field": "cases.submitter_id",
                                        "value": ["C3L-00006"]}},
    "fields": "file_id,file_name,data_type,cases.project.project_id,cases.primary_site",
    "size": "3",
})
