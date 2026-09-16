#!/usr/bin/env python3
"""Third probe: pull slide entities from the GDC *cases* (biospecimen) endpoint.
GDC no longer hosts CPTAC-3 Slide Image files, but the slide entities still live
under samples.portions.slides -- that is where section_location comes from."""
import json
import urllib.request

CASES = "https://api.gdc.cancer.gov/cases"


def q(payload, timeout=120):
    req = urllib.request.Request(
        CASES, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


payload = {
    "filters": {"op": "in", "content": {"field": "submitter_id",
                                        "value": ["C3L-00006", "C3L-00017", "C3N-00646"]}},
    "expand": "samples,samples.portions,samples.portions.slides",
    "fields": "submitter_id,primary_site,project.project_id",
    "size": "10",
}
d = q(payload)
hits = d["data"]["hits"]
print(f"total={d['data']['pagination']['total']} shown={len(hits)}")

for h in hits:
    print("\n" + "=" * 70)
    print(f"case {h.get('submitter_id')}  site={h.get('primary_site')}")
    for s in h.get("samples", []):
        print(f"  sample {s.get('submitter_id')}  type={s.get('sample_type')!r}")
        for p in s.get("portions", []):
            for sl in p.get("slides", []):
                print(f"      slide {sl.get('submitter_id')!r} "
                      f"section_location={sl.get('section_location')!r} "
                      f"pct_tumor={sl.get('percent_tumor_cells')!r}")

print("\n\n--- raw first hit (truncated) ---")
print(json.dumps(hits[0], indent=2)[:3000])
