#!/usr/bin/env python3
"""Find the correct PDC study ids for the UCEC and PDAC proteome studies, and
confirm the ones already used for GBM/PDAC/LUAD/ccRCC."""
import json
import urllib.request

PDC = "https://pdc.cancer.gov/graphql"


def q(query, timeout=120):
    req = urllib.request.Request(
        PDC, data=json.dumps({"query": query}).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


d = q('{ allPrograms { program_name projects { project_name studies { '
      'pdc_study_id study_id submitter_id_name analytical_fraction experiment_type } } } }')

rows = []
for p in d["data"]["allPrograms"]:
    for proj in p.get("projects") or []:
        for s in proj.get("studies") or []:
            rows.append((s.get("pdc_study_id"), s.get("submitter_id_name"),
                         s.get("analytical_fraction"), s.get("experiment_type"),
                         proj.get("project_name")))

keys = ["endometrial", "uterine", "ucec", "pancrea", "pda", "glioblastoma",
        "lung adeno", "renal"]
print(f"{len(rows)} studies total\n")
for k in keys:
    print(f"--- match {k!r} ---")
    for r in sorted(set(rows)):
        name = (r[1] or "").lower()
        if k in name and (r[2] or "") == "Proteome":
            print(f"  {r[0]:<12} {r[1]:<60} {r[2]:<10} {r[3]:<8} {r[4]}")
    print()
