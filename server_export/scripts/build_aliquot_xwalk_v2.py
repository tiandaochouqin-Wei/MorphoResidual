#!/usr/bin/env python3
"""
Build the tumor-only aliquot -> case crosswalk that residual_analysis.py uses to
pick protein columns.

Why this replaces build_aliquot_xwalk.py: that version kept every aliquot whose
sample_type did not literally contain "Normal". For GBM that quietly admitted one
"Not Reported" aliquot; for PDAC (PDC000270) it would admit 6 "Cell Lines", 5
"Tumor" and 3 "Not Reported" alongside the 145 Primary Tumor aliquots. Cell-line
channels are not patient tissue and must not enter a patient-level tumor cohort,
and load_protein_matrix() averages all aliquots that map to the same case, so a
stray channel silently contaminates that patient's protein vector.

This version keeps sample_type == "Primary Tumor" and nothing else, and reports
exactly what it dropped so the choice is auditable.

Verified study ids (protein-column coverage / WSI-case overlap measured live):
  ccrcc PDC000127 | luad PDC000153-family | gbm PDC000204 (110/110, 100/100)
  ucec  PDC000125 (149/149 cols, 102/102 cases)
  pdac  PDC000270 (222/224 cols, 140/140 cases)

Usage: python build_aliquot_xwalk_v2.py <PDC_ID> <cancer>
"""
import json
import sys
import urllib.request

KEEP = {"Primary Tumor"}

# scripts/ is shared with other people on this account and the legacy loose-filter
# build_aliquot_xwalk.py writes the same filenames; on 2026-09-01 a concurrent run
# of it put 6 Cell Lines, 3 Not Reported and 5 "Tumor" aliquots back into the PDAC
# crosswalk mid-session. MORPHO_MAP_DIR points this at a private directory so the
# pipeline's inputs cannot be replaced underneath it.
import os
OUTDIR = os.environ.get("MORPHO_MAP_DIR", "/public/home/fjhui/ZW/scripts")

pdc_id, cancer = sys.argv[1], sys.argv[2]
os.makedirs(OUTDIR, exist_ok=True)
out = f"{OUTDIR}/aliquot_to_case_tumor_{cancer}.tsv"

q = ('{ biospecimenPerStudy(pdc_study_id: "%s" acceptDUA: true) '
     '{ aliquot_submitter_id case_submitter_id sample_type } }' % pdc_id)
req = urllib.request.Request(
    "https://pdc.cancer.gov/graphql", data=json.dumps({"query": q}).encode(),
    headers={"Content-Type": "application/json"})
d = json.load(urllib.request.urlopen(req, timeout=180))
if d.get("errors"):
    sys.exit(f"PDC graphql error: {d['errors']}")
rows = d["data"]["biospecimenPerStudy"] or []

dropped = {}
seen = {}
for r in rows:
    a = r.get("aliquot_submitter_id")
    st = (r.get("sample_type") or "").strip()
    if not a:
        continue
    if st not in KEEP:
        dropped[st] = dropped.get(st, 0) + 1
        continue
    if a in seen:
        continue
    seen[a] = (r.get("case_submitter_id"), st)

cases = {}
for a, (c, _st) in seen.items():
    cases.setdefault(c, []).append(a)
multi = {c: v for c, v in cases.items() if len(v) > 1}

with open(out, "w") as f:
    f.write("aliquot_id\tcase_id\tsample_type\n")
    for a, (c, st) in sorted(seen.items()):
        f.write(f"{a}\t{c}\t{st}\n")

print(f"{cancer} ({pdc_id}): {len(rows)} biospecimen rows")
print(f"  kept {len(seen)} Primary Tumor aliquots over {len(cases)} cases -> {out}")
print(f"  dropped by sample_type: {dropped}")
if multi:
    print(f"  NOTE: {len(multi)} cases have >1 tumor aliquot "
          f"(load_protein_matrix averages them): "
          f"{dict(list(multi.items())[:5])}")
