import json, sys, urllib.request
pdc_id, cancer = sys.argv[1], sys.argv[2]
out = f"/public/home/fjhui/ZW/scripts/aliquot_to_case_tumor_{cancer}.tsv"
q = '{ biospecimenPerStudy(pdc_study_id: "%s" acceptDUA: true) { aliquot_submitter_id case_submitter_id sample_type } }' % pdc_id
req = urllib.request.Request("https://pdc.cancer.gov/graphql",
    data=json.dumps({"query": q}).encode(), headers={"Content-Type": "application/json"})
rows = json.load(urllib.request.urlopen(req, timeout=90))["data"]["biospecimenPerStudy"]
seen = {}
for r in rows:
    a = r["aliquot_submitter_id"]; st = (r["sample_type"] or "").strip()
    if not a or "Normal" in st or a in seen:
        continue
    seen[a] = (r["case_submitter_id"], st)
with open(out, "w") as f:
    f.write("aliquot_id\tcase_id\tsample_type\n")
    for a,(c,st) in seen.items(): f.write(f"{a}\t{c}\t{st}\n")
print(f"{cancer}: wrote {len(seen)} tumor aliquots -> {out}")
