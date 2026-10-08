#!/usr/bin/env python3
"""vial_match.py -- metadata only. For the 102 DX-arm cases: which vial / portion does each
modality sit on? PDC MS aliquot (PDC biospecimenPerStudy), GDC STAR-Counts RNA file sample,
GDC DX slides and TS/BS/MS slides. Reads inventory_ts_files.tsv + inventory_dx_cases.tsv
(written by gdc_inventory_ts.py). Writes vial_match.tsv + vial_match_summary.json.
No molecular / morphology statistic is computed."""
import csv, collections, json, os, re, time
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
S = requests.Session()


def vial(b):
    m = re.match(r"(TCGA-\w\w-\w{4})-(\d\d[A-Z])(?:-(\d\d))?", b or "")
    return (m.group(2), m.group(3)) if m else (None, None)


cases = [l.strip() for l in open(f"{HERE}/inventory_dx_cases.tsv", encoding="utf-8").read().split("\n")[1:] if l.strip()]
q = '{ biospecimenPerStudy(pdc_study_id: "PDC000173") { case_submitter_id sample_submitter_id sample_type aliquot_submitter_id } }'
pdc = S.post("https://pdc.cancer.gov/graphql", json={"query": q}, timeout=60).json()["data"]["biospecimenPerStudy"]
ms = collections.defaultdict(list)
for r in pdc:
    if "Tumor" in (r.get("sample_type") or ""):
        ms[r["case_submitter_id"]].append(r["aliquot_submitter_id"])

# RNA files for these cases
F = {"op": "and", "content": [
    {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]}},
    {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
    {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
    {"op": "in", "content": {"field": "cases.submitter_id", "value": cases}}]}
r = S.post("https://api.gdc.cancer.gov/files", json={"filters": F, "size": 5000, "format": "json",
           "fields": "file_name,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type"}, timeout=180).json()
rna = collections.defaultdict(list)
for h in r["data"]["hits"]:
    cs = h["cases"][0]
    for s in cs.get("samples", []):
        if "Normal" not in (s.get("sample_type") or ""):
            rna[cs["submitter_id"]].append(vial(s["submitter_id"])[0])

ts = collections.defaultdict(list); dx = collections.defaultdict(list)
for row in csv.DictReader(open(f"{HERE}/inventory_ts_files.tsv", encoding="utf-8"), delimiter="\t"):
    (ts if row["strategy"] == "Tissue Slide" else dx)[row["case"]].append(row)

out = []
cnt = collections.Counter()
for c in cases:
    ms_v = sorted({vial(a)[0] for a in ms[c]})
    rna_v = sorted(set(rna[c]))
    ts_prim = [x for x in ts[c] if x["sample_code"].startswith("01")]
    ts_v = sorted({x["sample_code"] for x in ts_prim})
    dx_v = sorted({re.match(r"TCGA-\w\w-\w{4}-(\d\d[A-Z])", x["file_name"]).group(1) for x in dx[c]})
    same = bool(set(ms_v) & set(ts_v))
    cnt["ms_vial_has_TS"] += same
    cnt["ms_vial_has_TS_and_BS"] += (same and {"T", "B"} <= {x["slide_label"][0] for x in ts_prim if x["sample_code"] in ms_v})
    cnt["n_ts_prim_ge1"] += bool(ts_prim)
    cnt["ms_aliquots_gt1"] += len(ms[c]) > 1
    out.append(dict(case=c, ms_aliquots="|".join(ms[c]), ms_vials=",".join(ms_v), rna_vials=",".join(rna_v),
                    dx_vials=",".join(dx_v), ts_primary_vials=",".join(ts_v),
                    ts_primary_labels=",".join(sorted(x["slide_label"] + "@" + x["sample_code"] for x in ts_prim)),
                    ts_nonprimary=",".join(sorted(x["slide_label"] + "@" + x["sample_code"] for x in ts[c] if not x["sample_code"].startswith("01"))),
                    ms_vial_has_ts=same))
with open(f"{HERE}/vial_match.tsv", "w", encoding="utf-8") as fh:
    fh.write("\t".join(out[0]) + "\n")
    for o in out:
        fh.write("\t".join(str(v) for v in o.values()) + "\n")
summ = dict(cnt)
summ["ms_vial_codes"] = dict(collections.Counter(o["ms_vials"] for o in out))
summ["rna_vial_codes"] = dict(collections.Counter(o["rna_vials"] for o in out))
summ["dx_vial_codes"] = dict(collections.Counter(o["dx_vials"] for o in out))
summ["ts_primary_vial_codes"] = dict(collections.Counter(o["ts_primary_vials"] for o in out))
summ["cases_with_06_ts"] = [o["case"] for o in out if "@06" in o["ts_nonprimary"]]
summ["cases_with_11_ts"] = sum("@11" in o["ts_nonprimary"] for o in out)
json.dump(summ, open(f"{HERE}/vial_match_summary.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print(json.dumps(summ, indent=1, ensure_ascii=False))
