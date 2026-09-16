#!/usr/bin/env python3
"""Pull CPTAC-3 case-level clinical metadata for all five cohorts from the public
GDC API, locally (no HPC round trip).

Three things depend on this:
  1. UCEC endometrioid-only sensitivity. FIGO grade is defined only for the
     endometrioid histotype; serous, clear-cell and carcinosarcoma cases are
     high-grade by definition. The manuscript says so but has never actually
     restricted to endometrioid, because the server-side clinical file carries no
     histology column. primary_diagnosis supplies it.
  2. Whether the age field really is unpopulated in all five cohorts. The
     manuscript asserts this (Methods, clinical association) on the basis of the
     server-derived files; asserting it from the authoritative source is better,
     and if it turns out to be populated the adjusted Cox models must be redone.
  3. Sex reporting. Nature Communications requires the sex/gender of human
     participants to be stated and disaggregated where relevant.

Writes figures/figdata/gdc_clinical_<cohort>.csv and prints a summary.
"""
import json
import time
import urllib.parse
import urllib.request

import pandas as pd

SITES = {                       # GDC primary_site strings, per the pan-cancer extension notes
    "ccrcc": "Kidney",
    "luad": "Bronchus and lung",
    "ucec": "Uterus, NOS",
    "gbm": "Brain",
    "pdac": "Pancreas",
}
FIELDS = [
    "submitter_id",
    "disease_type",
    "diagnoses.primary_diagnosis",
    "diagnoses.tumor_grade",
    "diagnoses.figo_stage",
    "diagnoses.ajcc_pathologic_stage",
    "demographic.sex_at_birth",
    "demographic.year_of_birth",
    "demographic.days_to_birth",
    "demographic.age_at_index",
    "demographic.race",
    "demographic.ethnicity",
    "demographic.vital_status",
]
API = "https://api.gdc.cancer.gov/cases"


def fetch(site):
    flt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "project.project_id", "value": ["CPTAC-3"]}},
        {"op": "in", "content": {"field": "primary_site", "value": [site]}}]}
    url = (f"{API}?filters={urllib.parse.quote(json.dumps(flt))}"
           f"&fields={','.join(FIELDS)}&size=2000&format=json")
    with urllib.request.urlopen(url, timeout=120) as fh:
        r = json.load(fh)
    if r.get("warnings"):
        print(f"  [warn] {r['warnings']}")
    rows = []
    for h in r["data"]["hits"]:
        dg = (h.get("diagnoses") or [{}])[0]
        dm = h.get("demographic") or {}
        rows.append(dict(
            case=h.get("submitter_id"),
            disease_type=h.get("disease_type"),
            primary_diagnosis=dg.get("primary_diagnosis"),
            tumor_grade=dg.get("tumor_grade"),
            figo_stage=dg.get("figo_stage"),
            ajcc_stage=dg.get("ajcc_pathologic_stage"),
            sex_at_birth=dm.get("sex_at_birth"),
            year_of_birth=dm.get("year_of_birth"),
            days_to_birth=dm.get("days_to_birth"),
            age_at_index=dm.get("age_at_index"),
            race=dm.get("race"),
            ethnicity=dm.get("ethnicity"),
            vital_status=dm.get("vital_status"),
        ))
    return pd.DataFrame(rows)


for c, site in SITES.items():
    print(f"\n=== {c.upper()}  (primary_site='{site}') ===")
    d = fetch(site)
    d.to_csv(f"figures/figdata/gdc_clinical_{c}.csv", index=False)

    # restrict to the cases this study actually analysed
    an = pd.read_csv(f"figures/figdata/clinical_link_{c}_clinical.csv")
    m = d[d.case.isin(an.case)]
    print(f"  GDC returned {len(d)} CPTAC-3 cases at this site; "
          f"{len(m)}/{len(an)} of the analysed cases matched")
    print(f"  age_at_index non-null: {m.age_at_index.notna().sum()}/{len(m)}")
    print(f"  sex_at_birth: {m.sex_at_birth.value_counts(dropna=False).to_dict()}")
    print(f"  year_of_birth non-null: {m.year_of_birth.notna().sum()}/{len(m)}; "
          f"days_to_birth non-null: {m.days_to_birth.notna().sum()}/{len(m)}")
    print(f"  race: {m.race.value_counts(dropna=False).to_dict()}")
    vc = m.primary_diagnosis.value_counts(dropna=False)
    print("  primary_diagnosis:")
    for k, v in vc.items():
        print(f"    {v:4d}  {k}")
    print(f"  wrote figures/figdata/gdc_clinical_{c}.csv")
    time.sleep(0.5)
