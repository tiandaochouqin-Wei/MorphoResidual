#!/usr/bin/env python3
"""ADAPTED COPY (2026-10-08) of the scoping crossval_route2.py. Changes: STUDIES = C1-UCEC
confirmatory PDC000439 only; additionally saves route 2's own case -> (aliquot, label,
experiment_number) table to route2_cases_<PDC>.json so that build_plex_maps.py can derive
case -> run a second, independent way (experiment_number) and compare it with route 1.
Original docstring: Independent second route: paginatedCasesSamplesAliquots ->
aliquots.aliquot_run_metadata{aliquot_run_metadata_id label}. For every
(run, channel, aliquot) of the per-channel design (route 1), check that the same
aliquot_run_metadata_id is attached to that aliquot in route 2 with
label == 'TMT_<channel>'. Pseudo-case aliquots (pools, QC) may be absent from
route 2 because they are not under a real case; they are reported separately."""
import json
import os
from collections import Counter

from pdc_gql import gql, HERE

STUDIES = {"ucec_c1": "PDC000439"}
TMT = ["tmt_126", "tmt_127n", "tmt_127c", "tmt_128n", "tmt_128c", "tmt_129n",
       "tmt_129c", "tmt_130n", "tmt_130c", "tmt_131", "tmt_131c"]
res = {}
for coh, pid in STUDIES.items():
    q = (f'{{ paginatedCasesSamplesAliquots(pdc_study_id: "{pid}", offset: 0, limit: 1000, acceptDUA: true) '
         f'{{ total casesSamplesAliquots {{ case_submitter_id samples {{ sample_type aliquots {{ aliquot_submitter_id '
         f'aliquot_run_metadata {{ aliquot_run_metadata_id label experiment_number }} }} }} }} }} }}')
    st, d = gql(q, f"r5_pcsa_arm_{pid}", f"paginatedCasesSamplesAliquots({pid}) [{coh}] aliquot_run_metadata id/label/experiment_number")
    if d.get("errors"):
        print(coh, "ERR", d["errors"][0]["message"][:300]); continue
    p = d["data"]["paginatedCasesSamplesAliquots"]
    arm = {}   # aliquot_run_metadata_id -> (aliquot_submitter_id, label)
    if len(p["casesSamplesAliquots"]) != p["total"]:
        raise SystemExit(f"{coh}: got {len(p['casesSamplesAliquots'])} of {p['total']} cases (pagination)")
    r2cases = []   # route 2's own case -> aliquot -> run table
    for c in p["casesSamplesAliquots"]:
        for s in c["samples"] or []:
            for a in s["aliquots"] or []:
                for m in a.get("aliquot_run_metadata") or []:
                    arm[m["aliquot_run_metadata_id"]] = (a["aliquot_submitter_id"], m["label"], c["case_submitter_id"])
                    r2cases.append(dict(case=c["case_submitter_id"], sample_type=s["sample_type"],
                                        aliquot=a["aliquot_submitter_id"], label=m["label"],
                                        experiment_number=m.get("experiment_number"),
                                        aliquot_run_metadata_id=m["aliquot_run_metadata_id"]))
    design = json.load(open(os.path.join(HERE, "raw", f"r2_sed_full_{pid}.json"), encoding="utf-8"))
    design = design["data"]["studyExperimentalDesign"]
    tally = Counter()
    mism = []
    for run in design:
        for ch in TMT:
            for e in run.get(ch) or []:
                hit = arm.get(e["aliquot_run_metadata_id"])
                if hit is None:
                    tally["absent_in_route2" + ("_patient" if e["aliquot_submitter_id"].startswith("CPT") else "_pseudo")] += 1
                elif hit[0] == e["aliquot_submitter_id"] and hit[1].lower() == ch.lower():
                    tally["agree"] += 1
                else:
                    tally["DISAGREE"] += 1
                    mism.append((run["plex_dataset_name"], ch, e["aliquot_submitter_id"], hit))
    with open(os.path.join(HERE, f"route2_cases_{pid}.json"), "w", encoding="utf-8") as f:
        json.dump(r2cases, f, indent=1)
    res[coh] = dict(total_cases=p["total"], tally=dict(tally), disagreements=mism[:20])
    print(coh, pid, "cases", p["total"], dict(tally), mism[:5])
with open(os.path.join(HERE, "route2_crossval.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=1)
