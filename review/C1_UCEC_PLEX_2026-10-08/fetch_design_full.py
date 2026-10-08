#!/usr/bin/env python3
"""ADAPTED COPY (2026-10-08) of review/BATCH_AXES_SCOPING_2026-10-02/fetch_design_full.py.
Only STUDIES changed: pulls the full per-run, per-channel TMT experimental design for the
C1-UCEC confirmatory proteome study PDC000439 (metadata only). Study ID source:
review/C1_UCEC_FROZEN_RULE_2026-09-15.md:82 and server_export/scripts/c1_pull_omics.py:233.
"""
import json
import os

from pdc_gql import gql, HERE

STUDIES = {"ucec_c1": "PDC000439"}
SUB = "{ aliquot_id aliquot_submitter_id aliquot_run_metadata_id }"
TMT = ["tmt_126", "tmt_127n", "tmt_127c", "tmt_128n", "tmt_128c", "tmt_129n",
       "tmt_129c", "tmt_130n", "tmt_130c", "tmt_131", "tmt_131c"]
RUN = ["study_id", "pdc_study_id", "study_submitter_id", "study_run_metadata_id",
       "study_run_metadata_submitter_id", "plex_dataset_name", "experiment_number",
       "experiment_type", "analyte", "acquisition_type", "number_of_fractions",
       "aliquot_is_ref"]


def main():
    fields = " ".join(RUN) + " " + " ".join(f"{c} {SUB}" for c in TMT)
    summary = {}
    for coh, pid in STUDIES.items():
        q = f'{{ studyExperimentalDesign(pdc_study_id: "{pid}", acceptDUA: true) {{ {fields} }} }}'
        st, d = gql(q, f"r2_sed_full_{pid}", f"studyExperimentalDesign({pid}) [{coh}] run fields + tmt_126..tmt_131c {SUB}")
        errs = d.get("errors")
        rows = (d.get("data") or {}).get("studyExperimentalDesign") or []
        print(coh, pid, st, errs[0]["message"][:200] if errs else f"rows={len(rows)}")
        if rows:
            types = sorted({str(r.get("experiment_type")) for r in rows})
            filled = {c: sum(1 for r in rows if r.get(c)) for c in TMT}
            summary[coh] = {"pdc_study_id": pid, "n_runs": len(rows), "experiment_type": types,
                            "runs_with_channel_filled": filled}
            print("   types", types, "filled", filled)
    with open(os.path.join(HERE, "design_fetch_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)


if __name__ == "__main__":
    main()
