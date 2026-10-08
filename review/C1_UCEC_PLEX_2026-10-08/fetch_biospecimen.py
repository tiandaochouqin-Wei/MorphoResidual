#!/usr/bin/env python3
"""ADAPTED COPY (2026-10-08) of the scoping fetch_biospecimen.py; only STUDIES changed
(C1-UCEC confirmatory PDC000439). Original: pull biospecimenPerStudy (aliquot -> case, sample_type, status) for the six
studies -- the same query build_aliquot_xwalk_v2.py uses to build the pinned
tumor crosswalks, plus aliquot_status (field proven valid in Step 0b)."""
from pdc_gql import gql

STUDIES = {"ucec_c1": "PDC000439"}
for coh, pid in STUDIES.items():
    q = ('{ biospecimenPerStudy(pdc_study_id: "%s" acceptDUA: true) '
         '{ aliquot_id aliquot_submitter_id case_submitter_id sample_type aliquot_status } }' % pid)
    st, d = gql(q, f"r3_biospecimen_{pid}", f"biospecimenPerStudy({pid}) [{coh}] aliquot/case/sample_type/status")
    errs = d.get("errors")
    rows = (d.get("data") or {}).get("biospecimenPerStudy") or []
    print(coh, pid, st, errs[0]["message"][:200] if errs else f"rows={len(rows)}")
