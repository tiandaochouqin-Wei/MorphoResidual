#!/usr/bin/env python3
"""Identify the right PDC study id for the UCEC and PDAC proteome tables.

Acceptance test is not "the query returned rows" but "the returned aliquots
actually cover the Log Ratio columns of the protein TSV we have on disk, and
the cases overlap the WSI cohort". A study id that returns data for the wrong
cohort would silently produce an empty 3-modality intersection later.
"""
import glob
import json
import os
import re
import sys
import urllib.request

PDC = "https://pdc.cancer.gov/graphql"
ROOT = "/public/home/fjhui/ZW"

CANDIDATES = {
    "ucec": ["PDC000125", "PDC000126", "PDC000153", "PDC000230"],
    "pdac": ["PDC000270", "PDC000271", "PDC000318"],
    # sanity: gbm's known-good id should pass the same test
    "gbm":  ["PDC000204"],
}

PROTEIN_GLOB = {
    "ucec": f"{ROOT}/ucec/omics/protein/*.tsv",
    "pdac": f"{ROOT}/pdac/omics/protein/*.tsv",
    "gbm":  f"{ROOT}/gbm/omics/protein/*.tsv",
}
EMB_DIR = {
    "ucec": f"{ROOT}/ucec/WSI/emb",
    "pdac": f"{ROOT}/pdac/WSI/emb",
    "gbm":  f"{ROOT}/gbm/WSI/emb",
}


def pdc(query, timeout=180):
    req = urllib.request.Request(
        PDC, data=json.dumps({"query": query}).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def protein_aliquots(cancer):
    fp = sorted(glob.glob(PROTEIN_GLOB[cancer]))[0]
    with open(fp) as f:
        header = f.readline().rstrip("\n").split("\t")
    cols = [c[:-len(" Log Ratio")] for c in header
            if c.endswith(" Log Ratio") and "Unshared" not in c]
    clean = [c for c in cols if not c.startswith("Withdrawn:") and not c.startswith("QC")]
    return os.path.basename(fp), set(cols), set(clean)


def wsi_cases(cancer):
    out = set()
    for p in glob.glob(EMB_DIR[cancer] + "/*.pt"):
        m = re.match(r"^(C3[A-Z]-\d{5})", os.path.basename(p))
        if m:
            out.add(m.group(1))
    return out


for cancer, ids in CANDIDATES.items():
    fname, all_cols, clean_cols = protein_aliquots(cancer)
    cases = wsi_cases(cancer)
    print(f"\n############ {cancer}  protein={fname}")
    print(f"  {len(all_cols)} Log Ratio cols ({len(clean_cols)} after dropping QC/Withdrawn), "
          f"{len(cases)} WSI cases")
    for pid in ids:
        q = ('{ biospecimenPerStudy(pdc_study_id: "%s" acceptDUA: true) '
             '{ aliquot_submitter_id case_submitter_id sample_type } }' % pid)
        try:
            d = pdc(q)
        except Exception as e:
            print(f"  {pid}: request failed {e}")
            continue
        if d.get("errors"):
            print(f"  {pid}: graphql error {str(d['errors'])[:160]}")
            continue
        rows = d["data"]["biospecimenPerStudy"] or []
        if not rows:
            print(f"  {pid}: 0 rows")
            continue
        aliq = {r["aliquot_submitter_id"] for r in rows if r.get("aliquot_submitter_id")}
        tumor = {r["aliquot_submitter_id"] for r in rows
                 if r.get("sample_type") and "Normal" not in (r["sample_type"] or "")}
        rcases = {r["case_submitter_id"] for r in rows if r.get("case_submitter_id")}
        cov = len(clean_cols & aliq) / max(1, len(clean_cols))
        tcov = len(clean_cols & tumor)
        print(f"  {pid}: {len(rows)} rows, {len(aliq)} aliquots, {len(rcases)} cases "
              f"| protein-col coverage {len(clean_cols & aliq)}/{len(clean_cols)} = {cov:.1%} "
              f"| tumor-only cols {tcov} "
              f"| WSI-case overlap {len(rcases & cases)}/{len(cases)}")
        types = {}
        for r in rows:
            types[r.get("sample_type")] = types.get(r.get("sample_type"), 0) + 1
        print(f"       sample_type tally: {types}")
