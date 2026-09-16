#!/usr/bin/env python3
"""
C1-UCEC confirmatory cohort: pull PDC000439 (TMT11) protein + GDC STAR-Counts RNA for
the 138 confirmatory cases, per the frozen rule (review/C1_UCEC_FROZEN_RULE_2026-09-15.md
§1). Protein follows cptac2_replicate.py's GraphQL quantDataMatrix method (frozen rule's
explicit instruction), not the discovery cohort's manual-portal-TSV method (that one exists
only because PDC's signed download URLs expire in 7 days -- quantDataMatrix has no such
expiry and is the simpler path when it works).

Unlike cptac2_replicate.py, the case set is NOT rediscovered here: it is the 138 cases
already fixed by the confirmatory slide download (pinned/c1_case_ids_ucec.txt). This script
only pulls the other two modalities for that fixed set -- it must not add or drop cases.

quantDataMatrix's column-header aliquot format is NOT assumed to match cptac2_replicate.py's
TCGA-barcode parsing (that parsing is specific to the TCGA-legacy CPTAC-2 cohorts it was
written for; PDC000439 is a native CPTAC-3 study using CPTAC aliquot_submitter_ids, e.g.
CPT0123456789, the same convention as the discovery protein TSV's "<aliquot> Log Ratio"
columns). The header is printed raw and joined against the aliquot->case crosswalk built
here; if fewer than half the crosswalk's tumour aliquots match, the script stops rather than
silently emitting a near-empty matrix.

Usage:
  python c1_pull_omics.py protein   [--pdc-study PDC000439] [--data-type log2_ratio]
  python c1_pull_omics.py rna
  python c1_pull_omics.py both
Outputs (relative to --out-dir, default /public/home/fjhui/ZW/ucec_c1):
  pinned_aliquot_to_case_tumor_ucec_c1.tsv   aliquot->case crosswalk (tumour only)
  protein_ucec_c1.csv                        case-indexed protein matrix
  manifest_rna_tumor_ucec_c1.tsv             RNA manifest in residual_analysis.py's format
  gdc_rna_download_manifest_ucec_c1.tsv      gdc_download.py-compatible download manifest
"""
import argparse
import io
import json
import os
import sys
import time
import urllib.request

import pandas as pd

PDC_URL = "https://pdc.cancer.gov/graphql"
GDC_FILES = "https://api.gdc.cancer.gov/files"


def pdc_query(q, tries=6):
    req = urllib.request.Request(
        PDC_URL, data=json.dumps({"query": q}).encode(),
        headers={"Content-Type": "application/json"})
    for a in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                d = json.load(r)
            if d.get("errors"):
                sys.exit(f"PDC graphql error: {str(d['errors'])[:300]}")
            return d["data"]
        except urllib.error.URLError as e:
            if a == tries - 1:
                raise
            print(f"  PDC attempt {a + 1}/{tries} failed: {e}", flush=True)
            time.sleep(4 * (a + 1))


def gdc_post(url, payload, tries=6):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
    for a in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except urllib.error.URLError as e:
            if a == tries - 1:
                raise
            print(f"  GDC attempt {a + 1}/{tries} failed: {e}", flush=True)
            time.sleep(4 * (a + 1))


def load_case_ids(path):
    with open(path, encoding="utf-8") as f:
        return sorted({line.strip() for line in f if line.strip()})


def build_aliquot_xwalk(pdc_study_id, case_ids, out_path):
    q = ('{ biospecimenPerStudy(pdc_study_id: "%s" acceptDUA: true) '
         '{ aliquot_submitter_id case_submitter_id sample_type } }' % pdc_study_id)
    rows = pdc_query(q)["biospecimenPerStudy"] or []
    print(f"[pdc] {pdc_study_id}: {len(rows)} biospecimen rows", flush=True)
    tally = {}
    for r in rows:
        tally[r.get("sample_type")] = tally.get(r.get("sample_type"), 0) + 1
    print(f"  sample_type tally: {tally}")

    seen = {}
    for r in rows:
        a = r.get("aliquot_submitter_id")
        st = (r.get("sample_type") or "").strip()
        if not a or "Normal" in st or a in seen:
            continue
        seen[a] = (r.get("case_submitter_id"), st)

    case_set = set(case_ids)
    covered = {c for _a, (c, _s) in seen.items() if c in case_set}
    print(f"[xwalk] {len(seen)} tumour aliquots, covering {len(covered)}/{len(case_ids)} "
          f"of the 138 confirmatory cases")
    missing = case_set - covered
    if missing:
        print(f"  {len(missing)} confirmatory cases have NO tumour aliquot in {pdc_study_id}: "
              f"{sorted(missing)[:10]}{' ...' if len(missing) > 10 else ''}")

    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("aliquot_id\tcase_id\tsample_type\n")
        for a, (c, st) in sorted(seen.items()):
            f.write(f"{a}\t{c}\t{st}\n")
    print(f"[xwalk] wrote {len(seen)} rows -> {out_path}")
    return seen


def pull_protein(args):
    case_ids = load_case_ids(args.case_ids)
    print(f"{len(case_ids)} confirmatory case ids loaded from {args.case_ids}")

    xwalk_path = os.path.join(args.out_dir, "pinned_aliquot_to_case_tumor_ucec_c1.tsv")
    xwalk = build_aliquot_xwalk(args.pdc_study, case_ids, xwalk_path)

    q = '{ quantDataMatrix(pdc_study_id: "%s", data_type: "%s") }' % (
        args.pdc_study, args.data_type)
    mat = pdc_query(q)["quantDataMatrix"]
    header = mat[0]
    print(f"[pdc] quantDataMatrix {args.pdc_study}/{args.data_type}: "
          f"{len(mat) - 1} genes x {len(header) - 1} columns")
    print(f"  raw header sample (first 5): {header[1:6]}")

    # observed header format: "<aliquot_uuid>:<aliquot_submitter_id>", e.g.
    # "e836eba1-f839-49f1-a708-1a103de8afcb:CPT0127610003" (verified 2026-09-16 against
    # PDC000439's actual quantDataMatrix response -- take the part after the colon).
    matched = sum(1 for c in header[1:] if c.split(":")[-1] in xwalk)
    n_tumour_aliquots = len(xwalk)
    print(f"  {matched}/{len(header) - 1} matrix columns match a tumour aliquot in the "
          f"crosswalk ({n_tumour_aliquots} tumour aliquots total)")
    if matched < n_tumour_aliquots / 2:
        sys.exit(
            "FATAL: fewer than half the crosswalk's tumour aliquots matched a matrix column "
            "-- the header format assumption (raw aliquot_submitter_id) is probably wrong for "
            "this study. Inspect the raw header above before adapting the parsing; do not "
            "guess a TCGA-barcode transform (that is for the CPTAC-2 legacy cohorts only).")

    cols, keep_idx = [], []
    for i, c in enumerate(header[1:]):
        case = xwalk.get(c.split(":")[-1], (None, None))[0]
        if case:
            cols.append(case)
            keep_idx.append(i)

    rows = {}
    for r in mat[1:]:
        gene = r[0]
        vals = r[1:]
        rows[gene] = [vals[i] for i in keep_idx]
    df = pd.DataFrame(rows, index=cols).apply(pd.to_numeric, errors="coerce")
    df = df.groupby(level=0).mean()  # >1 tumour aliquot per case -> average, matches discovery

    case_set = set(case_ids)
    overlap = len(set(df.index) & case_set)
    print(f"[protein] {df.shape[0]} cases x {df.shape[1]} genes; "
          f"{overlap}/{len(case_ids)} of the confirmatory case set present")

    out_path = os.path.join(args.out_dir, "protein_ucec_c1.csv")
    df.to_csv(out_path)
    print(f"[protein] wrote -> {out_path}")


def pull_rna(args):
    case_ids = load_case_ids(args.case_ids)
    print(f"{len(case_ids)} confirmatory case ids loaded from {args.case_ids}")

    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.submitter_id", "value": case_ids}},
            {"op": "in", "content": {"field": "data_type",
                                      "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type",
                                      "value": ["STAR - Counts"]}},
        ],
    }
    fields = "file_id,file_name,md5sum,file_size,cases.submitter_id,cases.samples.sample_type"
    payload = {"filters": filters, "size": 5000, "format": "json", "fields": fields}
    hits = gdc_post(GDC_FILES, payload)["data"]["hits"]
    print(f"[gdc] {len(hits)} STAR-Counts files across the 138 confirmatory cases "
          f"(before tumour/normal filtering)")

    kept, dropped_normal, dropped_nocase = [], 0, 0
    case_manifest_path = os.path.join(args.out_dir, "manifest_rna_tumor_ucec_c1.tsv")
    dl_manifest_path = os.path.join(args.out_dir, "gdc_rna_download_manifest_ucec_c1.tsv")
    with open(case_manifest_path, "w", encoding="utf-8", newline="\n") as fm, \
         open(dl_manifest_path, "w", encoding="utf-8", newline="\n") as fd:
        fm.write("file_id\tfile_name\tcase_submitter_id\tsample_type\n")
        fd.write("id\tfilename\tmd5\tsize\tstate\n")
        for h in hits:
            cs = (h.get("cases") or [{}])[0]
            case = cs.get("submitter_id")
            stype = (cs.get("samples", [{}])[0] or {}).get("sample_type", "")
            if not case:
                dropped_nocase += 1
                continue
            if "Normal" in stype:
                dropped_normal += 1
                continue
            fm.write(f"{h['file_id']}\t{h['file_name']}\t{case}\t{stype}\n")
            fd.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum', '')}\t"
                     f"{h.get('file_size', 0)}\treleased\n")
            kept.append(case)

    covered = len(set(kept) & set(case_ids))
    print(f"[rna] kept {len(kept)} tumour files (dropped {dropped_normal} normal, "
          f"{dropped_nocase} no-case), covering {covered}/{len(case_ids)} confirmatory cases")
    dup = len(kept) - len(set(kept))
    if dup:
        print(f"  note: {dup} cases have >1 tumour STAR file (multiple RNA aliquots); "
              f"load_rna_matrix() averages these per gene, matching discovery")
    missing = set(case_ids) - set(kept)
    if missing:
        print(f"  {len(missing)} confirmatory cases have no tumour STAR-Counts file: "
              f"{sorted(missing)[:10]}{' ...' if len(missing) > 10 else ''}")
    print(f"[rna] wrote {case_manifest_path}\n[rna] wrote {dl_manifest_path}")
    print(f"\nNEXT: nohup python gdc_download.py {dl_manifest_path} "
          f"{args.out_dir}/omics/rna 6 > {args.out_dir}/rna_dl.log 2>&1 &")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["protein", "rna", "both"])
    ap.add_argument("--pdc-study", default="PDC000439")
    ap.add_argument("--data-type", default="log2_ratio")
    ap.add_argument("--case-ids", default="pinned/c1_case_ids_ucec.txt")
    ap.add_argument("--out-dir", default="/public/home/fjhui/ZW/ucec_c1")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    if args.mode in ("protein", "both"):
        pull_protein(args)
    if args.mode in ("rna", "both"):
        pull_rna(args)


if __name__ == "__main__":
    main()
