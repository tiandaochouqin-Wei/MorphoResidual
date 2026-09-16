#!/usr/bin/env python3
"""cptac2_replicate.py — external validation on the ORIGINAL CPTAC-2 (2014/2016
Nature/Cell) mass-spec cohorts, which reprofiled TCGA-legacy cases: TCGA Colon
Cancer Proteome (Zhang 2014), TCGA Breast Cancer Proteome (Mertins 2016), TCGA
Ovarian JHU Proteome (Zhang 2016). Unlike the RPPA replication (tcga_replicate.py),
these are full mass-spec panels (~8.6k-10.6k proteins) that actually cover the
translation / ER-secretion / splicing / folding machinery the main CPTAC-3 analysis
converged on -- so this tests pathway-level convergence, not just "is morphology
predictable". No CPTAC-3 discovery cohort exists for these three organs, so there is
no gene-specific replication-rate step (see cptac2_residual.py / cptac2_enrichment.py
instead).

Studies (PDC, "CPTAC2 Retrospective" program):
  COAD  PDC000111  Label Free   90 MS cases  (data_type TBD -- quantDataMatrix
                                               rejected every guessed data_type;
                                               use --skip-protein and fetch the
                                               matrix by hand if pursuing COAD)
  BRCA  PDC000173  iTRAQ4      109 MS cases  data_type=log2_ratio  (RECOMMENDED:
                                               102 tri-modal, 10625 proteins, full
                                               coverage of all 4 signature families)
  OV    PDC000113  iTRAQ4      124 MS cases  data_type=log2_ratio  (16 tri-modal --
                                               too small, do not run)

Does, in one run (mn02, internet, minutes; slides/RNA download come after):
  1) PDC tumour case list for the study
  2) tri-modal cases (diagnostic slide + STAR RNA + PDC MS case)
  3) slide manifest for those cases        -> manifest_<short>_slides.txt
  4) tumour STAR-RNA manifest + crosswalk  -> manifest_<short>_rna.txt, rna_case_file.tsv
  5) PDC quant matrix, case-indexed        -> protein_<short>.csv
Run:  python cptac2_replicate.py --project TCGA-BRCA --pdc-study PDC000173 --data-type log2_ratio
Then (background, big):
  nohup python gdc_download.py <ROOT>/manifest_<short>_slides.txt <ROOT>/slides 4 > <ROOT>/sl.log 2>&1 &
  nohup python gdc_download.py <ROOT>/manifest_<short>_rna.txt    <ROOT>/rna    6 > <ROOT>/rna.log 2>&1 &
  (gpu02) python extract_phikon.py --raw-dir <ROOT>/slides --out-dir <ROOT>/emb_phikon --stride 512
  python cptac2_residual.py --project TCGA-BRCA
  python cptac2_enrichment.py TCGA-BRCA
"""
import argparse, io, os, sys, time
import pandas as pd
try:
    import requests
except Exception as e:
    sys.exit(f"requests unavailable: {e}")

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data/"
PDC_URL = "https://pdc.cancer.gov/graphql"
SESS = requests.Session()


def f_and(*c): return {"op": "and", "content": list(c)}
def f_in(f, v): return {"op": "in", "content": {"field": f, "value": v}}


def gdc_query(filters, fields, size=20000):
    for a in range(6):
        try:
            r = SESS.post(GDC_FILES, json={"filters": filters, "size": size, "format": "json", "fields": fields}, timeout=180)
            r.raise_for_status(); return r.json()["data"]["hits"]
        except Exception:
            if a == 5: raise
            time.sleep(3 * (a + 1))


def gdc_get_text(url, tries=6):
    for a in range(tries):
        try:
            r = SESS.get(url, timeout=180); r.raise_for_status(); return r.text
        except Exception:
            if a == tries - 1: raise
            time.sleep(3 * (a + 1))


def gdc_case_set(filters):
    return {c["submitter_id"] for h in gdc_query(filters, "cases.submitter_id")
            for c in (h.get("cases") or []) if c.get("submitter_id")}


def write_manifest(hits, path):
    tot = 0
    with io.open(path, "w") as fh:
        fh.write("id\tfilename\tmd5\tsize\tstate\n")
        for h in hits:
            fh.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum','')}\t{h.get('file_size',0)}\treleased\n")
            tot += h.get("file_size", 0) or 0
    return tot


def pdc_query(q, variables=None, tries=4):
    data = __import__("json").dumps({"query": q, "variables": variables or {}}).encode()
    req = __import__("urllib.request", fromlist=["Request"]).Request(
        PDC_URL, data=data, headers={"Content-Type": "application/json"})
    import urllib.request, json as _json
    for a in range(tries):
        try:
            return _json.load(urllib.request.urlopen(req, timeout=60))
        except Exception:
            if a == tries - 1: raise
            time.sleep(3 * (a + 1))


def pdc_tumour_cases(pdc_study_id):
    """Case-level TCGA barcodes with a Tumor sample in this PDC study."""
    q = '{ biospecimenPerStudy(pdc_study_id: "%s") { case_submitter_id sample_type } }' % pdc_study_id
    r = pdc_query(q)
    rows = r["data"]["biospecimenPerStudy"]
    tumour = {x["case_submitter_id"] for x in rows if "Tumor" in (x.get("sample_type") or "")}
    return tumour


def pdc_quant_matrix(pdc_study_id, data_type):
    q = '{ quantDataMatrix(pdc_study_id: "%s", data_type: "%s") }' % (pdc_study_id, data_type)
    r = pdc_query(q)
    if "errors" in r:
        sys.exit(f"PDC quantDataMatrix failed for {pdc_study_id}/{data_type}: {r['errors'][0]['message']}")
    return r["data"]["quantDataMatrix"]


def tcga_case(s):
    return "-".join(str(s).split("-")[:3])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True, help="e.g. TCGA-BRCA")
    ap.add_argument("--pdc-study", required=True, help="e.g. PDC000173")
    ap.add_argument("--data-type", default="log2_ratio", help="PDC quantDataMatrix data_type")
    ap.add_argument("--base", default="/public/home/fjhui/ZW")
    ap.add_argument("--skip-protein", action="store_true", help="skip the quant-matrix pull (e.g. COAD data_type unresolved)")
    a = ap.parse_args()
    proj = a.project; short = proj.split("-")[1].lower()
    ROOT = f"{a.base}/cptac2_{short}"; os.makedirs(ROOT, exist_ok=True)
    PROJ = f_in("cases.project.project_id", [proj])
    print(f"=== {proj} -> {ROOT}  (PDC study {a.pdc_study}, data_type={a.data_type}) ===")

    # 1) PDC tumour case list
    ms_cases = pdc_tumour_cases(a.pdc_study)
    print(f"[pdc] {len(ms_cases)} tumour cases in {a.pdc_study}")

    # 2) tri-modal
    slides_c = gdc_case_set(f_and(PROJ, f_in("data_type", ["Slide Image"]), f_in("experimental_strategy", ["Diagnostic Slide"])))
    rna_c = gdc_case_set(f_and(PROJ, f_in("data_type", ["Gene Expression Quantification"]), f_in("analysis.workflow_type", ["STAR - Counts"])))
    tri = ms_cases & slides_c & rna_c
    print(f"[cases] pdc_ms {len(ms_cases)}  slides {len(slides_c)}  rna {len(rna_c)}  tri-modal {len(tri)}")
    if len(tri) < 30:
        print("  WARNING: <30 tri-modal cases -- do not proceed, replication would be under-powered "
              "(this is what happened to TCGA-GBM/TCGA-PAAD in the RPPA line: n=85/112 already overfit "
              "the fixed 20-PC pipeline; do not go smaller)")

    # 3) slide manifest (tri-modal cases only)
    sl = gdc_query(f_and(PROJ, f_in("data_type", ["Slide Image"]), f_in("experimental_strategy", ["Diagnostic Slide"]),
                          f_in("cases.submitter_id", sorted(tri))), "file_id,file_name,md5sum,file_size,cases.submitter_id")
    tot = write_manifest(sl, f"{ROOT}/manifest_{short}_slides.txt")
    print(f"[slides] {len(sl)} slides, {tot/1e9:.1f} GB -> manifest_{short}_slides.txt")

    # 4) tumour STAR RNA manifest + crosswalk (tri-modal cases only)
    hits = gdc_query(f_and(PROJ, f_in("data_type", ["Gene Expression Quantification"]), f_in("analysis.workflow_type", ["STAR - Counts"]),
                            f_in("cases.submitter_id", sorted(tri))),
                      "file_id,file_name,md5sum,file_size,cases.submitter_id,cases.samples.sample_type")
    keep = []
    with io.open(f"{ROOT}/rna_case_file.tsv", "w") as fx:
        fx.write("case\tsample_type\tfilename\n")
        for h in hits:
            cs = (h.get("cases") or [{}])[0]; case = cs.get("submitter_id")
            stype = (cs.get("samples", [{}])[0] or {}).get("sample_type", "")
            if not case or "Normal" in stype:
                continue
            keep.append(h); fx.write(f"{case}\t{stype}\t{h['file_name']}\n")
    tot = write_manifest(keep, f"{ROOT}/manifest_{short}_rna.txt")
    print(f"[rna] {len(keep)} tumour STAR files, {tot/1e9:.2f} GB -> manifest_{short}_rna.txt + rna_case_file.tsv")

    # 5) PDC quant matrix -> case-indexed protein_<short>.csv
    if not a.skip_protein:
        mat = pdc_quant_matrix(a.pdc_study, a.data_type)
        header = mat[0]
        # column header format observed: "<aliquot_uuid>:<site>-<participant>-<sample>"
        # (no "TCGA-" prefix) -> prepend it, then truncate to case level.
        cols = []
        for c in header[1:]:
            suffix = c.split(":")[-1]
            cols.append(tcga_case("TCGA-" + suffix))
        rows = {}
        for r in mat[1:]:
            gene = r[0]
            vals = r[1:]
            rows[gene] = vals
        df = pd.DataFrame(rows, index=cols).apply(pd.to_numeric, errors="coerce")
        # multiple aliquots can map to the same case (dedupe by mean)
        df = df.groupby(level=0).mean()
        df.to_csv(f"{ROOT}/protein_{short}.csv")
        print(f"[protein] {df.shape[0]} cases x {df.shape[1]} proteins ({a.data_type}) -> protein_{short}.csv")
    else:
        print("[protein] skipped (--skip-protein)")

    print(f"\nNEXT:\n  nohup python gdc_download.py {ROOT}/manifest_{short}_slides.txt {ROOT}/slides 4 > {ROOT}/sl.log 2>&1 &"
          f"\n  nohup python gdc_download.py {ROOT}/manifest_{short}_rna.txt {ROOT}/rna 6 > {ROOT}/rna.log 2>&1 &"
          f"\n  (gpu02) python extract_phikon.py --raw-dir {ROOT}/slides --out-dir {ROOT}/emb_phikon --stride 512"
          f"\n  python cptac2_residual.py --project {proj}"
          f"\n  python cptac2_enrichment.py {proj}")


if __name__ == "__main__":
    main()
