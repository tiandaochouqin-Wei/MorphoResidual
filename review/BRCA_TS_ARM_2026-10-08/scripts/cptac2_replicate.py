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

PATCH 2026-10-08 (review/BRCA_TS_ARM_2026-10-08): --slide-type {dx,ts}.
  dx (DEFAULT) = the code path above, unchanged: the published Diagnostic-Slide arm.
  ts           = GDC 'Tissue Slide' arm (frozen TS/BS/MS sections cut either side of the analyte
                 portion). Uses the SAME case set as the DX arm (PDC tumour & GDC DX & GDC STAR),
                 asserts it equals the cases already in <ROOT>/rna_case_file.tsv, and writes ONLY
                 slide files (manifest_<short>_ts_slides.txt, md5_<short>_ts.txt, ts_slide_selection.tsv).
                 It never touches the RNA manifest, rna_case_file.tsv or protein_<short>.csv.
"""
import argparse, collections, io, os, re, sys, time
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


def slide_parts(fn):
    """TCGA-XX-YYYY-01A-01-TSA.<uuid>.svs -> ('01A', '01', 'TSA'); (None, None, None) if not a slide barcode"""
    m = re.match(r"TCGA-\w\w-\w{4}-(\d\d[A-Z])-(\d\d)-(\w+?)\.", fn)
    return m.groups() if m else (None, None, None)


def pdc_ms_vials(pdc_study_id):
    """case -> set of sample vials (e.g. '01A') that carry a PDC tumour aliquot."""
    q = '{ biospecimenPerStudy(pdc_study_id: "%s") { case_submitter_id sample_submitter_id sample_type } }' % pdc_study_id
    rows = pdc_query(q)["data"]["biospecimenPerStudy"]
    out = collections.defaultdict(set)
    for x in rows:
        if "Tumor" in (x.get("sample_type") or ""):
            out[x["case_submitter_id"]].add(x["sample_submitter_id"].split("-")[-1])
    return out


def build_ts_arm(a, proj, short, ROOT, PROJ, tri):
    """Tissue-Slide arm: manifest only. Same cases as the DX arm; omics are NOT re-pulled."""
    xw_path = f"{ROOT}/rna_case_file.tsv"
    if os.path.exists(xw_path):
        prev = set(pd.read_csv(xw_path, sep="\t")["case"].map(tcga_case))
        if prev != set(tri):
            print(f"  CASE-SET DRIFT vs {xw_path}: only_now={sorted(set(tri) - prev)} only_prev={sorted(prev - set(tri))}")
            if not a.allow_case_drift:
                sys.exit("abort: re-derived DX case set != published rna_case_file.tsv (--allow-case-drift to override; "
                         "the TS arm would then NOT be the same cases as the DX arm)")
    else:
        print(f"  WARNING: {xw_path} missing -> cannot assert that the case set equals the DX arm's")
    hits = gdc_query(f_and(PROJ, f_in("data_type", ["Slide Image"]), f_in("experimental_strategy", ["Tissue Slide"]),
                           f_in("cases.submitter_id", sorted(tri))),
                     "file_id,file_name,md5sum,file_size,cases.submitter_id,cases.samples.sample_type")
    vials = pdc_ms_vials(a.pdc_study) if a.ts_rule == "ms_vial" else {}
    keep, rows = [], []
    for h in hits:
        code, portion, label = slide_parts(h["file_name"])
        case = tcga_case(h["file_name"])
        prim = code is not None and code.startswith("01")        # 01 = Primary Solid Tumor; drops 06 metastatic, 11 normal
        ok = prim
        if ok and a.ts_rule == "primary_ts_only":
            ok = label.startswith("TS")
        if ok and a.ts_rule == "ms_vial":
            ok = code in vials.get(case, set())
        why = "keep" if ok else (f"non-primary sample code {code}" if not prim else f"rule {a.ts_rule}")
        rows.append((case, h["file_name"], code, portion, label, h.get("file_size", 0), "1" if ok else "0", why))
        if ok:
            keep.append(h)
    with io.open(f"{ROOT}/ts_slide_selection.tsv", "w") as fh:
        fh.write("case\tfilename\tsample_code\tportion\tslide_label\tsize\tkeep\treason\n")
        for r in sorted(rows):
            fh.write("\t".join(str(v) for v in r) + "\n")
    tot = write_manifest(keep, f"{ROOT}/manifest_{short}_ts_slides.txt")
    with io.open(f"{ROOT}/md5_{short}_ts.txt", "w") as fh:         # md5sum -c list, run inside slides_ts/
        for h in keep:
            fh.write(f"{h.get('md5sum', '')}  {h['file_name']}\n")
    cov = {tcga_case(h["file_name"]) for h in keep}
    print(f"[ts] rule={a.ts_rule}: {len(keep)}/{len(hits)} TS-arm files kept for {len(cov)}/{len(tri)} cases, "
          f"{tot/1e9:.1f} GB -> manifest_{short}_ts_slides.txt ; selection table ts_slide_selection.tsv")
    miss = sorted(set(tri) - cov)
    if miss:
        print(f"  cases with no kept TS-arm slide ({len(miss)}): {miss}")
    print(f"\nNEXT (ts arm):"
          f"\n  nohup python gdc_download.py {ROOT}/manifest_{short}_ts_slides.txt {ROOT}/slides_ts 4 > {ROOT}/sl_ts.log 2>&1 &"
          f"\n  (cd {ROOT}/slides_ts && md5sum -c {ROOT}/md5_{short}_ts.txt | grep -v ': OK'; echo md5-check-done)"
          f"\n  extract_phikon.py --raw-dir {ROOT}/slides_ts --out-dir {ROOT}/emb_phikon_ts --stride 0 --manifest {ROOT}/manifest_{short}_ts_slides.txt"
          f"\n  python cptac2_residual.py --project {proj} --slide-type ts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True, help="e.g. TCGA-BRCA")
    ap.add_argument("--pdc-study", required=True, help="e.g. PDC000173")
    ap.add_argument("--data-type", default="log2_ratio", help="PDC quantDataMatrix data_type")
    ap.add_argument("--base", default="/public/home/fjhui/ZW")
    ap.add_argument("--skip-protein", action="store_true", help="skip the quant-matrix pull (e.g. COAD data_type unresolved)")
    ap.add_argument("--slide-type", choices=["dx", "ts"], default="dx",
                    help="dx = Diagnostic Slide (published arm; default, unchanged); ts = Tissue Slide arm")
    ap.add_argument("--ts-rule", choices=["primary_all", "primary_ts_only", "ms_vial"], default="primary_all",
                    help="ts only. primary_all: every TS/BS/MS slide whose sample code is 01x (primary tumour); "
                         "primary_ts_only: primary_all restricted to top slides (TS*); "
                         "ms_vial: primary_all restricted to the vial (01A/01B) of the PDC MS aliquot")
    ap.add_argument("--allow-case-drift", action="store_true",
                    help="ts only: do not abort when the re-derived case set differs from rna_case_file.tsv")
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

    if a.slide_type == "ts":
        build_ts_arm(a, proj, short, ROOT, PROJ, tri)
        return

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
