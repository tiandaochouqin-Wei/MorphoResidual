#!/usr/bin/env python3
"""gdc_inventory_ts.py -- METADATA-ONLY inventory of TCGA-BRCA slide files on the GDC for the
CPTAC-2 BRCA (PDC000173) case set used by the published Diagnostic-Slide (DX) arm.

Reproduces the DX-arm case-set rule of cptac2_replicate.py:131-144 (PDC tumour cases
& GDC 'Diagnostic Slide' cases & GDC 'STAR - Counts' cases) against the CURRENT GDC release,
then lists every 'Tissue Slide' (TS, incl. 'Bottom Slide'/'Top Slide' are all
experimental_strategy 'Tissue Slide') file for those cases.

Reads public metadata only (PDC GraphQL, GDC /files). Downloads nothing. No protein/mRNA/morphology
statistic is computed.

Outputs (this folder):
  inventory_dx_cases.tsv        case  n_dx  ms=1 rna=1
  inventory_ts_files.tsv        one row per TS file (+ DX rows for reference, slide_kind column)
  inventory_ts_case_summary.tsv one row per case
  inventory_summary.json        counts / sizes / sample-type tables
"""
import json, os, re, sys, time, collections
import requests

GDC_FILES = "https://api.gdc.cancer.gov/files"
PDC_URL = "https://pdc.cancer.gov/graphql"
HERE = os.path.dirname(os.path.abspath(__file__))
S = requests.Session()
PROJ = "TCGA-BRCA"
PDC_STUDY = "PDC000173"


def f_and(*c): return {"op": "and", "content": list(c)}
def f_in(f, v): return {"op": "in", "content": {"field": f, "value": v}}


def gdc(filters, fields, size=20000):
    for a in range(6):
        try:
            r = S.post(GDC_FILES, json={"filters": filters, "size": size, "format": "json", "fields": fields}, timeout=180)
            r.raise_for_status()
            return r.json()["data"]["hits"]
        except Exception:
            if a == 5:
                raise
            time.sleep(3 * (a + 1))


def pdc_tumour_cases():
    q = '{ biospecimenPerStudy(pdc_study_id: "%s") { case_submitter_id sample_type } }' % PDC_STUDY
    for a in range(4):
        try:
            r = S.post(PDC_URL, json={"query": q}, timeout=60); r.raise_for_status()
            rows = r.json()["data"]["biospecimenPerStudy"]
            allc = {x["case_submitter_id"] for x in rows}
            tum = {x["case_submitter_id"] for x in rows if "Tumor" in (x.get("sample_type") or "")}
            return tum, allc, collections.Counter(x.get("sample_type") for x in rows)
        except Exception:
            if a == 3:
                raise
            time.sleep(3 * (a + 1))


def case_set(filters):
    return {c["submitter_id"] for h in gdc(filters, "cases.submitter_id")
            for c in (h.get("cases") or []) if c.get("submitter_id")}


def barcode_sample(fn):
    """TCGA-XX-YYYY-01A-01-TS1.<uuid>.svs -> ('01A', '01', 'TS1')"""
    m = re.match(r"TCGA-\w\w-\w{4}-(\d\d[A-Z])-(\d\d)-(\w+)\.", fn)
    return m.groups() if m else (None, None, None)


def main():
    PJ = f_in("cases.project.project_id", [PROJ])
    ms, ms_all, ms_types = pdc_tumour_cases()
    dx_c = case_set(f_and(PJ, f_in("data_type", ["Slide Image"]), f_in("experimental_strategy", ["Diagnostic Slide"])))
    rna_c = case_set(f_and(PJ, f_in("data_type", ["Gene Expression Quantification"]),
                           f_in("analysis.workflow_type", ["STAR - Counts"])))
    tri = sorted(ms & dx_c & rna_c)
    print(f"[cases] pdc tumour {len(ms)} (all {len(ms_all)}), DX {len(dx_c)}, RNA {len(rna_c)}, tri-modal {len(tri)}")

    fields = "file_id,file_name,md5sum,file_size,experimental_strategy,data_type,access,cases.submitter_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.portions.slides.submitter_id"
    rows = []
    for strat in ("Tissue Slide", "Diagnostic Slide"):
        hits = gdc(f_and(PJ, f_in("data_type", ["Slide Image"]), f_in("experimental_strategy", [strat]),
                         f_in("cases.submitter_id", tri)), fields)
        for h in hits:
            cs = (h.get("cases") or [{}])[0]
            samples = cs.get("samples") or []
            st = ";".join(sorted({(s or {}).get("sample_type", "") for s in samples}))
            code, portion, slide = barcode_sample(h["file_name"])
            rows.append(dict(case=cs.get("submitter_id"), strategy=strat, file_id=h["file_id"], file_name=h["file_name"],
                             size=h.get("file_size", 0) or 0, md5=h.get("md5sum", ""), access=h.get("access"),
                             sample_code=code, portion=portion, slide_label=slide, sample_type_gdc=st))
    # TS-only extras: tumour/normal flag from barcode sample code (01-09 tumour, 10-19 normal)
    for r in rows:
        c = r["sample_code"]
        r["tumour"] = (c is not None and 1 <= int(c[:2]) <= 9)
    with open(f"{HERE}/inventory_ts_files.tsv", "w", encoding="utf-8") as fh:
        cols = ["case", "strategy", "file_id", "file_name", "size", "md5", "access", "sample_code", "portion", "slide_label", "sample_type_gdc", "tumour"]
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["case"], x["strategy"], x["file_name"])):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        by[r["case"]][r["strategy"]].append(r)
    summ_rows = []
    for c in tri:
        ts = by[c]["Tissue Slide"]; dx = by[c]["Diagnostic Slide"]
        ts_t = [r for r in ts if r["tumour"]]
        summ_rows.append(dict(case=c, n_dx=len(dx), n_ts_all=len(ts), n_ts_tumour=len(ts_t),
                              n_ts_normal=len(ts) - len(ts_t),
                              ts_tumour_codes=",".join(sorted({r["sample_code"] for r in ts_t})),
                              ts_tumour_labels=",".join(sorted(r["slide_label"] for r in ts_t)),
                              ts_tumour_gb=sum(r["size"] for r in ts_t) / 1e9))
    with open(f"{HERE}/inventory_ts_case_summary.tsv", "w", encoding="utf-8") as fh:
        cols = list(summ_rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in summ_rows:
            fh.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    with open(f"{HERE}/inventory_dx_cases.tsv", "w", encoding="utf-8") as fh:
        fh.write("case\n" + "\n".join(tri) + "\n")

    ts_all = [r for r in rows if r["strategy"] == "Tissue Slide"]
    ts_t = [r for r in ts_all if r["tumour"]]
    dx_all = [r for r in rows if r["strategy"] == "Diagnostic Slide"]
    out = dict(
        query_date=time.strftime("%Y-%m-%d %H:%M:%S"),
        n_pdc_tumour=len(ms), n_pdc_all=len(ms_all), pdc_sample_types=dict(ms_types),
        n_dx_cases_gdc=len(dx_c), n_rna_cases_gdc=len(rna_c), n_trimodal=len(tri),
        ts_files_all=len(ts_all), ts_gb_all=sum(r["size"] for r in ts_all) / 1e9,
        ts_files_tumour=len(ts_t), ts_gb_tumour=sum(r["size"] for r in ts_t) / 1e9,
        ts_sample_codes=dict(collections.Counter(r["sample_code"] for r in ts_all)),
        ts_slide_labels=dict(collections.Counter(re.sub(r"\d+$", "", r["slide_label"] or "") + "#" + (r["slide_label"] or "")[-1:] for r in ts_all)),
        ts_cases_any=sum(1 for s in summ_rows if s["n_ts_all"] > 0),
        ts_cases_tumour=sum(1 for s in summ_rows if s["n_ts_tumour"] > 0),
        ts_cases_no_tumour=[s["case"] for s in summ_rows if s["n_ts_tumour"] == 0],
        ts_tumour_per_case=dict(collections.Counter(s["n_ts_tumour"] for s in summ_rows)),
        dx_files=len(dx_all), dx_gb=sum(r["size"] for r in dx_all) / 1e9,
        dx_per_case=dict(collections.Counter(s["n_dx"] for s in summ_rows)),
        ts_access=dict(collections.Counter(r["access"] for r in ts_all)),
        ts_ext=dict(collections.Counter(os.path.splitext(r["file_name"])[1] for r in ts_all)),
    )
    json.dump(out, open(f"{HERE}/inventory_summary.json", "w", encoding="utf-8"), indent=1, default=str, ensure_ascii=False)
    print(json.dumps(out, indent=1, default=str, ensure_ascii=False))


if __name__ == "__main__":
    main()
