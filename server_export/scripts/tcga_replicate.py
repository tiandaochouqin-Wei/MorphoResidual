#!/usr/bin/env python3
"""tcga_replicate.py — multi-cohort external replication (Tier-1 ①). Generalises the
TCGA-KIRC pipeline (tcga_kirc_pull + tcga_kirc_rna_manifest + tcga_kirc_omics) to ANY
TCGA project via --project, so each CPTAC organ gets an independent TCGA cohort on the
RPPA platform:  TCGA-LUAD (lung)  TCGA-UCEC (uterus)  TCGA-PAAD (pancreas)  TCGA-GBM (brain).

Does, in one run (mn02, internet, minutes; slides/RNA download themselves come after):
  1) tri-modal cases (diagnostic slide + STAR RNA + RPPA)
  2) slide manifest for those cases        -> manifest_<short>_slides.txt
  3) tumour STAR-RNA manifest + crosswalk  -> manifest_<short>_rna.txt, rna_case_file.tsv
  4) RPPA total-protein matrix             -> rppa_<short>.csv
  5) coverage vs the matching CPTAC cohort's significant residual set
Run:  python tcga_replicate.py --project TCGA-LUAD
Then (background, big):
  nohup python gdc_download.py <ROOT>/manifest_<short>_slides.txt <ROOT>/slides 4 > sl.log 2>&1 &
  nohup python gdc_download.py <ROOT>/manifest_<short>_rna.txt    <ROOT>/rna    6 > rna.log 2>&1 &
  (gpu02) python extract_phikon.py --raw-dir <ROOT>/slides --out-dir <ROOT>/emb_phikon --stride 512
  python tcga_residual.py --project TCGA-LUAD
"""
import argparse, io, os, re, sys, time
import pandas as pd
try:
    import requests
except Exception as e:
    sys.exit(f"requests unavailable: {e}")

FILES = "https://api.gdc.cancer.gov/files"
DATA = "https://api.gdc.cancer.gov/data/"
ANNOT = "https://api.gdc.cancer.gov/v0/data/62647302-b4d3-4a81-a7c0-d141f5dbd300"
PINNED = "/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/{c}/residual_results_tumoronly.csv"
CPTAC_OF = {"TCGA-KIRC": "ccrcc", "TCGA-LUAD": "luad", "TCGA-UCEC": "ucec", "TCGA-PAAD": "pdac", "TCGA-GBM": "gbm"}
PHOSPHO = re.compile(r"_p[STY]\d")
SESS = requests.Session()


def f_and(*c): return {"op": "and", "content": list(c)}
def f_in(f, v): return {"op": "in", "content": {"field": f, "value": v}}


def query(filters, fields, size=20000):
    for a in range(6):
        try:
            r = SESS.post(FILES, json={"filters": filters, "size": size, "format": "json", "fields": fields}, timeout=180)
            r.raise_for_status(); return r.json()["data"]["hits"]
        except Exception:
            if a == 5: raise
            time.sleep(3 * (a + 1))


def get_text(url, tries=6):
    for a in range(tries):
        try:
            r = SESS.get(url, timeout=180); r.raise_for_status(); return r.text
        except Exception:
            if a == tries - 1: raise
            time.sleep(3 * (a + 1))


def case_set(filters):
    return {c["submitter_id"] for h in query(filters, "cases.submitter_id")
            for c in (h.get("cases") or []) if c.get("submitter_id")}


def write_manifest(hits, path):
    tot = 0
    with io.open(path, "w") as fh:
        fh.write("id\tfilename\tmd5\tsize\tstate\n")
        for h in hits:
            fh.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum','')}\t{h.get('file_size',0)}\treleased\n")
            tot += h.get("file_size", 0) or 0
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True, help="e.g. TCGA-LUAD")
    ap.add_argument("--base", default="/public/home/fjhui/ZW")
    a = ap.parse_args()
    proj = a.project; short = proj.split("-")[1].lower()
    ROOT = f"{a.base}/tcga_{short}"; os.makedirs(ROOT, exist_ok=True)
    PROJ = f_in("cases.project.project_id", [proj])
    print(f"=== {proj} -> {ROOT}  (CPTAC match: {CPTAC_OF.get(proj, '?')}) ===")

    # 1) tri-modal cases
    slides_c = case_set(f_and(PROJ, f_in("data_type", ["Slide Image"]), f_in("experimental_strategy", ["Diagnostic Slide"])))
    rna_c = case_set(f_and(PROJ, f_in("data_type", ["Gene Expression Quantification"]), f_in("analysis.workflow_type", ["STAR - Counts"])))
    rppa_c = case_set(f_and(PROJ, f_in("data_type", ["Protein Expression Quantification"])))
    tri = slides_c & rna_c & rppa_c
    print(f"[cases] slides {len(slides_c)}  rna {len(rna_c)}  rppa {len(rppa_c)}  tri-modal {len(tri)}")
    if len(tri) < 30:
        print("  WARNING: <30 tri-modal cases -- replication would be under-powered")

    # 2) slide manifest (tri-modal cases only)
    sl = query(f_and(PROJ, f_in("data_type", ["Slide Image"]), f_in("experimental_strategy", ["Diagnostic Slide"]),
                     f_in("cases.submitter_id", sorted(tri))), "file_id,file_name,md5sum,file_size,cases.submitter_id")
    tot = write_manifest(sl, f"{ROOT}/manifest_{short}_slides.txt")
    print(f"[slides] {len(sl)} slides, {tot/1e9:.1f} GB -> manifest_{short}_slides.txt")

    # 3) tumour STAR RNA manifest + crosswalk
    hits = query(f_and(PROJ, f_in("data_type", ["Gene Expression Quantification"]), f_in("analysis.workflow_type", ["STAR - Counts"])),
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

    # 4) RPPA total-protein matrix
    ann = pd.read_csv(io.StringIO(get_text(ANNOT)), sep="\t"); ann.columns = [c.strip().lower() for c in ann.columns]
    gcol = next(c for c in ann.columns if "gene" in c)
    acol = next(c for c in ann.columns if c in ("agid", "ag_id"))
    pcol = next((c for c in ann.columns if "peptide" in c or "antibody" in c or "target" in c), None)
    ag2g = {}
    for _, r in ann.iterrows():
        ag = str(r[acol]).strip(); gene = str(r[gcol]).split()[0] if pd.notna(r[gcol]) else ""
        if ag and gene and gene != "nan":
            ag2g[ag] = (gene, not bool(PHOSPHO.search(str(r[pcol]) if pcol else "")))
    rh = query(f_and(PROJ, f_in("data_type", ["Protein Expression Quantification"])),
               "file_id,cases.submitter_id,cases.samples.sample_type", size=5000)
    rows = {}
    for i, h in enumerate(rh, 1):
        cs = h.get("cases", [{}])[0]; case = cs.get("submitter_id")
        stype = (cs.get("samples", [{}])[0] or {}).get("sample_type", "")
        if not case or "Normal" in stype:
            continue
        df = pd.read_csv(io.StringIO(get_text(DATA + h["file_id"])), sep="\t"); df.columns = [c.strip().lower() for c in df.columns]
        per = {}
        for _, rr in df.iterrows():
            ag = str(rr.get("agid", "")).strip()
            if ag in ag2g and ag2g[ag][1]:
                v = rr.get("protein_expression")
                if pd.notna(v):
                    per.setdefault(ag2g[ag][0], []).append(float(v))
        rows[case] = {g: sum(v) / len(v) for g, v in per.items()}
        if i % 100 == 0:
            print(f"  rppa {i}/{len(rh)} ...", flush=True)
    rppa = pd.DataFrame(rows).T; rppa.to_csv(f"{ROOT}/rppa_{short}.csv")
    print(f"[rppa] {rppa.shape[0]} cases x {rppa.shape[1]} total-protein genes -> rppa_{short}.csv")

    # 5) coverage vs matching CPTAC significant set
    cp = CPTAC_OF.get(proj)
    if cp and os.path.exists(PINNED.format(c=cp)):
        md = pd.read_csv(PINNED.format(c=cp))
        sig = set(md[(md.fdr < 0.05) & (md.incremental_r2 > 0)].gene.astype(str))
        ov = sig & set(rppa.columns)
        print(f"[coverage] CPTAC-{cp.upper()} significant {len(sig)}; on RPPA {len(ov)} -> gene-specific replication set")
    print(f"\nNEXT:\n  nohup python gdc_download.py {ROOT}/manifest_{short}_slides.txt {ROOT}/slides 4 > {ROOT}/sl.log 2>&1 &"
          f"\n  nohup python gdc_download.py {ROOT}/manifest_{short}_rna.txt {ROOT}/rna 6 > {ROOT}/rna.log 2>&1 &"
          f"\n  (gpu02) python extract_phikon.py --raw-dir {ROOT}/slides --out-dir {ROOT}/emb_phikon --stride 512"
          f"\n  python tcga_residual.py --project {proj}")


if __name__ == "__main__":
    main()
