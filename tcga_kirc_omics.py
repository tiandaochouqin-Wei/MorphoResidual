#!/usr/bin/env python3
"""tcga_kirc_omics.py — B-11 RPPA side. Downloads the GDC RPPA antibody annotation
(antibody->gene_symbol) + all TCGA-KIRC RPPA files, assembles a per-case total-protein
matrix (cases x genes), and reports how many genes overlap our CCRCC significant
residual set. Small/fast (RPPA files are tiny). Runs on mn02.  python tcga_kirc_omics.py"""
import io, os, re, sys, time
import pandas as pd
try:
    import requests
except Exception as e:
    sys.exit(f"requests unavailable: {e}")

SESS = requests.Session()


def get_text(url, tries=6):
    """GET with retries (robust to GDC connection resets while slides download)."""
    for a in range(tries):
        try:
            r = SESS.get(url, timeout=180)
            r.raise_for_status()
            return r.text
        except Exception as e:
            if a == tries - 1:
                raise
            time.sleep(3 * (a + 1))

OUT = "/public/home/fjhui/ZW/tcga_kirc"
os.makedirs(OUT, exist_ok=True)
FILES = "https://api.gdc.cancer.gov/files"
DATA = "https://api.gdc.cancer.gov/data/"
ANNOT = "https://api.gdc.cancer.gov/v0/data/62647302-b4d3-4a81-a7c0-d141f5dbd300"
PINNED_CCRCC = "/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/ccrcc/residual_results_tumoronly.csv"
PHOSPHO = re.compile(r"_p[STY]\d")


def f_and(*c): return {"op": "and", "content": list(c)}
def f_in(f, v): return {"op": "in", "content": {"field": f, "value": v}}
PROJ = f_in("cases.project.project_id", ["TCGA-KIRC"])

# 1) antibody -> gene annotation ------------------------------------------------
atxt = get_text(ANNOT)
ann = pd.read_csv(io.StringIO(atxt), sep="\t")
ann.columns = [c.strip().lower() for c in ann.columns]
print("[annot] columns:", list(ann.columns))
gcol = next((c for c in ann.columns if "gene" in c), None)
acol = next((c for c in ann.columns if c in ("agid", "ag_id")), None)
pcol = next((c for c in ann.columns if "peptide" in c or "antibody" in c or "target" in c), None)
print(f"[annot] {len(ann)} antibodies; gene col='{gcol}' agid col='{acol}' target col='{pcol}'")
# AGID -> (gene, is_total)
agid2gene = {}
for _, r in ann.iterrows():
    ag = str(r[acol]).strip()
    gene = str(r[gcol]).split()[0] if pd.notna(r[gcol]) else ""
    tgt = str(r[pcol]) if pcol else ""
    is_total = not bool(PHOSPHO.search(tgt))
    if ag and gene and gene != "nan":
        agid2gene[ag] = (gene, is_total)
n_total = sum(1 for v in agid2gene.values() if v[1])
print(f"[annot] {len(agid2gene)} AGID->gene; total-protein antibodies: {n_total}")

# 2) list + download all KIRC RPPA files ---------------------------------------
hits = SESS.post(FILES, json={"filters": f_and(PROJ, f_in("data_type", ["Protein Expression Quantification"])),
                              "size": 2000, "format": "json",
                              "fields": "file_id,cases.submitter_id,cases.samples.sample_type"}, timeout=120).json()["data"]["hits"]
print(f"[rppa] {len(hits)} files")
rows = {}
for i, h in enumerate(hits, 1):
    cs = h.get("cases", [{}])[0]
    case = cs.get("submitter_id")
    stype = (cs.get("samples", [{}])[0] or {}).get("sample_type", "")
    if not case or "Normal" in stype:   # tumour only
        continue
    txt = get_text(DATA + h["file_id"])
    df = pd.read_csv(io.StringIO(txt), sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]
    per = {}
    for _, rr in df.iterrows():
        ag = str(rr.get("agid", "")).strip()
        if ag in agid2gene and agid2gene[ag][1]:      # total protein only
            g = agid2gene[ag][0]
            v = rr.get("protein_expression")
            if pd.notna(v):
                per.setdefault(g, []).append(float(v))
    rows[case] = {g: sum(vs) / len(vs) for g, vs in per.items()}
    if i % 100 == 0:
        print(f"  {i}/{len(hits)} ...", flush=True)
rppa = pd.DataFrame(rows).T   # cases x genes
rppa.to_csv(f"{OUT}/rppa_kirc.csv")
print(f"[rppa] matrix {rppa.shape[0]} cases x {rppa.shape[1]} genes -> rppa_kirc.csv")

# 3) overlap with our CCRCC significant residual set ---------------------------
md = pd.read_csv(PINNED_CCRCC)
sig = set(md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]["gene"].astype(str))
rppa_genes = set(rppa.columns)
ov = sig & rppa_genes
print("\n==================== B-11 RPPA COVERAGE ====================")
print(f" RPPA total-protein genes (KIRC): {len(rppa_genes)}")
print(f" our CCRCC significant residual genes: {len(sig)}")
print(f" overlap (testable gene-specific replication): {len(ov)}")
print(f" -> phenomenon-level test uses all {len(rppa_genes)} RPPA genes;")
print(f"    gene-specific agreement uses the {len(ov)} overlapping genes.")
print(" sample overlap genes:", sorted(list(ov))[:25])
