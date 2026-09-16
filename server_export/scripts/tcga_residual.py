#!/usr/bin/env python3
"""tcga_residual.py — external replication for ANY TCGA project (generalises
kirc_residual.py). Tests morphology's incremental R^2 over an mRNA baseline for the
RPPA-protein residual with the SAME estimand as CPTAC (5-fold OOF, closed-form ridge
lambda=1, 1000 patient<->slide permutations, BH-FDR), and the gene-specific replication
rate against the matching CPTAC cohort's significant set.

Inputs under /public/home/fjhui/ZW/tcga_<short>/ (from tcga_replicate.py + downloads):
  rppa_<short>.csv   rna_case_file.tsv   rna/   emb_phikon/*.pt (or MORPHO_WSI_EMB_DIR)
Run:  python -u tcga_residual.py --project TCGA-LUAD
Outputs: results/<short>_rppa_results.csv + results/replication_summary_<short>.csv
"""
import argparse, os, glob, sys
import numpy as np, pandas as pd, torch
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

CV_FOLDS, N_PERM, RIDGE_ALPHA = 5, 1000, 1.0
PINNED = "/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/{c}/residual_results_tumoronly.csv"
CPTAC_OF = {"TCGA-KIRC": "ccrcc", "TCGA-LUAD": "luad", "TCGA-UCEC": "ucec", "TCGA-PAAD": "pdac", "TCGA-GBM": "gbm"}


def tcga_case(s): return "-".join(str(s).split("-")[:3])


def cv_r2(X, y, folds=CV_FOLDS):
    kf = KFold(n_splits=folds, shuffle=True, random_state=0)
    preds = np.zeros_like(y, dtype=float); I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        Xtr, Xte, ytr = X[tr], X[te], y[tr]
        xm = Xtr.mean(0); xs = Xtr.std(0); xs[xs == 0] = 1.0; ym = ytr.mean()
        A = (Xtr - xm) / xs
        w = np.linalg.solve(A.T @ A + I, A.T @ (ytr - ym))
        preds[te] = ((Xte - xm) / xs) @ w + ym
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - np.sum((y - preds) ** 2) / ss_tot if ss_tot > 0 else 0.0


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--base", default="/public/home/fjhui/ZW")
    a = ap.parse_args()
    proj = a.project; short = proj.split("-")[1].lower(); cp = CPTAC_OF.get(proj)
    ROOT = f"{a.base}/tcga_{short}"
    EMB = os.environ.get("MORPHO_WSI_EMB_DIR", f"{ROOT}/emb_phikon")
    OUT = os.environ.get("MORPHO_OUT", f"{ROOT}/results"); os.makedirs(OUT, exist_ok=True)

    rppa = pd.read_csv(f"{ROOT}/rppa_{short}.csv", index_col=0)
    rppa.index = [tcga_case(c) for c in rppa.index]; rppa = rppa[~rppa.index.duplicated()]

    xw = pd.read_csv(f"{ROOT}/rna_case_file.tsv", sep="\t"); rows = {}
    for _, r in xw.iterrows():
        fp = f"{ROOT}/rna/{r['filename']}"
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp, sep="\t", comment="#"); d.columns = [c.strip() for c in d.columns]
        d = d[["gene_name", "tpm_unstranded"]].dropna()
        d = d[~d["gene_name"].astype(str).str.startswith("N_")]
        rows.setdefault(tcga_case(r["case"]), []).append(d.groupby("gene_name")["tpm_unstranded"].mean())
    rna = np.log2(pd.DataFrame({c: pd.concat(v, axis=1).mean(axis=1) for c, v in rows.items()}).T.fillna(0) + 1)

    emb = {}
    for pt in glob.glob(f"{EMB}/*.pt"):
        e = torch.load(pt, map_location="cpu")["embeddings"].float().mean(0).numpy()
        emb.setdefault(tcga_case(os.path.basename(pt)), []).append(e)
    wsi = pd.DataFrame({c: np.mean(v, axis=0) for c, v in emb.items()}).T
    if wsi.empty:
        sys.exit(f"no embeddings in {EMB}")

    common = sorted(set(rppa.index) & set(rna.index) & set(wsi.index))
    print(f"[{proj}] cases: rppa {len(rppa)}, rna {len(rna)}, wsi {len(wsi)}, common {len(common)}")
    rppa, rna, wsi = rppa.loc[common], rna.loc[common], wsi.loc[common]
    pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
    perm = [np.random.RandomState(i).permutation(len(common)) for i in range(N_PERM)]
    genes = [g for g in rppa.columns if g in rna.columns]
    print(f"[{proj}] testing {len(genes)} RPPA genes present in RNA")

    res = []
    for g in genes:
        y = np.asarray(rppa[g].values, float); v = ~np.isnan(y)
        if v.sum() < 30:
            continue
        yv = y[v]; xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1); Wv = pcs[v]
        r2r = cv_r2(xr, yv); r2b = cv_r2(np.hstack([xr, Wv]), yv); incr = r2b - r2r
        null = np.array([cv_r2(np.hstack([xr, pcs[p][v]]), yv) - r2r for p in perm])
        res.append({"gene": g, "n": int(v.sum()), "r2_rna": r2r, "r2_both": r2b,
                    "incremental_r2": incr, "pval": (np.sum(null >= incr) + 1) / (len(perm) + 1)})
    res = pd.DataFrame(res); res["fdr"] = bh_fdr(res["pval"].values)
    res.to_csv(f"{OUT}/{short}_rppa_results.csv", index=False)
    sig = res[(res.fdr < 0.05) & (res.incremental_r2 > 0)]

    summ = {"project": proj, "cptac_match": cp, "n_cases": len(common), "tested": len(res),
            "significant": len(sig), "pct_sig": 100 * len(sig) / max(len(res), 1),
            "median_incr_sig": float(sig.incremental_r2.median()) if len(sig) else np.nan}
    if cp and os.path.exists(PINNED.format(c=cp)):
        md = pd.read_csv(PINNED.format(c=cp))
        csig = set(md[(md.fdr < 0.05) & (md.incremental_r2 > 0)].gene.astype(str))
        ov = [g for g in res.gene if g in csig]; rep = set(sig.gene) & set(ov)
        summ.update({"cptac_sig_on_rppa": len(ov), "replicated": len(rep),
                     "replication_rate": 100 * len(rep) / max(len(ov), 1)})
    pd.DataFrame([summ]).to_csv(f"{OUT}/replication_summary_{short}.csv", index=False)
    print(f"\n==================== EXTERNAL REPLICATION {proj} (RPPA) ====================")
    for k, v in summ.items():
        print(f" {k}: {v:.1f}" if isinstance(v, float) else f" {k}: {v}")
    print(f" -> {OUT}/{short}_rppa_results.csv , replication_summary_{short}.csv")


if __name__ == "__main__":
    main()
