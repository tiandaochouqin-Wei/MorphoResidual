#!/usr/bin/env python3
"""kirc_residual.py — B-11 external replication (TCGA-KIRC, independent consortium +
RPPA platform). Tests morphology's incremental R^2 over an mRNA baseline for the
RPPA-protein residual, with the SAME estimand as the CPTAC analysis (embedded cv_r2,
5-fold OOF, closed-form ridge lambda=1). Self-contained (no residual_analysis import).

Inputs (all under /public/home/fjhui/ZW/tcga_kirc):
  rppa_kirc.csv          cases x genes (RPPA level-4 total protein)
  rna_case_file.tsv      case<->STAR filename; STAR files in rna/
  emb_phikon/*.pt (or MORPHO_WSI_EMB_DIR)   Phikon slide embeddings for KIRC
                          slides -- deliberately NOT UNI: external replication
                          uses an encoder independent of the discovery pipeline
                          (see main.tex, external-replication Methods). Set
                          MORPHO_WSI_EMB_DIR to point at a UNI export instead if
                          you specifically want the (non-published) UNI variant.
Run:  MORPHO_WSI_EMB_DIR=.../emb_phikon MORPHO_N_WORKERS=8 python -u kirc_residual.py
"""
import os, glob, sys
import numpy as np, pandas as pd
import torch
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from scipy import stats as st

ROOT = "/public/home/fjhui/ZW/tcga_kirc"
EMB = os.environ.get("MORPHO_WSI_EMB_DIR", f"{ROOT}/emb_phikon")
OUT = os.environ.get("MORPHO_OUT", f"{ROOT}/results")
os.makedirs(OUT, exist_ok=True)
CV_FOLDS, N_PERM, RIDGE_ALPHA = 5, 1000, 1.0
PINNED_CCRCC = "/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/ccrcc/residual_results_tumoronly.csv"


def tcga_case(s):
    return "-".join(str(s).split("-")[:3])


def cv_r2(X, y, folds=CV_FOLDS):
    kf = KFold(n_splits=folds, shuffle=True, random_state=0)
    preds = np.zeros_like(y, dtype=float)
    ridge_I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        Xtr, Xte, ytr = X[tr], X[te], y[tr]
        xm = Xtr.mean(0); xs = Xtr.std(0); xs[xs == 0] = 1.0
        Xtr_s = (Xtr - xm) / xs; Xte_s = (Xte - xm) / xs; ym = ytr.mean()
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + ridge_I, Xtr_s.T @ (ytr - ym))
        preds[te] = Xte_s @ w + ym
    ss_res = np.sum((y - preds) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
    return q


def main():
    # ---- RPPA protein ----
    rppa = pd.read_csv(f"{ROOT}/rppa_kirc.csv", index_col=0)
    rppa.index = [tcga_case(c) for c in rppa.index]
    rppa = rppa[~rppa.index.duplicated()]

    # ---- RNA (STAR TPM -> log2(TPM+1)) ----
    xw = pd.read_csv(f"{ROOT}/rna_case_file.tsv", sep="\t")
    rows = {}
    for _, r in xw.iterrows():
        fp = f"{ROOT}/rna/{r['filename']}"
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp, sep="\t", comment="#")
        d.columns = [c.strip() for c in d.columns]
        d = d[["gene_name", "tpm_unstranded"]].dropna()
        d = d[~d["gene_name"].astype(str).str.startswith("N_")]
        s = d.groupby("gene_name")["tpm_unstranded"].mean()
        rows.setdefault(tcga_case(r["case"]), []).append(s)
    rna = pd.DataFrame({c: pd.concat(v, axis=1).mean(axis=1) for c, v in rows.items()}).T
    rna = np.log2(rna.fillna(0) + 1)

    # ---- WSI UNI embeddings -> per-case mean-pool ----
    emb = {}
    for pt in glob.glob(f"{EMB}/*.pt"):
        case = tcga_case(os.path.basename(pt))
        e = torch.load(pt, map_location="cpu")["embeddings"].float().mean(0).numpy()
        emb.setdefault(case, []).append(e)
    wsi = pd.DataFrame({c: np.mean(v, axis=0) for c, v in emb.items()}).T
    if wsi.empty:
        sys.exit(f"no embeddings in {EMB} -- extract UNI features from the KIRC slides first")

    common = sorted(set(rppa.index) & set(rna.index) & set(wsi.index))
    print(f"cases: rppa {len(rppa)}, rna {len(rna)}, wsi {len(wsi)}, common {len(common)}")
    rppa, rna, wsi = rppa.loc[common], rna.loc[common], wsi.loc[common]
    wsi_pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
    perm = [np.random.RandomState(i).permutation(len(common)) for i in range(N_PERM)]

    genes = [g for g in rppa.columns if g in rna.columns]
    print(f"testing {len(genes)} RPPA genes present in RNA")
    res = []
    for g in genes:
        y = np.asarray(rppa[g].values, float); v = ~np.isnan(y)
        if v.sum() < 30:
            continue
        yv = y[v]; xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1); Wv = wsi_pcs[v]
        r2r = cv_r2(xr, yv); r2b = cv_r2(np.hstack([xr, Wv]), yv); incr = r2b - r2r
        null = np.array([cv_r2(np.hstack([xr, wsi_pcs[p][v]]), yv) - r2r for p in perm])
        pval = (np.sum(null >= incr) + 1) / (len(perm) + 1)
        res.append({"gene": g, "n": int(v.sum()), "r2_rna": r2r, "r2_both": r2b,
                    "incremental_r2": incr, "pval": pval})
    res = pd.DataFrame(res)
    res["fdr"] = bh_fdr(res["pval"].values)
    res.to_csv(f"{OUT}/kirc_rppa_results.csv", index=False)

    sig = res[(res["fdr"] < 0.05) & (res["incremental_r2"] > 0)]
    print("\n==================== B-11 EXTERNAL REPLICATION (TCGA-KIRC / RPPA) ====================")
    print(f" tested {len(res)} RPPA proteins; SIGNIFICANT {len(sig)} (FDR<0.05, incr>0) "
          f"= {100*len(sig)/max(len(res),1):.0f}%")
    print(f" median incremental R2 (significant): {sig['incremental_r2'].median():+.4f}")
    # gene-specific agreement with CPTAC-CCRCC
    md = pd.read_csv(PINNED_CCRCC)
    ccrcc_sig = set(md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]["gene"].astype(str))
    ov = [g for g in res["gene"] if g in ccrcc_sig]
    kirc_sig_ov = set(sig["gene"]) & set(ov)
    print(f" of {len(ov)} CCRCC-significant genes also on RPPA: {len(kirc_sig_ov)} replicate in KIRC "
          f"({100*len(kirc_sig_ov)/max(len(ov),1):.0f}%)")
    print(f" -> {OUT}/kirc_rppa_results.csv")


if __name__ == "__main__":
    main()
