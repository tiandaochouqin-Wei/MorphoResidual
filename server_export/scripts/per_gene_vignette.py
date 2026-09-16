#!/usr/bin/env python3
"""per_gene_vignette.py (SERVER, compute-only) — dump per-patient data for a few named
proteins so we can draw a single-protein vignette (protein-vs-mRNA + morphology recovers
the residual). Same estimand as the main analysis (OOF ridge, KFold seed 0, alpha=1,
train-fold z-standardisation).  CPU only, no matplotlib.

Run with the main-analysis python:
    cd /public/home/fjhui/ZW/scripts && source morpho_env.sh ccrcc
    export MORPHO_GENES=RPL35A,SSR3
    export MORPHO_OUT=/public/home/fjhui/ZW/results
    /public/home/fjhui/miniconda3/bin/python -u per_gene_vignette.py
Then WinSCP  $MORPHO_OUT/vignette_data.csv  to figures/figdata/ .
"""
import os, numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
import residual_analysis as RA

RIDGE = getattr(RA, "RIDGE_ALPHA", 1.0)
FOLDS = getattr(RA, "CV_FOLDS", 5)
GENES = os.environ.get("MORPHO_GENES", "RPL35A,SSR3").split(",")
OUT = os.environ.get("MORPHO_OUT", ".")


def cv_predict(X, y):
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    y = np.asarray(y, float)
    kf = KFold(n_splits=FOLDS, shuffle=True, random_state=0)
    preds = np.full(len(y), np.nan); I = RIDGE * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        Xtr, Xte, ytr = X[tr], X[te], y[tr]
        xm = Xtr.mean(0); xs = Xtr.std(0); xs[xs == 0] = 1.0; ym = ytr.mean()
        w = np.linalg.solve(((Xtr - xm) / xs).T @ ((Xtr - xm) / xs) + I, ((Xtr - xm) / xs).T @ (ytr - ym))
        preds[te] = ((Xte - xm) / xs) @ w + ym
    return preds


def r2(a, b):
    a = np.asarray(a); b = np.asarray(b)
    return 1 - np.sum((a - b) ** 2) / np.sum((a - a.mean()) ** 2)


rna = RA.load_rna_matrix(); protein = RA.load_protein_matrix(); wsi = RA.load_wsi_embeddings()
common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]
pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
common = np.array(common)
print(f"[vig] {len(common)} patients")

rows = []
for g in GENES:
    g = g.strip()
    if g not in rna.columns or g not in protein.columns:
        print(f"[vig] skip {g}: not in matrices"); continue
    y = protein[g].values.astype(float); x = rna[g].values.astype(float)
    v = np.isfinite(y) & np.isfinite(x)
    yr, xr, P, idx = y[v], x[v], pcs[v], common[v]
    resid = yr - cv_predict(xr, yr)
    morph = cv_predict(P, resid)
    for i in range(int(v.sum())):
        rows.append({"gene": g, "patient": idx[i], "mrna": xr[i], "protein": yr[i],
                     "residual": resid[i], "morph_pred": morph[i]})
    print(f"[vig] {g}: n={v.sum()}  mRNA->protein R2={r2(yr, cv_predict(xr, yr)):+.3f}  "
          f"morphology incremental R2 (residual)={r2(resid, morph):+.3f}")

pd.DataFrame(rows).to_csv(f"{OUT}/vignette_data.csv", index=False)
print(f"[out] {OUT}/vignette_data.csv  ({len(rows)} rows)  -> WinSCP to figures/figdata/")
