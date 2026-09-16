#!/usr/bin/env python3
"""T1-7 step 3 (NC_ROADMAP.md): in-fold PCA+ridge refit of the H&E-only family
score across all five cohorts, using the exact main-analysis design (5-fold
KFold random_state=0, PCA(20, svd_solver='full') fit WITHIN each training fold,
standardised closed-form ridge lambda=1) applied to the full mean-pooled UNI
embedding in figures/figdata/loco_<cohort>.npz. This directly answers the P1-10
concern that the CCRCC-only ABMIL-vs-mean-pool comparison in Fig. 2j used a
mean-pool baseline computed on a subsampled tile bag (<=2000 tiles) rather than
the full-tile mean-pool the main analysis actually uses; this script uses the
main analysis's own already-mean-pooled embedding, so it is the right baseline.

Writes: figures/figdata/loco_infold_scores.csv, SuppTable_locoinfold.tex
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from scipy import stats as st

DD = "figures/figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
NICE = {"ccrcc": "CCRCC", "luad": "LUAD", "ucec": "UCEC", "gbm": "GBM", "pdac": "PDAC"}
FAMS = ["translation_meas", "secretion_ER_meas"]
FAM_NICE = {"translation_meas": "translation", "secretion_ER_meas": "ER-secretion"}


def infold_ridge_oof(X, y, n_splits=5, n_pcs=20, lam=1.0, seed=0):
    """Exact main-analysis design, but PCA fit strictly within each training
    fold (not once on all patients) -- the point of this script."""
    n = len(y)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof = np.full(n, np.nan)
    for tr, te in kf.split(X):
        pca = PCA(n_components=n_pcs, svd_solver="full", random_state=seed)
        Xtr = pca.fit_transform(X[tr])
        Xte = pca.transform(X[te])
        mu, sd = Xtr.mean(0), Xtr.std(0); sd[sd == 0] = 1.0
        Xtr = (Xtr - mu) / sd; Xte = (Xte - mu) / sd
        ytr_mu = y[tr].mean()
        w = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(n_pcs), Xtr.T @ (y[tr] - ytr_mu))
        oof[te] = Xte @ w + ytr_mu
    return oof


rows = []
for c in COH:
    d = np.load(f"{DD}/loco_{c}.npz")
    X = d["emb"].astype(float)
    for fam in FAMS:
        if fam not in d.files:
            continue
        y = d[fam].astype(float)
        m = np.isfinite(y)
        oof = infold_ridge_oof(X[m], y[m])
        r, p = st.pearsonr(y[m], oof)
        r2 = 1 - np.sum((y[m] - oof) ** 2) / np.sum((y[m] - y[m].mean()) ** 2)
        rows.append(dict(cohort=NICE[c], family=FAM_NICE[fam], n=int(m.sum()), r=r, R2=r2, p=p))
        print(f"{NICE[c]:6s} {FAM_NICE[fam]:12s} n={int(m.sum()):3d}  in-fold-PCA mean-pool r={r:.3f} R2={r2:.3f} p={p:.2e}")

df = pd.DataFrame(rows)
df.round(4).to_csv(f"{DD}/loco_infold_scores.csv", index=False)

# CCRCC translation vs the existing ABMIL numbers already reported in the paper
ccrcc_transl = df[(df.cohort == "CCRCC") & (df.family == "translation")].iloc[0]
print(f"\nCCRCC translation, in-fold-PCA full-tile mean-pool: r={ccrcc_transl.r:.3f}, R2={ccrcc_transl.R2:.3f}")
print("compare to main text: ABMIL r=0.59/R2=0.33 (Fig. 2j); previously reported subsampled-bag "
      "mean-pool r=0.41/R2=0.14")

lines = [r"\begin{table}[htbp]\centering\small",
         r"\caption{In-fold aggregation check, all five cohorts. Out-of-fold Pearson $r$ and $R^2$ for "
         r"predicting the measured translation and ER-secretion residual scores from the full mean-pooled "
         r"UNI tile embedding, with PCA(20) refit strictly within each training fold of an identical "
         r"5-fold split (random\_state=0) -- the same design as the main analysis, extended from the "
         r"CCRCC-only ABMIL comparison (Fig.~\ref{fig:master2}j, Table~\ref{tab:abmil}) to all five "
         r"cohorts, and using the main analysis's own full-tile mean-pooled embedding rather than a "
         r"subsampled tile bag.}",
         r"\label{tab:locoinfold}",
         r"\begin{tabular}{llrrrr}\toprule",
         r"Cohort & Family & $n$ & $r$ & $R^2$ & $p$ \\ \midrule"]
for _, r in df.iterrows():
    lines.append(f"{r['cohort']} & {r['family']} & {int(r['n'])} & {r['r']:.3f} & {r['R2']:.3f} & {r['p']:.1e} \\\\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
with open("SuppTable_locoinfold.tex", "w", encoding="utf8") as f:
    f.write("\n".join(lines) + "\n")
print(f"\nwrote SuppTable_locoinfold.tex ({len(df)} rows)")
