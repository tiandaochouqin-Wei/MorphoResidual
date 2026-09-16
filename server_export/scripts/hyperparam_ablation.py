#!/usr/bin/env python3
"""hyperparam_ablation.py (SERVER, per cohort, CPU) — robustness of the per-gene morphology
increment to the two analysis hyper-parameters: number of morphology PCs and ridge lambda.
Recomputes incremental R2 (no permutation) for the pinned significant genes under a grid and
reports median increment, fraction > 0, and per-gene correlation with the default (20 PCs,
lambda=1).  Same CV scheme as the main analysis.
Run (CCRCC first; loop cohorts if wanted):
    cd /public/home/fjhui/ZW/scripts && source morpho_env.sh ccrcc
    export MORPHO_COHORT=ccrcc
    export MORPHO_SIG_CSV=/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/ccrcc/residual_results_tumoronly.csv
    export MORPHO_OUT=/public/home/fjhui/ZW/results
    /public/home/fjhui/miniconda3/bin/python -u hyperparam_ablation.py
Then WinSCP  $MORPHO_OUT/ablation_<cohort>.csv  to figures/figdata/ .
"""
import os, numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
import residual_analysis as RA

C = os.environ.get("MORPHO_COHORT", "cohort")
SIG = os.environ["MORPHO_SIG_CSV"]
OUT = os.environ.get("MORPHO_OUT", ".")
GRID = [(5, 1.0), (10, 1.0), (20, 1.0), (40, 1.0), (20, 0.1), (20, 10.0)]

rna = RA.load_rna_matrix(); protein = RA.load_protein_matrix(); wsi = RA.load_wsi_embeddings()
common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]
md = pd.read_csv(SIG)
genes = [g for g in md[(md.fdr < 0.05) & (md.incremental_r2 > 0)]["gene"].astype(str)
         if g in rna.columns and g in protein.columns]
print(f"[abl] {C}: {len(common)} patients, {len(genes)} significant genes, grid {GRID}")


def cv_r2(X, y, alpha):
    kf = KFold(n_splits=5, shuffle=True, random_state=0)
    preds = np.zeros_like(y, dtype=float); I = alpha * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        Xtr, Xte, ytr = X[tr], X[te], y[tr]
        xm = Xtr.mean(0); xs = Xtr.std(0); xs[xs == 0] = 1.0; ym = ytr.mean()
        A = (Xtr - xm) / xs
        w = np.linalg.solve(A.T @ A + I, A.T @ (ytr - ym))
        preds[te] = ((Xte - xm) / xs) @ w + ym
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - np.sum((y - preds) ** 2) / ss_tot if ss_tot > 0 else 0.0


pc_cache = {}
def pcs(n):
    if n not in pc_cache:
        pc_cache[n] = PCA(n_components=n, svd_solver="full").fit_transform(wsi.values)
    return pc_cache[n]


incr = {}
for (n, a) in GRID:
    P = pcs(n); vals = []
    for g in genes:
        y = protein[g].values.astype(float); v = ~np.isnan(y)
        if v.sum() < 30:
            vals.append(np.nan); continue
        yv = y[v]; xr = rna[g].values.astype(float)[v].reshape(-1, 1)
        vals.append(cv_r2(np.hstack([xr, P[v]]), yv, a) - cv_r2(xr, yv, a))
    incr[(n, a)] = np.array(vals, float)
    print(f"[abl] pcs={n:2d} lambda={a:<4}: median incr {np.nanmedian(incr[(n, a)]):+.3f}")

ref = incr[(20, 1.0)]
rows = []
for (n, a), v in incr.items():
    m = np.isfinite(v) & np.isfinite(ref)
    rows.append(dict(cohort=C, n_pcs=n, ridge_lambda=a, n_genes=int(m.sum()),
                     median_incr=float(np.nanmedian(v)), frac_pos=float(np.nanmean(v > 0)),
                     corr_vs_default=float(np.corrcoef(v[m], ref[m])[0, 1])))
res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/ablation_{C}.csv", index=False)
print(res.round(3).to_string(index=False))
print(f"[out] {OUT}/ablation_{C}.csv -> WinSCP to figures/figdata/")
