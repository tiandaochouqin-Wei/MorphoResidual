#!/usr/bin/env python3
"""
n_fragility.py -- how much of the median incremental R^2 is a property of having
exactly these n patients, and how much is the cross-validation fold assignment?

WHY. Every "retained X%" figure in this project is a median of the same functional over
~100 patients: the composition controls (2.7-67.9%), the sex control (83.7-93.1%), the
transcriptome-wide 40%. If that median moves a lot when a handful of patients leave, all
of them are fragile at once. The trigger was an incidental observation: in the covariate
run, dropping the 5 LUAD patients with no recorded race moved the cohort's median
increment from +0.0704 to +0.0210, and dropping 3 of GBM's 99 moved it from +0.0580 to
+0.0068 -- a larger loss than the entire operator batch correction.

BUT THAT OBSERVATION CONFOUNDS TWO THINGS. residual_analysis.cv_r2 uses
KFold(n_splits=5, shuffle=True, random_state=0): change n and every fold boundary moves,
so a subsample is not only smaller, it is cross-validated differently. This script
separates them.

    arm A  full n, CV seeds 0..N_SEEDS-1                  fold assignment alone
    arm B  subsample to a fraction of n, seed 0, repeated  sampling (plus its folds)

If arm B's spread at 95% of n is no wider than arm A's at full n, the "five patients"
effect is fold noise and the paper's existing seed-reliability analysis already covers
it. If B is clearly wider, the estimate is genuinely n-fragile and every retained-%
figure in the paper inherits that.

The PCs are fitted once on the full cohort and then row-subset, which is what the
covariate and operator runs did, so this reproduces their regime exactly. A genuinely
smaller cohort would also refit the PCA; that is a different (larger) effect and is not
what is being isolated here.

Env:  source morpho_env.sh <cancer>, plus
      MORPHO_SIG_CSV, MORPHO_COHORT
      MORPHO_NFRAG_SEEDS   default 10
      MORPHO_NFRAG_DRAWS   default 10
      MORPHO_NFRAG_FRACS   default 0.95,0.90,0.85,0.80
      MORPHO_NFRAG_OUT     default n_fragility_<cohort>.csv
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

MIN_PATIENTS = 30
RIDGE_ALPHA = 1.0


def cv_r2_seed(X, y, folds, seed):
    """residual_analysis.cv_r2 with the fold seed exposed. The arithmetic below is copied
    from residual_analysis.py lines 184-199 and must stay identical to it; only
    random_state is parameterised."""
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    preds = np.zeros_like(y, dtype=float)
    ridge_I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        X_tr, X_te, y_tr = X[tr], X[te], y[tr]
        x_mean = X_tr.mean(axis=0)
        x_std = X_tr.std(axis=0)
        x_std[x_std == 0] = 1.0
        Xtr_s = (X_tr - x_mean) / x_std
        Xte_s = (X_te - x_mean) / x_std
        y_mean = y_tr.mean()
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + ridge_I, Xtr_s.T @ (y_tr - y_mean))
        preds[te] = Xte_s @ w + y_mean
    ss_res = np.sum((y - preds) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def median_increment(rna, protein, pcs, genes, rows, folds, seed):
    """Median incremental R^2 over `genes`, on the patients in `rows`, at fold seed `seed`."""
    vals = []
    for g in genes:
        y = np.asarray(protein[g].values, float)[rows]
        v = ~np.isnan(y)
        if v.sum() < MIN_PATIENTS:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[rows][v].reshape(-1, 1)
        Wv = pcs[rows][v]
        r = cv_r2_seed(xr, yv, folds, seed)
        vals.append(cv_r2_seed(np.hstack([xr, Wv]), yv, folds, seed) - r)
    return (float(np.median(vals)) if vals else float("nan"), len(vals))


def main():
    cohort = os.environ.get("MORPHO_COHORT", "").strip().lower()
    sig_csv = os.environ.get("MORPHO_SIG_CSV")
    if not cohort or not sig_csv:
        sys.exit("[nfrag] set MORPHO_COHORT and MORPHO_SIG_CSV")
    out = os.environ.get("MORPHO_NFRAG_OUT", f"n_fragility_{cohort}.csv")
    n_seeds = int(os.environ.get("MORPHO_NFRAG_SEEDS", "10"))
    n_draws = int(os.environ.get("MORPHO_NFRAG_DRAWS", "10"))
    fracs = [float(x) for x in os.environ.get("MORPHO_NFRAG_FRACS", "0.95,0.90,0.85,0.80").split(",")]

    import batch_leak_check as blc
    paths = blc.cohort_paths(cohort)
    os.environ.update(
        MORPHO_ROOT=paths["root"], MORPHO_RNA_MANIFEST=paths["rna_manifest"],
        MORPHO_ALIQUOT_XWALK=paths["xwalk"], MORPHO_SLIDE_MAP=paths["slide_map"],
        MORPHO_PROTEIN_TSV=paths["protein"], MORPHO_WSI_EMB_DIR=paths["emb"])
    import residual_analysis as ra

    rna, protein, wsi = ra.load_rna_matrix(), ra.load_protein_matrix(), ra.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    gcols = sorted(set(rna.columns) & set(protein.columns))
    rna, protein = rna.loc[common, gcols], protein.loc[common, gcols]
    n = len(common)
    CV = getattr(ra, "CV_FOLDS", 5)
    n_pcs = min(getattr(ra, "N_PCS", 20), n - 1, wsi.shape[1])
    pcs = PCA(n_components=n_pcs, svd_solver="full").fit_transform(wsi.loc[common].values)

    md = pd.read_csv(sig_csv)
    sig = (md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]["gene"].astype(str).tolist()
           if {"fdr", "incremental_r2"}.issubset(md.columns) else md["gene"].astype(str).tolist())
    sig = [g for g in sig if g in rna.columns and g in protein.columns]
    cap = int(os.environ.get("MORPHO_MAX_GENES", "0") or 0)
    if cap:
        sig = sig[:cap]
        print(f"[nfrag] *** SMOKE TEST: {cap} genes -- NOT the real run ***", flush=True)
    print(f"[nfrag] {cohort}: n={n}, {len(sig)} genes, PCs {pcs.shape}, "
          f"{n_seeds} seeds / {n_draws} draws / fracs {fracs}", flush=True)

    allrows = np.arange(n)
    rows = []
    ref, ng = median_increment(rna, protein, pcs, sig, allrows, CV, 0)
    print(f"[nfrag] reference (full n, seed 0): median {ref:+.4f} over {ng} genes", flush=True)
    rows.append(dict(arm="reference", frac=1.0, n=n, seed=0, draw=0, median=ref, genes=ng))

    for s in range(n_seeds):
        m, k = median_increment(rna, protein, pcs, sig, allrows, CV, s)
        rows.append(dict(arm="A_fold_seed", frac=1.0, n=n, seed=s, draw=0, median=m, genes=k))
        print(f"[nfrag]   A seed {s:2d}: {m:+.4f}", flush=True)

    rng = np.random.RandomState(12345)
    for f in fracs:
        k_n = max(MIN_PATIENTS, int(round(f * n)))
        for d in range(n_draws):
            sub = np.sort(rng.choice(n, size=k_n, replace=False))
            m, k = median_increment(rna, protein, pcs, sig, sub, CV, 0)
            rows.append(dict(arm="B_subsample", frac=f, n=k_n, seed=0, draw=d, median=m, genes=k))
        got = [r["median"] for r in rows if r["arm"] == "B_subsample" and r["frac"] == f]
        print(f"[nfrag]   B frac {f:.2f} (n={k_n}): median of draws {np.median(got):+.4f}, "
              f"range {min(got):+.4f}..{max(got):+.4f}", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)

    A = res[res.arm == "A_fold_seed"]["median"]
    print("\n==================== n-FRAGILITY ====================")
    print(f" cohort {cohort}: reference median {ref:+.4f} (full n={n})")
    print(f" arm A, fold seed only      : median {A.median():+.4f}, "
          f"range {A.min():+.4f}..{A.max():+.4f}, spread {A.max() - A.min():.4f}")
    for f in fracs:
        B = res[(res.arm == "B_subsample") & (res.frac == f)]["median"]
        print(f" arm B, {100*f:.0f}% of n (n={int(res[(res.arm=='B_subsample') & (res.frac==f)].n.iloc[0])})"
              f"  : median {B.median():+.4f}, range {B.min():+.4f}..{B.max():+.4f}, "
              f"spread {B.max() - B.min():.4f}  ({B.max() - B.min():.4f} / {A.max() - A.min():.4f} = "
              f"{(B.max() - B.min()) / max(A.max() - A.min(), 1e-9):.1f}x the fold-only spread)")
    print(" READ IT AS: if B's spread is close to A's, losing a few patients is fold noise,")
    print(" not n-fragility. If B is several times wider, every retained-% median in the")
    print(" paper inherits that instability.")
    print(f" -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
