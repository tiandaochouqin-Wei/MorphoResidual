#!/usr/bin/env python3
"""
Matched-capacity / random-projection control (plan.md Section 3: "gain != extra
capacity"). A skeptical reviewer will ask whether the incremental R2 of the 20
WSI PCs is just the free variance any 20 extra columns buy in a ridge fit. This
compares, per morphology-predictable gene, the incremental R2 of:

  (real)     mRNA + 20 real WSI PCs
  (gauss)    mRNA + 20 columns of N(0,1) noise        -> pure capacity control
  (randproj) mRNA + 20 random-projection features of the raw WSI embedding
             (1024-dim -> random 1024x20 Gaussian matrix; JL-style, keeps
             embedding info but scrambles the variance-optimal PCA directions)

Each control is averaged over K independent random draws (a single draw is noisy).
Headline: real incremental R2 must be strongly right-shifted vs BOTH controls.
gauss ~ 0 proves capacity alone doesn't inflate; real >> randproj shows the
structured top-variance morphology directions specifically carry the signal.

Everything is 5-fold CV closed-form ridge (identical to residual_analysis) so the
numbers are directly comparable. Thread limits set before numpy import.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import multiprocessing as mp
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent))
import residual_analysis as base
from residual_analysis import CV_FOLDS, N_PCS, OUT_DIR, cv_r2
from sklearn.decomposition import PCA

N_WORKERS = int(os.environ.get("MORPHO_N_WORKERS", 32))
SIG_CSV = Path(os.environ.get("MORPHO_SIG_CSV", str(OUT_DIR / "residual_results_wsibatch.csv")))
K_DRAWS = int(os.environ.get("MORPHO_K_DRAWS", 10))
RNG = np.random.RandomState(0)

_W = {}


def _init(rna, protein, wsi_pcs, gauss_list, randproj_list):
    _W["rna"] = rna
    _W["protein"] = protein
    _W["wsi"] = wsi_pcs
    _W["gauss"] = gauss_list          # list of K (n x 20) noise matrices
    _W["randproj"] = randproj_list    # list of K (n x 20) random-projection matrices


def _process_gene(gene):
    y = np.asarray(_W["protein"][gene])
    valid = ~np.isnan(y)
    if valid.sum() < 30:
        return None
    yv = y[valid]
    xr = np.asarray(_W["rna"][gene])[valid].reshape(-1, 1)
    r2_rna = cv_r2(xr, yv, CV_FOLDS)
    real = cv_r2(np.hstack([xr, _W["wsi"][valid]]), yv, CV_FOLDS) - r2_rna
    gauss = np.mean([cv_r2(np.hstack([xr, G[valid]]), yv, CV_FOLDS) - r2_rna
                     for G in _W["gauss"]])
    rproj = np.mean([cv_r2(np.hstack([xr, P[valid]]), yv, CV_FOLDS) - r2_rna
                     for P in _W["randproj"]])
    return {"gene": gene, "incr_real": real, "incr_gauss": gauss, "incr_randproj": rproj}


def main():
    rna = base.load_rna_matrix()
    protein = base.load_protein_matrix()
    wsi = base.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    print(f"\ncommon patients: {len(common)}")

    sig = pd.read_csv(SIG_CSV)
    sig_genes = sig[(sig["fdr"] < 0.05) & (sig["incremental_r2"] > 0)]["gene"].tolist()
    genes = [g for g in sig_genes if g in rna.columns and g in protein.columns]
    print(f"testing {len(genes)} morphology-predictable proteins (from {SIG_CSV.name})")

    rna_c = rna.loc[common]
    protein_c = protein.loc[common]
    wsi_raw = wsi.loc[common].values                       # (n, 1024) mean-pooled embedding
    n, d = wsi_raw.shape

    # svd_solver='full' is load-bearing: at this matrix shape sklearn's auto
    # solver is randomized SVD with random_state=None, which made this script
    # non-reproducible (GBM: 699/700/704/701 significant across four identical
    # runs, 90.4% set overlap). 'full' is exact and deterministic.
    pca = PCA(n_components=min(N_PCS, len(common) - 1, d), svd_solver="full")
    wsi_pcs = pca.fit_transform(wsi_raw)
    k = wsi_pcs.shape[1]
    print(f"WSI PCA: {k} comps; building {K_DRAWS} draws each of "
          f"N(0,1) noise and random-projection controls")

    # standardize the raw embedding before random projection (fair JL projection)
    wsi_std = (wsi_raw - wsi_raw.mean(0)) / (wsi_raw.std(0) + 1e-8)
    gauss_list = [RNG.standard_normal((n, k)) for _ in range(K_DRAWS)]
    randproj_list = [wsi_std @ RNG.standard_normal((d, k)) / np.sqrt(d) for _ in range(K_DRAWS)]

    rna_dict = {g: rna_c[g].values for g in genes}
    protein_dict = {g: protein_c[g].values for g in genes}

    print(f"parallelizing across {N_WORKERS} workers")
    rows = []
    done = 0
    with mp.Pool(N_WORKERS, initializer=_init,
                 initargs=(rna_dict, protein_dict, wsi_pcs, gauss_list, randproj_list)) as pool:
        for r in pool.imap_unordered(_process_gene, genes, chunksize=8):
            if r is not None:
                rows.append(r)
            done += 1
            if done % 400 == 0:
                print(f"  {done}/{len(genes)} done...")

    res = pd.DataFrame(rows)
    out = OUT_DIR / "randproj_control_results.csv"
    res.to_csv(out, index=False)

    def summ(col):
        return f"median={res[col].median():+.4f}  mean={res[col].mean():+.4f}  " \
               f">0 in {(res[col] > 0).mean():.0%}"
    print(f"\nincremental R2 across {len(res)} morphology-predictable proteins:")
    print(f"  real WSI PCs       : {summ('incr_real')}")
    print(f"  N(0,1) noise (cap) : {summ('incr_gauss')}")
    print(f"  random projection  : {summ('incr_randproj')}")

    # paired comparisons (sign test)
    for ctrl in ["incr_gauss", "incr_randproj"]:
        wins = (res["incr_real"] > res[ctrl]).mean()
        gap = (res["incr_real"] - res[ctrl]).median()
        print(f"\nreal vs {ctrl}: real larger in {wins:.0%} of genes, "
              f"median gap {gap:+.4f}")
    print(f"\nfull results: {out}")
    print("expected: noise ~0 (capacity does not inflate); real >> both controls")


if __name__ == "__main__":
    main()
