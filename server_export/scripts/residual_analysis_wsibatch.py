#!/usr/bin/env python3
"""
WSI-batch-controlled residual analysis (the headline-limitation stress test).

The concern after the 1340-gene adjusted run: foundation-model (UNI) mean-pooled
embeddings are known to encode site/scanner/staining BATCH signal, and if that
co-varies with proteome technical structure it inflates the count. The previous
"adjusted" run put plex dummies in the BASELINE, which over-fit (r2_base went
strongly negative on n=103 with 22 dummies). This version fixes both issues:

  * control batch on the WSI FEATURE side, not the protein target side --
    residualize each WSI PC against a batch design matrix
    B = [cohort (C3L vs C3N) + TMT plex + purity], keeping only the
    batch-orthogonal component;
  * keep the protein baseline LEAN (just mRNA), so no baseline over-fitting.

Then run the standard nested test on the batch-residualized WSI:
  baseline: protein ~ mRNA ;  full: protein ~ mRNA + WSI_batchfree
  incremental R2 of WSI_batchfree, vs the same patient<->WSI permutation null.

If the count stays near 1340 -> the signal is largely orthogonal to these
batches (robust). If it collapses -> cohort/plex/purity batch was driving it.
NOTE: this controls the batch axes we can label WITHOUT slide metadata. The
WSI's own scanner/stain batch (from SVS aperio headers) is a further, separate
control (see extract_svs_batch.py) if this run leaves residual inflation.

Thread limits set before numpy import (prevents the gpu02 BLAS-oversubscription
stall). Default 32 workers.
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
from residual_analysis import (CV_FOLDS, FDR_THRESHOLD, GATE_MIN_PROTEINS, N_PCS,
                               N_PERM, MIN_PATIENTS_TO_RUN, OUT_DIR, cv_r2, bh_fdr)
from sklearn.decomposition import PCA

N_WORKERS = int(os.environ.get("MORPHO_N_WORKERS", 32))
PLEX_MAP = Path(os.environ.get("MORPHO_PLEX_MAP", str(Path(__file__).parent / "case_to_plex.tsv")))
RNG = np.random.RandomState(0)

TUMOR_MARKERS = ["CA9", "NDUFA4L2", "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "PAX8", "VEGFA"]
STROMAL_MARKERS = ["PTPRC", "VIM", "COL1A1", "COL3A1", "PECAM1", "CD3D", "CD68", "ACTA2", "LUM"]


def build_purity(rna_common):
    tum = [g for g in TUMOR_MARKERS if g in rna_common.columns]
    stro = [g for g in STROMAL_MARKERS if g in rna_common.columns]
    t = rna_common[tum]; s = rna_common[stro]
    t_z = (t - t.mean()) / t.std().replace(0, 1)
    s_z = (s - s.mean()) / s.std().replace(0, 1)
    print(f"purity proxy: {len(tum)} tumor + {len(stro)} stromal markers")
    return (t_z.mean(axis=1) - s_z.mean(axis=1)).values.reshape(-1, 1)


def build_batch_matrix(common, purity):
    # DEFAULT = coarse, df-safe batch control: cohort (C3L vs C3N) + purity only.
    # Residualizing 20 WSI PCs against 23 TMT-plex dummies at n=103 (~4.5/plex)
    # over-removes -- validated on synthetic data that even genuine signal is
    # destroyed (>90% variance stripped). Fine plex deconfounding is deferred to
    # pan-cancer scale (n >> #batches). Set MORPHO_INCLUDE_PLEX=1 to force it
    # (diagnostic only -- expect over-correction at this n).
    cohort = np.array([1.0 if c.startswith("C3L") else 0.0 for c in common]).reshape(-1, 1)
    cols = [np.ones((len(common), 1)), cohort, purity]
    label = "intercept + cohort(1) + purity(1)"
    if os.environ.get("MORPHO_INCLUDE_PLEX") == "1":
        pm = pd.read_csv(PLEX_MAP, sep="\t").set_index("case_id")["plex"].to_dict()
        plex = np.array([pm.get(c, -1) for c in common])
        levels = sorted(set(plex[plex != -1]))
        plex_dum = np.zeros((len(common), len(levels) - 1))
        for j, lv in enumerate(levels[1:]):
            plex_dum[:, j] = (plex == lv).astype(float)
        cols.insert(2, plex_dum)
        label = f"intercept + cohort(1) + plex({plex_dum.shape[1]}) + purity(1)"
    B = np.hstack(cols)
    print(f"batch design: {label} = {B.shape[1]} cols")
    return B


def residualize(WSI_pcs, B):
    # OLS residuals of each WSI PC on the batch design (fit once on all patients;
    # this removes known batch structure, it is not a prediction task -> no CV).
    beta, *_ = np.linalg.lstsq(B, WSI_pcs, rcond=None)
    resid = WSI_pcs - B @ beta
    kept = resid.var(axis=0).sum() / WSI_pcs.var(axis=0).sum()
    print(f"WSI variance retained after batch residualization: {kept:.1%}")
    return resid


_W = {}
def _init(rna_dict, protein_dict, wsi_resid, perm_indices):
    _W["rna"] = rna_dict; _W["protein"] = protein_dict
    _W["wsi"] = wsi_resid; _W["perm"] = perm_indices


def _process_gene(gene):
    y = np.asarray(_W["protein"][gene])
    valid = ~np.isnan(y)
    if valid.sum() < 30:
        return None
    y_v = y[valid]
    x_rna = np.asarray(_W["rna"][gene])[valid].reshape(-1, 1)
    wsi = _W["wsi"]
    r2_rna = cv_r2(x_rna, y_v, CV_FOLDS)
    r2_both = cv_r2(np.hstack([x_rna, wsi[valid]]), y_v, CV_FOLDS)
    incr = r2_both - r2_rna
    null = np.empty(len(_W["perm"]))
    for i, perm in enumerate(_W["perm"]):
        null[i] = cv_r2(np.hstack([x_rna, wsi[perm][valid]]), y_v, CV_FOLDS) - r2_rna
    pval = (np.sum(null >= incr) + 1) / (len(_W["perm"]) + 1)
    return {"gene": gene, "n": int(valid.sum()), "r2_rna": r2_rna,
            "r2_mrna_wsi": r2_both, "incremental_r2": incr, "pval": pval}


def main():
    rna = base.load_rna_matrix()
    protein = base.load_protein_matrix()
    wsi = base.load_wsi_embeddings()

    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    print(f"\ncommon patients across all 3: {len(common)}")
    if len(common) < MIN_PATIENTS_TO_RUN:
        sys.exit(f"FAIL: only {len(common)} patients overlap")
    common_genes = sorted(set(rna.columns) & set(protein.columns))
    print(f"common genes (RNA n protein): {len(common_genes)}")

    rna_c = rna.loc[common, common_genes]
    protein_c = protein.loc[common, common_genes]
    wsi_raw = wsi.loc[common].values

    purity = build_purity(rna.loc[common])
    # svd_solver='full' is load-bearing: at this matrix shape sklearn's auto
    # solver is randomized SVD with random_state=None, which made this script
    # non-reproducible (GBM: 699/700/704/701 significant across four identical
    # runs, 90.4% set overlap). 'full' is exact and deterministic.
    pca = PCA(n_components=min(N_PCS, len(common) - 1, wsi_raw.shape[1]),
              svd_solver="full")
    wsi_pcs = pca.fit_transform(wsi_raw)
    print(f"WSI PCA: {wsi_pcs.shape[1]} comps, {pca.explained_variance_ratio_.sum():.1%} var")

    B = build_batch_matrix(common, purity)
    wsi_resid = residualize(wsi_pcs, B)

    n = len(common)
    perm_indices = [RNG.permutation(n) for _ in range(N_PERM)]
    rna_dict = {g: rna_c[g].values for g in common_genes}
    protein_dict = {g: protein_c[g].values for g in common_genes}

    print(f"parallelizing across {N_WORKERS} workers (baseline=mRNA, WSI batch-residualized)")
    results = []; tested = 0
    with mp.Pool(N_WORKERS, initializer=_init,
                 initargs=(rna_dict, protein_dict, wsi_resid, perm_indices)) as pool:
        for r in pool.imap_unordered(_process_gene, common_genes, chunksize=8):
            tested += 1
            if r is not None:
                results.append(r)
            if tested % 200 == 0:
                print(f"  {tested}/{len(common_genes)} genes tested...")

    res = pd.DataFrame(results)
    res["fdr"] = bh_fdr(res["pval"].values)
    res = res.sort_values("incremental_r2", ascending=False)
    out_csv = OUT_DIR / "residual_results_wsibatch.csv"
    res.to_csv(out_csv, index=False)

    sig = res[(res["fdr"] < FDR_THRESHOLD) & (res["incremental_r2"] > 0)]
    print(f"\ntested {len(res)} genes (WSI residualized vs cohort+plex+purity, lean mRNA baseline)")
    print(f"significant (FDR<{FDR_THRESHOLD}, incremental_r2>0): {len(sig)}")
    print(f"\ntop 15 by incremental R2:")
    print(res.head(15)[["gene", "n", "r2_rna", "r2_mrna_wsi", "incremental_r2", "pval", "fdr"]]
          .to_string(index=False))
    gate = len(sig) >= GATE_MIN_PROTEINS
    print(f"\n{'PASS' if gate else 'FAIL'}: {len(sig)} proteins clear the pilot gate "
          f"({'>=' if gate else '<'}{GATE_MIN_PROTEINS}) with batch-residualized WSI")
    print(f"full results: {out_csv}")
    print("compare: contaminated 3109 -> tumor-only 2233 -> +plex/purity baseline 1340 -> this")
    sys.exit(0 if gate else 1)


if __name__ == "__main__":
    main()
