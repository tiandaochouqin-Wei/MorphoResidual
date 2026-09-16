#!/usr/bin/env python3
"""
Confound-ADJUSTED residual analysis (plan.md Section 3 controls). The tumor-only
run still flagged 2233/9635 genes as significant with the tell-tale confound
signature (p-value floor, r2_rna~0 -> r2_mrna_wsi~0.5 for thousands of unrelated
genes). This asks the honest question: does H&E morphology explain protein
residual BEYOND mRNA *and the known technical/biological confounds* -- namely
  (1) TMT plex/batch (23 plexes; plexes 19,20 are all-tumor -> imbalanced design,
      the documented ComBat-corrected batch effect in this cohort), and
  (2) tumor purity (RNA marker-based epithelial-vs-stromal proxy).

Baseline model:  protein ~ mRNA + purity + plex_dummies
Full model:      protein ~ mRNA + purity + plex_dummies + WSI_PCs
Incremental R2 of WSI is measured over this richer baseline, vs the SAME
permutation null (shuffle patient<->WSI). If the WSI signal collapses toward the
null -> it was batch/purity confound. If it survives with a structured,
pathway-enriched gene set -> that is the real, defensible morphology signal.

Thread limits are set at the very top (before numpy import) to prevent the
BLAS-oversubscription stall seen on gpu02 (each worker spawning many BLAS
threads -> thousands of threads thrashing).
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
import residual_analysis as base  # reuse validated loaders + cv_r2 + bh_fdr
from residual_analysis import (CV_FOLDS, FDR_THRESHOLD, GATE_MIN_PROTEINS, N_PCS,
                               N_PERM, MIN_PATIENTS_TO_RUN, OUT_DIR, cv_r2, bh_fdr)
from sklearn.decomposition import PCA

N_WORKERS = int(os.environ.get("MORPHO_N_WORKERS", 32))
PLEX_MAP = Path(os.environ.get("MORPHO_PLEX_MAP", str(Path(__file__).parent / "case_to_plex.tsv")))
RNG = np.random.RandomState(0)

# CCRCC epithelial (tumor) vs stromal/immune markers -> a simple purity proxy.
TUMOR_MARKERS = ["CA9", "NDUFA4L2", "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "PAX8", "VEGFA"]
STROMAL_MARKERS = ["PTPRC", "VIM", "COL1A1", "COL3A1", "PECAM1", "CD3D", "CD68", "ACTA2", "LUM"]


def build_purity(rna_common):
    tum = [g for g in TUMOR_MARKERS if g in rna_common.columns]
    stro = [g for g in STROMAL_MARKERS if g in rna_common.columns]
    print(f"purity proxy: {len(tum)}/{len(TUMOR_MARKERS)} tumor + {len(stro)}/{len(STROMAL_MARKERS)} stromal markers found")
    t = rna_common[tum]; s = rna_common[stro]
    t_z = (t - t.mean()) / t.std().replace(0, 1)
    s_z = (s - s.mean()) / s.std().replace(0, 1)
    return (t_z.mean(axis=1) - s_z.mean(axis=1)).values.reshape(-1, 1)


def build_plex_dummies(common):
    pm = pd.read_csv(PLEX_MAP, sep="\t").set_index("case_id")["plex"].to_dict()
    plex = np.array([pm.get(c, -1) for c in common])
    missing = (plex == -1).sum()
    if missing:
        print(f"  WARNING: {missing} patients have no plex label")
    # one-hot, drop first level to avoid collinearity with intercept
    levels = sorted(set(plex[plex != -1]))
    dummies = np.zeros((len(common), len(levels) - 1))
    for j, lv in enumerate(levels[1:]):
        dummies[:, j] = (plex == lv).astype(float)
    print(f"  plex dummies: {len(levels)} plexes -> {dummies.shape[1]} dummy columns")
    return dummies


_W = {}


def _init(rna_dict, protein_dict, baseline_extra, wsi_pcs, perm_indices):
    _W["rna"] = rna_dict; _W["protein"] = protein_dict
    _W["base_extra"] = baseline_extra; _W["wsi_pcs"] = wsi_pcs
    _W["perm"] = perm_indices


def _process_gene(gene):
    y = np.asarray(_W["protein"][gene])
    valid = ~np.isnan(y)
    if valid.sum() < 30:
        return None
    y_v = y[valid]
    x_rna = np.asarray(_W["rna"][gene])[valid].reshape(-1, 1)
    base_extra = _W["base_extra"][valid]        # purity + plex dummies
    wsi = _W["wsi_pcs"]
    X_base = np.hstack([x_rna, base_extra])
    X_full = np.hstack([x_rna, base_extra, wsi[valid]])
    r2_base = cv_r2(X_base, y_v, CV_FOLDS)
    r2_full = cv_r2(X_full, y_v, CV_FOLDS)
    incr = r2_full - r2_base

    null = np.empty(len(_W["perm"]))
    for i, perm in enumerate(_W["perm"]):
        X_perm = np.hstack([x_rna, base_extra, wsi[perm][valid]])
        null[i] = cv_r2(X_perm, y_v, CV_FOLDS) - r2_base
    pval = (np.sum(null >= incr) + 1) / (len(_W["perm"]) + 1)
    return {"gene": gene, "n": int(valid.sum()), "r2_base": r2_base,
            "r2_full": r2_full, "incremental_r2": incr, "pval": pval}


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
    plex_dum = build_plex_dummies(common)
    baseline_extra = np.hstack([purity, plex_dum])
    print(f"baseline covariates beyond mRNA: purity(1) + plex({plex_dum.shape[1]}) "
          f"= {baseline_extra.shape[1]} columns")

    # svd_solver='full' is load-bearing: at this matrix shape sklearn's auto
    # solver is randomized SVD with random_state=None, which made this script
    # non-reproducible (GBM: 699/700/704/701 significant across four identical
    # runs, 90.4% set overlap). 'full' is exact and deterministic.
    pca = PCA(n_components=min(N_PCS, len(common) - 1, wsi_raw.shape[1]),
              svd_solver="full")
    wsi_pcs = pca.fit_transform(wsi_raw)
    print(f"WSI PCA: {wsi_pcs.shape[1]} comps, {pca.explained_variance_ratio_.sum():.1%} var")

    n = len(common)
    perm_indices = [RNG.permutation(n) for _ in range(N_PERM)]
    rna_dict = {g: rna_c[g].values for g in common_genes}
    protein_dict = {g: protein_c[g].values for g in common_genes}

    print(f"parallelizing across {N_WORKERS} workers (baseline = mRNA + purity + plex)")
    results = []
    tested = 0
    with mp.Pool(N_WORKERS, initializer=_init,
                 initargs=(rna_dict, protein_dict, baseline_extra, wsi_pcs, perm_indices)) as pool:
        for r in pool.imap_unordered(_process_gene, common_genes, chunksize=8):
            tested += 1
            if r is not None:
                results.append(r)
            if tested % 200 == 0:
                print(f"  {tested}/{len(common_genes)} genes tested...")

    res = pd.DataFrame(results)
    res["fdr"] = bh_fdr(res["pval"].values)
    res = res.sort_values("incremental_r2", ascending=False)
    out_csv = OUT_DIR / "residual_results_adjusted.csv"
    res.to_csv(out_csv, index=False)

    sig = res[(res["fdr"] < FDR_THRESHOLD) & (res["incremental_r2"] > 0)]
    print(f"\ntested {len(res)} genes (baseline ADJUSTED for mRNA + purity + plex batch)")
    print(f"significant (FDR<{FDR_THRESHOLD}, incremental_r2>0): {len(sig)}")
    print(f"\ntop 15 by incremental R2:")
    print(res.head(15)[["gene", "n", "r2_base", "r2_full", "incremental_r2", "pval", "fdr"]]
          .to_string(index=False))
    gate = len(sig) >= GATE_MIN_PROTEINS
    print(f"\n{'PASS' if gate else 'FAIL'}: {len(sig)} proteins clear the pilot gate "
          f"({'>=' if gate else '<'}{GATE_MIN_PROTEINS}) AFTER confound adjustment")
    print(f"full results: {out_csv}")
    print("\n(compare with residual_results_tumoronly.csv: how much of the 2233 survives batch+purity control)")
    sys.exit(0 if gate else 1)


if __name__ == "__main__":
    main()
