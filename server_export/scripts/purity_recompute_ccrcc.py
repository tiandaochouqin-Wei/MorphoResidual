#!/usr/bin/env python3
"""purity_recompute_ccrcc.py -- P0-3 fix (adversarial review, 2026-09-06), bypasses the
shared/mutable OUT_DIR entirely. Reads the PINNED CCRCC significant-gene results directly
(never from OUT_DIR, which is a shared path other jobs also write to -- that's what caused
the earlier 3,109-vs-2,191 confusion), and recomputes the tumour-purity-adjusted increment
for exactly those genes using the current residual_analysis loaders, restricted to the
protein-matrix patient set (confirmed = the correct 103-patient CCRCC tumour-only cohort).

Run on the server:
    cd /public/home/fjhui/ZW/scripts
    source morpho_env.sh ccrcc
    python3 -u purity_recompute_ccrcc.py

Sanity checks are hard asserts: if the pinned file or patient count doesn't match what the
manuscript already reports (9,635 genes / 2,191 significant / ~103 patients), the script
stops instead of silently producing wrong numbers.
"""
import numpy as np, pandas as pd
from residual_analysis import load_rna_matrix, load_protein_matrix, load_wsi_embeddings, cv_r2, CV_FOLDS
from sklearn.decomposition import PCA

PINNED = "/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/ccrcc/residual_results_tumoronly.csv"
OUT = "/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/ccrcc/confound_check_results_FIXED.csv"

TUMOR_MARKERS = ["CA9", "NDUFA4L2", "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "PAX8"]
STROMAL_MARKERS = ["PTPRC", "VIM", "COL1A1", "COL3A1", "PECAM1", "CD3D", "CD68", "ACTA2"]

results = pd.read_csv(PINNED)
sig_check = ((results.fdr < 0.05) & (results.incremental_r2 > 0)).sum()
print(f"[sanity] pinned file: shape={results.shape}, sig(fdr<0.05 & incr>0)={sig_check}  (must be 2191)")
assert results.shape[0] == 9635 and sig_check == 2191, "pinned file does not match the expected CCRCC snapshot -- STOP, do not proceed"

rna = load_rna_matrix()
protein = load_protein_matrix()
wsi = load_wsi_embeddings()
print(f"[loaders] rna {rna.shape}  protein {protein.shape}  wsi {wsi.shape}")

# anchor the patient set to the protein matrix (confirmed correctly CCRCC-scoped, n=103),
# not to whatever rna.index happens to contain (which was seen loading 254 patients once --
# likely a shared/generic RNA matrix that needs intersecting down, not a bug per se, but we
# don't trust it blindly)
common = sorted(set(protein.index) & set(rna.index) & set(wsi.index))
print(f"[loaders] common patients (anchored on protein): {len(common)}  (expect ~103)")
assert len(common) in (102, 103, 104), f"unexpected patient count {len(common)} -- STOP, investigate before trusting output"

common_genes = sorted(set(rna.columns) & set(protein.columns) & set(results["gene"]))
rna_c = rna.loc[common, common_genes]
protein_c = protein.loc[common, common_genes]
wsi_raw = wsi.loc[common].values
pca = PCA(n_components=min(20, len(common) - 1, wsi_raw.shape[1]), svd_solver="full")
wsi_pcs = pca.fit_transform(wsi_raw)
print(f"[loaders] genes tested: {len(common_genes)}  wsi PCs: {wsi_pcs.shape[1]}")

purity_tumor = [g for g in TUMOR_MARKERS if g in rna_c.columns]
purity_stroma = [g for g in STROMAL_MARKERS if g in rna_c.columns]
t = rna_c[purity_tumor]; s = rna_c[purity_stroma]
t_z = (t - t.mean()) / t.std().replace(0, 1); s_z = (s - s.mean()) / s.std().replace(0, 1)
purity = (t_z.mean(axis=1) - s_z.mean(axis=1)).values
print(f"[purity] proxy from {len(purity_tumor)}/{len(TUMOR_MARKERS)} tumour + "
      f"{len(purity_stroma)}/{len(STROMAL_MARKERS)} stromal markers")

rows = []
for j, g in enumerate(common_genes):
    y = protein_c[g].values
    valid = ~np.isnan(y)
    if valid.sum() < 30:
        continue
    y_v = y[valid]
    x_rna = rna_c[g].values[valid].reshape(-1, 1)
    x_purity = purity[valid].reshape(-1, 1)
    x_wsi = wsi_pcs[valid]
    r2_base = cv_r2(np.hstack([x_rna, x_purity]), y_v, CV_FOLDS)
    r2_full = cv_r2(np.hstack([x_rna, x_purity, x_wsi]), y_v, CV_FOLDS)
    rows.append({"gene": g, "incr_r2_purity_adjusted": r2_full - r2_base})
    if (j + 1) % 2000 == 0:
        print(f"  {j+1}/{len(common_genes)} genes done...")

adj = pd.DataFrame(rows)
merged = results.merge(adj, on="gene", how="inner")
sig = merged[(merged["fdr"] < 0.05) & (merged["incremental_r2"] > 0)]
print(f"\n[RESULT] significant genes with a purity-adjustment computed: {len(sig)}  (must be close to 2191)")
print(f"[RESULT] median incremental_r2 (original): {sig['incremental_r2'].median():.4f}")
print(f"[RESULT] median incremental_r2 (purity-adjusted): {sig['incr_r2_purity_adjusted'].median():.4f}")
print(f"[RESULT] median retained: {100*np.median(sig['incr_r2_purity_adjusted']/sig['incremental_r2']):.0f}%")
print(f"[RESULT] frac still positive: {100*(sig['incr_r2_purity_adjusted']>0).mean():.0f}%")
merged.to_csv(OUT, index=False)
print(f"\nwrote {OUT}")
