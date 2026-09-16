#!/usr/bin/env python3
"""
Confound diagnostic for the pilot residual-analysis result (plan.md Section 3
"Controls (rule out boring explanations)"). The main run found 3109/9675
genes "significant" (32%), which is suspiciously high -- many p-values sit at
the permutation floor (0.000999) and R2 jumps from ~0/negative (mRNA alone)
to 0.5-0.8 (mRNA+WSI) for thousands of unrelated genes. That pattern is the
signature of a shared confound (something that moves WSI appearance AND bulk
protein levels together for many genes at once), not gene-specific
post-transcriptional regulation. Runs three of plan.md's pre-registered
controls:

  A. mRNA-strength confound -- does incremental_r2 correlate with the
     baseline mRNA->protein R2? (plan.md: "not just easy proteins")
  B. Batch/site confound -- does TMT plex (derived from the RAW column order
     of the protein TSV, since consecutive columns between reference/QC
     channels form one plex block) predict the WSI PCs?
  C. Tumor-purity confound -- does an RNA-based epithelial-vs-stromal marker
     score predict the WSI PCs, and does incremental R2 shrink once that
     purity proxy is forced into the baseline model alongside mRNA?

Run on the same node/env as residual_analysis.py (needs its
residual_results.csv already written). No GPU needed.
"""
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent))
from residual_analysis import (  # noqa: E402
    ALIQUOT_XWALK, OUT_DIR, PROTEIN_TSV, RIDGE_ALPHA, CV_FOLDS,
    cv_r2, load_protein_matrix, load_rna_matrix, load_wsi_embeddings,
)
from sklearn.decomposition import PCA  # noqa: E402

N_PCS = 20
N_PERM_BATCH = 2000
RNG = np.random.RandomState(1)

TUMOR_MARKERS = ["CA9", "NDUFA4L2", "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "PAX8"]
STROMAL_MARKERS = ["PTPRC", "VIM", "COL1A1", "COL3A1", "PECAM1", "CD3D", "CD68", "ACTA2"]


def detect_plex_blocks():
    """Consecutive case-sample columns between reference/QC channels = one
    TMT plex. Reference channels are the ` Log Ratio` columns whose aliquot
    ID is NOT in aliquot_to_case.tsv (excluded there as QC/NCI7 channels)."""
    xwalk = pd.read_csv(ALIQUOT_XWALK, sep="\t").set_index("aliquot_id")["case_id"].to_dict()
    df = pd.read_csv(PROTEIN_TSV, sep="\t", index_col=0, nrows=1)
    log_ratio_cols = [c for c in df.columns if c.endswith(" Log Ratio") and "Unshared" not in c]

    blocks = []
    current = []
    for c in log_ratio_cols:
        aliquot = c.split(" ")[0]
        case = xwalk.get(aliquot)
        if case is None:  # reference/QC channel -> block boundary
            if current:
                blocks.append(current)
                current = []
            continue
        current.append(case)
    if current:
        blocks.append(current)

    plex_of_case = {}
    for i, block in enumerate(blocks):
        for case in block:
            plex_of_case[case] = i
    sizes = [len(b) for b in blocks]
    print(f"[Check B] detected {len(blocks)} plex blocks from column order, sizes: {sizes}")
    return plex_of_case


def groupmean_r2(y, labels):
    df = pd.DataFrame({"y": y, "g": labels})
    grand_mean = df["y"].mean()
    ss_tot = ((df["y"] - grand_mean) ** 2).sum()
    group_means = df.groupby("g")["y"].transform("mean")
    ss_between = ((group_means - grand_mean) ** 2).sum()
    return ss_between / ss_tot if ss_tot > 0 else 0.0


def permutation_p_batch(y, labels, n=N_PERM_BATCH):
    obs = groupmean_r2(y, labels)
    labels = np.asarray(labels)
    null = np.empty(n)
    for i in range(n):
        null[i] = groupmean_r2(y, RNG.permutation(labels))
    return obs, (np.sum(null >= obs) + 1) / (n + 1)


def compute_purity_proxy(rna_common):
    tumor = [g for g in TUMOR_MARKERS if g in rna_common.columns]
    stroma = [g for g in STROMAL_MARKERS if g in rna_common.columns]
    print(f"[Check C] purity proxy using {len(tumor)}/{len(TUMOR_MARKERS)} tumor markers, "
          f"{len(stroma)}/{len(STROMAL_MARKERS)} stromal markers found in RNA matrix")
    t = rna_common[tumor]
    s = rna_common[stroma]
    t_z = (t - t.mean()) / t.std().replace(0, 1)
    s_z = (s - s.mean()) / s.std().replace(0, 1)
    return (t_z.mean(axis=1) - s_z.mean(axis=1)).values  # high = tumor-rich, low = stroma-rich


def main():
    rna = load_rna_matrix()
    protein = load_protein_matrix()
    wsi = load_wsi_embeddings()
    results = pd.read_csv(OUT_DIR / "residual_results.csv")

    # ---- Check A: mRNA-strength confound ----
    corr = results[["r2_rna", "incremental_r2"]].corr(method="spearman").iloc[0, 1]
    print(f"\n[Check A] Spearman(incremental_r2, r2_rna) = {corr:.3f} "
          f"(near 0 expected; strongly negative would mean WSI is just "
          f"'rescuing' genes with bad mRNA baselines -- consistent with a "
          f"generic confound rather than gene-specific signal)")

    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    common_genes = sorted(set(rna.columns) & set(protein.columns))
    rna_c = rna.loc[common, common_genes]
    protein_c = protein.loc[common, common_genes]
    wsi_raw = wsi.loc[common].values

    # svd_solver='full' is load-bearing: at this matrix shape sklearn's auto
    # solver is randomized SVD with random_state=None, which made this script
    # non-reproducible (GBM: 699/700/704/701 significant across four identical
    # runs, 90.4% set overlap). 'full' is exact and deterministic.
    pca = PCA(n_components=min(N_PCS, len(common) - 1, wsi_raw.shape[1]),
              svd_solver="full")
    wsi_pcs = pca.fit_transform(wsi_raw)
    print(f"WSI PCA (recomputed, same as main run): {wsi_pcs.shape[1]} comps, "
          f"{pca.explained_variance_ratio_.sum():.1%} var")

    # ---- Check B: TMT batch/plex confound ----
    plex_of_case = detect_plex_blocks()
    plex_labels = np.array([plex_of_case.get(c, -1) for c in common])
    n_missing = (plex_labels == -1).sum()
    if n_missing:
        print(f"[Check B] WARNING: {n_missing}/{len(common)} common patients have no plex label, dropping them")
    keep = plex_labels != -1
    for i in range(min(5, wsi_pcs.shape[1])):
        obs_r2, p = permutation_p_batch(wsi_pcs[keep, i], plex_labels[keep])
        print(f"  WSI-PC{i+1} ~ plex_batch:  R2={obs_r2:.3f}  perm-p={p:.4f}")

    # ---- Check C: tumor-purity confound ----
    purity = compute_purity_proxy(rna_c)
    for i in range(min(5, wsi_pcs.shape[1])):
        r = np.corrcoef(purity, wsi_pcs[:, i])[0, 1]
        print(f"  purity_proxy vs WSI-PC{i+1}: r={r:.3f}")

    print("\n[Check C] recomputing incremental R2 with purity forced into the "
          "baseline model (protein ~ mRNA+purity  vs  protein ~ mRNA+purity+WSI)")
    print("this is point-estimate only (no permutation) -- fast, just checking "
          "whether the effect size survives, not re-deriving FDR")
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
        X_base = np.hstack([x_rna, x_purity])
        X_full = np.hstack([x_rna, x_purity, x_wsi])
        r2_base = cv_r2(X_base, y_v, CV_FOLDS)
        r2_full = cv_r2(X_full, y_v, CV_FOLDS)
        rows.append({"gene": g, "incr_r2_purity_adjusted": r2_full - r2_base})
        if (j + 1) % 2000 == 0:
            print(f"  {j+1}/{len(common_genes)} genes purity-adjusted...")

    adj = pd.DataFrame(rows)
    merged = results.merge(adj, on="gene", how="inner")
    orig_med = merged["incremental_r2"].median()
    adj_med = merged["incr_r2_purity_adjusted"].median()
    print(f"\n[Check C] median incremental_r2 (original, no purity control):     {orig_med:.4f}")
    print(f"[Check C] median incremental_r2 (purity-adjusted):                  {adj_med:.4f}")

    sig = merged[(merged["fdr"] < 0.05) & (merged["incremental_r2"] > 0)]
    survive = (sig["incr_r2_purity_adjusted"] > 0.5 * sig["incremental_r2"])
    print(f"[Check C] of the {len(sig)} originally-significant genes, "
          f"{survive.sum()} ({survive.mean():.0%}) retain >=50% of their "
          f"incremental_r2 after adjusting for the purity proxy")

    merged.to_csv(OUT_DIR / "confound_check_results.csv", index=False)
    print(f"\nfull table: {OUT_DIR / 'confound_check_results.csv'}")

    print("\ntop 30 significant genes (for manual pathway-coherence check "
          "against plan.md's post-transcriptional/degradation/secretion/cell-cycle expectation):")
    print(results[(results.fdr < 0.05) & (results.incremental_r2 > 0)]
          .sort_values("incremental_r2", ascending=False)
          .head(30)[["gene", "incremental_r2", "fdr"]].to_string(index=False))


if __name__ == "__main__":
    main()
