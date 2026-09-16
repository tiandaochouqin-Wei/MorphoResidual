#!/usr/bin/env python3
"""Motivating panel (Clark Fig3A / Buccitelli Fig2 form): distribution of the per-protein
mRNA->protein predictive R^2 across all five cohorts. Most proteins are poorly explained by
mRNA -> a large, pervasive 'beyond-mRNA' residual, which is what morphology then predicts."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import mrstyle as S

BASE = "D:/claude/pajinsen/MorphoResidual/pilot_ccrcc/results_server"
PATHS = {"CCRCC": f"{BASE}/residual_results_tumoronly.csv",
         "LUAD": f"{BASE}/luad/results/residual_results_tumoronly.csv",
         "UCEC": f"{BASE}/ucec/results/residual_results_tumoronly.csv",
         "GBM": f"{BASE}/gbm/results/residual_results_tumoronly.csv",
         "PDAC": f"{BASE}/pdac/results/residual_results_tumoronly.csv"}

r2 = {c: pd.read_csv(p)["r2_rna"].values for c, p in PATHS.items()}
pooled = np.concatenate(list(r2.values()))
pooled = pooled[np.isfinite(pooled)]
med = np.median(pooled)
frac_low = 100 * np.mean(pooled < 0.25)
frac_neg = 100 * np.mean(pooled <= 0)

fig, ax = plt.subplots(figsize=(3.5, 2.8))
ax.hist(np.clip(pooled, -0.1, 1.0), bins=62, range=(-0.1, 1.0),
        color="#9ecae1", edgecolor="white", linewidth=0.3)
ymax = ax.get_ylim()[1]
# shade the poorly-explained region
ax.axvspan(-0.1, 0.25, color=S.ORG["PDAC"], alpha=0.10, zorder=0)
ax.axvline(med, color=S.INK, lw=1.2, ls="--")
ax.text(med + 0.02, ymax * 0.92, f"median $R^2$ = {med:.2f}", fontsize=6.6, color=S.INK)
ax.text(0.07, ymax * 0.62, f"{frac_low:.0f}% of proteins\n$R^2$ < 0.25", fontsize=6.4,
        color=S.ORG["PDAC"], ha="left")
ax.set_xlabel(r"mRNA$\rightarrow$protein $R^2$  (per protein, cross-validated)", fontsize=7.5)
ax.set_ylabel("proteins", fontsize=7.5)
ax.set_xlim(-0.1, 1.0)
ax.set_title("Most protein variation is not explained by mRNA",
             fontsize=8, fontweight="bold", loc="left")
# per-cohort median ticks under the axis
for c in S.COH:
    m = np.median(r2[c][np.isfinite(r2[c])])
    ax.plot([m, m], [-ymax * 0.045, -ymax * 0.02], color=S.ORG[c], lw=1.4, clip_on=False)
ax.text(0.5, -0.24, "coloured ticks = per-cohort medians", transform=ax.transAxes,
        fontsize=5.6, color=S.GREY, ha="center")
S.save_pub(fig, "Fig_residual_motivation")
print(f"wrote Fig_residual_motivation  (n={len(pooled)}, median r2_rna={med:.3f}, "
      f"%<0.25={frac_low:.0f}, %<=0={frac_neg:.0f})")
print("per-cohort median r2_rna:", {c: round(float(np.median(r2[c][np.isfinite(r2[c])])), 3) for c in S.COH})
