#!/usr/bin/env python3
"""vignette_plot.py (LOCAL) — single-protein vignette from per_gene_vignette.py output.
Per protein: (left) protein vs mRNA (mRNA barely predicts it -> a large residual);
(right) the measured residual recovered by the morphology-only prediction."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as st
import mrstyle as S

CSV = "figdata/vignette_data.csv"
if not os.path.exists(CSV):
    raise SystemExit(f"missing {CSV} -- run per_gene_vignette.py on the server and WinSCP it back")

df = pd.read_csv(CSV)
genes = list(dict.fromkeys(df["gene"]))          # preserve order
col = [S.ORG["CCRCC"], S.ORG["UCEC"]]

fig, axes = plt.subplots(len(genes), 2, figsize=(5.4, 2.7 * len(genes)))
axes = np.atleast_2d(axes)


def r2(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return 1 - np.sum((a - b) ** 2) / np.sum((a - a.mean()) ** 2)


for gi, g in enumerate(genes):
    d = df[df["gene"] == g]
    c = col[gi % len(col)]
    # left: protein vs mRNA
    ax = axes[gi, 0]
    ax.scatter(d["mrna"], d["protein"], s=12, color=c, alpha=0.7, edgecolors="none")
    b1, b0 = np.polyfit(d["mrna"], d["protein"], 1)
    xr = np.array([d["mrna"].min(), d["mrna"].max()])
    ax.plot(xr, b1 * xr + b0, color=S.GREY, lw=1.0, ls="--")
    rr = st.pearsonr(d["mrna"], d["protein"])[0]
    ax.text(0.05, 0.93, f"{g}\nmRNA$\\rightarrow$protein $r$={rr:.2f}", transform=ax.transAxes,
            fontsize=6.8, va="top", fontweight="bold", color=c)
    ax.set_xlabel("mRNA (log$_2$ TPM)", fontsize=7); ax.set_ylabel("protein (measured)", fontsize=7)
    # right: residual recovered by morphology
    ax = axes[gi, 1]
    ax.scatter(d["residual"], d["morph_pred"], s=12, color=c, alpha=0.7, edgecolors="none")
    lo = min(d["residual"].min(), d["morph_pred"].min()); hi = max(d["residual"].max(), d["morph_pred"].max())
    ax.plot([lo, hi], [lo, hi], color=S.GREY, lw=0.8, ls=":")
    b1, b0 = np.polyfit(d["residual"], d["morph_pred"], 1)
    xr = np.array([d["residual"].min(), d["residual"].max()])
    ax.plot(xr, b1 * xr + b0, color=S.INK, lw=1.0)
    incr = r2(d["residual"], d["morph_pred"])
    rr = st.pearsonr(d["residual"], d["morph_pred"])[0]
    ax.text(0.05, 0.93, f"morphology recovers residual\n$r$={rr:.2f}  (morphology $R^2$={incr:.2f})",
            transform=ax.transAxes, fontsize=6.4, va="top", color=S.INK)
    ax.set_xlabel("measured residual\n(protein $-$ f(mRNA))", fontsize=7)
    ax.set_ylabel("morphology-predicted\nresidual", fontsize=7)

for i, letter in enumerate("abcdefgh"[:axes.size]):
    axes.flat[i].text(-0.22, 1.02, letter, transform=axes.flat[i].transAxes, fontsize=10,
                      fontweight="bold", va="bottom", ha="right")
fig.suptitle("Single-protein vignettes: mRNA misses it, morphology reads the residual",
             fontsize=8.6, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.96))
S.save_pub(fig, "Fig_vignette")
print(f"wrote Fig_vignette  (genes: {genes})")
