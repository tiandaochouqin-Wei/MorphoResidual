#!/usr/bin/env python3
"""Supplementary figure for the pan-organ gene-level core (T1-3, NC_ROADMAP.md).
Two panels: (a) k-of-5 histogram, observed vs Poisson-binomial and permutation
nulls; (b) the 50 machine-annotated k>=4 genes as a heat-strip of per-cohort
incremental R^2, grouped by machine. Not a headline figure -- Supplement only."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import mrstyle as S

COH = S.COH
COHL = [c.lower() for c in COH]
DD = "figdata"

hist = pd.read_csv(f"{DD}/pan_organ_k5hist.csv")
core = pd.read_csv(f"{DD}/pan_organ_core.csv", index_col=0)

fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.5, left=0.07, right=0.97, top=0.95, bottom=0.09)


def lab(ax, letter, x=-0.045, y=1.05):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


# ---- a: k-of-5 histogram ----
axa = fig.add_subplot(gs[0]); lab(axa, "a")
x = hist["k"].values
w = 0.36
axa.bar(x - w / 2, hist["observed"], w, color=S.INK, label="observed", zorder=3)
axa.bar(x + w / 2, hist["poisson_binomial_expected"], w, color=S.LGREY, edgecolor=S.GREY, lw=0.6,
        label="Poisson-binomial null", zorder=3)
for xi, obs, exp in zip(x, hist["observed"], hist["poisson_binomial_expected"]):
    if xi >= 4:
        # each value label centred on its own bar: at 7 pt a tick-centred label
        # straddles the bar boundary and the grey one lands on the dark bar
        axa.text(xi - w / 2, obs + 6, f"{obs:.0f}", ha="center", fontsize=7.0, color=S.INK, fontweight="bold")
        axa.text(xi + w / 2, exp + 6, f"{exp:.0f}", ha="center", fontsize=7.0, color=S.GREY)
axa.set_yscale("symlog", linthresh=10)
axa.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
axa.set_xticks(x)
axa.set_xlabel("number of cohorts (of 5) in which a gene is FDR-significant\nwith positive incremental R\u00b2 (k)")
axa.set_ylabel("genes (of 8,309 tested\nin all five cohorts)")
axa.legend(loc="upper right", fontsize=7.0, frameon=False)
axa.set_title("Most genes are significant in 0-3 cohorts; k$\\geq$4 exceeds both nulls by ~3.4$\\times$",
              fontsize=7.6, fontweight="bold", loc="left")
axa.text(0.98, 0.7, "independent within-cohort\npermutation null (N=1,000):\nmean k$\\geq$4 = 50.5, p<0.001",
          transform=axa.transAxes, ha="right", va="top", fontsize=7.0, color=S.GREY, linespacing=1.3)

# ---- b: per-machine median incremental R2, one row per machine (individual genes in SuppTable~11) ----
axb = fig.add_subplot(gs[1]); lab(axb, "b")
mach_nice = {"translocon_OST_SPC": "Translocon / OST / SPC", "spliceosome_mRNA_export": "Spliceosome / mRNA export",
             "ribosome_translation": "Ribosome / translation", "cohesin": "Cohesin", "other": "Unclassified"}
order = ["translocon_OST_SPC", "spliceosome_mRNA_export", "ribosome_translation", "cohesin", "other"]
summary = core.groupby("machine")[COHL].median().reindex(order)
counts = core["machine"].value_counts().reindex(order)
mat = summary.values.astype(float)
im = axb.imshow(mat, cmap="RdBu_r", vmin=-0.15, vmax=0.15, aspect="auto")
axb.set_xticks(range(5)); axb.set_xticklabels(COH, fontsize=7.0)
axb.set_yticks(range(len(order)))
axb.set_yticklabels([f"{mach_nice[m]}  (n={counts[m]})" for m in order], fontsize=7.0)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = mat[i, j]
        axb.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7.0,
                  color="white" if abs(v) > 0.09 else "#333")
cb = fig.colorbar(im, ax=axb, fraction=0.035, pad=0.03)
cb.set_label("median incremental R\u00b2\n(morphology over own mRNA, UNI)", fontsize=7.0)
cb.ax.tick_params(labelsize=7.0)
axb.set_title("Median incremental R\u00b2 of the k$\\geq$4 core, by obligate-complex machine "
              "(170 genes total, 50 classifiable; see Supplementary Table)",
              fontsize=7.6, fontweight="bold", loc="left")
for s in axb.spines.values():
    s.set_visible(False)
axb.tick_params(length=0)

S.save_pub(fig, "Fig_pan_organ_core", tiff=False)
print("wrote Fig_pan_organ_core")
