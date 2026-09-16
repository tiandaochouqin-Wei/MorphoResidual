#!/usr/bin/env python3
"""Dense study-overview figure (multi-panel workflow + real-data panels + ranked bars).
Built around OUR study's actual components; uses OUR data only."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import mrstyle as S

coh = ["UCEC", "LUAD", "CCRCC", "PDAC", "GBM"]
organ = {"UCEC": "Uterus", "LUAD": "Lung", "CCRCC": "Kidney", "PDAC": "Pancreas", "GBM": "Brain"}
sig = np.array([2566, 2310, 2191, 1572, 702])
pct = np.array([24.4, 21.5, 22.7, 15.9, 6.5])
med_incr = np.array([0.118, 0.070, 0.086, 0.054, 0.058])
col = [S.ORG[c] for c in coh]

fig = plt.figure(figsize=(7.4, 8.2))
gs = fig.add_gridspec(3, 3, height_ratios=[0.62, 0.95, 1.05], hspace=0.32, wspace=0.16,
                      left=0.085, right=0.975, top=0.955, bottom=0.075)


def imgpanel(ax, path, letter, title):
    ax.imshow(mpimg.imread(path)); ax.axis("off")
    ax.set_title(f"{letter}   {title}", loc="left", fontsize=8, fontweight="bold", pad=3)


# a: workflow schematic (full width)
axA = fig.add_subplot(gs[0, :])
imgpanel(axA, "Fig1_concept.png", "a", "Study workflow: H&E reads the beyond-mRNA protein residual")

# b,c,d: OUR "assess | interpret images | interpret omics" functional row (our data)
imgpanel(fig.add_subplot(gs[1, 0]), "Fig10_morph_recovers.png", "b", "Assess: H&E recovers the residual")
imgpanel(fig.add_subplot(gs[1, 1]), "Fig11_tiles_hilo.png", "c", "Interpret images: high vs low residual")
imgpanel(fig.add_subplot(gs[1, 2]), "Fig3_enrichment.png", "d", "Interpret omics: pathway enrichment")

# e: ranked bars (full width)
axE = fig.add_subplot(gs[2, :])
x = np.arange(5)
axE.bar(x, sig, color=col, width=0.64, edgecolor="none")
for i, (s, p) in enumerate(zip(sig, pct)):
    axE.text(i, s + 30, f"{s:,}", ha="center", va="bottom", fontsize=7.2, color="#B22222", fontweight="bold")
    axE.text(i, s + 175, f"{p:.1f}%", ha="center", va="bottom", fontsize=6, color="#777")
axE.set_xticks(x); axE.set_xticklabels([f"{c}\n({organ[c]})" for c in coh], fontsize=7)
axE.set_ylabel("Morphology-predictable residual\nproteins  (FDR < 0.05)", fontsize=7.5)
axE.set_ylim(0, 3050)
axE.set_title("e   Reproducible across five organs (red = count, grey = % of tested; line = median incremental $R^2$)",
              loc="left", fontsize=8, fontweight="bold")
axE2 = axE.twinx()
axE2.plot(x, med_incr, "o-", color=S.INK, lw=1.2, ms=5, zorder=5)
for i, m in enumerate(med_incr):
    axE2.text(i + 0.03, m + 0.004, f"{m:.2f}", fontsize=5.8, color=S.INK, va="bottom")
axE2.set_ylabel(r"median incr $R^2$", fontsize=7, color=S.INK); axE2.set_ylim(0, 0.16)
axE2.tick_params(axis="y", labelsize=6); axE2.spines["top"].set_visible(False)

S.save_pub(fig, "Fig_overview")
print("wrote Fig_overview (dense multi-panel, our data)")
