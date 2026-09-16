#!/usr/bin/env python3
"""MorphoResidual Figure 2 — incremental-R2 distribution with capacity controls as
reference lines (CCRCC). Single panel: the three capacity-control medians (real WSI PCs,
20-d Gaussian noise, random 20-d projection; all on the FDR-significant gene set) are
drawn as labelled vertical reference lines on the all-gene histogram."""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke

import mrstyle as S  # Arial, editable-text SVG, 600-dpi TIFF (submission grade)

# CCRCC incremental-R2 histogram (server, current deterministic run; n=9635)
lo = np.arange(-0.70, 0.80, 0.05)
cnt = np.array([84,112,149,172,257,317,366,474,573,663,807,965,974,916,699,
                533,374,263,140,88,54,24,9,5,3,1,0,1,0,0])
median = -0.1375
frac_pos = 0.228
# capacity control on the significant gene set (median incremental R^2)
cap = {  # key: (value, colour, linestyle, linewidth)
    "real":   (0.086,  S.FAM["translation"], "-",  1.5),
    "gauss":  (-0.293, S.GREY,               "--", 1.3),
    "rproj":  (-0.083, S.LGREY,              ":",  2.0),
}

fig, a = plt.subplots(figsize=(4.8, 3.0))

# --- distribution (all genes) ---
col = ["#D55E00" if x >= 0 else "#B9C0C7" for x in lo]
a.bar(lo, cnt, width=0.048, color=col, align="edge", edgecolor="none", zorder=2)
a.axvline(0, color="k", lw=0.7, zorder=3)
a.axvline(median, color="#333", lw=0.9, ls="--", zorder=3)
a.text(median - 0.015, 1010, f"median {median:.2f}".replace("-", "−"),
       ha="right", va="bottom", fontsize=6.2, color="#333", zorder=6,
       bbox=dict(fc="white", ec="none", pad=0.8))
a.text(0.30, 780, "22.8% of genes positive;\n2,191 (22.7%)\nFDR-significant", fontsize=6.6, color="#D55E00", va="top")
a.text(-0.695, 300, "612 genes\n< −0.70\nnot shown", fontsize=5.4, color="#8C8C8C", va="bottom", ha="left", linespacing=1.1)

# --- capacity controls as reference lines (significant genes; medians only) ---
YTOP = 1400
# each line stops just under its own label so that no line crosses another label
specs = {  # key: (line top, label text, label x-offset, ha)
    "real":  (1000, "real WSI PCs +0.09",        +0.012, "left"),
    "gauss": (880,  "Gaussian noise −0.29",   -0.010, "right"),
    "rproj": (1130, "random projection −0.08", -0.012, "right"),
}
for k, (v, c, ls, lw) in cap.items():
    ytop, lab, dx, ha = specs[k]
    ln = a.plot([v, v], [0, ytop], color=c, ls=ls, lw=lw, zorder=4, solid_capstyle="butt")[0]
    if k == "rproj":  # light dotted line needs a dark halo to read over the grey bars
        ln.set_path_effects([withStroke(linewidth=lw + 1.8, foreground=S.INK)])
    tcol = S.INK if k == "rproj" else c
    a.text(v + dx, ytop + 20, lab, ha=ha, va="bottom", fontsize=6.2, color=tcol,
           fontweight="bold" if k == "real" else "normal", zorder=6,
           bbox=dict(fc="white", ec="none", pad=0.8))   # white backing so no line crosses text
a.text(0.74, YTOP, "reference lines: median incremental $R^2$\non the FDR-significant set,\n"
       "per feature set (capacity control)", ha="right", va="top", fontsize=5.8,
       color=S.GREY, linespacing=1.05)

a.set_xlabel(r"Incremental $R^2$ of morphology over mRNA")
a.set_ylabel("Genes")
a.set_xlim(-0.72, 0.75)
a.set_ylim(0, YTOP)
a.set_yticks([0, 200, 400, 600, 800, 1000, 1200])
a.set_title("Most genes show no gain; the positive tail is a specific, non-capacity effect",
            loc="left", fontsize=7.6, fontweight="bold")

fig.tight_layout()
S.save_pub(fig, "Fig2_distribution")
print("wrote Fig2_distribution.{svg,pdf,png,tiff}")
