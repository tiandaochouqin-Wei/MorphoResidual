#!/usr/bin/env python3
"""Real-data figures 2: (1) morphology-only vs measured residual pathway score scatter;
(2) representative H&E tiles from high- vs low-residual patients. Uses figdata/."""
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy import stats as st
import mrstyle as S

COH = S.COH

# ---- load per-patient scores, pooled across cohorts ----
frames = []
for c in COH:
    d = pd.read_csv(f"figdata/scores_{c}.csv", index_col=0)
    d["cohort"] = c
    frames.append(d)
sc = pd.concat(frames)

# ================= FIG 10 — morphology recovers the measured residual programme =====
fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.9))
for ax, fam, lab in [(a, "translation", "translation"), (b, "secretion_ER", "ER-secretion")]:
    mcol, pcol = f"{fam}_meas", f"{fam}_morph"
    sub = sc.dropna(subset=[mcol, pcol])
    for c in COH:
        m = sub["cohort"].eq(c).values
        ax.scatter(sub[mcol].values[m], sub[pcol].values[m], s=9, c=S.ORG[c],
                   edgecolors="none", alpha=0.8, label=c if ax is a else None)
    r, p = st.pearsonr(sub[mcol], sub[pcol])
    # overall least-squares trend
    z = np.polyfit(sub[mcol], sub[pcol], 1)
    xs = np.array([sub[mcol].min(), sub[mcol].max()])
    ax.plot(xs, np.polyval(z, xs), lw=1.1, color=S.INK)
    ax.axhline(0, lw=0.4, color="#DDD"); ax.axvline(0, lw=0.4, color="#DDD")
    ax.text(0.04, 0.95, f"$r={r:.2f}$\n$P={p:.1e}$", transform=ax.transAxes,
            fontsize=6.4, va="top", fontweight="bold")
    ax.set_xlabel(f"measured {lab} residual score\n(from proteomics)")
    ax.set_ylabel(f"morphology-only score\n(from H&E)")
    ax.set_title(f"{'a' if ax is a else 'b'}  {lab}", loc="left", fontsize=8, fontweight="bold")
a.legend(loc="lower right", fontsize=5.4, handletextpad=0.2, labelspacing=0.25,
         markerscale=1.1, ncol=1)
fig.suptitle("A morphology-only score recovers the measured residual programme",
             fontsize=8.5, fontweight="bold", x=0.02, ha="left", y=1.02)
fig.tight_layout()
S.save_pub(fig, "Fig10_morph_recovers")
plt.close(fig)
print("wrote Fig10_morph_recovers")

# ================= FIG 11 — representative H&E tiles, high vs low residual ==========
# curated cellular tiles (drop blank/blood artefacts by index, from visual QA)
hi = [f"figdata/tile_hi_{i}_{j}.png" for i, j in [(1, 0), (2, 0), (2, 1), (3, 1)]]
lo = [f"figdata/tile_lo_{i}_{j}.png" for i, j in [(0, 1), (1, 1), (3, 0), (3, 1)]]
hi = [f for f in hi if glob.glob(f)]
lo = [f for f in lo if glob.glob(f)]

fig, axes = plt.subplots(2, 4, figsize=(5.4, 3.0))
for j, f in enumerate(hi[:4]):
    axes[0, j].imshow(mpimg.imread(f)); axes[0, j].axis("off")
for j, f in enumerate(lo[:4]):
    axes[1, j].imshow(mpimg.imread(f)); axes[1, j].axis("off")
axes[0, 0].set_ylabel("")
fig.text(0.01, 0.72, "highest\nresidual", fontsize=6.6, va="center", ha="left",
         fontweight="bold", color=S.ORG["CCRCC"])
fig.text(0.01, 0.28, "lowest\nresidual", fontsize=6.6, va="center", ha="left",
         fontweight="bold", color=S.GREY)
fig.suptitle("Representative H&E (CCRCC): highest vs lowest translation-residual patients",
             fontsize=7.6, fontweight="bold", x=0.02, ha="left", y=1.0)
fig.subplots_adjust(left=0.14, right=0.99, top=0.9, bottom=0.02, wspace=0.06, hspace=0.06)
S.save_pub(fig, "Fig11_tiles_hilo", tiff=False)
plt.close(fig)
print("wrote Fig11_tiles_hilo (curated)")
