#!/usr/bin/env python3
"""Fig 12 — interpretable nuclear-morphometry correlates of the residual programmes.
Spearman rho of 6 hand-crafted H&E nuclear features vs the measured translation and
ER-secretion residual scores, per cohort (server-verified pathomics_link output)."""
import numpy as np
import matplotlib.pyplot as plt
import mrstyle as S

COH = S.COH
FEAT = ["nuclear density", "nucleus count", "mean nucleus area",
        "nucleus area SD", "chromatin OD", "eosin fraction"]

# rows = FEAT, cols = COH (CCRCC LUAD UCEC GBM PDAC)
transl = np.array([
    [ 0.18, -0.09,  0.05, -0.03,  0.23],
    [ 0.13,  0.38, -0.43, -0.08,  0.13],
    [ 0.08, -0.33,  0.34, -0.01, -0.04],
    [ 0.13, -0.28,  0.30, -0.06, -0.02],
    [ 0.05, -0.21,  0.33, -0.04,  0.13],
    [ 0.15, -0.05, -0.15, -0.04,  0.08]])
secr = np.array([
    [-0.08, -0.01, -0.25,  0.04,  0.26],
    [ 0.35,  0.24,  0.38, -0.10,  0.23],
    [-0.27, -0.17, -0.42,  0.02, -0.11],
    [-0.19, -0.18, -0.45,  0.03, -0.02],
    [-0.22,  0.03, -0.37, -0.03,  0.08],
    [ 0.27, -0.21,  0.21, -0.03,  0.03]])

fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), gridspec_kw={"wspace": 0.12})
for ax, M, ttl in [(axes[0], transl, "a  translation residual"),
                   (axes[1], secr, "b  ER-secretion residual")]:
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_xticks(range(5)); ax.set_xticklabels(COH, rotation=45, ha="right", fontsize=6.2)
    if ax is axes[0]:
        ax.set_yticks(range(len(FEAT))); ax.set_yticklabels(FEAT, fontsize=6.2)
    else:
        ax.set_yticks([])
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=5.4,
                    color="white" if abs(v) > 0.28 else "#333")
    ax.set_title(ttl, loc="left", fontsize=8, fontweight="bold")
    for s in ax.spines.values():
        s.set_visible(False)
cb = fig.colorbar(im, ax=axes, fraction=0.028, pad=0.02)
cb.set_label(r"Spearman $\rho$", fontsize=6.2); cb.ax.tick_params(labelsize=5.6)
fig.suptitle("Interpretable nuclear correlates of the residual programmes",
             fontsize=8.5, fontweight="bold", x=0.02, ha="left", y=1.03)
S.save_pub(fig, "Fig12_pathomics")
print("wrote Fig12_pathomics")
