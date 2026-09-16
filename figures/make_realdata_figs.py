#!/usr/bin/env python3
"""Real-data figures for MorphoResidual: per-gene UNI-vs-Phikon scatter (Fig 7) and
a t-SNE of slide embeddings coloured by organ (Fig 9). Uses figdata/ exported from
the server. Style: mrstyle (Arial, editable SVG, 600-dpi TIFF)."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as st
from sklearn.manifold import TSNE
import mrstyle as S

COH = S.COH
x = np.arange(5)

# ==================================================== FIG 7 — second FM, REAL scatter
r_incr = np.array([0.808, 0.780, 0.845, 0.669, 0.854])
recall = np.array([0.737, 0.587, 0.742, 0.380, 0.642])

fig = plt.figure(figsize=(6.8, 2.7))
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.42)
a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1])

# a: hero hexbin — per-gene incremental R2, UNI vs Phikon (CCRCC)
d = pd.read_csv("figdata/merged_incr_ccrcc.csv")
xu, yp = d["incremental_r2_uni"].values, d["incremental_r2_phi"].values
lim = (-0.35, 0.75)
hb = a.hexbin(xu, yp, gridsize=45, extent=(*lim, *lim), cmap="Blues",
              bins="log", mincnt=1, linewidths=0)
a.plot(lim, lim, ls="--", lw=0.8, color=S.GREY)
a.axhline(0, lw=0.5, color="#CCC"); a.axvline(0, lw=0.5, color="#CCC")
r = np.corrcoef(xu, yp)[0, 1]
a.text(0.05, 0.93, f"$r = {r:.2f}$\n$n = {len(d):,}$ genes", transform=a.transAxes,
       fontsize=6.6, va="top", fontweight="bold")
a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel(r"UNI incremental $R^2$"); a.set_ylabel(r"Phikon incremental $R^2$")
a.set_title("a  Per-gene agreement (CCRCC)", loc="left", fontsize=8, fontweight="bold")
cb = fig.colorbar(hb, ax=a, fraction=0.046, pad=0.03); cb.set_label("genes (log)", fontsize=5.6)
cb.ax.tick_params(labelsize=5.4)

# b: cross-cohort agreement
wa = 0.36
b.bar(x - wa / 2, r_incr, wa, color="#0072B2", label=r"incr-$R^2$ corr $r$")
b.bar(x + wa / 2, recall, wa, color="#B9D4E8", label="UNI-significant recall")
for i, v in enumerate(r_incr):
    b.text(i - wa / 2, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=5.2)
b.set_xticks(x); b.set_xticklabels(COH, rotation=45, ha="right")
b.set_ylabel("Agreement (UNI vs Phikon)"); b.set_ylim(0, 1.15)
b.legend(loc="upper right", fontsize=5.6, handlelength=1.0)
b.set_title("b  Consistent across cohorts", loc="left", fontsize=8, fontweight="bold")
fig.suptitle("The signal is not specific to one foundation model",
             fontsize=8.5, fontweight="bold", x=0.02, ha="left", y=1.02)
S.save_pub(fig, "Fig7_second_fm")
plt.close(fig)
print("wrote Fig7_second_fm (real per-gene scatter)")

# ==================================================== FIG 9 — embedding t-SNE by organ
emb = pd.read_csv("figdata/emb_pca50.csv")
Z = TSNE(n_components=2, perplexity=30, init="pca", learning_rate="auto",
         random_state=0).fit_transform(emb.filter(like="pc").values)
fig, ax = plt.subplots(figsize=(3.6, 3.2))
name2coh = {"ccrcc": "CCRCC", "luad": "LUAD", "ucec": "UCEC", "gbm": "GBM", "pdac": "PDAC"}
for c in COH:
    key = c.lower()
    m = emb["cohort"].str.lower().eq(key).values
    ax.scatter(Z[m, 0], Z[m, 1], s=5, c=S.ORG[c], edgecolors="none", alpha=0.75,
               label=f"{c} ({S.ORGAN[c]})")
ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("t-SNE 1", fontsize=6.5); ax.set_ylabel("t-SNE 2", fontsize=6.5)
ax.set_title("UNI slide embeddings separate by organ",
             loc="left", fontsize=8, fontweight="bold")
ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=6, handletextpad=0.3,
          borderpad=0.3, labelspacing=0.35)
for s in ("top", "right", "left", "bottom"):
    ax.spines[s].set_visible(True); ax.spines[s].set_linewidth(0.6)
S.save_pub(fig, "Fig9_embedding_tsne")
plt.close(fig)
print(f"wrote Fig9_embedding_tsne ({len(emb)} slides)")
