#!/usr/bin/env python3
"""make_txome_fig.py (LOCAL) -- figure for the transcriptome-wide baseline analysis (P0-1).

Three panels, deliberately none of them bar charts:
  a  cohort-level scatter: morphology's noise-adjusted advantage against the own-mRNA baseline (x)
     vs against own mRNA + 20 RNA PCs (y). Diagonal = 100% retained; guide lines at 50% / 25%.
     Each cohort is one point labelled with its retained fraction.
  b  per-gene distribution of the paired excess (real morphology minus same-size noise) under the
     transcriptome-wide baseline, one violin per cohort, with the % of genes > 0 annotated.
     Needs figdata/transcriptome_baseline_<cohort>.csv (WinSCP from the server); the panel is
     drawn as a placeholder note if they are absent.
  c  dumbbell of median canonical correlation with the 20 RNA PCs: morphology block (filled) vs
     same-size noise block (open), per cohort -- the collinearity check that makes (a) an upper bound.

Inputs: figdata/transcriptome_baseline_summary.csv (cohort-level, from the server printout) and,
optionally, figdata/transcriptome_baseline_<cohort>.csv (per-gene)."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import mrstyle as S

DD = "figdata"
COH = ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]
d = pd.read_csv(f"{DD}/transcriptome_baseline_summary.csv").set_index("cohort").loc[COH]

fig = plt.figure(figsize=(7.2, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 0.95], wspace=0.55, left=0.07, right=0.985, top=0.86, bottom=0.25)

# ---- a: retained advantage, cohort-level ----
a = fig.add_subplot(gs[0])
lim = 0.48
a.plot([0, lim], [0, lim], color=S.GREY, lw=0.8, ls="--")
a.plot([0, lim], [0, 0.5 * lim], color=S.LGREY, lw=0.6, ls=":")
a.plot([0, lim], [0, 0.25 * lim], color=S.LGREY, lw=0.6, ls=":")
a.text(lim - 0.01, lim - 0.01, "100%", fontsize=5.4, color=S.GREY, ha="right", va="top")
a.text(0.13, 0.5 * 0.13 + 0.008, "50%", fontsize=5.4, color=S.GREY, ha="left", va="bottom")
a.text(0.13, 0.25 * 0.13 - 0.008, "25%", fontsize=5.4, color=S.GREY, ha="left", va="top")
for c in COH:
    r = d.loc[c]
    a.scatter(r.adv_own, r.adv_pcs20, s=34, color=S.ORG[c], edgecolors=S.INK, linewidths=0.4, zorder=4)
    off = {"UCEC": (0.014, 0.0, "left", "center"), "CCRCC": (0.014, 0.0, "left", "center"),
           "GBM": (0.014, 0.018, "left", "center"), "LUAD": (-0.014, -0.004, "right", "center"),
           "PDAC": (0.0, -0.028, "center", "top")}[c]
    a.text(r.adv_own + off[0], r.adv_pcs20 + off[1], f"{c} {r.retains_pct:.0f}%",
           fontsize=5.6, color=S.ORG[c], ha=off[2], va=off[3])
a.set_xlim(0, lim); a.set_ylim(-0.02, lim)
a.set_xlabel("advantage over noise,\nbaseline = own mRNA", fontsize=7)
a.set_ylabel("advantage over noise,\nbaseline = own mRNA + 20 RNA PCs", fontsize=7)
a.set_title("Retained against the transcriptome", fontsize=7.6, fontweight="bold", loc="left")
a.text(-0.32, 1.06, "a", transform=a.transAxes, fontsize=10, fontweight="bold")

# ---- b: per-gene excess distributions (if the per-gene files are present) ----
b = fig.add_subplot(gs[1])
per = {c: f"{DD}/transcriptome_baseline_{c.lower()}.csv" for c in COH}
have = all(os.path.exists(p) for p in per.values())
if have:
    data, pos = [], []
    for c in COH:
        g = pd.read_csv(per[c])
        ex = (g["incr_over_pcs20"] - g["incr_over_pcs20_randctrl"]).dropna().values
        data.append(ex); pos.append(100 * (ex > 0).mean())
    vp = b.violinplot(data, positions=range(5), widths=0.8, showextrema=False, showmedians=False)
    for body, c in zip(vp["bodies"], COH):
        body.set_facecolor(S.ORG[c]); body.set_edgecolor("none"); body.set_alpha(0.75)
    for i, (ex, c) in enumerate(zip(data, COH)):
        b.scatter(i, np.median(ex), marker="_", s=180, linewidths=1.5, color=S.INK, zorder=5)
        b.text(i, b.get_ylim()[1] if False else max(np.percentile(ex, 99), 0.55), f"{pos[i]:.0f}%\n>0",
               ha="center", va="bottom", fontsize=5.6, color=S.INK)
    b.axhline(0, color="k", lw=0.5)
    b.set_xticks(range(5)); b.set_xticklabels(COH, fontsize=6.2)
    b.set_ylabel("per-gene excess of real morphology\nover same-size noise (ΔR², txome baseline)", fontsize=6.6)
    b.set_ylim(-0.6, 0.85)
else:
    b.text(0.5, 0.5, "per-gene panel pending:\nWinSCP transcriptome_baseline_<cohort>.csv\ninto figdata/",
           ha="center", va="center", fontsize=6.4, color=S.LGREY, style="italic", transform=b.transAxes)
    b.set_xticks([]); b.set_yticks([])
b.set_title("Gene-by-gene excess over noise", fontsize=7.6, fontweight="bold", loc="left")
b.text(-0.25, 1.06, "b", transform=b.transAxes, fontsize=10, fontweight="bold")

# ---- c: collinearity with the RNA PCs, morphology vs noise (dumbbell) ----
c_ = fig.add_subplot(gs[2])
for i, c in enumerate(COH):
    r = d.loc[c]; y = 4 - i
    c_.plot([r.cc_noise_med, r.cc_morph_med], [y, y], color=S.LGREY, lw=1.4, zorder=1)
    c_.scatter(r.cc_noise_med, y, s=26, facecolors="white", edgecolors=S.ORG[c], linewidths=1.0, zorder=3)
    c_.scatter(r.cc_morph_med, y, s=30, color=S.ORG[c], edgecolors=S.INK, linewidths=0.4, zorder=4)
    c_.text(r.cc_morph_med + 0.02, y, f"{int(r.cc_morph_gt05)} vs {int(r.cc_noise_gt05)} of 20 > 0.5",
            fontsize=5.2, color=S.GREY, va="center")
c_.set_yticks(range(5)); c_.set_yticklabels(COH[::-1], fontsize=6.2)
c_.set_xlim(0.25, 0.78)
c_.set_xlabel("median canonical correlation\nwith the 20 RNA PCs", fontsize=7)
c_.set_title("Collinear with the RNA PCs\n(so a is an upper bound)", fontsize=7, fontweight="bold", loc="left")
c_.text(0.5, -0.36, "filled = morphology PCs;  open = same-size noise", transform=c_.transAxes,
        ha="center", fontsize=5.4, color=S.GREY)
c_.text(-0.3, 1.06, "c", transform=c_.transAxes, fontsize=10, fontweight="bold")

S.save_pub(fig, "Fig_txome")
print("wrote Fig_txome", "(per-gene panel b drawn)" if have else "(panel b placeholder; per-gene CSVs missing)")
