#!/usr/bin/env python3
r"""Supplementary figure Fig_tiles_supp: ALL 16 candidate H&E tiles behind Fig. 3f (Fig_biology panel f).

Provenance of the 16 tiles (figdata/tile_{hi,lo}_{i}_{j}.png), from fig_data_export2.py:
  * CCRCC patients with a measured translation-residual score (figdata/scores_ccrcc.csv, column
    translation_meas, n = 103) are sorted ascending; lo_cases = first four, hi_cases = last four.
    Tile file index i is the position in those lists, so hi_3 is the highest-scoring patient and
    lo_0 the lowest.
  * For each patient the first slide file (sorted by name) is scanned on a fixed raster
    (start 5 tiles in, step 8 tiles) at level 0; the first two 256x256 tiles with mean HSV saturation
    above 15 in more than 35 % of pixels are saved (j = 0, 1). No other selection is applied there.
  * make_fig_biology.py panel f then shows 4 + 4 of the 16, chosen by visual inspection
    (make_realdata_figs2.py comment: 'drop blank/blood artefacts by index, from visual QA').
The reason tags below are the authors' reading of each excluded tile image (this script has no
per-tile QA log to draw on; none was saved when the selection was made).
Run from this directory: python make_fig_tiles_supp.py"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe
import mrstyle as S

mpl.rcParams.update({"font.family": ["Arial", "DejaVu Sans"]})
DD = "figdata"

# ---- which (i, j) tiles Fig. 3f shows: copied verbatim from make_fig_biology.py lines 363-364 ----
SHOWN_HI = [(1, 0), (2, 0), (2, 1), (3, 1)]
SHOWN_LO = [(0, 1), (1, 1), (3, 0), (3, 1)]
SHOWN = {("hi", i, j) for i, j in SHOWN_HI} | {("lo", i, j) for i, j in SHOWN_LO}

# ---- reason tags for the 8 candidates not shown (visual reading of the tile files) ----
REASON = {
    ("hi", 0, 0): "glass, blurred",
    ("hi", 0, 1): "blood",
    ("hi", 1, 1): "pale stroma",
    ("hi", 3, 0): "sparse stroma",
    ("lo", 0, 0): "pale stroma",
    ("lo", 1, 0): "blood",
    ("lo", 2, 0): "dark, blurred",
    ("lo", 2, 1): "dark, saturated",
}

# ---- case IDs and measured score, replicating fig_data_export2.py ----
sc = pd.read_csv(f"{DD}/scores_ccrcc.csv", index_col=0)
col = "translation_meas"
srt = sc.dropna(subset=[col]).sort_values(col)
N = len(srt)
lo_cases = list(srt.index[:4])
hi_cases = list(srt.index[-4:])
score = srt[col].to_dict()
rank_hi = {c: N - k for k, c in enumerate(srt.index)}          # 1 = highest

# ---- geometry (inches) ----
FW = 7.2
S_IN, G_IN, P_IN = 0.76, 0.03, 0.14      # tile size, gap within patient, gap between patients
LEFT = 0.55
HEAD, TAG, VGAP = 0.36, 0.46, 0.06
BLOCK = HEAD + S_IN + TAG + VGAP
FH = 2 * BLOCK + 0.05

MPP_UM, TILE_PX, SCALEBAR_UM = 0.4942, 256, 50            # as in make_fig_biology.py
SCALEBAR_PX = SCALEBAR_UM / MPP_UM
BLUE = S.ORG["CCRCC"]


def add_scalebar(ax, margin=14):
    x1 = y = TILE_PX - margin
    x0 = x1 - SCALEBAR_PX
    ax.plot([x0, x1], [y, y], color="white", lw=2.4, solid_capstyle="butt", zorder=5,
            path_effects=[pe.Stroke(linewidth=3.6, foreground="black"), pe.Normal()])
    ax.text((x0 + x1) / 2.0, y - 7, f"{SCALEBAR_UM} µm", fontsize=6.0, color="white",
            ha="center", va="bottom", zorder=5,
            path_effects=[pe.Stroke(linewidth=1.4, foreground="black"), pe.Normal()])


fig = plt.figure(figsize=(FW, FH))


def fx(x_in):
    return x_in / FW


def fy(y_in):
    return y_in / FH


n_shown = 0
groups = [("hi", [hi_cases[i] for i in (3, 2, 1, 0)], (3, 2, 1, 0), "highest\nresidual", BLUE),
          ("lo", lo_cases, (0, 1, 2, 3), "lowest\nresidual", S.GREY)]
for g, (tag, cases, idxs, glabel, gcol) in enumerate(groups):
    top = FH - g * BLOCK
    y_tile_top = top - HEAD
    y_tile_bot = y_tile_top - S_IN
    fig.text(fx(LEFT - 0.07), fy(y_tile_top - S_IN / 2), glabel, fontsize=7.0, fontweight="bold",
             color=gcol, ha="right", va="center", linespacing=1.1)
    for p, (case, i) in enumerate(zip(cases, idxs)):
        x0 = LEFT + p * (2 * S_IN + G_IN + P_IN)
        rk_txt = f"rank {rank_hi[case]} of {N}"      # one ranking for both rows: 1 = highest score
        fig.text(fx(x0 + S_IN + G_IN / 2), fy(top - 0.02),
                 case, fontsize=7.0, fontweight="bold", color=S.INK, ha="center", va="top")
        fig.text(fx(x0 + S_IN + G_IN / 2), fy(top - 0.15),
                 f"score {score[case]:+.2f}, {rk_txt}".replace("-", "−"), fontsize=6.0, color=S.INK, ha="center", va="top")
        for j in range(2):
            xt = x0 + j * (S_IN + G_IN)
            ax = fig.add_axes([fx(xt), fy(y_tile_bot), fx(S_IN), fy(S_IN)])
            ax.imshow(mpimg.imread(f"{DD}/tile_{tag}_{i}_{j}.png"))
            ax.set_xticks([]); ax.set_yticks([])
            shown = (tag, i, j) in SHOWN
            n_shown += shown
            for sp in ax.spines.values():
                sp.set_visible(True)
                sp.set_edgecolor(BLUE if shown else S.LGREY)
                sp.set_linewidth(1.6 if shown else 0.5)
            add_scalebar(ax)
            tid = f"{tag}_{i}_{j}"
            if shown:
                fig.text(fx(xt + S_IN / 2), fy(y_tile_bot - 0.03), f"{tid}\nin Fig. 3f",
                         fontsize=6.0, color=BLUE, fontweight="bold", ha="center", va="top", linespacing=1.15)
            else:
                fig.text(fx(xt + S_IN / 2), fy(y_tile_bot - 0.03),
                         f"{tid}\nnot shown\n{REASON[(tag, i, j)]}",
                         fontsize=6.0, color=S.INK, ha="center", va="top", linespacing=1.15)

assert n_shown == 8, n_shown
S.save_pub(fig, "Fig_tiles_supp")
plt.close(fig)
print(f"Fig_tiles_supp: 16 candidates, {n_shown} shown in Fig. 3f; hi={hi_cases[::-1]} lo={lo_cases}; "
      f"scores hi {[round(score[c], 2) for c in hi_cases[::-1]]} lo {[round(score[c], 2) for c in lo_cases]}")
