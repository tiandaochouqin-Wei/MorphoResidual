#!/usr/bin/env python3
"""Flagship Figure 1 — top-journal study schematic for MorphoResidual.
Real H&E micrographs (CCRCC tiles) + vector organ icons + strict grid.
Data coords: 720 x 270 (100 units = 1 inch), equal aspect."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
from matplotlib.patches import (Ellipse, Circle, Rectangle, FancyBboxPatch,
                                FancyArrowPatch, Polygon, Arc)
import mrstyle as S

INK, GREY, LGREY = S.INK, S.GREY, S.LGREY
RESID = "#D55E00"        # residual accent (vermillion)
MRNA = "#C9C9C9"         # mRNA-explained grey
UNI = "#0B7A6E"          # foundation-model teal
GBOX_F, GBOX_E = "#E6F3EF", "#0B7A6E"   # predict box
OBOX_F, OBOX_E = "#FCEBE0", "#D55E00"   # enriched box
ORG = S.ORG
TILE = "figdata/tile_{:02d}.png"

fig = plt.figure(figsize=(7.2, 1.86))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 720); ax.set_ylim(78, 264)
ax.set_aspect("equal"); ax.axis("off")


def disc(x, y, n, r=8.5, c=INK):
    ax.add_patch(Circle((x, y), r, facecolor=c, ec="none", zorder=5))
    ax.text(x, y - 0.3, str(n), color="white", fontsize=8, fontweight="bold",
            ha="center", va="center", zorder=6)


def header(x, y, n, t):
    disc(x, y, n)
    ax.text(x + 14, y, t, fontsize=9, fontweight="bold", color=INK, va="center", ha="left")


def tile(path, x0, y0, s, border=INK, lw=0.8):
    ax.imshow(mpimg.imread(path), extent=[x0, x0 + s, y0, y0 + s], zorder=3, aspect="auto")
    ax.add_patch(Rectangle((x0, y0), s, s, facecolor="none", edgecolor=border, lw=lw, zorder=4))


def arrow(x1, y1, x2, y2, c=INK, lw=1.3, ms=8):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=ms,
                 lw=lw, color=c, shrinkA=0, shrinkB=0, zorder=4))


# ---------- organ vector icons ----------
def ic_kidney(cx, cy, r, c):
    ax.add_patch(Ellipse((cx, cy), 1.7 * r, 1.35 * r, facecolor=c, ec="none", zorder=3))
    ax.add_patch(Circle((cx, cy + 0.82 * r), 0.55 * r, facecolor="white", ec="none", zorder=4))


def ic_lung(cx, cy, r, c):
    for s in (-1, 1):
        ax.add_patch(Ellipse((cx + s * 0.5 * r, cy - 0.05 * r), 0.92 * r, 1.7 * r,
                     facecolor=c, ec="none", zorder=3))
    ax.add_patch(Rectangle((cx - 0.07 * r, cy - 0.1 * r), 0.14 * r, 1.05 * r, facecolor=c, ec="none", zorder=3))
    ax.add_patch(Rectangle((cx - 0.34 * r, cy + 0.5 * r), 0.68 * r, 0.12 * r, facecolor=c, ec="none", zorder=3))
    ax.add_patch(Rectangle((cx - 0.05 * r, cy - 0.95 * r), 0.10 * r, 1.05 * r, facecolor="white", ec="none", zorder=4))


def ic_brain(cx, cy, r, c):
    ax.add_patch(Ellipse((cx, cy), 1.75 * r, 1.4 * r, facecolor=c, ec="none", zorder=3))
    ax.add_patch(Rectangle((cx - 0.03 * r, cy - 0.62 * r), 0.06 * r, 1.24 * r, facecolor="white", ec="none", zorder=4))
    for dy in (0.35, -0.05, -0.45):
        xs = np.linspace(cx - 0.7 * r, cx - 0.08 * r, 20)
        ax.plot(xs, cy + dy * r + 0.06 * r * np.sin((xs - cx) / (0.13 * r)), color="white", lw=0.8, zorder=4)
        ax.plot(-xs + 2 * cx, cy + dy * r + 0.06 * r * np.sin((xs - cx) / (0.13 * r)), color="white", lw=0.8, zorder=4)


def ic_uterus(cx, cy, r, c):
    ax.add_patch(Polygon([(cx - 0.42 * r, cy + 0.28 * r), (cx + 0.42 * r, cy + 0.28 * r),
                          (cx + 0.2 * r, cy - 0.85 * r), (cx - 0.2 * r, cy - 0.85 * r)],
                         closed=True, facecolor=c, ec="none", zorder=3))
    for s in (-1, 1):
        ax.add_patch(FancyArrowPatch((cx + s * 0.35 * r, cy + 0.25 * r), (cx + s * 0.92 * r, cy + 0.8 * r),
                     arrowstyle="-", lw=2.0, color=c, zorder=3))
        ax.add_patch(Circle((cx + s * 0.95 * r, cy + 0.85 * r), 0.16 * r, facecolor=c, ec="none", zorder=3))


def ic_pancreas(cx, cy, r, c):
    ax.add_patch(Ellipse((cx + 0.15 * r, cy - 0.1 * r), 2.0 * r, 0.66 * r, angle=-18,
                 facecolor=c, ec="none", zorder=3))
    ax.add_patch(Circle((cx - 0.82 * r, cy + 0.24 * r), 0.4 * r, facecolor=c, ec="none", zorder=3))


ORGAN_IC = {"CCRCC": ic_kidney, "LUAD": ic_lung, "UCEC": ic_uterus, "GBM": ic_brain, "PDAC": ic_pancreas}

# separators
for xs in (244, 480):
    ax.plot([xs, xs], [90, 242], color=LGREY, lw=0.8, zorder=1)

# ==================== STAGE 1 : the beyond-mRNA proteome ====================
header(26, 250, 1, "The beyond-mRNA proteome")
ax.text(26, 208, "protein  =  f(mRNA)  +", fontsize=9, color=INK, va="center", ha="left")
ax.text(163, 208, "residual", fontsize=9, color=RESID, fontweight="bold", va="center", ha="left")
# stacked bar
# Segment widths are set by the 8 pt text each must hold. They fitted at the smaller
# font this schematic was first drawn with; at 8 pt "residual" overflows a 40-unit
# box and spills onto the grey segment, so the split moves to 106/58.
ax.add_patch(Rectangle((26, 168), 106, 22, facecolor=MRNA, ec="white", lw=1, zorder=3))
ax.add_patch(Rectangle((132, 168), 58, 22, facecolor=RESID, ec="white", lw=1, zorder=3))
ax.text(79, 179, "explained by mRNA", fontsize=8, color="white", ha="center", va="center", zorder=4)
ax.text(161, 179, "residual", fontsize=8, color="white", fontweight="bold", ha="center", va="center", zorder=4)
# elbow from residual down to caption
ax.plot([161, 161, 40], [168, 156, 156], color=RESID, lw=1.1, zorder=3)
arrow(40, 156, 40, 146, c=RESID, lw=1.1, ms=7)
ax.text(26, 138, "Post-transcriptional layer", fontsize=8, color=RESID, fontweight="bold", va="top", ha="left")
ax.text(26, 126, "translation · secretion · degradation", fontsize=8, color=INK, va="top", ha="left")
ax.text(26, 115, "— not read by the gene's own transcript", fontsize=8, color=GREY, style="italic", va="top", ha="left")

# ==================== STAGE 2 : read morphology from H&E ====================
header(258, 250, 2, "Read morphology from H&E")
# WSI (real tissue) with inset region
tile(TILE.format(4), 258, 150, 58)
ax.text(287, 142, "H&E whole-slide image", fontsize=8, color=INK, ha="center", va="top")
ax.add_patch(Rectangle((296, 178), 13, 13, facecolor="none", edgecolor=INK, lw=1.0, zorder=5))
# 2x2 extracted tiles
tx, ty, ts = [340, 366], [180, 154], 24
for i, xx in enumerate(tx):
    for j, yy in enumerate(ty):
        tile(TILE.format(i * 2 + j), xx, yy, ts, lw=0.7)
ax.text(365, 142, "tiles", fontsize=8, color=INK, ha="center", va="top")
# zoom lines from inset to tile grid
for (px, py) in [(309, 191), (309, 178)]:
    ax.plot([px, 340], [py, 204 if py > 184 else 154], color=GREY, lw=0.6, ls=(0, (3, 2)), zorder=2)
# UNI model
# 48x32 held these three lines at the old font size but clips "foundation" at 8 pt.
ax.add_patch(FancyBboxPatch((398, 162), 62, 36, boxstyle="round,pad=2,rounding_size=4",
             facecolor=UNI, ec="none", zorder=3))
ax.text(429, 180, "UNI\nfoundation\nmodel", fontsize=8, color="white", fontweight="bold",
        ha="center", va="center", linespacing=1.15, zorder=4)
arrow(386, 180, 398, 180)
# embedding cloud
rng = np.random.RandomState(3)
ecx, ecy = 469, 180
for _ in range(16):
    a, rr = rng.rand() * 2 * np.pi, rng.rand() ** 0.5 * 9
    cc = ORG[list(ORG)[rng.randint(5)]]
    ax.add_patch(Circle((ecx + rr * np.cos(a), ecy + rr * np.sin(a)), 1.5, facecolor=cc, ec="none", zorder=3))
arrow(462, 180, 466, 180)
ax.text(469, 152, "slide\nembedding", fontsize=8, color=INK, ha="center", va="top", linespacing=1.1)

# ==================== STAGE 3 : morphology predicts the residual ====================
header(494, 250, 3, "Morphology predicts it")
# Both labels overflowed a 204-unit box on the left at 8 pt (measured 215 and 211
# units of text). Widened to 220 and re-centred at 596; the first also drops "an".
ax.add_patch(FancyBboxPatch((486, 198), 220, 26, boxstyle="round,pad=1,rounding_size=3",
             facecolor=GBOX_F, ec=GBOX_E, lw=1.0, zorder=3))
ax.text(602, 211, "Incremental R² over an mRNA-only baseline", fontsize=8, color=INK,
        ha="center", va="center", zorder=4)
ax.add_patch(FancyBboxPatch((486, 166), 220, 26, boxstyle="round,pad=1,rounding_size=3",
             facecolor=OBOX_F, ec=OBOX_E, lw=1.0, zorder=3))
ax.text(602, 179, "Enriched for translation · secretion · splicing", fontsize=8, color=INK,
        ha="center", va="center", zorder=4)
ax.text(498, 150, "Reproducible across five organs", fontsize=8, color=INK, fontweight="bold",
        ha="left", va="center")
# organ icons + counts
counts = {"CCRCC": "2,191", "LUAD": "2,310", "UCEC": "2,566", "GBM": "702", "PDAC": "1,572"}
xs = np.linspace(524, 680, 5)
for xi, k in zip(xs, S.COH):
    ORGAN_IC[k](xi, 120, 11, ORG[k])
    ax.text(xi, 100, S.ORGAN[k], fontsize=8, color=INK, ha="center", va="top")
    ax.text(xi, 92, counts[k], fontsize=8, color=ORG[k], fontweight="bold", ha="center", va="top")

S.save_pub(fig, "Fig1_flagship")
print("wrote Fig1_flagship")
