#!/usr/bin/env python3
"""MorphoResidual Figure 1 -- concept + pipeline schematic with REAL H&E tiles.
Convention: UNI (Nat Med 2024 Fig 1d) / Wang (Cell Rep Med 2023 Fig 1a) -- rounded
stage boxes, green foundation-model box, real histology thumbnails, embedding dots."""
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
import mrstyle as S

UNIG = "#A8D5B5"; DOWN = "#F3C98B"; RES = "#D55E00"; GREY = "#D9D9D9"
TILES = sorted(glob.glob("figdata/tile_0*.png"))
IMG = [mpimg.imread(t) for t in TILES]

fig, ax = plt.subplots(figsize=(7.2, 2.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")


def rbox(x, y, w, h, fc, ec="none", lw=0, r=2.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, fc=fc, ec=ec, lw=lw,
                 boxstyle=f"round,pad=0,rounding_size={r}", mutation_aspect=1, zorder=2))


def arrow(x1, y1, x2, y2, c="#555", lw=1.5):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=10, lw=lw, color=c, shrinkA=2, shrinkB=2, zorder=6))


def head(x, t):
    ax.text(x, 38.5, t, fontsize=7.4, fontweight="bold", ha="left", color="#222")


def tile(img, x, y, w, h, ec="#9C6B8A", lw=0.6):
    ax.imshow(img, extent=[x, x + w, y, y + h], aspect="auto", zorder=3)
    ax.add_patch(Rectangle((x, y), w, h, fc="none", ec=ec, lw=lw, zorder=4))


# ===== Stage 1 : the beyond-mRNA proteome =====
head(1, "1  The beyond-mRNA proteome")
ax.text(1, 33.0, "protein  =  $f$(mRNA)  +  residual", fontsize=8.3, color="#222")
ax.add_patch(Rectangle((1, 25.5), 15, 4.2, fc=GREY, ec="none", zorder=2))
ax.add_patch(Rectangle((16, 25.5), 7.5, 4.2, fc=RES, ec="none", zorder=2))
ax.text(8.5, 27.6, "explained\nby mRNA", ha="center", va="center", fontsize=5.7, color="#555")
ax.text(19.7, 27.6, "residual", ha="center", va="center", fontsize=5.9, color="white", fontweight="bold")
ax.annotate("post-transcriptional layer\n(translation · secretion · degradation)\ninvisible to the transcriptome",
            xy=(19.7, 25.5), xytext=(2.5, 18.0), fontsize=5.9, color=RES, ha="left", va="top",
            arrowprops=dict(arrowstyle="-|>", color=RES, lw=1.0))
arrow(25, 27.6, 30.5, 27.6)

# ===== Stage 2 : read morphology from real H&E =====
head(31, "2  Read morphology from H&E")
tile(IMG[2], 31, 24.2, 8, 8.2, ec="#8A5070", lw=0.9)      # WSI region (real)
ax.text(35, 22.6, "H&E WSI", ha="center", fontsize=5.9)
arrow(39.5, 28.3, 42.3, 28.3)
grid = [0, 1, 3, 4, 5, 6]                                  # 3x2 real tiles
for idx, ti in enumerate(grid):
    r, c = divmod(idx, 3)
    tile(IMG[ti % len(IMG)], 42.6 + c * 2.35, 26.0 + (1 - r) * 2.35, 2.15, 2.15, lw=0.4)
ax.text(46.1, 22.6, "tiles", ha="center", fontsize=5.9)
arrow(49.9, 28.3, 52.5, 28.3)
rbox(53, 25.4, 10, 5.6, UNIG, r=1.5)                       # green UNI box
ax.text(58, 28.2, "UNI\nfoundation model", ha="center", va="center", fontsize=6.0,
        color="#20502F", fontweight="bold")
arrow(63.5, 28.3, 66.3, 28.3)
for i in range(4):                                         # embedding dots
    ax.add_patch(Circle((67.4 + i * 2.2, 28.3), 0.82, fc="white", ec="#6E9E7C", lw=0.8, zorder=4))
    ax.add_patch(Circle((67.4 + i * 2.2, 28.3), 0.34, fc="#4E7E5C", ec="none", zorder=5))
ax.text(70.7, 22.6, "slide embedding\n$\\rightarrow$ 20 PCs", ha="center", va="top", fontsize=5.6)
arrow(76.3, 28.3, 80.5, 24.3)

# ===== Stage 3 : does morphology predict the residual? =====
head(72.5, "3  Predict the residual?")
rbox(72.5, 18.3, 26, 5.2, "#EAF3EE", ec="#2E8B57", lw=1.0, r=1.3)
ax.text(85.5, 20.9, "incremental $R^2$ of morphology\nover an mRNA-only baseline",
        ha="center", va="center", fontsize=5.9, color="#1C5C38")
rbox(72.5, 13.3, 26, 3.4, "#FBEEE2", ec=DOWN, lw=1.0, r=1.1)
ax.text(85.5, 15.0, "enriched for translation, secretion & splicing",
        ha="center", va="center", fontsize=5.8, color=RES, fontweight="bold")
ax.text(85.5, 10.8, "reproducible across five organs", ha="center", fontsize=5.8, color="#333")
labs = ["kidney", "lung", "uterus", "brain", "pancreas"]
for i, (c, l) in enumerate(zip([S.ORG[x] for x in S.COH], labs)):
    cx = 74.5 + i * 5.6
    ax.add_patch(Circle((cx, 7.6), 1.3, fc=c, ec="none", zorder=3))
    ax.text(cx, 5.0, l, ha="center", fontsize=4.9, color="#333")

for sx in (29.5, 71.0):
    ax.plot([sx, sx], [3.5, 37], color="#EDEDED", lw=0.8, zorder=0)

fig.tight_layout(pad=0.3)
S.save_pub(fig, "Fig1_concept")
print("wrote Fig1_concept with real H&E tiles")
