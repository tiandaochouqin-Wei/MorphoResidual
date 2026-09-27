#!/usr/bin/env python3
"""Flagship Figure 1 panel -- study schematic for MorphoResidual.

Real H&E micrographs (CCRCC tiles, the only raster) + vector organ icons + vector text.
Every label is a live matplotlib text artist authored at >= 7 pt (headers 8 pt).

The schematic is drawn by draw_flagship(ax) so that make_master_composite.py can place it
NATIVELY as panel a of Fig_master (previously it pasted Fig1_flagship.png with
aspect='auto', which shrank the 8 pt labels to ~3 pt printed and stretched the icons).
Data coordinates are XMAX x (YTOP-YBOT) units; the axes must have an equal aspect with
UNIT_IN inches per unit (PANEL_W_IN / XMAX) so that 1 unit = 0.62 pt and the authored point
sizes below are what is printed at authored scale.  Run standalone to (re)write
Fig1_flagship.{svg,pdf,png,tiff}."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import (Ellipse, Circle, Rectangle, FancyBboxPatch,
                                FancyArrowPatch, Polygon)
import mrstyle as S

INK, GREY, LGREY = S.INK, S.GREY, S.LGREY
RESID = "#D55E00"        # residual accent (vermillion)
MRNA = "#C9C9C9"         # mRNA-explained grey
UNI = "#0B7A6E"          # foundation-model teal
GBOX_F, GBOX_E = "#E6F3EF", "#0B7A6E"   # predict box
OBOX_F, OBOX_E = "#FCEBE0", "#D55E00"   # enriched box
ORG = S.ORG
TILE = "figdata/tile_{:02d}.png"

FS = 7.0                 # body / label size (pt, authored)
FH = 8.0                 # stage-header size (pt, authored)
XMAX, YBOT, YTOP = 752, 98, 262
PANEL_W_IN = 6.3         # width the panel occupies in Fig_master (0.875 x 7.2 in)
UNIT_IN = PANEL_W_IN / XMAX
PANEL_H_IN = (YTOP - YBOT) * UNIT_IN


def draw_flagship(ax):
    ax.set_xlim(0, XMAX); ax.set_ylim(YBOT, YTOP)
    ax.set_aspect("equal"); ax.axis("off")

    def disc(x, y, n, r=8.5, c=INK):
        ax.add_patch(Circle((x, y), r, facecolor=c, ec="none", zorder=5))
        ax.text(x, y - 0.3, str(n), color="white", fontsize=FS, fontweight="bold",
                ha="center", va="center", zorder=6)

    def header(x, y, n, t):
        disc(x, y, n)
        ax.text(x + 14, y, t, fontsize=FH, fontweight="bold", color=INK, va="center", ha="left")

    def tile(path, x0, y0, s, border=INK, lw=0.8):
        ax.imshow(mpimg.imread(path), extent=[x0, x0 + s, y0, y0 + s], zorder=3, aspect="auto")
        ax.add_patch(Rectangle((x0, y0), s, s, facecolor="none", edgecolor=border, lw=lw, zorder=4))

    def arrow(x1, y1, x2, y2, c=INK, lw=1.0, ms=7):
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
            ax.plot(xs, cy + dy * r + 0.06 * r * np.sin((xs - cx) / (0.13 * r)), color="white", lw=0.7, zorder=4)
            ax.plot(-xs + 2 * cx, cy + dy * r + 0.06 * r * np.sin((xs - cx) / (0.13 * r)), color="white", lw=0.7, zorder=4)

    def ic_uterus(cx, cy, r, c):
        ax.add_patch(Polygon([(cx - 0.42 * r, cy + 0.28 * r), (cx + 0.42 * r, cy + 0.28 * r),
                              (cx + 0.2 * r, cy - 0.85 * r), (cx - 0.2 * r, cy - 0.85 * r)],
                             closed=True, facecolor=c, ec="none", zorder=3))
        for s in (-1, 1):
            ax.add_patch(FancyArrowPatch((cx + s * 0.35 * r, cy + 0.25 * r), (cx + s * 0.92 * r, cy + 0.8 * r),
                         arrowstyle="-", lw=1.8, color=c, zorder=3))
            ax.add_patch(Circle((cx + s * 0.95 * r, cy + 0.85 * r), 0.16 * r, facecolor=c, ec="none", zorder=3))

    def ic_pancreas(cx, cy, r, c):
        ax.add_patch(Ellipse((cx + 0.15 * r, cy - 0.1 * r), 2.0 * r, 0.66 * r, angle=-18,
                     facecolor=c, ec="none", zorder=3))
        ax.add_patch(Circle((cx - 0.82 * r, cy + 0.24 * r), 0.4 * r, facecolor=c, ec="none", zorder=3))

    ORGAN_IC = {"CCRCC": ic_kidney, "LUAD": ic_lung, "UCEC": ic_uterus, "GBM": ic_brain, "PDAC": ic_pancreas}

    # separators
    for xs in (238, 498):
        ax.plot([xs, xs], [104, 238], color=LGREY, lw=0.8, zorder=1)

    # ==================== STAGE 1 : the beyond-mRNA proteome ====================
    header(26, 250, 1, "The beyond-mRNA proteome")
    ax.text(26, 222, "protein  =  f(mRNA)  +", fontsize=FS, color=INK, va="center", ha="left")
    ax.text(146, 222, "residual", fontsize=FS, color=RESID, fontweight="bold", va="center", ha="left")
    # stacked bar; segment widths are set by the text each must hold
    ax.add_patch(Rectangle((26, 190), 110, 22, facecolor=MRNA, ec="white", lw=1, zorder=3))
    ax.add_patch(Rectangle((136, 190), 60, 22, facecolor=RESID, ec="white", lw=1, zorder=3))
    ax.text(81, 201, "explained by mRNA", fontsize=FS, color=INK, ha="center", va="center", zorder=4)
    ax.text(166, 201, "residual", fontsize=FS, color="white", fontweight="bold", ha="center", va="center", zorder=4)
    # elbow from residual down to caption
    ax.plot([166, 166, 40], [190, 178, 178], color=RESID, lw=1.0, zorder=3)
    arrow(40, 178, 40, 168, c=RESID, lw=1.0, ms=6)
    ax.text(26, 162, "Post-transcriptional layer", fontsize=FS, color=RESID, fontweight="bold", va="top", ha="left")
    ax.text(26, 148, "translation · secretion · degradation", fontsize=FS, color=INK, va="top", ha="left")
    ax.text(26, 134, "— not read by the gene's own transcript", fontsize=FS, color=GREY, style="italic", va="top", ha="left")

    # ==================== STAGE 2 : read morphology from H&E ====================
    ox = -8                                   # stage-2 origin offset (keeps the layout compact)
    header(258 + ox, 250, 2, "Read morphology from H&E")
    # WSI (real tissue) with inset region
    tile(TILE.format(4), 258 + ox, 164, 54)
    ax.text(285 + ox, 158, "H&E whole-\nslide image", fontsize=FS, color=INK, ha="center", va="top", linespacing=1.15)
    ax.add_patch(Rectangle((292 + ox, 190), 12, 12, facecolor="none", edgecolor=INK, lw=0.9, zorder=5))
    # 2x2 extracted tiles
    tx, ty, ts = [326 + ox, 352 + ox], [194, 168], 24
    for i, xx in enumerate(tx):
        for j, yy in enumerate(ty):
            tile(TILE.format(i * 2 + j), xx, yy, ts, lw=0.7)
    ax.text(351 + ox, 158, "tiles", fontsize=FS, color=INK, ha="center", va="top")
    # zoom lines from inset to tile grid
    for (px, py, qy) in [(304, 202, 218), (304, 190, 168)]:
        ax.plot([px + ox, 326 + ox], [py, qy], color=GREY, lw=0.6, ls=(0, (3, 2)), zorder=2)
    # UNI model
    ax.add_patch(FancyBboxPatch((392 + ox, 176), 66, 36, boxstyle="round,pad=2,rounding_size=4",
                 facecolor=UNI, ec="none", zorder=3))
    ax.text(425 + ox, 194, "UNI\nfoundation\nmodel", fontsize=FS, color="white", fontweight="bold",
            ha="center", va="center", linespacing=1.15, zorder=4)
    arrow(379 + ox, 194, 390 + ox, 194)
    # embedding cloud
    rng = np.random.RandomState(3)
    ecx, ecy = 471 + ox, 194
    for _ in range(16):
        a, rr = rng.rand() * 2 * np.pi, rng.rand() ** 0.5 * 6.0
        cc = ORG[list(ORG)[rng.randint(5)]]
        ax.add_patch(Circle((ecx + rr * np.cos(a), ecy + rr * np.sin(a)), 1.4, facecolor=cc, ec="none", zorder=3))
    arrow(461 + ox, 194, 464 + ox, 194)
    ax.text(ecx, 176, "slide\nembedding", fontsize=FS, color=INK, ha="center", va="top", linespacing=1.15)

    # ==================== STAGE 3 : morphology predicts the residual ====================
    X3 = 506
    header(X3 + 8, 250, 3, "Morphology predicts it")
    ax.add_patch(FancyBboxPatch((X3, 207), 240, 22, boxstyle="round,pad=1,rounding_size=3",
                 facecolor=GBOX_F, ec=GBOX_E, lw=1.0, zorder=3))
    ax.text(X3 + 120, 218, "Incremental R² over an mRNA-only baseline", fontsize=FS, color=INK,
            ha="center", va="center", zorder=4)
    ax.add_patch(FancyBboxPatch((X3, 181), 240, 22, boxstyle="round,pad=1,rounding_size=3",
                 facecolor=OBOX_F, ec=OBOX_E, lw=1.0, zorder=3))
    ax.text(X3 + 120, 192, "Enriched for translation · secretion · splicing", fontsize=FS, color=INK,
            ha="center", va="center", zorder=4)
    ax.text(X3 + 4, 166, "Reproducible across five organs", fontsize=FS, color=INK, fontweight="bold",
            ha="left", va="center")
    # organ icons + counts
    counts = {"CCRCC": "2,191", "LUAD": "2,310", "UCEC": "2,566", "GBM": "702", "PDAC": "1,572"}
    xs = np.linspace(X3 + 32, X3 + 208, 5)
    for xi, k in zip(xs, S.COH):
        ORGAN_IC[k](xi, 143, 10, ORG[k])
        ax.text(xi, 126, S.ORGAN[k], fontsize=FS, color=INK, ha="center", va="top")
        ax.text(xi, 113, counts[k], fontsize=FS, color=ORG[k], fontweight="bold", ha="center", va="top")


if __name__ == "__main__":
    fig = plt.figure(figsize=(PANEL_W_IN, PANEL_H_IN))
    ax = fig.add_axes([0, 0, 1, 1])
    draw_flagship(ax)
    S.save_pub(fig, "Fig1_flagship")
    print("wrote Fig1_flagship")
