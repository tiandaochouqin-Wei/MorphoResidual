#!/usr/bin/env python3
"""Methods pipeline schematic (new, for the top of the Methods section): data ->
morphology representation -> residual estimand -> significance test, with the
confound/robustness controls and the downstream analyses that consume the
significant set, laid out as a technical flowchart (boxes + arrows), not a data
plot. Style matches make_fig1_flagship.py (mrstyle colours, FancyArrowPatch).

NOTE: matplotlib here runs WITHOUT text.usetex, so plain (non-$...$) strings
must never contain LaTeX-only escapes like \\, or {,} -- they render literally.

TYPOGRAPHY CONTRACT (Nature Communications technical check): every authored
font size in this file is >= 7.0 pt, so that after \\includegraphics rescaling
nothing prints below 5 pt.  Mathtext is deliberately avoided for the two
symbols that used to need it -- the superscript in $R^2$ was typeset by
matplotlib at 0.7x the base size (5.9 pt -> 4.1 pt authored, 3.4 pt printed),
which was the single worst offender in the figure.  Both are now the literal
Unicode characters U+00B2 and U+2194, which render at the full 7 pt.  Box widths
below are sized from measured Arial 7 pt string widths (see _audit at the end,
which re-checks every label against its box on each build)."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import mrstyle as S

INK, GREY, LGREY = S.INK, S.GREY, S.LGREY
ORG = S.ORG
SPINE = "#0B7A6E"     # core-pipeline teal (matches Fig 1's UNI colour)
CTRL = "#5B6B7A"      # robustness-control slate
DOWN = "#B5541F"      # downstream-analysis burnt orange
DATA = "#6B4C9A"      # data-layer purple

FS = 7.0              # authored floor: nothing in this figure is smaller
SUP2 = "²"       # superscript two, full-size glyph (not mathtext)
LRARR = "↔"      # left-right arrow, full-size glyph (not mathtext)

W, Y0, Y1 = 720, 40, 500
fig = plt.figure(figsize=(7.2, 5.0))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(Y0, Y1)
ax.set_aspect("equal"); ax.axis("off")

_BOXES = []           # (text_artist, x, y, w, h, label) for the build-time fit audit


def box(x, y, w, h, text, fc, ec, fs=FS, fw="normal", tc=None, round_=6, lw=1.0, z=4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={round_}",
                 facecolor=fc, edgecolor=ec, lw=lw, zorder=z))
    t = ax.text(x + w / 2, y + h / 2, text, fontsize=fs, fontweight=fw, color=tc or INK,
                ha="center", va="center", zorder=z + 1, linespacing=1.3)
    _BOXES.append((t, x, y, w, h))
    return x, y, w, h


def arrow(x1, y1, x2, y2, c=INK, lw=1.1, ms=7, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=ms,
                 lw=lw, color=c, shrinkA=0, shrinkB=0, zorder=3, linestyle=ls))


def tag(x, y, text, c=INK, fs=7.2, fw="bold", ha="left"):
    ax.text(x, y, text, fontsize=fs, fontweight=fw, color=c, ha=ha, va="center")


# ============================================================ row 1: data
tag(8, 490, "DATA", c=DATA)
cohorts = [("CCRCC", "Kidney"), ("LUAD", "Lung"), ("UCEC", "Uterus"), ("GBM", "Brain"), ("PDAC", "Pancreas")]
cw, gap, x0 = 132, 10, 10
for i, (c, organ) in enumerate(cohorts):
    x = x0 + i * (cw + gap)
    box(x, 446, cw, 32, f"{c}\n({organ})", fc="#FFFFFF", ec=ORG[c], tc=ORG[c], fw="bold", fs=FS, lw=1.3)
    arrow(x + cw / 2, 446, x + cw / 2, 430)
box(10, 404, 700, 24, "Proteome (MS, TMT)          +          RNA-seq          +          H&E whole-slide image",
    fc="#F2F2F2", ec=GREY, fs=FS, fw="bold")

# ============================================================ row 2: representation spine
tag(8, 391, "MORPHOLOGY REPRESENTATION", c=SPINE)
arrow(360, 404, 360, 380)
seq = [("WSI tiled\n(256 px, tissue-filtered)", 160), ("UNI foundation\nmodel (frozen)", 130),
       ("mean-pool\nper slide/patient", 132), ("PCA, 20 components\n(deterministic full SVD)", 152)]
x = 20
for i, (t, w) in enumerate(seq):
    box(x, 346, w, 32, t, fc="#E6F3EF", ec=SPINE, tc=SPINE, fs=FS)
    if i < len(seq) - 1:
        arrow(x + w, 362, x + w + 12, 362, c=SPINE)
    x += w + 12

# ============================================================ row 3: residual estimand
tag(8, 320, "RESIDUAL ESTIMAND", c=SPINE)
arrow(600, 346, 600, 328, c=SPINE)
ax.text(607, 337, "morph.\nPCs", fontsize=FS, color=SPINE, va="center", ha="left", linespacing=1.1)
box(20, 260, 224, 46,
    "protein = f(own mRNA) + residual\nnested ridge, closed form (λ=1)\nbaseline vs. baseline+morphology",
    fc="#FFFFFF", ec=SPINE, fs=FS)
arrow(244, 283, 266, 283, c=SPINE)
box(266, 260, 200, 46,
    f"5-fold out-of-fold R{SUP2}\n(paired, identical folds)\nincremental R{SUP2} = morph. gain",
    fc="#FFFFFF", ec=SPINE, fs=FS)
arrow(466, 283, 488, 283, c=SPINE)
box(488, 260, 192, 46,
    f"patient{LRARR}slide permutation\n(N=1,000), BH-FDR<0.05\nand incremental R{SUP2}>0",
    fc="#FFFFFF", ec=SPINE, fs=FS)
arrow(584, 260, 584, 244, c=SPINE, lw=1.4, ms=9)
box(200, 206, 440, 36,
    "significant morphology-predictable residual proteins\n(702-2,566 per cohort)",
    fc=SPINE, ec=SPINE, tc="white", fw="bold", fs=FS)

# ============================================================ row 4: confound / robustness controls
tag(8, 192, "CONFOUND AND ROBUSTNESS CONTROLS", c=CTRL)
arrow(420, 206, 420, 178, c=GREY, lw=1.0, ms=7)
ctrl_items = [("Capacity\n(noise, rand. proj.)", 94), ("Scanner and batch\n(operator, plex)", 97),
              ("Tumour purity", 89), ("Tissue composition\n(stroma, immune,\nproliferation)", 99),
              ("Transcriptome-wide\n(own mRNA +\n20 RNA PCs)", 101), ("Encoder\n(Phikon, ResNet-50)", 104),
              ("Aggregation\n(ABMIL vs.\nmean-pool)", 80)]
gap2 = 6; x = 10
for t, w in ctrl_items:
    box(x, 132, w, 44, t, fc="#EEF1F4", ec=CTRL, tc=CTRL, fs=FS)
    x += w + gap2

# ============================================================ row 5: downstream analyses
tag(8, 118, "DOWNSTREAM ANALYSES", c=DOWN)
down_items = [("Pathway enrichment\n(Reactome, tested bkg.)", 119), ("Leave-one-organ-out\ntransfer", 105),
              ("Molecular subtypes", 99), ("mTORC1 phospho\nmechanism anchor", 98),
              ("Clinical association\n(grade, stage, survival)", 114), ("External replication\n(TCGA-KIRC, RPPA)", 108)]
gap3 = 11; x = 11
for t, w in down_items:
    box(x, 72, w, 32, t, fc="#FBEDE3", ec=DOWN, tc=DOWN, fs=FS)
    x += w + gap3

# ============================================================ legend
ax.text(W / 2, 50, "Full procedures for every stage are given in the corresponding "
        "Methods subsections below.", fontsize=FS, color=GREY, ha="center", style="italic")


def _audit():
    """Warn if any label has outgrown its box (data units == 1/100 inch here)."""
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bad = 0
    for t, x, y, w, h in _BOXES:
        bb = inv.transform_bbox(t.get_window_extent(renderer=r))
        if bb.width > w - 8 or bb.height > h - 6:
            bad += 1
            print(f"  TIGHT {t.get_text()[:28]!r}: text {bb.width:.0f}x{bb.height:.0f} "
                  f"in box {w}x{h}")
    print(f"fit audit: {bad} tight box(es) of {len(_BOXES)}")


fig.canvas.draw()
_audit()
S.save_pub(fig, "Fig_methods_pipeline", tiff=False)
print("wrote Fig_methods_pipeline")
