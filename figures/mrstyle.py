"""Shared submission-grade style for MorphoResidual figures (Genome Medicine).
Arial 7pt, editable-text SVG, 600-dpi TIFF, Okabe-Ito colourblind-safe palette."""
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    # DejaVu is last-resort glyph cover, not a design choice: Arial has the
    # superscript digits but not U+207B / U+2080-2089, and mathtext is not an
    # option because matplotlib sets sub/superscripts at 0.7x, which drops a
    # 7 pt label to 4.9 pt and under the 5 pt print floor. Latin text still
    # resolves to Arial.
    "svg.fonttype": "none",      # editable text in SVG
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "legend.frameon": False, "savefig.dpi": 600, "figure.dpi": 150,
})

COH = ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]
ORGAN = {"CCRCC": "Kidney", "LUAD": "Lung", "UCEC": "Uterus", "GBM": "Brain", "PDAC": "Pancreas"}
ORG = {"CCRCC": "#0072B2", "LUAD": "#009E73", "UCEC": "#CC79A7", "GBM": "#E69F00", "PDAC": "#D55E00"}
FAM = {"translation": "#009E73", "secretion": "#E69F00", "splicing": "#0072B2",
       "ptm": "#CC79A7", "ecm": "#8C8C8C", "folding": "#56B4E9"}
INK, GREY, LGREY = "#2F2F2F", "#8C8C8C", "#D9D9D9"


def save_pub(fig, name, tiff=True):
    fig.savefig(f"{name}.svg", bbox_inches="tight")
    fig.savefig(f"{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{name}.png", dpi=300, bbox_inches="tight")
    if tiff:
        fig.savefig(f"{name}.tiff", dpi=600, bbox_inches="tight")


def panel_label(ax, letter, x=-0.06, y=1.04, size=9):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="right")


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else ""
