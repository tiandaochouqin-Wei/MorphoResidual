#!/usr/bin/env python3
"""Supplementary composite: "Reproducibility across organs and encoders"
(Fig_supp1_reproducibility).  Five panels, two rows, no bar charts.

ROW 1  a  morphology-predictable residual proteins per organ (lollipop dot plot;
          server-verified counts / % of tested genes)
       b  leading enriched Reactome programme per organ (dot plot coloured by family;
          server-verified -log10 adj. p, tested-proteome background)
       c  t-SNE of UNI slide embeddings (figdata/emb_pca50.csv), coloured by cohort
ROW 2  d  per-gene agreement between encoders, CCRCC (hexbin of
          figdata/merged_incr_ccrcc.csv, UNI vs Phikon incremental R2)
       e  per-cohort agreement dumbbell: per-gene r (filled) and recall of
          UNI-significant proteins by Phikon (open); figdata/encoder_comparison.csv
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")                     # headless; letters are placed after a draw pass
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.manifold import TSNE
import mrstyle as S

ORG, INK, GREY, LGREY, FAM = S.ORG, S.INK, S.GREY, S.LGREY, S.FAM
COH, ORGAN = S.COH, S.ORGAN
DD = "figdata"
ROWS = np.arange(5)[::-1]          # CCRCC on top (y = 4) ... PDAC at the bottom (y = 0)


def place_letters(fig, items, gap_pt=5.0, size=10):
    """Panel letters after a draw pass: baseline on the title's first-line baseline,
    right edge `gap_pt` left of the panel's tight bbox (tick labels / axis labels)."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    _, lp_h, lp_d = rend.get_text_width_height_descent("lp", items[0][1].get_fontproperties(), ismath=False)
    placed = {}
    for ax, title, letter in items:
        tb = ax.get_tightbbox(rend)
        tt = title.get_window_extent(rend)
        line0 = title.get_text().split("\n")[0]
        _, h, d = rend.get_text_width_height_descent(line0, title.get_fontproperties(), ismath=False)
        h, d = max(h, lp_h), max(d, lp_d)          # same line metrics matplotlib uses for layout
        x = tb.x0 - gap_pt * fig.dpi / 72.0
        fx, fy = inv.transform((x, tt.y1 - (h - d)))
        fig.text(fx, fy, letter, fontsize=size, fontweight="bold", ha="right", va="baseline")
        placed[letter] = (tb.x0 / fig.dpi, fx, fy)
    return placed


fig = plt.figure(figsize=(7.2, 6.4))
# two rows of equal height (0.955 -> 0.10 with the equivalent of hspace = 0.42); row 2 starts
# further left so that panel d's y-axis label is flush with panel a's tick labels, keeping
# the left edge of the composite straight.  Panel d keeps its width; e absorbs the extra.
g1 = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.15], wspace=0.42,
                      left=0.13, right=0.975, top=0.955, bottom=0.6017)
g2 = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.036], wspace=0.42,
                      left=0.0935, right=0.975, top=0.4483, bottom=0.10)

# =====================================================================================
# a : significant residual proteins per organ
# (figdata/merged_incr_<cohort>.csv, same numbers as make_master_composite.py's
# counts/tested dicts -- FDR<0.05 & incremental_r2_uni>0, checked 2026-09-17 to match
# the previously hand-typed values exactly)
# =====================================================================================
axa = fig.add_subplot(g1[0])
_incr_a = {c: pd.read_csv(f"{DD}/merged_incr_{c.lower()}.csv") for c in COH}
_sig_a = {c: d[(d.fdr_uni < 0.05) & (d.incremental_r2_uni > 0)] for c, d in _incr_a.items()}
sig_n = {c: len(_sig_a[c]) for c in COH}
sig_pct = {c: 100.0 * len(_sig_a[c]) / int(_incr_a[c]["incremental_r2_uni"].notna().sum())
           for c in COH}
for c, y in zip(COH, ROWS):
    n = sig_n[c]
    axa.plot([0, n], [y, y], color=LGREY, lw=1.2, zorder=1, solid_capstyle="butt")
    axa.scatter(n, y, s=40, color=ORG[c], edgecolors=INK, linewidths=0.4, zorder=3)
    axa.text(n + 150, y, f"{n:,}\n({sig_pct[c]:.1f}%)", fontsize=7.0, color=INK,
             ha="left", va="center", linespacing=1.0)
axa.set_yticks(ROWS); axa.set_yticklabels([f"{c} ({ORGAN[c]})" for c in COH], fontsize=7.0)
axa.set_xlim(0, 3300); axa.set_ylim(-0.6, 4.6)
axa.set_xticks([0, 1000, 2000, 3000]); axa.set_xticklabels(["0", "1,000", "2,000", "3,000"], fontsize=7.0)
axa.tick_params(axis="y", length=0)
axa.spines["left"].set_visible(False); axa.spines["bottom"].set_bounds(0, 3000)
axa.set_xlabel("significant residual proteins\n(FDR < 0.05; % of tested proteins)", fontsize=7)
ta = axa.set_title("Morphology-predictable\nresidual proteins per organ", loc="left", fontsize=7.6, fontweight="bold")

# =====================================================================================
# b : leading enriched Reactome programme per organ (server-verified topq / topfam)
# =====================================================================================
axb = fig.add_subplot(g1[1])
# Same representative term per organ, same tested-proteome-background ORA (P0-5) and source
# file as make_composite2.py panel a; the earlier topq/topfam values here were leftover from
# the superseded whole-genome-background analysis (mathematically impossible under the
# current methodology -- CCRCC's best available -log10 q is ~2.1, not 16.3 -- and topfam
# contradicted panel a's family assignments). Fixed 2026-09-17.
_lead_term = {"CCRCC": "Chaperonin-mediated Protein Folding",
              "LUAD": "Cap-dependent Translation Initiation",
              "UCEC": "Asparagine N-linked Glycosylation",
              "GBM": "mRNA Splicing - Major Pathway",
              "PDAC": "Extracellular Matrix Organization"}
_lead_fam = {"CCRCC": "folding", "LUAD": "translation", "UCEC": "secretion",
             "GBM": "splicing", "PDAC": "ecm"}
_efull_b = pd.read_csv(f"{DD}/enrichment_tested_background_full.csv")
topq, topfam = {}, {}
for _c, _term in _lead_term.items():
    _row = _efull_b[(_efull_b.cohort == _c) & _efull_b.Term.str.startswith(_term)]
    assert len(_row) == 1, f"expected exactly one match for {_c}/{_term}, got {len(_row)}"
    topq[_c] = -np.log10(float(_row["Adjusted P-value"].iloc[0]))
    topfam[_c] = _lead_fam[_c]
XB = 40
for c, y in zip(COH, ROWS):
    q, f = topq[c], topfam[c]
    axb.plot([0, q], [y, y], color=LGREY, lw=1.2, zorder=1, solid_capstyle="butt")
    axb.scatter(q, y, s=40, color=FAM[f], edgecolors=INK, linewidths=0.4, zorder=3)
    axb.text(q + 1.5, y, f, fontsize=7.0, color=FAM[f], fontweight="bold", ha="left", va="center")
axb.set_yticks(ROWS); axb.set_yticklabels(COH, fontsize=7.0)
axb.set_xlim(0, XB); axb.set_ylim(-0.6, 4.6)
axb.set_xticks([0, 10, 20, 30]); axb.tick_params(axis="x", labelsize=7.0)
axb.tick_params(axis="y", length=0)
axb.spines["left"].set_visible(False); axb.spines["bottom"].set_bounds(0, 30)
axb.set_xlabel("top enriched Reactome term\n(\u2212log10 adj. p, tested-background ORA)", fontsize=7)
tb = axb.set_title("Leading enriched programme\nper organ", loc="left", fontsize=7.6, fontweight="bold")
fams = ["translation", "secretion", "splicing", "folding", "ecm"]
hb_ = [Line2D([], [], marker="o", ls="none", ms=4.2, mfc=FAM[f], mec=INK, mew=0.4) for f in fams]
# upper right: the only region of the panel with no row label (CCRCC 'splicing' ends at x ~ 37;
# the LUAD 'translation' label sits a full row lower)
axb.legend(hb_, fams, loc="upper right", fontsize=7.0, handlelength=0.9, handletextpad=0.4,
           labelspacing=0.25, borderaxespad=0.0, title="family", title_fontsize=7.0, ncol=2)

# =====================================================================================
# c : t-SNE of UNI slide embeddings, coloured by cohort
# =====================================================================================
axc = fig.add_subplot(g1[2])
emb = pd.read_csv(f"{DD}/emb_pca50.csv")
X = emb.filter(like="pc").values
coh = emb["cohort"].str.lower().map({c.lower(): c for c in COH}).values
n_slides = len(emb)
ts = TSNE(n_components=2, perplexity=30, init="pca", learning_rate="auto", random_state=0).fit_transform(X)
n_per = {}
for c in COH:
    m = coh == c; n_per[c] = int(m.sum())
    axc.scatter(ts[m, 0], ts[m, 1], s=5, color=ORG[c], alpha=0.75, edgecolors="none", label=c, zorder=3)
axc.set_xticks([]); axc.set_yticks([])
for sp in axc.spines.values():
    sp.set_visible(True); sp.set_linewidth(0.6)
axc.set_xlabel("t-SNE 1", fontsize=7); axc.set_ylabel("t-SNE 2", fontsize=7)
pad = 0.05 * np.ptp(ts, axis=0)
axc.set_xlim(ts[:, 0].min() - pad[0], ts[:, 0].max() + pad[0])
axc.set_ylim(ts[:, 1].min() - pad[1], ts[:, 1].max() + pad[1] * 3.2)   # head-room for the n text
axc.text(0.03, 0.97, f"n = {n_slides:,} slides", transform=axc.transAxes, fontsize=7.0,
         color=GREY, ha="left", va="top")
# 7 pt labels are wider than the 5.6 pt ones this legend was tuned for; the handles and the
# inter-column gaps absorb the difference so the strip still fits the panel it belongs to.
axc.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=5, fontsize=7.0, markerscale=1.7,
           handlelength=0.5, handletextpad=0.2, columnspacing=0.38, borderaxespad=0.0)
tc = axc.set_title("UNI slide embeddings\nseparate by organ", loc="left", fontsize=7.6, fontweight="bold")

# =====================================================================================
# d : per-gene agreement between encoders (CCRCC hexbin)
# =====================================================================================
axd = fig.add_subplot(g2[0])
mg = pd.read_csv(f"{DD}/merged_incr_ccrcc.csv")
dd = mg[["incremental_r2_uni", "incremental_r2_phi"]].replace([np.inf, -np.inf], np.nan).dropna()
xu, yp = dd["incremental_r2_uni"].values, dd["incremental_r2_phi"].values
r_ccrcc = float(np.corrcoef(xu, yp)[0, 1]); n_genes = len(dd)
EXT = (-0.35, 0.75)
n_out = int(((xu < EXT[0]) | (xu > EXT[1]) | (yp < EXT[0]) | (yp > EXT[1])).sum())
axd.axhline(0, color=LGREY, lw=0.5, zorder=1); axd.axvline(0, color=LGREY, lw=0.5, zorder=1)
hb = axd.hexbin(xu, yp, gridsize=45, extent=EXT + EXT, cmap="Blues", bins="log", mincnt=1,
                linewidths=0, zorder=2)
axd.plot(EXT, EXT, ls="--", color=GREY, lw=0.8, zorder=3)
axd.text(EXT[1] - 0.02, EXT[1] - 0.17, "identity", fontsize=7.0, color=GREY, ha="right", va="top", zorder=4)
axd.text(0.04, 0.96, f"r = {r_ccrcc:.2f}\nn = {n_genes:,} genes", transform=axd.transAxes,
         fontsize=7.0, ha="left", va="top", zorder=5)
if n_out:
    axd.text(0.96, 0.04, f"{n_out:,} genes outside\nthe plotted range", transform=axd.transAxes,
             fontsize=7.0, color=GREY, ha="right", va="bottom", linespacing=1.05)
axd.set_xlim(EXT); axd.set_ylim(EXT)
axd.set_xticks([-0.25, 0, 0.25, 0.5, 0.75]); axd.set_yticks([-0.25, 0, 0.25, 0.5, 0.75])
axd.tick_params(labelsize=7.0)
axd.set_xlabel("UNI incremental R\u00b2", fontsize=7)
axd.set_ylabel("Phikon incremental R\u00b2", fontsize=7)
cb = fig.colorbar(hb, ax=axd, fraction=0.046, pad=0.03)
cb.set_label("genes per hexagon", fontsize=7.0); cb.ax.tick_params(labelsize=7.0, length=2, width=0.5)
cb.outline.set_linewidth(0.5)
# plain integers instead of the default 10^n sci-notation: mathtext exponents render at
# ~0.7x the span size, which would put a 4-pt glyph in a figure floored at 7 pt.
cb.ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
cb.ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
td = axd.set_title("Per-gene agreement between encoders (CCRCC)", loc="left", fontsize=7.6, fontweight="bold")

# =====================================================================================
# e : agreement across cohorts (dumbbell: filled = per-gene r, open = recall)
# =====================================================================================
axe = fig.add_subplot(g2[1])
r_hard = {"CCRCC": 0.808, "LUAD": 0.780, "UCEC": 0.845, "GBM": 0.669, "PDAC": 0.854}
rec_hard = {"CCRCC": 0.737, "LUAD": 0.587, "UCEC": 0.742, "GBM": 0.380, "PDAC": 0.642}
ec = pd.read_csv(f"{DD}/encoder_comparison.csv")
ec["cohort"] = ec["cohort"].str.upper()
ph = ec[ec["encoder"] == "Phikon"].set_index("cohort")
have_csv = set(COH) <= set(ph.index) and {"r_incr_vs_uni", "recall_of_uni_sig"} <= set(ph.columns)
if have_csv:
    r_use = {c: float(ph.loc[c, "r_incr_vs_uni"]) for c in COH}
    rec_use = {c: float(ph.loc[c, "recall_of_uni_sig"]) for c in COH}
    src = "encoder_comparison.csv"
else:
    r_use, rec_use, src = r_hard, rec_hard, "hard-coded"
for c, y in zip(COH, ROWS):
    r, rc = r_use[c], rec_use[c]
    axe.plot([rc, r], [y, y], color=LGREY, lw=1.4, zorder=1, solid_capstyle="butt")
    axe.scatter(rc, y, s=30, facecolors="white", edgecolors=ORG[c], linewidths=1.0, zorder=3)
    axe.scatter(r, y, s=36, color=ORG[c], edgecolors=INK, linewidths=0.4, zorder=4)
    axe.text(rc - 0.03, y, f"{rc:.2f}", fontsize=7.0, color=GREY, ha="right", va="center")
    axe.text(r + 0.03, y, f"r {r:.2f}", fontsize=7.0, color=ORG[c], ha="left", va="center")
axe.set_yticks(ROWS); axe.set_yticklabels(COH, fontsize=7.0)
axe.set_xlim(0, 1.06); axe.set_ylim(-0.6, 4.6)
axe.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); axe.tick_params(axis="x", labelsize=7.0)
axe.tick_params(axis="y", length=0)
axe.spines["left"].set_visible(False); axe.spines["bottom"].set_bounds(0, 1)
axe.set_xlabel("agreement, UNI vs Phikon", fontsize=7)
axe.text(0.5, -0.175, "filled = per-gene r of incremental R\u00b2 (UNI vs Phikon)\n"
         "open = recall of UNI-significant proteins by Phikon",
         transform=axe.transAxes, ha="center", va="top", fontsize=7.0, color=GREY, linespacing=1.15)
te = axe.set_title("Agreement is consistent across cohorts", loc="left", fontsize=7.6, fontweight="bold")

# --- panel letters (after a draw pass) and layout diagnostics ---------------------------------
placed = place_letters(fig, [(axa, ta, "a"), (axb, tb, "b"), (axc, tc, "c"), (axd, td, "d"), (axe, te, "e")])
left_misalign_in = placed["a"][0] - placed["d"][0]          # >0: panel d content starts right of a's
bbox = fig.get_tightbbox(fig.canvas.get_renderer())          # inches; save_pub adds 0.1 in pad each side
out_w, out_h = bbox.width + 0.2, bbox.height + 0.2

S.save_pub(fig, "Fig_supp1_reproducibility")
agree = all(abs(r_use[c] - r_hard[c]) < 5e-4 and abs(rec_use[c] - rec_hard[c]) < 5e-4 for c in COH)
print(f"layout: a/d left-edge misalignment {left_misalign_in:+.3f} in; output canvas {out_w:.2f} x {out_h:.2f} in")
print(f"wrote Fig_supp1_reproducibility | a: sig counts {[sig_n[c] for c in COH]} | "
      f"b: top -log10q {[topq[c] for c in COH]} | c: t-SNE n={n_slides} slides "
      f"({', '.join(f'{c} {n_per[c]}' for c in COH)}) | d: CCRCC UNI-vs-Phikon r={r_ccrcc:.3f} "
      f"n={n_genes} genes, {n_out} outside hexbin extent | e: source={src}; "
      f"r csv={[round(r_use[c], 3) for c in COH]} hard={[r_hard[c] for c in COH]}; "
      f"recall csv={[round(rec_use[c], 3) for c in COH]} hard={[rec_hard[c] for c in COH]}; "
      f"agree={agree}")
