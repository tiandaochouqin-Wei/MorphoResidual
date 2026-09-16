#!/usr/bin/env python3
"""MorphoResidual figures 4 & 5 from summary numbers (server-verified 2026-09-01)."""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import mrstyle as S  # Arial, editable-text SVG, 600-dpi TIFF (submission grade)

cohorts = ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]
organs  = ["Kidney", "Lung", "Uterus", "Brain", "Pancreas"]
# --- counts, scanner census and enrichment survival: read, never hardcode ---------
# Regenerate with recount_batch_levels.py. The baseline row is asserted against the
# values Table 1 quotes, so a re-run that moves it breaks the build rather than
# silently drawing a figure the text contradicts.
import pandas as _pd
_BL = _pd.read_csv("figdata/batch_levels.csv").set_index("cohort").loc[
    ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]]
assert list(_BL.index) == ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]
assert list(_BL.baseline) == [2191, 2310, 2566, 702, 1572], (
    f"baseline counts moved: {list(_BL.baseline)} -- update Table 1 before the figures")
sig       = _BL.baseline.to_numpy(int)
pct       = _BL.pct_of_tested.to_numpy(float)
batchcorr = _BL.batchcorr.to_numpy(int)
strictest = _BL.strictest.to_numpy(int)
scan_pct  = _BL.scan_pct.to_numpy(float)
n_scan    = _BL.n_scanners.tolist()
# signature terms (of top-15) at baseline / batch-corrected / strictest
sigterms  = _BL[["sigterms_baseline", "sigterms_batchcorr",
                 "sigterms_strictest"]].to_numpy(int)
# top signature-family Reactome term at baseline: (-log10 q, family)
topq   = np.array([16.3, 50.7, 41.5, 14.6, 24.3])
topfam = ["splicing", "translation", "secretion", "splicing", "secretion"]

# Okabe-Ito colourblind-safe
ORG = {"CCRCC": "#0072B2", "LUAD": "#009E73", "UCEC": "#CC79A7",
       "GBM": "#E69F00", "PDAC": "#D55E00"}
FAM = {"splicing": "#0072B2", "translation": "#009E73",
       "secretion": "#E69F00", "ptm": "#CC79A7"}
ocol = [ORG[c] for c in cohorts]

# ===================== FIGURE 4 — five-organ reproducibility =====================
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.9))
x = np.arange(5)

a.bar(x, sig, color=ocol, width=0.68, edgecolor="none")
for i, (s, p) in enumerate(zip(sig, pct)):
    a.text(i, s + 40, f"{s:,}\n{p:.1f}%", ha="center", va="bottom", fontsize=6.6, linespacing=0.95)
a.set_xticks(x); a.set_xticklabels([f"{c}\n({o})" for c, o in zip(cohorts, organs)], fontsize=6.8)
a.set_ylabel("Morphology-predictable\nresidual proteins (FDR < 0.05)")
a.set_ylim(0, 3050)
a.set_title("a  Reproducible across five organs", loc="left", fontsize=9, fontweight="bold")

famcol = [FAM[f] for f in topfam]
b.barh(x, topq, color=famcol, height=0.66, edgecolor="none")
for i, (q, f) in enumerate(zip(topq, topfam)):
    b.text(q - 1.2, i, f, ha="right", va="center", color="white", fontsize=6.6, fontweight="bold")
b.set_yticks(x); b.set_yticklabels(cohorts, fontsize=7)
b.invert_yaxis()
b.set_xlabel(r"Top enriched Reactome term  ($-\log_{10}$ adj. $p$)")
b.set_title("b  Convergent programme", loc="left", fontsize=9, fontweight="bold")
leg = [Patch(fc=FAM[k], ec="none", label=k) for k in ["translation", "secretion", "splicing", "ptm"]]
b.legend(handles=leg, fontsize=6.2, frameon=False, loc="lower right", handlelength=1.0, ncol=2)
fig.tight_layout()
S.save_pub(fig, "Fig4_reproducibility")
plt.close(fig)

# ===================== FIGURE 5 — acquisition robustness =====================
fig = plt.figure(figsize=(7.4, 3.1))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.18], wspace=0.55,
                      left=0.075, right=0.955, top=0.82, bottom=0.18)
a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1]); c = fig.add_subplot(gs[2])
TT = dict(loc="left", fontsize=8, fontweight="bold", pad=6)

# a — scanner census as a dot strip: one dot per cohort on a single 80-100 % axis,
# cohorts staggered on their own row so 99/100/100/100 do not pile up
rows = np.arange(5)[::-1]                       # CCRCC on top ... PDAC at bottom
for i, ch in enumerate(cohorts):
    v, y, col = scan_pct[i], rows[i], ORG[ch]
    a.plot([80, 100], [y, y], ls=":", lw=0.5, color="#D9D9D9", zorder=1)   # faint row guide
    if ch == "PDAC":     # the only genuine two-scanner contrast -> ringed dot
        a.plot(v, y, "o", ms=9, mfc="none", mec=col, mew=0.9, zorder=3)
        note = f"two scanners ({100 - v:.1f}% on the second)"
    elif n_scan[i] == 2:  # CCRCC: nominally two scanners, but a 0.8% minority
        note = f"two scanners ({100 - v:.1f}% minority)"
    else:
        note = "single scanner"
    a.plot(v, y, "o", ms=4.6, color=col, mec="white", mew=0.5, zorder=4)
    a.text(v - 1.1, y + 0.12, f"{ch}  {v:.0f}%", ha="right", va="bottom", fontsize=6.4,
           color=col, fontweight="bold")
    a.text(v - 1.1, y - 0.12, note, ha="right", va="top", fontsize=5.6, color="#8C8C8C")
a.set_xlim(80, 101.5); a.set_ylim(-0.75, 4.75)
a.set_xticks([80, 85, 90, 95, 100]); a.tick_params(axis="x", labelsize=6.2)
a.set_yticks([]); a.spines["left"].set_visible(False); a.spines["bottom"].set_bounds(80, 100)
a.set_xlabel("Slides on dominant scanner (%)")
a.set_title("a  Scanner $\\neq$ driver", **TT)

# b — slope chart: significant proteins across three levels of batch control
# (symlog: log above 10, linear 0-10, so PDAC's collapse to 0 is a real axis position)
stages = ["baseline", "batch-\ncorrected", "strictest\n(grouped CV)"]
xs = np.array([0, 1, 2])
counts = np.c_[sig, batchcorr, strictest]          # rows = cohorts, cols = stages
b.set_yscale("symlog", linthresh=10, linscale=0.35)
b.set_ylim(0, 4200); b.set_xlim(-0.12, 3.05)
b.set_yticks([0, 10, 100, 1000]); b.set_yticklabels(["0", "10", "100", "1000"])
b.yaxis.set_minor_locator(mpl.ticker.FixedLocator(
    [k * 10 ** e for e in (1, 2, 3) for k in range(2, 10)]))
b.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
b.plot([-0.12, 2.0], [50, 50], ls="--", lw=0.7, color="k", zorder=1)   # retention threshold
b.text(2.06, 50, "50", fontsize=6, color="k", ha="left", va="center")
b.spines["bottom"].set_bounds(0, 2)
for i, ch in enumerate(cohorts):
    yv = counts[i].astype(float); col = ORG[ch]
    b.plot(xs, yv, "-", lw=1.1, color=col, zorder=3)
    filled = yv > 0
    b.plot(xs[filled], yv[filled], "o", ms=3.6, color=col, mec="white", mew=0.5, zorder=4)
    if (~filled).any():   # open marker at the true 0 (end-label spells out "0")
        b.plot(xs[~filled], yv[~filled], "o", ms=3.8, mfc="white", mec=col, mew=0.9,
               zorder=5, clip_on=False)

def _repel(ax, ys, min_pt=7.2):
    """Nudge label y-positions apart (display space) so end-labels never collide."""
    ys = np.asarray(ys, float); order = np.argsort(ys)
    disp = ax.transData.transform(np.c_[np.zeros(len(ys)), ys])[:, 1]
    gap = min_pt / 72 * ax.figure.dpi
    d = disp[order].copy()
    for k in range(1, len(d)):
        d[k] = max(d[k], d[k - 1] + gap)
    out = np.empty_like(ys)
    out[order] = ax.transData.inverted().transform(np.c_[np.zeros(len(d)), d])[:, 1]
    return out

yl = _repel(b, strictest)
for i, ch in enumerate(cohorts):
    b.text(2.09, yl[i], f"{ch}  {strictest[i]:,}", fontsize=6.4, color=ORG[ch],
           ha="left", va="center", fontweight="bold")
b.set_xticks(xs); b.set_xticklabels(stages, fontsize=6.2, linespacing=0.95)
b.set_ylabel("Significant proteins\n(log scale; linear below 10)")
b.set_title("b  Signal survives batch", **TT)

# c — enrichment-family survival heatmap (signature terms of top 15)
im = c.imshow(sigterms, cmap="YlGnBu", vmin=0, vmax=15, aspect="auto")
c.set_xticks([0, 1, 2]); c.set_xticklabels(["base", "batch", "strict"], fontsize=6.4)
c.set_yticks(x); c.set_yticklabels(cohorts, fontsize=6.8)
for i in range(5):
    for j in range(3):
        v = sigterms[i, j]
        c.text(j, i, str(v), ha="center", va="center",
               color="white" if v > 8 else "black", fontsize=7)
c.set_title("c  Biology survives", **TT)
cb = fig.colorbar(im, ax=c, fraction=0.05, pad=0.05)
cb.set_label("enriched signature\nterms (of top 15)", fontsize=5.6)
cb.ax.tick_params(labelsize=5.6)
S.save_pub(fig, "Fig5_robustness")
plt.close(fig)
print("wrote Fig4_reproducibility and Fig5_robustness")

# ===================== FIGURE 3 — per-cohort Reactome enrichment =====================
ECMc = "#8C8C8C"
FAM2 = dict(FAM); FAM2["ecm"] = ECMc
enr = {  # cohort: (terms, -log10 q, families)  baseline top-5 signature/leading terms
 "CCRCC": (["Capped-intron pre-mRNA proc.", "Membrane trafficking",
            "mRNA splicing (major)", "Post-transl. modification", "Protein localization"],
           [16.33, 15.89, 13.44, 13.32, 12.48],
           ["splicing", "secretion", "splicing", "ptm", "secretion"]),
 "LUAD":  (["rRNA processing", "Translation", "L13a transl. silencing",
            "60S subunit joining", "SRP co-transl. targeting"],
           [50.69, 48.28, 47.22, 45.56, 35.72],
           ["translation"] * 5),
 "UCEC":  (["Membrane trafficking", "N-linked glycosylation",
            "Capped-intron pre-mRNA proc.", "Vesicle-mediated transport", "mRNA splicing (major)"],
           [41.54, 39.35, 36.88, 31.68, 31.50],
           ["secretion", "secretion", "splicing", "secretion", "splicing"]),
 "GBM":   (["mRNA splicing (major)", "Capped-intron pre-mRNA proc.", "SUMO E3 ligases",
            "SUMOylation", "mRNA splicing"],
           [14.61, 14.61, 11.70, 11.42, 10.99],
           ["splicing", "splicing", "ptm", "ptm", "splicing"]),
 "PDAC":  (["Neutrophil degranulation", "Membrane trafficking", "ECM organization",
            "Vesicle-mediated transport", "ECM proteoglycans"],
           [24.31, 23.27, 22.19, 21.87, 15.39],
           ["secretion", "secretion", "ecm", "secretion", "ecm"]),
}
fig, axes = plt.subplots(5, 1, figsize=(4.9, 6.6), sharex=False)
for ax, ch in zip(axes, cohorts):
    terms, q, fams = enr[ch]
    yy = np.arange(len(terms))[::-1]
    ax.barh(yy, q, color=[FAM2[f] for f in fams], height=0.66, edgecolor="none")
    ax.set_yticks(yy); ax.set_yticklabels(terms, fontsize=6.2)
    ax.set_xlim(0, max(q) * 1.06)
    ax.tick_params(axis="x", labelsize=6)
    ax.set_title(f"{ch} ({dict(zip(cohorts, organs))[ch]})", loc="left",
                 fontsize=8, fontweight="bold", pad=2, color=ORG[ch])
    ax.margins(y=0.08)
axes[-1].set_xlabel(r"Reactome enrichment  ($-\log_{10}$ adjusted $p$)", fontsize=7.5)
leg = [Patch(fc=FAM2[k], label=k) for k in ["translation", "secretion", "splicing", "ptm", "ecm"]]
fig.legend(handles=leg, fontsize=6.4, frameon=False, ncol=5,
           loc="lower center", bbox_to_anchor=(0.5, -0.005), handlelength=1.0, columnspacing=1.1)
fig.suptitle("Enrichment of the post-transcriptional programme",
             fontsize=8.5, fontweight="bold", y=0.996)
fig.tight_layout(rect=(0, 0.03, 1, 0.975))
S.save_pub(fig, "Fig3_enrichment_legacy_bars")  # superseded by make_fig3_dotmatrix.py (2026-09-07); kept only for provenance
plt.close(fig)
print("wrote Fig3_enrichment")
