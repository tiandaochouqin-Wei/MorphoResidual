#!/usr/bin/env python3
"""Composite: "alternative explanations excluded" (Fig_controls).  Seven panels (eight
when the per-gene transcriptome files are present), all native, server-verified numbers,
no bar charts.

ROW 1  a  CCRCC incremental-R2 histogram with capacity-control reference lines
          (computed from figdata/merged_incr_ccrcc.csv; encoding of make_fig2.py)
       b  scanner census dot strip (make_figs.py Fig5 a)
ROW 2  c  batch-retention slope chart, symlog (make_figs.py Fig5 b)
       d  enriched-family survival heatmap (make_figs.py Fig5 c)
       e  tumour-purity scatter (make_purity_fig.py; figdata/confound_check_results.csv)
ROW 3  f  retained advantage vs transcriptome-wide baseline (make_txome_fig.py a)
       g  canonical-correlation dumbbell (make_txome_fig.py c)
       -- only if all five figdata/transcriptome_baseline_<cohort>.csv exist, row 3 is
          three panels instead: f scatter, g per-gene excess violins (make_txome_fig.py b),
          h canonical-correlation dumbbell.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke
from scipy import stats as st
import mrstyle as S

ORG, INK, GREY, LGREY, FAM = S.ORG, S.INK, S.GREY, S.LGREY, S.FAM
COH, ORGAN = S.COH, S.ORGAN
DD = "figdata"


def lab(ax, s, x=-0.18, y=1.05):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.92, 0.92], hspace=0.40,
                      left=0.075, right=0.975, top=0.975, bottom=0.05)
g1 = gs[0].subgridspec(1, 2, width_ratios=[58, 42], wspace=0.28)
g2 = gs[1].subgridspec(1, 3, width_ratios=[1.12, 1.0, 1.0], wspace=0.62)
# g3 (row 3) is built below: two panels, or three when the per-gene CSVs are present.

# =====================================================================================
# a : CCRCC incremental-R2 distribution with capacity controls (make_fig2.py encoding,
#     counts computed from figdata/merged_incr_ccrcc.csv instead of hardcoded)
# =====================================================================================
axa = fig.add_subplot(g1[0]); lab(axa, "a", x=-0.12)
mg = pd.read_csv(f"{DD}/merged_incr_ccrcc.csv")
v_all = mg["incremental_r2_uni"].replace([np.inf, -np.inf], np.nan).dropna().values
edges = np.round(np.arange(-0.70, 0.80 + 1e-9, 0.05), 4)           # 30 bins of 0.05
lo = edges[:-1]
cnt, _ = np.histogram(v_all, bins=edges)
n_below = int((v_all < -0.70).sum())
median = float(np.median(v_all))
frac_pos = float((v_all > 0).mean())
n_sig = int(((mg["fdr_uni"] < 0.05) & (mg["incremental_r2_uni"] > 0)).sum())
pct_sig = 100.0 * n_sig / len(v_all)
# capacity control on the significant gene set (median incremental R^2 across genes, from
# server_export/pinned/ccrcc/randproj_control_results.csv -- per-gene incr_real/incr_gauss/
# incr_randproj; median confirmed 2026-09-17 to match the previously hand-typed values here)
_cap_raw = pd.read_csv("../server_export/pinned/ccrcc/randproj_control_results.csv")
_cap_med = _cap_raw[["incr_real", "incr_gauss", "incr_randproj"]].median()
cap = {  # key: (value, colour, linestyle, linewidth)
    "real":   (float(_cap_med["incr_real"]),     FAM["translation"], "-",  1.5),
    "gauss":  (float(_cap_med["incr_gauss"]),    GREY,               "--", 1.3),
    "rproj":  (float(_cap_med["incr_randproj"]), LGREY,              ":",  2.0),
}
col = ["#D55E00" if x >= 0 else "#B9C0C7" for x in lo]
axa.bar(lo, cnt, width=0.048, color=col, align="edge", edgecolor="none", zorder=2)   # histogram
axa.axvline(0, color="k", lw=0.7, zorder=3)
axa.axvline(median, color="#333", lw=0.9, ls="--", zorder=3)
axa.text(median - 0.015, 1010, f"median {median:.2f}".replace("-", "−"),
         ha="right", va="bottom", fontsize=7.0, color="#333", zorder=6,
         bbox=dict(fc="white", ec="none", pad=0.8))
axa.text(0.26, 780, f"{100*frac_pos:.1f}% of genes positive;\n{n_sig:,} ({pct_sig:.1f}%)\nFDR-significant",
         fontsize=7.0, color="#D55E00", va="top")
axa.text(-0.695, 300, f"{n_below} genes\n< −0.70\nnot shown", fontsize=7.0, color="#8C8C8C",
         va="bottom", ha="left", linespacing=1.1)
YTOP = 1400
specs = {  # key: (line top, label text, label x-offset, ha)
    "real":  (1000, "real WSI PCs +0.09",        +0.012, "left"),
    "gauss": (880,  "Gaussian noise −0.29",   -0.003, "right"),
    "rproj": (1130, "random projection −0.08", -0.012, "right"),
}
for k, (v, c, ls, lw) in cap.items():
    ytop, labt, dx, ha = specs[k]
    ln = axa.plot([v, v], [0, ytop], color=c, ls=ls, lw=lw, zorder=4, solid_capstyle="butt")[0]
    if k == "rproj":
        ln.set_path_effects([withStroke(linewidth=lw + 1.8, foreground=INK)])
    tcol = INK if k == "rproj" else c
    axa.text(v + dx, ytop + 20, labt, ha=ha, va="bottom", fontsize=7.0, color=tcol,
             fontweight="bold" if k == "real" else "normal", zorder=6,
             bbox=dict(fc="white", ec="none", pad=0.8))
axa.text(0.74, YTOP, "reference lines: median incremental R²\non the FDR-significant set,\n"
         "per feature set (capacity control)", ha="right", va="top", fontsize=7.0,
         color=GREY, linespacing=1.05)
axa.set_xlabel(r"Incremental R² of morphology over mRNA", fontsize=7)
axa.set_ylabel("Genes", fontsize=7)
axa.set_xlim(-0.72, 0.75); axa.set_ylim(0, YTOP)
axa.set_yticks([0, 200, 400, 600, 800, 1000, 1200]); axa.tick_params(labelsize=7.0)
axa.set_title("Most genes show no gain; the positive tail is\na specific, non-capacity effect",
              loc="left", fontsize=7.6, fontweight="bold")

# =====================================================================================
# b : scanner census dot strip (make_figs.py Fig5 a)
# =====================================================================================
axb = fig.add_subplot(g1[1]); lab(axb, "b", x=-0.05, y=1.03)   # single-line title: letter baseline on the title baseline
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
scan_pct = _BL.scan_pct.to_numpy(float)
n_scan = _BL.n_scanners.tolist()
rows = np.arange(5)[::-1]
for i, ch in enumerate(COH):
    v, y, c = scan_pct[i], rows[i], ORG[ch]
    axb.plot([80, 100], [y, y], ls=":", lw=0.5, color="#D9D9D9", zorder=1)
    if ch == "PDAC":
        axb.plot(v, y, "o", ms=9, mfc="none", mec=c, mew=0.9, zorder=3)
        note = f"two scanners ({100 - v:.1f}% on the second)"
    elif n_scan[i] == 2:
        note = f"two scanners ({100 - v:.1f}% minority)"
    else:
        note = "single scanner"
    axb.plot(v, y, "o", ms=4.6, color=c, mec="white", mew=0.5, zorder=4)
    axb.text(v - 1.1, y + 0.12, f"{ch}  {v:.0f}%", ha="right", va="bottom", fontsize=7.0,
             color=c, fontweight="bold")
    axb.text(v - 1.1, y - 0.12, note, ha="right", va="top", fontsize=7.0, color="#8C8C8C")
axb.set_xlim(80, 101.5); axb.set_ylim(-0.75, 4.75)
axb.set_xticks([80, 85, 90, 95, 100]); axb.tick_params(axis="x", labelsize=7.0)
axb.set_yticks([]); axb.spines["left"].set_visible(False); axb.spines["bottom"].set_bounds(80, 100)
axb.set_xlabel("Slides on dominant scanner (%)", fontsize=7)
axb.set_title("Scanner $\\neq$ driver", loc="left", fontsize=7.6, fontweight="bold")

# =====================================================================================
# c : batch-retention slope chart (make_figs.py Fig5 b)
# =====================================================================================
axc = fig.add_subplot(g2[0]); lab(axc, "c", x=-0.30)
sig_n = _BL.baseline.to_numpy(int)
batchcorr = _BL.batchcorr.to_numpy(int)
strictest = _BL.strictest.to_numpy(int)
stages = ["baseline", "batch-\ncorrected", "strictest\n(grouped\nCV)"]
xs = np.array([0, 1, 2])
counts = np.c_[sig_n, batchcorr, strictest]
axc.set_yscale("symlog", linthresh=10, linscale=0.35)
axc.set_ylim(0, 4200); axc.set_xlim(-0.12, 3.05)
axc.set_yticks([0, 10, 100, 1000]); axc.set_yticklabels(["0", "10", "100", "1000"], fontsize=7.0)
axc.yaxis.set_minor_locator(mpl.ticker.FixedLocator(
    [k * 10 ** e for e in (1, 2, 3) for k in range(2, 10)]))
axc.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
axc.plot([-0.12, 2.0], [50, 50], ls="--", lw=0.7, color="k", zorder=1)
axc.text(2.06, 50, "50", fontsize=7.0, color="k", ha="left", va="center")
axc.spines["bottom"].set_bounds(0, 2)
for i, ch in enumerate(COH):
    yv = counts[i].astype(float); c = ORG[ch]
    axc.plot(xs, yv, "-", lw=1.1, color=c, zorder=3)
    filled = yv > 0
    axc.plot(xs[filled], yv[filled], "o", ms=3.6, color=c, mec="white", mew=0.5, zorder=4)
    if (~filled).any():
        axc.plot(xs[~filled], yv[~filled], "o", ms=3.8, mfc="white", mec=c, mew=0.9,
                 zorder=5, clip_on=False)


def _repel(ax, ys, min_pt=8.0):
    """Nudge label y-positions apart (display space) so end-labels never collide."""
    ys = np.asarray(ys, float); order = np.argsort(ys)
    disp = ax.transData.transform(np.c_[np.zeros(len(ys)), ys])[:, 1]
    gap = min_pt / 72 * ax.figure.dpi
    dd_ = disp[order].copy()
    for k in range(1, len(dd_)):
        dd_[k] = max(dd_[k], dd_[k - 1] + gap)
    out = np.empty_like(ys)
    out[order] = ax.transData.inverted().transform(np.c_[np.zeros(len(dd_)), dd_])[:, 1]
    return out


yl = _repel(axc, strictest)
for i, ch in enumerate(COH):
    axc.text(2.09, yl[i], f"{ch}  {strictest[i]:,}", fontsize=7.0, color=ORG[ch],
             ha="left", va="center", fontweight="bold")
axc.set_xticks(xs); axc.set_xticklabels(stages, fontsize=7.0, linespacing=0.95)
axc.set_ylabel("Significant proteins\n(log scale; linear below 10)", fontsize=7)
axc.set_title("Signal survives batch", loc="left", fontsize=7.6, fontweight="bold")

# =====================================================================================
# d : enrichment-family survival heatmap (make_figs.py Fig5 c)
# =====================================================================================
axd = fig.add_subplot(g2[1]); lab(axd, "d", x=-0.36)
sigterms = _BL[["sigterms_baseline", "sigterms_batchcorr",
                "sigterms_strictest"]].to_numpy(int)
im = axd.imshow(sigterms, cmap="YlGnBu", vmin=0, vmax=15, aspect="auto")
axd.set_xticks([0, 1, 2]); axd.set_xticklabels(["base", "batch", "strict"], fontsize=7.0)
axd.set_yticks(range(5)); axd.set_yticklabels(COH, fontsize=7.0)
for i in range(5):
    for j in range(3):
        vv = sigterms[i, j]
        axd.text(j, i, str(vv), ha="center", va="center",
                 color="white" if vv > 8 else "black", fontsize=7)
axd.set_title("Biology survives", loc="left", fontsize=7.6, fontweight="bold")
cb = fig.colorbar(im, ax=axd, fraction=0.05, pad=0.05)
cb.set_ticks([0, 5, 10, 15])          # four integer ticks: legible at 7 pt, no decimal clutter
cb.set_label("enriched terms\n(of top 15)", fontsize=7.0, labelpad=2)
cb.ax.tick_params(labelsize=7.0)

# =====================================================================================
# e : tumour-purity control scatter (make_purity_fig.py)
# =====================================================================================
axe = fig.add_subplot(g2[2]); lab(axe, "e", x=-0.36)
pc = pd.read_csv(f"{DD}/confound_check_results.csv")
pc = pc.replace([np.inf, -np.inf], np.nan).dropna(subset=["incremental_r2", "incr_r2_purity_adjusted", "fdr"])
psig = pc[(pc.fdr < 0.05) & (pc.incremental_r2 > 0) & (pc.incremental_r2 < 1.5)
          & (pc.incr_r2_purity_adjusted.abs() < 1.5)].copy()
px = psig.incremental_r2.values; py = psig.incr_r2_purity_adjusted.values
pr = st.pearsonr(px, py)[0]
ratio = np.median(py / px); pfrac_pos = np.mean(py > 0)
axe.scatter(px, py, s=5, color=ORG["CCRCC"], alpha=0.45, edgecolors="none")
plim = [0, max(0.72, max(px.max(), py.max()) * 1.03)]      # axis and identity line both reach ~0.72
axe.plot(plim, plim, ls="--", color=GREY, lw=0.9)
# stats block in the empty lower-right region (clear of the identity line and the outliers)
axe.text(0.97, 0.04, f"n = {len(psig):,} proteins\nr = {pr:.2f}\n"
         f"median retained {100*ratio:.0f}%\n{100*pfrac_pos:.0f}% remain > 0",
         transform=axe.transAxes, fontsize=7.0, ha="right", va="bottom", linespacing=1.25)
axe.set_xlabel(r"morphology incremental R² (mRNA baseline)", fontsize=7.0)
axe.set_ylabel(r"incremental R² after purity adjustment", fontsize=7.0)
axe.set_xlim(plim); axe.set_ylim(min(0, py.min()) - 0.02, plim[1])
axe.set_xticks([0, 0.2, 0.4, 0.6])
axe.tick_params(labelsize=7.0)
axe.set_title("Not a tumour-purity artefact (CCRCC)", fontsize=7.6, fontweight="bold", loc="left")

# =====================================================================================
# f, g (, h) : transcriptome-wide baseline (make_txome_fig.py a, c; b only when the five
# per-gene CSVs figdata/transcriptome_baseline_<cohort>.csv are all present)
# =====================================================================================
tx = pd.read_csv(f"{DD}/transcriptome_baseline_summary.csv").set_index("cohort").loc[COH]
per = {c: f"{DD}/transcriptome_baseline_{c.lower()}.csv" for c in COH}
have = all(os.path.exists(p) for p in per.values())
if have:   # three panels: f scatter, g per-gene violins, h dumbbell
    g3 = gs[2].subgridspec(1, 3, width_ratios=[1.0, 1.15, 0.95], wspace=0.60)
    axf = fig.add_subplot(g3[0]); lab(axf, "f", x=-0.32)
    axg = fig.add_subplot(g3[1]); lab(axg, "g", x=-0.25)
    axcc = fig.add_subplot(g3[2]); lab(axcc, "h", x=-0.30)
else:      # two panels: f scatter (~52 % of the row), g dumbbell (~48 %)
    g3 = gs[2].subgridspec(1, 2, width_ratios=[52, 48], wspace=0.38)
    axf = fig.add_subplot(g3[0]); lab(axf, "f", x=-0.16)
    axg = None
    axcc = fig.add_subplot(g3[1]); lab(axcc, "g", x=-0.17)

# ---- f: retained advantage, cohort-level ----
lim = 0.48
axf.plot([0, lim], [0, lim], color=GREY, lw=0.8, ls="--")
axf.plot([0, lim], [0, 0.5 * lim], color=LGREY, lw=0.6, ls=":")
axf.plot([0, lim], [0, 0.25 * lim], color=LGREY, lw=0.6, ls=":")
axf.text(lim - 0.01, lim - 0.01, "100%", fontsize=7.0, color=GREY, ha="right", va="top")
axf.text(0.13, 0.5 * 0.13 + 0.008, "50%", fontsize=7.0, color=GREY, ha="left", va="bottom")
axf.text(0.13, 0.25 * 0.13 - 0.008, "25%", fontsize=7.0, color=GREY, ha="left", va="top")
for c in COH:
    r = tx.loc[c]
    axf.scatter(r.adv_own, r.adv_pcs20, s=34, color=ORG[c], edgecolors=INK, linewidths=0.4, zorder=4)
    off = {"UCEC": (0.014, 0.0, "left", "center"), "CCRCC": (0.014, 0.0, "left", "center"),
           "GBM": (0.014, 0.018, "left", "center"), "LUAD": (-0.014, -0.004, "right", "center"),
           "PDAC": (0.0, -0.028, "center", "top")}[c]
    axf.text(r.adv_own + off[0], r.adv_pcs20 + off[1], f"{c} {r.retains_pct:.0f}%",
             fontsize=7.0, color=ORG[c], ha=off[2], va=off[3])
axf.set_xlim(0, lim); axf.set_ylim(-0.02, lim)
axf.tick_params(labelsize=7.0)
axf.set_xlabel("advantage over noise,\nbaseline = own mRNA", fontsize=7)
axf.set_ylabel("advantage over noise,\nbaseline = own mRNA + 20 RNA PCs", fontsize=7)
axf.set_title("Retained against the transcriptome", fontsize=7.6, fontweight="bold", loc="left")

# ---- g (three-panel branch only): per-gene excess distributions ----
if have:
    data, pos = [], []
    for c in COH:
        gg = pd.read_csv(per[c])
        ex = (gg["incr_over_pcs20"] - gg["incr_over_pcs20_randctrl"]).dropna().values
        pos.append(100 * (ex > 0).mean())          # share >0 on ALL genes
        lo, hi = np.percentile(ex, [1, 99])         # violin drawn on the central 98 % so no tail is clipped
        data.append(ex[(ex >= lo) & (ex <= hi)])
    vp = axg.violinplot(data, positions=range(5), widths=0.8, showextrema=False, showmedians=False)
    for body, c in zip(vp["bodies"], COH):
        body.set_facecolor(ORG[c]); body.set_edgecolor("none"); body.set_alpha(0.75)
    for i, ex in enumerate(data):
        axg.scatter(i, np.median(ex), marker="_", s=180, linewidths=1.5, color=INK, zorder=5)
    axg.axhline(0, color="k", lw=0.5)
    # share of genes beating noise goes into the tick label (a free-floating text at the 99th
    # percentile escaped the axes into panel d when the tail exceeded the y-limit)
    axg.set_xticks(range(5))
    # 7 pt is wider than one violin slot, so the cohort and its share go on one 45-degree
    # line (the ">0" is spelled out in the note below) instead of two colliding rows
    axg.set_xticklabels([f"{c} {p:.0f}%" for c, p in zip(COH, pos)], fontsize=7.0,
                        rotation=30, ha="right", rotation_mode="anchor")
    axg.set_ylabel("per-gene excess of real morphology\nover same-size noise (ΔR², txome baseline)", fontsize=7.0)
    ymin = min(d.min() for d in data) - 0.03; ymax = max(d.max() for d in data) + 0.03
    axg.set_ylim(ymin, ymax)
    axg.text(0.99, 0.02, "violins: central 98% of genes\n% = genes with excess > 0",
             transform=axg.transAxes, ha="right", va="bottom", fontsize=7.0,
             color=GREY, linespacing=1.25)
    axg.set_title("Gene-by-gene excess over noise", fontsize=7.6, fontweight="bold", loc="left")

# ---- g (or h): collinearity with the RNA PCs, morphology vs noise (dumbbell) ----
for i, c in enumerate(COH):
    r = tx.loc[c]; y = 4 - i
    axcc.plot([r.cc_noise_med, r.cc_morph_med], [y, y], color=LGREY, lw=1.4, zorder=1)
    axcc.scatter(r.cc_noise_med, y, s=26, facecolors="white", edgecolors=ORG[c], linewidths=1.0, zorder=3)
    axcc.scatter(r.cc_morph_med, y, s=30, color=ORG[c], edgecolors=INK, linewidths=0.4, zorder=4)
    # at 7 pt this count no longer fits to the right of the filled marker without running
    # off the canvas, so it sits above its own dumbbell instead
    axcc.text(0.285, y + 0.28, f"{int(r.cc_morph_gt05)} vs {int(r.cc_noise_gt05)} of 20 > 0.5",
              fontsize=7.0, color=GREY, va="bottom", ha="left")
axcc.set_yticks(range(5)); axcc.set_yticklabels(COH[::-1], fontsize=7.0)
axcc.set_xlim(0.25, 0.78); axcc.set_ylim(-1.05, 4.65); axcc.tick_params(axis="x", labelsize=7.0)
axcc.set_xlabel("median canonical correlation\nwith the 20 RNA PCs", fontsize=7)
axcc.set_title("Collinear with the RNA PCs\n(not a directional bound)", fontsize=7.6, fontweight="bold", loc="left")
# the marker key sits inside the axes (below the PDAC row): at 7 pt a line under the
# two-line x-label ran off the bottom of the canvas
axcc.text(0.5, 0.015, "filled = morphology PCs\nopen = same-size noise", transform=axcc.transAxes,
          ha="center", va="bottom", fontsize=7.0, color=GREY, linespacing=1.25)

S.save_pub(fig, "Fig_controls")
print(f"wrote Fig_controls | hist n={len(v_all)} median={median:.4f} pos={100*frac_pos:.1f}% "
      f"sig={n_sig} ({pct_sig:.1f}%) below-0.70={n_below} | purity n={len(psig)} r={pr:.3f} "
      f"retained={100*ratio:.0f}% >0={100*pfrac_pos:.0f}% | row 3 = "
      f"{'three panels f-g-h (per-gene violins drawn)' if have else 'two panels f-g (per-gene CSVs absent)'}")
