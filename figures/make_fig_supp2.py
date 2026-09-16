#!/usr/bin/env python3
"""Supplementary composite Fig_supp2_scores: composition controls, morphology-only residual
scores, nuclear correlates and clinical correlates.  Seven panels, three rows, no bar charts.

ROW 1  a  composition-score control, Cleveland dot plot (server-verified retained-% numbers)
       b  morphology-only vs measured translation residual, pooled scatter (figdata/scores_*.csv)
       c  same for the ER-secretion residual
ROW 2  d  nuclear-feature Spearman rho heatmap, translation residual (server-verified pathomics_link)
       e  same for the ER-secretion residual (shared colourbar)
ROW 3  f  Spearman rho with tumour grade, dot-and-interval (Fisher-z 95% CI) for the four gradeable cohorts
       g  CCRCC overall-survival forest: univariable vs grade+stage-adjusted translation score, and stage

Layout notes: panel letters are placed from the rendered bounding boxes (left column a/d/f share one x,
each row shares one baseline), and the outer margins are fitted so that the tight-bbox export is exactly
7.2 x 9.0 in (183 x 229 mm) -- nothing spills outside the canvas, so no post-hoc down-scaling.
"""
import json as _json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats as st
import mrstyle as S

# mathtext in Arial so italic statistical symbols (r, P, n, q, rho, z) match the body font
mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                     "mathtext.bf": "Arial:bold"})

ORG, INK, GREY, LGREY, FAM = S.ORG, S.INK, S.GREY, S.LGREY, S.FAM
COH = S.COH
DD = "figdata"
FS_TITLE = 7.6
FS_TICK = 7.0
FS_LEG = 7.0      # every legend / note / footnote
FS_STAT = 7.0     # r/P/n corner text, HR/P forest text
W_IN, H_IN = 7.2, 9.0
PAD_IN = 0.10     # content margin inside the canvas == savefig pad_inches, so the tight export is W_IN x H_IN

PANELS = []       # (letter, axes, row index, column group) -- letters are placed after layout, see place_letters()


def lab(ax, s, row, colgrp=None):
    PANELS.append((s, ax, row, colgrp))


def fmt_p(p):
    """P value: 3 dp when >= 1e-3 (0.007, 0.269); otherwise two-significant-digit mantissa in e-notation
    (4.93e-49 -> '4.9e-49', 7.1e-4 -> '7.1e-4').  Used identically in panels b, c and g."""
    if p >= 1e-3:
        return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e")
    return f"{m}e{int(e)}"


fig = plt.figure(figsize=(W_IN, H_IN))
# One top-level gridspec per row so that each row's left/right margins can be fitted to its own label
# widths (row 2 carries the wide feature labels and the colourbar label; rows 1 and 3 do not).
H_RATIOS, HSPACE = [1.0, 0.88, 0.98], 0.50
TOP, BOTTOM = 0.955, 0.06


def row_bounds(top, bottom):
    """Row (top, bottom) pairs in figure fraction, same arithmetic as a 3-row GridSpec."""
    cell = (top - bottom) / (len(H_RATIOS) + HSPACE * (len(H_RATIOS) - 1))
    hs = [cell * len(H_RATIOS) * r / sum(H_RATIOS) for r in H_RATIOS]
    out, t = [], top
    for h in hs:
        out.append((t, t - h)); t = t - h - HSPACE * cell
    return out


(_t1, _b1), (_t2, _b2), (_t3, _b3) = row_bounds(TOP, BOTTOM)
g1 = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.05, 1.05], wspace=0.42, left=0.10, right=0.97, top=_t1, bottom=_b1)
g2 = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.05], wspace=0.22, left=0.14, right=0.94, top=_t2, bottom=_b2)
g3 = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.48, left=0.10, right=0.97, top=_t3, bottom=_b3)

# =====================================================================================
# a : composition controls -- % of the morphology increment retained over mRNA + score(s)
# =====================================================================================
axa = fig.add_subplot(g1[0]); lab(axa, "a", 0, "L")
# Series colours deliberately avoid every cohort (S.ORG) and family (S.FAM) hue used elsewhere in the
# figure, and each series also gets its own marker shape (colour-blind safe).
# Computed from the per-gene composition_controls tables (Source Data), same definition and
# source files as make_composite2.py panel b: retained % = median adjusted increment / median
# original increment over the significant set. Verified 2026-09-17 to reproduce the previous
# hand-typed values exactly, cohort order CCRCC LUAD UCEC GBM PDAC.
_META = {"stroma": ("incr_W_over_stroma", "#9CC3E0", "o", 4.8, +0.15),
         "immune": ("incr_W_over_immune", "#6A51A3", "s", 4.2, 0.0),
         "proliferation": ("incr_W_over_proliferation", "#8C564B", "^", 5.0, -0.15),
         "all three": ("incr_W_over_ALL", INK, "D", 4.6, 0.0)}
ret = {}
for _name, (_col, _c1, _mk, _ms, _dy) in _META.items():
    _vals = []
    for _c in COH:
        _d = pd.read_csv(f"{DD}/composition_controls_{_c.lower()}.csv").dropna(subset=["incr_orig", _col])
        _vals.append(int(round(100 * _d[_col].median() / _d["incr_orig"].median())))
    ret[_name] = (_vals, _c1, _mk, _ms, _dy)
ypos = np.arange(5)[::-1]                          # CCRCC on top ... PDAC at bottom
axa.axvspan(95, 105, color="#EFEFEF", zorder=0, lw=0)
axa.axvline(100, color=GREY, lw=0.7, ls=(0, (3, 2)), zorder=1)
for y in ypos:
    axa.plot([0, 100], [y, y], color=LGREY, lw=0.8, zorder=1)
for name, (vals, col, mk, ms, dy) in ret.items():
    for i, v in enumerate(vals):
        axa.plot(v, ypos[i] + dy, mk, ms=ms, color=col, mec="white" if mk != "D" else INK,
                 mew=0.5, zorder=4 if mk != "D" else 5)
        if name == "all three":   # white box lifts the bold value off the grey guide line
            axa.text(v + 3.2, ypos[i] + 0.02, f"{v}", fontsize=7.0, fontweight="bold", color=INK,
                     ha="left", va="center", zorder=6,
                     bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none"))
axa.set_yticks(ypos); axa.set_yticklabels(COH, fontsize=FS_TICK)
axa.tick_params(axis="y", length=0); axa.tick_params(axis="x", labelsize=FS_TICK)
axa.set_xlim(-3, 110); axa.set_ylim(-2.0, 4.6)
axa.set_xticks([0, 25, 50, 75, 100])
axa.spines["left"].set_visible(False)
axa.set_xlabel("% of morphology increment retained\nover mRNA + composition score(s)", fontsize=7)
handles = [Line2D([], [], marker=mk, ls="none", ms=ms, color=col, mec="white" if mk != "D" else INK,
                  mew=0.5, label=name) for name, (_, col, mk, ms, _dy) in ret.items()]
axa.legend(handles=handles, fontsize=FS_LEG, frameon=False, loc="lower left", ncol=2,
           bbox_to_anchor=(-0.02, 0.035), handletextpad=0.3, columnspacing=0.9, labelspacing=0.3,
           borderaxespad=0.0)
axa.set_title("Not explained by\nmicroenvironment composition", loc="left", fontsize=FS_TITLE, fontweight="bold")

# =====================================================================================
# b, c : morphology-only out-of-fold score vs measured residual score (pooled over cohorts)
# =====================================================================================
frames = []
for c in COH:
    d = pd.read_csv(f"{DD}/scores_{c.lower()}.csv", index_col=0)
    d["cohort"] = c
    frames.append(d)
sc = pd.concat(frames)
stats_bc = {}
specs = [("b", "translation", "Morphology-only score recovers\nthe translation residual",
          "measured translation residual score\n(from proteomics)"),
         ("c", "secretion_ER", "Morphology-only score recovers\nthe ER-secretion residual",
          "measured ER-secretion residual score\n(from proteomics)")]
for k, (letter, key, title, xlab) in enumerate(specs):
    ax = fig.add_subplot(g1[k + 1]); lab(ax, letter, 0)
    sub = sc[[f"{key}_meas", f"{key}_morph", "cohort"]].dropna()
    x = sub[f"{key}_meas"].values; y = sub[f"{key}_morph"].values
    ax.axhline(0, color=LGREY, lw=0.6, zorder=0); ax.axvline(0, color=LGREY, lw=0.6, zorder=0)
    for c in COH:
        m = sub.cohort.values == c
        ax.scatter(x[m], y[m], s=9, color=ORG[c], alpha=0.8, edgecolors="none", zorder=3, label=c)
    r, p = st.pearsonr(x, y)
    slope, icpt = np.polyfit(x, y, 1)
    xx = np.array([x.min(), x.max()])
    ax.plot(xx, slope * xx + icpt, color=INK, lw=1.1, zorder=4)
    ax.text(0.03, 0.97, f"$r$ = {r:.2f}\n$P$ = {fmt_p(p)}\n$n$ = {len(sub)}", transform=ax.transAxes,
            fontsize=FS_STAT, ha="left", va="top", linespacing=1.15, zorder=6)
    ax.set_xlabel(xlab, fontsize=7)
    ax.set_ylabel("morphology-only score (from H&E)", fontsize=7)
    ax.tick_params(labelsize=FS_TICK)
    ax.set_title(title, loc="left", fontsize=FS_TITLE, fontweight="bold")
    stats_bc[key] = (len(sub), r, p, slope)
    if letter == "b":
        ax.legend(loc="lower right", fontsize=FS_LEG, frameon=False, handletextpad=0.2, labelspacing=0.25,
                  markerscale=1.1, borderaxespad=0.3)
        axb = ax
    else:
        axc = ax

# =====================================================================================
# d, e : nuclear correlates (Spearman rho, pathomics features x cohorts), shared colourbar
# =====================================================================================
feats = ["nuclear density", "nucleus count", "mean nucleus area", "nucleus area SD",
         "chromatin OD", "eosin fraction"]
_fkey = ["nuclear_density", "nucleus_count", "mean_nuc_area", "nuc_area_std",
         "chromatin_od", "eosin_fraction"]
# server_export/results/results__clinical_link_<cohort>_pathomics_corr.csv, same source as
# make_composite2.py panel c; verified 2026-09-17 to reproduce the previous hand-typed
# matrices to 2dp for both translation_meas and secretion_ER_meas.
_pc5 = {c: pd.read_csv(f"../server_export/results/results__clinical_link_{c.lower()}_pathomics_corr.csv")
        for c in COH}
rho_tr = np.array(
    [[float(_pc5[c][(_pc5[c].feature == k) & (_pc5[c].score == "translation_meas")]["rho"].iloc[0])
      for c in COH] for k in _fkey])
rho_se = np.array(
    [[float(_pc5[c][(_pc5[c].feature == k) & (_pc5[c].score == "secretion_ER_meas")]["rho"].iloc[0])
      for c in COH] for k in _fkey])
axd = fig.add_subplot(g2[0]); lab(axd, "d", 1, "L")
axe = fig.add_subplot(g2[1]); lab(axe, "e", 1)
cax = fig.add_subplot(g2[2])
for ax, M, title in [(axd, rho_tr, "Nuclear correlates of the translation residual"),
                     (axe, rho_se, "Nuclear correlates of the ER-secretion residual")]:
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            ax.text(j, i, f"{v:+.2f}".replace("-", "−"), ha="center", va="center", fontsize=7.0,
                    color="white" if abs(v) > 0.28 else "#333")
    ax.set_xticks(range(5)); ax.set_xticklabels(COH, rotation=45, ha="right", fontsize=FS_TICK)
    ax.set_yticks(range(6))
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, loc="left", fontsize=FS_TITLE, fontweight="bold")
axd.set_yticklabels(feats, fontsize=FS_TICK)
axe.set_yticklabels([])
cb = fig.colorbar(im, cax=cax)
cb.set_label(r"Spearman $\rho$", fontsize=FS_TICK, labelpad=2)
cb.ax.tick_params(labelsize=7.0, length=2)
cb.outline.set_visible(False)

# =====================================================================================
# f : Spearman rho with tumour grade, Fisher-z 95% CI (four gradeable cohorts)
# =====================================================================================
axf = fig.add_subplot(g3[0]); lab(axf, "f", 2, "L")
GCOH = ["CCRCC", "LUAD", "UCEC", "PDAC"]
# rho, n and BH q within each cohort's own clinical test family -- the primary
# correction (figdata/clinical_family_all.csv, 108-test enumeration). Identical
# source and rule as Fig 2d (make_composite2.py) and Table 3. Read live rather than
# hand-typed: corrects the same n transcription error found 2026-09-17 in Fig_clinical
# panel a (n was [103,103,97,137], file gives [103,102,100,135] for CCRCC/LUAD/UCEC/PDAC).
_cfa_f = pd.read_csv(f"{DD}/clinical_family_all.csv")
_cfa_f = _cfa_f[(_cfa_f.outcome == "grade") & (_cfa_f["mode"] == "H&E-only")]
def _famf(fam):
    r = _cfa_f[_cfa_f.family == fam].set_index("cohort").loc[GCOH]
    return r["n"].values, r["stat"].values, r["q_cohort"].values
n_tr, rho_tr_f, q_tr_f = _famf("Transl.")
n_se, rho_se_f, q_se_f = _famf("ER-secr.")
assert (n_tr == n_se).all(), "translation/ER-secretion n mismatch within a cohort"
n_grade = dict(zip(GCOH, n_tr))
grade = {
    "translation": {"rho": list(rho_tr_f), "q": list(q_tr_f)},
    "secretion":   {"rho": list(rho_se_f), "q": list(q_se_f)},
}
yf = np.arange(4)[::-1]                             # CCRCC top ... PDAC bottom
JIT = {"translation": +0.17, "secretion": -0.17}
axf.plot([0, 0], [-0.55, 3.55], color=INK, lw=0.6, ls=(0, (3, 2)), zorder=0)   # stops clear of legend/footnote
ci_log = []
for fam in ["translation", "secretion"]:
    col = FAM[fam]
    for i, c in enumerate(GCOH):
        rho, q, n = grade[fam]["rho"][i], grade[fam]["q"][i], n_grade[c]
        z = np.arctanh(rho); se = 1.0 / np.sqrt(n - 3)
        lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
        ci_log.append((fam, c, n, rho, lo, hi, q))
        y = yf[i] + JIT[fam]
        axf.plot([lo, hi], [y, y], color=col, lw=1.2, solid_capstyle="round", zorder=2)
        label = S.stars(q)
        filled = label != ""
        axf.plot(rho, y, "o", ms=5.0, mfc=col if filled else "white", mec=col, mew=1.0, zorder=4)
        if label:
            axf.text(hi + 0.02, y, label, fontsize=7.0, color=col,
                     ha="left", va="center", zorder=5)
axf.set_yticks(yf); axf.set_yticklabels(GCOH, fontsize=FS_TICK)
axf.tick_params(axis="y", length=0); axf.tick_params(axis="x", labelsize=FS_TICK)
axf.set_xlim(-0.35, 0.75); axf.set_ylim(-1.55, 4.75)
axf.set_xticks([-0.2, 0, 0.2, 0.4, 0.6])
axf.spines["left"].set_visible(False)
axf.set_xlabel(r"Spearman $\rho$ with tumour grade (95% CI, Fisher $z$)", fontsize=7)
hf = [Line2D([], [], marker="o", ls="-", lw=1.2, ms=5, color=FAM["translation"], mec=FAM["translation"],
             label="translation score"),
      Line2D([], [], marker="o", ls="-", lw=1.2, ms=5, color=FAM["secretion"], mec=FAM["secretion"],
             label="ER-secretion score"),
      Line2D([], [], marker="o", ls="none", ms=5, mfc=INK, mec=INK, label="filled = BH $q$ < 0.05, cohort family"),
      Line2D([], [], marker="o", ls="none", ms=5, mfc="white", mec=INK, mew=1.0, label="open = not significant")]
# two-column strip in the extra head-room above the CCRCC row (the zero line stops below it)
axf.legend(handles=hf, loc="upper left", bbox_to_anchor=(-0.01, 1.0), fontsize=FS_LEG, frameon=False, ncol=2,
           handletextpad=0.4, labelspacing=0.35, handlelength=1.6, columnspacing=1.2, borderaxespad=0.0)
# re-wrapped into four short lines at 7 pt so the block clears the UCEC translation CI, which reaches rho = 0.35
axf.text(0.74, 1.0, "GBM not shown: the\nGDC records neither\ngrade nor stage\n(uniformly grade IV)",
         fontsize=FS_LEG, color=GREY, ha="right", va="center", linespacing=1.1)
axf.text(-0.35, -1.45, "LUAD ER-secretion is the closest test to the threshold\n($p$ = 0.011, $q$ = 0.054); see Supplementary Table",
         fontsize=FS_LEG, color=GREY, ha="left", va="bottom", linespacing=1.1)
axf.set_title("Residual score tracks tumour grade", loc="left", fontsize=FS_TITLE, fontweight="bold")

# =====================================================================================
# g : CCRCC overall-survival forest (Cox, per SD) -- translation score keeps its panel-f colour
# =====================================================================================
axg = fig.add_subplot(g3[1]); lab(axg, "g", 2)
# review/recalc/recalc_survival_raw.json ["cox"] -- same source as Fig_clinical panel c,
# additionally pulling the "stage" term from the grade+stage-adjusted model. Verified to
# reproduce the panel's previous hand-typed values exactly (HR 1.708/1.335/3.022).
_cox_g = _json.load(open("../review/recalc/recalc_survival_raw.json", encoding="utf-8"))["cox"]
def _termg(block, name):
    t = next(t for t in _cox_g[block]["terms"] if t["term"] == name)
    return t["HR"], t["CI95"][0], t["CI95"][1], t["p"]
forest = [  # label, HR, lo, hi, P, face colour, edge/line colour
    ("translation score,\nunivariable",       *_termg("uni_score", "translation_morph"), FAM["translation"], FAM["translation"]),
    ("translation score,\nadj. grade + stage", *_termg("adj_grade_stage", "translation_morph"), GREY, GREY),
    ("stage (same model,\nper SD)",           *_termg("adj_grade_stage", "stage"), LGREY, INK),
]
yg = np.arange(3)[::-1]
axg.plot([1, 1], [-0.7, 2.7], color=INK, lw=0.6, ls=(0, (3, 2)), zorder=0)
for y, (labt, hr, lo, hi, p, fc, ec) in zip(yg, forest):
    axg.plot([lo, hi], [y, y], color=ec, lw=1.4, solid_capstyle="round", zorder=2)
    axg.plot(hr, y, "o", ms=6, mfc=fc, mec=ec, mew=0.8, zorder=4)
    axg.text(hi * 1.07, y, f"HR {hr:.2f}\n$P$ = {fmt_p(p)}", fontsize=FS_STAT, ha="left", va="center", color=INK,
             linespacing=1.05, zorder=5)
axg.set_xscale("log")
axg.set_xlim(0.6, 9.5)
axg.set_xticks([0.7, 1, 2, 3, 5, 7]); axg.set_xticklabels(["0.7", "1", "2", "3", "5", "7"], fontsize=FS_TICK)
axg.xaxis.set_minor_locator(mpl.ticker.NullLocator())
axg.set_yticks(yg); axg.set_yticklabels([f[0] for f in forest], fontsize=FS_TICK, linespacing=1.0)
axg.tick_params(axis="y", length=0)
axg.set_ylim(-0.7, 2.7)
axg.spines["left"].set_visible(False)
axg.set_xlabel("overall-survival hazard ratio (per SD)\nCCRCC, $n$ = 98, 21 events", fontsize=7)
axg.set_title("Stage absorbs the survival\neffect (CCRCC)", loc="left", fontsize=FS_TITLE, fontweight="bold")

# =====================================================================================
# panel letters from rendered bounding boxes + canvas fit
# =====================================================================================
LETTER_GAP_IN = 0.05
ASCENT = 0.716          # Arial cap/ascender height in em -> first-line baseline = bbox top - ASCENT*fontsize
DPI = fig.dpi
fig.canvas.draw()
_t = fig.text(0, 0, "g", fontsize=10, fontweight="bold")
LETTER_W = _t.get_window_extent(fig.canvas.get_renderer()).width / DPI; _t.remove()
LEFT_CONTENT = PAD_IN + LETTER_W + LETTER_GAP_IN     # inches: where each row's own content (labels) starts
ROWS = [(g1, [axa, axb, axc]), (g2, [axd, axe, cax]), (g3, [axf, axg])]


def place_letters():
    """Left-column letters (a/d/f) at one fixed x; the others 0.05 in left of their panel's label block.
    Each row's letters share the baseline of the first title line of the tallest title in that row."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ttop = {s: ax._left_title.get_window_extent(r).y1 / DPI for s, ax, _, _ in PANELS}
    rows = {}
    for s, ax, row, cg in PANELS:
        rows.setdefault(row, []).append(s)
    row_base = {row: max(ttop[s] for s in ss) - ASCENT * FS_TITLE / 72 for row, ss in rows.items()}
    out = []
    for s, ax, row, cg in PANELS:
        xl = (PAD_IN + LETTER_W) if cg == "L" else ax.get_tightbbox(r).x0 / DPI - LETTER_GAP_IN
        out.append(fig.text(xl / W_IN, row_base[row] / H_IN, s, fontsize=10, fontweight="bold",
                            ha="right", va="baseline"))
    return out


letters = []
for _ in range(4):                       # margins -> letters -> tight bboxes -> margins; converges in 2-3 passes
    for t in letters:
        t.remove()
    letters = place_letters()
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = fig.get_tightbbox(r)                                                # inches, letters included
    TOP += ((H_IN - PAD_IN) - bb.y1) / H_IN
    BOTTOM += (PAD_IN - bb.y0) / H_IN
    for (g, axs), (t_, b_) in zip(ROWS, row_bounds(TOP, BOTTOM)):
        x0 = min(ax.get_tightbbox(r).x0 for ax in axs) / DPI                 # row content, letters excluded
        x1 = max(ax.get_tightbbox(r).x1 for ax in axs) / DPI
        sp = g.get_subplot_params(fig)
        g.update(left=sp.left + (LEFT_CONTENT - x0) / W_IN, right=sp.right + ((W_IN - PAD_IN) - x1) / W_IN,
                 top=t_, bottom=b_)
for t in letters:
    t.remove()
letters = place_letters()
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
out_w, out_h = bb.width + 0.2, bb.height + 0.2                    # savefig pad_inches = 0.1 on each side

S.save_pub(fig, "Fig_supp2_scores")
nb, rb, pb, sb = stats_bc["translation"]; nc, rc, pc_, scc = stats_bc["secretion_ER"]
ci_txt = "; ".join(f"{f[:3]}/{c} n={n} rho={r:.2f} CI[{lo:.2f},{hi:.2f}] q={q:.2g}"
                   for f, c, n, r, lo, hi, q in ci_log)
print(f"wrote Fig_supp2_scores (canvas {W_IN} x {H_IN} in, tight export {out_w:.2f} x {out_h:.2f} in, 7 panels a-g) "
      f"| a: composition retained-% dot plot (5 cohorts x 4 series) "
      f"| b: translation pooled n={nb} r={rb:.3f} P={pb:.2e} slope={sb:.2f} | c: ER-secretion pooled n={nc} r={rc:.3f} "
      f"P={pc_:.2e} slope={scc:.2f} | d,e: 6x5 Spearman heatmaps | f: Fisher-z CIs {ci_txt} | g: CCRCC forest 3 rows "
      f"HR {forest[0][1]:.2f}/{forest[1][1]:.2f}/{forest[2][1]:.2f}")
