#!/usr/bin/env python3
"""Fig_clinical -- the clinical block, split out of the old ten-panel Fig_master2.

WHY THIS IS ITS OWN FIGURE
Fig_master2 carried ten panels on a 7.2 x 9.6 in canvas. To fit a page with its
caption it had to be scaled to 60.7%, which put its smallest labels at 2.2 pt --
far below the 5 pt a journal accepts at technical check, and unfixable by raising
fonts alone because there was no room left on the canvas. The panels are therefore
split by message: Fig_master2 keeps the robustness evidence, this figure takes the
clinical block, and the two WSI spatial maps move to the Supplement (they were
always a per-patient illustration rather than a result; see the T1-8 note in
review/NC_ROADMAP.md). Each resulting figure prints at close to its authored size.

Panels (formerly Fig_master2 d, e, f):
  a  morphology-only residual score versus tumour grade, four gradeable cohorts
  b  CCRCC overall survival by a median split of the H&E-only translation score
  c  the same association univariable versus grade+stage-adjusted

Every number and encoding is carried over unchanged from make_composite2.py.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats as st
from scipy.stats import chi2
import mrstyle as S

ORG, INK, GREY, LGREY, FAM = S.ORG, S.INK, S.GREY, S.LGREY, S.FAM
DD = "figdata"


def lab(ax, s, x=-0.18, y=1.05):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


# B7: height 2.7 -> 3.2 in so panel b can carry a numbers-at-risk table under its
# x-axis; panels a and c keep their widths and gain the extra height.
FH = 3.2
fig = plt.figure(figsize=(7.2, FH))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.42,
                      left=0.085, right=0.985, top=(FH - 0.38) / FH, bottom=0.60 / FH)

# ===== a : residual score tracks grade =====
# title carries pad=18 to clear the in-axes legend, so the panel letter must rise
# above it or it reads as belonging to the panel below
axd = fig.add_subplot(gs[0, 0]); lab(axd, "a", x=-0.14, y=1.155)
gc = ["CCRCC", "LUAD", "UCEC", "PDAC"]
# filled = BH q<0.05 within the cohort's own clinical test family -- the primary
# correction declared in Methods and used identically in Table 3 and Supp Fig 2f.
# rho, p and the BH q_cohort<0.05 flags are read from figdata/clinical_family_all.csv
# (108-test enumeration), mode=H&E-only, outcome=grade. The CSV's "n" column is NOT the
# grade n: it is the survival-evaluable n (LUAD 102, UCEC 100, PDAC 135). The grade n
# is the number of analysed patients with a G1-G4 grade (GX / Unknown / missing are not
# gradeable), counted below from the same clinical files and checked against the CSV's
# rho. (An earlier revision of this script replaced the correct hand-typed
# [103,103,97,137] with that survival n; the hand-typed values were right.)
_cfa = pd.read_csv(f"{DD}/clinical_family_all.csv")
_cfa = _cfa[(_cfa.outcome == "grade") & (_cfa["mode"] == "H&E-only")]
def _fam(fam):
    r = _cfa[_cfa.family == fam].set_index("cohort").loc[gc]
    return r["stat"].values, r["p"].values, (r["q_cohort"].values < 0.05)
trho, tp, tsig = _fam("Transl.")
srho, sp, ssig = _fam("ER-secr.")
from scipy.stats import spearmanr as _spr
ng = []
for _i, _c in enumerate(gc):
    _sc = pd.read_csv(f"{DD}/scores_{_c.lower()}.csv", index_col=0)
    _cl = pd.read_csv(f"{DD}/gdc_clinical_{_c.lower()}.csv").set_index("case").reindex(_sc.index)
    _g = _cl["tumor_grade"].astype(str).map({"G1": 1, "G2": 2, "G3": 3, "G4": 4})
    _ok = _g.notna()
    ng.append(int(_ok.sum()))
    assert abs(_spr(_sc.loc[_ok, "translation_morph"], _g[_ok])[0] - trho[_i]) < 1e-6, \
        f"{_c}: rho on the counted graded set does not reproduce the enumeration's rho"
ng = np.array(ng)
assert list(ng) == [103, 103, 97, 137], f"graded n moved: {list(ng)}"
# dot-and-interval: 95% CI from the Fisher z-transform, se = 1/sqrt(n-3)
yrow = np.arange(4)[::-1]; JIT = 0.17
for off, rho, sig, col, l in [(+JIT, trho, tsig, FAM["translation"], "translation"),
                              (-JIT, srho, ssig, FAM["secretion"], "ER-secretion")]:
    z = np.arctanh(rho); se = 1.0 / np.sqrt(ng - 3)
    lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    for y_, r_, lo_, hi_, s_ in zip(yrow + off, rho, lo, hi, sig):
        axd.plot([lo_, hi_], [y_, y_], color=col, lw=1.0, solid_capstyle="butt", zorder=2)
        axd.plot(r_, y_, "o", ms=4.6, mfc=(col if s_ else "white"), mec=col, mew=1.0, zorder=3)
axd.axvline(0, color="k", lw=0.6, zorder=1)
axd.set_yticks(yrow); axd.set_yticklabels(gc, fontsize=7)
axd.set_ylim(-0.6, 4.9)
axd.set_xlim(-0.35, 0.72)
axd.set_xlabel(r"Spearman $\rho$ with grade (95% CI, Fisher $z$)", fontsize=7)
axd.tick_params(axis="y", length=0)
hd = [Line2D([], [], ls="none", marker="o", ms=4.2, mfc=FAM["translation"], mec=FAM["translation"], label="translation"),
      Line2D([], [], ls="none", marker="o", ms=4.2, mfc=FAM["secretion"], mec=FAM["secretion"], label="ER-secretion"),
      Line2D([], [], ls="none", marker="o", ms=4.2, mfc=INK, mec=INK, label="BH $q$<0.05, cohort family"),
      Line2D([], [], ls="none", marker="o", ms=4.2, mfc="white", mec=INK, label="not significant")]
axd.legend(handles=hd, loc="upper left", fontsize=7, handlelength=1.0, handletextpad=0.4,
           labelspacing=0.25, ncol=2, columnspacing=0.8, bbox_to_anchor=(-0.02, 1.03), frameon=False)
axd.set_title("Residual score tracks tumour grade", fontsize=7.6, fontweight="bold", loc="left", pad=18)

# ===== b : KM survival (descriptive) =====
axe = fig.add_subplot(gs[0, 1]); lab(axe, "b", x=-0.20)
sc = pd.read_csv(f"{DD}/scores_ccrcc.csv", index_col=0)
cl = pd.read_csv(f"{DD}/clinical_link_ccrcc_clinical.csv", index_col=0)
df = sc.join(cl[["time", "event"]]).dropna(subset=["translation_morph", "time", "event"])
df = df[df["time"] >= 0]
t = df["time"].values.astype(float); ev = df["event"].values.astype(int)
g = (df["translation_morph"].values >= np.median(df["translation_morph"].values)).astype(int)


def km(t, e):
    xs, ys, Sv = [0.0], [1.0], 1.0
    for tt in np.unique(t):
        n = (t >= tt).sum(); d = int(((t == tt) & (e == 1)).sum())
        Sv *= (1 - d / n) if n > 0 else 1
        xs.append(tt); ys.append(Sv)
    return np.array(xs), np.array(ys)


def logrank(t, e, g):
    O1 = E1 = V = 0.0
    for tt in np.unique(t[e == 1]):
        risk = t >= tt; n = risk.sum(); n1 = (risk & (g == 1)).sum()
        d = (t == tt) & (e == 1); dd = int(d.sum()); d1 = int((d & (g == 1)).sum())
        O1 += d1
        if n > 1:
            E1 += dd * n1 / n; V += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    return float(1 - chi2.cdf((O1 - E1) ** 2 / V, 1)) if V > 0 else 1.0


pv = logrank(t, ev, g)
GROUPS = [(1, ORG["CCRCC"], "high score"), (0, GREY, "low score")]
for gi, col, l in GROUPS:
    mm = g == gi; xs, ys = km(t[mm], ev[mm])
    axe.step(xs, ys, where="post", color=col, lw=1.4, label=f"{l} (n={mm.sum()})")
    # censoring ticks: one vertical mark per censored patient, on the curve at the
    # survival probability in force at that patient's censoring time
    tc = np.sort(t[mm & (ev == 0)])
    yc = ys[np.searchsorted(xs, tc, side="right") - 1]
    axe.plot(tc, yc, ls="none", marker="|", ms=4.5, mew=0.8, color=col, zorder=3)
axe.set_xlabel("overall survival (days)", fontsize=7); axe.set_ylabel("survival probability", fontsize=7)
axe.set_ylim(0, 1.02); axe.text(0.04, 0.10, f"log-rank P = {pv:.3f}", transform=axe.transAxes, fontsize=7)
# frameon=False put the legend samples directly on the KM curves, so the blue
# high-score swatch was hidden under the grey low-score curve
# B7: with censoring ticks the upper-right corner is no longer clear (the low-score
# curve runs along 0.8-0.98 there), so the legend sits in the empty band below the
# curves, above the log-rank text
axe.legend(fontsize=7, loc="lower left", bbox_to_anchor=(0.02, 0.19), frameon=True,
           framealpha=0.92, edgecolor="none", borderpad=0.3)
axe.set_title("Survival (CCRCC, descriptive)", fontsize=7.6, fontweight="bold", loc="left")

# numbers at risk (n with time >= t) at evenly spaced times, computed from the same
# t / ev / g arrays the curves use; observed follow-up runs to ~2006 days
TT = np.arange(0, 2001, 500)
axe.set_xticks(TT); axe.set_xlim(-60, float(t.max()) * 1.02)
_p = axe.get_position()
axe.set_position([_p.x0, 0.98 / FH, _p.width, _p.y1 - 0.98 / FH])
axr = fig.add_axes([_p.x0, 0.06 / FH, _p.width, 0.46 / FH], sharex=axe)
axr.set_ylim(-0.5, 2.5)
for sp_ in axr.spines.values():
    sp_.set_visible(False)
axr.set_yticks([]); axr.tick_params(axis="x", length=0, labelbottom=False)
axr.patch.set_visible(False)
axr.text(-0.02, 2.0, "No. at risk", transform=axr.get_yaxis_transform(), fontsize=7,
         fontweight="bold", ha="left", va="center")
_atrisk = {}
for row, (gi, col, l) in zip((1, 0), GROUPS):
    axr.text(-0.04, row - 0.15, l.split()[0], transform=axr.get_yaxis_transform(),
             fontsize=7, color=col, ha="right", va="center")
    _atrisk[gi] = [int(((g == gi) & (t >= tt)).sum()) for tt in TT]
    for tt, nn in zip(TT, _atrisk[gi]):
        axr.text(tt, row - 0.15, str(nn), fontsize=7, color=col, ha="center", va="center")
assert _atrisk[1][0] + _atrisk[0][0] == len(df)

# ===== c : honest forest (uni vs grade/stage-adjusted) =====
axf = fig.add_subplot(gs[0, 2]); lab(axf, "c", x=-0.28)

# review/recalc/recalc_survival_raw.json ["cox"]["uni_score"/"adj_grade_stage"]["terms"][0]
# (term="translation_morph"): true model CI, not the HR+p log-normal approximation this
# panel previously used. Corrects a transcription error found 2026-09-17: adjusted p was
# hand-typed as 0.270, the model gives 0.26855 -> 0.269 (matches panel g, which already
# had this right).
import json as _json
_cox = _json.load(open("../review/recalc/recalc_survival_raw.json", encoding="utf-8"))["cox"]
def _term(block, name="translation_morph"):
    t = next(t for t in _cox[block]["terms"] if t["term"] == name)
    return t["HR"], t["p"], t["CI95"]
_hr_u, _p_u, _ci_u = _term("uni_score")
_hr_a, _p_a, _ci_a = _term("adj_grade_stage")
rows = [("univariable", _hr_u, _p_u, _ci_u, ORG["CCRCC"]),
        ("adj. grade+stage", _hr_a, _p_a, _ci_a, GREY)]
for i, (name, hr, p, civals, col) in enumerate(rows):
    lo, hi = civals; yv = 1 - i
    axf.plot([lo, hi], [yv, yv], color=col, lw=1.6); axf.plot(hr, yv, "o", color=col, ms=7)
    axf.text(hi + 0.10, yv, f"HR {hr:.2f}\np={p:.3f}", fontsize=7, va="center")
axf.axvline(1, ls="--", lw=0.8, color="k")
axf.set_yticks([0, 1]); axf.set_yticklabels(["adj.\ngrade+stage", "univariable"], fontsize=7)
axf.set_xscale("log"); axf.minorticks_off()
axf.set_xticks([0.7, 1, 2, 3]); axf.set_xticklabels(["0.7", "1", "2", "3"], fontsize=7)
axf.set_xlim(0.62, 5.2); axf.set_ylim(-0.6, 1.6)
axf.set_xlabel("OS hazard ratio (per SD)", fontsize=7)
axf.set_title("Not independent of grade", fontsize=7.6, fontweight="bold", loc="left")

S.save_pub(fig, "Fig_clinical")
print(f"wrote Fig_clinical | KM log-rank P={pv:.3f}, n={len(df)}, at-risk high={_atrisk[1]} low={_atrisk[0]}, censor ticks={int((ev==0).sum())}")
