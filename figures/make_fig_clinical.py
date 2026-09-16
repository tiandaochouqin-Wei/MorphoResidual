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


fig = plt.figure(figsize=(7.2, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.42,
                      left=0.085, right=0.985, top=0.86, bottom=0.20)

# ===== a : residual score tracks grade =====
# title carries pad=18 to clear the in-axes legend, so the panel letter must rise
# above it or it reads as belonging to the panel below
axd = fig.add_subplot(gs[0, 0]); lab(axd, "a", x=-0.14, y=1.19)
gc = ["CCRCC", "LUAD", "UCEC", "PDAC"]
# filled = BH q<0.05 within the cohort's own clinical test family -- the primary
# correction declared in Methods and used identically in Table 3 and Supp Fig 2f.
# Read from figdata/clinical_family_all.csv (108-test enumeration), mode=H&E-only,
# outcome=grade. Corrects a transcription error found 2026-09-17: n was previously
# hand-typed as [103,103,97,137] for CCRCC/LUAD/UCEC/PDAC; the file gives
# [103,102,100,135] (LUAD/UCEC/PDAC off by 1/3/2). rho, p and the q_cohort<0.05 flags
# below were already correct.
_cfa = pd.read_csv(f"{DD}/clinical_family_all.csv")
_cfa = _cfa[(_cfa.outcome == "grade") & (_cfa["mode"] == "H&E-only")]
def _fam(fam):
    r = _cfa[_cfa.family == fam].set_index("cohort").loc[gc]
    return r["n"].values, r["stat"].values, r["p"].values, (r["q_cohort"].values < 0.05)
ng_t, trho, tp, tsig = _fam("Transl.")
ng_s, srho, sp, ssig = _fam("ER-secr.")
assert (ng_t == ng_s).all(), "translation/ER-secretion n mismatch within a cohort"
ng = ng_t
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
for gi, col, l in [(1, ORG["CCRCC"], "high score"), (0, GREY, "low score")]:
    mm = g == gi; xs, ys = km(t[mm], ev[mm])
    axe.step(xs, ys, where="post", color=col, lw=1.4, label=f"{l} (n={mm.sum()})")
axe.set_xlabel("overall survival (days)", fontsize=7); axe.set_ylabel("survival probability", fontsize=7)
axe.set_ylim(0, 1.02); axe.text(0.04, 0.10, f"log-rank P = {pv:.3f}", transform=axe.transAxes, fontsize=7)
# frameon=False put the legend samples directly on the KM curves, so the blue
# high-score swatch was hidden under the grey low-score curve
axe.legend(fontsize=7, loc="upper right", frameon=True, framealpha=0.92,
           edgecolor="none", borderpad=0.3)
axe.set_title("Survival (CCRCC, descriptive)", fontsize=7.6, fontweight="bold", loc="left")

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
print(f"wrote Fig_clinical | KM log-rank P={pv:.3f}, n={len(df)}")
