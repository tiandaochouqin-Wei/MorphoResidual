#!/usr/bin/env python3
"""New MorphoResidual figures (server-verified numbers): composition controls,
second foundation-model replication, clinical association. Submission-grade (mrstyle)."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats as st
import mrstyle as S

COH = S.COH
x = np.arange(5)

# ============================================================ FIG 6 — composition controls
# WSI increment retained (% of the original mRNA->morphology increment) after adjusting for
# each transcriptomic composition/state score, and for all three jointly.
retain = {  # cohort order: CCRCC LUAD UCEC GBM PDAC
    "stroma":        np.array([87, 93, 70, 86, 85]),
    "immune":        np.array([99, 78, 92, 77, 90]),
    "proliferation": np.array([86, 32, 83, 44, 48]),
    "all three":     np.array([68,  6, 36,  3, 34]),
}
ctrl_col = {"stroma": "#9CC3E0", "immune": "#5AA0B8",
            "proliferation": "#E69F00", "all three": S.INK}

fig, ax = plt.subplots(figsize=(5.4, 2.7))
ax.axhspan(95, 105, color="#EFEFEF", zorder=0)
ax.axhline(100, ls="--", lw=0.6, color=S.GREY, zorder=1)
w = 0.2
for k, (name, vals) in enumerate(retain.items()):
    ax.bar(x + (k - 1.5) * w, vals, w, color=ctrl_col[name], edgecolor="none",
           label=name, zorder=2)
# annotate the strict "all three" bar per cohort
for i, v in enumerate(retain["all three"]):
    ax.text(i + 1.5 * w, v + 2, f"{v}", ha="center", va="bottom",
            fontsize=5.6, color=S.INK, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels([f"{c}" for c in COH])
ax.set_ylabel("Morphology increment retained\nover mRNA + composition (%)")
ax.set_ylim(0, 112)
ax.set_title("Not microenvironment composition; partly a proliferation programme",
             loc="left", fontsize=8, fontweight="bold")
ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.34),
          handlelength=1.1, columnspacing=1.2, fontsize=6)
S.save_pub(fig, "Fig6_composition")
plt.close(fig)

# NOTE: Fig 7 (second foundation model) is built with the REAL per-gene hexbin
# scatter in make_realdata_figs.py -- do not rebuild the bar version here.

# ============================================================ FIG 8 — clinical association
gc = ["CCRCC", "LUAD", "UCEC", "PDAC"]  # gradeable (GBM is uniformly grade IV)
xg = np.arange(4)
transl_rho = np.array([0.37, 0.40, 0.16, 0.21]); transl_p = np.array([1.3e-4, 3.0e-5, 0.124, 0.013])
secr_rho   = np.array([0.44, 0.25, -0.07, 0.19]); secr_p   = np.array([3.5e-6, 0.011, 0.52, 0.028])
# BH q over the 17 reported clinical tests (adversarial review, 2026-09-06; review/recalc/clinical.json);
# PDAC's q<0.05 here does not survive correction within PDAC's own larger test family (q=0.052/0.067,
# SuppTable_clinical_family) -- it is annotated "nominal" rather than starred.
transl_q = np.array([0.00073, 0.00026, 0.2218, 0.0368]); secr_q = np.array([6.0e-5, 0.0368, 0.6314, 0.0667])
NOMINAL = {("transl", "PDAC")}  # (score, cohort) pairs to annotate as nominal, not starred

fig = plt.figure(figsize=(6.8, 2.7))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.42)
a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1])

# a: grade Spearman rho, two morphology-only scores
wc = 0.36
for key, xi, rho, q, col, lab in [("transl", xg - wc / 2, transl_rho, transl_q, S.FAM["translation"], "translation score"),
                                   ("secr", xg + wc / 2, secr_rho, secr_q, S.FAM["secretion"], "ER-secretion score")]:
    a.bar(xi, rho, wc, color=col, label=lab)
    for c_, xx, r, qq in zip(gc, xi, rho, q):
        s = "n.s." if (key, c_) in NOMINAL else S.stars(qq)
        if s:
            a.text(xx, r + (0.01 if r >= 0 else -0.03), s, ha="center",
                   va="bottom" if r >= 0 else "top", fontsize=6.4 if s == "n.s." else 7, color=S.INK)
a.axhline(0, color="k", lw=0.6)
a.set_xticks(xg); a.set_xticklabels(gc)
a.set_ylabel(r"Spearman $\rho$ with tumour grade" + "\n(morphology-only score)")
a.set_ylim(-0.15, 0.55)
a.legend(loc="upper right", fontsize=5.8, handlelength=1.0)
a.set_title("a  H&E residual score tracks grade in CCRCC, LUAD",
            loc="left", fontsize=8, fontweight="bold")
a.text(0.02, 0.03, "stars = BH $q$<0.05 (17 reported tests); PDAC translation ($q$=0.037) is nominal only"
       "\nwithin PDAC's own larger test family (Supplementary Table)", transform=a.transAxes, fontsize=4.6, color=S.GREY)

# b: CCRCC survival forest -- univariable vs grade+stage-adjusted (Wald 95% CI, lifelines; n=98, 21 events)
# exact CIs from clinical_family_ccrcc.py / review/recalc/clinical.json ccrcc_survival, not back-computed from p
rows = [("univariable", 1.708, 1.157, 2.521, 0.0070, S.ORG["CCRCC"]),
        ("adj. grade+stage", 1.335, 0.800, 2.226, 0.269, S.GREY),
        ("adj. stage (per SD)", 3.022, 1.593, 5.733, 7.1e-4, S.LGREY)]
for i, (lab, hr, lo, hi, p, col) in enumerate(rows):
    yy = len(rows) - 1 - i
    b.plot([lo, hi], [yy, yy], color=col, lw=1.4, solid_capstyle="round")
    b.plot(hr, yy, "o", color=col, ms=6)
    b.text(hi + 0.05, yy, f"HR {hr:.2f}\n$p$={p:.3g}", va="center", fontsize=5.8, color=S.INK)
b.axvline(1, ls="--", lw=0.7, color="k")
b.set_yticks(range(len(rows))); b.set_yticklabels(["stage\n(same model)", "translation score\nadj. grade+stage",
                                                    "translation score\nunivariable"], fontsize=6.0)
b.set_ylim(-0.6, len(rows) - 0.4)
b.set_xscale("log")
b.set_xlim(0.6, 7.5)
b.xaxis.set_major_locator(plt.FixedLocator([0.7, 1, 2, 3, 5, 7]))
b.xaxis.set_minor_locator(plt.NullLocator())
b.xaxis.set_major_formatter(plt.FixedFormatter(["0.7", "1", "2", "3", "5", "7"]))
b.set_xlabel("Overall-survival hazard ratio (per SD)\nCCRCC, $n$=98, 21 events")
b.set_title("b  Stage absorbs the effect", loc="left", fontsize=8, fontweight="bold")
b.spines["left"].set_visible(False); b.tick_params(axis="y", length=0)
fig.tight_layout()
S.save_pub(fig, "Fig8_clinical")
plt.close(fig)

print("wrote Fig6_composition, Fig7_second_fm, Fig8_clinical {svg,pdf,png,tiff}")
