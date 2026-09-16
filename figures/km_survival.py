#!/usr/bin/env python3
"""KM survival panel (CCRCC), morphology-only translation-residual score split at median.
DESCRIPTIVE ONLY — paired with the grade/stage-adjusted forest plot (Fig8b), which shows
the association is NOT independent of grade. Kaplan-Meier + log-rank in numpy (no lifelines).
Needs (WinSCP from server /public/home/fjhui/ZW/results/ into figdata/):
  clinical_link_ccrcc_scores.csv   clinical_link_ccrcc_clinical.csv"""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2
import mrstyle as S

DD = "figdata"
SCORE = "translation_morph"
sp = f"{DD}/scores_ccrcc.csv"                              # already local (morph scores)
cp = f"{DD}/clinical_link_ccrcc_clinical.csv"              # WinSCP this one (time/event)
if not os.path.exists(cp):
    raise SystemExit(f"missing {cp} -- WinSCP it from server /public/home/fjhui/ZW/results/ into figdata/")

sc = pd.read_csv(sp, index_col=0); cl = pd.read_csv(cp, index_col=0)
df = sc.join(cl[["time", "event"]]).dropna(subset=[SCORE, "time", "event"])
df = df[df["time"] >= 0]
t = df["time"].values.astype(float); e = df["event"].values.astype(int)
g = (df[SCORE].values >= np.median(df[SCORE].values)).astype(int)   # 1 = high score


def km(t, e):
    xs, ys, Sv = [0.0], [1.0], 1.0
    for tt in np.unique(t):
        n = (t >= tt).sum(); d = int(((t == tt) & (e == 1)).sum())
        if n > 0:
            Sv *= (1 - d / n)
        xs.append(tt); ys.append(Sv)
    return np.array(xs), np.array(ys)


def logrank(t, e, g):
    O1 = E1 = V = 0.0
    for tt in np.unique(t[e == 1]):
        risk = t >= tt; n = risk.sum(); n1 = (risk & (g == 1)).sum()
        d = (t == tt) & (e == 1); dd = int(d.sum()); d1 = int((d & (g == 1)).sum())
        O1 += d1
        if n > 1:
            E1 += dd * n1 / n
            V += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    x2 = (O1 - E1) ** 2 / V if V > 0 else 0.0
    return float(1 - chi2.cdf(x2, 1))


p = logrank(t, e, g)
fig, ax = plt.subplots(figsize=(3.2, 2.7))
for gi, col, lab in [(1, S.ORG["CCRCC"], "high residual score"), (0, S.GREY, "low residual score")]:
    m = g == gi
    xs, ys = km(t[m], e[m])
    ax.step(xs, ys, where="post", color=col, lw=1.5, label=f"{lab} (n={m.sum()})")
ax.set_xlabel("overall survival (days)"); ax.set_ylabel("survival probability")
ax.set_ylim(0, 1.02)
ax.text(0.04, 0.06, f"log-rank P = {p:.3f}", transform=ax.transAxes, fontsize=7)
ax.text(0.04, 0.14, "descriptive; not independent of grade (see forest)", transform=ax.transAxes,
        fontsize=5.6, color=S.GREY, style="italic")
ax.legend(fontsize=6, loc="upper right")
S.save_pub(fig, "Fig_km_ccrcc")
print(f"wrote Fig_km_ccrcc  (n={len(df)}, deaths={e.sum()}, log-rank P={p:.4f})")
