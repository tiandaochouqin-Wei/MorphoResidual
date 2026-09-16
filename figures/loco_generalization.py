#!/usr/bin/env python3
"""loco_generalization.py (LOCAL) — leave-one-cancer-out generalisation of the
morphology->residual map. For each held-out organ: train PCA(20)+ridge on the other four
organs pooled (UNI mean-embedding -> measured translation residual score), predict the held-out
organ zero-shot; compare with the within-organ 5-fold CV reference.
Needs figdata/loco_<cohort>.npz (from export_cohort_features.py)."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
import mrstyle as S

COH, DD, NPC, ALPHA = S.COH, "figdata", 20, 1.0
TARGETS = ["translation_meas", "secretion_ER_meas"]
data = {}
for c in COH:
    p = f"{DD}/loco_{c.lower()}.npz"
    if not os.path.exists(p):
        raise SystemExit(f"missing {p} -- run export_cohort_features.py on the server and WinSCP loco_*.npz back")
    data[c] = np.load(p, allow_pickle=True)


def fit_predict(Xtr, ytr, Xte):
    pca = PCA(NPC, svd_solver="full").fit(Xtr)
    Ptr, Pte = pca.transform(Xtr), pca.transform(Xte)
    xm, xs = Ptr.mean(0), Ptr.std(0); xs[xs == 0] = 1.0; ym = ytr.mean()
    A = (Ptr - xm) / xs
    w = np.linalg.solve(A.T @ A + ALPHA * np.eye(NPC), A.T @ (ytr - ym))
    return ((Pte - xm) / xs) @ w + ym


def r2(y, p):
    return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)


rows = []
for tgt in TARGETS:
    if not all(tgt in data[c] for c in COH):
        continue
    for H in COH:
        Xh, yh = data[H]["emb"].astype(float), data[H][tgt].astype(float)
        m = np.isfinite(yh); Xh, yh = Xh[m], yh[m]
        Xtr = np.vstack([data[c]["emb"].astype(float)[np.isfinite(data[c][tgt])] for c in COH if c != H])
        ytr = np.concatenate([data[c][tgt].astype(float)[np.isfinite(data[c][tgt])] for c in COH if c != H])
        pz = fit_predict(Xtr, ytr, Xh)
        po = np.full(len(yh), np.nan)
        for tr, te in KFold(5, shuffle=True, random_state=0).split(Xh):
            po[te] = fit_predict(Xh[tr], yh[tr], Xh[te])
        rows.append(dict(target=tgt.replace("_meas", ""), held_out=H, n=len(yh),
                         zero_shot_r=np.corrcoef(yh, pz)[0, 1], zero_shot_R2=r2(yh, pz),
                         within_r=np.corrcoef(yh, po)[0, 1], within_R2=r2(yh, po)))
res = pd.DataFrame(rows)
res.to_csv("loco_results.csv", index=False)
print(res.round(3).to_string(index=False))

# ---- dumbbell (paired-dot) chart: within-organ CV r (filled) -> zero-shot r (open), per organ ----
from matplotlib.lines import Line2D

tgts = [t for t in ["translation", "secretion_ER"] if t in set(res.target)]
PANEL = {"translation": "Translation residual programme", "secretion_ER": "ER-secretion residual programme"}
fig, axes = plt.subplots(1, len(tgts), figsize=(6.4, 2.6), sharey=True)
axes = np.atleast_1d(axes)
fig.subplots_adjust(left=0.115, right=0.985, bottom=0.25, top=0.79, wspace=0.10)
ypos = np.arange(len(COH))[::-1]                      # Kidney on top ... Pancreas at bottom
XLIM, XTICKS, OFF = (-0.22, 0.85), [-0.2, 0, 0.2, 0.4, 0.6, 0.8], 0.046
for ax, tgt, letter in zip(axes, tgts, "ab"):
    tr = res[res.target == tgt].set_index("held_out").loc[COH]
    ax.axvline(0, color=S.INK, lw=0.6, ls=(0, (3, 2)), zorder=0)
    for y, c in zip(ypos, COH):
        wi, zs, col = float(tr.within_r[c]), float(tr.zero_shot_r[c]), S.ORG[c]
        ax.plot([wi, zs], [y, y], color=S.GREY, lw=1.0, solid_capstyle="round", zorder=1)   # connector = drop
        ax.scatter([zs], [y], s=30, facecolor="white", edgecolor=col, linewidths=1.1, zorder=2)  # open = zero-shot
        ax.scatter([wi], [y], s=30, facecolor=col, edgecolor=col, linewidths=0.8, zorder=3)    # filled = within-organ
        # if the two dots coincide (< 1 marker width apart) add a white halo so the ring stays visible on the disc
        if abs(zs - wi) * ax.get_position().width * fig.get_figwidth() * 72 / (XLIM[1] - XLIM[0]) < 6:
            ax.scatter([zs], [y], s=30, facecolor="none", edgecolor="white", linewidths=3.2, zorder=4)
        ax.scatter([zs], [y], s=30, facecolor="none", edgecolor=col, linewidths=1.1, zorder=5)   # ring drawn on top
        side = -1 if zs <= wi else 1                     # label on the side away from the filled dot
        ax.text(zs + side * OFF, y, f"{zs:.2f}".replace("-", "−"), ha="right" if side < 0 else "left", va="center",
                fontsize=6.2, color=S.INK, zorder=4,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
    ax.set_xlim(*XLIM); ax.set_xticks(XTICKS)
    ax.set_ylim(-0.65, len(COH) - 0.35)
    ax.set_yticks(ypos); ax.set_yticklabels([S.ORGAN[c] for c in COH])
    ax.tick_params(axis="y", length=0)
    ax.set_title(PANEL[tgt], fontsize=7.5, fontweight="bold", loc="left", pad=4)
for ax, letter, dx in zip(axes, "ab", (0.0, -0.035)):
    bb = ax.get_position()
    fig.text(0.005 if letter == "a" else bb.x0 + dx, bb.y1 + 0.03, letter,
             fontsize=9, fontweight="bold", ha="left", va="bottom")
fig.text(0.005, 0.965, "Morphology$\\rightarrow$residual map: shared direction, organ-specific parts",
         fontsize=8, fontweight="bold", ha="left", va="top")
fig.text((axes[0].get_position().x0 + axes[-1].get_position().x1) / 2, 0.115,
         "Pearson $r$, predicted vs measured residual score", fontsize=7, ha="center", va="center")
handles = [Line2D([], [], marker="o", ls="none", ms=5, mfc=S.INK, mec=S.INK,
                  label="filled = within-organ 5-fold CV"),
           Line2D([], [], marker="o", ls="none", ms=5, mfc="white", mec=S.INK, mew=1.0,
                  label="open = zero-shot (trained on the other 4 organs)")]
fig.legend(handles=handles, fontsize=6, frameon=False, ncol=2, loc="lower center",
           bbox_to_anchor=(0.5, 0.0), handletextpad=0.4, columnspacing=1.8)
S.save_pub(fig, "Fig_loco")
print("wrote Fig_loco + loco_results.csv")
