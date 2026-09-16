#!/usr/bin/env python3
"""Fig_wsi_spatial -- the two whole-slide residual maps, moved to the Supplement.

These were panels g and h of the old ten-panel Fig_master2. Two reasons for the
move, and the second is the substantive one:

  * Fig_master2 had to be scaled to 60.7% to fit a page, putting its smallest text
    at 2.2 pt. Splitting it is the only way to print any of it at a legible size.
  * These maps were a per-patient illustration, not a result, because the direction
    they project was fitted on the full sample and the two slides were the extremes
    of that same fit. T1-8 has now redone both: PCA and the ridge direction are
    refitted inside each of the 5 folds, so a patient's tiles are scored by a model
    that never saw them, and the extremes are chosen on those out-of-fold scores.
    The high-residual slide changes (C3L-00011-21 -> C3N-01651-21); the low one does
    not.

    The old caption also claimed the within-slide share of score variance was
    roughly 10%, i.e. that the signal was overwhelmingly between patients rather
    than spatial. That number was computed on the two displayed slides ONLY, which
    were selected as the score extremes and therefore maximise between-slide
    variance by construction; it reproduces at 9.8% on that pair and is an artefact
    of the selection. Over all 365 CCRCC slides the three-level decomposition is
    within-slide 45.0% / between-slide-within-patient 16.4% / between-patient 38.6%
    in sample, and 49.0 / 12.4 / 38.6 out of fold. The spatial component is the
    largest of the three, not a tenth.

Encoding is carried over unchanged from make_composite2.py.
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import mrstyle as S

INK, GREY, LGREY = S.INK, S.GREY, S.LGREY
DD = "figdata"


def lab(ax, s, x=-0.06, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


npz = f"{DD}/wsi_spatial_data_oof.npz"
if not os.path.exists(npz):
    raise SystemExit(f"{npz} not found -- run wsi_spatial_map_oof.py on the server first")

dd = np.load(npz, allow_pickle=True)
assert bool(dd["oof"][0]), "this npz is the in-sample version; the figure must use the out-of-fold one"
vc = float(dd["vc"][0]); vspan = float(dd["vspan"][0]); fam = str(dd["family"][0])


def paint(ax, coords, ts):
    xx, yy = coords[:, 0], coords[:, 1]
    step = np.median([np.median(np.diff(np.unique(xx))), np.median(np.diff(np.unique(yy)))])
    if not np.isfinite(step) or step <= 0:
        step = 1.0
    ix = np.round((xx - xx.min()) / step).astype(int)
    iy = np.round((yy - yy.min()) / step).astype(int)
    G = np.full((iy.max() + 1, ix.max() + 1), np.nan)
    G[iy, ix] = ts
    im = ax.imshow(np.ma.masked_invalid(G), cmap="RdBu_r",
                   vmin=vc - vspan, vmax=vc + vspan, origin="upper")
    ax.axis("off")
    return im


fig = plt.figure(figsize=(7.2, 3.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.0, 0.045], wspace=0.08,
                      left=0.03, right=0.90, top=0.88, bottom=0.04)

axg = fig.add_subplot(gs[0]); lab(axg, "a", x=-0.02)
img = paint(axg, dd["coords0"], dd["scores0"])
axg.set_title("Spatial map: high-residual slide", fontsize=7.6, fontweight="bold", loc="left")

axh = fig.add_subplot(gs[1]); lab(axh, "b", x=-0.03)
paint(axh, dd["coords1"], dd["scores1"])
axh.set_title("low-residual slide", fontsize=7.6, fontweight="bold", loc="left")

cax = fig.add_subplot(gs[2])
cb = fig.colorbar(img, cax=cax)
cb.set_label(f"predicted {fam} residual (per tile)", fontsize=7)
cb.ax.tick_params(labelsize=7)

S.save_pub(fig, "Fig_wsi_spatial_supp")
print(f"wrote Fig_wsi_spatial_supp | family={fam}, "
      f"{len(dd['scores0'])} and {len(dd['scores1'])} tiles")
