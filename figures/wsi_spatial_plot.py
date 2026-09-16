#!/usr/bin/env python3
"""wsi_spatial_plot.py (LOCAL) — render the WSI spatial residual map from the .npz
that wsi_spatial_map.py dumps on the server. Paints per-tile predicted residual at tile
coordinates: high- vs low-residual slide, shared diverging colour scale."""
import os, numpy as np
import matplotlib.pyplot as plt
import mrstyle as S

NPZ = "figdata/wsi_spatial_data.npz"
if not os.path.exists(NPZ):
    raise SystemExit(f"missing {NPZ} -- run wsi_spatial_map.py on the server and WinSCP it back")

d = np.load(NPZ, allow_pickle=True)
files = list(d["files"]); fam = str(d["family"][0]); vc = float(d["vc"][0]); vspan = float(d["vspan"][0])
n = len(files)
labels = ["high-residual slide", "low-residual slide"] if n == 2 else [str(f) for f in files]


def paint(ax, coords, ts, title):
    x, y = coords[:, 0], coords[:, 1]
    step = np.median([np.median(np.diff(np.unique(x))), np.median(np.diff(np.unique(y)))])
    if not np.isfinite(step) or step <= 0:
        step = 1.0
    ix = np.round((x - x.min()) / step).astype(int); iy = np.round((y - y.min()) / step).astype(int)
    G = np.full((iy.max() + 1, ix.max() + 1), np.nan); G[iy, ix] = ts
    im = ax.imshow(np.ma.masked_invalid(G), cmap="RdBu_r", vmin=vc - vspan, vmax=vc + vspan, origin="upper")
    ax.axis("off"); ax.set_title(title, fontsize=7.5, fontweight="bold")
    return im

fig, axes = plt.subplots(1, n, figsize=(3.1 * n, 3.3))
axes = np.atleast_1d(axes)
for i, (ax, lab) in enumerate(zip(axes, labels)):
    im = paint(ax, d[f"coords{i}"], d[f"scores{i}"], lab)
cb = fig.colorbar(im, ax=list(axes), fraction=0.035, pad=0.02)
cb.set_label(f"predicted {fam} residual (per tile)", fontsize=6.5); cb.ax.tick_params(labelsize=5.6)
fig.suptitle("Where high-residual morphology sits (CCRCC)", fontsize=8.6, fontweight="bold")
S.save_pub(fig, "Fig_wsi_spatial")
print(f"wrote Fig_wsi_spatial  ({n} slides: {files})")
