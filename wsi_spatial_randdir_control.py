#!/usr/bin/env python3
r"""wsi_spatial_randdir_control.py (SERVER, compute-only) -- the control T1-8 left open.

THE QUESTION
The out-of-fold whole-slide maps show strongly spatially autocorrelated tile scores
(Moran's I = +0.76 and +0.43 on the two displayed CCRCC slides, z = +46 and +22
against a tile-label permutation null). That null is the wrong one for the claim we
would like to make. Permuting tile labels destroys ALL spatial structure, so it only
asks "is this map non-random", and the answer is trivially yes: adjacent tiles hold
adjacent tissue, so ANY quantity derived from morphology inherits the spatial
autocorrelation of the section. Ribosome content, stroma fraction and a coin flip
weighted by tissue type would all clear that bar.

The claim that would actually upgrade the figure from illustration to result is
narrower: that the morphology-to-residual direction is MORE spatially organised than
an arbitrary direction in the same morphology subspace. That needs a different null
-- random directions, not random tile positions.

THE CONTROL
For every slide, tiles are projected onto the same fold-specific PCA and standardised
exactly as in wsi_spatial_map_oof.py, giving a (tiles x 20) matrix. The true
out-of-fold ridge direction and K random unit directions in that same 20-dimensional
space are then applied to it, and Moran's I is computed for each on the slide's own
rook-adjacency graph. Random directions in the 20-PC space, not in the raw 1024-d
embedding: the PCA step is held constant, so the comparison isolates whether the
FITTED direction is special rather than whether 20 UNI components are.

Reported per slide as the percentile of the true I within its own random
distribution, then aggregated. If the residual direction carries no spatial
organisation beyond generic morphology, those percentiles are uniform on [0,1] and
their median sits at 0.5; a one-sided Wilcoxon signed-rank test against 0.5 is the
summary statistic. PC1 is scored alongside as a reference for what a strong but
non-residual morphology axis looks like.

The random directions are drawn once and shared across slides, so a direction that
happens to be spatially smooth is smooth everywhere and cannot manufacture a
per-slide advantage for the true direction by resampling.

WHAT A POSITIVE RESULT WOULD AND WOULD NOT LICENSE
It would license "the residual direction is more spatially organised than arbitrary
directions in the same subspace", which is enough to describe the maps as a result
rather than an illustration, and to rewrite Limitations (xiv) in the positive. It
would NOT license any claim about what the spatial organisation corresponds to
biologically: no expert annotation of these slides exists (Limitations (xii)).

USAGE
    cd /public/home/fjhui/ZW/scripts && source morpho_env.sh ccrcc
    export MORPHO_SCORES_CSV=/public/home/fjhui/ZW/results/clinical_link_ccrcc_scores.csv
    export MORPHO_OUT=/public/home/fjhui/ZW/results
    bsub -q smp -n 4 -o wsi_randdir.out "/public/home/fjhui/miniconda3/bin/python -u wsi_spatial_randdir_control.py"

    Optional: MORPHO_NDIR (default 50 random directions), MORPHO_MAXSLIDES (default
    all). CPU only, one slide in memory at a time. Moran's I is vectorised over
    directions, so cost is dominated by reading every tile file once, as in T1-8
    (about two minutes there).

    Then WinSCP  $MORPHO_OUT/wsi_randdir_control.csv  back to figures/figdata/ .

This imports wsi_spatial_map_oof.py rather than reimplementing it, so the true
direction tested here is bit-for-bit the one the published maps are drawn with.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as st

from wsi_spatial_map_oof import CASE_RE, EMB, build_models, load_slide

NDIR = int(os.environ.get("MORPHO_NDIR", "50"))
MAXSLIDES = int(os.environ.get("MORPHO_MAXSLIDES", "0")) or None
OUT = os.environ.get("MORPHO_OUT", ".")
RNG = np.random.default_rng(0)


def rook_pairs(coords):
    """Index arrays (i, j) for tiles adjacent on the tiling grid, each pair listed once."""
    xx, yy = coords[:, 0], coords[:, 1]
    ux, uy = np.unique(xx), np.unique(yy)
    steps = []
    if len(ux) > 1:
        steps.append(np.diff(ux))
    if len(uy) > 1:
        steps.append(np.diff(uy))
    step = np.median(np.concatenate(steps)) if steps else 1.0
    if not np.isfinite(step) or step <= 0:
        step = 1.0
    ix = np.round((xx - xx.min()) / step).astype(np.int64)
    iy = np.round((yy - yy.min()) / step).astype(np.int64)
    key = {(a, b): k for k, (a, b) in enumerate(zip(ix, iy))}
    I, J = [], []
    for (a, b), k in key.items():
        for da, db in ((1, 0), (0, 1)):          # right and down only: each pair once
            j = key.get((a + da, b + db))
            if j is not None:
                I.append(k)
                J.append(j)
    return np.asarray(I, np.int64), np.asarray(J, np.int64)


def morans_I_many(S, I, J):
    """Moran's I for every column of S (tiles x directions) on one adjacency graph."""
    Z = S - S.mean(0, keepdims=True)
    denom = (Z ** 2).sum(0)
    num = (Z[I] * Z[J]).sum(0) * 2.0          # symmetrise: each pair counts both ways
    W = 2.0 * len(I)
    n = S.shape[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        return (n / W) * num / denom


def main():
    per_case, fold_of, insample, cases, perf = build_models()
    case_set = set(cases)

    R = RNG.normal(size=(20, NDIR))
    R /= np.linalg.norm(R, axis=0, keepdims=True)
    e1 = np.zeros(20)
    e1[0] = 1.0                                # PC1 reference axis

    files = sorted(f for f in os.listdir(EMB) if f.endswith(".pt"))
    rows = []
    for fn in files:
        stem = fn[:-3]
        cm = CASE_RE.match(stem)
        if not cm or cm.group(1) not in case_set:
            continue
        case = cm.group(1)
        try:
            emb, coords = load_slide(os.path.join(EMB, fn))
        except Exception as e:                                    # noqa: BLE001
            print(f"  !! {stem}: {str(e)[:70]}")
            continue
        if len(emb) < 100:
            continue
        I, J = rook_pairs(coords)
        if len(I) < 50:
            print(f"  .. {stem}: only {len(I)} adjacent tile pairs, skipped")
            continue

        pca, xm, xs, ym, w = per_case[case]
        pc = (emb - pca.mean_) @ pca.components_.T
        Zs = (pc - xm) / xs                                        # tiles x 20

        S = Zs @ np.column_stack([w / np.linalg.norm(w), e1, R])   # true, PC1, randoms
        Is = morans_I_many(S, I, J)
        i_true, i_pc1 = Is[0], Is[1]
        i_rand = Is[2:][np.isfinite(Is[2:])]
        if not np.isfinite(i_true) or not len(i_rand):
            continue
        rows.append(dict(slide=stem, case=case, fold=fold_of[case], n_tiles=len(emb),
                         n_pairs=len(I), I_true=float(i_true), I_pc1=float(i_pc1),
                         I_rand_median=float(np.median(i_rand)),
                         I_rand_p95=float(np.percentile(i_rand, 95)),
                         percentile=float((i_rand < i_true).mean()),
                         n_rand=int(len(i_rand))))
        if len(rows) % 50 == 0:
            print(f"  {len(rows)} slides done", flush=True)
        if MAXSLIDES and len(rows) >= MAXSLIDES:
            break

    if not rows:
        sys.exit("no slides scored -- check MORPHO_WSI_EMB_DIR")
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/wsi_randdir_control.csv", index=False)

    print(f"\n[result] {len(df)} slides, {NDIR} random directions each")
    print(f"  Moran's I  true direction    median {df.I_true.median():+.3f}   "
          f"IQR {df.I_true.quantile(.25):+.3f} .. {df.I_true.quantile(.75):+.3f}")
    print(f"  Moran's I  random directions median {df.I_rand_median.median():+.3f}   "
          f"(their 95th pct, median over slides {df.I_rand_p95.median():+.3f})")
    print(f"  Moran's I  PC1 reference     median {df.I_pc1.median():+.3f}")
    print("\n  percentile of the true direction within its own random set")
    print(f"    median {df.percentile.median():.3f}    "
          f"(0.5 = indistinguishable from an arbitrary direction)")
    print(f"    slides beating every random direction: "
          f"{int((df.percentile >= 1.0).sum())}/{len(df)}")
    print(f"    slides above the random 95th percentile: "
          f"{int((df.percentile >= 0.95).sum())}/{len(df)}")
    stat, p = st.wilcoxon(df.percentile - 0.5, alternative="greater")
    print(f"    one-sided Wilcoxon vs 0.5: W={stat:.0f}, p={p:.3e}")
    if p < 0.05 and df.percentile.median() > 0.5:
        print("\n  -> RESULT: the residual direction is more spatially organised than")
        print("     arbitrary directions in the same subspace. Fig. can be described")
        print("     as a result; rewrite Limitations (xiv) in the positive.")
    else:
        print("\n  -> NULL: not distinguishable from an arbitrary direction. The map")
        print("     stays an illustration and Limitations (xiv) stands as written.")
    print(f"\nwrote {OUT}/wsi_randdir_control.csv")


if __name__ == "__main__":
    main()
