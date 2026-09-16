#!/usr/bin/env python3
"""wsi_spatial_map.py (SERVER, compute-only, NO matplotlib) — WSI spatial residual map.
Back-projects the SAME morphology->residual direction used in the main analysis onto every
tile of a slide, and dumps per-tile coords+scores to a small .npz for local plotting.

Method (identical estimand to clinical_link.py): RA.load_wsi_embeddings() -> per-patient mean
embeddings; PCA(20); ridge(alpha=1) fit on (wsi_pcs -> measured family residual) -> direction w
in PC space. Per tile: project onto SAME PCA, standardise with fitted stats, dot with w.
Auto-picks the highest- and lowest-residual slide. CPU only.

Run with the SAME python that runs the main analysis (has numpy/pandas/torch/sklearn):
    cd /public/home/fjhui/ZW/scripts && source morpho_env.sh ccrcc
    export MORPHO_SCORES_CSV=/public/home/fjhui/ZW/results/clinical_link_ccrcc_scores.csv
    export MORPHO_OUT=/public/home/fjhui/ZW/results
    /public/home/fjhui/miniconda3/bin/python -u wsi_spatial_map.py
Then WinSCP  $MORPHO_OUT/wsi_spatial_data.npz  back to figures/figdata/ .
"""
import os, sys, glob
import numpy as np, pandas as pd, torch
from sklearn.decomposition import PCA
import residual_analysis as RA

RIDGE_ALPHA = getattr(RA, "RIDGE_ALPHA", 1.0)
EMB = os.environ.get("MORPHO_WSI_EMB_DIR") or os.environ.get("WSI_EMB_DIR")
SCORES = os.environ.get("MORPHO_SCORES_CSV", "clinical_link_ccrcc_scores.csv")
FAMILY = os.environ.get("MORPHO_FAMILY", "translation")
OUT = os.environ.get("MORPHO_OUT", ".")
if not EMB:
    sys.exit("set MORPHO_WSI_EMB_DIR (or WSI_EMB_DIR) to the raw per-slide .pt directory")


def fit_direction():
    wsi = RA.load_wsi_embeddings()
    pca = PCA(n_components=20, svd_solver="full").fit(wsi.values)
    pcs = pca.transform(wsi.values)
    sc = pd.read_csv(SCORES, index_col=0)
    meas = pd.to_numeric(sc[f"{FAMILY}_meas"], errors="coerce").reindex(wsi.index).values
    m = np.isfinite(meas)
    X, y = pcs[m], meas[m]
    xm, xs = X.mean(0), X.std(0); xs[xs == 0] = 1.0; ym = y.mean()
    Xs = (X - xm) / xs
    w = np.linalg.solve(Xs.T @ Xs + RIDGE_ALPHA * np.eye(X.shape[1]), Xs.T @ (y - ym))
    print(f"[fit] direction from {m.sum()} patients")
    return pca, xm, xs, ym, w


def tile_scores(emb, pca, xm, xs, ym, w):
    pc = (emb - pca.mean_) @ pca.components_.T
    return ((pc - xm) / xs) @ w + ym


def load_slide(fp):
    d = torch.load(fp, map_location="cpu")
    return d["embeddings"].float().numpy(), np.asarray(d["coords"], float)


def main():
    pca, xm, xs, ym, w = fit_direction()
    args = [a if os.path.isabs(a) else f"{EMB}/{a}" for a in sys.argv[1:]]
    if args:
        chosen = args
    else:
        files = sorted(glob.glob(f"{EMB}/*.pt"))
        means = []
        for fp in files:
            try:
                emb, _ = load_slide(fp)
                means.append((float(tile_scores(emb, pca, xm, xs, ym, w).mean()), fp))
            except Exception:
                continue
        means.sort()
        chosen = [means[-1][1], means[0][1]]
        print(f"[pick] high={os.path.basename(chosen[0])} ({means[-1][0]:+.2f})  "
              f"low={os.path.basename(chosen[1])} ({means[0][0]:+.2f})  of {len(means)} slides")

    data, allv = {}, []
    files_out = []
    for i, fp in enumerate(chosen):
        emb, coords = load_slide(fp)
        ts = tile_scores(emb, pca, xm, xs, ym, w)
        data[f"coords{i}"] = coords.astype(np.float32)
        data[f"scores{i}"] = ts.astype(np.float32)
        files_out.append(os.path.basename(fp)); allv.append(ts)
    allv = np.concatenate(allv); vc = float(np.median(allv))
    data["files"] = np.array(files_out); data["family"] = np.array([FAMILY])
    data["vc"] = np.array([vc]); data["vspan"] = np.array([float(np.percentile(np.abs(allv - vc), 96))])
    path = f"{OUT}/wsi_spatial_data.npz"
    np.savez(path, **data)
    print(f"[out] {path}  ({len(chosen)} slides, {sum(len(data[f'scores{i}']) for i in range(len(chosen)))} tiles)")
    print("      WinSCP it to figures/figdata/ ; then run wsi_spatial_plot.py locally")


if __name__ == "__main__":
    main()
