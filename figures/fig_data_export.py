#!/usr/bin/env python3
"""fig_data_export.py — pull RAW data from the server so the figures can carry real
data visuals (H&E tiles, per-gene scatter, embedding UMAP) like the reference papers.
Run once on mn02:  export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib ; python fig_data_export.py
Then WinSCP the whole /public/home/fjhui/ZW/figdata/ folder to your local machine."""
import os, glob
import numpy as np, pandas as pd
import openslide  # before torch
import torch
from sklearn.decomposition import PCA

ROOT = "/public/home/fjhui/ZW"
OUT = f"{ROOT}/figdata"
os.makedirs(OUT, exist_ok=True)
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
EMB = {c: (f"{ROOT}/WSI/emb/ccrcc_real" if c == "ccrcc" else f"{ROOT}/{c}/WSI/emb") for c in COH}

# 1) merged per-gene incremental R^2, UNI (pinned) vs Phikon --------------------
UNI = ROOT + "/scripts/pinned/results_snapshot_20260901/{c}/residual_results_tumoronly.csv"
PHI = ROOT + "/results_phikon/{c}/residual_results_tumoronly.csv"
for c in COH:
    fu, fp = UNI.format(c=c), PHI.format(c=c)
    if os.path.exists(fu) and os.path.exists(fp):
        u, p = pd.read_csv(fu), pd.read_csv(fp)
        m = u[["gene", "incremental_r2", "fdr"]].merge(
            p[["gene", "incremental_r2", "fdr"]], on="gene", suffixes=("_uni", "_phi"))
        m.to_csv(f"{OUT}/merged_incr_{c}.csv", index=False)
        print(f"[csv] {c}: {len(m)} genes -> merged_incr_{c}.csv")
    else:
        print(f"[csv] {c}: missing UNI/Phikon result, skipped")

# 2) real H&E tiles (variety, tissue only) -------------------------------------
TILE, SAT, FRAC = 256, 15, 0.35
def is_tissue(t):
    hsv = np.array(t.convert("HSV")); return (hsv[:, :, 1] > SAT).mean() > FRAC
tile_slides = sorted(glob.glob(f"{ROOT}/WSI/raw/ccrcc_real/*.svs"))[:2] + \
              sorted(glob.glob(f"{ROOT}/luad/WSI/raw/*.svs"))[:1] + \
              sorted(glob.glob(f"{ROOT}/pdac/WSI/raw/*.svs"))[:1]
k = 0
for sp in tile_slides:
    try:
        s = openslide.OpenSlide(sp); w, h = s.level_dimensions[0]; got = 0
        for y in range(TILE * 4, h - TILE, TILE * 9):
            for x in range(TILE * 4, w - TILE, TILE * 9):
                t = s.read_region((x, y), 0, (TILE, TILE)).convert("RGB")
                if is_tissue(t):
                    t.save(f"{OUT}/tile_{k:02d}.png"); k += 1; got += 1
                    if got >= 4:
                        break
            if got >= 4:
                break
    except Exception as e:
        print(f"[tile] {sp}: {e}")
print(f"[tile] wrote {k} tiles")

# 3) per-slide embedding, PCA-50, tagged by cohort (for a UMAP) -----------------
rows, labs = [], []
for c in COH:
    d = EMB[c]
    if not os.path.isdir(d):
        print(f"[emb] {c}: no dir {d}"); continue
    n = 0
    for pt in sorted(glob.glob(f"{d}/*.pt")):
        try:
            e = torch.load(pt, map_location="cpu")["embeddings"].float().mean(0).numpy()
            rows.append(e); labs.append(c); n += 1
        except Exception:
            pass
    print(f"[emb] {c}: {n} slides")
if rows:
    X = np.vstack(rows)
    Z = PCA(n_components=50, svd_solver="full", random_state=0).fit_transform(X)
    df = pd.DataFrame(Z, columns=[f"pc{i+1}" for i in range(Z.shape[1])])
    df.insert(0, "cohort", labs)
    df.to_csv(f"{OUT}/emb_pca50.csv", index=False)
    print(f"[emb] {X.shape[0]} slides x {X.shape[1]} -> PCA50 -> emb_pca50.csv")
print(f"\nDONE. WinSCP the folder to your local machine:\n  {OUT}")
