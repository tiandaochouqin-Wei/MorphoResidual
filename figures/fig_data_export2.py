#!/usr/bin/env python3
"""fig_data_export2.py — second raw-data pull for the real-data figures:
(1) per-patient measured vs morphology-predicted residual pathway scores (all cohorts);
(2) H&E tiles from the highest- and lowest-residual CCRCC patients (Wang-style plate).
Run once on mn02:  export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib ; python fig_data_export2.py
Then WinSCP the new files in /public/home/fjhui/ZW/figdata/ back to your local figures/figdata/."""
import os, glob, shutil
import numpy as np, pandas as pd
import openslide

ROOT = "/public/home/fjhui/ZW"
OUT = f"{ROOT}/figdata"
os.makedirs(OUT, exist_ok=True)
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]

# 1) per-patient scores (measured *_meas + morphology-only *_morph) ------------
for c in COH:
    src = f"{ROOT}/results/clinical_link_{c}_scores.csv"
    if os.path.exists(src):
        shutil.copy(src, f"{OUT}/scores_{c}.csv")
        print(f"[scores] {c} -> scores_{c}.csv")
    else:
        print(f"[scores] {c}: missing {src}")

# 2) top/bottom measured-residual CCRCC tiles ---------------------------------
TILE, SAT, FRAC = 256, 15, 0.35
def is_tissue(t):
    hsv = np.array(t.convert("HSV")); return (hsv[:, :, 1] > SAT).mean() > FRAC

sc = pd.read_csv(f"{ROOT}/results/clinical_link_ccrcc_scores.csv", index_col=0)
col = "translation_meas" if "translation_meas" in sc.columns else sc.columns[0]
sc = sc.dropna(subset=[col]).sort_values(col)
lo_cases = list(sc.index[:4]); hi_cases = list(sc.index[-4:])
print(f"[tiles] rank by {col}; hi={hi_cases}  lo={lo_cases}")

def grab(case, tag, idx, ntiles=2):
    svs = sorted(glob.glob(f"{ROOT}/WSI/raw/ccrcc_real/{case}-*.svs"))
    if not svs:
        print(f"  no slide for {case}"); return
    s = openslide.OpenSlide(svs[0]); w, h = s.level_dimensions[0]; got = 0
    for y in range(TILE * 5, h - TILE, TILE * 8):
        for x in range(TILE * 5, w - TILE, TILE * 8):
            t = s.read_region((x, y), 0, (TILE, TILE)).convert("RGB")
            if is_tissue(t):
                t.save(f"{OUT}/tile_{tag}_{idx}_{got}.png"); got += 1
                if got >= ntiles:
                    return
        if got >= ntiles:
            return

for i, case in enumerate(hi_cases):
    grab(case, "hi", i)
for i, case in enumerate(lo_cases):
    grab(case, "lo", i)
print(f"[tiles] wrote hi/lo tiles for {col}")
print(f"\nDONE. WinSCP new files from {OUT} back to local figures/figdata/")
