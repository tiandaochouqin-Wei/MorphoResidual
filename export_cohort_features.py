#!/usr/bin/env python3
"""export_cohort_features.py (SERVER, per cohort, compute-only) — dump per-patient UNI
mean-embedding + measured family residual scores, for the leave-one-cancer-out (LOCO)
generalisation analysis done locally.  Run once per cohort:
    cd /public/home/fjhui/ZW/scripts
    for c in ccrcc luad ucec gbm pdac; do
      source morpho_env.sh $c
      MORPHO_COHORT=$c MORPHO_SCORES_CSV=/public/home/fjhui/ZW/results/clinical_link_${c}_scores.csv \
      MORPHO_OUT=/public/home/fjhui/ZW/results /public/home/fjhui/miniconda3/bin/python -u export_cohort_features.py
    done
Then WinSCP  /public/home/fjhui/ZW/results/loco_*.npz  (5 small files) to figures/figdata/ .
"""
import os, numpy as np, pandas as pd
import residual_analysis as RA

c = os.environ["MORPHO_COHORT"]
OUT = os.environ.get("MORPHO_OUT", ".")
sc = pd.read_csv(os.environ["MORPHO_SCORES_CSV"], index_col=0)
wsi = RA.load_wsi_embeddings()
common = [p for p in wsi.index if p in sc.index]
meas = {k: pd.to_numeric(sc.loc[common, k], errors="coerce").values.astype(float)
        for k in sc.columns if k.endswith("_meas")}
np.savez(f"{OUT}/loco_{c}.npz", patients=np.array(common),
         emb=wsi.loc[common].values.astype(np.float32), **meas)
print(f"[export] {c}: {len(common)} patients, emb {wsi.shape[1]}d, scores {list(meas)} -> {OUT}/loco_{c}.npz")
