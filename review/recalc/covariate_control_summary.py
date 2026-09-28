#!/usr/bin/env python3
"""Read back the cluster demographic controls and produce the numbers the manuscript needs.

Run after scp-ing covariate_control_<cohort>.csv into figures/figdata/.

VALIDATE   Blocks whose patient set is the full cohort must reproduce the pinned per-gene
           incremental_r2 exactly; that is what proves the run is the published estimand.
           Blocks that drop patients (race is unrecorded for some) cannot be checked this
           way and are reported as such.

SUMMARISE  For each cohort and block: how much of the increment survives, and how much a
           SAME-WIDTH block of pure Gaussian noise leaves. Adding any column costs
           cross-validated degrees of freedom, so only the gap between the two is
           attributable to the covariate. The comparison is made per gene and paired
           (Wilcoxon signed-rank), not just on the medians.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FD = os.path.join(ROOT, "figures", "figdata")
PIN = os.path.join(ROOT, "server_export", "pinned")
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
LABEL = {"ccrcc": "CCRCC", "luad": "LUAD", "ucec": "UCEC", "gbm": "GBM", "pdac": "PDAC"}
BLOCKS = ["sex", "age", "race", "all"]

rows, val = [], []
for c in COHORTS:
    p = os.path.join(FD, f"covariate_control_{c}.csv")
    if not os.path.exists(p):
        print(f"MISSING: {p}")
        continue
    d = pd.read_csv(p).drop_duplicates("gene")
    pin = pd.read_csv(os.path.join(PIN, c, "residual_results_tumoronly.csv")).drop_duplicates("gene")
    m = d.merge(pin[["gene", "incremental_r2", "n"]].rename(columns={"n": "n_pub"}), on="gene")
    for k in BLOCKS:
        col = f"incr_W_over_{k}"
        if col not in d.columns:
            continue
        sub = m[m[col].notna()].copy()
        if sub.empty:
            continue
        full = sub[f"n_{k}"] == sub["n_pub"]          # same patient set as the published run
        dev = (sub.loc[full, f"incr_orig_{k}"] - sub.loc[full, "incremental_r2"]).abs()
        val.append(dict(cohort=c, block=k, genes=len(sub), full_set_genes=int(full.sum()),
                        max_dev=float(dev.max()) if len(dev) else np.nan))
        med = sub[f"incr_orig_{k}"].median()
        kept = 100 * sub[col].median() / med
        nk = 100 * sub[f"incr_W_over_noise_{k}"].median() / med
        paired = sub[col] - sub[f"incr_W_over_noise_{k}"]
        w = stats.wilcoxon(paired, alternative="less") if len(paired) > 20 else None
        rows.append(dict(cohort=c, block=k, width=1 if k != "all" else np.nan,
                         patients=int(sub[f"n_{k}"].max()), genes=len(sub),
                         med_incr=med, block_alone=float(sub[f"{k}_gain"].median()),
                         noise_alone=float(sub[f"noise_gain_{k}"].median()),
                         kept_pct=kept, noise_kept_pct=nk, gap=nk - kept,
                         med_paired_diff=float(paired.median()),
                         pct_below_noise=100 * float((paired < 0).mean()),
                         wilcoxon_p=float(w.pvalue) if w else np.nan,
                         pct_genes_pos=100 * float((sub[col] > 0).mean())))

v = pd.DataFrame(val)
t = pd.DataFrame(rows)
print("=== VALIDATION: does incr_orig reproduce the published per-gene increments? ===")
print("(only genes whose block patient set equals the published set can be checked)")
for _, r in v.iterrows():
    if r.full_set_genes == 0:
        print(f"  {LABEL[r.cohort]:6s} {r.block:5s} {int(r.genes):5d} genes | patient set reduced, not checkable")
    else:
        flag = "OK" if r.max_dev < 1e-4 else "*** MISMATCH ***"
        print(f"  {LABEL[r.cohort]:6s} {r.block:5s} {int(r.genes):5d} genes | "
              f"{int(r.full_set_genes):5d} on the full set, max|diff| {r.max_dev:.2e}  {flag}")

print("\n=== COVARIATE vs SAME-WIDTH NOISE ===")
print(f"{'cohort':7s} {'block':5s} {'pts':>4s} {'genes':>6s} {'blk alone':>10s} {'noise alone':>12s} "
      f"{'kept %':>8s} {'noise %':>8s} {'gap':>7s} {'<noise':>7s} {'Wilcoxon p':>11s}")
for _, r in t.iterrows():
    print(f"{LABEL[r.cohort]:7s} {r.block:5s} {int(r.patients):4d} {int(r.genes):6d} "
          f"{r.block_alone:+10.4f} {r.noise_alone:+12.4f} {r.kept_pct:7.1f}% {r.noise_kept_pct:7.1f}% "
          f"{r.gap:+6.1f} {r.pct_below_noise:6.0f}% {r.wilcoxon_p:11.2e}")

print("\n=== per block, across cohorts ===")
for k in BLOCKS:
    s = t[t.block == k]
    if s.empty:
        continue
    print(f"  {k:5s} kept {s.kept_pct.min():.0f}-{s.kept_pct.max():.0f}% | "
          f"noise {s.noise_kept_pct.min():.0f}-{s.noise_kept_pct.max():.0f}% | "
          f"gap {s.gap.min():+.1f} to {s.gap.max():+.1f} pts | "
          f"cohorts where the covariate is clearly below noise (gap > 5 pts): "
          f"{', '.join(LABEL[x] for x in s[s.gap > 5].cohort) or 'none'}")
t.to_csv(os.path.join(ROOT, "review", "recalc", "covariate_control_summary.csv"), index=False)
print("\n-> review/recalc/covariate_control_summary.csv")
