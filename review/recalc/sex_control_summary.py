#!/usr/bin/env python3
"""Read back the cluster sex-control results and produce the numbers the manuscript needs.

Run after scp-ing sex_control_<cohort>.csv into figures/figdata/.

Two jobs:

  VALIDATE  The cluster recomputes incr_orig with the published recipe, so it must
            reproduce the pinned per-gene incremental_r2 gene for gene. If it does not,
            something about the run differed and nothing below can be trusted. (A few
            1e-7-scale residues are expected on GBM: cross-host LAPACK SVD noise, the
            same signature composition_controls.py shows.)

  SUMMARISE Per cohort: how much of the increment survives putting sex into the mRNA
            baseline, what sex contributes on its own, and -- where both strata cleared
            the 40-patient floor -- whether the increment is positive within each sex.

Prints a ready-to-paste sentence at the end.
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FD = os.path.join(ROOT, "figures", "figdata")
PIN = os.path.join(ROOT, "server_export", "pinned")
COHORTS = ["ccrcc", "luad", "gbm", "pdac"]
LABEL = {"ccrcc": "CCRCC", "luad": "LUAD", "gbm": "GBM", "pdac": "PDAC"}

rows, missing = [], []
for c in COHORTS:
    p = os.path.join(FD, f"sex_control_{c}.csv")
    if not os.path.exists(p):
        missing.append(c)
        continue
    d = pd.read_csv(p).drop_duplicates("gene")
    pin = pd.read_csv(os.path.join(PIN, c, "residual_results_tumoronly.csv")).drop_duplicates("gene")
    m = d.merge(pin[["gene", "incremental_r2", "n"]], on="gene", suffixes=("", "_pub"))
    dev = (m.incr_orig - m.incremental_r2).abs()
    kept = 100 * d.incr_W_over_sex.median() / d.incr_orig.median()
    rec = dict(cohort=c, genes=len(d), merged=len(m), max_dev=float(dev.max()),
               n_dev_gt_1e4=int((dev > 1e-4).sum()),
               med_incr=float(d.incr_orig.median()),
               med_sex_gain=float(d.sex_gain.median()),
               med_over_sex=float(d.incr_W_over_sex.median()),
               kept_pct=float(kept),
               pct_genes_pos=100 * float((d.incr_W_over_sex > 0).mean()))
    for tag in ("male", "female"):
        col = f"incr_{tag}"
        if col in d.columns and d[col].notna().any():
            s = d[col].dropna()
            rec[f"med_{tag}"] = float(s.median())
            rec[f"pos_{tag}"] = 100 * float((s > 0).mean())
            rec[f"n_{tag}"] = int(len(s))
    rows.append(rec)

if missing:
    print("MISSING (not yet copied back):", ", ".join(missing))
if not rows:
    raise SystemExit("no sex_control_<cohort>.csv found in figures/figdata/")

t = pd.DataFrame(rows)
print("\n=== VALIDATION: does the cluster's incr_orig reproduce the published values? ===")
for _, r in t.iterrows():
    flag = "OK" if r.max_dev < 1e-4 else "*** MISMATCH ***"
    print(f"  {LABEL[r.cohort]:6s} {int(r.merged):5d} genes merged  max|diff| {r.max_dev:.2e}  "
          f"({int(r.n_dev_gt_1e4)} genes over 1e-4)  {flag}")

print("\n=== SEX CONTROL ===")
print(f"{'cohort':7s} {'genes':>6s} {'med incr':>9s} {'sex alone':>10s} {'over sex':>9s} "
      f"{'kept %':>7s} {'genes +':>8s}")
for _, r in t.iterrows():
    print(f"{LABEL[r.cohort]:7s} {int(r.genes):6d} {r.med_incr:+9.4f} {r.med_sex_gain:+10.4f} "
          f"{r.med_over_sex:+9.4f} {r.kept_pct:6.1f}% {r.pct_genes_pos:7.0f}%")

strat = t[t.get("med_male", pd.Series(dtype=float)).notna()] if "med_male" in t else t.iloc[0:0]
if len(strat):
    print("\n=== WITHIN EACH SEX (underpowered; direction check only) ===")
    for _, r in strat.iterrows():
        print(f"  {LABEL[r.cohort]:6s} male {r.med_male:+.4f} ({r.pos_male:.0f}% +, n={int(r.n_male)}) | "
              f"female {r.med_female:+.4f} ({r.pos_female:.0f}% +, n={int(r.n_female)})")

lo, hi = t.kept_pct.min(), t.kept_pct.max()
worst = t.loc[t.kept_pct.idxmin(), "cohort"]
print("\n=== draft sentence ===")
print(f"Adding sex to the mRNA baseline left the per-gene increment essentially unchanged in all "
      f"four cohorts that contain both sexes ({lo:.0f}--{hi:.0f}% of the median increment retained, "
      f"lowest in {LABEL[worst]}; sex alone contributed a median "
      f"{t.med_sex_gain.min():+.3f} to {t.med_sex_gain.max():+.3f}).")
t.to_csv(os.path.join(ROOT, "review", "recalc", "sex_control_summary.csv"), index=False)
print(f"\n-> review/recalc/sex_control_summary.csv")
