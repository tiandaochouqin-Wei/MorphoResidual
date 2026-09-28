#!/usr/bin/env python3
"""Read back the operator-design width sweep and say what it means for the batch claim.

Run after scp-ing operator_width_<cohort>.csv into figures/figdata/.

VALIDATE   The 'none' design applies no correction, so incr_op_none must reproduce the
           pinned per-gene incremental_r2 exactly. If it does not, the run is not the
           published estimand and nothing below can be read.

SUMMARISE  Per cohort and design: the median increment, the median under a SAME-WIDTH
           Gaussian design, and the gap. Residualising against w columns removes w
           dimensions' worth of variance whatever those columns are, so only the gap is
           attributable to the operator axis. Also paired per gene (Wilcoxon signed-rank
           and the fraction of genes whose operator arm falls below its own noise arm),
           because a gap between two medians can hide a split population.

The designs are CAP:MIN as in residual_analysis_sitepack.py. 12:3 is the published one.
"""
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FD = os.path.join(ROOT, "figures", "figdata")
PIN = os.path.join(ROOT, "server_export", "pinned")
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
LABEL = {"ccrcc": "CCRCC", "luad": "LUAD", "ucec": "UCEC", "gbm": "GBM", "pdac": "PDAC"}
# published batch-corrected significant counts, for context
PUB = {"ccrcc": (2191, 815), "luad": (2310, 726), "ucec": (2566, 201),
       "gbm": (702, 680), "pdac": (1572, 682)}

rows = []
for c in COHORTS:
    p = os.path.join(FD, f"operator_width_{c}.csv")
    if not os.path.exists(p):
        print(f"MISSING: operator_width_{c}.csv")
        continue
    d = pd.read_csv(p).drop_duplicates("gene")
    pin = pd.read_csv(os.path.join(PIN, c, "residual_results_tumoronly.csv")).drop_duplicates("gene")
    m = d.merge(pin[["gene", "incremental_r2"]], on="gene")
    dev = (m["incr_op_none"] - m["incremental_r2"]).abs()
    print(f"VALIDATE {LABEL[c]:6s} {len(m):5d} genes  max|incr_op_none - published| = {dev.max():.2e}  "
          f"{'OK' if dev.max() < 1e-4 else '*** MISMATCH ***'}")
    ref = d["incr_op_none"].median()
    designs = [re.sub(r"^incr_op_", "", col) for col in d.columns
               if col.startswith("incr_op_") and col != "incr_op_none"]
    for g in designs:
        o, n = d[f"incr_op_{g}"], d.get(f"incr_noise_{g}")
        rec = dict(cohort=c, design=g, genes=len(d), ref=ref,
                   med_op=float(o.median()), pct_of_ref=100 * float(o.median()) / ref,
                   pct_genes_pos=100 * float((o > 0).mean()))
        if n is not None:
            paired = o - n
            rec.update(med_noise=float(n.median()), noise_pct=100 * float(n.median()) / ref,
                       gap=100 * float(n.median() - o.median()) / ref,
                       pct_below_noise=100 * float((paired < 0).mean()),
                       wilcoxon_p=float(stats.wilcoxon(paired, alternative="less").pvalue)
                       if len(paired) > 20 else np.nan,
                       med_paired_diff=float(paired.median()))
        rows.append(rec)

t = pd.DataFrame(rows)
if t.empty:
    raise SystemExit("no operator_width_*.csv found")

print("\n=== OPERATOR AXIS vs SAME-WIDTH NOISE ===")
print("(% of ref = median increment as a percentage of the uncorrected median)")
print(f"{'cohort':7s} {'design':>7s} {'genes':>6s} {'med incr':>10s} {'% ref':>7s} {'genes+':>7s} "
      f"{'noise':>9s} {'noise%':>7s} {'gap':>7s} {'<noise':>7s} {'Wilcoxon p':>11s}")
for _, r in t.iterrows():
    print(f"{LABEL[r.cohort]:7s} {r.design:>7s} {int(r.genes):6d} {r.med_op:+10.4f} "
          f"{r.pct_of_ref:6.1f}% {r.pct_genes_pos:6.0f}% {r.get('med_noise', np.nan):+9.4f} "
          f"{r.get('noise_pct', np.nan):6.1f}% {r.get('gap', np.nan):+6.1f} "
          f"{r.get('pct_below_noise', np.nan):6.0f}% {r.get('wilcoxon_p', np.nan):11.2e}")

print("\n=== the published design (12:3), against the paper's own counts ===")
pubd = t[t.design == "12:3"]
for _, r in pubd.iterrows():
    base, corr = PUB[r.cohort]
    print(f"  {LABEL[r.cohort]:6s} increment retained {r.pct_of_ref:5.1f}% of uncorrected, "
          f"noise arm {r.get('noise_pct', float('nan')):5.1f}%, gap {r.get('gap', float('nan')):+5.1f} pts | "
          f"paper reports {base} -> {corr} significant after batch correction")

print("\n=== reading ===")
for _, r in pubd.iterrows():
    gap = r.get("gap", np.nan)
    if gap != gap:
        continue
    verdict = ("the drop is the price of the columns, not the operator axis" if gap < 5 else
               "the operator axis removes signal beyond what its width costs")
    print(f"  {LABEL[r.cohort]:6s} gap {gap:+5.1f} pts -> {verdict}")
t.to_csv(os.path.join(ROOT, "review", "recalc", "operator_width_summary.csv"), index=False)
print("\n-> review/recalc/operator_width_summary.csv")
