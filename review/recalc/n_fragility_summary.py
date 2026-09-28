#!/usr/bin/env python3
"""Read back the n-fragility runs: is the median increment fragile to n, or to the folds?

Run after scp-ing n_fragility_<cohort>.csv into figures/figdata/.

VALIDATE   The reference row is the full cohort at fold seed 0, which is the published
           configuration, so its median must equal the median of the pinned per-gene
           incremental_r2 over that cohort's significant set.

READ       arm A varies only the cross-validation fold seed at full n.
           arm B subsamples patients at a fixed seed.
           Arm B's spread includes arm A's, because changing n also moves every fold
           boundary. So the question is not whether B has spread but whether it has MORE
           than A: the excess is what is attributable to losing patients.
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FD = os.path.join(ROOT, "figures", "figdata")
PIN = os.path.join(ROOT, "server_export", "pinned")
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
LABEL = {"ccrcc": "CCRCC", "luad": "LUAD", "ucec": "UCEC", "gbm": "GBM", "pdac": "PDAC"}


def iqr(s):
    return float(np.percentile(s, 75) - np.percentile(s, 25))


rows = []
print("=== VALIDATION: reference (full n, seed 0) vs the pinned median ===")
for c in COHORTS:
    p = os.path.join(FD, f"n_fragility_{c}.csv")
    if not os.path.exists(p):
        print(f"  MISSING n_fragility_{c}.csv")
        continue
    d = pd.read_csv(p)
    ref = float(d[d.arm == "reference"]["median"].iloc[0])
    pin = pd.read_csv(os.path.join(PIN, c, "residual_results_tumoronly.csv")).drop_duplicates("gene")
    sig = pin[(pin.fdr < 0.05) & (pin.incremental_r2 > 0)]
    pub = float(sig.incremental_r2.median())
    print(f"  {LABEL[c]:6s} reference {ref:+.4f}  published {pub:+.4f}  "
          f"diff {abs(ref - pub):.2e}  {'OK' if abs(ref - pub) < 1e-4 else '*** MISMATCH ***'}")
    A = d[d.arm == "A_fold_seed"]["median"]
    rec = dict(cohort=c, ref=ref, n=int(d[d.arm == "reference"]["n"].iloc[0]),
               A_median=float(A.median()), A_min=float(A.min()), A_max=float(A.max()),
               A_spread=float(A.max() - A.min()), A_iqr=iqr(A), A_n=len(A))
    for f in sorted(d[d.arm == "B_subsample"].frac.unique()):
        B = d[(d.arm == "B_subsample") & (d.frac == f)]["median"]
        k = f"B{int(round(100*f))}"
        rec[f"{k}_n"] = int(d[(d.arm == "B_subsample") & (d.frac == f)].n.iloc[0])
        rec[f"{k}_median"] = float(B.median())
        rec[f"{k}_spread"] = float(B.max() - B.min())
        rec[f"{k}_iqr"] = iqr(B)
        rec[f"{k}_ratio"] = float(B.max() - B.min()) / max(float(A.max() - A.min()), 1e-12)
        rec[f"{k}_iqr_ratio"] = iqr(B) / max(iqr(A), 1e-12)
    rows.append(rec)

t = pd.DataFrame(rows)
if t.empty:
    raise SystemExit("no n_fragility_*.csv found")

print("\n=== FOLD NOISE ALONE (arm A: full n, only the CV seed changes) ===")
print(f"{'cohort':7s} {'n':>4s} {'reference':>10s} {'A median':>9s} {'A range':>22s} "
      f"{'spread':>8s} {'as % of ref':>12s}")
for _, r in t.iterrows():
    print(f"{LABEL[r.cohort]:7s} {int(r.n):4d} {r.ref:+10.4f} {r.A_median:+9.4f} "
          f"{('[%+.4f, %+.4f]' % (r.A_min, r.A_max)):>22s} {r.A_spread:8.4f} "
          f"{100*r.A_spread/abs(r.ref):11.0f}%")

print("\n=== SUBSAMPLING (arm B), against the fold-only spread ===")
fr = [c[:-7] for c in t.columns if c.endswith("_spread") and c.startswith("B")]
print(f"{'cohort':7s} " + " ".join(f"{k+' spread/ratio':>20s}" for k in fr))
for _, r in t.iterrows():
    cells = " ".join(f"{r[k+'_spread']:8.4f} /{r[k+'_ratio']:6.1f}x" for k in fr)
    print(f"{LABEL[r.cohort]:7s} {cells}")
print("\n(IQR-based, less sensitive to a single extreme draw)")
for _, r in t.iterrows():
    cells = " ".join(f"{k}:{r[k+'_iqr_ratio']:5.1f}x" for k in fr)
    print(f"  {LABEL[r.cohort]:7s} A IQR {r.A_iqr:.4f} | {cells}")

print("\n=== verdict per cohort (on the 95% arm, IQR ratio) ===")
for _, r in t.iterrows():
    k = "B95"
    if f"{k}_iqr_ratio" not in r:
        continue
    q = r[f"{k}_iqr_ratio"]
    v = ("losing 5% of patients adds nothing beyond fold noise" if q < 1.3 else
         "losing 5% of patients adds modestly to fold noise" if q < 2 else
         "genuinely n-fragile: subsampling spread clearly exceeds fold noise")
    print(f"  {LABEL[r.cohort]:6s} IQR ratio {q:4.1f}x -> {v}")
t.to_csv(os.path.join(ROOT, "review", "recalc", "n_fragility_summary.csv"), index=False)
print("\n-> review/recalc/n_fragility_summary.csv")
