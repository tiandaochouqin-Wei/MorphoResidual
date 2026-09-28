#!/usr/bin/env python3
"""Read back the held-out-effect-size-on-corrected-set runs and validate + summarise.

Run after scp-ing held_out_corrected_<cohort>.csv into figures/figdata/.

VALIDATE   seed-0 must reproduce the published incremental_r2_batchresid exactly for
           every gene (this is the same estimand as the batch-corrected column already
           in Table 1; only seeds 1-9 are new).

AGGREGATE  the verified convention (see held_out_corrected.py docstring): pool every
           (gene, held-out-seed) pair for seeds 1-9 and take ONE grand median. This is
           the rule that reproduces Table 1's existing Held-out R^2 column (0.026-0.113,
           computed on the uncorrected sets) to 3 d.p., so the new number is computed
           the same way and is directly comparable.
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FD = os.path.join(ROOT, "figures", "figdata")
RES = os.path.join(ROOT, "server_export", "results")
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
LABEL = {"ccrcc": "CCRCC", "luad": "LUAD", "ucec": "UCEC", "gbm": "GBM", "pdac": "PDAC"}
SITEPACK = {"ccrcc": "results__sitepack_operator.csv", "luad": "luad__sitepack_operator.csv",
           "ucec": "ucec__sitepack_operator.csv", "gbm": "gbm__sitepack_operator.csv",
           "pdac": "pdac__sitepack_operator.csv"}
PUB_UNCORRECTED_HELDOUT = {"ccrcc": 0.061, "luad": 0.032, "ucec": 0.113, "gbm": 0.026, "pdac": 0.040}
PUB_CORRECTED_COUNT = {"ccrcc": 815, "luad": 726, "ucec": 201, "gbm": 680, "pdac": 682}

rows = []
print("=== VALIDATION: seed-0 vs published incremental_r2_batchresid ===")
for c in COHORTS:
    p = os.path.join(FD, f"held_out_corrected_{c}.csv")
    if not os.path.exists(p):
        print(f"  MISSING held_out_corrected_{c}.csv")
        continue
    d = pd.read_csv(p).drop_duplicates("gene")
    site = pd.read_csv(os.path.join(RES, SITEPACK[c])).drop_duplicates("gene")
    m = d.merge(site[["gene", "incremental_r2_batchresid"]], on="gene")
    dev = (m["incr_s0"] - m["incremental_r2_batchresid"]).abs()
    ok = dev.max() < 1e-4
    print(f"  {LABEL[c]:6s} {len(m):4d} genes  max|diff| {dev.max():.2e}  {'OK' if ok else '*** MISMATCH ***'} "
          f"| genes here {len(d)} vs published corrected-significant count {PUB_CORRECTED_COUNT[c]}")
    held_cols = [f"incr_s{s}" for s in range(1, 10)]
    pooled = d[held_cols].values.ravel()
    rows.append(dict(cohort=c, n_genes=len(d), n_pub_corrected=PUB_CORRECTED_COUNT[c],
                     seed0_median=float(d["incr_s0"].median()),
                     held_out_median=float(np.nanmedian(pooled)),
                     frac_positive=100 * float((pooled > 0).mean()),
                     n_pooled_obs=len(pooled), validated=ok))

t = pd.DataFrame(rows)
if t.empty:
    raise SystemExit("no held_out_corrected_*.csv found")

print("\n=== HELD-OUT EFFECT SIZE ON THE BATCH-CORRECTED SIGNIFICANT SET ===")
print(f"{'cohort':7s} {'genes':>6s} {'seed-0':>9s} {'held-out':>10s} {'% pos':>7s} "
      f"| {'uncorrected held-out (Table 1)':>32s}")
for _, r in t.iterrows():
    print(f"{LABEL[r.cohort]:7s} {int(r.n_genes):6d} {r.seed0_median:+9.4f} {r.held_out_median:+10.4f} "
          f"{r.frac_positive:6.1f}% | {PUB_UNCORRECTED_HELDOUT[r.cohort]:+32.3f}")

print("\n=== draft sentence for the manuscript ===")
vals = ", ".join(f"{r.held_out_median:.3f}" for _, r in t.iterrows())
print(f"Re-estimated under nine held-out seeds, the batch-corrected significant proteins "
      f"retain a median incremental R^2 of {vals} for CCRCC/LUAD/UCEC/GBM/PDAC "
      f"(cf. {'/'.join(f'{v:.3f}' for v in PUB_UNCORRECTED_HELDOUT.values())} on the uncorrected sets).")
t.to_csv(os.path.join(ROOT, "review", "recalc", "held_out_corrected_summary.csv"), index=False)
print("\n-> review/recalc/held_out_corrected_summary.csv")
