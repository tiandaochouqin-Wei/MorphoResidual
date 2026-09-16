#!/usr/bin/env python3
"""
Acceptance test for the PCA fix, plus the cost of the fix in each cohort.

Part 1 is the one that matters: two independent LSF jobs run the same script on
the same GBM inputs. With svd_solver="full" their outputs must be identical, not
merely similar. Asserting determinism without testing it would be exactly the
kind of unverified claim that produced the original defect.

Part 2 compares each cohort's new deterministic result against the randomised
draw it replaces, so the atlas can state how much the published counts move.
"""
import hashlib
import os

import numpy as np
import pandas as pd

ROOT = "/public/home/fjhui/ZW"
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
FDR = 0.05


def rdir(c):
    return ROOT + ("/results" if c == "ccrcc" else f"/{c}/results")


def sig(df):
    return set(df.loc[(df["fdr"] < FDR) & (df["incremental_r2"] > 0), "gene"])


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


print("=" * 74)
print("PART 1 - determinism: two independent jobs, same inputs, same script")
print("=" * 74)
a = f"{ROOT}/gbm/results/residual_results_tumoronly.csv"
b = f"{ROOT}/gbm/results_detverify/residual_results_tumoronly.csv"
if not (os.path.exists(a) and os.path.exists(b)):
    print("  one of the two runs is missing; cannot verify")
else:
    ha, hb = sha(a), sha(b)
    print(f"  run A sha256 {ha[:32]}...")
    print(f"  run B sha256 {hb[:32]}...")
    if ha == hb:
        print("  RESULT: byte-identical. Determinism confirmed.")
    else:
        da, db = pd.read_csv(a), pd.read_csv(b)
        m = da.merge(db, on="gene", suffixes=("_a", "_b"))
        d = (m["incremental_r2_a"] - m["incremental_r2_b"]).abs()
        print(f"  files differ (row order or content). max |d incremental_r2| = {d.max():.3e}")
        print(f"  significant: A={len(sig(da))} B={len(sig(db))}  "
              f"identical sets: {sig(da) == sig(db)}")
        print("  RESULT: NOT byte-identical -- investigate before trusting the rerun.")

print()
print("=" * 74)
print("PART 2 - deterministic result vs the randomised draw it replaces")
print("=" * 74)
rows = []
for c in COHORTS:
    new_fp = os.path.join(rdir(c), "residual_results_tumoronly.csv")
    old_fp = new_fp + ".pre_pcafix.bak"
    if not (os.path.exists(new_fp) and os.path.exists(old_fp)):
        print(f"  {c}: missing new or backup file, skipped")
        continue
    new, old = pd.read_csv(new_fp), pd.read_csv(old_fp)
    sn, so = sig(new), sig(old)
    inter, union = len(sn & so), len(sn | so)
    m = new[["gene", "incremental_r2"]].merge(
        old[["gene", "incremental_r2"]], on="gene", suffixes=("_new", "_old"))
    r = np.corrcoef(m["incremental_r2_new"], m["incremental_r2_old"])[0, 1]
    rows.append(dict(cohort=c, tested=len(new), sig_deterministic=len(sn),
                     sig_randomised=len(so), delta=len(sn) - len(so),
                     jaccard=inter / union if union else np.nan,
                     shared=inter, pearson_incr=r))

df = pd.DataFrame(rows)
if len(df):
    out = df.copy()
    out["jaccard"] = out["jaccard"].map(lambda v: f"{v:.3f}")
    out["pearson_incr"] = out["pearson_incr"].map(lambda v: f"{v:.4f}")
    out["delta"] = out["delta"].map(lambda v: f"{v:+d}")
    print(out.to_string(index=False))
    df.to_csv(f"{ROOT}/scripts/determinism_check.tsv", sep="\t", index=False)
    print(f"\nwrote {ROOT}/scripts/determinism_check.tsv")
    print("\nFor scale: four randomised GBM draws spanned 699-704 (SD 2.2) with "
          "Jaccard 0.93-0.95,\nso a delta inside that band and a Jaccard near 0.95 "
          "means the fix changed nothing beyond\nremoving the sampling noise itself.")
