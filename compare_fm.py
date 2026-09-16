#!/usr/bin/env python3
"""
compare_fm.py — UNI vs Phikon second-foundation-model robustness (B-8).

Reads each cohort's UNI (pinned snapshot) and Phikon residual_results_tumoronly.csv
and asks whether the morphology-residual signal is UNI-specific:
  - significant-set agreement (Jaccard, and recall of the UNI-significant set)
  - per-gene incremental-R^2 correlation (Pearson) between the two encoders
No residual_analysis import needed -- just pandas. Run anywhere the CSVs are visible:
    python compare_fm.py
"""
import os, numpy as np, pandas as pd

UNI = "/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/{c}/residual_results_tumoronly.csv"
PHI = "/public/home/fjhui/ZW/results_phikon/{c}/residual_results_tumoronly.csv"
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]


def sig_set(df):
    return set(df[(df["fdr"] < 0.05) & (df["incremental_r2"] > 0)]["gene"].astype(str))


def main():
    rows = []
    for c in COHORTS:
        fu, fp = UNI.format(c=c), PHI.format(c=c)
        if not os.path.exists(fu) or not os.path.exists(fp):
            miss = ("UNI " if not os.path.exists(fu) else "") + ("PHI" if not os.path.exists(fp) else "")
            print(f"{c}: missing {miss}")
            continue
        u, p = pd.read_csv(fu), pd.read_csv(fp)
        su, sp = sig_set(u), sig_set(p)
        inter, union = su & sp, su | sp
        m = u[["gene", "incremental_r2"]].merge(
            p[["gene", "incremental_r2"]], on="gene", suffixes=("_uni", "_phi"))
        r_all = m["incremental_r2_uni"].corr(m["incremental_r2_phi"])
        ms = m[m["gene"].isin(union)]
        r_sig = ms["incremental_r2_uni"].corr(ms["incremental_r2_phi"]) if len(ms) > 2 else float("nan")
        rows.append({
            "cohort": c, "n_uni": len(su), "n_phi": len(sp), "overlap": len(inter),
            "jaccard": len(inter) / len(union) if union else float("nan"),
            "recall_uni": len(inter) / len(su) if su else float("nan"),
            "med_incr_uni": u.loc[u["gene"].isin(su), "incremental_r2"].median(),
            "med_incr_phi": p.loc[p["gene"].isin(sp), "incremental_r2"].median(),
            "r_incr_all": r_all, "r_incr_sig": r_sig,
        })
    if not rows:
        print("no cohorts ready yet")
        return
    res = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    res.to_csv("/public/home/fjhui/ZW/results_phikon/uni_vs_phikon.csv", index=False)
    print("\nREAD: jaccard = significant-set agreement; recall_uni = fraction of UNI-significant genes")
    print("      also significant with Phikon; r_incr_* = Pearson corr of per-gene incremental R^2.")
    print("      High recall/jaccard + high r_incr => signal is NOT UNI-specific (replicates on Phikon).")
    print(" -> /public/home/fjhui/ZW/results_phikon/uni_vs_phikon.csv")


if __name__ == "__main__":
    main()
