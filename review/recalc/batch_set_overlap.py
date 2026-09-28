#!/usr/bin/env python3
"""Are the batch-corrected significant sets subsets of the uncorrected ones?

Table 1's "Batch-corr." column sits beside the uncorrected "Sig." count, so a reader
takes 702 -> 680 in GBM as 97% survival. It is not: the two columns are separately
selected gene sets of comparable size, and membership turns over substantially. This
script measures the turnover so the manuscript can say so.

Sources: server_export/pinned/<cohort>/residual_results_tumoronly.csv   (uncorrected)
         server_export/results/<prefix>sitepack_operator.csv            (batch-corrected)
Both are selected by the same rule: FDR < 0.05 and a positive increment.
"""
import os

import pandas as pd

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
PIN = os.path.join(ROOT, "server_export", "pinned")
RES = os.path.join(ROOT, "server_export", "results")
COH = {"ccrcc": "results__sitepack_operator.csv", "luad": "luad__sitepack_operator.csv",
       "ucec": "ucec__sitepack_operator.csv", "gbm": "gbm__sitepack_operator.csv",
       "pdac": "pdac__sitepack_operator.csv"}
LABEL = {"ccrcc": "CCRCC", "luad": "LUAD", "ucec": "UCEC", "gbm": "GBM", "pdac": "PDAC"}
PUB = {"ccrcc": (2191, 815), "luad": (2310, 726), "ucec": (2566, 201),
       "gbm": (702, 680), "pdac": (1572, 682)}


def sig_set(df, fdr_col, incr_col):
    d = df.drop_duplicates("gene")
    return set(d.loc[(d[fdr_col] < 0.05) & (d[incr_col] > 0), "gene"].astype(str))


rows = []
for c, fn in COH.items():
    base = pd.read_csv(os.path.join(PIN, c, "residual_results_tumoronly.csv"))
    corr = pd.read_csv(os.path.join(RES, fn))
    # The sitepack output carries several selection variants. Table 1's "Batch-corr."
    # column is the batch-RESIDUALISED pair; pairing fdr_batchresid with
    # incremental_r2_plain (which startswith-matching does) silently mixes a corrected
    # FDR with an uncorrected increment and reproduces nothing.
    fdr, inc = "fdr_batchresid", "incremental_r2_batchresid"
    if fdr not in corr.columns or inc not in corr.columns:
        print(f"{LABEL[c]}: expected {fdr}/{inc}; have {list(corr.columns)}")
        continue
    B = sig_set(base, "fdr", "incremental_r2")
    C = sig_set(corr, fdr, inc)
    both, only_c = B & C, C - B
    pb, pc = PUB[c]
    rows.append(dict(cohort=c, baseline=len(B), corrected=len(C), shared=len(both),
                     new_in_corrected=len(only_c), lost=len(B - both),
                     retained_pct=100 * len(both) / len(B),
                     new_pct=100 * len(only_c) / len(C) if C else float("nan"),
                     jaccard=len(both) / len(B | C) if (B | C) else float("nan"),
                     pub_baseline=pb, pub_corrected=pc,
                     matches_pub=(len(B) == pb and len(C) == pc),
                     corr_fdr_col=fdr, corr_incr_col=inc))

t = pd.DataFrame(rows)
print("=== do the recomputed counts match the published Table 1 columns? ===")
for _, r in t.iterrows():
    print(f"  {LABEL[r.cohort]:6s} baseline {int(r.baseline):5d} (published {int(r.pub_baseline):5d})  "
          f"corrected {int(r.corrected):5d} (published {int(r.pub_corrected):5d})  "
          f"{'OK' if r.matches_pub else '*** MISMATCH ***'}")

print("\n=== set turnover between the two columns ===")
print(f"{'cohort':7s} {'baseline':>9s} {'corrected':>10s} {'shared':>7s} {'retained':>9s} "
      f"{'new':>6s} {'new %':>7s} {'Jaccard':>8s}")
for _, r in t.iterrows():
    print(f"{LABEL[r.cohort]:7s} {int(r.baseline):9d} {int(r.corrected):10d} {int(r.shared):7d} "
          f"{r.retained_pct:8.1f}% {int(r.new_in_corrected):6d} {r.new_pct:6.1f}% {r.jaccard:8.2f}")

g = t[t.cohort == "gbm"].iloc[0]
print(f"\nGBM reads as {int(g.pub_baseline)} -> {int(g.pub_corrected)}, a 3% change, but only "
      f"{int(g.shared)} of its {int(g.baseline)} baseline proteins survive ({g.retained_pct:.1f}%) "
      f"and {int(g.new_in_corrected)} of the {int(g.corrected)} corrected proteins "
      f"({g.new_pct:.0f}%) were not baseline-significant.")
t.to_csv(os.path.join(ROOT, "review", "recalc", "batch_set_overlap.csv"), index=False)
print("\n-> review/recalc/batch_set_overlap.csv")
