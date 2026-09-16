#!/usr/bin/env python3
"""
Assemble the five-cancer summary: baseline counts, what the acquisition-batch
pack does to them, the negative control, and the scanner census -- i.e. the
table that answers the reviewer's fatal objection in one place.

Usage: atlas_summary.py [axis]
Output: scripts/atlas_summary.tsv + a readable table on stdout
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = "/public/home/fjhui/ZW"
SCRIPTS = f"{ROOT}/scripts"
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
FDR = 0.05


def rdir(c):
    return ROOT + ("/results" if c == "ccrcc" else f"/{c}/results")


def main():
    axis = sys.argv[1] if len(sys.argv) > 1 else "operator"
    rows = []

    meta_fp = f"{SCRIPTS}/slide_acquisition_meta.tsv"
    meta = pd.read_csv(meta_fp, sep="\t", dtype=str).fillna("") if os.path.exists(meta_fp) else None

    for c in COHORTS:
        row = {"cohort": c}

        base_fp = os.path.join(rdir(c), "residual_results_tumoronly.csv")
        if os.path.exists(base_fp):
            b = pd.read_csv(base_fp)
            row["genes_tested"] = len(b)
            row["n_patients"] = int(b["n"].max()) if "n" in b else np.nan
            row["sig_baseline"] = int(((b["fdr"] < FDR) & (b["incremental_r2"] > 0)).sum())
            row["rate_baseline"] = row["sig_baseline"] / max(len(b), 1)
            row["median_incr_baseline"] = float(b["incremental_r2"].median())

        sp_fp = os.path.join(rdir(c), f"sitepack_{axis}.csv")
        if os.path.exists(sp_fp):
            s = pd.read_csv(sp_fp)
            need = {"fdr_batchresid", "fdr_both", "incremental_r2_batchresid",
                    "incremental_r2_both", "incremental_r2_batchonly"}
            if need.issubset(s.columns):
                sel_b = s[(s["fdr_batchresid"] < FDR) & (s["incremental_r2_batchresid"] > 0)]
                sel_o = s[(s["fdr_both"] < FDR) & (s["incremental_r2_both"] > 0)]
                row["sig_batchresid"] = len(sel_b)
                row["sig_both"] = len(sel_o)
                # Survival must be measured as an intersection, not a ratio of
                # counts. The two significant sets are NOT nested -- Jaccard runs
                # 0.15-0.49 -- so sig_batchresid/sig_baseline overstates survival in
                # every cohort and exceeded 100% in GBM, where batch removal adds
                # proteins rather than only removing them.
                base_set = set(b.loc[(b["fdr"] < FDR) & (b["incremental_r2"] > 0), "gene"])
                bres_set = set(sel_b["gene"])
                both_set = set(sel_o["gene"])
                inter = len(base_set & bres_set)
                row["shared_base_batchresid"] = inter
                row["survival_batchresid"] = inter / max(len(base_set), 1)
                row["jaccard_base_batchresid"] = (
                    inter / max(len(base_set | bres_set), 1))
                row["new_in_batchresid"] = len(bres_set - base_set)
                row["survival_both"] = (len(base_set & both_set) /
                                        max(len(base_set), 1))
                for col, name in [("incremental_r2_plain", "median_incr_plain"),
                                  ("incremental_r2_batchresid", "median_incr_batchresid"),
                                  ("incremental_r2_groupcv", "median_incr_groupcv"),
                                  ("incremental_r2_both", "median_incr_both"),
                                  ("incremental_r2_batchonly", "median_incr_batchonly")]:
                    row[name] = float(s[col].median()) if col in s else np.nan
                    row[name.replace("median_incr", "fracpos")] = (
                        float((s[col] > 0).mean()) if col in s else np.nan)
            else:
                row["sitepack_status"] = "stale_schema"

        if meta is not None:
            m = meta[meta["cohort"] == c]
            sc = m["scanner_id"].replace("", np.nan).dropna()
            row["slides"] = len(m)
            row["n_scanners"] = int(sc.nunique())
            if len(sc):
                top = sc.value_counts()
                row["top_scanner"] = top.index[0]
                row["top_scanner_frac"] = float(top.iloc[0] / len(sc))

        rows.append(row)

    df = pd.DataFrame(rows)
    out = f"{SCRIPTS}/atlas_summary.tsv"
    df.to_csv(out, sep="\t", index=False)

    def show(cols, title, fmt=None):
        have = [c for c in cols if c in df.columns]
        if len(have) <= 1:
            return
        print(f"\n=== {title} ===")
        sub = df[have].copy()
        for c in have:
            if fmt and c in fmt:
                sub[c] = sub[c].map(lambda v: fmt[c](v) if pd.notna(v) else "-")
        print(sub.to_string(index=False))

    pct = lambda v: f"{v:.1%}"
    num = lambda v: f"{v:+.4f}"

    show(["cohort", "n_patients", "genes_tested", "sig_baseline", "rate_baseline"],
         "C1 baseline, five cancers", {"rate_baseline": pct})
    show(["cohort", "slides", "n_scanners", "top_scanner", "top_scanner_frac"],
         "scanner census (the confound's precondition)", {"top_scanner_frac": pct})
    show(["cohort", "sig_baseline", "sig_batchresid", "sig_both",
          "shared_base_batchresid", "survival_batchresid", "new_in_batchresid",
          "jaccard_base_batchresid", "survival_both"],
         f"C1 under the acquisition-batch pack (axis={axis})",
         {"survival_batchresid": pct, "jaccard_base_batchresid":
          lambda v: f"{v:.2f}", "survival_both": pct})
    show(["cohort", "median_incr_plain", "median_incr_batchresid",
          "median_incr_groupcv", "median_incr_both", "median_incr_batchonly"],
         "median incremental R2 by estimator",
         {k: num for k in ["median_incr_plain", "median_incr_batchresid",
                           "median_incr_groupcv", "median_incr_both",
                           "median_incr_batchonly"]})
    show(["cohort", "fracpos_plain", "fracpos_batchresid", "fracpos_groupcv",
          "fracpos_both", "fracpos_batchonly"],
         "fraction of proteins with a positive increment",
         {k: pct for k in ["fracpos_plain", "fracpos_batchresid", "fracpos_groupcv",
                           "fracpos_both", "fracpos_batchonly"]})

    print(f"\nwrote {out}")
    if "sig_batchresid" in df.columns:
        print("\nREAD THIS AS: sig_batchresid is C1 after removing acquisition-batch "
              "structure from the embedding.\nThe pack passes where it stays > 0 and "
              "fracpos_batchonly (the negative control) stays far\nbelow "
              "fracpos_plain. If batchonly approaches plain, the signal is batch and "
              "the run is void.")


if __name__ == "__main__":
    main()
