#!/usr/bin/env python3
"""
posthoc_summarise_arms.py -- read the per-arm outputs that lsf_posthoc_residual.sh
+ residual_analysis.py (unmodified) wrote for one cohort's variant-input arms
(R, F, D01..D19, Q -- built by posthoc_build_arms.py) and report, per arm:
  - genes tested
  - Sig. count (fdr < 0.05 & incremental_r2 > 0)
  - overlap (count) and Jaccard with arm R's Sig. set
  - median incremental_r2, restricted to the genes that are in R's Sig. set
and, if an F result and at least one D-arm result are present, where F's Sig.
count empirically ranks among the D-arms' Sig. counts.

POST HOC. This script computes no new statistic of its own (no permutation, no
residualisation, no CV refit) -- it only reads residual_analysis.py's own
output table and summarises it. It is governed by, and does not override,
review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md (being written now). It cannot
read, and never touches, any C1-UCEC or C1-LUAD confirmatory file.

Inputs: for each requested arm, <--posthoc-root>/<arm>/residual_results_tumoronly.csv
-- the exact, fixed filename residual_analysis.py's own main() writes
(`out_csv = OUT_DIR / "residual_results_tumoronly.csv"`; OUT_DIR is that run's
MORPHO_OUT, which lsf_posthoc_residual.sh set to <posthoc-root>/<arm>/). An arm
whose file is absent is skipped with a printed note, not a FATAL -- a cohort may
legitimately have no D arms (k == 0) or a job still queued.

sha256 note: every OTHER script in this task's input chain sha256-asserts its
inputs against a constant or a required --*-sha256 argument. That is
deliberately NOT done here for the per-arm result CSVs themselves: they are
this run's own freshly produced output, generated moments before this script
reads them, so there is no PRIOR expected hash to assert against (unlike a
pinned/frozen reference file, whose wrong vintage would corrupt an analysis
silently). Each file's sha256 is still computed and printed/recorded in the
output table below, for audit traceability, which is the sha256 rule's spirit
even where a FATAL-on-mismatch check does not apply.

Output: a table printed to stdout (arm, genes_tested, n_sig, overlap_with_R_sig,
jaccard_with_R_sig, median_incr_r2_on_R_sig, sha256, path), optionally also
written to --out-csv, and (if applicable) one line reporting F's empirical rank
among the D-arm Sig. counts.

Run (HPC or anywhere the arm CSVs are reachable -- this is plain pandas, no
torch, no MORPHO_* env needed):
  python3 -u posthoc_summarise_arms.py --cohort ucec \\
    --posthoc-root /public/home/fjhui/ZW/ucec/results/posthoc \\
    --out-csv /public/home/fjhui/ZW/ucec/results/posthoc/SUMMARY_ucec.csv
"""
import argparse
import hashlib
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

RESULT_FILENAME = "residual_results_tumoronly.csv"  # fixed name residual_analysis.py writes
REQUIRED_COLUMNS = {"gene", "incremental_r2", "fdr"}
D_ARM_RE_PREFIX = "D"


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_arm_result(posthoc_root, arm):
    path = Path(posthoc_root) / arm / RESULT_FILENAME
    if not path.is_file():
        return None, path
    df = pd.read_csv(path)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        sys.exit(f"FATAL: {path} is missing column(s) {sorted(missing)} -- not a "
                 f"residual_analysis.py output (or wrong file).")
    if not df["gene"].is_unique:
        sys.exit(f"FATAL: {path} has duplicate gene rows -- the Sig. set would be ambiguous.")
    return df, path


def sig_genes(df):
    sig = df[(df["fdr"] < 0.05) & (df["incremental_r2"] > 0)]
    return set(sig["gene"])


def jaccard(a, b):
    u = a | b
    return (len(a & b) / len(u)) if u else 1.0


def median_incr_on_set(df, gene_set):
    present = sorted(gene_set & set(df["gene"]))
    if not present:
        return float("nan"), 0
    vals = df.set_index("gene").loc[present, "incremental_r2"]
    return float(vals.median()), len(present)


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("--cohort", required=True, help="label only, for the printed header; not used to resolve paths")
    ap.add_argument("--posthoc-root", required=True,
                    help="e.g. /public/home/fjhui/ZW/<cohort_root>/results/posthoc "
                         "(what lsf_posthoc_residual.sh set MORPHO_OUT's parent to)")
    ap.add_argument("--arms", nargs="*", default=None,
                    help="defaults to R, F, Q, D01..D19 (arms whose file is absent are skipped, not an error)")
    ap.add_argument("--out-csv", default=None)
    args = ap.parse_args()

    arms = args.arms or (["R", "F", "Q"] + [f"{D_ARM_RE_PREFIX}{s:02d}" for s in range(1, 20)])

    print(f"=== posthoc_summarise_arms: cohort={args.cohort} root={args.posthoc_root} ===")
    loaded = {}
    for arm in arms:
        df, path = load_arm_result(args.posthoc_root, arm)
        if df is None:
            print(f"  (skip) {arm}: no {RESULT_FILENAME} under {path.parent}")
            continue
        loaded[arm] = (df, path)

    if "R" not in loaded:
        sys.exit(f"FATAL: arm R's {RESULT_FILENAME} was not found under {args.posthoc_root}/R -- "
                 f"every comparison below is relative to R's Sig. set, so there is nothing to "
                 f"compare against. Run lsf_posthoc_residual.sh <cohort> R first.")

    r_df, r_path = loaded["R"]
    r_sig = sig_genes(r_df)
    print(f"\narm R: {len(r_df)} genes tested, {len(r_sig)} significant "
          f"(fdr<0.05 & incremental_r2>0)  <- {r_path}")

    rows = []
    for arm, (df, path) in loaded.items():
        sig = sig_genes(df)
        med, n_covered = median_incr_on_set(df, r_sig)
        rows.append(dict(
            arm=arm,
            genes_tested=len(df),
            n_sig=len(sig),
            overlap_with_R_sig=len(sig & r_sig),
            jaccard_with_R_sig=(1.0 if arm == "R" else jaccard(sig, r_sig)),
            median_incr_r2_on_R_sig=med,
            r_sig_genes_covered=n_covered,
            sha256=sha256_of(path),
            path=str(path),
        ))
    out = pd.DataFrame(rows)

    def arm_sort_key(a):
        return (0, "") if a == "R" else (1, "") if a == "F" else (2, "") if a == "Q" else (3, a)
    out = out.iloc[sorted(range(len(out)), key=lambda i: arm_sort_key(out.loc[i, "arm"]))].reset_index(drop=True)

    print()
    print(out.drop(columns=["path"]).to_string(index=False))

    d_rows = out[out["arm"].str.match(r"^D\d\d$")]
    if "F" in out["arm"].values and len(d_rows) > 0:
        f_sig = int(out.loc[out["arm"] == "F", "n_sig"].iloc[0])
        d_sigs = sorted(int(v) for v in d_rows["n_sig"].tolist())
        n_le = sum(1 for v in d_sigs if v <= f_sig)
        n_ge = sum(1 for v in d_sigs if v >= f_sig)
        rank = sum(1 for v in d_sigs if v < f_sig) + 1  # 1 = F is the smallest of {F} u D
        print(f"\nF vs D empirical rank (descriptive; see "
              f"review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md for how this is read):")
        print(f"  F n_sig = {f_sig}")
        print(f"  D arms (n={len(d_sigs)}): min={min(d_sigs)} median={np.median(d_sigs):.1f} "
              f"max={max(d_sigs)} values={d_sigs}")
        print(f"  F ranks {rank} of {len(d_sigs) + 1} ({{F}} union D, ascending; 1 = smallest); "
              f"{n_le}/{len(d_sigs)} D-arms have n_sig<=F; {n_ge}/{len(d_sigs)} have n_sig>=F.")
        if len(d_sigs) < 19:
            print(f"  NOTE: only {len(d_sigs)}/19 D-arm results were available -- this rank is "
                  f"NOT the frozen 19-arm reading until all 19 are present.")
    elif "F" not in out["arm"].values:
        print("\n(no F arm result present -- F/D empirical-rank comparison skipped)")
    elif len(d_rows) == 0:
        print("\n(no D-arm results present -- consistent with k==0 at build time, or jobs not yet run; "
              "F/D empirical-rank comparison skipped)")

    if args.out_csv:
        out.to_csv(args.out_csv, index=False)
        print(f"\nwrote {args.out_csv}")


if __name__ == "__main__":
    main()
