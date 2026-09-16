#!/usr/bin/env python3
"""discovery_subsample_diag.py — F11 diagnostic (fast, NO permutation), run on the
MAIN DISCOVERY cohort's own data (UNI embeddings, dense tiling, same
collection protocol as the CPTAC-3 confirmatory cohorts) instead of the
TCGA/CPTAC-2 external data tcga_subsample_diag.py used.

Question: the CPTAC-3 confirmatory cohort (C1) that could rescue the NC
external-validation story turns out to have a real tri-modal n far below the
raw PDC case count once TCIA slide availability is checked -- UCEC lands at
n=60, LUAD at n=25, CCRCC at n=29 (GBM has no TCIA slide collection at all).
tcga_subsample_diag.py already showed that TCGA-external data (Phikon,
stride-512, FFPE diagnostic slides) collapses (negative_frac_all approx 0.90-0.99)
when subsampled to n=85-112 -- but that is a DIFFERENT sample-origin story
(adjacent-section FFPE vs the omics aliquot) from CPTAC-3 confirmatory data,
which is collected under the SAME protocol as the discovery cohorts that
worked. This script asks the cheap, decisive question before committing to a
~1-week C1 download: does the DISCOVERY cohort's own data (same protocol,
same encoder, dense tiling) ALSO collapse when subsampled down to n=60, or
does same-protocol data tolerate small n better than the TCGA-external data
did? If discovery-quality data collapses too, C1 is pointless at n=60. If it
stays healthy, that is real evidence the sample-correspondence quality (not
n alone) was driving the TCGA-external nulls, and C1 is worth the week.

Reuses residual_analysis.py's own load_rna_matrix/load_protein_matrix/
load_wsi_embeddings/cv_r2 verbatim (imported, not reimplemented) so this is
reading the EXACT data the main analysis and hyperparam_ablation.py use --
same env-var contract, same tumour-only filtering, same PCA(20, full SVD).

IMPORTANT -- run exactly like hyperparam_ablation.py, in the SAME shell,
because residual_analysis.py reads MORPHO_* env vars at import time:
    cd /public/home/fjhui/ZW/scripts && source morpho_env.sh luad
    /public/home/fjhui/miniconda3/bin/python -u discovery_subsample_diag.py
This prints the resolved MORPHO_* paths first (same fields as morpho_env.sh's
own morpho_print_env) so a stale/contaminated env var is visible immediately,
before any compute -- check the printed root/wsi_emb path is really the
cohort you meant to source before trusting the numbers below it.

Output: <MORPHO_OUT>/subsample_diag_discovery_<cohort>.csv
"""
import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import residual_analysis as RA  # noqa: E402  (reuses its load_* / cv_r2 verbatim)

N_REPEATS = 20
CANDIDATE_N = [60, 85, 103, 112]  # 60 = real CPTAC-3 UCEC-confirmatory tri-modal n;
                                   # 85/103/112 = same grid as tcga_subsample_diag.py
                                   # (GBM / BRCA / PDAC external n) for cross-comparison


def run_one(rna, protein, genes, wsi_raw, idx, min_n=30):
    sub_wsi = wsi_raw[idx]
    n_pcs = min(20, len(idx) - 1, sub_wsi.shape[1])
    pcs = PCA(n_components=n_pcs, svd_solver="full").fit_transform(sub_wsi)
    incrs = []
    for g in genes:
        y = np.asarray(protein[g].values, float)[idx]
        v = ~np.isnan(y)
        if v.sum() < min_n:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[idx][v].reshape(-1, 1)
        Wv = pcs[v]
        r2r = RA.cv_r2(xr, yv, RA.CV_FOLDS)
        r2b = RA.cv_r2(np.hstack([xr, Wv]), yv, RA.CV_FOLDS)
        incrs.append(r2b - r2r)
    incrs = np.array(incrs)
    return dict(
        n=len(idx), n_pcs=n_pcs, n_tested=len(incrs),
        median_incr_all=float(np.median(incrs)) if len(incrs) else float("nan"),
        negative_frac_all=float((incrs < 0).mean()) if len(incrs) else float("nan"),
    )


def main():
    cohort = os.environ.get("MORPHO_CANCER", "?")
    print("=== resolved MORPHO_* paths (verify this is the cohort you sourced) ===")
    print(f"  MORPHO_CANCER   : {cohort}")
    print(f"  MORPHO_ROOT     : {RA.ROOT}")
    print(f"  MORPHO_PROTEIN  : {RA.PROTEIN_TSV}")
    print(f"  MORPHO_WSI_EMB  : {RA.WSI_EMB_DIR}")
    print(f"  MORPHO_OUT      : {RA.OUT_DIR}")
    print("=" * 70, flush=True)

    rna = RA.load_rna_matrix()
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()

    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    N = len(common)
    genes = sorted(set(rna.columns) & set(protein.columns))
    print(f"\n[{cohort}] full n={N}, {len(genes)} common genes, "
          f"subsample sizes {CANDIDATE_N} x {N_REPEATS} repeats", flush=True)

    rna, protein = rna.loc[common, genes], protein.loc[common, genes]
    wsi_raw = wsi.loc[common].values

    ns = [n for n in CANDIDATE_N if n < N]
    if not ns:
        sys.exit(f"full n={N} is not larger than any requested subsample size {CANDIDATE_N}")

    rows = []
    full = run_one(rna, protein, genes, wsi_raw, np.arange(N))
    full["repeat"] = -1
    full["label"] = "full_n"
    rows.append(full)
    print(f"  full-n baseline: n={full['n']} n_pcs={full['n_pcs']} n_tested={full['n_tested']} "
          f"median_incr_all={full['median_incr_all']:.4f} negative_frac_all={full['negative_frac_all']:.4f}",
          flush=True)

    for n in ns:
        for rep in range(N_REPEATS):
            rng = np.random.RandomState(1000 * n + rep)
            idx = rng.choice(N, size=n, replace=False)
            res = run_one(rna, protein, genes, wsi_raw, idx)
            res["repeat"] = rep
            res["label"] = f"n{n}"
            rows.append(res)
        sub = [r for r in rows if r["label"] == f"n{n}"]
        med = float(np.median([r["median_incr_all"] for r in sub]))
        negf = float(np.median([r["negative_frac_all"] for r in sub]))
        negf_lo = float(np.min([r["negative_frac_all"] for r in sub]))
        negf_hi = float(np.max([r["negative_frac_all"] for r in sub]))
        print(f"  n={n}: across {N_REPEATS} draws -> median(median_incr_all)={med:.4f}  "
              f"median(negative_frac_all)={negf:.4f}  range=[{negf_lo:.4f}, {negf_hi:.4f}]", flush=True)

    df = pd.DataFrame(rows)
    outfp = RA.OUT_DIR / f"subsample_diag_discovery_{cohort}.csv"
    df.to_csv(outfp, index=False)
    print(f"\n-> {outfp}")
    print("\nCompare against tcga_subsample_diag.py's TCGA-external numbers at the same n "
          "(LUAD/KIRC subsampled: negative_frac_all approx 0.90-0.96 at n=85-112).")
    print("If this discovery-protocol data stays well below that at n=60-112, sample "
          "correspondence quality (not n alone) explains the external nulls -- C1 is worth "
          "the week. If it collapses similarly, C1 at n=60 is predictably uninformative too.")


if __name__ == "__main__":
    main()
