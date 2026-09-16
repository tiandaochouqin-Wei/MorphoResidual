#!/usr/bin/env python3
"""
export_gene_covariates.py -- MorphoResidual T1-6 (NC_ROADMAP.md).

TMT-proteomics abundance/missingness confound check. Reviewers will ask: is the
morphology-predictable residual set just an artefact of which proteins are
abundantly and completely quantified by mass spectrometry, rather than a
biological pattern? This script exports, for EVERY tested gene (not only the
pinned significant set), the quantities needed to answer that locally:
    mean_log_ratio   mean protein TMT log-ratio across non-missing patients
                     (a relative-abundance proxy; these data have no absolute
                     abundance axis, only log-ratios to a pooled reference channel)
    pct_missing      % of the common rna/protein/wsi patient set with no protein value
    n                non-missing patients used for that gene's ridge fits
    r2_rna           cv_r2(own mRNA only)               -- the manuscript's baseline
    r2_wsi_only      cv_r2(WSI 20 PCs only)              -- morphology alone, no mRNA
    r2_both          cv_r2(own mRNA + WSI 20 PCs)        -- the manuscript's full model
    incremental_r2   r2_both - r2_rna                    -- must match the pinned
                     residual_results_tumoronly.csv for the significant subset (sanity
                     check printed at the end); small numerical drift from a different
                     random CV split is expected, large drift is not.

Reuses residual_analysis (RA) loaders + RA.cv_r2 -> identical estimand to the main
pipeline and to transcriptome_baseline.py (same patient-intersection, same PCA
convention, MORPHO_SIG_CSV read directly from the pinned snapshot, never via the
shared/mutable OUT_DIR).

Usage (per cohort, smp queue, same env as transcriptome_baseline.py):
  source morpho_env.sh <cohort>
  export MORPHO_SIG_CSV=/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/<cohort>/residual_results_tumoronly.csv
  export MORPHO_COV_OUT=gene_covariates_<cohort>.csv
  python -u export_gene_covariates.py

Run on all five cohorts (CCRCC/LUAD/UCEC/GBM/PDAC). Cost: ~3 cv_r2 calls per
tested gene (9,600-10,800 genes/cohort), same order of magnitude as the main
discovery run and transcriptome_baseline.py's default (fast) mode -- no
permutation loop, so this should finish well under an hour per cohort.
"""
import os, re, sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import residual_analysis as RA  # loaders + constants (env read at import)

CV_FOLDS = getattr(RA, "CV_FOLDS", 5)

EXPECTED = {  # (n_patients, n_sig_genes), from the pinned manuscript numbers (Table 1).
    "ccrcc": (103, 2191), "luad": (105, 2310), "ucec": (100, 2566),
    "gbm": (99, 702), "pdac": (137, 1572),
}


def infer_cohort(sig_csv):
    m = re.search(r"(ccrcc|luad|ucec|gbm|pdac)", sig_csv.lower())
    return m.group(1) if m else (os.environ.get("MORPHO_COHORT", "").lower() or None)


def main():
    out = os.environ.get("MORPHO_COV_OUT", "gene_covariates.csv")
    sig_csv = os.environ.get("MORPHO_SIG_CSV") or os.environ.get("MORPHO_MAIN_CSV")
    if not sig_csv or not os.path.exists(sig_csv):
        sys.exit("set MORPHO_SIG_CSV to the pinned residual_results_tumoronly.csv for this cohort")

    cohort = infer_cohort(sig_csv)
    print(f"[sanity] MORPHO_SIG_CSV: {sig_csv}")
    print(f"[sanity] cohort inferred as: '{cohort}'")

    md0 = pd.read_csv(sig_csv)
    if "gene" not in md0.columns:
        sys.exit(f"{sig_csv} has no 'gene' column -- wrong file")
    sig0 = None
    if {"fdr", "incremental_r2"}.issubset(md0.columns):
        sig0 = int(((md0["fdr"] < 0.05) & (md0["incremental_r2"] > 0)).sum())
        print(f"[sanity] pinned significant genes in this file: {sig0}")
        if cohort in EXPECTED:
            _, exp_sig = EXPECTED[cohort]
            if sig0 != exp_sig:
                sys.exit(f"[sanity] STOP: {sig_csv} has {sig0} significant genes, expected {exp_sig} for "
                          f"cohort '{cohort}' (Table 1). Check MORPHO_SIG_CSV points at the pinned snapshot "
                          f"for the RIGHT cohort, not a stale/shared file.")

    rna = RA.load_rna_matrix()
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()
    print(f"[sanity] raw loader shapes -- rna {rna.shape}, protein {protein.shape}, wsi {wsi.shape}")
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]
    n = len(common)
    if cohort in EXPECTED:
        exp_n, _ = EXPECTED[cohort]
        if abs(n - exp_n) > 2:
            sys.exit(f"[sanity] STOP: {n} common patients after intersecting rna/protein/wsi, expected "
                      f"~{exp_n} for cohort '{cohort}'. The loaders are likely returning the wrong "
                      f"cohort's data -- do not trust downstream numbers.")

    wsi_pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
    print(f"[cov] {n} patients; wsi PCs {wsi_pcs.shape}")

    tested = sorted(set(rna.columns) & set(protein.columns))
    print(f"[cov] tested genes (rna and protein overlap): {len(tested)}")

    sig_lookup = {}
    if {"fdr", "incremental_r2"}.issubset(md0.columns):
        for _, r in md0.iterrows():
            sig_lookup[str(r["gene"])] = (bool(r["fdr"] < 0.05 and r["incremental_r2"] > 0), float(r["incremental_r2"]))

    rows = []
    for gi, g in enumerate(tested):
        y = np.asarray(protein[g].values, float)
        v = ~np.isnan(y)
        n_g = int(v.sum())
        pct_missing = 100.0 * (1 - n_g / n)
        if n_g < 30:
            rows.append(dict(gene=g, n=n_g, pct_missing=round(pct_missing, 2),
                              mean_log_ratio=np.nan, r2_rna=np.nan, r2_wsi_only=np.nan,
                              r2_both=np.nan, incremental_r2=np.nan,
                              pinned_significant=sig_lookup.get(g, (np.nan, np.nan))[0],
                              pinned_incremental_r2=sig_lookup.get(g, (np.nan, np.nan))[1]))
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1)
        Wv = wsi_pcs[v]
        r2_rna = RA.cv_r2(xr, yv, CV_FOLDS)
        r2_wsi_only = RA.cv_r2(Wv, yv, CV_FOLDS)
        r2_both = RA.cv_r2(np.hstack([xr, Wv]), yv, CV_FOLDS)
        pinned_sig, pinned_incr = sig_lookup.get(g, (np.nan, np.nan))
        rows.append(dict(gene=g, n=n_g, pct_missing=round(pct_missing, 2),
                          mean_log_ratio=round(float(np.mean(yv)), 4),
                          r2_rna=round(r2_rna, 5), r2_wsi_only=round(r2_wsi_only, 5),
                          r2_both=round(r2_both, 5), incremental_r2=round(r2_both - r2_rna, 5),
                          pinned_significant=pinned_sig, pinned_incremental_r2=pinned_incr))
        if gi % 1000 == 0:
            print(f"[cov] {gi}/{len(tested)} genes done", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)

    # sanity check: recomputed incremental_r2 on the pinned-significant subset should track the
    # pinned numbers closely (same estimand, small drift expected from a fresh 5-fold split)
    chk = res.dropna(subset=["incremental_r2", "pinned_incremental_r2"])
    chk = chk[chk["pinned_significant"] == True]  # noqa: E712
    if len(chk):
        corr = np.corrcoef(chk["incremental_r2"], chk["pinned_incremental_r2"])[0, 1]
        print(f"[sanity] recomputed vs pinned incremental_r2 on the {len(chk)} pinned-significant genes: "
              f"Pearson r={corr:.3f} (expect close to 1.0; large drift means something upstream changed)")

    print(f"\n[cov] wrote {out}: {len(res)} genes x [gene,n,pct_missing,mean_log_ratio,"
          f"r2_rna,r2_wsi_only,r2_both,incremental_r2,pinned_significant,pinned_incremental_r2]")
    print("DECISIVE SUMMARY:")
    d = res.dropna(subset=["incremental_r2"])
    med_missing_all = d["pct_missing"].median()
    med_missing_sig = d.loc[d["pinned_significant"] == True, "pct_missing"].median()  # noqa: E712
    med_missing_nonsig = d.loc[d["pinned_significant"] != True, "pct_missing"].median()  # noqa: E712
    print(f"  median %missing: all tested={med_missing_all:.1f}, "
          f"pinned-significant={med_missing_sig:.1f}, non-significant={med_missing_nonsig:.1f}")
    med_ab_all = d["mean_log_ratio"].median()
    med_ab_sig = d.loc[d["pinned_significant"] == True, "mean_log_ratio"].median()  # noqa: E712
    print(f"  median mean_log_ratio (abundance proxy): all tested={med_ab_all:.3f}, "
          f"pinned-significant={med_ab_sig:.3f}")


if __name__ == "__main__":
    main()
