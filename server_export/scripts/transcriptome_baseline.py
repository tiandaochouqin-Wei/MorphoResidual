#!/usr/bin/env python3
"""
transcriptome_baseline.py -- MorphoResidual P0-1 fix (adversarial review, 2026-09-06/07).

THE decisive experiment: the manuscript's baseline is the gene's OWN mRNA, a single
linear feature. Anything morphology captures that is carried by OTHER transcripts
(trans-programmes: proliferation, secretory differentiation, purity) gets counted as
"post-transcriptional" / "transcriptome-invisible" even though a transcriptome-wide
baseline would have caught it. This script re-estimates the morphology increment over
a stronger baseline and reports how many of the pinned significant genes survive.

DEFAULT (fast, matches composition_controls.py's own convention exactly -- point
estimates, no permutation, no FDR re-derivation):
  (A) own mRNA + top-K RNA-seq PCA components (K=10, 20), fit on the FULL transcriptome
      once (not per-gene). Reports median incremental R2 retained and % of genes still
      positive, same style as composition_controls.py's stroma/immune/proliferation
      numbers already in the manuscript (Results sec:ablres, "retained... % genes +").
      This is cheap: ~2,000 genes x 5 cv_r2 calls, same order of magnitude as
      composition_controls.py, which already ran fine on this cluster.

OPTIONAL, more expensive, OFF by default (set the env var to turn on):
  MORPHO_TXOME_PERM=1        also computes a 1000-permutation null + BH-FDR for the
                             K=20 increment (honest FDR<0.05 headline number, but ~1000x
                             the compute of the default -- budget accordingly, or run it
                             on a random subsample via MORPHO_TXOME_PERM_N=<n_genes>).
  MORPHO_TXOME_RIDGE=1       also fits a transcriptome-wide ridge (nested CV, all ~9,500
                             genes as features) as an exploratory, more powerful but more
                             fragile baseline (B). Restricted to the top
                             MORPHO_TXOME_RIDGE_N (default 200) genes by original
                             incremental R2, because a nested ridge on the full
                             transcriptome for every significant gene would be far too
                             slow to run for the whole significant set.

Reuses residual_analysis (RA) loaders + cv_r2 -> identical estimand to the main
pipeline (same pattern as composition_controls.py; MORPHO_SIG_CSV read directly, never
via the shared/mutable OUT_DIR that caused the confound_check.py mixup on 2026-09-06).

Usage (per cohort, smp queue, same env as composition_controls.py):
  source morpho_env.sh <cohort>
  export MORPHO_SIG_CSV=/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/<cohort>/residual_results_tumoronly.csv
  export MORPHO_TXOME_OUT=transcriptome_baseline_<cohort>.csv
  python -u transcriptome_baseline.py

Run the default (fast) version on all five cohorts first (CCRCC/LUAD/UCEC/GBM/PDAC).
The headline number the manuscript needs: for each cohort, the median incremental R2
retained and % of pinned significant genes still positive over own-mRNA+20 RNA PCs.
If that stays high (as composition_controls.py found for stroma/immune, 70-99%
retained) the "beyond-transcriptome" framing survives largely intact. If it collapses
toward the 3-6% seen when adding just 3 marker scores in LUAD/GBM, the title needs to
change to "beyond the gene's own mRNA" -- see review/REVIEW_REPORT.md P0-1 for the
exact rewrite either way. Only turn on MORPHO_TXOME_PERM/RIDGE afterwards, if the fast
result is borderline and you need the extra rigor to decide.
"""
import os, re, sys, numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

import residual_analysis as RA

CV_FOLDS = getattr(RA, "CV_FOLDS", 5)
DO_PERM = os.environ.get("MORPHO_TXOME_PERM", "0") == "1"
N_PERM = int(os.environ.get("MORPHO_TXOME_NPERM", "1000"))
PERM_N_GENES = os.environ.get("MORPHO_TXOME_PERM_N")  # optional: only perm-test a subsample
DO_RIDGE = os.environ.get("MORPHO_TXOME_RIDGE", "0") == "1"
RIDGE_N_GENES = int(os.environ.get("MORPHO_TXOME_RIDGE_N", "200"))
RIDGE_ALPHAS = np.logspace(-1, 3, 9)

EXPECTED = {  # (n_patients, n_sig_genes), from the pinned manuscript numbers (Table 1).
    "ccrcc": (103, 2191), "luad": (105, 2310), "ucec": (100, 2566),
    "gbm": (99, 702), "pdac": (137, 1572),
}


def infer_cohort(sig_csv):
    m = re.search(r"(ccrcc|luad|ucec|gbm|pdac)", sig_csv.lower())
    return m.group(1) if m else (os.environ.get("MORPHO_COHORT", "").lower() or None)


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
    return q


def txome_ridge_cv_r2(Xfull, y, folds=CV_FOLDS, seed=0):
    """Nested CV: outer loop for the honest R2 (like RA.cv_r2), inner RidgeCV picks
    alpha on the training fold only, to avoid leakage on a p>>n design matrix."""
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    preds = np.zeros_like(y, dtype=float)
    for tr, te in kf.split(Xfull):
        Xtr, Xte, ytr = Xfull[tr], Xfull[te], y[tr]
        xm = Xtr.mean(0); xs = Xtr.std(0); xs[xs == 0] = 1.0
        Xtr_s = (Xtr - xm) / xs; Xte_s = (Xte - xm) / xs; ym = ytr.mean()
        inner = RidgeCV(alphas=RIDGE_ALPHAS, cv=min(5, len(tr) - 1))
        inner.fit(Xtr_s, ytr - ym)
        preds[te] = inner.predict(Xte_s) + ym
    ss_res = np.sum((y - preds) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def main():
    out = os.environ.get("MORPHO_TXOME_OUT", "transcriptome_baseline.csv")
    sig_csv = os.environ.get("MORPHO_SIG_CSV") or os.environ.get("MORPHO_MAIN_CSV")
    if not sig_csv:
        sys.exit("set MORPHO_SIG_CSV to the pinned residual_results_tumoronly.csv for this cohort")

    # ---- sanity checks (2026-09-07: confound_check.py's shared OUT_DIR silently returned
    # another job's data once; this reads MORPHO_SIG_CSV directly so should be safe, but the
    # RA loaders are shared code, so verify before spending compute) ----
    cohort = infer_cohort(sig_csv)
    md0 = pd.read_csv(sig_csv)
    sig0 = ((md0["fdr"] < 0.05) & (md0["incremental_r2"] > 0)).sum() if {"fdr", "incremental_r2"}.issubset(md0.columns) else None
    print(f"[sanity] MORPHO_SIG_CSV: {sig_csv}")
    print(f"[sanity] cohort inferred as: {cohort!r}; sig_csv shape={md0.shape}, significant genes in it={sig0}")
    if cohort in EXPECTED:
        exp_n, exp_sig = EXPECTED[cohort]
        if sig0 != exp_sig:
            sys.exit(f"[sanity] STOP: {sig_csv} has {sig0} significant genes, expected {exp_sig} for "
                      f"cohort '{cohort}' (Table 1). Check MORPHO_SIG_CSV points at the pinned snapshot "
                      f"for the RIGHT cohort, not a stale/shared file.")
    else:
        print("[sanity] WARNING: could not identify cohort from the path/MORPHO_COHORT -- "
              "skipping the significant-gene-count cross-check. Proceed with extra caution.")

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
    print(f"[txome] {n} patients; wsi PCs {wsi_pcs.shape}; RNA matrix {rna.shape}")

    # ---- baseline (A): transcriptome-wide PCA (fit ONCE on the full RNA matrix, not per-gene) ----
    rna_vals = rna.values.astype(float)
    rna_z = (rna_vals - rna_vals.mean(0)) / np.where(rna_vals.std(0) == 0, 1.0, rna_vals.std(0))
    K_MAX = min(20, n - 10, rna_z.shape[1])
    if K_MAX < 10:
        sys.exit(f"too few patients ({n}) or genes for a 10-20-component RNA PCA")
    rna_pca = PCA(n_components=K_MAX, svd_solver="full").fit_transform(rna_z)
    rna_pcs10 = rna_pca[:, :min(10, K_MAX)]
    rna_pcs20 = rna_pca[:, :K_MAX]
    print(f"[txome] transcriptome-wide PCA: using top {rna_pcs10.shape[1]} / {rna_pcs20.shape[1]} components (K_MAX={K_MAX})")

    # 3 independent noise draws, averaged, so the comparator is not one lucky/unlucky matrix
    N_NOISE_DRAWS = int(os.environ.get("MORPHO_TXOME_NOISE_DRAWS", "3"))
    noise_mats = [np.random.RandomState(42 + i).randn(n, rna_pcs20.shape[1]) for i in range(N_NOISE_DRAWS)]

    # ---- canonical correlations between the morphology block and the RNA-PC block ----
    # Decisive for whether the "excess over noise" column can be read QUANTITATIVELY. An isotropic
    # noise block is independent of the enlarged baseline, whereas morphology PCs may be partly
    # collinear with the RNA PCs; a near-collinear block costs fewer effective degrees of freedom
    # under ridge, so it can beat independent noise even carrying no extra information. High
    # canonical correlations => the excess column is an upper bound only.
    def _cancorr(A, B):
        Az = (A - A.mean(0)) / np.where(A.std(0) == 0, 1.0, A.std(0))
        Bz = (B - B.mean(0)) / np.where(B.std(0) == 0, 1.0, B.std(0))
        Qa, _ = np.linalg.qr(Az); Qb, _ = np.linalg.qr(Bz)
        return np.clip(np.linalg.svd(Qa.T @ Qb, compute_uv=False), 0, 1)
    cc = _cancorr(wsi_pcs, rna_pcs20)
    cc_noise = _cancorr(noise_mats[0], rna_pcs20)
    print(f"[cancorr] morphology PCs vs RNA PCs : max {cc[0]:.3f}, median {np.median(cc):.3f}, "
          f"mean {cc.mean():.3f}, #>0.5: {(cc > 0.5).sum()}/{len(cc)}")
    print(f"[cancorr] SAME-SIZE NOISE vs RNA PCs: max {cc_noise[0]:.3f}, median {np.median(cc_noise):.3f}, "
          f"mean {cc_noise.mean():.3f}, #>0.5: {(cc_noise > 0.5).sum()}/{len(cc_noise)}   <-- the honest reference")
    print("[cancorr] if morphology's canonical correlations are FAR above the noise reference, the")
    print("[cancorr] 'EXCESS over noise' column is inflated by collinearity and is an UPPER BOUND only.")

    if {"fdr", "incremental_r2"}.issubset(md0.columns):
        sig_df = md0[(md0["fdr"] < 0.05) & (md0["incremental_r2"] > 0)][["gene", "incremental_r2"]]
    else:
        sig_df = md0[["gene"]].assign(incremental_r2=np.nan)
    sig_df = sig_df[sig_df["gene"].isin(rna.columns) & sig_df["gene"].isin(protein.columns)]
    sig = sig_df["gene"].astype(str).tolist()
    print(f"[txome] re-testing {len(sig)} pinned significant genes over transcriptome-wide baselines "
          f"(permutation FDR: {'ON, N=' + str(N_PERM) if DO_PERM else 'OFF (fast mode)'}; "
          f"exploratory ridge baseline: {'ON, top ' + str(RIDGE_N_GENES) + ' genes' if DO_RIDGE else 'OFF'})")
    if sig0 is not None and len(sig) < 0.95 * sig0:
        sys.exit(f"[sanity] STOP: only {len(sig)}/{sig0} pinned significant genes survived intersection "
                  f"with the freshly-loaded rna/protein columns -- bigger drop than expected, suggests a "
                  f"gene-naming or cohort mismatch between MORPHO_SIG_CSV and the loaders.")

    perm = None
    if DO_PERM:
        perm = [np.random.RandomState(i).permutation(n) for i in range(N_PERM)]
    ridge_genes = set()
    if DO_RIDGE:
        ridge_genes = set(sig_df.sort_values("incremental_r2", ascending=False)["gene"].head(RIDGE_N_GENES))

    # cheap capacity control (2026-09-07): adding K more dimensions to an already n~100-patient
    # regression costs OOF R2 by itself (the manuscript's own PC-count ablation found a 40-PC
    # baseline -> +20 more overfits to median -0.43 with ZERO real signal). Comparing the real
    # WSI-PC increment against this SAME-DIMENSION random-noise increment (fixed seed, generated
    # once, NOT per-gene) isolates genuine information loss to the transcriptome baseline from
    # plain capacity cost -- cheap (one extra cv_r2 call per gene, no permutation loop needed for
    # this aggregate-level check).

    rows = []
    for gi, g in enumerate(sig):
        y = np.asarray(protein[g].values, float)
        v = ~np.isnan(y)
        if v.sum() < 30:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1)
        Wv = wsi_pcs[v]
        P10v = rna_pcs10[v]; P20v = rna_pcs20[v]
        r2_rna = RA.cv_r2(xr, yv, CV_FOLDS)
        r2_p10 = RA.cv_r2(np.hstack([xr, P10v]), yv, CV_FOLDS)
        r2_p20 = RA.cv_r2(np.hstack([xr, P20v]), yv, CV_FOLDS)
        r2_p10w = RA.cv_r2(np.hstack([xr, P10v, Wv]), yv, CV_FOLDS)
        r2_p20w = RA.cv_r2(np.hstack([xr, P20v, Wv]), yv, CV_FOLDS)
        # WSI added to the ORIGINAL (own-mRNA-only) baseline -- this is the manuscript's own estimand
        r2_rnaw = RA.cv_r2(np.hstack([xr, Wv]), yv, CV_FOLDS)
        # same-size noise added to BOTH baselines, so morphology's noise-adjusted advantage can be
        # compared like-for-like before vs after conditioning on the transcriptome
        r2_p20noise = np.mean([RA.cv_r2(np.hstack([xr, P20v, Nm[v]]), yv, CV_FOLDS) for Nm in noise_mats])
        r2_rnanoise = np.mean([RA.cv_r2(np.hstack([xr, Nm[v]]), yv, CV_FOLDS) for Nm in noise_mats])
        incr_p10 = r2_p10w - r2_p10
        incr_p20 = r2_p20w - r2_p20
        incr_p20_randctrl = r2_p20noise - r2_p20
        incr_own = r2_rnaw - r2_rna                     # manuscript's published increment
        incr_own_randctrl = r2_rnanoise - r2_rna        # its capacity cost (cf. Fig 1's -0.29 in CCRCC)

        row = {"gene": g, "n": int(v.sum()), "r2_rna_only": r2_rna,
               "r2_rna_pcs10": r2_p10, "r2_rna_pcs20": r2_p20,
               "incr_over_pcs10": incr_p10, "incr_over_pcs20": incr_p20,
               "incr_over_pcs20_randctrl": incr_p20_randctrl,
               "incr_over_own_mrna": incr_own, "incr_over_own_mrna_randctrl": incr_own_randctrl}

        if DO_PERM and (PERM_N_GENES is None or gi < int(PERM_N_GENES)):
            null20 = np.array([RA.cv_r2(np.hstack([xr, P20v, wsi_pcs[p][v]]), yv, CV_FOLDS) - r2_p20 for p in perm])
            row["pval_over_pcs20"] = (np.sum(null20 >= incr_p20) + 1) / (len(perm) + 1)

        if DO_RIDGE and g in ridge_genes:
            Xtx = rna_z[v]  # full transcriptome, standardised
            r2_tx = txome_ridge_cv_r2(Xtx, yv)
            r2_txw = txome_ridge_cv_r2(np.hstack([Xtx, Wv]), yv)
            row["r2_txome_ridge"] = r2_tx
            row["incr_over_txome_ridge"] = r2_txw - r2_tx

        rows.append(row)
        if (gi + 1) % 200 == 0:
            print(f"[txome]  ... {gi+1}/{len(sig)} genes done")

    res = pd.DataFrame(rows)
    if "pval_over_pcs20" in res.columns:
        m = res["pval_over_pcs20"].notna()
        res.loc[m, "fdr_over_pcs20"] = bh_fdr(res.loc[m, "pval_over_pcs20"].values)
    res.to_csv(out, index=False)

    print("\n==================== P0-1 TRANSCRIPTOME-WIDE BASELINE ====================")
    print(f" pinned significant genes re-tested         : {len(res)}")
    print(f" median incr. over own-mRNA+10 RNA PCs      : {res['incr_over_pcs10'].median():+.4f}  "
          f"({100*(res['incr_over_pcs10']>0).mean():.0f}% still positive)")
    print(f" median incr. over own-mRNA+20 RNA PCs      : {res['incr_over_pcs20'].median():+.4f}  "
          f"({100*(res['incr_over_pcs20']>0).mean():.0f}% still positive)  <-- HEADLINE NUMBER (raw, capacity-CONFOUNDED)")
    print(f" median incr. over own-mRNA+20 SAME-SIZE RANDOM NOISE: {res['incr_over_pcs20_randctrl'].median():+.4f}  "
          f"({100*(res['incr_over_pcs20_randctrl']>0).mean():.0f}% still positive)  <-- CAPACITY CONTROL")
    excess = res["incr_over_pcs20"] - res["incr_over_pcs20_randctrl"]
    print(f" median EXCESS of real WSI over random noise (capacity-corrected): {excess.median():+.4f}  "
          f"({100*(excess > 0).mean():.0f}% of genes where real WSI beats random noise)  <-- READ THIS, NOT THE RAW HEADLINE")
    print(" (if the raw headline % is low but the capacity control is ALSO similarly low/negative,")
    print("  the collapse is mostly capacity/overfitting cost, not real information loss -- read the")
    print("  'EXCESS' line instead, which nets that cost out and is the fair comparison.)")

    # ---- THE decisive comparison: morphology's noise-adjusted advantage, before vs after
    # conditioning on the transcriptome. Same quantity, same estimator, two baselines. ----
    adv_own = (res["incr_over_own_mrna"] - res["incr_over_own_mrna_randctrl"]).median()
    adv_txome = excess.median()
    print("\n --- P0-1 DECISIVE COMPARISON (morphology's advantage over same-size noise) ---")
    print(f" baseline = own mRNA only            : median incr {res['incr_over_own_mrna'].median():+.4f}, "
          f"noise {res['incr_over_own_mrna_randctrl'].median():+.4f}  -> advantage {adv_own:+.4f}")
    print(f" baseline = own mRNA + 20 RNA PCs    : median incr {res['incr_over_pcs20'].median():+.4f}, "
          f"noise {res['incr_over_pcs20_randctrl'].median():+.4f}  -> advantage {adv_txome:+.4f}")
    if adv_own > 0:
        print(f" ==> morphology RETAINS {100*adv_txome/adv_own:.0f}% of its noise-adjusted advantage "
              f"after conditioning on the transcriptome.")
        print("     (high % -> the signal is genuinely not reducible to transcriptome-wide axes;")
        print("      low % -> most of what was called 'transcriptome-invisible' is carried by other transcripts.)")
    if "fdr_over_pcs20" in res.columns:
        sig20 = res[(res["fdr_over_pcs20"] < 0.05) & (res["incr_over_pcs20"] > 0)]
        print(f" STILL FDR<0.05 & positive over 20 RNA PCs  : {len(sig20)} / {len(res)} = "
              f"{100*len(sig20)/max(len(res),1):.0f}%")
    if "incr_over_txome_ridge" in res.columns:
        sub = res[res["incr_over_txome_ridge"].notna()]
        print(f" [exploratory, top {len(sub)} genes] median incr. over transcriptome-wide ridge: "
              f"{sub['incr_over_txome_ridge'].median():+.4f}  ({100*(sub['incr_over_txome_ridge']>0).mean():.0f}% still positive)")
    print(f" -> {out}")
    print("\nVERDICT GUIDE: if the headline % is high (as composition_controls.py found for stroma/immune,")
    print("70-99% retained), the 'beyond-transcriptome' framing survives largely intact.")
    print("If it collapses toward the 3-6% seen when adding just 3 marker scores in LUAD/GBM,")
    print("retitle to 'beyond the gene's own mRNA' and rewrite per review/REVIEW_REPORT.md P0-1.")


if __name__ == "__main__":
    main()
