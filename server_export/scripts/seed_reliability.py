#!/usr/bin/env python3
"""
seed_reliability.py -- MorphoResidual T1-7 item 2 (NC_ROADMAP.md).

WHAT THIS IS FOR
The manuscript reports several per-gene correlations of incremental R^2 between
things that ought to agree if the estimate were noiseless:
    CPTAC-CCRCC vs TCGA-KIRC        r = -0.066   (external replication)
    UNI vs Phikon / Phikon-v2 / RN50  r = 0.10-0.25 (encoder robustness)
Those numbers are currently interpreted as if 1.0 were the target. It is not. The
per-gene incremental R^2 is itself a noisy statistic -- it is a difference of two
cross-validated R^2 values estimated from ~100 patients -- and no correlation
between two such estimates can exceed the reliability of the estimate itself.
This script measures that ceiling by re-running the SAME estimand under 10
different cross-validation seeds on the SAME data, and reporting how well the
seeds agree with each other.

If the single-seed reliability turns out to be, say, 0.45, then two perfectly
concordant datasets could not correlate above ~0.45 either, an observed 0.10-0.25
is a substantial fraction of what is attainable, and -0.066 is still clearly a
null. If the reliability is near 1.0, the observed values are genuinely low and
the manuscript should say so. Either way the number is needed before those
correlations can be interpreted at all.

WHAT VARIES AND WHAT DOES NOT
Only the KFold split varies with the seed. The 20-component morphology PCA is
cohort-level and label-free in the pinned pipeline, so it is fitted ONCE here
too -- this measures the reliability of the pinned estimand, not of a different,
fold-internal one. (The fold-internal PCA variant is a separate analysis, in
loco_infold_refit.py, and answers a different question.)

No permutations. Point estimates only. This is deliberately cheap.

REPRODUCTION GATE
The ridge used here must be the same ridge the main pipeline uses, or the
reliability measured would be that of a different estimator. The script first
checks whether RA.cv_r2 accepts a seed argument and uses it if so; otherwise it
reimplements the closed-form ridge and verifies on a sample of genes that at
seed 0 it reproduces RA.cv_r2 to within 1e-6. If that check fails the script
stops rather than reporting a number that describes the wrong estimator.

USAGE (per cohort, smp queue, same env as export_gene_covariates.py)
  source morpho_env.sh <cohort>
  export MORPHO_SIG_CSV=/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/<cohort>/residual_results_tumoronly.csv
  export MORPHO_SEED_OUT=seed_reliability_<cohort>.csv
  python -u seed_reliability.py

Cost: 10 seeds x 2 ridge fits x ~10,000 genes = ~200k tiny solves per cohort;
each is a 5-fold loop over a <=21-feature, ~100-sample system. Expect minutes,
not hours -- far cheaper than export_gene_covariates.py, which also did the
20-feature WSI-only fit. Run on all five cohorts.
"""
import inspect
import os
import re
import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

import residual_analysis as RA

CV_FOLDS = getattr(RA, "CV_FOLDS", 5)
N_SEEDS = 10
LAM = 1.0
MIN_N = 30

EXPECTED = {   # (n_patients, n_sig_genes) from the pinned manuscript numbers, Table 1
    "ccrcc": (103, 2191), "luad": (105, 2310), "ucec": (100, 2566),
    "gbm": (99, 702), "pdac": (137, 1572),
}


# ----------------------------------------------------------------- the estimator
def cv_r2_seed(X, y, n_splits, seed):
    """Closed-form ridge, out-of-fold R^2, KFold(shuffle=True, random_state=seed).

    Matches the main pipeline's convention: standardise X on the training fold,
    centre y on the training fold, l2 penalty LAM on the standardised design, no
    penalty on the intercept.
    """
    X = np.asarray(X, float).reshape(len(y), -1)
    y = np.asarray(y, float)
    oof = np.empty_like(y)
    for tr, te in KFold(n_splits=n_splits, shuffle=True, random_state=seed).split(X):
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1.0
        Xtr, Xte = (X[tr] - mu) / sd, (X[te] - mu) / sd
        ym = y[tr].mean()
        w = np.linalg.solve(Xtr.T @ Xtr + LAM * np.eye(Xtr.shape[1]), Xtr.T @ (y[tr] - ym))
        oof[te] = Xte @ w + ym
    ss_res = float(((y - oof) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def resolve_scorer():
    """Prefer RA.cv_r2 with a seed kwarg; else use cv_r2_seed after a gate check."""
    try:
        params = inspect.signature(RA.cv_r2).parameters
    except (TypeError, ValueError):
        params = {}
    for kw in ("seed", "random_state", "cv_seed"):
        if kw in params:
            print(f"[gate] RA.cv_r2 exposes '{kw}'; using the pipeline's own scorer directly")
            return (lambda X, y, s, _kw=kw: RA.cv_r2(X, y, CV_FOLDS, **{_kw: s})), kw
    print("[gate] RA.cv_r2 takes no seed argument; using the local closed-form ridge, "
          "gated against RA.cv_r2 at seed 0")
    return (lambda X, y, s: cv_r2_seed(X, y, CV_FOLDS, s)), None


def infer_cohort(p):
    m = re.search(r"(ccrcc|luad|ucec|gbm|pdac)", (p or "").lower())
    return m.group(1) if m else (os.environ.get("MORPHO_COHORT", "").lower() or None)


def main():
    out = os.environ.get("MORPHO_SEED_OUT", "seed_reliability.csv")
    sig_csv = os.environ.get("MORPHO_SIG_CSV") or os.environ.get("MORPHO_MAIN_CSV")
    if not sig_csv or not os.path.exists(sig_csv):
        sys.exit("set MORPHO_SIG_CSV to the pinned residual_results_tumoronly.csv for this cohort")
    cohort = infer_cohort(sig_csv)
    print(f"[sanity] MORPHO_SIG_CSV: {sig_csv}\n[sanity] cohort inferred as: '{cohort}'")

    md = pd.read_csv(sig_csv)
    if "gene" not in md.columns:
        sys.exit(f"{sig_csv} has no 'gene' column -- wrong file")
    sig_set = set()
    if {"fdr", "incremental_r2"}.issubset(md.columns):
        sig_set = set(md.loc[(md.fdr < 0.05) & (md.incremental_r2 > 0), "gene"].astype(str))
        print(f"[sanity] pinned significant genes: {len(sig_set)}")
        if cohort in EXPECTED and len(sig_set) != EXPECTED[cohort][1]:
            sys.exit(f"[sanity] STOP: {len(sig_set)} significant genes, expected "
                     f"{EXPECTED[cohort][1]} for '{cohort}'. Wrong or stale MORPHO_SIG_CSV.")

    rna, protein, wsi = RA.load_rna_matrix(), RA.load_protein_matrix(), RA.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]
    n = len(common)
    if cohort in EXPECTED and abs(n - EXPECTED[cohort][0]) > 2:
        sys.exit(f"[sanity] STOP: {n} common patients, expected ~{EXPECTED[cohort][0]} "
                 f"for '{cohort}'. Loaders are returning the wrong cohort.")
    W = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
    print(f"[seed] {n} patients; WSI PCs {W.shape}; {N_SEEDS} seeds x {CV_FOLDS}-fold")

    score, kw = resolve_scorer()
    tested = sorted(set(rna.columns) & set(protein.columns))
    print(f"[seed] tested genes: {len(tested)}")

    # ---- gate: at seed 0 the scorer must agree with the pipeline's own cv_r2
    if kw is None:
        diffs = []
        for g in tested[::max(1, len(tested) // 200)][:200]:
            y = np.asarray(protein[g].values, float)
            v = ~np.isnan(y)
            if v.sum() < MIN_N:
                continue
            x = np.asarray(rna[g].values, float)[v].reshape(-1, 1)
            diffs.append(abs(cv_r2_seed(x, y[v], CV_FOLDS, 0) - RA.cv_r2(x, y[v], CV_FOLDS)))
        diffs = np.array(diffs)
        print(f"[gate] {len(diffs)} genes: max|local - RA.cv_r2| = {diffs.max():.3e}, "
              f"median {np.median(diffs):.3e}")
        if diffs.max() > 1e-6:
            print("\n[gate] STOP: the local ridge does not reproduce RA.cv_r2 at seed 0, so the "
                  "reliability measured here would describe a DIFFERENT estimator than the one the "
                  "manuscript reports. Nothing was written.")
            print("[gate] The source of the pipeline's own scorer follows -- send it back and the "
                  "local ridge can be aligned to it in one step, without another round trip:\n")
            try:
                print("-" * 78)
                print(inspect.getsource(RA.cv_r2))
                print("-" * 78)
                for name in ("CV_FOLDS", "RIDGE_LAMBDA", "LAMBDA", "ALPHA", "N_PCS", "RANDOM_STATE", "SEED"):
                    if hasattr(RA, name):
                        print(f"RA.{name} = {getattr(RA, name)!r}")
            except (OSError, TypeError) as exc:
                print(f"  (could not read the source: {exc}; send residual_analysis.py instead)")
            sys.exit(1)
        print("[gate] PASS -- local ridge reproduces the pipeline's scorer at seed 0")

    rows = []
    for gi, g in enumerate(tested):
        y = np.asarray(protein[g].values, float)
        v = ~np.isnan(y)
        if v.sum() < MIN_N:
            continue
        yv, xr, Wv = y[v], np.asarray(rna[g].values, float)[v].reshape(-1, 1), W[v]
        Xb = np.hstack([xr, Wv])
        rec = dict(gene=g, n=int(v.sum()), pinned_significant=g in sig_set)
        for s in range(N_SEEDS):
            rec[f"incr_s{s}"] = round(score(Xb, yv, s) - score(xr, yv, s), 6)
        rows.append(rec)
        if gi % 2000 == 0:
            print(f"[seed] {gi}/{len(tested)} genes", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)
    cols = [f"incr_s{s}" for s in range(N_SEEDS)]

    def report(d, lab):
        if len(d) < 50:
            print(f"\n  {lab}: only {len(d)} genes, skipped")
            return
        M = d[cols].values
        P = np.corrcoef(M.T)
        S = pd.DataFrame(M).corr(method="spearman").values
        iu = np.triu_indices(N_SEEDS, 1)
        rp, rs = P[iu], S[iu]
        sb = N_SEEDS * rp.mean() / (1 + (N_SEEDS - 1) * rp.mean())   # Spearman-Brown, 10-seed mean
        print(f"\n  {lab} (n={len(d)} genes)")
        print(f"    pairwise Pearson  r : mean {rp.mean():.3f}  median {np.median(rp):.3f}  "
              f"range {rp.min():.3f}-{rp.max():.3f}")
        print(f"    pairwise Spearman rho: mean {rs.mean():.3f}  median {np.median(rs):.3f}")
        print(f"    reliability of ONE seed          = {rp.mean():.3f}   <-- the ceiling on any "
              f"correlation with an independent estimate")
        print(f"    reliability of the {N_SEEDS}-seed mean = {sb:.3f}   (Spearman-Brown)")
        print(f"    max attainable |r| vs an equally noisy estimate = {rp.mean():.3f}; "
              f"vs a noiseless one = {np.sqrt(max(rp.mean(),0)):.3f}")
        med = np.nanmedian(M, 0)
        pos = (M > 0).sum(0)
        print(f"    median incremental R^2 per seed : {med.mean():.4f} +/- {med.std():.4f} "
              f"(range {med.min():.4f}-{med.max():.4f})")
        print(f"    genes with incr > 0 per seed    : {pos.mean():.0f} +/- {pos.std():.0f} "
              f"of {len(d)} ({100*pos.mean()/len(d):.1f}%)")

    print(f"\n[seed] wrote {out}: {len(res)} genes x {N_SEEDS} seeds")
    print("=" * 78)
    print("SEED-TO-SEED RELIABILITY OF THE PER-GENE INCREMENTAL R^2")
    report(res, "all tested genes")
    report(res[res.pinned_significant], "pinned-significant genes only")
    print("\n  How to use this: an observed correlation r_obs between two independent")
    print("  estimates of the same per-gene quantity is attenuated by the reliability of")
    print("  each. Disattenuated r = r_obs / sqrt(rel_1 * rel_2). Report both the observed")
    print("  value and the ceiling; do not present r_obs against an implicit target of 1.0.")


if __name__ == "__main__":
    main()
