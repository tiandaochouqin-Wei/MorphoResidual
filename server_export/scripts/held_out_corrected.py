#!/usr/bin/env python3
"""
held_out_corrected.py -- held-out effect size on the BATCH-CORRECTED significant set.

WHY. Table 1's "Held-out R^2" column (0.026-0.113) is the held-out median for the
UNCORRECTED significant set (2,191/2,310/2,566/702/1,572 genes) -- it says nothing
about the 815/726/201/680/682 genes that survive batch correction, which is what a
reader actually cares about after §confmeth argues those are the genes that matter.
No effect-size number for the corrected set has ever been reported.

AGGREGATION CONVENTION, VERIFIED, NOT ASSUMED. Fig. 1c's caption states seed 0 selects
the significant genes ("the split used to select it, which inflates it") and nine
HELD-OUT seeds (1-9) re-estimate them. Which of several plausible ways to collapse
9 seeds x N genes into Table 1's single published number was tested by reproducing all
five published values from figures/figdata/seed_reliability_<cohort>.csv:
    mean of 9 per-seed medians     : 0.0607 / 0.0271 / 0.1131 / 0.0257 / 0.0402
    median of 9 per-seed medians   : 0.0619 / 0.0289 / 0.1142 / 0.0262 / 0.0397
    POOL all (gene, seed) pairs,
      ONE grand median             : 0.0605 / 0.0316 / 0.1129 / 0.0258 / 0.0405   <- matches
    published (Table 1)            : 0.061  / 0.032  / 0.113  / 0.026  / 0.040
"Pool all (gene, seed) pairs, one grand median" is the only rule that reproduces every
published value to 3 d.p.; this script uses it.

WHAT IT DOES. Reuses residual_analysis_sitepack.py's exact operator design (CAP=12,
MIN=3, the published one) to residualise the WSI PCs ONCE (deterministic), then
re-estimates the per-gene increment under 10 KFold seeds (0 = selection split, 1-9 =
held-out), restricted to the batch-CORRECTED significant set (fdr_batchresid<0.05 and
incremental_r2_batchresid>0, from server_export/results/<prefix>sitepack_operator.csv).

VALIDATION. At seed 0, the per-gene increment on the corrected set must reproduce
incremental_r2_batchresid exactly (max|diff| ~ 0) -- this is the same estimand as the
published batch-corrected result, just re-derived here so seeds 1-9 can be added.

Env:  MORPHO_COHORT, MORPHO_SIG_CSV (unused here but kept for interface parity),
      MORPHO_SITEPACK_CSV   the sitepack_operator.csv for this cohort
      MORPHO_HELDOUT_OUT    default held_out_corrected_<cohort>.csv
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

MIN_PATIENTS = 30
RIDGE_ALPHA = 1.0
N_SEEDS = 10                  # seed 0 = selection split, 1-9 = held-out
CAP, MIN_GROUP = 12, 3        # the published operator design (residual_analysis_sitepack.py)


def cv_r2_seed(X, y, folds, seed):
    """Bit-identical to residual_analysis.cv_r2 / n_fragility.cv_r2_seed, seed exposed."""
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    preds = np.zeros_like(y, dtype=float)
    ridge_I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        X_tr, X_te, y_tr = X[tr], X[te], y[tr]
        mu, sd = X_tr.mean(axis=0), X_tr.std(axis=0)
        sd[sd == 0] = 1.0
        Xtr_s, Xte_s = (X_tr - mu) / sd, (X_te - mu) / sd
        ym = y_tr.mean()
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + ridge_I, Xtr_s.T @ (y_tr - ym))
        preds[te] = Xte_s @ w + ym
    ss_res = np.sum((y - preds) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def operator_design(raw, counts):
    """Published rule: dummy-code the CAP largest batches with >= MIN_GROUP cases."""
    levels = sorted(counts)
    dummy_levels = [lv for lv, k in counts.most_common(CAP) if k >= MIN_GROUP]
    if len(dummy_levels) == len(levels) and dummy_levels:
        dummy_levels = dummy_levels[:-1]
    D = np.zeros((len(raw), len(dummy_levels)))
    for j, lv in enumerate(dummy_levels):
        D[:, j] = (raw == lv).astype(float)
    return D


def residualise(pcs, D):
    B = np.hstack([np.ones((len(pcs), 1)), D])
    coef, *_ = np.linalg.lstsq(B, pcs, rcond=None)
    return pcs - B @ coef


def main():
    cohort = os.environ.get("MORPHO_COHORT", "").strip().lower()
    site_csv = os.environ.get("MORPHO_SITEPACK_CSV")
    if not cohort or not site_csv:
        sys.exit("[hoc] set MORPHO_COHORT and MORPHO_SITEPACK_CSV")
    out = os.environ.get("MORPHO_HELDOUT_OUT", f"held_out_corrected_{cohort}.csv")

    import batch_leak_check as blc
    paths = blc.cohort_paths(cohort)
    os.environ.update(
        MORPHO_ROOT=paths["root"], MORPHO_RNA_MANIFEST=paths["rna_manifest"],
        MORPHO_ALIQUOT_XWALK=paths["xwalk"], MORPHO_SLIDE_MAP=paths["slide_map"],
        MORPHO_PROTEIN_TSV=paths["protein"], MORPHO_WSI_EMB_DIR=paths["emb"])
    import residual_analysis as ra

    rna, protein, wsi = ra.load_rna_matrix(), ra.load_protein_matrix(), ra.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    gcols = sorted(set(rna.columns) & set(protein.columns))
    rna, protein = rna.loc[common, gcols], protein.loc[common, gcols]
    n = len(common)
    n_pcs = min(getattr(ra, "N_PCS", 20), n - 1, wsi.shape[1])
    pcs = PCA(n_components=n_pcs, svd_solver="full").fit_transform(wsi.loc[common].values)

    axes = blc.case_batch_labels(cohort, paths, common)
    lab, _ = axes["operator"]
    raw = np.array([lab.get(c, "") or "_missing" for c in common])
    counts = Counter(raw.tolist())
    D = operator_design(raw, counts)
    pcs_corr = residualise(pcs, D)
    print(f"[hoc] {cohort}: n={n}, {D.shape[1]} operator dummies (cap {CAP}, min {MIN_GROUP})", flush=True)

    site = pd.read_csv(site_csv).drop_duplicates("gene")
    corrected_sig = set(site.loc[(site.fdr_batchresid < 0.05) & (site.incremental_r2_batchresid > 0),
                                 "gene"].astype(str))
    genes = [g for g in corrected_sig if g in rna.columns and g in protein.columns]
    cap_genes = int(os.environ.get("MORPHO_MAX_GENES", "0") or 0)
    if cap_genes:
        genes = genes[:cap_genes]
        print(f"[hoc] *** SMOKE TEST: {cap_genes} genes -- NOT the real run ***", flush=True)
    print(f"[hoc] {len(genes)} batch-corrected significant genes (published count check: "
          f"{len(corrected_sig)})", flush=True)
    CV = getattr(ra, "CV_FOLDS", 5)

    rows = []
    for g in genes:
        y = np.asarray(protein[g].values, float)
        v = ~np.isnan(y)
        if v.sum() < MIN_PATIENTS:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1)
        Wv = pcs_corr[v]
        row = {"gene": g, "n": int(v.sum())}
        for s in range(N_SEEDS):
            r0 = cv_r2_seed(xr, yv, CV, s)
            row[f"incr_s{s}"] = cv_r2_seed(np.hstack([xr, Wv]), yv, CV, s) - r0
        rows.append(row)

    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)
    if res.empty:
        print(f"[hoc] nothing to summarise -> {out}")
        return 0

    # ---- validation: seed 0 must reproduce incremental_r2_batchresid exactly ----
    m = res.merge(site[["gene", "incremental_r2_batchresid"]], on="gene")
    dev = (m["incr_s0"] - m["incremental_r2_batchresid"]).abs()
    print(f"\n[hoc] VALIDATION: seed-0 vs published incremental_r2_batchresid, "
          f"{len(m)} genes, max|diff| = {dev.max():.2e} "
          f"{'OK' if dev.max() < 1e-4 else '*** MISMATCH ***'}", flush=True)

    # ---- the verified aggregation: pool all (gene, seed 1-9) pairs, one grand median ----
    held_cols = [f"incr_s{s}" for s in range(1, N_SEEDS)]
    pooled = res[held_cols].values.ravel()
    grand_median = float(np.nanmedian(pooled))
    per_seed_median = res[held_cols].median()
    frac_pos = float((pooled > 0).mean())

    print("\n==================== HELD-OUT EFFECT SIZE, BATCH-CORRECTED SET ====================")
    print(f" cohort {cohort}: {len(res)} genes (of {len(corrected_sig)} batch-corrected significant)")
    print(f" seed-0 (= published corrected median, validation) : {res['incr_s0'].median():+.4f}")
    print(f" per-seed medians, seeds 1-9                        : "
          + ", ".join(f"{v:+.4f}" for v in per_seed_median))
    print(f" HELD-OUT MEDIAN (pooled gene x seed, seeds 1-9)    : {grand_median:+.4f}")
    print(f" fraction of (gene, held-out-seed) pairs positive   : {100*frac_pos:.1f}%")
    print(f" -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
