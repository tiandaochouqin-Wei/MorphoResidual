#!/usr/bin/env python3
"""
Acquisition-robustness pack for C1 (morphology explains the beyond-mRNA protein
residual).

The reviewer's fatal objection: the random-projection capacity control shows the
WSI gain is not just extra degrees of freedom, but it cannot show the embedding
reads morphology rather than a scanner / centre fingerprint; and cross-organ
replication does not help, because a multi-centre artefact could appear in every
organ.

Two facts established before this script, which shape what it must do:

  1. A census of all 2007 SVS headers we hold (build_acquisition_meta.py) found
     1956/1996 slides = 98.0% on a SINGLE Aperio scanner (SS1553) across all five
     cancers; GBM is 247/247 on it. Where a second scanner does exist in ccRCC,
     batch_leak_check.py measured eta^2 = 0.038 with 0/20 PCs significant. So the
     scanner confound is refuted -- by absence of variance in LUAD/UCEC/GBM and
     empirically in ccRCC. PDAC is the exception (432 vs 37 slides, eta^2 = 0.66,
     2/20 PCs) and is the one cohort where scanner must actually be corrected.

  2. The acquisition SESSION does structure the embedding everywhere: the
     scanning-operator axis reaches eta^2 = 0.51-0.76 with 2-7/20 PCs significant
     in all five cohorts, and scan quarter is significant too. Operator and
     quarter are plausible proxies for the tissue-source site GDC does not expose
     (tissue_source_site is _missing for all 1866 CPTAC-3 cases), since slides
     shipped together are scanned together.

So the pack corrects the axis that varies (acquisition batch), not the one that
does not (scanner).

Estimators reported per protein -- kept SEPARATE rather than stacked, because
the reviewer asked for four independent checks and stacking them produces an
estimator so stringent it could kill the signal for statistical rather than
scientific reasons:

  incremental_r2_plain      reference; reproduces residual_analysis.py
  incremental_r2_batchresid A. WSI PCs projected off the batch design, plain CV
  incremental_r2_groupcv    B. raw WSI PCs, GroupKFold BY BATCH (test folds hold
                               acquisition batches never seen in training -- a
                               fingerprint cannot generalise to an unseen batch)
  incremental_r2_both       A+B together, the strictest reading
  incremental_r2_batchonly  D. batch dummies in place of the WSI, PLAIN CV.
                               Plain CV is essential here: under GroupKFold the
                               held-out batch's dummies are all-zero in test, so
                               the control could not use batch information at all
                               and would look falsely clean.

C. Permutation nulls are WITHIN-BATCH: the case<->WSI correspondence is shuffled
only inside a batch, destroying morphology-protein pairing while preserving batch
structure exactly. Nulls are computed for batchresid and for both.

Numerical guard: cross-validated R^2 is unbounded below, and at n~100 with a
small valid subset a single degenerate fold can produce absurd values (a smoke
run produced incremental R^2 = +35.7 for one gene with n=44). R^2 is therefore
clipped to [-5, 1]. Clipping is monotone and applied identically to observed and
permuted statistics, so the permutation p-values stay valid.

Pass criterion, fixed in advance (PAPER_PLAN): the significant count may fall a
lot, but C1 must stay > 0 and the C2 translation/secretion enrichment must
survive; otherwise the claim is downgraded honestly.

Usage:  residual_analysis_sitepack.py <cancer>
Env:    MORPHO_BATCH_AXIS  operator|quarter|cohort_stream|scanner|plex (default operator)
        MORPHO_N_PERM      default 1000   (must be >=1000: with N perms the
                                           smallest p is 1/(N+1), and BH-FDR over
                                           ~10k genes cannot reach 0.05 unless
                                           that floor is well below it)
        MORPHO_MAX_LEVELS  cap on pooled batch levels (default 12)
        MORPHO_N_WORKERS   default 8
Output: <cohort>/results/sitepack_<axis>.csv
"""
import multiprocessing as mp
import os
import sys
import warnings
from collections import Counter

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

SCRIPTS = "/public/home/fjhui/ZW/scripts"
sys.path.insert(0, SCRIPTS)

AXIS = os.environ.get("MORPHO_BATCH_AXIS", "operator")
N_PERM = int(os.environ.get("MORPHO_N_PERM", "1000"))
MAX_LEVELS = int(os.environ.get("MORPHO_MAX_LEVELS", "12"))
N_WORKERS = int(os.environ.get("MORPHO_N_WORKERS", "8"))
N_PCS = 20
CV_FOLDS = 5
RIDGE_ALPHA = 1.0
FDR = 0.05
MIN_VALID = 30
MIN_TRAIN_FOLD = 15
MIN_TEST_FOLD = 3
MIN_GROUP = 3          # a batch must reach this size to earn its own dummy column
MIN_PERMUTABLE = 0.50  # below this, the within-batch null has too little freedom
R2_FLOOR = -5.0
RNG = np.random.RandomState(0)


def _fit(X, y, folds):
    preds = np.zeros_like(y, dtype=float)
    ridge_I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in folds:
        X_tr, X_te, y_tr = X[tr], X[te], y[tr]
        mu = X_tr.mean(axis=0)
        sd = X_tr.std(axis=0)
        sd[sd == 0] = 1.0
        Xtr = (X_tr - mu) / sd
        Xte = (X_te - mu) / sd
        ym = y_tr.mean()
        w = np.linalg.solve(Xtr.T @ Xtr + ridge_I, Xtr.T @ (y_tr - ym))
        preds[te] = Xte @ w + ym
    ss_res = np.sum((y - preds) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    # clip: monotone, applied to observed and null alike, so permutation p-values
    # remain valid while one degenerate fold cannot dominate the whole statistic
    return float(np.clip(r2, R2_FLOOR, 1.0))


def plain_folds(n):
    from sklearn.model_selection import KFold
    return list(KFold(n_splits=CV_FOLDS, shuffle=True,
                      random_state=0).split(np.zeros((n, 1))))


def group_folds(groups):
    """GroupKFold by acquisition batch, with a usability check.

    Falls back to plain CV when the grouping would leave a fold too small to fit
    or to score -- a silently degenerate fold is what produced the absurd R^2
    values in the smoke run.
    """
    from sklearn.model_selection import GroupKFold
    ngroups = len(set(groups))
    if ngroups < 2:
        return plain_folds(len(groups)), False
    k = min(CV_FOLDS, ngroups)
    folds = list(GroupKFold(n_splits=k).split(np.zeros((len(groups), 1)), groups=groups))
    for tr, te in folds:
        if len(tr) < MIN_TRAIN_FOLD or len(te) < MIN_TEST_FOLD:
            return plain_folds(len(groups)), False
    return folds, True


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


_W = {}


def _init(rna, prot, pcs, pcs_resid, dummies, codes, perms, noise_draws):
    _W.update(rna=rna, prot=prot, pcs=pcs, pcs_resid=pcs_resid,
              dummies=dummies, codes=codes, perms=perms, noise=noise_draws)


def _gene(g):
    y = np.asarray(_W["prot"][g])
    valid = ~np.isnan(y)
    nv = int(valid.sum())
    if nv < MIN_VALID:
        return None
    yv = y[valid]
    x_rna = np.asarray(_W["rna"][g])[valid].reshape(-1, 1)

    pf = plain_folds(nv)
    gf, grouped = group_folds(_W["codes"][valid])

    pcs = _W["pcs"][valid]
    pcs_r = _W["pcs_resid"][valid]
    dum = _W["dummies"][valid]

    r2_rna_p = _fit(x_rna, yv, pf)
    r2_rna_g = _fit(x_rna, yv, gf)

    r2_plain = _fit(np.hstack([x_rna, pcs]), yv, pf)
    r2_bres = _fit(np.hstack([x_rna, pcs_r]), yv, pf)
    r2_gcv = _fit(np.hstack([x_rna, pcs]), yv, gf)
    r2_both = _fit(np.hstack([x_rna, pcs_r]), yv, gf)
    r2_bonly = _fit(np.hstack([x_rna, dum]), yv, pf) if dum.shape[1] else np.nan
    # Capacity-matched partner for the batch-only control: the SAME number of
    # columns, but pure noise. Comparing batch dummies against the 20-PC WSI
    # block is not a fair contest -- a handful of dummies pays a far smaller
    # cross-validation penalty, which is why the batch-only control could read
    # 22.9% positive against the real signal's 16.0% on PDAC's scanner axis and
    # look alarming for no good reason. Batch carries real signal only insofar
    # as it beats noise of its own width.
    if dum.shape[1]:
        r2_nm = np.mean([_fit(np.hstack([x_rna, nd[valid]]), yv, pf)
                         for nd in _W["noise"]])
    else:
        r2_nm = np.nan

    incr_plain = r2_plain - r2_rna_p
    incr_bres = r2_bres - r2_rna_p
    incr_gcv = r2_gcv - r2_rna_g
    incr_both = r2_both - r2_rna_g
    incr_bonly = (r2_bonly - r2_rna_p) if dum.shape[1] else np.nan
    incr_nm = (r2_nm - r2_rna_p) if dum.shape[1] else np.nan

    null_bres = np.empty(len(_W["perms"]))
    null_both = np.empty(len(_W["perms"]))
    for i, perm in enumerate(_W["perms"]):
        xr = _W["pcs_resid"][perm][valid]
        null_bres[i] = _fit(np.hstack([x_rna, xr]), yv, pf) - r2_rna_p
        null_both[i] = _fit(np.hstack([x_rna, xr]), yv, gf) - r2_rna_g
    n_p = len(_W["perms"])
    p_bres = (np.sum(null_bres >= incr_bres) + 1) / (n_p + 1)
    p_both = (np.sum(null_both >= incr_both) + 1) / (n_p + 1)

    return dict(gene=g, n=nv, grouped_cv=bool(grouped),
                r2_rna=r2_rna_p, r2_rna_groupcv=r2_rna_g,
                incremental_r2_plain=incr_plain,
                incremental_r2_batchresid=incr_bres,
                incremental_r2_groupcv=incr_gcv,
                incremental_r2_both=incr_both,
                incremental_r2_batchonly=incr_bonly,
                incremental_r2_noisematched=incr_nm,
                pval_batchresid=p_bres, pval_both=p_both)


def main():
    cancer = sys.argv[1]
    if N_PERM < 1000:
        print(f"WARNING: N_PERM={N_PERM} < 1000. The smallest attainable p is "
              f"{1 / (N_PERM + 1):.4f}; after BH over ~10k genes nothing can reach "
              f"FDR<{FDR}, so a zero count would be an artefact of this setting, "
              f"not a result.", flush=True)

    import batch_leak_check as blc
    paths = blc.cohort_paths(cancer)
    os.environ.update(
        MORPHO_ROOT=paths["root"], MORPHO_RNA_MANIFEST=paths["rna_manifest"],
        MORPHO_ALIQUOT_XWALK=paths["xwalk"], MORPHO_SLIDE_MAP=paths["slide_map"],
        MORPHO_PROTEIN_TSV=paths["protein"], MORPHO_WSI_EMB_DIR=paths["emb"],
        MORPHO_OUT=paths["out"])
    import residual_analysis as ra

    print(f"=== sitepack: cancer={cancer} axis={AXIS} n_perm={N_PERM} ===", flush=True)
    rna = ra.load_rna_matrix()
    protein = ra.load_protein_matrix()
    wsi = ra.load_wsi_embeddings()

    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    genes = sorted(set(rna.columns) & set(protein.columns))
    print(f"3-way patients: {len(common)}; common genes: {len(genes)}", flush=True)
    if len(common) < 20:
        sys.exit(f"FAIL: only {len(common)} patients")
    rna = rna.loc[common, genes]
    protein = protein.loc[common, genes]

    from sklearn.decomposition import PCA
    n_pcs = min(N_PCS, len(common) - 1, wsi.shape[1])
    # svd_solver="full" is REQUIRED, not cosmetic. On a 99x1024 matrix sklearn's
    # auto solver picks randomized SVD with random_state=None, so residual_analysis.py
    # is not reproducible run to run: two identical PCA(n_components=20) calls on the
    # same data were measured to differ by up to 24.4 in the scores, with different
    # explained_variance_ratio_. Every count this project has recorded inherits that
    # noise. "full" makes the decomposition exact and deterministic.
    pcs = PCA(n_components=n_pcs, svd_solver="full").fit_transform(wsi.loc[common].values)
    print(f"WSI PCs: {pcs.shape} (svd_solver=full, deterministic)", flush=True)

    axes = blc.case_batch_labels(cancer, paths, common)
    if AXIS not in axes:
        sys.exit(f"axis {AXIS!r} not available; have {sorted(axes)}")
    lab, _purity = axes[AXIS]
    raw = np.array([lab.get(c, "") or "_missing" for c in common])
    counts = Counter(raw.tolist())
    levels = sorted(counts)
    if len(levels) < 2:
        sys.exit(f"axis {AXIS} is degenerate for {cancer} ({levels}) -- no batch "
                 f"variation exists to correct for, which is itself the answer")

    # Grouping for CV and for the permutation null uses the RAW batch and never
    # merges anything. An earlier version pooled every batch outside the top
    # MAX_LEVELS into one "_other" level, which held 45% of PDAC's patients
    # across 54 distinct operators. That single bucket corrupted all three
    # checks at once: residualisation subtracted one shared mean that described
    # none of its members, the within-batch permutation was free to swap cases
    # between genuinely different operators, and GroupKFold treated the whole
    # heterogeneous mass as one group that a test fold could never generalise to
    # -- which is why PDAC's strictest estimator collapsed to zero.
    gidx = {lv: i for i, lv in enumerate(levels)}
    codes = np.array([gidx[v] for v in raw])

    # Dummies are a separate question from grouping: their count must stay small
    # relative to n (23 plex dummies at n=103 drove the baseline R2 negative in
    # residual_analysis_adjusted.py). So dummy-code only the largest batches;
    # cases in the remaining small batches get an all-zero row and are absorbed
    # by the intercept. They then receive NO batch correction, which is honest,
    # instead of a shared correction that is simply wrong for them.
    dummy_levels = [lv for lv, k in counts.most_common(MAX_LEVELS) if k >= MIN_GROUP]
    if len(dummy_levels) == len(levels) and dummy_levels:
        dummy_levels = dummy_levels[:-1]      # keep the design full-rank
    D = np.zeros((len(common), len(dummy_levels)))
    for j, lv in enumerate(dummy_levels):
        D[:, j] = (raw == lv).astype(float)
    B = np.hstack([np.ones((len(common), 1)), D])

    n_undummied = int(sum(counts[lv] for lv in levels if lv not in set(dummy_levels)))
    print(f"axis {AXIS}: {len(levels)} raw batches, kept as {len(levels)} CV/permutation "
          f"groups (no pooling)", flush=True)
    print(f"  dummy-coded batches: {len(dummy_levels)} (cap {MAX_LEVELS}, min size "
          f"{MIN_GROUP}); {n_undummied}/{len(common)} cases left to the intercept",
          flush=True)
    print(f"  batch sizes: {dict(counts.most_common(12))}"
          f"{' ...' if len(levels) > 12 else ''}", flush=True)

    coef, *_ = np.linalg.lstsq(B, pcs, rcond=None)
    pcs_resid = pcs - B @ coef
    kept = pcs_resid.var(axis=0).sum() / max(pcs.var(axis=0).sum(), 1e-12)
    print(f"batch residualisation retains {kept:.1%} of total PC variance "
          f"({D.shape[1]} dummies at n={len(common)})", flush=True)

    n = len(common)
    members = {lv: np.where(codes == gidx[lv])[0] for lv in levels}
    perms = []
    for _ in range(N_PERM):
        p = np.arange(n)
        for _lv, m in members.items():
            if len(m) > 1:
                p[m] = RNG.permutation(m)
        perms.append(p)
    fixed = sum(len(m) for m in members.values() if len(m) <= 1)
    permutable = 1.0 - fixed / n
    print(f"{N_PERM} within-batch permutations built; {fixed}/{n} cases sit in "
          f"singleton batches and are necessarily fixed "
          f"({permutable:.0%} of cases permutable)", flush=True)
    if permutable < MIN_PERMUTABLE:
        # Not merging them into a fake batch to buy freedom -- that is exactly the
        # bug this version removes. A fragmented axis simply cannot support this
        # test, and a p-value computed from a nearly-frozen null would be
        # meaningless rather than conservative.
        print(f"  WARNING: only {permutable:.0%} of cases can move under the "
              f"within-batch null (threshold {MIN_PERMUTABLE:.0%}). The permutation "
              f"test on axis '{AXIS}' is under-powered for {cancer}; read the "
              f"point estimates and use a coarser axis (e.g. quarter) for "
              f"significance.", flush=True)

    rna_d = {g: rna[g].values for g in genes}
    prot_d = {g: protein[g].values for g in genes}
    # noise blocks matched to the dummy count, averaged over draws so a single
    # lucky draw cannot set the bar
    noise_draws = [RNG.standard_normal((n, D.shape[1])) for _ in range(3)] \
        if D.shape[1] else []

    results = []
    with mp.Pool(N_WORKERS, initializer=_init,
                 initargs=(rna_d, prot_d, pcs, pcs_resid, D, codes, perms,
                           noise_draws)) as pool:
        for i, r in enumerate(pool.imap_unordered(_gene, genes, chunksize=8), 1):
            if r is not None:
                results.append(r)
            if i % 1000 == 0:
                print(f"  {i}/{len(genes)} genes...", flush=True)

    res = pd.DataFrame(results)
    res["fdr_batchresid"] = bh_fdr(res["pval_batchresid"].values)
    res["fdr_both"] = bh_fdr(res["pval_both"].values)
    res = res.sort_values("incremental_r2_batchresid", ascending=False)
    res.insert(0, "cohort", cancer)
    res.insert(1, "batch_axis", AXIS)
    os.makedirs(paths["out"], exist_ok=True)
    out = os.path.join(paths["out"], f"sitepack_{AXIS}.csv")
    res.to_csv(out, index=False)

    n_bres = int(((res["fdr_batchresid"] < FDR) & (res["incremental_r2_batchresid"] > 0)).sum())
    n_both = int(((res["fdr_both"] < FDR) & (res["incremental_r2_both"] > 0)).sum())
    ngrouped = int(res["grouped_cv"].sum())

    print(f"\n=== {cancer} / axis={AXIS} : {len(res)} genes tested ===")
    print(f"  A. batch-residualised, plain CV : {n_bres} significant (FDR<{FDR})")
    print(f"  A+B. residualised + GroupKFold  : {n_both} significant (FDR<{FDR})")
    print(f"  genes where GroupKFold was usable: {ngrouped}/{len(res)}")
    print()
    for col in ["incremental_r2_plain", "incremental_r2_batchresid",
                "incremental_r2_groupcv", "incremental_r2_both",
                "incremental_r2_batchonly", "incremental_r2_noisematched"]:
        v = res[col].dropna()
        if len(v):
            print(f"  {col:<32} median {v.median():+.4f}  frac>0 {np.mean(v > 0):6.1%}")
    bo = res["incremental_r2_batchonly"].dropna()
    nm = res["incremental_r2_noisematched"].dropna()
    if len(bo) and len(nm):
        print(f"\n  D. negative control, read against its capacity-matched partner:")
        print(f"     batch dummies ({D.shape[1]} cols) frac>0 = {np.mean(bo > 0):.1%}")
        print(f"     noise, same {D.shape[1]} cols        frac>0 = {np.mean(nm > 0):.1%}")
        print(f"     excess attributable to batch STRUCTURE rather than width: "
              f"{np.mean(bo > 0) - np.mean(nm > 0):+.1%}")
        print(f"     (comparing batch dummies directly against the 20-PC WSI block "
              f"is not a fair\n      contest -- fewer columns, smaller CV penalty. "
              f"This is the fair one.)")
    print(f"\n  top 12 by batch-residualised incremental R2:")
    cols = ["gene", "n", "incremental_r2_plain", "incremental_r2_batchresid",
            "incremental_r2_both", "incremental_r2_batchonly",
            "pval_batchresid", "fdr_batchresid"]
    print(res.head(12)[cols].to_string(index=False))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
