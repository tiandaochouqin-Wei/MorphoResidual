#!/usr/bin/env python3
"""
theta_kernel.py -- shared numerics for the own-mRNA measurement-error statistics
(MorphoResidual, Limitations (v); theory/DERIVATION.md sections 2.1, 2.3, 2.4, 4.1).
Pure numpy; no project data, no I/O.  Imported by local_theta.py (the local analysis) and,
for the OPTIONAL on-cluster permutation path, by measurement_error_checks.py.

The estimator is the published one (residual_analysis.py = "RA"):
    5-fold KFold(shuffle=True, random_state=0) on the gene's non-missing-protein patients,
    ridge alpha = 1 on columns standardised on the TRAINING fold (ddof 0, zero sd -> 1),
    response centred on the training fold, closed form (X'X + alpha I) w = X'(y - mean).
`oof_ridge` is RA.cv_r2's body returning the out-of-fold predictions; it is the reference
implementation.  `wonly_oof_weighted` is a batched/weighted version of the same ridge used for
the permutation and bootstrap loops; it is asserted equal to `oof_ridge` (1e-9) wherever both
are available (local_theta.py runs that assertion on every gene before it does anything else).

Statistics (per gene; M, P standardised over the gene's fit patients so beta = delta = rho):
    r_p = P - rho M ,  r_m = M - rho P                       (in-sample OLS residuals)
    rp_hat = W-ONLY out-of-fold ridge prediction of r_p      (NOT the nested increment d_p)
    N   = rho * cov(rp_hat, M)         own-transcript alignment numerator
    Nc  = rho * corr(rp_hat, M)        scale-free form (the one that is TESTED)
    D   = cov(rp_hat, r_p)             denominator, ~ permutation-centred W-only increment
    theta = N / D ,  lambda_hat_noise = theta / (1 + theta)    (set level: sum N / sum D)
    S_cross = sign(rho) * corr(rp_hat, r_m)                    (F18 sign test, cross form)
Secondary (optional):  S_cross_m = sign(rho) corr(rm_hat, r_p),  S_pp = sign(rho) corr(rp_hat, rm_hat).
All covariances/correlations use ddof 0 (means of centred products), exactly as in
theory/sim_signtest.py, which this module is cross-checked against in test_synthetic.py.
"""
import numpy as np
from sklearn.model_selection import KFold

CV_FOLDS = 5
RIDGE_ALPHA = 1.0
STAT_KEYS = ["N", "Nc", "D", "S_cross"]
SEC_KEYS = ["S_cross_m", "S_pp"]

# --------------------------------------------------------------------------------------
# folds (identical to RA.cv_r2: KFold(5, shuffle, random_state=0) on the gene's nv patients)
# --------------------------------------------------------------------------------------
_FOLDS = {}


def folds_for(nv):
    f = _FOLDS.get(nv)
    if f is None:
        kf = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=0)
        f = [(tr, te) for tr, te in kf.split(np.arange(nv))]
        _FOLDS[nv] = f
    return f


def fold_ids_for(nv):
    fid = np.empty(nv, dtype=np.int64)
    for k, (_, te) in enumerate(folds_for(nv)):
        fid[te] = k
    return fid


# --------------------------------------------------------------------------------------
# reference ridge (line-for-line RA.cv_r2, returning predictions)
# --------------------------------------------------------------------------------------
def oof_ridge(X, y, folds, alpha=RIDGE_ALPHA):
    preds = np.zeros_like(y, dtype=float)
    ridge_I = alpha * np.eye(X.shape[1])
    for tr, te in folds:
        X_tr, X_te, y_tr = X[tr], X[te], y[tr]
        x_mean = X_tr.mean(axis=0)
        x_std = X_tr.std(axis=0)
        x_std[x_std == 0] = 1.0
        Xtr_s = (X_tr - x_mean) / x_std
        Xte_s = (X_te - x_mean) / x_std
        y_mean = y_tr.mean()
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + ridge_I, Xtr_s.T @ (y_tr - y_mean))
        preds[te] = Xte_s @ w + y_mean
    return preds


def r2_of(y, p):
    ss_res = np.sum((y - p) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def published_increment(M, P, pcs, v=None):
    """Published statistic for one gene: (n_fit, r2_rna, r2_both, incremental_r2) with the gene's
    fit patients v = ~isnan(P).  Same arithmetic as RA._process_gene (hstack keeps pcs' dtype
    rules: float64 mRNA column upcasts a float32 PC block exactly)."""
    if v is None:
        v = ~np.isnan(P)
    nv = int(v.sum())
    y = P[v]
    xa = M[v].reshape(-1, 1)
    folds = folds_for(nv)
    r2_rna = r2_of(y, oof_ridge(xa, y, folds))
    r2_both = r2_of(y, oof_ridge(np.hstack([xa, pcs[v]]), y, folds))
    return nv, r2_rna, r2_both, r2_both - r2_rna


def published_null_increments(M, P, pcs, perms, v=None):
    """RA's null for one gene: incr under wsi_pcs[perm][valid] for each permutation in perms."""
    if v is None:
        v = ~np.isnan(P)
    nv = int(v.sum())
    y = P[v]
    xa = M[v].reshape(-1, 1)
    folds = folds_for(nv)
    r2_rna = r2_of(y, oof_ridge(xa, y, folds))
    out = np.empty(len(perms))
    for i, p in enumerate(perms):
        out[i] = r2_of(y, oof_ridge(np.hstack([xa, pcs[p][v]]), y, folds)) - r2_rna
    return out


# --------------------------------------------------------------------------------------
# literal (reference) per-gene statistics: unbatched, uses oof_ridge
# --------------------------------------------------------------------------------------
def _cov(a, b):
    return np.mean((a - a.mean()) * (b - b.mean()))


def _corr(a, b):
    with np.errstate(all="ignore"):
        return _cov(a, b) / np.sqrt(_cov(a, a) * _cov(b, b))


def prep_gene(M, P):
    """Standardise a gene's (M, P) over its fit patients and form the OLS residuals.
    Returns None for a degenerate gene (zero variance)."""
    sM, sP = M.std(), P.std()
    if not (sM > 0 and sP > 0):
        return None
    Mz = (M - M.mean()) / sM
    Pz = (P - P.mean()) / sP
    rho = float(np.mean(Mz * Pz))
    return dict(Mz=Mz, Pz=Pz, rho=rho, rp=Pz - rho * Mz, rm=Mz - rho * Pz)


def gene_stats_literal(W, M, P, secondary=False, folds=None):
    """Reference statistics for ONE gene.  W (nv,K) is the morphology block as seen by this
    gene's fit patients (observed or permuted); M, P its raw mRNA / protein on those patients."""
    g = prep_gene(M, P)
    if g is None:
        return None
    W = np.asarray(W, dtype=np.float64)     # exact upcast of a float32 PC block (as RA's hstack does)
    folds = folds if folds is not None else folds_for(len(M))
    rho, rp, rm, Mz = g["rho"], g["rp"], g["rm"], g["Mz"]
    rph = oof_ridge(W, rp, folds)
    sg = np.sign(rho)
    out = dict(rho=rho,
               N=rho * _cov(rph, Mz), Nc=rho * _corr(rph, Mz),
               D=_cov(rph, rp), S_cross=sg * _corr(rph, rm))
    if secondary:
        rmh = oof_ridge(W, rm, folds)
        out["S_cross_m"] = sg * _corr(rmh, rp)
        out["S_pp"] = sg * _corr(rph, rmh)
    return out


# --------------------------------------------------------------------------------------
# batched / weighted version (permutations: weights = 1, W varies; bootstrap: W fixed, weights vary)
# --------------------------------------------------------------------------------------
def wonly_oof_weighted(Wst, Y, fold_id, wts, alpha=RIDGE_ALPHA):
    """Out-of-fold ridge prediction of Y from W only, batched over b.

    Wst    (Bw, nv, K)   Bw in {1, B}
    Y      (B, T, nv)    T targets per draw (r_p [, r_m])
    fold_id (nv,)        fold of each patient (from folds_for(nv))
    wts    (Bt, nv)      Bt in {1, B}; integer multiplicities (bootstrap) or ones (permutation).
                         A weight-c row is exactly c duplicated rows (ddof-0 weighted moments).
    Returns (B, T, nv); entries are valid for every row (rows of weight 0 are simply ignored
    downstream)."""
    Wst = np.asarray(Wst, dtype=np.float64)
    B, T, nv = Y.shape
    K = Wst.shape[2]
    pred = np.empty((B, T, nv))
    I = alpha * np.eye(K)
    wts = np.asarray(wts, float)
    for k in range(CV_FOLDS):
        te = fold_id == k
        wt = np.broadcast_to(wts * (~te)[None, :], (B, nv))         # (B, nv) train weights
        Ntr = wt.sum(1)                                             # (B,)
        mu = (wt[:, None, :] @ Wst)[:, 0, :] / Ntr[:, None]         # (B, K)
        Wc = Wst - mu[:, None, :]                                   # (B, nv, K)
        var = (wt[:, None, :] @ (Wc * Wc))[:, 0, :] / Ntr[:, None]
        sd = np.sqrt(var)
        sd = np.where(sd == 0, 1.0, sd)
        As = Wc / sd[:, None, :]                                    # (B, nv, K)
        ym = (Y * wt[:, None, :]).sum(-1) / Ntr[:, None]            # (B, T)
        Yc = Y - ym[..., None]
        Aw = As * wt[..., None]
        AwT = np.swapaxes(Aw, 1, 2)
        Gm = AwT @ As + I                                           # (B, K, K)
        rhs = AwT @ np.swapaxes(Yc, 1, 2)                           # (B, K, T)
        w = np.linalg.solve(Gm, rhs)                                # (B, K, T)
        pk = As[:, te, :] @ w                                       # (B, nte, T)
        pred[:, :, te] = np.swapaxes(pk, 1, 2) + ym[..., None]
    return pred


def _wmean(a, w, N):
    return (a * w).sum(-1) / N


def _wcov(a, b, w, N):
    am = _wmean(a, w, N)[..., None]
    bm = _wmean(b, w, N)[..., None]
    return (w * (a - am) * (b - bm)).sum(-1) / N


def _wcorr(a, b, w, N):
    with np.errstate(all="ignore"):
        return _wcov(a, b, w, N) / np.sqrt(_wcov(a, a, w, N) * _wcov(b, b, w, N))


def gene_stats_batch(Wst, M, P, fold_id, wts, secondary=False):
    """Statistics for ONE gene over a batch of B draws.
    Wst (Bw, nv, K); M, P raw (nv,); wts (Bt, nv).  M, P are (weighted-)standardised PER DRAW so a
    bootstrap draw is a complete re-analysis of the resampled patients; for permutation draws
    (wts = 1) the standardisation is identical across draws.  Returns dict of (B,) arrays or None
    for a degenerate gene."""
    B = max(Wst.shape[0], wts.shape[0])
    wts = np.asarray(wts, float)
    Nn = wts.sum(1)                                     # (Bt,)
    mM = (wts * M).sum(1) / Nn
    mP = (wts * P).sum(1) / Nn
    sM = np.sqrt((wts * (M[None] - mM[:, None]) ** 2).sum(1) / Nn)
    sP = np.sqrt((wts * (P[None] - mP[:, None]) ** 2).sum(1) / Nn)
    if not (np.all(sM > 0) and np.all(sP > 0)):
        return None
    Mz = (M[None] - mM[:, None]) / sM[:, None]          # (Bt, nv)
    Pz = (P[None] - mP[:, None]) / sP[:, None]
    rho = (wts * Mz * Pz).sum(1) / Nn                   # (Bt,)
    rp = Pz - rho[:, None] * Mz
    rm = Mz - rho[:, None] * Pz
    rp, rm, Mz, wts_b = (np.broadcast_to(a, (B, a.shape[1])) for a in (rp, rm, Mz, wts))
    rho_b = np.broadcast_to(rho, (B,))
    Nn_b = np.broadcast_to(Nn, (B,))
    Y = np.stack([rp, rm], axis=1) if secondary else rp[:, None, :]
    pred = wonly_oof_weighted(Wst, np.ascontiguousarray(Y), fold_id, wts, RIDGE_ALPHA)
    rph = pred[:, 0]
    sg = np.sign(rho_b)
    out = dict(rho=np.array(rho_b),
               N=rho_b * _wcov(rph, Mz, wts_b, Nn_b),
               Nc=rho_b * _wcorr(rph, Mz, wts_b, Nn_b),
               D=_wcov(rph, rp, wts_b, Nn_b),
               S_cross=sg * _wcorr(rph, rm, wts_b, Nn_b))
    if secondary:
        rmh = pred[:, 1]
        out["S_cross_m"] = sg * _wcorr(rmh, rp, wts_b, Nn_b)
        out["S_pp"] = sg * _wcorr(rph, rmh, wts_b, Nn_b)
    return out


# --------------------------------------------------------------------------------------
# permutation / bootstrap draws
# --------------------------------------------------------------------------------------
def build_perms(n, B, scheme="ra", seed=0):
    """Unstratified patient<->PC-row permutations.  'ra': RandomState(seed), sequential draws, the
    very permutations of residual_analysis.py (RNG = RandomState(0); first draw first).
    'indexed': RandomState(seed + i) per draw (transcriptome_baseline.py)."""
    if B == 0:
        return np.zeros((0, n), dtype=np.int64)
    if scheme == "ra":
        rng = np.random.RandomState(seed)
        return np.asarray([rng.permutation(n) for _ in range(B)], dtype=np.int64)
    if scheme == "indexed":
        return np.asarray([np.random.RandomState(seed + i).permutation(n) for i in range(B)],
                          dtype=np.int64)
    raise ValueError(scheme)


def labels_to_codes(labels):
    """Raw batch labels -> (levels, codes), the sitepack convention: empty label -> '_missing',
    levels sorted, no pooling of small batches."""
    raw = np.array([(str(x) if str(x) != "" else "_missing") for x in labels])
    levels = sorted(set(raw.tolist()))
    idx = {lv: i for i, lv in enumerate(levels)}
    return levels, np.array([idx[x] for x in raw]), raw


def build_perms_stratified(labels, B, seed=0):
    """Within-stratum permutations, same construction as residual_analysis_sitepack.py:
    p = arange(n); for each level (sorted) with >1 member, p[members] = RNG.permutation(members),
    RNG = RandomState(seed), draws taken sequentially.  Row i of the permuted PC block is
    pcs[p[i]].  With seed 0 and the same labels these are sitepack's permutations."""
    levels, codes, _ = labels_to_codes(labels)
    n = len(codes)
    members = {lv: np.where(codes == i)[0] for i, lv in enumerate(levels)}
    rng = np.random.RandomState(seed)
    perms = []
    for _ in range(B):
        p = np.arange(n)
        for _lv, m in members.items():
            if len(m) > 1:
                p[m] = rng.permutation(m)
        perms.append(p)
    fixed = int(sum(len(m) for m in members.values() if len(m) <= 1))
    return np.asarray(perms, dtype=np.int64).reshape(B, n), dict(
        n_levels=len(levels), n_fixed_singletons=fixed, frac_permutable=1.0 - fixed / n,
        level_sizes=sorted((len(m) for m in members.values()), reverse=True)[:12])


def boot_counts(n, B, seed=0):
    """(B, n) patient multiplicities of B resamples of n patients with replacement."""
    rng = np.random.RandomState(seed)
    idx = rng.randint(0, n, size=(B, n))
    return np.stack([np.bincount(r, minlength=n) for r in idx]).astype(float)


# --------------------------------------------------------------------------------------
# depth-side helpers shared with the export
# --------------------------------------------------------------------------------------
def icc_oneway(groups):
    """One-way random-effects ICC(1) for unbalanced groups; vectorised over genes.
    groups: list of (k_i, G) arrays (k_i >= 2).  Returns (icc, within_sd, n_groups, mean_files)."""
    a = len(groups)
    ns = np.array([g.shape[0] for g in groups], float)
    N = ns.sum()
    means = np.stack([g.mean(0) for g in groups])
    grand = np.concatenate(groups).mean(0)
    msb = (ns[:, None] * (means - grand) ** 2).sum(0) / (a - 1)
    msw = sum(((g - g.mean(0)) ** 2).sum(0) for g in groups) / (N - a)
    k0 = (N - (ns ** 2).sum() / N) / (a - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        icc = (msb - msw) / (msb + (k0 - 1) * msw)
    return icc, np.sqrt(msw), a, float(N / a)


def lambda_poisson(counts_fit, pseudo=0.5):
    """Per-gene Poisson-implied UPPER bound on single-library reliability (DERIVATION 3, item 5):
        lambda_P = 1 - mean_i(1/(c_i+pseudo)) / var_i(ln(c_i+pseudo)),  clipped to [0,1].
    counts_fit (n_fit, G) with NaN outside the gene's fit patients.  The pseudo-count and the
    natural-log scale are CONVENTIONS that the pre-specification must freeze."""
    c = counts_fit + pseudo
    with np.errstate(all="ignore"):
        pois = np.nanmean(1.0 / c, axis=0)
        var = np.nanvar(np.log(c), axis=0, ddof=1)
        lam = 1.0 - pois / var
    return np.clip(lam, 0.0, 1.0)
