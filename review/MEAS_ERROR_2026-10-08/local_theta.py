#!/usr/bin/env python3
"""
local_theta.py -- LOCAL computation of the own-transcript-alignment statistics of
theory/DERIVATION.md (N, N_c, D, theta, lambda_hat_noise, S_cross) from the files exported by
measurement_error_checks.py v2 (inputs_<c>.npz).  Drafted 2026-10-08; revised 2026-10-08 to the FINAL
reading rules of DERIVATION.md section 9 / "Final reading rules".  NOT YET RUN on project data.

STATISTICS (per gene, same folds / standardisation / ridge alpha as the published estimator --
theta_kernel.py, which is asserted equal to RA.cv_r2's ridge):
    r_p = P - rho M ,  r_m = M - rho P  (M, P standardised over the gene's fit patients)
    rp_hat = W-ONLY 5-fold OOF ridge prediction of r_p  (20 morphology PCs; NOT the nested d_p)
    N  = rho cov(rp_hat, M),   D = cov(rp_hat, r_p),   S_cross = sign(rho) corr(rp_hat, r_m)
    N_c = sign(rho) corr(rp_hat, M)      <-- THE TESTED STATISTIC (set level: sum over genes).
          Computed as (kernel Nc = rho corr)/|rho|, for the observed value and for every null draw alike
          (rho depends on (M, P) only, which the null preserves), so theta_kernel.py needs no change.
          The old rho-weighted sum is kept ONLY as the diagnostic column Nc_rhow_* (z of the biased form).
    theta_raw = sum N / sum D ;  theta_pc = (sum N - E0 sum N)/(sum D - E0 sum D)  (E0 = mean over the null
    draws of the arm; THE theta of record), lambda_hat_noise = theta_pc/(1+theta_pc).
    theta_pc CI: patient bootstrap (resample the n patients, shared across genes; each patient keeps the
    fold it has in the gene's published KFold; integer weights == duplicated rows) of
    (sum N_b - E0 sum N)/(sum D_b - E0 sum D) with E0 FIXED at the observed-sample null mean (convention).

ELIGIBILITY AND STRATA (DERIVATION 3 item 1, section 9).  A gene's fit set is the patients with non-missing
protein.  n_fit = n (complete gene) is the PRIMARY eligibility; every base set (all, pinned_sig = published
FDR<0.05 & incr>0, --sets-file families, pinned_sig&<family>) is built on complete genes only.  Strata of a base
set <S> (reported only with >= 15 genes; theta fields are flagged theta_reading_ok only with >= 20 genes):
    <S>|depthT1..3        tertiles of the raw-count median (of the primary set's genes)
    <S>|miss0-5           0 < missing fraction <= 5%      <S>|miss5-20   5% < f <= 20%     <S>|miss20+  f > 20%
    <S>|miss<=5           f <= 5% (complete + 0-5%): the stratum a pinned set with < 20 complete genes is read on
                          (flag mnar_exposed = fewer than 20 complete genes in the base set)
Selection strength s = n_sel / G_tested of the stratum itself (pinned genes in the stratum over published-
tested genes satisfying the same stratum definition).  kappa(s, lambda_lo) = table of DERIVATION 4.4
(lower-inclusive row boundaries); unselected sets (all, plain families) and --confirmatory: kappa = 1.

NULL ARMS (--strata).  Joint patient<->PC-row permutation, ONE draw shared by all genes and all targets, folds /
standardisation / alpha fixed.  none = unrestricted; operator / plex / scanner / operator_scanner = within batch
levels (residual_analysis_sitepack.py construction, levels sorted, no pooling, singletons fixed).  PRIMARY arm =
--primary-arm (default operator; sensitivity = none).  Draws are FRESH: default --seed 20261007 is NOT the
RandomState(0) sequence that produced the published selection (--seed 0 reproduces it).  Nothing is residualised
on operator.

PER SET AND ARM (sets.csv):  z/p of sum N, sum N_c (sign-weighted), sum D; theta_raw, theta_pc (+ bootstrap CI,
SE), lambda_hat; S_cross; diagnostic Wm = corr(W-only OOF prediction of M, M) within the fit sets (set mean, z
against the same null); per-gene proportions outside the 2.5/97.5 % N_c band of their own null; depth, ICC,
lambda_Poisson; rho2bar.  READING NUMBERS (reading.csv; numbers only, NEVER a verdict): z_D, z_N under the primary
and the unrestricted arm, theta_pc with CI, s, kappa, lambda_lo (+source), lambda_hi, theta_noise_min =
kappa lam_lo/(1-lam_lo), margins theta_pc - theta_noise_min and CI_hi - theta_noise_min, the R1 band limit and
|theta_pc| - band, f_X range, post-transcriptional-share lower bound, phi_op (D and nested increment), and the
sensitivities below.  hetero.csv: pairwise stratum differences of theta_pc in units of the larger bootstrap SE.
SENSITIVITIES (primary-arm permutations; variants.csv and columns of reading.csv):
    spline   baseline P~M replaced by a 3-df natural cubic spline fitted IN-SAMPLE on the gene's fit patients (as
             the linear r_p is): r_p^s = P - fhat(M), N^s = cov(rp_hat, fhat), theta^s = N^s/D^s, and the nested
             increment Delta_p^spline (published ridge on [spline basis] vs [spline basis, PCs]; knots at the fit
             patients' min / terciles / max -- a convention) beside the published Delta_p
    rnapc    baseline P ~ M + the 20 transcriptome-wide RNA PCs of the export (in-sample OLS); own-transcript
             component = b_M * M_perp (M residualised on the RNA PCs); N, D defined as in the linear case
    purity   M and P residualised in-sample on --purity-file (patient<TAB>purity) before rho, r_p, rp_hat
    (the linear 'raw theta_set' and the unrestricted arm are the other two sensitivities)

FAIL-CLOSED (same idea as server_export/scripts/c1_ucec_posthoc_sets.py).  On any inputs file that
does not carry is_synthetic=True the statistics run ONLY if the sidecar named in SIDECAR exists and
every file listed in it (the frozen pre-specification, this script, theta_kernel.py, every inputs
file used, and any --sets-file / --purity-file given) still has its listed sha256 -- and this script,
theta_kernel.py, each inputs file and those side files must themselves be listed.  Exceptions that never
compute a discriminator:
    --repro-only    published quantities only: re-derives the published per-gene incremental R^2
                    from the exported inputs (optionally --repro-pvals K: the published
                    permutation p-values of K genes) and checks them to 1e-10.
    --synthetic     allowed only for files flagged is_synthetic=True (test data).
Before any statistic the full run ALSO re-derives the published increments for all genes and
refuses if they differ from the published results by > 1e-10, and asserts batched == reference
numerics gene by gene.

Usage:
   python local_theta.py --inputs export/ucec/inputs_ucec.npz --repro-only
   python local_theta.py --inputs export/ucec/inputs_ucec.npz --strata operator,none --B 200 --boot 200
   python local_theta.py --inputs syn.npz --synthetic --B 100        (synthetic test data only)
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import theta_kernel as tk  # noqa: E402

PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))          # MorphoResidual_paper/
PRESPEC_REL = "review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md"
SIDECAR = os.path.join(PAPER, "review", "POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md.sha256")
REPRO_TOL = 1e-10
BATCH_TOL = 1e-9
MIN_SET = 5          # minimum set size for any z
MIN_STRATUM = 15     # a stratum is reported only with >= 15 genes
MIN_THETA = 20       # theta (and its CI) is flagged readable only with >= 20 genes
SEED_DEFAULT = 20261007
BOOT_SEED_DEFAULT = 20261008
# constants of the final reading rules (DERIVATION section 9); printed beside the numbers, never applied as a verdict
Z_GATE_D = 3.0
Z_ANTI = -3.0
Z_ALIGN = 2.0
Z_STRONG = 3.0
R1_BAND_LARGE = 0.05      # n_genes >= 100
R1_BAND_SMALL = 0.10      # 20 <= n_genes < 100
R4_CUT = -0.3
# kappa(s, lambda_lo): rows = lower bound of the strength band (lower-inclusive), columns = lambda_lo <=0.55 / <=0.65 / >0.65
KAPPA_TABLE = [(0.50, (0.65, 0.65, 0.55)),
               (0.20, (0.55, 0.45, 0.35)),
               (0.10, (0.45, 0.40, 0.25)),
               (0.05, (0.40, 0.35, 0.25)),
               (0.02, (0.35, 0.35, 0.25)),
               (0.00, (0.30, 0.30, 0.25))]
VARIANTS_ALL = ("spline", "rnapc", "purity")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def kappa_of(s, lam_lo):
    """Winner's-curse floor kappa(s, lambda_lo) of DERIVATION 4.4 (lower-inclusive strength bands)."""
    if not (np.isfinite(s) and np.isfinite(lam_lo)):
        return np.nan
    col = 0 if lam_lo <= 0.55 else (1 if lam_lo <= 0.65 else 2)
    for lo, vals in KAPPA_TABLE:
        if s >= lo:
            return vals[col]
    return KAPPA_TABLE[-1][1][col]


# --------------------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------------------
def check_sidecar(sidecar, root, required_files, prespec_rel=PRESPEC_REL):
    """Fail-closed.  Every line '<sha256> *<path relative to root>' is verified; the files in
    `required_files` (script, kernel, inputs, side files) and the pre-specification must all be listed."""
    if not os.path.exists(sidecar):
        sys.exit(f"REFUSED: sidecar {sidecar} not found -- the pre-specification is not frozen.")
    listed = {}
    for line in open(sidecar, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        want, rel = line.split(None, 1)
        rel = rel.strip().lstrip("*")
        path = os.path.normcase(os.path.abspath(os.path.join(root, rel)))
        if not os.path.exists(path) or sha256(path) != want.lower():
            sys.exit(f"REFUSED: {rel} missing or changed since freezing")
        listed[path] = want
    pre = os.path.normcase(os.path.abspath(os.path.join(root, prespec_rel)))
    if pre not in listed:
        sys.exit(f"REFUSED: the pre-specification {prespec_rel} is not listed in the sidecar")
    for f in required_files:
        if os.path.normcase(os.path.abspath(f)) not in listed:
            sys.exit(f"REFUSED: {f} is not listed in the sidecar (not frozen)")
    log("sidecar verified")


def is_synthetic(path):
    z = np.load(path, allow_pickle=False)
    return bool(z["is_synthetic"]) if "is_synthetic" in z.files else False


def side_files(args):
    return [os.path.abspath(p) for p in (args.sets_file, args.purity_file) if p]


def gate(args, sidecar, root):
    """Returns the mode: 'synthetic' | 'repro' | 'real'."""
    if args.repro_only:
        return "repro"
    flags = [is_synthetic(p) for p in args.inputs]
    if args.synthetic:
        if not all(flags):
            sys.exit("REFUSED: --synthetic is only for files flagged is_synthetic=True; "
                     "real data requires the frozen sidecar.")
        return "synthetic"
    required = [os.path.abspath(__file__), os.path.abspath(tk.__file__)] + [os.path.abspath(p) for p in args.inputs] \
        + side_files(args)
    check_sidecar(sidecar, root, required)
    return "real"


# --------------------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------------------
def load_inputs(path):
    z = np.load(path, allow_pickle=False)
    d = {k: z[k] for k in z.files}
    for k in ("patients", "genes"):
        d[k] = d[k].astype(str)
    # exact upcast: the published pipeline's hstack([float64 mRNA, float32 PCs]) already computes in float64
    d["pcs"] = d["pcs"].astype(np.float64)
    d["path"] = path
    d["tag"] = str(d["cohort"]) if "cohort" in d else os.path.basename(path)
    return d


def gene_arrays(d, g):
    v = ~np.isnan(d["protein"][:, g])
    return v, int(v.sum()), d["mrna_log2tpm1"][v, g], d["protein"][v, g]


def repro_published(d, tol=REPRO_TOL):
    """Re-derive the published (n, r2_rna, r2_both, incremental_r2) per gene from the exported
    inputs and compare with the published results stored beside them."""
    G = d["protein"].shape[1]
    if "pub_incremental_r2" not in d or not np.isfinite(d["pub_incremental_r2"]).any():
        return dict(verdict="NO PUBLISHED RESULTS in the inputs file", worst=None, ok=False)
    got = np.full((G, 4), np.nan)
    for g in range(G):
        v, nv, M, P = gene_arrays(d, g)
        got[g] = tk.published_increment(d["mrna_log2tpm1"][:, g], d["protein"][:, g], d["pcs"], v)
    pub = np.column_stack([d["pub_n"], d["pub_r2_rna"], d["pub_r2_both"], d["pub_incremental_r2"]])
    ok = np.isfinite(pub[:, 3])
    diffs = {nm: float(np.max(np.abs(got[ok, j] - pub[ok, j])))
             for j, nm in enumerate(["n", "r2_rna", "r2_both", "incremental_r2"])}
    worst = max(diffs["r2_rna"], diffs["r2_both"], diffs["incremental_r2"])
    good = worst <= tol and diffs["n"] == 0
    return dict(n_compared=int(ok.sum()), n_missing_from_published=int((~ok).sum()), max_abs_diff=diffs,
                worst=worst, ok=bool(good),
                verdict=("REPRODUCES the published per-gene results (<=%.0e)" % tol) if good
                else "DOES NOT REPRODUCE the published per-gene results: statistics refused")


def repro_pvalues(d, k):
    """Published permutation p-values of the first k fit genes (RA: RandomState(0), sequential,
    1000 draws, wsi_pcs[perm][valid], p = (1 + #null >= obs) / 1001)."""
    n = d["pcs"].shape[0]
    nperm = int(d["ra_n_perm"]) if "ra_n_perm" in d else 1000
    perms = tk.build_perms(n, nperm, "ra", int(d["ra_seed"]) if "ra_seed" in d else 0)
    idx = [g for g in range(d["protein"].shape[1]) if np.isfinite(d["pub_pval"][g])][:k]
    worst, mism = 0.0, 0
    for g in idx:
        v, nv, M, P = gene_arrays(d, g)
        inc = tk.published_increment(d["mrna_log2tpm1"][:, g], d["protein"][:, g], d["pcs"], v)[3]
        null = tk.published_null_increments(d["mrna_log2tpm1"][:, g], d["protein"][:, g], d["pcs"], perms, v)
        p = (np.sum(null >= inc) + 1) / (nperm + 1)
        diff = abs(p - d["pub_pval"][g])
        worst = max(worst, diff)
        mism += diff > 1e-12
    return dict(n_genes=len(idx), n_perm=nperm, max_abs_diff_pval=float(worst), n_mismatch=int(mism))


def check_pca(d):
    from sklearn.decomposition import PCA
    k = d["pcs"].shape[1]
    new = PCA(n_components=k, svd_solver="full").fit_transform(d["wsi_pooled"])
    return float(np.max(np.abs(np.abs(new) - np.abs(d["pcs"]))))


def check_icc(d, path):
    rp = os.path.join(os.path.dirname(os.path.abspath(path)), f"replicate_files_{d['tag']}.npz")
    if not os.path.exists(rp):
        return None
    r = np.load(rp, allow_pickle=False)
    cases = r["cases"]
    groups = [r["logtpm"][cases == c] for c in dict.fromkeys(cases.tolist())]
    groups = [g for g in groups if g.shape[0] >= 2]
    icc, _, a, kbar = tk.icc_oneway(groups)
    ok = np.isfinite(d["icc_log2tpm1"]) & np.isfinite(icc)
    return dict(n_cases=a, mean_files=kbar, max_abs_diff_icc=float(np.max(np.abs(icc[ok] - d["icc_log2tpm1"][ok]))) if ok.any() else None)


# --------------------------------------------------------------------------------------
# small numerics for the variants / diagnostics (row-wise ddof-0 moments, uniform weights)
# --------------------------------------------------------------------------------------
def _rcov(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return np.mean((a - a.mean(-1, keepdims=True)) * (b - b.mean(-1, keepdims=True)), axis=-1)


def _rcorr(a, b):
    with np.errstate(all="ignore"):
        return _rcov(a, b) / np.sqrt(_rcov(a, a) * _rcov(b, b))


def ns_basis(x, knots):
    """Natural cubic spline basis (Hastie-Tibshirani-Friedman 5.2.1) without the intercept: K knots -> K-1 columns
    (x, d_k - d_{K-1}); 4 knots = 3 df."""
    x = np.asarray(x, float)
    K = len(knots)

    def dk(k):
        return (np.maximum(x - knots[k], 0) ** 3 - np.maximum(x - knots[K - 1], 0) ** 3) / (knots[K - 1] - knots[k])
    cols = [x] + [dk(k) - dk(K - 2) for k in range(K - 2)]
    return np.column_stack(cols)


def _spline_design(Mz):
    q = np.quantile(Mz, [1 / 3, 2 / 3])
    knots = np.array([Mz.min(), q[0], q[1], Mz.max()])
    if not np.all(np.diff(knots) > 0):
        return None
    return ns_basis(Mz, knots)


def _ols_fit(X, y):
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return coef, X @ coef


def variant_parts(M, P, kind, R=None):
    """Baseline fit for a sensitivity variant on ONE gene's fit patients.  Returns (fM, r_p) where fM is the
    own-transcript component of the in-sample baseline fit and r_p the baseline residual (both of standardised P);
    N = cov(rp_hat, fM), N_c = corr(rp_hat, fM), D = cov(rp_hat, r_p).  'linear' reproduces the primary statistics
    exactly (fM = rho M)."""
    sM, sP = M.std(), P.std()
    if not (sM > 0 and sP > 0):
        return None
    Mz = (M - M.mean()) / sM
    Pz = (P - P.mean()) / sP
    if kind == "linear":
        rho = float(np.mean(Mz * Pz))
        return rho * Mz, Pz - rho * Mz
    one = np.ones((len(Mz), 1))
    if kind == "spline":
        B = _spline_design(Mz)
        if B is None:
            return None
        _, fit = _ols_fit(np.hstack([one, B]), Pz)
        return fit - fit.mean(), Pz - fit
    if kind == "rnapc":
        if R is None or len(Mz) < R.shape[1] + 12:
            return None
        Rc = R - R.mean(0)
        coef, fit = _ols_fit(np.hstack([one, Mz[:, None], Rc]), Pz)
        _, mfit = _ols_fit(np.hstack([one, Rc]), Mz)
        return coef[1] * (Mz - mfit), Pz - fit
    raise ValueError(kind)


# --------------------------------------------------------------------------------------
# compute
# --------------------------------------------------------------------------------------
_S = {}


def _init(mrna, prot, pcs, PP, wts_all, secondary, kind, diag, rnapcs):
    _S.update(mrna=mrna, prot=prot, pcs=pcs, PP=PP, wts=wts_all, sec=secondary, kind=kind, diag=diag, rna=rnapcs)


def _job(gi):
    """All draws of the current run for gene gi -> (gi, dict stat -> (B,)) or (gi, None)."""
    kind = _S["kind"]
    v = ~np.isnan(_S["prot"][:, gi])
    nv = int(v.sum())
    M, P = _S["mrna"][v, gi], _S["prot"][v, gi]
    fid = tk.fold_ids_for(nv)
    if kind == "stats":
        out = tk.gene_stats_batch(_S["PP"][:, v], M, P, fid, np.ones((1, nv)), _S["sec"])
        if out is not None and _S["diag"]:
            B = _S["PP"].shape[0]
            Mz = (M - M.mean()) / M.std()
            pred = tk.wonly_oof_weighted(_S["PP"][:, v], np.broadcast_to(Mz, (B, 1, nv)).copy(), fid, np.ones((1, nv)))[:, 0]
            out = dict(out)
            out["Wm"] = _rcorr(pred, Mz[None, :])
        return gi, out
    if kind == "boot":
        return gi, tk.gene_stats_batch(_S["pcs"][v][None], M, P, fid, _S["wts"][:, v], False)
    if kind == "delta":                                     # nested increment under the supplied W (published ridge on [M, PCs])
        PPv = _S["PP"][:, v]
        B = PPv.shape[0]
        Wst = np.concatenate([np.broadcast_to(M[None, :, None], (B, nv, 1)), PPv], axis=2)
        pred = tk.wonly_oof_weighted(Wst, np.broadcast_to(P, (B, 1, nv)).copy(), fid, np.ones((1, nv)))[:, 0]
        ss_tot = np.sum((P - P.mean()) ** 2)
        if not ss_tot > 0:
            return gi, None
        r2_both = 1 - np.sum((P[None, :] - pred) ** 2, axis=1) / ss_tot
        r2_rna = tk.r2_of(P, tk.oof_ridge(M.reshape(-1, 1), P, tk.folds_for(nv)))
        return gi, dict(dp=r2_both - r2_rna)
    if kind == "dpspl":                                     # spline-baseline nested increment, observed only
        sM = M.std()
        if not (sM > 0):
            return gi, None
        Mz = (M - M.mean()) / sM
        Bs = _spline_design(Mz)
        if Bs is None:
            return gi, None
        folds = tk.folds_for(nv)
        r2_b = tk.r2_of(P, tk.oof_ridge(Bs, P, folds))
        r2_w = tk.r2_of(P, tk.oof_ridge(np.hstack([Bs, _S["pcs"][v]]), P, folds))
        return gi, dict(dp=np.array([r2_w - r2_b]))
    if kind.startswith("var:"):
        R = _S["rna"][v] if _S["rna"] is not None else None
        parts = variant_parts(M, P, kind[4:], R)
        if parts is None:
            return gi, None
        fM, rp = parts
        PPv = _S["PP"][:, v]
        B = PPv.shape[0]
        pred = tk.wonly_oof_weighted(PPv, np.broadcast_to(rp, (B, 1, nv)).copy(), fid, np.ones((1, nv)))[:, 0]
        return gi, dict(N=_rcov(pred, fM[None, :]), Nc=_rcorr(pred, fM[None, :]), D=_rcov(pred, rp[None, :]))
    raise ValueError(kind)


def run_draws(d, PP, wts, keys, secondary, workers, kind="stats", genes=None, diag=False, rnapcs=None):
    """Run `kind` for the genes of `genes` (default all) over the draws PP (B,n,K) [or bootstrap weights wts].
    Returns dict key -> (B, G) float array (NaN for genes not run / degenerate)."""
    G = d["protein"].shape[1]
    B = (PP.shape[0] if PP is not None else wts.shape[0])
    res = {k: np.full((B, G), np.nan) for k in keys}
    initargs = (d["mrna_log2tpm1"], d["protein"], d["pcs"], PP, wts, secondary, kind, diag, rnapcs)
    gl = list(range(G)) if genes is None else [int(g) for g in genes]
    t0 = time.time()

    def take(gi, out):
        if out is not None:
            for k in keys:
                res[k][:, gi] = out[k]
    if workers <= 1:
        _init(*initargs)
        for i, gi in enumerate(gl, 1):
            take(*_job(gi))
            if i % 500 == 0:
                log(f"    {i}/{len(gl)} genes  {time.time()-t0:.0f}s")
    else:
        with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=initargs) as pool:
            for i, (gi, out) in enumerate(pool.imap_unordered(_job, gl, chunksize=8), 1):
                take(gi, out)
                if i % 1000 == 0:
                    log(f"    {i}/{len(gl)} genes  {time.time()-t0:.0f}s")
    return res


def observed(d, keys, secondary):
    """Observed statistics: REFERENCE (unbatched, published-style ridge) path, with the batched path
    asserted equal gene by gene.  Adds Ncs = Nc/|rho| = sign(rho) corr(rp_hat, M), the tested statistic."""
    G = d["protein"].shape[1]
    obs = {k: np.full(G, np.nan) for k in keys + ["rho"]}
    worst = 0.0
    for g in range(G):
        v, nv, M, P = gene_arrays(d, g)
        lit = tk.gene_stats_literal(d["pcs"][v], M, P, secondary)
        if lit is None:
            continue
        bat = tk.gene_stats_batch(d["pcs"][v][None], M, P, tk.fold_ids_for(nv), np.ones((1, nv)), secondary)
        obs["rho"][g] = lit["rho"]
        for k in keys:
            obs[k][g] = lit[k]
            worst = max(worst, abs(lit[k] - bat[k][0]))
    if worst > BATCH_TOL:
        sys.exit(f"FATAL: batched statistics differ from the reference implementation by {worst:.2e}")
    add_sign_weighted(obs)
    return obs, worst


def add_sign_weighted(st):
    """st['Ncs'] = st['Nc']/|rho|  (works for the (G,) observed dict and for a (B,G) null dict given st['rho'] (G,))."""
    rho = np.abs(st["rho"])
    with np.errstate(all="ignore"):
        st["Ncs"] = np.where(rho > 0, st["Nc"] / rho, np.nan)
    return st


def p_two(null_sum, obs_sum):
    B = len(null_sum)
    ge = (1 + np.sum(null_sum >= obs_sum)) / (B + 1)
    le = (1 + np.sum(null_sum <= obs_sum)) / (B + 1)
    return float(min(1.0, 2 * min(ge, le)))


def p_up(null_sum, obs_sum):
    return float((1 + np.sum(null_sum >= obs_sum)) / (len(null_sum) + 1))


def zsc(null_sum, obs_sum):
    sd = null_sum.std(ddof=1)
    return float((obs_sum - null_sum.mean()) / sd) if sd > 0 else np.nan


def _pct(x, q):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, q)) if len(x) else np.nan


def lam_of(theta):
    return theta / (1 + theta) if np.isfinite(theta) and theta > -1 else np.nan


# --------------------------------------------------------------------------------------
# sets and strata
# --------------------------------------------------------------------------------------
def miss_stratum_masks(nfit, n):
    nm = n - nfit
    f = nm / float(n)
    e = 1e-12
    return {"miss0": nm == 0,
            "miss0-5": (nm > 0) & (f <= 0.05 + e),
            "miss5-20": (f > 0.05 + e) & (f <= 0.20 + e),
            "miss20+": f > 0.20 + e,
            "miss<=5": f <= 0.05 + e}


def build_sets(d, ok, nfit, depth, sets_file):
    """name -> meta dict (idx, parent, kind, stratum, selected, n_sel, G_tested, s, n_complete, mnar_exposed).
    Base sets are built on COMPLETE genes (n_fit = n); strata as documented in the module header."""
    n = d["pcs"].shape[0]
    genes = d["genes"]
    G = len(genes)
    have_pub = "pub_fdr" in d and np.isfinite(d["pub_fdr"]).any() and "pub_incremental_r2" in d
    if have_pub:
        tested = np.isfinite(d["pub_fdr"]) & np.isfinite(d["pub_incremental_r2"])
        with np.errstate(invalid="ignore"):
            pinned = tested & (d["pub_fdr"] < 0.05) & (d["pub_incremental_r2"] > 0)
    else:
        tested = np.zeros(G, bool)
        pinned = np.zeros(G, bool)
    members = [("all", np.ones(G, bool), False)]
    if have_pub:
        members.append(("pinned_sig", np.ones(G, bool), True))
    if sets_file:
        sf = pd.read_csv(sets_file, sep="\t", dtype=str)
        pos = {g: i for i, g in enumerate(genes)}
        for name, grp in sf.groupby(sf.columns[0], sort=False):
            m = np.zeros(G, bool)
            for g in grp[sf.columns[1]]:
                if g in pos:
                    m[pos[g]] = True
            members.append((str(name), m, False))
            if have_pub:
                members.append((f"pinned_sig&{name}", m, True))
    ms = miss_stratum_masks(nfit, n)
    finite_depth = np.isfinite(depth)
    out = {}

    def add(name, parent, kind, stratum, selected, scope, member, n_complete, mnar):
        mask = scope & member & ok & (pinned if selected else True)
        idx = np.flatnonzero(mask)
        if len(idx) < (MIN_SET if kind == "primary" else MIN_STRATUM):
            return None
        uni = scope & member & tested
        out[name] = dict(idx=idx, parent=parent, kind=kind, stratum=stratum, selected=bool(selected),
                         n_sel=float(len(idx)) if selected else np.nan,
                         G_tested=float(uni.sum()) if have_pub else np.nan,
                         s=(len(idx) / uni.sum()) if (selected and have_pub and uni.sum() > 0) else np.nan,
                         n_complete=int(n_complete), mnar_exposed=bool(mnar))
        return idx

    for pname, member, selected in members:
        n_complete = int((ms["miss0"] & member & ok & (pinned if selected else True)).sum())
        mnar = n_complete < MIN_THETA
        prim = add(pname, pname, "primary", "miss0", selected, ms["miss0"], member, n_complete, mnar)
        if prim is not None and len(prim) >= MIN_STRATUM:
            cut = np.quantile(depth[prim], [1 / 3, 2 / 3])
            tert = np.full(G, -1)
            tert[finite_depth] = np.digitize(depth[finite_depth], cut)
            for j in range(3):
                add(f"{pname}|depthT{j+1}", pname, "depth", f"depthT{j+1}", selected, ms["miss0"] & (tert == j), member,
                    n_complete, mnar)
        for k in ("miss0-5", "miss5-20", "miss20+", "miss<=5"):
            add(f"{pname}|{k}", pname, "miss", k, selected, ms[k], member, n_complete, mnar)
    return out


def sets_depth(d):
    prot = d["protein"]
    cnt = np.where(~np.isnan(prot), d["counts"], np.nan)
    with np.errstate(all="ignore"):
        depth = np.nanmedian(cnt, axis=0)
    return depth


# --------------------------------------------------------------------------------------
# set-level rows
# --------------------------------------------------------------------------------------
def _lam_lo(r, args):
    """lambda_lo: set-level (median) replicate ICC; else the transferred value; else rho2bar.  Returns (value, source)."""
    icc = r.get("median_icc", np.nan)
    if np.isfinite(icc):
        return float(np.clip(icc, 1e-3, 0.99)), "icc"
    tr = None
    if args.lambda_lo_transfer_tertiles and r["stratum"].startswith("depthT"):
        parts = [float(x) for x in args.lambda_lo_transfer_tertiles.split(",")]
        tr = parts[int(r["stratum"][-1]) - 1]
    elif args.lambda_lo_transfer is not None:
        tr = args.lambda_lo_transfer
    if tr is not None:
        return float(np.clip(tr, 1e-3, 0.99)), "transferred"
    return float(np.clip(r["rho2bar"], 1e-3, 0.99)), "rho2bar"


def reading_numbers(r, args):
    """Fill the numbers the reading rules read from this row (no verdict)."""
    lam_lo, src = _lam_lo(r, args)
    r["lambda_lo"], r["lambda_lo_source"] = lam_lo, src
    lam_hi = r["median_lambda_poisson_ub"]
    r["lambda_hi"] = lam_hi
    if args.confirmatory or not r["selected"]:
        kap = 1.0
        r["kappa_source"] = "confirmatory" if args.confirmatory else "unselected"
    else:
        kap = kappa_of(r["s"], lam_lo)
        r["kappa_source"] = "table"
    r["kappa"] = kap
    tnm = kap * lam_lo / (1 - lam_lo) if np.isfinite(kap) else np.nan
    r["theta_noise_min"] = tnm
    th = r["theta_pc"]
    ok20 = r["n_genes"] >= MIN_THETA
    r["theta_reading_ok"] = bool(ok20)
    r["margin_theta_pc_minus_tnm"] = th - tnm
    r["margin_ci_hi_minus_tnm"] = r.get("theta_pc_ci_hi", np.nan) - tnm
    band = R1_BAND_LARGE if r["n_genes"] >= 100 else R1_BAND_SMALL
    r["r1_band_limit"] = band
    r["margin_abs_theta_pc_minus_band"] = abs(th) - band
    r["z_D_minus_gate"] = r["D_z"] - Z_GATE_D
    if np.isfinite(lam_hi) and lam_hi > 0 and np.isfinite(kap) and kap > 0:
        fx_lo = (1 - lam_hi) * th / lam_hi
        fx_hi = (1 - lam_lo) * th / (kap * lam_lo)
        r["fX_lo"], r["fX_hi"] = float(np.clip(fx_lo, 0, 1)), float(np.clip(fx_hi, 0, 1))
        r["post_share_lb"] = float(min(1.0, 1 - fx_hi))
    else:
        r["fX_lo"] = r["fX_hi"] = r["post_share_lb"] = np.nan
    return r


def set_rows(d, arm, sets, obs, null, boot, depth, lam_pois, with_sec, args, diag):
    rows = []
    Nc_null = null["Ncs"]
    lo = np.nanpercentile(Nc_null, 2.5, axis=0)
    hi = np.nanpercentile(Nc_null, 97.5, axis=0)
    icc = d["icc_log2tpm1"] if "icc_log2tpm1" in d else np.full(len(d["genes"]), np.nan)
    nfit = (~np.isnan(d["protein"])).sum(0)
    n = d["pcs"].shape[0]
    for name, meta in sets.items():
        idx = meta["idx"]
        if len(idx) < MIN_SET:
            continue
        r = dict(arm=arm, set=name, parent=meta["parent"], kind=meta["kind"], stratum=meta["stratum"],
                 selected=meta["selected"], n_genes=len(idx), n_complete=meta["n_complete"],
                 mnar_exposed=meta["mnar_exposed"], n_sel=meta["n_sel"], G_tested=meta["G_tested"], s=meta["s"],
                 median_n_fit=float(np.median(nfit[idx])), median_miss_frac=float(np.median(1 - nfit[idx] / n)))
        for label, k in (("N", "N"), ("Nc", "Ncs"), ("D", "D"), ("Nc_rhow", "Nc")):
            o = float(np.sum(obs[k][idx]))
            ns = null[k][:, idx].sum(1)
            r[f"{label}_set"] = o
            r[f"{label}_null_mean"] = float(ns.mean())
            r[f"{label}_null_sd"] = float(ns.std(ddof=1))
            r[f"{label}_z"] = zsc(ns, o)
            if label != "Nc_rhow":
                r[f"{label}_p2"] = p_two(ns, o)
                r[f"{label}_p_up"] = p_up(ns, o)
        th_raw = r["N_set"] / r["D_set"] if r["D_set"] != 0 else np.nan
        r["theta_raw"] = th_raw
        r["theta_set"] = th_raw                                    # kept name (raw, sensitivity)
        dN, dD = r["N_set"] - r["N_null_mean"], r["D_set"] - r["D_null_mean"]
        r["theta_pc"] = dN / dD if dD != 0 else np.nan
        r["lambda_hat_noise"] = lam_of(r["theta_pc"])
        r["lambda_hat_noise_raw"] = lam_of(th_raw)
        if boot is not None:
            bN = np.nansum(boot["N"][:, idx], axis=1)
            bD = np.nansum(boot["D"][:, idx], axis=1)
            with np.errstate(all="ignore"):
                tb = (bN - r["N_null_mean"]) / (bD - r["D_null_mean"])
                tr = bN / bD
            r["theta_pc_ci_lo"], r["theta_pc_ci_hi"] = _pct(tb, 2.5), _pct(tb, 97.5)
            tbf = tb[np.isfinite(tb)]
            r["theta_pc_boot_se"] = float(tbf.std(ddof=1)) if len(tbf) > 2 else np.nan
            r["theta_raw_ci_lo"], r["theta_raw_ci_hi"] = _pct(tr, 2.5), _pct(tr, 97.5)
            r["theta_ci_lo"], r["theta_ci_hi"] = r["theta_raw_ci_lo"], r["theta_raw_ci_hi"]
            lb = np.array([lam_of(x) for x in tb])
            r["lambda_ci_lo"], r["lambda_ci_hi"] = _pct(lb, 2.5), _pct(lb, 97.5)   # NaN when theta <= -1
            r["n_boot"] = len(tb)
        so = float(np.mean(obs["S_cross"][idx]))
        sn = null["S_cross"][:, idx].mean(1)
        r["S_cross_mean"], r["S_cross_null_mean"], r["S_cross_z"] = so, float(sn.mean()), zsc(sn, so)
        if with_sec:
            for k in ("S_cross_m", "S_pp"):
                o = float(np.mean(obs[k][idx]))
                nn = null[k][:, idx].mean(1)
                r[f"{k}_mean"], r[f"{k}_null_mean"], r[f"{k}_z"] = o, float(nn.mean()), zsc(nn, o)
        if diag:
            o = float(np.mean(obs["Wm"][idx]))
            nn = null["Wm"][:, idx].mean(1)
            r["Wm_mean"], r["Wm_null_mean"], r["Wm_z"] = o, float(nn.mean()), zsc(nn, o)
        r["frac_Nc_below_band"] = float(np.mean(obs["Ncs"][idx] < lo[idx]))
        r["frac_Nc_above_band"] = float(np.mean(obs["Ncs"][idx] > hi[idx]))
        r["rho2bar"] = float(np.mean(obs["rho"][idx] ** 2))
        r["median_cnt_depth"] = float(np.nanmedian(depth[idx]))
        r["median_icc"] = float(np.nanmedian(icc[idx])) if np.isfinite(icc[idx]).any() else np.nan
        r["mean_icc"] = float(np.nanmean(icc[idx])) if np.isfinite(icc[idx]).any() else np.nan
        r["median_lambda_poisson_ub"] = float(np.nanmedian(lam_pois[idx]))
        # implied post-transcriptional share at the depth-implied reliability bounds (descriptive)
        for tag, lam in (("icc", r["median_icc"]), ("poisson", r["median_lambda_poisson_ub"])):
            fx = (1 - lam) * r["theta_pc"] / lam if np.isfinite(lam) and lam > 0 else np.nan
            r[f"fX_if_lambda_{tag}"] = fx
        reading_numbers(r, args)
        rows.append(r)
    return rows


def pergene_frame(d, arm, obs, null, depth, lam_pois, keys, nfit, diag):
    n = d["pcs"].shape[0]
    df = pd.DataFrame({"arm": arm, "gene": d["genes"], "n_fit": nfit, "miss_frac": 1 - nfit / n, "rho": obs["rho"]})
    for k in keys + ["Ncs"]:
        df[k] = obs[k]
        df[f"{k}_null_mean"] = null[k].mean(0)
        df[f"{k}_null_sd"] = null[k].std(0, ddof=1)
    if diag:
        df["Wm"] = obs["Wm"]
        df["Wm_null_mean"] = null["Wm"].mean(0)
        df["Wm_null_sd"] = null["Wm"].std(0, ddof=1)
    with np.errstate(all="ignore"):
        df["theta_gene"] = obs["N"] / obs["D"]
    B = null["Ncs"].shape[0]
    ge = (1 + (null["Ncs"] >= obs["Ncs"][None]).sum(0)) / (B + 1)
    le = (1 + (null["Ncs"] <= obs["Ncs"][None]).sum(0)) / (B + 1)
    df["Nc_p2"] = np.minimum(1.0, 2 * np.minimum(ge, le))
    df["D_p_up"] = (1 + (null["D"] >= obs["D"][None]).sum(0)) / (B + 1)
    df["cnt_median"] = depth
    df["icc_log2tpm1"] = d["icc_log2tpm1"] if "icc_log2tpm1" in d else np.nan
    df["lambda_poisson_ub"] = lam_pois
    for c in ("pub_incremental_r2", "pub_fdr"):
        if c in d:
            df[c] = d[c]
    return df


def arm_perms(d, arm, B, seed):
    n = d["pcs"].shape[0]
    if arm == "none":
        return tk.build_perms(n, B, "ra", seed), dict(frac_permutable=1.0, n_levels=1)
    if arm == "operator_scanner":                           # PDAC: operator x scanner stratum
        for k in ("label_operator", "label_scanner"):
            if k not in d:
                sys.exit(f"REFUSED: inputs have no {k}")
        lab = np.array([f"{a}|{b}" for a, b in zip(d["label_operator"].astype(str), d["label_scanner"].astype(str))])
        key = "label_operator_scanner"
    else:
        key = f"label_{arm}"
        if key not in d:
            sys.exit(f"REFUSED: inputs have no {key}")
        lab = d[key].astype(str)
    st = json.loads(str(d["label_status"])) if "label_status" in d else {}
    if "_error" in st:
        sys.exit(f"REFUSED: label export failed on the cluster ({st['_error']}); stratified arm '{arm}' impossible")
    if (lab == "").all():
        sys.exit(f"REFUSED: all '{arm}' labels are empty")
    if (lab == "").any():
        log(f"[arm {arm}] WARNING: {(lab == '').sum()}/{n} patients unlabelled -> '_missing' stratum (sitepack convention)")
    perms, info = tk.build_perms_stratified(lab, B, seed)
    if info["n_levels"] < 2:
        sys.exit(f"REFUSED: axis '{arm}' has a single level; no within-batch null exists")
    return perms, info


# --------------------------------------------------------------------------------------
# purity
# --------------------------------------------------------------------------------------
def purity_residualised(d, purity_file):
    """Copy of d with M and P residualised (in-sample OLS with intercept, per gene, on the gene's fit patients that
    have a purity value) on the supplied purity.  Patients without purity leave the gene's fit set."""
    pf = pd.read_csv(purity_file, sep="\t", dtype={0: str})
    pur_map = dict(zip(pf.iloc[:, 0].astype(str), pd.to_numeric(pf.iloc[:, 1], errors="coerce")))
    pur = np.array([pur_map.get(p, np.nan) for p in d["patients"]], float)
    if not np.isfinite(pur).any():
        sys.exit("REFUSED: no patient of the inputs file has a purity value in --purity-file")
    dp = dict(d)
    M, P = d["mrna_log2tpm1"], d["protein"]
    M2 = np.full_like(M, np.nan)
    P2 = np.full_like(P, np.nan)
    for g in range(P.shape[1]):
        v = ~np.isnan(P[:, g]) & np.isfinite(pur)
        if v.sum() < 30:
            continue
        X = np.column_stack([np.ones(v.sum()), pur[v]])
        for src, dst in ((M, M2), (P, P2)):
            coef, fit = _ols_fit(X, src[v, g])
            dst[v, g] = src[v, g] - fit
    dp["mrna_log2tpm1"], dp["protein"] = M2, P2
    return dp, int(np.isfinite(pur).sum())


# --------------------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------------------
def variant_rows(kind, arm, sets, vobs, vnull, extra=None):
    """Set rows of a sensitivity variant.  The variants already return corr(rp_hat, fM), i.e. the sign-weighted N_c
    (no rho involved); theta_pc is centred on the variant's own null mean."""
    rows = []
    for name, meta in sets.items():
        idx = meta["idx"]
        idx = idx[np.isfinite(vobs["N"][idx]) & np.isfinite(vnull["N"][:, idx]).all(0) & np.isfinite(vnull["D"][:, idx]).all(0)]
        if len(idx) < MIN_SET:
            continue
        r = dict(variant=kind, arm=arm, set=name, n_genes=len(idx))
        for label, k in (("N", "N"), ("Nc", "Nc"), ("D", "D")):
            o = float(np.sum(vobs[k][idx]))
            ns = vnull[k][:, idx].sum(1)
            r[f"{label}_set"], r[f"{label}_null_mean"], r[f"{label}_z"] = o, float(ns.mean()), zsc(ns, o)
        r["theta_raw"] = r["N_set"] / r["D_set"] if r["D_set"] != 0 else np.nan
        dD = r["D_set"] - r["D_null_mean"]
        r["theta_pc"] = (r["N_set"] - r["N_null_mean"]) / dD if dD != 0 else np.nan
        r["lambda_hat"] = lam_of(r["theta_pc"])
        if extra is not None:
            r.update(extra(idx))
        rows.append(r)
    return rows


def analyse(args, mode, out_dir):
    d = load_inputs(args.inputs[0]) if len(args.inputs) == 1 else None
    if d is None:
        sys.exit("one cohort per call: pass a single --inputs file")
    tag = d["tag"]
    os.makedirs(out_dir, exist_ok=True)
    if args.max_genes:
        for k in ("genes", "mrna_log2tpm1", "protein", "counts", "tpm", "icc_log2tpm1",
                  "pub_n", "pub_r2_rna", "pub_r2_both", "pub_incremental_r2", "pub_pval", "pub_fdr"):
            if k in d:
                d[k] = d[k][..., :args.max_genes] if d[k].ndim == 2 else d[k][:args.max_genes]
    arms = [a.strip() for a in args.strata.split(",") if a.strip()]
    primary = args.primary_arm if args.primary_arm in arms else arms[0]
    variants = [] if args.variants in ("", "none") else [v.strip() for v in args.variants.split(",")]
    for v in variants:
        if v not in VARIANTS_ALL:
            sys.exit(f"unknown variant {v!r}")
    manifest = dict(mode=mode, cohort=tag, started=time.strftime("%Y-%m-%d %H:%M:%S"),
                    args={k: v for k, v in vars(args).items()}, primary_arm=primary, arms=arms,
                    sha256_inputs={p: sha256(p) for p in args.inputs},
                    sha256_script=sha256(__file__), sha256_kernel=sha256(tk.__file__),
                    sha256_side={p: sha256(p) for p in side_files(args)},
                    n=int(d["pcs"].shape[0]), G=int(d["protein"].shape[1]))
    n = d["pcs"].shape[0]

    # ---- published reproduction first ----
    rep = repro_published(d)
    manifest["repro_published"] = rep
    log(f"[repro] {rep['verdict']}  {rep.get('max_abs_diff', '')}")
    if not rep["ok"]:
        if mode == "synthetic" and rep["worst"] is None:
            log("[repro] (synthetic file without published results: skipped)")
        else:
            json.dump(manifest, open(os.path.join(out_dir, f"theta_{tag}_manifest.json"), "w"), indent=1, default=str)
            sys.exit("REFUSED: exported inputs do not reproduce the published per-gene results.")

    secondary = args.secondary
    diag = not args.no_diag
    keys = tk.STAT_KEYS + (tk.SEC_KEYS if secondary else [])
    log("[observed] reference path + batched assertion ...")
    obs, worst = observed(d, keys, secondary)
    manifest["batched_vs_reference_max_abs_diff"] = worst
    ok = np.isfinite(obs["N"]) & np.isfinite(obs["D"])
    nfit = (~np.isnan(d["protein"])).sum(0)
    log(f"[observed] {int(ok.sum())}/{len(ok)} genes with finite statistics; batched==reference to {worst:.1e}; "
        f"complete (n_fit = n): {int((ok & (nfit == n)).sum())}")

    prot = d["protein"]
    cnt_fit = np.where(~np.isnan(prot), d["counts"], np.nan)
    lam_pois = tk.lambda_poisson(cnt_fit)
    depth = sets_depth(d)
    sets = build_sets(d, ok, nfit, depth, args.sets_file)
    log("[sets] " + ", ".join(f"{k}:{len(v['idx'])}" for k, v in sets.items() if v["kind"] == "primary")
        + f" | strata reported: {sum(1 for v in sets.values() if v['kind'] != 'primary')}")
    sel_prim = [k for k, v in sets.items() if v["kind"] == "primary" and v["selected"]]
    for k in sel_prim:
        m = sets[k]
        log(f"[sets] {k}: complete n={len(m['idx'])}  s={m['s']:.3f}  mnar_exposed={m['mnar_exposed']}")
    union = np.unique(np.concatenate([v["idx"] for v in sets.values()])) if sets else np.array([], int)

    boot = None
    if args.boot > 0:
        log(f"[bootstrap] {args.boot} patient resamples over {len(union)} genes ...")
        wts = tk.boot_counts(n, args.boot, args.boot_seed)
        boot = run_draws(d, None, wts, ["N", "D"], False, args.workers, kind="boot", genes=union)

    rows, pergene, arms_info, nulls = [], [], {}, {}
    perms_by_arm = {}
    rho = obs["rho"]
    if diag:
        # observed W-only prediction of M: one "draw" with the real PCs, same code path as the null
        ob = run_draws(d, d["pcs"][None], None, ["Wm"], False, args.workers, kind="stats", diag=True)
        obs["Wm"] = ob["Wm"][0]
    for arm in arms:
        log(f"[arm {arm}] building {args.B} permutations ...")
        perms, info = arm_perms(d, arm, args.B, args.seed)
        arms_info[arm] = info
        perms_by_arm[arm] = perms
        PP = d["pcs"][perms]
        log(f"[arm {arm}] running null ({info})")
        null = run_draws(d, PP, None, keys + (["Wm"] if diag else []), secondary, args.workers, kind="stats", diag=diag)
        null["rho"] = rho
        add_sign_weighted(null)
        nulls[arm] = null
    for arm in arms:
        null = nulls[arm]
        # a gene with any non-finite null draw is dropped from set sums in this arm
        okarm = ok & np.isfinite(null["N"]).all(0) & np.isfinite(null["Nc"]).all(0) & np.isfinite(null["D"]).all(0) \
            & np.isfinite(null["Ncs"]).all(0)
        arm_sets = {}
        for k, v in sets.items():
            v2 = dict(v)
            v2["idx"] = v["idx"][okarm[v["idx"]]]
            arm_sets[k] = v2
        rows += set_rows(d, arm, arm_sets, obs, null, boot, depth, lam_pois, secondary, args, diag)
        pergene.append(pergene_frame(d, arm, obs, null, depth, lam_pois, keys, nfit, diag))
        if args.save_null:
            np.savez_compressed(os.path.join(out_dir, f"theta_{tag}_null_{arm}.npz"),
                                **{k: v.astype(np.float32) for k, v in null.items()}, perms=perms_by_arm[arm].astype(np.int16))
    sets_df = pd.DataFrame(rows)

    # ---- phi_op: between-operator share of the permutation-centred W-only signal ----
    phi_cols = {}
    if primary != "none" and "none" in nulls:
        log("[phi_op] D (and nested increment Delta_p) ...")
        reading_genes = np.unique(np.concatenate([v["idx"] for v in sets.values() if v["selected"]])) \
            if any(v["selected"] for v in sets.values()) else union
        dnull = {}
        dobs = None
        if args.phi_delta and len(reading_genes):
            for arm in (primary, "none"):
                PPa = d["pcs"][perms_by_arm[arm]]
                dnull[arm] = run_draws(d, PPa, None, ["dp"], False, args.workers, kind="delta", genes=reading_genes)["dp"]
            dobs = run_draws(d, d["pcs"][None], None, ["dp"], False, 1, kind="delta", genes=reading_genes)["dp"][0]
            manifest["nested_increment_obs_vs_published"] = float(
                np.nanmax(np.abs(dobs[reading_genes] - d["pub_incremental_r2"][reading_genes]))) \
                if "pub_incremental_r2" in d and np.isfinite(d["pub_incremental_r2"][reading_genes]).any() else None
        for name, v in sets.items():
            idx = v["idx"]
            idx = idx[np.isfinite(nulls[primary]["D"][:, idx]).all(0) & np.isfinite(nulls["none"]["D"][:, idx]).all(0)]
            if len(idx) < MIN_SET:
                continue
            e_w = float(nulls[primary]["D"][:, idx].sum(1).mean())
            e_u = float(nulls["none"]["D"][:, idx].sum(1).mean())
            den = float(obs["D"][idx].sum()) - e_u
            c = dict(phi_op_D=(e_w - e_u) / den if den != 0 else np.nan, E0_D_within=e_w, E0_D_unres=e_u)
            if dobs is not None and name in sets and np.isfinite(dobs[idx]).all() and np.isfinite(dnull[primary][:, idx]).all():
                ew = float(dnull[primary][:, idx].sum(1).mean())
                eu = float(dnull["none"][:, idx].sum(1).mean())
                dd = float(dobs[idx].sum()) - eu
                c.update(phi_op_dp=(ew - eu) / dd if dd != 0 else np.nan, dp_obs_set=float(dobs[idx].sum()),
                         E0_dp_within=ew, E0_dp_unres=eu)
            phi_cols[name] = c

    # ---- sensitivities (primary-arm permutations) ----
    var_dfs, var_by_set = [], {}
    pperm = perms_by_arm[primary]
    PPp = d["pcs"][pperm]
    genes_var = union
    for kind in variants:
        dd = d
        rn = None
        if kind == "purity":
            if not args.purity_file:
                log("[variant purity] NOT RUN: no --purity-file")
                manifest.setdefault("variants_not_run", []).append("purity (no --purity-file)")
                continue
            dd, npur = purity_residualised(d, args.purity_file)
            log(f"[variant purity] purity available for {npur}/{n} patients")
            kk = "var:linear"
        elif kind == "rnapc":
            if "rna_pcs20" not in d:
                log("[variant rnapc] NOT RUN: no rna_pcs20 in the inputs")
                manifest.setdefault("variants_not_run", []).append("rnapc (no rna_pcs20)")
                continue
            rn = d["rna_pcs20"].astype(np.float64)
            kk = "var:rnapc"
        else:
            kk = "var:spline"
        log(f"[variant {kind}] observed + {args.B} draws of the {primary} arm over {len(genes_var)} genes ...")
        vo = run_draws(dd, dd["pcs"][None], None, ["N", "Nc", "D"], False, args.workers, kind=kk, genes=genes_var, rnapcs=rn)
        vo = {k: x[0] for k, x in vo.items()}
        vn = run_draws(dd, PPp, None, ["N", "Nc", "D"], False, args.workers, kind=kk, genes=genes_var, rnapcs=rn)
        extra = None
        if kind == "spline":
            dsp = run_draws(d, d["pcs"][None], None, ["dp"], False, args.workers, kind="dpspl", genes=genes_var)["dp"][0]
            pubdp = d["pub_incremental_r2"] if "pub_incremental_r2" in d else np.full(len(d["genes"]), np.nan)

            def extra(idx, dsp=dsp, pubdp=pubdp):
                return dict(dp_spline_mean=float(np.nanmean(dsp[idx])),
                            dp_linear_published_mean=float(np.nanmean(pubdp[idx])) if np.isfinite(pubdp[idx]).any() else np.nan)
        vr = variant_rows(kind, primary, sets, vo, vn, extra)
        var_dfs.append(pd.DataFrame(vr))
        for r in vr:
            ent = {f"{kind}_{k}": r[k] for k in ("theta_pc", "theta_raw", "N_z", "Nc_z", "D_z", "lambda_hat", "n_genes")}
            if kind == "spline":
                ent["spline_dp"] = r["dp_spline_mean"]
                ent["dp_linear_published_mean"] = r["dp_linear_published_mean"]
            var_by_set.setdefault(r["set"], {}).update(ent)
    variants_df = pd.concat(var_dfs) if var_dfs else pd.DataFrame()

    # ---- reading table: primary arm numbers + unrestricted + phi_op + sensitivities ----
    reading = pd.DataFrame()
    if len(sets_df):
        pr = sets_df[sets_df["arm"] == primary].copy().set_index("set")
        keep_u = ["N_z", "Nc_z", "D_z", "theta_pc", "theta_raw", "S_cross_z"]
        if "none" in arms and primary != "none":
            un = sets_df[sets_df["arm"] == "none"].set_index("set")[keep_u].add_suffix("_unres")
            pr = pr.join(un)
        else:
            for k in keep_u:
                pr[f"{k}_unres"] = np.nan
        for nm in list(pr.index):
            for k, val in phi_cols.get(nm, {}).items():
                pr.loc[nm, k] = val
            for k, val in var_by_set.get(nm, {}).items():
                pr.loc[nm, k] = val
        # theta numbers are masked below the theta threshold in the reading table
        mask = pr["n_genes"] < MIN_THETA
        for c in [c for c in pr.columns if c.startswith("theta") or c.startswith("lambda_hat") or c.startswith("margin_")
                  or c in ("fX_lo", "fX_hi", "post_share_lb")
                  or c.endswith("_theta_pc") or c.endswith("_theta_raw")]:
            if c not in ("theta_reading_ok", "theta_noise_min"):
                pr.loc[mask, c] = np.nan
        reading = pr.reset_index()
        # R5 numbers: pairwise differences of theta_pc between strata of one parent and kind
        het = []
        for (par, kind), g in reading[reading["kind"].isin(["depth", "miss"])].groupby(["parent", "kind"]):
            g = g[np.isfinite(g["theta_pc"])]
            gl = g.to_dict("records")
            for i in range(len(gl)):
                for j in range(i + 1, len(gl)):
                    a, b = gl[i], gl[j]
                    se = max(a.get("theta_pc_boot_se", np.nan), b.get("theta_pc_boot_se", np.nan))
                    diff = a["theta_pc"] - b["theta_pc"]
                    het.append(dict(parent=par, kind=kind, stratum_a=a["stratum"], stratum_b=b["stratum"],
                                    theta_pc_a=a["theta_pc"], theta_pc_b=b["theta_pc"], diff=diff,
                                    se_a=a.get("theta_pc_boot_se", np.nan), se_b=b.get("theta_pc_boot_se", np.nan),
                                    abs_diff_over_larger_se=abs(diff) / se if se and np.isfinite(se) and se > 0 else np.nan))
        hetero = pd.DataFrame(het)
    else:
        hetero = pd.DataFrame()

    sets_df.to_csv(os.path.join(out_dir, f"theta_{tag}_sets.csv"), index=False)
    reading.to_csv(os.path.join(out_dir, f"theta_{tag}_reading.csv"), index=False)
    hetero.to_csv(os.path.join(out_dir, f"theta_{tag}_hetero.csv"), index=False)
    if len(variants_df):
        variants_df.to_csv(os.path.join(out_dir, f"theta_{tag}_variants.csv"), index=False)
    pd.concat(pergene).to_csv(os.path.join(out_dir, f"theta_{tag}_pergene.csv"), index=False)
    if boot is not None:
        np.savez_compressed(os.path.join(out_dir, f"theta_{tag}_boot.npz"), N=boot["N"].astype(np.float32),
                            D=boot["D"].astype(np.float32))
    manifest.update(arms_info=arms_info, finished=time.strftime("%Y-%m-%d %H:%M:%S"))
    json.dump(manifest, open(os.path.join(out_dir, f"theta_{tag}_manifest.json"), "w"), indent=1, default=str)
    print_reading(reading, primary)
    log(f"[done] wrote {out_dir}")
    sets_df.attrs["reading"] = reading
    sets_df.attrs["variants"] = variants_df
    sets_df.attrs["hetero"] = hetero
    return sets_df


def print_reading(reading, primary):
    """The numbers the final reading rules read (DERIVATION section 9), for the base sets; no rule is applied."""
    if not len(reading):
        return
    cols = ["set", "n_genes", "mnar_exposed", "s", "kappa", "lambda_lo", "lambda_lo_source", "lambda_hi", "z_D_minus_gate",
            "D_z", "Nc_z", "Nc_z_unres", "theta_pc", "theta_pc_ci_lo", "theta_pc_ci_hi", "theta_noise_min",
            "margin_theta_pc_minus_tnm", "margin_ci_hi_minus_tnm", "r1_band_limit", "margin_abs_theta_pc_minus_band",
            "phi_op_D", "Wm_z"]
    cols = [c for c in cols if c in reading.columns]
    base = reading[reading["kind"] == "primary"][cols]
    log(f"[reading numbers; primary arm = {primary}; numbers only, no verdict]")
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.float_format", "{:.3f}".format):
        print(base.to_string(index=False), flush=True)


def main(argv=None, sidecar=SIDECAR, root=PAPER):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs", nargs="+", required=True, help="inputs_<c>.npz from measurement_error_checks.py v2")
    ap.add_argument("--out", default=None, help="output dir (default: <inputs dir>/local_theta)")
    ap.add_argument("--repro-only", action="store_true", help="published quantities only; never computes a discriminator")
    ap.add_argument("--repro-pvals", type=int, default=0, help="with --repro-only: also re-derive K published permutation p-values")
    ap.add_argument("--check-pca", action="store_true")
    ap.add_argument("--check-icc", action="store_true")
    ap.add_argument("--synthetic", action="store_true", help="only for is_synthetic=True test files")
    ap.add_argument("--strata", default="operator,none", help="comma list of null arms: operator | none | plex | scanner | operator_scanner")
    ap.add_argument("--primary-arm", default="operator", help="arm the reading numbers are taken from (falls back to the first arm)")
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--seed", type=int, default=SEED_DEFAULT,
                    help="permutation RandomState seed; default 20261007 = FRESH draws (0 would be the published selection sequence)")
    ap.add_argument("--boot", type=int, default=200, help="patient bootstrap resamples for the theta CI (0 = off)")
    ap.add_argument("--boot-seed", type=int, default=BOOT_SEED_DEFAULT)
    ap.add_argument("--secondary", action="store_true", help="also S_cross_m and S_pp (with their nulls)")
    ap.add_argument("--sets-file", default=None, help="TSV set<TAB>gene (e.g. Reactome families)")
    ap.add_argument("--variants", default="spline,rnapc,purity", help="sensitivities: comma list of spline,rnapc,purity or 'none'")
    ap.add_argument("--purity-file", default=None, help="TSV patient<TAB>purity (for the purity variant)")
    ap.add_argument("--no-diag", action="store_true", help="skip the diagnostic corr(W-only prediction of M, M)")
    ap.add_argument("--no-phi-delta", dest="phi_delta", action="store_false", help="skip phi_op for the nested increment")
    ap.add_argument("--confirmatory", action="store_true", help="unselected confirmatory cohort: kappa = 1")
    ap.add_argument("--lambda-lo-transfer", type=float, default=None,
                    help="transferred lambda_lo (ICC of the matched CCRCC/UCEC depth stratum) when the cohort has no replicates")
    ap.add_argument("--lambda-lo-transfer-tertiles", default=None, help="a,b,c transferred lambda_lo for depth tertiles 1..3")
    ap.add_argument("--save-null", action="store_true")
    ap.add_argument("--max-genes", type=int, default=0, help="smoke test")
    ap.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)))
    args = ap.parse_args(argv)

    mode = gate(args, sidecar, root)
    out_dir = args.out or os.path.join(os.path.dirname(os.path.abspath(args.inputs[0])), "local_theta")
    if mode == "repro":
        d = load_inputs(args.inputs[0])
        rep = repro_published(d)
        log(f"[repro-only] {rep}")
        res = dict(repro=rep)
        if args.repro_pvals:
            res["pvals"] = repro_pvalues(d, args.repro_pvals)
            log(f"[repro-only] p-values: {res['pvals']}")
        if args.check_pca:
            res["pca_max_abs_diff_abs_scores"] = check_pca(d)
            log(f"[repro-only] PCA re-derived from wsi_pooled: max ||score|| diff {res['pca_max_abs_diff_abs_scores']:.2e} (informational)")
        if args.check_icc:
            res["icc"] = check_icc(d, args.inputs[0])
            log(f"[repro-only] ICC re-derived from replicate_files: {res['icc']}")
        os.makedirs(out_dir, exist_ok=True)
        json.dump(res, open(os.path.join(out_dir, f"repro_only_{d['tag']}.json"), "w"), indent=1, default=str)
        sys.exit(0 if rep["ok"] else 2)
    sets_df = analyse(args, mode, out_dir)
    if args.check_pca or args.check_icc:
        d = load_inputs(args.inputs[0])
        if args.check_pca:
            log(f"[check] PCA diff {check_pca(d):.2e}")
        if args.check_icc:
            log(f"[check] ICC {check_icc(d, args.inputs[0])}")
    return sets_df


if __name__ == "__main__":
    main()
