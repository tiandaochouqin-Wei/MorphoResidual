#!/usr/bin/env python3
"""
c1_luad_pergene.py -- (b1)/(b2) per-gene rediscovery, DESCRIPTIVE ONLY, per
review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 3. Neither version can enter
the section 5 reading; both are reported in one table whatever the outcome (section
3, "Descriptive only").

Reuses c1_run_test_luad.py's prepare() (same load/independence-check/gene-set path,
same required CLI arguments) for the usable universe, and its FROZEN_DUMMY_OPERATORS /
build_operator_dummy_design / load_operator_labels / build_within_operator_perms for
(b2)'s residualisation and null -- so "same frozen labels, same frozen dummy design"
is enforced by construction, not by re-typing the list a second time.

(b1) Uncorrected discovery estimator.
  K = 20 fixed (residual_analysis.N_PCS, NOT the K=round(20n/100) formula the primary
  family uses -- this reproduces the discovery cohort's own per-gene detection rule,
  which always used K=20).
  residual_analysis._process_gene / bh_fdr, UNCHANGED (imported, not reimplemented).
  B = 1,000 permutations shared across every gene (one permutation-index list, drawn
  once, exactly as discovery did it), from RandomState(crc32(b"c1_luad|pergene")).

(b2) Batch-corrected estimator.
  K = 20 fixed (residual_analysis_sitepack.N_PCS).
  PCs residualised on the FROZEN dummy set (c1_run_test_luad.build_operator_dummy_design).
  Estimator = sitepack's `_fit` + `plain_folds`, REIMPLEMENTED here (not imported --
  residual_analysis_sitepack.py's logic is inline in its main(), which cannot be
  pointed at confirmatory data; see c1_run_test_luad.py's own docstring for the same
  reasoning): per gene, on that gene's valid cases, 5-fold KFold(shuffle,
  random_state=0), ridge alpha=1.0 on features standardised within each training fold,
  R^2 clipped to [-5, 1], increment = R^2(RNA + residualised PCs) - R^2(RNA). This is
  the estimator that defined `incremental_r2_batchresid` and the 726-gene set.
  B = 1,000 WITHIN-OPERATOR permutations (c1_run_test_luad.build_within_operator_perms),
  from RandomState(crc32(b"c1_luad|pergene_op")); unlabelled cases held fixed.

For each version, reports: m (genes tested), BH rejections before and after the
increment>0 filter, overlap of the final count with the 2,310 and 726 sets against the
hypergeometric expectation, whether the count reaches 50 (the discovery positivity
criterion, itself set post hoc), and a censoring diagnostic (genes at the permutation p
floor, genes at <= 2/(B+1), and k* = ceil(m / (0.05*(B+1))), the minimum non-zero BH
rejection count at this B -- reaching 50 mostly records whether BH rejects anything at
all, not a finer distinction).

Run (smp; env from lsf_c1_luad_test.sh, same required args as c1_run_test_luad.py):
  bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_pergene_smoke.out \
    lsf_c1_luad_test.sh descriptive --perms 5 --tag smoke --blind [... same required paths]
  bsub -q smp -n 16 -R "span[hosts=1]" -o c1_luad_pergene.out \
    lsf_c1_luad_test.sh descriptive --perms 1000 [... same required paths]
(lsf_c1_luad_test.sh dispatches to this script when mode=descriptive, since
c1_run_test_luad.py's own "descriptive" mode is only a pointer to this file.)
"""
import argparse
import json
import math
import os
import sys
import time
import zlib

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import residual_analysis as RA  # noqa: E402
import c1_run_test_luad as C1  # noqa: E402  (prepare, sha256, frozen constants, operator design)

B1_SEED = b"c1_luad|pergene"
B2_SEED = b"c1_luad|pergene_op"
RIDGE_ALPHA = 1.0
R2_FLOOR = -5.0
CV_FOLDS = 5
K_FIXED = 20


def emp_p(obs, null):
    null = np.asarray(null, dtype=float)
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))


def hypergeom_expected(n_hit_pool, n_universe, n_draws):
    """Expected overlap count E[|A cap B|] for two subsets of size n_hit_pool and
    n_draws drawn from a universe of size n_universe (informational only)."""
    if n_universe == 0:
        return 0.0
    return n_hit_pool * n_draws / n_universe


def run_b1(common, tested, R, P, W, disc_incr, disc_fdr, workers, perms, log):
    log("\n=== (b1): uncorrected discovery estimator, K=20 fixed ===")
    n = len(common)
    n_pcs = min(K_FIXED, n - 1, W.shape[1])
    pcs = PCA(n_components=n_pcs, svd_solver="full").fit_transform(W)
    log(f"  PCA: {n_pcs} components (K fixed at {K_FIXED}, not the primary family's "
        f"round(20n/100) formula)")

    rng = np.random.RandomState(zlib.crc32(B1_SEED))
    perm_indices = [rng.permutation(n) for _ in range(perms)]
    rna_dict = {g: R[:, j] for j, g in enumerate(tested)}
    prot_dict = {g: P[:, j] for j, g in enumerate(tested)}

    import multiprocessing as mp
    results = []
    with mp.Pool(workers, initializer=RA._init_worker,
                 initargs=(rna_dict, prot_dict, pcs, perm_indices)) as pool:
        for r in pool.imap_unordered(RA._process_gene, tested, chunksize=8):
            if r is not None:
                results.append(r)
    if not results:
        sys.exit("FATAL (b1): no gene produced a result (every gene had <30 confirmatory values?)")
    res = pd.DataFrame(results)
    res["fdr"] = RA.bh_fdr(res["pval"].values)
    m = len(res)
    before_filter = int((res.fdr < 0.05).sum())
    sig = res[(res.fdr < 0.05) & (res.incremental_r2 > 0)]
    after_filter = len(sig)
    log(f"  (b1): m={m} genes tested; BH<0.05 before increment>0 filter: "
        f"{before_filter}; after: {after_filter}")
    return res, sig, m, before_filter, after_filter


def _fit_sitepack(X, y, folds):
    """Verbatim reimplementation of residual_analysis_sitepack._fit (module docstring
    and c1_run_test_luad.py's own docstring explain why this is copied, not imported).
    """
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
    return float(np.clip(r2, R2_FLOOR, 1.0))


def _plain_folds_sitepack(n):
    return list(KFold(n_splits=CV_FOLDS, shuffle=True, random_state=0).split(np.zeros((n, 1))))


# Module-level worker state for (b2)'s multiprocessing pool. MUST be module-level
# (not a closure) -- a nested function cannot be pickled by multiprocessing.Pool, on
# either fork or spawn start methods, since pickle locates a function by its
# module-qualified name. This mirrors residual_analysis._WORKER /
# split_replication_power._W, the project's existing pattern for exactly this reason.
_B2W = {}


def _b2_init(rna_arr, prot_arr, pcs_resid, perm_list, tested_genes, min_n):
    _B2W.update(R=rna_arr, P=prot_arr, PCS=pcs_resid, PERMS=perm_list,
               GENES=tested_genes, MIN_N=min_n)


def _b2_worker(j):
    y = _B2W["P"][:, j]
    v = ~np.isnan(y)
    nv = int(v.sum())
    if nv < _B2W["MIN_N"]:
        return None
    yv = y[v]
    xr = _B2W["R"][v, j].reshape(-1, 1)
    pf = _plain_folds_sitepack(nv)
    r2_rna = _fit_sitepack(xr, yv, pf)
    r2_full = _fit_sitepack(np.hstack([xr, _B2W["PCS"][v]]), yv, pf)
    incr = r2_full - r2_rna
    perms = _B2W["PERMS"]
    null = np.empty(len(perms))
    for i, p in enumerate(perms):
        null[i] = _fit_sitepack(np.hstack([xr, _B2W["PCS"][p][v]]), yv, pf) - r2_rna
    pval = emp_p(incr, null)
    return dict(gene=_B2W["GENES"][j], n=nv, r2_rna=r2_rna, r2_batchresid=r2_full,
               incremental_r2_batchresid=incr, pval=pval)


def run_b2(common, tested, R, P, W, disc_incr_bres, workers, perms, log):
    log("\n=== (b2): batch-corrected estimator, K=20 fixed, frozen operator dummies ===")
    n = len(common)
    n_pcs = min(K_FIXED, n - 1, W.shape[1])
    pcs = PCA(n_components=n_pcs, svd_solver="full").fit_transform(W)
    log(f"  PCA: {n_pcs} components (K fixed at {K_FIXED})")

    raw_op = C1.load_operator_labels(common, args_batch_labels_path, log)
    pcs_resid, dummy_levels = C1.build_operator_dummy_design(raw_op, pcs, log)
    perms_op, qualifies, permutable, levels, members = C1.build_within_operator_perms(
        raw_op, perms, log)
    if not qualifies:
        log("  operator stratum not testable on the realised analysed set -- (b2) is "
            "reported as 'not computable' (no other axis substituted).")
        return None, None, 0, 0, 0

    import multiprocessing as mp
    results = []
    with mp.Pool(workers, initializer=_b2_init,
                 initargs=(R, P, pcs_resid, perms_op, tested, C1.MIN_N)) as pool:
        for r in pool.imap_unordered(_b2_worker, range(len(tested)), chunksize=8):
            if r is not None:
                results.append(r)
    if not results:
        sys.exit("FATAL (b2): no gene produced a result.")
    res = pd.DataFrame(results)
    res["fdr"] = RA.bh_fdr(res["pval"].values)
    m = len(res)
    before_filter = int((res.fdr < 0.05).sum())
    sig = res[(res.fdr < 0.05) & (res.incremental_r2_batchresid > 0)]
    after_filter = len(sig)
    log(f"  (b2): m={m} genes tested; BH<0.05 before increment>0 filter: "
        f"{before_filter}; after: {after_filter}")
    return res, sig, m, before_filter, after_filter


def censoring_diag(res, m, perms):
    p_floor = 1.0 / (perms + 1)
    p_2nd = 2.0 / (perms + 1)
    n_at_floor = int((res.pval <= p_floor + 1e-12).sum())
    n_at_2nd = int((res.pval <= p_2nd + 1e-12).sum())
    k_star = math.ceil(m / (0.05 * (perms + 1))) if m else 0
    return dict(p_floor=p_floor, n_at_floor=n_at_floor, n_at_2nd_floor=n_at_2nd, k_star=k_star)


def summarise(name, res, sig, m, before, after, h1_set, h2_set, perms, log):
    if res is None:
        return dict(status="not_computable", m=0, before_filter=0, after_filter=0)
    n_sig = after
    overlap_h1 = len(set(sig.gene) & h1_set)
    overlap_h2 = len(set(sig.gene) & h2_set)
    exp_h1 = hypergeom_expected(len(h1_set & set(res.gene)), len(res), n_sig)
    exp_h2 = hypergeom_expected(len(h2_set & set(res.gene)), len(res), n_sig)
    diag = censoring_diag(res, m, perms)
    censored = (n_sig == 0 and diag["n_at_floor"] > 0)
    log(f"\n[{name}] count={n_sig} (before-filter BH<0.05: {before}); "
        f">=50 (discovery positivity criterion): {n_sig >= 50}")
    log(f"[{name}] overlap with 2,310 H1 set: {overlap_h1} (hypergeometric expectation "
        f"{exp_h1:.1f}); overlap with 726 H2 set: {overlap_h2} (expectation {exp_h2:.1f})")
    log(f"[{name}] censoring: {diag['n_at_floor']} genes at p-floor "
        f"{diag['p_floor']:.4g}, {diag['n_at_2nd_floor']} at <= 2x floor, "
        f"k*={diag['k_star']} (min non-zero BH count at this B) -> "
        + ("ZERO AT B={0}, LIMITED BY THE PERMUTATION FLOOR; NOT EVIDENCE OF ABSENCE".format(perms)
           if censored else
           ("resolution-limited at B={0}".format(perms) if diag["n_at_floor"] > 0 else "not censored")))
    return dict(status="computed", m=m, before_filter=before, after_filter=after,
               overlap_h1=overlap_h1, overlap_h2=overlap_h2,
               hypergeometric_expected_h1=exp_h1, hypergeometric_expected_h2=exp_h2,
               reaches_50=bool(n_sig >= 50), censoring=diag, censored=bool(censored))


args_batch_labels_path = None  # set in main(), used by run_b2 via module global for brevity


def main():
    global args_batch_labels_path
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=RA.N_WORKERS)
    ap.add_argument("--tag", default="")
    ap.add_argument("--blind", action="store_true")
    ap.add_argument("--which", choices=["b1", "b2", "both"], default="both")
    # same required paths as c1_run_test_luad.py
    ap.add_argument("--rule", required=True)
    ap.add_argument("--rule-sha256", required=True)
    ap.add_argument("--case-ids", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--discovery-results", required=True)
    ap.add_argument("--discovery-sitepack", required=True)
    ap.add_argument("--discovery-slide-map", required=True)
    ap.add_argument("--discovery-xwalk", required=True)
    ap.add_argument("--discovery-rna-manifest", required=True)
    ap.add_argument("--batch-labels", required=True)
    args = ap.parse_args()
    args_batch_labels_path = args.batch_labels

    def log(msg):
        print(msg, flush=True)

    got = C1.sha256(args.rule)
    if got != args.rule_sha256:
        sys.exit(f"FATAL: --rule sha256 mismatch: got={got} want={args.rule_sha256}")

    t0 = time.time()
    common, tested, sel_h1, sel_h2, disc_incr, disc_fdr, disc_incr_bres, R, P, W, K = \
        C1.prepare(args, log)
    n = len(common)

    if args.blind:
        log(f"\nBLIND SMOKE (pergene): n={n} n_tested={len(tested)}")
        return

    out = {}
    if args.which in ("b1", "both"):
        res1, sig1, m1, b1_before, b1_after = run_b1(
            common, tested, R, P, W, disc_incr, disc_fdr, args.workers, args.perms, log)
        out["b1"] = summarise("b1", res1, sig1, m1, b1_before, b1_after, sel_h1, sel_h2,
                              args.perms, log)
        suffix = f"_{args.tag}" if args.tag else ""
        res1.to_csv(RA.OUT_DIR / f"c1_luad_pergene_b1{suffix}.csv", index=False)
    if args.which in ("b2", "both"):
        res2, sig2, m2, b2_before, b2_after = run_b2(
            common, tested, R, P, W, disc_incr_bres, args.workers, args.perms, log)
        out["b2"] = summarise("b2", res2, sig2, m2, b2_before, b2_after, sel_h1, sel_h2,
                              args.perms, log)
        if res2 is not None:
            suffix = f"_{args.tag}" if args.tag else ""
            res2.to_csv(RA.OUT_DIR / f"c1_luad_pergene_b2{suffix}.csv", index=False)

    result = dict(cohort="luad_c1", n_cases=n, perms=args.perms, **out,
                 minutes=round((time.time() - t0) / 60.0, 1),
                 sha256_this_script=C1.sha256(os.path.abspath(__file__)))
    suffix = f"_{args.tag}" if args.tag else ""
    fp = RA.OUT_DIR / f"c1_luad_pergene_result{suffix}.json"
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)
    log(f"\n  -> {fp}\n  elapsed {result['minutes']} min")


if __name__ == "__main__":
    main()
