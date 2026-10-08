#!/usr/bin/env python3
"""
local_theta.py -- LOCAL computation of the own-transcript-alignment statistics of
theory/DERIVATION.md (N, N_c, D, theta, lambda_hat_noise, S_cross) from the files exported by
measurement_error_checks.py v2 (inputs_<c>.npz).  Drafted 2026-10-08.  NOT YET RUN on project data.

STATISTICS (per gene, same folds / standardisation / ridge alpha as the published estimator --
theta_kernel.py, which is asserted equal to RA.cv_r2's ridge):
    r_p = P - rho M ,  r_m = M - rho P  (M, P standardised over the gene's fit patients)
    rp_hat = W-ONLY 5-fold OOF ridge prediction of r_p  (20 morphology PCs; NOT the nested d_p)
    N = rho cov(rp_hat, M),  N_c = rho corr(rp_hat, M),  D = cov(rp_hat, r_p),
    theta = sum N / sum D over a gene set,  lambda_hat_noise = theta / (1 + theta),
    S_cross = sign(rho) corr(rp_hat, r_m).      [--secondary adds S_cross_m, S_pp]
NULL: joint patient<->PC-row permutation, ONE draw shared by all genes and all targets, folds /
standardisation / alpha fixed.  Arms (--strata): none = the published unstratified
patient<->slide permutation (RandomState(0), sequential, the very draws of residual_analysis.py);
operator / plex = permutation WITHIN batch levels (residual_analysis_sitepack.py construction,
RandomState(0), levels sorted, no pooling).  Set level: sum N_c (the TESTED quantity), sum N, sum D
against the summed null; theta_set with a patient bootstrap CI (resample the n patients with
replacement, shared across genes; each patient keeps the fold it has in the gene's published
KFold, so copies never straddle train/test; implemented by integer weights == duplicated rows).
Per gene: the fraction of genes whose N_c lies below / above the 2.5 / 97.5 % band of their OWN
null (null expectation 2.5 % each) -- proportions, not gene calls.

FAIL-CLOSED (same idea as server_export/scripts/c1_ucec_posthoc_sets.py).  On any inputs file that
does not carry is_synthetic=True the statistics run ONLY if the sidecar named in SIDECAR exists and
every file listed in it (the frozen pre-specification, this script, theta_kernel.py, every inputs
file used) still has its listed sha256 -- and this script, theta_kernel.py and each inputs file used
must themselves be listed.  Exceptions that never compute a discriminator:
    --repro-only    published quantities only: re-derives the published per-gene incremental R^2
                    from the exported inputs (optionally --repro-pvals K: the published
                    permutation p-values of K genes) and checks them to 1e-10.
    --synthetic     allowed only for files flagged is_synthetic=True (test data).
Before any statistic the full run ALSO re-derives the published increments for all genes and
refuses if they differ from the published results by > 1e-10, and asserts batched == reference
numerics gene by gene.

Usage:
   python local_theta.py --inputs export/ucec/inputs_ucec.npz --repro-only
   python local_theta.py --inputs export/ucec/inputs_ucec.npz --strata none,operator,plex --B 200 --boot 200
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
MIN_SET = 5
MIN_SPLIT = 15


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------------------
def check_sidecar(sidecar, root, required_files, prespec_rel=PRESPEC_REL):
    """Fail-closed.  Every line '<sha256> *<path relative to root>' is verified; the files in
    `required_files` (script, kernel, inputs) and the pre-specification must all be listed."""
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
    required = [os.path.abspath(__file__), os.path.abspath(tk.__file__)] + [os.path.abspath(p) for p in args.inputs]
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
# compute
# --------------------------------------------------------------------------------------
_S = {}


def _init(mrna, prot, pcs, PP, wts_all, secondary):
    _S.update(mrna=mrna, prot=prot, pcs=pcs, PP=PP, wts=wts_all, sec=secondary)


def _job(gi):
    """All draws of the current arm for gene gi -> (gi, dict stat -> (B,)) or (gi, None)."""
    v = ~np.isnan(_S["prot"][:, gi])
    nv = int(v.sum())
    M, P = _S["mrna"][v, gi], _S["prot"][v, gi]
    fid = tk.fold_ids_for(nv)
    if _S["PP"] is not None:                                   # permutation arm
        out = tk.gene_stats_batch(_S["PP"][:, v], M, P, fid, np.ones((1, nv)), _S["sec"])
    else:                                                      # bootstrap
        out = tk.gene_stats_batch(_S["pcs"][v][None], M, P, fid, _S["wts"][:, v], _S["sec"])
    return gi, out


def run_draws(d, PP, wts, keys, secondary, workers):
    G = d["protein"].shape[1]
    B = (PP.shape[0] if PP is not None else wts.shape[0])
    res = {k: np.full((B, G), np.nan) for k in keys}
    initargs = (d["mrna_log2tpm1"], d["protein"], d["pcs"], PP, wts, secondary)
    t0 = time.time()

    def take(gi, out):
        if out is not None:
            for k in keys:
                res[k][:, gi] = out[k]
    if workers <= 1:
        _init(*initargs)
        for gi in range(G):
            take(*_job(gi))
            if (gi + 1) % 500 == 0:
                log(f"    {gi+1}/{G} genes  {time.time()-t0:.0f}s")
    else:
        with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=initargs) as pool:
            for i, (gi, out) in enumerate(pool.imap_unordered(_job, range(G), chunksize=8), 1):
                take(gi, out)
                if i % 1000 == 0:
                    log(f"    {i}/{G} genes  {time.time()-t0:.0f}s")
    return res


def observed(d, keys, secondary):
    """Observed statistics: REFERENCE (unbatched, published-style ridge) path, with the batched path
    asserted equal gene by gene."""
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
    return obs, worst


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


def build_sets(d, ok, sets_file):
    """name -> index array.  Built-in: all, pinned_sig (published selection FDR<0.05 & incr>0);
    --sets-file (TSV set<TAB>gene); each set also split into depth tertiles (raw-count median)."""
    genes = d["genes"]
    base = {"all": np.flatnonzero(ok)}
    if "pub_fdr" in d and np.isfinite(d["pub_fdr"]).any():
        sig = (d["pub_fdr"] < 0.05) & (d["pub_incremental_r2"] > 0)
        base["pinned_sig"] = np.flatnonzero(sig & ok)
    if sets_file:
        sf = pd.read_csv(sets_file, sep="\t", dtype=str)
        pos = {g: i for i, g in enumerate(genes)}
        for name, grp in sf.groupby(sf.columns[0], sort=False):
            idx = np.array(sorted(pos[g] for g in grp[sf.columns[1]] if g in pos and ok[pos[g]]), dtype=int)
            base[str(name)] = idx
    prot = d["protein"]
    cnt = np.where(~np.isnan(prot), d["counts"], np.nan)
    with np.errstate(all="ignore"):
        depth = np.nanmedian(cnt, axis=0)
    out = dict(base)
    for name, idx in base.items():
        if len(idx) >= MIN_SPLIT:
            q = np.quantile(depth[idx], [1 / 3, 2 / 3])
            t = np.digitize(depth[idx], q)               # 0,1,2
            for j in range(3):
                out[f"{name}|depthT{j+1}"] = idx[t == j]
    return out, depth


def set_rows(d, arm, sets, obs, null, boot, depth, lam_pois, with_sec):
    rows = []
    B = null["N"].shape[0]
    Nc_null = null["Nc"]
    lo = np.nanpercentile(Nc_null, 2.5, axis=0)
    hi = np.nanpercentile(Nc_null, 97.5, axis=0)
    icc = d["icc_log2tpm1"] if "icc_log2tpm1" in d else np.full(len(d["genes"]), np.nan)
    for name, idx in sets.items():
        if len(idx) < MIN_SET:
            continue
        r = dict(arm=arm, set=name, n_genes=len(idx))
        for k, label in (("N", "N"), ("Nc", "Nc"), ("D", "D")):
            o = float(np.sum(obs[k][idx]))
            ns = null[k][:, idx].sum(1)
            r[f"{label}_set"] = o
            r[f"{label}_null_mean"] = float(ns.mean())
            r[f"{label}_null_sd"] = float(ns.std(ddof=1))
            r[f"{label}_z"] = zsc(ns, o)
            r[f"{label}_p2"] = p_two(ns, o)
            r[f"{label}_p_up"] = p_up(ns, o)
        th = r["N_set"] / r["D_set"] if r["D_set"] != 0 else np.nan
        r["theta_set"] = th
        r["lambda_hat_noise"] = lam_of(th)
        if boot is not None:
            bN = np.nansum(boot["N"][:, idx], axis=1)
            bD = np.nansum(boot["D"][:, idx], axis=1)
            with np.errstate(all="ignore"):
                tb = bN / bD
            r["theta_ci_lo"], r["theta_ci_hi"] = _pct(tb, 2.5), _pct(tb, 97.5)
            lb = np.array([lam_of(x) for x in tb])
            r["lambda_ci_lo"], r["lambda_ci_hi"] = _pct(lb, 2.5), _pct(lb, 97.5)   # NaN when theta <= -1 (H_err signature)
            r["n_boot"] = len(tb)
        so = float(np.mean(obs["S_cross"][idx]))
        sn = null["S_cross"][:, idx].mean(1)
        r["S_cross_mean"], r["S_cross_null_mean"], r["S_cross_z"] = so, float(sn.mean()), zsc(sn, so)
        if with_sec:
            for k in ("S_cross_m", "S_pp"):
                o = float(np.mean(obs[k][idx]))
                nn = null[k][:, idx].mean(1)
                r[f"{k}_mean"], r[f"{k}_null_mean"], r[f"{k}_z"] = o, float(nn.mean()), zsc(nn, o)
        r["frac_Nc_below_band"] = float(np.mean(obs["Nc"][idx] < lo[idx]))
        r["frac_Nc_above_band"] = float(np.mean(obs["Nc"][idx] > hi[idx]))
        r["median_cnt_depth"] = float(np.nanmedian(depth[idx]))
        r["median_icc"] = float(np.nanmedian(icc[idx])) if np.isfinite(icc[idx]).any() else np.nan
        r["median_lambda_poisson_ub"] = float(np.nanmedian(lam_pois[idx]))
        # implied post-transcriptional share at the depth-implied reliability bounds (descriptive)
        for tag, lam in (("icc", r["median_icc"]), ("poisson", r["median_lambda_poisson_ub"])):
            fx = (1 - lam) * th / lam if np.isfinite(lam) and lam > 0 else np.nan
            r[f"fX_if_lambda_{tag}"] = fx
        rows.append(r)
    return rows


def pergene_frame(d, arm, obs, null, depth, lam_pois, keys):
    G = len(d["genes"])
    df = pd.DataFrame({"arm": arm, "gene": d["genes"], "n_fit": (~np.isnan(d["protein"])).sum(0),
                       "rho": obs["rho"]})
    for k in keys:
        df[k] = obs[k]
        df[f"{k}_null_mean"] = null[k].mean(0)
        df[f"{k}_null_sd"] = null[k].std(0, ddof=1)
    with np.errstate(all="ignore"):
        df["theta_gene"] = obs["N"] / obs["D"]
    B = null["Nc"].shape[0]
    ge = (1 + (null["Nc"] >= obs["Nc"][None]).sum(0)) / (B + 1)
    le = (1 + (null["Nc"] <= obs["Nc"][None]).sum(0)) / (B + 1)
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
    key = f"label_{arm}"
    if key not in d:
        sys.exit(f"REFUSED: inputs have no {key}")
    st = json.loads(str(d["label_status"])) if "label_status" in d else {}
    if "_error" in st:
        sys.exit(f"REFUSED: label export failed on the cluster ({st['_error']}); stratified arm '{arm}' impossible")
    lab = d[key].astype(str)
    if (lab == "").all():
        sys.exit(f"REFUSED: all '{arm}' labels are empty")
    if (lab == "").any():
        log(f"[arm {arm}] WARNING: {(lab == '').sum()}/{n} patients unlabelled -> '_missing' stratum (sitepack convention)")
    perms, info = tk.build_perms_stratified(lab, B, seed)
    if info["n_levels"] < 2:
        sys.exit(f"REFUSED: axis '{arm}' has a single level; no within-batch null exists")
    return perms, info


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
    manifest = dict(mode=mode, cohort=tag, started=time.strftime("%Y-%m-%d %H:%M:%S"),
                    args={k: v for k, v in vars(args).items()},
                    sha256_inputs={p: sha256(p) for p in args.inputs},
                    sha256_script=sha256(__file__), sha256_kernel=sha256(tk.__file__),
                    n=int(d["pcs"].shape[0]), G=int(d["protein"].shape[1]))

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
    keys = tk.STAT_KEYS + (tk.SEC_KEYS if secondary else [])
    log("[observed] reference path + batched assertion ...")
    obs, worst = observed(d, keys, secondary)
    manifest["batched_vs_reference_max_abs_diff"] = worst
    ok = np.isfinite(obs["N"]) & np.isfinite(obs["D"])
    log(f"[observed] {int(ok.sum())}/{len(ok)} genes with finite statistics; batched==reference to {worst:.1e}")

    prot = d["protein"]
    cnt_fit = np.where(~np.isnan(prot), d["counts"], np.nan)
    lam_pois = tk.lambda_poisson(cnt_fit)
    sets, depth = build_sets(d, ok, args.sets_file)
    log(f"[sets] " + ", ".join(f"{k}:{len(v)}" for k, v in sets.items() if "|" not in k))

    boot = None
    if args.boot > 0:
        log(f"[bootstrap] {args.boot} patient resamples ...")
        wts = tk.boot_counts(d["pcs"].shape[0], args.boot, args.boot_seed)
        boot = run_draws(d, None, wts, ["N", "D"], False, args.workers)

    rows, pergene, arms_info = [], [], {}
    for arm in args.strata.split(","):
        arm = arm.strip()
        log(f"[arm {arm}] building {args.B} permutations ...")
        perms, info = arm_perms(d, arm, args.B, args.seed)
        arms_info[arm] = info
        PP = d["pcs"][perms]
        log(f"[arm {arm}] running null ({info})")
        null = run_draws(d, PP, None, keys, secondary, args.workers)
        # a gene with any non-finite null draw is dropped from set sums in this arm
        okarm = ok & np.isfinite(null["N"]).all(0) & np.isfinite(null["Nc"]).all(0) & np.isfinite(null["D"]).all(0)
        arm_sets = {k: v[okarm[v]] for k, v in sets.items()}
        rows += set_rows(d, arm, arm_sets, obs, null, boot, depth, lam_pois, secondary)
        pergene.append(pergene_frame(d, arm, obs, null, depth, lam_pois, keys))
        if args.save_null:
            np.savez_compressed(os.path.join(out_dir, f"theta_{tag}_null_{arm}.npz"),
                                **{k: v.astype(np.float32) for k, v in null.items()}, perms=perms.astype(np.int16))
    sets_df = pd.DataFrame(rows)
    sets_df.to_csv(os.path.join(out_dir, f"theta_{tag}_sets.csv"), index=False)
    pd.concat(pergene).to_csv(os.path.join(out_dir, f"theta_{tag}_pergene.csv"), index=False)
    if boot is not None:
        np.savez_compressed(os.path.join(out_dir, f"theta_{tag}_boot.npz"), N=boot["N"].astype(np.float32),
                            D=boot["D"].astype(np.float32))
    manifest.update(arms=arms_info, finished=time.strftime("%Y-%m-%d %H:%M:%S"))
    json.dump(manifest, open(os.path.join(out_dir, f"theta_{tag}_manifest.json"), "w"), indent=1, default=str)
    log(f"[done] wrote {out_dir}")
    return sets_df


def main(argv=None, sidecar=SIDECAR, root=PAPER):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs", nargs="+", required=True, help="inputs_<c>.npz from measurement_error_checks.py v2")
    ap.add_argument("--out", default=None, help="output dir (default: <inputs dir>/local_theta)")
    ap.add_argument("--repro-only", action="store_true", help="published quantities only; never computes a discriminator")
    ap.add_argument("--repro-pvals", type=int, default=0, help="with --repro-only: also re-derive K published permutation p-values")
    ap.add_argument("--check-pca", action="store_true")
    ap.add_argument("--check-icc", action="store_true")
    ap.add_argument("--synthetic", action="store_true", help="only for is_synthetic=True test files")
    ap.add_argument("--strata", default="none,operator,plex", help="comma list of null arms: none | operator | plex")
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0, help="permutation RandomState seed (0 = the published sequence)")
    ap.add_argument("--boot", type=int, default=200, help="patient bootstrap resamples for theta_set CI (0 = off)")
    ap.add_argument("--boot-seed", type=int, default=20261008)
    ap.add_argument("--secondary", action="store_true", help="also S_cross_m and S_pp (with their nulls)")
    ap.add_argument("--sets-file", default=None, help="TSV set<TAB>gene (e.g. Reactome families)")
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
