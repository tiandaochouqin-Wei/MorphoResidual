"""
sim_batchnull.py -- which permutation null for the own-transcript alignment test N_c:
unrestricted patient<->slide permutation (as published) vs operator(batch)-stratified.
Synthetic only; NO project data.

Generative model = sim_signtest.py plus an acquisition-batch (operator) structure:
    operator op(i) in {1..n_ops} (uneven sizes) + n_single singleton patients
    W_i  = sqrt(nu_W) Z_{op(i)} + sqrt(1 - nu_W) W0_i          Z_k ~ N(0, I_K): operator offset in PC space
    X_g  = sqrt(nu_X) (Z_{op(i)} v_g) + sqrt(1 - nu_X) [ s (W0 v_g)/sqrt(1-nu_W) + sqrt(1-s^2) xi ]
           (nu_X > 0: "accrual" variant, between-operator transcript differences that morphology can read)
    e_g  = sqrt((1-lam)/lam) [ sqrt(nu_e) zeta_{g,op(i)} + sqrt(1-nu_e) (q (W0 z_g)/sqrt(1-nu_W) + sqrt(1-q^2) eta) ]
           (nu_e > 0: gene-specific RNA batch effect shared by the operator's patients)
    u_g  = sigma_u [ t (W0 w_g)/sqrt(1-nu_W) + sqrt(1-t^2) zeta ]
    eps_g = sigma_eps [ sqrt(nu_P) psi_{g,op(i)} + sqrt(1-nu_P) eps0 ]   (nu_P > 0: protein-side site effect)
    M = X + e,  P = b X + u + eps.
The biological signal lives in the idiosyncratic part W0 of the morphology block (so that the
operator offset carries no gene signal unless nu_X > 0); the signal strengths are scaled so that
the population increment equals the target after the batch dilution (rho_W = s^2 (1 - nu_W)).

For every replicate: observed statistics on the RAW W (as the pinned sets were selected),
selection by the published rule under the UNRESTRICTED null, then the discriminators' joint null
(a) unrestricted and (b) within-operator (singletons fixed), and theta on operator-residualised W.
Usage:  python sim_batchnull.py --out sim_batchnull.csv [--reps 8 --G 200 --B 200]
"""
import argparse
import csv
import math
import time

import numpy as np
from scipy.special import erfc

from sim_signtest import make_folds, gene_stats, nested_oof, r2cols, bh_select, unit
from sim_mixture import strengths, population, classify


def make_operators(rng, n, n_ops, n_single):
    n_grp = n - n_single
    w = rng.dirichlet(np.full(n_ops, 1.5))
    sizes = np.maximum(3, np.round(w * n_grp).astype(int))
    while sizes.sum() > n_grp:
        sizes[np.argmax(sizes)] -= 1
    while sizes.sum() < n_grp:
        sizes[np.argmin(sizes)] += 1
    op = np.concatenate([np.full(s, k) for k, s in enumerate(sizes)] + [np.arange(n_ops, n_ops + n_single)])
    return rng.permutation(op)


def within_perm(rng, op):
    pm = np.arange(len(op))
    for k in np.unique(op):
        idx = np.flatnonzero(op == k)
        if len(idx) > 1:
            pm[idx] = idx[rng.permutation(len(idx))]
    return pm


def residualise(W, op, max_levels=12):
    """Residualise W on dummies of the `max_levels` largest operators (others: no column), as the
    published batch-corrected estimator does."""
    lev, cnt = np.unique(op, return_counts=True)
    keep = lev[np.argsort(-cnt)[:max_levels]]
    Dm = np.column_stack([np.ones(len(op))] + [(op == k).astype(float) for k in keep])
    beta, *_ = np.linalg.lstsq(Dm, W, rcond=None)
    return W - Dm @ beta


def simulate(rng, n, K, G, lam, R2, frac_u, hyp, delta, fX, n_ops, n_single, nu_W, nu_e, nu_P, nu_X):
    b = np.sqrt(R2 / lam)
    sig_u2 = frac_u * (1 - b ** 2)
    sig_e2 = 1 - b ** 2 - sig_u2
    s2, t2, q2, feas = strengths(hyp, lam, b, sig_u2, delta, fX)
    dil = 1.0 - nu_W
    s2d, t2d, q2d = s2 / dil, t2 / dil, q2 / dil       # undo the batch dilution of rho_W, rho_U, rho_E
    if max(s2d, t2d, q2d) > 0.95:
        feas = False
    s2d, t2d, q2d = min(s2d, 0.95), min(t2d, 0.95), min(q2d, 0.95)
    s, t, q = np.sqrt(s2d), np.sqrt(t2d), np.sqrt(q2d)
    op = make_operators(rng, n, n_ops, n_single)
    n_lev = op.max() + 1
    Z = rng.normal(size=(n_lev, K))
    W0 = rng.normal(size=(n, K))
    W = np.sqrt(nu_W) * Z[op] + np.sqrt(1 - nu_W) * W0
    Wn = W0  # unit-variance idiosyncratic block (W0 ~ N(0, I))
    M = np.empty((n, G))
    P = np.empty((n, G))
    for g in range(G):
        va = unit(rng.normal(size=K))
        wc = rng.normal(size=K)
        wc = unit(wc - (wc @ va) * va)
        ze = unit(rng.normal(size=K))
        Xb = s * (Wn @ va) + np.sqrt(1 - s2d) * rng.normal(size=n)
        X = np.sqrt(nu_X) * (Z[op] @ va) + np.sqrt(1 - nu_X) * Xb
        e0 = q * (Wn @ ze) + np.sqrt(1 - q2d) * rng.normal(size=n)
        e = np.sqrt((1 - lam) / lam) * (np.sqrt(nu_e) * rng.normal(size=n_lev)[op] + np.sqrt(1 - nu_e) * e0)
        u = np.sqrt(sig_u2) * (t * (Wn @ wc) + np.sqrt(1 - t2d) * rng.normal(size=n))
        eps = np.sqrt(sig_e2) * (np.sqrt(nu_P) * rng.normal(size=n_lev)[op] + np.sqrt(1 - nu_P) * rng.normal(size=n))
        M[:, g] = X + e
        P[:, g] = b * X + u + eps
    pop = population(hyp, lam, R2, b, sig_u2, s2, t2, q2)
    pop["feasible"] = feas
    return W, M, P, op, pop


def set_stats(obs, null, mask):
    out = {}
    for kk, lab in (("N_corr", "Nc"), ("N", "N"), ("D", "D"), ("S_cross", "S")):
        o = obs[kk][mask].sum()
        nl = null[kk][:, mask].sum(1)
        out[f"{lab}_z"] = float((o - nl.mean()) / (nl.std() + 1e-12))
        out[f"{lab}_nullmean"] = float(nl.mean())
    return out


def run_one(rng, n, K, G, lam, R2, frac_u, hyp, delta, fX, cfg, B, n_ops=12, n_single=10, lam_ridge=1.0):
    nu_W, nu_e, nu_P, nu_X = cfg
    W, M, P, op, pop = simulate(rng, n, K, G, lam, R2, frac_u, hyp, delta, fX, n_ops, n_single, nu_W, nu_e, nu_P, nu_X)
    Mz = (M - M.mean(0)) / M.std(0)
    Pz = (P - P.mean(0)) / P.std(0)
    folds = make_folds(n, rng)
    obs = gene_stats(W, Mz, Pz, folds, lam_ridge)
    # selection: published rule, unrestricted null, normal-tail p of the permutation z (see sim_mixture.py)
    inc_null = np.empty((B, G))
    for b in range(B):
        pm = rng.permutation(n)
        p0, p1 = nested_oof(Mz, W[pm], Pz, folds, lam_ridge)
        inc_null[b] = r2cols(Pz, p1) - r2cols(Pz, p0)
    z_inc = (obs["incr_p"] - inc_null.mean(0)) / (inc_null.std(0) + 1e-12)
    sel = bh_select(0.5 * erfc(z_inc / math.sqrt(2)), obs["incr_p"])
    keys = ["S_cross", "N", "N_corr", "D"]
    nullU = {k: np.empty((B, G)) for k in keys}
    nullS = {k: np.empty((B, G)) for k in keys}
    for b in range(B):
        st = gene_stats(W[rng.permutation(n)], Mz, Pz, folds, lam_ridge)
        for k in keys:
            nullU[k][b] = st[k]
        st = gene_stats(W[within_perm(rng, op)], Mz, Pz, folds, lam_ridge)
        for k in keys:
            nullS[k][b] = st[k]
    # operator-residualised W: observed theta only (the published batch-corrected estimator residualises
    # the PCs on the 12 largest operators)
    Wr = residualise(W, op, 12)
    obsR = gene_stats(Wr, Mz, Pz, folds, lam_ridge)
    res = dict(hyp=hyp, lam=lam, R2rna=R2, delta_target=delta, fX=fX if hyp == "mix" else np.nan,
               nu_W=nu_W, nu_e=nu_e, nu_P=nu_P, nu_X=nu_X, n=n, G=G, B=B, n_ops=n_ops, n_single=n_single,
               feasible=pop["feasible"], theta_pop=pop["theta_pop"], incr_p_pop=pop["incr_p_pop"],
               n_sel=int(sel.sum()), incr_p_mean_sel=float(obs["incr_p"][sel].mean()) if sel.any() else np.nan)
    if sel.sum() >= 5:
        res["theta_set_sel"] = float(obs["N"][sel].sum() / obs["D"][sel].sum())
        res["theta_set_sel_resid"] = float(obsR["N"][sel].sum() / obsR["D"][sel].sum())
        res["theta_set_all"] = float(obs["N"].sum() / obs["D"].sum())
        res["theta_set_all_resid"] = float(obsR["N"].sum() / obsR["D"].sum())
        for tag, nl in (("unres", nullU), ("strat", nullS)):
            for k2, v in set_stats(obs, nl, sel).items():
                res[f"{k2}_{tag}"] = v
            res[f"class_{tag}"] = classify(res[f"Nc_z_{tag}"], res["theta_set_sel"], res[f"D_z_{tag}"])
            # permutation-centred theta under each null
            Nn = nl["N"][:, sel].sum(1).mean()
            Dn = nl["D"][:, sel].sum(1).mean()
            res[f"theta_pc_{tag}"] = float((obs["N"][sel].sum() - Nn) / (obs["D"][sel].sum() - Dn))
        res["S_mean_sel"] = float(obs["S_cross"][sel].mean())
    return res


CONFIGS = {              # (nu_W, nu_e, nu_P, nu_X)
    "nobatch":   (0.0, 0.0, 0.0, 0.0),
    "W_only":    (0.3, 0.0, 0.0, 0.0),
    "rna":       (0.3, 0.3, 0.0, 0.0),
    "rna_strong": (0.3, 0.5, 0.0, 0.0),
    "prot":      (0.3, 0.0, 0.3, 0.0),
    "both":      (0.3, 0.3, 0.3, 0.0),
    "accrualX":  (0.3, 0.0, 0.0, 0.3),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--G", type=int, default=200)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--K", type=int, default=20)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261008)
    ap.add_argument("--configs", default=",".join(CONFIGS))
    ap.add_argument("--part", default="0/1")
    args = ap.parse_args()
    k, m = (int(x) for x in args.part.split("/"))
    hyps = [("post", 0.5, 0.10, np.nan), ("noise", 0.5, 0.10, np.nan), ("err", 0.5, 0.10, np.nan), ("mix", 0.5, 0.10, 0.3)]
    jobs = []
    i = 0
    for cname in args.configs.split(","):
        for hyp, lam, d, fX in hyps:
            for r in range(args.reps):
                if i % m == k:
                    jobs.append((i, cname, hyp, lam, d, fX, r))
                i += 1
    rows = []
    t0 = time.time()
    for j, (i, cname, hyp, lam, d, fX, r) in enumerate(jobs):
        rng = np.random.default_rng([args.seed, i])
        t1 = time.time()
        res = run_one(rng, args.n, args.K, args.G, lam, 0.3, 0.5, hyp, d, fX if np.isfinite(fX) else 0.5,
                      CONFIGS[cname], args.B)
        res.update(config=cname, rep=r, seconds=time.time() - t1)
        rows.append(res)
        print(f"[{j+1}/{len(jobs)}] {cname:10s} {hyp:5s} rep={r} feas={res['feasible']} sel={res['n_sel']:3d} "
              f"th={res.get('theta_set_sel', np.nan):+.2f} th_res={res.get('theta_set_sel_resid', np.nan):+.2f} (pop {res['theta_pop']:+.2f}) "
              f"Ncz unres={res.get('Nc_z_unres', np.nan):+.1f} strat={res.get('Nc_z_strat', np.nan):+.1f} | "
              f"Dz unres={res.get('D_z_unres', np.nan):+.1f} strat={res.get('D_z_strat', np.nan):+.1f} | "
              f"{res.get('class_unres','-')} / {res.get('class_strat','-')} [{res['seconds']:.0f}s]", flush=True)
        keys = sorted({kk for rr in rows for kk in rr})
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for rr in rows:
                w.writerow({kk: rr.get(kk, "") for kk in keys})
    print(f"wrote {args.out} ({len(rows)} rows, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
