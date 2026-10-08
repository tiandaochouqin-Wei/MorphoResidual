"""
sim_mixture.py -- null + signal gene MIXTURE at realistic gene counts, with replicate seeds
(synthetic only; NO project data).  Companion of sim_signtest.py: same generative model, same
estimators (5-fold OOF ridge alpha=1, in-fold standardisation), same joint permutation null.

What is new relative to sim_signtest.py
  * only a fraction pi1 of the G genes carries the hypothesis signal; the rest are pure null
    (s = t = q = 0), so selection operates against a real null background and the selected set
    contains false discoveries, as in the published analysis;
  * replicate seeds per scenario, so that the DISTRIBUTION of the set-level statistics
    (theta_set, lambda_hat, Nc_set z, S_cross z) is estimated and false-classification rates of
    the reading rules can be stated;
  * three estimates of theta are compared: theta_set_sel (selected set, winner's curse),
    theta_set_all (all tested genes, selection-free but diluted by null genes) and the
    oracle theta_set_signal (all signal genes, no selection), plus the permutation-centred
    theta_set_sel_pc = (sum N - null mean) / (sum D - null mean);
  * selection p-values: the published rule is an upper-tail permutation p with B draws and BH
    FDR < 0.05.  With G = 2,000 genes and B = 200 the permutation p cannot resolve below
    1/(B+1) = 0.005, so BH would pass either nothing or everything (the paper notes the same
    granularity at B = 1,000).  To keep the selection pressure realistic without B = 1,000 per
    draw, the per-gene p used for BH is the normal-tail p of the permutation z-score
    (obs - null mean)/null sd; the permutation p is also stored.  Selection = BH(p) < 0.05 and
    incr_p > 0.

Usage:
  python sim_mixture.py --scen main --out sim_mixture_main.csv [--reps 5 --G 2000 --pi1 0.2 --part 0/2]
  python sim_mixture.py --scen size --out sim_mixture_size.csv
"""
import argparse
import csv
import math
import sys
import time

import numpy as np
from scipy.special import erfc

from sim_signtest import make_folds, gene_stats, nested_oof, r2cols, bh_select, unit


# ----------------------------------------------------------------------------------------
# generative model with a signal mask
# ----------------------------------------------------------------------------------------

def strengths(hyp, lam, b, sig_u2, delta, fX):
    """Population signal strengths (s2, t2, q2) that deliver population increment `delta`
    under hypothesis `hyp`; identical bookkeeping to sim_signtest.simulate_cohort."""
    s2 = t2 = q2 = 0.0
    feas = True
    if hyp == "noise":
        s2 = delta / (b ** 2 * (1 - lam) ** 2 + delta * lam)
    elif hyp == "post":
        t2 = delta / sig_u2
    elif hyp == "err":
        q2 = delta / (b ** 2 * lam * (1 - lam))
    elif hyp == "mix":
        s2 = fX * delta / (b ** 2 * (1 - lam) ** 2)
        t2 = (1 - fX) * delta / sig_u2
    elif hyp == "null":
        pass
    else:
        raise ValueError(hyp)
    feas = max(s2, t2, q2) <= 0.95
    s2, t2, q2 = min(s2, 0.95), min(t2, 0.95), min(q2, 0.95)
    return s2, t2, q2, feas


def population(hyp, lam, R2rna, b, sig_u2, s2, t2, q2):
    s, t, q = np.sqrt(s2), np.sqrt(t2), np.sqrt(q2)
    cp_x = b * (1 - lam) * s
    cp_u = np.sqrt(sig_u2) * t
    cp_e = -b * lam * np.sqrt((1 - lam) / lam) * q
    cov2 = cp_x ** 2 + cp_u ** 2 + cp_e ** 2
    h2 = lam * (s2 + (1 - lam) / lam * q2)
    hC = np.sqrt(lam) * (s * cp_x + np.sqrt((1 - lam) / lam) * q * cp_e)
    incr = cov2 + (hC ** 2 / (1 - h2) if h2 < 1 else 0.0)
    fx = cp_x ** 2 / cov2 if cov2 > 0 else 0.0
    fe = cp_e ** 2 / cov2 if cov2 > 0 else 0.0
    return dict(incr_p_pop=incr, theta_pop=fx * lam / (1 - lam) - fe, fX_pop=fx)


def simulate_mixture(rng, n, K, G, pi1, lam, R2rna, frac_u, hyp, delta, fX):
    b = np.sqrt(R2rna / lam)
    sig_u2 = frac_u * (1.0 - b ** 2)
    sig_e2 = 1.0 - b ** 2 - sig_u2
    if sig_e2 <= 0:
        raise ValueError("sigma_eps^2 <= 0")
    s2, t2, q2, feas = strengths(hyp, lam, b, sig_u2, delta, fX)
    s, t, q = np.sqrt(s2), np.sqrt(t2), np.sqrt(q2)
    n_sig = int(round(pi1 * G))
    signal = np.zeros(G, bool)
    signal[rng.choice(G, n_sig, replace=False)] = True
    W = rng.normal(size=(n, K))
    M = np.empty((n, G))
    P = np.empty((n, G))
    for g in range(G):
        sg, tg, qg = (s, t, q) if signal[g] else (0.0, 0.0, 0.0)
        va = unit(rng.normal(size=K))
        wc = rng.normal(size=K)
        wc = unit(wc - (wc @ va) * va)
        ze = unit(rng.normal(size=K))
        X = sg * (W @ va) + np.sqrt(1 - sg ** 2) * rng.normal(size=n)
        u = np.sqrt(sig_u2) * (tg * (W @ wc) + np.sqrt(1 - tg ** 2) * rng.normal(size=n))
        e = np.sqrt((1 - lam) / lam) * (qg * (W @ ze) + np.sqrt(1 - qg ** 2) * rng.normal(size=n))
        M[:, g] = X + e
        P[:, g] = b * X + u + np.sqrt(sig_e2) * rng.normal(size=n)
    pop = population(hyp, lam, R2rna, b, sig_u2, s2, t2, q2)
    pop.update(feasible=feas, s2=s2, t2=t2, q2=q2, b=b, sig_u2=sig_u2)
    return W, M, P, signal, pop


# ----------------------------------------------------------------------------------------
# reading rule (section 4.4 / 7.4 of DERIVATION.md), applied to set-level statistics
# ----------------------------------------------------------------------------------------

def classify(zN, theta, zD, z_lo=-3.0, z_hi=2.0, theta_noise=0.5):
    if not np.isfinite(zD) or zD < 2.0:
        return "R6_nosignal"
    if zN <= z_lo:
        return "R4_err"
    if zN < z_hi:
        return "R1_post"
    if theta < theta_noise:
        return "R3_mixed"
    return "R2_noise_compatible"


# ----------------------------------------------------------------------------------------
# one replicate
# ----------------------------------------------------------------------------------------

def run_one(rng, n, K, G, pi1, lam, R2rna, frac_u, hyp, delta, fX, B_sel, B_null, lam_ridge=1.0):
    W, M, P, signal, pop = simulate_mixture(rng, n, K, G, pi1, lam, R2rna, frac_u, hyp, delta, fX)
    Mz = (M - M.mean(0)) / M.std(0)
    Pz = (P - P.mean(0)) / P.std(0)
    folds = make_folds(n, rng)
    obs = gene_stats(W, Mz, Pz, folds, lam_ridge)
    # ---- selection null (separate draws) ----
    cnt = np.zeros(G)
    inc_null = np.empty((B_sel, G))
    for b in range(B_sel):
        pm = rng.permutation(n)
        p0, p1 = nested_oof(Mz, W[pm], Pz, folds, lam_ridge)
        inc_null[b] = r2cols(Pz, p1) - r2cols(Pz, p0)
    cnt = (inc_null >= obs["incr_p"]).sum(0)
    p_perm = (1 + cnt) / (1 + B_sel)
    z_inc = (obs["incr_p"] - inc_null.mean(0)) / (inc_null.std(0) + 1e-12)
    p_norm = 0.5 * erfc(z_inc / math.sqrt(2))
    sel = bh_select(p_norm, obs["incr_p"])
    sel_perm = bh_select(p_perm, obs["incr_p"])
    # ---- joint null for the discriminators ----
    keys = ["S_cross", "N", "N_corr", "D"]
    null = {k: np.empty((B_null, G)) for k in keys}
    for b in range(B_null):
        pm = rng.permutation(n)
        st = gene_stats(W[pm], Mz, Pz, folds, lam_ridge)
        for k in keys:
            null[k][b] = st[k]
    res = dict(hyp=hyp, lam=lam, R2rna=R2rna, delta_target=delta, fX=fX if hyp == "mix" else np.nan,
               pi1=pi1, n=n, K=K, G=G, B_sel=B_sel, B_null=B_null, n_signal=int(signal.sum()),
               feasible=pop["feasible"], incr_p_pop=pop["incr_p_pop"], theta_pop=pop["theta_pop"],
               fX_pop=pop["fX_pop"], n_sel=int(sel.sum()), n_sel_perm_p=int(sel_perm.sum()),
               n_sel_null=int((sel & ~signal).sum()), sel_rate_signal=float(sel[signal].mean()),
               fdr_realised=float((sel & ~signal).sum() / max(1, sel.sum())),
               incr_p_mean_sel=float(obs["incr_p"][sel].mean()) if sel.any() else np.nan,
               incr_p_mean_signal=float(obs["incr_p"][signal].mean()),
               rho2_mean_sel=float((obs["rho"][sel] ** 2).mean()) if sel.any() else np.nan)
    N, D, Nc = obs["N"], obs["D"], obs["N_corr"]

    def theta_of(mask):
        return float(N[mask].sum() / D[mask].sum()) if mask.any() and D[mask].sum() != 0 else np.nan

    res["theta_set_signal"] = theta_of(signal)                 # oracle: all signal genes, no selection
    res["theta_set_sel"] = theta_of(sel)
    res["theta_set_sel_signalonly"] = theta_of(sel & signal)
    res["theta_set_all"] = theta_of(np.ones(G, bool))
    ts = res["theta_set_sel"]
    res["lamhat_noise_sel"] = ts / (1 + ts) if np.isfinite(ts) and ts > -1 else np.nan
    # set-level statistics on the selected set vs the joint null; the same on all genes
    for tag, mask in (("sel", sel), ("all", np.ones(G, bool)), ("signal", signal)):
        if mask.sum() < 5:
            continue
        for kk, lab in (("N", "N"), ("N_corr", "Nc"), ("D", "D"), ("S_cross", "S")):
            o = obs[kk][mask].sum()
            nl = null[kk][:, mask].sum(1)
            res[f"{lab}_{tag}_z"] = float((o - nl.mean()) / (nl.std() + 1e-12))
            res[f"{lab}_{tag}_nullmean"] = float(nl.mean())
            res[f"{lab}_{tag}_nullsd"] = float(nl.std())
        # permutation-centred theta
        Nn = null["N"][:, mask].sum(1).mean()
        Dn = null["D"][:, mask].sum(1).mean()
        res[f"theta_pc_{tag}"] = float((N[mask].sum() - Nn) / (D[mask].sum() - Dn))
        res[f"S_mean_{tag}"] = float(obs["S_cross"][mask].mean())
        res[f"S_nullmean_{tag}"] = float(null["S_cross"][:, mask].mean())
    # per-gene two-sided rejection at 5% of the covariance and correlation forms (selected genes)
    for kk, lab in (("N", "N"), ("N_corr", "Nc")):
        o, nl = obs[kk], null[kk]
        ge = (nl >= o).mean(0)
        le = (nl <= o).mean(0)
        p2 = np.minimum(1.0, 2 * np.minimum((1 + ge * B_null) / (1 + B_null), (1 + le * B_null) / (1 + B_null)))
        res[f"{lab}_pg_rej_sel"] = float((p2 < 0.05)[sel].mean()) if sel.any() else np.nan
        res[f"{lab}_pg_rej_null"] = float((p2 < 0.05)[~signal].mean())
        hi = (o > np.quantile(nl, 0.975, axis=0))
        lo = (o < np.quantile(nl, 0.025, axis=0))
        res[f"{lab}_pg_fracabove_sel"] = float(hi[sel].mean()) if sel.any() else np.nan
        res[f"{lab}_pg_fracbelow_sel"] = float(lo[sel].mean()) if sel.any() else np.nan
    # gene-jackknife SE of theta_set_sel (approximate: genes share W)
    if sel.sum() >= 5:
        Ns, Ds = N[sel], D[sel]
        jk = np.array([(Ns.sum() - Ns[i]) / (Ds.sum() - Ds[i]) for i in range(len(Ns))])
        res["theta_set_sel_jkse"] = float(np.sqrt((len(jk) - 1) / len(jk) * ((jk - jk.mean()) ** 2).sum()))
    res["class_sel"] = classify(res.get("Nc_sel_z", np.nan), res["theta_set_sel"], res.get("D_sel_z", np.nan))
    res["class_all"] = classify(res.get("Nc_all_z", np.nan), res["theta_set_all"], res.get("D_all_z", np.nan))
    return res


def scenarios(which):
    S = []
    if which == "main":
        for lam in (0.6, 0.8):
            for d in (0.10, 0.20):
                S.append(("post", lam, d, np.nan, 0.3))
        for lam, d in ((0.5, 0.10), (0.6, 0.10), (0.7, 0.10), (0.5, 0.20)):
            S.append(("noise", lam, d, np.nan, 0.3))
        for lam in (0.5, 0.6):
            S.append(("err", lam, 0.10, np.nan, 0.3))
        for fX in (0.1, 0.3, 0.5):
            for lam in (0.6, 0.8):
                S.append(("mix", lam, 0.20, fX, 0.3))
        S.append(("post", 0.6, 0.20, np.nan, 0.1))
        S.append(("mix", 0.6, 0.20, 0.1, 0.1))
        S.append(("mix", 0.6, 0.20, 0.3, 0.1))
    elif which == "size":
        for hyp, lam, d, fX in (("post", 0.6, 0.20, np.nan), ("noise", 0.5, 0.20, np.nan), ("mix", 0.6, 0.20, 0.3)):
            for G, pi1 in ((2000, 0.05), (2000, 0.5), (5000, 0.5)):
                S.append((hyp, lam, d, fX, 0.3, G, pi1))
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scen", default="main", choices=["main", "size"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--G", type=int, default=2000)
    ap.add_argument("--pi1", type=float, default=0.2)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--K", type=int, default=20)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--frac_u", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=20261008)
    ap.add_argument("--part", default="0/1", help="k/m: run every m-th scenario starting at k")
    ap.add_argument("--resume", action="store_true", help="keep rows already in --out, skip their (scenario, rep)")
    args = ap.parse_args()
    k, m = (int(x) for x in args.part.split("/"))
    S = scenarios(args.scen)
    rows = []
    done = set()
    if args.resume:
        import os
        if os.path.exists(args.out):
            with open(args.out, newline="", encoding="utf-8") as f:
                for rr in csv.DictReader(f):
                    rows.append(rr)
                    done.add((int(rr["scenario"]), int(rr["rep"])))
            print(f"resume: {len(rows)} rows already in {args.out}", flush=True)
    t0 = time.time()
    jobs = []
    for i, sc in enumerate(S):
        if i % m != k:
            continue
        for r in range(args.reps):
            if (i, r) not in done:
                jobs.append((i, r, sc))
    for j, (i, r, sc) in enumerate(jobs):
        hyp, lam, d, fX, R2 = sc[:5]
        G = sc[5] if len(sc) > 5 else args.G
        pi1 = sc[6] if len(sc) > 6 else args.pi1
        rng = np.random.default_rng([args.seed, i, r])
        t1 = time.time()
        res = run_one(rng, args.n, args.K, G, pi1, lam, R2, args.frac_u, hyp, d,
                      fX if np.isfinite(fX) else 0.5, args.B, args.B)
        res["rep"] = r
        res["scenario"] = i
        res["seconds"] = time.time() - t1
        rows.append(res)
        print(f"[{j+1}/{len(jobs)}] {hyp:5s} R2={R2} lam={lam} d={d} fX={fX} G={G} pi1={pi1} rep={r} "
              f"feas={res['feasible']} nsig={res['n_signal']} sel={res['n_sel']} (null {res['n_sel_null']}, permBH {res['n_sel_perm_p']}) "
              f"th_sel={res['theta_set_sel']:+.3f} th_all={res['theta_set_all']:+.3f} th_sig={res['theta_set_signal']:+.3f} "
              f"(pop {res['theta_pop']:+.2f}) Ncz_sel={res.get('Nc_sel_z', np.nan):+.1f} Ncz_all={res.get('Nc_all_z', np.nan):+.1f} "
              f"Dz={res.get('D_sel_z', np.nan):+.1f} S_z={res.get('S_sel_z', np.nan):+.1f} -> {res['class_sel']} [{res['seconds']:.0f}s]",
              flush=True)
        keys = sorted({kk for rr in rows for kk in rr})
        tmp = args.out + ".tmp"
        with open(tmp, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for rr in rows:
                w.writerow({kk: rr.get(kk, "") for kk in keys})
        import os
        os.replace(tmp, args.out)
    print(f"wrote {args.out} ({len(rows)} rows, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
