"""
sim_signtest.py -- synthetic-data check of the measurement-error discriminators
(MorphoResidual, Limitations (v)).  NO project data are touched.

Generative model (per gene g, patients i = 1..n, morphology PCs W shared across genes):
    W      ~ N(0, I_K)                                   (the K morphology PCs, standardised)
    X_g    = s_g (W v_g) + sqrt(1 - s_g^2) xi_g          true expression, var 1, rho_W = s_g^2
    u_g    = sigma_u [ t_g (W w_g) + sqrt(1 - t_g^2) zeta_g ]   post-transcriptional, rho_U = t_g^2
    e_g    = sqrt((1-lam)/lam) [ q_g (W z_g) + sqrt(1 - q_g^2) eta_g ]  mRNA error, rho_E = q_g^2
    M_g    = X_g + e_g                                    measured mRNA, reliability lam
    P_g    = b X_g + u_g + eps_g,  var(P) = 1,  R2_rna = b^2 lam
Hypotheses: H_noise (t=q=0), H_post (s=q=0), H_err (s=t=0), H_mix (s,t>0, v_g _|_ w_g,
X-share f_X of the population increment).

Estimators (all 5-fold OOF, ridge alpha=1 on in-fold standardised features, exactly as the
published pipeline):
    incr_p   nested [M] vs [M, W] OOF R2 increment (published statistic)
    incr_m   mirror: nested [P] vs [P, W] for target M
    rp_hat   W-ONLY OOF ridge prediction of the OLS residual r_p = P - beta M
    rm_hat   W-ONLY OOF ridge prediction of the OLS residual r_m = M - delta P
    S_cross  sign(rho) * corr(rp_hat, r_m)            (F18 mirror sign test, cross form)
    S_pp     sign(rho) * corr(rp_hat, rm_hat)         (prediction-prediction form; mechanically biased)
    N        rho * cov(rp_hat, M~)                    (own-transcript alignment numerator)
    D        cov(rp_hat, r_p)                         (denominator; ~ increment)
    theta    N / D                                     population: f_X lam/(1-lam) | 0 | -1
Null: joint patient<->slide permutation of the rows of W (same permutation for every gene and
every target), B_null draws; selection of 'significant' genes by a separate B_sel-draw
permutation p-value on incr_p with BH FDR<0.05 and incr_p>0 (published rule).

Usage:  python sim_signtest.py --out results.csv [--genes 200 --n 100 --K 20 --B 200 --quick]
"""
import argparse
import json
import sys
import time

import numpy as np

# ----------------------------------------------------------------------------------------
# estimators
# ----------------------------------------------------------------------------------------

def make_folds(n, rng, k=5):
    perm = rng.permutation(n)
    return [(np.setdiff1d(np.arange(n), te), te) for te in np.array_split(perm, k)]


def _std_fit(A):
    mu = A.mean(0)
    sd = A.std(0)
    sd[sd == 0] = 1.0
    return mu, sd


def wonly_oof(W, Y, folds, lam=1.0):
    """OOF ridge prediction of every column of Y (n,G) from the shared block W (n,K)."""
    n, G = Y.shape
    K = W.shape[1]
    pred = np.empty((n, G))
    I = lam * np.eye(K)
    for tr, te in folds:
        mu, sd = _std_fit(W[tr])
        A = (W[tr] - mu) / sd
        Bm = (W[te] - mu) / sd
        ym = Y[tr].mean(0)
        H = Bm @ np.linalg.solve(A.T @ A + I, A.T)      # (n_te, n_tr)
        pred[te] = H @ (Y[tr] - ym) + ym
    return pred


def nested_oof(F, W, Y, folds, lam=1.0):
    """Per-gene scalar baseline F (n,G) vs [F_g, W]; returns (pred_base, pred_full), each (n,G)."""
    n, G = Y.shape
    K = W.shape[1]
    p0 = np.empty((n, G))
    p1 = np.empty((n, G))
    for tr, te in folds:
        mu, sd = _std_fit(W[tr])
        A = (W[tr] - mu) / sd
        Bm = (W[te] - mu) / sd
        fm, fs = _std_fit(F[tr])
        Ftr = (F[tr] - fm) / fs
        Fte = (F[te] - fm) / fs
        ym = Y[tr].mean(0)
        Yc = Y[tr] - ym
        fy = (Ftr * Yc).sum(0)                          # (G,)
        ff = (Ftr ** 2).sum(0) + lam                    # (G,)
        p0[te] = Fte * (fy / ff) + ym
        WtW = A.T @ A + lam * np.eye(K)
        fW = Ftr.T @ A                                  # (G,K)
        Gm = np.empty((G, K + 1, K + 1))
        Gm[:, 0, 0] = ff
        Gm[:, 0, 1:] = fW
        Gm[:, 1:, 0] = fW
        Gm[:, 1:, 1:] = WtW
        rhs = np.empty((G, K + 1))
        rhs[:, 0] = fy
        rhs[:, 1:] = (A.T @ Yc).T
        w = np.linalg.solve(Gm, rhs[..., None])[..., 0]  # (G,K+1)
        p1[te] = Fte * w[:, 0] + Bm @ w[:, 1:].T + ym
    return p0, p1


def r2cols(Y, P):
    return 1.0 - ((Y - P) ** 2).sum(0) / ((Y - Y.mean(0)) ** 2).sum(0)


def covcols(A, B):
    return ((A - A.mean(0)) * (B - B.mean(0))).mean(0)


def corrcols(A, B):
    return covcols(A, B) / np.sqrt(covcols(A, A) * covcols(B, B))


def gene_stats(W, Mz, Pz, folds, lam=1.0):
    """All per-gene statistics for one morphology block W (observed or permuted)."""
    rho = covcols(Mz, Pz)                       # Mz, Pz standardised -> beta = delta = rho
    rp = Pz - rho * Mz
    rm = Mz - rho * Pz
    rp_hat = wonly_oof(W, rp, folds, lam)
    rm_hat = wonly_oof(W, rm, folds, lam)
    sg = np.sign(rho)
    out = {}
    out["S_cross"] = sg * corrcols(rp_hat, rm)
    out["S_cross_m"] = sg * corrcols(rm_hat, rp)
    out["S_pp"] = sg * corrcols(rp_hat, rm_hat)
    out["N"] = rho * covcols(rp_hat, Mz)
    out["N_corr"] = rho * corrcols(rp_hat, Mz)   # scale-free form: its permutation null stays valid when W has a non-X signal
    out["D"] = covcols(rp_hat, rp)
    b0, b1 = nested_oof(Mz, W, Pz, folds, lam)
    out["incr_p"] = r2cols(Pz, b1) - r2cols(Pz, b0)
    m0, m1 = nested_oof(Pz, W, Mz, folds, lam)
    out["incr_m"] = r2cols(Mz, m1) - r2cols(Mz, m0)
    d_p = b1 - b0
    d_m = m1 - m0
    out["S_dd"] = sg * corrcols(d_p, d_m)       # F18 as literally specified (nested increments)
    out["S_drm"] = sg * corrcols(d_p, rm)       # nested increment vs mirror residual
    out["rho"] = rho
    return out


# ----------------------------------------------------------------------------------------
# generative model
# ----------------------------------------------------------------------------------------

def unit(v):
    return v / np.linalg.norm(v)


def simulate_cohort(rng, n, K, G, lam, R2rna, sig_u2, hyp, delta_target, fX=0.5):
    """Returns W (n,K), M (n,G), P (n,G), and the population parameters used."""
    b = np.sqrt(R2rna / lam)
    if sig_u2 is None or sig_u2 < 0:            # negative value = fraction of the non-mRNA variance
        frac = 0.5 if sig_u2 is None else -sig_u2
        sig_u2 = frac * (1.0 - b ** 2)
    sig_e2 = 1.0 - b ** 2 - sig_u2
    if sig_e2 <= 0:
        raise ValueError("sigma_eps^2 <= 0: lower sig_u2 or R2rna")
    W = rng.normal(size=(n, K))
    M = np.empty((n, G))
    P = np.empty((n, G))
    # population signal strengths
    s2 = t2 = q2 = 0.0
    feas = True
    if hyp == "noise":
        s2 = delta_target / (b ** 2 * (1 - lam) ** 2 + delta_target * lam)
        if s2 > 0.95:
            s2, feas = 0.95, False
    elif hyp == "post":
        t2 = delta_target / sig_u2
        if t2 > 0.95:
            t2, feas = 0.95, False
    elif hyp == "err":
        q2 = delta_target / (b ** 2 * lam * (1 - lam))
        if q2 > 0.95:
            q2, feas = 0.95, False
    elif hyp == "mix":
        s2 = fX * delta_target / (b ** 2 * (1 - lam) ** 2)
        t2 = (1 - fX) * delta_target / sig_u2
        if s2 > 0.95:
            s2, feas = 0.95, False
        if t2 > 0.95:
            t2, feas = 0.95, False
    else:
        raise ValueError(hyp)
    s, t, q = np.sqrt(s2), np.sqrt(t2), np.sqrt(q2)
    for g in range(G):
        va = unit(rng.normal(size=K))
        wc = rng.normal(size=K)
        wc = unit(wc - (wc @ va) * va)             # orthogonal to va (exact X-share bookkeeping)
        ze = unit(rng.normal(size=K))
        X = s * (W @ va) + np.sqrt(1 - s2) * rng.normal(size=n)
        u = np.sqrt(sig_u2) * (t * (W @ wc) + np.sqrt(1 - t2) * rng.normal(size=n))
        e = np.sqrt((1 - lam) / lam) * (q * (W @ ze) + np.sqrt(1 - q2) * rng.normal(size=n))
        M[:, g] = X + e
        P[:, g] = b * X + u + np.sqrt(sig_e2) * rng.normal(size=n)
    # population values
    pop = {}
    cp_x = b * (1 - lam) * s                    # X-channel covariance amplitude (along va)
    cp_u = sig_u2 * t / np.sqrt(sig_u2)         # u-channel: cov(W,u)=sigma_u t  -> enters r_p with coef 1
    cp_u = np.sqrt(sig_u2) * t
    cp_e = -b * lam * np.sqrt((1 - lam) / lam) * q
    pop["cov2"] = cp_x ** 2 + cp_u ** 2 + cp_e ** 2          # |C_p|^2 (orthogonal directions)
    # exact population increment (var P = 1): Delta_p = C_p' (I - h h')^{-1} C_p with
    # h = cov(W,M)/sd(M) = sqrt(lam) (s v_a + sqrt((1-lam)/lam) q z_e); Sherman-Morrison.
    # v_a, w_c orthogonal by construction; z_e is a random direction (treated as orthogonal).
    h2 = lam * (s2 + (1 - lam) / lam * q2)
    hC = np.sqrt(lam) * (s * cp_x + np.sqrt((1 - lam) / lam) * q * cp_e)
    pop["incr_p_pop"] = pop["cov2"] + hC ** 2 / (1 - h2)
    # theta = f_X lam/(1-lam) + f_E * (-1)   (shares of |C_p|^2)
    fx = cp_x ** 2 / pop["cov2"] if pop["cov2"] > 0 else 0.0
    fe = cp_e ** 2 / pop["cov2"] if pop["cov2"] > 0 else 0.0
    pop["theta_pop"] = fx * lam / (1 - lam) - fe
    pop["fX_pop"] = fx
    pop["feasible"] = feas
    pop["s2"], pop["t2"], pop["q2"], pop["b"], pop["sig_u2"] = s2, t2, q2, b, sig_u2
    # mirror increment, population, orthogonal channels, var(M)=1/lam:
    cm_x = s * (1 - R2rna / lam)                 # X-channel into r_m: coef (1 - delta b) = 1 - rho^2/lam
    cm_u = -b * np.sqrt(sig_u2) * t              # u-channel into r_m: coef -delta = -b (var P = 1)
    cm_e = np.sqrt((1 - lam) / lam) * q          # e-channel: coef 1
    # exact: Delta_m = lam * C_m' (I - g g')^{-1} C_m with g = cov(W,P) = b s v_a + sigma_u t w_c
    g2 = b ** 2 * s2 + sig_u2 * t2
    gC = b * s * cm_x + np.sqrt(sig_u2) * t * cm_u
    pop["incr_m_pop"] = lam * (cm_x ** 2 + cm_u ** 2 + cm_e ** 2 + gC ** 2 / (1 - g2))
    return W, M, P, pop


def bh_select(pvals, incr, alpha=0.05):
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * m / (np.arange(m) + 1)
    passed = np.zeros(m, bool)
    thr = np.where(ranked <= alpha)[0]
    if len(thr):
        passed[order[: thr.max() + 1]] = True
    return passed & (incr > 0)


# ----------------------------------------------------------------------------------------
# one scenario
# ----------------------------------------------------------------------------------------

def run_scenario(rng, n, K, G, lam, R2rna, sig_u2, hyp, delta, fX, B_sel, B_null, lam_ridge=1.0):
    W, M, P, pop = simulate_cohort(rng, n, K, G, lam, R2rna, sig_u2, hyp, delta, fX)
    Mz = (M - M.mean(0)) / M.std(0)
    Pz = (P - P.mean(0)) / P.std(0)
    folds = make_folds(n, rng)
    obs = gene_stats(W, Mz, Pz, folds, lam_ridge)
    # selection null (separate draws)
    cnt = np.zeros(G)
    for b in range(B_sel):
        pm = rng.permutation(n)
        p0, p1 = nested_oof(Mz, W[pm], Pz, folds, lam_ridge)
        inc = r2cols(Pz, p1) - r2cols(Pz, p0)
        cnt += inc >= obs["incr_p"]
    pv = (1 + cnt) / (1 + B_sel)
    sel = bh_select(pv, obs["incr_p"])
    # joint null for the discriminators
    keys = ["S_cross", "S_cross_m", "S_pp", "N", "N_corr", "D", "S_dd", "S_drm", "incr_m"]
    null = {k: np.empty((B_null, G)) for k in keys}
    for b in range(B_null):
        pm = rng.permutation(n)
        st = gene_stats(W[pm], Mz, Pz, folds, lam_ridge)
        for k in keys:
            null[k][b] = st[k]
    res = {"hyp": hyp, "lam": lam, "delta_target": delta, "fX": fX if hyp == "mix" else np.nan,
           "R2rna": R2rna, "sig_u2": pop["sig_u2"], "n": n, "K": K, "G": G, "B_sel": B_sel, "B_null": B_null,
           "feasible": pop["feasible"], "s2": pop["s2"], "t2": pop["t2"], "q2": pop["q2"],
           "incr_p_pop": pop["incr_p_pop"], "incr_m_pop": pop["incr_m_pop"], "theta_pop": pop["theta_pop"],
           "fX_pop": pop["fX_pop"], "n_sel": int(sel.sum()), "sel_rate": sel.mean(),
           "incr_p_mean_all": obs["incr_p"].mean(), "incr_p_mean_sel": obs["incr_p"][sel].mean() if sel.any() else np.nan,
           "incr_m_mean_all": obs["incr_m"].mean(), "incr_m_mean_sel": obs["incr_m"][sel].mean() if sel.any() else np.nan,
           "ratio_m_over_p_sel": (obs["incr_m"][sel].mean() / obs["incr_p"][sel].mean()) if sel.any() else np.nan,
           "rho_mean": obs["rho"].mean()}

    def pergene(k, sign_expected):
        o = obs[k]
        nl = null[k]
        # two-sided per-gene p
        ge = (nl >= o).mean(0)
        le = (nl <= o).mean(0)
        p2 = np.minimum(1.0, 2 * np.minimum((1 + ge * B_null) / (1 + B_null), (1 + le * B_null) / (1 + B_null)))
        z = (o - nl.mean(0)) / (nl.std(0) + 1e-12)
        d = {f"{k}_mean_all": o.mean(), f"{k}_null_mean_all": nl.mean(),
             f"{k}_mean_sel": o[sel].mean() if sel.any() else np.nan,
             f"{k}_null_mean_sel": nl[:, sel].mean() if sel.any() else np.nan,
             f"{k}_pg_power_sel": ((p2 < 0.05) & (np.sign(z) == sign_expected))[sel].mean() if sel.any() and sign_expected != 0 else np.nan,
             f"{k}_pg_rej_sel": (p2 < 0.05)[sel].mean() if sel.any() else np.nan,
             f"{k}_pg_rej_all": (p2 < 0.05).mean(),
             f"{k}_pg_fracpos_sel": (z > 0)[sel].mean() if sel.any() else np.nan}
        # aggregate over the selected set: median vs null-of-median (shared permutations)
        if sel.sum() >= 5:
            om = np.median(o[sel])
            nm = np.median(nl[:, sel], axis=1)
            d[f"{k}_agg_med_sel"] = om
            d[f"{k}_agg_null_mean"] = nm.mean()
            d[f"{k}_agg_null_sd"] = nm.std()
            d[f"{k}_agg_z"] = (om - nm.mean()) / (nm.std() + 1e-12)
            d[f"{k}_agg_p2"] = min(1.0, 2 * min((1 + (nm >= om).sum()) / (1 + B_null), (1 + (nm <= om).sum()) / (1 + B_null)))
        return d

    exp_sign = {"noise": +1, "post": -1, "err": -1, "mix": int(np.sign(res["theta_pop"] - 0))}
    for k in ["S_cross", "S_cross_m", "S_pp", "S_dd", "S_drm"]:
        # expected sign of the sign tests: + iff theta_pop > rho^2/(1-rho^2) (DERIVATION.md section 2.2).
        # NOTE: the sim_grid_v2_*.csv files were produced with the earlier (incorrect) threshold rho^2 here;
        # it only affects the S_*_pg_power_sel columns, which are not used in DERIVATION.md section 7.
        thr = R2rna / (1.0 - R2rna)
        res.update(pergene(k, +1 if res["theta_pop"] > thr else -1))
    nsign = +1 if hyp in ("noise",) or (hyp == "mix" and fX > 0) else (-1 if hyp == "err" else 0)
    res.update(pergene("N", nsign))
    res.update(pergene("N_corr", nsign))
    res.update(pergene("D", +1))
    # theta: per-gene median and ratio-of-sums on the selected set, with patient-free null for N
    th = obs["N"] / obs["D"]
    res["theta_pg_median_all"] = np.median(th)
    res["theta_pg_median_sel"] = np.median(th[sel]) if sel.any() else np.nan
    res["theta_set_all"] = obs["N"].sum() / obs["D"].sum()
    res["theta_set_sel"] = obs["N"][sel].sum() / obs["D"][sel].sum() if sel.any() else np.nan
    # lambda_hat_noise = theta/(1+theta)
    ts = res["theta_set_sel"]
    res["lamhat_noise_sel"] = ts / (1 + ts) if np.isfinite(ts) and ts > -1 else np.nan
    # aggregate N over selected genes vs joint null
    if sel.sum() >= 5:
        for kk, tag in (("N", "N_set"), ("N_corr", "Nc_set")):
            Nsum = obs[kk][sel].sum()
            Nn = null[kk][:, sel].sum(1)
            res[f"{tag}_z"] = (Nsum - Nn.mean()) / (Nn.std() + 1e-12)
            res[f"{tag}_p2"] = min(1.0, 2 * min((1 + (Nn >= Nsum).sum()) / (1 + B_null), (1 + (Nn <= Nsum).sum()) / (1 + B_null)))
        # all-genes version of the scale-free aggregate (no selection)
        Nsum = obs["N_corr"].sum()
        Nn = null["N_corr"].sum(1)
        res["Nc_all_z"] = (Nsum - Nn.mean()) / (Nn.std() + 1e-12)
        # bootstrap-free CI proxy for theta_set: gene-level jackknife (genes share W -> approximate)
        Ns, Ds = obs["N"][sel], obs["D"][sel]
        jk = np.array([(Ns.sum() - Ns[i]) / (Ds.sum() - Ds[i]) for i in range(len(Ns))])
        res["theta_set_sel_jkse"] = np.sqrt((len(jk) - 1) / len(jk) * ((jk - jk.mean()) ** 2).sum())
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="sim_signtest_results.csv")
    ap.add_argument("--genes", type=int, default=200)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--K", type=int, default=20)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261008)
    ap.add_argument("--R2rna", type=str, default="0.1,0.3", help="comma list of own-mRNA R2 values")
    ap.add_argument("--frac_u", type=float, default=0.5, help="sigma_u^2 as a fraction of the non-mRNA protein variance")
    ap.add_argument("--lams", type=str, default="0.5,0.6,0.7,0.8,0.9")
    ap.add_argument("--deltas", type=str, default="0.10,0.20", help="population increments (var P = 1)")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    lams = [float(x) for x in args.lams.split(",")]
    deltas = [float(x) for x in args.deltas.split(",")]
    r2s = [float(x) for x in args.R2rna.split(",")]
    scen = []
    for R2 in r2s:
        for lam in lams:
            for d in deltas:
                scen.append(("noise", lam, d, np.nan, R2))
                scen.append(("post", lam, d, np.nan, R2))
                scen.append(("err", lam, d, np.nan, R2))
                for fX in (0.1, 0.3, 0.5):
                    scen.append(("mix", lam, d, fX, R2))
    if args.quick:
        scen = [s for s in scen if s[1] in (0.6, 0.9) and s[2] == 0.20 and (s[0] != "mix" or s[3] == 0.3)]
    import csv
    rows = []
    t0 = time.time()
    for i, (hyp, lam, d, fX, R2) in enumerate(scen):
        t1 = time.time()
        r = run_scenario(rng, args.n, args.K, args.genes, lam, R2, -args.frac_u, hyp, d,
                         fX if np.isfinite(fX) else 0.5, args.B, args.B)
        r["seconds"] = time.time() - t1
        rows.append(r)
        print(f"[{i+1}/{len(scen)}] {hyp:5s} R2={R2:.1f} lam={lam:.2f} delta={d:.2f} fX={fX} feas={r['feasible']} "
              f"sel={r['n_sel']:3d} incr_p={r['incr_p_mean_all']:+.4f}(pop {r['incr_p_pop']:.4f}) "
              f"theta_set_sel={r['theta_set_sel']:+.2f}(pop {r['theta_pop']:+.2f}) "
              f"S_cross_sel={r['S_cross_mean_sel']:+.3f} null={r['S_cross_null_mean_sel']:+.3f} "
              f"aggz={r.get('S_cross_agg_z', np.nan):+.1f} Nz={r.get('N_set_z', np.nan):+.1f} Ncz={r.get('Nc_set_z', np.nan):+.1f} NcRej={r.get('N_corr_pg_rej_sel', np.nan):.2f} "
              f"[{r['seconds']:.0f}s]", flush=True)
    keys = sorted({k for r in rows for k in r})
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})
    print(f"wrote {args.out}  ({time.time()-t0:.0f}s total)")


if __name__ == "__main__":
    main()
