"""
finite_sample.py -- independent finite-sample Monte Carlo (n = 100, K = 20, 5-fold OOF ridge
alpha = 1 on in-fold standardised features) of the numerical claims in DERIVATION.md.  Own
implementation (not imported from ../sim_signtest.py).  Synthetic data only.

Checks
  1. OOF capacity offset of the nested increment, by hypothesis, at the SAME population Delta_p
     and R2_rna  (claim: -0.15..-0.20; and 'H_noise penalised more: -0.33 vs -0.18')
  2. null means of S_pp, S_cross, D and the nested-increment forms (claims: -0.46..-0.31, +0.01, ...)
  3. N_c (correlation form) per-gene rejection rate under H_post vs the permutation null
       (a) independent genes            (claim: nominal 5%)
       (b) co-expressed family sharing one morphology axis   (my adversarial case)
     and the set-level z of sum N_c under H_post for both
  4. theta_set with and without selection vs population (winner's curse)
"""
import sys
import time
import numpy as np

rng = np.random.default_rng(99)


def folds_of(n, k=5):
    perm = rng.permutation(n)
    return [(np.setdiff1d(np.arange(n), te), te) for te in np.array_split(perm, k)]


def std_fit(A):
    mu = A.mean(0); sd = A.std(0); sd[sd == 0] = 1
    return mu, sd


def oof_ridge(W, Y, folds, alpha=1.0):
    """W-only OOF ridge prediction of each column of Y."""
    n, G = Y.shape; K = W.shape[1]
    out = np.empty((n, G))
    for tr, te in folds:
        mu, sd = std_fit(W[tr]); A = (W[tr] - mu) / sd; B = (W[te] - mu) / sd
        ym = Y[tr].mean(0)
        H = B @ np.linalg.solve(A.T @ A + alpha * np.eye(K), A.T)
        out[te] = H @ (Y[tr] - ym) + ym
    return out


def oof_nested(F, W, Y, folds, alpha=1.0):
    """baseline [F_g] vs full [F_g, W] per gene; returns OOF preds (base, full)."""
    n, G = Y.shape; K = W.shape[1]
    p0 = np.empty((n, G)); p1 = np.empty((n, G))
    for tr, te in folds:
        mu, sd = std_fit(W[tr]); A = (W[tr] - mu) / sd; B = (W[te] - mu) / sd
        fm, fs = std_fit(F[tr]); Ftr = (F[tr] - fm) / fs; Fte = (F[te] - fm) / fs
        ym = Y[tr].mean(0); Yc = Y[tr] - ym
        for g in range(G):
            x0 = Ftr[:, g:g + 1]
            c0 = np.linalg.solve(x0.T @ x0 + alpha, x0.T @ Yc[:, g])
            p0[te, g] = Fte[:, g] * c0[0] + ym[g]
            X1 = np.column_stack([x0, A])
            c1 = np.linalg.solve(X1.T @ X1 + alpha * np.eye(K + 1), X1.T @ Yc[:, g])
            p1[te, g] = np.column_stack([Fte[:, g:g + 1], B]) @ c1 + ym[g]
    return p0, p1


def r2(Y, Pr):
    return 1 - ((Y - Pr) ** 2).sum(0) / ((Y - Y.mean(0)) ** 2).sum(0)


def ccov(A, B):
    return ((A - A.mean(0)) * (B - B.mean(0))).mean(0)


def ccorr(A, B):
    return ccov(A, B) / np.sqrt(ccov(A, A) * ccov(B, B))


def gene_stats(W, Mz, Pz, folds):
    rho = ccov(Mz, Pz)
    rp = Pz - rho * Mz; rm = Mz - rho * Pz
    rph = oof_ridge(W, rp, folds); rmh = oof_ridge(W, rm, folds)
    sg = np.sign(rho)
    b0, b1 = oof_nested(Mz, W, Pz, folds)
    m0, m1 = oof_nested(Pz, W, Mz, folds)
    dp = b1 - b0; dm = m1 - m0
    return dict(N=rho * ccov(rph, Mz), Nc=rho * ccorr(rph, Mz), D=ccov(rph, rp),
                S_cross=sg * ccorr(rph, rm), S_pp=sg * ccorr(rph, rmh), S_dd=sg * ccorr(dp, dm),
                S_drm=sg * ccorr(dp, rm), incr=r2(Pz, b1) - r2(Pz, b0), incr_m=r2(Mz, m1) - r2(Mz, m0), rho=rho)


def cohort(n, K, G, hyp, lam, R2, Dp, shared_axis=False, coexpr=0.0):
    """One synthetic cohort. sigma_P^2 = 1; sigma_u^2 = half the non-mRNA protein variance.
    shared_axis: all genes' u-channel points along the same W direction; coexpr: shared X factor
    across genes with that loading (co-expression)."""
    b = np.sqrt(R2 / lam); su2 = 0.5 * (1 - b * b); sg2 = 1 - b * b - su2
    W = rng.normal(size=(n, K))
    M = np.empty((n, G)); P = np.empty((n, G))
    if hyp == "noise":
        s2 = Dp / (b * b * (1 - lam) ** 2 + Dp * lam); t2 = q2 = 0.0
        assert s2 <= 1, "infeasible"
    elif hyp == "post":
        t2 = Dp / su2; s2 = q2 = 0.0
    elif hyp == "err":
        q2 = Dp / (b * b * lam * (1 - lam) + Dp * (1 - lam)); s2 = t2 = 0.0   # exact inverse incl. denominator
    s, t, q = np.sqrt(s2), np.sqrt(t2), np.sqrt(q2)
    Xshared = rng.normal(size=n); ushared = rng.normal(size=n)
    wc_shared = rng.normal(size=K); wc_shared /= np.linalg.norm(wc_shared)
    for g in range(G):
        va = rng.normal(size=K); va /= np.linalg.norm(va)
        wc = wc_shared if shared_axis else rng.normal(size=K)
        wc = wc - (wc @ va) * va; wc /= np.linalg.norm(wc)
        ze = rng.normal(size=K); ze /= np.linalg.norm(ze)
        xi = np.sqrt(1 - coexpr ** 2) * rng.normal(size=n) + coexpr * Xshared
        X = s * (W @ va) + np.sqrt(1 - s2) * xi
        zeta = np.sqrt(1 - coexpr ** 2) * rng.normal(size=n) + coexpr * ushared
        u = np.sqrt(su2) * (t * (W @ wc) + np.sqrt(1 - t2) * zeta)
        e = np.sqrt((1 - lam) / lam) * (q * (W @ ze) + np.sqrt(1 - q2) * rng.normal(size=n))
        M[:, g] = X + e
        P[:, g] = b * X + u + np.sqrt(sg2) * rng.normal(size=n)
    theta_pop = {"noise": lam / (1 - lam), "post": 0.0, "err": -1.0}[hyp]
    return W, M, P, theta_pop


def bh(p, alpha=0.05):
    m = len(p); o = np.argsort(p); r = p[o] * m / (np.arange(m) + 1)
    ok = np.zeros(m, bool); idx = np.where(r <= alpha)[0]
    if len(idx): ok[o[:idx.max() + 1]] = True
    return ok


def run_scenario(hyp, lam, R2, Dp, n=100, K=20, G=300, B=200, shared_axis=False, coexpr=0.0, tag=""):
    t0 = time.time()
    W, M, P, th_pop = cohort(n, K, G, hyp, lam, R2, Dp, shared_axis, coexpr)
    Mz = (M - M.mean(0)) / M.std(0); Pz = (P - P.mean(0)) / P.std(0)
    folds = folds_of(n)
    obs = gene_stats(W, Mz, Pz, folds)
    keys = list(obs)
    null = {k: np.empty((B, G)) for k in keys}
    for bb in range(B):
        st = gene_stats(W[rng.permutation(n)], Mz, Pz, folds)
        for k in keys: null[k][bb] = st[k]
    # selection: permutation p on incr (same draws; fine for a rate check)
    pv = (1 + (null["incr"] >= obs["incr"]).sum(0)) / (1 + B)
    sel = bh(pv) & (obs["incr"] > 0)
    R2full = R2 + Dp
    off_pred = -(K / (0.8 * n - K - 1)) * (1 - R2full)
    print(f"\n[{hyp} lam={lam} R2={R2} Dp={Dp} {tag}]  G={G} B={B} n_sel={sel.sum()}  ({time.time()-t0:.0f}s)")
    print(f"  nested incr: obs mean={obs['incr'].mean():+.3f}  pop={Dp:.3f}  offset={obs['incr'].mean()-Dp:+.3f}  "
          f"pred K/(ntr-K-1)(1-R2full)={off_pred:+.3f}  null mean={null['incr'].mean():+.3f}")
    print(f"  null means: S_pp={null['S_pp'].mean():+.3f}  S_cross={null['S_cross'].mean():+.3f}  D={null['D'].mean():+.4f}  "
          f"S_drm={null['S_drm'].mean():+.3f}  S_dd={null['S_dd'].mean():+.3f}  Nc={null['Nc'].mean():+.4f}")
    print(f"  observed means: S_cross={obs['S_cross'].mean():+.3f} S_pp={obs['S_pp'].mean():+.3f} S_drm={obs['S_drm'].mean():+.3f} "
          f"S_dd={obs['S_dd'].mean():+.3f}  incr_m={obs['incr_m'].mean():+.3f}")
    for k in ("N", "Nc"):
        o = obs[k]; nl = null[k]
        ge = (1 + (nl >= o).sum(0)) / (1 + B); le = (1 + (nl <= o).sum(0)) / (1 + B)
        p2 = np.minimum(1, 2 * np.minimum(ge, le))
        print(f"  {k}: per-gene two-sided rejection at 5%: all genes={np.mean(p2 < 0.05):.3f}  selected={np.mean(p2[sel]<0.05) if sel.any() else np.nan:.3f}")
        Ssum = o.sum(); Sn = nl.sum(1)
        zs = (Ssum - Sn.mean()) / Sn.std()
        if sel.sum() >= 5:
            zsel = (o[sel].sum() - nl[:, sel].sum(1).mean()) / nl[:, sel].sum(1).std()
        else:
            zsel = np.nan
        print(f"  {k}: set-level z (all genes)={zs:+.2f}   (selected)={zsel:+.2f}")
    th_all = obs["N"].sum() / obs["D"].sum()
    th_sel = obs["N"][sel].sum() / obs["D"][sel].sum() if sel.any() else np.nan
    print(f"  theta_pop={th_pop:+.3f}  theta_set(all)={th_all:+.3f}  theta_set(selected)={th_sel:+.3f}  "
          f"lamhat(all)={th_all/(1+th_all):+.3f} lamhat(sel)={th_sel/(1+th_sel) if np.isfinite(th_sel) else np.nan:+.3f}")
    return obs, null, sel


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "offset"):
        # 1-2: same population Delta_p and R2 under the three hypotheses
        for hyp in ("noise", "post", "err"):
            run_scenario(hyp, 0.6, 0.3, 0.10, G=300, B=100)
        for hyp in ("noise", "post"):
            run_scenario(hyp, 0.5, 0.3, 0.20, G=300, B=100)
    if which in ("all", "post"):
        # 3: N_c null under H_post: independent genes vs co-expressed family with a shared morphology axis
        run_scenario("post", 0.8, 0.3, 0.15, G=400, B=200, tag="independent genes")
        run_scenario("post", 0.8, 0.3, 0.15, G=400, B=200, shared_axis=True, coexpr=0.0, tag="shared u-axis, no co-expression")
        run_scenario("post", 0.8, 0.3, 0.15, G=400, B=200, shared_axis=True, coexpr=0.7, tag="shared u-axis + co-expression 0.7")
        run_scenario("post", 0.8, 0.3, 0.15, G=400, B=200, shared_axis=False, coexpr=0.7, tag="co-expression 0.7, independent axes")
    if which in ("all", "repeat"):
        # replicate the set-level z under H_post over independent cohorts (type-I at set level)
        zs_ind, zs_fam = [], []
        for rep in range(12):
            o, nl, sel = run_scenario("post", 0.8, 0.3, 0.15, G=60, B=150, tag=f"rep{rep} independent")
            Sn = nl["Nc"].sum(1); zs_ind.append((o["Nc"].sum() - Sn.mean()) / Sn.std())
            o, nl, sel = run_scenario("post", 0.8, 0.3, 0.15, G=60, B=150, shared_axis=True, coexpr=0.7, tag=f"rep{rep} family")
            Sn = nl["Nc"].sum(1); zs_fam.append((o["Nc"].sum() - Sn.mean()) / Sn.std())
        zs_ind, zs_fam = np.array(zs_ind), np.array(zs_fam)
        print("\nSET-LEVEL Nc z under H_post across cohorts:")
        print(f"  independent genes: z = {np.round(zs_ind,2)}  sd={zs_ind.std():.2f}  |z|>1.96 in {np.mean(np.abs(zs_ind)>1.96):.2f}")
        print(f"  co-expressed family: z = {np.round(zs_fam,2)}  sd={zs_fam.std():.2f}  |z|>1.96 in {np.mean(np.abs(zs_fam)>1.96):.2f}")
