"""vt_finite.py -- independent finite-sample checks (own implementation, synthetic only).
5-fold OOF, ridge alpha=1 on in-fold standardised features, as the published pipeline.
E2: H_post n=100 K=20 G=200 with selection: per-gene rejection of N (cov) vs N_c (corr), set-level z.
E3: H_post large selected set (~1500 genes): z_N with in-sample vs in-fold (OOF) baseline residuals.
E4: H_noise winner's curse vs selection intensity (theta_sel/theta_pop).
E5: H_post with MNAR protein truncation at n=100 (per-gene fit set), set-level z_N and theta.
"""
import sys, time
import numpy as np
from scipy.stats import norm

def folds_of(rng, n, k=5):
    perm = rng.permutation(n)
    return [(np.setdiff1d(np.arange(n), te), te) for te in np.array_split(perm, k)]

def stdz(A, mu, sd):
    return (A - mu) / sd

def fit_std(A):
    mu = A.mean(0); sd = A.std(0); sd = np.where(sd == 0, 1.0, sd); return mu, sd

def wonly_hat(W, folds, alpha=1.0):
    """Per fold: (te, H) with H (n_te, n_tr) mapping centred training targets to test predictions."""
    K = W.shape[1]; out = []
    for tr, te in folds:
        mu, sd = fit_std(W[tr]); A = stdz(W[tr], mu, sd); B = stdz(W[te], mu, sd)
        H = B @ np.linalg.solve(A.T @ A + alpha * np.eye(K), A.T)
        out.append((tr, te, H))
    return out

def oof_wonly(hats, Y):
    pred = np.empty_like(Y)
    for tr, te, H in hats:
        ym = Y[tr].mean(0); pred[te] = H @ (Y[tr] - ym) + ym
    return pred

def nested_incr(W, M, P, folds, alpha=1.0):
    """Published statistic: OOF R2(P|M,W) - OOF R2(P|M); ridge alpha on in-fold standardised [M_g, W]."""
    n, G = P.shape; K = W.shape[1]
    p0 = np.empty((n, G)); p1 = np.empty((n, G))
    for tr, te in folds:
        mu, sd = fit_std(W[tr]); A = stdz(W[tr], mu, sd); B = stdz(W[te], mu, sd)
        mm, ms = fit_std(M[tr]); Mtr = stdz(M[tr], mm, ms); Mte = stdz(M[te], mm, ms)
        ym = P[tr].mean(0); Yc = P[tr] - ym
        fy = (Mtr * Yc).sum(0); ff = (Mtr ** 2).sum(0) + alpha
        p0[te] = Mte * (fy / ff) + ym
        AtA = A.T @ A + alpha * np.eye(K)
        fW = Mtr.T @ A                                   # (G,K)
        Gm = np.empty((G, K + 1, K + 1))
        Gm[:, 0, 0] = ff; Gm[:, 0, 1:] = fW; Gm[:, 1:, 0] = fW; Gm[:, 1:, 1:] = AtA
        rhs = np.empty((G, K + 1)); rhs[:, 0] = fy; rhs[:, 1:] = (A.T @ Yc).T
        w = np.linalg.solve(Gm, rhs[..., None])[..., 0]
        p1[te] = Mte * w[:, 0] + B @ w[:, 1:].T + ym
    sst = ((P - P.mean(0)) ** 2).sum(0)
    return (((P - p0) ** 2).sum(0) - ((P - p1) ** 2).sum(0)) / sst

def cov0(A, B): return ((A - A.mean(0)) * (B - B.mean(0))).mean(0)
def corr0(A, B): return cov0(A, B) / np.sqrt(cov0(A, A) * cov0(B, B))

def discrim(W, M, P, folds, mode="insample", alpha=1.0):
    """N (cov form), N_c (corr form), D, S_cross for all genes; mode = insample | infold residuals."""
    hats = wonly_hat(W, folds, alpha)
    n, G = P.shape
    if mode == "insample":
        Mz = (M - M.mean(0)) / M.std(0); Pz = (P - P.mean(0)) / P.std(0)
        rho = cov0(Mz, Pz); rp = Pz - rho * Mz; rm = Mz - rho * Pz
        rp_hat = oof_wonly(hats, rp)
    else:
        rp = np.empty((n, G)); rm = np.empty((n, G)); rp_hat = np.empty((n, G)); rhos = np.empty((len(hats), G))
        Mz = np.empty((n, G))
        for f, (tr, te, H) in enumerate(hats):
            mm, ms = fit_std(M[tr]); pm, ps = fit_std(P[tr])
            Mtr = stdz(M[tr], mm, ms); Ptr = stdz(P[tr], pm, ps)
            Mte = stdz(M[te], mm, ms); Pte = stdz(P[te], pm, ps)
            r = cov0(Mtr, Ptr); rhos[f] = r
            rtr = Ptr - r * Mtr; rp[te] = Pte - r * Mte; rm[te] = Mte - r * Pte; Mz[te] = Mte
            rp_hat[te] = H @ (rtr - rtr.mean(0)) + rtr.mean(0)
        rho = rhos.mean(0)
    sg = np.sign(rho)
    return dict(N=rho * cov0(rp_hat, Mz), Nc=rho * corr0(rp_hat, Mz), D=cov0(rp_hat, rp),
                S=sg * corr0(rp_hat, rm), rho=rho)

def simulate(rng, n, K, G, lam, R2, hyp, delta, frac_u=0.5, fX=0.5):
    """Same generative model as the document (W ~ N(0,I_K); orthogonal X/u directions)."""
    b = np.sqrt(R2 / lam); su2 = frac_u * (1 - b ** 2); se2 = 1 - b ** 2 - su2
    W = rng.normal(size=(n, K)); M = np.empty((n, G)); P = np.empty((n, G))
    s2 = t2 = q2 = 0.0
    if hyp == "noise": s2 = min(0.95, delta / (b ** 2 * (1 - lam) ** 2 + delta * lam))
    elif hyp == "post": t2 = min(0.95, delta / su2)
    elif hyp == "err": q2 = min(0.95, delta / (b ** 2 * lam * (1 - lam)))
    elif hyp == "mix":
        s2 = min(0.95, fX * delta / (b ** 2 * (1 - lam) ** 2)); t2 = min(0.95, (1 - fX) * delta / su2)
    elif hyp == "null": pass
    s, t, q = np.sqrt(s2), np.sqrt(t2), np.sqrt(q2)
    for g in range(G):
        va = rng.normal(size=K); va /= np.linalg.norm(va)
        wc = rng.normal(size=K); wc -= (wc @ va) * va; wc /= np.linalg.norm(wc)
        ze = rng.normal(size=K); ze /= np.linalg.norm(ze)
        X = s * (W @ va) + np.sqrt(1 - s2) * rng.normal(size=n)
        u = np.sqrt(su2) * (t * (W @ wc) + np.sqrt(1 - t2) * rng.normal(size=n))
        e = np.sqrt((1 - lam) / lam) * (q * (W @ ze) + np.sqrt(1 - q2) * rng.normal(size=n))
        M[:, g] = X + e; P[:, g] = b * X + u + np.sqrt(se2) * rng.normal(size=n)
    cpx2 = (b * (1 - lam) * s) ** 2; cpu2 = su2 * t2; cpe2 = (b * lam * np.sqrt((1 - lam) / lam) * q) ** 2
    tot = cpx2 + cpu2 + cpe2
    theta_pop = (cpx2 / tot) * lam / (1 - lam) - cpe2 / tot if tot > 0 else np.nan
    return W, M, P, theta_pop

def bh(p, alpha=0.05):
    m = len(p); o = np.argsort(p); r = p[o] * m / (np.arange(m) + 1)
    ok = np.zeros(m, bool); hit = np.where(r <= alpha)[0]
    if len(hit): ok[o[:hit.max() + 1]] = True
    return ok

def select(rng, W, M, P, folds, B_sel, use_normal_tail=False):
    obs = nested_incr(W, M, P, folds)
    nulls = np.empty((B_sel, M.shape[1]))
    for b in range(B_sel):
        nulls[b] = nested_incr(W[rng.permutation(len(W))], M, P, folds)
    if use_normal_tail:
        z = (obs - nulls.mean(0)) / nulls.std(0); p = norm.sf(z)
    else:
        p = (1 + (nulls >= obs).sum(0)) / (1 + B_sel)
    return bh(p) & (obs > 0), obs, nulls

def joint_null(rng, W, M, P, folds, B, mode):
    keys = ["N", "Nc", "D", "S"]; out = {k: np.empty((B, M.shape[1])) for k in keys}
    for b in range(B):
        st = discrim(W[rng.permutation(len(W))], M, P, folds, mode)
        for k in keys: out[k][b] = st[k]
    return out

def setz(obs, null, sel):
    o = obs[sel].sum(); nn = null[:, sel].sum(1); return (o - nn.mean()) / nn.std()

def pergene_rej(obs, null):
    B = null.shape[0]
    ge = (null >= obs).sum(0); le = (null <= obs).sum(0)
    p2 = np.minimum(1, 2 * np.minimum((1 + ge) / (1 + B), (1 + le) / (1 + B)))
    return p2 < 0.05

def E2(seeds=6):
    print("\n=== E2: H_post, n=100, K=20, G=200, R2=0.3, lam=0.7, selection + joint null (B=200) ===")
    print("seed delta n_sel | per-gene rej (sel): N_cov N_corr | (all genes) N_cov N_corr | z_N(Nc) z_D theta_set theta_pc | z_N infold")
    rows = []
    for delta in (0.1, 0.2):
        for sd in range(seeds):
            rng = np.random.default_rng(1000 + sd)
            W, M, P, _ = simulate(rng, 100, 20, 200, 0.7, 0.3, "post", delta)
            folds = folds_of(rng, 100)
            sel, inc, _ = select(rng, W, M, P, folds, 200)
            obs = discrim(W, M, P, folds, "insample"); null = joint_null(rng, W, M, P, folds, 200, "insample")
            obs2 = discrim(W, M, P, folds, "infold"); null2 = joint_null(rng, W, M, P, folds, 200, "infold")
            rN = pergene_rej(obs["N"], null["N"]); rNc = pergene_rej(obs["Nc"], null["Nc"])
            if sel.sum() >= 5:
                zN = setz(obs["Nc"], null["Nc"], sel); zD = setz(obs["D"], null["D"], sel)
                th = obs["N"][sel].sum() / obs["D"][sel].sum()
                thpc = (obs["N"][sel].sum() - null["N"][:, sel].sum(1).mean()) / (obs["D"][sel].sum() - null["D"][:, sel].sum(1).mean())
                zN2 = setz(obs2["Nc"], null2["Nc"], sel)
            else:
                zN = zD = th = thpc = zN2 = np.nan
            rows.append((delta, sel.sum(), rN[sel].mean() if sel.any() else np.nan, rNc[sel].mean() if sel.any() else np.nan, rN.mean(), rNc.mean(), zN, zD, th, thpc, zN2))
            print(f"{sd:4d} {delta:.1f} {sel.sum():5d} | {rows[-1][2]:.3f} {rows[-1][3]:.3f} | {rN.mean():.3f} {rNc.mean():.3f} | {zN:+.2f} {zD:+.1f} {th:+.3f} {thpc:+.3f} | {zN2:+.2f}", flush=True)
    r = np.array(rows, float)
    for delta in (0.1, 0.2):
        m = r[r[:, 0] == delta]
        print(f"delta={delta}: mean n_sel={np.nanmean(m[:,1]):.0f}; per-gene rej on selected: cov {np.nanmean(m[:,2]):.3f} corr {np.nanmean(m[:,3]):.3f}; "
              f"all genes: cov {np.nanmean(m[:,4]):.3f} corr {np.nanmean(m[:,5]):.3f}; z_N mean {np.nanmean(m[:,6]):+.2f} max {np.nanmax(m[:,6]):+.2f}; "
              f"theta_set mean {np.nanmean(m[:,8]):+.3f} max|.| {np.nanmax(np.abs(m[:,8])):.3f}; z_N infold mean {np.nanmean(m[:,10]):+.2f} max {np.nanmax(m[:,10]):+.2f}")

def E3(seeds=4, G=3000):
    print(f"\n=== E3: H_post large set: n=100, K=20, G={G} all signal, R2=0.3, lam=0.6, delta=0.2; normal-tail BH selection ===")
    print("seed n_sel | z_N insample (sel, all) | z_N infold (sel, all) | theta_pc insample infold | per-gene mean excess Nc (sel) insample infold")
    for sd in range(seeds):
        t0 = time.time(); rng = np.random.default_rng(2000 + sd)
        W, M, P, _ = simulate(rng, 100, 20, G, 0.6, 0.3, "post", 0.2)
        folds = folds_of(rng, 100)
        sel, inc, _ = select(rng, W, M, P, folds, 200, use_normal_tail=True)
        res = {}
        for mode in ("insample", "infold"):
            obs = discrim(W, M, P, folds, mode); null = joint_null(rng, W, M, P, folds, 200, mode)
            allg = np.ones(G, bool)
            thpc = (obs["N"][sel].sum() - null["N"][:, sel].sum(1).mean()) / (obs["D"][sel].sum() - null["D"][:, sel].sum(1).mean())
            res[mode] = (setz(obs["Nc"], null["Nc"], sel), setz(obs["Nc"], null["Nc"], allg), thpc,
                         (obs["Nc"][sel] - null["Nc"][:, sel].mean(0)).mean())
        a, b_ = res["insample"], res["infold"]
        print(f"{sd:4d} {sel.sum():5d} | {a[0]:+.2f} {a[1]:+.2f} | {b_[0]:+.2f} {b_[1]:+.2f} | {a[2]:+.4f} {b_[2]:+.4f} | {a[3]:+.5f} {b_[3]:+.5f}  [{time.time()-t0:.0f}s]", flush=True)

def E4(seeds=3):
    print("\n=== E4: H_noise winner's curse vs selection intensity (n=100, K=20, R2=0.3, lam=0.5 -> theta_pop=1) ===")
    print("(a) G=2000 all-signal genes at delta=0.1, ranked by published increment; theta_set (raw) of the top fraction")
    fr = (1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005)
    for sd in range(seeds):
        rng = np.random.default_rng(3000 + sd)
        W, M, P, thp = simulate(rng, 100, 20, 2000, 0.5, 0.3, "noise", 0.1)
        folds = folds_of(rng, 100)
        inc = nested_incr(W, M, P, folds); obs = discrim(W, M, P, folds, "insample")
        order = np.argsort(-inc)
        out = []
        for f in fr:
            k = max(5, int(round(f * 2000))); idx = order[:k]
            out.append(obs["N"][idx].sum() / obs["D"][idx].sum())
        print(f"  seed {sd}: " + " ".join(f"top{f*100:g}%={o/thp:.2f}" for f, o in zip(fr, out)))
    print("(b) mixture G=2000 with 20% signal genes; BH(normal-tail p) selection; theta_set,sel/theta_pop vs delta")
    for delta in (0.08, 0.1, 0.15, 0.2):
        for sd in range(2):
            rng = np.random.default_rng(4000 + sd)
            W, M, P, thp = simulate(rng, 100, 20, 1600, 0.5, 0.3, "null", 0.0)
            # share the same W: regenerate signal genes on W
            rng2 = np.random.default_rng(6000 + sd)
            _, M2, P2, thp = simulate_onW(rng2, W, 400, 0.5, 0.3, "noise", delta)
            Mall = np.hstack([M, M2]); Pall = np.hstack([P, P2]); sig = np.r_[np.zeros(1600, bool), np.ones(400, bool)]
            folds = folds_of(rng, 100)
            sel, inc, _ = select(rng, W, Mall, Pall, folds, 200, use_normal_tail=True)
            obs = discrim(W, Mall, Pall, folds, "insample")
            if sel.sum() >= 5:
                th = obs["N"][sel].sum() / obs["D"][sel].sum()
                print(f"  delta={delta} seed={sd}: n_sel={sel.sum()} (false {(sel & ~sig).sum()}), signal selected {sel[sig].mean()*100:.1f}%, theta_set/theta_pop={th/thp:.2f}")
            else:
                print(f"  delta={delta} seed={sd}: n_sel={sel.sum()} < 5")

def simulate_onW(rng, W, G, lam, R2, hyp, delta, frac_u=0.5):
    n, K = W.shape
    b = np.sqrt(R2 / lam); su2 = frac_u * (1 - b ** 2); se2 = 1 - b ** 2 - su2
    s2 = min(0.95, delta / (b ** 2 * (1 - lam) ** 2 + delta * lam)) if hyp == "noise" else 0.0
    t2 = min(0.95, delta / su2) if hyp == "post" else 0.0
    s, t = np.sqrt(s2), np.sqrt(t2)
    M = np.empty((n, G)); P = np.empty((n, G))
    for g in range(G):
        va = rng.normal(size=K); va /= np.linalg.norm(va)
        wc = rng.normal(size=K); wc -= (wc @ va) * va; wc /= np.linalg.norm(wc)
        X = s * (W @ va) + np.sqrt(1 - s2) * rng.normal(size=n)
        u = np.sqrt(su2) * (t * (W @ wc) + np.sqrt(1 - t2) * rng.normal(size=n))
        e = np.sqrt((1 - lam) / lam) * rng.normal(size=n)
        M[:, g] = X + e; P[:, g] = b * X + u + np.sqrt(se2) * rng.normal(size=n)
    cpx2 = (b * (1 - lam) * s) ** 2; cpu2 = su2 * t2; tot = cpx2 + cpu2
    return W, M, P, (cpx2 / tot) * lam / (1 - lam) if tot > 0 else np.nan

def E5(seeds=2):
    print("\n=== E5: H_post with MNAR protein (per-gene: drop patients whose P is in the bottom q), n=130 so ~100 remain at q=0.25 ===")
    print("q seed n_sel | z_N(Nc) theta_set theta_pc | z_N on all genes")
    for q in (0.1, 0.25):
        for sd in range(seeds):
            rng = np.random.default_rng(7000 + sd)
            n = 100 if q == 0 else int(round(100 / (1 - q)))
            W, M, P, _ = simulate(rng, n, 20, 200, 0.7, 0.3, "post", 0.2)
            # per gene, different fit sets -> loop genes (small G), common folds on the kept patients differ per gene;
            # to keep the shared-W joint null exact we use a common patient permutation and per-gene subsetting.
            G = 200
            keep = P >= np.quantile(P, q, axis=0) if q > 0 else np.ones_like(P, bool)
            perms = [rng.permutation(n) for _ in range(120)]
            fold_seed = rng.integers(1 << 30)
            def stats_for(Wx):
                N = np.empty(G); Nc = np.empty(G); D = np.empty(G); inc = np.empty(G)
                for g in range(G):
                    k = keep[:, g]; Wg = Wx[k]; Mg = M[k, g:g + 1]; Pg = P[k, g:g + 1]
                    fl = folds_of(np.random.default_rng(fold_seed + g), k.sum())
                    st = discrim(Wg, Mg, Pg, fl, "insample")
                    N[g], Nc[g], D[g] = st["N"][0], st["Nc"][0], st["D"][0]
                    inc[g] = nested_incr(Wg, Mg, Pg, fl)[0]
                return N, Nc, D, inc
            N, Nc, D, inc = stats_for(W)
            nullN = np.empty((120, G)); nullNc = np.empty((120, G)); nullD = np.empty((120, G)); nullI = np.empty((120, G))
            for b, pm in enumerate(perms):
                nullN[b], nullNc[b], nullD[b], nullI[b] = stats_for(W[pm])
            z = (inc - nullI.mean(0)) / nullI.std(0); sel = bh(norm.sf(z)) & (inc > 0)
            allg = np.ones(G, bool)
            if sel.sum() >= 5:
                zN = setz(Nc, nullNc, sel); th = N[sel].sum() / D[sel].sum()
                thpc = (N[sel].sum() - nullN[:, sel].sum(1).mean()) / (D[sel].sum() - nullD[:, sel].sum(1).mean())
            else:
                zN = th = thpc = np.nan
            print(f"{q:.2f} {sd:4d} {sel.sum():5d} | {zN:+.2f} {th:+.3f} {thpc:+.3f} | {setz(Nc, nullNc, allg):+.2f}", flush=True)

if __name__ == "__main__":
    which = sys.argv[1:] or ["E2", "E3", "E4", "E5"]
    for w in which:
        globals()[w]()
