"""vt_final.py -- synthetic checks for the post-review DERIVATION.md (final reading rules).
Synthetic data only.  Reuses the independent harness vt_finite.py (5-fold OOF, ridge alpha=1,
in-fold standardisation, joint W-row permutation null).

F1  sign(rho_hat)-weighted N_c,set under H_post for a ~1,500-gene selected set (vs rho_hat-weighted)
F2  closed-form corrections: Delta_p^err denominator; Delta_m/Delta_p factors (population limit)
F3  false-classification rates of the FINAL numeric rules (grids + mixtures), sign-weighted statistic
F4  winner's-curse kappa vs selection strength at a second (lambda, Delta) point
F5  MNAR protein: per-gene fit sets with 0 / 5 / 10 % truncation; rule outcome; corr(W-pred M, M) diagnostic
"""
import sys, time, os
import numpy as np
from scipy.stats import norm

VT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verify-theory")
sys.path.insert(0, VT)
from vt_finite import (simulate, simulate_onW, folds_of, discrim, joint_null, select, setz,
                       nested_incr, wonly_hat, oof_wonly, cov0, corr0, bh, fit_std, stdz)


# ---------------------------------------------------------------- final rules
def kappa_of(strength):
    """winner's-curse floor by observable selection strength n_sel / G_tested (review curve)."""
    if strength >= 0.20: return 0.60
    if strength >= 0.10: return 0.50
    if strength >= 0.05: return 0.45
    if strength >= 0.005: return 0.40
    return 0.35


def classify(zD, zN, thpc, n_sel, lam_lo, strength):
    if n_sel < 5: return "R6"
    if zD < 3: return "R6"
    if zN <= -3:
        return "R4" if thpc <= -0.3 else "R4w"           # R4w = weak anti-alignment (undecided)
    band = 0.05 if n_sel >= 100 else 0.10
    if zN < 2:
        return "R1" if abs(thpc) <= band else "R7"       # R7 = undecidable (D noise)
    thmin = kappa_of(strength) * lam_lo / (1 - lam_lo)
    tag = "R2" if thpc >= thmin else "R3"
    return tag + ("w" if zN < 3 else "")


def set_stats(obs, null, sel):
    """sign-weighted N_c (= N_c/|rho|, rho is W-free so the null divides by the same |rho|)."""
    a = np.abs(obs["rho"])
    zNs = setz(obs["Nc"] / a, null["Nc"] / a, sel)
    zNr = setz(obs["Nc"], null["Nc"], sel)
    zD = setz(obs["D"], null["D"], sel)
    thpc = (obs["N"][sel].sum() - null["N"][:, sel].sum(1).mean()) / (obs["D"][sel].sum() - null["D"][:, sel].sum(1).mean())
    return zNs, zNr, zD, thpc


# ---------------------------------------------------------------- F1
def F1(seeds=4, B=150):
    print("\n=== F1: H_post, n=100, K=20, G=2500 all-signal, R2=0.3, Delta=0.2, normal-tail BH -> n_sel ~ 1,500 ===")
    print("lam seed n_sel | z[sign*corr] z[rho*corr] | theta_pc | time")
    out = {}
    for lam in (0.6, 0.8):
        zs = []
        for sd in range(seeds):
            t0 = time.time(); rng = np.random.default_rng(9100 + sd)
            W, M, P, _ = simulate(rng, 100, 20, 2500, lam, 0.3, "post", 0.2)
            folds = folds_of(rng, 100)
            sel, inc, _ = select(rng, W, M, P, folds, B, use_normal_tail=True)
            obs = discrim(W, M, P, folds, "insample"); null = joint_null(rng, W, M, P, folds, B, "insample")
            zNs, zNr, zD, thpc = set_stats(obs, null, sel)
            zs.append(zNs)
            print(f"{lam:.1f} {sd:4d} {sel.sum():5d} | {zNs:+.2f} {zNr:+.2f} | {thpc:+.4f} | {time.time()-t0:.0f}s", flush=True)
        out[lam] = zs
        print(f"  lam={lam}: z[sign*corr] mean {np.mean(zs):+.2f}, max {np.max(zs):+.2f}, n(z>=2)={sum(z >= 2 for z in zs)}")
    return out


# ---------------------------------------------------------------- F2
def F2(n=400000):
    print("\n=== F2: population-limit closed forms (one signal direction, K=20, R2_rna=0.3, var P=1, rho_W=rho_E=0.6) ===")
    rng = np.random.default_rng(55)
    R2 = 0.3; rW = 0.6
    def ols_r2(Y, Xs):
        Xc = Xs - Xs.mean(0); Yc = Y - Y.mean()
        beta = np.linalg.lstsq(Xc, Yc, rcond=None)[0]
        return 1 - ((Yc - Xc @ beta) ** 2).mean() / (Yc ** 2).mean()
    for lam in (0.5, 0.7, 0.9):
        b = np.sqrt(R2 / lam); su2 = 0.5 * (1 - b ** 2); se2 = 1 - b ** 2 - su2; sig_e = np.sqrt((1 - lam) / lam)
        W = rng.normal(size=(n, 20)); Z = W[:, 0]
        for hyp in ("noise", "err"):
            if hyp == "noise":
                X = np.sqrt(rW) * Z + np.sqrt(1 - rW) * rng.normal(size=n); e = sig_e * rng.normal(size=n)
            else:
                X = rng.normal(size=n); e = sig_e * (np.sqrt(rW) * Z + np.sqrt(1 - rW) * rng.normal(size=n))
            u = np.sqrt(su2) * rng.normal(size=n); M = X + e; P = b * X + u + np.sqrt(se2) * rng.normal(size=n)
            Dp = ols_r2(P, np.column_stack([M, W])) - ols_r2(P, M[:, None])
            Dm = ols_r2(M, np.column_stack([P, W])) - ols_r2(M, P[:, None])
            if hyp == "noise":
                Dp_pred = R2 * (1 - lam) ** 2 * rW / (lam * (1 - rW * lam))
                ratio_old = (lam - R2) ** 2 / ((1 - lam) ** 2 * R2)
                ratio_new = ratio_old * (1 - lam * rW) / (1 - rW * R2 / lam)
            else:
                Dp_old = rW * R2 * (1 - lam)
                Dp_pred = rW * R2 * (1 - lam) / (1 - (1 - lam) * rW)
                ratio_old = 1 / R2
                ratio_new = (1 - (1 - lam) * rW) / R2
            extra = f" (old text Dp={Dp_old:.4f})" if hyp == "err" else ""
            print(f"{hyp:5s} lam={lam}: Dp sim={Dp:.4f} pred={Dp_pred:.4f}{extra} | Dm/Dp sim={Dm/Dp:.3f} old={ratio_old:.3f} corrected={ratio_new:.3f}")


# ---------------------------------------------------------------- F3
def run_cohort(rng, n, G, lam, R2, hyp, delta, fX, B, normal_tail=False):
    W, M, P, thp = simulate(rng, n, 20, G, lam, R2, hyp, delta, fX=fX)
    folds = folds_of(rng, n)
    sel, inc, _ = select(rng, W, M, P, folds, B, use_normal_tail=normal_tail)
    obs = discrim(W, M, P, folds, "insample"); null = joint_null(rng, W, M, P, folds, B, "insample")
    return W, M, P, thp, folds, sel, obs, null


def F3(seeds=3, B=150, part="ab"):
    print("\n=== F3: final rules on synthetic sets (sign-weighted z_N, z>=2 everywhere, kappa by n_sel/G_tested, lam_lo = true lam) ===")
    cells = []
    for lam in (0.6, 0.8):
        for delta in (0.1, 0.2):
            cells.append(("post", lam, delta, 0.5))
    for lam in (0.5, 0.6, 0.7): cells.append(("noise", lam, 0.1, 0.5))
    cells.append(("noise", 0.5, 0.2, 0.5))
    for lam in (0.5, 0.6): cells.append(("err", lam, 0.1, 0.5))
    for fX in (0.1, 0.3, 0.5):
        for lam in (0.6, 0.8): cells.append(("mix", lam, 0.2, fX))
    cells.append(("mix", 0.6, 0.1, 0.1)); cells.append(("mix", 0.8, 0.1, 0.05))
    print("hyp lam delta fX | seed n_sel strength | z_D z_N(sign) theta_pc theta_pop | class")
    tally = {}
    if "a" not in part: cells = []
    for ci, (hyp, lam, delta, fX) in enumerate(cells):
        for sd in range(seeds):
            rng = np.random.default_rng(10000 + 97 * sd + ci)
            W, M, P, thp, folds, sel, obs, null = run_cohort(rng, 100, 200, lam, 0.3, hyp, delta, fX, B)
            ns = int(sel.sum())
            if ns >= 5:
                zNs, zNr, zD, thpc = set_stats(obs, null, sel)
            else:
                zNs = zNr = zD = thpc = np.nan
            cl = classify(zD, zNs, thpc, ns, lam, ns / 200)
            key = (hyp if hyp != "mix" else f"mix{fX}")
            tally.setdefault(key, []).append(cl)
            print(f"{hyp:5s} {lam:.1f} {delta:.1f} {fX:.2f} | {sd} {ns:4d} {ns/200:.2f} | {zD:+5.1f} {zNs:+5.2f} {thpc:+.3f} {thp:+.2f} | {cl}", flush=True)
    print("\n-- mixtures: G=2000, 20% signal genes, n=100, Delta=0.2, R2=0.3, normal-tail BH (strength ~ 0.1) --")
    mixcells = (("post", 0.6, 0.5), ("post", 0.8, 0.5), ("noise", 0.5, 0.5), ("mix", 0.6, 0.1), ("mix", 0.6, 0.3), ("err", 0.5, 0.5)) if "b" in part else ()
    for (hyp, lam, fX) in mixcells:
        for sd in range(seeds):
            rng = np.random.default_rng(20000 + sd)
            W, M0, P0, _ = simulate(rng, 100, 20, 1600, lam, 0.3, "null", 0.0)
            rng2 = np.random.default_rng(30000 + sd + int(100 * fX))
            if hyp in ("noise", "post"):
                _, M1, P1, thp = simulate_onW(rng2, W, 400, lam, 0.3, hyp, 0.2)
            else:
                # mix / err on the shared W: regenerate via simulate with the same W by re-seeding simulate's draws
                W1, M1, P1, thp = simulate(rng2, 100, 20, 400, lam, 0.3, hyp, 0.2, fX=fX)
                # simulate() draws its own W; replace by projecting onto the shared W is not possible, so use its W for all
                W = W1; _, M0, P0, _ = simulate_onW(np.random.default_rng(40000 + sd), W, 1600, lam, 0.3, "null", 0.0)
            M = np.hstack([M0, M1]); P = np.hstack([P0, P1]); sig = np.r_[np.zeros(1600, bool), np.ones(400, bool)]
            folds = folds_of(rng, 100)
            sel, inc, _ = select(rng, W, M, P, folds, B, use_normal_tail=True)
            obs = discrim(W, M, P, folds, "insample"); null = joint_null(rng, W, M, P, folds, B, "insample")
            ns = int(sel.sum())
            if ns >= 5:
                zNs, zNr, zD, thpc = set_stats(obs, null, sel)
            else:
                zNs = zNr = zD = thpc = np.nan
            cl = classify(zD, zNs, thpc, ns, lam, ns / 2000)
            key = "MIXCOHORT-" + (hyp if hyp != "mix" else f"mix{fX}")
            tally.setdefault(key, []).append(cl)
            print(f"{hyp:5s} {lam:.1f} fX={fX:.2f} | {sd} n_sel={ns:4d} (false {int((sel & ~sig).sum())}) strength={ns/2000:.3f} | {zD:+5.1f} {zNs:+5.2f} {thpc:+.3f} {thp:+.2f} | {cl}", flush=True)
    print("\n-- tally --")
    for k, v in tally.items():
        from collections import Counter
        print(f"{k:18s} n={len(v):2d}  {dict(Counter(v))}")


# ---------------------------------------------------------------- F4
def F4(seeds=2):
    print("\n=== F4: kappa curve at a second point: G=2000 all-signal, ranked by published increment; theta_set/theta_pop of the top fraction ===")
    fr = (1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005)
    for (lam, delta) in ((0.6, 0.2), (0.6, 0.1), (0.7, 0.1)):
        for sd in range(seeds):
            rng = np.random.default_rng(5000 + sd)
            W, M, P, thp = simulate(rng, 100, 20, 2000, lam, 0.3, "noise", delta)
            folds = folds_of(rng, 100)
            inc = nested_incr(W, M, P, folds); obs = discrim(W, M, P, folds, "insample")
            order = np.argsort(-inc); out = []
            for f in fr:
                k = max(5, int(round(f * 2000))); idx = order[:k]
                out.append(obs["N"][idx].sum() / obs["D"][idx].sum())
            print(f"  lam={lam} delta={delta} theta_pop={thp:.2f} seed {sd}: " + " ".join(f"top{f*100:g}%={o/thp:.2f}" for f, o in zip(fr, out)), flush=True)


# ---------------------------------------------------------------- F5
def F5(seeds=2, G=150, B=100, qs=(0.0, 0.05, 0.10)):
    print("\n=== F5: H_post with per-gene MNAR protein truncation q (n chosen so ~100 fit patients remain), lam=0.7, R2=0.3, Delta=0.2 ===")
    print("q seed n_sel | z_D z_N(sign) theta_pc | class | diag: z of sum corr(W-pred M, M) on selected set / all genes")
    for q in qs:
        for sd in range(seeds):
            rng = np.random.default_rng(7000 + sd)
            n = int(round(100 / (1 - q))) if q > 0 else 100
            W, M, P, _ = simulate(rng, n, 20, G, 0.7, 0.3, "post", 0.2)
            keep = P >= np.quantile(P, q, axis=0) if q > 0 else np.ones_like(P, bool)
            perms = [rng.permutation(n) for _ in range(B)]
            fold_seed = int(rng.integers(1 << 30))

            def stats_for(Wx):
                N = np.empty(G); Nc = np.empty(G); D = np.empty(G); inc = np.empty(G); dg = np.empty(G); rho = np.empty(G)
                for g in range(G):
                    k = keep[:, g]; Wg = Wx[k]; Mg = M[k, g:g + 1]; Pg = P[k, g:g + 1]
                    fl = folds_of(np.random.default_rng(fold_seed + g), k.sum())
                    st = discrim(Wg, Mg, Pg, fl, "insample")
                    N[g], Nc[g], D[g], rho[g] = st["N"][0], st["Nc"][0], st["D"][0], st["rho"][0]
                    inc[g] = nested_incr(Wg, Mg, Pg, fl)[0]
                    Mz = (Mg - Mg.mean(0)) / Mg.std(0)
                    Mh = oof_wonly(wonly_hat(Wg, fl), Mz)
                    dg[g] = corr0(Mh, Mz)[0]
                return N, Nc, D, inc, dg, rho
            N, Nc, D, inc, dg, rho = stats_for(W)
            nN = np.empty((B, G)); nNc = np.empty((B, G)); nD = np.empty((B, G)); nI = np.empty((B, G)); ndg = np.empty((B, G))
            for b, pm in enumerate(perms):
                nN[b], nNc[b], nD[b], nI[b], ndg[b], _ = stats_for(W[pm])
            z = (inc - nI.mean(0)) / nI.std(0); sel = bh(norm.sf(z)) & (inc > 0)
            allg = np.ones(G, bool); a = np.abs(rho)
            ns = int(sel.sum())
            if ns >= 5:
                zN = setz(Nc / a, nNc / a, sel); zD = setz(D, nD, sel)
                thpc = (N[sel].sum() - nN[:, sel].sum(1).mean()) / (D[sel].sum() - nD[:, sel].sum(1).mean())
            else:
                zN = zD = thpc = np.nan
            cl = classify(zD, zN, thpc, ns, 0.7, ns / G)
            print(f"{q:.2f} {sd:4d} {ns:5d} | {zD:+5.1f} {zN:+5.2f} {thpc:+.3f} | {cl:4s} | {setz(dg, ndg, sel) if ns>=5 else np.nan:+.2f} / {setz(dg, ndg, allg):+.2f}", flush=True)


def F3a(): F3(part="a")
def F3b(): F3(part="b")
def F5b(): F5(seeds=2, G=150, B=100, qs=(0.02,))


if __name__ == "__main__":
    which = sys.argv[1:] or ["F2", "F1", "F4", "F3", "F5"]
    for w in which:
        t0 = time.time(); globals()[w](); print(f"[{w} done in {time.time()-t0:.0f}s]", flush=True)

