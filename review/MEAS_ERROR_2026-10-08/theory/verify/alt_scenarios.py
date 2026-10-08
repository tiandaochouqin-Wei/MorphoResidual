"""
alt_scenarios.py -- large-n Monte Carlo (n = 40,000, synthetic only) of alternatives that are
NOT in the H_noise / H_post / H_err family, to see what theta and the sign test would report.
In-sample OLS residuals and in-sample OLS W-regressions are used (population limit; no folds),
so every number is the population target of the proposed statistics.

Scenarios
  A  aligned channels: morphology reads X and u through the SAME morphological axis (a || c), lam = 1
  B  purity confounder: pi drives M and P with different loadings and W reads pi, lam = 1
  C  nonlinear transcript->protein map P = bX + g(X^2-1) + eps, W reads X and X^2, lam = 1
  D  MNAR protein: P missing when below a quantile, fit on observed, W reads X, lam = 1
  E  unshared transcript component: M = X1 + X2 + e, P = b X1 + eps, W reads X1  (biological 'dilution')
  F  heteroscedastic noise (var eps and var e depend on X), H_post truth -> theta?
  G  replicate averaging: M = mean of k libraries -> theta under H_noise
"""
import numpy as np

rng = np.random.default_rng(2026)
n = 40000


def ols_resid(y, x):
    x = x - x.mean(); y = y - y.mean()
    return y - x * (x @ y) / (x @ x)


def stats(W, M, P, label):
    """Population-limit versions of the proposed statistics (in-sample OLS everywhere)."""
    M = M - M.mean(); P = P - P.mean(); W = W - W.mean(0)
    beta = (M @ P) / (M @ M); delta = (M @ P) / (P @ P); rho2 = beta * delta
    rp = P - beta * M; rm = M - delta * P
    coef = np.linalg.lstsq(W, rp, rcond=None)[0]
    rp_hat = W @ coef
    N = beta * np.cov(rp_hat, M)[0, 1]; D = np.cov(rp_hat, rp)[0, 1]
    theta = N / D
    S = np.sign(beta) * np.corrcoef(rp_hat, rm)[0, 1]
    Dp = D / P.var()                               # W-only increment (= nested increment when W _|_ M)
    # nested increment proper
    X1 = np.column_stack([M, W]); full = X1 @ np.linalg.lstsq(X1, P, rcond=None)[0]
    Dp_nested = 1 - ((P - full) ** 2).mean() / P.var() - rho2
    lamhat = theta / (1 + theta) if theta > -1 else np.nan
    print(f"  {label:58s} R2rna={rho2:.3f} Dp_nested={Dp_nested:+.4f} theta={theta:+.3f} lamhat={lamhat:+.3f} S_cross={S:+.3f}")
    return theta


print("A. aligned channels, NO measurement error (lam = 1): W = a X + c u + eta with a || c")
X = rng.normal(size=n); u = rng.normal(size=n); eps = rng.normal(size=n)
b = 0.7
for cos_ac in (-1.0, -0.5, 0.0, 0.5, 1.0):
    K = 4
    a = np.array([1.0, 0, 0, 0]); c_perp = np.array([0, 1.0, 0, 0])
    c = cos_ac * a + np.sqrt(1 - cos_ac ** 2) * c_perp
    W = np.outer(X, 1.0 * a) + np.outer(u, 0.8 * c) + 0.5 * rng.normal(size=(n, K))
    M = X
    P = b * X + 0.8 * u + 0.6 * eps
    stats(W, M, P, f"cos(a,c)={cos_ac:+.1f}, lam=1 (true theta should read 0 = H_post)")

print("\nB. purity confounder pi -> M and P (different loadings) and morphology reads pi; lam = 1")
pi = rng.normal(size=n)
for (gM, gP) in ((0.6, 0.3), (0.6, 0.6 * b), (0.6, 0.8), (0.6, -0.3)):
    X0 = rng.normal(size=n); u0 = rng.normal(size=n); eps = rng.normal(size=n)
    Xtot = np.sqrt(1 - gM ** 2) * X0 + gM * pi            # true transcript includes purity effect
    M = Xtot                                              # no measurement error
    P = b * np.sqrt(1 - gM ** 2) * X0 + gP * pi + 0.6 * u0 + 0.6 * eps   # purity acts on P with loading gP (not b*gM)
    W = np.column_stack([pi + 0.5 * rng.normal(size=n), 0.5 * rng.normal(size=(n, 3))])
    tag = "purity on P exactly b*gM (pure X-channel)" if abs(gP - b * gM) < 1e-9 else f"purity loading on P = {gP:+.2f} (b*gM={b*gM:.2f})"
    stats(W, M, P, "W reads purity only; " + tag)

print("\nC. nonlinear transcript->protein map, lam = 1, morphology reads X and X^2")
for skew in (0.0, 1.0):
    for g in (0.3, 0.6):
        Z = rng.normal(size=n)
        X = Z if skew == 0 else (np.exp(0.6 * Z) - np.exp(0.18)) / np.sqrt((np.exp(0.36) - 1) * np.exp(0.36))
        X = (X - X.mean()) / X.std()
        Q = X ** 2 - 1; Q = (Q - Q.mean()) / Q.std()
        eps = rng.normal(size=n)
        P = b * X + g * Q + 0.6 * eps
        W = np.column_stack([X + 0.5 * rng.normal(size=n), Q + 0.5 * rng.normal(size=n), 0.5 * rng.normal(size=(n, 2))])
        stats(W, X, P, f"P = bX + {g}(X^2-1) + eps, skew={skew}, lam=1 (no u at all)")

print("\nD. MNAR protein (observed only if P above its q-quantile), lam = 1, W reads X only")
for q in (0.0, 0.2, 0.4):
    X = rng.normal(size=n); eps = rng.normal(size=n)
    P = b * X + 0.7 * eps
    W = np.column_stack([X + 0.4 * rng.normal(size=n), 0.5 * rng.normal(size=(n, 3))])
    keep = P > np.quantile(P, q)
    stats(W[keep], X[keep], P[keep], f"missing fraction={q}, W reads X (no u, no e)")

print("\nE. unshared transcript component: M = X1 + X2 + e (lam_tech), P = b X1 + u + eps, W reads X1")
for (v2, lam_t) in ((0.0, 0.9), (0.5, 0.9), (0.5, 1.0), (1.0, 0.9)):
    X1 = rng.normal(size=n); X2 = np.sqrt(v2) * rng.normal(size=n)
    e = np.sqrt((1 - lam_t) / lam_t * (1 + v2)) * rng.normal(size=n)
    M = X1 + X2 + e
    P = b * X1 + 0.6 * rng.normal(size=n)
    W = np.column_stack([X1 + 0.4 * rng.normal(size=n), 0.5 * rng.normal(size=(n, 3))])
    lam_eff = 1.0 / (1 + v2 + e.var())
    th = stats(W, M, P, f"var X2={v2}, lam_tech={lam_t}: lam_eff={lam_eff:.3f}, lam_eff/(1-lam_eff)={lam_eff/(1-lam_eff):.2f}")

print("\nF. heteroscedastic noise (sd of e and eps grow with X), H_post truth")
X = rng.normal(size=n); u = rng.normal(size=n)
for het in (0.0, 0.8):
    e = np.sqrt(0.4) * rng.normal(size=n) * np.exp(het * X / 2) / np.sqrt(np.exp(het ** 2 / 2))
    eps = 0.6 * rng.normal(size=n) * np.exp(het * X / 2) / np.sqrt(np.exp(het ** 2 / 2))
    M = X + e
    P = b * X + 0.8 * u + eps
    W = np.column_stack([u + 0.5 * rng.normal(size=n), 0.5 * rng.normal(size=(n, 3))])
    stats(W, M, P, f"het={het}: H_post truth, lam_marginal={1/(1+e.var()):.2f}")

print("\nG. replicate averaging under H_noise: single-library lam=0.6, M = mean of k libraries")
for k in (1, 2, 3):
    X = rng.normal(size=n)
    libs = [X + np.sqrt(0.4 / 0.6) * rng.normal(size=n) for _ in range(k)]
    M = np.mean(libs, axis=0)
    P = b * X + 0.8 * rng.normal(size=n)
    W = np.column_stack([X + 0.3 * rng.normal(size=n), 0.5 * rng.normal(size=(n, 3))])
    lam_avg = k * 0.6 / (1 + (k - 1) * 0.6)
    stats(W, M, P, f"k={k}: Spearman-Brown lam_avg={lam_avg:.3f} -> expected theta={lam_avg/(1-lam_avg):.2f}")
