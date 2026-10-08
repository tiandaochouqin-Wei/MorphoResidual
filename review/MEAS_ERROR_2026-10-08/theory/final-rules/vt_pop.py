"""vt_pop.py -- independent population-limit checks of DERIVATION.md (synthetic only).
In-sample OLS everywhere (n large), own implementation. Checks:
  I1/I2, Delta_p closed forms (incl. corrected H_err denominator), Delta_m/Delta_p ratio,
  theta under H_noise/H_post/H_err/mix, sign-test threshold rho^2/(1-rho^2),
  nested-increment negativity, and the alternative scenarios (unshared X component, nonlinearity,
  purity mixing, MNAR truncation, heteroscedastic noise, replicate averaging, anti-aligned a,c).
"""
import numpy as np
rng = np.random.default_rng(12345)
n = 400_000

def c(a, b):
    return np.mean((a - a.mean()) * (b - b.mean()))

def stats(W, M, P):
    W = W - W.mean(0); M = M - M.mean(); P = P - P.mean()
    beta = c(M, P) / c(M, M); delta = c(M, P) / c(P, P); rho2 = beta * delta
    rp = P - beta * M; rm = M - delta * P
    # W-only OLS prediction of rp
    w = np.linalg.lstsq(W, rp, rcond=None)[0]; rp_hat = W @ w
    N = beta * c(rp_hat, M); D = c(rp_hat, rp); theta = N / D
    S = np.sign(beta) * c(rp_hat, rm) / np.sqrt(c(rp_hat, rp_hat) * c(rm, rm))
    # nested increments
    XM = np.column_stack([M, W]); fullP = XM @ np.linalg.lstsq(XM, P, rcond=None)[0]
    baseP = M * beta
    Dp = (c(fullP, fullP) - c(baseP, baseP)) / c(P, P)
    XP = np.column_stack([P, W]); fullM = XP @ np.linalg.lstsq(XP, M, rcond=None)[0]
    baseM = P * delta
    Dm = (c(fullM, fullM) - c(baseM, baseM)) / c(M, M)
    dp = fullP - baseP
    S_drm = np.sign(beta) * c(dp, rm) / np.sqrt(c(dp, dp) * c(rm, rm))
    I2 = c(rp, rm) / np.sqrt(c(rp, rp) * c(rm, rm))
    I1 = np.max(np.abs(rm - ((1 - rho2) * M - delta * rp)))
    return dict(beta=beta, rho2=rho2, theta=theta, N=N, D=D, S=S, Dp=Dp, Dm=Dm, S_drm=S_drm,
                I2=I2, I1err=I1, thr=rho2 / (1 - rho2))

def gen(lam, b, su2, se2, hyp, K=5, fX=0.5, a_norm=0.8, c_norm=0.8, align=None):
    """W = a X + c u + d e + eta (eta ~ N(0,I)). Returns W,M,P plus params."""
    X = rng.normal(size=n)
    u = np.sqrt(su2) * rng.normal(size=n)
    e = np.sqrt((1 - lam) / lam) * rng.normal(size=n)
    eps = np.sqrt(se2) * rng.normal(size=n)
    M = X + e; P = b * X + u + eps
    a = np.zeros(K); cc = np.zeros(K); d = np.zeros(K)
    if hyp == "noise":
        a[0] = a_norm
    elif hyp == "post":
        cc[1] = c_norm
    elif hyp == "err":
        d[2] = a_norm
    elif hyp == "mix":
        a[0] = a_norm * np.sqrt(fX); cc[1] = c_norm * np.sqrt(1 - fX)
    elif hyp == "mix_aligned":      # a and c on the same axis, same sign
        a[0] = a_norm * np.sqrt(fX); cc[0] = c_norm * np.sqrt(1 - fX)
    elif hyp == "mix_anti":         # a and c on the same axis, opposite sign
        a[0] = a_norm * np.sqrt(fX); cc[0] = -c_norm * np.sqrt(1 - fX)
    W = np.outer(X, a) + np.outer(u, cc) + np.outer(e, d) + rng.normal(size=(n, K))
    return W, M, P, dict(a=a, c=cc, d=d, X=X, u=u, e=e)

print("=== closed-form checks (R2_rna=0.3 at lam, var P = 1) ===")
for lam in (0.5, 0.7, 0.9):
    R2 = 0.3; b = np.sqrt(R2 / lam); su2 = 0.5 * (1 - b**2); se2 = 1 - b**2 - su2
    for hyp in ("noise", "post", "err"):
        W, M, P, pr = gen(lam, b, su2, se2, hyp)
        st = stats(W, M, P)
        # population quantities
        a, cc, d = pr["a"], pr["c"], pr["d"]
        SigW = np.eye(5) + np.outer(a, a) + np.outer(cc, cc) * su2 + np.outer(d, d) * (1 - lam) / lam
        rhoW = a @ np.linalg.solve(SigW, a)
        rhoU = su2 * cc @ np.linalg.solve(SigW, cc)
        rhoE = (1 - lam) / lam * d @ np.linalg.solve(SigW, d)
        sP2 = 1.0
        if hyp == "noise":
            Dp_pred = R2 * (1 - lam)**2 * rhoW / (lam * (1 - rhoW * lam))
            Dm_pred = lam * rhoW * (1 - R2 / lam)**2 / (1 - rhoW * R2 / lam)
            th_pred = lam / (1 - lam)
            ratio_text = (lam - R2)**2 / ((1 - lam)**2 * R2)
            ratio_mine = ratio_text * (1 - lam * rhoW) / (1 - rhoW * R2 / lam)
        elif hyp == "post":
            Dp_pred = rhoU * su2 / sP2
            Dm_pred = lam * rhoU * su2 * b**2 / (1 - rhoU * su2)
            th_pred = 0.0; ratio_text = R2; ratio_mine = R2 / (1 - Dp_pred)
        else:
            Dp_text = rhoE * R2 * (1 - lam)
            Dp_pred = rhoE * R2 * (1 - lam) / (1 - (1 - lam) * rhoE)
            Dm_pred = rhoE * (1 - lam)
            th_pred = -1.0; ratio_text = 1 / R2; ratio_mine = (1 - (1 - lam) * rhoE) / R2
        print(f"{hyp:5s} lam={lam:.1f} | I1err={st['I1err']:.1e} I2={st['I2']:+.3f}(-rho={-np.sqrt(R2):+.3f}) "
              f"| Dp sim={st['Dp']:.4f} pred={Dp_pred:.4f}" + (f" text={Dp_text:.4f}" if hyp == "err" else "") +
              f" | Dm sim={st['Dm']:.4f} pred={Dm_pred:.4f} | Dm/Dp sim={st['Dm']/st['Dp']:.3f} text={ratio_text:.3f} mine={ratio_mine:.3f}"
              f" | theta sim={st['theta']:+.3f} pred={th_pred:+.3f} | S={st['S']:+.3f} S_drm={st['S_drm']:+.3f}")

print("\n=== mixes: theta = f_X lam/(1-lam) where f_X = X-share of D (not of Delta_p); sign test threshold ===")
for lam, R2 in ((0.8, 0.2), (0.6, 0.3)):
    b = np.sqrt(R2 / lam); su2 = 0.5 * (1 - b**2); se2 = 1 - b**2 - su2
    thr_fX = (1 - lam) * R2 / (lam * (1 - R2))
    print(f"lam={lam} R2={R2}: S flips at f_X = {thr_fX:.3f}; theta threshold rho2/(1-rho2) = {R2/(1-R2):.3f}")
    for fX in (0.02, 0.05, 0.1, 0.3, 0.5, 0.8):
        W, M, P, pr = gen(lam, b, su2, se2, "mix", fX=fX)
        st = stats(W, M, P)
        # exact f_X of D: w = SigW^-1 C_p, C_p = a b(1-lam) + c su2
        a, cc = pr["a"], pr["c"]
        SigW = np.eye(5) + np.outer(a, a) + np.outer(cc, cc) * su2
        Cp = a * b * (1 - lam) + cc * su2
        w = np.linalg.solve(SigW, Cp)
        fX_D = (w @ a) * b * (1 - lam) / (w @ Cp)
        print(f"  fX_in={fX:.2f} fX_D={fX_D:.3f} theta sim={st['theta']:+.3f} pred={fX_D*lam/(1-lam):+.3f} "
              f"S={st['S']:+.4f} (pred sign {'+' if st['theta'] > st['thr'] else '-'})  Dp={st['Dp']:.4f} S_drm={st['S_drm']:+.3f}")

print("\n=== aligned / anti-aligned a and c (same morphological axis) ===")
lam, R2 = 0.7, 0.3; b = np.sqrt(R2 / lam); su2 = 0.5 * (1 - b**2); se2 = 1 - b**2 - su2
for hyp in ("mix_aligned", "mix_anti"):
    for fX in (0.1, 0.3, 0.5):
        W, M, P, pr = gen(lam, b, su2, se2, hyp, fX=fX)
        st = stats(W, M, P)
        print(f"{hyp:12s} fX={fX} theta={st['theta']:+.3f} D={st['D']:+.4f} Dp={st['Dp']:+.4f} S={st['S']:+.3f}")

print("\n=== alternatives ===")
# A. unshared transcript component: M = X1 + X2 + e, P = b X1 + u + eps, W reads X2 (or X1+X2)
lam = 0.8
for read in ("X2", "X1+X2", "X1"):
    X1 = rng.normal(size=n); X2 = 0.6 * rng.normal(size=n); e = np.sqrt((1 - lam) / lam) * rng.normal(size=n)
    u = 0.5 * rng.normal(size=n); eps = 0.5 * rng.normal(size=n)
    M = X1 + X2 + e; P = 0.6 * X1 + u + eps
    Z = {"X2": X2, "X1+X2": X1 + X2, "X1": X1}[read]
    W = np.column_stack([0.8 * Z + rng.normal(size=n), rng.normal(size=(n, 4))])
    st = stats(W, M, P)
    print(f"A unshared component, W reads {read:6s}: theta={st['theta']:+.3f} Dp={st['Dp']:.4f} S={st['S']:+.3f} rho2={st['rho2']:.3f}")

# B. nonlinearity with lam = 1: P = X + g (X^2-1), W reads X and X^2
for g in (0.2, 0.4):
    X = rng.normal(size=n); M = X.copy(); P = X + g * (X**2 - 1) + 0.6 * rng.normal(size=n)
    W = np.column_stack([X + rng.normal(size=n), (X**2 - 1) / np.sqrt(2) + rng.normal(size=n), rng.normal(size=(n, 3))])
    st = stats(W, M, P)
    print(f"B nonlinear (lam=1, g={g}): theta={st['theta']:+.3f} Dp={st['Dp']:.4f} S={st['S']:+.3f}  (reads 'H_post')")
    # same with lam = 0.8 and W reading X and X^2: mixed reading
    e = 0.5 * rng.normal(size=n); M2 = X + e
    st = stats(W, M2, P)
    print(f"B nonlinear (lam=0.8, g={g}): theta={st['theta']:+.3f} (pure H_noise would be 4.0) Dp={st['Dp']:.4f}")

# C. purity mixing: M = pi*Xt + (1-pi)*Xs + e ; P = pi*bt*Xt + (1-pi)*bs*Xs + eps ; W reads pi
for (mt, ms, bt, bs, lab) in ((2.0, 0.0, 1.0, 1.0, "tumour-enriched, equal translation"),
                              (0.0, 2.0, 1.0, 1.0, "stroma-enriched, equal translation"),
                              (1.0, 1.0, 1.0, 0.3, "equal expr, stroma translates less"),
                              (1.0, 1.0, 1.0, 1.0, "no composition effect (control)")):
    pi = np.clip(0.6 + 0.2 * rng.normal(size=n), 0.05, 0.98)
    Xt = mt + rng.normal(size=n); Xs = ms + rng.normal(size=n)
    M = pi * Xt + (1 - pi) * Xs + 0.3 * rng.normal(size=n)
    P = pi * bt * Xt + (1 - pi) * bs * Xs + 0.3 * rng.normal(size=n)
    W = np.column_stack([pi / pi.std() + rng.normal(size=n), rng.normal(size=(n, 4))])
    st = stats(W, M, P)
    print(f"C purity [{lab}]: theta={st['theta']:+.3f} Dp={st['Dp']:.4f} S={st['S']:+.3f} rho2={st['rho2']:.3f}")

# D. MNAR protein truncation under H_post (fit only patients with P above quantile q)
lam, R2 = 0.7, 0.3; b = np.sqrt(R2 / lam); su2 = 0.5 * (1 - b**2); se2 = 1 - b**2 - su2
for q in (0.0, 0.1, 0.3, 0.5):
    W, M, P, pr = gen(lam, b, su2, se2, "post")
    keep = P >= np.quantile(P, q)
    st = stats(W[keep], M[keep], P[keep])
    print(f"D MNAR H_post, drop bottom {int(q*100):2d}% of P: theta={st['theta']:+.3f} corr(W1,M)={np.corrcoef(W[keep,1], M[keep])[0,1]:+.3f} Dp={st['Dp']:.4f} S={st['S']:+.3f}")
for q in (0.3,):
    W, M, P, pr = gen(lam, b, su2, se2, "noise")
    keep = P >= np.quantile(P, q)
    st = stats(W[keep], M[keep], P[keep])
    print(f"D MNAR H_noise, drop bottom {int(q*100):2d}% of P: theta={st['theta']:+.3f} (pop {lam/(1-lam):.2f}) Dp={st['Dp']:.4f}")

# E. heteroscedastic noise (var e depends on W-readable quality variable), H_post and H_noise
for hyp in ("post", "noise"):
    X = rng.normal(size=n); u = np.sqrt(su2) * rng.normal(size=n); Q = rng.normal(size=n)
    sig_e = np.sqrt((1 - lam) / lam) * np.exp(0.5 * Q) / np.sqrt(np.mean(np.exp(Q)))
    e = sig_e * rng.normal(size=n); M = X + e; P = b * X + u + np.sqrt(se2) * rng.normal(size=n)
    Z = u if hyp == "post" else X
    W = np.column_stack([0.8 * Z / Z.std() + rng.normal(size=n), Q + rng.normal(size=n), rng.normal(size=(n, 3))])
    st = stats(W, M, P)
    lam_eff = 1 / (1 + np.mean(sig_e**2))
    print(f"E heteroscedastic e read by W, H_{hyp}: theta={st['theta']:+.3f} (pred {0 if hyp=='post' else lam_eff/(1-lam_eff):+.3f}, lam_eff={lam_eff:.3f}) Dp={st['Dp']:.4f}")

# F. replicate averaging: half the patients have k=3 libraries averaged, H_noise
X = rng.normal(size=n); k = np.where(rng.random(n) < 0.5, 3, 1)
e = np.sqrt((1 - lam) / lam) * rng.normal(size=n) / np.sqrt(k); M = X + e
P = b * X + np.sqrt(su2) * rng.normal(size=n) + np.sqrt(se2) * rng.normal(size=n)
W = np.column_stack([0.8 * X + rng.normal(size=n), rng.normal(size=(n, 4))])
st = stats(W, M, P); lam_eff = 1 / (1 + np.mean(e**2))
print(f"F replicate averaging H_noise: theta={st['theta']:+.3f} pred lam_eff/(1-lam_eff)={lam_eff/(1-lam_eff):+.3f} (single-library {lam/(1-lam):.2f})")

# G. ICC as lower bound: RNA/protein from same portion vs different portions
vc, vw, s2e = 1.0, 0.5, 0.4
Xc = rng.normal(size=n); Xw1 = np.sqrt(vw) * rng.normal(size=n); Xw2 = np.sqrt(vw) * rng.normal(size=n)
e = np.sqrt(s2e) * rng.normal(size=n); M = Xc + Xw1 + e
ICC = vc / (vc + vw + s2e); lam_port = (vc + vw) / (vc + vw + s2e)
W = np.column_stack([0.8 * Xc + rng.normal(size=n), rng.normal(size=(n, 4))])
for lab, P in (("same portion", 0.6 * (Xc + Xw1) + 0.6 * rng.normal(size=n)),
               ("different portion", 0.6 * (Xc + Xw2) + 0.6 * rng.normal(size=n))):
    st = stats(W, M, P)
    print(f"G ICC bound [{lab}]: theta={st['theta']:+.3f} -> lam_hat={st['theta']/(1+st['theta']):.3f}; ICC={ICC:.3f} lam_portion={lam_port:.3f}")
