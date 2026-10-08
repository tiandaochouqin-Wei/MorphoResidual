"""
pop_algebra.py -- independent, exact (no Monte Carlo) population-level re-derivation of the
quantities in ../DERIVATION.md.  Synthetic covariance algebra only; no project data.

Latent independent components  z = (X, e, u, eps, eta_1..eta_K)  with variances
(1, se2, su2, sg2, Seta).  Observables are linear in z:
    M = X + e,   P = b X + u + eps,   W = a X + c u + d e + eta        (a, c, d in R^K)
Every population quantity (OLS slopes, residuals, W-only prediction, nested projections) is a
linear functional of z, so covariances are exact quadratic forms.  Each check prints the value
obtained from the generic algebra next to the closed form claimed in DERIVATION.md.
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


class Model:
    def __init__(self, b, se2, su2, sg2, a, c, d, Seta):
        a, c, d = map(np.asarray, (a, c, d))
        K = len(a)
        self.K = K
        self.b, self.se2, self.su2, self.sg2 = b, se2, su2, sg2
        nz = 4 + K
        V = np.zeros(nz)
        V[0], V[1], V[2], V[3] = 1.0, se2, su2, sg2
        Sz = np.diag(V)
        Sz[4:, 4:] = np.asarray(Seta)
        self.Sz = Sz
        # loading vectors (rows) on z
        self.M = np.zeros(nz); self.M[0] = 1; self.M[1] = 1
        self.P = np.zeros(nz); self.P[0] = b; self.P[2] = 1; self.P[3] = 1
        Wl = np.zeros((K, nz))
        Wl[:, 0] = a; Wl[:, 2] = c; Wl[:, 1] = d; Wl[:, 4:] = np.eye(K)
        self.W = Wl
        self.lam = 1.0 / (1.0 + se2)
        self.sP2 = self.cov(self.P, self.P)
        self.rho2 = self.cov(self.M, self.P) ** 2 / (self.cov(self.M, self.M) * self.sP2)
        self.beta = self.cov(self.P, self.M) / self.cov(self.M, self.M)
        self.delta = self.cov(self.M, self.P) / self.sP2
        self.rp = self.P - self.beta * self.M
        self.rm = self.M - self.delta * self.P
        self.SW = Wl @ Sz @ Wl.T
        # W-only population regression of r_p (and r_m, M) on W
        self.Cp = Wl @ Sz @ self.rp
        self.Cm = Wl @ Sz @ self.rm
        self.CM = Wl @ Sz @ self.M
        self.w = np.linalg.solve(self.SW, self.Cp)
        self.rp_hat = self.w @ Wl
        self.rm_hat = np.linalg.solve(self.SW, self.Cm) @ Wl
        self.M_hat = np.linalg.solve(self.SW, self.CM) @ Wl
        # nested increments: projection of r_p on W residualised on M (and r_m on W residualised on P)
        WperpM = Wl - np.outer(self.CM / self.cov(self.M, self.M), self.M)
        SWM = WperpM @ Sz @ WperpM.T
        self.d_p = np.linalg.solve(SWM, WperpM @ Sz @ self.rp) @ WperpM
        CWP = Wl @ Sz @ self.P
        WperpP = Wl - np.outer(CWP / self.sP2, self.P)
        SWP = WperpP @ Sz @ WperpP.T
        self.d_m = np.linalg.solve(SWP, WperpP @ Sz @ self.rm) @ WperpP
        self.Delta_p = self.cov(self.d_p, self.d_p) / self.sP2
        self.Delta_m = self.cov(self.d_m, self.d_m) / self.cov(self.M, self.M)
        # discriminators
        self.N = self.beta * self.cov(self.rp_hat, self.M)
        self.D = self.cov(self.rp_hat, self.rp)
        self.theta = self.N / self.D if self.D > 0 else np.nan
        sg = np.sign(self.beta)
        self.S_cross = sg * self.corr(self.rp_hat, self.rm)
        self.S_pp = sg * self.corr(self.rp_hat, self.rm_hat)
        self.S_dd = sg * self.corr(self.d_p, self.d_m)
        self.S_drm = sg * self.corr(self.d_p, self.rm)
        # oracle R2 of X, u, e from W
        self.rhoW = (Wl @ Sz @ self.X()) @ np.linalg.solve(self.SW, Wl @ Sz @ self.X())
        uu = np.zeros(nz); uu[2] = 1
        ee = np.zeros(nz); ee[1] = 1
        self.rhoU = (Wl @ Sz @ uu) @ np.linalg.solve(self.SW, Wl @ Sz @ uu) / su2 if su2 > 0 else 0.0
        self.rhoE = (Wl @ Sz @ ee) @ np.linalg.solve(self.SW, Wl @ Sz @ ee) / se2 if se2 > 0 else 0.0

    def X(self):
        v = np.zeros(4 + self.K); v[0] = 1; return v

    def cov(self, u, v):
        return float(u @ self.Sz @ v)

    def corr(self, u, v):
        return self.cov(u, v) / np.sqrt(self.cov(u, u) * self.cov(v, v))


def line(tag, got, claim):
    flag = "ok" if abs(got - claim) < 1e-9 else ("~" if abs(got - claim) < 0.02 else "DIFF")
    print(f"  {tag:46s} algebra={got:+.5f}  claimed={claim:+.5f}  [{flag}]")


def run():
    K = 6
    rng = np.random.default_rng(1)
    # anisotropic eta so that Sigma_W is a general PD matrix (the derivation's 'general case')
    A = rng.normal(size=(K, K)); Seta = A @ A.T / K + 0.3 * np.eye(K)
    b, su2, sg2 = 0.8, 0.5, 0.4           # sigma_P^2 = 0.64+0.5+0.4 = 1.54
    for lam in (0.6, 0.8, 0.95):
        se2 = (1 - lam) / lam
        a0 = rng.normal(size=K); c0 = rng.normal(size=K); d0 = rng.normal(size=K)
        print(f"\n===== lam={lam}  (b={b}, su2={su2}, sg2={sg2}) =====")
        # --- identities I1, I2 (any model)
        m = Model(b, se2, su2, sg2, 0.7 * a0, 0.5 * c0, 0.3 * d0, Seta)
        print("I1/I2 (general model):")
        lhs = m.rm; rhs = (1 - m.rho2) * m.M - m.delta * m.rp
        print(f"  I1 max|r_m - [(1-rho2)M - delta r_p]| loadings = {np.abs(lhs - rhs).max():.2e}")
        line("I2 corr(r_p, r_m) = -rho", m.corr(m.rp, m.rm), -np.sqrt(m.rho2))
        line("rho^2 = b^2 lam / sP2", m.rho2, b * b * lam / m.sP2)
        line("beta = b lam", m.beta, b * lam)
        # --- H_noise
        m = Model(b, se2, su2, sg2, a0, 0 * c0, 0 * d0, Seta)
        rW, r2 = m.rhoW, m.rho2
        print("H_noise:")
        line("Delta_p = rho2 (1-lam)^2 rhoW/(lam(1-rhoW lam))", m.Delta_p, r2 * (1 - lam) ** 2 * rW / (lam * (1 - rW * lam)))
        line("  <= dilution bound rho2 (1-lam)/lam", m.Delta_p, r2 * (1 - lam) / lam)
        line("Delta_m = lam rhoW (1-rho2/lam)^2/(1-rhoW rho2/lam)", m.Delta_m, lam * rW * (1 - r2 / lam) ** 2 / (1 - rW * r2 / lam))
        line("theta = lam/(1-lam)", m.theta, lam / (1 - lam))
        line("Dm/Dp claimed (lam-rho2)^2/((1-lam)^2 rho2)", m.Delta_m / m.Delta_p, (lam - r2) ** 2 / ((1 - lam) ** 2 * r2))
        line("Dm/Dp corrected x(1-lam rhoW)/(1-rhoW rho2/lam)", m.Delta_m / m.Delta_p,
             (lam - r2) ** 2 / ((1 - lam) ** 2 * r2) * (1 - lam * rW) / (1 - rW * r2 / lam))
        line("S_cross = (lam-rho2) sqrt(rhoW)/sqrt(lam(1-rho2))", m.S_cross, (lam - r2) * np.sqrt(rW) / np.sqrt(lam * (1 - r2)))
        print(f"  S_pp={m.S_pp:+.3f}  S_dd={m.S_dd:+.3f}  S_drm={m.S_drm:+.3f}  rhoW={rW:.3f}")
        # --- H_post
        m = Model(b, se2, su2, sg2, 0 * a0, c0, 0 * d0, Seta)
        rU = m.rhoU
        print("H_post:")
        line("Delta_p = rhoU su2/sP2", m.Delta_p, rU * su2 / m.sP2)
        line("Delta_m = lam rhoU su2 b^2/sP2^2/(1-rhoU su2/sP2)", m.Delta_m, lam * rU * su2 * b * b / m.sP2 ** 2 / (1 - rU * su2 / m.sP2))
        line("theta = 0", m.theta, 0.0)
        print(f"  S_cross={m.S_cross:+.3f}  S_pp={m.S_pp:+.3f}  S_dd={m.S_dd:+.3f}  S_drm={m.S_drm:+.3f}  Dm/Dp={m.Delta_m/m.Delta_p:.3f} vs rho2={m.rho2:.3f}")
        # --- H_err
        m = Model(b, se2, su2, sg2, 0 * a0, 0 * c0, d0, Seta)
        rE = m.rhoE
        print("H_err:")
        line("Delta_p claimed rhoE rho2 (1-lam)", m.Delta_p, rE * m.rho2 * (1 - lam))
        line("Delta_p corrected rhoE rho2 (1-lam)/(1-(1-lam)rhoE)", m.Delta_p, rE * m.rho2 * (1 - lam) / (1 - (1 - lam) * rE))
        line("Delta_m = rhoE (1-lam)", m.Delta_m, rE * (1 - lam))
        line("theta = -1", m.theta, -1.0)
        line("Dm/Dp claimed 1/rho2", m.Delta_m / m.Delta_p, 1 / m.rho2)
        print(f"  S_cross={m.S_cross:+.3f}  S_pp={m.S_pp:+.3f}  S_dd={m.S_dd:+.3f}  S_drm={m.S_drm:+.3f}")
        # --- H_err with the X channel too (realistic: shared RNA artefact + morphology reads X)
        m = Model(b, se2, su2, sg2, a0, 0 * c0, d0, Seta)
        print(f"H_err+H_noise (a,d both): theta={m.theta:+.3f}  (between -1 and lam/(1-lam)={lam/(1-lam):.2f})")

    # --- sign-test threshold: theta < rho2  or  theta < rho2/(1-rho2) ?
    print("\n===== S_cross threshold (mix, a _|_ c, isotropic eta): scan f_X =====")
    K = 4; Seta = np.eye(K)
    a_dir = np.array([1, 0, 0, 0.]); c_dir = np.array([0, 1, 0, 0.])
    for lam, b in ((0.8, 0.8), (0.6, 1.0)):
        se2 = (1 - lam) / lam
        rows = []
        for amp in np.linspace(0.0, 3.0, 601):
            m = Model(b, se2, su2, sg2, amp * a_dir, 1.0 * c_dir, 0 * a_dir, Seta)
            rows.append((amp, m.theta, m.S_cross, m.S_pp, m.S_dd, m.rho2))
        rows = np.array(rows)
        r2 = rows[0, 5]
        for name, col in (("S_cross", 2), ("S_pp", 3), ("S_dd", 4)):
            i = np.where(np.diff(np.sign(rows[:, col])) != 0)[0]
            th_cross = rows[i[0] + 1, 1] if len(i) else np.nan
            print(f"  lam={lam} rho2={r2:.3f}: {name} changes sign at theta={th_cross:.4f}; "
                  f"rho2={r2:.4f}, rho2/(1-rho2)={r2/(1-r2):.4f}")

    # --- scale invariance of theta, ridge-direction (non-)invariance
    print("\n===== theta invariances =====")
    K = 6; rng = np.random.default_rng(3)
    A = rng.normal(size=(K, K)); Seta = A @ A.T / K + 0.3 * np.eye(K)
    lam = 0.7; se2 = (1 - lam) / lam
    a0 = rng.normal(size=K); c0 = rng.normal(size=K)
    for tag, (a, c) in {"noise": (a0, 0 * c0), "post": (0 * a0, c0), "mix": (a0, c0),
                        "mix a_|_c (Sigma_W metric)": (a0, None)}.items():
        if c is None:
            m0 = Model(1.0, se2, su2, sg2, a0, c0, 0 * a0, Seta)
            # make c Sigma_W-orthogonal to a: c <- c - a (a' SW^-1 c)/(a' SW^-1 a)  (approximately, SW depends on c)
            SWi = np.linalg.inv(m0.SW)
            c = c0 - a0 * (a0 @ SWi @ c0) / (a0 @ SWi @ a0)
        m = Model(1.0, se2, su2, sg2, a, c, 0 * a, Seta)
        # scaled M and P: multiply loadings (same as scaling the variables)
        vals = [m.theta]
        for alpha_ridge in (0.0, 0.5, 2.0, 10.0):
            w = np.linalg.solve(m.SW + alpha_ridge * np.eye(K), m.Cp)
            rp_hat = w @ m.W
            vals.append(m.beta * m.cov(rp_hat, m.M) / m.cov(rp_hat, m.rp))
        print(f"  {tag:28s} theta(OLS w)={vals[0]:+.4f}  ridge alpha=0/0.5/2/10 -> {vals[1]:+.4f} {vals[2]:+.4f} {vals[3]:+.4f} {vals[4]:+.4f}   lam/(1-lam)={lam/(1-lam):.3f}")


if __name__ == "__main__":
    run()
