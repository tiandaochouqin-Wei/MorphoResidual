"""vt_E7.py -- sign(rho_hat)-weighted vs rho_hat-weighted set statistic: H_post large selected set, and
H_noise / mix power check (synthetic only)."""
import numpy as np
from scipy.stats import norm
from vt_finite import simulate, folds_of, discrim, joint_null, select, setz

print("H_post, G=3000, lam=0.6, R2=0.3, delta=0.2, normal-tail BH selection; z on selected set")
for sd in range(3):
    rng = np.random.default_rng(8200 + sd)
    W, M, P, _ = simulate(rng, 100, 20, 3000, 0.6, 0.3, "post", 0.2)
    folds = folds_of(rng, 100)
    sel, inc, _ = select(rng, W, M, P, folds, 150, use_normal_tail=True)
    obs = discrim(W, M, P, folds, "insample"); null = joint_null(rng, W, M, P, folds, 150, "insample")
    rho = obs["rho"]; sgn = np.sign(rho)
    z1 = setz(obs["Nc"], null["Nc"], sel); z2 = setz(sgn * obs["Nc"] / rho, sgn * null["Nc"] / rho, sel)
    thpc = (obs["N"][sel].sum() - null["N"][:, sel].sum(1).mean()) / (obs["D"][sel].sum() - null["D"][:, sel].sum(1).mean())
    print(f"  seed {sd}: n_sel={sel.sum()}  z[rho*corr]={z1:+.2f}  z[sign*corr]={z2:+.2f}  theta_pc={thpc:+.4f}", flush=True)
print("power check, G=200, n=100: H_noise lam=0.5 delta=0.1 and mix fX=0.1 lam=0.6 delta=0.2 (R2=0.3)")
for hyp, lam, delta, fX in (("noise", 0.5, 0.1, 0.5), ("mix", 0.6, 0.2, 0.1), ("err", 0.5, 0.1, 0.5)):
    for sd in range(2):
        rng = np.random.default_rng(8300 + sd)
        W, M, P, thp = simulate(rng, 100, 20, 200, lam, 0.3, hyp, delta, fX=fX)
        folds = folds_of(rng, 100)
        sel, inc, _ = select(rng, W, M, P, folds, 150)
        obs = discrim(W, M, P, folds, "insample"); null = joint_null(rng, W, M, P, folds, 150, "insample")
        rho = obs["rho"]; sgn = np.sign(rho)
        if sel.sum() >= 5:
            z1 = setz(obs["Nc"], null["Nc"], sel); z2 = setz(sgn * obs["Nc"] / rho, sgn * null["Nc"] / rho, sel)
            th = obs["N"][sel].sum() / obs["D"][sel].sum()
            print(f"  {hyp:5s} seed {sd}: n_sel={sel.sum()} z[rho*corr]={z1:+.1f} z[sign*corr]={z2:+.1f} theta_set={th:+.2f} (pop {thp:+.2f})", flush=True)
        else:
            print(f"  {hyp:5s} seed {sd}: n_sel={sel.sum()} < 5")
