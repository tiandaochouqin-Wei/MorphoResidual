"""vt_E6.py -- large H_post sets: is the positive z of sum(rho_hat * corr) removed by sign(rho_hat) weighting
or by rho_pop-free weighting? Also check how the excess depends on rho_hat (synthetic only)."""
import numpy as np
from vt_finite import simulate, folds_of, discrim, joint_null

for sd in range(4):
    rng = np.random.default_rng(8100 + sd)
    W, M, P, _ = simulate(rng, 100, 20, 3000, 0.6, 0.3, "post", 0.2)
    folds = folds_of(rng, 100)
    obs = discrim(W, M, P, folds, "insample"); null = joint_null(rng, W, M, P, folds, 120, "insample")
    rho = obs["rho"]; corr_o = obs["Nc"] / rho; corr_n = null["Nc"] / rho
    def z(wt):
        o = (wt * corr_o).sum(); nn = (wt * corr_n).sum(1); return (o - nn.mean()) / nn.std()
    # excess by rho_hat tertile
    q = np.quantile(rho, [1 / 3, 2 / 3]); ex = corr_o - corr_n.mean(0)
    lo, mid, hi = ex[rho <= q[0]].mean(), ex[(rho > q[0]) & (rho <= q[1])].mean(), ex[rho > q[1]].mean()
    print(f"seed {sd}: z[rho*corr]={z(rho):+.2f}  z[sign*corr]={z(np.sign(rho)):+.2f}  z[corr]={z(np.ones_like(rho)):+.2f}  "
          f"z[(rho-mean rho)*corr]={z(rho - rho.mean()):+.2f} | mean excess corr by rho tertile: lo {lo:+.5f} mid {mid:+.5f} hi {hi:+.5f} | rho mean {rho.mean():.3f} sd {rho.std():.3f}", flush=True)
