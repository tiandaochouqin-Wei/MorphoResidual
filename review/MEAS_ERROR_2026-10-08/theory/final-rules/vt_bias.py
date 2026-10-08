"""vt_bias.py -- where does the small positive excess of N_c under H_post come from? (synthetic)
Compare per-gene mean of corr(rp_hat, M) and cov(rp_hat, M) under: no-signal cohort, H_post cohort,
each against its own W-permutation null; in-sample vs in-fold residuals; and a 'signal-only' decomposition.
"""
import numpy as np
from vt_finite import simulate, folds_of, discrim, joint_null, wonly_hat, oof_wonly, cov0, corr0

for hyp in ("null", "post"):
    for sd in range(2):
        rng = np.random.default_rng(8000 + sd)
        W, M, P, _ = simulate(rng, 100, 20, 3000, 0.6, 0.3, hyp, 0.2)
        folds = folds_of(rng, 100)
        for mode in ("insample", "infold"):
            obs = discrim(W, M, P, folds, mode); null = joint_null(rng, W, M, P, folds, 100, mode)
            c_obs = (obs["Nc"] / obs["rho"]).mean(); c_null = (null["Nc"] / obs["rho"]).mean()
            v_obs = (obs["N"] / obs["rho"]).mean(); v_null = (null["N"] / obs["rho"]).mean()
            z = (obs["Nc"].sum() - null["Nc"].sum(1).mean()) / null["Nc"].sum(1).std()
            print(f"{hyp:4s} seed{sd} {mode:8s}: mean corr(rp_hat,M) obs={c_obs:+.5f} null={c_null:+.5f} diff={c_obs-c_null:+.5f} | "
                  f"mean cov obs={v_obs:+.5f} null={v_null:+.5f} | sd(rp_hat) obs={np.sqrt(cov0(obs['D'],obs['D'])).mean() if False else 0:.0f} | z_all={z:+.2f}", flush=True)
# decomposition under H_post, in-sample: rp_hat = H(P) - rho H(M)  (both OOF); correlate pieces with M
rng = np.random.default_rng(9000)
W, M, P, _ = simulate(rng, 100, 20, 3000, 0.6, 0.3, "post", 0.2)
folds = folds_of(rng, 100); hats = wonly_hat(W, folds)
Mz = (M - M.mean(0)) / M.std(0); Pz = (P - P.mean(0)) / P.std(0); rho = cov0(Mz, Pz)
HP = oof_wonly(hats, Pz); HM = oof_wonly(hats, Mz); rp_hat = HP - rho * HM
print("H_post obs: mean corr(HP,M)=%+.5f  mean corr(HM,M)=%+.5f  mean corr(rp_hat,M)=%+.5f  mean sd(HP)=%.3f sd(HM)=%.3f sd(rp_hat)=%.3f" % (
    corr0(HP, Mz).mean(), corr0(HM, Mz).mean(), corr0(rp_hat, Mz).mean(), HP.std(0).mean(), HM.std(0).mean(), rp_hat.std(0).mean()))
pm = rng.permutation(100); hats2 = wonly_hat(W[pm], folds)
HP2 = oof_wonly(hats2, Pz); HM2 = oof_wonly(hats2, Mz); rp2 = HP2 - rho * HM2
print("perm null : mean corr(HP,M)=%+.5f  mean corr(HM,M)=%+.5f  mean corr(rp_hat,M)=%+.5f  mean sd(HP)=%.3f sd(HM)=%.3f sd(rp_hat)=%.3f" % (
    corr0(HP2, Mz).mean(), corr0(HM2, Mz).mean(), corr0(rp2, Mz).mean(), HP2.std(0).mean(), HM2.std(0).mean(), rp2.std(0).mean()))
# cov pieces
print("H_post obs: mean cov(HP,M)=%+.5f cov(HM,M)=%+.5f ; perm: cov(HP,M)=%+.5f cov(HM,M)=%+.5f" % (
    cov0(HP, Mz).mean(), cov0(HM, Mz).mean(), cov0(HP2, Mz).mean(), cov0(HM2, Mz).mean()))
