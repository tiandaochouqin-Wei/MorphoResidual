"""
sim_population_check.py -- large-n check (no permutations) that the closed forms in
DERIVATION.md are right: incremental R2 (own and mirror), theta = f_X lam/(1-lam) - f_E,
and the sign of the mirror sign test S_cross = sign(rho) corr(rp_hat, r_m).
Synthetic data only.   Usage: python sim_population_check.py
"""
import numpy as np
from sim_signtest import simulate_cohort, gene_stats, make_folds

rng = np.random.default_rng(7)
n, K, G = 20000, 20, 12
R2, su2 = 0.20, -0.5
print(f"n={n} K={K} genes/scenario={G}  R2_rna={R2} sigma_u^2={su2}")
print(f"{'hyp':5s} {'lam':>4s} {'fX':>4s} | {'incr_p':>7s} {'pop':>7s} | {'incr_m':>7s} {'pop':>7s} | "
      f"{'theta':>7s} {'pop':>7s} | {'S_cross':>8s} {'S_pp':>7s} {'S_dd':>7s} {'S_drm':>7s} | thr(rho^2)={R2}")
for hyp, fXs in (("noise", [None]), ("post", [None]), ("err", [None]), ("mix", [0.05, 0.1, 0.3, 0.5])):
    for lam in (0.6, 0.8, 0.9):
        for fX in fXs:
            d = 0.04
            W, M, P, pop = simulate_cohort(rng, n, K, G, lam, R2, su2, hyp, d, fX if fX else 0.5)
            if not pop["feasible"]:
                print(f"{hyp:5s} {lam:4.2f} {str(fX):>4s} | infeasible at delta={d} (needs rho>0.95)")
                continue
            Mz = (M - M.mean(0)) / M.std(0)
            Pz = (P - P.mean(0)) / P.std(0)
            st = gene_stats(W, Mz, Pz, make_folds(n, rng), 1.0)
            th = st["N"].sum() / st["D"].sum()
            print(f"{hyp:5s} {lam:4.2f} {str(fX):>4s} | {st['incr_p'].mean():7.4f} {pop['incr_p_pop']:7.4f} | "
                  f"{st['incr_m'].mean():7.4f} {pop['incr_m_pop']:7.4f} | {th:7.3f} {pop['theta_pop']:7.3f} | "
                  f"{st['S_cross'].mean():+8.4f} {st['S_pp'].mean():+7.4f} {st['S_dd'].mean():+7.4f} {st['S_drm'].mean():+7.4f}")
