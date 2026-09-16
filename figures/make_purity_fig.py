#!/usr/bin/env python3
"""Fig -- tumour-purity control (CCRCC). Per-gene morphology increment before vs after
additionally adjusting the mRNA baseline for tumour purity (server-verified
confound_check output). If the increment survives purity adjustment, the signal is not
a purity/cellularity artefact.

P0-3 FIX (adversarial review, 2026-09-06): this used to read a STALE pilot-stage file
(pilot_ccrcc/results_server/confound_check_results.csv; ~110-patient pilot cohort,
pre-deterministic-PCA-fix, giving the wrong "n=3,109" instead of the pinned 2,191) via
a hardcoded absolute path outside this project. It now reads figdata/confound_check_results.csv,
which must be produced by RE-RUNNING pilot_ccrcc/scripts_server/confound_check.py on the
server against the CURRENT pinned tumour-only residual_results.csv (n=103; the script
already has the svd_solver="full" determinism fix -- see review/REVIEW_REPORT.md P0-3
and review/recalc/ for the exact re-run instructions) and WinSCPing the output here."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as st
import mrstyle as S

_p = "figdata/confound_check_results.csv"
if not os.path.exists(_p):
    raise SystemExit(f"missing {_p} -- re-run confound_check.py on the server against the PINNED "
                      "tumour-only residual_results.csv (not the stale pilot file) and WinSCP it here; "
                      "see review/REVIEW_REPORT.md P0-3.")
d = pd.read_csv(_p)
d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=["incremental_r2", "incr_r2_purity_adjusted", "fdr"])
sig = d[(d.fdr < 0.05) & (d.incremental_r2 > 0) & (d.incremental_r2 < 1.5)
        & (d.incr_r2_purity_adjusted.abs() < 1.5)].copy()
x = sig.incremental_r2.values; y = sig.incr_r2_purity_adjusted.values
r = st.pearsonr(x, y)[0]
ratio = np.median(y / x); frac_pos = np.mean(y > 0); frac_half = np.mean(y > 0.5 * x)

fig, ax = plt.subplots(figsize=(3.4, 3.2))
ax.scatter(x, y, s=5, color=S.ORG["CCRCC"], alpha=0.45, edgecolors="none")
lim = [0, max(x.max(), y.max()) * 1.03]
ax.plot(lim, lim, ls="--", color=S.GREY, lw=0.9)
ax.text(0.04, 0.95, f"n = {len(sig):,} significant proteins\nr = {r:.2f}\n"
        f"median increment retained = {100*ratio:.0f}%\n{100*frac_pos:.0f}% remain > 0",
        transform=ax.transAxes, fontsize=6.5, va="top")
ax.set_xlabel(r"morphology incremental $R^2$ (mRNA baseline)", fontsize=7)
ax.set_ylabel(r"incremental $R^2$ after purity adjustment", fontsize=7)
ax.set_xlim(lim); ax.set_ylim(min(0, y.min()) - 0.02, lim[1])
ax.set_title("Not a tumour-purity artefact (CCRCC)", fontsize=8, fontweight="bold", loc="left")
S.save_pub(fig, "Fig_purity")
print(f"wrote Fig_purity  n={len(sig)} r={r:.3f} median_retained={100*ratio:.0f}% "
      f"frac>0={100*frac_pos:.0f}% frac>half={100*frac_half:.0f}%")
