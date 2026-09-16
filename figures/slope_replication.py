#!/usr/bin/env python3
"""External-replication slope chart (Fu Fig6c form): morphology incremental R^2 for the
SAME proteins in the discovery cohort (CPTAC-CCRCC, UNI) vs the independent cohort
(TCGA-KIRC, RPPA platform, Phikon). Each line = one protein; green = replicates in KIRC.
Needs (WinSCP from server): figdata/kirc_rppa_results.csv  (output of kirc_residual.py)
plus local figdata/merged_incr_ccrcc.csv."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import mrstyle as S

DD = "figdata"
kp = f"{DD}/kirc_rppa_results.csv"
if not os.path.exists(kp):
    raise SystemExit(f"missing {kp} -- run kirc_residual.py on the server and WinSCP it back")

kirc = pd.read_csv(kp)                                    # gene, incremental_r2, fdr, ...
ccrcc = pd.read_csv(f"{DD}/merged_incr_ccrcc.csv")        # gene, incremental_r2_uni, fdr_uni
m = ccrcc.merge(kirc, on="gene", suffixes=("_ccrcc", "_kirc"))
m = m.dropna(subset=["incremental_r2_uni", "incremental_r2"])
m = m[(m["fdr_uni"] < 0.05) & (m["incremental_r2_uni"] > 0)]      # CPTAC-CCRCC significant proteins only
rep = (m["fdr"] < 0.05) & (m["incremental_r2"] > 0)       # replicates in KIRC
rate = 100 * rep.mean()
# base rate: replication among ALL 360 tested RPPA proteins (the null a gene-specific claim must beat)
kirc_sig = (kirc["fdr"] < 0.05) & (kirc["incremental_r2"] > 0)
base_rate = 100 * kirc_sig.mean()
from scipy.stats import fisher_exact
non_m = kirc[~kirc["gene"].isin(m["gene"])]
a_, b_ = int(rep.sum()), len(m) - int(rep.sum())
c_, d_ = int(((non_m["fdr"] < 0.05) & (non_m["incremental_r2"] > 0)).sum()), len(non_m) - int(((non_m["fdr"] < 0.05) & (non_m["incremental_r2"] > 0)).sum())
_, fisher_p = fisher_exact([[a_, b_], [c_, d_]])

fig, ax = plt.subplots(figsize=(3.0, 3.0))
for _, row in m.iterrows():
    ok = (row["fdr"] < 0.05) and (row["incremental_r2"] > 0)
    ax.plot([0, 1], [row["incremental_r2_uni"], row["incremental_r2"]],
            color=(S.ORG["CCRCC"] if ok else S.LGREY), lw=0.7, alpha=0.75 if ok else 0.5,
            zorder=3 if ok else 1)
# medians
med = [m["incremental_r2_uni"].median(), m["incremental_r2"].median()]
ax.plot([0, 1], med, "o-", color=S.INK, lw=1.8, ms=5, zorder=5)
ax.set_xticks([0, 1]); ax.set_xticklabels(["CPTAC-CCRCC\n(UNI)", "TCGA-KIRC\n(RPPA · Phikon)"], fontsize=6.6)
ax.set_ylabel(r"morphology incremental $R^2$")
ax.set_xlim(-0.25, 1.25)
ax.set_title(f"{len(m)} CPTAC-significant proteins on RPPA: {rate:.0f}% replicate\n"
             f"vs.\\ {base_rate:.0f}% panel base rate (Fisher $P$={fisher_p:.2f}) — phenomenon-level",
             fontsize=6.8, fontweight="bold")
ax.axhline(np.median(kirc.loc[kirc_sig, "incremental_r2"]) if kirc_sig.any() else np.nan,
           color=S.GREY, lw=0.8, ls="--", zorder=0)
ax.text(0.5, 0.02, "blue = FDR<0.05 & incr.>0 in KIRC (replicates)", transform=ax.transAxes, fontsize=5.6,
        color=S.ORG["CCRCC"], ha="center")
S.save_pub(fig, "Fig_replication_slope")
print(f"wrote Fig_replication_slope  (shared={len(m)}, replicate={int(rep.sum())} = {rate:.0f}%, "
      f"base_rate={base_rate:.1f}%, fisher_p={fisher_p:.3f})")
