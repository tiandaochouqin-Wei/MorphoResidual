#!/usr/bin/env python3
"""Fig_replication -- all four RPPA replication attempts, not KIRC alone.

WHY THIS IS ITS OWN FIGURE
Fig_master2 panel d showed only TCGA-KIRC, because that panel's slope-plot design
was built before LUAD/GBM/PDAC were run. Adding three more panels to Fig_master2
would repeat the print-scale problem that made Fig_clinical and
Fig_wsi_spatial_supp their own figures in the first place (see make_fig_clinical.py's
docstring): Fig_master2 is already at its post-split 5-panel size. This figure
carries the full four-cohort replication result instead, at its own authored scale;
Fig_master2 panel d is kept as a single-cohort discrimination example and now
cross-references this figure and Table~tab:replication in its caption.

Each panel is the same slope-plot encoding as the original KIRC panel: every gene
significant in the matching CPTAC discovery cohort (FDR<0.05, positive incremental
R^2, UNI encoder) is drawn as a line from its CPTAC-side incremental R^2 to its
value in the RPPA cohort (Phikon encoder), coloured blue if it stays significant
there and grey if it does not, with the bold line marking the two medians. Titles
and the significant/base-rate/P numbers are read directly from
review/recalc/rppa_baserate_all.csv (the narrow, corrected-denominator frame used
throughout the manuscript -- see SuppTable_replication.tex) rather than recomputed
here, so the figure cannot silently drift from the reported statistics.

BRCA (the fifth attempt, CPTAC-2 mass spectrometry) is not a targeted-panel
platform and has no matching "CPTAC-significant subset on the panel" concept, so
it is reported in the main text only and is not a panel here.
"""
import pandas as pd
import matplotlib.pyplot as plt
import mrstyle as S

ORG, INK, LGREY = S.ORG, S.INK, S.LGREY
DD = "figdata"

COHORTS = [
    ("kirc",  "CCRCC", "TCGA-KIRC"),
    ("luad",  "LUAD",  "TCGA-LUAD"),
    ("gbm",   "GBM",   "TCGA-GBM"),
    ("pdac",  "PDAC",  "TCGA-PDAC"),
]
RPPA_FILE = {"kirc": "kirc_rppa_results.csv", "luad": "luad_rppa_results.csv",
             "gbm": "gbm_rppa_results.csv", "pdac": "paad_rppa_results.csv"}

stats = pd.read_csv("../review/recalc/rppa_baserate_all.csv").set_index("cohort")


def lab(ax, s, x=-0.24, y=1.10):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


fig = plt.figure(figsize=(7.2, 2.6))
gs = fig.add_gridspec(1, 4, wspace=0.55, left=0.075, right=0.99, top=0.72, bottom=0.30)

for i, (key, org, tcga) in enumerate(COHORTS):
    ccr = pd.read_csv(f"{DD}/merged_incr_{'ccrcc' if key == 'kirc' else key}.csv")
    ext = pd.read_csv(f"{DD}/{RPPA_FILE[key]}")
    mm = ccr.merge(ext, on="gene").dropna(subset=["incremental_r2_uni", "incremental_r2"])
    mm = mm[(mm.fdr_uni < 0.05) & (mm.incremental_r2_uni > 0)]
    okk = (mm.fdr < 0.05) & (mm.incremental_r2 > 0)

    ax = fig.add_subplot(gs[0, i])
    lab(ax, chr(ord("a") + i))
    for _, r in mm.iterrows():
        o = (r.fdr < 0.05) and (r.incremental_r2 > 0)
        ax.plot([0, 1], [r.incremental_r2_uni, r.incremental_r2],
                color=(ORG[org] if o else LGREY), lw=0.5, alpha=0.7 if o else 0.4,
                zorder=3 if o else 1)
    if len(mm):
        ax.plot([0, 1], [mm.incremental_r2_uni.median(), mm.incremental_r2.median()],
                "o-", color=INK, lw=1.6, ms=4, zorder=5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CPTAC\n(MS·UNI)", f"{tcga}\n(RPPA·Phikon)"], fontsize=6.3)
    ax.set_xlim(-0.35, 1.35)
    if i == 0:
        ax.set_ylabel("morphology incremental $R^2$", fontsize=7)
    s = stats.loc[key]
    ax.set_title(f"{org}: {s.rate_num:.0%} of {int(s.cptac_sig_on_rppa)}\n"
                 f"vs. {s.narrow_rate_base:.0%} base, $P$={s.narrow_p_two:.2f}",
                 fontsize=6.6, fontweight="bold", loc="left")

fig.text(0.5, 0.965,
         "Gene-specific replication does not exceed base rate in any of four independent RPPA cohorts",
         ha="center", fontsize=8, fontweight="bold")
fig.text(0.5, 0.015, "blue = FDR<0.05 & incr.>0 in the external cohort (colour = discovery organ)",
         ha="center", fontsize=6.3, color=INK, style="italic")

S.save_pub(fig, "Fig_replication")
print("wrote Fig_replication | per-panel N (CPTAC-sig ∩ RPPA-tested):",
      {k: int(stats.loc[k].cptac_sig_on_rppa) for k, _, _ in COHORTS})
