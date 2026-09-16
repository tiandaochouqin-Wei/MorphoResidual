#!/usr/bin/env python3
"""Composite Figure 2 — robustness & interpretation.  All native, server-verified numbers.
a residual programme (enrichment) | b composition controls | c interpretable nuclear
correlates | d residual score tracks grade | e survival (KM, descriptive) | f survival
adjusted for grade/stage (honest forest).  [slots for slope-replication + WSI heatmap]."""
import os          # was imported inside the removed g/h block
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats as st
from scipy.stats import chi2
import mrstyle as S

ORG, INK, GREY, LGREY, FAM = S.ORG, S.INK, S.GREY, S.LGREY, S.FAM
COH, ORGAN = S.COH, S.ORGAN
DD = "figdata"


def lab(ax, s, x=-0.18, y=1.05):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


# Split 2026-09-09: the clinical block (old d,e,f) moved to make_fig_clinical.py and the
# WSI spatial maps (old g,h) to make_fig_wsi_spatial_supp.py. Ten panels on one canvas
# forced a 60.7% print scale and 2.2 pt labels; five panels print at ~90%.
fig = plt.figure(figsize=(7.2, 7.0))
gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 1, 1], hspace=0.68, wspace=0.34,
                      left=0.10, right=0.955, top=0.955, bottom=0.075)

# ===== a : residual programme (leading enrichment per organ) =====
axa = fig.add_subplot(gs[0, :]); lab(axa, "a", x=-0.075)
# tested-proteome-background ORA (P0-5, adversarial review 2026-09-06);
# see enrichment_tested_background.py / figdata/enrichment_tested_background_full.csv.
# The representative term per organ is a manual pick of the clearest name from among that
# organ's top handful of terms in the same family (e.g. CCRCC's smallest-q term is a
# TRiC/CCT-cooperation term at q=0.0083 vs this term's q=0.0113 -- both are the folding
# chaperonin complex, this name is just more legible); the q VALUE plotted is read fresh
# from the file below, not typed in, so it cannot silently drift from a re-run.
_lead_term = {"CCRCC": "Chaperonin-mediated Protein Folding",
              "LUAD": "Cap-dependent Translation Initiation",
              "UCEC": "Asparagine N-linked Glycosylation",
              "GBM": "mRNA Splicing - Major Pathway",
              "PDAC": "Extracellular Matrix Organization"}
_lead_label = {"CCRCC": "Chaperonin protein folding", "LUAD": "Cap-dependent translation",
               "UCEC": "N-linked glycosylation", "GBM": "mRNA splicing (major)",
               "PDAC": "ECM organisation"}
_lead_fam = {"CCRCC": "folding", "LUAD": "translation", "UCEC": "secretion",
             "GBM": "splicing", "PDAC": "ecm"}
_efull = pd.read_csv(f"{DD}/enrichment_tested_background_full.csv")
lead = {}
for _c, _term in _lead_term.items():
    _row = _efull[(_efull.cohort == _c) & _efull.Term.str.startswith(_term)]
    assert len(_row) == 1, f"expected exactly one match for {_c}/{_term}, got {len(_row)}"
    _q = float(_row["Adjusted P-value"].iloc[0])
    lead[_c] = (_lead_label[_c], -np.log10(_q), _lead_fam[_c])
# single-row dot matrix (same encoding as the Figure-10 term x cohort matrix, make_fig3_dotmatrix.py):
# one dot per organ, AREA = -log10 q (capped at 30), colour = programme family.
import textwrap
CAP_A = 30.0


def dot_area(v):                     # v = -log10 q ; marker area in pt^2
    return 14.0 + 200.0 * min(max(v, 0.0), CAP_A) / CAP_A


Y_DOT, Y_ORG, Y_TERM = 0.60, 0.84, 0.44
for j, c in enumerate(COH):
    term, v, famk = lead[c]
    axa.scatter(j, Y_DOT, s=dot_area(v), facecolor=FAM[famk], edgecolor="white", linewidth=0.4, zorder=3)
    axa.text(j, Y_ORG, ORGAN[c], ha="center", va="bottom", fontsize=7, fontweight="bold", color=INK)
    axa.text(j, Y_TERM, textwrap.fill(term, 15), ha="center", va="top", fontsize=7, color=INK,
             linespacing=1.05)
    axa.text(j, Y_TERM - 0.27, f"{v:.1f}", ha="center", va="top", fontsize=7, color=GREY)
axa.set_xlim(-0.55, 4.55); axa.set_ylim(-0.32, 1.0)
for s_ in axa.spines.values():
    s_.set_visible(False)
axa.set_xticks([]); axa.set_yticks([])
axa.set_title("Organ-specific post-transcriptional programme", fontsize=7.6, fontweight="bold", loc="left")
# Unicode rather than mathtext: matplotlib sets sub/superscripts at 0.7x the base
# size, so a 7 pt "$\log_{10}$" ships 4.9 pt glyphs and fails the 5 pt print floor.
axa.text(-0.5, Y_TERM - 0.27, "−log10 q:", ha="right", va="top", fontsize=7, color=GREY, clip_on=False)
from matplotlib.lines import Line2D
leg_fam = [Line2D([], [], ls="none", marker="o", ms=4.2, mfc=FAM[k], mec="white", mew=0.3, label=k)
           for k in ["translation", "secretion", "splicing", "folding", "ecm"]]
lg1 = axa.legend(handles=leg_fam, fontsize=7, frameon=False, ncol=3, loc="lower left",
                 handlelength=0.8, handletextpad=0.3, columnspacing=0.8, labelspacing=0.2,
                 bbox_to_anchor=(-0.02, -0.20), title="family", title_fontsize=7, alignment="left")
axa.add_artist(lg1)
leg_sz = [Line2D([], [], ls="none", marker="o", ms=np.sqrt(dot_area(v)), mfc=GREY, mec="white", mew=0.3, label=l)
          for v, l in [(2.0, "q = 0.01"), (20.0, "q = 1e−20")]]
axa.legend(handles=leg_sz, fontsize=7, frameon=False, ncol=1, loc="lower right",
           handlelength=2.4, handletextpad=1.0, labelspacing=0.7,
           bbox_to_anchor=(1.03, -0.20), title="dot area = −log10 BH q\n(tested-background ORA, cap 30)",
           title_fontsize=7, alignment="left")

# ===== b : composition controls =====
axb = fig.add_subplot(gs[1, 0]); lab(axb, "b")
# Computed from the per-gene composition_controls tables (Source Data) rather than
# hardcoded: retained % = median adjusted increment / median original increment over
# the significant set. The previously hardcoded values are kept as an assertion so a
# silent drift in the data or the definition fails loudly instead of re-plotting.
_HARD = {"stroma": [87, 93, 70, 86, 85], "immune": [99, 78, 92, 77, 90],
         "proliferation": [86, 32, 83, 44, 48], "all three": [68, 6, 36, 3, 34]}
_KEY = {"stroma": "incr_W_over_stroma", "immune": "incr_W_over_immune",
        "proliferation": "incr_W_over_proliferation", "all three": "incr_W_over_ALL"}
retain = {}
for _ax, _col in _KEY.items():
    vals = []
    for _c in COH:
        _d = pd.read_csv(f"{DD}/composition_controls_{_c.lower()}.csv").dropna(subset=["incr_orig", _col])
        vals.append(100 * _d[_col].median() / _d["incr_orig"].median())
    retain[_ax] = vals
    for _v, _h in zip(vals, _HARD[_ax]):
        assert abs(_v - _h) <= 1.0, f"composition {_ax}: data gives {_v:.1f}, figure text says {_h}"
retain = {k: [int(round(v)) for v in vs] for k, vs in retain.items()}
# annotated heatmap (rows = adjustment, cols = cohort) instead of 20 grouped bars: the quantity
# is a 4x5 matrix of "% of the morphology increment retained", and the cohort x adjustment
# structure (LUAD/GBM collapse only under proliferation) reads directly off the grid.
M = np.array([retain[k] for k in retain])
im = axb.imshow(M, cmap="Blues", vmin=0, vmax=100, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        axb.text(j, i, f"{v:.0f}%", ha="center", va="center", fontsize=7,
                 color="white" if v > 60 else INK, fontweight="bold" if v < 10 else "normal")
axb.set_xticks(range(5)); axb.set_xticklabels(COH, fontsize=7)
axb.set_yticks(range(4)); axb.set_yticklabels(["+ stroma", "+ immune", "+ proliferation", "+ all three"], fontsize=7)
for s_ in axb.spines.values():
    s_.set_visible(False)
axb.tick_params(length=0)
axb.set_title("Not explained by\nmicroenvironment composition", fontsize=7.6, fontweight="bold", loc="left")
axb.text(0.5, -0.20, "% of morphology increment retained after adding the score(s) to the mRNA baseline",
         transform=axb.transAxes, ha="center", fontsize=7, color=GREY)

# ===== c : interpretable nuclear correlates (pathomics, translation) =====
axc = fig.add_subplot(gs[1, 1]); lab(axc, "c")
FEAT = ["nuclear density", "nucleus count", "mean nucleus area", "nucleus area SD",
        "chromatin OD", "eosin fraction"]
_FKEY = ["nuclear_density", "nucleus_count", "mean_nuc_area", "nuc_area_std",
         "chromatin_od", "eosin_fraction"]
# server_export/results/results__clinical_link_<cohort>_pathomics_corr.csv (long format:
# feature, score, rho, p); score=="translation_meas", pivoted feature x cohort. Verified
# 2026-09-17 to reproduce the previously hand-typed matrix to 2dp for every value checked.
_transl_cols = []
for _c in COH:
    _pc = pd.read_csv(f"../server_export/results/results__clinical_link_{_c.lower()}_pathomics_corr.csv")
    _pc = _pc[_pc.score == "translation_meas"].set_index("feature")["rho"]
    _transl_cols.append([float(_pc[k]) for k in _FKEY])
transl = np.array(_transl_cols).T
im = axc.imshow(transl, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
axc.set_xticks(range(5)); axc.set_xticklabels(COH, rotation=45, ha="right", fontsize=7)
axc.set_yticks(range(len(FEAT))); axc.set_yticklabels(FEAT, fontsize=7)
for i in range(transl.shape[0]):
    for j in range(transl.shape[1]):
        axc.text(j, i, f"{transl[i,j]:+.2f}", ha="center", va="center", fontsize=7,
                 color="white" if abs(transl[i, j]) > 0.28 else "#333")
for s in axc.spines.values():
    s.set_visible(False)
axc.set_title("Nuclear correlates (translation residual)", fontsize=7.6, fontweight="bold", loc="left")
cb = fig.colorbar(im, ax=axc, fraction=0.045, pad=0.03); cb.set_label(r"Spearman $\rho$", fontsize=7)
cb.ax.tick_params(labelsize=7)

# ===== d : external replication slope (CPTAC-CCRCC -> TCGA-KIRC; RPPA platform, Phikon) =====
kp = f"{DD}/kirc_rppa_results.csv"
if os.path.exists(kp):
    kirc = pd.read_csv(kp); ccr = pd.read_csv(f"{DD}/merged_incr_ccrcc.csv")
    mm = ccr.merge(kirc, on="gene").dropna(subset=["incremental_r2_uni", "incremental_r2"])
    mm = mm[(mm.fdr_uni < 0.05) & (mm.incremental_r2_uni > 0)]
    okk = (mm.fdr < 0.05) & (mm.incremental_r2 > 0)
    kirc_sig_all = (kirc.fdr < 0.05) & (kirc.incremental_r2 > 0)
    base_rate = 100 * kirc_sig_all.mean()
    non_m = kirc[~kirc["gene"].isin(mm["gene"])]
    non_sig = int(((non_m.fdr < 0.05) & (non_m.incremental_r2 > 0)).sum())
    from scipy.stats import fisher_exact
    _, fisher_p = fisher_exact([[int(okk.sum()), len(mm) - int(okk.sum())],
                                 [non_sig, len(non_m) - non_sig]])
    axi = fig.add_subplot(gs[2, 0]); lab(axi, "d")
    for _, r in mm.iterrows():
        o = (r.fdr < 0.05) and (r.incremental_r2 > 0)
        axi.plot([0, 1], [r.incremental_r2_uni, r.incremental_r2], color=(ORG["CCRCC"] if o else LGREY),
                 lw=0.6, alpha=0.75 if o else 0.5, zorder=3 if o else 1)
    axi.plot([0, 1], [mm.incremental_r2_uni.median(), mm.incremental_r2.median()], "o-", color=INK,
             lw=1.6, ms=4.5, zorder=5)
    axi.set_xticks([0, 1]); axi.set_xticklabels(["CPTAC-CCRCC\n(MS · UNI)", "TCGA-KIRC\n(RPPA · Phikon)"], fontsize=7)
    axi.set_xlim(-0.3, 1.3); axi.set_ylabel("morphology incremental R²", fontsize=7)
    axi.set_title(f"Phenomenon-level replication: {100*okk.mean():.0f}% of {len(mm)} vs.\n"
                  f"{base_rate:.0f}% panel base rate ($P$={fisher_p:.2f})",
                  fontsize=7, fontweight="bold", loc="left")
    axi.text(0.5, 0.03, "blue = FDR<0.05 & incr.>0 in KIRC", transform=axi.transAxes, fontsize=7,
             color=ORG["CCRCC"], ha="center")
else:
    fig.text(0.5, 0.01, "external-replication slope pending (kirc_rppa_results.csv)", ha="center",
             fontsize=7, color=LGREY, style="italic")

# ===== e : aggregation ablation (ABMIL vs mean-pool) =====
ap_ = f"{DD}/abmil_ccrcc.csv"
if os.path.exists(ap_):
    ab = pd.read_csv(ap_)
    mp_ = ab[ab.method.str.startswith("mean")].iloc[0]; am_ = ab[ab.method.str.startswith("ABMIL")].iloc[0]
    axj = fig.add_subplot(gs[2, 1]); lab(axj, "e")
    vals = [float(mp_.r), float(am_.r)]; r2s = [float(mp_.R2), float(am_.R2)]
    # horizontal dumbbell on a 0-0.8 r axis (no per-fold values exist locally: abmil_ccrcc.csv holds r/R2 only)
    axj.axvline(0, color=LGREY, lw=0.6, zorder=0)
    axj.plot(vals, [0, 0], color=GREY, lw=1.4, zorder=1, solid_capstyle="butt")
    axj.plot(vals[0], 0, "o", ms=7, mfc="white", mec=GREY, mew=1.3, zorder=3)
    axj.plot(vals[1], 0, "o", ms=7, mfc=ORG["CCRCC"], mec=ORG["CCRCC"], zorder=3)
    axj.text(vals[0], -0.22, f"mean-pool + PCA + ridge (main)\n$r$ = {vals[0]:.2f}, R² = {r2s[0]:.2f}",
             ha="center", va="top", fontsize=7, color=INK, linespacing=1.25)
    axj.text(vals[1], 0.22, f"attention MIL (ABMIL)\n$r$ = {vals[1]:.2f}, R² = {r2s[1]:.2f}",
             ha="center", va="bottom", fontsize=7, color=ORG["CCRCC"], linespacing=1.25)
    axj.set_xlim(-0.06, 0.8); axj.set_ylim(-1.0, 1.0); axj.set_yticks([])
    axj.spines["left"].set_visible(False)
    axj.set_xticks([0, 0.2, 0.4, 0.6, 0.8]); axj.tick_params(axis="x", labelsize=7)
    axj.set_xlabel("OOF Pearson $r$ (translation residual, same patient folds)", fontsize=7)
    axj.set_title("Not an artefact of mean-pooling (CCRCC)", fontsize=7.6, fontweight="bold", loc="left")
S.save_pub(fig, "Fig_master2")
print("wrote Fig_master2 (5 panels: enrichment, composition, nuclear correlates, "
      "KIRC replication, ABMIL). Clinical block -> make_fig_clinical.py; "
      "spatial maps -> make_fig_wsi_spatial_supp.py")
