#!/usr/bin/env python3
"""Composite figure "what the residual is" (Fig_biology).
a  term-by-cohort Reactome enrichment dot matrix (plotting code copied from make_fig3_dotmatrix.py,
   data figdata/enrichment_tested_background_full.csv; legends moved below the matrix)
b-e single-protein vignettes RPL35A / SSR3 (vignette_plot.py, figdata/vignette_data.csv)
f  representative H&E tiles, highest vs lowest translation-residual patients
   (make_realdata_figs2.py Fig11_tiles_hilo, curated tile indices copied verbatim)
No source script is imported (they all compute at module level); all numbers are recomputed from
the saved figdata files with the same code."""
import glob
import textwrap
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory
from scipy import stats as st
import mrstyle as S

mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     # Sub/superscripts are set as real Unicode glyphs rather than mathtext, because
                     # mathtext shrinks them to 0.7x the body size (4.9 pt out of 7 pt) and the
                     # journal's technical check measures every span. Arial carries the superscript
                     # digits but not U+2080-2089 / U+207B, so DejaVu Sans (shipped with matplotlib)
                     # is appended purely as a glyph fallback; all Latin text stays Arial.
                     "font.family": ["Arial", "DejaVu Sans"]})

COH, ORGAN, ORG, FAM = S.COH, S.ORGAN, S.ORG, S.FAM
INK, GREY, LGREY = S.INK, S.GREY, S.LGREY
DD = "figdata"


def lab(ax, s, x=-0.18, y=1.05):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


# =====================================================================================
# a : enrichment dot matrix -- term bookkeeping copied verbatim from make_fig3_dotmatrix.py
# =====================================================================================
Q_SIG, CAP, TOP_N, EXT_N, EXT_Q = 0.05, 30.0, 6, 10, 1e-8
CANON = {
    "Cap-dependent Translation Initiation": [
        "L13a-mediated Translational Silencing Of Ceruloplasmin Expression",
        "GTP Hydrolysis And Joining Of 60S Ribosomal Subunit",
        "Formation Of A Pool Of Free 40S Subunits", "Translation",
        "Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency",
        "Translation Initiation Complex Formation",
        "mRNA Activation Upon Binding Of Cap-Binding Complex And eIFs, Subsequent Binding To 43S",
        "Ribosomal Scanning And Start Codon Recognition",
        "Formation Of Ternary Complex, And Subsequently, 43S Complex",
        "Selenocysteine Synthesis", "Selenoamino Acid Metabolism", "Viral mRNA Translation",
        "Nonsense Mediated Decay (NMD) Independent Of Exon Junction Complex (EJC)",
        "Nonsense Mediated Decay (NMD) Enhanced By Exon Junction Complex (EJC)",
        "Nonsense-Mediated Decay (NMD)", "Eukaryotic Translation Termination",
        "Cellular Response To Starvation", "SARS-CoV-2 Modulates Host Translation Machinery",
        "Influenza Viral RNA Transcription And Replication", "Influenza Infection",
        "Regulation Of Expression Of SLITs And ROBOs", "Signaling By ROBO Receptors"],
    "Eukaryotic Translation Elongation": ["Peptide Chain Elongation"],
    "Major Pathway Of rRNA Processing In Nucleolus And Cytosol": [
        "rRNA Processing In Nucleus And Cytosol", "rRNA Processing",
        "rRNA Modification In Nucleus And Cytosol"],
    "mRNA Splicing - Major Pathway": ["mRNA Splicing"],
    "ER To Golgi Anterograde Transport": ["Transport To Golgi And Subsequent Modification",
                                          "COPI-mediated Anterograde Transport",
                                          "COPII-mediated Vesicle Transport"],
    "Membrane Trafficking": ["Vesicle-mediated Transport"],
    "Transport Of Small Molecules": ["SLC-mediated Transmembrane Transport",
                                     "Disorders Of Transmembrane Transporters"],
    "Chaperonin-mediated Protein Folding": ["Protein Folding"],
    "Prefoldin Mediated Transfer Of Substrate To CCT/TriC": [
        "Cooperation Of Prefoldin And TriC/CCT In Actin And Tubulin Folding",
        "Folding Of Actin By CCT/TriC"],
    "Formation Of Tubulin Folding Intermediates By CCT/TriC": ["Post-chaperonin Tubulin Folding Pathway"],
    "SUMOylation": ["SUMO E3 Ligases SUMOylate Target Proteins"],
    "RHO GTPase Effectors": ["Signaling By Rho GTPases",
                             "Signaling By Rho GTPases, Miro GTPases And RHOBTB3", "RHO GTPase Cycle"],
    "Chromatin Modifying Enzymes": ["HDACs Deacetylate Histones", "Chromatin Organization"],
    "Positive Epigenetic Regulation Of rRNA Expression": [
        "ERCC6 (CSB) And EHMT2 (G9a) Positively Regulate rRNA Expression",
        "Epigenetic Regulation Of Gene Expression", "RNA Polymerase I Transcription",
        "RNA Polymerase I Promoter Clearance", "RNA Polymerase I Transcription Initiation"],
    "Glucose Metabolism": ["Gluconeogenesis"],
    "Assembly Of Collagen Fibrils And Other Multimeric Structures": ["Collagen Formation"],
}
ALIAS = {a: k for k, v in CANON.items() for a in v}
EXCLUDE = {
    "Disease", "Infectious Disease", "Metabolism Of RNA", "Metabolism Of Proteins",
    "Post-translational Protein Modification", "Gene Expression (Transcription)",
    "RNA Polymerase II Transcription", "Cellular Responses To Stress", "Cellular Responses To Stimuli",
    "Innate Immune System", "Immune System", "Metabolism", "Signal Transduction",
    "Protein Localization", "Metabolism Of Amino Acids And Derivatives",
    "Axon Guidance", "Nervous System Development", "Developmental Biology",
    "Interactions Of Vpr With Host Cellular Proteins", "Leishmania Infection", "SARS-CoV Infections",
    "SARS-CoV-2-host Interactions", "Potential Therapeutics For SARS", "Late SARS-CoV-2 Infection Events",
    "Viral Messenger RNA Synthesis", "HIV Infection", "Host Interactions Of HIV Factors",
}
FAMILY = {
    "Cap-dependent Translation Initiation": "translation",
    "Eukaryotic Translation Elongation": "translation",
    "SRP-dependent Cotranslational Protein Targeting To Membrane": "translation",
    "Major Pathway Of rRNA Processing In Nucleolus And Cytosol": "translation",
    "Asparagine N-linked Glycosylation": "secretion",
    "ER To Golgi Anterograde Transport": "secretion",
    "Membrane Trafficking": "secretion",
    "mRNA Splicing - Major Pathway": "splicing",
    "Processing Of Capped Intron-Containing Pre-mRNA": "splicing",
    "Chaperonin-mediated Protein Folding": "folding",
    "Prefoldin Mediated Transfer Of Substrate To CCT/TriC": "folding",
    "Cooperation Of PDCL (PhLP1) And TRiC/CCT In G-protein Beta Folding": "folding",
    "Formation Of Tubulin Folding Intermediates By CCT/TriC": "folding",
    "Extracellular Matrix Organization": "ecm",
    "ECM Proteoglycans": "ecm",
    "Assembly Of Collagen Fibrils And Other Multimeric Structures": "ecm",
    "SUMOylation": "ptm",
    "SUMOylation Of Transcription Cofactors": "ptm",
    "SUMOylation Of DNA Damage Response And Repair Proteins": "ptm",
}
# Row labels are abbreviated so that the widest one fits the ~110 pt label gutter at 7 pt (they
# were tuned for 5.8 pt); the caption spells the Reactome terms out in full.
SHORT = {
    "Cap-dependent Translation Initiation": "Cap-dependent transl. initiation",
    "Eukaryotic Translation Elongation": "Translation elongation",
    "SRP-dependent Cotranslational Protein Targeting To Membrane": "SRP cotranslational targeting",
    "Major Pathway Of rRNA Processing In Nucleolus And Cytosol": "rRNA processing (major)",
    "Asparagine N-linked Glycosylation": "N-linked glycosylation",
    "ER To Golgi Anterograde Transport": "ER-to-Golgi transport",
    "Membrane Trafficking": "Membrane trafficking",
    "mRNA Splicing - Major Pathway": "mRNA splicing (major)",
    "Processing Of Capped Intron-Containing Pre-mRNA": "Processing of capped pre-mRNA",
    "Chaperonin-mediated Protein Folding": "Chaperonin-mediated folding",
    "Prefoldin Mediated Transfer Of Substrate To CCT/TriC": "Prefoldin–CCT/TriC transfer",
    "Cooperation Of PDCL (PhLP1) And TRiC/CCT In G-protein Beta Folding": "PDCL–TRiC/CCT β folding",
    "Formation Of Tubulin Folding Intermediates By CCT/TriC": "Tubulin folding (CCT/TriC)",
    "Extracellular Matrix Organization": "ECM organisation",
    "ECM Proteoglycans": "ECM proteoglycans",
    "Assembly Of Collagen Fibrils And Other Multimeric Structures": "Collagen fibril assembly",
    "SUMOylation": "SUMOylation",
    "SUMOylation Of Transcription Cofactors": "SUMOylation of TF cofactors",
    "SUMOylation Of DNA Damage Response And Repair Proteins": "SUMOylation, DNA repair",
    "Transport Of Small Molecules": "Transport of small molecules",
    "Glucose Metabolism": "Glucose metabolism",
    "Class I Peroxisomal Membrane Protein Import": "Peroxisomal protein import",
    "Neutrophil Degranulation": "Neutrophil degranulation",
    "RHO GTPase Effectors": "RHO GTPase effectors",
    "Chromatin Modifying Enzymes": "Chromatin-modifying enzymes",
    "Positive Epigenetic Regulation Of rRNA Expression": "Epigenetic rRNA activation",
}
FAMS = ["translation", "secretion", "splicing", "folding", "ecm", "ptm", "other"]
FAMC = dict(FAM); FAMC["other"] = INK
FAMTAG = {"translation": "translation", "secretion": "secretion", "splicing": "splicing",
          "folding": "folding", "ecm": "ECM", "ptm": "PTM", "other": "other"}
FAMLEG = {"translation": "translation", "secretion": "secretion / ER–Golgi", "splicing": "mRNA splicing",
          "folding": "folding / chaperonin", "ecm": "ECM / matrisome", "ptm": "PTM (SUMOylation)",
          "other": "other"}


def short(t):
    s = SHORT.get(t, t)
    return s if len(s) <= 32 else s[:31] + "…"


df = pd.read_csv(f"{DD}/enrichment_tested_background_full.csv")
df["term"] = df["Term"].str.replace(r"\s+R-HSA-\d+$", "", regex=True)
df["q"] = df["Adjusted P-value"].astype(float)
df["canon"] = df["term"].map(lambda t: ALIAS.get(t, t))
qtab = (df[df["term"] == df["canon"]]
        .pivot_table(index="term", columns="cohort", values="q", aggfunc="min")
        .reindex(columns=COH))


def q_of(term, c):
    return float(qtab.loc[term, c]) if term in qtab.index and pd.notna(qtab.loc[term, c]) else np.nan


picked = {}
for c in COH:
    seen = []
    for _, r in df[df["cohort"] == c].sort_values("q").iterrows():
        k = r["canon"]
        if r["term"] in EXCLUDE or k in EXCLUDE or k in seen:
            continue
        seen.append(k)
        if len(seen) >= EXT_N:
            break
    picked[c] = seen[:TOP_N] + [k for k in seen[TOP_N:EXT_N] if q_of(k, c) < EXT_Q]

# union in first-appearance order (deterministic; the standalone uses set(), whose order for exact
# q ties -- e.g. rRNA processing vs cap-dependent initiation in LUAD, both q=1.4e-26 -- varies by hash seed)
rows = list(dict.fromkeys(sum(picked.values(), [])))
rows = sorted(rows, key=lambda t: np.nanmin([q_of(t, c) for c in COH]))
fam_of = {t: FAMILY.get(t, "other") for t in rows}

GAP = 0.55
ordered, ys, spans = [], [], {}
y = 0.0
for f in FAMS:
    rf = [t for t in rows if fam_of[t] == f]
    if not rf:
        continue
    top = y
    for t in rf:
        ordered.append(t); ys.append(y); y -= 1.0
    spans[f] = (top, y + 1.0)
    y -= GAP
n = len(ordered)
ymin, ymax = min(ys) - 0.5, 0.5
fam_present = [f for f in FAMS if f in spans]
print(f"a: {n} rows; picked per cohort: " + "; ".join(f"{c}={len(picked[c])}" for c in COH))


def dot_size(q):
    v = min(max(-np.log10(q), 0.0), CAP)
    return 8.0 + 92.0 * v / CAP                  # marker area (pt^2): 12 at q=0.05 ... 100 at cap


# =====================================================================================
# figure
# =====================================================================================
FW, FH = 7.2, 8.4
fig = plt.figure(figsize=(FW, FH))

# ---------- a : matrix (left column, ~46 % width, full height) ----------
# at 7 pt the widest row label is ~105 pt and the widest column header ("Pancreas") ~30 pt, so the
# label gutter needs ~110 pt (0.213 of the 518 pt canvas) and the five columns need >=26 pt of
# pitch; the matrix therefore starts further left and ends further right than at 5.8 pt, and the
# bottom edge rises to give the legend block below room for 7 pt rows.
AX_L, AX_R, AX_B, AX_T = 0.243, 0.494, 0.180, 0.892
ax = fig.add_axes([AX_L, AX_B, AX_R - AX_L, AX_T - AX_B])
ax.set_xlim(-0.5, len(COH) - 0.5); ax.set_ylim(ymin, ymax)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(axis="both", length=0)
ax.set_xticks([])
ax.set_yticks(ys); ax.set_yticklabels([short(t) for t in ordered], fontsize=7.0, color=INK)
ax.tick_params(axis="y", pad=3)
for yy in ys:
    ax.plot([-0.5, len(COH) - 0.5], [yy, yy], color=LGREY, lw=0.35, ls=(0, (1, 2)), zorder=0)
for f_hi, f_lo in zip(fam_present[:-1], fam_present[1:]):
    ysep = (spans[f_hi][1] + spans[f_lo][0]) / 2.0
    ax.plot([-0.5, len(COH) - 0.5], [ysep, ysep], color=LGREY, lw=0.6, zorder=0, clip_on=False)
for j, c in enumerate(COH):
    for t, yy in zip(ordered, ys):
        q = q_of(t, c)
        if np.isnan(q) or q >= Q_SIG:
            ax.scatter(j, yy, s=7, facecolor="white", edgecolor=GREY, linewidth=0.5, zorder=2)
        else:
            ax.scatter(j, yy, s=dot_size(q), facecolor=FAMC[fam_of[t]], edgecolor="white",
                       linewidth=0.35, zorder=3)
bt = blended_transform_factory(ax.transData, ax.transAxes)
for j, c in enumerate(COH):
    ax.annotate(c, xy=(j, 1.0), xycoords=bt, xytext=(0, 11), textcoords="offset points",
                ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=ORG[c])
    ax.annotate(ORGAN[c], xy=(j, 1.0), xycoords=bt, xytext=(0, 2.5), textcoords="offset points",
                ha="center", va="bottom", fontsize=7.0, color=GREY)

# family tags (coloured band with rotated white label) left of the term labels
axb = fig.add_axes([0.004, AX_B, 0.019, AX_T - AX_B]); axb.set_axis_off()
axb.set_xlim(0, 1); axb.set_ylim(ymin, ymax)
for f, (top, bot) in spans.items():
    axb.add_patch(Rectangle((0.15, bot - 0.42), 0.7, top - bot + 0.84, facecolor=FAMC[f], edgecolor="none"))
    axb.text(0.5, (top + bot) / 2.0, FAMTAG[f], rotation=90, ha="center", va="center",
             fontsize=7.0, fontweight="bold", color="white")

# legends moved BELOW the matrix: dot-size (left) and family colour (right).
# The block is taller and the row pitch larger than at 5.8 pt: seven 7 pt rows need >= 8.4 pt of
# pitch each, and the family column now steps at the full pitch instead of 0.86 of it.
axl = fig.add_axes([0.03, 0.012, 0.46, 0.145]); axl.set_axis_off()
axl.set_xlim(0, 1); axl.set_ylim(0, 1)
yl, step = 0.97, 0.118
axl.text(0.0, yl, "adjusted P (q)", fontsize=7.0, fontweight="bold", color=INK, va="top"); yl -= step * 1.05
# Unicode superscripts, not mathtext: mathtext would set the exponent at 0.7 x 7 = 4.9 pt.
for q, lq in [(0.05, "0.05"), (1e-5, "10⁻⁵"), (1e-10, "10⁻¹⁰"),
              (1e-20, "10⁻²⁰")]:
    axl.scatter(0.055, yl, s=dot_size(q), facecolor=GREY, edgecolor="white", linewidth=0.35, clip_on=False)
    axl.text(0.135, yl, lq, fontsize=7.0, color=INK, va="center"); yl -= step
axl.scatter(0.055, yl, s=7, facecolor="white", edgecolor=GREY, linewidth=0.5)
axl.text(0.135, yl, "≥ 0.05 (n.s.)", fontsize=7.0, color=INK, va="center"); yl -= step * 0.95
axl.text(0.0, yl, "dot area: −log₁₀ q (capped at 30)", fontsize=7.0, color=GREY, va="center")
yl = 0.97
axl.text(0.50, yl, "family", fontsize=7.0, fontweight="bold", color=INK, va="top"); yl -= step * 1.05
for f in fam_present:
    axl.scatter(0.555, yl, s=26, facecolor=FAMC[f], edgecolor="white", linewidth=0.35)
    axl.text(0.625, yl, FAMLEG[f], fontsize=7.0, color=INK, va="center"); yl -= step

fig.text(0.012, 0.990, "a", fontsize=10, fontweight="bold", va="top", ha="left")
fig.text(0.055, 0.990, "Organ-specific enrichment of the residual proteins\n(tested-proteome background)",
         fontsize=7.6, fontweight="bold", color=INK, ha="left", va="top", linespacing=1.15)
# wrapped so the last (lowest) line, which sits at the height of the column headers, is short
# enough to stop ~30 pt left of "CCRCC"
fig.text(0.055, 0.958, "Reactome ORA of morphology-predictable residual proteins;\n"
         "hypergeometric test, BH-adjusted; background = each cohort's\nMS-detected proteome",
         fontsize=7.0, color=GREY, ha="left", va="top", linespacing=1.15)

# ---------- b-e : single-protein vignettes (right column, top ~58 %) ----------
vd = pd.read_csv(f"{DD}/vignette_data.csv")
genes = list(dict.fromkeys(vd["gene"]))          # preserve order: RPL35A, SSR3
vcol = [S.ORG["CCRCC"], S.ORG["UCEC"]]


def r2(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return 1 - np.sum((a - b) ** 2) / np.sum((a - a.mean()) ** 2)


gsv = fig.add_gridspec(2, 2, left=0.605, right=0.985, top=0.935, bottom=0.420, hspace=0.50, wspace=0.55)
axes_v = np.array([[fig.add_subplot(gsv[i, j]) for j in range(2)] for i in range(2)])
letters = [["b", "c"], ["d", "e"]]
vig_numbers = {}
for gi, g in enumerate(genes):
    d = vd[vd["gene"] == g]
    c = vcol[gi % len(vcol)]
    # left: protein vs mRNA
    axv = axes_v[gi, 0]; lab(axv, letters[gi][0], x=-0.42, y=1.02)
    axv.scatter(d["mrna"], d["protein"], s=12, color=c, alpha=0.7, edgecolors="none")
    b1, b0 = np.polyfit(d["mrna"], d["protein"], 1)
    xr = np.array([d["mrna"].min(), d["mrna"].max()])
    axv.plot(xr, b1 * xr + b0, color=S.GREY, lw=1.0, ls="--")
    rr = st.pearsonr(d["mrna"], d["protein"])[0]
    # gene label: the 2-line block (0.91 in wide) cannot sit anywhere inside the 1.07-in axes without
    # covering points (b's cloud fills the panel), so the name + relation go on the panel-letter
    # baseline above the axes and only the coefficient sits in the empty lower-left corner
    # (b: nothing below y=-0.5 for x<7.7; d: nothing below y=-2 for x<6.3)
    # x=-0.24: starts ~0.2 in right of the panel letter (same letter-title gap as panels a and f)
    axv.text(-0.24, 1.02, f"{g}, mRNA$\\rightarrow$protein", transform=axv.transAxes,
             fontsize=7.0, ha="left", va="bottom", fontweight="bold", color=c)
    axv.text(0.04, 0.04, f"$r$={rr:.2f}", transform=axv.transAxes,
             fontsize=7.0, ha="left", va="bottom", fontweight="bold", color=c)
    axv.set_xlabel("mRNA (log₂ TPM)", fontsize=7.0); axv.set_ylabel("protein (measured)", fontsize=7.0)
    axv.tick_params(labelsize=7.0)
    vig_numbers[f"{g} mRNA->protein r"] = round(float(rr), 2)
    # right: residual recovered by morphology
    axv = axes_v[gi, 1]; lab(axv, letters[gi][1], x=-0.42, y=1.02)
    axv.scatter(d["residual"], d["morph_pred"], s=12, color=c, alpha=0.7, edgecolors="none")
    lo = min(d["residual"].min(), d["morph_pred"].min()); hi = max(d["residual"].max(), d["morph_pred"].max())
    axv.plot([lo, hi], [lo, hi], color=S.GREY, lw=0.8, ls=":")
    b1, b0 = np.polyfit(d["residual"], d["morph_pred"], 1)
    xr = np.array([d["residual"].min(), d["residual"].max()])
    axv.plot(xr, b1 * xr + b0, color=S.INK, lw=1.0)
    incr = r2(d["residual"], d["morph_pred"])
    rr = st.pearsonr(d["residual"], d["morph_pred"])[0]
    # "morphology recovers residual" is 1.13 in wide at 6 pt (wider than the 1.07-in axes), so any
    # in-axes placement crosses the identity line; it goes on the panel-letter baseline above the
    # axes and the two statistics sit as a narrow block in the empty lower-right corner
    axv.text(-0.24, 1.02, "morphology recovers residual", transform=axv.transAxes,
             fontsize=7.0, ha="left", va="bottom", color=S.INK)
    # (0.98, 0.02) rather than (0.97, 0.04): at (0.97, 0.04) one RPL35A point sits 1 pt from the block
    axv.text(0.98, 0.02, f"$r$={rr:.2f}\n(morphology\nR²={incr:.2f})",
             transform=axv.transAxes, fontsize=7.0, ha="right", va="bottom", color=S.INK, linespacing=1.15)
    axv.set_xlabel("measured residual\n(protein $-$ f(mRNA))", fontsize=7.0)
    axv.set_ylabel("morphology-predicted\nresidual", fontsize=7.0)
    axv.tick_params(labelsize=7.0)
    vig_numbers[f"{g} residual r"] = round(float(rr), 2)
    vig_numbers[f"{g} morphology R2"] = round(float(incr), 2)
fig.text(0.535, 0.990, "Single-protein vignettes: mRNA misses it,\nmorphology reads the residual",
         fontsize=7.6, fontweight="bold", color=INK, ha="left", va="top", linespacing=1.15)

# ---------- f : H&E tile montage (right column, bottom ~42 %) ----------
# curated cellular tiles, indices copied verbatim from make_realdata_figs2.py (Fig11_tiles_hilo)
hi = [f"{DD}/tile_hi_{i}_{j}.png" for i, j in [(1, 0), (2, 0), (2, 1), (3, 1)]]
lo = [f"{DD}/tile_lo_{i}_{j}.png" for i, j in [(0, 1), (1, 1), (3, 0), (3, 1)]]
hi = [f for f in hi if glob.glob(f)]
lo = [f for f in lo if glob.glob(f)]
# square cells: 2 rows x 4 columns of 256x256 tiles, so the rows touch with only a hairline gap.
# gridspec w/hspace are fractions of the average cell size: total width = (4 + 3*wspace) * cell,
# total height = (2 + hspace) * cell, cell = column width / (4 + 3*wspace)
T_L, T_R, T_TOP, T_SP = 0.605, 0.985, 0.290, 0.06
T_CELL_IN = (T_R - T_L) * FW / (4 + 3 * T_SP)                 # ~0.655 in per tile
T_BOT = T_TOP - (2 + T_SP) * T_CELL_IN / FH
gst = fig.add_gridspec(2, 4, left=T_L, right=T_R, top=T_TOP, bottom=T_BOT, hspace=T_SP, wspace=T_SP)
print(f"f: tile cell = {T_CELL_IN:.3f} in, grid top={T_TOP:.3f} bottom={T_BOT:.3f}")
axes_t = np.array([[fig.add_subplot(gst[i, j]) for j in range(4)] for i in range(2)])
for j, f in enumerate(hi[:4]):
    axes_t[0, j].imshow(mpimg.imread(f)); axes_t[0, j].axis("off")
for j, f in enumerate(lo[:4]):
    axes_t[1, j].imshow(mpimg.imread(f)); axes_t[1, j].axis("off")
for i in range(2):
    for j in range(4):
        axes_t[i, j].axis("off")
axes_t[0, 0].text(-0.08, 0.5, "highest\nresidual", transform=axes_t[0, 0].transAxes, fontsize=7.0,
                  va="center", ha="right", fontweight="bold", color=S.ORG["CCRCC"])
axes_t[1, 0].text(-0.08, 0.5, "lowest\nresidual", transform=axes_t[1, 0].transAxes, fontsize=7.0,
                  va="center", ha="right", fontweight="bold", color=S.GREY)
fig.text(0.535, 0.335, "f", fontsize=10, fontweight="bold", va="top", ha="left")
fig.text(0.575, 0.335, "Representative H&E (CCRCC): highest vs lowest\ntranslation-residual patients",
         fontsize=7.6, fontweight="bold", color=INK, ha="left", va="top", linespacing=1.15)

S.save_pub(fig, "Fig_biology")
plt.close(fig)
print("vignette numbers:", vig_numbers)
print("tiles hi:", hi); print("tiles lo:", lo)
print("wrote Fig_biology.{svg,pdf,png,tiff}")
