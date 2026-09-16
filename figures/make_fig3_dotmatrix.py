#!/usr/bin/env python3
"""make_fig3_dotmatrix.py -- Figure 3 as a term-by-cohort dot matrix.

Reactome over-representation of the morphology-predictable residual proteins,
tested against each cohort's own MS-detected proteome (P0-5 fix; produced by
enrichment_tested_background.py -> figdata/enrichment_tested_background_full.csv).
Replaces the five stacked horizontal-bar panels of Fig3_enrichment.

Row selection (reproducible, see CANON / EXCLUDE below):
  1. nested or near-duplicate Reactome terms (parent/child pairs, or terms whose
     hit is the same ribosomal-protein core) are collapsed onto one canonical,
     recognisable member; the dot always shows the canonical term's OWN q;
  2. Reactome top-level umbrellas, promiscuous developmental sets and
     host-pathogen terms are dropped as uninformative rows;
  3. per cohort: the first six collapsed terms by BH q, extended to any term at
     collapsed rank 7-10 that is still overwhelming (q < 1e-8);
  4. rows = union across cohorts, grouped by family, strongest first.
Writes Fig3_enrichment.{svg,pdf,png,tiff} via mrstyle.save_pub.
"""
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory
import mrstyle as S

mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})

COH, ORGAN, ORG, FAM = S.COH, S.ORGAN, S.ORG, S.FAM
INK, GREY, LGREY = S.INK, S.GREY, S.LGREY
DD = "figdata"
Q_SIG, CAP, TOP_N, EXT_N, EXT_Q = 0.05, 30.0, 6, 10, 1e-8

# ---------------------------------------------------------------- term bookkeeping
# canonical term -> nested / near-duplicate members collapsed onto it
CANON = {
    "Cap-dependent Translation Initiation": [        # the ~70-protein ribosomal-subunit core
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

# Reactome top-level umbrellas, promiscuous developmental sets, host-pathogen terms
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
SHORT = {
    "Cap-dependent Translation Initiation": "Cap-dependent translation initiation",
    "Eukaryotic Translation Elongation": "Eukaryotic translation elongation",
    "SRP-dependent Cotranslational Protein Targeting To Membrane": "SRP-dependent cotranslational targeting",
    "Major Pathway Of rRNA Processing In Nucleolus And Cytosol": "rRNA processing (major pathway)",
    "Asparagine N-linked Glycosylation": "Asparagine N-linked glycosylation",
    "ER To Golgi Anterograde Transport": "ER-to-Golgi anterograde transport",
    "Membrane Trafficking": "Membrane trafficking",
    "mRNA Splicing - Major Pathway": "mRNA splicing (major pathway)",
    "Processing Of Capped Intron-Containing Pre-mRNA": "Processing of capped pre-mRNA",
    "Chaperonin-mediated Protein Folding": "Chaperonin-mediated protein folding",
    "Prefoldin Mediated Transfer Of Substrate To CCT/TriC": "Prefoldin substrate transfer to CCT/TriC",
    "Cooperation Of PDCL (PhLP1) And TRiC/CCT In G-protein Beta Folding": "PDCL–TRiC/CCT G-protein β folding",
    "Formation Of Tubulin Folding Intermediates By CCT/TriC": "Tubulin folding intermediates (CCT/TriC)",
    "Extracellular Matrix Organization": "Extracellular matrix organisation",
    "ECM Proteoglycans": "ECM proteoglycans",
    "Assembly Of Collagen Fibrils And Other Multimeric Structures": "Collagen fibril assembly",
    "SUMOylation": "SUMOylation",
    "SUMOylation Of Transcription Cofactors": "SUMOylation of transcription cofactors",
    "SUMOylation Of DNA Damage Response And Repair Proteins": "SUMOylation of DNA-repair proteins",
    "Transport Of Small Molecules": "Transport of small molecules",
    "Glucose Metabolism": "Glucose metabolism",
    "Class I Peroxisomal Membrane Protein Import": "Peroxisomal membrane protein import",
    "Neutrophil Degranulation": "Neutrophil degranulation",
    "RHO GTPase Effectors": "RHO GTPase effectors",
    "Chromatin Modifying Enzymes": "Chromatin-modifying enzymes",
    "Positive Epigenetic Regulation Of rRNA Expression": "Epigenetic activation of rRNA expression",
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
    return s if len(s) <= 40 else s[:39] + "…"


# ---------------------------------------------------------------- data + row selection
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

rows = sorted(set(sum(picked.values(), [])), key=lambda t: np.nanmin([q_of(t, c) for c in COH]))
fam_of = {t: FAMILY.get(t, "other") for t in rows}
for t in rows:
    if t not in FAMILY:
        print(f"  [note] '{t}' -> family 'other'")

# ---------------------------------------------------------------- layout bookkeeping
GAP = 0.55                                   # extra vertical space between families (row units)
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
print(f"{n} rows; picked per cohort: " + "; ".join(f"{c}={len(picked[c])}" for c in COH))


def dot_size(q):
    v = min(max(-np.log10(q), 0.0), CAP)
    return 8.0 + 92.0 * v / CAP                  # marker area (pt^2): 12 at q=0.05 ... 100 at cap


# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(5.2, 6.0))
AX_L, AX_R, AX_B, AX_T = 0.415, 0.775, 0.025, 0.885
ax = fig.add_axes([AX_L, AX_B, AX_R - AX_L, AX_T - AX_B])
ax.set_xlim(-0.5, len(COH) - 0.5); ax.set_ylim(ymin, ymax)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(axis="both", length=0)
ax.set_xticks([])
ax.set_yticks(ys); ax.set_yticklabels([short(t) for t in ordered], fontsize=6.1, color=INK)
ax.tick_params(axis="y", pad=3)

# row guides + family separators
for yy in ys:
    ax.plot([-0.5, len(COH) - 0.5], [yy, yy], color=LGREY, lw=0.35, ls=(0, (1, 2)), zorder=0)
fam_present = [f for f in FAMS if f in spans]
for f_hi, f_lo in zip(fam_present[:-1], fam_present[1:]):
    ysep = (spans[f_hi][1] + spans[f_lo][0]) / 2.0
    ax.plot([-0.5, len(COH) - 0.5], [ysep, ysep], color=LGREY, lw=0.6, zorder=0, clip_on=False)

# dots
for j, c in enumerate(COH):
    for t, yy in zip(ordered, ys):
        q = q_of(t, c)
        if np.isnan(q) or q >= Q_SIG:
            ax.scatter(j, yy, s=7, facecolor="white", edgecolor=GREY, linewidth=0.5, zorder=2)
        else:
            ax.scatter(j, yy, s=dot_size(q), facecolor=FAMC[fam_of[t]], edgecolor="white",
                       linewidth=0.35, zorder=3)

# column headers: cohort (bold, organ colour) over organ (grey)
bt = blended_transform_factory(ax.transData, ax.transAxes)
for j, c in enumerate(COH):
    ax.annotate(c, xy=(j, 1.0), xycoords=bt, xytext=(0, 11), textcoords="offset points",
                ha="center", va="bottom", fontsize=6.8, fontweight="bold", color=ORG[c])
    ax.annotate(ORGAN[c], xy=(j, 1.0), xycoords=bt, xytext=(0, 2.5), textcoords="offset points",
                ha="center", va="bottom", fontsize=6.0, color=GREY)

# family tags (coloured band with rotated white label) left of the term labels
axb = fig.add_axes([0.012, AX_B, 0.026, AX_T - AX_B]); axb.set_axis_off()
axb.set_xlim(0, 1); axb.set_ylim(ymin, ymax)
for f, (top, bot) in spans.items():
    axb.add_patch(Rectangle((0.15, bot - 0.42), 0.7, top - bot + 0.84, facecolor=FAMC[f],
                            edgecolor="none"))
    axb.text(0.5, (top + bot) / 2.0, FAMTAG[f], rotation=90, ha="center", va="center",
             fontsize=5.6, fontweight="bold", color="white")

# legends (right margin): dot size, then family colour
axl = fig.add_axes([0.79, AX_B, 0.205, AX_T - AX_B]); axl.set_axis_off()
axl.set_xlim(0, 1); axl.set_ylim(0, 1)
yl, step = 0.985, 0.047
axl.text(0.0, yl, "adjusted P (q)", fontsize=6.4, fontweight="bold", color=INK, va="top"); yl -= step * 0.9
for q, lab in [(0.05, "0.05"), (1e-5, r"$10^{-5}$"), (1e-10, r"$10^{-10}$"), (1e-20, r"$10^{-20}$")]:
    axl.scatter(0.16, yl, s=dot_size(q), facecolor=GREY, edgecolor="white", linewidth=0.35, clip_on=False)
    axl.text(0.36, yl, lab, fontsize=6.2, color=INK, va="center"); yl -= step
axl.scatter(0.16, yl, s=7, facecolor="white", edgecolor=GREY, linewidth=0.5)
axl.text(0.36, yl, "≥ 0.05 (n.s.)", fontsize=6.2, color=INK, va="center"); yl -= step * 0.85
axl.text(0.0, yl, "dot area: −log$_{10}$ q\n(capped at 30)", fontsize=5.4, color=GREY, va="top",
         linespacing=1.15); yl -= step * 1.9
axl.text(0.0, yl, "family", fontsize=6.4, fontweight="bold", color=INK, va="top"); yl -= step * 0.9
for f in fam_present:
    axl.scatter(0.16, yl, s=26, facecolor=FAMC[f], edgecolor="white", linewidth=0.35)
    axl.text(0.36, yl, FAMLEG[f], fontsize=6.0, color=INK, va="center"); yl -= step * 0.9

# title + method line
fig.text(0.012, 0.992, "Organ-specific enrichment of the residual proteins (tested-proteome background)",
         fontsize=8, fontweight="bold", color=INK, ha="left", va="top")
fig.text(0.012, 0.962, "Reactome ORA of morphology-predictable residual proteins; hypergeometric test, "
         "BH-adjusted; background = each cohort's MS-detected proteome",
         fontsize=5.8, color=GREY, ha="left", va="top")

S.save_pub(fig, "Fig3_enrichment")
plt.close(fig)

# ---------------------------------------------------------------- console report
print("\nrows by family (q per cohort; '.' = q >= 0.05 or absent):")
for f in fam_present:
    print(f"  [{f}]")
    for t in ordered:
        if fam_of[t] != f:
            continue
        cells = "  ".join(f"{c}={q_of(t, c):.0e}" if q_of(t, c) < Q_SIG else f"{c}=." for c in COH)
        print(f"    {short(t):<42} {cells}")
print("\npicked per cohort (collapsed ranking):")
for c in COH:
    print(f"  {c}: " + " | ".join(short(t) for t in picked[c]))
print("wrote Fig3_enrichment.{svg,pdf,png,tiff}")
