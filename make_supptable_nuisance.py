#!/usr/bin/env python3
"""SuppTable_nuisance.tex -- the acquisition/stain nuisance block, and the UCEC
site-effect diagnosis it turned up. Written with the Write tool, not a heredoc."""
import numpy as np
import pandas as pd

DD = "figures/figdata"
ORDER = {"CCRCC": 0, "LUAD": 1, "UCEC": 2, "GBM": 3, "PDAC": 4}
NICE = {"translation_meas": "Transl.", "secretion_ER_meas": "ER-secr."}

blk = pd.read_csv(f"{DD}/nuisance_block.csv")
cca = pd.read_csv(f"{DD}/nuisance_block_cca.csv").set_index("cohort")
wch = pd.read_csv(f"{DD}/nuisance_which.csv")

rows, prev = [], None
for _, r in blk.sort_values(["cohort", "score"],
                            key=lambda s: s.map(ORDER) if s.name == "cohort" else s).iterrows():
    if prev is not None and r.cohort != prev:
        rows.append(r"\midrule")
    prev = r.cohort
    c = cca.loc[r.cohort]
    ratio = f"{100*r.r2_nuisance/r.r2_morphology:.0f}" if r.r2_morphology > 0.02 and r.r2_nuisance > 0 else "--"
    rows.append(f"{r.cohort} & {NICE[r.score]} & {int(r.n)} & "
                f"{c.cca1:.2f} & {c.cca1_null_median:.2f} & {100*c.r2_weighted_20pc:.0f} & "
                f"{r.r2_nuisance:+.3f} & {r.r2_morphology:+.3f} & {ratio} " + r"\\")

# which variable carries UCEC
u = wch[(wch.cohort == "UCEC") & (wch.score == "translation_meas")].sort_values("r2_alone", ascending=False)
urows = [f"{r.variable.replace('_', chr(92)+'_')} & {r.r2_alone:+.3f} & {r.r2_without:+.3f} & "
         f"{r.spearman_with_outcome:+.2f} " + r"\\" for _, r in u.iterrows()]

caption = (
    r"An acquisition-and-stain nuisance block does not reproduce the morphology "
    r"signal (\S\ref{sec:encoder}). Because an ImageNet ResNet-50 recovers 79\% of the "
    r"UNI-significant set, a natural objection is that the embedding keys on gross stain "
    r"and acquisition statistics rather than morphology. Per patient we formed a five-term "
    r"block carrying no tumour morphology---$\log_{10}$ total tiles, number of slides, and "
    r"the slide-mean nuclear density, chromatin optical density and eosin fraction from the "
    r"hand-crafted morphometry---and asked what it can do. CCA$_1$, first canonical "
    r"correlation with the 20 morphology principal components, against the median of a "
    r"label-permutation null (200 permutations; the observed value exceeds the null in "
    r"every cohort, $p=0.005$, so the alignment is real and not a consequence of comparing "
    r"5 variables with 20 components at $n\approx100$). \% PC var., the share of 20-PC "
    r"variance the block predicts, cross-validated. The last three columns are the point: "
    r"out-of-fold $R^2$ for the measured residual score from the nuisance block alone, from "
    r"the 20 morphology components, and their ratio. The block is strongly aligned with the "
    r"embedding yet predicts the residual in only one cohort; in eight of ten "
    r"cohort-by-score combinations it reaches $R^2\leq0.08$ and five are negative. Stain "
    r"statistics are therefore present in the embedding but are not what carries the "
    r"residual, and the ResNet-50 result cannot be attributed to them. The lower block "
    r"resolves the exception: in UCEC the predictive term is the number of tumour "
    r"slides embedded per patient, not any stain variable."
)

tex = [
    r"\begin{table}[htbp]\centering\small",
    r"\caption{" + caption + r"}",
    r"\label{tab:nuisance}",
    r"\resizebox{\textwidth}{!}{\begin{tabular}{llrrrrrrr}",
    r"\toprule",
    r"& & & \multicolumn{2}{c}{CCA$_1$} & & \multicolumn{3}{c}{OOF $R^2$ for the residual score} \\",
    r"\cmidrule(lr){4-5}\cmidrule(lr){7-9}",
    r"Cohort & Score & $n$ & obs.\ & null & \% PC var.\ & nuisance & morphology & ratio \% \\",
    r"\midrule",
]
tex += rows
tex += [
    r"\midrule",
    r"\multicolumn{9}{l}{\emph{UCEC translation: which term of the block carries it}} \\",
    r"\multicolumn{2}{l}{variable} & \multicolumn{2}{r}{alone} & \multicolumn{2}{r}{block without it} & "
    r"\multicolumn{3}{r}{Spearman with residual} \\",
]
for ur in urows:
    v, alone, without, rho = ur.replace(r" \\", "").split(" & ")
    tex.append(rf"\multicolumn{{2}}{{l}}{{{v}}} & \multicolumn{{2}}{{r}}{{{alone}}} & "
               rf"\multicolumn{{2}}{{r}}{{{without}}} & \multicolumn{{3}}{{r}}{{{rho}}} \\")
tex += [r"\bottomrule", r"\end{tabular}}", r"\end{table}"]

with open("SuppTable_nuisance.tex", "w", encoding="utf8") as f:
    f.write("\n".join(tex) + "\n")
print(f"wrote SuppTable_nuisance.tex ({len(rows)} block rows, {len(urows)} UCEC rows)")
