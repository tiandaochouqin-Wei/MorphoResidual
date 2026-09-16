#!/usr/bin/env python3
"""SuppTable_composition_path.tex -- composition covariates not derived from the mRNA
baseline: pathologist-read (PDAC) and published deconvolution (CCRCC/UCEC/GBM/PDAC).
Written with the Write tool, never a heredoc."""
import numpy as np
import pandas as pd

DD = "figures/figdata"
d = pd.read_csv(f"{DD}/composition_pathologist.csv")
d = d[d.n.fillna(0) > 0]
ag = pd.read_csv(f"{DD}/composition_proxy_agreement.csv")

NICE = {"translation": "Transl.", "secretion_ER": "ER-secr."}
ORDER = {"CCRCC": 0, "LUAD": 1, "UCEC": 2, "GBM": 3, "PDAC": 4}
d = d.sort_values(["cohort", "covariates", "score"],
                  key=lambda s: s.map(ORDER) if s.name == "cohort" else s)

rows, prev = [], None
for _, r in d.iterrows():
    if prev is not None and r.cohort != prev:
        rows.append(r"\midrule")
    prev = r.cohort
    cov = r.covariates.replace("&", r"\&").replace("+", "${+}$")
    rows.append(f"{r.cohort} & {NICE[r.score]} & {cov} & {int(r.n)} & "
                f"{r.rho_morph_cov:+.2f} & {r.rho_meas_cov:+.2f} & "
                f"{r.r_raw:.2f} & {r.r_partial:.2f} & {r.retained_pct:.0f} " + r"\\")

agree_rows = []
NICE_T = {"StromalScore_ESTIMATE": "ESTIMATE stromal score",
          "stromal_deconv": "CPTAC deconvolution, stroma",
          "ImmuneScore_ESTIMATE": "ESTIMATE immune score",
          "immune_deconv": "CPTAC deconvolution, immune",
          "epithelial_cancer_deconv": "CPTAC deconvolution, epithelial"}
NICE_P = {"Stromal_fraction": "stromal fraction", "Inflammation_fraction": "inflammation",
          "Neoplastic_cellularity": "neoplastic cellularity"}
for _, r in ag.iterrows():
    agree_rows.append(f"{NICE_T[r.transcriptomic]} & {NICE_P[r.pathologist]} & "
                      f"{r.rho:+.2f} & {int(r.n)} " + r"\\")

caption = (
    r"Composition covariates not derived from the mRNA baseline (\S\ref{sec:stroma}). "
    r"The transcriptomic stromal, immune and proliferation scores used in "
    r"Fig.~\ref{fig:master2}b are computed from the same RNA that defines the baseline. "
    r"Here the morphology--residual link is instead tested against a pathologist's own "
    r"reading of the slide (PDAC: stromal fraction, neoplastic cellularity and inflammation "
    r"recorded for 77 of 137 cases) and against published deconvolution scores external to "
    r"this pipeline (ESTIMATE in CCRCC, UCEC and PDAC; xCell in GBM; the CPTAC cellular "
    r"deconvolution fractions in PDAC). Per cohort and score: $\rho_{\mathrm{morph}}$ and "
    r"$\rho_{\mathrm{meas}}$, Spearman correlation of the H\&E-only and of the measured "
    r"residual score with the first-listed covariate; $r$ raw, Pearson correlation of the "
    r"H\&E-only with the measured score; $r$ partial, the same after linear residualisation "
    r"of both on the covariate block; retained, partial as a percentage of raw. Only PDAC "
    r"has a pathologist-read composition, and there it removes more of the link (72\% "
    r"retained for translation) than any transcriptomic estimate does in the same patients "
    r"(97--100\%). The lower block shows why: transcriptomic composition estimates track "
    r"what the pathologist saw only weakly (Spearman, PDAC, $n$ cases with both). The LUAD "
    r"covariates are the ESTIMATE and purity estimates published with Gillette et al.\ 2020 "
    r"and redistributed by cBioPortal (study \texttt{luad\_cptac\_2020}), covering 104 of its "
    r"105 cases; LUAD is the cohort whose morphology score is least related to composition "
    r"($|\rho|\leq0.13$) and the only one that retains the link in full under every covariate "
    r"set. The retained fractions of "
    r"Fig.~\ref{fig:master2}b should therefore be read as upper bounds."
)

tex = [
    r"\begin{table}[htbp]\centering\small",
    r"\caption{" + caption + r"}",
    r"\label{tab:composition_path}",
    # nine columns overflowed the text width by ~77 pt: the retained column landed
    # off the page edge with no overfull warning, because the float is centred.
    # resizebox is the house style for the other wide supplementary tables.
    r"\resizebox{\textwidth}{!}{\begin{tabular}{lllrrrrrr}",
    r"\toprule",
    r"Cohort & Score & Covariate block & $n$ & $\rho_{\mathrm{morph}}$ & $\rho_{\mathrm{meas}}$ & "
    r"$r$ raw & $r$ partial & retained \% \\",
    r"\midrule",
]
tex += rows
tex += [
    r"\midrule",
    r"\multicolumn{9}{l}{\emph{PDAC: agreement of transcriptomic composition estimates with the pathologist}} \\",
    r"\multicolumn{3}{l}{transcriptomic estimate} & \multicolumn{3}{l}{pathologist read} & "
    r"\multicolumn{2}{r}{Spearman $\rho$} & $n$ \\",
]
for ar in agree_rows:
    t, p, rho, n = ar.replace(r" \\", "").split(" & ")
    tex.append(rf"\multicolumn{{3}}{{l}}{{{t}}} & \multicolumn{{3}}{{l}}{{{p}}} & "
               rf"\multicolumn{{2}}{{r}}{{{rho}}} & {n} \\")
tex += [r"\bottomrule", r"\end{tabular}}", r"\end{table}"]

with open("SuppTable_composition_path.tex", "w", encoding="utf8") as f:
    f.write("\n".join(tex) + "\n")
print(f"wrote SuppTable_composition_path.tex ({len(rows)} covariate rows, {len(agree_rows)} agreement rows)")
