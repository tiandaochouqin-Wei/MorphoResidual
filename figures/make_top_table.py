#!/usr/bin/env python3
"""Supplementary Table — top morphology-predictable residual proteins per cohort
(gene, n, incremental R2, FDR, functional family). Writes CSV + LaTeX (\\hline, no booktabs)."""
import re, numpy as np, pandas as pd
import mrstyle as S

BASE = "D:/claude/pajinsen/MorphoResidual/pilot_ccrcc/results_server"
PATHS = {"CCRCC": f"{BASE}/residual_results_tumoronly.csv",
         "LUAD": f"{BASE}/luad/results/residual_results_tumoronly.csv",
         "UCEC": f"{BASE}/ucec/results/residual_results_tumoronly.csv",
         "GBM": f"{BASE}/gbm/results/residual_results_tumoronly.csv",
         "PDAC": f"{BASE}/pdac/results/residual_results_tumoronly.csv"}
TOPN = 12
SECRETION_ER = {"SSR1", "SSR2", "SSR3", "SSR4", "SEC11A", "SEC11C", "SEC61A1", "SEC61B", "SEC61G",
                "SRP9", "SRP14", "SRP54", "SRP68", "SRP72", "SRPRA", "SRPRB", "SPCS1", "SPCS2", "SPCS3",
                "RPN1", "RPN2", "DDOST", "STT3A", "STT3B", "OST4", "PIGS", "TANGO2", "SEC13", "SEC23A",
                "SEC24C", "SEC31A", "LMAN1", "SURF4", "TMED2", "TMED9"}
MATRISOME = {"LUM", "ASPN", "DCN", "POSTN", "COL1A1", "COL1A2", "COL3A1", "HAPLN3", "FN1", "SPARC",
             "THBS2", "FBN1", "COL5A1", "COL6A3", "OGN", "PRELP", "GPX3"}


def family(g):
    if re.match(r"^(RPS|RPL|EEF|EIF)\d", g) or g in {"RACK1", "FAU"}:
        return "translation"
    if g in SECRETION_ER:
        return "secretion-ER"
    if g in MATRISOME:
        return "matrisome"
    if re.match(r"^(SF3|SNRP|SRSF|HNRNP|PRPF|U2AF|RBM|DDX|DHX|LSM|SNRNP|SNU|CWC)", g):
        return "splicing"
    return "other"


rows = []
for c in S.COH:
    d = pd.read_csv(PATHS[c])
    s = d[(d.fdr < 0.05) & (d.incremental_r2 > 0)].sort_values("incremental_r2", ascending=False).head(TOPN)
    for _, r in s.iterrows():
        rows.append({"cohort": c, "gene": r.gene, "n": int(r.n), "incr_r2": float(r.incremental_r2),
                     "fdr": float(r.fdr), "family": family(str(r.gene))})
t = pd.DataFrame(rows)
t.to_csv("../SuppTable_top_proteins.csv", index=False)

# longtable (non-floating, page-spanning) so a 60-row table cannot orphan the section heading
L = [r"\begin{small}",
     r"\begin{longtable}{llrrrl}",
     r"\caption{Top morphology-predictable residual proteins per cohort. For each cohort the "
     + str(TOPN) + r" FDR-significant proteins with the largest incremental $R^2$ of morphology over "
     r"the mRNA-only baseline; $n$, patients with a measurement; family, functional annotation.}"
     r"\label{tab:topproteins}\\",
     r"\toprule",
     r"Cohort & Protein & $n$ & incr.\ $R^2$ & FDR & family \\ \midrule",
     r"\endfirsthead",
     r"\multicolumn{6}{l}{\textit{(continued)}}\\ \toprule",
     r"Cohort & Protein & $n$ & incr.\ $R^2$ & FDR & family \\ \midrule",
     r"\endhead",
     r"\endfoot",
     r"\bottomrule",
     r"\endlastfoot"]
for i, c in enumerate(S.COH):
    sub = t[t.cohort == c]
    for j, (_, r) in enumerate(sub.iterrows()):
        L.append(f"{c if j == 0 else ''} & {r.gene} & {r.n} & {r.incr_r2:.2f} & {r.fdr:.1e} & {r.family} \\\\")
    if i < len(S.COH) - 1:
        L.append(r"\midrule")
L += [r"\end{longtable}", r"\end{small}"]
open("../SuppTable_top_proteins.tex", "w", encoding="utf-8").write("\n".join(L) + "\n")
print(f"wrote SuppTable_top_proteins.csv/.tex  ({len(t)} rows)")
print(t.groupby("cohort")["family"].value_counts().unstack(fill_value=0))
