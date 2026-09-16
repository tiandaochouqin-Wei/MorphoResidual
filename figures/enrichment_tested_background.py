#!/usr/bin/env python3
"""enrichment_tested_background.py (LOCAL) -- P0-5 fix (adversarial review, 2026-09-06).

The originally reported Reactome over-representation used Enrichr's default
whole-genome background, but only 9,635-10,785 MS-detected proteins were
actually tested per cohort -- a universe already enriched for abundant
ribosomal/ER/secretory proteins. This script re-runs the over-representation
analysis (ORA) with the cohort's own tested-protein set as the background
(hypergeometric test via gseapy.enrich, no network dependency beyond a
one-time Reactome_2022.gmt download cached in figures/figdata/), and
separately quantifies RPS/RPL ribosomal-protein-specific enrichment (the
reviewer's most concrete falsifiable claim).

Requires: pip install gseapy (network access once, to cache Reactome_2022).
Writes:
  figdata/enrichment_tested_background.csv   (top-8 corrected terms x 5 cohorts)
  figdata/ribosomal_enrichment.csv           (RPS/RPL fold-enrichment x 5 cohorts)
  ../SuppTable_enrichment.tex                (formatted supplementary table)
"""
import os, re, json, warnings
import numpy as np, pandas as pd
from scipy import stats
import gseapy as gp

warnings.filterwarnings("ignore")
DD = "figdata"
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
NICE = {"ccrcc": "CCRCC", "luad": "LUAD", "ucec": "UCEC", "gbm": "GBM", "pdac": "PDAC"}
GMT = f"{DD}/Reactome_2022.gmt"
RIBO = re.compile(r"^RP[SL]\d+[A-Z]?$|^RPLP[012]$|^RPSA$")

if not os.path.exists(GMT):
    lib = gp.get_library(name="Reactome_2022", organism="Human")
    with open(GMT, "w", encoding="utf-8") as f:
        for term, genes in lib.items():
            f.write(term + "\t\t" + "\t".join(genes) + "\n")
    print(f"downloaded and cached Reactome_2022 ({len(lib)} terms) -> {GMT}")

table_rows, ribo_rows, full_rows = [], [], []
for c in COHORTS:
    d = pd.read_csv(f"{DD}/merged_incr_{c}.csv")
    tested = d["gene"].dropna().unique().tolist()
    sig = d[(d["fdr_uni"] < 0.05) & (d["incremental_r2_uni"] > 0)]["gene"].dropna().unique().tolist()
    enr = gp.enrich(gene_list=sig, gene_sets=GMT, background=tested, outdir=None)
    full_rows.append(enr.results.assign(cohort=NICE[c]))          # keep everything, for the dot matrix
    top = enr.results.sort_values("Adjusted P-value").head(3)
    for _, row in top.iterrows():
        table_rows.append(dict(cohort=NICE[c], term=row["Term"].split(" R-HSA")[0][:48],
                                overlap=row["Overlap"], q=row["Adjusted P-value"]))
    # ribosomal-protein-specific fold enrichment (Fisher exact, tested-proteome background)
    ribo_tested = [g for g in tested if RIBO.match(str(g))]
    ribo_sig = [g for g in sig if RIBO.match(str(g))]
    a = len(ribo_sig); b = len(ribo_tested) - a
    cc = len(sig) - a; dd = len(tested) - len(sig) - b
    odds, p = stats.fisher_exact([[a, b], [cc, dd]])
    frac_ribo = a / len(ribo_tested) if ribo_tested else np.nan
    frac_all = len(sig) / len(tested)
    ribo_rows.append(dict(cohort=NICE[c], n_ribo_tested=len(ribo_tested), n_ribo_sig=a,
                           frac_ribo_sig=round(frac_ribo, 3), frac_all_sig=round(frac_all, 3),
                           fold=round(frac_ribo / frac_all, 2), OR=round(odds, 2), p=p))

tab = pd.DataFrame(table_rows)
ribo = pd.DataFrame(ribo_rows)
tab.to_csv(f"{DD}/enrichment_tested_background.csv", index=False)
pd.concat(full_rows, ignore_index=True)[["cohort", "Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]]     .to_csv(f"{DD}/enrichment_tested_background_full.csv", index=False)
ribo.to_csv(f"{DD}/ribosomal_enrichment.csv", index=False)
print(tab.to_string(index=False))
print()
print(ribo.round(3).to_string(index=False))

L = [r"\begin{table}[htbp]\centering\small",
     r"\caption{Reactome over-representation of morphology-predictable residual proteins against "
     r"the cohort's own tested (MS-detected) protein set as background (hypergeometric test; "
     r"top three terms per cohort by Benjamini--Hochberg $q$), correcting the whole-genome-background "
     r"analysis reported in the main text (\S\ref{sec:repro}). Overlap, significant genes in term / "
     r"term size within the tested background.}",
     r"\label{tab:enrichment}",
     r"\resizebox{\textwidth}{!}{\begin{tabular}{llrr}\toprule",
     r"Cohort & Reactome term & Overlap & $q$ \\ \midrule"]
for _, r in tab.iterrows():
    L.append(f"{r['cohort']} & {r['term'].replace('_', r'{\\_}')} & {r['overlap']} & {r['q']:.2e} \\\\")
L += [r"\bottomrule", r"\end{tabular}}",
      r"\vspace{2pt}", r"{\footnotesize Ribosomal-protein (RPS/RPL) family fold-enrichment over the "
      r"cohort's overall significant fraction, same tested background (Fisher exact): " +
      "; ".join(f"{r['cohort']} {r['fold']:.1f}$\\times$ ($p={r['p']:.1e}$)" for _, r in ribo.iterrows()) +
      r".}", r"\end{table}"]
open("../SuppTable_enrichment.tex", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\nwrote SuppTable_enrichment.tex")
