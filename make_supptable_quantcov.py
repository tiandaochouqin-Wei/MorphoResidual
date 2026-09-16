#!/usr/bin/env python3
"""Build SuppTable_quantcov.tex -- the TMT quantification-completeness confound check (T1-6).

Written as a .py file rather than a shell heredoc because this project's Bash tool
eats backslashes, silently corrupting LaTeX macros.
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

DD = "figures/figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]

lg = pd.read_csv(f"{DD}/gene_covariate_logistic.csv")
en = pd.read_csv(f"{DD}/enrichment_complete_only.csv").set_index("cohort")
rb = pd.read_csv(f"{DD}/ribo_enrichment_matched.csv").set_index("cohort")

rows, summary = [], []
for c in COH:
    C = c.upper()
    d = pd.read_csv(f"{DD}/gene_covariates_{c}.csv")
    d = d[d.pinned_significant.notna()].dropna(subset=["incremental_r2"])
    sig = d.pinned_significant == True  # noqa: E712
    comp = d.pct_missing == 0
    a, b = int((sig & comp).sum()), int((~sig & comp).sum())
    e, f = int((sig & ~comp).sum()), int((~sig & ~comp).sum())
    rc, ri = 100 * a / (a + b), 100 * e / (e + f)
    _, pf = fisher_exact([[a, b], [e, f]])

    def g(term):
        r = lg[(lg.cohort == C) & (lg.term == term)].iloc[0]
        return r["or_"], r["ci_lo"], r["ci_hi"], r["p"]

    om, lo, hi, pm = g("pct_missing")
    ou, _, _, pu = g("ribo_unadjusted")
    oa, _, _, pa = g("ribo")
    ee = en.loc[C]

    rows.append(
        f"{C} & {100*comp.mean():.0f} & {rc:.1f} & {ri:.1f} & {rc/ri:.2f}$\\times$ & "
        f"{om:.2f} [{lo:.2f}, {hi:.2f}] & {ou:.2f} & {oa:.2f}{'' if pa<0.05 else '$^{\\dagger}$'} & "
        f"{rb.loc[C,'fold_reported']:.2f}$\\times$ & {rb.loc[C,'fold_matched']:.2f}$\\times$ & "
        f"{int(ee.top10_overlap)}/10 & {int(ee.n_sig_terms_shared)}/{int(ee.n_sig_terms_all)} \\\\"
    )
    summary.append((C, rc / ri, pf, om, pm, ou, pu, oa, pa,
                    rb.loc[C, "fold_matched"], rb.loc[C, "p_matched"],
                    int(ee.top10_overlap), int(ee.n_sig_terms_shared), int(ee.n_sig_terms_all),
                    ee.top_term_all, ee.top_term_complete))

caption = (
    r"TMT quantification completeness as a confounder of morphology-predictability "
    r"(\S\ref{sec:quantcov}). \% complete, tested proteins with no missing value across "
    r"the cohort's patients. Sig.\ rate, percentage of tested proteins that are "
    r"morphology-predictable (FDR${<}0.05$, positive increment) among fully quantified "
    r"versus incompletely quantified proteins, and their ratio (Fisher $p<10^{-20}$ in "
    r"every cohort). Missingness OR, adjusted odds ratio per 10 percentage points of "
    r"missing values from the logistic model "
    r"$\mathrm{logit}\,P(\text{significant}) \sim \text{\%\,missing} + \text{mean log-ratio} "
    r"+ R^2_{\mathrm{mRNA}} + \mathbf{1}[\text{ribosomal}]$ fitted over all tested genes; "
    r"values below 1 mean incompletely quantified proteins are less often significant. "
    r"Ribosomal OR, the cytosolic ribosomal-protein (RPS/RPL) family effect before and "
    r"after adjustment in that model. Ribosomal fold, the family's enrichment over the "
    r"cohort hit rate as reported in \S\ref{sec:enrich}, and after matching the background "
    r"on quantification completeness and abundance decile. Top-10 overlap and shared terms, "
    r"agreement of the Reactome over-representation analysis when both the significant set "
    r"and the background are restricted to fully quantified proteins, against the same "
    r"analysis on all tested proteins. $^{\dagger}$adjusted ribosomal effect not significant "
    r"at $p<0.05$."
)

tex = [
    r"\begin{table}[htbp]\centering\small",
    r"\caption{" + caption + r"}",
    r"\label{tab:quantcov}",
    r"\resizebox{\textwidth}{!}{\begin{tabular}{lrrrrlrrrrcc}\toprule",
    r"& \% & \multicolumn{3}{c}{Sig.\ rate (\%)} & Missingness OR & "
    r"\multicolumn{2}{c}{Ribosomal OR} & \multicolumn{2}{c}{Ribosomal fold} & "
    r"Top-10 & Shared \\",
    r"\cmidrule(lr){3-5}\cmidrule(lr){7-8}\cmidrule(lr){9-10}",
    r"Cohort & complete & complete & incompl. & ratio & (per 10\% missing) & unadj. & adj. & "
    r"reported & matched & overlap & terms \\ \midrule",
]
tex += rows
tex += [r"\bottomrule", r"\end{tabular}}", r"\end{table}"]

with open("SuppTable_quantcov.tex", "w", encoding="utf8") as fh:
    fh.write("\n".join(tex) + "\n")
print("wrote SuppTable_quantcov.tex\n")

print(f"{'coh':6s} {'ratio':>6s} {'missOR':>7s} {'riboU':>7s} {'riboA':>7s} {'p_adj':>9s} "
      f"{'matchFold':>9s} {'p_match':>9s} {'top10':>5s} {'shared':>9s}")
for (C, ratio, pf, om, pm, ou, pu, oa, pa, mf, pmt, t10, sh, na, ta, tc) in summary:
    print(f"{C:6s} {ratio:6.2f} {om:7.2f} {ou:7.2f} {oa:7.2f} {pa:9.1e} "
          f"{mf:9.2f} {pmt:9.1e} {t10:4d}/10 {sh:4d}/{na:<4d}")
    print(f"       all: {ta}\n       cmp: {tc}")
