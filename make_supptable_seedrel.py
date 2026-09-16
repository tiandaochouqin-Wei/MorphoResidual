#!/usr/bin/env python3
"""SuppTable_seedrel.tex -- seed-reliability of the per-gene incremental R^2 (T1-7).

Written with the Write tool, never via a shell heredoc: this project's Bash tool
eats backslashes, which silently turns a LaTeX row terminator into a single
backslash and produces a "Misplaced \\noalign" a hundred lines later.
"""
import pandas as pd

S = pd.read_csv("figures/figdata/seed_reliability_summary.csv")
ENC = {"CCRCC": 0.808, "LUAD": 0.780, "UCEC": 0.845, "GBM": 0.669, "PDAC": 0.854}
ORDER = ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]

rows = []
for c in ORDER:
    a = S[(S.cohort == c) & (S.gene_set == "all tested")].iloc[0]
    g = S[(S.cohort == c) & (S.gene_set == "significant")].iloc[0]
    rows.append(
        f"{c} & {int(a.n_genes)} & {a.rel_single_seed:.3f} & {a.icc21:.3f} & "
        f"{a.rel_10seed_mean:.3f} & {int(g.n_genes)} & {g.rel_single_seed:.3f} & "
        f"{100 * g.frac_incr_pos:.0f} & {ENC[c]:.3f} & "
        f"{min(ENC[c] / a.rel_single_seed, 1.0) * 100:.0f} " + r"\\"
    )

caption = (
    r"Reliability of the per-gene incremental $R^2$ under ten cross-validation seeds "
    r"(\S\ref{sec:seedrel}). Only the fold assignment varies; patients, the cohort-level "
    r"morphology PCA basis and the ridge are held fixed, and the seed-0 estimates reproduce "
    r"the pinned per-gene values at Pearson $r=1.00000$ (maximum absolute difference "
    r"$5\times10^{-7}$). $\rho_1$, reliability of a single-seed estimate, taken as the mean "
    r"correlation over the 45 seed pairs; ICC, the more conservative intraclass correlation "
    r"ICC$(2,1)$, which also penalises systematic between-seed shifts; $\rho_{10}$, "
    r"Spearman--Brown reliability of the ten-seed mean. \% incr.${>}0$, the percentage of "
    r"each cohort's significant genes that retain a positive point estimate when the folds "
    r"are redrawn---a per-gene stability statement, not an error rate. $r$ vs.\ UNI, the "
    r"per-gene incremental-$R^2$ correlation between UNI and Phikon "
    r"(SuppTable~\ref{tab:encoders}); \% of ceiling, that correlation as a percentage of "
    r"$\rho_1$, the highest value attainable between two equally noisy estimates of the same "
    r"quantity. Two encoders therefore agree about as well as one encoder agrees with itself "
    r"under a different fold split. Because seed resampling isolates only the "
    r"fold-assignment variance component, $\rho_1$ is an upper bound on true test--retest "
    r"reliability and these percentages are lower bounds."
)

tex = [
    r"\begin{table}[htbp]\centering\small",
    r"\caption{" + caption + r"}",
    r"\label{tab:seedrel}",
    r"\begin{tabular}{lrrrrrrrrr}",
    r"\toprule",
    r"& \multicolumn{4}{c}{All tested genes} & \multicolumn{3}{c}{Significant set} & "
    r"\multicolumn{2}{c}{Encoder agreement} \\",
    r"\cmidrule(lr){2-5}\cmidrule(lr){6-8}\cmidrule(lr){9-10}",
    r"Cohort & genes & $\rho_1$ & ICC & $\rho_{10}$ & genes & $\rho_1$ & \% incr.${>}0$ & "
    r"$r$ vs.\ UNI & \% of ceiling \\",
    r"\midrule",
]
tex += rows
tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

with open("SuppTable_seedrel.tex", "w", encoding="utf8") as f:
    f.write("\n".join(tex) + "\n")
print("wrote SuppTable_seedrel.tex")
