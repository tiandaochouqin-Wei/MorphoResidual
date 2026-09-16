#!/usr/bin/env python3
"""Build SuppTable_txomeperm.tex from the permutation-FDR summary.
Written as a file (not a Bash heredoc) because this project's Bash tool eats
backslashes, which silently corrupts LaTeX macros such as \\times and \\\\."""
import pandas as pd

s = pd.read_csv("figures/figdata/transcriptome_baseline_perm_summary.csv")

rows = []
for _, r in s.iterrows():
    cens = "yes" if r["censored"] else "no"
    pct = 100 * r["n_p05"] / r["m"]
    rows.append(
        f"{r['cohort']} & {int(r['m'])} & {int(r['n_p05'])} ({pct:.1f}\\%) & "
        f"{int(r['exp_p05'])} & {r['fold']:.1f}$\\times$ & {int(r['n_at_floor'])} & "
        f"{int(r['need_at_floor'])} & {int(r['n_fdr_sig_pos'])} & {cens} \\\\"
    )

caption = (
    r"Permutation null and FDR under the transcriptome-wide baseline "
    r"(\S\ref{sec:txomemeth}). For each cohort's pinned significant genes, the morphology "
    r"block was re-tested over the own-mRNA${+}$20 RNA-component baseline under a "
    r"patient$\leftrightarrow$slide permutation null ($N=1{,}000$) with Benjamini--Hochberg "
    r"FDR. $p{<}0.05$, genes clearing an uncorrected permutation $p$ (the "
    r"resolution-independent quantity); expected, the 5\% null expectation; fold, their ratio "
    r"(binomial $p<10^{-119}$ in every cohort). At floor, genes at the smallest attainable "
    r"$p=1/1001$; needed, the number that must sit at that floor before Benjamini--Hochberg "
    r"can certify any gene at $\alpha=0.05$ in a set of this size. FDR-sig., genes passing "
    r"FDR${<}0.05$ with a positive increment. Censored, whether the cohort has fewer floor "
    r"genes than the procedure requires, in which case a count of zero reflects the "
    r"permutation budget rather than absence of signal."
)

tex = [
    r"\begin{table}[htbp]\centering\small",
    r"\caption{" + caption + r"}",
    r"\label{tab:txomeperm}",
    r"\begin{tabular}{lrrrrrrrc}\toprule",
    r"Cohort & Tested & $p{<}0.05$ & Expected & Fold & At floor & Needed & FDR-sig. & Censored \\ \midrule",
]
tex += rows
tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

with open("SuppTable_txomeperm.tex", "w", encoding="utf8") as f:
    f.write("\n".join(tex) + "\n")
print("wrote SuppTable_txomeperm.tex")
