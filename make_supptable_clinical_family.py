#!/usr/bin/env python3
"""Rebuild SuppTable_clinical_family.tex over the complete five-cohort family (T1-7).

108 tests. q printed to two significant figures throughout rather than censored at
"<0.001", so that every value the main text quotes can be checked against this
table -- it is the source of record for the multiplicity accounting and a referee
who cannot reproduce a quoted q from it has a legitimate complaint.

longtable, not table: 108 rows will not fit on a page.
Written as a .py file because this project's Bash tool eats backslashes.
"""
import pandas as pd

d = pd.read_csv("figures/figdata/clinical_family_all.csv")
for col in ["family", "mode", "outcome"]:
    d[col] = d[col].str.replace("&", r"\&", regex=False)

ORDER = {"CCRCC": 0, "LUAD": 1, "UCEC": 2, "GBM": 3, "PDAC": 4}


def fmt_stat(r):
    if r["outcome"].startswith("Cox"):
        return f"HR {r['stat']:.2f}"
    if "log-rank" in r["outcome"]:
        return f"{r['stat']:.2f}"
    return f"{r['stat']:+.3f}"


def q(v):
    return f"{v:.2e}" if v < 0.01 else f"{v:.3f}"


MARK = {"": "", "missing-adjusted": r"$^{\ddagger}$", "age": r"$^{\S}$"}

rows, prev = [], None
for _, r in d.sort_values(["cohort", "family", "mode", "outcome"],
                          key=lambda s: s.map(ORDER) if s.name == "cohort" else s).iterrows():
    if prev is not None and r["cohort"] != prev:
        rows.append(r"\midrule")
    prev = r["cohort"]
    rows.append(f"{r['cohort']} & {r['family']} & {r['mode']} & "
                f"{r['outcome']}{MARK[r['added'] if isinstance(r['added'], str) else '']} & "
                f"{fmt_stat(r)} & {r['p']:.2e} & {q(r['q_cohort'])} & {q(r['q_all'])} \\\\")

n = len(d)
per = d.groupby("cohort").size().reindex(list(ORDER))
n_age = int((d.added == "age").sum())
caption = (
    r"Complete five-cohort clinical test family (\S\ref{sec:clinical}). Every clinical "
    r"test run in this study is enumerated: each residual pathway score (translation, "
    r"ER-secretion and, where a matrisome score was formed, matrisome), as the measured "
    r"proteomic score and as the H\&E-only morphology-predicted score, against tumour "
    r"grade and pathologic stage (Spearman) and overall survival (median-split log-rank), "
    r"plus the Cox models---univariable, grade${+}$stage-adjusted, age-adjusted and "
    r"grade${+}$stage${+}$age-adjusted---giving "
    f"{n}"
    r" tests in total ("
    + ", ".join(f"{c} {int(v)}" for c, v in per.items()) + r"). "
    r"$^{\ddagger}$grade${+}$stage-adjusted Cox models for the CCRCC ER-secretion and "
    r"matrisome scores, which the original analysis omitted and which the enumeration rule "
    r"requires; refitted here. "
    r"$^{\S}$the "
    f"{n_age}"
    r" age-adjusted Cox models, possible because \texttt{days\_to\_birth} is populated "
    r"even though \texttt{age\_at\_index} is not (\S\ref{sec:clinical}). They are counted "
    r"in the family rather than exempted as sensitivity analyses, which is the more "
    r"conservative choice: excluding them would move exactly three borderline tests from "
    r"just above to just below $q{=}0.05$ (LUAD ER-secretion vs.\ grade, CCRCC measured "
    r"translation vs.\ survival, LUAD measured translation vs.\ stage) and would change no "
    r"claim made in the text. "
    r"GBM contributes survival tests only: the CPTAC GBM clinical file records neither "
    r"grade (the cohort is uniformly WHO grade~IV) nor pathologic stage, so no "
    r"grade${+}$stage adjustment is possible there either. "
    r"$q_{\mathrm{cohort}}$, Benjamini--Hochberg $q$ within that cohort's own family; this "
    r"is the primary correction, and the one the main text, Table~\ref{tab:clinical}, "
    r"Fig.~\ref{fig:master2}d and Supplementary Fig.~\ref{fig:supp2}f all use. "
    r"$q_{\mathrm{all}}$, over all "
    f"{n}"
    r" tests. Pooling is not uniformly more conservative: because it changes ranks, "
    r"$q_{\mathrm{all}}$ is smaller than $q_{\mathrm{cohort}}$ in "
    f"{int((d.q_all < d.q_cohort).sum())}"
    r" of the "
    f"{n}"
    r" rows. The two levels reach different conclusions at $\alpha{=}0.05$ for only "
    f"{int((((d.q_cohort < 0.05) & (d.q_all >= 0.05)) | ((d.q_cohort >= 0.05) & (d.q_all < 0.05))).sum())}"
    r" tests, both of them the CCRCC univariable and age-adjusted Cox models, which the "
    r"text already reports as not independent of stage."
)

tex = [
    r"\begingroup\small",
    r"\begin{longtable}{lllllrrr}",
    r"\caption{" + caption + r"}\label{tab:clinical_family}\\",
    r"\toprule",
    r"Cohort & Score & Input & Outcome & Statistic & $p$ & $q_{\mathrm{cohort}}$ & $q_{\mathrm{all}}$ \\ \midrule",
    r"\endfirsthead",
    r"\multicolumn{8}{l}{\small\itshape Table~\ref{tab:clinical_family} continued}\\",
    r"\toprule",
    r"Cohort & Score & Input & Outcome & Statistic & $p$ & $q_{\mathrm{cohort}}$ & $q_{\mathrm{all}}$ \\ \midrule",
    r"\endhead",
    r"\bottomrule",
    r"\endlastfoot",
]
tex += rows
tex += [r"\end{longtable}", r"\endgroup"]

with open("SuppTable_clinical_family.tex", "w", encoding="utf8") as f:
    f.write("\n".join(tex) + "\n")
print(f"wrote SuppTable_clinical_family.tex: {n} tests "
      f"({int((d.added=='').sum())} original, {int((d.added=='missing-adjusted').sum())} restored, "
      f"{n_age} age), {len(rows)} body lines")
