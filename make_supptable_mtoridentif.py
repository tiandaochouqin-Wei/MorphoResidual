#!/usr/bin/env python3
"""SuppTable_mtoridentif.tex -- the protein-normalised phospho-S6 test and why it is
not identifiable in this design (T1-4). Written with the Write tool, not a heredoc."""
import pandas as pd

d = pd.read_csv("figures/figdata/mtor_identifiability.csv")
ORDER = ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]
d = d.set_index("cohort").loc[ORDER].reset_index()

rows = []
for _, r in d.iterrows():
    forced = 100 * (r.rho_raw - r.placebo_denom_mean) / (r.rho_raw - r.rho_norm)
    lvl = "gene" if r.gene_level else "site"
    rows.append(
        f"{r.cohort} & {lvl} & {int(r.n)} & {int(r.k_arm)} & "
        f"{r.rho_outcome_vs_ribo_panel:+.2f} & {r.rho_raw:+.2f} & {r.rho_norm:+.2f} & "
        f"{r.placebo_denom_mean:+.2f} & {r.pct_in_placebo_denom:.0f} & {forced:.0f} & "
        f"{r.perm_null_median:+.2f} & {int(r.bg_n)} & {r.bg_median_norm:+.2f} & "
        f"{r.pct_bg_size_matched:.0f} " + r"\\"
    )

caption = (
    r"The protein-normalised phospho-S6 test is not identifiable in this design "
    r"(\S\ref{sec:mtormeth}). For the S6K1 arm (phospho-RPS6 S235/236/240/244; gene-level "
    r"RPS6 in CCRCC), Spearman $\rho$ with the measured translation residual before (raw) "
    r"and after normalising each site to its host protein (norm.). $\rho_{\mathrm{ribo}}$, "
    r"correlation of the mean abundance of the cytosolic ribosomal proteins with the "
    r"translation residual, i.e.\ how nearly the denominator is a copy of the outcome. "
    r"Placebo, the mean $\rho$ obtained when each of the other ribosomal proteins is "
    r"substituted for RPS6 as the denominator, and the percentile at which the true "
    r"denominator falls in that distribution; forced, the fraction of the raw-to-normalised "
    r"attenuation that the placebo denominator reproduces without any host-specific "
    r"information. Perm.\ null, median normalised $\rho$ after permuting the phosphosite "
    r"across patients ($N{=}1{,}000$), the correct null for the normalised statistic. "
    r"Background, non-mTORC1 ribosomal phosphosites each normalised to their own host; "
    r"percentile is against a size-matched null (means of $k$ background sites, $2{,}000$ "
    r"draws). Three things are seen in every cohort: the denominator is collinear with the "
    r"outcome, most of the attenuation is reproduced by an unrelated denominator, and the "
    r"arm's rank in an identically normalised background is not reproducible across organs. "
    r"Neither the attenuation nor the sign of the normalised statistic can therefore be read "
    r"as evidence about S6K1 activity, and the manuscript's conclusion that these data cannot "
    r"separate mTORC1 signalling from ribosomal-protein abundance stands for a structural "
    r"reason rather than for want of normalisation."
)

tex = [
    r"\begin{table}[htbp]\centering\small",
    r"\caption{" + caption + r"}",
    r"\label{tab:mtoridentif}",
    r"\resizebox{\textwidth}{!}{\begin{tabular}{llrrrrrrrrrrrr}",
    r"\toprule",
    r"& & & & & \multicolumn{2}{c}{S6K1 arm $\rho$} & \multicolumn{3}{c}{Placebo denominator} & "
    r"& \multicolumn{3}{c}{Normalised background} \\",
    r"\cmidrule(lr){6-7}\cmidrule(lr){8-10}\cmidrule(lr){12-14}",
    r"Cohort & level & $n$ & $k$ & $\rho_{\mathrm{ribo}}$ & raw & norm.\ & mean & own pctile & "
    r"forced \% & perm.\ null & $n$ & median & pctile \\",
    r"\midrule",
]
tex += rows
tex += [r"\bottomrule", r"\end{tabular}}", r"\end{table}"]

with open("SuppTable_mtoridentif.tex", "w", encoding="utf8") as f:
    f.write("\n".join(tex) + "\n")
print("wrote SuppTable_mtoridentif.tex")
