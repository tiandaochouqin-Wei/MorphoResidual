#!/usr/bin/env python3
"""make_ablation_tables.py (LOCAL) — supplementary tables from server CSVs in figdata/:
  ablation_<cohort>.csv      (⑤ hyper-parameter robustness)
  abmil_<cohort>.csv         (⑥ aggregation ablation: ABMIL vs mean-pool)
  encoder_comparison.csv     (③ UNI / Phikon / H-optimus-0 / ResNet50-ImageNet)
Emits one LaTeX table per available file (\\hline, no booktabs) into ../SuppTable_*.tex ."""
import glob, os, pandas as pd

DD = "figdata"


def tex_table(df, caption, label, fmt, resize=False):
    cols = list(df.columns)
    L = [r"\begin{table}[htbp]\centering\small", rf"\caption{{{caption}}}", rf"\label{{{label}}}",
         (r"\resizebox{\textwidth}{!}{" if resize else "") + r"\begin{tabular}{" + "l" * len(cols) + r"}\hline",
         " & ".join((c if "$" in c else c.replace("_", r"\_").replace("R2", "$R^2$")) for c in cols) + r" \\ \hline"]
    for _, r in df.iterrows():
        L.append(" & ".join(fmt(c, r[c]) for c in cols) + r" \\")
    L += [r"\hline", r"\end{tabular}" + ("}" if resize else ""), r"\end{table}"]
    return "\n".join(L) + "\n"


def f(c, v):
    if isinstance(v, float):
        if v != v:  # NaN: the reference encoder has no self-agreement entry
            return "--"
        return f"{v:.3f}" if abs(v) < 10 else f"{v:.1f}"
    return str(v).replace("_", r"\_")


made = []
# ⑤ hyper-parameter ablation
abl = sorted(glob.glob(f"{DD}/ablation_*.csv"))
if abl:
    d = pd.concat([pd.read_csv(p) for p in abl])
    d = d[["cohort", "n_pcs", "ridge_lambda", "n_genes", "median_incr", "frac_pos", "corr_vs_default"]]
    open("../SuppTable_ablation.tex", "w", encoding="utf-8").write(tex_table(
        d, "Robustness of the per-gene morphology increment to the analysis hyper-parameters "
           "(number of morphology principal components; ridge $\\lambda$). For the FDR-significant "
           "proteins: median incremental $R^2$, fraction positive, and per-gene correlation with the "
           "default setting (20 PCs, $\\lambda=1$).", "tab:ablation", f))
    made.append("SuppTable_ablation.tex"); print(d.round(3).to_string(index=False))
# ⑥ ABMIL vs mean-pool
ab = sorted(glob.glob(f"{DD}/abmil_*.csv"))
if ab:
    d = pd.concat([pd.read_csv(p) for p in ab])[["cohort", "method", "n", "r", "R2"]]
    open("../SuppTable_abmil.tex", "w", encoding="utf-8").write(tex_table(
        d, "Aggregation ablation. Out-of-fold Pearson $r$ and $R^2$ for predicting the measured "
           "translation-residual score from tile embeddings with gated attention MIL (ABMIL) versus the "
           "main-analysis aggregation (mean-pool, 20 PCs, ridge), same patient-level folds.", "tab:abmil", f))
    made.append("SuppTable_abmil.tex"); print(d.round(3).to_string(index=False))
# ③ encoder comparison
if os.path.exists(f"{DD}/encoder_comparison.csv"):
    d = pd.read_csv(f"{DD}/encoder_comparison.csv")
    d["cohort"] = d["cohort"].str.upper()
    d["encoder"] = d["encoder"].replace({"ResNet50-ImageNet": "ResNet-50 (ImageNet)"})
    d = d.rename(columns={"n_sig": "significant proteins", "med_incr_sig": "median incr. $R^2$",
                          "jaccard_vs_uni": "Jaccard vs UNI", "recall_of_uni_sig": "recall of UNI set",
                          "r_incr_vs_uni": "per-gene $r$ vs UNI"})
    open("../SuppTable_encoders.tex", "w", encoding="utf-8").write(tex_table(
        d, "Encoder robustness. Per cohort and tile encoder: number of FDR-significant morphology-"
           "predictable residual proteins, median incremental $R^2$, and agreement with UNI (Jaccard of "
           "significant sets, recall of the UNI-significant set, per-gene incremental-$R^2$ correlation). "
           "ResNet-50 ImageNet is a non-pathology baseline.", "tab:encoders", f, resize=True))
    made.append("SuppTable_encoders.tex"); print(d.round(3).to_string(index=False))
print("wrote:", made if made else "nothing yet (no server CSVs in figdata/)")
