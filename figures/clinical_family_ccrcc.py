#!/usr/bin/env python3
"""clinical_family_ccrcc.py (LOCAL) -- P1-12 fix (adversarial review, 2026-09-06).

Full clinical test family for CCRCC (the only cohort with a local clinical file):
3 scores (translation, ER-secretion, matrisome) x 2 versions (measured,
morphology-only H&E) x 3 outcomes (grade Spearman, stage Spearman, survival
median-split log-rank) = 18 association tests, plus 3 morphology-only Cox models
(univariable) = 21, plus the grade+stage-adjusted Cox for translation_morph = the
24-test family the main text's Table~tab:clinical draws 4 cells from. Also runs
the full CCRCC survival analysis (KM, Cox uni/adjusted, Schoenfeld) with lifelines.

Requires: pip install lifelines. Writes ../SuppTable_clinical_family.tex and
figdata/clinical_family_ccrcc.csv .
"""
import numpy as np, pandas as pd
from scipy import stats as st
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

DD = "figdata"
sc = pd.read_csv(f"{DD}/scores_ccrcc.csv", index_col=0)
cl = pd.read_csv(f"{DD}/clinical_link_ccrcc_clinical.csv", index_col=0)
d = sc.join(cl, how="inner")
print(f"n={len(d)} (grade/stage); survival n={int((d['time'] >= 0).sum())} after excluding negative follow-up")

SCORES = ["translation", "secretion_ER", "matrisome"]
rows = []
for score in SCORES:
    for ver in ["meas", "morph"]:
        col = f"{score}_{ver}"
        for var in ["grade", "stage"]:
            m = d[col].notna() & d[var].notna()
            rho, p = st.spearmanr(d.loc[m, col], d.loc[m, var])
            rows.append(dict(score=col, var=var, stat="spearman_rho", value=rho, p=p, n=int(m.sum())))
        dv = d[d["time"] >= 0].copy()
        med = dv[col].median()
        hi, lo = dv[dv[col] >= med], dv[dv[col] < med]
        lr = logrank_test(hi["time"], lo["time"], hi["event"], lo["event"])
        rows.append(dict(score=col, var="survival_logrank", stat="chi2", value=lr.test_statistic,
                          p=lr.p_value, n=len(dv)))

# Cox models for the three morphology-only scores (univariable) + adjusted translation_morph
dv = d[d["time"] >= 0].copy()
for zc in ["translation_morph", "secretion_ER_morph", "matrisome_morph", "grade", "stage"]:
    dv[zc + "_z"] = (dv[zc] - dv[zc].mean()) / dv[zc].std(ddof=0)
for score in ["translation_morph", "secretion_ER_morph", "matrisome_morph"]:
    cph = CoxPHFitter(penalizer=1e-3)
    cph.fit(dv[[score + "_z", "time", "event"]].rename(columns={score + "_z": "score"}),
            duration_col="time", event_col="event")
    s = cph.summary.loc["score"]
    rows.append(dict(score=score, var="cox_uni", stat="HR", value=np.exp(s["coef"]), p=s["p"], n=len(dv)))

cph_adj = CoxPHFitter(penalizer=1e-3)
cph_adj.fit(dv[["translation_morph_z", "grade_z", "stage_z", "time", "event"]],
            duration_col="time", event_col="event")
sa = cph_adj.summary
rows.append(dict(score="translation_morph", var="cox_adj", stat="HR", value=np.exp(sa.loc["translation_morph_z", "coef"]),
                  p=sa.loc["translation_morph_z", "p"], n=len(dv)))
rows.append(dict(score="stage(adj)", var="cox_adj", stat="HR", value=np.exp(sa.loc["stage_z", "coef"]),
                  p=sa.loc["stage_z", "p"], n=len(dv)))
rows.append(dict(score="grade(adj)", var="cox_adj", stat="HR", value=np.exp(sa.loc["grade_z", "coef"]),
                  p=sa.loc["grade_z", "p"], n=len(dv)))
cph_uni_t = CoxPHFitter(penalizer=1e-3); cph_uni_t.fit(dv[["translation_morph_z", "time", "event"]], "time", "event")
s_uni = cph_uni_t.summary.loc["translation_morph_z"]
rows.append(dict(score="translation_morph", var="cox_uni_ci_lo", stat="HR", value=np.exp(s_uni["coef lower 95%"]), p=np.nan, n=len(dv)))
rows.append(dict(score="translation_morph", var="cox_uni_ci_hi", stat="HR", value=np.exp(s_uni["coef upper 95%"]), p=np.nan, n=len(dv)))

full = pd.DataFrame(rows)
# BH correction: (a) the 18 association tests only, (b) all 24 (assoc + Cox)
assoc = full[full["var"].isin(["grade", "stage", "survival_logrank"])].copy()
assoc["q_ccrcc18"] = multipletests(assoc["p"], method="fdr_bh")[1]
cox_only = full[full["var"].isin(["cox_uni", "cox_adj"]) & full["score"].str.endswith("_morph")].copy()
allp = pd.concat([assoc[["p"]], cox_only[["p"]]])
q_all24 = multipletests(allp["p"], method="fdr_bh")[1]
assoc["q_ccrcc24"] = q_all24[:len(assoc)]
cox_only["q_ccrcc24"] = q_all24[len(assoc):]
out = pd.concat([assoc, cox_only], ignore_index=True)
out.to_csv(f"{DD}/clinical_family_ccrcc.csv", index=False)
pd.set_option("display.width", 200); print(out.round(4).to_string(index=False))

# ---- survival headline numbers (printed for the manuscript text; not re-typed by hand) ----
km = KaplanMeierFitter(); km.fit(dv["time"], dv["event"])
med_fu = KaplanMeierFitter().fit(dv["time"], 1 - dv["event"]).median_survival_time_
cph0 = CoxPHFitter(penalizer=1e-3); cph0.fit(dv[["translation_morph_z", "time", "event"]], "time", "event")
sph = cph0.check_assumptions(dv[["translation_morph_z", "time", "event"]], show_plots=False, p_value_threshold=1.0)
print(f"\nn_survival={len(dv)} events={int(dv['event'].sum())} median_followup_days(reverseKM)={med_fu:.0f}")
print(f"Cox uni HR={np.exp(cph0.summary.loc['translation_morph_z','coef']):.3f} "
      f"CI=[{np.exp(cph0.summary.loc['translation_morph_z','coef lower 95%']):.3f},"
      f"{np.exp(cph0.summary.loc['translation_morph_z','coef upper 95%']):.3f}] "
      f"p={cph0.summary.loc['translation_morph_z','p']:.4f}")
print(f"Cox adj HR={np.exp(sa.loc['translation_morph_z','coef']):.3f} "
      f"CI=[{np.exp(sa.loc['translation_morph_z','coef lower 95%']):.3f},"
      f"{np.exp(sa.loc['translation_morph_z','coef upper 95%']):.3f}] "
      f"p={sa.loc['translation_morph_z','p']:.4f}  "
      f"| stage HR={np.exp(sa.loc['stage_z','coef']):.3f} p={sa.loc['stage_z','p']:.4f}  "
      f"grade HR={np.exp(sa.loc['grade_z','coef']):.3f} p={sa.loc['grade_z','p']:.4f}")

# ---- SuppTable ----
BS = chr(92); EOL = BS + BS
NICE_SCORE = {"translation_meas": "Transl.\\ (meas.)", "translation_morph": "Transl.\\ (H\\&E-only)",
              "secretion_ER_meas": "ER-secr.\\ (meas.)", "secretion_ER_morph": "ER-secr.\\ (H\\&E-only)",
              "matrisome_meas": "Matrisome (meas.)", "matrisome_morph": "Matrisome (H\\&E-only)"}
NICE_VAR = {"grade": "grade", "stage": "stage", "survival_logrank": "survival (log-rank)"}
L = [BS + "begin{table}[htbp]" + BS + "centering" + BS + "small",
     BS + "caption{Full CCRCC clinical test family (P1-12, adversarial review): all 18 score-by-outcome "
          "association tests (3 residual-pathway scores $" + BS + "times$ measured or H" + BS + "&E-only "
          "$" + BS + "times$ grade/stage Spearman correlation or median-split log-rank survival test), plus "
          "the 3 univariable and 1 grade+stage-adjusted Cox models on the morphology-only scores reported in "
          "the text (24 tests total). $q_{18}$, Benjamini--Hochberg $q$ over the 18 association tests; "
          "$q_{24}$, over all 24. This is the only cohort with a locally available clinical file; the "
          "equivalent family for LUAD, UCEC, GBM and PDAC requires the server-side "
          "results/clinical\\_link\\_$" + BS + "langle$cohort$" + BS + "rangle$\\_\\{assoc,cox\\}.csv files.}",
     BS + "label{tab:clinical_family}",
     BS + "resizebox{" + BS + "textwidth}{!}{" + BS + "begin{tabular}{llrrrr}" + BS + "toprule",
     "Score & Outcome & statistic & $p$ & $q_{18}$ & $q_{24}$ " + EOL + " " + BS + "midrule"]
for _, r in assoc.sort_values(["score", "var"]).iterrows():
    stat_s = f"{r.value:+.3f}" if r["stat"] == "spearman_rho" else f"{r.value:.2f}"
    L.append(f"{NICE_SCORE.get(r.score, r.score)} & {NICE_VAR.get(r['var'], r['var'])} & {stat_s} & {r.p:.2e} & {r.q_ccrcc18:.3f} & {r.q_ccrcc24:.3f} " + EOL)
L.append(BS + "midrule")
for _, r in cox_only.iterrows():
    L.append(f"{NICE_SCORE.get(r.score, r.score)} & {r['var'].replace('_', ' ')} & HR {r.value:.2f} & {r.p:.2e} & --- & {r.q_ccrcc24:.3f} " + EOL)
L += [BS + "bottomrule", BS + "end{tabular}}", BS + "end{table}"]
open("../SuppTable_clinical_family.tex", "w", encoding="utf-8").write(chr(10).join(L) + chr(10))
print("\nwrote SuppTable_clinical_family.tex")
