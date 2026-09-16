#!/usr/bin/env python3
r"""T1-5: does the morphology->residual link survive composition covariates that
are NOT derived from the mRNA baseline?

WHY A SECOND COMPOSITION CONTROL
The manuscript's composition control (sec:stroma) adds transcriptomic stromal,
immune and proliferation scores to the mRNA baseline and re-measures the
morphology increment. A referee's fair objection: those scores are computed from
the same RNA that defines the baseline, so they are the baseline's own view of
composition, not an independent measurement of it. Two independent views exist:

  PATHOLOGIST-READ  In PDAC, a pathologist recorded stromal fraction, neoplastic
                    cellularity, inflammation and other tissue fractions for 77
                    of the 137 analysed cases (Cao et al. 2021 Table S1). This is
                    the gold standard: what a human sees on the slide.
  DECONVOLUTION     ESTIMATE (CCRCC, UCEC, PDAC), xCell (GBM) and the CPTAC
                    cellular-deconvolution fractions (PDAC) are transcriptomic
                    but published, external to this pipeline, and not the
                    24-gene proxy the manuscript used.

LUAD's covariates are the ESTIMATE stromal/immune scores and the two purity
estimates published with Gillette et al. 2020 (Cell 182:200) and redistributed by
cBioPortal as study luad_cptac_2020, which is the same Table S1 content reached
without the publisher paywall. 104 of the 105 analysed LUAD cases are covered.

WHAT IS TESTED (per cohort, per score: translation and ER-secretion)
  rho(H&E-only score, covariate)      does the morphology score read composition?
  rho(measured score, covariate)      does the residual itself track composition?
  r(morph, measured) raw              the link the paper reports
  r(morph, measured | covariates)     the link after removing composition
  retained %                          partial / raw

and, for PDAC where both exist,
  rho(transcriptomic estimate, pathologist read)   how well a deconvolution score
                                                   tracks what the pathologist saw,
                                                   which bounds how much any
                                                   transcriptomic adjustment can
                                                   have missed.

Writes figures/figdata/composition_pathologist.csv
       figures/figdata/composition_proxy_agreement.csv
"""
import numpy as np
import pandas as pd
from scipy import stats as st

DD = "figures/figdata"
SCORES = ["translation", "secretion_ER"]


def norm_id(s):
    s = str(s).strip().replace(".", "-")
    return "-".join(s.split("-")[:2]) if s.startswith("C3") else s


def load_scores(c):
    d = pd.read_csv(f"{DD}/scores_{c}.csv").rename(columns={"Unnamed: 0": "case"})
    d["case"] = d["case"].map(norm_id)
    return d.set_index("case")


def partial_r(x, y, C):
    """Pearson r between x and y after linear residualisation on covariate block C."""
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(C), axis=1)
    x, y, C = x[m], y[m], C[m]
    if m.sum() < 20:
        return np.nan, int(m.sum())
    Z = np.column_stack([np.ones(len(x)), C])
    rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    return st.pearsonr(rx, ry)[0], int(m.sum())


def spearman(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return (st.spearmanr(a[m], b[m])[0], int(m.sum())) if m.sum() >= 20 else (np.nan, int(m.sum()))


# ----------------------------------------------------------- covariates per cohort
def cov_pdac():
    cl = pd.read_excel(f"{DD}/subtype_pdac_raw.xlsx", sheet_name="Clinical_data")
    cl["case"] = cl["case_id"].map(norm_id)
    path = cl.set_index("case")[["Stromal_fraction", "Neoplastic_cellularity",
                                 "Inflammation_fraction"]].apply(pd.to_numeric, errors="coerce")
    mp = pd.read_excel(f"{DD}/subtype_pdac_raw.xlsx", sheet_name="Molecular_phenotype_data")
    mp["case"] = mp["case_id"].map(norm_id)
    dec = mp.set_index("case")[["stromal_deconv", "immune_deconv", "epithelial_cancer_deconv",
                                "StromalScore_ESTIMATE", "ImmuneScore_ESTIMATE"]] \
            .apply(pd.to_numeric, errors="coerce")
    return {
        "pathologist (stroma+cellularity+inflammation)":
            path[["Stromal_fraction", "Neoplastic_cellularity", "Inflammation_fraction"]],
        "pathologist stromal fraction only": path[["Stromal_fraction"]],
        "CPTAC deconvolution (stroma+immune)": dec[["stromal_deconv", "immune_deconv"]],
        "ESTIMATE (stroma+immune)": dec[["StromalScore_ESTIMATE", "ImmuneScore_ESTIMATE"]],
    }, path, dec


def cov_ccrcc():
    raw = pd.read_excel(f"{DD}/subtype_ccrcc_raw.xlsx", sheet_name="ESTIMATE scores",
                        header=3, index_col=0)
    est = raw.T                                          # aliquots x scores
    est = est.apply(pd.to_numeric, errors="coerce")
    xw = pd.read_csv("server_export/scripts/aliquot_to_case_tumor.tsv", sep="\t")
    xw = xw[xw["sample_type"].str.contains("Tumor", case=False, na=False)]
    a2c = dict(zip(xw["aliquot_id"].astype(str), xw["case_id"].map(norm_id)))
    est["case"] = [a2c.get(str(i)) for i in est.index]
    est = est.dropna(subset=["case"]).groupby("case").mean()
    return {"ESTIMATE (stroma+immune)": est[["StromalScore", "ImmuneScore"]],
            "ESTIMATE tumour purity only": est[["TumorPurity"]]}


def cov_ucec():
    d = pd.read_csv(f"{DD}/subtype_ucec.txt", sep="\t", encoding="latin-1", low_memory=False)
    if "Proteomics_Tumor_Normal" in d.columns:
        d = d[d["Proteomics_Tumor_Normal"].astype(str).str.lower().str.startswith("tumor")]
    d["case"] = d["Proteomics_Participant_ID"].map(norm_id)
    d = d.set_index("case")
    # Tumor_purity is NOT a purity: it holds the string "Normal"/"Tumor". The
    # numeric cancer-cell fraction is Purity_Cancer. Pointing at the wrong column
    # silently produced an all-NaN covariate set (n=0) in every earlier run.
    cols = ["ESTIMATE_StromalScore", "ESTIMATE_ImmuneScore", "Purity_Cancer",
            "Purity_Stroma", "Purity_Immune"]
    d = d[cols].apply(pd.to_numeric, errors="coerce")
    d = d[~d.index.duplicated()]
    return {"ESTIMATE (stroma+immune)": d[["ESTIMATE_StromalScore", "ESTIMATE_ImmuneScore"]],
            "tumour purity only": d[["Purity_Cancer"]],
            "purity-derived stroma+immune": d[["Purity_Stroma", "Purity_Immune"]]}


def cov_gbm():
    d = pd.read_excel(f"{DD}/subtype_gbm_raw.xlsx", sheet_name="additional_annotations")
    if "sample_type" in d.columns:
        t = d[d["sample_type"].astype(str).str.lower().str.contains("tumor")]
        if len(t):
            d = t
    d["case"] = d["case"].map(norm_id)
    d = d.set_index("case")[["xcell_stroma_score", "xcell_immune_score"]] \
         .apply(pd.to_numeric, errors="coerce")
    d = d[~d.index.duplicated()]
    return {"xCell (stroma+immune)": d}


def cov_luad():
    """Gillette et al. 2020 Table S1 composition, via cBioPortal luad_cptac_2020.

    Same deconvolution tier as CCRCC and UCEC -- ESTIMATE stromal/immune plus a
    purity estimate -- not the pathologist read that only PDAC has. patientId in
    that table is already the C3L/C3N case identifier.
    """
    d = pd.read_csv(f"{DD}/composition_luad_cbioportal.csv")
    d["case"] = d["patientId"].map(norm_id)
    d = d.set_index("case")
    d = d[~d.index.duplicated()]

    def num(cols):
        return d[cols].apply(pd.to_numeric, errors="coerce")

    return {"ESTIMATE (stroma+immune)": num(["ESTIMATE STROMALSCORE",
                                             "ESTIMATE IMMUNESCORE"]),
            "ESTIMATE tumour purity only": num(["TUMOR_PURITY_BYESTIMATE_RNASEQ"]),
            "TSNet purity only": num(["TSNET PURITY"])}


COV = {"pdac": None, "ccrcc": cov_ccrcc, "ucec": cov_ucec, "gbm": cov_gbm,
       "luad": cov_luad}

# ----------------------------------------------------------------------- run
rows, agree = [], []
pdac_sets, pdac_path, pdac_dec = cov_pdac()
COV["pdac"] = lambda: pdac_sets

for c in ["ccrcc", "luad", "ucec", "gbm", "pdac"]:
    sc = load_scores(c)
    sets = COV[c]()
    for _lab, _C in sets.items():
        assert _C.notna().any(axis=None), (
            f"{c}/{_lab}: covariate block is entirely missing -- a wrong column name "
            f"would otherwise be written out as n=0 NaN rows")
    for label, C in sets.items():
        j = sc.join(C, how="inner")
        cols = list(C.columns)
        for s in SCORES:
            x = j[f"{s}_morph"].astype(float).values
            y = j[f"{s}_meas"].astype(float).values
            Cm = j[cols].astype(float).values
            r_raw, n0 = partial_r(x, y, np.zeros((len(x), 0)))
            r_par, n = partial_r(x, y, Cm)
            # single-covariate correlations use the first covariate as the headline
            rho_mc, _ = spearman(x, Cm[:, 0])
            rho_yc, _ = spearman(y, Cm[:, 0])
            rows.append(dict(cohort=c.upper(), score=s, covariates=label, n=n,
                             rho_morph_cov=rho_mc, rho_meas_cov=rho_yc,
                             r_raw=r_raw, r_partial=r_par,
                             retained_pct=100 * r_par / r_raw if r_raw else np.nan))
            print(f"{c.upper():6s} {s:13s} {label:44s} n={n:3d}  "
                  f"morph~cov {rho_mc:+.2f}  meas~cov {rho_yc:+.2f}  "
                  f"r {r_raw:.2f} -> {r_par:.2f} ({100*r_par/r_raw:.0f}% retained)")

# ---------------------------------------------- PDAC: transcriptomic vs pathologist
print("\nPDAC: how well do transcriptomic composition estimates track the pathologist?")
pairs = [("StromalScore_ESTIMATE", "Stromal_fraction"),
         ("stromal_deconv", "Stromal_fraction"),
         ("ImmuneScore_ESTIMATE", "Inflammation_fraction"),
         ("immune_deconv", "Inflammation_fraction"),
         ("epithelial_cancer_deconv", "Neoplastic_cellularity")]
jj = pdac_dec.join(pdac_path, how="inner")
for t, p in pairs:
    r, n = spearman(jj[t].values, jj[p].values)
    agree.append(dict(transcriptomic=t, pathologist=p, rho=r, n=n))
    print(f"  {t:28s} vs {p:24s} rho={r:+.2f} (n={n})")

pd.DataFrame(rows).to_csv(f"{DD}/composition_pathologist.csv", index=False)
pd.DataFrame(agree).to_csv(f"{DD}/composition_proxy_agreement.csv", index=False)
print(f"\nwrote {DD}/composition_pathologist.csv and composition_proxy_agreement.csv")
