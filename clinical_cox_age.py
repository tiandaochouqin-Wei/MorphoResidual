#!/usr/bin/env python3
"""T1-7 follow-up: age-adjusted Cox models for the morphology-only residual scores.

The manuscript states that the GDC clinical record carries an age field that is
entirely unpopulated, and therefore adjusts only for grade and stage. That is true
of demographic.age_at_index (null for all 544 analysed cases) but NOT of
demographic.days_to_birth, which is populated for 542/544 and encodes the same
quantity (age at index = -days_to_birth/365.25). Age is available after all, so
"did you adjust for age?" -- the first question any reviewer asks of a survival
association in oncology -- can now be answered.

Two case sets are used, and they are kept distinct on purpose:

  GATE set   time >= 0 (CCRCC alone has 5 strictly negative follow-up intervals,
             -9 to -6 days, all censored; the manuscript documents excluding them,
             n=98), complete on grade and stage where the cohort records them.
             This reproduces the original univariable and grade+stage models and
             is the trustworthiness check: if these do not match, the age models
             are not reported.

  ANALYSIS   the GATE set further restricted to cases with a derivable age, so
  set        that all four models are fitted on identical rows and the hazard
             ratios are directly comparable across adjustments. Only PDAC loses
             anything to this (2 cases without days_to_birth).

Models: univariable; + grade + stage; + age; + grade + stage + age. GBM admits
only the age adjustment, since GDC records neither grade nor stage for any of its
99 cases.

Conventions follow the original scripts: hazard ratio per SD of the score, l2
penalizer 1e-3, Wald 95% CI.

Writes figures/figdata/clinical_cox_age.csv
"""
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")
DD = "figures/figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
PEN = 1e-3
TOL = 0.02          # HR agreement required of the reproduction gate


def fit(df, covs):
    d = df[["time", "event", "score"] + covs].dropna()
    if len(d) < 20 or d.event.sum() < 5:
        return None
    cph = CoxPHFitter(penalizer=PEN)
    cph.fit(d, duration_col="time", event_col="event")
    ci = cph.confidence_intervals_.loc["score"]
    return dict(n=len(d), events=int(d.event.sum()),
                HR=float(np.exp(cph.params_["score"])),
                lo=float(np.exp(ci.iloc[0])), hi=float(np.exp(ci.iloc[1])),
                p=float(cph.summary.loc["score", "p"]))


def server_table(c):
    if c == "ccrcc":       # CCRCC's Cox models were run locally, into the family file
        f0 = pd.read_csv(f"{DD}/clinical_family_ccrcc.csv")
        t = pd.DataFrame({"uni_HR": f0[f0["var"] == "cox_uni"].set_index("score")["value"],
                          "adj_HR": f0[f0["var"] == "cox_adj"].set_index("score")["value"]})
        t["adj_for"] = np.where(t.adj_HR.notna(), "grade+stage", "-")
        return t
    return pd.read_csv(f"{DD}/clinical_link_{c}_cox.csv").set_index("score")


rows, gate = [], []
for c in COH:
    sc = pd.read_csv(f"{DD}/scores_{c}.csv").rename(columns={"Unnamed: 0": "case"})
    cl = pd.read_csv(f"{DD}/clinical_link_{c}_clinical.csv")
    gd = pd.read_csv(f"{DD}/gdc_clinical_{c}.csv")
    gd["age"] = -gd.days_to_birth / 365.25
    d0 = (sc.merge(cl[["case", "grade", "stage", "time", "event"]], on="case")
            .merge(gd[["case", "age"]], on="case", how="left"))

    has_gs = d0.grade.notna().any() and d0.stage.notna().any()
    ADJ = ["grade", "stage"] if has_gs else []

    base = d0[d0.time.notna() & (d0.time >= 0)]
    gate_set = base.dropna(subset=ADJ) if ADJ else base
    ana = gate_set[gate_set.age.notna()]
    srv = server_table(c)

    for s in [x for x in sc.columns if x.endswith("_morph")]:
        def prep(df):
            d = df.copy()
            d["score"] = (d[s] - d[s].mean()) / d[s].std()
            return d

        # ---- gate: reproduce the published models on the published case set
        gu = fit(prep(gate_set), [])
        ga = fit(prep(gate_set), ADJ) if ADJ else None
        if s in srv.index:
            du = abs(gu["HR"] - float(srv.loc[s, "uni_HR"]))
            ok = du < TOL
            msg = f"{c}:{s} uni {gu['HR']:.3f} vs {srv.loc[s,'uni_HR']:.3f} (d={du:.3f})"
            if ga is not None and str(srv.loc[s, "adj_for"]).strip() == "grade+stage":
                da = abs(ga["HR"] - float(srv.loc[s, "adj_HR"]))
                ok = ok and da < TOL
                msg += f"; adj {ga['HR']:.3f} vs {srv.loc[s,'adj_HR']:.3f} (d={da:.3f})"
            gate.append((ok, msg))

        # ---- analysis: all models on identical rows, so adjustments are comparable
        d = prep(ana)
        for lab, covs in [("univariable", []), ("grade+stage", ADJ),
                          ("age", ["age"]), ("grade+stage+age", ADJ + ["age"])]:
            if lab in ("grade+stage", "grade+stage+age") and not ADJ:
                continue
            r = fit(d, covs)
            if r:
                rows.append(dict(cohort=c.upper(), score=s, model=lab, **r))

res = pd.DataFrame(rows)
res.to_csv(f"{DD}/clinical_cox_age.csv", index=False)

print("=== reproduction gate: local refit vs the published models ===")
bad = [m for ok, m in gate if not ok]
for ok, m in gate:
    print(("  PASS  " if ok else "  FAIL  ") + m)
print(f"  -> {len(gate)-len(bad)}/{len(gate)} pass at |dHR| < {TOL}")
if bad:
    print("  !! age models NOT trustworthy; investigate before using them")

print("\n=== age availability (days_to_birth) among analysed cases ===")
for c in COH:
    gd = pd.read_csv(f"{DD}/gdc_clinical_{c}.csv")
    an = pd.read_csv(f"{DD}/clinical_link_{c}_clinical.csv")
    m = gd[gd.case.isin(an.case)]
    a = -m.days_to_birth / 365.25
    print(f"  {c.upper():6s} {a.notna().sum():3d}/{len(m):3d}  median {a.median():.1f} y "
          f"(range {a.min():.0f}-{a.max():.0f}); age_at_index non-null "
          f"{m.age_at_index.notna().sum()}/{len(m)}")

print("\n=== all models on identical rows (HR per SD of the H&E-only score) ===")
for c in COH:
    g = res[res.cohort == c.upper()]
    if not len(g):
        continue
    print(f"\n{c.upper()}")
    for s, gg in g.groupby("score", sort=False):
        print(f"  {s}")
        for _, r in gg.iterrows():
            print(f"    {r.model:16s} n={int(r.n):3d} ev={int(r.events):3d} "
                  f"HR={r.HR:5.3f} [{r.lo:.2f}, {r.hi:.2f}] p={r.p:.4f}"
                  + (" *" if r.p < 0.05 else ""))
print(f"\nwrote {DD}/clinical_cox_age.csv")
