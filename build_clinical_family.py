#!/usr/bin/env python3
"""T1-7: assemble the complete five-cohort clinical test family and run
Benjamini--Hochberg over it.

Revision history of this file matters, because the family definition is exactly
what an adversarial referee will check:
  v1  CCRCC only (the one cohort whose clinical file was locally available)
  v2  all five cohorts, 84 tests
  v3  (this) two corrections applied by the audit:
      - CCRCC was missing the grade+stage-adjusted Cox models for the
        ER-secretion and matrisome scores, which every other gradeable cohort
        has. By the stated enumeration rule CCRCC should carry one adjusted Cox
        per score, so the "complete" family was not complete. Added.
      - age-adjusted Cox models are now possible (demographic.days_to_birth is
        populated for 542/544 cases, though demographic.age_at_index is null for
        all of them) and are included, so that a robustness check the manuscript
        reports is also carried in the multiplicity accounting rather than
        exempted from it. The script prints the family both with and without
        them, so the effect of that choice is visible rather than assumed.

Two corrections are reported per test:
  q_cohort  BH within that cohort's own family -- PRIMARY, and the one the main
            text, Table 3, Fig 2d and Supp Fig 2f all use
  q_all     BH over every clinical test in the study -- a second, pooled reading.
            Note it is NOT uniformly more conservative per test: pooling changes
            ranks, so q_all is smaller than q_cohort for a substantial minority
            of rows. The script counts them so the manuscript can say so.

Writes figures/figdata/clinical_family_all.csv
"""
import numpy as np
import pandas as pd

DD = "figures/figdata"
NICE = {"translation": "Transl.", "secretion_ER": "ER-secr.", "matrisome": "Matrisome"}
OUT = {"grade": "grade", "stage": "stage", "survival_logrank": "survival (log-rank)",
       "cox_uni": "Cox univariable", "cox_adj": "Cox grade+stage-adj."}
COX_AGE = {"age": "Cox age-adj.", "grade+stage+age": "Cox grade+stage+age-adj."}
ORDER = {"CCRCC": 0, "LUAD": 1, "UCEC": 2, "GBM": 3, "PDAC": 4}


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def split(score):
    base, kind = score.rsplit("_", 1)
    return NICE[base], "meas." if kind == "meas" else "H&E-only"


rows = []

# ------------------------------------------------------------ CCRCC (run locally)
cc = pd.read_csv(f"{DD}/clinical_family_ccrcc.csv")
for _, r in cc.iterrows():
    fam, mode = split(r["score"])
    rows.append(dict(cohort="CCRCC", family=fam, mode=mode, outcome=OUT[r["var"]],
                     stat=r["value"], p=r["p"], n=int(r["n"]), added=""))

# --------------------------------------------- LUAD / UCEC / GBM / PDAC (server)
for c in ["luad", "ucec", "gbm", "pdac"]:
    C = c.upper()
    a = pd.read_csv(f"{DD}/clinical_link_{c}_assoc.csv")
    for _, r in a.iterrows():
        fam, mode = split(r["score"])
        for var, st, pv in [("grade", "grade_rho", "grade_p"),
                            ("stage", "stage_rho", "stage_p"),
                            ("survival_logrank", "surv_logrank_z", "surv_logrank_p")]:
            if pd.isna(r[pv]):
                continue          # GBM records neither grade nor stage
            rows.append(dict(cohort=C, family=fam, mode=mode, outcome=OUT[var],
                             stat=r[st], p=r[pv], n=int(r["n"]), added=""))
    x = pd.read_csv(f"{DD}/clinical_link_{c}_cox.csv")
    for _, r in x.iterrows():
        fam, mode = split(r["score"])
        rows.append(dict(cohort=C, family=fam, mode=mode, outcome=OUT["cox_uni"],
                         stat=r["uni_HR"], p=r["uni_p"], n=int(r["n"]), added=""))
        if str(r["adj_for"]).strip() not in ("-", "", "nan"):
            rows.append(dict(cohort=C, family=fam, mode=mode, outcome=OUT["cox_adj"],
                             stat=r["adj_HR"], p=r["adj_p"], n=int(r["n"]), added=""))

d = pd.DataFrame(rows)

# ---------------- the Cox models the original runs omitted, plus the age models
# All refitted locally on one case-selection rule; clinical_cox_age.py gates this
# refit against every published model and reproduces them to |dHR| < 0.015.
ca = pd.read_csv(f"{DD}/clinical_cox_age.csv")
have = set(zip(d.cohort, d.family, d["mode"], d.outcome))
for _, r in ca.iterrows():
    fam, mode = split(r["score"])
    if r["model"] == "univariable":
        continue                                   # already carried, from the published run
    if r["model"] == "grade+stage":
        key = (r["cohort"], fam, mode, OUT["cox_adj"])
        if key in have:
            continue                               # published value kept; do not double-count
        rows.append(dict(cohort=r["cohort"], family=fam, mode=mode, outcome=OUT["cox_adj"],
                         stat=r["HR"], p=r["p"], n=int(r["n"]), added="missing-adjusted"))
    else:
        rows.append(dict(cohort=r["cohort"], family=fam, mode=mode, outcome=COX_AGE[r["model"]],
                         stat=r["HR"], p=r["p"], n=int(r["n"]), added="age"))

d = pd.DataFrame(rows)
d["q_all"] = bh(d.p.values)
d["q_cohort"] = np.nan
for C, g in d.groupby("cohort"):
    d.loc[g.index, "q_cohort"] = bh(g.p.values)
d = d.sort_values(["cohort", "family", "mode", "outcome"],
                  key=lambda s: s.map(ORDER) if s.name == "cohort" else s).reset_index(drop=True)
d.to_csv(f"{DD}/clinical_family_all.csv", index=False)

# ------------------------------------------------------------------- reporting
print(f"family size: {len(d)} tests  "
      f"({(d.added=='').sum()} original, {(d.added=='missing-adjusted').sum()} restored "
      f"adjusted-Cox, {(d.added=='age').sum()} age models)")
print(d.groupby("cohort", sort=False).size().reindex(list(ORDER)).to_string(), "\n")

# what would the q values be WITHOUT the age models? (is including them load-bearing?)
noage = d[d.added != "age"].copy()
noage["q_all_noage"] = bh(noage.p.values)
noage["q_cohort_noage"] = np.nan
for C, g in noage.groupby("cohort"):
    noage.loc[g.index, "q_cohort_noage"] = bh(g.p.values)
chk = noage.assign(flip=((noage.q_cohort < 0.05) != (noage.q_cohort_noage < 0.05)))
print(f"[age models] family without them = {len(noage)}; conclusions that flip when they are "
      f"included: {int(chk.flip.sum())}")
if chk.flip.any():
    print(chk[chk.flip][["cohort", "family", "mode", "outcome", "p",
                         "q_cohort_noage", "q_cohort"]].to_string(index=False))

dis = d[((d.q_cohort < 0.05) & (d.q_all >= 0.05)) | ((d.q_cohort >= 0.05) & (d.q_all < 0.05))]
print(f"\n[two levels] tests where the two corrections disagree: {len(dis)}/{len(d)}")
print(dis[["cohort", "family", "mode", "outcome", "p", "q_cohort", "q_all"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print(f"\n[two levels] q_all < q_cohort in {(d.q_all < d.q_cohort).sum()}/{len(d)} rows, "
      f"> in {(d.q_all > d.q_cohort).sum()} -- pooling is NOT uniformly more conservative")

print("\n[survives q_cohort<0.05 -- the primary correction]")
print(d[d.q_cohort < 0.05].sort_values("p")[["cohort", "family", "mode", "outcome", "stat", "p",
                                             "q_cohort", "q_all"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))

print("\n[grade tests, H&E-only -- what Fig 2d and Supp Fig 2f plot]")
g = d[(d["mode"] == "H&E-only") & (d.outcome == "grade")]
print(g[["cohort", "family", "stat", "p", "q_cohort", "q_all"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
