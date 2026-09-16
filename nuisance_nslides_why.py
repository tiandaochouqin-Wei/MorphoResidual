#!/usr/bin/env python3
r"""Why does the number of slides per patient predict the UCEC residual?

nuisance_which.py found that n_slides alone reaches OOF R^2 = 0.300 for the UCEC
translation residual -- higher than the whole nuisance block (0.292) and 76% of
the 20 morphology PCs (0.393). Two readings, with opposite consequences:

  A  DISEASE BURDEN. A larger or higher-stage tumour is sampled more, so more
     slides are banked; burden also shifts the proteome. n_slides is then a
     legitimate correlate of disease, not an artefact, and the morphology link is
     unaffected. Signature: n_slides tracks stage/grade, and conditioning on
     stage/grade removes its predictive power.

  B  DENOISING ARTEFACT. A patient with more slides gets a mean-pooled embedding
     averaged over more tiles, hence a less noisy morphology score, hence a
     higher correlation with the residual. Signature: an INTERACTION, not a main
     effect -- morphology predicts the residual better in patients with many
     slides -- and n_slides should not predict the residual on its own.

Note before running: the observed pattern is already a large main effect
(rho = +0.48), which is A's signature and not B's. This tests both explicitly,
and reports every cohort so that UCEC can be seen in context.

Writes figures/figdata/nuisance_nslides_why.csv
"""
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

EX = Path("server_export")
DD = "figures/figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
CASE_RE = re.compile(r"(C3[A-Z]-\d{5})")


def case_of(s):
    m = CASE_RE.match(str(s))
    return m.group(1) if m else None


def cv_r2(X, y, seed=0, lam=1.0):
    X = np.asarray(X, float).reshape(len(y), -1)
    y = np.asarray(y, float)
    oof = np.empty_like(y)
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1.0
        Xtr, Xte = (X[tr] - mu) / sd, (X[te] - mu) / sd
        ym = y[tr].mean()
        w = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Xtr.T @ (y[tr] - ym))
        oof[te] = Xte @ w + ym
    ss = ((y - y.mean()) ** 2).sum()
    return 1 - ((y - oof) ** 2).sum() / ss if ss > 0 else np.nan


def partial_r(x, y, C):
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(C), axis=1)
    if m.sum() < 25:
        return np.nan, int(m.sum())
    x, y, C = x[m], y[m], C[m]
    Z = np.column_stack([np.ones(len(x)), C])
    rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    return st.pearsonr(rx, ry)[0], int(m.sum())


rows = []
for c in COH:
    ep = EX / "embeddings" / (f"slide_mean_{c}_emb.npz" if c != "ccrcc"
                              else "slide_mean_ccrcc_emb_ccrcc_real.npz")
    if not ep.exists():
        continue
    d = np.load(ep, allow_pickle=True)
    sl = pd.DataFrame({"slide": [str(x) for x in d["slide_ids"]],
                       "n_tiles": d["n_tiles"].astype(float)})
    sl["case"] = sl["slide"].map(case_of)
    E = pd.DataFrame(d["matrix"], index=sl["slide"].values)
    Epat = E.groupby(sl["case"].values).mean().dropna()
    ns = sl.dropna(subset=["case"]).groupby("case").size().rename("n_slides")

    sc = pd.read_csv(f"{DD}/scores_{c}.csv").rename(columns={"Unnamed: 0": "case"})
    sc["case"] = sc["case"].map(lambda s: "-".join(str(s).replace(".", "-").split("-")[:2]))
    sc = sc.set_index("case")
    cl = pd.read_csv(f"{DD}/clinical_link_{c}_clinical.csv")
    cl["case"] = cl["case"].map(lambda s: "-".join(str(s).replace(".", "-").split("-")[:2]))
    cl = cl.set_index("case")[["grade", "stage"]]

    common = sorted(set(Epat.index) & set(ns.index) & set(sc.index))
    if len(common) < 40:
        continue
    nsv = ns.loc[common].values.astype(float)
    PCs = PCA(n_components=20, svd_solver="full", random_state=0).fit_transform(Epat.loc[common].values)
    G = cl.reindex(common)[["grade", "stage"]].values.astype(float)

    for s in ["translation_meas", "secretion_ER_meas"]:
        if s not in sc.columns:
            continue
        y = sc.loc[common, s].astype(float).values
        m = np.isfinite(y)
        if m.sum() < 40:
            continue
        rho_y = st.spearmanr(nsv[m], y[m])[0]
        rho_g = st.spearmanr(nsv, G[:, 0], nan_policy="omit")[0]
        rho_s = st.spearmanr(nsv, G[:, 1], nan_policy="omit")[0]
        r_raw, _ = partial_r(nsv, y, np.zeros((len(y), 0)))
        r_adj, n_adj = partial_r(nsv, y, G)          # A: does stage+grade absorb it?
        # B: interaction -- is morphology a better predictor when there are more slides?
        med = np.median(nsv[m])
        hi, lo = m & (nsv >= med), m & (nsv < med)
        r2_hi = cv_r2(PCs[hi], y[hi]) if hi.sum() >= 35 else np.nan
        r2_lo = cv_r2(PCs[lo], y[lo]) if lo.sum() >= 35 else np.nan
        rows.append(dict(cohort=c.upper(), score=s, n=int(m.sum()),
                         rho_nslides_outcome=rho_y, rho_nslides_grade=rho_g,
                         rho_nslides_stage=rho_s,
                         r_nslides_outcome=r_raw, r_nslides_outcome_adj_gradestage=r_adj,
                         n_adj=n_adj, r2_morph_many_slides=r2_hi, r2_morph_few_slides=r2_lo,
                         n_hi=int(hi.sum()), n_lo=int(lo.sum())))

res = pd.DataFrame(rows)
res.to_csv(f"{DD}/nuisance_nslides_why.csv", index=False)

print(f"{'cohort':7s} {'score':18s} {'rho(ns,y)':>10s} {'rho(ns,grade)':>13s} "
      f"{'rho(ns,stage)':>13s} {'r(ns,y)':>8s} {'adj g+s':>8s}")
for _, r in res.iterrows():
    print(f"{r.cohort:7s} {r.score:18s} {r.rho_nslides_outcome:+10.2f} "
          f"{r.rho_nslides_grade:+13.2f} {r.rho_nslides_stage:+13.2f} "
          f"{r.r_nslides_outcome:+8.2f} {r.r_nslides_outcome_adj_gradestage:+8.2f}")

print(f"\nReading B (denoising) predicts morphology works better with more slides:")
print(f"{'cohort':7s} {'score':18s} {'R2 many':>9s} {'R2 few':>9s} {'n hi/lo':>10s}")
for _, r in res.iterrows():
    print(f"{r.cohort:7s} {r.score:18s} {r.r2_morph_many_slides:+9.3f} "
          f"{r.r2_morph_few_slides:+9.3f} {int(r.n_hi):5d}/{int(r.n_lo):<4d}")
print(f"\nwrote {DD}/nuisance_nslides_why.csv")
