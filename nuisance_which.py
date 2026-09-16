#!/usr/bin/env python3
r"""Which nuisance variable drives the UCEC result?

nuisance_block.py found that the acquisition/stain block predicts UCEC's
translation residual at OOF R^2 = 0.292, against 0.393 for the full 20 morphology
PCs -- 74%. That number means different things depending on which variable
carries it:

  log_tiles / n_slides    an ACQUISITION artefact: how much tissue was scanned.
                          Worrying in a different way from stain, and fixable in
                          principle by conditioning on scan extent.
  chromatin_od / eosin_fraction / nuclear_density
                          a STAIN or CELLULARITY effect. Cellularity is not a
                          nuisance at all -- nuclear density is a real
                          morphological property that the paper already reports as
                          an interpretable correlate of the residual
                          (sec:pathomics), so if that is the driver the "nuisance"
                          framing is wrong for it.

Reports, per cohort and score: each variable alone, each variable dropped, and
the raw Spearman correlation with the outcome.

Writes figures/figdata/nuisance_which.csv
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
STAIN = ["nuclear_density", "chromatin_od", "eosin_fraction"]
SCORES = ["translation_meas", "secretion_ER_meas"]


def case_of(s):
    m = CASE_RE.match(str(s))
    return m.group(1) if m else None


def cv_r2(X, y, seed=0, lam=1.0, n_splits=5):
    X = np.asarray(X, float).reshape(len(y), -1)
    y = np.asarray(y, float)
    oof = np.empty_like(y)
    for tr, te in KFold(n_splits=n_splits, shuffle=True, random_state=seed).split(X):
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1.0
        Xtr, Xte = (X[tr] - mu) / sd, (X[te] - mu) / sd
        ym = y[tr].mean()
        w = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Xtr.T @ (y[tr] - ym))
        oof[te] = Xte @ w + ym
    ss_res = float(((y - oof) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


rows = []
for c in COH:
    ep = EX / "embeddings" / (f"slide_mean_{c}_emb.npz" if c != "ccrcc"
                              else "slide_mean_ccrcc_emb_ccrcc_real.npz")
    hits = glob.glob(str(EX / "pathomics" / f"slide_features_{c}_*.csv"))
    if not ep.exists() or not hits:
        continue
    d = np.load(ep, allow_pickle=True)
    sl = pd.DataFrame({"slide": [str(x) for x in d["slide_ids"]],
                       "n_tiles": d["n_tiles"].astype(float)})
    sl["case"] = sl["slide"].map(case_of)
    E = pd.DataFrame(d["matrix"], index=sl["slide"].values)
    Epat = E.groupby(sl["case"].values).mean().dropna()
    nui = sl.dropna(subset=["case"]).groupby("case").agg(
        log_tiles=("n_tiles", lambda v: np.log10(v.sum())), n_slides=("slide", "size"))
    ph = pd.read_csv(hits[0])
    ph["case"] = ph["slide_id"].map(case_of)
    have = [k for k in STAIN if k in ph.columns]
    nui = nui.join(ph.dropna(subset=["case"]).groupby("case")[have].mean())
    nui = nui.dropna()

    sc = pd.read_csv(f"{DD}/scores_{c}.csv").rename(columns={"Unnamed: 0": "case"})
    sc["case"] = sc["case"].map(lambda s: "-".join(str(s).replace(".", "-").split("-")[:2]))
    sc = sc.set_index("case")
    common = sorted(set(Epat.index) & set(nui.index) & set(sc.index))
    if len(common) < 40:
        continue
    cols = list(nui.columns)
    N = nui.loc[common].values.astype(float)
    PCs = PCA(n_components=20, svd_solver="full", random_state=0).fit_transform(Epat.loc[common].values)

    for s in SCORES:
        if s not in sc.columns:
            continue
        y = sc.loc[common, s].astype(float).values
        m = np.isfinite(y)
        if m.sum() < 40:
            continue
        full = cv_r2(N[m], y[m])
        morph = cv_r2(PCs[m], y[m])
        for j, v in enumerate(cols):
            alone = cv_r2(N[m][:, [j]], y[m])
            drop = cv_r2(np.delete(N[m], j, axis=1), y[m])
            rho = st.spearmanr(N[m][:, j], y[m])[0]
            rows.append(dict(cohort=c.upper(), score=s, variable=v, n=int(m.sum()),
                             r2_full_block=full, r2_morphology=morph,
                             r2_alone=alone, r2_without=drop,
                             loss_when_dropped=full - drop, spearman_with_outcome=rho))

res = pd.DataFrame(rows)
res.to_csv(f"{DD}/nuisance_which.csv", index=False)

for (c, s), g in res.groupby(["cohort", "score"], sort=False):
    full, morph = g.r2_full_block.iloc[0], g.r2_morphology.iloc[0]
    if full < 0.05:
        continue
    print(f"\n{c} {s}   block R2={full:+.3f}  (morphology 20PC {morph:+.3f})")
    print(f"    {'variable':18s} {'alone':>8s} {'without':>8s} {'loss':>7s} {'rho w/ outcome':>15s}")
    for _, r in g.sort_values("r2_alone", ascending=False).iterrows():
        print(f"    {r.variable:18s} {r.r2_alone:+8.3f} {r.r2_without:+8.3f} "
              f"{r.loss_when_dropped:+7.3f} {r.spearman_with_outcome:+15.2f}")
print(f"\nwrote {DD}/nuisance_which.csv")
