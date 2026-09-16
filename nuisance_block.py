#!/usr/bin/env python3
r"""T1-5 nuisance block: is the morphology embedding just acquisition and stain?

THE OBJECTION THIS ANSWERS
An ImageNet-pretrained ResNet-50 recovers 79% of UNI's significant set
(SuppTable_encoders). A referee's reading: if a network never trained on
pathology does almost as well, the embedding may be keying on gross stain and
acquisition statistics rather than on morphology. The manuscript currently
answers this only by analogy. This measures it.

THE NUISANCE BLOCK
Per patient, a small block of quantities that carry no morphological information
about the tumour but do describe how the slide was made and stained:
    log tile count      how much tissue was tiled (scan area / magnification)
    n slides            how many slides that patient contributed
    chromatin_od        mean haematoxylin optical density  -- global stain intensity
    eosin_fraction      fraction of eosin-positive pixels  -- global stain balance
    nuclear_density     fraction of haematoxylin-positive pixels
(the last three are the stain-level summaries from pathomics.py, mean-pooled per
slide then averaged per patient; they are exactly the kind of statistic a
non-pathology network could pick up.)

WHAT IS MEASURED
  1. CCA between the nuisance block and the 20 morphology principal components
     that the whole analysis is built on. If the embedding were mostly stain and
     acquisition, the leading canonical correlation would be near 1.
  2. How much of the 20-PC variance a linear model on the nuisance block explains
     (per-PC R^2, and variance-weighted total).
  3. The nuisance block as a predictor of the residual pathway scores, out of
     fold, against the same 20 PCs -- i.e. can stain and acquisition alone
     reproduce what morphology predicts?

Writes figures/figdata/nuisance_block.csv and nuisance_block_cca.csv
"""
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

EX = Path("server_export")
DD = "figures/figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
CASE_RE = re.compile(r"(C3[A-Z]-\d{5})")
STAIN = ["nuclear_density", "chromatin_od", "eosin_fraction"]
SCORES = ["translation_meas", "secretion_ER_meas"]


def case_of(slide_id):
    m = CASE_RE.match(str(slide_id))
    return m.group(1) if m else None


def emb_path(c):
    p = EX / "embeddings" / (f"slide_mean_{c}_emb.npz" if c != "ccrcc"
                             else "slide_mean_ccrcc_emb_ccrcc_real.npz")
    return p if p.exists() else None


def path_path(c):
    hits = glob.glob(str(EX / "pathomics" / f"slide_features_{c}_*.csv"))
    return Path(hits[0]) if hits else None


def cv_r2(X, y, seed=0, lam=1.0, n_splits=5):
    """Same closed-form ridge and folds as the main pipeline."""
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


rows, cca_rows = [], []
for c in COH:
    ep, pp = emb_path(c), path_path(c)
    if ep is None:
        print(f"{c.upper()}: no embedding, skipped")
        continue
    d = np.load(ep, allow_pickle=True)
    sl = pd.DataFrame({"slide": [str(x) for x in d["slide_ids"]],
                       "n_tiles": d["n_tiles"].astype(float)})
    sl["case"] = sl["slide"].map(case_of)
    E = pd.DataFrame(d["matrix"], index=sl["slide"].values)

    # patient-level morphology matrix, exactly as the pipeline builds it
    Epat = E.groupby(sl["case"].values).mean().dropna()
    # nuisance: acquisition
    nui = sl.dropna(subset=["case"]).groupby("case").agg(
        log_tiles=("n_tiles", lambda v: np.log10(v.sum())),
        n_slides=("slide", "size"))
    # nuisance: stain
    if pp is not None:
        ph = pd.read_csv(pp)
        ph["case"] = ph["slide_id"].map(case_of)
        have = [k for k in STAIN if k in ph.columns]
        nui = nui.join(ph.dropna(subset=["case"]).groupby("case")[have].mean())
    common = sorted(set(Epat.index) & set(nui.dropna().index))
    if len(common) < 40:
        print(f"{c.upper()}: only {len(common)} patients with both, skipped")
        continue
    W = PCA(n_components=20, svd_solver="full", random_state=0).fit(Epat.loc[common].values)
    PCs = W.transform(Epat.loc[common].values)
    N = nui.loc[common].values.astype(float)
    ncol = list(nui.columns)

    # ---- 1. CCA between the nuisance block and the 20 morphology PCs.
    # A canonical correlation between 5 variables and 20 PCs on ~100 patients is
    # large by construction, so the observed value is uninterpretable without a
    # null: permute the patient labels of the nuisance block and refit.
    k = min(N.shape[1], 20)

    def cca1_of(Nm, P):
        cc_ = CCA(n_components=1, max_iter=2000).fit(Nm, P)
        u, v = cc_.transform(Nm, P)
        return abs(np.corrcoef(u[:, 0], v[:, 0])[0, 1])

    cc = CCA(n_components=k, max_iter=2000).fit(N, PCs)
    U, V = cc.transform(N, PCs)
    cors = [abs(np.corrcoef(U[:, i], V[:, i])[0, 1]) for i in range(k)]
    rng = np.random.default_rng(0)
    null_cca = np.array([cca1_of(N[rng.permutation(len(N))], PCs) for _ in range(200)])
    p_cca = (1 + (null_cca >= cors[0]).sum()) / (1 + len(null_cca))

    # ---- 2. variance of each PC explained by the nuisance block (in-sample R^2)
    # cross-validated, not in-sample: five predictors on ~100 patients makes an
    # in-sample R^2 optimistic for every PC
    r2 = np.array([cv_r2(N, PCs[:, j]) for j in range(PCs.shape[1])])
    wv = W.explained_variance_ratio_[:20]
    r2_weighted = float((r2 * wv).sum() / wv.sum())

    # ---- 3. can the nuisance block predict the residual scores, out of fold?
    sc = pd.read_csv(f"{DD}/scores_{c}.csv").rename(columns={"Unnamed: 0": "case"})
    sc["case"] = sc["case"].map(lambda s: "-".join(str(s).replace(".", "-").split("-")[:2]))
    sc = sc.set_index("case")
    both = [x for x in common if x in sc.index]
    for s in SCORES:
        if s not in sc.columns or len(both) < 40:
            continue
        y = sc.loc[both, s].astype(float).values
        m = np.isfinite(y)
        idx = [common.index(x) for x in both]
        r2_nui = cv_r2(N[idx][m], y[m])
        r2_pc = cv_r2(PCs[idx][m], y[m])
        rows.append(dict(cohort=c.upper(), score=s, n=int(m.sum()),
                         r2_nuisance=r2_nui, r2_morphology=r2_pc,
                         ratio=r2_nui / r2_pc if r2_pc > 0 else np.nan))
        print(f"{c.upper():6s} {s:18s} n={int(m.sum()):3d}  "
              f"nuisance OOF R2={r2_nui:+.3f}  morphology 20PC R2={r2_pc:+.3f}")

    cca_rows.append(dict(cohort=c.upper(), n=len(common), n_nuisance=N.shape[1],
                         nuisance_cols=";".join(ncol),
                         cca1=cors[0], cca1_null_median=float(np.median(null_cca)),
                         cca1_null_p95=float(np.percentile(null_cca, 95)), p_cca=p_cca,
                         r2_pc1=r2[0], r2_pc_max=float(r2.max()),
                         r2_weighted_20pc=r2_weighted))
    print(f"{c.upper():6s} CCA1={cors[0]:.3f} (null median {np.median(null_cca):.3f}, "
          f"95th {np.percentile(null_cca,95):.3f}, p={p_cca:.3f})  "
          f"nuisance explains {100*r2_weighted:.1f}% of 20-PC variance (CV; "
          f"max single PC {100*r2.max():.1f}%)\n")

pd.DataFrame(rows).to_csv(f"{DD}/nuisance_block.csv", index=False)
pd.DataFrame(cca_rows).to_csv(f"{DD}/nuisance_block_cca.csv", index=False)
print(f"wrote {DD}/nuisance_block.csv and nuisance_block_cca.csv")
