#!/usr/bin/env python3
r"""What in the UCEC clinical table travels with the number of slides per patient?

n_slides predicts the UCEC translation residual at rho = +0.48 (and ER-secretion
at -0.30, opposite sign). It is not a disease-burden proxy -- it is uncorrelated
with grade and stage (rho 0.09, 0.08) and conditioning on both leaves the
association untouched (0.59 -> 0.58). It is not a mean-pooling denoising artefact
either -- morphology predicts WORSE, not better, in the patients with more slides.

This screens all 179 columns of the CPTAC UCEC clinical table for association
with n_slides, then asks whether any hit also accounts for the residual link.

This is a SCREEN over 179 columns on ~100 patients: Benjamini-Hochberg is applied
and nothing here is confirmatory. Its purpose is to name a candidate or to
establish that none is visible, so that the manuscript can say which.

Writes figures/figdata/nuisance_ucec_screen.csv
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

EX = Path("server_export")
DD = "figures/figdata"
CASE_RE = re.compile(r"(C3[A-Z]-\d{5})")


def case_of(s):
    m = CASE_RE.match(str(s))
    return m.group(1) if m else None


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


d = np.load(EX / "embeddings" / "slide_mean_ucec_emb.npz", allow_pickle=True)
sl = pd.DataFrame({"slide": [str(x) for x in d["slide_ids"]]})
sl["case"] = sl["slide"].map(case_of)
ns = sl.dropna(subset=["case"]).groupby("case").size().rename("n_slides")

sc = pd.read_csv(f"{DD}/scores_ucec.csv").rename(columns={"Unnamed: 0": "case"})
sc["case"] = sc["case"].map(lambda s: "-".join(str(s).replace(".", "-").split("-")[:2]))
sc = sc.set_index("case")

cl = pd.read_csv(f"{DD}/subtype_ucec.txt", sep="\t", encoding="latin-1", low_memory=False)
if "Proteomics_Tumor_Normal" in cl.columns:
    cl = cl[cl["Proteomics_Tumor_Normal"].astype(str).str.lower().str.startswith("tumor")]
cl["case"] = cl["Proteomics_Participant_ID"].map(case_of)
cl = cl.dropna(subset=["case"]).drop_duplicates("case").set_index("case")

common = sorted(set(ns.index) & set(sc.index) & set(cl.index))
print(f"UCEC: {len(common)} patients with slides + scores + clinical")
y_ns = ns.loc[common].values.astype(float)
y_tr = sc.loc[common, "translation_meas"].astype(float).values

rows = []
for col in cl.columns:
    v = cl.loc[common, col]
    if col in ("case", "Proteomics_Participant_ID"):
        continue
    num = pd.to_numeric(v, errors="coerce")
    if num.notna().sum() >= 40 and num.nunique() > 2:
        m = np.isfinite(num.values) & np.isfinite(y_ns)
        if m.sum() < 40:
            continue
        rho, p = st.spearmanr(num.values[m], y_ns[m])
        kind, n = "numeric", int(m.sum())
    else:
        s = v.astype(str).str.strip()
        s = s[~s.isin(["nan", "NA", "", "Unknown", "unknown"])]
        vc = s.value_counts()
        lev = vc[vc >= 8].index
        if len(lev) < 2:
            continue
        groups = [y_ns[[common.index(i) for i in s[s == L].index]] for L in lev]
        groups = [g for g in groups if len(g) >= 8]
        if len(groups) < 2:
            continue
        try:
            _, p = st.kruskal(*groups)
        except ValueError:
            continue
        rho, kind, n = np.nan, f"categorical({len(groups)})", int(sum(len(g) for g in groups))
    rows.append(dict(column=col, kind=kind, n=n, rho_with_nslides=rho, p=p))

res = pd.DataFrame(rows).dropna(subset=["p"])
res["q"] = bh(res.p.values)
res = res.sort_values("p")
res.to_csv(f"{DD}/nuisance_ucec_screen.csv", index=False)

print(f"\nscreened {len(res)} columns; {int((res.q < 0.05).sum())} pass BH q<0.05\n")
print(f"{'column':46s} {'kind':16s} {'n':>4s} {'rho':>7s} {'p':>9s} {'q':>8s}")
for _, r in res.head(12).iterrows():
    rr = f"{r.rho_with_nslides:+7.2f}" if np.isfinite(r.rho_with_nslides) else "      -"
    print(f"{str(r.column)[:46]:46s} {r.kind:16s} {int(r.n):4d} {rr} {r.p:9.2e} {r.q:8.3f}")

# do the top hits account for the residual association?
print("\nDoes the top hit explain n_slides -> translation residual?")
top = res[res.q < 0.05].head(5)
if not len(top):
    print("  no column survives BH; nothing to condition on")
for _, r in top.iterrows():
    v = pd.to_numeric(cl.loc[common, r.column], errors="coerce").values
    m = np.isfinite(v) & np.isfinite(y_tr) & np.isfinite(y_ns)
    if m.sum() < 40:
        print(f"  {r.column}: too few complete cases")
        continue
    Z = np.column_stack([np.ones(m.sum()), v[m]])
    rx = y_ns[m] - Z @ np.linalg.lstsq(Z, y_ns[m], rcond=None)[0]
    ry = y_tr[m] - Z @ np.linalg.lstsq(Z, y_tr[m], rcond=None)[0]
    raw = st.pearsonr(y_ns[m], y_tr[m])[0]
    adj = st.pearsonr(rx, ry)[0]
    print(f"  {str(r.column)[:44]:44s} r(ns,resid) {raw:+.2f} -> {adj:+.2f} after adjustment (n={m.sum()})")
print(f"\nwrote {DD}/nuisance_ucec_screen.csv")
