#!/usr/bin/env python3
"""subtype_classify.py (LOCAL) — Tier-2.5 ⑦: does the H&E-ONLY residual score track the
published CPTAC molecular subtypes?  Per cohort: merge figdata/scores_<c>.csv (morphology-only
family residual scores) with a subtype label table, then
  (1) Kruskal-Wallis of each morph score across subtypes,
  (2) 5-fold CV multinomial logistic regression on the morph scores -> macro OvR AUROC,
      with a label-permutation null (P),
  (3) box + strip figure of the most separating score by subtype.
Label files go in figdata/ (see CONFIG); loader is defensive about LinkedOmics .tsi
(transposed) layout and dotted IDs (C3L.00001 -> C3L-00001)."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as st
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import mrstyle as S

DD = "figdata"
CONFIG = {  # cohort: (label file, id column or None=index, subtype column, scheme name)
    "LUAD": (f"{DD}/subtype_luad.tsi", None, "NMF.consensus", "multi-omic NMF C1-C4"),
    "UCEC": (f"{DD}/subtype_ucec.txt", "Proteomics_Participant_ID", "Genomics_subtype", "TCGA genomic subtype"),
    "CCRCC": (f"{DD}/subtype_ccrcc.csv", None, "immune_subtype", "immune subtype (Clark 2019)"),
    "GBM": (f"{DD}/subtype_gbm.csv", None, "nmf_subtype", "nmf1/2/3 (Wang 2021)"),
    "PDAC": (f"{DD}/subtype_pdac.csv", None, "subtype", "Moffitt classical/basal-like (Cao 2021 annotation)"),
}
MORPH = ["translation_morph", "secretion_ER_morph", "matrisome_morph"]
N_PERM = 300


def norm_id(s):
    s = str(s).strip().replace(".", "-")
    return "-".join(s.split("-")[:2]) if s.startswith("C3") else s


def load_labels(path, id_col, sub_col):
    sep = "\t" if path.endswith((".tsi", ".txt", ".tsv")) else ","
    df = pd.read_csv(path, sep=sep, index_col=0 if id_col is None else None, dtype=str, encoding="latin-1")
    if sub_col not in df.columns and sub_col in df.index:      # LinkedOmics transposed layout
        df = df.T
    if sub_col not in df.columns:
        raise SystemExit(f"{path}: subtype column '{sub_col}' not found. columns={list(df.columns)[:15]} index={list(df.index)[:8]}")
    ids = df.index if id_col is None else df[id_col]
    return pd.Series(df[sub_col].values, index=[norm_id(i) for i in ids]).dropna()


rows = []; panels = []
for c, (path, id_col, sub_col, scheme) in CONFIG.items():
    if not os.path.exists(path):
        print(f"[{c}] label file missing: {path} (skipped)"); continue
    lab = load_labels(path, id_col, sub_col)
    sc = pd.read_csv(f"{DD}/scores_{c.lower()}.csv", index_col=0)
    sc.index = [norm_id(i) for i in sc.index]
    df = sc.join(lab.rename("subtype"), how="inner").dropna(subset=["subtype"])
    df = df[df["subtype"].astype(str).str.lower().isin(["nan", "na", ""]) == False]
    vc = df["subtype"].value_counts(); keep = vc[vc >= 5].index
    df = df[df["subtype"].isin(keep)]
    mcols = [m for m in MORPH if m in df.columns and df[m].notna().sum() > 10]
    print(f"[{c}] {len(df)} patients with labels, subtypes {dict(vc[keep])}")
    if len(keep) < 2:
        continue
    # (1) Kruskal-Wallis per morph score
    kw = {}
    for m in mcols:
        groups = [df.loc[df.subtype == s, m].dropna().values for s in keep]
        kw[m] = st.kruskal(*groups).pvalue
    best = min(kw, key=kw.get)
    # (2) CV multinomial LR on all morph scores
    X = df[mcols].fillna(df[mcols].median()).values; y = df["subtype"].astype(str).values
    classes = sorted(set(y)); yi = np.array([classes.index(v) for v in y])

    def cv_auc(yi, seed=0):
        P = np.zeros((len(yi), len(classes)))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, yi):
            sc_ = StandardScaler().fit(X[tr])
            clf = LogisticRegression(max_iter=2000, C=1.0).fit(sc_.transform(X[tr]), yi[tr])
            pr = clf.predict_proba(sc_.transform(X[te]))
            for j, k in enumerate(clf.classes_):
                P[te, k] = pr[:, j]
        if len(classes) == 2:
            return roc_auc_score(yi, P[:, 1])
        return roc_auc_score(yi, P, multi_class="ovr", average="macro")
    auc = cv_auc(yi)
    rng = np.random.RandomState(0)
    null = np.array([cv_auc(rng.permutation(yi)) for _ in range(N_PERM)])
    p_auc = (np.sum(null >= auc) + 1) / (N_PERM + 1)
    rows.append(dict(cohort=c, scheme=scheme, n=len(df), n_subtypes=len(classes),
                     best_score=best.replace("_morph", ""), kw_p_best=kw[best],
                     cv_auroc=auc, auroc_perm_p=p_auc))
    # keep what the combined figure needs
    order = sorted(keep, key=str)
    panels.append(dict(c=c, scheme=scheme, best=best, order=order, kw=kw[best], auc=auc, p=p_auc,
                       data=[df.loc[df.subtype == s, best].dropna().values for s in order]))

res = pd.DataFrame(rows)
if len(res):
    res.to_csv("subtype_results.csv", index=False); print(res.round(3).to_string(index=False))

# ---- combined 1 x N supplementary figure ----
if panels:
    fig, axes = plt.subplots(1, len(panels), figsize=(2.2 * len(panels), 2.7))
    axes = np.atleast_1d(axes)
    for ax, P, letter in zip(axes, panels, "abcdefgh"):
        c = P["c"]
        bp = ax.boxplot(P["data"], widths=0.55, showfliers=False, patch_artist=True,
                        medianprops=dict(color=S.INK, lw=1.2))
        for b in bp["boxes"]:
            b.set_facecolor(S.ORG[c]); b.set_alpha(0.35); b.set_edgecolor(S.ORG[c])
        for i, d in enumerate(P["data"]):
            ax.scatter(np.random.RandomState(i).normal(i + 1, 0.06, len(d)), d, s=6, color=S.ORG[c],
                       alpha=0.65, edgecolors="none")
        labs = [s.replace(" Immune ", " Imm. ").replace("Metabolic", "Metab.") for s in P["order"]]
        ax.set_xticks(range(1, len(labs) + 1)); ax.set_xticklabels(labs, fontsize=5.4, rotation=30, ha="right")
        ax.set_ylabel(f"H&E-only {P['best'].replace('_morph', '')} residual", fontsize=6.4)
        sig = "*" if P["p"] < 0.05 else ""
        ax.set_title(f"{c}: {P['scheme'].split(' (')[0]}\nKW P={P['kw']:.1g} · AUROC={P['auc']:.2f} (P={P['p']:.3f}){sig}",
                     fontsize=6.2, fontweight="bold", loc="left")
        ax.text(-0.25, 1.12, letter, transform=ax.transAxes, fontsize=9, fontweight="bold")
    fig.suptitle("Does the H&E-only residual score track published molecular subtypes?",
                 fontsize=8, fontweight="bold", y=1.06)
    fig.tight_layout()
    S.save_pub(fig, "Fig_subtype"); plt.close(fig)
    print("wrote Fig_subtype (combined)")
print("READ: cv_auroc>0.5 with auroc_perm_p<0.05 => H&E residual score carries molecular-subtype information.")
