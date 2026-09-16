#!/usr/bin/env python3
r"""Is the UCEC site effect really explained by the scanning operator, or by degrees of freedom?

WHY THIS SCRIPT EXISTS
scanner_site_test.py found that adjusting for aperio.User collapses UCEC's
n_slides -> translation-residual association from r = +0.585 to +0.042 (p = 0.68).
That is complete attenuation, where accrual country managed only 0.585 -> 0.31 and
OCT embedding 0.585 -> 0.25. Taken at face value it would say the site effect is
the scanning operator.

It cannot be taken at face value. aperio.User has 28 levels among 102 UCEC cases,
so the adjustment spends 27 degrees of freedom on n = 100. A 27-df categorical
covariate will attenuate a lot of things. Before any of this reaches the manuscript
the attenuation has to be separated from the arithmetic, and the way to do that is a
placebo: permute the operator labels across cases, preserving the number of levels
and the size of each, and re-run the identical adjustment. If a random partition of
the same shape also drives 0.585 to near zero, the finding is an artefact. If the
placebo null stays near 0.585, the operator is really carrying it.

Three further things have to be true before the operator can be called the mechanism,
and each is tested here:
  * the operator must be associated with BOTH the slide count and the residual --
    a confounder needs both legs, and a variable associated with neither cannot
    attenuate anything except by chance;
  * it must not simply be a relabelling of accrual country and OCT embedding, the
    two variables already in the manuscript. If operator is nested within country,
    this is a finer measurement of a known effect and should be written that way,
    not as a new mechanism;
  * and critically for the paper: adjusting for the operator must NOT destroy the
    morphology signal itself. If it does, UCEC's headline result is an operator
    artefact and the cohort has a much bigger problem than Limitations (vi) admits.

PDAC gets the same treatment for its own finding: the single genuine scanner
contrast in the dataset (SS7559, 37 slides over 12 cases, and the only slides at
MPP 0.5031) is associated with the ER-secretion residual at d = -0.77, p = 0.013.
Scanner, MPP and possibly scan date are collinear there by construction, so the
question is whether anything survives separating them, and again whether the
morphology signal is touched.

Console output is ascii-folded: the Windows console is GBK.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

EX = Path("server_export")
DD = Path("figures/figdata")
CASE_RE = re.compile(r"(C3[A-Z]-\d{5})")
NPERM = 4000
RNG = np.random.default_rng(20260909)


def case_of(s):
    m = CASE_RE.search(str(s))
    return m.group(1) if m else None


def load_slides():
    s = pd.read_csv(EX / "wsi_scanner_map.tsv", sep="\t", dtype=str, keep_default_na=False)
    s = s[s.error == ""].copy()
    s["cohort"] = s["cohort"].str.upper()
    s["mpp"] = pd.to_numeric(s["mpp_x"], errors="coerce")
    d = pd.to_datetime(s["scan_date"].str.strip(), format="%m/%d/%y", errors="coerce")
    s["date"] = d
    return s


def per_case(s, cohort):
    g = s[(s.cohort == cohort) & (s.case_id != "")].copy()
    has_t = g.groupby("case_id")["sample_type"].transform(lambda x: (x == "Primary Tumor").any())
    g = g[(~has_t) | (g.sample_type == "Primary Tumor")]
    out = g.groupby("case_id").agg(
        n_banked=("file", "size"),
        n_users=("user", "nunique"),
        user=("user", lambda x: x.mode().iloc[0] if len(x.mode()) else ""),
        scanner=("scanscope_id", lambda x: x.mode().iloc[0] if len(x.mode()) else ""),
        mpp=("mpp", "median"),
        date=("date", "median"),
    )
    return out


def scores(c):
    s = pd.read_csv(DD / f"scores_{c.lower()}.csv").rename(columns={"Unnamed: 0": "case"})
    s["case"] = s["case"].map(lambda x: "-".join(str(x).replace(".", "-").split("-")[:2]))
    return s.set_index("case")


def emb_nslides(c, fn):
    d = np.load(EX / "embeddings" / fn, allow_pickle=True)
    ids = pd.Series([str(x) for x in d["slide_ids"]]).map(case_of).dropna()
    return ids.value_counts().rename("n_emb")


def ucec_clinical(idx):
    cl = pd.read_csv(DD / "subtype_ucec.txt", sep="\t", encoding="latin-1", low_memory=False)
    if "Proteomics_Tumor_Normal" in cl.columns:
        cl = cl[cl["Proteomics_Tumor_Normal"].astype(str).str.lower().str.startswith("tumor")]
    cl["case"] = cl["Proteomics_Participant_ID"].map(case_of)
    cl = cl.dropna(subset=["case"]).drop_duplicates("case").set_index("case")
    return cl.reindex(idx)


# --------------------------------------------------------------------- machinery
def dummies(labels):
    return pd.get_dummies(pd.Series(labels).astype(str), drop_first=True).values.astype(float)


def partial_r(x, y, Z=None):
    """Pearson r between x and y after linearly removing Z from both."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if Z is None or (hasattr(Z, "shape") and Z.shape[1] == 0):
        return st.pearsonr(x, y)
    A = np.column_stack([np.ones(len(x)), Z])
    rx = x - A @ np.linalg.lstsq(A, x, rcond=None)[0]
    ry = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
    return st.pearsonr(rx, ry)


def placebo_null(x, y, labels, nperm=NPERM):
    """Partial r after adjusting for a RANDOM partition of the same shape.

    This is the whole point. Permuting the labels keeps the level count and every
    level's size, so the adjustment costs exactly the same degrees of freedom, but
    destroys any real link to x and y. The spread of this null is how much
    attenuation the arithmetic alone can buy.
    """
    lab = np.array(labels, dtype=object)
    out = np.empty(nperm)
    for i in range(nperm):
        out[i] = partial_r(x, y, dummies(RNG.permutation(lab)))[0]
    return out


def cramers_v(a, b):
    t = pd.crosstab(pd.Series(a).astype(str), pd.Series(b).astype(str))
    if t.shape[0] < 2 or t.shape[1] < 2:
        return np.nan
    chi2 = st.chi2_contingency(t, correction=False)[0]
    n = t.values.sum()
    return np.sqrt((chi2 / n) / (min(t.shape) - 1))


def nestedness(inner, outer):
    """Fraction of inner-level members that fall in that level's modal outer value.
    1.00 means inner is perfectly nested within outer (a finer relabelling of it)."""
    df = pd.DataFrame({"i": pd.Series(inner).astype(str), "o": pd.Series(outer).astype(str)})
    df = df[df.o.str.lower().isin(["nan", "", "unknown"]) == False]  # noqa: E712
    if not len(df):
        return np.nan
    return df.groupby("i")["o"].apply(lambda s: s.value_counts().iloc[0] / len(s)).mean()


def kw(y, labels, minn=5):
    lab = pd.Series(labels).astype(str)
    vc = lab.value_counts()
    lev = vc[vc >= minn].index
    gs = [np.asarray(y)[(lab == L).values] for L in lev]
    gs = [g[np.isfinite(g)] for g in gs]
    gs = [g for g in gs if len(g) >= minn]
    if len(gs) < 2:
        return np.nan, np.nan, 0
    s, p = st.kruskal(*gs)
    return s, p, len(gs)


def banner(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


# -------------------------------------------------------------------- UCEC block
def ucec(sl):
    banner("UCEC -- is the operator the mechanism, or 27 degrees of freedom?")
    h = per_case(sl, "UCEC")
    sc = scores("UCEC")
    ne = emb_nslides("UCEC", "slide_mean_ucec_emb.npz")
    idx = sorted(set(h.index) & set(sc.index) & set(ne.index))
    h, sc, ne = h.loc[idx], sc.loc[idx], ne.loc[idx]
    cl = ucec_clinical(idx)

    x = ne.values.astype(float)
    y = pd.to_numeric(sc["translation_meas"], errors="coerce").values
    U = h["user"].astype(str).values
    nlev = pd.Series(U).nunique()

    print(f"n = {len(idx)} cases; operator levels = {nlev}; "
          f"df spent = {nlev - 1}; residual df = {len(idx) - nlev - 1}")
    vc = pd.Series(U).value_counts()
    print(f"  cases per operator: median {vc.median():.0f}, max {vc.max()}, "
          f"singletons {int((vc == 1).sum())}/{nlev}")
    print(f"  cases whose slides span >1 operator: {int((h.n_users > 1).sum())}/{len(h)} "
          f"(the per-case operator is a mode, so this is the summary's error rate)")
    print(f"  banked vs embedded slide count: rho = "
          f"{st.spearmanr(h.n_banked.values, x)[0]:+.3f}  "
          f"(low rho means the header table counts a different thing than the analysis did)")

    r0 = partial_r(x, y, None)
    rU = partial_r(x, y, dummies(U))
    print(f"\n  raw               r = {r0[0]:+.3f}  p = {r0[1]:.2e}")
    print(f"  adjusted for user r = {rU[0]:+.3f}  p = {rU[1]:.2e}")

    null = placebo_null(x, y, U)
    p_lo = (np.sum(null <= rU[0]) + 1) / (len(null) + 1)
    print(f"\n  PLACEBO ({NPERM} random partitions of identical shape):")
    print(f"    null partial r: median {np.median(null):+.3f}, "
          f"5th pct {np.percentile(null, 5):+.3f}, min {null.min():+.3f}")
    print(f"    observed {rU[0]:+.3f} sits at p = {p_lo:.4f} of that null")
    verdict = ("REAL: a random partition of the same shape does not do this"
               if p_lo < 0.05 else
               "ARTEFACT: the arithmetic alone reproduces the attenuation")
    print(f"    -> {verdict}")

    print("\n  both legs of the confounder (an operator that predicts neither cannot confound):")
    s1, p1, k1 = kw(x, U)
    s2, p2, k2 = kw(y, U)
    print(f"    operator -> n_slides          H = {s1:7.2f}  p = {p1:.2e}  ({k1} levels used)")
    print(f"    operator -> translation resid H = {s2:7.2f}  p = {p2:.2e}  ({k2} levels used)")

    print("\n  is the operator just a relabelling of what is already in the manuscript?")
    rows = []
    for col in ["Country", "Proteomics_OCT"]:
        if col not in cl.columns:
            print(f"    {col}: not in clinical table")
            continue
        v = cl[col].astype(str).values
        print(f"    operator vs {col:16s} Cramers V = {cramers_v(U, v):.3f}   "
              f"nested-in-{col} = {nestedness(U, v):.3f}   "
              f"({col} levels = {pd.Series(v).nunique()})")
        rows.append((col, v))

    print("\n  who explains whom (sequential adjustment):")
    for col, v in rows:
        rc = partial_r(x, y, dummies(v))
        rcu = partial_r(x, y, np.column_stack([dummies(v), dummies(U)]))
        print(f"    adjust {col:16s} alone r = {rc[0]:+.3f} ; then add operator r = {rcu[0]:+.3f}")
    if rows:
        Zc = np.column_stack([dummies(v) for _, v in rows])
        print(f"    adjust country+OCT together r = {partial_r(x, y, Zc)[0]:+.3f} ; "
              f"then add operator r = "
              f"{partial_r(x, y, np.column_stack([Zc, dummies(U)]))[0]:+.3f}")

    print("\n  THE PAPER-CRITICAL TEST: does adjusting for operator damage the morphology signal?")
    for sname in [c for c in sc.columns if c.endswith("_meas")]:
        mm = sname.replace("_meas", "_morph")
        if mm not in sc.columns:
            continue
        ym = pd.to_numeric(sc[sname], errors="coerce").values
        xm = pd.to_numeric(sc[mm], errors="coerce").values
        m = np.isfinite(xm) & np.isfinite(ym)
        rr = partial_r(xm[m], ym[m], None)
        ru = partial_r(xm[m], ym[m], dummies(U[m]))
        nl = placebo_null(xm[m], ym[m], U[m], nperm=1500)
        pp = (np.sum(nl <= ru[0]) + 1) / (len(nl) + 1)
        print(f"    {sname:20s} morph vs meas  raw r = {rr[0]:+.3f} (p={rr[1]:.1e})  "
              f"| operator r = {ru[0]:+.3f} (p={ru[1]:.1e})  placebo-null median "
              f"{np.median(nl):+.3f}, p = {pp:.3f}")
    return dict(n=len(idx), levels=nlev, r_raw=r0[0], r_user=rU[0],
                placebo_median=float(np.median(null)), placebo_p=p_lo)


# -------------------------------------------------------------------- PDAC block
def pdac(sl):
    banner("PDAC -- the only genuine scanner contrast in the dataset (SS7559, 12 cases)")
    h = per_case(sl, "PDAC")
    sc = scores("PDAC")
    idx = sorted(set(h.index) & set(sc.index))
    h, sc = h.loc[idx], sc.loc[idx]
    S = (h["scanner"].astype(str) == "SS7559").values
    print(f"n = {len(idx)} cases; SS7559 = {int(S.sum())}, SS1553 = {int((~S).sum())}")
    print(f"  scanner vs MPP: SS7559 median {np.nanmedian(h.mpp.values[S]):.4f}, "
          f"SS1553 median {np.nanmedian(h.mpp.values[~S]):.4f}  "
          f"(collinear by construction -- they cannot be separated)")
    d = h["date"]
    if d.notna().any():
        dd = (d - d.min()).dt.days.astype(float).values
        print(f"  scanner vs scan date: SS7559 median day {np.nanmedian(dd[S]):.0f}, "
              f"SS1553 {np.nanmedian(dd[~S]):.0f}  "
              f"(Mann-Whitney p = {st.mannwhitneyu(dd[S], dd[~S])[1]:.2e})")
    for sname in [c for c in sc.columns if c.endswith("_meas")]:
        y = pd.to_numeric(sc[sname], errors="coerce").values
        m = np.isfinite(y)
        a, b = y[m & S], y[m & ~S]
        if len(a) < 5:
            continue
        u, p = st.mannwhitneyu(a, b)
        pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                         / (len(a) + len(b) - 2))
        dz = (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan
        # does it survive removing scan date, the one thing not perfectly collinear?
        extra = ""
        if d.notna().all():
            rr = partial_r(S[m].astype(float), y[m], dd[m][:, None])
            extra = f"  | date r = {rr[0]:+.3f} (p={rr[1]:.2e})"
        print(f"  {sname:20s} d = {dz:+.3f}  MannWhitney p = {p:.2e}{extra}")
        mm = sname.replace("_meas", "_morph")
        if mm in sc.columns:
            xm = pd.to_numeric(sc[mm], errors="coerce").values
            k = m & np.isfinite(xm)
            r1 = partial_r(xm[k], y[k], None)
            r2 = partial_r(xm[k], y[k], S[k].astype(float)[:, None])
            print(f"    {'morph vs meas':20s} raw r = {r1[0]:+.3f} -> "
                  f"adjusting for scanner r = {r2[0]:+.3f}")


def main():
    sl = load_slides()
    u = ucec(sl)
    pdac(sl)
    banner("bottom line")
    if u["placebo_p"] < 0.05:
        print(f"UCEC: operator adjustment takes r {u['r_raw']:+.3f} -> {u['r_user']:+.3f} where a")
        print(f"      random partition of the same shape leaves it at {u['placebo_median']:+.3f}")
        print(f"      (placebo p = {u['placebo_p']:.4f}). The attenuation is not a df artefact.")
    else:
        print(f"UCEC: the attenuation is within what a random {u['levels']}-level partition")
        print(f"      achieves (placebo median {u['placebo_median']:+.3f}, p = {u['placebo_p']:.3f}).")
        print("      Do NOT write the operator into the manuscript as the mechanism.")


if __name__ == "__main__":
    main()
