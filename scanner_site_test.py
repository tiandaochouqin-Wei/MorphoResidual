#!/usr/bin/env python3
r"""Can the physical acquisition recorded in the SVS headers explain the UCEC site effect?

BACKGROUND
The nuisance-block analysis left one loose end. In UCEC, and only in UCEC, a term
with no tumour morphology in it -- the number of slides banked per patient --
predicts the translation residual at out-of-fold R^2 = 0.30. It is not disease
burden and not a mean-pooling artefact; it tracks accrual country (p = 2.7e-12)
and OCT embedding (p = 3.7e-10), and adjusting for either halves the association
(r = 0.59 -> 0.31 and 0.25) where grade and stage do not touch it. Country and
embedding medium sit outside the operator / quarter / scanner axes the manuscript's
batch suite corrects, so we went to the SVS headers for a per-slide site proxy.

WHAT THE HEADERS ACTUALLY SHOWED (extract_svs_batch5.py, 2007 slides, 552 cases)
Aperio ScanScope ID survived GDC repackaging -- it is NOT blank -- but it is very
nearly constant: 387/390 ccRCC, 247/247 GBM, 501/508 LUAD, 389/393 UCEC and
432/469 PDAC report SS1553, and 551 of 552 cases have every tumour slide on a
single scanner. MPP-x is 0.4942 in every cohort. So the headline variable is
degenerate, and in UCEC it is degenerate to a single level.

That is a logical exclusion, not merely an underpowered test: a covariate that is
constant within a cohort cannot generate a within-cohort association with anything.
Whatever the UCEC site effect is, it is not the scanner.

WHAT IS STILL WORTH TESTING
Three header fields do carry variance and were never examined:
  * scan_date -- an operator/quarter proxy that is measured rather than inferred,
    and the only header field with real spread. If UCEC's high-slide-count patients
    were also scanned in a distinct window, that is a concrete accrual mechanism.
  * aperio.User -- the operator string, if the vendor wrote one.
  * PDAC's SS7559, 37 slides over 12 cases -- the single genuine scanner contrast
    in the whole dataset. Underpowered, but it is the only direct test available,
    and a null there is worth stating with its detectable effect size attached.

So this script does not assume the answer. It first quantifies how degenerate each
header variable is per cohort, then tests only what has variance, and finally
re-runs the UCEC n_slides -> residual association adjusting for each header field
in turn, in the same before/after form as the existing country and OCT numbers so
the results drop straight into the manuscript beside them.

INPUTS   server_export/wsi_scanner_map.tsv      (per slide)
         server_export/wsi_scanner_by_case.tsv  (per case)
         figures/figdata/scores_<cohort>.csv
         server_export/embeddings/slide_mean_<cohort>_emb*.npz
OUTPUTS  figures/figdata/scanner_degeneracy.csv
         figures/figdata/scanner_assoc.csv
         figures/figdata/scanner_ucec_adjust.csv

Console output is ascii-folded: the Windows console is GBK and non-ascii aborts it.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

EX = Path("server_export")
DD = Path("figures/figdata")
COHORTS = ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]
EMB = {"CCRCC": "slide_mean_ccrcc_emb_ccrcc_real.npz", "LUAD": "slide_mean_luad_emb.npz",
       "UCEC": "slide_mean_ucec_emb.npz", "GBM": "slide_mean_gbm_emb.npz",
       "PDAC": "slide_mean_pdac_emb.npz"}
CASE_RE = re.compile(r"(C3[A-Z]-\d{5})")
RNG = np.random.default_rng(0)
NPERM = 20000


def ascii_(x):
    return str(x).encode("ascii", "replace").decode("ascii")


def case_of(s):
    m = CASE_RE.search(str(s))
    return m.group(1) if m else None


# --------------------------------------------------------------------------- load
def load_headers():
    ps, pc = EX / "wsi_scanner_map.tsv", EX / "wsi_scanner_by_case.tsv"
    for p in (ps, pc):
        if not p.exists():
            sys.exit(f"missing {p}\n  transfer it from /public/home/fjhui/ZW/results/")
    sl = pd.read_csv(ps, sep="\t", dtype=str, keep_default_na=False)
    ca = pd.read_csv(pc, sep="\t")
    sl["cohort"] = sl["cohort"].str.upper()
    ca["cohort"] = ca["cohort"].str.upper()
    # Aperio writes MM/DD/YY. Parse that strictly and only then fall back, keeping the
    # column datetime64 throughout: a generic parse of the whole column returns mixed
    # tz-aware objects, which silently degrades the dtype to object and kills .dt.
    raw = sl["scan_date"].astype(str).str.strip()
    d = pd.to_datetime(raw, format="%m/%d/%y", errors="coerce")
    miss = d.isna() & (raw != "") & (raw.str.lower() != "nan")
    if miss.any():
        alt = pd.to_datetime(raw[miss], errors="coerce", utc=True)
        d.loc[miss] = alt.dt.tz_localize(None)
        print(f"  note: {int(miss.sum())} scan_date values not in MM/DD/YY, "
              f"{int(d[miss].notna().sum())} recovered by fallback parse")
    sl["date"] = pd.to_datetime(d)
    sl["mpp"] = pd.to_numeric(sl["mpp_x"], errors="coerce")
    return sl, ca


def load_scores(c):
    s = pd.read_csv(DD / f"scores_{c.lower()}.csv").rename(columns={"Unnamed: 0": "case"})
    s["case"] = s["case"].map(lambda x: "-".join(str(x).replace(".", "-").split("-")[:2]))
    return s.set_index("case")


def embedded_nslides(c):
    """n_slides as the nuisance analysis measured it -- embedded, not banked."""
    p = EX / "embeddings" / EMB[c]
    if not p.exists():
        return None
    d = np.load(p, allow_pickle=True)
    ids = pd.Series([str(x) for x in d["slide_ids"]]).map(case_of).dropna()
    return ids.value_counts().rename("n_emb")


# ------------------------------------------------------------------- degeneracy
def degeneracy(sl):
    """How much variance does each header field actually have, per cohort?

    A field whose modal level covers ~everything cannot explain within-cohort
    variation in anything, however it is tested. Reporting this first keeps the
    subsequent nulls honest: most of them are absence of variance, not absence of
    effect, and the manuscript must not conflate the two.
    """
    fields = ["scanscope_id", "user", "appmag", "make", "model", "vendor"]
    rows = []
    for c in COHORTS:
        g = sl[(sl.cohort == c) & (sl.error == "")]
        if not len(g):
            continue
        for f in fields:
            if f not in g.columns:
                continue
            v = g[f].replace("", "<blank>")
            vc = v.value_counts()
            rows.append(dict(cohort=c, field=f, n=len(v), levels=int(vc.size),
                             modal=ascii_(vc.index[0]), modal_frac=vc.iloc[0] / len(v),
                             second=ascii_(vc.index[1]) if vc.size > 1 else "",
                             second_n=int(vc.iloc[1]) if vc.size > 1 else 0))
        m = g["mpp"].dropna().round(4)
        vc = m.value_counts()
        rows.append(dict(cohort=c, field="mpp_x", n=int(len(m)),
                         levels=int(vc.size), modal=f"{vc.index[0]:.4f}" if vc.size else "",
                         modal_frac=vc.iloc[0] / len(m) if len(m) else np.nan,
                         second=f"{vc.index[1]:.4f}" if vc.size > 1 else "",
                         second_n=int(vc.iloc[1]) if vc.size > 1 else 0))
        d = g["date"].dropna()
        rows.append(dict(cohort=c, field="scan_date", n=int(len(d)),
                         levels=int(d.dt.to_period("M").nunique()),
                         modal=str(d.min().date()) if len(d) else "",
                         modal_frac=np.nan,
                         second=str(d.max().date()) if len(d) else "", second_n=0))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------ tests
def perm_p(a, b, nperm=NPERM):
    """Two-sided permutation test on the difference in means. Small-n safe."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    obs = abs(a.mean() - b.mean())
    pool = np.concatenate([a, b])
    n = len(a)
    hits = 0
    for _ in range(nperm):
        RNG.shuffle(pool)
        hits += abs(pool[:n].mean() - pool[n:].mean()) >= obs - 1e-12
    return (hits + 1) / (nperm + 1)


def cohen_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    s = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                / max(len(a) + len(b) - 2, 1))
    return (a.mean() - b.mean()) / s if s > 0 else np.nan


def mde(n1, n2, alpha=0.05, power=0.80):
    """Smallest |d| this comparison could have caught. A null with no MDE beside it
    is not interpretable when one arm has twelve cases."""
    z = st.norm.ppf(1 - alpha / 2) + st.norm.ppf(power)
    return z * np.sqrt(1 / n1 + 1 / n2)


def per_case_headers(sl):
    """Collapse tumour slides to one row per case, carrying the fields with variance."""
    g = sl[sl.error == ""].copy()
    g = g[g.case_id != ""]
    tum = g.groupby(["cohort", "case_id"])["sample_type"].transform(
        lambda s: (s == "Primary Tumor").any())
    g = g[(~tum) | (g.sample_type == "Primary Tumor")]
    out = g.groupby(["cohort", "case_id"]).agg(
        n_banked=("file", "size"),
        scanner=("scanscope_id", lambda s: s.replace("", np.nan).mode().iloc[0]
                 if s.replace("", np.nan).notna().any() else "<blank>"),
        user=("user", lambda s: s.replace("", np.nan).mode().iloc[0]
              if s.replace("", np.nan).notna().any() else "<blank>"),
        appmag=("appmag", lambda s: s.mode().iloc[0] if len(s.mode()) else ""),
        mpp=("mpp", "median"),
        date=("date", "median"),
    ).reset_index()
    return out


def assoc_tests(hdr, deg):
    """Test only the header fields that have variance, against every residual score."""
    rows = []
    for c in COHORTS:
        h = hdr[hdr.cohort == c].set_index("case_id")
        try:
            sc = load_scores(c)
        except FileNotFoundError:
            continue
        meas = [k for k in sc.columns if k.endswith("_meas")]
        common = sorted(set(h.index) & set(sc.index))
        if len(common) < 30:
            continue
        h, sc = h.loc[common], sc.loc[common]
        for score in meas:
            y = pd.to_numeric(sc[score], errors="coerce").values
            # continuous: scan date, mpp
            for f in ["date", "mpp"]:
                v = h[f]
                x = (v - v.min()).dt.days.astype(float).values if f == "date" \
                    else pd.to_numeric(v, errors="coerce").values
                m = np.isfinite(x) & np.isfinite(y)
                if m.sum() < 30 or np.nanstd(x[m]) == 0:
                    rows.append(dict(cohort=c, score=score, field=f, kind="continuous",
                                     n=int(m.sum()), levels=1, stat=np.nan, p=np.nan,
                                     note="no variance" if m.sum() >= 30 else "too few"))
                    continue
                rho, p = st.spearmanr(x[m], y[m])
                rows.append(dict(cohort=c, score=score, field=f, kind="continuous",
                                 n=int(m.sum()), levels=int(pd.Series(x[m]).nunique()),
                                 stat=rho, p=p, note=""))
            # categorical: scanner, user, appmag
            for f in ["scanner", "user", "appmag"]:
                v = h[f].astype(str)
                vc = v.value_counts()
                lev = vc[vc >= 8].index
                if len(lev) < 2:
                    rows.append(dict(cohort=c, score=score, field=f, kind="categorical",
                                     n=len(v), levels=int(vc.size), stat=np.nan, p=np.nan,
                                     note=f"degenerate: {vc.iloc[0]}/{len(v)} in one level"))
                    continue
                groups = [y[(v == L).values & np.isfinite(y)] for L in lev]
                groups = [g for g in groups if len(g) >= 8]
                if len(groups) < 2:
                    continue
                if len(groups) == 2:
                    stat = cohen_d(groups[0], groups[1])
                    p = perm_p(groups[0], groups[1])
                    note = f"MDE(80%)={mde(len(groups[0]), len(groups[1])):.2f}"
                else:
                    stat, p = st.kruskal(*groups)
                    note = "kruskal"
                rows.append(dict(cohort=c, score=score, field=f, kind="categorical",
                                 n=int(sum(len(g) for g in groups)), levels=len(groups),
                                 stat=stat, p=p, note=note))
    return pd.DataFrame(rows)


def resid(y, Z):
    Z = np.column_stack([np.ones(len(y)), Z])
    return y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]


def ucec_adjust(hdr, sl):
    """Re-run n_slides -> translation residual, adjusting for each header field.

    Reported in the same before/after form as the existing country (0.59 -> 0.31)
    and OCT (-> 0.25) numbers so the three sit in one row of the manuscript. The
    n_slides used is the embedded count the nuisance analysis used, not the banked
    count from the raw directory -- and the two are cross-checked, because if they
    disagree the header table is describing a different quantity.
    """
    c = "UCEC"
    h = hdr[hdr.cohort == c].set_index("case_id")
    sc = load_scores(c)
    ne = embedded_nslides(c)
    common = sorted(set(h.index) & set(sc.index) & set(ne.index))
    h, sc, ne = h.loc[common], sc.loc[common], ne.loc[common]
    nb = h["n_banked"].astype(float).values
    nemb = ne.values.astype(float)
    agree = st.spearmanr(nb, nemb)
    y = pd.to_numeric(sc["translation_meas"], errors="coerce").values
    x = nemb
    base = st.pearsonr(x, y)
    rows = [dict(adjust_for="(none)", n=len(common), r=base[0], p=base[1],
                 note=f"banked vs embedded slide count rho={agree[0]:.3f}")]
    for f, kind in [("date", "continuous"), ("mpp", "continuous"),
                    ("scanner", "categorical"), ("user", "categorical"),
                    ("appmag", "categorical")]:
        v = h[f]
        if kind == "continuous":
            z = ((v - v.min()).dt.days.astype(float).values if f == "date"
                 else pd.to_numeric(v, errors="coerce").values)
            m = np.isfinite(z)
            if m.sum() < 30 or np.nanstd(z[m]) == 0:
                rows.append(dict(adjust_for=f, n=int(m.sum()), r=np.nan, p=np.nan,
                                 note="no variance in UCEC"))
                continue
            Z = z[m][:, None]
        else:
            s = v.astype(str)
            d = pd.get_dummies(s, drop_first=True)
            if d.shape[1] == 0:
                rows.append(dict(adjust_for=f, n=len(s), r=np.nan, p=np.nan,
                                 note=f"constant in UCEC ({s.value_counts().iloc[0]}/{len(s)})"))
                continue
            m = np.ones(len(s), bool)
            Z = d.values.astype(float)
        rx, ry = resid(x[m], Z), resid(y[m], Z)
        r, p = st.pearsonr(rx, ry)
        rows.append(dict(adjust_for=f, n=int(m.sum()), r=r, p=p, note=""))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------- main
def main():
    sl, ca = load_headers()
    print(f"per-slide {len(sl)} rows, per-case {len(ca)} rows")
    bad = (sl.error != "").sum()
    if bad:
        print(f"  !! {bad} slides failed to open:")
        print(sl[sl.error != ""][["cohort", "file", "error"]].head(5).to_string(index=False))

    deg = degeneracy(sl)
    deg.to_csv(DD / "scanner_degeneracy.csv", index=False)
    print("\n=== 1. how much variance does each header field have? ===")
    print("(modal_frac near 1.00 means the field is constant and can explain nothing)")
    show = deg[deg.field.isin(["scanscope_id", "user", "appmag", "mpp_x", "scan_date"])]
    for c in COHORTS:
        g = show[show.cohort == c]
        if not len(g):
            continue
        print(f"\n  {c}")
        for _, r in g.iterrows():
            mf = f"{r.modal_frac:.3f}" if np.isfinite(r.modal_frac) else "   - "
            extra = f"  2nd={r.second}({r.second_n})" if r.second_n else ""
            if r.field == "scan_date":
                print(f"    {r.field:13s} n={r.n:4d} months={r.levels:3d} "
                      f"range {r.modal} .. {r.second}")
            else:
                print(f"    {r.field:13s} n={r.n:4d} levels={r.levels:2d} "
                      f"modal={r.modal[:16]:16s} frac={mf}{extra}")

    hdr = per_case_headers(sl)
    print(f"\nper-case header table: {len(hdr)} cases "
          f"({hdr.groupby('cohort').size().to_dict()})")

    ass = assoc_tests(hdr, deg)
    ass.to_csv(DD / "scanner_assoc.csv", index=False)
    print("\n=== 2. does any header field with variance predict a residual score? ===")
    live = ass[ass.p.notna()].sort_values("p")
    if not len(live):
        print("  nothing testable: every header field is degenerate in every cohort")
    else:
        print(f"{'cohort':7s} {'score':18s} {'field':9s} {'n':>4s} {'stat':>7s} "
              f"{'p':>9s}  note")
        for _, r in live.head(20).iterrows():
            print(f"{r.cohort:7s} {r.score[:18]:18s} {r.field:9s} {int(r.n):4d} "
                  f"{r.stat:+7.3f} {r.p:9.2e}  {r.note}")
        dead = ass[ass.p.isna()]
        if len(dead):
            print(f"\n  {len(dead)} field-by-score tests skipped for lack of variance; "
                  f"reasons: {sorted(set(dead.note))[:3]}")

    print("\n=== 3. UCEC: does any header field account for n_slides -> residual? ===")
    adj = ucec_adjust(hdr, sl)
    adj.to_csv(DD / "scanner_ucec_adjust.csv", index=False)
    print(f"{'adjust for':12s} {'n':>4s} {'r':>7s} {'p':>9s}  note")
    for _, r in adj.iterrows():
        rr = f"{r.r:+7.3f}" if np.isfinite(r.r) else "      -"
        pp = f"{r.p:9.2e}" if np.isfinite(r.p) else "        -"
        print(f"{r.adjust_for:12s} {int(r.n):4d} {rr} {pp}  {r.note}")
    print("\n  for comparison, the clinical-table result already in the manuscript:")
    print("    accrual country  0.59 -> 0.31      OCT embedding  0.59 -> 0.25")
    print("    grade and stage  0.59 -> 0.58  (i.e. no attenuation)")

    print(f"\nwrote {DD}/scanner_degeneracy.csv, scanner_assoc.csv, scanner_ucec_adjust.csv")


if __name__ == "__main__":
    main()
