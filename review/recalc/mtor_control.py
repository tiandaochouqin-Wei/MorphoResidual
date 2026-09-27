#!/usr/bin/env python3
"""mtor_control.py (review recalc; READ-ONLY w.r.t. the manuscript, writes only into review/recalc/)

Ribosomal-phosphosite control for the claim "the translation residual is an mTORC1-S6K1 output".
LinkedOmics phosphosite matrices are NOT normalised to host-protein abundance, so a phosphosite on
any ribosomal protein can track ribosomal-protein abundance rather than kinase activity.

Per site-level cohort (LUAD, UCEC, GBM, PDAC):
  1. Spearman rho(measured translation residual, every phosphosite) for sites with >=50 % non-missing
     values among the joined patients; distributions for (a) non-RPS6 ribosomal-protein sites,
     (b) the whole phosphoproteome; percentile rank of each canonical mTORC1 site, the S6 arm score,
     the composite and the 4E-BP1 arm inside each distribution.
  2. Partial Spearman of each canonical site / arm score with the translation residual controlling for
     (A) mean z of non-canonical RPS6 sites (same host protein, where present),
     (B) mean z of all non-RPS6 ribosomal-protein phosphosites (proxy for ribosomal-protein abundance),
     (C) a phospho-derived proliferation proxy (mean z of phosphosites on the 24 proliferation markers
         used by composition_controls.py; the transcriptomic score is not available locally),
     and the reverse partial rho(ribosomal proxy, residual | canonical site).
  3. BH across the 5 cohorts x 6 quantities of SuppTable_mtor (from figures/mtor_results.csv).
CCRCC (gene-level only) gets the gene-level analogue of 1-2.
Loader / parser / scoring functions are copied verbatim from figures/mtor_mechanism.py.
"""
import os, re, json
import numpy as np, pandas as pd
from scipy import stats as st
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
DD = os.path.join(PAPER, "figures", "figdata")
FIGDIR = os.path.join(PAPER, "figures")

# ---------------- copied verbatim from figures/mtor_mechanism.py ----------------
MTOR = {"RPS6": ["S235", "S236", "S240", "S244"], "EIF4EBP1": ["T37", "T46", "S65", "T70"],
        "EIF4B": ["S422"], "RPS6KB1": ["T389", "T390"]}   # T390 = T389 under alternative isoform numbering (GBM)
ERK = {"MAPK1": ["T185", "Y187"], "MAPK3": ["T202", "Y204"]}


def norm_id(s):
    s = str(s).strip().replace(".", "-")
    return "-".join(s.split("-")[:2]) if s.startswith("C3") else s


def parse_site(rid, gene_hint=None):
    """Return (gene, set of residues like 'S235') from a row id; None if unparsable."""
    rid = str(rid)
    residues = set(re.findall(r"[STY]\d{1,5}", rid.upper()))
    gene = None
    for tok in re.split(r"[:_\-\s|/]+", rid):
        if re.fullmatch(r"[A-Z][A-Z0-9]{1,9}", tok) and not tok.startswith(("NP", "NM", "XP")) and not re.fullmatch(r"[STY]\d+", tok):
            gene = tok; break
    if gene is None and gene_hint is not None:
        gene = str(gene_hint)
    return (gene, residues) if gene and residues else None


def load_phospho(path):
    """Return DataFrame sites x samples (float), with a parsed (gene, residues) map."""
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    if "Gene" in df.columns:                                   # PDAC layout: Index | Gene | Peptide | samples
        sites = pd.Series(df.index.astype(str)).str.extract(r"((?:[STY]\d+)+)$")[0].fillna("").values
        df.index = df["Gene"].astype(str).values + "_" + sites
        df = df.drop(columns=[c for c in ["Gene", "Peptide"] if c in df.columns])
    # transpose if samples are rows (index looks like C3L-/C3N-)
    if pd.Series(df.index.astype(str)).str.match(r"^C3[LN]").mean() > 0.5:
        df = df.T
    # drop non-numeric annotation rows/cols
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.loc[df.notna().sum(axis=1) > 0, df.notna().sum(axis=0) > 0]
    df.columns = [norm_id(c) for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    return df


def site_score(df, spec):
    """mean z across available sites in spec; returns (score Series, per-site z DataFrame)."""
    rows = {}
    for rid in df.index:
        p = parse_site(rid)
        if not p:
            continue
        g, res = p
        if g in spec and res & set(spec[g]):
            rows[f"{g} {'/'.join(sorted(res & set(spec[g])))}"] = rid
    if not rows:                                               # gene-level matrix (e.g. CCRCC LinkedOmics)
        for g in spec:
            if g in df.index:
                rows[f"{g} (all sites)"] = g
    if not rows:
        return None, None
    z = df.loc[list(rows.values())].T
    z.columns = list(rows.keys())
    z = (z - z.mean()) / z.std(ddof=0)
    return z.mean(axis=1), z
# ---------------------------------------------------------------------------------

RIBO_RE = re.compile(r"^RP[SL]\d+[A-Z]?$")
RIBO_EXTRA = {"RPLP0", "RPLP1", "RPLP2", "RPSA"}
RPS6_S6K_SITES = {"S235", "S236", "S240", "S244", "S247"}     # S247 is also an S6K site: excluded from "other" RPS6 sites
# proliferation markers used by composition_controls.py (MARKERS["proliferation"])
PROLIF = ["MKI67", "PCNA", "TOP2A", "CCNB1", "CCNB2", "CDK1", "BUB1", "BUB1B", "AURKA", "AURKB", "FOXM1", "TYMS",
          "RRM2", "CENPA", "PLK1", "CCNA2", "CDC20", "KIF11", "TPX2", "UBE2C", "BIRC5", "CENPF", "MCM2", "MCM6"]
MIN_FRAC, MIN_N, NBOOT = 0.5, 10, 2000
SCORE_COL = "translation_meas"   # the column mtor_mechanism.py correlates against


def is_ribo(g):
    return bool(RIBO_RE.match(g)) or g in RIBO_EXTRA


def spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < MIN_N:
        return np.nan, np.nan, int(m.sum())
    r, p = st.spearmanr(a[m], b[m])
    return float(r), float(p), int(m.sum())


def partial_spearman(x, y, C):
    """Partial Spearman: rank x, y and each covariate on complete cases, residualise x and y on the
    covariates (+intercept) by OLS, Pearson of the residuals; t-test with df = n - 2 - k."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    C = np.asarray(C, float)
    if C.ndim == 1:
        C = C[:, None]
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(C), axis=1)
    n, k = int(m.sum()), C.shape[1]
    if n < MIN_N + k:
        return np.nan, np.nan, n
    rx, ry = st.rankdata(x[m]), st.rankdata(y[m])
    RC = np.column_stack([np.ones(n)] + [st.rankdata(C[m, j]) for j in range(k)])
    ex = rx - RC @ np.linalg.lstsq(RC, rx, rcond=None)[0]
    ey = ry - RC @ np.linalg.lstsq(RC, ry, rcond=None)[0]
    r = float(np.corrcoef(ex, ey)[0, 1])
    df = n - 2 - k
    t = r * np.sqrt(df / max(1e-12, 1 - r * r))
    return r, float(2 * st.t.sf(abs(t), df)), n


def pct_rank(bg, v):
    bg = np.asarray(bg, float); bg = bg[np.isfinite(bg)]
    if not np.isfinite(v) or len(bg) == 0:
        return np.nan
    return float(100 * ((bg < v).sum() + 0.5 * (bg == v).sum()) / len(bg))


def zmean(P, ids):
    """per-sample mean z over the site rows `ids` of P (sites x samples); NaN if no ids."""
    if len(ids) == 0:
        return pd.Series(np.nan, index=P.columns)
    Z = P.loc[list(ids)].T.astype(float)
    Z = (Z - Z.mean()) / Z.std(ddof=0)
    return Z.mean(axis=1)


def boot_diff(a, b, y, seed=0):
    """paired bootstrap 95 % CI for rho(a,y) - rho(b,y) over complete cases."""
    a = np.asarray(a, float); b = np.asarray(b, float); y = np.asarray(y, float)
    m = np.isfinite(a) & np.isfinite(b) & np.isfinite(y)
    a, b, y = a[m], b[m], y[m]; n = len(y)
    rng = np.random.default_rng(seed); d = []
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        d.append(st.spearmanr(a[i], y[i])[0] - st.spearmanr(b[i], y[i])[0])
    d = np.asarray(d)
    return dict(diff=float(st.spearmanr(a, y)[0] - st.spearmanr(b, y)[0]),
                ci95=[float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))], n=int(n))


def r(x, nd=3):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), nd)


def dist_summary(rhos, label):
    v = np.asarray(rhos, float); v = v[np.isfinite(v)]
    if len(v) == 0:
        return dict(label=label, n=0)
    return dict(label=label, n=int(len(v)), median=r(np.median(v)), p90=r(np.percentile(v, 90)),
                p95=r(np.percentile(v, 95)), max=r(v.max()), min=r(v.min()),
                frac_gt_0p3=r((v > 0.3).mean()), frac_gt_0p5=r((v > 0.5).mean()))


# =====================================================================================
out = dict(script=os.path.relpath(__file__, PAPER), score_column=SCORE_COL, min_nonmissing_frac=MIN_FRAC,
           min_n=MIN_N, cohorts={}, notes=[])
site_tables, top_tables = [], []
res_prev = pd.read_csv(os.path.join(FIGDIR, "mtor_results.csv"))

for c in ["LUAD", "UCEC", "GBM", "PDAC"]:
    ph = load_phospho(os.path.join(DD, f"phospho_{c.lower()}.txt"))
    sc = pd.read_csv(os.path.join(DD, f"scores_{c.lower()}.csv"), index_col=0)
    sc.index = [norm_id(i) for i in sc.index]
    # composite / arm scores exactly as mtor_mechanism.py (z over all phospho samples, then join)
    mt, mtz = site_score(ph, MTOR)
    s6, _ = site_score(ph, {"RPS6": MTOR["RPS6"]})
    bp, _ = site_score(ph, {"EIF4EBP1": MTOR["EIF4EBP1"]})
    common = [s for s in sc.index if s in ph.columns]
    y = sc.loc[common, SCORE_COL].values.astype(float)
    P = ph[common].copy()
    if P.index.duplicated().any():                     # PDAC: same gene+site from several RefSeq isoforms
        idx = pd.Series(P.index.astype(str)); dup = idx.duplicated(keep=False).values; k = idx.groupby(idx).cumcount().values
        P.index = [f"{a}#{j + 1}" if d else a for a, d, j in zip(idx, dup, k)]
        out["notes"].append(f"{c}: {int(dup.sum())} duplicated gene_site ids (isoforms) were suffixed #1/#2 in this analysis")
    frac = P.notna().mean(axis=1)
    n_prev = int(res_prev.loc[res_prev.cohort == c, "n"].iloc[0])
    print(f"[{c}] joined n={len(common)} (mtor_results.csv n={n_prev}); sites={P.shape[0]}; "
          f"sites>=50% non-missing={(frac >= MIN_FRAC).sum()}")

    # ---- classify + per-site rho ----------------------------------------------------
    recs = []
    for rid in P.index:
        p = parse_site(rid)
        if not p:
            continue
        g, res = p
        canon = g in MTOR and bool(res & set(MTOR[g]))
        if canon:
            cls = "canonical_mTORC1"
        elif g == "RPS6":
            cls = "RPS6_S6K_noncanonical" if res & RPS6_S6K_SITES else "RPS6_other"
        elif is_ribo(g):
            cls = "ribosomal_other"
        elif g in PROLIF:
            cls = "proliferation_marker"
        else:
            cls = "other"
        passes = bool(frac[rid] >= MIN_FRAC)
        if not (passes or cls != "other"):
            continue                       # skip missing-heavy 'other' sites (speed); all special classes kept
        rho, pv, n = spearman(P.loc[rid].values, y)
        recs.append(dict(cohort=c, site_id=rid, gene=g, residues="/".join(sorted(res, key=lambda s: int(s[1:]))),
                         cls=cls, n=n, frac_nonmissing=round(float(frac[rid]), 3), pass50=passes, rho=rho, p=pv))
    T = pd.DataFrame(recs)
    site_tables.append(T)

    bg_ribo = T[(T.cls == "ribosomal_other") & T.pass50].dropna(subset=["rho"])
    bg_all = T[T.pass50].dropna(subset=["rho"])
    canon = T[T.cls == "canonical_mTORC1"].copy()
    canon["pct_rank_in_ribosomal"] = [pct_rank(bg_ribo.rho, v) for v in canon.rho]
    canon["pct_rank_in_proteome"] = [pct_rank(bg_all.rho, v) for v in canon.rho]
    # ribosomal sites that beat each canonical pRPS6 site
    top = bg_ribo.sort_values("rho", ascending=False).head(15).copy()
    top_tables.append(top)

    # ---- arm scores on the joined patients ---------------------------------------------
    S6 = s6.reindex(common).values if s6 is not None else np.full(len(common), np.nan)
    MT = mt.reindex(common).values
    BP = bp.reindex(common).values if bp is not None else np.full(len(common), np.nan)
    rho_s6, p_s6, n_s6 = spearman(S6, y)
    rho_mt, p_mt, n_mt = spearman(MT, y)
    rho_bp, p_bp, n_bp = spearman(BP, y)

    # ---- proxies -----------------------------------------------------------------------
    ribo_ids = bg_ribo.site_id.tolist()
    host_ids = T[(T.cls == "RPS6_other") & T.pass50].site_id.tolist()
    host_ids_any = T[T.cls == "RPS6_other"].site_id.tolist()
    prol_ids = T[(T.cls == "proliferation_marker") & T.pass50].site_id.tolist()
    riboBG = zmean(P, ribo_ids).values
    hostBG = zmean(P, host_ids).values
    prolBG = zmean(P, prol_ids).values
    rho_ribo, p_ribo, n_ribo = spearman(riboBG, y)
    rho_host, p_host, n_host = spearman(hostBG, y)
    rho_prol, p_prol, n_prol = spearman(prolBG, y)
    rho_ribo_s6 = spearman(riboBG, S6)
    rho_prol_ribo = spearman(prolBG, riboBG)

    def block(name, x):
        raw = spearman(x, y)
        pr_ribo = partial_spearman(x, y, riboBG)
        pr_host = partial_spearman(x, y, hostBG) if len(host_ids) else (np.nan, np.nan, 0)
        pr_prol = partial_spearman(x, y, prolBG) if len(prol_ids) else (np.nan, np.nan, 0)
        pr_both = partial_spearman(x, y, np.column_stack([riboBG, prolBG])) if len(prol_ids) else (np.nan, np.nan, 0)
        rev = partial_spearman(riboBG, y, x)
        return dict(target=name, rho_raw=r(raw[0]), p_raw=raw[1], n=raw[2],
                    partial_given_ribosomal_proxy=r(pr_ribo[0]), p_partial_ribosomal=pr_ribo[1], n_partial_ribosomal=pr_ribo[2],
                    partial_given_RPS6_other_sites=r(pr_host[0]), p_partial_host=pr_host[1], n_partial_host=pr_host[2],
                    partial_given_proliferation_proxy=r(pr_prol[0]), p_partial_prolif=pr_prol[1],
                    partial_given_ribosomal_and_proliferation=r(pr_both[0]), p_partial_both=pr_both[1], n_partial_both=pr_both[2],
                    reverse_partial_ribosomal_given_target=r(rev[0]), p_reverse=rev[1])

    partials = [block("S6 arm score (pRPS6 mean z)", S6), block("mTORC1 composite", MT), block("4E-BP1 arm score", BP)]
    for _, row in canon.iterrows():
        partials.append(block(f"site {row.site_id}", P.loc[row.site_id].values))

    bd_s6 = boot_diff(S6, riboBG, y)
    bd_mt = boot_diff(MT, riboBG, y)

    out["cohorts"][c] = dict(
        n_joined=len(common), n_in_mtor_results=n_prev, n_sites_total=int(P.shape[0]),
        n_sites_pass50=int((frac >= MIN_FRAC).sum()),
        arm_scores=dict(S6_arm=dict(rho=r(rho_s6), p=p_s6, n=n_s6,
                                    pct_rank_in_ribosomal=r(pct_rank(bg_ribo.rho, rho_s6), 1),
                                    pct_rank_in_proteome=r(pct_rank(bg_all.rho, rho_s6), 1)),
                        composite=dict(rho=r(rho_mt), p=p_mt, n=n_mt,
                                       pct_rank_in_ribosomal=r(pct_rank(bg_ribo.rho, rho_mt), 1),
                                       pct_rank_in_proteome=r(pct_rank(bg_all.rho, rho_mt), 1)),
                        EBP1_arm=dict(rho=r(rho_bp), p=p_bp, n=n_bp,
                                      pct_rank_in_ribosomal=r(pct_rank(bg_ribo.rho, rho_bp), 1),
                                      pct_rank_in_proteome=r(pct_rank(bg_all.rho, rho_bp), 1))),
        ribosomal_nonRPS6_sites=dict(**dist_summary(bg_ribo.rho, "non-RPS6 ribosomal-protein phosphosites, >=50% non-missing"),
                                     max_site=(bg_ribo.sort_values("rho").iloc[-1].site_id if len(bg_ribo) else None),
                                     n_genes=int(bg_ribo.gene.nunique()),
                                     top10=[dict(site=t.site_id, rho=r(t.rho), n=int(t.n)) for _, t in top.head(10).iterrows()]),
        proteome_all_sites=dist_summary(bg_all.rho, "all phosphosites, >=50% non-missing"),
        canonical_sites=[dict(site=row.site_id, residues=row.residues, n=int(row.n), frac_nonmissing=row.frac_nonmissing,
                              pass50=bool(row.pass50), rho=r(row.rho), p=row.p,
                              pct_rank_in_ribosomal=r(row.pct_rank_in_ribosomal, 1),
                              pct_rank_in_proteome=r(row.pct_rank_in_proteome, 1),
                              n_ribosomal_sites_with_higher_rho=int((bg_ribo.rho > row.rho).sum()) if np.isfinite(row.rho) else None)
                         for _, row in canon.iterrows()],
        proxies=dict(ribosomal_proxy=dict(definition="mean z of non-RPS6 ribosomal-protein phosphosites (>=50% non-missing)",
                                          n_sites=len(ribo_ids), rho_with_translation=r(rho_ribo), p=p_ribo, n=n_ribo,
                                          rho_with_S6_arm=r(rho_ribo_s6[0])),
                     RPS6_other_sites_proxy=dict(definition="mean z of RPS6 sites outside S235/236/240/244/247",
                                                 sites_pass50=host_ids, sites_any=host_ids_any,
                                                 rho_with_translation=r(rho_host), p=p_host, n=n_host),
                     proliferation_proxy=dict(definition="mean z of phosphosites on the 24 composition_controls.py proliferation markers (>=50% non-missing); phospho-derived, NOT the transcriptomic score",
                                              n_sites=len(prol_ids), genes=sorted(T[(T.cls == 'proliferation_marker') & T.pass50].gene.unique().tolist()),
                                              rho_with_translation=r(rho_prol), p=p_prol, n=n_prol,
                                              rho_with_ribosomal_proxy=r(rho_prol_ribo[0]))),
        partial_correlations=partials,
        bootstrap_diff=dict(S6arm_minus_ribosomal_proxy=bd_s6, composite_minus_ribosomal_proxy=bd_mt),
    )
    print(f"   S6 arm rho={rho_s6:+.2f} (pct in ribosomal bg {pct_rank(bg_ribo.rho, rho_s6):.0f}, proteome {pct_rank(bg_all.rho, rho_s6):.0f}); "
          f"ribosomal proxy rho={rho_ribo:+.2f}; partial S6|ribo={partials[0]['partial_given_ribosomal_proxy']}; "
          f"reverse ribo|S6={partials[0]['reverse_partial_ribosomal_given_target']}; "
          f"ribosomal bg n={len(bg_ribo)} median={bg_ribo.rho.median():+.2f} p90={bg_ribo.rho.quantile(.9):+.2f} max={bg_ribo.rho.max():+.2f} ({out['cohorts'][c]['ribosomal_nonRPS6_sites']['max_site']})")

# ---------------- CCRCC gene-level analogue -----------------------------------------------
ph = load_phospho(os.path.join(DD, "phospho_ccrcc.txt"))
sc = pd.read_csv(os.path.join(DD, "scores_ccrcc.csv"), index_col=0); sc.index = [norm_id(i) for i in sc.index]
common = [s for s in sc.index if s in ph.columns]
y = sc.loc[common, SCORE_COL].values.astype(float); P = ph[common]; frac = P.notna().mean(axis=1)
recs = []
for g in P.index:
    if frac[g] < MIN_FRAC:
        continue
    rho, pv, n = spearman(P.loc[g].values, y)
    recs.append(dict(cohort="CCRCC", site_id=g, gene=g, residues="(gene-level)", n=n, frac_nonmissing=round(float(frac[g]), 3), pass50=True,
                     cls="canonical_mTORC1" if g in MTOR else ("ribosomal_other" if is_ribo(g) else ("proliferation_marker" if g in PROLIF else "other")),
                     rho=rho, p=pv))
Tc = pd.DataFrame(recs); site_tables.append(Tc)
bg_ribo = Tc[Tc.cls == "ribosomal_other"].dropna(subset=["rho"]); bg_all = Tc.dropna(subset=["rho"])
rps6 = Tc[Tc.gene == "RPS6"].iloc[0]
riboBG = zmean(P, bg_ribo.site_id.tolist()).values
prolBG = zmean(P, Tc[Tc.cls == "proliferation_marker"].site_id.tolist()).values
x = P.loc["RPS6"].values.astype(float)
pr = partial_spearman(x, y, riboBG); rev = partial_spearman(riboBG, y, x); prb = partial_spearman(x, y, np.column_stack([riboBG, prolBG]))
out["cohorts"]["CCRCC"] = dict(
    n_joined=len(common), level="gene-level phospho aggregate (no site resolution)",
    RPS6_gene=dict(rho=r(rps6.rho), p=rps6.p, n=int(rps6.n), pct_rank_in_ribosomal=r(pct_rank(bg_ribo.rho, rps6.rho), 1),
                   pct_rank_in_proteome=r(pct_rank(bg_all.rho, rps6.rho), 1), n_ribosomal_genes_with_higher_rho=int((bg_ribo.rho > rps6.rho).sum())),
    ribosomal_nonRPS6_genes=dict(**dist_summary(bg_ribo.rho, "non-RPS6 ribosomal-protein gene-level phospho"),
                                 max_gene=bg_ribo.sort_values("rho").iloc[-1].gene,
                                 top10=[dict(gene=t.gene, rho=r(t.rho)) for _, t in bg_ribo.sort_values("rho", ascending=False).head(10).iterrows()]),
    proteome_all_genes=dist_summary(bg_all.rho, "all gene-level phospho aggregates"),
    proxies=dict(ribosomal_proxy=dict(n_genes=len(bg_ribo), rho_with_translation=r(spearman(riboBG, y)[0])),
                 proliferation_proxy=dict(n_genes=int((Tc.cls == "proliferation_marker").sum()), rho_with_translation=r(spearman(prolBG, y)[0]))),
    partial_correlations=[dict(target="RPS6 gene-level phospho", rho_raw=r(rps6.rho), partial_given_ribosomal_proxy=r(pr[0]), p_partial_ribosomal=pr[1],
                               partial_given_ribosomal_and_proliferation=r(prb[0]), p_partial_both=prb[1],
                               reverse_partial_ribosomal_given_target=r(rev[0]), p_reverse=rev[1], n=pr[2])],
    bootstrap_diff=dict(RPS6_minus_ribosomal_proxy=boot_diff(x, riboBG, y)))
print(f"[CCRCC] RPS6 gene rho={rps6.rho:+.2f} pct in ribosomal genes {pct_rank(bg_ribo.rho, rps6.rho):.0f} (n={len(bg_ribo)}, median {bg_ribo.rho.median():+.2f}, max {bg_ribo.rho.max():+.2f} {out['cohorts']['CCRCC']['ribosomal_nonRPS6_genes']['max_gene']}); partial|ribo={pr[0]:+.2f}; reverse={rev[0]:+.2f}")

# ---------------- skeptic's specific numbers ------------------------------------------------
T_all = pd.concat(site_tables, ignore_index=True)
def look(cohort, pat):
    s = T_all[(T_all.cohort == cohort) & T_all.site_id.str.contains(pat, regex=True)]
    return [dict(site=a.site_id, rho=r(a.rho), n=int(a.n), frac=a.frac_nonmissing) for _, a in s.iterrows()]
out["skeptic_check"] = dict(PDAC_RPL13_S140=look("PDAC", r"^RPL13_S140(#\d)?$"), PDAC_RPS6_S240=look("PDAC", r"^RPS6_S240(#\d)?$"),
                            LUAD_RPL4_S295=look("LUAD", r"^RPL4:.*S295"), LUAD_RPS6_S240=look("LUAD", r"^RPS6:.*S240"))

# ---------------- BH over the SuppTable_mtor grid (5 cohorts x 6 quantities) --------------
grid = [("composite (mTORC1 vs measured translation residual)", "rho_mtor_vs_translation_meas", "p_meas"),
        ("S6K1 arm (pRPS6) vs translation residual", "rho_S6arm_vs_translation_meas", "p_s6"),
        ("4E-BP1 arm vs translation residual", "rho_4EBP1arm_vs_translation_meas", "p_bp"),
        ("composite vs H&E-only translation score", "rho_mtor_vs_translation_HE", "p_HE"),
        ("composite vs ER-secretion residual", "rho_mtor_vs_secretion_meas", "p_sec"),
        ("ERK control vs translation residual", "rho_ERK_vs_translation_meas", "p_erk")]
bh = []
for _, row in res_prev.iterrows():
    for lab, rc, pc in grid:
        bh.append(dict(cohort=row.cohort, quantity=lab, rho=float(row[rc]), p=float(row[pc])))
bh = pd.DataFrame(bh)
bh["q_BH"] = multipletests(bh.p.values, method="fdr_bh")[1]
bh["survives_q05"] = bh.q_BH < 0.05
bh = bh.sort_values("p").reset_index(drop=True)
out["BH_supp_table_grid"] = dict(n_tests=int(len(bh)), n_survive_q05=int(bh.survives_q05.sum()),
                                 survivors=[dict(cohort=a.cohort, quantity=a.quantity, rho=r(a.rho), p=a.p, q=r(a.q_BH, 4)) for _, a in bh[bh.survives_q05].iterrows()],
                                 non_survivors=[dict(cohort=a.cohort, quantity=a.quantity, rho=r(a.rho), p=a.p, q=r(a.q_BH, 4)) for _, a in bh[~bh.survives_q05].iterrows()],
                                 note="ERK 'control' rows are negative correlations; a survivor there means a significant NEGATIVE association.")
bh.to_csv(os.path.join(HERE, "mtor_control_BH.csv"), index=False)
print(bh.round(4).to_string(index=False))

out["notes"] += [
    "Phosphosite intensities are LinkedOmics site-level log-ratios NOT normalised to host-protein abundance; the ribosomal proxy (mean z of non-RPS6 ribosomal phosphosites) is a surrogate for ribosomal-protein abundance, not a true protein-level normalisation.",
    "No transcriptomic proliferation score exists in figures/figdata or locally (RNA/protein matrices live on the HPC); the proliferation proxy here is phospho-derived from the same 24 markers used by composition_controls.py.",
    "Canonical sites with <50% non-missing values are still reported (n shown) but are excluded from the background distributions.",
    "Percentile rank = % of background sites with lower rho (ties half-weighted).",
]

# ---------------- write ------------------------------------------------------------------
keep = T_all[T_all.cls != "other"].copy()
top50 = pd.concat([g.sort_values("rho", ascending=False).head(50) for _, g in T_all[T_all.cls == "other"].groupby("cohort")])
keep = pd.concat([keep, top50], ignore_index=True).sort_values(["cohort", "cls", "rho"], ascending=[True, True, False])
keep.to_csv(os.path.join(HERE, "mtor_control.csv"), index=False)
T_all.to_csv(os.path.join(HERE, "mtor_control_allsites.csv.gz"), index=False, compression="gzip")


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


json.dump(clean(out), open(os.path.join(HERE, "mtor_control.json"), "w"), indent=1)
print("wrote mtor_control.json / mtor_control.csv / mtor_control_allsites.csv.gz / mtor_control_BH.csv")
