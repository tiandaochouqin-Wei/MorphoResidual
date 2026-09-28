#!/usr/bin/env python3
"""Reproduce V(operator,race) and V(scanner,race), and build the permutation null."""
import numpy as np, pandas as pd, itertools, sys

BASE = "D:/claude/pajinsen/MorphoResidual_paper/"
w = pd.read_csv(BASE + "server_export/wsi_scanner_by_case.tsv", sep="\t")
c = pd.read_csv(BASE + "server_export/scripts/covariates_by_case.tsv", sep="\t")
c["cohort"] = c["cohort"].str.upper()
c["case"] = c["case"].astype(str)

def cramer_v(a, b):
    tab = pd.crosstab(a, b).values.astype(float)
    n = tab.sum()
    if n == 0: return np.nan, np.nan, tab.shape
    exp = np.outer(tab.sum(1), tab.sum(0)) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = np.nansum(np.where(exp > 0, (tab - exp) ** 2 / exp, 0.0))
    r, k = tab.shape
    denom = n * (min(r, k) - 1)
    v = np.sqrt(chi2 / denom) if denom > 0 else np.nan
    # Bergsma bias-corrected V
    phi2 = chi2 / n
    phi2c = max(0.0, phi2 - (r - 1) * (k - 1) / (n - 1))
    rc = r - (r - 1) ** 2 / (n - 1)
    kc = k - (k - 1) ** 2 / (n - 1)
    dc = min(rc, kc) - 1
    vc = np.sqrt(phi2c / dc) if dc > 0 else np.nan
    return v, vc, tab.shape

def cramer_v_fast(codes_a, labels_b, na):
    """codes_a: int array of group ids 0..na-1 ; labels_b: 0/1 array."""
    n = len(labels_b)
    n1 = np.bincount(codes_a, weights=labels_b, minlength=na)
    n0 = np.bincount(codes_a, weights=1 - labels_b, minlength=na)
    tot = n0 + n1
    p1 = labels_b.sum() / n
    exp1 = tot * p1
    exp0 = tot * (1 - p1)
    chi2 = 0.0
    m = exp1 > 0
    chi2 += np.sum((n1[m] - exp1[m]) ** 2 / exp1[m])
    m = exp0 > 0
    chi2 += np.sum((n0[m] - exp0[m]) ** 2 / exp0[m])
    r = int((tot > 0).sum())
    denom = n * (min(r, 2) - 1)
    return np.sqrt(chi2 / denom)

RNG = np.random.RandomState(12345)
NPERM = 20000
out = []
for coh in ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]:
    ww = w[w.cohort == coh].drop_duplicates("case_id").set_index("case_id")
    cc = c[c.cohort == coh].drop_duplicates("case").set_index("case")
    common = sorted(set(ww.index) & set(cc.index))
    d = pd.DataFrame({"user": ww.loc[common, "user"].fillna("__NA_USER__"),
                      "scanner": ww.loc[common, "scanner_id"].fillna("__NA_SCAN__"),
                      "race_raw": cc.loc[common, "race"].astype(str).str.lower().str.strip(),
                      "sex": cc.loc[common, "sex"].astype(str).str.lower().str.strip()})
    d["known"] = ~d.race_raw.isin(["not reported", "unknown", "nan", ""])
    dr = d[d.known].copy()
    dr["nonwhite"] = (dr.race_raw != "white").astype(int)
    nmin = int(min(dr.nonwhite.sum(), (1 - dr.nonwhite).sum()))
    n = len(dr)
    vu, vuc, shp_u = cramer_v(dr.user, dr.nonwhite)
    vs, vsc, shp_s = cramer_v(dr.scanner, dr.nonwhite)
    # sex, for comparison
    ds = d[d.sex.isin(["male", "female"])].copy()
    ds["m"] = (ds.sex == "male").astype(int)
    vu_sex, vuc_sex, shp_usex = cramer_v(ds.user, ds.m)

    # how many operators see exactly one race group
    g = dr.groupby("user").nonwhite.agg(["size", "sum"])
    pure = int(((g["sum"] == 0) | (g["sum"] == g["size"])).sum())
    nops = len(g)

    # permutation null for V(operator, X) with same marginal
    codes, uniq = pd.factorize(dr.user)
    na = len(uniq)
    lab = dr.nonwhite.values.astype(float)
    k1 = int(lab.sum())
    null = np.empty(NPERM)
    pure_null = np.empty(NPERM)
    idx = np.arange(n)
    for i in range(NPERM):
        RNG.shuffle(idx)
        lp = np.zeros(n); lp[idx[:k1]] = 1.0
        null[i] = cramer_v_fast(codes, lp, na)
        n1 = np.bincount(codes, weights=lp, minlength=na)
        tot = np.bincount(codes, minlength=na)
        pure_null[i] = np.sum((tot > 0) & ((n1 == 0) | (n1 == tot)))
    p = (np.sum(null >= vu - 1e-12) + 1) / (NPERM + 1)
    out.append(dict(cohort=coh, n=n, n_minority=nmin, n_operators=nops,
                    med_pts_per_op=float(np.median(g["size"])),
                    V_op_race=vu, V_op_race_corrected=vuc,
                    V_scan_race=vs, V_scan_race_corr=vsc,
                    V_op_sex=vu_sex, V_op_sex_corr=vuc_sex, n_sex=len(ds),
                    pure_ops_obs=pure,
                    null_mean=null.mean(), null_med=np.median(null),
                    null_p05=np.percentile(null, 5), null_p50=np.percentile(null, 50),
                    null_p95=np.percentile(null, 95), null_max=null.max(),
                    null_frac_ge_0p9=float((null >= 0.9).mean()),
                    null_frac_eq_1=float((null >= 0.999).mean()),
                    pure_null_med=float(np.median(pure_null)),
                    pure_null_p95=float(np.percentile(pure_null, 95)),
                    perm_p=p))
    print(f"{coh}: n={n} minority={nmin} ops={nops} V_op_race={vu:.3f} (corr {vuc:.3f}) "
          f"V_scan_race={vs:.3f} | null med {np.median(null):.3f} p95 {np.percentile(null,95):.3f} "
          f"max {null.max():.3f} | perm p={p:.4f} | pure obs {pure}/{nops} null med {np.median(pure_null):.0f}")
    sys.stdout.flush()

df = pd.DataFrame(out)
df.to_csv("C:/Users/wei/AppData/Local/Temp/claude/D--claude-pajinsen/f15c3c55-39e5-4d4f-bc57-c97aa8275075/scratchpad/cramer_null.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
print(df.round(4).to_string())
