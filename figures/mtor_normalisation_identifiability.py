#!/usr/bin/env python3
r"""
mtor_normalisation_identifiability.py -- MorphoResidual T1-4, corrected.

WHAT HAPPENED
A first pass at this analysis (mtor_protein_normalised.py) normalised each
phosphosite by its host protein, saw the correlation with the translation
residual fall in every cohort, and concluded that the phospho-S6 association is
"substantially an abundance effect". An adversarial audit killed that conclusion,
and it was right. The outcome variable is the mean z-scored protein residual over
the translation family, which is overwhelmingly ribosomal proteins; the
denominator we subtract is a ribosomal protein. Denominator and outcome are
therefore collinear, and attenuation follows by arithmetic rather than by
biology.

WHAT THIS SCRIPT DOES INSTEAD
It stops asking "does the association survive normalisation?" -- which is not
identifiable here -- and asks whether the test can discriminate at all, using two
nulls that the first pass lacked.

  PLACEBO DENOMINATOR. Subtract, instead of the site's own host, each of the ~75
  OTHER cytosolic ribosomal proteins, which bear no relation to the site. If
  normalising by the site's own host is doing something host-specific, the real
  denominator should sit in the tail of that placebo distribution. If it sits in
  the middle, the normalisation is not removing host abundance; it is removing
  outcome variance, and the test is uninformative.

  PERMUTATION PLACEBO. Permute the phosphosite across patients, destroying any
  relation to the outcome, then apply the identical subtraction. Whatever
  correlation survives is what the arithmetic manufactures from noise, and IS the
  correct null for the normalised statistic. It is not zero.

It also fixes three defects the audit found in the first pass:
  * CCRCC's ribosomal background was declared impossible. It is not: at gene level
    the construction (gene-level phospho minus that gene's own protein) is exactly
    the same one used for the CCRCC arm, so the background is built here.
  * The gene-level branch reached EIF4EBP1/EIF4B only if their symbols happened to
    contain an [STY]<digit> substring, which they do not, so those arms were
    silently dropped for CCRCC. Symbols are now matched directly.
  * The percentile compared a MEAN of k arm sites against a null of INDIVIDUAL
    background sites, which is wider by roughly sqrt(k). The null is now
    size-matched: means of k randomly chosen background sites.

Writes figures/figdata/mtor_identifiability.csv
"""
import numpy as np
import pandas as pd
from scipy import stats as st

from mtor_common import MTOR, RIBO_PAT, load_phospho, norm_id, parse_site

DD = "figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
MIN_N = 25
N_PERM = 1000
RNG = np.random.default_rng(0)


def z(a):
    a = np.asarray(a, float)
    sd = np.nanstd(a)
    return (a - np.nanmean(a)) / (sd if sd else 1.0)


def rho(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return st.spearmanr(a[m], b[m])[0] if m.sum() >= MIN_N else np.nan


def arm_score(site_mat, denom=None):
    """Mean z across the arm's sites, optionally after subtracting a denominator."""
    cols = []
    for v in site_mat:
        cols.append(z(v if denom is None else v - denom))
    return np.nanmean(np.vstack(cols), axis=0)


rows = []
for c in COH:
    ph = load_phospho(f"{DD}/phospho_{c}.txt")
    hosts = pd.read_csv(f"{DD}/phospho_hosts_{c}.csv", index_col=0)
    hosts.index = [norm_id(i) for i in hosts.index]
    sc = pd.read_csv(f"{DD}/scores_{c}.csv").rename(columns={"Unnamed: 0": "case"}).set_index("case")
    sc.index = [norm_id(i) for i in sc.index]
    common = sorted(set(ph.columns) & set(hosts.index) & set(sc.index))
    y = sc.loc[common, "translation_meas"].astype(float).values
    gene_level = bool(set(MTOR) & set(ph.index.astype(str)))

    # --- how strongly is the OUTCOME itself a ribosomal-abundance proxy?
    ribo_hosts = [g for g in hosts.columns if RIBO_PAT.match(g)]
    panel = hosts.loc[common, ribo_hosts].mean(axis=1).values
    r_panel = rho(panel, y)
    r_rps6 = rho(hosts.loc[common, "RPS6"].astype(float).values, y) if "RPS6" in hosts else np.nan

    # --- collect the S6K1 arm sites and the ribosomal background sites
    arm, bg = [], {}
    if gene_level:
        if "RPS6" in ph.index and "RPS6" in hosts.columns:
            arm.append(ph.loc["RPS6", common].astype(float).values)
        for rid in ph.index.astype(str):
            if RIBO_PAT.match(rid) and rid != "RPS6" and rid in hosts.columns:
                bg[rid] = ph.loc[rid, common].astype(float).values
    else:
        for rid in ph.index:
            p = parse_site(rid)
            if not p:
                continue
            g, res = p
            if g == "RPS6" and (res & set(MTOR["RPS6"])) and "RPS6" in hosts.columns:
                arm.append(ph.loc[rid, common].astype(float).values)
            elif RIBO_PAT.match(g) and g != "RPS6" and g in hosts.columns:
                bg.setdefault(f"{g}|{'/'.join(sorted(res))}", (g, ph.loc[rid, common].astype(float).values))
    if not arm:
        print(f"{c.upper()}: no S6K1 arm sites; skipped")
        continue

    own = hosts.loc[common, "RPS6"].astype(float).values
    r_raw = rho(arm_score(arm), y)
    r_norm = rho(arm_score(arm, own), y)

    # --- PLACEBO DENOMINATOR: subtract an unrelated ribosomal protein instead
    others = [g for g in ribo_hosts if g != "RPS6"]
    plac = np.array([rho(arm_score(arm, hosts.loc[common, g].astype(float).values), y)
                     for g in others])
    plac = plac[np.isfinite(plac)]
    pct_in_placebo = 100 * (plac < r_norm).mean() if len(plac) else np.nan

    # --- PERMUTATION PLACEBO: the correct null for the normalised statistic
    null = []
    for _ in range(N_PERM):
        perm = [RNG.permutation(v) for v in arm]
        null.append(rho(arm_score(perm, own), y))
    null = np.array([v for v in null if np.isfinite(v)])
    p_vs_null = (1 + (null >= r_norm).sum()) / (1 + len(null))

    # --- background, normalised each to its OWN host, with a size-matched null
    bgn = []
    for k, v in bg.items():
        g, vec = (k, v) if gene_level else v
        h = hosts.loc[common, g].astype(float).values
        r = rho(z(vec - h), y)
        if np.isfinite(r):
            bgn.append(r)
    bgn = np.array(bgn)
    kk = len(arm)
    if len(bgn) >= kk:
        means = np.array([bgn[RNG.choice(len(bgn), kk, replace=False)].mean() for _ in range(2000)])
        pct_size_matched = 100 * (means < r_norm).mean()
    else:
        pct_size_matched = np.nan
    pct_single = 100 * (bgn < r_norm).mean() if len(bgn) else np.nan

    rows.append(dict(cohort=c.upper(), gene_level=gene_level, n=len(common), k_arm=kk,
                     rho_outcome_vs_ribo_panel=r_panel, rho_outcome_vs_rps6_protein=r_rps6,
                     rho_raw=r_raw, rho_norm=r_norm,
                     placebo_denom_mean=plac.mean() if len(plac) else np.nan,
                     placebo_denom_n=len(plac), pct_in_placebo_denom=pct_in_placebo,
                     perm_null_median=np.median(null) if len(null) else np.nan,
                     p_vs_perm_null=p_vs_null,
                     bg_n=len(bgn), bg_median_norm=np.median(bgn) if len(bgn) else np.nan,
                     pct_bg_single=pct_single, pct_bg_size_matched=pct_size_matched))
    print(f"{c.upper():6s} outcome~ribo panel r={r_panel:+.3f} | raw {r_raw:+.3f} -> norm {r_norm:+.3f}"
          f" | placebo denom mean {plac.mean():+.3f} (own at {pct_in_placebo:.0f}th pct)"
          f" | perm null median {np.median(null):+.3f}, p={p_vs_null:.3f}"
          f" | bg n={len(bgn)} median {np.median(bgn) if len(bgn) else float('nan'):+.3f},"
          f" arm at {pct_size_matched:.0f}th pct (size-matched)")

res = pd.DataFrame(rows)
res.to_csv(f"{DD}/mtor_identifiability.csv", index=False)
print(f"\nwrote {DD}/mtor_identifiability.csv")
print("\nREAD THIS AS: if 'own at Xth pct' of the placebo-denominator distribution is")
print("near 50, normalising by the site's own host is doing nothing host-specific, and")
print("the attenuation is the arithmetic of subtracting anything correlated with the")
print("outcome. In that case neither the drop nor the sign of rho_norm is interpretable.")
