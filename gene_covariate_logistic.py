#!/usr/bin/env python3
"""T1-6 logistic model: is morphology-predictability explained by TMT technical
covariates, or does the biology survive adjusting for them?

For every tested gene in each cohort, fit

    logit P(pinned_significant) = b0 + b1*pct_missing + b2*mean_log_ratio
                                     + b3*r2_rna + b4*1[cytosolic ribosomal]

pct_missing and mean_log_ratio are the two technical covariates a reviewer will
name (which proteins mass spectrometry quantifies completely, and how abundant
they are relative to the pooled reference channel); r2_rna is the gene's own
baseline predictability, which bounds how much room a morphology increment has.
b4 is then the adjusted association of the ribosomal family, i.e. the enrichment
that is NOT attributable to ribosomal proteins being the most completely
quantified and most abundant proteins on the panel.

n is omitted: n = N*(1 - pct_missing/100) within a cohort, so it is collinear
with pct_missing by construction.

Writes figures/figdata/gene_covariate_logistic.csv
"""
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore")
DD = "figures/figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]

# cytosolic ribosomal proteins: RPL*/RPS*/RPLP*/RPSA, excluding the mitochondrial
# (MRPL/MRPS) and the pseudogene-like RPL*P suffixes
def is_ribo(g):
    g = str(g).upper()
    if g.startswith(("MRPL", "MRPS")):
        return False
    return g == "RPSA" or (g.startswith(("RPL", "RPS")) and not g.endswith("P"))


rows = []
for c in COH:
    d = pd.read_csv(f"{DD}/gene_covariates_{c}.csv")
    d = d.dropna(subset=["incremental_r2", "mean_log_ratio"])
    d = d[d.pinned_significant.notna()].copy()
    d["ribo"] = d.gene.map(is_ribo).astype(float)
    y = (d.pinned_significant == True).astype(int).values  # noqa: E712

    X = d[["pct_missing", "mean_log_ratio", "r2_rna", "ribo"]].astype(float)
    X = sm.add_constant(X)
    m = sm.Logit(y, X).fit(disp=0)
    ci = m.conf_int()

    print(f"\n--- {c.upper()} (n={len(d)} genes, {y.sum()} significant, "
          f"{int(d.ribo.sum())} ribosomal) ---")
    for v in ["pct_missing", "mean_log_ratio", "r2_rna", "ribo"]:
        b, p = m.params[v], m.pvalues[v]
        lo, hi = ci.loc[v]
        unit = " (per 10%)" if v == "pct_missing" else ""
        sc = 10.0 if v == "pct_missing" else 1.0
        print(f"  {v:16s} OR={np.exp(b*sc):6.3f} "
              f"[{np.exp(lo*sc):.3f}, {np.exp(hi*sc):.3f}]  p={p:.2e}{unit}")
        rows.append(dict(cohort=c.upper(), term=v, n_genes=len(d), n_sig=int(y.sum()),
                         n_ribo=int(d.ribo.sum()),
                         beta=b, or_=np.exp(b*sc), ci_lo=np.exp(lo*sc), ci_hi=np.exp(hi*sc),
                         p=p, pseudo_r2=m.prsquared))
    print(f"  McFadden pseudo-R2 = {m.prsquared:.4f}")

    # unadjusted ribosomal OR, for the contrast
    Xu = sm.add_constant(d[["ribo"]].astype(float))
    mu = sm.Logit(y, Xu).fit(disp=0)
    print(f"  ribo unadjusted OR={np.exp(mu.params['ribo']):.3f} "
          f"p={mu.pvalues['ribo']:.2e}  ->  adjusted OR={np.exp(m.params['ribo']):.3f} "
          f"p={m.pvalues['ribo']:.2e}")
    rows.append(dict(cohort=c.upper(), term="ribo_unadjusted", n_genes=len(d),
                     n_sig=int(y.sum()), n_ribo=int(d.ribo.sum()),
                     beta=mu.params["ribo"], or_=np.exp(mu.params["ribo"]),
                     ci_lo=np.exp(mu.conf_int().loc["ribo", 0]),
                     ci_hi=np.exp(mu.conf_int().loc["ribo", 1]),
                     p=mu.pvalues["ribo"], pseudo_r2=mu.prsquared))

pd.DataFrame(rows).to_csv(f"{DD}/gene_covariate_logistic.csv", index=False)
print(f"\nwrote {DD}/gene_covariate_logistic.csv")
