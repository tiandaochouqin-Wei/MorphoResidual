#!/usr/bin/env python3
"""T1-6 decisive check: does the Reactome enrichment of the morphology-residual
proteins survive restriction to fully quantified proteins?

Ribosomal, ER and secretory proteins are the most completely quantified proteins
in any TMT experiment, and significance is 2.3-7.9x more likely among proteins
with no missing values (gene_covariates_*.csv). A reviewer will therefore ask
whether the enrichment reported in the paper is a quantification artefact. This
re-runs the same tested-proteome-background ORA on the subset of proteins with
zero missing values, so that both the significant set and the background are
drawn from the same, uniformly measured universe.

Writes figures/figdata/enrichment_complete_only.csv
"""
import warnings
import numpy as np
import pandas as pd
import gseapy as gp

warnings.filterwarnings("ignore")
DD = "figures/figdata"
GMT = f"{DD}/Reactome_2022.gmt"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]

rows = []
for c in COH:
    d = pd.read_csv(f"{DD}/gene_covariates_{c}.csv")
    d = d[d.pinned_significant.notna()]

    tested_all = d.gene.astype(str).tolist()
    sig_all = d.loc[d.pinned_significant == True, "gene"].astype(str).tolist()  # noqa: E712
    e_all = gp.enrich(gene_list=sig_all, gene_sets=GMT, background=tested_all, outdir=None).results

    cp = d[d.pct_missing == 0]
    tested_c = cp.gene.astype(str).tolist()
    sig_c = cp.loc[cp.pinned_significant == True, "gene"].astype(str).tolist()  # noqa: E712
    e_c = gp.enrich(gene_list=sig_c, gene_sets=GMT, background=tested_c, outdir=None).results

    a = e_all.sort_values("Adjusted P-value")
    b = e_c.sort_values("Adjusted P-value")
    top10a, top10b = set(a.head(10).Term), set(b.head(10).Term)
    sig_a = set(a[a["Adjusted P-value"] < 0.05].Term)
    sig_b = set(b[b["Adjusted P-value"] < 0.05].Term)

    print(f"\n--- {c.upper()}: {len(sig_all)}/{len(tested_all)} sig (all) -> "
          f"{len(sig_c)}/{len(tested_c)} sig (fully quantified only) ---")
    print("  all-tested background   :",
          " | ".join(f"{r['Term'].split(' R-HSA')[0][:36]} q={r['Adjusted P-value']:.1e}"
                     for _, r in a.head(3).iterrows()))
    print("  complete-only background:",
          " | ".join(f"{r['Term'].split(' R-HSA')[0][:36]} q={r['Adjusted P-value']:.1e}"
                     for _, r in b.head(3).iterrows()))
    print(f"  top-10 overlap {len(top10a & top10b)}/10 | q<0.05 terms {len(sig_a)} -> {len(sig_b)} "
          f"({len(sig_a & sig_b)} shared)")

    rows.append(dict(
        cohort=c.upper(), n_sig_all=len(sig_all), n_tested_all=len(tested_all),
        n_sig_complete=len(sig_c), n_tested_complete=len(tested_c),
        top_term_all=a.iloc[0]["Term"].split(" R-HSA")[0], q_all=a.iloc[0]["Adjusted P-value"],
        top_term_complete=b.iloc[0]["Term"].split(" R-HSA")[0], q_complete=b.iloc[0]["Adjusted P-value"],
        top10_overlap=len(top10a & top10b),
        n_sig_terms_all=len(sig_a), n_sig_terms_complete=len(sig_b),
        n_sig_terms_shared=len(sig_a & sig_b),
    ))

pd.DataFrame(rows).to_csv(f"{DD}/enrichment_complete_only.csv", index=False)
print(f"\nwrote {DD}/enrichment_complete_only.csv")
