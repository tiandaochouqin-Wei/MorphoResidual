#!/usr/bin/env python3
"""enrichment_c2_levels.py -- regenerate Table 2 (tab:c2) under the
tested-proteome background instead of the superseded whole-genome Enrichr
background.

Table 2 in main.tex compares Reactome enrichment-family survival across three
levels of acquisition-batch control (baseline / batch-corrected / strictest)
for each of the five CPTAC cohorts. The version currently in the manuscript
used the original whole-genome Enrichr background (server_export/scripts/
enrichment_check.py, see server_export/scripts/enrichment_survival.tsv for
that run's raw output). Limitations (ii) (main.tex ~L1229-1231) flags that a
tested-proteome-background re-run of this comparison was "a direct extension
not yet done". This script does that extension, reusing:

  - the gene-set library and background definition from
    figures/enrichment_tested_background.py (gseapy.enrich, hypergeometric
    test, background = the cohort's own MS-tested gene universe, library =
    figures/figdata/Reactome_2022.gmt, cached from gseapy.get_library);
  - the SIGNATURE_FAMILIES keyword map from
    server_export/scripts/enrichment_check.py (translation / secretion /
    splicing / ptm), applied to the top-15 terms by Adjusted P-value exactly
    as enrichment_check.py applied it to the top-15 Enrichr terms.

Gene sets per cohort/level (verified to reproduce the exact counts already
reported in main.tex and server_export/scripts/enrichment_survival.tsv):
  baseline       figures/figdata/merged_incr_<cohort>.csv,
                 fdr_uni<0.05 & incremental_r2_uni>0
                 (== server_export/pinned/<cohort>/residual_results_tumoronly.csv,
                 fdr<0.05 & incremental_r2>0 -- same universe, verified identical)
  batch-corrected server_export/results/<cohort>__sitepack_operator.csv
                 (results__sitepack_operator.csv for ccrcc),
                 fdr_batchresid<0.05 & incremental_r2_batchresid>0
  strictest      same sitepack file, fdr_both<0.05 & incremental_r2_both>0

Background for every level = the full tested-protein universe for that
cohort (the 'gene' column of merged_incr_<cohort>.csv, verified identical to
the sitepack file's own gene universe for every cohort).

Writes figdata/enrichment_c2_levels.csv (top-15 terms x level x cohort, with
family classification) and prints the Table-2-shaped summary to stdout.
"""
import os
import re
import pandas as pd
import gseapy as gp

DD = os.path.join(os.path.dirname(__file__), "figdata")
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "server_export", "results")
GMT = os.path.join(DD, "Reactome_2022.gmt")

COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
SITEPACK_PREFIX = {"ccrcc": "results", "luad": "luad", "ucec": "ucec",
                    "gbm": "gbm", "pdac": "pdac"}

# Verbatim from server_export/scripts/enrichment_check.py (SIGNATURE_FAMILIES),
# the family keyword map the project already uses to score C2 survival.
SIGNATURE_FAMILIES = {
    "translation": ["translation", "ribosom", "peptide chain", "40s", "60s",
                    "nonsense-mediated", "trna", "rrna processing"],
    "secretion":   ["srp", "cotranslational", "co-translational", "secret",
                    "golgi", "vesicle", "traffick", "copi", "copii",
                    "glycosylation", "endoplasmic reticulum", "er to ",
                    "protein localization", "protein targeting", "degranulation"],
    "splicing":    ["splic", "pre-mrna", "spliceosom"],
    "ptm":         ["post-translational", "sumoylat", "ubiquitin", "proteasom"],
}
# The four families above are enrichment_check.py's, left unchanged. The paper's
# enrichment figure (Fig. 4a, make_composite2.py; family palette mrstyle.FAM) also
# defines a folding and an ECM family -- CCRCC's and PDAC's leading programmes under
# the tested background. Without them those two cohorts' leading terms would score as
# "no family", so they are added here with keywords that name only those two families.
SIGNATURE_FAMILIES["folding"] = ["folding", "chaperon", "prefoldin"]
SIGNATURE_FAMILIES["ecm"] = ["extracellular matrix", "collagen", "laminin", "elastic fibre"]
_ECM_WORD = re.compile(r"\becm\b")
ABBR = {"translation": "tr", "secretion": "sec", "splicing": "spl", "ptm": "ptm",
        "folding": "fold", "ecm": "ecm"}
FAMILY_ORDER = ["splicing", "secretion", "translation", "folding", "ecm", "ptm"]


def families_of(term):
    t = term.lower()
    fams = [f for f, keys in SIGNATURE_FAMILIES.items() if any(k in t for k in keys)]
    if "ecm" not in fams and _ECM_WORD.search(t):
        fams.append("ecm")
    return fams


def run_level(cohort, level, sig_genes, tested_genes, rows):
    sig_genes = [g for g in sig_genes if isinstance(g, str) and g]
    if len(sig_genes) < 10:
        print(f"{cohort}/{level}: too few genes ({len(sig_genes)}), skipping")
        return
    enr = gp.enrich(gene_list=sig_genes, gene_sets=GMT, background=tested_genes,
                     outdir=None)
    res = enr.results.sort_values("Adjusted P-value").reset_index(drop=True)
    n_sig = 0
    fams_seen = set()
    for rank, r in res.head(15).iterrows():
        term = r["Term"].split(" R-HSA")[0]
        fams = families_of(term)
        n_sig += int(bool(fams))
        fams_seen.update(fams)
        rows.append(dict(cohort=cohort, level=level, n_genes=len(sig_genes),
                          rank=rank + 1, term=term, qvalue=r["Adjusted P-value"],
                          overlap=r["Overlap"], signature_hit=bool(fams),
                          families="+".join(fams)))
    abbr_fams = [ABBR[f] for f in FAMILY_ORDER if f in fams_seen]
    print(f"{cohort:6s} {level:16s} n={len(sig_genes):5d}  sig_hits/15={n_sig:2d}  "
          f"families={','.join(abbr_fams) or 'NONE'}")


def main():
    rows = []
    for c in COHORTS:
        merged = pd.read_csv(os.path.join(DD, f"merged_incr_{c}.csv"))
        tested = merged["gene"].dropna().unique().tolist()
        baseline_sig = merged[(merged["fdr_uni"] < 0.05) &
                               (merged["incremental_r2_uni"] > 0)]["gene"].tolist()
        run_level(c, "baseline", baseline_sig, tested, rows)

        sp_fp = os.path.join(RESULTS, f"{SITEPACK_PREFIX[c]}__sitepack_operator.csv")
        sp = pd.read_csv(sp_fp)
        assert set(sp["gene"]) == set(tested), f"{c}: sitepack gene universe mismatch"
        batchresid_sig = sp[(sp["fdr_batchresid"] < 0.05) &
                             (sp["incremental_r2_batchresid"] > 0)]["gene"].tolist()
        run_level(c, "batch-corrected", batchresid_sig, tested, rows)

        both_sig = sp[(sp["fdr_both"] < 0.05) &
                       (sp["incremental_r2_both"] > 0)]["gene"].tolist()
        run_level(c, "strictest", both_sig, tested, rows)

    out = os.path.join(DD, "enrichment_c2_levels.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
