#!/usr/bin/env python3
"""brca_tested_genes.py -- list the genes cptac2_residual.py will test, from LABELS and COUNTS ONLY (no morphology, no statistic).

Replicates the gene-selection steps of cptac2_residual.py exactly: protein_<short>.csv columns, rows restricted to the case list,
non-missing count n >= min_n (30), and presence among the RNA gene_name symbols of the STAR files in rna_case_file.tsv
(the 'g in rna.columns' step). Prints / writes symbols and n only. Expected for the published DX arm: 10,031 genes.
The RNA presence step reads the gene_name column only (no expression values are used).

  python brca_tested_genes.py --root /public/home/fjhui/ZW/cptac2_brca --case-list /public/home/fjhui/ZW/cptac2_brca/dx_case_list.tsv \
         --out /public/home/fjhui/ZW/cptac2_brca/brca_tested_genes.tsv
Bring brca_tested_genes.tsv back; feed its 'gene' column to brca_selection_set_inventory.py --brca-genes (one symbol per line:
the script also writes <out>.txt).
"""
import argparse
import os
import sys

import pandas as pd


def tcga_case(s):
    return "-".join(str(s).split("-")[:3])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/public/home/fjhui/ZW/cptac2_brca")
    ap.add_argument("--case-list", default="")
    ap.add_argument("--short", default="brca")
    ap.add_argument("--min-n", type=int, default=30)
    ap.add_argument("--out", required=True)
    ap.add_argument("--expect", type=int, default=10031, help="published tested-gene count (paper); a mismatch is printed, not fatal")
    a = ap.parse_args()
    prot = pd.read_csv(f"{a.root}/protein_{a.short}.csv", index_col=0)
    prot.index = [tcga_case(c) for c in prot.index]; prot = prot[~prot.index.duplicated()]
    xw = pd.read_csv(f"{a.root}/rna_case_file.tsv", sep="\t")
    rna_cases = {tcga_case(c) for c in xw["case"]}
    cases = set(prot.index) & rna_cases
    if a.case_list:
        want = {l.strip() for l in open(a.case_list) if l.strip().startswith("TCGA-")}
        cases &= want
    prot = prot.loc[sorted(cases)]
    names = set()
    for fn in xw["filename"]:
        fp = f"{a.root}/rna/{fn}"
        if os.path.exists(fp):
            d = pd.read_csv(fp, sep="\t", comment="#", usecols=lambda c: c.strip() == "gene_name")
            g = d.iloc[:, 0].dropna().astype(str)
            names |= set(g[~g.str.startswith("N_")])
    n = prot.notna().sum(axis=0)
    keep = [g for g in prot.columns if g in names and n[g] >= a.min_n]
    out = pd.DataFrame({"gene": keep, "n": [int(n[g]) for g in keep]})
    out.to_csv(a.out, sep="\t", index=False)
    with open(a.out.rsplit(".", 1)[0] + ".txt", "w") as fh:
        fh.write("\n".join(keep) + "\n")
    print(f"cases (protein & rna & case-list): {len(cases)}; protein columns {prot.shape[1]}; RNA symbols {len(names)}; "
          f"tested (n>={a.min_n}, in RNA): {len(keep)}  [published: {a.expect}]"
          + ("" if len(keep) == a.expect else "   <-- DIFFERS: investigate before using the gene list"))
    if prot.columns.duplicated().any():
        print("WARNING: duplicated protein column labels", file=sys.stderr)


if __name__ == "__main__":
    main()
