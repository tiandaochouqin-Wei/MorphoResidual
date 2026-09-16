#!/usr/bin/env python3
r"""
export_phospho_hosts.py -- MorphoResidual T1-4 (NC_ROADMAP.md).

WHY
The manuscript's mTORC1 section rests on phosphosite intensities that are NOT
normalised to the abundance of the protein carrying the site. Because the
translation residual is itself dominated by ribosomal proteins, a phospho-S6
signal could be reporting RPS6 protein abundance rather than S6K1 activity, and
\S\ref{sec:mtorres} currently says so and calls protein-normalised ratios "the
required direct test", which we have not done. This exports the one thing needed
to do it locally: the TOTAL protein abundance of the host genes, from the same
protein matrix the main pipeline uses.

Using the pipeline's own matrix rather than a fresh LinkedOmics download matters:
aliquot-to-case mapping, the tumour-only filter and the normalisation are already
settled and identical to everything else in the paper, so a ratio computed from
it is directly comparable to the residual scores. A separate download would
reintroduce all three as alignment risks for no benefit.

WHAT IT EXPORTS (per cohort, patients x genes, tiny -- a few hundred kB)
  the four mTORC1 phospho-site host proteins:
      RPS6, EIF4EBP1, EIF4B, RPS6KB1
  plus every cytosolic ribosomal protein (RPS*/RPL*/RPLP0-2/RPSA), which is the
  background the non-mTORC1 control needs: those sites must be normalised to
  THEIR own host proteins too, or the control is asymmetric and the comparison
  is meaningless.

USAGE (per cohort; seconds each, no GPU, no permutations)
  cd /public/home/fjhui/ZW/scripts
  for c in ccrcc luad ucec gbm pdac; do
    source morpho_env.sh $c
    export MORPHO_HOST_OUT=phospho_hosts_$c.csv
    python -u export_phospho_hosts.py
  done

  Then transfer the five phospho_hosts_*.csv back to figures/figdata/.

  NOTE: this script deliberately does NOT read MORPHO_ROOT or any other variable
  that morpho_env.sh sets for its own purposes -- the export_package.py run hit
  exactly that collision. The only variable it reads is MORPHO_HOST_OUT, which
  nothing else uses.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

import residual_analysis as RA  # loaders only; env is read at import, as usual

MTOR_HOSTS = ["RPS6", "EIF4EBP1", "EIF4B", "RPS6KB1"]
RIBO_PAT = re.compile(r"^RP[SL]\d+[A-Z]?$|^RPLP[012]$|^RPSA$")   # matches mtor_mechanism.py

EXPECTED_N = {"ccrcc": 103, "luad": 105, "ucec": 100, "gbm": 99, "pdac": 137}


def infer_cohort():
    for src in (os.environ.get("MORPHO_HOST_OUT", ""),
                os.environ.get("MORPHO_SIG_CSV", ""),
                os.environ.get("MORPHO_COHORT", "")):
        m = re.search(r"(ccrcc|luad|ucec|gbm|pdac)", src.lower())
        if m:
            return m.group(1)
    return None


def main():
    out = os.environ.get("MORPHO_HOST_OUT", "phospho_hosts.csv")
    cohort = infer_cohort()
    print(f"[hosts] output: {out}")
    print(f"[hosts] cohort inferred as: '{cohort}'")

    protein = RA.load_protein_matrix()
    print(f"[hosts] protein matrix: {protein.shape[0]} patients x {protein.shape[1]} genes")

    # sanity: the loaders must be pointing at the cohort we think they are
    if cohort in EXPECTED_N and abs(protein.shape[0] - EXPECTED_N[cohort]) > 6:
        sys.exit(f"[hosts] STOP: {protein.shape[0]} patients in the protein matrix, expected "
                 f"~{EXPECTED_N[cohort]} for '{cohort}'. morpho_env.sh is probably set to a "
                 f"different cohort than MORPHO_HOST_OUT names -- do not trust this export.")

    genes = list(protein.columns)
    ribo = sorted(g for g in genes if RIBO_PAT.match(str(g)))
    hosts = [g for g in MTOR_HOSTS if g in genes]
    missing_hosts = [g for g in MTOR_HOSTS if g not in genes]

    keep = sorted(set(hosts) | set(ribo))
    if not keep:
        sys.exit("[hosts] STOP: none of the requested genes are in this protein matrix.")

    sub = protein[keep].copy()
    sub.index.name = "case"
    sub.to_csv(out)

    print(f"[hosts] mTORC1 host proteins found: {hosts}")
    if missing_hosts:
        print(f"[hosts] mTORC1 host proteins MISSING from this cohort's matrix: {missing_hosts} "
              f"(the corresponding arm cannot be protein-normalised here; report as such)")
    print(f"[hosts] ribosomal background proteins: {len(ribo)}")
    print(f"[hosts] wrote {out}: {sub.shape[0]} patients x {sub.shape[1]} genes")

    nn = sub.notna().sum()
    print(f"[hosts] non-missing per gene: median {int(nn.median())}, "
          f"min {int(nn.min())}, max {int(nn.max())} of {sub.shape[0]} patients")
    for g in hosts:
        v = sub[g]
        print(f"    {g:10s} n={int(v.notna().sum()):3d}  "
              f"mean log-ratio {np.nanmean(v.values):+.4f}  sd {np.nanstd(v.values):.4f}")


if __name__ == "__main__":
    main()
