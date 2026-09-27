"""Single definition of the RPPA antibody rule used by every figure and table script.

The Methods exclude phospho- and other modification-specific antibodies from the TCGA
RPPA replication panels. The original pipeline (tcga_residual.py / kirc_residual.py)
removed only phospho-specific antibodies. Eight non-phospho modification-specific
antibodies (cleaved caspase-7 and -8, cleaved PARP, cleaved Notch1, acetyl-alpha-tubulin
Lys40, three histone H3/H2B modification-specific antibodies) therefore fed seven
gene-level proteins in every panel; they are removed here, whole gene by whole gene.

`load_rppa(path)` returns the per-gene results with those seven genes dropped and the
Benjamini-Hochberg FDR recomputed over the genes that remain (the per-gene permutation
p-values are unchanged, so nothing is re-permuted). `load_rppa(path, exclude=False)`
returns the panel exactly as originally analysed, for the sensitivity analysis.
"""
import numpy as np
import pandas as pd

MODIFIED_GENES = ["CASP7", "CASP8", "PARP1", "NOTCH1", "TUBA1B", "H3C1", "H2BC3"]


def bh(p):
    """Benjamini-Hochberg q-values; identical to bh_fdr in tcga_residual.py."""
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1))
        q[o[i]] = prev
    return q


def load_rppa(path, exclude=True):
    d = pd.read_csv(path).drop_duplicates("gene")
    # the stored FDR must be BH over the stored p-values, otherwise re-running BH after
    # dropping genes would not be the same procedure
    assert np.allclose(bh(d.pval.values), d.fdr.values, atol=1e-12), f"{path}: stored fdr is not BH(pval)"
    if not exclude:
        return d
    assert set(MODIFIED_GENES) <= set(d.gene), f"{path}: a modification-antibody gene is missing"
    d = d[~d.gene.isin(MODIFIED_GENES)].copy()
    d["fdr"] = bh(d.pval.values)
    return d
