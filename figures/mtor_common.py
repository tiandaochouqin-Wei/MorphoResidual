#!/usr/bin/env python3
r"""Shared phosphoproteomics helpers for the mTORC1 analyses.

Extracted verbatim from mtor_mechanism.py so that mtor_protein_normalised.py
(T1-4: protein-normalised phosphosite ratios) uses byte-identical parsing rather
than a second, drifting copy. mtor_mechanism.py now imports from here; its
outputs were diffed before and after the extraction and are unchanged.

The site-id parsing is deliberately liberal because the five CPTAC cohorts do not
share a layout: LUAD/UCEC/GBM carry site-level rows keyed like 'RPS6_S235s',
PDAC keys on Index | Gene | Peptide, and CCRCC is gene-level only (which is why
its mTORC1 result is reported as gene-level in the manuscript).
"""
import re

import numpy as np
import pandas as pd
from scipy import stats as st

# the canonical mTORC1 phospho-outputs, and the ERK specificity control
MTOR = {"RPS6": ["S235", "S236", "S240", "S244"], "EIF4EBP1": ["T37", "T46", "S65", "T70"],
        "EIF4B": ["S422"], "RPS6KB1": ["T389", "T390"]}   # T390 = T389 under alternative isoform numbering (GBM)
ERK = {"MAPK1": ["T185", "Y187"], "MAPK3": ["T202", "Y204"]}
# ribosomal-protein genes, for the non-mTORC1 control: is phospho-S6 special
# among ribosomal phosphosites, or typical of the class?
RIBO_PAT = re.compile(r"^RP[SL]\d+[A-Z]?$|^RPLP[012]$|^RPSA$")


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


def partial_r(x, y, c):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
    x, y, c = x[m], y[m], c[m]
    rx = x - np.polyval(np.polyfit(c, x, 1), c); ry = y - np.polyval(np.polyfit(c, y, 1), c)
    return st.pearsonr(rx, ry)[0], int(m.sum())
