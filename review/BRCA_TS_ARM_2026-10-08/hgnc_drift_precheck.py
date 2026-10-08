#!/usr/bin/env python
"""Symbol-drift exposure of the candidate selection sets, from PUBLIC METADATA ONLY (HGNC REST API, symbol-change dates).

Question: if the CPTAC-2 iTRAQ table (PDC000173) names genes with an older symbol than the CPTAC-3 TMT tables, which candidate-set
genes could silently fail an exact-symbol match? One HGNC query per cut-off year returns every approved symbol whose
date_symbol_changed is on or after that date (no per-gene calls, no data download). A discovery gene whose symbol was last
changed AFTER the BRCA table's annotation date is the at-risk class; a symbol changed before it is safe.
The annotation date of PDC000173 is not known locally, so exposure is reported for several cut-offs (the user may pin one from the
PDC study metadata; this script does not choose it).

Needs network on first run (cache: setlevel/hgnc_changed_<YYYY>.json). Reads setlevel/universe_common.tsv and candidate_sets/*.tsv
written by brca_selection_set_inventory.py. NO BRCA DATA.
"""
import json
import os
import sys
import urllib.parse
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "setlevel")
CUTOFFS = [2012, 2014, 2016, 2018, 2020]


def changed_since(year):
    cache = os.path.join(OUT, f"hgnc_changed_{year}.json")
    if os.path.exists(cache):
        return set(json.load(open(cache)))
    q = f"date_symbol_changed:[{year}-01-01T00:00:00Z TO *]"
    url = "https://rest.genenames.org/search/" + urllib.parse.quote(q, safe=":[]*") + "?rows=20000"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        d = json.load(r)["response"]
    syms = sorted(x["symbol"] for x in d["docs"])
    if d["numFound"] != len(syms):
        sys.exit(f"FATAL: HGNC returned {len(syms)} of {d['numFound']} for {year}")
    json.dump(syms, open(cache, "w"))
    return set(syms)


def main():
    uni = pd.read_csv(os.path.join(OUT, "universe_common.tsv"), sep="\t")
    ch = {y: changed_since(y) for y in CUTOFFS}
    rows = []
    for sid in ["pub_k1", "pub_k2", "pub_k3", "pub_k4", "op_k1", "op_k2", "plex_k1", "plex_k2", "phi_k2", "phi_k3", "phi_k4"]:
        d = pd.read_csv(os.path.join(OUT, "candidate_sets", f"{sid}.tsv"), sep="\t")
        sel = set(d.loc[d.selected == 1, "gene"])
        r = dict(set_id=sid, n_selected=len(sel))
        for y in CUTOFFS:
            r[f"n_changed_since_{y}"] = len(sel & ch[y])
            r[f"frac_{y}"] = round(len(sel & ch[y]) / len(sel), 4)
        rows.append(r)
    r = dict(set_id="universe_common(8309)", n_selected=len(uni))
    for y in CUTOFFS:
        r[f"n_changed_since_{y}"] = int(uni.gene.isin(ch[y]).sum()); r[f"frac_{y}"] = round(r[f"n_changed_since_{y}"] / len(uni), 4)
    rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "hgnc_drift_summary.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(out.to_string(index=False))
    # the 'changed 2014-2020' window = renamed between a plausible CPTAC-2 annotation date and the CPTAC-3 tables
    win = ch[2014] - ch[2020]
    print(f"symbols changed in [2014, 2020): {len(win)} HGNC-wide; in the common universe: {int(uni.gene.isin(win).sum())}")
    pd.DataFrame({"gene": sorted(set(uni.gene) & ch[2014])}).to_csv(os.path.join(OUT, "at_risk_changed_since_2014.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
