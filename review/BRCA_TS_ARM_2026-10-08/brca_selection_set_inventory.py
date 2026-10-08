#!/usr/bin/env python
"""Inventory of candidate discovery 'selection sets' for a set-level BRCA read-out. NO BRCA DATA, NO MORPHOLOGY-PROTEIN STATISTIC.

Only gene lists and counts from the already-published, pinned CPTAC-3 discovery tables (fdr < 0.05 & increment > 0 -- the
paper's own definition of 'significant') are used here. The only inputs are:
  server_export/pinned/<c>/residual_results_tumoronly.csv      published discovery (UNI), c in ccrcc luad ucec gbm pdac
  server_export/results/{results,<c>}__sitepack_operator.csv    operator batch axis: fdr_batchresid / fdr_both
  server_export/pinned/posthoc/<c>_plex.csv                     TMT-plex batch axis
  figures/figdata/merged_incr_<c>.csv                           Phikon columns (fdr_phi) -- the BRCA arm's encoder
  server_export/pinned/ucec_c1/c1_ucec_gene_level.csv           PDC000439 quantDataMatrix symbols (PDC-symbol PROXY only)
Optional: --brca-genes FILE (one symbol per line = the genes the BRCA arm tests; produced by brca_tested_genes.py from
labels / non-missing counts only) -> overlap columns. Without it the overlap columns say PENDING.

Writes into setlevel/ : candidate_sets/<ID>.tsv (gene, selected, score, k over the set's universe), candidate_sets_inventory.tsv/.json,
symbol_qc.json, universe_common.tsv.
Definitions (all fixed here, none tuned):
  universe  = genes tested in ALL five discovery cohorts (8,309; the universe of the paper's 170-gene core)
  k_fam(g)  = number of the five cohorts in which g is selected under family fam
  families  : pub   published discovery (UNI)                         op   operator batch-corrected (Table 1 'Batch-corr.')
              opb   operator 'strictest' (fdr_both)                  plex TMT-plex batch-corrected
              phi   published discovery, Phikon encoder
  sets      : <fam>_k{1,2,3,4,5} = k_fam >= m  (k1 = union of the five cohorts' sets; pub_k4 = the 170-gene core)
              single-organ sets <fam>_<cohort> over that cohort's own tested genes
  score     : mean over cohorts of the family's increment (mean of the cohorts that tested g) -- used only for the Spearman read-out
"""
import argparse
import hashlib
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
SE = os.path.join(PAPER, "server_export")
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
OUT = os.path.join(HERE, "setlevel")
SETDIR = os.path.join(OUT, "candidate_sets")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path):
    return os.path.relpath(path, PAPER).replace("\\", "/")


def table(fam, c):
    """-> (path, DataFrame indexed by gene with columns fdr, inc)."""
    if fam == "pub":
        p = os.path.join(SE, "pinned", c, "residual_results_tumoronly.csv"); f, i = "fdr", "incremental_r2"
    elif fam in ("op", "opb"):
        p = os.path.join(SE, "results", ("results" if c == "ccrcc" else c) + "__sitepack_operator.csv")
        f, i = ("fdr_batchresid", "incremental_r2_batchresid") if fam == "op" else ("fdr_both", "incremental_r2_both")
    elif fam == "plex":
        p = os.path.join(SE, "pinned", "posthoc", f"{c}_plex.csv"); f, i = "fdr_batchresid", "incremental_r2_batchresid"
    elif fam == "phi":
        p = os.path.join(PAPER, "figures", "figdata", f"merged_incr_{c}.csv"); f, i = "fdr_phi", "incremental_r2_phi"
    else:
        raise ValueError(fam)
    d = pd.read_csv(p)
    if "cohort" in d.columns:
        d = d[d["cohort"].str.lower() == c]
    if "batch_axis" in d.columns:
        d = d[d["batch_axis"].isin(["operator", "plex"])]
    if d["gene"].duplicated().any():
        sys.exit(f"FATAL: duplicated genes in {p}")
    d = d.set_index("gene")[[f, i]].rename(columns={f: "fdr", i: "inc"})
    return p, d.dropna()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--brca-genes", default="", help="file with the BRCA-tested gene symbols (one per line); optional")
    a = ap.parse_args()
    os.makedirs(SETDIR, exist_ok=True)
    fams = ["pub", "op", "opb", "plex", "phi"]
    T, src = {}, {}
    for fam in fams:
        for c in COH:
            p, d = table(fam, c)
            T[fam, c] = d; src[fam, c] = dict(path=rel(p), sha256=sha256(p), n_rows=int(len(d)))
    # lineage check: the 'uni' columns of merged_incr == the pinned published tables
    chk = {}
    for c in COH:
        m = pd.read_csv(os.path.join(PAPER, "figures", "figdata", f"merged_incr_{c}.csv")).set_index("gene")
        d = T["pub", c]
        same_genes = set(m.index) == set(d.index)
        sel_m = set(m.index[(m.fdr_uni < 0.05) & (m.incremental_r2_uni > 0)])
        sel_d = set(d.index[(d.fdr < 0.05) & (d.inc > 0)])
        chk[c] = dict(same_gene_set=bool(same_genes), n_sel_pinned=len(sel_d), n_sel_merged_uni=len(sel_m),
                      same_selection=bool(sel_m == sel_d))
    print("lineage pinned-vs-merged_incr(uni):", json.dumps(chk))

    tested = {c: set(T["pub", c].index) for c in COH}
    common = sorted(set.intersection(*tested.values()))
    union_tested = set.union(*tested.values())
    print(f"genes tested per cohort: {({c: len(tested[c]) for c in COH})}; common {len(common)}; union {len(union_tested)}")

    sel, K, SC = {}, {}, {}
    for fam in fams:
        for c in COH:
            d = T[fam, c]
            sel[fam, c] = set(d.index[(d.fdr < 0.05) & (d.inc > 0)])
        K[fam] = pd.Series(0, index=common, dtype=int)
        incs = pd.DataFrame({c: T[fam, c]["inc"].reindex(common) for c in COH})
        for c in COH:
            K[fam] += pd.Series(common, index=common).isin(sel[fam, c]).astype(int)
        SC[fam] = incs.mean(axis=1, skipna=True)
    # universe file
    uni = pd.DataFrame({"gene": common})
    for fam in fams:
        uni[f"k_{fam}"] = K[fam].values; uni[f"score_{fam}"] = SC[fam].values
    uni.to_csv(os.path.join(OUT, "universe_common.tsv"), sep="\t", index=False)

    # symbol QC (no BRCA input): non-standard-looking symbols in the discovery universe
    allsyms = sorted(union_tested)
    pats = {
        "date_like": re.compile(r"^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$|^(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-\d{1,2}$", re.I),
        "has_dot_suffix": re.compile(r"\.\d+$"),
        "has_space_or_semicolon": re.compile(r"[\s;,]"),
        "lowercase_start": re.compile(r"^[a-z]"),
        "orf_style_old_name": re.compile(r"^C\d+orf\d+$|^C[0-9XY]+orf\d+$"),
        "lncrna_or_clone_style": re.compile(r"^(AC|AL|AP|LINC|RP\d+-)"),
    }
    qc = {k: sorted(s for s in allsyms if r.search(s)) for k, r in pats.items()}
    symbol_qc = {k: dict(n=len(v), examples=v[:15]) for k, v in qc.items()}

    # PDC-symbol proxy: UCEC-C1 (PDC000439 quantDataMatrix, same cptac2_replicate.py row[0] symbol path)
    c1 = pd.read_csv(os.path.join(SE, "pinned", "ucec_c1", "c1_ucec_gene_level.csv"))
    c1_sym = set(c1["gene"])
    ucec_t = tested["ucec"]
    lost = sorted(ucec_t - c1_sym)
    symbol_qc["pdc_proxy_ucec_c1"] = dict(
        discovery_ucec_tested=len(ucec_t), in_c1_pdc_matrix_usable=len(ucec_t & c1_sym),
        frac=round(len(ucec_t & c1_sym) / len(ucec_t), 4), n_lost=len(lost), lost_examples=lost[:40],
        note="lost = tested in discovery UCEC but absent from the usable (n>=30) PDC000439 matrix; mixes symbol drift with "
             "missing-protein/missing-RNA/low-n; an upper bound on symbol drift for a same-API PDC study")
    common_in_c1 = sum(g in c1_sym for g in common)
    symbol_qc["pdc_proxy_common_universe_in_c1"] = dict(n=common_in_c1, of=len(common), frac=round(common_in_c1 / len(common), 4))
    json.dump(symbol_qc, open(os.path.join(OUT, "symbol_qc.json"), "w"), indent=1)

    brca = None
    if a.brca_genes:
        brca = set(l.strip() for l in open(a.brca_genes, encoding="utf-8") if l.strip())
        print(f"BRCA gene list: {len(brca)} symbols <- {a.brca_genes}")

    rows = []

    def add(sid, desc, fam, genes_universe, selected, scorev, kv=None):
        selected = set(selected)
        df = pd.DataFrame({"gene": genes_universe})
        df["selected"] = df["gene"].isin(selected).astype(int)
        df["score"] = np.asarray(scorev, float)
        df["k"] = (np.asarray(kv, int) if kv is not None else -1)
        df.to_csv(os.path.join(SETDIR, f"{sid}.tsv"), sep="\t", index=False)
        u = set(genes_universe)
        r = dict(set_id=sid, family=fam, definition=desc, n_selected=len(selected), n_universe=len(u),
                 frac_selected=round(len(selected) / len(u), 4),
                 n_selected_in_pdc_proxy=int(sum(g in c1_sym for g in selected)),
                 file=rel(os.path.join(SETDIR, f"{sid}.tsv")))
        if brca is not None:
            r.update(n_selected_in_brca=len(selected & brca), n_universe_in_brca=len(u & brca),
                     frac_selected_in_brca=round(len(selected & brca) / max(len(selected), 1), 4))
        else:
            r.update(n_selected_in_brca="PENDING", n_universe_in_brca="PENDING", frac_selected_in_brca="PENDING")
        rows.append(r)

    FAMNAME = {"pub": "published discovery (UNI)", "op": "operator batch-corrected", "opb": "operator strictest (fdr_both)",
               "plex": "TMT-plex batch-corrected", "phi": "published discovery, Phikon encoder"}
    for fam in fams:
        for m in range(1, 6):
            nm = {1: "union of 5 cohorts' sets", 5: "all five"}.get(m, f">= {m} of 5 cohorts")
            sid = f"{fam}_k{m}"
            extra = " [= the paper's 170-gene core]" if (fam == "pub" and m == 4) else ""
            add(sid, f"{FAMNAME[fam]}: significant in {nm}{extra}", fam, common, K[fam].index[K[fam] >= m], SC[fam].values, K[fam].values)
    for fam in ("pub", "op", "plex"):
        for c in COH:
            g_u = sorted(tested[c])
            add(f"{fam}_{c}", f"{FAMNAME[fam]}: {c.upper()} only (universe = genes tested in {c.upper()})", fam, g_u,
                sel[fam, c], T[fam, c]["inc"].reindex(g_u).values)
    inv = pd.DataFrame(rows)
    inv.to_csv(os.path.join(OUT, "candidate_sets_inventory.tsv"), sep="\t", index=False)

    # nestedness / overlap among cross-organ published sets (set membership only)
    K_pub = K["pub"]
    ov = {}
    for fam in ("op", "opb", "plex", "phi"):
        for m in (1, 2, 3, 4):
            s_f = set(K[fam].index[K[fam] >= m]); s_p = set(K_pub.index[K_pub >= m])
            ov[f"{fam}_k{m}_vs_pub_k{m}"] = dict(n_fam=len(s_f), n_pub=len(s_p), shared=len(s_f & s_p),
                                                 jaccard=round(len(s_f & s_p) / max(len(s_f | s_p), 1), 4))
    # genes selected somewhere but outside the common universe (what the common-universe rule drops)
    outside = {c: len(sel["pub", c] - set(common)) for c in COH}
    meta = dict(lineage=chk, n_tested_per_cohort={c: len(tested[c]) for c in COH}, n_common=len(common),
                n_union_tested=len(union_tested), n_selected_published_outside_common=outside,
                n_selected_per_cohort={fam: {c: len(sel[fam, c]) for c in COH} for fam in fams},
                overlap_vs_published_cross_organ=ov, sources={f"{k[0]}:{k[1]}": v for k, v in src.items()},
                brca_gene_file=a.brca_genes or None)
    json.dump(meta, open(os.path.join(OUT, "candidate_sets_inventory.json"), "w"), indent=1)
    pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 70)
    print(inv.drop(columns=["file"]).to_string(index=False))
    print(json.dumps(meta["n_selected_per_cohort"]))
    print("selected outside common universe (published):", outside)
    print(json.dumps(ov, indent=0)[:2500])
    print("symbol QC:", json.dumps({k: (v if k.startswith("pdc") else v["n"]) for k, v in symbol_qc.items()})[:1500])


if __name__ == "__main__":
    main()
