#!/usr/bin/env python3
r"""Recompute Table 1 and Fig_controls panels b--d from the CURRENT analysis outputs.

WHY
The three count columns of Table 1 and panels c and d of Fig_controls are hardcoded
literals (make_fig_controls.py:104,131-133,181; make_figs.py:16-19). The baseline
row reproduces exactly from pinned/*/residual_results_tumoronly.csv. The
batch-corrected and strictest rows reproduced from nothing until the pre-fix
backups were recovered from the server, and they are stale: every sitepack output
exists twice, a *.pre_otherfix.bak written 2026-09-01 01:44-16:34 and a current
file written 17:45-18:27 the same day. residual_analysis_sitepack.py documents what
changed in between -- an earlier version pooled every batch outside the top
MAX_LEVELS into one "_other" level holding 45% of PDAC's patients across 54 distinct
operators, corrupting residualisation, the within-batch permutation and GroupKFold
at once, "which is why PDAC's strictest estimator collapsed to zero". The manuscript
reports exactly that zero. It quotes the pre-fix run.

THE CRITERION IS NOT GUESSED
It was recovered from the pre-fix backups by exhaustive search over both batch axes,
every incremental_r2 column, every fdr/pval column, thresholds 0.05/0.01/0.10, with
and without a positivity filter, and with and without nesting inside the baseline
significant set. Exactly one combination reproduces all five cohorts at once for
each row, and it is the same family the baseline uses:

    batch-corrected  fdr_batchresid < 0.05 and incremental_r2_batchresid > 0
    strictest        fdr_both       < 0.05 and incremental_r2_both       > 0
    axis             operator, and the sets are NOT nested inside the baseline

So the definition is held fixed and only the data changes. The one thing that would
have confounded the update -- switching to a nested "retention" definition, which
was considered because a slope chart reads as nesting -- is deliberately NOT done.
Note that the counts are therefore not monotone across stages in every cohort: BH
runs separately per estimator, so a gene can fail one stage and pass the next. That
is a property of the manuscript's own definition and is stated in the caption
rather than engineered away.

An independent check that the re-run is the correct one: under the pre-fix numbers
GBM has 847 batch-corrected genes against a baseline of 702, i.e. more significant
after correction than before. Under the fixed run it is 680 < 702. The impossible
number was itself a symptom of the "_other" bucket.

ENRICHMENT (panel d)
Recomputed, because PDAC's strictest cell showed 0 signature terms only because it
had 0 genes, and it now has 342. This calls the project's own enrichment_check.py
functions directly rather than reimplementing them -- post_list, enrich, families_of
and LIBS are imported, not copied, after a hand-transcribed copy was found to have
corrupted two of the four keyword lists. That path uses online Enrichr with the
standard whole-genome background and the Reactome_Pathways_2024 -> 2022 -> 2016
fallback chain, which is what panel d has always used and is consistent with
Limitations (ii); it is NOT the tested-proteome background of \S\ref{sec:repro},
which governs the enrichment table instead. This sends the gene symbol lists to
Enrichr, as the project's pipeline already does.

Writes figures/figdata/batch_levels.csv, which the figure scripts then read instead
of carrying literals.
"""
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, "server_export/scripts")
import enrichment_check as EC          # noqa: E402  post_list / enrich / families_of / LIBS

EX = Path("server_export")
DD = Path("figures/figdata")
COH = ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]
SITEPACK = {"CCRCC": "results__sitepack_operator.csv", "LUAD": "luad__sitepack_operator.csv",
            "UCEC": "ucec__sitepack_operator.csv", "GBM": "gbm__sitepack_operator.csv",
            "PDAC": "pdac__sitepack_operator.csv"}
# what the manuscript currently prints, for the change report
OLD = {"CCRCC": (2191, 1010, 913), "LUAD": (2310, 893, 332), "UCEC": (2566, 481, 573),
       "GBM": (702, 847, 223), "PDAC": (1572, 1055, 0)}
OLD_SIGTERMS = {"CCRCC": (7, 7, 5), "LUAD": (12, 10, 11), "UCEC": (11, 13, 11),
                "GBM": (6, 5, 4), "PDAC": (3, 2, 0)}


def sig_terms(cohort, label, genes):
    """Top-15 Reactome terms via the project's own Enrichr path; count signature hits."""
    genes = sorted(str(g) for g in genes if isinstance(g, str) and g)
    if len(genes) < 10:
        print(f"  {cohort:6s} {label:10s} n={len(genes):5d}  too few to enrich -> 0")
        return 0, "", 0
    uid = EC.post_list(genes[:3000])
    if uid is None:
        raise SystemExit(f"Enrichr unreachable for {cohort}/{label}")
    res = []
    for lib in EC.LIBS:
        res = EC.enrich(uid, lib)
        if res:
            break
    else:
        raise SystemExit(f"no Reactome library returned results for {cohort}/{label}")
    top = res[:15]
    hits = sum(1 for t in top if EC.families_of(t[1]))
    print(f"  {cohort:6s} {label:10s} n={len(genes):5d}  lib={lib:24s} "
          f"terms={len(top):2d}  signature={hits}")
    return hits, lib, len(top)


def main():
    rows = []
    for c in COH:
        base = pd.read_csv(EX / "pinned" / c.lower() / "residual_results_tumoronly.csv",
                           low_memory=False)
        sp = pd.read_csv(EX / "results" / SITEPACK[c], low_memory=False)

        g0 = base[(base.fdr < 0.05) & (base.incremental_r2 > 0)].gene.astype(str).tolist()
        g1 = sp[(sp.fdr_batchresid < 0.05)
                & (sp.incremental_r2_batchresid > 0)].gene.astype(str).tolist()
        g2 = sp[(sp.fdr_both < 0.05)
                & (sp.incremental_r2_both > 0)].gene.astype(str).tolist()

        ob, obc, ost = OLD[c]
        assert len(g0) == ob, f"{c} baseline {len(g0)} != published {ob}"

        st, lib, _ = {}, "", 0
        for label, g in [("baseline", g0), ("batchcorr", g1), ("strictest", g2)]:
            st[label], lib, _ = sig_terms(c, label, g)
            time.sleep(1.0)
        rows.append(dict(cohort=c, n_tested=len(base), baseline=len(g0),
                         batchcorr=len(g1), strictest=len(g2),
                         sigterms_baseline=st["baseline"],
                         sigterms_batchcorr=st["batchcorr"],
                         sigterms_strictest=st["strictest"], library=lib))

    df = pd.DataFrame(rows)

    # scanner percentage over ALL slides, a blank ScanScope ID counted as unknown
    sc = pd.read_csv(EX / "wsi_scanner_map.tsv", sep="\t", dtype=str, keep_default_na=False)
    sc = sc[sc.error == ""]
    pct, nsc = [], []
    for c in COH:
        g = sc[sc.cohort.str.upper() == c]
        named = g.scanscope_id.value_counts()
        named = named[named.index != ""]
        pct.append(round(100 * named.iloc[0] / len(g), 1))
        nsc.append(int(len(named)))
    df["scan_pct"] = pct
    df["n_scanners"] = nsc
    df["pct_of_tested"] = (100 * df.baseline / df.n_tested).round(1)

    DD.mkdir(parents=True, exist_ok=True)
    df.to_csv(DD / "batch_levels.csv", index=False)
    pd.set_option("display.width", 220)
    print("\n" + df.to_string(index=False))

    print("\nchange against the manuscript's current literals (definition held fixed):")
    print(f"  {'cohort':6s} {'batch-corrected':>22s} {'strictest':>18s} "
          f"{'sig.terms b/c/s':>26s}")
    for _, r in df.iterrows():
        ob, obc, ost = OLD[r.cohort]
        s0, s1, s2 = OLD_SIGTERMS[r.cohort]
        print(f"  {r.cohort:6s} {obc:8d} -> {int(r.batchcorr):-8d} "
              f"{ost:8d} -> {int(r.strictest):-6d}   "
              f"({s0},{s1},{s2}) -> ({int(r.sigterms_baseline)},"
              f"{int(r.sigterms_batchcorr)},{int(r.sigterms_strictest)})")
    print(f"\nwrote {DD/'batch_levels.csv'}")


if __name__ == "__main__":
    main()
