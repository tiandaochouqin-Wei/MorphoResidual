#!/usr/bin/env python3
"""
Does the C2 enrichment survive the acquisition-batch correction?

The pass criterion fixed in PAPER_PLAN is two-part: after the robustness pack the
significant count may fall a lot, but (i) C1 must stay > 0 and (ii) the C2
translation / secretion enrichment must still be there. This script tests (ii).

For each cohort it enriches three gene sets against Reactome via Enrichr:
  baseline   significant in residual_results_tumoronly.csv        (uncorrected)
  batchresid significant in sitepack_<axis>.csv, batch-residualised
  both       significant in sitepack_<axis>.csv, residualised + GroupKFold

and reports whether the translation / SRP / secretion / splicing families that
converged across ccRCC, LUAD and GBM are still on top.

Usage: enrichment_check.py [axis] [cohort ...]
Output: scripts/enrichment_survival.tsv  + a readable summary on stdout
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

ROOT = "/public/home/fjhui/ZW"
SCRIPTS = f"{ROOT}/scripts"
ADDLIST = "https://maayanlab.cloud/Enrichr/addList"
ENRICH = "https://maayanlab.cloud/Enrichr/enrich"
LIBS = ["Reactome_Pathways_2024", "Reactome_2022", "Reactome_2016"]

# The families the project reported converging across ccRCC / LUAD / GBM, grouped
# so the report can say WHICH one survived rather than just counting hits.
#
# The secretion family deliberately includes the vesicle/Golgi/glycosylation
# vocabulary: a first pass keyed only on the literal word "secretion" scored UCEC
# at 3/15 while its actual top terms were Membrane Trafficking, Vesicle-mediated
# Transport, Intra-Golgi and Retrograde Golgi-to-ER Traffic, ER to Golgi
# Anterograde Transport and COPI-mediated Anterograde Transport -- i.e. the
# secretory pathway under every name except that one. Keying on the narrow word
# understates C2 survival.
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


def families_of(term):
    t = term.lower()
    return [f for f, keys in SIGNATURE_FAMILIES.items() if any(k in t for k in keys)]


def results_dir(cohort):
    return ROOT + ("/results" if cohort == "ccrcc" else f"/{cohort}/results")


def post_list(genes):
    boundary = "----MorphoResidualBoundary"
    parts = []
    for name, val in (("list", "\n".join(genes)), ("description", "morphoresidual")):
        parts.append(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\"\r\n\r\n{val}\r\n")
    body = ("".join(parts) + f"--{boundary}--\r\n").encode()
    req = urllib.request.Request(
        ADDLIST, data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)["userListId"]
        except Exception as e:
            print(f"    addList attempt {attempt}/4 failed: {e}", flush=True)
            time.sleep(3 * attempt)
    return None


def enrich(user_id, lib):
    url = f"{ENRICH}?{urllib.parse.urlencode({'userListId': user_id, 'backgroundType': lib})}"
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r).get(lib, [])
        except Exception as e:
            print(f"    enrich attempt {attempt}/4 failed: {e}", flush=True)
            time.sleep(3 * attempt)
    return []


def run_set(label, cohort, genes, rows):
    genes = [g for g in genes if isinstance(g, str) and g]
    print(f"\n--- {cohort} / {label}: {len(genes)} genes", flush=True)
    if len(genes) < 10:
        print("    too few genes to enrich")
        rows.append(dict(cohort=cohort, gene_set=label, n_genes=len(genes),
                         library="", rank=-1, term="", pvalue=None, qvalue=None,
                         signature_hit=None))
        return
    uid = post_list(genes[:3000])
    if uid is None:
        print("    Enrichr unreachable")
        return
    for lib in LIBS:
        res = enrich(uid, lib)
        if res:
            break
    else:
        print("    no Reactome library returned results")
        return
    print(f"    library {lib}, {len(res)} terms")
    seen_families = set()
    sig_hits = 0
    for rank, t in enumerate(res[:15], 1):
        term, pval, qval = t[1], t[2], t[6]
        fams = families_of(term)
        seen_families.update(fams)
        sig_hits += int(bool(fams))
        mark = f"  <== {'+'.join(fams)}" if fams else ""
        print(f"    {rank:>2}. q={qval:.3g}  {term[:66]}{mark}")
        rows.append(dict(cohort=cohort, gene_set=label, n_genes=len(genes),
                         library=lib, rank=rank, term=term, pvalue=pval,
                         qvalue=qval, signature_hit=bool(fams),
                         families="+".join(fams)))
    print(f"    signature terms in top 15: {sig_hits}; "
          f"families present: {sorted(seen_families) or 'NONE'}")


def main():
    axis = sys.argv[1] if len(sys.argv) > 1 else "operator"
    cohorts = sys.argv[2:] or ["ccrcc", "luad", "ucec", "gbm", "pdac"]
    rows = []
    for c in cohorts:
        rd = results_dir(c)
        base_fp = os.path.join(rd, "residual_results_tumoronly.csv")
        sp_fp = os.path.join(rd, f"sitepack_{axis}.csv")
        if os.path.exists(base_fp):
            b = pd.read_csv(base_fp)
            sel = b[(b["fdr"] < 0.05) & (b["incremental_r2"] > 0)]
            run_set("baseline", c, sel["gene"].tolist(), rows)
        else:
            print(f"\n({c}: no {base_fp})")
        if os.path.exists(sp_fp):
            s = pd.read_csv(sp_fp)
            for label, fdr_col, incr_col in [
                    ("batchresid", "fdr_batchresid", "incremental_r2_batchresid"),
                    ("both", "fdr_both", "incremental_r2_both")]:
                if fdr_col not in s.columns or incr_col not in s.columns:
                    # an older sitepack CSV (e.g. left over from a smoke run) has a
                    # different schema; skip rather than crash the whole sweep
                    print(f"\n({c}: {os.path.basename(sp_fp)} predates the current "
                          f"schema -- missing {fdr_col}; rerun the sitepack)")
                    continue
                sel = s[(s[fdr_col] < 0.05) & (s[incr_col] > 0)]
                run_set(label, c, sel["gene"].tolist(), rows)
        else:
            print(f"\n({c}: no {sp_fp} yet)")

    if rows:
        out = f"{SCRIPTS}/enrichment_survival.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        print(f"\nwrote {out}")
        df = pd.DataFrame(rows)
        df = df[df["rank"] > 0]
        print("\n=== C2 SURVIVAL SUMMARY ===")
        print(f"  {'cohort':<7} {'set':<11} {'n_genes':>7} {'sig/15':>7}  families")
        for (c, gs), g in df.groupby(["cohort", "gene_set"]):
            fams = sorted({f for s in g["families"].fillna("") for f in s.split("+") if f})
            print(f"  {c:<7} {gs:<11} {g['n_genes'].iloc[0]:>7} "
                  f"{int(g['signature_hit'].sum()):>7}  {', '.join(fams) or 'NONE'}")
        print("\n  C2 passes for a cohort if the corrected sets keep the same "
              "families as its baseline set, even with far fewer genes.")


if __name__ == "__main__":
    main()
