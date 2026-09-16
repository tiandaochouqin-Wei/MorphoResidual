#!/usr/bin/env python3
"""cptac2_enrichment.py — pathway enrichment (Reactome via Enrichr) for the
morphology-predictable gene set found by cptac2_residual.py in a CPTAC-2 organ
(BRCA/COAD/OV), using the identical method and family taxonomy as
server_export/scripts/enrichment_check.py (the one used for the five CPTAC-3
discovery organs). The point of this script is a single yes/no: does the
significant set converge on the same translation / ER-secretion / splicing /
folding families that the five discovery organs converged on, or on something
else / nothing at all?

Run:  python cptac2_enrichment.py TCGA-BRCA
Output: cptac2_<short>/results/enrichment_<short>.tsv + a printed summary
"""
import json, os, sys, time, urllib.parse, urllib.request
import pandas as pd

ROOT = "/public/home/fjhui/ZW"
ADDLIST = "https://maayanlab.cloud/Enrichr/addList"
ENRICH = "https://maayanlab.cloud/Enrichr/enrich"
LIBS = ["Reactome_Pathways_2024", "Reactome_2022", "Reactome_2016"]

# identical taxonomy to enrichment_check.py -- keep in sync if that file changes
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


def post_list(genes):
    boundary = "----MorphoResidualBoundary"
    parts = []
    for name, val in (("list", "\n".join(genes)), ("description", "morphoresidual_cptac2")):
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


def run_set(cohort, genes, rows):
    genes = [g for g in genes if isinstance(g, str) and g]
    print(f"\n--- {cohort}: {len(genes)} morphology-predictable genes", flush=True)
    if len(genes) < 10:
        print("    too few genes to enrich (this itself is the finding)")
        rows.append(dict(cohort=cohort, n_genes=len(genes), library="", rank=-1,
                          term="", pvalue=None, qvalue=None, signature_hit=None))
        return
    uid = post_list(genes[:3000])
    if uid is None:
        print("    Enrichr unreachable"); return
    for lib in LIBS:
        res = enrich(uid, lib)
        if res:
            break
    else:
        print("    no Reactome library returned results"); return
    print(f"    library {lib}, {len(res)} terms")
    seen_families = set(); sig_hits = 0
    for rank, t in enumerate(res[:15], 1):
        term, pval, qval = t[1], t[2], t[6]
        fams = families_of(term)
        seen_families.update(fams); sig_hits += int(bool(fams))
        mark = f"  <== {'+'.join(fams)}" if fams else ""
        print(f"    {rank:>2}. q={qval:.3g}  {term[:66]}{mark}")
        rows.append(dict(cohort=cohort, n_genes=len(genes), library=lib, rank=rank,
                          term=term, pvalue=pval, qvalue=qval, signature_hit=bool(fams),
                          families="+".join(fams)))
    print(f"    signature terms in top 15: {sig_hits}; families present: {sorted(seen_families) or 'NONE'}")


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: cptac2_enrichment.py TCGA-BRCA")
    proj = sys.argv[1]; short = proj.split("-")[1].lower()
    RD = f"{ROOT}/cptac2_{short}/results"
    fp = f"{RD}/{short}_protein_results.csv"
    if not os.path.exists(fp):
        sys.exit(f"missing {fp} -- run cptac2_residual.py --project {proj} first")
    r = pd.read_csv(fp)
    sig = r[(r["fdr"] < 0.05) & (r["incremental_r2"] > 0)]
    rows = []
    run_set(short, sig["gene"].tolist(), rows)
    if rows:
        out = f"{RD}/enrichment_{short}.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        print(f"\nwrote {out}")
        print("\n  Compare against the five CPTAC-3 discovery organs: translation in the "
              "lung, matrisome/secretion in the pancreas, splicing in the brain, folding/"
              "glycosylation in the kidney and uterus. This organ converging on ANY of the "
              "same families is the finding; converging on none of them is also a finding.")


if __name__ == "__main__":
    main()
