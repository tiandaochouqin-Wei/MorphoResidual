#!/usr/bin/env python3
"""build_side_files.py -- F2 of review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md (section 10).

Builds, in the folder of this script (M/side/):
  families_reactome2022.tsv   columns set, gene (tab-separated, header line): the six signature families of
                              figures/enrichment_c2_levels.py applied to the term names of
                              figures/figdata/Reactome_2022.gmt (section 4)
  purity_ccrcc.tsv, purity_luad.tsv, purity_ucec.tsv, purity_pdac.tsv   columns patient, purity (section 6 item 2)
  build_side_files.log        the console output of this script (written by this script itself)

It reads ONLY: the GMT, the keyword map (imported from figures/enrichment_c2_levels.py), the covariate tables of
section 6 item 2 (the same sources and ID normalisation as composition_pathologist.py), the aliquot-to-case table,
and the 'patients' array of each M/export/<c>/inputs_<c>.npz (for coverage).  It loads no morphology, protein or
RNA array and computes no statistic that relates morphology to protein or mRNA.

Families: families_of / SIGNATURE_FAMILIES are IMPORTED from figures/enrichment_c2_levels.py (not re-implemented);
that module has no import-time side effects (its plotting-free main() is guarded by __main__; module level only
defines constants and functions; it imports gseapy, which must be installed).  Following run_level() of that
module (line 99: term = r["Term"].split(" R-HSA")[0]) the Reactome accession suffix is removed from the GMT term
name before matching.

Purity ID normalisation = norm_id() of composition_pathologist.py (lines 51-53): dots -> dashes, and C3* IDs are cut
to their first two dash parts.  Purity reader of local_theta.py purity_residualised(): TSV, first column patient
(str), second column numeric.
"""
import hashlib
import importlib.util
import io
import os
import sys
import datetime

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(M))      # MorphoResidual_paper/
FIG = os.path.join(ROOT, "figures")
DD = os.path.join(FIG, "figdata")
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
PURITY_COHORTS = ["ccrcc", "luad", "ucec", "pdac"]      # GBM has no purity estimate (section 6 item 2)

LOG = io.StringIO()


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def write_lf(path, text):
    assert "\r" not in text
    with open(path, "wb") as f:
        f.write(text.encode("utf-8"))


def norm_id(s):          # composition_pathologist.py lines 51-53
    s = str(s).strip().replace(".", "-")
    return "-".join(s.split("-")[:2]) if s.startswith("C3") else s


def inputs_patients(c):
    z = np.load(os.path.join(M, "export", c, f"inputs_{c}.npz"), allow_pickle=False)
    return [str(p) for p in z["patients"]]        # only this key is read


# ------------------------------------------------------------------------------------ families
def build_families():
    spec = importlib.util.spec_from_file_location("enrichment_c2_levels", os.path.join(FIG, "enrichment_c2_levels.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    fam_names = list(mod.SIGNATURE_FAMILIES.keys())
    say("families imported from figures/enrichment_c2_levels.py:", fam_names)
    gmt = os.path.join(DD, "Reactome_2022.gmt")
    genes = {f: set() for f in fam_names}
    terms = {f: 0 for f in fam_names}
    n_terms = 0
    n_nomatch = 0
    with open(gmt, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            parts = line.split("\t")
            name = parts[0].split(" R-HSA")[0]
            gl = [g.strip() for g in parts[2:] if g.strip()]
            n_terms += 1
            fams = mod.families_of(name)
            if not fams:
                n_nomatch += 1
            for fam in fams:
                terms[fam] += 1
                genes[fam].update(gl)
    rows = sorted((fam, g) for fam in fam_names for g in genes[fam])
    out = os.path.join(HERE, "families_reactome2022.tsv")
    write_lf(out, "set\tgene\n" + "".join(f"{a}\t{b}\n" for a, b in rows))
    say(f"\n[families] GMT terms read: {n_terms}; terms matching no family: {n_nomatch}")
    say(f"[families] wrote {out}")
    say(f"  sha256 {sha256_file(out)}  rows {len(rows)}")
    say("  family        terms_matched  genes")
    for fam in fam_names:
        say(f"  {fam:12s}  {terms[fam]:13d}  {len(genes[fam]):5d}")
    allg = set().union(*genes.values())
    say(f"  union of the six families: {len(allg)} distinct genes")
    return out


# ------------------------------------------------------------------------------------ purity
def purity_ccrcc():
    raw = pd.read_excel(os.path.join(DD, "subtype_ccrcc_raw.xlsx"), sheet_name="ESTIMATE scores", header=3, index_col=0)
    # The sheet holds TWO ESTIMATE blocks, each with a TumorPurity row: rows 1-8 "based on RNASeq Data" (first) and
    # rows 12-19 "based on Global Proteomics Data" (second).  composition_pathologist.py's est[["TumorPurity"]] therefore
    # returns both rows as two columns.  Here the FIRST (RNA-seq ESTIMATE) row is used, like the RNA-based purity of
    # the other cohorts; asserted below so a changed sheet cannot silently pick the other block.
    ix = [str(i) for i in raw.index]
    pos = [k for k, i in enumerate(ix) if i == "TumorPurity"]
    assert len(pos) == 2, f"expected two TumorPurity rows, found {len(pos)}"
    est = raw.iloc[[pos[0]]].T.apply(pd.to_numeric, errors="coerce")
    est.columns = ["TumorPurity"]
    say("[purity ccrcc] sheet has 2 TumorPurity rows (RNASeq-based first, Global-Proteomics-based second); first used; "
        f"proteomics-block row ignored (its mean {raw.iloc[pos[1]].apply(pd.to_numeric, errors='coerce').mean():.4g} "
        f"vs RNASeq-block {est['TumorPurity'].mean():.4g})")
    xw = pd.read_csv(os.path.join(ROOT, "server_export", "scripts", "aliquot_to_case_tumor.tsv"), sep="\t")
    xw = xw[xw["sample_type"].str.contains("Tumor", case=False, na=False)]
    a2c = dict(zip(xw["aliquot_id"].astype(str), xw["case_id"].map(norm_id)))
    est["case"] = [a2c.get(str(i)) for i in est.index]
    n_alq = int(est["case"].notna().sum())
    est = est.dropna(subset=["case"])
    s = est.groupby("case")["TumorPurity"].mean()
    return s, f"{n_alq} tumour aliquots mapped to {est['case'].nunique()} cases, mean per case"


def purity_ucec():
    d = pd.read_csv(os.path.join(DD, "subtype_ucec.txt"), sep="\t", encoding="latin-1", low_memory=False)
    n0 = len(d)
    if "Proteomics_Tumor_Normal" in d.columns:
        d = d[d["Proteomics_Tumor_Normal"].astype(str).str.lower().str.startswith("tumor")]
    d = d.copy()
    d["case"] = d["Proteomics_Participant_ID"].map(norm_id)
    d = d.set_index("case")
    v = pd.to_numeric(d["Purity_Cancer"], errors="coerce")
    ndup = int(v.index.duplicated().sum())
    v = v[~v.index.duplicated()]
    return v, f"{len(d)} tumour rows of {n0}; duplicate cases dropped (first kept): {ndup}"


def purity_luad():
    d = pd.read_csv(os.path.join(DD, "composition_luad_cbioportal.csv"))
    d["case"] = d["patientId"].map(norm_id)
    d = d.set_index("case")
    v = pd.to_numeric(d["TUMOR_PURITY_BYESTIMATE_RNASEQ"], errors="coerce")
    ndup = int(v.index.duplicated().sum())
    v = v[~v.index.duplicated()]
    return v, f"{len(d)} rows; duplicate cases dropped (first kept): {ndup}"


def purity_pdac():
    mp = pd.read_excel(os.path.join(DD, "subtype_pdac_raw.xlsx"), sheet_name="Molecular_phenotype_data")
    mp["case"] = mp["case_id"].map(norm_id)
    mp = mp.set_index("case")
    v = pd.to_numeric(mp["epithelial_cancer_deconv"], errors="coerce")
    ndup = int(v.index.duplicated().sum())
    v = v[~v.index.duplicated()]
    return v, f"{len(mp)} rows; duplicate cases dropped (first kept): {ndup}"


def build_purity(c, fn):
    s, note = fn()
    s = s.astype(float)
    s = s[np.isfinite(s.values)]
    out = os.path.join(HERE, f"purity_{c}.tsv")
    items = sorted(s.items())
    write_lf(out, "patient\tpurity\n" + "".join(f"{k}\t{repr(float(v))}\n" for k, v in items))
    pats = inputs_patients(c)
    have = [p for p in pats if p in s.index]
    lack = [p for p in pats if p not in s.index]
    say(f"\n[purity {c}] source: {note}")
    say(f"[purity {c}] wrote {out}")
    say(f"  sha256 {sha256_file(out)}  rows {len(items)}  (cases with a finite value in the source)")
    say(f"  coverage: {len(have)} of {len(pats)} inputs patients have a value; {len(lack)} lack one"
        + (f" ({', '.join(lack)})" if lack else ""))
    say(f"  purity range {s.min():.4g} .. {s.max():.4g}, mean {s.mean():.4g}; "
        f"file cases not in inputs: {len(set(s.index) - set(pats))}")
    return out


def main():
    say("build_side_files.py (F2 of POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md, section 10)")
    say("run time (UTC):", datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
    say("python", sys.version.split()[0], "pandas", pd.__version__, "numpy", np.__version__)
    say("reads: GMT, keyword map, covariate tables, aliquot map, 'patients' array of each inputs file; "
        "reads no morphology, protein or RNA values")
    build_families()
    for c, fn in (("ccrcc", purity_ccrcc), ("luad", purity_luad), ("ucec", purity_ucec), ("pdac", purity_pdac)):
        build_purity(c, fn)
    say("\nGBM: no purity estimate (only xCell stroma/immune); no purity file is built (section 6 item 2).")
    say("done.")
    write_lf(os.path.join(HERE, "build_side_files.log"), LOG.getvalue())


if __name__ == "__main__":
    main()
