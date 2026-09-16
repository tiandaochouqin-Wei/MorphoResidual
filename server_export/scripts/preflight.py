#!/usr/bin/env python3
"""
Cheap preflight for residual_analysis.py on a new cohort.

Loads the same three matrices through the analysis module's own loaders and
reports the 3-way patient intersection and gene overlap, WITHOUT running the
1000-permutation test. A misconfigured cohort (wrong crosswalk, wrong slide map,
wrong protein table) shows up here in a couple of minutes instead of after an
hour of permutation testing.

Usage: preflight.py <cancer>
"""
import os
import sys

cancer = sys.argv[1]
ROOT = "/public/home/fjhui/ZW"
SCRIPTS = f"{ROOT}/scripts"

PROTEIN = {
    "ucec": f"{ROOT}/ucec/omics/protein/CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Proteome.tmt10.tsv",
    "pdac": f"{ROOT}/pdac/omics/protein/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv",
    "gbm":  f"{ROOT}/gbm/omics/protein/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv",
}

os.environ["MORPHO_ROOT"] = f"{ROOT}/{cancer}"
os.environ["MORPHO_RNA_MANIFEST"] = f"{SCRIPTS}/manifest_rna_tumor_{cancer}.tsv"
os.environ["MORPHO_ALIQUOT_XWALK"] = f"{SCRIPTS}/aliquot_to_case_tumor_{cancer}.tsv"
os.environ["MORPHO_SLIDE_MAP"] = f"{SCRIPTS}/slide_type_map_{cancer}.tsv"
os.environ["MORPHO_PROTEIN_TSV"] = PROTEIN[cancer]
os.environ["MORPHO_WSI_EMB_DIR"] = f"{ROOT}/{cancer}/WSI/emb"
os.environ["MORPHO_OUT"] = f"{ROOT}/{cancer}/results"

sys.path.insert(0, SCRIPTS)
import residual_analysis as ra  # noqa: E402

print(f"\n### preflight {cancer}")
print(f"  RNA manifest : {ra.RNA_MANIFEST}")
print(f"  xwalk        : {ra.ALIQUOT_XWALK}")
print(f"  slide map    : {ra.SLIDE_TYPE_MAP}")
print(f"  protein      : {ra.PROTEIN_TSV}")
print(f"  wsi emb      : {ra.WSI_EMB_DIR}")
print(f"  out          : {ra.OUT_DIR}")

# guard: the hardcoded ccRCC exclusion list must not silently bite another cohort
import pandas as pd  # noqa: E402
xw = pd.read_csv(ra.ALIQUOT_XWALK, sep="\t")
collide = set(xw["case_id"]) & ra.NON_CCRCC_CASES
print(f"  NON_CCRCC_CASES collision with this cohort: {sorted(collide) or 'none'}")

rna = ra.load_rna_matrix()
protein = ra.load_protein_matrix()
wsi = ra.load_wsi_embeddings()

common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
genes = sorted(set(rna.columns) & set(protein.columns))
print(f"\n  patients: rna={len(rna.index)} protein={len(protein.index)} wsi={len(wsi.index)}")
print(f"  3-way common patients : {len(common)}")
print(f"  common genes          : {len(genes)}")
print(f"  pairwise: rna&prot={len(set(rna.index)&set(protein.index))} "
      f"rna&wsi={len(set(rna.index)&set(wsi.index))} "
      f"prot&wsi={len(set(protein.index)&set(wsi.index))}")
print(f"\n  VERDICT: {'OK to run' if len(common) >= 40 else 'TOO FEW PATIENTS - investigate'}")
