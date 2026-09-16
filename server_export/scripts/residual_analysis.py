#!/usr/bin/env python3
"""
Core pilot analysis (plan.md Phase A): does H&E morphology explain the part
of protein abundance that mRNA alone cannot -- and does it clear the pilot
go/no-go gate (>=50 proteins, FDR<0.05, positive incremental R²)?

Design (n=110 patients << D=1024 embedding dims, so this MUST be
cross-validated + regularized, not naive in-sample R² -- an unregularized
fit would trivially "explain" the residual with 1024 free parameters):
  1. Per-patient WSI embedding = mean-pool ALL tiles across ALL of that
     patient's slides (mean-pooling, not ABMIL -- first-pass choice already
     agreed with the user), then PCA-reduced to 20 PCs (keeps the WSI design
     matrix small relative to n, standard practice at this sample size).
  2. Per protein g: 5-fold CV R² for (i) protein~mRNA, (ii) protein~mRNA+WSI,
     (iii) protein~WSI-only. Incremental R² = CV_R²(mRNA+WSI) - CV_R²(mRNA).
  3. Permutation null: permute the patient<->WSI-embedding correspondence
     (shared across all proteins per draw -- this is what makes testing
     thousands of proteins computationally tractable), recompute incremental
     R² for every protein under each permutation, use the resulting
     per-protein null distribution for an empirical p-value.
  4. BH-FDR across all tested proteins.
  5. Gate: count proteins with incremental R² > 0 AND FDR < 0.05.

Run on mn02 (or any CPU node) -- no GPU needed for this step.
"""
import multiprocessing as mp
import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")

N_WORKERS = int(os.environ.get("MORPHO_N_WORKERS", max(1, (os.cpu_count() or 4) - 1)))

ROOT = Path(os.environ.get("MORPHO_ROOT", "/public/home/fjhui/ZW"))
RNA_DIR = ROOT / "omics" / "rna"
# TUMOR-ONLY manifest (sample_type column, Primary Tumor rows only). The old
# manifest_rna.tsv dropped the sample_type field entirely, so load_rna_matrix's
# "first file per case wins" was silently picking tumor OR normal at random for
# the 168/261 cases that have both -- a tumor/normal contamination fixed here.
RNA_MANIFEST = Path(os.environ.get("MORPHO_RNA_MANIFEST",
                                   str(Path(__file__).parent / "manifest_rna_tumor.tsv")))
PROTEIN_TSV = Path(os.environ.get("MORPHO_PROTEIN_TSV",
                                  str(ROOT / "omics" / "protein" / "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv")))
# TUMOR-ONLY crosswalk: 103 ccRCC Primary Tumor aliquots (110 tumor - 7 official
# non-ccRCC exclusions), 1:1 aliquot->case. The old aliquot_to_case.tsv held all
# 194 tumor+normal aliquots and load_protein_matrix's groupby.mean AVERAGED each
# patient's tumor and normal protein values together -- fixed here.
ALIQUOT_XWALK = Path(os.environ.get("MORPHO_ALIQUOT_XWALK",
                                    str(Path(__file__).parent / "aliquot_to_case_tumor.tsv")))
# Per-slide sample_type (authoritative, from GDC slide entities). Needed because
# the -NN slide suffix does NOT encode tumor/normal by a fixed rule (e.g. -21..-26
# each appear in BOTH classes) -- must look up each slide individually.
SLIDE_TYPE_MAP = Path(os.environ.get("MORPHO_SLIDE_MAP",
                                     str(Path(__file__).parent / "slide_type_map.tsv")))
WSI_EMB_DIR = Path(os.environ.get("MORPHO_WSI_EMB_DIR", str(ROOT / "WSI" / "emb" / "ccrcc_real")))
OUT_DIR = Path(os.environ.get("MORPHO_OUT", str(ROOT / "results")))
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 7 tumor samples flagged by CPTAC as molecularly non-ccRCC and excluded from all
# official ccRCC downstream analyses (per PDC000127 study note). Dropped from every
# modality so the pilot cohort is a homogeneous ccRCC set.
NON_CCRCC_CASES = {"C3L-00359", "C3N-00313", "C3N-00435", "C3N-00492",
                   "C3N-00832", "C3N-01175", "C3N-01180"}

N_PCS = 20
N_PERM = 1000
CV_FOLDS = 5
RIDGE_ALPHA = 1.0
FDR_THRESHOLD = 0.05
GATE_MIN_PROTEINS = 50   # pilot go/no-go: >= this many significant proteins
MIN_PATIENTS_TO_RUN = 20  # sanity floor for CV to be meaningful at all
RNG = np.random.RandomState(0)


def load_rna_matrix():
    manifest = pd.read_csv(RNA_MANIFEST, sep="\t")
    # Defensive: keep only Primary Tumor rows even if a mixed manifest is passed.
    if "sample_type" in manifest.columns:
        manifest = manifest[manifest["sample_type"] == "Primary Tumor"]
    manifest = manifest[~manifest["case_submitter_id"].isin(NON_CCRCC_CASES)]
    # A case may have several tumor RNA files (40 cases do) -> average them per
    # gene rather than arbitrarily keeping the first.
    per_case = {}  # case -> list of per-gene Series
    for _, r in manifest.iterrows():
        case = r["case_submitter_id"]
        fp = RNA_DIR / r["file_name"]
        if not fp.exists():
            continue
        try:
            df = pd.read_csv(fp, sep="\t", skiprows=1, usecols=["gene_name", "tpm_unstranded"])
        except Exception:
            continue
        df = df.dropna(subset=["gene_name"]).groupby("gene_name")["tpm_unstranded"].mean()
        per_case.setdefault(case, []).append(df)
    rows = []
    for case, series_list in per_case.items():
        merged = pd.concat(series_list, axis=1).mean(axis=1)  # avg tumor replicates
        rows.append(merged.rename(case))
    mat = pd.concat(rows, axis=1).T  # patients x genes
    mat = np.log2(mat.fillna(0) + 1)
    print(f"RNA matrix (tumor-only): {mat.shape[0]} patients x {mat.shape[1]} genes")
    return mat


def load_protein_matrix():
    # Crosswalk is tumor-only (Primary Tumor aliquots), so only tumor protein
    # columns are selected -- normal-tissue channels are never mapped in.
    xwalk = pd.read_csv(ALIQUOT_XWALK, sep="\t").set_index("aliquot_id")["case_id"].to_dict()
    df = pd.read_csv(PROTEIN_TSV, sep="\t", index_col=0)
    log_ratio_cols = [c for c in df.columns if c.endswith(" Log Ratio") and "Unshared" not in c]
    keep = {}
    for c in log_ratio_cols:
        aliquot = c.split(" ")[0]
        case = xwalk.get(aliquot)
        if case and case not in NON_CCRCC_CASES:
            keep[c] = case
    sub = df[list(keep.keys())].rename(columns=keep)
    mat = sub.T.groupby(level=0).mean()  # tumor crosswalk is 1:1, so this is a no-op safeguard
    print(f"Protein matrix (tumor-only): {mat.shape[0]} patients x {mat.shape[1]} genes")
    return mat


def load_wsi_embeddings():
    # Authoritative per-slide sample_type map (slide_submitter_id -> type). We
    # mean-pool ONLY Primary Tumor slides per case; normal-adjacent slides (180 of
    # 561 in this cohort) previously got averaged into the same per-patient vector,
    # so "WSI predicts protein residual" could have been reading tumor-vs-normal.
    slide_df = pd.read_csv(SLIDE_TYPE_MAP, sep="\t")
    slide_type = slide_df.set_index("slide_submitter_id")["sample_type"].to_dict()

    files = sorted(WSI_EMB_DIR.glob("*.pt"))
    slide_re = re.compile(r"^(C3[A-Z]-\d{5}-\d+)")   # full slide id: case + section code
    case_re = re.compile(r"^(C3[A-Z]-\d{5})")
    by_case = {}
    n_tumor = n_normal = n_unmatched = 0
    for f in files:
        sm = slide_re.match(f.stem)
        cm = case_re.match(f.stem)
        if not cm:
            continue
        case = cm.group(1)
        if case in NON_CCRCC_CASES:
            continue
        slide_id = sm.group(1) if sm else None
        stype = slide_type.get(slide_id)
        if stype == "Primary Tumor":
            n_tumor += 1
        elif stype is None:
            n_unmatched += 1
            continue  # unknown sample type -> exclude to stay conservative
        else:
            n_normal += 1
            continue  # Solid Tissue Normal (or other) -> exclude
        d = torch.load(f, map_location="cpu")
        by_case.setdefault(case, []).append(d["embeddings"])
    print(f"WSI slides: {n_tumor} tumor kept, {n_normal} normal dropped, "
          f"{n_unmatched} unmatched dropped")
    rows = {}
    for case, tensors in by_case.items():
        allt = torch.cat(tensors, dim=0)
        rows[case] = allt.mean(dim=0).numpy()
    mat = pd.DataFrame(rows).T
    print(f"WSI matrix (tumor-only, mean-pooled): {mat.shape[0]} patients x {mat.shape[1]} dims")
    return mat


def cv_r2(X, y, folds):
    # Closed-form ridge (no sklearn Ridge/StandardScaler object creation) --
    # confirmed live 2026-07-17 that per-call sklearn overhead, not the
    # actual linear algebra, was the bottleneck: 127 parallel workers ran
    # for ~2 hours (each burning real CPU time, not deadlocked) without
    # finishing 200 genes. At this problem size (~90 train rows x <=21
    # features) the matrix math is trivial; the cost was Ridge/StandardScaler
    # Python-level overhead x ~5000 fits/gene x ~9800 genes.
    kf = KFold(n_splits=folds, shuffle=True, random_state=0)
    preds = np.zeros_like(y, dtype=float)
    ridge_I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        X_tr, X_te, y_tr = X[tr], X[te], y[tr]
        x_mean = X_tr.mean(axis=0)
        x_std = X_tr.std(axis=0)
        x_std[x_std == 0] = 1.0
        Xtr_s = (X_tr - x_mean) / x_std
        Xte_s = (X_te - x_mean) / x_std
        y_mean = y_tr.mean()
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + ridge_I, Xtr_s.T @ (y_tr - y_mean))
        preds[te] = Xte_s @ w + y_mean
    ss_res = np.sum((y - preds) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


_WORKER = {}  # populated once per worker process via _init_worker, avoids re-pickling per gene


def _init_worker(rna_dict, protein_dict, wsi_pcs, perm_indices):
    _WORKER["rna"] = rna_dict
    _WORKER["protein"] = protein_dict
    _WORKER["wsi_pcs"] = wsi_pcs
    _WORKER["perm_indices"] = perm_indices


def _process_gene(gene):
    protein_col = _WORKER["protein"][gene]
    y = np.asarray(protein_col)
    valid = ~np.isnan(y)
    if valid.sum() < 30:
        return None
    y_v = y[valid]
    x_rna = np.asarray(_WORKER["rna"][gene])[valid].reshape(-1, 1)
    wsi_pcs = _WORKER["wsi_pcs"]
    x_both = np.hstack([x_rna, wsi_pcs[valid]])
    r2_rna = cv_r2(x_rna, y_v, CV_FOLDS)
    r2_both = cv_r2(x_both, y_v, CV_FOLDS)
    incr = r2_both - r2_rna

    perm_indices = _WORKER["perm_indices"]
    null_incr = np.empty(len(perm_indices))
    for i, perm in enumerate(perm_indices):
        x_wsi_perm = wsi_pcs[perm][valid]
        x_both_perm = np.hstack([x_rna, x_wsi_perm])
        null_incr[i] = cv_r2(x_both_perm, y_v, CV_FOLDS) - r2_rna
    pval = (np.sum(null_incr >= incr) + 1) / (len(perm_indices) + 1)

    return {"gene": gene, "n": int(valid.sum()), "r2_rna": r2_rna,
            "r2_mrna_wsi": r2_both, "incremental_r2": incr, "pval": pval}


def bh_fdr(pvals):
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    fdr = np.empty(n)
    fdr[order] = np.clip(ranked, 0, 1)
    return fdr


def main():
    rna = load_rna_matrix()
    protein = load_protein_matrix()
    wsi = load_wsi_embeddings()

    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    print(f"\ncommon patients across all 3: {len(common)}")
    if len(common) < MIN_PATIENTS_TO_RUN:
        sys.exit(f"FAIL: only {len(common)} patients overlap -- cannot proceed")

    common_genes = sorted(set(rna.columns) & set(protein.columns))
    print(f"common genes (RNA ∩ protein): {len(common_genes)}")

    rna = rna.loc[common, common_genes]
    protein = protein.loc[common, common_genes]
    wsi_raw = wsi.loc[common].values

    # svd_solver='full' is load-bearing: at this matrix shape sklearn's auto
    # solver is randomized SVD with random_state=None, which made this script
    # non-reproducible (GBM: 699/700/704/701 significant across four identical
    # runs, 90.4% set overlap). 'full' is exact and deterministic.
    pca = PCA(n_components=min(N_PCS, len(common) - 1, wsi_raw.shape[1]),
              svd_solver="full")
    wsi_pcs = pca.fit_transform(wsi_raw)
    print(f"WSI PCA: {wsi_pcs.shape[1]} components, "
          f"{pca.explained_variance_ratio_.sum():.1%} variance retained")

    n = len(common)

    # IMPORTANT: every gene with enough data gets the FULL permutation test.
    # An earlier version of this script pre-filtered to genes with a
    # positive point-estimate incremental R² before permutation-testing only
    # those -- that is circular analysis (selecting on the same statistic
    # being tested biases almost everything to look "significant", since the
    # selected set is, by construction, the noise-favored tail). Confirmed
    # live 2026-07-16: that version returned 3102/3103 "significant" genes,
    # which is implausibly close to 100% and was the tell that something was
    # wrong. No pre-filter this time -- BH-FDR is computed across the true
    # full test universe.
    perm_indices = [RNG.permutation(n) for _ in range(N_PERM)]
    rna_dict = {g: rna[g].values for g in common_genes}
    protein_dict = {g: protein[g].values for g in common_genes}

    print(f"parallelizing across {N_WORKERS} worker processes "
          f"(override with MORPHO_N_WORKERS env var)")
    results = []
    tested = 0
    with mp.Pool(N_WORKERS, initializer=_init_worker,
                 initargs=(rna_dict, protein_dict, wsi_pcs, perm_indices)) as pool:
        for r in pool.imap_unordered(_process_gene, common_genes, chunksize=8):
            tested += 1
            if r is not None:
                results.append(r)
            if tested % 200 == 0:
                print(f"  {tested}/{len(common_genes)} genes permutation-tested...")

    res = pd.DataFrame(results)
    res["fdr"] = bh_fdr(res["pval"].values)
    res = res.sort_values("incremental_r2", ascending=False)
    out_csv = OUT_DIR / "residual_results_tumoronly.csv"  # keep the old contaminated result for A/B comparison
    res.to_csv(out_csv, index=False)

    sig = res[(res["fdr"] < FDR_THRESHOLD) & (res["incremental_r2"] > 0)]
    print(f"\ntested {len(res)} genes (full permutation test, no pre-filter, TUMOR-ONLY)")
    print(f"significant (FDR<{FDR_THRESHOLD}, incremental_r2>0): {len(sig)}")
    print(f"\ntop 15 by incremental R2:")
    print(res.head(15)[["gene", "n", "r2_rna", "r2_mrna_wsi", "incremental_r2", "pval", "fdr"]]
          .to_string(index=False))

    gate = len(sig) >= GATE_MIN_PROTEINS
    print(f"\n{'PASS' if gate else 'FAIL'}: {len(sig)} proteins clear the pilot gate "
          f"({'>=' if gate else '<'}{GATE_MIN_PROTEINS})")
    print(f"full results: {out_csv}")
    sys.exit(0 if gate else 1)


if __name__ == "__main__":
    main()
