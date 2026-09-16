#!/usr/bin/env python3
"""
stroma_control.py — Devil's-Advocate control for MorphoResidual.

Does the WSI-morphology residual survive adjusting the baseline for a
transcriptome-derived STROMAL score?  Reuses residual_analysis.cv_r2 + loaders
so the estimand is byte-for-byte the same as the main analysis.

Per significant gene (same CV_FOLDS / closed-form ridge as the main run):
    r2_rna   = cv_r2([mRNA])                         # reproduce baseline
    r2_both  = cv_r2([mRNA, WSI])                    # reproduce main model
    r2_mS    = cv_r2([mRNA, stroma])                 # mRNA + stroma
    r2_mSW   = cv_r2([mRNA, stroma, WSI])            # + WSI on top of stroma
    incr_orig      = r2_both - r2_rna                # == reported incremental R^2
    stroma_gain    = r2_mS   - r2_rna                # how much stroma alone adds
    incr_W_over_S  = r2_mSW  - r2_mS                 # WSI gain ON TOP of stroma  <-- key

VERDICT: signal is NOT stromal if incr_W_over_S stays positive for most
significant genes and tracks incr_orig; it IS (partly) stromal if it collapses.

Run from the scripts/ dir with the SAME MORPHO_* env vars as residual_analysis.py,
plus MORPHO_MAIN_CSV pointing at that cohort's residual_results_tumoronly.csv:
    cd /public/home/fjhui/ZW/scripts
    MORPHO_ROOT=... MORPHO_RNA_MANIFEST=... MORPHO_ALIQUOT_XWALK=... MORPHO_SLIDE_MAP=... \
    MORPHO_PROTEIN_TSV=... MORPHO_WSI_EMB_DIR=... MORPHO_MAIN_CSV=.../residual_results_tumoronly.csv \
    MORPHO_OUT=.../stroma_control_<cohort>.csv  python -u stroma_control.py
"""
import os, sys, numpy as np, pandas as pd
from sklearn.decomposition import PCA

import residual_analysis as RA  # reuse loaders, cv_r2, constants (env read at import)

CV_FOLDS = getattr(RA, "CV_FOLDS", 5)

# Compact ESTIMATE-style stromal signature (canonical CAF / ECM markers).
STROMAL_MARKERS = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A3", "FN1", "SPARC",
    "FAP", "PDGFRA", "PDGFRB", "ACTA2", "TAGLN", "DCN", "LUM", "POSTN",
    "THBS2", "FBN1", "MMP2", "TIMP1", "VIM", "COL4A1", "COL4A2", "SERPINH1",
]


def main():
    # dedicated names so `source morpho_env.sh <cancer>` (which sets MORPHO_OUT to a
    # DIRECTORY) does not collide with our file output.
    out = os.environ.get("MORPHO_STROMA_OUT", "stroma_control_results.csv")
    main_csv = os.environ.get("MORPHO_SIG_CSV") or os.environ.get("MORPHO_MAIN_CSV")
    if not main_csv:
        sys.exit("set MORPHO_SIG_CSV to the cohort's residual_results_tumoronly.csv")

    rna = RA.load_rna_matrix()            # patients x genes, log2(TPM+1), symbols
    protein = RA.load_protein_matrix()    # patients x genes
    wsi = RA.load_wsi_embeddings()        # patients x dims
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]
    wsi_pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
    print(f"[stroma] {len(common)} patients; wsi PCs {wsi_pcs.shape}")

    # per-patient stromal score = mean z-scored stromal-marker expression
    mk = [m for m in STROMAL_MARKERS if m in rna.columns]
    if len(mk) < 8:
        sys.exit(f"only {len(mk)} stromal markers found in RNA — check symbols")
    z = rna[mk].sub(rna[mk].mean()).div(rna[mk].std().replace(0, np.nan))
    S = z.mean(axis=1).values.reshape(-1, 1)
    print(f"[stroma] stromal score from {len(mk)}/{len(STROMAL_MARKERS)} markers")

    # originally-significant genes.  Accept either a full residual-results CSV
    # (filter FDR<0.05 & incremental_r2>0) or a pre-filtered gene list.
    md = pd.read_csv(main_csv)
    if {"fdr", "incremental_r2"}.issubset(md.columns):
        sig_df = md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]
    else:
        sig_df = md
    sig = sig_df["gene"].astype(str).tolist()
    sig = [g for g in sig if g in rna.columns and g in protein.columns]
    print(f"[stroma] testing {len(sig)} significant genes")

    rows = []
    for g in sig:
        y = np.asarray(protein[g].values, dtype=float)
        v = ~np.isnan(y)
        if v.sum() < 30:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, dtype=float)[v].reshape(-1, 1)
        Sv = S[v]
        Wv = wsi_pcs[v]
        r2_rna = RA.cv_r2(xr, yv, CV_FOLDS)
        r2_both = RA.cv_r2(np.hstack([xr, Wv]), yv, CV_FOLDS)
        r2_mS = RA.cv_r2(np.hstack([xr, Sv]), yv, CV_FOLDS)
        r2_mSW = RA.cv_r2(np.hstack([xr, Sv, Wv]), yv, CV_FOLDS)
        rows.append({
            "gene": g, "n": int(v.sum()),
            "incr_orig": r2_both - r2_rna,
            "stroma_gain": r2_mS - r2_rna,
            "incr_W_over_S": r2_mSW - r2_mS,
        })
    res = pd.DataFrame(rows)
    # self-check: our reproduced incr_orig must track the PUBLISHED incremental_r2
    if "incremental_r2" in md.columns:
        rep = md.set_index("gene")["incremental_r2"]
        res["incr_reported"] = res["gene"].map(rep)
        ok = res.dropna(subset=["incr_reported"])
        rho = ok["incr_orig"].corr(ok["incr_reported"]) if len(ok) > 2 else float("nan")
        mad = (ok["incr_orig"] - ok["incr_reported"]).abs().median() if len(ok) else float("nan")
    else:
        rho = mad = float("nan")
    res.to_csv(out, index=False)

    n = len(res)
    surv = int((res["incr_W_over_S"] > 0).sum())
    keep_half = int((res["incr_W_over_S"] > 0.5 * res["incr_orig"]).sum())
    print("\n==================== STROMA CONTROL ====================")
    print(f" estimand self-check: corr(reproduced, published incr) = {rho:.3f}")
    print(f"                      median|reproduced - published|   = {mad:.4f}")
    print(f"   (corr ~1 & MAD tiny => cv_r2/PCA match the main run)")
    print(f" significant genes tested          : {n}")
    print(f" median incremental R2 (WSI|mRNA)  : {res['incr_orig'].median():+.4f}")
    print(f" median stroma-alone gain          : {res['stroma_gain'].median():+.4f}")
    print(f" median WSI gain OVER (mRNA+stroma): {res['incr_W_over_S'].median():+.4f}")
    print(f" genes keeping WSI gain > 0        : {surv}/{n} ({100*surv/max(n,1):.1f}%)")
    print(f" genes keeping >=50% of orig gain  : {keep_half}/{n} ({100*keep_half/max(n,1):.1f}%)")
    print(" VERDICT: NOT stromal if WSI gain over (mRNA+stroma) stays positive for")
    print("          most genes and near incr_orig; stromal if it collapses to ~0.")
    print(f" -> {out}")


if __name__ == "__main__":
    main()
