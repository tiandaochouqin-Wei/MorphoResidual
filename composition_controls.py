#!/usr/bin/env python3
"""
composition_controls.py — MorphoResidual tissue-composition controls (Batch A-4).

Generalises the stromal control to THREE microenvironment/composition axes and a
combined test: does the WSI-morphology residual survive adjusting the baseline for
stromal, immune AND proliferation content simultaneously?

Reuses residual_analysis (RA) cv_r2 + loaders + PCA -> identical estimand.

Per significant gene, over mRNA baselines augmented with composition score(s) C:
    incr_orig          = R2(mRNA,WSI)          - R2(mRNA)              # published increment
    <C>_gain           = R2(mRNA,C)            - R2(mRNA)              # composition alone
    incr_W_over_<C>    = R2(mRNA,C,WSI)        - R2(mRNA,C)            # WSI on top of C
    incr_W_over_ALL    = R2(mRNA,stroma,immune,prolif,WSI) - R2(mRNA,stroma,immune,prolif)

Run like stroma_control (smp batch), MORPHO_SIG_CSV = pinned snapshot,
MORPHO_STROMA_OUT reused as the output-file path.
"""
import os, sys, numpy as np, pandas as pd
from sklearn.decomposition import PCA

import residual_analysis as RA

CV_FOLDS = getattr(RA, "CV_FOLDS", 5)

MARKERS = {
    "stroma": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A3", "FN1", "SPARC",
               "FAP", "PDGFRA", "PDGFRB", "ACTA2", "TAGLN", "DCN", "LUM", "POSTN",
               "THBS2", "FBN1", "MMP2", "TIMP1", "VIM", "COL4A1", "COL4A2", "SERPINH1"],
    "immune": ["PTPRC", "CD3D", "CD3E", "CD3G", "CD8A", "CD4", "CD2", "CD5", "CD68", "CD163",
               "MS4A1", "CD19", "NKG7", "GZMB", "GZMA", "PRF1", "FOXP3", "IL2RA", "CD14",
               "ITGAM", "LYZ", "HLA-DRA", "CD74", "IRF8", "CCL5"],
    "proliferation": ["MKI67", "PCNA", "TOP2A", "CCNB1", "CCNB2", "CDK1", "BUB1", "BUB1B",
                      "AURKA", "AURKB", "FOXM1", "TYMS", "RRM2", "CENPA", "PLK1", "CCNA2",
                      "CDC20", "KIF11", "TPX2", "UBE2C", "BIRC5", "CENPF", "MCM2", "MCM6"],
}


def comp_score(rna, markers):
    mk = [m for m in markers if m in rna.columns]
    if len(mk) < 6:
        return None, mk
    z = rna[mk].sub(rna[mk].mean()).div(rna[mk].std().replace(0, np.nan))
    return z.mean(axis=1).values.reshape(-1, 1), mk


def main():
    out = os.environ.get("MORPHO_STROMA_OUT", "composition_controls.csv")
    sig_csv = os.environ.get("MORPHO_SIG_CSV") or os.environ.get("MORPHO_MAIN_CSV")
    if not sig_csv:
        sys.exit("set MORPHO_SIG_CSV to the pinned residual_results_tumoronly.csv")

    rna = RA.load_rna_matrix()
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]
    wsi_pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
    print(f"[comp] {len(common)} patients; wsi PCs {wsi_pcs.shape}")

    scores = {}
    for fam, mk in MARKERS.items():
        s, used = comp_score(rna, mk)
        if s is None:
            print(f"[comp] {fam}: only {len(used)} markers -> skipped")
            continue
        scores[fam] = s
        print(f"[comp] {fam}: score from {len(used)}/{len(mk)} markers")
    fams = list(scores.keys())
    ALL = np.hstack([scores[f] for f in fams]) if fams else None

    md = pd.read_csv(sig_csv)
    if {"fdr", "incremental_r2"}.issubset(md.columns):
        sig = md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]["gene"].astype(str).tolist()
    else:
        sig = md["gene"].astype(str).tolist()
    sig = [g for g in sig if g in rna.columns and g in protein.columns]
    print(f"[comp] testing {len(sig)} significant genes over {len(fams)} composition axes")

    rows = []
    for g in sig:
        y = np.asarray(protein[g].values, float)
        v = ~np.isnan(y)
        if v.sum() < 30:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1)
        Wv = wsi_pcs[v]
        r2_rna = RA.cv_r2(xr, yv, CV_FOLDS)
        r2_both = RA.cv_r2(np.hstack([xr, Wv]), yv, CV_FOLDS)
        row = {"gene": g, "n": int(v.sum()), "incr_orig": r2_both - r2_rna}
        for f in fams:
            Cv = scores[f][v]
            r2_c = RA.cv_r2(np.hstack([xr, Cv]), yv, CV_FOLDS)
            r2_cw = RA.cv_r2(np.hstack([xr, Cv, Wv]), yv, CV_FOLDS)
            row[f"{f}_gain"] = r2_c - r2_rna
            row[f"incr_W_over_{f}"] = r2_cw - r2_c
        if ALL is not None:
            Av = ALL[v]
            r2_a = RA.cv_r2(np.hstack([xr, Av]), yv, CV_FOLDS)
            r2_aw = RA.cv_r2(np.hstack([xr, Av, Wv]), yv, CV_FOLDS)
            row["ALL_gain"] = r2_a - r2_rna
            row["incr_W_over_ALL"] = r2_aw - r2_a
        rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)

    n = len(res)
    print("\n==================== COMPOSITION CONTROLS ====================")
    print(f" significant genes tested            : {n}")
    print(f" median incremental R2 (WSI|mRNA)    : {res['incr_orig'].median():+.4f}")
    for f in fams:
        gain = res[f"{f}_gain"].median()
        wov = res[f"incr_W_over_{f}"].median()
        pos = 100 * (res[f"incr_W_over_{f}"] > 0).mean()
        print(f" {f:<14}: alone {gain:+.4f} | WSI over mRNA+{f:<9} {wov:+.4f} "
              f"({wov/res['incr_orig'].median()*100:.0f}% kept, {pos:.0f}% genes +)")
    if "incr_W_over_ALL" in res:
        wov = res["incr_W_over_ALL"].median()
        pos = 100 * (res["incr_W_over_ALL"] > 0).mean()
        print(f" {'ALL 3 axes':<14}: WSI over mRNA+stroma+immune+prolif {wov:+.4f} "
              f"({wov/res['incr_orig'].median()*100:.0f}% kept, {pos:.0f}% genes +)")
    print(" VERDICT: not composition-driven if WSI stays positive over mRNA+ALL.")
    print(f" -> {out}")


if __name__ == "__main__":
    main()
