#!/usr/bin/env python3
"""
pathomics_link.py — which nuclear feature encodes which residual programme (B-9).

Reuses residual_analysis.load_wsi_embeddings (point MORPHO_WSI_EMB_DIR at the
pathomics output dir) to get a per-patient morphometry matrix, then correlates each
interpretable feature with the MEASURED residual pathway scores written by
clinical_link.py. Spearman, so scale-invariant.

Run with the cohort's env (source morpho_env.sh <c>), then:
    export MORPHO_WSI_EMB_DIR=<.../WSI/pathomics/...>
    export MORPHO_CLIN_OUT=/public/home/fjhui/ZW/results/clinical_link_<c>
    python pathomics_link.py
"""
import os, numpy as np, pandas as pd
from scipy import stats
import residual_analysis as RA

FEAT_NAMES = ["nuclear_density", "nucleus_count", "mean_nuc_area",
              "nuc_area_std", "chromatin_od", "eosin_fraction"]


def main():
    prefix = os.environ["MORPHO_CLIN_OUT"]
    cohort = os.path.basename(prefix).replace("clinical_link_", "")
    morph = RA.load_wsi_embeddings()  # patients x 6 (mean-pooled morphometry)
    M = pd.DataFrame(np.asarray(morph), index=morph.index, columns=FEAT_NAMES)

    scores = pd.read_csv(f"{prefix}_scores.csv", index_col=0)
    meas = [c for c in scores.columns if c.endswith("_meas")]
    df = M.join(scores[meas], how="inner")
    print(f"\n===== {cohort}: {len(df)} patients =====")
    print(f"[patho] features {FEAT_NAMES}")
    print(f"[patho] scores   {meas}")

    rows = []
    for f in FEAT_NAMES:
        for s in meas:
            r, p = stats.spearmanr(df[f], df[s], nan_policy="omit")
            rows.append({"feature": f, "score": s, "rho": r, "p": p})
    res = pd.DataFrame(rows)
    piv = res.pivot(index="feature", columns="score", values="rho")
    print("\nSpearman rho (nuclear feature x residual pathway score):")
    print(piv.to_string(float_format=lambda x: f"{x:+.2f}"))
    res.to_csv(f"{prefix}_pathomics_corr.csv", index=False)
    print(f"\n -> {prefix}_pathomics_corr.csv")
    print(" READ: e.g. a positive nuclear_density x translation_meas rho means denser")
    print("       nuclei mark higher translation-residual -> that is what morphology reads.")


if __name__ == "__main__":
    main()
