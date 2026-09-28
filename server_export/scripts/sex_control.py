#!/usr/bin/env python3
"""
sex_control.py -- per-gene sex control for MorphoResidual.

Answers the reviewer question the score-level analysis cannot: does the published
per-gene morphology increment survive putting sex into the mRNA baseline, and is the
phenomenon present in each sex separately?

Reuses residual_analysis (RA) cv_r2 + loaders + PCA, so the estimand is identical to
the published one -- the only change is which columns enter the baseline.

Per gene in the cohort's own significant set, on the patients with a known sex and a
measured protein value (one patient set per gene, so every quantity below is computed
on the same rows and is directly comparable):

    incr_orig        = R2(mRNA, WSI)      - R2(mRNA)          # the published increment
    sex_gain         = R2(mRNA, sex)      - R2(mRNA)          # sex alone
    incr_W_over_sex  = R2(mRNA, sex, WSI) - R2(mRNA, sex)     # morphology on top of sex

and, only where BOTH strata have at least MIN_STRATUM patients for that gene, the same
increment fitted within each sex separately (underpowered by construction -- CV with
~21 features on 40-70 rows -- so it is reported as a direction check, not an estimate):

    incr_male / incr_female

Both or neither: a male-only column with no female counterpart cannot be read. On the
real cohorts this means the stratified block runs for GBM (55/44) and PDAC (71/66) and
is skipped for CCRCC (77/26) and LUAD (69/36).

UCEC is skipped: it is uniformly female, so sex is constant and cannot be adjusted for.

Inputs
    the usual MORPHO_* environment (source morpho_env.sh <cancer>), plus
    MORPHO_SIG_CSV   pinned residual_results_tumoronly.csv for this cohort
    MORPHO_SEX_TSV   sex_by_case.tsv with columns cohort, case, sex (default: next to this
                     script). Built from the GDC clinical records for the five cohorts --
                     one row per case, sex taken verbatim from demographic.gender -- and
                     not shipped with the code, because this repository carries no
                     patient-level data. It is a three-column projection of the same GDC
                     pull the manuscript's Data availability section points to.
    MORPHO_COHORT    cohort key, lowercase (default: taken from MORPHO_SIG_CSV's path)
Output
    MORPHO_SEX_OUT   default sex_control_<cohort>.csv

Run exactly like the other controls:
    source ${ZW}/scripts/morpho_env.sh <cancer>
    MORPHO_SIG_CSV=... python -u sex_control.py
"""
import os

# Single-thread BLAS before numpy is imported, as randproj_control.py does: the only
# multithreaded call here is one (n x 1024) SVD, and unpinned threads on a shared node
# are a nuisance to everyone else.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import residual_analysis as RA

CV_FOLDS = getattr(RA, "CV_FOLDS", 5)
N_PCS = getattr(RA, "N_PCS", 20)
MIN_PATIENTS = 30         # the published per-gene floor: residual_analysis._process_gene
                          # and composition_controls.py both use < 30, not MIN_PATIENTS_TO_RUN
MIN_STRATUM = 40          # floor for the per-sex fits; below this the CV is not worth reporting
HERE = os.path.dirname(os.path.abspath(__file__))


def guess_cohort():
    c = os.environ.get("MORPHO_COHORT", "").strip().lower()
    if c:
        return c
    sig = os.environ.get("MORPHO_SIG_CSV", "")
    for k in ("ccrcc", "luad", "ucec", "gbm", "pdac"):
        if f"/{k}/" in sig.replace("\\", "/"):
            return k
    # ccRCC's results sit at ${ZW}/results, with no cohort directory in the path
    return "ccrcc" if sig else ""


def load_sex(cohort):
    path = os.environ.get("MORPHO_SEX_TSV", os.path.join(HERE, "sex_by_case.tsv"))
    if not os.path.exists(path):
        sys.exit(f"[sex] missing sex table: {path} (set MORPHO_SEX_TSV)")
    t = pd.read_csv(path, sep="\t")
    need = {"cohort", "case", "sex"}
    if not need.issubset(t.columns):
        sys.exit(f"[sex] {path} must have columns {sorted(need)}; found {list(t.columns)}")
    t = t[t["cohort"].astype(str).str.lower() == cohort]
    t = t[t["sex"].astype(str).str.lower().isin(["male", "female"])]
    if t.empty:
        sys.exit(f"[sex] no male/female rows for cohort '{cohort}' in {path}")
    return dict(zip(t["case"].astype(str), t["sex"].astype(str).str.lower().eq("male").astype(float)))


def main():
    cohort = guess_cohort()
    if not cohort:
        sys.exit("[sex] set MORPHO_COHORT (ccrcc|luad|gbm|pdac)")
    if cohort == "ucec":
        print("[sex] UCEC is uniformly female: sex is constant and cannot be adjusted for. Nothing to do.")
        return 0
    out = os.environ.get("MORPHO_SEX_OUT", f"sex_control_{cohort}.csv")
    sig_csv = os.environ.get("MORPHO_SIG_CSV") or os.environ.get("MORPHO_MAIN_CSV")
    if not sig_csv:
        sys.exit("[sex] set MORPHO_SIG_CSV to the pinned residual_results_tumoronly.csv")

    sex_of = load_sex(cohort)
    rna = RA.load_rna_matrix()
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]

    sex_vec = np.array([sex_of.get(str(c), np.nan) for c in common], float)
    known = ~np.isnan(sex_vec)
    n_m, n_f = int(np.nansum(sex_vec == 1)), int(np.nansum(sex_vec == 0))
    print(f"[sex] cohort {cohort}: {len(common)} patients; sex known for {int(known.sum())} "
          f"({n_m} male / {n_f} female)")
    if known.sum() < MIN_PATIENTS or min(n_m, n_f) < 5:
        sys.exit(f"[sex] not enough patients with both sexes to run ({n_m}m/{n_f}f)")

    wsi_pcs = PCA(n_components=N_PCS, svd_solver="full").fit_transform(wsi.values)
    print(f"[sex] wsi PCs {wsi_pcs.shape}; CV folds {CV_FOLDS}")

    md = pd.read_csv(sig_csv)
    if {"fdr", "incremental_r2"}.issubset(md.columns):
        sig = md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]["gene"].astype(str).tolist()
    else:
        sig = md["gene"].astype(str).tolist()
    sig = [g for g in sig if g in rna.columns and g in protein.columns]
    cap = int(os.environ.get("MORPHO_MAX_GENES", "0") or 0)
    if cap > 0:
        sig = sig[:cap]
        print(f"[sex] *** SMOKE TEST: capped at {cap} genes -- results are NOT the real run ***")
    print(f"[sex] testing {len(sig)} significant genes")

    rows = []
    for g in sig:
        y = np.asarray(protein[g].values, float)
        v = (~np.isnan(y)) & known
        if v.sum() < MIN_PATIENTS:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1)
        Wv = wsi_pcs[v]
        Sv = sex_vec[v].reshape(-1, 1)
        r2_rna = RA.cv_r2(xr, yv, CV_FOLDS)
        r2_both = RA.cv_r2(np.hstack([xr, Wv]), yv, CV_FOLDS)
        r2_sex = RA.cv_r2(np.hstack([xr, Sv]), yv, CV_FOLDS)
        r2_sexw = RA.cv_r2(np.hstack([xr, Sv, Wv]), yv, CV_FOLDS)
        row = {"gene": g, "n": int(v.sum()),
               "n_male": int((Sv == 1).sum()), "n_female": int((Sv == 0).sum()),
               "incr_orig": r2_both - r2_rna,
               "sex_gain": r2_sex - r2_rna,
               "incr_W_over_sex": r2_sexw - r2_sex}
        # Direction check within each sex. All-or-nothing: a male-only column with no
        # female counterpart is not interpretable, so both strata must clear the floor.
        # Never allowed to break the main result.
        masks = {"male": Sv[:, 0] == 1, "female": Sv[:, 0] == 0}
        if min(int(m.sum()) for m in masks.values()) >= MIN_STRATUM:
            try:
                for tag, mask in masks.items():
                    a = RA.cv_r2(xr[mask], yv[mask], CV_FOLDS)
                    b = RA.cv_r2(np.hstack([xr[mask], Wv[mask]]), yv[mask], CV_FOLDS)
                    row[f"incr_{tag}"] = b - a
            except Exception as e:                                # noqa: BLE001
                print(f"[sex] stratified fit skipped for {g}: {e}")
                row.pop("incr_male", None)
                row.pop("incr_female", None)
        rows.append(row)

    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)
    if res.empty:
        print(f"[sex] no gene cleared the {MIN_PATIENTS}-patient floor; nothing to summarise -> {out}")
        return 0

    med = res["incr_orig"].median()
    kept = res["incr_W_over_sex"].median()
    print("\n==================== SEX CONTROL ====================")
    print(f" cohort                                : {cohort}  ({n_m} male / {n_f} female)")
    print(f" significant genes tested              : {len(res)}")
    print(f" median incremental R2 (WSI|mRNA)      : {med:+.4f}")
    print(f" median sex alone over mRNA            : {res['sex_gain'].median():+.4f}")
    print(f" median WSI over mRNA+sex              : {kept:+.4f}"
          f"   ({100 * kept / med:.1f}% kept, {100 * (res['incr_W_over_sex'] > 0).mean():.0f}% genes +)")
    for tag in ("male", "female"):
        col = f"incr_{tag}"
        if col in res and res[col].notna().any():
            s = res[col].dropna()
            print(f" median increment, {tag+' only':<20}: {s.median():+.4f}   "
                  f"({100 * (s > 0).mean():.0f}% genes +, n_genes {len(s)}; underpowered, direction only)")
    print(" VERDICT: not sex-driven if WSI stays positive over mRNA+sex and the kept % is near 100.")
    print(f" -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
