#!/usr/bin/env python3
"""
covariate_control.py -- per-gene demographic controls for MorphoResidual.

Supersedes sex_control.py: same estimand, but it runs sex, age and race in one pass and
adds the capacity control that makes the numbers interpretable.

WHY THE CAPACITY CONTROL IS NOT OPTIONAL. Every covariate column added to the baseline
costs cross-validated degrees of freedom, so the "retained %" falls even when the
covariate carries no information at all. A one-column block already costs something; a
three-column block costs more. Reporting "morphology retains 85% over mRNA+sex" without
saying what a *useless* one-column block would have left invites a reviewer to read the
15% as confounding. So for every covariate block of width w this script also fits the
same thing with w columns of standard Gaussian noise (three independent draws, averaged),
exactly as the manuscript's dimension-matched noise control does. The comparison that
matters is covariate versus same-width noise, not covariate versus 100%.

Per gene, on that cohort's analysis set:

    incr_orig            = R2(mRNA, WSI)      - R2(mRNA)
    <C>_gain             = R2(mRNA, C)        - R2(mRNA)          # what C adds on its own
    incr_W_over_<C>      = R2(mRNA, C, WSI)   - R2(mRNA, C)       # morphology on top of C
    noise<w>_gain        = R2(mRNA, N_w)      - R2(mRNA)          # same width, pure noise
    incr_W_over_noise<w> = R2(mRNA, N_w, WSI) - R2(mRNA, N_w)

Blocks, and when each is testable:
    sex    1 col, male/female            all cohorts except UCEC (uniformly female)
    age    1 col, years from days_to_birth   all cohorts (age_at_index is empty in the
                                             GDC records for every cohort, so it is
                                             derived from days_to_birth)
    race   1 col, white vs non-white     only where the minority has >= MIN_RACE_MINORITY
                                         patients; CCRCC (6 non-white) and UCEC (4) are
                                         refused rather than reported as a token control
    all    the testable blocks together

EACH BLOCK USES ITS OWN PATIENT SET -- the patients for whom that block's covariates are
recorded -- and incr_orig, the block fit and the noise fit are all recomputed on it. A
single shared set would throw away patients for no benefit (race is unrecorded for 11 of
137 PDAC cases, 5 of 105 LUAD) and, worse, would stop incr_orig reproducing the published
per-gene increments, which is the check that proves the run is the published estimand.
So the output carries n_<block> and incr_orig_<block> per block, and every ratio below is
formed within one block.

Inputs
    the usual MORPHO_* environment (source morpho_env.sh <cancer>), plus
    MORPHO_SIG_CSV    pinned residual_results_tumoronly.csv for this cohort
    MORPHO_COVAR_TSV  covariates_by_case.tsv: cohort, case, sex, age_years, race
                      Built from the GDC clinical records (demographic.gender,
                      demographic.days_to_birth, demographic.race); not shipped with the
                      code, because this repository carries no patient-level data.
    MORPHO_COHORT     cohort key, lowercase
Output
    MORPHO_COVAR_OUT  default covariate_control_<cohort>.csv
"""
import os

# single-thread BLAS before numpy, as randproj_control.py does
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
MIN_PATIENTS = 30              # the published per-gene floor (residual_analysis._process_gene)
MIN_RACE_MINORITY = 20         # below this a race adjustment is not worth reporting
N_NOISE_DRAWS = 3              # as in the manuscript's dimension-matched noise control
NOISE_SEED = 0
HERE = os.path.dirname(os.path.abspath(__file__))


def guess_cohort():
    c = os.environ.get("MORPHO_COHORT", "").strip().lower()
    if c:
        return c
    sig = os.environ.get("MORPHO_SIG_CSV", "").replace("\\", "/")
    for k in ("ccrcc", "luad", "ucec", "gbm", "pdac"):
        if f"/{k}/" in sig:
            return k
    return "ccrcc" if sig else ""


def load_covariates(cohort):
    path = os.environ.get("MORPHO_COVAR_TSV", os.path.join(HERE, "covariates_by_case.tsv"))
    if not os.path.exists(path):
        sys.exit(f"[cov] missing covariate table: {path} (set MORPHO_COVAR_TSV)")
    t = pd.read_csv(path, sep="\t")
    need = {"cohort", "case", "sex", "age_years", "race"}
    if not need.issubset(t.columns):
        sys.exit(f"[cov] {path} must have columns {sorted(need)}; found {list(t.columns)}")
    t = t[t["cohort"].astype(str).str.lower() == cohort].copy()
    if t.empty:
        sys.exit(f"[cov] no rows for cohort '{cohort}' in {path}")
    t["case"] = t["case"].astype(str)
    return t.set_index("case")


def main():
    cohort = guess_cohort()
    if not cohort:
        sys.exit("[cov] set MORPHO_COHORT (ccrcc|luad|ucec|gbm|pdac)")
    out = os.environ.get("MORPHO_COVAR_OUT", f"covariate_control_{cohort}.csv")
    sig_csv = os.environ.get("MORPHO_SIG_CSV") or os.environ.get("MORPHO_MAIN_CSV")
    if not sig_csv:
        sys.exit("[cov] set MORPHO_SIG_CSV to the pinned residual_results_tumoronly.csv")

    cov = load_covariates(cohort)
    rna = RA.load_rna_matrix()
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]

    c = cov.reindex([str(x) for x in common])
    male = c["sex"].astype(str).str.lower().map({"male": 1.0, "female": 0.0})
    age = pd.to_numeric(c["age_years"], errors="coerce")
    race_raw = c["race"].astype(str).str.lower().str.strip()
    known_race = ~race_raw.isin(["not reported", "unknown", "nan", ""])
    nonwhite = np.where(known_race, (race_raw != "white").astype(float), np.nan)
    nonwhite = pd.Series(nonwhite, index=c.index)

    # which blocks are testable in this cohort?
    blocks, why = {}, []
    n_sex = int(male.notna().sum())
    if male.dropna().nunique() > 1 and min((male == 1).sum(), (male == 0).sum()) >= 5:
        blocks["sex"] = male
    else:
        why.append(f"sex: only one level present ({int((male == 1).sum())}m/{int((male == 0).sum())}f)")
    if age.notna().sum() >= MIN_PATIENTS and age.dropna().nunique() > 5:
        blocks["age"] = age
    else:
        why.append(f"age: only {int(age.notna().sum())} patients with a usable value")
    minority = int(min((nonwhite == 1).sum(), (nonwhite == 0).sum()))
    if minority >= MIN_RACE_MINORITY:
        blocks["race"] = nonwhite
    else:
        why.append(f"race: smaller group has {minority} patients (< {MIN_RACE_MINORITY}); "
                   f"a race adjustment here would not be interpretable")
    if not blocks:
        print(f"[cov] {cohort}: no testable covariate. " + "; ".join(why))
        return 0

    print(f"[cov] cohort {cohort}: {len(common)} patients; testable blocks {list(blocks)}")
    for w in why:
        print(f"[cov]   skipped {w}")

    # PCA is fitted once on all common patients, exactly as the published run does; the
    # per-block patient masks are applied to the resulting components, not before the fit.
    wsi_pcs = PCA(n_components=N_PCS, svd_solver="full").fit_transform(wsi.values)

    cols = {k: blocks[k].values.astype(float).reshape(-1, 1) for k in blocks}
    if "age" in cols:                       # standardise so its scale does not interact
        a = cols["age"]                     # with the ridge penalty on a 0/1 column
        fin = np.isfinite(a)
        cols["age"] = (a - np.nanmean(a)) / (np.nanstd(a) if np.nanstd(a) else 1.0)
        cols["age"][~fin] = np.nan
    C = {k: cols[k] for k in blocks}
    if len(blocks) > 1:
        C["all"] = np.hstack([cols[k] for k in blocks])

    rng = np.random.RandomState(NOISE_SEED)
    MASK, NOISE = {}, {}
    for k, M in C.items():
        MASK[k] = np.isfinite(M).all(axis=1)
        NOISE[k] = [rng.normal(size=(len(common), M.shape[1])) for _ in range(N_NOISE_DRAWS)]
        if MASK[k].sum() < MIN_PATIENTS:
            sys.exit(f"[cov] block {k}: only {int(MASK[k].sum())} patients; refusing to run")
    print(f"[cov] wsi PCs {wsi_pcs.shape}; blocks "
          f"{[(k, C[k].shape[1], int(MASK[k].sum())) for k in C]} (name, width, patients); "
          f"noise {N_NOISE_DRAWS} draws each; CV folds {CV_FOLDS}")

    md = pd.read_csv(sig_csv)
    if {"fdr", "incremental_r2"}.issubset(md.columns):
        sig = md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]["gene"].astype(str).tolist()
    else:
        sig = md["gene"].astype(str).tolist()
    sig = [g for g in sig if g in rna.columns and g in protein.columns]
    cap = int(os.environ.get("MORPHO_MAX_GENES", "0") or 0)
    if cap > 0:
        sig = sig[:cap]
        print(f"[cov] *** SMOKE TEST: capped at {cap} genes -- results are NOT the real run ***")
    print(f"[cov] testing {len(sig)} significant genes")

    rows = []
    for g in sig:
        y = np.asarray(protein[g].values, float)
        xg = np.asarray(rna[g].values, float)
        row, base = {"gene": g}, {}
        for k, M in C.items():
            v = (~np.isnan(y)) & MASK[k]
            if v.sum() < MIN_PATIENTS:
                continue
            yv, xr, Wv, Mv = y[v], xg[v].reshape(-1, 1), wsi_pcs[v], M[v]
            key = v.tobytes()
            if key not in base:          # incr_orig depends only on the mask -> cache per mask
                r = RA.cv_r2(xr, yv, CV_FOLDS)
                base[key] = (r, RA.cv_r2(np.hstack([xr, Wv]), yv, CV_FOLDS) - r)
            r2_rna, incr = base[key]
            r2_c = RA.cv_r2(np.hstack([xr, Mv]), yv, CV_FOLDS)
            gains, overs = [], []
            for Nz in NOISE[k]:
                Nv = Nz[v]
                r2_n = RA.cv_r2(np.hstack([xr, Nv]), yv, CV_FOLDS)
                gains.append(r2_n - r2_rna)
                overs.append(RA.cv_r2(np.hstack([xr, Nv, Wv]), yv, CV_FOLDS) - r2_n)
            row.update({f"n_{k}": int(v.sum()), f"incr_orig_{k}": incr,
                        f"{k}_gain": r2_c - r2_rna,
                        f"incr_W_over_{k}": RA.cv_r2(np.hstack([xr, Mv, Wv]), yv, CV_FOLDS) - r2_c,
                        f"noise_gain_{k}": float(np.mean(gains)),
                        f"incr_W_over_noise_{k}": float(np.mean(overs))})
        if len(row) > 1:
            rows.append(row)

    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)
    if res.empty:
        print(f"[cov] no gene cleared the {MIN_PATIENTS}-patient floor -> {out}")
        return 0

    print("\n==================== COVARIATE CONTROLS ====================")
    print(f" cohort : {cohort}   ({len(res)} genes)")
    print(f" {'block':<8s} {'w':>2s} {'pts':>5s} {'genes':>6s} {'med incr':>10s} "
          f"{'blk alone':>10s} {'WSI over':>10s} {'kept %':>8s} {'genes +':>8s} | {'noise kept %':>13s}")
    for k, M in C.items():
        col = f"incr_W_over_{k}"
        if col not in res:
            continue
        d = res[res[col].notna()]
        med = d[f"incr_orig_{k}"].median()
        kept = 100 * d[col].median() / med
        nk = 100 * d[f"incr_W_over_noise_{k}"].median() / med
        print(f" {k:<8s} {M.shape[1]:2d} {int(d[f'n_{k}'].max()):5d} {len(d):6d} {med:+10.4f} "
              f"{d[f'{k}_gain'].median():+10.4f} {d[col].median():+10.4f} {kept:7.1f}% "
              f"{100 * (d[col] > 0).mean():7.0f}% | {nk:12.1f}%")
    print(" READ IT AS: a covariate explains the increment only if its kept % is clearly")
    print(" BELOW the same-width noise kept %. Equal means the drop is the price of the")
    print(" extra column, not confounding. Each row is computed on its own patient set.")
    print(f" -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
