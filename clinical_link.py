#!/usr/bin/env python3
"""
clinical_link.py — MorphoResidual translational hook (pilot).

Question: does a MORPHOLOGY-derived residual pathway score stratify patients by
clinical grade / stage / survival?  If yes, H&E alone yields a clinically
meaningful readout of the post-transcriptional programme.

Reuses residual_analysis (RA) loaders + PCA + the EXACT ridge CV scheme, so the
per-patient scores are on the same estimand as the main analysis.  Clinical
variables are pulled from the GDC API (needs internet -> run on the login node).

Per patient, for each pathway family f (translation / secretion-ER / matrisome):
    residual_g        = protein_g - OOF_ridge(mRNA_g)          # measured residual
    score_meas_f      = mean over family-f significant genes of z(residual_g)
    score_morph_f     = OOF ridge prediction of score_meas_f from the 20 WSI PCs
Then associate score_meas_f and score_morph_f with:
    grade / stage  -> Spearman rho
    survival       -> median-split log-rank (manual) + optional Cox (lifelines)

Env (same as stroma_control plus a clinical cache path):
    MORPHO_ROOT/RNA_MANIFEST/ALIQUOT_XWALK/SLIDE_MAP/PROTEIN_TSV/WSI_EMB_DIR  (via morpho_env.sh)
    MORPHO_SIG_CSV        cohort's pinned residual_results_tumoronly.csv
    MORPHO_CLIN_OUT       output prefix (default clinical_link_<cohort>)
    MORPHO_CLIN_PROJECT   GDC project_id filter (default CPTAC-3)
Run on the login node (internet):
    cd /public/home/fjhui/ZW/scripts
    source morpho_env.sh ccrcc
    export MORPHO_SIG_CSV=/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/ccrcc/residual_results_tumoronly.csv
    export MORPHO_CLIN_OUT=/public/home/fjhui/ZW/results/clinical_link_ccrcc
    python -u clinical_link.py
"""
import os, re, sys, numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

import residual_analysis as RA  # loaders + constants (env read at import)

RIDGE_ALPHA = getattr(RA, "RIDGE_ALPHA", 1.0)
CV_FOLDS = getattr(RA, "CV_FOLDS", 5)

# ---- pathway families: symbol rule, intersected with significant genes present in proteome ----
SECRETION_ER = {
    "SSR1", "SSR2", "SSR3", "SSR4", "SEC11A", "SEC11C", "SEC61A1", "SEC61B", "SEC61G",
    "SRP9", "SRP14", "SRP54", "SRP68", "SRP72", "SRPRA", "SRPRB", "SPCS1", "SPCS2", "SPCS3",
    "RPN1", "RPN2", "DDOST", "STT3A", "STT3B", "OST4", "PIGS", "TANGO2", "SEC13", "SEC23A",
    "SEC24C", "SEC31A", "LMAN1", "SURF4", "TMED2", "TMED9", "SSR3",
}
MATRISOME = {
    "LUM", "ASPN", "DCN", "POSTN", "COL1A1", "COL1A2", "COL3A1", "HAPLN3", "FN1", "SPARC",
    "THBS2", "FBN1", "COL5A1", "COL6A3", "OGN", "PRELP", "GPX3",
}
FAMILIES = {
    "translation": lambda g: bool(re.match(r"^(RPS|RPL|EEF|EIF)\d", g)) or g in {"RACK1", "FAU"},
    "secretion_ER": lambda g: g in SECRETION_ER,
    "matrisome": lambda g: g in MATRISOME,
}


def cv_predict(X, y):
    """Out-of-fold predictions, identical scheme to RA.cv_r2 (KFold shuffle seed 0,
    train-fold z-standardisation, closed-form ridge with lambda=RIDGE_ALPHA)."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    kf = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=0)
    preds = np.full(len(y), np.nan)
    ridge_I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        Xtr, Xte, ytr = X[tr], X[te], y[tr]
        xm = Xtr.mean(0); xs = Xtr.std(0); xs[xs == 0] = 1.0
        Xtr_s = (Xtr - xm) / xs; Xte_s = (Xte - xm) / xs; ym = ytr.mean()
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + ridge_I, Xtr_s.T @ (ytr - ym))
        preds[te] = Xte_s @ w + ym
    return preds


def logrank(time, event, group):
    """Two-group log-rank; group is a boolean array. Returns (z, p)."""
    time = np.asarray(time, float); event = np.asarray(event, int); group = np.asarray(group, bool)
    ev_times = np.unique(time[event == 1])
    O1 = E1 = V = 0.0
    for t in ev_times:
        at_risk = time >= t
        n = at_risk.sum(); n1 = (at_risk & group).sum()
        d = ((time == t) & (event == 1)).sum()
        d1 = ((time == t) & (event == 1) & group).sum()
        if n > 1:
            E1 += d * n1 / n
            V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
        O1 += d1
    if V <= 0:
        return float("nan"), float("nan")
    from scipy import stats as _st
    z = (O1 - E1) / np.sqrt(V)
    return z, 2 * _st.norm.sf(abs(z))


def spearman(a, b):
    from scipy import stats as _st
    m = ~(pd.isna(a) | pd.isna(b))
    if m.sum() < 8:
        return float("nan"), float("nan"), int(m.sum())
    r, p = _st.spearmanr(np.asarray(a)[m], np.asarray(b)[m])
    return r, p, int(m.sum())


def parse_grade(g):
    if not isinstance(g, str):
        return np.nan
    m = re.search(r"G\s*([1-4])", g)
    return float(m.group(1)) if m else np.nan


ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4}


def parse_stage(s):
    if not isinstance(s, str):
        return np.nan
    m = re.search(r"Stage\s+(IV|III|II|I)", s)
    if m:
        return float(ROMAN[m.group(1)])
    m = re.search(r"\b(IV|III|II|I)\b", s)  # FIGO etc.
    return float(ROMAN[m.group(1)]) if m else np.nan


def fetch_gdc_clinical(cases, project):
    """Pull grade/stage/survival for the cohort cases from the GDC cases endpoint."""
    try:
        import requests
    except Exception as e:
        print(f"[clin] requests unavailable ({e}); skipping clinical pull")
        return None
    filt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "submitter_id", "value": list(cases)}}]}
    if project:
        filt["content"].append({"op": "in", "content": {"field": "project.project_id", "value": [project]}})
    body = {"filters": filt, "expand": "diagnoses,demographic", "size": len(cases) + 50, "format": "json"}
    try:
        r = requests.post("https://api.gdc.cancer.gov/cases", json=body, timeout=90)
        r.raise_for_status()
        hits = r.json()["data"]["hits"]
    except Exception as e:
        print(f"[clin] GDC query failed ({e}); skipping clinical association")
        return None
    rows = []
    for h in hits:
        demo = h.get("demographic") or {}
        dxs = h.get("diagnoses") or [{}]
        dx = dxs[0]
        rows.append({
            "case": h.get("submitter_id"),
            "grade": parse_grade(dx.get("tumor_grade")),
            "stage": parse_stage(dx.get("ajcc_pathologic_stage") or dx.get("figo_stage")
                                 or dx.get("uicc_pathologic_stage")),
            "vital": (demo.get("vital_status") or "").strip(),
            "dtd": demo.get("days_to_death"),
            "dtlf": demo.get("days_to_last_follow_up") or dx.get("days_to_last_follow_up"),
            "age": demo.get("age_at_index"),
        })
    df = pd.DataFrame(rows).dropna(subset=["case"]).set_index("case")
    return df


def main():
    out = os.environ.get("MORPHO_CLIN_OUT", "clinical_link")
    sig_csv = os.environ.get("MORPHO_SIG_CSV") or os.environ.get("MORPHO_MAIN_CSV")
    project = os.environ.get("MORPHO_CLIN_PROJECT", "CPTAC-3")
    if not sig_csv:
        sys.exit("set MORPHO_SIG_CSV to the pinned residual_results_tumoronly.csv")

    rna = RA.load_rna_matrix()
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    rna, protein, wsi = rna.loc[common], protein.loc[common], wsi.loc[common]
    wsi_pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
    print(f"[clin] {len(common)} patients; wsi PCs {wsi_pcs.shape}")

    md = pd.read_csv(sig_csv)
    if {"fdr", "incremental_r2"}.issubset(md.columns):
        sig = md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]["gene"].astype(str).tolist()
    else:
        sig = md["gene"].astype(str).tolist()
    sig = [g for g in sig if g in rna.columns and g in protein.columns]
    print(f"[clin] {len(sig)} significant residual genes")

    # per-patient residual matrix (protein - OOF ridge on mRNA), z-scored per gene
    resid = pd.DataFrame(index=common, columns=sig, dtype=float)
    for g in sig:
        y = protein[g].values.astype(float)
        v = ~np.isnan(y)
        if v.sum() < 30:
            continue
        pred = cv_predict(rna[g].values[v], y[v])
        col = np.full(len(common), np.nan)
        col[np.where(v)[0]] = y[v] - pred
        resid[g] = col
    resid = (resid - resid.mean()) / resid.std(ddof=0)

    # family scores (measured) and their morphology-predicted versions
    scores = pd.DataFrame(index=common)
    fam_genes = {}
    for fam, rule in FAMILIES.items():
        cols = [g for g in sig if rule(g)]
        fam_genes[fam] = cols
        if len(cols) < 5:
            print(f"[clin] family {fam}: only {len(cols)} genes -> skipped")
            continue
        meas = resid[cols].mean(axis=1)  # nan-aware
        scores[f"{fam}_meas"] = meas
        vv = ~meas.isna().values
        morph = np.full(len(common), np.nan)
        morph[vv] = cv_predict(wsi_pcs[vv], meas.values[vv])
        scores[f"{fam}_morph"] = morph
        # how well morphology predicts the measured pathway score (sanity)
        mm = ~np.isnan(morph) & ~meas.isna().values
        r2 = 1 - np.sum((meas.values[mm] - morph[mm]) ** 2) / np.sum(
            (meas.values[mm] - meas.values[mm].mean()) ** 2)
        print(f"[clin] family {fam}: {len(cols)} genes; morphology CV-R2 of measured score = {r2:+.3f}")

    clin = fetch_gdc_clinical(common, project)
    scores.to_csv(f"{out}_scores.csv")

    print("\n==================== CLINICAL LINK ====================")
    if clin is None or clin.empty:
        print(" clinical pull unavailable; wrote per-patient scores only ->", f"{out}_scores.csv")
        return
    clin = clin.reindex(common)
    # survival time / event
    dead = clin["vital"].str.lower().eq("dead")
    time = np.where(dead, pd.to_numeric(clin["dtd"], errors="coerce"),
                    pd.to_numeric(clin["dtlf"], errors="coerce"))
    event = dead.astype(int).values
    cov_grade = int(clin["grade"].notna().sum())
    cov_stage = int(clin["stage"].notna().sum())
    cov_surv = int((~np.isnan(time) & (np.asarray(time) >= 0)).sum())
    n_death = int(np.nansum(event))
    print(f" clinical coverage: grade {cov_grade}/{len(common)}, stage {cov_stage}/{len(common)}, "
          f"survival {cov_surv}/{len(common)} ({n_death} deaths)")
    clin.assign(time=time, event=event).to_csv(f"{out}_clinical.csv")

    rows = []
    for col in scores.columns:
        s = scores[col].values
        rg, pg, ng = spearman(s, clin["grade"].values)
        rs, ps, ns = spearman(s, clin["stage"].values)
        # survival: median split (high vs low score)
        m = ~np.isnan(s) & ~np.isnan(time) & (np.asarray(time) >= 0)
        lr_z = lr_p = float("nan"); nlr = int(m.sum())
        if m.sum() >= 20 and n_death >= 5:
            med = np.median(s[m])
            grp = s[m] > med
            lr_z, lr_p = logrank(np.asarray(time)[m], event[m], grp)
        rows.append({"score": col, "grade_rho": rg, "grade_p": pg,
                     "stage_rho": rs, "stage_p": ps,
                     "surv_logrank_z": lr_z, "surv_logrank_p": lr_p, "n": nlr})
    res = pd.DataFrame(rows)
    pd.set_option("display.width", 160); pd.set_option("display.max_columns", 20)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    res.to_csv(f"{out}_assoc.csv", index=False)
    print(f"\n -> {out}_scores.csv / {out}_clinical.csv / {out}_assoc.csv")
    print(" READ: for each family, *_morph is the H&E-only score; a significant grade/stage rho")
    print("       or survival log-rank p on *_morph is the translational payload.")


if __name__ == "__main__":
    main()
