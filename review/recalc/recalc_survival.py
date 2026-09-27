#!/usr/bin/env python3
"""Recompute CCRCC clinical associations + survival from local files (read-only on inputs).
Mirrors clinical_link.py (Spearman, median-split log-rank) and clinical_cox.py (numpy Breslow Cox,
z-standardised covariates, ridge l2=1e-3), and adds Wald CIs, Schoenfeld (Grambsch-Therneau) PH test,
reverse-KM median follow-up, EPV, grade/stage univariable and adjusted HRs.
Run:  python review/recalc/recalc_survival.py   (from the paper root)
"""
import json
import numpy as np
import os
import pandas as pd
from scipy import stats
from scipy.stats import chi2

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sc = pd.read_csv(f"{ROOT}/figures/figdata/scores_ccrcc.csv", index_col=0)
cl = pd.read_csv(f"{ROOT}/figures/figdata/clinical_link_ccrcc_clinical.csv", index_col=0)
out = {}


# ---------------- Spearman (full cohort, mirrors clinical_link.spearman) ----------------
def spearman(a, b):
    m = ~(pd.isna(a) | pd.isna(b))
    r, p = stats.spearmanr(np.asarray(a)[m], np.asarray(b)[m])
    return float(r), float(p), int(m.sum())


clin = cl.reindex(sc.index)
time_all = clin["time"].values.astype(float)
event_all = clin["event"].values.astype(int)


def logrank_link(time, event, group):
    """clinical_link.logrank (z-based, two-sided normal). Returns z, p, O1, E1."""
    time = np.asarray(time, float); event = np.asarray(event, int); group = np.asarray(group, bool)
    O1 = E1 = V = 0.0
    for t in np.unique(time[event == 1]):
        at_risk = time >= t
        n = at_risk.sum(); n1 = (at_risk & group).sum()
        d = ((time == t) & (event == 1)).sum()
        d1 = ((time == t) & (event == 1) & group).sum()
        if n > 1:
            E1 += d * n1 / n
            V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
        O1 += d1
    z = (O1 - E1) / np.sqrt(V)
    return float(z), float(2 * stats.norm.sf(abs(z))), float(O1), float(E1)


fam_rows = []
for col in sc.columns:
    s = sc[col].values
    rg, pg, ng = spearman(s, clin["grade"].values)
    rs, ps, ns = spearman(s, clin["stage"].values)
    m = ~np.isnan(s) & ~np.isnan(time_all) & (time_all >= 0)
    med = np.median(s[m])
    z_gt, p_gt, O1, E1 = logrank_link(time_all[m], event_all[m], s[m] > med)   # clinical_link uses '>'
    z_ge, p_ge, _, _ = logrank_link(time_all[m], event_all[m], s[m] >= med)    # km_survival uses '>='
    fam_rows.append(dict(score=col, grade_rho=rg, grade_p=pg, n_grade=ng,
                         stage_rho=rs, stage_p=ps, n_stage=ns,
                         surv_n=int(m.sum()), surv_logrank_z=z_gt, surv_logrank_p=p_gt,
                         surv_logrank_p_ge=p_ge,
                         n_high=int((s[m] > med).sum()), n_low=int((s[m] <= med).sum()),
                         obs_high=O1, exp_high=E1))
assoc = pd.DataFrame(fam_rows)
print("=== CCRCC full clinical_link family (6 scores x 3 tests) ===")
print(assoc.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
out["ccrcc_assoc_full_family"] = assoc.to_dict(orient="records")

# ---------------- survival subset ----------------
SCORE = "translation_morph"
df = sc.join(cl[["grade", "stage", "time", "event"]])
df = df[df["time"].notna() & (df["time"] >= 0)].dropna(subset=[SCORE, "grade", "stage"])
t = df["time"].values.astype(float)
e = df["event"].values.astype(int)
n = len(df); d = int(e.sum())
n_neg = int((cl["time"] < 0).sum())
print(f"\n=== survival subset: n={n}, deaths={d}, excluded (time<0)={n_neg}, "
      f"tied death times={d - len(np.unique(t[e == 1]))}")
out["survival"] = dict(
    n=n, deaths=d, censored=n - d, excluded_negative_followup=n_neg, n_cohort=len(sc),
    grade_missing=int(cl["grade"].isna().sum()), stage_missing=int(cl["stage"].isna().sum()),
    age_available=int(cl["age"].notna().sum()),
    tied_death_times=int(d - len(np.unique(t[e == 1]))),
    grade_coding="ordinal 1-4 (GDC tumor_grade G1-G4), z-standardised in Cox (clinical_cox.zstd)",
    stage_coding="ordinal 1-4 (AJCC pathologic stage I-IV), z-standardised in Cox",
    grade_counts={int(k): int(v) for k, v in df["grade"].value_counts().sort_index().items()},
    stage_counts={int(k): int(v) for k, v in df["stage"].value_counts().sort_index().items()},
)


# median follow-up
def km_curve(tt, ee):
    ts = np.unique(tt); S = 1.0; xs = []; ys = []
    for x in ts:
        nr = (tt >= x).sum(); dd = ((tt == x) & (ee == 1)).sum()
        S *= (1 - dd / nr) if nr > 0 else 1
        xs.append(x); ys.append(S)
    return np.array(xs), np.array(ys)


def km_median(tt, ee):
    xs, ys = km_curve(tt, ee)
    idx = np.where(ys <= 0.5)[0]
    return float(xs[idx[0]]) if len(idx) else float("nan")


med_fu_revkm = km_median(t, 1 - e)   # reverse KM: censoring as the event
med_fu_censored = float(np.median(t[e == 0]))
med_fu_all = float(np.median(t))
med_os = km_median(t, e)
out["survival"].update(
    median_followup_reverse_KM_days=med_fu_revkm,
    median_followup_censored_only_days=med_fu_censored,
    median_time_all_days=med_fu_all, median_OS_days=med_os,
    followup_range_days=[float(t.min()), float(t.max())],
    time_definition="days_to_death if vital_status==Dead else days_to_last_follow_up (GDC); "
                    "event = death from any cause",
)
print(f"median follow-up: reverse-KM {med_fu_revkm:.0f} d ({med_fu_revkm / 365.25:.2f} y); "
      f"median of censored times {med_fu_censored:.0f} d; median of all times {med_fu_all:.0f} d; "
      f"median OS {'not reached' if np.isnan(med_os) else med_os}")

# median-split log-rank (km_survival.py rule: >= median -> high)
s = df[SCORE].values; med = np.median(s); hi = s >= med
z_lr, p_lr, O1, E1 = logrank_link(t, e, hi)
p_lr_chi2 = float(1 - chi2.cdf(z_lr ** 2, 1))
out["survival"]["logrank_median_split"] = dict(
    n_high=int(hi.sum()), n_low=int((~hi).sum()),
    deaths_high=int(e[hi].sum()), deaths_low=int(e[~hi].sum()),
    O_high=O1, E_high=E1, chi2=float(z_lr ** 2), p=p_lr_chi2,
    split_rule="score >= median -> high (km_survival.py); '>' gives the identical split",
)
print(f"log-rank median split: high n={hi.sum()} ({e[hi].sum()} deaths), low n={(~hi).sum()} "
      f"({e[~hi].sum()} deaths), chi2={z_lr ** 2:.3f}, p={p_lr_chi2:.4f}")
os_land = {}
for lab, m in [("high", hi), ("low", ~hi)]:
    xs, ys = km_curve(t[m], e[m])
    for yr in (3, 5):
        k = np.where(xs <= yr * 365.25)[0]
        v = float(ys[k[-1]]) if len(k) else 1.0
        os_land[f"{lab}_{yr}y"] = v
        print(f"  {lab}: {yr}-y OS = {v:.3f}")
out["survival"]["landmark_OS"] = os_land


# ---------------- Cox (mirrors clinical_cox.cox_fit) ----------------
def cox_fit(X, time, event, l2=1e-3, max_iter=100, tol=1e-9):
    X = np.asarray(X, float); time = np.asarray(time, float); event = np.asarray(event, int)
    n_, p_ = X.shape; beta = np.zeros(p_)
    ev_times = np.unique(time[event == 1])
    risk_masks = [(time >= tt) for tt in ev_times]
    death_masks = [((time == tt) & (event == 1)) for tt in ev_times]
    d_counts = [int(m.sum()) for m in death_masks]
    Xdeath_sum = [X[m].sum(0) for m in death_masks]
    H = None
    for _ in range(max_iter):
        eta = X @ beta; eta -= eta.max(); w = np.exp(eta)
        grad = -l2 * beta; H = l2 * np.eye(p_)
        for rm, dd, xds in zip(risk_masks, d_counts, Xdeath_sum):
            Xr = X[rm]; wr = w[rm]; S0 = wr.sum()
            S1 = (Xr * wr[:, None]).sum(0); xbar = S1 / S0
            S2 = (Xr[:, :, None] * Xr[:, None, :] * wr[:, None, None]).sum(0)
            grad += xds - dd * xbar
            H += dd * (S2 / S0 - np.outer(xbar, xbar))
        step = np.linalg.solve(H, grad); beta += step
        if np.max(np.abs(step)) < tol:
            break
    cov = np.linalg.inv(H); se = np.sqrt(np.diag(cov)); z = beta / se
    return beta, se, 2 * stats.norm.sf(np.abs(z)), cov


def zstd(a):
    a = np.asarray(a, float); sd = np.nanstd(a)
    return (a - np.nanmean(a)) / (sd if sd > 0 else 1.0)


def summarise(names, beta, se, p):
    return [dict(term=nm, beta=float(b), se=float(s_), HR=float(np.exp(b)),
                 CI95=[float(np.exp(b - 1.96 * s_)), float(np.exp(b + 1.96 * s_))], p=float(pp))
            for nm, b, s_, pp in zip(names, beta, se, p)]


def partial_loglik(beta, X):
    eta = X @ beta; ll = 0.0
    for tt in np.unique(t[e == 1]):
        rm = t >= tt; dm = (t == tt) & (e == 1)
        ll += eta[dm].sum() - dm.sum() * np.log(np.exp(eta[rm]).sum())
    return ll


def lrt(X):
    b, _, _, _ = cox_fit(X, t, e, l2=0.0)
    return float(2 * (partial_loglik(b, X) - partial_loglik(np.zeros(X.shape[1]), X)))


models = {}
Zs = {c: zstd(df[c].values) for c in [SCORE, "grade", "stage", "translation_meas",
                                       "secretion_ER_morph", "secretion_ER_meas", "matrisome_morph"]}
specs = {
    "uni_score": [SCORE],
    "adj_grade_stage": [SCORE, "grade", "stage"],
    "adj_grade_only": [SCORE, "grade"],
    "adj_stage_only": [SCORE, "stage"],
    "uni_grade": ["grade"],
    "uni_stage": ["stage"],
    "grade_stage_only": ["grade", "stage"],
    "uni_translation_meas": ["translation_meas"],
    "uni_secretionER_morph": ["secretion_ER_morph"],
    "uni_secretionER_meas": ["secretion_ER_meas"],
    "uni_matrisome_morph": ["matrisome_morph"],
    "adj_secretionER_morph": ["secretion_ER_morph", "grade", "stage"],
    "adj_matrisome_morph": ["matrisome_morph", "grade", "stage"],
}
print("\n=== Cox PH (Breslow, z-standardised covariates; l2=1e-3 as in clinical_cox.py) ===")
for nm, cols in specs.items():
    X = np.column_stack([Zs[c] for c in cols])
    b, s_, p, cov = cox_fit(X, t, e, l2=1e-3)
    b0, s0, p0, _ = cox_fit(X, t, e, l2=0.0)
    models[nm] = dict(terms=summarise(cols, b, s_, p), terms_l2_0=summarise(cols, b0, s0, p0),
                      n=n, events=d, epv=d / len(cols), lrt_chi2=lrt(X), lrt_df=len(cols))
    models[nm]["lrt_p"] = float(1 - chi2.cdf(models[nm]["lrt_chi2"], len(cols)))
    print(f"[{nm}] n={n} events={d} EPV={d / len(cols):.1f} LRT chi2={models[nm]['lrt_chi2']:.2f} "
          f"p={models[nm]['lrt_p']:.4f}")
    for tm in models[nm]["terms"]:
        print(f"   {tm['term']:<22} HR={tm['HR']:.3f} 95%CI [{tm['CI95'][0]:.3f}, {tm['CI95'][1]:.3f}] "
              f"p={tm['p']:.4f}")

# per-unit (unstandardised) HR for grade & stage in the adjusted model
X = np.column_stack([Zs[SCORE], df["grade"].values, df["stage"].values])
b, s_, p, cov = cox_fit(X, t, e, l2=1e-3)
models["adj_grade_stage_perunit"] = dict(
    terms=summarise([SCORE + "(per SD)", "grade(per level)", "stage(per level)"], b, s_, p))
print("[adj per-unit grade/stage]",
      [(tm["term"], round(tm["HR"], 3), [round(x, 3) for x in tm["CI95"]], round(tm["p"], 4))
       for tm in models["adj_grade_stage_perunit"]["terms"]])
out["cox"] = models
out["cox_meta"] = dict(
    sd_grade=float(np.nanstd(df["grade"].values)), sd_stage=float(np.nanstd(df["stage"].values)),
    corr_score_grade_spearman=float(stats.spearmanr(df[SCORE], df["grade"])[0]),
    corr_score_stage_spearman=float(stats.spearmanr(df[SCORE], df["stage"])[0]),
    corr_grade_stage_spearman=float(stats.spearmanr(df["grade"], df["stage"])[0]),
)


# ---------------- Schoenfeld / Grambsch-Therneau PH test (cox.zph, survival 2.x algorithm) ----------------
def schoenfeld_test(X, names):
    b, s_, p, cov = cox_fit(X, t, e, l2=0.0)   # unpenalised fit for the test
    eta = X @ b; w = np.exp(eta - eta.max())
    dead_idx = np.where(e == 1)[0]
    dead_idx = dead_idx[np.argsort(t[dead_idx], kind="stable")]
    res = []; tk = []
    for i in dead_idx:
        rm = t >= t[i]
        xbar = (X[rm] * w[rm, None]).sum(0) / w[rm].sum()
        res.append(X[i] - xbar); tk.append(t[i])
    res = np.array(res); tk = np.array(tk); nd = len(tk)
    xs, ys = km_curve(t, e)
    km_left = np.array([ys[xs < tt][-1] if (xs < tt).any() else 1.0 for tt in tk])   # KM(t-)
    transforms = {"km": 1 - km_left, "rank": stats.rankdata(tk), "identity": tk, "log": np.log(tk)}
    result = {}
    for tr, g in transforms.items():
        xx = g - g.mean()
        r2 = res @ cov * nd                                 # scaled Schoenfeld residuals (minus beta)
        test = xx @ r2
        zstat = test ** 2 / (np.diag(cov) * nd * (xx ** 2).sum())
        per = {nm: dict(rho=float(np.corrcoef(xx, r2[:, j])[0, 1]), chi2=float(zstat[j]),
                        p=float(1 - chi2.cdf(zstat[j], 1)))
               for j, nm in enumerate(names)}
        tg = xx @ res
        gchi = float(tg @ cov @ tg * nd / (xx ** 2).sum())
        per["GLOBAL"] = dict(chi2=gchi, df=len(names), p=float(1 - chi2.cdf(gchi, len(names))))
        result[tr] = per
    return result


print("\n=== Schoenfeld PH test (Grambsch-Therneau; transforms km / rank / identity / log) ===")
ph = {}
for nm, cols in [("uni_score", [SCORE]), ("adj_grade_stage", [SCORE, "grade", "stage"])]:
    X = np.column_stack([Zs[c] for c in cols])
    ph[nm] = schoenfeld_test(X, cols)
    for tr in ("km", "rank", "log"):
        print(f"[{nm}] transform={tr}: " +
              "; ".join(f"{k}: chi2={v['chi2']:.2f} p={v['p']:.3f}" for k, v in ph[nm][tr].items()))
out["schoenfeld"] = ph


# ---------------- BH within the local CCRCC families ----------------
def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    ranked = p[o] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    outq = np.empty(m); outq[o] = np.minimum(q, 1)
    return outq


fam = []
for r in fam_rows:
    fam += [dict(score=r["score"], var="grade", p=r["grade_p"]),
            dict(score=r["score"], var="stage", p=r["stage_p"]),
            dict(score=r["score"], var="survival_logrank", p=r["surv_logrank_p"])]
famdf = pd.DataFrame(fam)
famdf["q_ccrcc_all18"] = bh(famdf["p"].values)
morph = famdf[famdf.score.str.endswith("_morph")].copy()
morph["q_ccrcc_morph9"] = bh(morph["p"].values)
famdf = famdf.merge(morph[["score", "var", "q_ccrcc_morph9"]], how="left")
print("\n=== BH within CCRCC computed families ===")
print(famdf.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
out["ccrcc_bh_local"] = famdf.to_dict(orient="records")

with open(f"{ROOT}/review/recalc/recalc_survival_raw.json", "w") as fh:
    json.dump(out, fh, indent=1)
print("\nwrote recalc_survival_raw.json")
