#!/usr/bin/env python3
"""Sex as a covariate: does it explain the morphology-residual link?

The cohorts were assembled by CPTAC and their sex composition is set by disease
epidemiology, not by design (77/26 male/female in CCRCC, 69/36 LUAD, 0/100 UCEC,
55/44 GBM, 71/66 PDAC). UCEC is anatomically sex-specific and is therefore excluded:
sex is constant there and cannot be adjusted for.

Two questions, both at the level of the per-patient residual pathway scores:

  Q1  Does either score differ by sex? (Mann-Whitney U, two-sided, plus a
      rank-biserial effect size.) A morphology score that tracks sex would make
      sex a candidate confounder of everything downstream.

  Q2  Does the morphology-to-measured-residual link survive adjusting for sex?
      Partial correlation of the H&E-only score with the measured score given sex,
      as a percentage of the unadjusted correlation -- the same frame as the
      pathologist-read composition control (Supplementary Table, composition_path),
      where "retained %" is r_partial / r_raw.

  Q3  Is morphology reading sex directly? If it were, Y-linked proteins -- whose
      abundance is close to a deterministic function of sex -- should be
      morphology-predictable more often than the cohort's own base rate, because
      morphology could recover what their own mRNA leaves simply by reading sex off
      the tissue. Tested against each cohort's base rate (exact binomial per cohort;
      pooled observed versus expected under the per-cohort rates). UCEC is the
      built-in negative control: no Y-linked protein is detected there at all.

SCOPE. This is a score-level control. The per-gene estimates cannot be re-run here:
that needs the protein and RNA matrices, which live on the compute cluster, not in
this repository. What it does establish is whether the patient-level quantity the
clinical claims rest on is confounded by sex.

Inputs (both already in the repository):
  figures/figdata/scores_<cohort>.csv        per-patient measured and morphology-only scores
  figures/figdata/gdc_clinical_<cohort>.csv  sex_at_birth per patient
Output:
  review/recalc/sex_sensitivity.csv
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FD = os.path.join(ROOT, "figures", "figdata")
COHORTS = ["ccrcc", "luad", "gbm", "pdac"]          # UCEC excluded: uniformly female
FAMILIES = ["translation", "secretion_ER", "matrisome"]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1))
        q[o[i]] = prev
    return q


def partial_r(x, y, z):
    """Pearson correlation of x and y after regressing both on z (an intercept is included)."""
    Z = np.column_stack([np.ones(len(z)), z])
    bx = np.linalg.lstsq(Z, x, rcond=None)[0]
    by = np.linalg.lstsq(Z, y, rcond=None)[0]
    return stats.pearsonr(x - Z @ bx, y - Z @ by)


rows = []
for c in COHORTS:
    sc = pd.read_csv(f"{FD}/scores_{c}.csv", index_col=0)
    cl = pd.read_csv(f"{FD}/gdc_clinical_{c}.csv").set_index("case")
    sex = cl["sex_at_birth"].reindex(sc.index)
    male = sex.astype(str).str.lower().str.startswith("male")
    keep = sex.astype(str).str.lower().isin(["male", "female"])
    for fam in FAMILIES:
        mcol, pcol = f"{fam}_meas", f"{fam}_morph"
        if mcol not in sc.columns or pcol not in sc.columns:
            continue
        ok = keep & sc[mcol].notna() & sc[pcol].notna()
        if ok.sum() < 20:
            continue
        m = male[ok].values.astype(float)
        meas, morph = sc.loc[ok, mcol].values, sc.loc[ok, pcol].values
        n_m, n_f = int(m.sum()), int((1 - m).sum())
        if min(n_m, n_f) < 5:
            continue
        rec = dict(cohort=c, family=fam, n=int(ok.sum()), n_male=n_m, n_female=n_f)
        # Q1: does each score differ by sex?
        for tag, v in (("meas", meas), ("morph", morph)):
            u = stats.mannwhitneyu(v[m == 1], v[m == 0], alternative="two-sided")
            rec[f"{tag}_p_sex"] = float(u.pvalue)
            # rank-biserial: +1 = male scores uniformly higher
            rec[f"{tag}_rank_biserial"] = float(2 * u.statistic / (n_m * n_f) - 1)
            rec[f"{tag}_median_male"] = float(np.median(v[m == 1]))
            rec[f"{tag}_median_female"] = float(np.median(v[m == 0]))
        # Q2: does the morphology-to-measured link survive adjusting for sex?
        raw = stats.pearsonr(morph, meas)
        par = partial_r(morph, meas, m)
        rec.update(r_raw=float(raw.statistic), p_raw=float(raw.pvalue),
                   r_partial_sex=float(par.statistic), p_partial_sex=float(par.pvalue),
                   retained_pct=100.0 * float(par.statistic) / float(raw.statistic))
        rows.append(rec)

d = pd.DataFrame(rows)
for tag in ("meas", "morph"):
    d[f"{tag}_q_sex"] = bh(d[f"{tag}_p_sex"].values)
d.to_csv(os.path.join(ROOT, "review", "recalc", "sex_sensitivity.csv"), index=False)

pd.set_option("display.width", 200)
print("Q1  Does the score differ by sex? (Mann-Whitney, two-sided; BH q over the "
      f"{len(d)} tests of each type)")
print(d[["cohort", "family", "n", "n_male", "n_female", "meas_p_sex", "meas_q_sex",
         "meas_rank_biserial", "morph_p_sex", "morph_q_sex", "morph_rank_biserial"]]
      .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("\nQ2  Morphology-to-measured link, unadjusted vs adjusted for sex")
print(d[["cohort", "family", "n", "r_raw", "r_partial_sex", "retained_pct"]]
      .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print(f"\nretained %: min {d.retained_pct.min():.1f}, median {d.retained_pct.median():.1f}, "
      f"max {d.retained_pct.max():.1f}")
print(f"smallest BH q for a sex difference: measured {d.meas_q_sex.min():.3f}, "
      f"morphology-only {d.morph_q_sex.min():.3f}")

# ---- Q3: are Y-linked proteins morphology-predictable above the cohort base rate? ----
YLINKED = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "USP9Y", "EIF1AY", "NLGN4Y", "TXLNGY", "ZFY", "TMSB4Y"]
print("\nQ3  Y-linked proteins versus each cohort's base rate")
obs_tot, exp_tot, yrows = 0, 0.0, []
for c in COHORTS + ["ucec"]:
    g = pd.read_csv(os.path.join(ROOT, "server_export", "pinned", c,
                                 "residual_results_tumoronly.csv")).drop_duplicates("gene")
    sig = (g.fdr < 0.05) & (g.incremental_r2 > 0)
    base = float(sig.mean())
    sub = g[g.gene.isin(YLINKED)]
    k = int(((sub.fdr < 0.05) & (sub.incremental_r2 > 0)).sum())
    p = (stats.binomtest(k, len(sub), base, alternative="greater").pvalue
         if len(sub) else float("nan"))
    hits = sorted(sub.gene[(sub.fdr < 0.05) & (sub.incremental_r2 > 0)])
    print(f"   {c:6s} base {100*base:5.1f}%  Y-linked tested {len(sub)}  significant {k} "
          f"{hits}  exact binomial p {p:.3f}" if len(sub) else
          f"   {c:6s} base {100*base:5.1f}%  Y-linked tested 0 (negative control: none detected)")
    yrows.append(dict(cohort=c, base_rate=base, y_tested=len(sub), y_significant=k, binom_p=p))
    if c != "ucec":
        obs_tot += k
        exp_tot += len(sub) * base
print(f"   pooled over the four mixed-sex cohorts: {obs_tot} significant of "
      f"{sum(r['y_tested'] for r in yrows if r['cohort'] != 'ucec')} Y-linked proteins, "
      f"{exp_tot:.2f} expected at the cohorts' own base rates")
pd.DataFrame(yrows).to_csv(os.path.join(ROOT, "review", "recalc", "sex_sensitivity_ylinked.csv"), index=False)

print("\nUCEC excluded from Q1/Q2 (uniformly female). Per-gene increments are not re-estimated "
      "here; that needs the protein/RNA matrices on the cluster.")
