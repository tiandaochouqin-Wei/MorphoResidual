#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig. 6 (v7) "confirmation and boundaries" -- v6 (built to
review/QUALITY_REVIEW_2026-11-20/fig6_v6/spec.md) plus, in panel a, the post hoc rows of the
covariate re-reads (pre-specification 2026-11-21, section 7.9): one row for each of the four primary cells and one
row per cohort for the placebo margins. Panels c and d are v6's, shifted down by DY.

Reads only figdata/fig6v6_*.csv and figdata/fig6v7_a_covariates.csv. Every number printed on the figure is
formatted from a CSV value and first compared with the published string
(asserts below). No random numbers. Canvas 183 mm wide, drawn at final size.

Layout: all text and marks are placed in millimetres on one full-canvas
overlay axes (origin top left); only panel b (ECDF) and the tick axes use
real data axes.
"""
import sys, io, itertools
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import mrstyle as S

S.use_v2()
mpl.rcParams.update({"font.family": ["Arial", "DejaVu Sans"]})

INK2, GREY_D, GREY_M, LGREY, ORG = S.INK2, S.GREY_D, S.GREY_M, S.LGREY, S.ORG
DD = "figdata/fig6v6_"
W = print
MINUS = "−"


def fm(v, nd):
    """format with a true minus sign"""
    s = f"{v:.{nd}f}"
    return s.replace("-", MINUS)


# ---------------------------------------------------------------- data
ta = pd.read_csv(DD + "a_tests.csv").set_index("analysis_id")
tb = pd.read_csv(DD + "a_bootstrap.csv").set_index("analysis_id")
bg = pd.read_csv(DD + "b_genes.csv")
bs = pd.read_csv(DD + "b_summary.csv")
rg = pd.read_csv(DD + "c_rppa_groups.csv")
rs = pd.read_csv(DD + "c_rppa_stats.csv")
bgene = pd.read_csv(DD + "c_brca_gene.csv").set_index("section")
bset = pd.read_csv(DD + "c_brca_set.csv").set_index("section")
dth = pd.read_csv(DD + "d_theta.csv")
drd = pd.read_csv(DD + "d_reading.csv").set_index("cohort")
cv = pd.read_csv("figdata/fig6v7_a_covariates.csv")


def summ(cohort, quantity):
    r = bs[(bs.cohort == cohort) & (bs.quantity == quantity)]
    assert len(r) == 1, (cohort, quantity, len(r))
    return r.iloc[0]


def commas(n):
    return f"{int(round(n)):,}"


# =============================================================== PART 1 a
W("=" * 70, "\nPART 1 -- panel a asserts\n", "=" * 70, sep="")
# analysis_id -> (status printed, label, AUC, null, p, dagger)
A_ROWS = {
    "U": [("U_unrestricted", "pre-specified", "unrestricted", "0.597", "0.498 ± 0.036", "0.004"),
          ("U_within_operator", "post hoc", "within operator", "0.597", "0.506 ± 0.032", "0.004"),
          ("U_operator_residualised", "post hoc", "within operator, residualised§", "0.567", "0.489 ± 0.031", "0.006"),
          ("U_compl_discovery", "post hoc", "completeness-stratified", "0.561", "0.467 ± 0.037", "0.007"),
          ("U_within_plex", "post hoc", "within TMT plex", "0.597", "0.503 ± 0.036", "0.004"),
          ("U_plex_residualised", "post hoc", "within TMT plex, residualised§", "0.608", "0.545 ± 0.033", "0.026")],
    "L": [("L_step1_unrestricted", "registered", "step 1, unrestricted", "0.644", "0.524 ± 0.019", "0.001†"),
          ("L_step2_within_operator", "registered", "step 2, within operator", "0.644", "0.547 ± 0.018", "0.001†"),
          ("L_step3_operator_residualised", "registered", "step 3, 726 set, within operator‡", "0.654", "0.527 ± 0.026", "0.001†"),
          ("L_compl_discovery", "registered", "completeness-stratified§", "0.587", "0.470 ± 0.017", "0.001†"),
          ("L_within_plex", "post hoc", "within TMT plex", "0.644", "0.546 ± 0.019", "0.001†"),
          ("L_plex_residualised", "post hoc", "726 set, within plex, residualised§", "0.651", "0.519 ± 0.024", "0.001†")],
}
STATUS_OK = {"pre-specified": "pre-specified", "post hoc": "post hoc",
             "registered": "registered", "registered (descriptive item)": "registered"}
RNG_EXPECT = {"U_unrestricted": "0.48–0.65", "L_step1_unrestricted": "0.55–0.65"}
A_OUT = {}
for k, rows in A_ROWS.items():
    out = []
    for aid, stat, lab, auc, nul, p in rows:
        r = ta.loc[aid]
        a_s = f"{r.obs_auc:.3f}"
        n_s = f"{r.null_mean:.3f} ± {r.null_sd:.3f}"
        floor = abs(r.p_one_sided - 1 / 1001) < 1e-9
        p_s = f"{r.p_one_sided:.3f}" + ("†" if floor else "")
        assert (a_s, n_s, p_s) == (auc, nul, p), (aid, (a_s, n_s, p_s), (auc, nul, p))
        assert STATUS_OK[r.status] == stat, (aid, r.status)
        assert int(r.n_perm) == 1000
        if "†" in p:
            assert int(r.n_exceed) == 0
        if aid in RNG_EXPECT:
            b = tb.loc[aid]
            rng = f"{b.ci95_lo:.2f}–{b.ci95_hi:.2f}"
            assert rng == RNG_EXPECT[aid], (aid, rng)
            lab = lab.format(rng=rng)
        out.append((aid, stat, lab, r.obs_auc, r.null_mean, r.null_sd, a_s, n_s, p_s))
    A_OUT[k] = out
assert int(ta.loc["U_unrestricted", "cohort_n"]) == 138 and int(ta.loc["L_step1_unrestricted", "cohort_n"]) == 112
assert int(ta.loc["U_unrestricted", "n_selected"]) == 2536 and int(ta.loc["L_step1_unrestricted", "n_selected"]) == 2287
assert int(ta.loc["L_step3_operator_residualised", "n_selected"]) == 721
assert int(ta.loc["U_unrestricted", "n_background"]) == 7610 and int(ta.loc["L_step1_unrestricted", "n_background"]) == 8119
# ---- post hoc covariate rows (pre-specification 2026-11-21; results record 2026-11-23): four primary cells + placebos
COV_EXPECT = {("U", "C3"): ("0.551", "0.515 ± 0.035", "0.161"), ("U", "R20"): ("0.546", "0.508 ± 0.041", "0.199"),
              ("L", "C3"): ("0.623", "0.556 ± 0.016", "0.001†"), ("L", "R20"): ("0.598", "0.564 ± 0.017", "0.029")}
COV_LAB = {"C3": "+ 3 composition scores", "R20": "+ 20 expression components"}
PLC = {}
for k, coh in (("U", "ucec"), ("L", "luad")):
    for cov in ("C3", "R20"):
        r = cv[(cv.cohort == coh) & (cv.covariate == cov)]
        assert len(r) == 1, (coh, cov)
        r = r.iloc[0]
        a_s = f"{r.obs_auc:.3f}"
        n_s = f"{r.null_mean:.3f} ± {r.null_sd:.3f}"
        floor = abs(r.p_one_sided - 1 / 1001) < 1e-9
        p_s = f"{r.p_one_sided:.3f}" + ("†" if floor else "")
        assert (a_s, n_s, p_s) == COV_EXPECT[(k, cov)], (coh, cov, (a_s, n_s, p_s))
        assert int(r.n_perm) == 1000
        if floor:
            assert int(r.n_exceed) == 0
        A_OUT[k].append((f"{k}_cov_{cov}", "post hoc", COV_LAB[cov], r.obs_auc, r.null_mean, r.null_sd, a_s, n_s, p_s))
        if cov == "R20":
            mj = np.array([float(x) for x in r.m_j.split(";")])
            assert len(mj) == 119 and int(r.e) == int((mj <= float(r.m_R)).sum())
            assert (r.null_mean + mj).max() < 0.70 and (r.null_mean + mj).min() > 0.40   # all on the 0.40-0.70 axis
            PLC[k] = dict(mj=mj, nm=float(r.null_mean), mR=float(r.m_R), e=int(r.e), p_abs=float(r.p_abs))
assert (PLC["U"]["e"], f"{PLC['U']['p_abs']:.3f}", PLC["L"]["e"], f"{PLC['L']['p_abs']:.3f}") == (46, "0.392", 0, "0.008")
assert all(len(A_OUT[k]) == 8 for k in A_OUT)
W("[assert OK] panel a: four post hoc covariate rows, two placebo rows")
# bootstrap brackets
b3 = tb.loc["L_step3_operator_residualised"]
assert (f"{b3.ci95_lo:.3f}", f"{b3.ci95_hi:.3f}") == ("0.536", "0.660")
assert all(int(tb.loc[i, "boots"]) == 1000 for i in tb.index)
# footnote numbers
fn_lo = ta.loc["U_compl_joint", "obs_auc"]; fn_hi = ta.loc["U_compl_discovery", "obs_auc"]
fn_p = ta.loc["U_compl_joint", "p_one_sided"]
assert (f"{fn_lo:.3f}", f"{fn_hi:.3f}", f"{fn_p:.3f}") == ("0.554", "0.561", "0.010")
others = ["U_compl_discovery", "U_compl_confirmatory", "U_compl_joint", "U_compl_confirmatory_exact"]
assert max(ta.loc[o, "p_one_sided"] for o in others) <= 0.0100 + 1e-4
assert min(ta.loc[o, "obs_auc"] for o in others) == fn_lo
assert max(ta.loc[o, "obs_auc"] for o in others) == fn_hi   # all four strata, as main.tex (0.554--0.561)
assert len(others) == 4
n_sel_u = summ("UCEC", "discovery-selected genes (uncorrected set)").stored_value
n_sel_l = summ("LUAD", "discovery-selected genes (uncorrected set)").stored_value
assert (int(n_sel_u), int(n_sel_l)) == (2566, 2310)
assert int(summ("UCEC", "discovery-selected genes testable in confirmatory").stored_value) == 2536
assert float(bg[(bg.cohort == "LUAD")].confirmatory_increment.max()) <= 0.45  # no hidden right truncation in b
assert int(summ("LUAD", "discovery-selected genes testable in confirmatory").stored_value) == 2287
GT_U = (f"UCEC (uterus), n = {int(ta.loc['U_unrestricted','cohort_n'])}; {commas(ta.loc['U_unrestricted','n_selected'])} of "
        f"{commas(n_sel_u)} selected genes testable, {commas(ta.loc['U_unrestricted','n_background'])} background")
GT_L = (f"LUAD (lung), n = {int(ta.loc['L_step1_unrestricted','cohort_n'])}; {commas(ta.loc['L_step1_unrestricted','n_selected'])} of "
        f"{commas(n_sel_l)} selected genes testable, {commas(ta.loc['L_step1_unrestricted','n_background'])} background")
W("[assert OK] panel a: 12 published rows (+4 post hoc covariate rows, +2 placebo rows), brackets, footnote numbers")

# =============================================================== PART 1 b
W("PART 1 -- panel b asserts")
B = {}
for c, nsel, nbg, tot, neg, rs_, rb_ in (("UCEC", 2536, 7610, 10146, "93.2", -0.053, 0.331),
                                       ("LUAD", 2287, 8119, 10406, "85.8", 0.070, 0.419)):
    d = bg[bg.cohort == c]
    sel = d[d.group == "selected"].confirmatory_increment.values
    bgd = d[d.group == "background"].confirmatory_increment.values
    assert len(sel) == nsel and len(bgd) == nbg and len(d) == tot
    assert int(summ(c, "genes in both cohorts (usable)" if c == "UCEC" else "testable genes (usable universe)").stored_value) == tot
    share = 100 * d.confirmatory_negative.mean()
    pub = summ(c, "share of genes with negative confirmatory increment" if c == "UCEC" else "share of genes with negative confirmatory increment")
    assert f"{share:.1f}" == neg == f"{pub.published_value:.1f}", (c, share, pub.published_value)
    r_sel = summ(c, "Spearman rho within selected genes"); r_bg = summ(c, "Spearman rho within background genes")
    assert (f"{r_sel.stored_value:.3f}", f"{r_bg.stored_value:.3f}") == (f"{rs_:.3f}", f"{rb_:.3f}"), c
    assert round(float(r_sel.published_value), 3) == round(rs_, 3) and round(float(r_bg.published_value), 3) == round(rb_, 3)
    B[c] = dict(sel=np.sort(sel), bg=np.sort(bgd), neg=neg, tot=tot, rs=float(r_sel.stored_value),
                rb=float(r_bg.stored_value), nsel=nsel, nbg=nbg)
W("[assert OK] panel b: group sizes, negative shares, Spearman rho")

# =============================================================== PART 1 c
W("PART 1 -- panel c asserts")
RP = {  # cohort -> (left label, line1, line2)
    "KIRC": ("TCGA-KIRC (n = 451)", "35/62 (56.5%) vs 47.0%", "OR 1.46, p = 0.201"),
    "LUAD": ("TCGA-LUAD (n = 328)", "45/79 (57.0%) vs 55.6%", "OR 1.06, p = 0.896"),
    "GBM": ("TCGA-GBM (n = 85)", "0/23 (0.0%) vs 0.0%", "OR not defined, p = 1.000"),
    "PDAC": ("TCGA-PDAC (n = 112)", "7/42 (16.7%) vs 15.7%", "OR 1.08, p = 0.822"),
}
WILSON_PUB = {}  # published Wilson CIs checked against SuppTable_replication copy
pc = pd.read_csv(DD + "c_published_check.csv")
RPPA = {}
for coh, (lab, l1, l2) in RP.items():
    g = rg[rg.cohort == coh]
    sig = g[g.group.str.startswith("CPTAC-significant")].iloc[0]
    oth = g[g.group.str.startswith("other tested")].iloc[0]
    st = rs[(rs.cohort == coh) & (rs.block == "rppa")].iloc[0]
    n_s = f"{int(sig.n_predictable)}/{int(sig.n_proteins)} ({100*sig.share:.1f}%) vs {100*oth.share:.1f}%"
    if np.isnan(st.odds_ratio):
        o_s = f"OR not defined, p = {st.p_two_sided:.3f}"
    else:
        o_s = f"OR {st.odds_ratio:.2f}, p = {st.p_two_sided:.3f}"
    lab_s = f"TCGA-{coh} (n = {int(st.n_cases)})"
    assert (lab_s, n_s, o_s) == (lab, l1, l2), (coh, (lab_s, n_s, o_s))
    RPPA[coh] = dict(sig=sig, oth=oth, lab=lab, l1=n_s, l2=o_s, n=int(st.n_cases))
pool = rs[rs.block == "rppa_pooled"].iloc[0]
POOL = f"pooled CMH OR {pool.odds_ratio:.2f} (95% CI {pool.or_ci_lo:.2f}–{pool.or_ci_hi:.2f})"
assert POOL == "pooled CMH OR 1.20 (95% CI 0.85–1.70)"
# pooled CMH p: stored 0.3447; main.tex l.689 and SuppTable_replication.tex print 0.34 (corrected 2026-11-20)
POOL_P = f"p = {pool.p_two_sided:.2f}"
assert POOL_P == "p = 0.34" and f"{pool.p_two_sided:.4f}" == "0.3447", pool.p_two_sided
# Wilson interval agreement with the check table (KIRC CPTAC-significant 44.1-68.1)
k_sig = RPPA["KIRC"]["sig"]
assert (f"{100*k_sig.ci_lo:.1f}", f"{100*k_sig.ci_hi:.1f}") == ("44.1", "68.1")
_mm = pc[pc.status.str.lower() != "match"]
assert _mm.empty, list(_mm.quantity)  # every published number of panel c matches, pooled p included
_pp = pc[pc.quantity == "pooled CMH continuity-corrected p"].iloc[0]
assert float(_pp.published) == 0.34 and _pp.status == "match"
# BRCA
assert all(int(bgene.loc[s, "n_patients"]) == 102 and int(bgene.loc[s, "n_tested"]) == 10031 and
           int(bgene.loc[s, "n_significant"]) == 0 for s in ("DX", "TS"))
dx_neg = f"{100*bgene.loc['DX','share_negative']:.1f}"; ts_neg = f"{100*bgene.loc['TS','share_negative']:.1f}"
assert (dx_neg, ts_neg) == ("99.1", "98.9")
assert all(int(bset.loc[s, "n_set_file"]) == 777 and int(bset.loc[s, "n_selected_tested"]) == 766 and
           int(bset.loc[s, "n_background"]) == 7307 and int(bset.loc[s, "n_perm"]) == 1000 for s in ("DX", "TS"))
BRCA = []
for s, expect in (("DX", ("0.517", "0.498 ± 0.027", "0.234")), ("TS", ("0.473", "0.498 ± 0.027", "0.819"))):
    r = bset.loc[s]
    got = (f"{r.auc_obs:.3f}", f"{r.auc_null_mean:.3f} ± {r.auc_null_sd:.3f}", f"{r.p_one_sided:.3f}")
    assert got == expect, (s, got)
    BRCA.append((s, r.auc_obs, r.auc_null_mean, r.auc_null_sd, *got))
W("[assert OK] panel c: RPPA rows, pooled, BRCA gene/set")

# =============================================================== PART 1 d
W("PART 1 -- panel d asserts")
D_ORDER = ["CCRCC", "LUAD", "UCEC", "GBM", "PDAC"]
D_EXPECT = {  # cohort: (n, lambda, status, class, dilution minimum, label2, reading)
    "CCRCC": (1946, "0.69", "measured", "R3", "0.79"),
    "LUAD": (2205, "0.17", "transferred", "R2", "0.11"),
    "UCEC": (2418, "0.17", "measured", "R2", "0.11"),
    "GBM": (684, "0.17", "transferred", "R2", "0.08"),
    "PDAC": (1454, "0.17", "transferred", "R2", "0.11"),
}
D = {}
for c in D_ORDER:
    r = drd.loc[c]
    t = dth[(dth.cohort == c) & (dth.element == "theta")].iloc[0]
    m = dth[(dth.cohort == c) & (dth.element == "dilution_minimum")].iloc[0]
    n, lam, stt, cls, dmin = D_EXPECT[c]
    assert int(r.n_genes) == n and f"{r.lambda_lo:.2f}" == lam and r.lambda_lo_status == stt
    assert r.reading_class == cls and f"{m.value:.2f}" == dmin
    # theta and interval equal the published values at the published precision
    assert round(t.value, 2) == round(t.published_value, 2), c
    assert round(t.ci_lo, 2) == round(t.published_ci_lo, 2) and round(t.ci_hi, 2) == round(t.published_ci_hi, 2), c
    assert t.value == r.theta and t.ci_lo == r.ci_lo and t.ci_hi == r.ci_hi
    stated = bool(r.dilution_bound_stated_in_manuscript)
    assert stated == (c == "CCRCC")
    if stated:
        assert f"{r.dilution_bound_printed:.2f}" == "0.08"
    D[c] = dict(n=n, lam=lam, stt=stt, theta=t.value, lo=t.ci_lo, hi=t.ci_hi, dmin=m.value,
                reading=("excluded up to a bound" if cls == "R3" else "not excluded"), bound=stated)
W("[assert OK] panel d: theta, intervals, minima, reliabilities, readings")

# =============================================================== PART 3 figure
MM = 25.4
FWmm = 183.0
DY = 12.4               # panels c and d are v6's, shifted down by DY for the post hoc rows of panel a
FHmm = 151.5 + DY       # spec section 1 allows at most 165 mm
assert FHmm <= 165.0
fig = plt.figure(figsize=(FWmm / MM, FHmm / MM))
FS, FS_T, FS_L = 5.5, 6.0, 8.0
assert min(FS, FS_T) >= 5.0

ov = fig.add_axes([0, 0, 1, 1], zorder=5)
ov.set_xlim(0, FWmm); ov.set_ylim(FHmm, 0); ov.axis("off"); ov.patch.set_alpha(0)
TXT = []
YOFF = 0.0              # vertical offset of everything drawn through the helpers below (0 in a and b, DY in c and d)


def T(x, y, s, ha="left", va="center", fs=FS, bold=False, color=INK2, rot=0):
    t = ov.text(x, y + YOFF, s, ha=ha, va=va, fontsize=fs, fontweight="bold" if bold else "normal",
                color=color, rotation=rot)
    TXT.append(t)
    return t


# panel letter and panel title share one top line (va="top" at the same y), so the
# letter's top sits on the title's cap/ascender line
def letter(ch, x, y):
    fig.text(x / FWmm, 1 - (y + YOFF + 2.1) / FHmm, ch, fontsize=FS_L, fontweight="bold", va="baseline", ha="left", color="black")


def ptitle(x, y, s):
    T(x, y + 2.1, s, va="baseline", fs=FS_T, bold=False)


def tick_axis(x0, x1, ybot, lim, ticks, labels=None, log=False):
    ax = fig.add_axes([x0 / FWmm, 1 - (ybot + YOFF) / FHmm, (x1 - x0) / FWmm, 0.4 / FHmm])
    if log:                      # panel d only; major ticks only, as on the linear axes
        ax.set_xscale("log")
        ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_xlim(*lim); ax.set_xticks(ticks)
    ax.set_xticklabels(labels if labels else [f"{t:g}" for t in ticks], color=INK2)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.patch.set_alpha(0)
    ax.tick_params(axis="x", colors=INK2, labelsize=FS)
    return ax


def mapper(x0, x1, lo, hi):
    return lambda v: x0 + (v - lo) / (hi - lo) * (x1 - x0)


def logmapper(x0, x1, lo, hi):
    assert lo > 0
    return lambda v: x0 + (np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo)) * (x1 - x0)


def dot(x, y, color, ms=3.0, z=6, mec="white", mfc=None):
    ov.plot([x], [y + YOFF], "o", ms=ms, mfc=color if mfc is None else mfc, mec=mec, mew=0.3 if mfc is None else 0.7,
            zorder=z, clip_on=False)


def hline(xa, xb, y, color, lw=0.6, z=4):
    ov.plot([xa, xb], [y + YOFF, y + YOFF], color=color, lw=lw, solid_capstyle="butt", zorder=z, clip_on=False)


def vtick(x, ya, yb, color, lw=0.6, z=5):
    ov.plot([x, x], [ya + YOFF, yb + YOFF], color=color, lw=lw, solid_capstyle="butt", zorder=z, clip_on=False)


def null_bar(y, nm, sd):
    ov.add_patch(Rectangle((mA(nm - sd), y + YOFF - 0.75), mA(nm + sd) - mA(nm - sd), 1.5, fc=LGREY, ec="none", zorder=3))
    vtick(mA(nm), y - 0.95, y + 0.95, GREY_D)


def stat_headers(yb):
    """AUC / null mean ± s.d. / one-sided p, all bottom lines on one baseline yb"""
    T(X_AUC, yb, "AUC", ha="right", va="baseline")
    T(X_NUL, yb, "null", ha="right", va="baseline")
    T(X_P, yb - 2.3, "one-sided", ha="right", va="baseline"); T(X_P, yb, "p", ha="right", va="baseline")


def p_cell(y, p_s):
    """decimal-aligned p: digits right-aligned at X_P, dagger hangs outside the column"""
    T(X_P, y, p_s.rstrip("†"), ha="right")
    if p_s.endswith("†"):
        T(X_P + 0.15, y, "†", ha="left")


# ----------------------------------------------------------------- columns
AX0, AX1 = 49.0, 77.0            # a / c-BRCA AUC axis (shared x mapping)
mA = mapper(AX0, AX1, 0.40, 0.70)
X_ST, X_LAB = 1.0, 14.0
X_AUC, X_NUL, X_P = 82.8, 96.0, 106.0
LET_X = 1.0                      # left-column letters inset like the content

# ================================================================= panel a
letter("a", LET_X, 1.0)
ptitle(5.0, 1.0, "Set-level ordering in two later CPTAC-3 cohorts")
stat_headers(8.2)
PITCH = 2.9


def a_block(key, y_title, y0, title, color, brackets):
    T(X_ST, y_title, title)
    for i, (aid, stat, lab, obs, nm, sd, a_s, n_s, p_s) in enumerate(A_OUT[key]):
        y = y0 + PITCH * i
        T(X_ST, y, stat); T(X_LAB, y, lab)
        T(X_AUC, y, a_s, ha="right"); T(X_NUL, y, n_s, ha="right"); p_cell(y, p_s)
        null_bar(y, nm, sd)
        dot(mA(obs), y, color)
        if aid in brackets:
            b = tb.loc[aid]
            ya = y                  # on the row's own line, through the dot (as in c and d)
            bc = color
            hline(mA(b.ci95_lo), mA(b.ci95_hi), ya, bc, lw=0.6, z=4)
            vtick(mA(b.ci95_lo), ya - 0.55, ya + 0.55, bc, lw=0.6); vtick(mA(b.ci95_hi), ya - 0.55, ya + 0.55, bc, lw=0.6)
    # placebo row (post hoc): 119 random reassignments of the 20 components; each grey dot is the margin of one placebo
    # placed at the real cell's null mean plus that margin; the coloured dot is the real components' margin (columns 0, 1-25)
    P_ = PLC[key]
    y = y0 + PITCH * len(A_OUT[key])
    T(X_ST, y, "post hoc"); T(X_LAB, y, f"placebos: {P_['e']} of 119 ≤ real")
    for m in P_["mj"]:
        ov.plot([mA(P_["nm"] + m)], [y + YOFF], "o", ms=1.3, mfc=GREY_M, mec="none", alpha=0.55, zorder=4, clip_on=False)
    dot(mA(P_["nm"] + P_["mR"]), y, color)
    T(X_P, y, f"{P_['p_abs']:.3f}", ha="right"); T(X_P + 0.15, y, "‖", ha="left")



BR = {"U_unrestricted", "L_step1_unrestricted", "L_step3_operator_residualised"}
a_block("U", 10.6, 13.4, GT_U, ORG["UCEC"], BR)
a_block("L", 40.0, 42.8, GT_L, ORG["LUAD"], BR)
tick_axis(AX0, AX1, 68.4, (0.40, 0.70), [0.40, 0.50, 0.60, 0.70], ["0.40", "0.50", "0.60", "0.70"])
T((AX0 + AX1) / 2, 72.6, "set-level AUC, selected vs background genes", ha="center")
KY = 75.4   # key for panel a, below the axis
dot(2.6, KY, INK2)
T(4.4, KY, "observed AUC (cohort colour)")
ov.add_patch(Rectangle((34.0, KY - 0.75), 3.6, 1.5, fc=LGREY, ec="none", zorder=3)); vtick(35.8, KY - 0.95, KY + 0.95, GREY_D)
T(39.0, KY, "null mean ± s.d.")
_kb = INK2
vtick(57.6, KY - 0.55, KY + 0.55, _kb, lw=0.6); hline(57.6, 61.9, KY, _kb, lw=0.6); vtick(61.9, KY - 0.55, KY + 0.55, _kb, lw=0.6)
dot(59.75, KY, INK2)
T(63.0, KY, "95% patient-bootstrap interval, 3 rows only")
T(1, 78.1, "† floor at 1,000 permutations.  ‡ operator-batch-corrected set, 721 of 726 testable.")
T(1, 80.5, f"§ descriptive item, not read.  All four UCEC completeness strata: AUC {fn_lo:.3f}–{fn_hi:.3f}, p ≤ {fn_p:.3f}.")
T(1, 82.9, "‖ randomisation p of the capacity comparison, (1 + e)/120; grey dots, 119 reassignments of the 20 components.")

# ================================================================= panel b
letter("b", 110, 1.0)
ptitle(114, 1.0, "Per-gene increment in the same cohorts")
BX0, BX1 = 123.0, 181.0
XL = (-0.8, 0.45)
axes_b = {}
for c, (yt, yb) in (("UCEC", (6.5, 35.0)), ("LUAD", (39.9, 68.4))):
    ax = fig.add_axes([BX0 / FWmm, 1 - yb / FHmm, (BX1 - BX0) / FWmm, (yb - yt) / FHmm])
    d = B[c]
    for arr, col in ((d["bg"], GREY_D), (d["sel"], ORG[c])):
        y = np.arange(1, len(arr) + 1) / len(arr)
        ax.step(np.r_[arr[0], arr], np.r_[0, y], where="post", color=col, lw=0.9, solid_joinstyle="round")
    ax.axvline(0, color=GREY_M, lw=0.4, ls=(0, (3, 2)), zorder=0)
    # where each curve crosses x = 0 (share of that group below zero), as data on the dashed line
    for arr, col in ((d["bg"], GREY_D), (d["sel"], ORG[c])):
        f0 = np.searchsorted(arr, 0.0, side="left") / len(arr)
        ax.plot([0], [f0], "o", ms=2.6, color=col, mec="white", mew=0.3, zorder=6, clip_on=False)
    ax.set_xlim(*XL); ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(["0", "0.5", "1"], color=INK2)
    ax.set_xticks([-0.8, -0.4, 0.0, 0.4])
    ax.tick_params(axis="both", colors=INK2, labelsize=FS)
    ax.patch.set_alpha(0)
    axes_b[c] = ax
    if c == "UCEC":
        ax.set_xticklabels([])
    else:
        ax.set_xticklabels(["≤ " + MINUS + "0.8", MINUS + "0.4", "0", "0.4"], color=INK2)
    x_t = BX0 + 1.2
    T(x_t, yt + 1.6, c, bold=True)
    T(x_t, yt + 4.2, f"all genes ({commas(d['tot'])}): {d['neg']}% < 0")
    T(x_t, yt + 6.8, "Spearman ρ with discovery")
    T(x_t, yt + 9.2, f"increment: selected {fm(d['rs'], 3)},")
    T(x_t, yt + 11.6, f"background {fm(d['rb'], 3)}")
    # direct curve labels in the empty lower-right region, two lines each
    xs = BX0 + 40.5
    for j, (nm, col, n_) in enumerate((("selected", ORG[c], d["nsel"]), ("background", GREY_D, d["nbg"]))):
        y1 = yb - 11.4 + 5.8 * j
        hline(xs, xs + 2.6, y1, col, lw=0.9)
        T(xs + 3.6, y1, nm)
        T(xs + 3.6, y1 + 2.4, f"({commas(n_)})")
    if c == "UCEC":
        dot(xs + 0.5, yb - 18.0, ORG[c], ms=2.6); dot(xs + 2.3, yb - 18.0, GREY_D, ms=2.6)
        T(xs + 4.2, yb - 19.2, "dot: share of"); T(xs + 4.2, yb - 16.8, "each group < 0")
T(116.3, 37.45, "cumulative fraction of genes", ha="center", rot=90)
T((BX0 + BX1) / 2, 73.4, "confirmatory incremental R² (morphology over own mRNA)", ha="center")

# ================================================================= panel c
YOFF = DY
letter("c", LET_X, 73.3)
ptitle(5.0, 73.3, "Outside CPTAC-3")
T(1, 78.2, "TCGA, RPPA (Phikon): CPTAC-significant vs other tested proteins")
RAX0 = AX0
mC = mapper(RAX0, AX1, -3, 75)
XSTAT = 81.0
T(XSTAT, 81.5, "two-sided Fisher", va="baseline")
rows_y = [86.5, 94.0, 101.5, 109.0]
dot(3.0, 81.6, INK2)
T(4.8, 81.6, "CPTAC-significant (cohort colour)")
dot(41.0, 81.6, GREY_M)
T(42.8, 81.6, "other tested (grey)")
OFF = 1.3
for (coh, d), yc in zip(RPPA.items(), rows_y):
    col = ORG["CCRCC" if coh == "KIRC" else coh]
    # left-aligned on the panel a / BRCA grid: cohort in the X_ST column, n in the X_LAB column
    assert d["lab"] == f"TCGA-{coh} (n = {d['n']})"
    T(X_ST, yc, f"TCGA-{coh}"); T(X_LAB, yc, f"n = {d['n']}")
    T(XSTAT, yc - 1.15, d["l1"]); T(XSTAT, yc + 1.15, d["l2"])
    for g, cc, off in ((d["sig"], col, -OFF), (d["oth"], GREY_M, OFF)):
        y = yc + off
        hline(mC(100 * g.ci_lo), mC(100 * g.ci_hi), y, cc, lw=0.6)
        dot(mC(100 * g.share), y, cc)
# GBM other-tested: Wilson upper bound is 1.4%, narrower than the marker (checked, drawn above)
_gb = RPPA["GBM"]["oth"]
assert abs(100 * _gb.ci_hi - 1.4) < 0.05 and _gb.ci_lo == 0
# pooled row: its own text line below the axis title, not on the tick-label line
tick_axis(RAX0, AX1, 113.0, (-3, 75), [0, 25, 50, 75])
T((RAX0 + AX1) / 2, 117.6, "morphology-predictable proteins (%), Wilson 95% CI", ha="center")
_p1 = f"pooled CMH OR {pool.odds_ratio:.2f}"
_p2 = f"(95% CI {pool.or_ci_lo:.2f}–{pool.or_ci_hi:.2f})"
assert f"{_p1} {_p2}" == POOL
T(X_ST, 121.6, f"Four cohorts pooled: CMH OR {pool.odds_ratio:.2f} {_p2}, {POOL_P}")
# BRCA
T(1, 126.8, f"CPTAC-2 TCGA-BRCA, iTRAQ proteome, 102 patients; gene level: {int(bgene.loc['DX','n_significant'])} of {commas(bgene.loc['DX','n_tested'])} proteins at FDR < 0.05")
T(1, 129.5, f"set: {commas(bset.loc['DX','n_set_file'])} genes significant in ≥ 3 CPTAC-3 cohorts, "
            f"{commas(bset.loc['DX','n_selected_tested'])} tested, vs {commas(bset.loc['DX','n_background'])} background")
dot(3.0, 133.4, INK2); T(4.8, 133.4, "BRCA observed AUC")
ov.add_patch(Rectangle((30.0, 133.4 + YOFF - 0.75), 3.6, 1.5, fc=LGREY, ec="none", zorder=3)); vtick(31.8, 133.4 - 0.95, 133.4 + 0.95, GREY_D)
T(35.0, 133.4, "null mean ± s.d.")
stat_headers(134.6)
for (s_, obs, nm, sd, a_s, n_s, p_s), y in zip(BRCA, (137.8, 141.2)):
    T(X_ST, y, "post hoc"); T(X_LAB, y, "diagnostic slides" if s_ == "DX" else "frozen sections")
    T(X_AUC, y, a_s, ha="right"); T(X_NUL, y, n_s, ha="right"); p_cell(y, p_s)
    null_bar(y, nm, sd)
    dot(mA(obs), y, INK2)
tick_axis(AX0, AX1, 144.0, (0.40, 0.70), [0.40, 0.50, 0.60, 0.70], ["0.40", "0.50", "0.60", "0.70"])
T((AX0 + AX1) / 2, 148.7, "set-level AUC, selected vs background genes", ha="center")

# ================================================================= panel d
letter("d", 110, 73.3)
ptitle(114, 73.3, "Own-mRNA measurement error (post hoc)")
DX0, DX1 = 130.0, 156.5
# logarithmic theta axis (brief 1b-i): on the linear 0-0.85 axis four cohorts sat in 0.03-0.24 and theta,
# range and dilution minimum overlapped; the log axis gives every row the same treatment
D_LIM = (0.02, 1.0)
D_TICKS = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
assert all(D_LIM[0] < min(D[c]["lo"], D[c]["dmin"]) and max(D[c]["hi"], D[c]["dmin"]) < D_LIM[1] for c in D_ORDER)
mD = logmapper(DX0, DX1, *D_LIM)
XR = 160.0
# key (defined once, applies to every row): one line, symbols on the text's centre line
DTICK = S.GREY          # mid grey (#8C8C8C), distinct from the dark-grey key dot and the cohort colours
KYD = 78.4
vtick(112.0, KYD - 1.0, KYD + 1.0, DTICK, lw=0.6, z=7)
T(113.4, KYD, "dilution minimum")
hline(132.0, 135.0, KYD, INK2, lw=0.8); dot(133.5, KYD, INK2, z=6, mec=INK2)
T(136.4, KYD, "θ with 95% patient-bootstrap interval")
T(111, 84.2, "cohort, proteins", va="baseline"); T(111, 86.5, "reliability λ (lower bound)", va="baseline")
T(XR, 84.2, "dilution by", va="baseline"); T(XR, 86.5, "RNA-seq noise", va="baseline")
ys = [93.5, 103.0, 112.5, 122.0, 131.5]
# printed theta vs dilution minimum, as main.tex prints them ("theta=0.07 ... below the dilution
# minimum 0.79"; "theta=0.10; 0.08; 0.13; 0.09" ... "(0.11; 0.11; 0.08; 0.11)")
D_PRINTED = {"CCRCC": "θ 0.07 vs min. 0.79", "LUAD": "θ 0.10 vs min. 0.11", "UCEC": "θ 0.08 vs min. 0.11",
             "GBM": "θ 0.13 vs min. 0.08", "PDAC": "θ 0.09 vs min. 0.11"}
for c, yc in zip(D_ORDER, ys):
    d = D[c]
    T(111, yc - 1.2, f"{c}, {commas(d['n'])}")
    T(111, yc + 1.2, f"λ {d['lam']}, {d['stt']}")
    col = ORG[c]
    xm = mD(d["dmin"])
    hline(mD(d["lo"]), mD(d["hi"]), yc, col, lw=0.8, z=4)
    dot(mD(d["theta"]), yc, col, ms=2.6, z=6, mec=col)
    vtick(xm, yc - 1.9, yc + 1.9, DTICK, lw=0.6, z=7)   # same grey/height/width on all five rows, above the dots
    # theta and minimum printed under the reading (stage-2 audit): in LUAD the dot sits 0.34 mm left of
    # the tick, too close to see at print size; both values are printed in main.tex (measurement-error
    # paragraph) at this precision
    th_s = f"θ {d['theta']:.2f} vs min. {d['dmin']:.2f}"
    assert th_s == D_PRINTED[c], (c, th_s)
    if d["bound"]:
        T(XR, yc - 2.4, "excluded up to a bound"); T(XR, yc, "(at most 0.08)"); T(XR, yc + 2.4, th_s)
    else:
        T(XR, yc - 1.2, "not excluded"); T(XR, yc + 1.2, th_s)
tick_axis(DX0, DX1, 136.8, D_LIM, D_TICKS, ["0.02", "0.05", "0.1", "0.2", "0.5", "1"], log=True)
T((DX0 + DX1) / 2, 141.5, "alignment θ with measured transcript (log scale)", ha="center")
# which side of its dilution minimum each row's theta and range lie (drawn; printed only as the reading)
for c in D_ORDER:
    d = D[c]
    W(f"  d {c}: theta {d['theta']:.4f} [{d['lo']:.4f}, {d['hi']:.4f}] vs min {d['dmin']:.4f}; "
      f"dot-to-tick {mD(d['theta']) - mD(d['dmin']):+.2f} mm, range end-to-tick {mD(d['hi']) - mD(d['dmin']):+.2f} mm")
    assert (d["theta"] < d["dmin"] and d["hi"] < d["dmin"]) == d["bound"], c   # reading rule R3 vs R2
T(111, 146.6, "linear contribution of the gene's own measured transcript only;")
T(111, 149.1, "batch-uncorrected significant sets")

# ------------------------------------------------------------------ checks
fig.canvas.draw()
rend = fig.canvas.get_renderer()
W_PX, H_PX = fig.bbox.width, fig.bbox.height
allt = list(TXT)
for ax in fig.axes:
    if ax is ov:
        continue
    allt += [t for t in ax.get_xticklabels() + ax.get_yticklabels() if t.get_text().strip()]
boxes = []
CLIP = []
for t in allt:
    if not t.get_text().strip():
        continue
    bb = t.get_window_extent(rend)
    boxes.append((t.get_text()[:40], bb))
    if not (bb.x0 >= -1 and bb.y0 >= -1 and bb.x1 <= W_PX + 1 and bb.y1 <= H_PX + 1):
        W("[clipped]", t.get_text(), bb.x1 / W_PX * FWmm)
        CLIP.append(t.get_text())
bad = []
for (s1, b1), (s2, b2) in itertools.combinations(boxes, 2):
    if b1.x0 < b2.x1 - 0.5 and b2.x0 < b1.x1 - 0.5 and b1.y0 < b2.y1 - 0.5 and b2.y0 < b1.y1 - 0.5:
        bad.append((s1, s2))
for p in bad:
    W("[overlap]", p)
W(f"text overlaps: {len(bad)}")
assert not bad and not CLIP, "overlapping or clipped text"

fig.savefig("Fig6_v7.svg", bbox_inches=None)
fig.savefig("Fig6_v7.pdf", bbox_inches=None)
fig.savefig("Fig6_v7.png", dpi=300, bbox_inches=None)
W(f"wrote Fig6_v7.{{svg,pdf,png}} -- canvas {FWmm:.1f} x {FHmm:.1f} mm")
