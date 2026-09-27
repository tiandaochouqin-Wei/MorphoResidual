#!/usr/bin/env python3
"""Corrected-baseline external-replication recalculation for LUAD/GBM/PDAC, extending
review/recalc/kirc_baserate.py's method (verified against its published numbers below)
to the three RPPA cohorts that were run after KIRC but never had this correction applied
(audit F8/F9: "base rate" must exclude the CPTAC-significant proteins from the
denominator, and both one- and two-sided Fisher p plus 95% CI must be reported).

Per cohort:
  - "replicated" (numerator) = of the RPPA proteins that are CPTAC-<cohort> significant
    (fdr<0.05 & incremental_r2>0 in the discovery result) AND tested on RPPA, the count
    that are ALSO significant in the RPPA cohort (same two-condition definition).
  - "base rate" (denominator frame) = among the tested-on-RPPA proteins OUTSIDE the
    CPTAC-significant set, the fraction significant in the RPPA cohort. This is the
    frame KIRC used for its headline Fisher test (the 297-of-360 frame); the
    all-shared-tested frame is reported as a sensitivity check, matching kirc_baserate.md.
  - Wilson score 95% CI for both proportions.
  - Fisher exact, both one- and two-sided p (F9's explicit requirement).
  - Pearson r / Spearman rho of incremental_r2 across the shared-tested proteins.
Then a Cochran-Mantel-Haenszel test stratified by cohort across all four RPPA cohorts
(KIRC/LUAD/GBM/PDAC), continuity-corrected, plus Breslow-Day for homogeneity of the OR
across strata (F9).
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable

FD = "figures/figdata"
COHORTS = {
    "kirc":  dict(rppa=f"{FD}/kirc_rppa_results.csv",  cptac="server_export/pinned/ccrcc/residual_results_tumoronly.csv", n=451, label="TCGA-KIRC"),
    "luad":  dict(rppa=f"{FD}/luad_rppa_results.csv",  cptac="server_export/pinned/luad/residual_results_tumoronly.csv",  n=328, label="TCGA-LUAD"),
    "gbm":   dict(rppa=f"{FD}/gbm_rppa_results.csv",   cptac="server_export/pinned/gbm/residual_results_tumoronly.csv",   n=85,  label="TCGA-GBM"),
    "pdac":  dict(rppa=f"{FD}/paad_rppa_results.csv",  cptac="server_export/pinned/pdac/residual_results_tumoronly.csv",  n=112, label="TCGA-PDAC"),
}


def wilson_ci(k, n, z=1.959963984540054):
    if n == 0:
        return (float("nan"),) * 2
    phat = k / n
    denom = 1 + z**2 / n
    centre = phat + z**2 / (2 * n)
    adj = z * np.sqrt(phat * (1 - phat) / n + z**2 / (4 * n**2))
    return ((centre - adj) / denom, (centre + adj) / denom)


def sig(df):
    return (df.fdr < 0.05) & (df.incremental_r2 > 0)


def one_cohort(key, cfg):
    rppa = pd.read_csv(cfg["rppa"]).drop_duplicates("gene")
    cptac = pd.read_csv(cfg["cptac"]).drop_duplicates("gene")
    rppa_sig = set(rppa.gene[sig(rppa)])
    cptac_sig = set(cptac.gene[sig(cptac)])

    tested = set(rppa.gene)
    shared_tested = tested & set(cptac.gene)  # RPPA-tested proteins also tested in CPTAC discovery
    cptac_sig_on_rppa = shared_tested & cptac_sig
    replicated = cptac_sig_on_rppa & rppa_sig
    a, b = len(replicated), len(cptac_sig_on_rppa) - len(replicated)

    # HEADLINE frame (kirc_baserate.md's published Results number, e.g. KIRC's
    # [[35,28],[133,164]] OR=1.54): denominator is ALL RPPA-tested proteins minus the
    # CPTAC-significant-on-RPPA ones, whether or not each was itself tested in the CPTAC
    # discovery cohort -- so it also includes RPPA proteins CPTAC never measured.
    base_all = tested - cptac_sig_on_rppa
    hits_all = base_all & rppa_sig
    c, d = len(hits_all), len(base_all) - len(hits_all)
    table = [[a, b], [c, d]]
    or_, p_two = stats.fisher_exact(table, alternative="two-sided")
    _, p_one = stats.fisher_exact(table, alternative="greater")

    # SENSITIVITY frame (kirc_baserate.md's "301 shared-tested frame"): denominator
    # restricted to proteins tested in BOTH RPPA and the CPTAC discovery cohort.
    non_sig_on_rppa = shared_tested - cptac_sig_on_rppa
    baseline_hits = non_sig_on_rppa & rppa_sig
    c_s, d_s = len(baseline_hits), len(non_sig_on_rppa) - len(baseline_hits)
    or_s, p_two_s = stats.fisher_exact([[a, b], [c_s, d_s]], alternative="two-sided")
    _, p_one_s = stats.fisher_exact([[a, b], [c_s, d_s]], alternative="greater")
    rate_base_narrow = c_s / len(non_sig_on_rppa) if non_sig_on_rppa else float("nan")
    ci_base_narrow = wilson_ci(c_s, len(non_sig_on_rppa))

    m = rppa.merge(cptac, on="gene", suffixes=("_rppa", "_cptac"))
    m = m[m.gene.isin(shared_tested)]
    pear = stats.pearsonr(m.incremental_r2_cptac, m.incremental_r2_rppa)
    spear = stats.spearmanr(m.incremental_r2_cptac, m.incremental_r2_rppa)

    rate_num = a / len(cptac_sig_on_rppa) if cptac_sig_on_rppa else float("nan")
    rate_base = c / len(base_all) if base_all else float("nan")
    ci_num = wilson_ci(a, len(cptac_sig_on_rppa))
    ci_base = wilson_ci(c, len(base_all))
    tot_sig = len(rppa_sig & tested)

    print(f"\n=== {cfg['label']} (n={cfg['n']}) ===")
    print(f"  RPPA tested: {len(tested)} | phenomenon-level significant: {tot_sig} "
          f"({tot_sig/len(tested):.1%})")
    print(f"  shared-tested (RPPA ∩ CPTAC-discovery-tested): {len(shared_tested)}")
    print(f"  CPTAC-{key}-significant on RPPA: {len(cptac_sig_on_rppa)}  -> replicated: {a} "
          f"({rate_num:.1%}, 95% CI {ci_num[0]:.1%}-{ci_num[1]:.1%})")
    print(f"  WIDE base rate -- all tested minus CPTAC-sig, {len(base_all)} proteins "
          f"(superseded, kept only for the manuscript's disclosed-uncorrected-frame "
          f"comparison): {c} ({rate_base:.1%}, 95% CI {ci_base[0]:.1%}-{ci_base[1]:.1%})")
    print(f"  Fisher table [[{a},{b}],[{c},{d}]]  OR={or_:.3f}  two-sided p={p_two:.4f}  "
          f"one-sided p={p_one:.4f}   <- WIDE (uncorrected)")
    print(f"  NARROW base rate -- shared-tested-only frame, {len(non_sig_on_rppa)} proteins "
          f"(this is the manuscript's primary number): {c_s} ({rate_base_narrow:.1%}, "
          f"95% CI {ci_base_narrow[0]:.1%}-{ci_base_narrow[1]:.1%})")
    print(f"  Fisher table [[{a},{b}],[{c_s},{d_s}]]  OR={or_s:.3f}  two-sided p={p_two_s:.4f}  "
          f"one-sided p={p_one_s:.4f}   <- NARROW (manuscript primary)")
    print(f"  per-gene agreement (n={len(m)}): Pearson r={pear.statistic:.4f} (p={pear.pvalue:.3f})  "
          f"Spearman rho={spear.statistic:.4f} (p={spear.pvalue:.3f})")

    return dict(cohort=key, label=cfg["label"], n_cases=cfg["n"], tested=len(tested),
                phenomenon_sig=tot_sig, shared_tested=len(shared_tested),
                cptac_sig_on_rppa=len(cptac_sig_on_rppa), replicated=a,
                rate_num=rate_num, ci_num_lo=ci_num[0], ci_num_hi=ci_num[1],
                a=a, b=b,
                # NARROW frame = the manuscript's primary/reported numbers.
                narrow_base_hits=c_s, narrow_base_n=len(non_sig_on_rppa),
                narrow_rate_base=rate_base_narrow,
                narrow_ci_base_lo=ci_base_narrow[0], narrow_ci_base_hi=ci_base_narrow[1],
                narrow_OR=or_s, narrow_p_two=p_two_s, narrow_p_one=p_one_s,
                # WIDE frame = superseded/uncorrected, kept for the disclosed comparison only.
                wide_base_hits=c, wide_base_n=len(base_all), wide_rate_base=rate_base,
                wide_ci_base_lo=ci_base[0], wide_ci_base_hi=ci_base[1],
                wide_OR=or_, wide_p_two=p_two, wide_p_one=p_one,
                pearson_r=pear.statistic, pearson_p=pear.pvalue,
                spearman_rho=spear.statistic, spearman_p=spear.pvalue)


def main():
    rows = [one_cohort(k, cfg) for k, cfg in COHORTS.items()]
    df = pd.DataFrame(rows)
    out = "review/recalc/rppa_baserate_all.csv"
    df.to_csv(out, index=False)

    print("\n" + "=" * 78)
    print("CROSS-CHECK vs kirc_baserate.md published numbers (must match to validate the script)")
    print("=" * 78)
    k = df[df.cohort == "kirc"].iloc[0]
    print(f"  expected (wide/uncorrected frame): a=35 b=28 c=133 d=164 OR=1.54 two-sided p=0.128 one-sided p=0.078")
    print(f"  got     : a={k.a} b={k.b} c={k.wide_base_hits} d={k.wide_base_n - k.wide_base_hits} "
          f"OR={k.wide_OR:.2f} two-sided p={k.wide_p_two:.3f} one-sided p={k.wide_p_one:.3f}")

    def cmh_block(frame, hits_col, n_col, label):
        print("\n" + "=" * 78)
        print(f"STRATIFIED (CMH) ACROSS THE FOUR RPPA COHORTS -- {label} FRAME")
        print("=" * 78)
        tabs = [(r.a, r.b, getattr(r, hits_col), getattr(r, n_col) - getattr(r, hits_col))
                for r in df.itertuples()]
        tables = np.array([[[a, b], [c, d]] for a, b, c, d in tabs]).transpose(1, 2, 0)
        st = StratifiedTable(tables)
        cmh = st.test_null_odds(correction=True)
        pooled_or = st.oddsratio_pooled
        ci_lo, ci_hi = st.oddsratio_pooled_confint()
        print(f"  per-stratum tables: {tabs}")
        print(f"  Mantel-Haenszel pooled OR = {pooled_or:.3f}  95% CI [{ci_lo:.3f}, {ci_hi:.3f}]")
        print(f"  CMH (continuity-corrected): statistic={cmh.statistic:.4f}  p={cmh.pvalue:.4f}")
        # Breslow-Day over the non-degenerate cohorts only (GBM's zero cells make it undefined).
        keep = [i for i, (a, b, c, d) in enumerate(tabs) if not (a == 0 and c == 0)]
        tabs_nd = [tabs[i] for i in keep]
        tables_nd = np.array([[[a, b], [c, d]] for a, b, c, d in tabs_nd]).transpose(1, 2, 0)
        bd = StratifiedTable(tables_nd).test_equal_odds()
        print(f"  Breslow-Day (excl. degenerate strata, n={len(tabs_nd)}): "
              f"statistic={bd.statistic:.4f}  p={bd.pvalue:.4f}  "
              f"(p<0.05 => OR is NOT homogeneous across cohorts)")
        return pooled_or, ci_lo, ci_hi, cmh.pvalue, bd.pvalue

    print("\nNARROW frame is the manuscript's primary CMH result "
          "(SuppTable_replication.tex, main.tex sec:kircres).")
    cmh_block(df, "narrow_base_hits", "narrow_base_n", "NARROW (manuscript primary)")
    print("\nWIDE frame is reported only as the disclosed superseded/uncorrected comparison.")
    cmh_block(df, "wide_base_hits", "wide_base_n", "WIDE (superseded)")

    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
