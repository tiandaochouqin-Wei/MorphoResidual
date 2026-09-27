#!/usr/bin/env python3
"""External-replication statistics under the stated RPPA antibody rule.

WHY THIS EXISTS
The Methods say RPPA antibodies were mapped to genes "excluding phospho- and other
modification-specific antibodies". The code (tcga_replicate.py / tcga_kirc_omics.py)
excluded only phospho-specific antibodies (regex _p[STY]\\d). Eight non-phospho
modification-specific antibodies (cleaved caspase-7/-8, PARP, Notch1; acetyl-alpha-
tubulin Lys40; three histone H3/H2B modification antibodies) therefore fed seven
gene-level proteins in every panel: CASP7, CASP8, PARP1, NOTCH1, TUBA1B, H3C1, H2BC3.
CASP7, TUBA1B and H2BC3 are represented ONLY by a modified antibody; CASP8, PARP1,
NOTCH1 and H3C1 average a modified antibody with the total-protein antibody.

WHAT IT DOES
Two definitions, computed side by side from the checked-in per-gene results
(figures/figdata/{kirc,luad,gbm,paad}_rppa_results.csv; columns gene, incremental_r2,
pval, fdr):
  as_run    -- the panels exactly as analysed (reproduces every number the manuscript
               carried before this correction; asserted below)
  excluded  -- the stated rule: all seven genes removed, then Benjamini-Hochberg re-run
               on the remaining genes' permutation p-values (the p-values are per gene;
               nothing else changes). Removing the whole gene is the conservative reading
               for the four averaged genes, whose total-protein antibody would otherwise
               remain; keeping them would need antibody-level data that are not stored
               locally, and adds at most four proteins per panel.
The BH implementation is asserted to reproduce the stored FDR exactly. No permutation is
repeated.

OUTPUT  review/recalc/rppa_exclusion_results.json  (every statistic, both definitions)
        review/recalc/rppa_baserate_all.csv        ('excluded': primary analysis, read by
                                                    the figure and Source Data scripts)
        review/recalc/rppa_baserate_all_asrun.csv  ('as_run': the sensitivity analysis)
This script supersedes the earlier rppa_baserate_all.py, whose as-run output it reproduces
exactly (asserted below).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures"))
from rppa_filter import MODIFIED_GENES, bh, load_rppa  # noqa: E402

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FD = os.path.join(ROOT, "figures", "figdata")
COHORTS = {
    "kirc": dict(rppa="kirc_rppa_results.csv", cptac="ccrcc", n=451, label="TCGA-KIRC"),
    "luad": dict(rppa="luad_rppa_results.csv", cptac="luad", n=328, label="TCGA-LUAD"),
    "gbm": dict(rppa="gbm_rppa_results.csv", cptac="gbm", n=85, label="TCGA-GBM"),
    "pdac": dict(rppa="paad_rppa_results.csv", cptac="pdac", n=112, label="TCGA-PDAC"),
}


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return [float("nan")] * 2
    ph = k / n; den = 1 + z ** 2 / n; ctr = ph + z ** 2 / (2 * n)
    adj = z * np.sqrt(ph * (1 - ph) / n + z ** 2 / (4 * n ** 2))
    return [float((ctr - adj) / den), float((ctr + adj) / den)]


def load(key, mode):
    return load_rppa(os.path.join(FD, COHORTS[key]["rppa"]), exclude=(mode == "excluded"))


def sig(df):
    return (df.fdr < 0.05) & (df.incremental_r2 > 0)


def cohort_stats(key, mode):
    cfg = COHORTS[key]
    rppa = load(key, mode)
    cptac = pd.read_csv(os.path.join(ROOT, "server_export", "pinned", cfg["cptac"], "residual_results_tumoronly.csv")).drop_duplicates("gene")
    rppa_sig = set(rppa.gene[sig(rppa)]); cptac_sig = set(cptac.gene[sig(cptac)])
    tested = set(rppa.gene); shared = tested & set(cptac.gene)
    cs_on = shared & cptac_sig; rep = cs_on & rppa_sig
    a, b = len(rep), len(cs_on) - len(rep)
    nonsig = shared - cs_on
    c_s = len(nonsig & rppa_sig); d_s = len(nonsig) - c_s
    or_s, p2 = stats.fisher_exact([[a, b], [c_s, d_s]], alternative="two-sided")
    _, p1 = stats.fisher_exact([[a, b], [c_s, d_s]], alternative="greater")
    base_all = tested - cs_on; c_w = len(base_all & rppa_sig); d_w = len(base_all) - c_w
    or_w, p2w = stats.fisher_exact([[a, b], [c_w, d_w]], alternative="two-sided")
    _, p1w = stats.fisher_exact([[a, b], [c_w, d_w]], alternative="greater")
    m = rppa.merge(cptac, on="gene", suffixes=("_rppa", "_cptac")); m = m[m.gene.isin(shared)]
    pr = stats.pearsonr(m.incremental_r2_cptac, m.incremental_r2_rppa)
    sp = stats.spearmanr(m.incremental_r2_cptac, m.incremental_r2_rppa)
    sigdf = rppa[sig(rppa)]
    out = dict(
        cohort=key, label=cfg["label"], n_cases=cfg["n"], tested=len(tested), significant=int(len(rppa_sig)),
        pct_sig=len(rppa_sig) / len(tested), median_incr_sig=float(sigdf.incremental_r2.median()) if len(sigdf) else None,
        any_sign_fdr005=int((rppa.fdr < 0.05).sum()), shared_tested=len(shared),
        cptac_sig_on_rppa=len(cs_on), replicated=a, a=a, b=b, rate_num=a / len(cs_on) if cs_on else float("nan"),
        ci_num=wilson(a, len(cs_on)),
        narrow_base_hits=c_s, narrow_base_n=len(nonsig), narrow_rate_base=c_s / len(nonsig), narrow_ci_base=wilson(c_s, len(nonsig)),
        narrow_OR=float(or_s), narrow_p_two=float(p2), narrow_p_one=float(p1),
        wide_base_hits=c_w, wide_base_n=len(base_all), wide_rate_base=c_w / len(base_all), wide_ci_base=wilson(c_w, len(base_all)),
        wide_OR=float(or_w), wide_p_two=float(p2w), wide_p_one=float(p1w),
        pearson_r=float(pr.statistic), pearson_p=float(pr.pvalue), spearman_rho=float(sp.statistic), spearman_p=float(sp.pvalue),
        n_pairs=len(m))
    if key == "kirc":
        # discrimination detail (Supplementary Table, 'KIRC discrimination detail')
        others = m[~m.gene.isin(cs_on)].incremental_r2_rppa.values     # shared-tested, not CPTAC-sig
        mine = m[m.gene.isin(cs_on)].incremental_r2_rppa.values
        u = stats.mannwhitneyu(mine, others, alternative="greater")
        allother = rppa[~rppa.gene.isin(cs_on)].incremental_r2.values  # every other tested panel protein
        u2 = stats.mannwhitneyu(mine, allother, alternative="greater")
        lr = stats.linregress(m.incremental_r2_cptac, m.incremental_r2_rppa)
        n_by = rppa.n.value_counts().to_dict()
        out.update(mw_greater_p_vs_shared_nonsig=float(u.pvalue), mw_greater_p_vs_all_other_tested=float(u2.pvalue),
                   auroc_vs_shared_nonsig=float(u.statistic / (len(mine) * len(others))),
                   auroc_vs_all_other_tested=float(u2.statistic / (len(mine) * len(allother))),
                   median_mine=float(np.median(mine)), median_others_shared=float(np.median(others)),
                   median_all_other_tested=float(np.median(allother)),
                   calib_slope=float(lr.slope), calib_se=float(lr.stderr), calib_r=float(lr.rvalue), calib_p=float(lr.pvalue),
                   n_genes_calib=len(m), n_value_counts={int(k): int(v) for k, v in n_by.items()})
    return out


def cmh(rows, hits, ncol):
    tabs = [(r["a"], r["b"], r[hits], r[ncol] - r[hits]) for r in rows]
    T = np.array([[[a, b], [c, d]] for a, b, c, d in tabs]).transpose(1, 2, 0)
    st = StratifiedTable(T); t = st.test_null_odds(correction=True); lo, hi = st.oddsratio_pooled_confint()
    keep = [i for i, (a, b, c, d) in enumerate(tabs) if not (a == 0 and c == 0)]
    Tn = np.array([[[tabs[i][0], tabs[i][1]], [tabs[i][2], tabs[i][3]]] for i in keep]).transpose(1, 2, 0)
    bd = StratifiedTable(Tn).test_equal_odds()
    return dict(pooled_OR=float(st.oddsratio_pooled), ci=[float(lo), float(hi)], cmh_p=float(t.pvalue), cmh_stat=float(t.statistic),
                breslow_day_p=float(bd.pvalue), tables=[list(map(int, t_)) for t_ in tabs])


def main():
    res = {}
    for mode in ("as_run", "excluded"):
        rows = [cohort_stats(k, mode) for k in COHORTS]
        res[mode] = dict(cohorts={r["cohort"]: r for r in rows}, cmh_narrow=cmh(rows, "narrow_base_hits", "narrow_base_n"),
                         cmh_wide=cmh(rows, "wide_base_hits", "wide_base_n"))
    # ---- validate the 'as_run' arm against the numbers the manuscript carried before this correction ----
    A = res["as_run"]["cohorts"]
    expect = {  # (tested, significant, cptac_sig_on_rppa, replicated, narrow base n, narrow base hits)
        "kirc": (360, 168, 63, 35, 238, 111), "luad": (360, 193, 80, 45, None, None),
        "gbm": (355, 0, 24, 0, None, None), "pdac": (350, 64, 43, 7, None, None)}
    for k, e in expect.items():
        got = (A[k]["tested"], A[k]["significant"], A[k]["cptac_sig_on_rppa"], A[k]["replicated"])
        assert got == e[:4], (k, got, e[:4])
    assert (A["kirc"]["narrow_base_n"], A["kirc"]["narrow_base_hits"]) == (238, 111)
    assert round(A["kirc"]["narrow_OR"], 2) == 1.43 and round(A["kirc"]["narrow_p_two"], 3) == 0.257
    assert round(A["kirc"]["pearson_r"], 3) == -0.066 and round(A["kirc"]["spearman_rho"], 3) == 0.048
    assert round(A["kirc"]["wide_OR"], 2) == 1.54 and round(A["kirc"]["wide_p_two"], 3) == 0.128
    assert A["kirc"]["any_sign_fdr005"] == 223 and A["kirc"]["n_genes_calib"] == 301
    assert (round(A["luad"]["narrow_OR"], 2), round(A["luad"]["narrow_p_two"], 3)) == (1.02, 1.0)
    assert (round(A["pdac"]["narrow_OR"], 2), round(A["pdac"]["narrow_p_two"], 3)) == (1.07, 0.823)
    c = res["as_run"]["cmh_narrow"]
    assert round(c["pooled_OR"], 2) == 1.17 and [round(x, 2) for x in c["ci"]] == [0.83, 1.66]
    assert round(c["cmh_p"], 2) == 0.42 and round(c["breslow_day_p"], 2) == 0.67
    print("as_run arm reproduces every pre-correction manuscript number checked (asserted).")
    # ---- outputs ----
    with open(os.path.join(ROOT, "review", "recalc", "rppa_exclusion_results.json"), "w", encoding="utf-8") as f:
        json.dump(dict(modified_genes=MODIFIED_GENES, **res), f, indent=1)
    for mode, fn in (("excluded", "rppa_baserate_all.csv"), ("as_run", "rppa_baserate_all_asrun.csv")):
        rows = []
        for k, r in res[mode]["cohorts"].items():
            rows.append(dict(cohort=k, label=r["label"], n_cases=r["n_cases"], tested=r["tested"], phenomenon_sig=r["significant"],
                             shared_tested=r["shared_tested"], cptac_sig_on_rppa=r["cptac_sig_on_rppa"], replicated=r["replicated"],
                             rate_num=r["rate_num"], ci_num_lo=r["ci_num"][0], ci_num_hi=r["ci_num"][1], a=r["a"], b=r["b"],
                             narrow_base_hits=r["narrow_base_hits"], narrow_base_n=r["narrow_base_n"], narrow_rate_base=r["narrow_rate_base"],
                             narrow_ci_base_lo=r["narrow_ci_base"][0], narrow_ci_base_hi=r["narrow_ci_base"][1],
                             narrow_OR=r["narrow_OR"], narrow_p_two=r["narrow_p_two"], narrow_p_one=r["narrow_p_one"],
                             wide_base_hits=r["wide_base_hits"], wide_base_n=r["wide_base_n"], wide_rate_base=r["wide_rate_base"],
                             wide_ci_base_lo=r["wide_ci_base"][0], wide_ci_base_hi=r["wide_ci_base"][1],
                             wide_OR=r["wide_OR"], wide_p_two=r["wide_p_two"], wide_p_one=r["wide_p_one"],
                             pearson_r=r["pearson_r"], pearson_p=r["pearson_p"], spearman_rho=r["spearman_rho"], spearman_p=r["spearman_p"]))
        pd.DataFrame(rows).to_csv(os.path.join(ROOT, "review", "recalc", fn), index=False)
    # ---- print side-by-side ----
    def line(mode, k):
        r = res[mode]["cohorts"][k]
        return (f"{mode:8s} {k}: tested {r['tested']} sig {r['significant']} ({r['pct_sig']:.1%}) med {r['median_incr_sig']} | "
                f"CPTAC-sig-on-RPPA {r['cptac_sig_on_rppa']} rep {r['replicated']} ({r['rate_num']:.1%}) base {r['narrow_base_hits']}/{r['narrow_base_n']} "
                f"({r['narrow_rate_base']:.1%}) OR {r['narrow_OR']:.2f} p2 {r['narrow_p_two']:.3f} p1 {r['narrow_p_one']:.3f} | r {r['pearson_r']:.3f} rho {r['spearman_rho']:.3f}")
    for k in COHORTS:
        print(line("as_run", k)); print(line("excluded", k))
    for mode in ("as_run", "excluded"):
        for fr in ("cmh_narrow", "cmh_wide"):
            c = res[mode][fr]; print(f"{mode:8s} {fr}: OR {c['pooled_OR']:.3f} CI [{c['ci'][0]:.3f},{c['ci'][1]:.3f}] CMH p {c['cmh_p']:.3f} BD p {c['breslow_day_p']:.3f}")
    for mode in ("as_run", "excluded"):
        r = res[mode]["cohorts"]["kirc"]
        print(f"{mode:8s} KIRC detail: any-sign {r['any_sign_fdr005']} ({r['any_sign_fdr005']/r['tested']:.1%}); MW greater p vs shared-nonsig {r['mw_greater_p_vs_shared_nonsig']:.3f}, vs all other tested {r['mw_greater_p_vs_all_other_tested']:.3f}; "
              f"AUROC {r['auroc_vs_shared_nonsig']:.3f}/{r['auroc_vs_all_other_tested']:.3f}; medians {r['median_mine']:.4f} vs {r['median_others_shared']:.4f}/{r['median_all_other_tested']:.4f}; "
              f"calib slope {r['calib_slope']:.3f} se {r['calib_se']:.3f} r {r['calib_r']:.3f} p {r['calib_p']:.3f} n {r['n_genes_calib']}; n counts {r['n_value_counts']}")


if __name__ == "__main__":
    sys.exit(main())
