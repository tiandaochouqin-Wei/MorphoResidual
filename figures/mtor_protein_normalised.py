#!/usr/bin/env python3
r"""
mtor_protein_normalised.py -- MorphoResidual T1-4 (NC_ROADMAP.md).

THE QUESTION
The manuscript reports that phospho-S6 tracks the translation residual in all five
cancers, and then honestly refuses to call it mTORC1 activity, because the
phosphosite intensities are not normalised to the abundance of the protein that
carries the site. Since the residual is itself dominated by ribosomal proteins, a
phospho-S6 signal could be reporting RPS6 protein abundance rather than S6K1
activity. Section 3.14 currently calls protein-normalised ratios "the required
direct test" and leaves it for future work. This does it.

THE TEST
For every phosphosite we form, in log space and per patient,

    normalised site  =  log phosphosite intensity  -  log host-protein abundance

using the same protein matrix the main pipeline uses (exported by
export_phospho_hosts.py). We then ask two things:

  1. does the normalised S6K1 arm still track the translation residual?
  2. where does it sit in the distribution of the NORMALISED non-mTORC1 ribosomal
     phosphosites?

Point 2 is the whole point, and it is why the background sites are normalised to
THEIR OWN host proteins too. Comparing a normalised pS6 against unnormalised
ribosomal sites would be an asymmetric control and worse than the honest null the
paper currently reports.

INTERPRETATION, PRE-COMMITTED BEFORE LOOKING
  * normalised pS6 still tracks the residual AND now sits above the normalised
    ribosomal background  -> mTORC1-specific reading is supported; the manuscript
    can be upgraded from "cannot establish" to "consistent with", still short of
    proof.
  * normalised pS6 sits inside the background -> the existing conclusion stands,
    but is now TESTED rather than deferred, which is the point.

A CAVEAT THAT MUST BE REPORTED EITHER WAY
The translation score is the mean z-scored residual over the significant genes of
the translation family. If RPS6 is itself in that set, then normalising pS6 by
RPS6 while RPS6 contributes to the score creates a ratio-versus-component
artefact biased towards a negative correlation. Checking gene_covariates_*.csv:
RPS6 is in the significant set only in LUAD, where it is 1 of 81 genes (1.2%);
EIF4EBP1 is in no cohort's set. So the artefact is confined to one cohort at ~1%
leverage, and CCRCC, UCEC, GBM and PDAC give the clean test. The script reports
which cohorts are clean rather than assuming it.

Writes figures/figdata/mtor_protein_normalised.csv (per-arm summary)
       figures/figdata/mtor_protein_normalised_sites.csv (per-site detail)
"""
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats as st

from mtor_common import MTOR, RIBO_PAT, load_phospho, norm_id, parse_site

DD = "figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
MIN_N = 25                       # patients needed for a correlation to be reported

# whether RPS6 / EIF4EBP1 are in that cohort's translation-family significant set
# (recomputed here rather than hardcoded, so the caveat cannot go stale)
TRANS_EXTRA = {"EEF1A1", "EEF1B2", "EEF1G", "EEF2", "EIF3A", "EIF3B", "EIF3C", "EIF3D",
               "EIF3E", "EIF3F", "EIF3G", "EIF3H", "EIF3I", "EIF3L", "EIF4A1", "EIF4A2",
               "EIF4B", "EIF4E", "EIF4G1", "EIF4G2", "EIF2S1", "EIF2S2", "EIF2S3",
               "EIF5", "EIF5A", "EIF6", "PABPC1"}


def zscore(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(ddof=0)
    return (s - s.mean()) / (sd if sd and np.isfinite(sd) else 1.0)


def family_membership(cohort):
    """Is RPS6 (or EIF4EBP1) inside this cohort's translation score?"""
    try:
        d = pd.read_csv(f"{DD}/gene_covariates_{cohort}.csv")
    except FileNotFoundError:
        return None, None, None
    sig = set(d.loc[d.pinned_significant == True, "gene"].astype(str))  # noqa: E712
    fam = {g for g in sig if RIBO_PAT.match(g) or g in TRANS_EXTRA}
    return ("RPS6" in sig), ("EIF4EBP1" in sig), len(fam)


rows, site_rows = [], []
for c in COH:
    print(f"\n=== {c.upper()} ===")
    try:
        ph = load_phospho(f"{DD}/phospho_{c}.txt")
        hosts = pd.read_csv(f"{DD}/phospho_hosts_{c}.csv", index_col=0)
        scores = pd.read_csv(f"{DD}/scores_{c}.csv").rename(columns={"Unnamed: 0": "case"}).set_index("case")
    except FileNotFoundError as exc:
        print(f"  skipped: {exc}")
        continue
    hosts.index = [norm_id(i) for i in hosts.index]
    scores.index = [norm_id(i) for i in scores.index]

    common = sorted(set(ph.columns) & set(hosts.index) & set(scores.index))
    print(f"  patients with phospho + host protein + score: {len(common)}")
    if len(common) < MIN_N:
        print("  too few; skipped")
        continue
    y = scores.loc[common, "translation_meas"].astype(float)
    y_morph = scores.loc[common, "translation_morph"].astype(float)

    r6_in, e4_in, n_fam = family_membership(c)
    clean = not r6_in
    print(f"  translation family significant genes: {n_fam}; RPS6 in score: {r6_in}; "
          f"EIF4EBP1 in score: {e4_in}  -> {'CLEAN test' if clean else 'ratio-vs-component caveat applies'}")

    # Gene-level detection must NOT use "does the id contain [STY]<digits>": CCRCC's
    # row ids are bare gene symbols and RPS6 / RPS10 then read as sites S6 / S10.
    # That exact trap is documented in mtor_mechanism.py. Detect instead by asking
    # how many row ids ARE gene symbols we hold protein for.
    # precise test: in a gene-level matrix the bare mTORC1 symbols ARE row ids
    # ("RPS6"); in a site-level one they never are ("RPS6_S235s").
    gene_level = bool(set(MTOR) & set(ph.index.astype(str)))
    if gene_level:
        print("  gene-level matrix: matching on bare gene symbols, and no ribosomal "
              "background is possible (symbol digits are not phosphosites)")

    def site_series(rid):
        return ph.loc[rid, common].astype(float)

    # ---------------- collect sites: mTORC1 arms + non-RPS6 ribosomal background
    collected = {}          # sid -> (host_gene, raw series)
    for rid in ph.index:
        p = parse_site(rid, gene_hint=rid if gene_level else None)
        if not p:
            continue
        g, res = p
        if gene_level:
            # exact gene-symbol rows only; residues parsed out of a symbol are noise
            is_mtor = str(rid) in MTOR
            is_ribo_bg = False          # cannot be built from a gene-level matrix
            res = set()
        else:
            is_mtor = g in MTOR and bool(res & set(MTOR[g]))
            is_ribo_bg = bool(RIBO_PAT.match(g)) and g != "RPS6"
        if not (is_mtor or is_ribo_bg):
            continue
        if g not in hosts.columns:          # cannot normalise without the host protein
            continue
        tag = f"{g} {'/'.join(sorted(res))}" if res else f"{g} (gene-level)"
        if tag in collected:
            continue
        # Residues must be re-extracted with the gene symbol removed, or RPS6 donates
        # a phantom "S6" and no peptide is ever a subset of the canonical set.
        canonical = set(MTOR.get(g, []))
        body = re.sub(re.escape(g), "", str(rid), count=1, flags=re.I)
        res_clean = set(re.findall(r"[STY]\d{1,5}", body.upper()))
        strict = bool(res_clean) and res_clean.issubset(canonical) if is_mtor else False
        collected[tag] = (g, site_series(rid),
                          "mTORC1" if is_mtor else "ribosomal background", strict)

    print(f"  normalisable sites: {len(collected)} "
          f"({sum(1 for v in collected.values() if v[2]=='mTORC1')} mTORC1)")

    for tag, (g, raw, kind, strict) in collected.items():
        host = hosts.loc[common, g].astype(float)
        m = raw.notna() & host.notna() & y.notna()
        if m.sum() < MIN_N:
            continue
        raw_z = zscore(raw[m])
        norm_z = zscore(raw[m].values - host[m].values)      # log ratio: site minus host
        rho_raw = st.spearmanr(raw_z, y[m])[0]
        rho_norm = st.spearmanr(norm_z, y[m])[0]
        rho_host = st.spearmanr(host[m], y[m])[0]
        site_rows.append(dict(cohort=c.upper(), site=tag, host_gene=g, kind=kind,
                              canonical_only=strict, n=int(m.sum()),
                              rho_raw=rho_raw, rho_norm=rho_norm,
                              rho_host_protein=rho_host))

    # ---------------- arm-level summary
    sub = pd.DataFrame([r for r in site_rows if r["cohort"] == c.upper()])
    if sub.empty:
        continue
    bg = sub[sub.kind == "ribosomal background"]
    for arm, genes in [("S6K1 arm (pRPS6)", ["RPS6"]),
                       ("4E-BP1 arm (p4E-BP1)", ["EIF4EBP1"]),
                       ("other mTORC1 (EIF4B, S6K1)", ["EIF4B", "RPS6KB1"])]:
        a = sub[(sub.kind == "mTORC1") & (sub.host_gene.isin(genes))]
        if a.empty:
            continue
        pct_raw = 100 * (bg.rho_raw < a.rho_raw.mean()).mean() if len(bg) else np.nan
        pct_norm = 100 * (bg.rho_norm < a.rho_norm.mean()).mean() if len(bg) else np.nan
        rows.append(dict(cohort=c.upper(), arm=arm, n_sites=len(a),
                         rho_raw=a.rho_raw.mean(), rho_norm=a.rho_norm.mean(),
                         bg_n=len(bg), bg_median_raw=bg.rho_raw.median() if len(bg) else np.nan,
                         bg_median_norm=bg.rho_norm.median() if len(bg) else np.nan,
                         pctile_raw=pct_raw, pctile_norm=pct_norm,
                         rps6_in_score=r6_in, clean_test=clean))
        print(f"    {arm:28s} raw rho={a.rho_raw.mean():+.3f} (pctile {pct_raw:5.1f}) "
              f"-> normalised rho={a.rho_norm.mean():+.3f} (pctile {pct_norm:5.1f}) "
              f"| bg median {bg.rho_norm.median():+.3f} of {len(bg)}")

res = pd.DataFrame(rows)
sites = pd.DataFrame(site_rows)
res.to_csv(f"{DD}/mtor_protein_normalised.csv", index=False)
sites.to_csv(f"{DD}/mtor_protein_normalised_sites.csv", index=False)

print("\n" + "=" * 78)
print("SUMMARY -- does the S6K1 arm survive host-protein normalisation, and does it")
print("rise above the equally-normalised ribosomal background?")
s6 = res[res.arm.str.startswith("S6K1")]
if not s6.empty:
    print(s6[["cohort", "rho_raw", "rho_norm", "bg_median_norm", "pctile_norm",
              "bg_n", "clean_test"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    above = s6[s6.pctile_norm > 90]
    print(f"\n  cohorts where the normalised S6K1 arm is above the 90th percentile of the "
          f"normalised ribosomal background: {len(above)}/{len(s6)}"
          + (f" ({', '.join(above.cohort)})" if len(above) else ""))
print(f"\nwrote {DD}/mtor_protein_normalised.csv and _sites.csv")
