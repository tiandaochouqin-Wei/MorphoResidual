#!/usr/bin/env python3
"""T1-7 item 2: summarise the five-cohort seed-reliability run and use it to
disattenuate the per-gene correlations the manuscript reports.

The per-gene incremental R^2 is a difference of two cross-validated R^2 values on
~100 patients. Before any correlation between two such estimates can be read as
"agreement" or "disagreement", we need to know how well the estimate agrees with
ITSELF -- i.e. its reliability. seed_reliability.py re-ran the identical estimand
under 10 cross-validation seeds (the local ridge reproduced the pipeline's own
scorer exactly, max|diff| = 0), and this script turns those into:

  rel_1   reliability of a single-seed estimate = mean pairwise correlation over
          the 45 seed pairs. This is the ceiling on any correlation between two
          independent, equally noisy estimates of the same quantity.
  rel_10  Spearman-Brown reliability of the 10-seed mean.

and then applies r_true = r_obs / sqrt(rel_1 * rel_2) to the two comparisons the
manuscript makes:
  - encoder agreement (SuppTable_encoders): SAME patients, SAME platform,
    different tile encoder. Both sides are estimated with the same estimator on
    the same data, so rel_1 applies to both.
  - external replication in TCGA-KIRC (SuppTable_replication): DIFFERENT patients,
    different protein platform (RPPA). rel_1 is not available for the KIRC side,
    so we report the bound obtained by assuming the KIRC estimate is noiseless,
    which is the most generous possible correction.

Caveat carried into the manuscript: seed resampling varies only the fold
assignment. Patients, PCA basis, platform and encoder are all held fixed, so this
isolates ONE variance component and is an UPPER bound on true test-retest
reliability. Disattenuated correlations computed from it are correspondingly
LOWER bounds.

Writes figures/figdata/seed_reliability_summary.csv
"""
import numpy as np
import pandas as pd

DD = "figures/figdata"
COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
K = 10
COLS = [f"incr_s{i}" for i in range(K)]

# per-gene incremental-R^2 correlation with UNI, from SuppTable_encoders.tex
ENC = {"ccrcc": [("Phikon", 0.808), ("Phikon-v2", 0.585), ("ResNet-50", 0.681)],
       "luad": [("Phikon", 0.780)], "ucec": [("Phikon", 0.845)],
       "gbm": [("Phikon", 0.669)], "pdac": [("Phikon", 0.854)]}
KIRC_R = -0.066          # CPTAC-CCRCC vs TCGA-KIRC, all genes tested in both


def rel(M):
    """Single-measurement reliability = mean pairwise Pearson r over seed pairs."""
    P = np.corrcoef(M.T)
    iu = np.triu_indices(M.shape[1], 1)
    return float(P[iu].mean()), float(np.median(P[iu])), float(P[iu].min()), float(P[iu].max())


def icc21(M):
    """ICC(2,1) cross-check: two-way random, absolute agreement, single measure."""
    n, k = M.shape
    gm = M.mean()
    ms_r = k * ((M.mean(1) - gm) ** 2).sum() / (n - 1)
    ms_c = n * ((M.mean(0) - gm) ** 2).sum() / (k - 1)
    ms_e = ((M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + gm) ** 2).sum() / ((n - 1) * (k - 1))
    return float((ms_r - ms_e) / (ms_r + (k - 1) * ms_e + k * (ms_c - ms_e) / n))


rows = []
print(f"{'cohort':7s} {'set':10s} {'genes':>6s} {'rel_1':>7s} {'rel_10':>7s} {'ICC':>6s} "
      f"{'med incr':>16s} {'incr>0 per seed':>18s}")
for c in COH:
    d = pd.read_csv(f"{DD}/seed_reliability_{c}.csv")
    for lab, sub in [("all tested", d), ("significant", d[d.pinned_significant])]:
        M = sub[COLS].values
        r1, rmed, rlo, rhi = rel(M)
        r10 = K * r1 / (1 + (K - 1) * r1)
        med = np.nanmedian(M, 0)
        pos = (M > 0).sum(0)
        print(f"{c.upper():7s} {lab:10s} {len(sub):6d} {r1:7.3f} {r10:7.3f} {icc21(M):6.3f} "
              f"{med.mean():8.4f}+/-{med.std():.4f} {pos.mean():8.0f}+/-{pos.std():3.0f} "
              f"({100*pos.mean()/len(sub):4.1f}%)")
        rows.append(dict(cohort=c.upper(), gene_set=lab, n_genes=len(sub),
                         rel_single_seed=r1, rel_pairwise_median=rmed,
                         rel_pairwise_min=rlo, rel_pairwise_max=rhi,
                         rel_10seed_mean=r10, icc21=icc21(M),
                         median_incr_mean=med.mean(), median_incr_sd=med.std(),
                         n_incr_pos_mean=pos.mean(), n_incr_pos_sd=pos.std(),
                         frac_incr_pos=pos.mean() / len(sub)))

res = pd.DataFrame(rows)
res.to_csv(f"{DD}/seed_reliability_summary.csv", index=False)

A = res[res.gene_set == "all tested"].set_index("cohort")
print("\n" + "=" * 78)
print("DISATTENUATION 1 -- encoder agreement (same patients, same platform)")
print("  both sides carry the same estimator noise, so r_true = r_obs / rel_1")
for c in COH:
    r1 = A.loc[c.upper(), "rel_single_seed"]
    for enc, r in ENC[c]:
        dis = r / r1
        flag = "  <- at the ceiling" if dis > 0.95 else ""
        print(f"  {c.upper():6s} UNI vs {enc:12s} r_obs={r:.3f}  rel_1={r1:.3f}  "
              f"r_true>={min(dis,1.0):.3f}{' (>1, i.e. at ceiling)' if dis>1 else ''}{flag}")

print("\nDISATTENUATION 2 -- external replication (different patients, RPPA platform)")
r1 = A.loc["CCRCC", "rel_single_seed"]
print(f"  CPTAC-CCRCC vs TCGA-KIRC  r_obs={KIRC_R:+.3f}")
print(f"    assuming the KIRC estimate is NOISELESS (most generous): "
      f"r_true = {KIRC_R:+.3f}/sqrt({r1:.3f}) = {KIRC_R/np.sqrt(r1):+.3f}")
print(f"    assuming it is equally noisy:                          "
      f"r_true = {KIRC_R:+.3f}/{r1:.3f} = {KIRC_R/r1:+.3f}")
print("    -> the null survives every correction; it is not an artefact of estimator noise.")

print("\nRANGES FOR THE MANUSCRIPT")
for lab in ["all tested", "significant"]:
    g = res[res.gene_set == lab]
    print(f"  {lab:11s}: rel_1 {g.rel_single_seed.min():.2f}-{g.rel_single_seed.max():.2f}, "
          f"rel_10 {g.rel_10seed_mean.min():.3f}-{g.rel_10seed_mean.max():.3f}, "
          f"median incr {g.median_incr_mean.min():+.4f}..{g.median_incr_mean.max():+.4f} "
          f"(seed SD {g.median_incr_sd.min():.4f}-{g.median_incr_sd.max():.4f})")
g = res[res.gene_set == "significant"]
print(f"  fraction of pinned-significant genes with a positive increment under a random "
      f"reseed: {100*g.frac_incr_pos.min():.0f}-{100*g.frac_incr_pos.max():.0f}%")
print(f"\nwrote {DD}/seed_reliability_summary.csv")
