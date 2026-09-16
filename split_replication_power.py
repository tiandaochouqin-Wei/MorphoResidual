#!/usr/bin/env python3
"""split_replication_power.py -- feasibility pre-check for a SET-LEVEL replication
design, run before spending about a week on the CPTAC-3 confirmatory UCEC cohort
(n=60 tri-modal cases).

Why. Per-gene re-discovery in a new cohort (per-gene permutation + BH-FDR) has no
power at n of about 60-110: discovery_subsample_diag.py and tcga_subsample_diag.py
showed the increment distribution collapsing there. A replication does not need
per-gene significance. It needs one answer: do the genes that morphology predicts
in one set of patients also carry higher increments, AS A SET, in independent
patients? This script estimates the power of that set-level test by splitting a
single discovery cohort into two disjoint halves, many times:

  selection half (n_sel = N - n_test)
      PCA(K_sel) on this half's own WSI embeddings; raw incremental R^2 for every
      gene; the top --top-frac genes form the "selected" set.
  test half (n_test, default 60 = real C1-UCEC n)
      PCA(K_test) on this half's own WSI embeddings; incremental R^2 for every gene.
      Primary statistic   : AUC = P(increment of a selected gene > increment of a
                            background gene).
      Secondary statistic : Spearman rho between selection-half and test-half
                            increments over all usable genes.
      Null                : permute the patient<->slide rows of the test-half PCs B
                            times and recompute EVERY gene's increment and both
                            statistics. Genes stay together, so gene-gene
                            correlation is preserved and the p-value is not
                            anti-conservative. One-sided (upper tail).
  Model capacity scales with n exactly like the discovery design (20 PCs at n=100):
      K = round(pcs_per_100 * n / 100), so n=60 -> 12, n=77 -> 15.

Conservative on purpose. The real C1 selection would be the permutation-FDR set from
the FULL discovery cohort (n=100); here selection is a raw-increment top-k on fewer
patients, which is noisier. PDAC is the weakest-signal discovery cohort, so running
on PDAC under-states rather than over-states the power UCEC would have.

Reading, fixed before running:
  >= 80% of splits with p_auc < 0.05  -> set-level design powered at n_test: run C1-UCEC
  50-80%                              -> marginal
  <  50%                              -> under-powered: do not run C1

Run inside an LSF job so the MORPHO_* variables come from a clean shell
(residual_analysis.py reads them at import time):

  smoke test (minutes; checks it runs and gives a timing):
    bsub -q smp -n 16 -R "span[hosts=1]" -o splitpow_pdac_smoke.out \
      "bash -c 'source morpho_env.sh pdac && python -u split_replication_power.py --workers 16 --splits 1 --perms 10 --tag smoke'"
  full run:
    bsub -q smp -n 16 -R "span[hosts=1]" -o splitpow_pdac.out \
      "bash -c 'source morpho_env.sh pdac && python -u split_replication_power.py --workers 16'"

Output: <MORPHO_OUT>/split_replication_power_<cohort>[_<tag>].csv, one row per split,
plus a summary on stdout.
"""
import argparse
import math
import multiprocessing as mp
import os
import sys
import time
import zlib

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import residual_analysis as RA  # noqa: E402  (same loaders and cv_r2 as the main analysis)

_W = {}


def _init(rna, prot, pcs, perms, min_n):
    _W.update(rna=rna, prot=prot, pcs=pcs, perms=perms, min_n=min_n)


def _gene_sel(j):
    y = _W["prot"][:, j]
    v = ~np.isnan(y)
    if v.sum() < _W["min_n"]:
        return j, np.nan
    yv = y[v]
    xr = _W["rna"][v, j].reshape(-1, 1)
    pcs = _W["pcs"]
    return j, RA.cv_r2(np.hstack([xr, pcs[v]]), yv, RA.CV_FOLDS) - RA.cv_r2(xr, yv, RA.CV_FOLDS)


def _gene_test(j):
    y = _W["prot"][:, j]
    v = ~np.isnan(y)
    if v.sum() < _W["min_n"]:
        return j, None
    yv = y[v]
    xr = _W["rna"][v, j].reshape(-1, 1)
    pcs = _W["pcs"]
    r2r = RA.cv_r2(xr, yv, RA.CV_FOLDS)
    out = np.empty(len(_W["perms"]) + 1)
    out[0] = RA.cv_r2(np.hstack([xr, pcs[v]]), yv, RA.CV_FOLDS) - r2r
    for b, p in enumerate(_W["perms"], 1):
        out[b] = RA.cv_r2(np.hstack([xr, pcs[p][v]]), yv, RA.CV_FOLDS) - r2r
    return j, out


def run_pool(fn, n_genes, workers, rna, prot, pcs, perms, min_n):
    with mp.Pool(workers, initializer=_init, initargs=(rna, prot, pcs, perms, min_n)) as pool:
        return dict(pool.imap_unordered(fn, range(n_genes), chunksize=64))


def k_for(n, per100):
    return max(2, int(round(per100 * n / 100.0)))


def cohort_name():
    name = os.path.basename(str(RA.ROOT).rstrip("/"))
    return "ccrcc" if name == "ZW" else name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-test", type=int, default=60)
    ap.add_argument("--splits", type=int, default=10)
    ap.add_argument("--perms", type=int, default=200)
    ap.add_argument("--top-frac", type=float, default=0.15)
    ap.add_argument("--pcs-per-100", type=float, default=20.0)
    ap.add_argument("--min-n", type=int, default=30)
    ap.add_argument("--workers", type=int, default=RA.N_WORKERS)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="")
    ap.add_argument("--group-csv", default="",
                     help="CSV with a case_submitter_id column plus one or more group "
                          "columns valued 'test'/'sel'. When given (with --group-col), "
                          "replaces the random patient split with this FIXED, "
                          "batch/operator-disjoint assignment -- one deterministic split "
                          "per invocation, ignoring --splits and --n-test.")
    ap.add_argument("--group-col", default="",
                     help="which column of --group-csv to use as the test/sel assignment "
                          "for this run (e.g. q_early_vs_late, op_balanced_0)")
    a = ap.parse_args()

    cohort = cohort_name()
    print("=== resolved MORPHO_* paths (verify this is the cohort you sourced) ===")
    print(f"  cohort (from root) : {cohort}")
    print(f"  MORPHO_ROOT        : {RA.ROOT}")
    print(f"  MORPHO_PROTEIN     : {RA.PROTEIN_TSV}")
    print(f"  MORPHO_WSI_EMB     : {RA.WSI_EMB_DIR}")
    print(f"  MORPHO_OUT         : {RA.OUT_DIR}")
    print("=" * 70, flush=True)

    rna = RA.load_rna_matrix()
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    genes = sorted(set(rna.columns) & set(protein.columns))

    grouped = bool(a.group_csv and a.group_col)
    if grouped:
        gdf = pd.read_csv(a.group_csv)
        gdf = gdf.set_index("case_submitter_id")[a.group_col]
        gdf = gdf[gdf.index.isin(common)]
        dropped = sorted(set(common) - set(gdf.index))
        if dropped:
            print(f"  NOTE: {len(dropped)} of {len(common)} analysis cases have no group label "
                  f"in {a.group_csv}:{a.group_col} and are dropped: {dropped[:10]}"
                  f"{' ...' if len(dropped) > 10 else ''}", flush=True)
        common = sorted(gdf.index)
        split_specs = [(a.group_col, gdf)]
    else:
        split_specs = None  # random splits computed per-iteration below

    R = rna.loc[common, genes].values.astype(float)
    P = protein.loc[common, genes].values.astype(float)
    Wraw = wsi.loc[common].values.astype(float)
    N, G = len(common), len(genes)

    if grouped:
        label, gser = split_specs[0]
        case_pos = {c: i for i, c in enumerate(common)}
        test_idx = np.array([case_pos[c] for c in gser.index if gser[c] == "test"])
        sel_idx = np.array([case_pos[c] for c in gser.index if gser[c] == "sel"])
        n_test, n_sel = len(test_idx), len(sel_idx)
        if n_sel < 40 or n_test < 30:
            sys.exit(f"grouped split '{label}': n_sel={n_sel}, n_test={n_test} -- too small")
        splits_iter = [(label, test_idx, sel_idx)]
        print(f"\n[{cohort}] GROUPED split '{label}': N={N}, genes={G}; "
              f"n_sel={n_sel}, n_test={n_test}; perms={a.perms}, top_frac={a.top_frac}, "
              f"workers={a.workers} (--splits and --n-test ignored)", flush=True)
    else:
        n_test = a.n_test
        n_sel = N - n_test
        if n_sel < 40:
            sys.exit(f"N={N}: a test half of {n_test} leaves only {n_sel} patients for selection (need >= 40)")
        print(f"\n[{cohort}] N={N}, genes={G}; n_sel={n_sel}, n_test={n_test}; "
              f"splits={a.splits}, perms={a.perms}, top_frac={a.top_frac}, workers={a.workers}", flush=True)
        splits_iter = []
        for s in range(a.splits):
            rng0 = np.random.RandomState(a.seed * 1000 + s)
            order = rng0.permutation(N)
            splits_iter.append((s, order[:n_test], order[n_test:]))

    rows = []
    for label, test_idx, sel_idx in splits_iter:
        s = label
        t0 = time.time()
        # zlib.crc32 (not Python's hash()) so this seed is stable across processes/runs --
        # hash() on a str is randomised per-process (PYTHONHASHSEED) and would make the
        # permutation null irreproducible.
        rng = np.random.RandomState(zlib.crc32(f"perm|{cohort}|{label}".encode()))
        n_test, n_sel = len(test_idx), len(sel_idx)
        k_sel, k_test = k_for(n_sel, a.pcs_per_100), k_for(n_test, a.pcs_per_100)

        pcs_sel = PCA(n_components=k_sel, svd_solver="full").fit_transform(Wraw[sel_idx])
        res_sel = run_pool(_gene_sel, G, a.workers, R[sel_idx], P[sel_idx], pcs_sel, [], a.min_n)
        incr_sel = np.array([res_sel[j] for j in range(G)])

        pcs_test = PCA(n_components=k_test, svd_solver="full").fit_transform(Wraw[test_idx])
        perms = [rng.permutation(n_test) for _ in range(a.perms)]
        res_test = run_pool(_gene_test, G, a.workers, R[test_idx], P[test_idx], pcs_test, perms, a.min_n)

        usable = [j for j in range(G) if np.isfinite(incr_sel[j]) and res_test[j] is not None]
        if len(usable) < 100:
            print(f"  split {s}: only {len(usable)} usable genes, skipped", flush=True)
            continue
        xs = incr_sel[usable]
        M = np.vstack([res_test[j] for j in usable])  # usable genes x (1 + perms); column 0 observed

        n_selected = int(math.ceil(a.top_frac * len(usable)))
        sel_mask = np.zeros(len(usable), dtype=bool)
        sel_mask[np.argsort(-xs)[:n_selected]] = True
        n1, n0 = sel_mask.sum(), (~sel_mask).sum()

        ranks = rankdata(M, axis=0)
        auc = (ranks[sel_mask].sum(axis=0) - n1 * (n1 + 1) / 2.0) / (n1 * n0)

        rs = rankdata(xs)
        rs = (rs - rs.mean()) / rs.std()
        rc = (ranks - ranks.mean(axis=0)) / ranks.std(axis=0)
        rho = (rs[:, None] * rc).mean(axis=0)

        auc_obs, auc_null = auc[0], auc[1:]
        rho_obs, rho_null = rho[0], rho[1:]
        p_auc = (np.sum(auc_null >= auc_obs) + 1) / (len(auc_null) + 1)
        p_rho = (np.sum(rho_null >= rho_obs) + 1) / (len(rho_null) + 1)
        sd_auc = auc_null.std()
        row = dict(split=s, n_sel=n_sel, k_sel=k_sel, n_test=n_test, k_test=k_test,
                   n_usable=len(usable), n_selected=int(n1),
                   auc_obs=float(auc_obs), auc_null_mean=float(auc_null.mean()), auc_null_sd=float(sd_auc),
                   z_auc=float((auc_obs - auc_null.mean()) / sd_auc) if sd_auc > 0 else float("nan"),
                   p_auc=float(p_auc),
                   rho_obs=float(rho_obs), rho_null_mean=float(rho_null.mean()), p_rho=float(p_rho),
                   test_negative_frac=float((M[:, 0] < 0).mean()),
                   minutes=round((time.time() - t0) / 60.0, 1))
        rows.append(row)
        print(f"  split {s}: AUC={auc_obs:.4f} (null {auc_null.mean():.4f}+/-{sd_auc:.4f}) p={p_auc:.4f} | "
              f"rho={rho_obs:.4f} p={p_rho:.4f} | usable={len(usable)} selected={n1} | "
              f"test neg_frac={row['test_negative_frac']:.3f} | {row['minutes']} min", flush=True)

    if not rows:
        sys.exit("no split produced a result")
    df = pd.DataFrame(rows)
    suffix = f"_{a.tag}" if a.tag else ""
    outfp = RA.OUT_DIR / f"split_replication_power_{cohort}{suffix}.csv"
    df.to_csv(outfp, index=False)

    power_auc = float((df.p_auc < 0.05).mean())
    power_rho = float((df.p_rho < 0.05).mean())
    print(f"\n==================== SET-LEVEL REPLICATION POWER [{cohort}] ====================")
    print(f"  splits completed        : {len(df)}")
    print(f"  n_sel / n_test          : {n_sel} / {n_test}   (K {k_sel} / {k_test})")
    print(f"  median AUC obs / null   : {df.auc_obs.median():.4f} / {df.auc_null_mean.median():.4f}")
    print(f"  median z (AUC)          : {df.z_auc.median():.2f}")
    print(f"  power (p_auc < 0.05)    : {power_auc:.0%}   <- primary")
    print(f"  power (p_rho < 0.05)    : {power_rho:.0%}   (secondary)")
    if power_auc >= 0.80:
        verdict = "POWERED: the set-level design detects replication at this n_test -> C1-UCEC worth running"
    elif power_auc >= 0.50:
        verdict = "MARGINAL: detectable in some splits only -> C1-UCEC is a gamble"
    else:
        verdict = "UNDER-POWERED: do not run C1-UCEC on this design"
    print(f"  pre-stated reading      : {verdict}")
    if len(df) < a.splits:
        print(f"  NOTE: only {len(df)} of {a.splits} splits completed")
    print(f"  -> {outfp}")


if __name__ == "__main__":
    main()
