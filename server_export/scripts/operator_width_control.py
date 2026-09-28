#!/usr/bin/env python3
"""
operator_width_control.py -- does the morphology increment survive WIDER operator
correction, and if it falls, is that confounding or capacity?

WHY. residual_analysis_sitepack.py dummy-codes only the largest MAX_LEVELS=12 operator
batches (lines 298-310); cases in the smaller batches get an all-zero row and receive no
batch correction at all. That is a deliberate, documented trade-off -- 23 plex dummies at
n=103 once drove the baseline R2 negative -- but it means the published correction covers
12 of 34 operators in LUAD and GBM, 28 in CCRCC/UCEC and 54 in PDAC, leaving a large
minority of cases uncorrected. Separately, recorded race is aliased with the operator
field far beyond chance (permutation p <= 0.001 in CCRCC/LUAD/GBM/PDAC, at the null in
UCEC, and NOT so for sex), and adjusting for race removes about two thirds of the LUAD
increment. Those two facts together raise a question the published run cannot answer:
would a wider operator design remove the increment as well?

WHAT IT DOES. Exactly the published operation -- residualise the WSI principal components
against the batch design, then re-measure the per-gene increment -- repeated at several
design widths, each paired with a control design of the SAME width made of Gaussian noise.
Residualising against w columns removes w dimensions' worth of variance whatever those
columns are, so the noise arm is what the increment would lose from width alone. Only the
gap between the two arms is attributable to the operator axis.

A design is specified as CAP:MIN, the two parameters sitepack uses: dummy-code the CAP
largest batches that have at least MIN cases. RAISING CAP ALONE DOES ALMOST NOTHING --
that was measured, not assumed: with MIN=3 only 13 of LUAD's 34 operators, 11 of GBM's 34
and 18 of PDAC's 54 batches qualify at all, so GBM's published CAP of 12 never even binds.
The parameter that actually controls coverage is MIN, and lowering it to 1 is what puts
every operator in the design -- 34 dummy columns on 105 patients in LUAD, which is exactly
the regime where the noise arm earns its place.

    0        no correction; reproduces the published per-gene increment
    12:3     the published design
    99:3     cap removed, MIN unchanged -- shows how little the cap was binding
    99:2     batches of two admitted
    99:1     every operator gets a column (the widest design the data allow)

    per width w, per gene:
        incr_op_<w>     = R2(mRNA, PCs residualised on the operator design) - R2(mRNA)
        incr_noise_<w>  = R2(mRNA, PCs residualised on w Gaussian columns) - R2(mRNA)
                          averaged over N_NOISE_DRAWS draws

READ IT AS: if incr_op_w tracks incr_noise_w as w grows, the loss is the price of
removing w dimensions and the operator axis carries no extra confounding. If incr_op_w
falls below incr_noise_w, the operator axis is doing real work and the published width
was too narrow.

This does NOT re-run the permutation null, so it produces no new FDR and no new counts;
it is a sensitivity analysis on the estimate, run on the cohort's published significant
set. Re-running sitepack at a larger MAX_LEVELS is the follow-up if this shows a gap.

Env:  the usual MORPHO_* (source morpho_env.sh <cancer>), plus
      MORPHO_SIG_CSV   pinned residual_results_tumoronly.csv
      MORPHO_COHORT    cohort key
      MORPHO_WIDTHS    comma-separated, default "0,12,20,full"
      MORPHO_OPW_OUT   default operator_width_<cohort>.csv
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
from collections import Counter

import numpy as np
import pandas as pd

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

MIN_GROUP_DEFAULT = 3    # as in residual_analysis_sitepack.py
N_NOISE_DRAWS = 3
NOISE_SEED = 0
MIN_PATIENTS = 30
AXIS = os.environ.get("MORPHO_BATCH_AXIS", "operator")


def label(spec):
    return "none" if spec[0] == 0 else f"{spec[0]}:{spec[1]}"


def parse_spec(s):
    """'0' -> no design; 'CAP:MIN' -> (cap, min_group)."""
    s = s.strip()
    if s == "0":
        return (0, 0)
    cap, _, mg = s.partition(":")
    return (int(cap), int(mg) if mg else MIN_GROUP_DEFAULT)


def design(raw, counts, spec):
    """The published dummy rule (sitepack lines 304-310), with both parameters exposed."""
    cap, min_group = spec
    levels = sorted(counts)
    if cap == 0:
        return np.zeros((len(raw), 0)), 0, len(raw)
    dummy_levels = [lv for lv, k in counts.most_common(cap) if k >= min_group]
    if len(dummy_levels) == len(levels) and dummy_levels:
        dummy_levels = dummy_levels[:-1]          # keep the design full-rank
    D = np.zeros((len(raw), len(dummy_levels)))
    for j, lv in enumerate(dummy_levels):
        D[:, j] = (raw == lv).astype(float)
    undummied = int(sum(counts[lv] for lv in levels if lv not in set(dummy_levels)))
    return D, len(dummy_levels), undummied


def residualise(pcs, D):
    if D.shape[1] == 0:
        return pcs
    B = np.hstack([np.ones((len(pcs), 1)), D])
    coef, *_ = np.linalg.lstsq(B, pcs, rcond=None)
    return pcs - B @ coef


def main():
    cohort = os.environ.get("MORPHO_COHORT", "").strip().lower()
    if not cohort:
        sys.exit("[opw] set MORPHO_COHORT")
    sig_csv = os.environ.get("MORPHO_SIG_CSV")
    if not sig_csv:
        sys.exit("[opw] set MORPHO_SIG_CSV")
    out = os.environ.get("MORPHO_OPW_OUT", f"operator_width_{cohort}.csv")
    specs = [parse_spec(s) for s in
             os.environ.get("MORPHO_WIDTHS", "0,12:3,99:3,99:2,99:1").split(",")]

    import batch_leak_check as blc
    paths = blc.cohort_paths(cohort)
    os.environ.update(
        MORPHO_ROOT=paths["root"], MORPHO_RNA_MANIFEST=paths["rna_manifest"],
        MORPHO_ALIQUOT_XWALK=paths["xwalk"], MORPHO_SLIDE_MAP=paths["slide_map"],
        MORPHO_PROTEIN_TSV=paths["protein"], MORPHO_WSI_EMB_DIR=paths["emb"])
    import residual_analysis as ra
    from sklearn.decomposition import PCA

    rna, protein, wsi = ra.load_rna_matrix(), ra.load_protein_matrix(), ra.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    genes_all = sorted(set(rna.columns) & set(protein.columns))
    rna, protein = rna.loc[common, genes_all], protein.loc[common, genes_all]
    n_pcs = min(getattr(ra, "N_PCS", 20), len(common) - 1, wsi.shape[1])
    pcs = PCA(n_components=n_pcs, svd_solver="full").fit_transform(wsi.loc[common].values)
    print(f"[opw] {cohort}: {len(common)} patients; PCs {pcs.shape}", flush=True)

    axes = blc.case_batch_labels(cohort, paths, common)
    if AXIS not in axes:
        sys.exit(f"[opw] axis {AXIS!r} unavailable; have {sorted(axes)}")
    lab, _ = axes[AXIS]
    raw = np.array([lab.get(c, "") or "_missing" for c in common])
    counts = Counter(raw.tolist())
    print(f"[opw] axis {AXIS}: {len(counts)} raw batches, "
          f"sizes {dict(counts.most_common(8))}{' ...' if len(counts) > 8 else ''}", flush=True)

    rng = np.random.RandomState(NOISE_SEED)
    PCSETS, INFO = {}, {}
    for w in specs:
        D, ndum, undum = design(raw, counts, w)
        PCSETS[("op", w)] = residualise(pcs, D)
        kept = (PCSETS[("op", w)].var(axis=0).sum() / max(pcs.var(axis=0).sum(), 1e-12))
        INFO[w] = dict(dummies=ndum, undummied=undum, pc_var_kept=float(kept))
        if ndum > 0:
            PCSETS[("noise", w)] = [residualise(pcs, rng.normal(size=(len(common), ndum)))
                                    for _ in range(N_NOISE_DRAWS)]
        print(f"[opw]   design {label(w):>6s}: {ndum:3d} dummies, {undum:4d}/{len(common)} cases "
              f"left to the intercept, PC variance retained {100 * kept:5.1f}%", flush=True)

    md = pd.read_csv(sig_csv)
    sig = (md[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0)]["gene"].astype(str).tolist()
           if {"fdr", "incremental_r2"}.issubset(md.columns) else md["gene"].astype(str).tolist())
    sig = [g for g in sig if g in rna.columns and g in protein.columns]
    cap = int(os.environ.get("MORPHO_MAX_GENES", "0") or 0)
    if cap:
        sig = sig[:cap]
        print(f"[opw] *** SMOKE TEST: {cap} genes -- NOT the real run ***", flush=True)
    print(f"[opw] testing {len(sig)} genes over designs {[label(w) for w in specs]}", flush=True)

    CV = getattr(ra, "CV_FOLDS", 5)
    rows = []
    for g in sig:
        y = np.asarray(protein[g].values, float)
        v = ~np.isnan(y)
        if v.sum() < MIN_PATIENTS:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1)
        r2_rna = ra.cv_r2(xr, yv, CV)
        row = {"gene": g, "n": int(v.sum())}
        for w in specs:
            row[f"incr_op_{label(w)}"] = ra.cv_r2(np.hstack([xr, PCSETS[("op", w)][v]]), yv, CV) - r2_rna
            if ("noise", w) in PCSETS:
                row[f"incr_noise_{label(w)}"] = float(np.mean(
                    [ra.cv_r2(np.hstack([xr, P[v]]), yv, CV) - r2_rna for P in PCSETS[("noise", w)]]))
        rows.append(row)

    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)
    if res.empty:
        print(f"[opw] nothing to summarise -> {out}")
        return 0

    ref = res[f"incr_op_{label(specs[0])}"].median() if specs[0][0] == 0 else np.nan
    print("\n==================== OPERATOR WIDTH SWEEP ====================")
    print(f" cohort {cohort}: {len(res)} genes, reference (uncorrected) median increment {ref:+.4f}")
    print(f" {'design':>6s} {'dummies':>8s} {'uncorr. cases':>14s} {'PC var':>7s} "
          f"{'median incr':>12s} {'% of ref':>9s} | {'noise incr':>11s} {'noise %':>8s} {'gap':>7s}")
    for w in specs:
        i = INFO[w]
        m = res[f"incr_op_{label(w)}"].median()
        pct = 100 * m / ref if ref == ref and ref != 0 else float("nan")
        if f"incr_noise_{label(w)}" in res:
            nm = res[f"incr_noise_{label(w)}"].median()
            npct = 100 * nm / ref
            print(f" {label(w):>6s} {i['dummies']:8d} {i['undummied']:14d} {100*i['pc_var_kept']:6.1f}% "
                  f"{m:+12.4f} {pct:8.1f}% | {nm:+11.4f} {npct:7.1f}% {npct - pct:+6.1f}")
        else:
            print(f" {label(w):>6s} {i['dummies']:8d} {i['undummied']:14d} {100*i['pc_var_kept']:6.1f}% "
                  f"{m:+12.4f} {pct:8.1f}% | {'--':>11s} {'--':>8s} {'--':>6s}")
    print(" READ IT AS: the operator axis carries confounding only where its median sits")
    print(" clearly BELOW the same-width noise arm. Tracking it means the loss is width.")
    print(f" -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
