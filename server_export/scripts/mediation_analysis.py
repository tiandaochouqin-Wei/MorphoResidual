#!/usr/bin/env python3
"""
Mediation arm C (plan.md contribution C), CONSERVATIVE rebuild. For proteins that
are morphology-predictable AND downstream of a recurrent ccRCC somatic driver,
test whether tumor morphology statistically MEDIATES the mutation -> protein
(post-transcriptional residual) effect.

  X = somatic driver mutation status (0/1), from WES MAF (VHL/PBRM1/BAP1/SETD2)
  M = morphology score = out-of-fold WSI prediction of protein g's mRNA-residual
  Y = protein g's mRNA-residual (post-transcriptional layer; avoids mut->mRNA path)

Mediation (Baron-Kenny): a: M~X ; b,c': Y~X+M ; ACME=a*b ; ADE=c' ;
prop_mediated=ACME/(ACME+ADE).

WHY THE PRIMARY TEST IS A PERMUTATION ON MUTATION LABELS, not the ACME bootstrap:
we test mediation ONLY on morphology-predictable proteins, so the M->Y link (b)
is strong BY SELECTION and any M-Y-confounding sensitivity (incl. Imai) trivially
looks robust. The genuinely uncertain, testable path is a (mutation -> morphology).
Permuting which patients carry the mutation gives an exact, small-n-robust null
for ACME that isolates path a. Bootstrap CI and (corrected) Imai rho* are reported
as secondary/context only.

Conservative choices: MIN_MUT=10 (skip fragile mutation groups); prop_mediated
outside [0,1] flagged 'unstable' (inconsistent mediation); BH-FDR on permutation p.

Thread limits set before numpy import. Parallelized over genes. CPU-only.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import glob
import gzip
import multiprocessing as mp
import sys
import warnings
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent))
import residual_analysis as base
from residual_analysis import ROOT, CV_FOLDS, N_PCS, OUT_DIR
from sklearn.decomposition import PCA

DRIVERS = ["VHL", "PBRM1", "BAP1", "SETD2"]
WES_DIR = ROOT / "omics" / "mutation"
WES_MANIFEST = ROOT / "omics" / "manifest_wes.tsv"
SIG_CSV = Path(os.environ.get("MORPHO_SIG_CSV", str(OUT_DIR / "residual_results_adjusted.csv")))
N_PERM = int(os.environ.get("MORPHO_N_PERM", 2000))
N_BOOT = int(os.environ.get("MORPHO_N_BOOT", 1000))
N_WORKERS = int(os.environ.get("MORPHO_N_WORKERS", 32))
MIN_MUT = int(os.environ.get("MORPHO_MIN_MUT", 10))   # >= this many mutated AND wild-type
FDR_THRESHOLD = 0.05
NONSILENT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
             "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
             "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins"}


def load_driver_mutations(common_cases):
    manifest = pd.read_csv(WES_MANIFEST, sep="\t")
    fname_to_case = dict(zip(manifest["file_name"], manifest["case_submitter_id"]))
    mut = {c: set() for c in common_cases}
    files = glob.glob(str(WES_DIR / "**" / "*.maf.gz"), recursive=True) + \
            glob.glob(str(WES_DIR / "**" / "*.maf"), recursive=True)
    parsed = 0
    for fp in files:
        case = fname_to_case.get(Path(fp).name)
        if case is None or case not in mut:
            continue
        opener = gzip.open if fp.endswith(".gz") else open
        try:
            with opener(fp, "rt") as fh:
                header = None
                for line in fh:
                    if line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if header is None:
                        header = {c: i for i, c in enumerate(cols)}
                        continue
                    if cols[header["Hugo_Symbol"]] in DRIVERS and \
                       cols[header["Variant_Classification"]] in NONSILENT:
                        mut[case].add(cols[header["Hugo_Symbol"]])
            parsed += 1
        except Exception as e:
            print(f"  WARN parse {Path(fp).name}: {e}")
    print(f"parsed {parsed} MAF files")
    for d in DRIVERS:
        print(f"  {d}: mutated in {sum(1 for c in mut if d in mut[c])}/{len(mut)} patients")
    return mut


def oof_predict(y, X, folds=CV_FOLDS):
    X = np.atleast_2d(X)
    if X.shape[0] != len(y):
        X = X.T
    n = len(y)
    idx = np.arange(n)
    rs = np.random.RandomState(0)
    rs.shuffle(idx)
    pred = np.zeros_like(y, dtype=float)
    aI = 1.0 * np.eye(X.shape[1])
    for k in range(folds):
        te = idx[k::folds]
        tr = np.setdiff1d(idx, te)
        Xtr, Xte, ytr = X[tr], X[te], y[tr]
        mu, sd = Xtr.mean(0), Xtr.std(0)
        sd[sd == 0] = 1
        Xs, Xes = (Xtr - mu) / sd, (Xte - mu) / sd
        w = np.linalg.solve(Xs.T @ Xs + aI, Xs.T @ (ytr - ytr.mean()))
        pred[te] = Xes @ w + ytr.mean()
    return pred


def ols(X, y):
    X1 = np.hstack([np.ones((len(X), 1)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return beta


def mediation(x, m, y):
    a = ols(x.reshape(-1, 1), m)[1]
    byx = ols(np.column_stack([x, m]), y)
    b, cprime = byx[2], byx[1]
    acme, ade = a * b, cprime
    total = acme + ade
    prop = acme / total if abs(total) > 1e-9 else np.nan
    return acme, ade, total, prop


def imai_rho_star(x, m, y):
    # Correct linear-model sensitivity: rho* = |b * sigma_M / sigma_Y| (clamped),
    # the M-Y error correlation that would drive ACME to 0. NOTE: because we test
    # only morphology-predictable proteins, b is large by selection so rho* is
    # usually near 1 -- report but do not over-interpret (see module docstring).
    am = ols(x.reshape(-1, 1), m)
    e_m = m - (am[0] + am[1] * x)
    sig_m = e_m.std()
    yb = ols(np.column_stack([x, m]), y)
    b = yb[2]
    e_y = y - (yb[0] + yb[1] * x + yb[2] * m)
    sig_y = e_y.std()
    if sig_y < 1e-9:
        return 1.0
    return min(abs(b * sig_m / sig_y), 1.0)


def bh_fdr(pvals):
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    ranked = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    fdr = np.empty(n)
    fdr[order] = np.clip(ranked, 0, 1)
    return fdr


_W = {}


def _init(protein, rna, wsi_pcs, driver_vec):
    _W["protein"] = protein
    _W["rna"] = rna
    _W["wsi"] = wsi_pcs
    _W["drv"] = driver_vec


def _seed(s):
    return zlib.crc32(s.encode()) & 0x7fffffff   # deterministic per-string seed


def _process_gene(gene):
    y = np.asarray(_W["protein"][gene])
    valid = ~np.isnan(y)
    if valid.sum() < 30:
        return []
    yv = y[valid]
    xr = _W["rna"][gene][valid].reshape(-1, 1)
    resid = yv - oof_predict(yv, xr, CV_FOLDS)               # Y = mRNA-residual
    M = oof_predict(resid, _W["wsi"][valid], CV_FOLDS)       # morphology score
    nv = len(resid)
    # per-gene permutations of the CORRECT length (valid patients may be < n)
    rperm = np.random.RandomState(_seed(gene))
    perms = [rperm.permutation(nv) for _ in range(N_PERM)]
    out = []
    for d in DRIVERS:
        x = _W["drv"][d][valid]
        nmut = int(x.sum())
        if nmut < MIN_MUT or (nv - nmut) < MIN_MUT:
            continue
        acme, ade, total, prop = mediation(x, M, resid)
        # PRIMARY: permute which patients carry the mutation -> null for path a
        null = np.array([mediation(x[p], M, resid)[0] for p in perms])
        pval = (np.sum(np.abs(null) >= abs(acme)) + 1) / (N_PERM + 1)
        # SECONDARY: bootstrap CI + corrected Imai
        rb = np.random.RandomState(_seed(gene + d))
        boot = np.array([mediation(*(lambda s: (x[s], M[s], resid[s]))(
            rb.choice(nv, nv, replace=True)))[0] for _ in range(N_BOOT)])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        rho0 = imai_rho_star(x, M, resid)
        unstable = not (0.0 <= prop <= 1.0) if np.isfinite(prop) else True
        out.append({"protein": gene, "driver": d, "n_mut": nmut, "ACME": acme,
                    "ADE": ade, "prop_mediated": prop, "acme_lo": lo, "acme_hi": hi,
                    "perm_pval": pval, "imai_rho_star": rho0, "unstable": unstable})
    return out


def main():
    rna = base.load_rna_matrix()
    protein = base.load_protein_matrix()
    wsi = base.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    print(f"\ncommon patients: {len(common)}")

    sig = pd.read_csv(SIG_CSV)
    sig_genes = sig[(sig["fdr"] < 0.05) & (sig["incremental_r2"] > 0)]["gene"].tolist()
    genes = [g for g in sig_genes if g in rna.columns and g in protein.columns]
    print(f"morphology-predictable proteins to test: {len(genes)} (from {SIG_CSV.name})")

    rna_c = rna.loc[common]
    protein_c = protein.loc[common]
    wsi_raw = wsi.loc[common].values
    # svd_solver='full' is load-bearing: at this matrix shape sklearn's auto
    # solver is randomized SVD with random_state=None, which made this script
    # non-reproducible (GBM: 699/700/704/701 significant across four identical
    # runs, 90.4% set overlap). 'full' is exact and deterministic.
    pca = PCA(n_components=min(N_PCS, len(common) - 1, wsi_raw.shape[1]),
              svd_solver="full")
    wsi_pcs = pca.fit_transform(wsi_raw)

    mut = load_driver_mutations(common)
    driver_vec = {d: np.array([1.0 if d in mut[c] else 0.0 for c in common]) for d in DRIVERS}
    protein_dict = {g: protein_c[g].values for g in genes}
    rna_dict = {g: rna_c[g].values for g in genes}

    print(f"parallelizing across {N_WORKERS} workers "
          f"(primary=permutation on mutation labels, N_PERM={N_PERM}, MIN_MUT={MIN_MUT})")
    rows = []
    done = 0
    with mp.Pool(N_WORKERS, initializer=_init,
                 initargs=(protein_dict, rna_dict, wsi_pcs, driver_vec)) as pool:
        for r in pool.imap_unordered(_process_gene, genes, chunksize=8):
            rows.extend(r)
            done += 1
            if done % 200 == 0:
                print(f"  {done}/{len(genes)} proteins done...")

    res = pd.DataFrame(rows)
    if len(res) == 0:
        print("no (driver, protein) pairs met MIN_MUT -- check WES parse")
        sys.exit(1)
    res["fdr"] = bh_fdr(res["perm_pval"].values)
    res = res.sort_values("perm_pval")
    out = OUT_DIR / "mediation_results.csv"
    res.to_csv(out, index=False)

    # conservative "significant": permutation FDR<0.05, ACME CI excludes 0, AND stable
    sigm = res[(res["fdr"] < FDR_THRESHOLD) &
               (np.sign(res["acme_lo"]) == np.sign(res["acme_hi"])) &
               (~res["unstable"])]
    print(f"\ntested {len(res)} (driver,protein) pairs")
    print(f"CONSERVATIVE significant mediations "
          f"(permFDR<{FDR_THRESHOLD} AND ACME CI excl 0 AND stable prop_mediated): {len(sigm)}")
    print(f"  (of which flagged unstable and excluded: "
          f"{((res['fdr'] < FDR_THRESHOLD) & (res['unstable'])).sum()})")
    print(f"\ntop 15 by permutation p:")
    show = ["protein", "driver", "n_mut", "ACME", "prop_mediated",
            "perm_pval", "fdr", "imai_rho_star", "unstable"]
    print(res.head(15)[show].to_string(index=False))
    print(f"\nby driver (conservative significant):")
    for d in DRIVERS:
        print(f"  {d}: {len(sigm[sigm.driver == d])}")
    print(f"full results: {out}")


if __name__ == "__main__":
    main()
