#!/usr/bin/env python3
"""tcga_subsample_diag.py — F11 diagnostic (fast, NO permutation): does the fixed
20-PC / ridge=1 pipeline on stride-512 Phikon embeddings collapse (median
incremental R2 far below 0, nearly all genes negative) purely because n is small,
independent of organ or platform? This subsamples the ALREADY-EXTRACTED large-n
external cohorts (TCGA-KIRC, TCGA-LUAD) down to n matching TCGA-GBM (85),
CPTAC2-BRCA (102) and TCGA-PAAD (112), refits PCA+ridge on each subsample and
reports the raw incremental-R2 distribution over ALL tested genes -- the same
diagnostic (median_incr_all, negative_frac_all) that showed BRCA at
median approx -? , negative_frac 0.991 and TCGA-GBM at negative_frac 0.99.
No permutation is needed for this: the diagnostic is the raw increment, not its
significance, so this is cheap (seconds per draw) and repeated many times per n
for a stable read.

Interpretation:
  - If LUAD/KIRC collapse the same way at matched n: it is an n-effect of this
    exact pipeline (stride-512 Phikon + fixed 20 PCs), independent of organ or
    platform, and BRCA's null is uninformative rather than biological -- exactly
    like TCGA-GBM.
  - If they do NOT collapse at matched n: BRCA's null needs a different
    explanation (Mertins 2016 QC filtering, iTRAQ ratio compression, a data
    problem specific to the PDC quant matrix -- see review item F31), not a
    generic capacity ceiling.

Inputs: reuses tcga_<short>/{rppa_<short>.csv, rna_case_file.tsv, rna/,
emb_phikon/*.pt} exactly as written by tcga_replicate.py + downloads +
extract_phikon.py. Unlike tcga_residual.py this does NOT read
MORPHO_WSI_EMB_DIR/MORPHO_OUT from the environment -- paths are fixed to the
named project's own directory, to avoid the stale-env-var contamination this
project has hit twice already.

Run:  python tcga_subsample_diag.py --project TCGA-LUAD
      python tcga_subsample_diag.py --project TCGA-KIRC
Output: tcga_<short>/results/subsample_diag_<short>.csv (one row per n x repeat,
plus a full-n baseline row), and a running summary printed to stdout.
"""
import argparse, os, glob, sys
import numpy as np, pandas as pd, torch
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

CV_FOLDS, RIDGE_ALPHA = 5, 1.0
N_REPEATS = 20
CANDIDATE_N = [85, 103, 112]  # match TCGA-GBM / CPTAC2-BRCA / TCGA-PAAD


def tcga_case(s):
    return "-".join(str(s).split("-")[:3])


def cv_r2(X, y, folds=CV_FOLDS, seed=0):
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    preds = np.zeros_like(y, dtype=float)
    I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        Xtr, Xte, ytr = X[tr], X[te], y[tr]
        xm = Xtr.mean(0); xs = Xtr.std(0); xs[xs == 0] = 1.0; ym = ytr.mean()
        A = (Xtr - xm) / xs
        w = np.linalg.solve(A.T @ A + I, A.T @ (ytr - ym))
        preds[te] = ((Xte - xm) / xs) @ w + ym
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - np.sum((y - preds) ** 2) / ss_tot if ss_tot > 0 else 0.0


def load_data(proj, base):
    short = proj.split("-")[1].lower()
    ROOT = f"{base}/tcga_{short}"
    EMB = f"{ROOT}/emb_phikon"

    rppa = pd.read_csv(f"{ROOT}/rppa_{short}.csv", index_col=0)
    rppa.index = [tcga_case(c) for c in rppa.index]
    rppa = rppa[~rppa.index.duplicated()]

    xw = pd.read_csv(f"{ROOT}/rna_case_file.tsv", sep="\t")
    rows = {}
    for _, r in xw.iterrows():
        fp = f"{ROOT}/rna/{r['filename']}"
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp, sep="\t", comment="#")
        d.columns = [c.strip() for c in d.columns]
        d = d[["gene_name", "tpm_unstranded"]].dropna()
        d = d[~d["gene_name"].astype(str).str.startswith("N_")]
        rows.setdefault(tcga_case(r["case"]), []).append(d.groupby("gene_name")["tpm_unstranded"].mean())
    rna = np.log2(pd.DataFrame({c: pd.concat(v, axis=1).mean(axis=1) for c, v in rows.items()}).T.fillna(0) + 1)

    emb = {}
    for pt in glob.glob(f"{EMB}/*.pt"):
        e = torch.load(pt, map_location="cpu")["embeddings"].float().mean(0).numpy()
        emb.setdefault(tcga_case(os.path.basename(pt)), []).append(e)
    wsi = pd.DataFrame({c: np.mean(v, axis=0) for c, v in emb.items()}).T
    if wsi.empty:
        sys.exit(f"no embeddings in {EMB}")

    common = sorted(set(rppa.index) & set(rna.index) & set(wsi.index))
    rppa, rna, wsi = rppa.loc[common], rna.loc[common], wsi.loc[common]
    genes = [g for g in rppa.columns if g in rna.columns]
    return rppa, rna, wsi, genes, common


def run_one(rppa, rna, wsi, genes, idx, min_n=30):
    sub_wsi = wsi.iloc[idx]
    pcs = PCA(n_components=20, svd_solver="full").fit_transform(sub_wsi.values)
    incrs = []
    for g in genes:
        y = np.asarray(rppa[g].values, float)[idx]
        v = ~np.isnan(y)
        if v.sum() < min_n:
            continue
        yv = y[v]
        xr = np.asarray(rna[g].values, float)[idx][v].reshape(-1, 1)
        Wv = pcs[v]
        r2r = cv_r2(xr, yv)
        r2b = cv_r2(np.hstack([xr, Wv]), yv)
        incrs.append(r2b - r2r)
    incrs = np.array(incrs)
    return dict(
        n=len(idx), n_tested=len(incrs),
        median_incr_all=float(np.median(incrs)) if len(incrs) else float("nan"),
        negative_frac_all=float((incrs < 0).mean()) if len(incrs) else float("nan"),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True, help="e.g. TCGA-LUAD or TCGA-KIRC (needs a large n cohort)")
    ap.add_argument("--base", default="/public/home/fjhui/ZW")
    ap.add_argument("--repeats", type=int, default=N_REPEATS)
    ap.add_argument("--ns", default=",".join(map(str, CANDIDATE_N)))
    a = ap.parse_args()
    short = a.project.split("-")[1].lower()
    OUT = f"{a.base}/tcga_{short}/results"
    os.makedirs(OUT, exist_ok=True)
    ns = [int(x) for x in a.ns.split(",")]

    rppa, rna, wsi, genes, common = load_data(a.project, a.base)
    N = len(common)
    print(f"[{a.project}] full n={N}, testing {len(genes)} genes, subsample sizes {ns} x {a.repeats} repeats", flush=True)
    if min(ns) >= N:
        sys.exit(f"full n={N} is not larger than the smallest requested subsample size {min(ns)} -- "
                  f"this project has no headroom to subsample from; use a larger cohort (TCGA-LUAD or TCGA-KIRC)")

    rows = []
    full = run_one(rppa, rna, wsi, genes, np.arange(N))
    full["repeat"] = -1
    full["label"] = "full_n"
    rows.append(full)
    print(f"  full-n baseline: n={full['n']} n_tested={full['n_tested']} "
          f"median_incr_all={full['median_incr_all']:.4f} negative_frac_all={full['negative_frac_all']:.4f}", flush=True)

    for n in ns:
        if n >= N:
            print(f"  skip n={n} (>= full n={N})", flush=True)
            continue
        for rep in range(a.repeats):
            rng = np.random.RandomState(1000 * n + rep)
            idx = rng.choice(N, size=n, replace=False)
            res = run_one(rppa, rna, wsi, genes, idx)
            res["repeat"] = rep
            res["label"] = f"n{n}"
            rows.append(res)
        sub = [r for r in rows if r["label"] == f"n{n}"]
        med = float(np.median([r["median_incr_all"] for r in sub]))
        negf = float(np.median([r["negative_frac_all"] for r in sub]))
        negf_lo = float(np.min([r["negative_frac_all"] for r in sub]))
        negf_hi = float(np.max([r["negative_frac_all"] for r in sub]))
        print(f"  n={n}: across {a.repeats} draws -> median(median_incr_all)={med:.4f}  "
              f"median(negative_frac_all)={negf:.4f}  range=[{negf_lo:.4f}, {negf_hi:.4f}]", flush=True)

    df = pd.DataFrame(rows)
    outfp = f"{OUT}/subsample_diag_{short}.csv"
    df.to_csv(outfp, index=False)
    print(f"\n-> {outfp}")
    print("\nCompare median(negative_frac_all) at each n against: TCGA-GBM (n=85) actual = 0.99, "
          "CPTAC2-BRCA (n=102) actual = 0.9912.")


if __name__ == "__main__":
    main()
