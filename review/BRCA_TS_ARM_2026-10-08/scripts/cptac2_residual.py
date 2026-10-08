#!/usr/bin/env python3
"""cptac2_residual.py — same estimand as tcga_residual.py / the CPTAC-3 discovery
analysis (5-fold OOF, closed-form ridge lambda=1, 1000 patient<->slide permutations,
BH-FDR), applied to the CPTAC-2 mass-spec protein matrix instead of RPPA. There is
no CPTAC-3 discovery cohort for these organs (COAD/BRCA/OV are new), so this reports
the base morphology-predictable rate only; run cptac2_enrichment.py afterwards to
test whether the significant genes converge on the same translation / ER-secretion /
splicing / folding families as the five CPTAC-3 discovery organs.

Inputs under /public/home/fjhui/ZW/cptac2_<short>/ (from cptac2_replicate.py + downloads):
  protein_<short>.csv   rna_case_file.tsv   rna/   emb_phikon/*.pt (or MORPHO_WSI_EMB_DIR)
Run:  python -u cptac2_residual.py --project TCGA-BRCA
Outputs: results/<short>_protein_results.csv + results/summary_<short>.csv
PATCH v2 (2026-10-08): also results/<short>_incr_matrix[_ts][_tilepool].npz = genes x (1+B) incremental-R2 matrix
(col 0 observed, cols 1..B the published permutations RandomState(0..B-1)) + per-gene table, for local set-level AUC.
PATCH v3 (2026-10-08, adversarial review M5/M6): ts arm asserts every kept embedding is a TS/BS/MS slide, prints the kept
count and, with --expect-kept N, exits unless exactly N are kept; the npz also records n_slides. Default dx path unchanged.

WARNING: the CPTAC-2 protein matrix covers ~8.6k-10.6k genes vs RPPA's ~360. At
~2s/gene (1000 permutations x closed-form ridge) this is a multi-hour CPU job --
do NOT run inline on the mn02 login node. Submit it (bsub, or a queued batch
job) or nohup it into the background; do not block a shared login shell.
"""
import argparse, os, glob, sys
import numpy as np, pandas as pd, torch
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

CV_FOLDS, N_PERM, RIDGE_ALPHA = 5, 1000, 1.0


def tcga_case(s): return "-".join(str(s).split("-")[:3])


def cv_r2(X, y, folds=CV_FOLDS):
    kf = KFold(n_splits=folds, shuffle=True, random_state=0)
    preds = np.zeros_like(y, dtype=float); I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in kf.split(X):
        Xtr, Xte, ytr = X[tr], X[te], y[tr]
        xm = Xtr.mean(0); xs = Xtr.std(0); xs[xs == 0] = 1.0; ym = ytr.mean()
        A = (Xtr - xm) / xs
        w = np.linalg.solve(A.T @ A + I, A.T @ (ytr - ym))
        preds[te] = ((Xte - xm) / xs) @ w + ym
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - np.sum((y - preds) ** 2) / ss_tot if ss_tot > 0 else 0.0


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--base", default="/public/home/fjhui/ZW")
    ap.add_argument("--min-n", type=int, default=30)
    # PATCH 2026-10-08 (review/BRCA_TS_ARM_2026-10-08). Defaults reproduce the published DX arm exactly.
    ap.add_argument("--slide-type", choices=["dx", "ts"], default="dx",
                    help="dx (default) = emb_phikon -> results/, no filename filter [published arm]; "
                         "ts = emb_phikon_ts -> results_ts/, only embeddings whose barcode sample code is 01x")
    ap.add_argument("--pool", choices=["slide", "tile"], default="slide",
                    help="slide (default, as published) = mean over tiles within a slide, then unweighted mean over a "
                         "case's slides [lines 77-80 of the original]; tile = mean over all tiles of all of a case's "
                         "slides (the CPTAC-3 discovery convention, residual_analysis.py:168-170)")
    ap.add_argument("--case-list", default="",
                    help="optional file, one TCGA case per line (first line may be a header 'case'); restricts `common` "
                         "to these cases, e.g. inventory_dx_cases.tsv to pin the DX-arm 102")
    ap.add_argument("--write-pool-manifest", action="store_true",
                    help="write pool_manifest_<short>[_ts].csv (case, n_slides, n_tiles, slide stems); implied by --slide-type ts")
    # PATCH v2 2026-10-08 (set-level read-out). Saving the matrix changes no number: the permutation loop is the same loop.
    ap.add_argument("--no-save-matrix", action="store_true",
                    help="do NOT write the per-gene increment matrix (genes x (1+B)) npz [default: write it]")
    ap.add_argument("--perms", type=int, default=N_PERM,
                    help="SMOKE TESTS ONLY. Number of permutations; the default (1000, seeds RandomState(0..999)) is the "
                         "published setting. Any other value tags every output file with _smoke.")
    ap.add_argument("--max-genes", type=int, default=0,
                    help="SMOKE TESTS ONLY. Test only the first N genes; tags every output file with _smoke.")
    # PATCH v3 2026-10-08 (adversarial review M5/M6): the TS arm must not silently read DX embeddings, and the number of
    # TS slides it keeps must equal the frozen list minus the recorded skips (the wrapper computes that number).
    ap.add_argument("--expect-kept", type=int, default=0,
                    help="ts arm: exit unless exactly this many slide embeddings pass the sample-code filter "
                         "(0 = no check; the LSF wrapper passes 189 minus the recorded skips)")
    a = ap.parse_args()
    smoke = (a.perms != N_PERM) or (a.max_genes > 0)
    proj = a.project; short = proj.split("-")[1].lower()
    ROOT = f"{a.base}/cptac2_{short}"
    ts = a.slide_type == "ts"; sfx = "_ts" if ts else ""
    EMB = os.environ.get("MORPHO_WSI_EMB_DIR", f"{ROOT}/emb_phikon{sfx}")
    OUT = os.environ.get("MORPHO_OUT", f"{ROOT}/results{sfx}"); os.makedirs(OUT, exist_ok=True)
    print(f"[{proj}] slide-type={a.slide_type} pool={a.pool} EMB={EMB} OUT={OUT}", flush=True)

    protein = pd.read_csv(f"{ROOT}/protein_{short}.csv", index_col=0)
    protein.index = [tcga_case(c) for c in protein.index]; protein = protein[~protein.index.duplicated()]

    xw = pd.read_csv(f"{ROOT}/rna_case_file.tsv", sep="\t"); rows = {}
    for _, r in xw.iterrows():
        fp = f"{ROOT}/rna/{r['filename']}"
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp, sep="\t", comment="#"); d.columns = [c.strip() for c in d.columns]
        d = d[["gene_name", "tpm_unstranded"]].dropna()
        d = d[~d["gene_name"].astype(str).str.startswith("N_")]
        rows.setdefault(tcga_case(r["case"]), []).append(d.groupby("gene_name")["tpm_unstranded"].mean())
    rna = np.log2(pd.DataFrame({c: pd.concat(v, axis=1).mean(axis=1) for c, v in rows.items()}).T.fillna(0) + 1)

    emb = {}
    if a.slide_type == "dx" and a.pool == "slide" and not a.write_pool_manifest:
        # ---- original code path, unchanged ----
        for pt in glob.glob(f"{EMB}/*.pt"):
            e = torch.load(pt, map_location="cpu")["embeddings"].float().mean(0).numpy()
            emb.setdefault(tcga_case(os.path.basename(pt)), []).append(e)
        wsi = pd.DataFrame({c: np.mean(v, axis=0) for c, v in emb.items()}).T
        n_kept = sum(len(v) for v in emb.values())            # v3: bookkeeping only (npz metadata)
    else:
        # ---- patched path: sample-code filter (ts), slide- or tile-weighted pooling, pooling manifest ----
        import re
        code_re = re.compile(r"^TCGA-\w\w-\w{4}-(\d\d)[A-Z]")
        label_re = re.compile(r"^TCGA-\w\w-\w{4}-\d\d[A-Z]-\d\d-(TS|BS|MS)")   # frozen Tissue Slides: top / bottom / middle
        per_case, n_skip = {}, 0
        for pt in sorted(glob.glob(f"{EMB}/*.pt")):
            stem = os.path.basename(pt)[:-3]
            m = code_re.match(stem)
            if ts and not (m and m.group(1) == "01"):          # primary solid tumour only; mirrors extraction-time rule
                n_skip += 1; continue
            if ts and not label_re.match(stem):                # M5: a DX slide (label DX1, code 01Z) would pass the 01 filter
                sys.exit(f"FATAL: ts arm kept {stem}, whose slide label is not TS/BS/MS -- EMB={EMB} is not the TS embedding "
                         f"folder (is MORPHO_WSI_EMB_DIR set?)")
            E = torch.load(pt, map_location="cpu")["embeddings"].float()
            d = per_case.setdefault(tcga_case(stem), {"slide_means": [], "tile_sum": 0.0, "n_tiles": 0, "stems": []})
            d["slide_means"].append(E.mean(0).numpy()); d["tile_sum"] = d["tile_sum"] + E.sum(0).numpy()
            d["n_tiles"] += int(E.shape[0]); d["stems"].append(stem)
        n_kept = sum(len(d['stems']) for d in per_case.values())
        print(f"[{proj}] embeddings: {n_kept} kept, {n_skip} skipped by sample-code filter, "
              f"{len(per_case)} cases", flush=True)
        if a.expect_kept and n_kept != a.expect_kept:
            sys.exit(f"FATAL: {n_kept} slide embeddings kept, expected {a.expect_kept} (frozen list minus recorded skips)")
        if a.pool == "slide":
            wsi = pd.DataFrame({c: np.mean(d["slide_means"], axis=0) for c, d in per_case.items()}).T
        else:
            wsi = pd.DataFrame({c: d["tile_sum"] / d["n_tiles"] for c, d in per_case.items()}).T
        if ts or a.write_pool_manifest:
            pm = pd.DataFrame([{"case": c, "n_slides": len(d["stems"]), "n_tiles": d["n_tiles"], "slides": "|".join(d["stems"])}
                               for c, d in sorted(per_case.items())])
            pm.to_csv(f"{OUT}/pool_manifest_{short}{sfx}.csv", index=False)
    if wsi.empty:
        sys.exit(f"no embeddings in {EMB}")

    common = sorted(set(protein.index) & set(rna.index) & set(wsi.index))
    if a.case_list:
        want = [l.strip() for l in open(a.case_list) if l.strip().startswith("TCGA-")]
        missing = sorted(set(want) - set(common))
        print(f"[{proj}] --case-list {a.case_list}: {len(want)} cases requested, {len(missing)} not in common: {missing}", flush=True)
        common = [c for c in common if c in set(want)]
    print(f"[{proj}] cases: protein {len(protein)}, rna {len(rna)}, wsi {len(wsi)}, common {len(common)}", flush=True)
    protein, rna, wsi = protein.loc[common], rna.loc[common], wsi.loc[common]
    pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.values)
    perm = [np.random.RandomState(i).permutation(len(common)) for i in range(a.perms)]   # a.perms == N_PERM unless --perms (smoke)
    genes = [g for g in protein.columns if g in rna.columns]
    if a.max_genes > 0:
        genes = genes[:a.max_genes]
    print(f"[{proj}] testing {len(genes)} proteins present in RNA -- this is a multi-hour CPU job, "
          f"progress prints every 200 genes", flush=True)

    res, null_rows = [], []
    for i, g in enumerate(genes, 1):
        y = np.asarray(protein[g].values, float); v = ~np.isnan(y)
        if v.sum() < a.min_n:
            continue
        yv = y[v]; xr = np.asarray(rna[g].values, float)[v].reshape(-1, 1); Wv = pcs[v]
        r2r = cv_r2(xr, yv); r2b = cv_r2(np.hstack([xr, Wv]), yv); incr = r2b - r2r
        null = np.array([cv_r2(np.hstack([xr, pcs[p][v]]), yv) - r2r for p in perm])
        res.append({"gene": g, "n": int(v.sum()), "r2_rna": r2r, "r2_both": r2b,
                    "incremental_r2": incr, "pval": (np.sum(null >= incr) + 1) / (len(perm) + 1)})
        null_rows.append(null)                                   # PATCH v2: keep the null draws (same array that gave pval)
        if i % 200 == 0:
            print(f"  [{i}/{len(genes)}] ...", flush=True)
    res = pd.DataFrame(res); res["fdr"] = bh_fdr(res["pval"].values)
    tag = sfx + ("_tilepool" if a.pool == "tile" else "")      # "" for the published default -> original file names
    if smoke:
        tag += "_smoke"
    res.to_csv(f"{OUT}/{short}_protein_results{tag}.csv", index=False)

    # ---- PATCH v2: per-arm increment matrix for local set-level read-outs (local_brca_setlevel.py) ----
    # M[g, 0] = observed incremental R^2 of gene g; M[g, 1..B] = the same statistic under permutation b = 1..B
    # (RandomState(b-1).permutation(len(common)) of the patient<->WSI-PC-row correspondence, every gene on the SAME permutations,
    # so gene-gene correlation is preserved column-wise). Row order = row order of the results csv.
    if not a.no_save_matrix:
        M = np.vstack(null_rows) if null_rows else np.zeros((0, len(perm)))
        M = np.hstack([res["incremental_r2"].values.reshape(-1, 1).astype(float), M.astype(float)])
        assert M.shape == (len(res), 1 + len(perm)), M.shape
        assert np.array_equal(M[:, 0], res["incremental_r2"].values)
        p_chk = (np.sum(M[:, 1:] >= M[:, [0]], axis=1) + 1) / (M.shape[1])           # = (#null>=obs + 1)/(B + 1)
        assert np.allclose(p_chk, res["pval"].values, rtol=0, atol=1e-15), "saved matrix does not reproduce pval"
        mpath = f"{OUT}/{short}_incr_matrix{tag}.npz"
        np.savez_compressed(
            mpath,
            genes=np.array(res["gene"].astype(str).values, dtype=str), M=M,
            n=res["n"].values.astype(np.int64), r2_rna=res["r2_rna"].values.astype(float),
            r2_both=res["r2_both"].values.astype(float), pval=res["pval"].values.astype(float),
            fdr=res["fdr"].values.astype(float), cases=np.array(common, dtype=str),
            perm_seeds=np.arange(len(perm), dtype=np.int64), B=np.int64(len(perm)),
            n_pcs=np.int64(20), ridge_alpha=np.float64(RIDGE_ALPHA), cv_folds=np.int64(CV_FOLDS), n_slides=np.int64(n_kept),
            arm=np.array(a.slide_type), pool=np.array(a.pool), project=np.array(proj), smoke=np.bool_(smoke),
            columns=np.array("0=observed; 1..B=permutation b-1 (RandomState(b-1)); rows follow genes"))
        print(f"[{proj}] increment matrix {M.shape[0]} genes x {M.shape[1]} -> {mpath}", flush=True)
    sig = res[(res.fdr < 0.05) & (res.incremental_r2 > 0)]

    summ = {"project": proj, "n_cases": len(common), "tested": len(res),
            "significant": len(sig), "pct_sig": 100 * len(sig) / max(len(res), 1),
            "median_incr_sig": float(sig.incremental_r2.median()) if len(sig) else np.nan,
            "negative_frac_all": float((res.incremental_r2 < 0).mean())}
    if tag:
        summ.update({"slide_type": a.slide_type, "pool": a.pool})
    pd.DataFrame([summ]).to_csv(f"{OUT}/summary_{short}{tag}.csv", index=False)
    print(f"\n==================== CPTAC-2 MS EXTERNAL TEST {proj} ====================")
    for k, v in summ.items():
        print(f" {k}: {v:.2f}" if isinstance(v, float) else f" {k}: {v}")
    print(f" -> {OUT}/{short}_protein_results{tag}.csv , summary_{short}{tag}.csv")
    print(f" next: python cptac2_enrichment.py {proj}")


if __name__ == "__main__":
    main()
