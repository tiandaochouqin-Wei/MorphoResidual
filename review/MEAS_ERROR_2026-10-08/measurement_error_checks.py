#!/usr/bin/env python3
"""
measurement_error_checks.py (v2, 2026-10-08) -- EXPORT job for the own-mRNA measurement-error
analysis of the five discovery cohorts (CCRCC, LUAD, UCEC, GBM, PDAC).  NOT YET RUN.

v2 changes the job's purpose.  v1 stored nested-increment sufficient statistics; the derivation
(theory/DERIVATION.md 2.3 and 4.5) shows those cannot give the discriminator theta, which needs the
W-ONLY out-of-fold ridge prediction of the own-mRNA protein residual.  Rather than compute any
discriminator on the cluster, the job now EXPORTS everything needed to compute every statistic
LOCALLY (local_theta.py) after the pre-specification has been frozen.  By default it computes NO
statistic that compares morphology with protein/mRNA except the PUBLISHED incremental R^2 (used
only to prove that the exported inputs reproduce the published per-gene results).  v1 files are
kept as *_v1.*.

Reuses, UNCHANGED, the published pipeline (residual_analysis.py = "RA"):
  * loaders           RA.load_rna_matrix / RA.load_protein_matrix / RA.load_wsi_embeddings
  * patient / gene    common = sorted(rna & protein & wsi patients);
    alignment         common_genes = sorted(rna & protein genes)   (RA.main)
  * morphology        PCA(n_components=min(RA.N_PCS, n-1, dim), svd_solver="full") fitted ONCE
                      on ALL common patients (RA.main) -- nowhere in the published pipeline is a
                      PCA fitted in-fold, so exporting the 20 PCs is exact; the pooled embeddings
                      are exported as well so the PCA can be re-derived and checked locally.
  * estimator         closed-form ridge alpha=1, columns standardised on the training fold,
                      KFold(5, shuffle, random_state=0) on the gene's non-missing-protein
                      patients (theta_kernel.oof_ridge == RA.cv_r2 body; asserted to 1e-10)
  * permutation null  (only if --B > 0) patient<->PC-row correspondence permuted, one draw
                      shared by all genes, RandomState(0) sequence as RA.

FILES WRITTEN (under --out; tag = cohort)
  inputs_<c>.npz          schema_version=2; ALL numeric arrays float64 EXCEPT pcs / wsi_pooled which
                          keep the dtype RA produced (float32 embeddings -> PCA float32) so the
                          published fit is reproduced bit-for-bit:
      patients (n,)  genes (G,)                        fit genes = >=30 protein-complete patients
      mrna_log2tpm1 (n,G)  protein (n,G, NaN=missing)  counts (n,G)  tpm (n,G)   [case-mean of
                          replicate files, exactly as RA averages TPM; NaN -> 0 as RA]
      pcs (n,K)   wsi_pooled (n,D)   pca_evr (K,)     the published 20 morphology PCs
      rna_pcs20 (n,Kr)                                 transcriptome-wide RNA PCs, exactly as
                          transcriptome_baseline.py builds them (z-scored full RNA matrix,
                          PCA full SVD, K=min(20,n-10)) -- robustness variant for H_err
      label_<axis> (n,) str   axes operator, plex, quarter, scanner, mpp, medium, cohort_stream
                          from batch_leak_check.case_batch_labels (the labels the sitepack /
                          plex analyses used); "" = no label.  label_status records failures.
      n_rna_files (n,)  icc_log2tpm1 (G,)  (NaN where no replicate design)
      pub_n, pub_r2_rna, pub_r2_both, pub_incremental_r2, pub_pval, pub_fdr (G,)
                          PUBLISHED per-gene results aligned to `genes` (MORPHO_SIG_CSV snapshot)
      is_synthetic=False, cohort, ridge_alpha, cv_folds, ra_seed=0, ra_n_perm=1000
  gene_table_<c>.csv      one row per tested gene: depth (raw counts, TPM, log2TPM sd, detection),
                          protein completeness, published values, and the export's recomputation
                          of the published increment (incr_p_export etc.).
  patient_table_<c>.csv   per patient: library size, STAR N_* counters, n files, genes >= 1 count.
  replicate_icc_<c>.csv   per-gene one-way ICC(1) of log2(TPM+1) across a case's RNA files
                          (cohorts with >= 8 replicate cases: CCRCC, UCEC).
  replicate_files_<c>.npz per-FILE log2(TPM+1) and raw counts of those cases (files, cases, genes,
                          logtpm, counts) -> the ICC is recomputable locally (local_theta --check-icc).
  repro_check_<c>.json    exported inputs (re-read from disk) vs the published per-gene results,
                          all fit genes: max |diff| for n, r2_rna, r2_mrna_wsi, incremental_r2.
                          Verdict EXACT needs <= 1e-10; the job exits 3 otherwise.
  export_sha256_<c>.txt   sha256 of every exported file ("<hash> *name") -- check after transfer
                          and use for the freeze sidecar.
  run_manifest_<c>.json   args, versions, sha256 of RA / this script / theta_kernel, shapes, timings.
OPTIONAL (default OFF; --B N and --prespec-sha256 <hash of the frozen pre-spec>):
  perms_<c>.npy, wonly_obs_<c>.npy (G,5) [rho,N,Nc,D,S_cross], wonly_perm_<c>.npy (G,B,4)
  [N,Nc,D,S_cross].  This computes the discriminator ON THE CLUSTER; it is a speed fallback only,
  and refuses to start unless the pre-spec hash is supplied (it is recorded in the manifest).

Run (one cohort; see lsf_measurement_error.sh):  python -u measurement_error_checks.py --cohort ucec
Smoke:   --max-genes 40            Timing of the optional path (synthetic noise):  --bench
Resume (optional path only): re-run the same command with --resume.
Blindness note: nothing here was run on project data.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")   # N workers x multi-threaded BLAS = oversubscription

import argparse
import hashlib
import json
import multiprocessing as mp
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
SCRIPTS_DIR = os.environ.get("MORPHO_SCRIPTS_DIR", "/public/home/fjhui/ZW/scripts")
for _p in (str(HERE), SCRIPTS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import theta_kernel as tk   # noqa: E402  (shared numerics; must sit next to this file)

CV_FOLDS = tk.CV_FOLDS
RIDGE_ALPHA = tk.RIDGE_ALPHA
MIN_FIT_PATIENTS = 30            # RA._process_gene: valid.sum() < 30 -> gene skipped
REPRO_TOL = 1e-10
SCHEMA_VERSION = 2
LABEL_AXES = ["operator", "plex", "quarter", "scanner", "mpp", "medium", "cohort_stream"]
EXPECTED = {                     # (n_patients, n_sig_genes): Table 1, same as peer scripts
    "ccrcc": (103, 2191), "luad": (105, 2310), "ucec": (100, 2566),
    "gbm": (99, 702), "pdac": (137, 1572),
}
WONLY_OBS_COLS = ["rho", "N", "Nc", "D", "S_cross"]
WONLY_PERM_COLS = ["N", "Nc", "D", "S_cross"]

RA = None  # the published module, imported lazily so --bench needs no cluster paths


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def _load_ra():
    global RA
    import residual_analysis as _ra   # reads MORPHO_* env at import; OUT_DIR.mkdir only
    RA = _ra
    assert RA.CV_FOLDS == CV_FOLDS and RA.RIDGE_ALPHA == RIDGE_ALPHA, \
        "published constants changed: update this script deliberately"
    return RA


def sha256(path, nbytes=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read() if nbytes is None else f.read(nbytes))
    return h.hexdigest()


# --------------------------------------------------------------------------------------
# data: published loaders + the extra raw-count pass
# --------------------------------------------------------------------------------------
def load_rna_depth():
    """Second pass over the SAME files RA.load_rna_matrix reads (same manifest filters, same
    replicate averaging), keeping raw `unstranded` counts and the per-file TPM / counts.
    Returns dict(cnt_case, tpm_case, file_logtpm, file_cnt, file_names, lib)."""
    manifest = pd.read_csv(RA.RNA_MANIFEST, sep="\t")
    if "sample_type" in manifest.columns:
        manifest = manifest[manifest["sample_type"] == "Primary Tumor"]
    manifest = manifest[~manifest["case_submitter_id"].isin(RA.NON_CCRCC_CASES)]
    want = {"gene_id", "gene_name", "unstranded", "tpm_unstranded"}
    per_cnt, per_tpm, per_logtpm, per_names, lib = {}, {}, {}, {}, {}
    n_missing = n_failed = 0
    first = True
    for _, r in manifest.iterrows():
        case = r["case_submitter_id"]
        fp = RA.RNA_DIR / r["file_name"]
        if not fp.exists():
            n_missing += 1
            continue
        try:
            df = pd.read_csv(fp, sep="\t", skiprows=1, usecols=lambda c: c in want)
        except Exception:
            n_failed += 1
            continue
        if first:
            log(f"[depth] STAR columns present in first file: {sorted(df.columns)} "
                "(RA uses tpm_unstranded; counts used here = unstranded)")
            first = False
        if not {"unstranded", "tpm_unstranded", "gene_name"}.issubset(df.columns):
            sys.exit(f"[depth] FATAL: {fp.name} lacks unstranded/tpm_unstranded/gene_name")
        summ = {}
        if "gene_id" in df.columns:
            sr = df[df["gene_id"].astype(str).str.startswith("N_")]
            summ = {str(k): float(x) for k, x in zip(sr["gene_id"], sr["unstranded"])}
        g = df.dropna(subset=["gene_name"])
        cnt = g.groupby("gene_name")["unstranded"].mean()
        tpm = g.groupby("gene_name")["tpm_unstranded"].mean()
        per_cnt.setdefault(case, []).append(cnt)
        per_tpm.setdefault(case, []).append(tpm)
        per_logtpm.setdefault(case, []).append(np.log2(tpm + 1.0))
        per_names.setdefault(case, []).append(str(r["file_name"]))
        lib.setdefault(case, []).append(dict(
            lib_size=float(g["unstranded"].sum()), n_genes_detected=int((cnt > 0).sum()),
            **{k: summ.get(k, np.nan) for k in
               ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")}))
    log(f"[depth] files read for {len(per_cnt)} cases; manifest files missing={n_missing}, "
        f"unreadable={n_failed}")

    def case_mat(d):                       # mirrors RA.load_rna_matrix
        rows = [pd.concat(sl, axis=1).mean(axis=1).rename(c) for c, sl in d.items()]
        return pd.concat(rows, axis=1).T
    cnt_case = case_mat(per_cnt)
    tpm_case = case_mat(per_tpm)
    libdf = pd.DataFrame({c: pd.DataFrame(l).mean() for c, l in lib.items()}).T
    libdf["n_rna_files"] = pd.Series({c: len(l) for c, l in lib.items()})
    return dict(cnt_case=cnt_case, tpm_case=tpm_case, file_logtpm=per_logtpm, file_cnt=per_cnt,
                file_names=per_names, lib=libdf)


def depth_table(genes, common, prot_df, rna_log, depth):
    """Depth + completeness for every common gene, on each gene's fit patients."""
    cnt = depth["cnt_case"].reindex(index=common, columns=genes).to_numpy(float)
    tpm = depth["tpm_case"].reindex(index=common, columns=genes).to_numpy(float)
    cnt = np.nan_to_num(cnt, nan=0.0)      # RA: fillna(0) after the case average
    tpm = np.nan_to_num(tpm, nan=0.0)
    P = prot_df.loc[common, genes].to_numpy(float)
    L = rna_log.loc[common, genes].to_numpy(float)
    V = ~np.isnan(P)
    nv = V.sum(0)

    def m(A):  # masked to the gene's fit patients
        return np.where(V, A, np.nan)
    with np.errstate(all="ignore"):
        out = pd.DataFrame({
            "gene": genes,
            "n_common": len(common),
            "n_fit": nv,
            "pct_missing_protein": 100.0 * (1 - nv / len(common)),
            "protein_mean_log_ratio": np.nanmean(P, 0),
            "protein_sd": np.nanstd(P, 0, ddof=1),
            "cnt_mean": np.nanmean(m(cnt), 0),
            "cnt_median": np.nanmedian(m(cnt), 0),
            "tpm_median": np.nanmedian(m(tpm), 0),
            "mrna_mean_log2tpm1": np.nanmean(m(L), 0),
            "mrna_sd_log2tpm1": np.nanstd(m(L), 0, ddof=1),
            "frac_tpm_lt1": np.nansum(np.where(V, tpm < 1, 0), 0) / np.maximum(nv, 1),
            "frac_cnt_zero": np.nansum(np.where(V, cnt == 0, 0), 0) / np.maximum(nv, 1),
        })
    out.loc[out["n_fit"] == 0, ["frac_tpm_lt1", "frac_cnt_zero"]] = np.nan
    return out, cnt, tpm


def export_labels(cancer, common, labels_file):
    """Per-patient batch labels, one array per axis ("" = none).  Source = the very function the
    sitepack / plex analyses call (batch_leak_check.case_batch_labels); failures are RECORDED and
    produce empty labels (a stratified local run then refuses), they do not abort the export."""
    out = {a: np.array([""] * len(common), dtype=object) for a in LABEL_AXES}
    status = {}
    if labels_file:
        df = pd.read_csv(labels_file, sep="\t", dtype=str).fillna("")
        key = "case_id" if "case_id" in df.columns else df.columns[0]
        mp_ = df.set_index(key)
        for a in LABEL_AXES:
            if a in mp_.columns:
                out[a] = np.array([mp_[a].get(c, "") for c in common], dtype=object)
                status[a] = f"from --labels-file ({int((out[a] != '').sum())}/{len(common)} labelled)"
        return out, status
    try:
        import batch_leak_check as blc
        paths = blc.cohort_paths(cancer)
        axes = blc.case_batch_labels(cancer, paths, common, strict_plex=False)
        for a in LABEL_AXES:
            if a in axes:
                lab = axes[a][0]
                out[a] = np.array([str(lab.get(c, "")) for c in common], dtype=object)
                status[a] = f"{int((out[a] != '').sum())}/{len(common)} labelled"
            else:
                status[a] = "axis not available for this cohort"
    except BaseException as e:   # SystemExit from a sha mismatch included
        status["_error"] = f"{type(e).__name__}: {e}"
        log(f"[labels] FAILED (continuing, labels left empty): {status['_error']}")
    return out, status


def load_published(sig_csv, genes_fit):
    """Published per-gene results aligned to the fit genes (NaN if absent)."""
    names = ["n", "r2_rna", "r2_mrna_wsi", "incremental_r2", "pval", "fdr"]
    out = {k: np.full(len(genes_fit), np.nan) for k in names}
    if not (sig_csv and os.path.exists(sig_csv)):
        return out, False
    pin = pd.read_csv(sig_csv, float_precision="round_trip")   # exact float parsing
    if pin["gene"].duplicated().any():
        sys.exit(f"FATAL: duplicated genes in {sig_csv}")
    pin = pin.set_index("gene").reindex(genes_fit)
    for k in names:
        if k in pin.columns:
            out[k] = pin[k].to_numpy(float)
    log(f"[pub] {int(np.isfinite(out['incremental_r2']).sum())}/{len(genes_fit)} fit genes found in {sig_csv}")
    return out, True


# --------------------------------------------------------------------------------------
# repro of the published increments from the exported file (re-read from disk)
# --------------------------------------------------------------------------------------
def repro_from_disk(npz_path, have_pub):
    z = np.load(npz_path, allow_pickle=False)
    mr, pr, pcs = z["mrna_log2tpm1"], z["protein"], z["pcs"]
    G = pr.shape[1]
    got = np.full((G, 4), np.nan)
    for g in range(G):
        v = ~np.isnan(pr[:, g])
        got[g] = tk.published_increment(mr[:, g], pr[:, g], pcs, v)
    rep = dict(n_fit_genes=G)
    if have_pub:
        pub = np.column_stack([z["pub_n"], z["pub_r2_rna"], z["pub_r2_both"], z["pub_incremental_r2"]])
        ok = np.isfinite(pub[:, 3])
        rep["n_compared"] = int(ok.sum())
        rep["n_fit_genes_missing_from_published"] = int((~ok).sum())
        for j, nm in enumerate(["n", "r2_rna", "r2_both", "incremental_r2"]):
            rep[f"max_abs_diff_{nm}"] = float(np.max(np.abs(got[ok, j] - pub[ok, j]))) if ok.any() else None
        worst = max(rep[f"max_abs_diff_{nm}"] for nm in ["r2_rna", "r2_both", "incremental_r2"])
        rep["verdict"] = ("EXACT (<=1e-10)" if worst <= REPRO_TOL and rep["max_abs_diff_n"] == 0
                          else ("CLOSE (<=1e-4) but not exact: do NOT freeze until explained"
                                if worst <= 1e-4 else "DRIFT vs published: do NOT use"))
    else:
        rep["verdict"] = "NO PUBLISHED CSV (MORPHO_SIG_CSV unset): reproduction unproven"
    return got, rep


# --------------------------------------------------------------------------------------
# OPTIONAL on-cluster W-only permutation path
# --------------------------------------------------------------------------------------
_W = {}


def _init_worker(mrna, prot, pcs, PP, perm_mask, n_perm):
    _W.update(mrna=mrna, prot=prot, pcs=pcs, PP=PP, perm_mask=perm_mask, B=n_perm)


def _gene_job(gi):
    prot_col = _W["prot"][:, gi]
    v = ~np.isnan(prot_col)
    nv = int(v.sum())
    P = prot_col[v]
    M = _W["mrna"][v, gi]
    fid = tk.fold_ids_for(nv)
    ones = np.ones((1, nv))
    o = tk.gene_stats_batch(_W["pcs"][v][None], M, P, fid, ones)
    obs = np.full(len(WONLY_OBS_COLS), np.nan)
    if o is not None:
        obs = np.array([o[k][0] for k in WONLY_OBS_COLS])
    perm = None
    if _W["B"] > 0 and _W["perm_mask"][gi] and o is not None:
        pb = tk.gene_stats_batch(_W["PP"][:, v], M, P, fid, ones)
        perm = np.column_stack([pb[k] for k in WONLY_PERM_COLS])
    return gi, obs, perm


def _open(path, shape, dtype, resume, fill=np.nan):
    path = Path(path)
    if resume and path.exists():
        return np.lib.format.open_memmap(path, mode="r+")
    mm = np.lib.format.open_memmap(path, mode="w+", dtype=dtype, shape=shape)
    step = max(1, 2_000_000 // max(1, int(np.prod(shape[1:]))))
    for i in range(0, shape[0], step):
        mm[i:i + step] = fill
    mm.flush()
    return mm


def build_perms(n, B, scheme, seed):
    return tk.build_perms(n, B, scheme, seed).astype(np.int16)


def optional_wonly(args, out, tag, genes_fit, mrna_f, prot_f, pcs, sig_set, t0):
    n, G = prot_f.shape
    B = args.B
    if args.perm_genes == "all":
        perm_mask = np.ones(G, bool)
    else:
        if not sig_set:
            sys.exit("--perm-genes sig needs MORPHO_SIG_CSV")
        perm_mask = np.array([g in sig_set for g in genes_fit])
        if args.perm_bg > 0:
            bg = np.flatnonzero(~perm_mask)
            pick = np.random.RandomState(args.bg_seed).choice(bg, size=min(args.perm_bg, len(bg)),
                                                              replace=False)
            perm_mask[pick] = True
    log(f"[optional] W-only permutation null for {int(perm_mask.sum())}/{G} genes, B={B}, "
        f"scheme={args.perm_scheme}  (pre-spec sha256 {args.prespec_sha256[:12]}...)")
    perms = tk.build_perms(n, B, args.perm_scheme, args.seed)
    np.save(out / f"perms_{tag}.npy", perms.astype(np.int16))
    pcs = pcs.astype(np.float64)                       # exact upcast (RA's hstack with float64 mRNA does the same)
    PP = pcs[perms]                                    # (B, n, K)
    cfg = dict(cohort=tag, B=B, scheme=args.perm_scheme, seed=args.seed, n=n, G=G,
               perm_genes=args.perm_genes, perm_bg=args.perm_bg, bg_seed=args.bg_seed,
               genes_sha=hashlib.sha256("|".join(genes_fit).encode()).hexdigest())
    cfg_fp = out / f"config_{tag}.json"
    if args.resume and cfg_fp.exists() and json.loads(cfg_fp.read_text()) != cfg:
        sys.exit("--resume: stored config differs from this run")
    cfg_fp.write_text(json.dumps(cfg, indent=1))
    obs_mm = _open(out / f"wonly_obs_{tag}.npy", (G, len(WONLY_OBS_COLS)), np.float64, args.resume)
    perm_mm = _open(out / f"wonly_perm_{tag}.npy", (G, B, len(WONLY_PERM_COLS)), np.float32, args.resume)
    done_fp = out / f"done_{tag}.npy"
    done = np.load(done_fp) if (args.resume and done_fp.exists()) else np.zeros(G, bool)
    todo = [i for i in range(G) if not done[i]]
    t1 = time.time()
    nd = 0
    with mp.get_context().Pool(args.workers, initializer=_init_worker,
                               initargs=(mrna_f, prot_f, pcs, PP, perm_mask, B)) as pool:
        for gi, obs, perm in pool.imap_unordered(_gene_job, todo, chunksize=4):
            obs_mm[gi] = obs
            if perm is not None:
                perm_mm[gi] = perm.astype(np.float32)
            done[gi] = True
            nd += 1
            if nd % 200 == 0 or nd == len(todo):
                obs_mm.flush()
                perm_mm.flush()
                np.save(str(done_fp) + ".tmp.npy", done)
                os.replace(str(done_fp) + ".tmp.npy", done_fp)
                el = time.time() - t1
                log(f"  {nd}/{len(todo)} genes  elapsed {el/60:.1f} min  ETA {el/nd*(len(todo)-nd)/60:.1f} min")
    obs_mm.flush()
    perm_mm.flush()
    np.save(done_fp, done)


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def main(args):
    t0 = time.time()
    if args.B > 0 and not args.prespec_sha256:
        sys.exit("REFUSED: --B > 0 computes the discriminator on the cluster; pass "
                 "--prespec-sha256 <sha256 of the frozen pre-specification> (recorded in the manifest). "
                 "The default export job computes no morphology-vs-protein statistic except the published one.")
    _load_ra()
    tag = args.cohort
    out = Path(args.out) if args.out else Path(RA.OUT_DIR)
    out.mkdir(parents=True, exist_ok=True)
    sig_csv = os.environ.get("MORPHO_SIG_CSV")

    # ---- published loaders, unchanged ----
    rna = RA.load_rna_matrix()
    protein = RA.load_protein_matrix()
    wsi = RA.load_wsi_embeddings()
    common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
    log(f"common patients across all 3: {len(common)}")
    if tag in EXPECTED and abs(len(common) - EXPECTED[tag][0]) > 2:
        sys.exit(f"[sanity] STOP: {len(common)} common patients, expected ~{EXPECTED[tag][0]} "
                 f"for {tag}: loaders point at the wrong cohort")
    common_genes = sorted(set(rna.columns) & set(protein.columns))
    rna_all = rna.loc[common]                       # ALL RNA genes (transcriptome PCs, as transcriptome_baseline.py)
    rna, protein = rna.loc[common, common_genes], protein.loc[common, common_genes]
    wsi_raw = wsi.loc[common].values
    pca = PCA(n_components=min(RA.N_PCS, len(common) - 1, wsi_raw.shape[1]), svd_solver="full")
    pcs = pca.fit_transform(wsi_raw)
    n = len(common)
    log(f"n={n} patients, {len(common_genes)} common genes, PCA {pcs.shape[1]} comps "
        f"({pca.explained_variance_ratio_.sum():.1%} var), dtypes: wsi {wsi_raw.dtype}, pcs {pcs.dtype}")

    # ---- transcriptome-wide RNA PCs, exactly as transcriptome_baseline.py ----
    rv = rna_all.values.astype(float)
    sdv = rv.std(0)
    rz = (rv - rv.mean(0)) / np.where(sdv == 0, 1.0, sdv)
    k_max = min(20, n - 10, rz.shape[1])
    rna_pcs20 = PCA(n_components=k_max, svd_solver="full").fit_transform(rz)
    log(f"RNA PCs: {rna_pcs20.shape} from {rv.shape[1]} genes")

    # ---- raw-count pass; must reproduce RA's RNA matrix exactly ----
    depth = load_rna_depth()
    chk = np.log2(depth["tpm_case"].fillna(0) + 1).reindex(index=common, columns=common_genes)
    dmax = float(np.nanmax(np.abs(chk.to_numpy() - rna.to_numpy())))
    log(f"[depth] max |log2(TPM+1) re-read - RA matrix| = {dmax:.2e}")
    if not dmax < 1e-9:
        sys.exit("[depth] FATAL: raw-count pass selected different files/cases than RA.load_rna_matrix")
    depth_df, cnt, tpm = depth_table(common_genes, common, protein, rna, depth)
    libdf = depth["lib"].reindex(common)
    libdf.index.name = "patient"
    libdf.to_csv(out / f"patient_table_{tag}.csv")

    # ---- replicate reliability (cohorts with >1 RNA file per case) ----
    icc_all = np.full(len(common_genes), np.nan)
    rep_cases = [c for c in common if len(depth["file_logtpm"].get(c, [])) >= 2]
    if len(rep_cases) >= 8:
        groups = [pd.concat(depth["file_logtpm"][c], axis=1).reindex(common_genes).fillna(0).to_numpy().T
                  for c in rep_cases]
        icc_all, sdw, a, kbar = tk.icc_oneway(groups)
        pd.DataFrame({"gene": common_genes, "icc1_log2tpm1": icc_all, "within_sd_log2tpm1": sdw,
                      "n_cases_with_replicates": a, "mean_files_per_case": kbar}
                     ).to_csv(out / f"replicate_icc_{tag}.csv", index=False)
        log(f"[icc] {a} cases with replicates, mean {kbar:.2f} files/case")
    else:
        log(f"[icc] only {len(rep_cases)} cases with >=2 RNA files in the analysis set: no ICC file")

    # ---- fit-gene set ----
    prot = protein.to_numpy(float)
    mrna = rna.to_numpy(float)
    nvalid = (~np.isnan(prot)).sum(0)
    fit_cols = np.flatnonzero(nvalid >= MIN_FIT_PATIENTS)
    genes_all = np.array(common_genes)
    if args.genes_file:
        keep = set(Path(args.genes_file).read_text().split())
        fit_cols = np.array([c for c in fit_cols if genes_all[c] in keep], dtype=int)
    if args.max_genes:
        fit_cols = fit_cols[:args.max_genes]
    genes_fit = genes_all[fit_cols].tolist()
    G = len(genes_fit)
    log(f"fit genes: {G} (>= {MIN_FIT_PATIENTS} protein-complete patients)")
    prot_f, mrna_f = prot[:, fit_cols], mrna[:, fit_cols]
    cnt_f, tpm_f = cnt[:, fit_cols], tpm[:, fit_cols]

    sig_set = set()
    if sig_csv and os.path.exists(sig_csv):
        md = pd.read_csv(sig_csv)
        if {"fdr", "incremental_r2"}.issubset(md.columns):
            sig_set = set(md.loc[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0), "gene"].astype(str))
            if tag in EXPECTED and len(sig_set) != EXPECTED[tag][1] and not args.skip_expected_check:
                sys.exit(f"[sanity] STOP: {sig_csv} has {len(sig_set)} significant genes, "
                         f"expected {EXPECTED[tag][1]} for {tag}")
    else:
        log("[sanity] WARNING: MORPHO_SIG_CSV not set/found: published results cannot be exported")
    pub, have_pub = load_published(sig_csv, genes_fit)

    # ---- labels ----
    labels, label_status = export_labels(tag, common, args.labels_file)
    log(f"[labels] {label_status}")

    # ---- self-test: oof_ridge == RA.cv_r2 on the first genes ----
    worst = 0.0
    for gi in range(min(5, G)):
        v = ~np.isnan(prot_f[:, gi])
        y = prot_f[v, gi]
        xa = mrna_f[v, gi].reshape(-1, 1)
        for X in (xa, np.hstack([xa, pcs[v]])):
            a_ = tk.r2_of(y, tk.oof_ridge(X, y, tk.folds_for(int(v.sum()))))
            worst = max(worst, abs(a_ - RA.cv_r2(X, y, CV_FOLDS)))
    log(f"[selftest] theta_kernel.oof_ridge vs RA.cv_r2: max |diff| = {worst:.2e}")
    if worst > 1e-10:
        sys.exit("[selftest] FATAL: oof_ridge diverges from the published cv_r2")

    # ---- the inputs file ----
    z = dict(schema_version=np.int64(SCHEMA_VERSION), is_synthetic=np.bool_(False), cohort=np.array(tag),
             patients=np.array(common), genes=np.array(genes_fit),
             mrna_log2tpm1=mrna_f.astype(np.float64), protein=prot_f.astype(np.float64),
             counts=cnt_f.astype(np.float64), tpm=tpm_f.astype(np.float64),
             pcs=pcs, wsi_pooled=wsi_raw, pca_evr=pca.explained_variance_ratio_.astype(np.float64),
             rna_pcs20=rna_pcs20.astype(np.float64),
             n_rna_files=libdf["n_rna_files"].to_numpy(float),
             icc_log2tpm1=icc_all[fit_cols].astype(np.float64),
             ridge_alpha=np.float64(RIDGE_ALPHA), cv_folds=np.int64(CV_FOLDS),
             ra_seed=np.int64(0), ra_n_perm=np.int64(RA.N_PERM),
             label_status=np.array(json.dumps(label_status)))
    for a in LABEL_AXES:
        z[f"label_{a}"] = np.array([str(x) for x in labels[a]])
    z.update(pub_n=pub["n"], pub_r2_rna=pub["r2_rna"], pub_r2_both=pub["r2_mrna_wsi"],
             pub_incremental_r2=pub["incremental_r2"], pub_pval=pub["pval"], pub_fdr=pub["fdr"])
    npz_path = out / f"inputs_{tag}.npz"
    np.savez_compressed(npz_path, **z)
    log(f"[export] {npz_path.name}: {npz_path.stat().st_size/1e6:.1f} MB")

    # per-file replicate values (ICC recomputable locally)
    if len(rep_cases) >= 8:
        fn, fc, lt, cn = [], [], [], []
        fit_genes_arr = np.array(genes_fit)
        for c in rep_cases:
            for nm, ls, cs in zip(depth["file_names"][c], depth["file_logtpm"][c], depth["file_cnt"][c]):
                fn.append(nm)
                fc.append(c)
                lt.append(ls.reindex(fit_genes_arr).fillna(0).to_numpy(float))
                cn.append(cs.reindex(fit_genes_arr).fillna(0).to_numpy(float))
        np.savez_compressed(out / f"replicate_files_{tag}.npz", files=np.array(fn), cases=np.array(fc),
                            genes=fit_genes_arr, logtpm=np.stack(lt), counts=np.stack(cn))

    # ---- reproduction of the published increments FROM THE FILE ON DISK ----
    log("[repro] recomputing the published increments from the exported file ...")
    got, rep = repro_from_disk(npz_path, have_pub)
    log(f"[repro] {rep}")
    (out / f"repro_check_{tag}.json").write_text(json.dumps(rep, indent=1))

    # ---- gene table (depth + completeness + published + export recomputation) ----
    gt = depth_df[depth_df["gene"].isin(genes_fit)].copy()
    gt = gt.set_index("gene").loc[genes_fit].reset_index()
    gt["icc1_log2tpm1"] = icc_all[fit_cols]
    gt["lambda_poisson_ub"] = tk.lambda_poisson(np.where(~np.isnan(prot_f), cnt_f, np.nan))
    gt["pub_incremental_r2"] = pub["incremental_r2"]
    gt["pub_r2_rna"] = pub["r2_rna"]
    gt["pub_pval"] = pub["pval"]
    gt["pub_fdr"] = pub["fdr"]
    gt["pub_significant"] = (pub["fdr"] < 0.05) & (pub["incremental_r2"] > 0)
    gt["r2_rna_export"] = got[:, 1]
    gt["r2_both_export"] = got[:, 2]
    gt["incr_p_export"] = got[:, 3]
    gt.to_csv(out / f"gene_table_{tag}.csv", index=False)

    # ---- optional on-cluster discriminator ----
    if args.B > 0:
        optional_wonly(args, out, tag, genes_fit, mrna_f, prot_f, pcs, sig_set, t0)

    # ---- manifest + checksums ----
    man = dict(cohort=tag, n_patients=n, n_common_genes=len(common_genes), n_fit_genes=G,
               n_pcs=int(pcs.shape[1]), pcs_dtype=str(pcs.dtype), wsi_dtype=str(wsi_raw.dtype),
               rna_pcs=int(rna_pcs20.shape[1]), args=vars(args), repro=rep, label_status=label_status,
               wall_seconds=round(time.time() - t0, 1), python=sys.version.split()[0],
               numpy=np.__version__, pandas=pd.__version__, host=platform.node(),
               finished=time.strftime("%Y-%m-%d %H:%M:%S"), sha256_script=sha256(__file__),
               sha256_theta_kernel=sha256(tk.__file__), sha256_residual_analysis=sha256(RA.__file__),
               sig_csv=sig_csv, sig_csv_sha256=(sha256(sig_csv) if sig_csv and os.path.exists(sig_csv) else None),
               schema_version=SCHEMA_VERSION)
    (out / f"run_manifest_{tag}.json").write_text(json.dumps(man, indent=1, default=str))
    keep_names = {f"perms_{tag}.npy", f"wonly_obs_{tag}.npy", f"wonly_perm_{tag}.npy"}
    names = [p for p in sorted(out.iterdir())
             if p.is_file() and (p.name in keep_names or p.name.endswith((f"_{tag}.npz", f"_{tag}.csv"))
                                 or (p.name.endswith(f"_{tag}.json") and not p.name.startswith("config_")))]
    with open(out / f"export_sha256_{tag}.txt", "w") as f:
        for p in names:
            f.write(f"{sha256(p)} *{p.name}\n")
    log(f"[done] {tag}: {len(names)} files checksummed; verdict: {rep['verdict']}")
    if not rep["verdict"].startswith("EXACT"):
        sys.exit(3)


def bench(args):
    """Time the OPTIONAL on-cluster path on SYNTHETIC noise (no project data)."""
    rs = np.random.RandomState(1)
    n, K = args.bench_n, 20
    W = rs.randn(n, K)
    M, P = rs.randn(n), rs.randn(n)
    fid = tk.fold_ids_for(n)
    for B in (50, 200):
        PP = W[tk.build_perms(n, B)]
        t = time.perf_counter()
        reps = 20
        for _ in range(reps):
            tk.gene_stats_batch(PP, M, P, fid, np.ones((1, n)))
        per = (time.perf_counter() - t) / reps
        print(f"n={n}, B={B}: {per*1e3:.2f} ms per gene (all B draws, W-only statistics)")
        for G in (2000, 10000):
            cpu_h = G * per * (1000 / B) / 3600   # scaled to B=1000
            print(f"  G={G:>6} B=1000: {cpu_h:6.2f} CPU-h  -> {cpu_h/args.workers:6.2f} h wall on {args.workers} workers")


def parse():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cohort", choices=sorted(EXPECTED), help="tag used in output file names")
    ap.add_argument("--out", default=None, help="output dir (default: RA.OUT_DIR = MORPHO_OUT)")
    ap.add_argument("--labels-file", default=None,
                    help="TSV (case_id + axis columns) replacing batch_leak_check labels (tests/fallback)")
    ap.add_argument("--skip-expected-check", action="store_true",
                    help="TEST ONLY: skip the Table-1 significant-gene-count assertion")
    ap.add_argument("--genes-file", default=None, help="restrict fit genes (whitespace-separated)")
    ap.add_argument("--max-genes", type=int, default=0, help="smoke test: first N fit genes")
    # optional on-cluster discriminator (default OFF)
    ap.add_argument("--B", type=int, default=0,
                    help="OPTIONAL on-cluster W-only permutation null (0 = off, the default)")
    ap.add_argument("--prespec-sha256", default="", help="required with --B > 0")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--perm-scheme", choices=["ra", "indexed"], default="ra")
    ap.add_argument("--perm-genes", choices=["all", "sig"], default="all")
    ap.add_argument("--perm-bg", type=int, default=0)
    ap.add_argument("--bg-seed", type=int, default=20261008)
    ap.add_argument("--workers", type=int,
                    default=int(os.environ.get("MORPHO_N_WORKERS", max(1, (os.cpu_count() or 4) - 1))))
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--bench-n", type=int, default=105)
    a = ap.parse_args()
    if not a.bench and not a.cohort:
        ap.error("--cohort is required")
    return a


if __name__ == "__main__":
    args = parse()
    if args.bench:
        bench(args)
    else:
        main(args)
