#!/usr/bin/env python3
"""
measurement_error_checks.py -- raw material for a measurement-error check of the
"beyond the gene's own mRNA" claim, for the five discovery cohorts (CCRCC, LUAD,
UCEC, GBM, PDAC).  Drafted 2026-10-08.  NOT YET RUN.  The exact test statistic is
deliberately NOT fixed here: this script writes per-gene, per-cohort raw material
(depth metrics, out-of-fold predictions, and per-permutation sufficient statistics)
so that the pre-specification can choose the statistic afterwards, locally.

Reuses, UNCHANGED, the published pipeline (residual_analysis.py = "RA"):
  * loaders           RA.load_rna_matrix / RA.load_protein_matrix / RA.load_wsi_embeddings
  * patient / gene    common = sorted(rna & protein & wsi patients);
    alignment         common_genes = sorted(rna & protein genes)   (RA.main lines 254-265)
  * morphology        PCA(n_components=min(RA.N_PCS, n-1, dim), svd_solver="full")
                      fitted on ALL common patients (RA.main lines 270-272)
  * estimator         closed-form ridge, alpha=RA.RIDGE_ALPHA, columns standardised on
                      the training fold, y centred, KFold(5, shuffle, random_state=0)
                      on the gene's non-missing-protein patients (RA.cv_r2)
  * permutation null  patient<->WSI-PC-row correspondence permuted, ONE permutation
                      shared by all genes per draw (RA.main lines 288-289); the default
                      permutation scheme reproduces RA's RandomState(0) sequence, so the
                      first B draws are the very draws of the published null.
The only new numerical code is oof_ridge(), a line-for-line copy of RA.cv_r2's body that
returns the out-of-fold predictions instead of R^2; _selftest() asserts that
1 - SSE/SST from oof_ridge equals RA.cv_r2 to 1e-10 on real genes before anything is
written.

WHAT IS WRITTEN (all under --out, tag = cohort):
  gene_table_<c>.csv      one row per tested gene (RA's common_genes; those with >=30
                          protein-complete patients are the "fit" genes):
        (a) depth / completeness (computed on the gene's fit patients = non-missing
            protein; raw counts = GDC STAR `unstranded`, replicate files of a case
            averaged exactly as RA averages TPM):
            n_common, n_fit, pct_missing_protein, protein_mean_log_ratio, protein_sd,
            cnt_mean, cnt_median, tpm_median, mrna_mean_log2tpm1, mrna_sd_log2tpm1,
            frac_tpm_lt1, frac_cnt_zero
        (b) observed OOF summaries (protein-given-mRNA = "P direction",
            mRNA-given-protein = "M direction"): n, sst_p/m, r2_rna_p, r2_both_p,
            incr_p (== RA incremental R^2), r2_rna_m, r2_both_m, incr_m, plus the raw
            sums listed under SUFFICIENT STATISTICS; descriptive permutation p-values
            pval_p / pval_m (upper tail, (1+#null>=obs)/(1+B)) -- convenience only.
  oof_pred_<c>.npy       float64 (G_fit, 4, n_common): OOF predictions on the gene's fit
                          patients, NaN elsewhere.  axis1 = [P_rna, P_both, M_rna, M_both]
                          P_rna  = OOF pred of protein from own mRNA           (baseline)
                          P_both = OOF pred of protein from own mRNA + 20 PCs  (RA model)
                          M_rna  = OOF pred of mRNA from own protein
                          M_both = OOF pred of mRNA from own protein + 20 PCs
                          Same folds in both directions; fold ids in oof_fold_<c>.npy.
  oof_fold_<c>.npy       int8 (G_fit, n_common): fold id (0..4) of each fit patient, -1 else.
  inputs_<c>.npz         patients, genes (fit genes), mrna_log2tpm1 (n x G_fit),
                          protein (n x G_fit), counts (n x G_fit, case-mean unstranded),
                          pcs (n x K), so every y / x vector can be rebuilt locally.
  perm_stats_<c>.npy     (G_fit, B, 7) per-permutation SUFFICIENT STATISTICS (dtype
                          --perm-dtype, default float32; computed in float64).  NaN rows =
                          gene not permutation-tested (--perm-genes).
  perms_<c>.npy          int16 (B, n_common): the permutations used (row b: rows of the PC
                          matrix are reordered as pcs[perm_b]), so the null is replayable.
  patient_table_<c>.csv  per-patient library size (sum of gene-row unstranded counts),
                          STAR N_unmapped/N_multimapping/N_noFeature/N_ambiguous, number of
                          files averaged, genes with >=1 count.
  replicate_icc_<c>.csv  per-gene one-way ICC of log2(TPM+1) across a case's RNA files
                          (only cohorts with replicate files: CCRCC, UCEC), NaN otherwise.
  repro_check_<c>.json   observed incr_p vs the pinned residual_results_tumoronly.csv
                          (max |diff|; p-value agreement when B=1000 and scheme=ra).
  run_manifest_<c>.json  args, versions, sha256 of RA / this script, shapes, timings.

SUFFICIENT STATISTICS.  For each gene, with r = y - p_rna (OOF residual of the baseline)
and d = p_both - p_rna (OOF increment due to morphology), per direction
(p: y=protein, m: y=mRNA), over the n_fit patients:
      s_dp   = sum d_p          s_dm   = sum d_m
      s_dpdp = sum d_p^2        s_dmdm = sum d_m^2
      s_dpdm = sum d_p d_m
      s_dprp = sum d_p r_p      s_dmrm = sum d_m r_m
Permutation-invariant (observed only, stored in the CSV): n, sst_p, sst_m (centred total
SS of y), sum_r_p, sum_r_m, ss_r_p, ss_r_m, sum_rp_rm.  Everything below follows:
      incremental R^2 (published)   incr = (2 s_dr - s_dd) / sst
      corr(d, r)                    = (s_dr - s_d*sum_r/n) / sqrt((s_dd - s_d^2/n)(ss_r - sum_r^2/n))
      corr(d_p, d_m)                = (s_dpdm - s_dp*s_dm/n) / sqrt(...)
Under permutation only p_both changes (the mRNA-only/protein-only baselines do not involve
the PCs), so d and these seven sums are the whole permutation-dependent state.

Run (one cohort; see lsf_measurement_error.sh):  python -u measurement_error_checks.py --cohort ucec --B 200
Smoke:   --max-genes 40 --B 3 --workers 2          Timing on the target node:  --bench
Resume:  re-run the same command with --resume (checkpoint every 200 genes).
Blindness note for the author of this script: nothing here was run on project data.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")   # N workers x multi-threaded BLAS = oversubscription

import argparse
import hashlib
import json
import multiprocessing as mp
import platform
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
SCRIPTS_DIR = os.environ.get("MORPHO_SCRIPTS_DIR", "/public/home/fjhui/ZW/scripts")
for _p in (str(HERE), SCRIPTS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Published constants (re-asserted against RA at run time in _load_ra()).
CV_FOLDS = 5
RIDGE_ALPHA = 1.0
MIN_FIT_PATIENTS = 30            # RA._process_gene: valid.sum() < 30 -> gene skipped
EXPECTED = {                     # (n_patients, n_sig_genes): Table 1, same as peer scripts
    "ccrcc": (103, 2191), "luad": (105, 2310), "ucec": (100, 2566),
    "gbm": (99, 702), "pdac": (137, 1572),
}
PERM_STATS = ["s_dp", "s_dm", "s_dpdp", "s_dmdm", "s_dpdm", "s_dprp", "s_dmrm"]
INV_COLS = ["n", "sst_p", "sst_m", "sum_r_p", "sum_r_m", "ss_r_p", "ss_r_m", "sum_rp_rm"]
OBS_COLS = INV_COLS + PERM_STATS          # 15 columns of per-gene observed scalars
OOF_AXES = ["P_rna", "P_both", "M_rna", "M_both"]

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


# --------------------------------------------------------------------------------------
# numerical kernels  (oof_ridge is RA.cv_r2's body, returning predictions)
# --------------------------------------------------------------------------------------
_FOLDS = {}


def folds_for(nv):
    """KFold(5, shuffle, random_state=0) split of the gene's nv fit patients. RA builds
    the splitter per call and splits X; the split depends only on nv, so caching is exact."""
    f = _FOLDS.get(nv)
    if f is None:
        kf = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=0)
        f = [(tr, te) for tr, te in kf.split(np.arange(nv))]
        _FOLDS[nv] = f
    return f


def oof_ridge(X, y, folds):
    preds = np.zeros_like(y, dtype=float)
    ridge_I = RIDGE_ALPHA * np.eye(X.shape[1])
    for tr, te in folds:
        X_tr, X_te, y_tr = X[tr], X[te], y[tr]
        x_mean = X_tr.mean(axis=0)
        x_std = X_tr.std(axis=0)
        x_std[x_std == 0] = 1.0
        Xtr_s = (X_tr - x_mean) / x_std
        Xte_s = (X_te - x_mean) / x_std
        y_mean = y_tr.mean()
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + ridge_I, Xtr_s.T @ (y_tr - y_mean))
        preds[te] = Xte_s @ w + y_mean
    return preds


def r2_of(y, p):
    ss_res = np.sum((y - p) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def inc_stats(d_p, d_m, r_p, r_m):
    """The seven permutation-dependent sufficient statistics (order = PERM_STATS)."""
    return np.array([d_p.sum(), d_m.sum(), d_p @ d_p, d_m @ d_m, d_p @ d_m,
                     d_p @ r_p, d_m @ r_m])


# --------------------------------------------------------------------------------------
# worker
# --------------------------------------------------------------------------------------
_W = {}


def _init_worker(mrna, prot, pcs, PP, perm_mask, n_perm):
    _W.update(mrna=mrna, prot=prot, pcs=pcs, PP=PP, perm_mask=perm_mask, B=n_perm)


def _gene_job(gi):
    """All work for fit gene index gi.  Returns (gi, oof(4,n), foldid(n), obs(15), perm(B,7)|None)."""
    prot_col = _W["prot"][:, gi]
    v = ~np.isnan(prot_col)
    nv = int(v.sum())
    yP = prot_col[v]                       # protein, fit patients
    xM = _W["mrna"][v, gi]                 # own mRNA log2(TPM+1), fit patients
    Wv = _W["pcs"][v]
    folds = folds_for(nv)

    xa = xM.reshape(-1, 1)                 # P direction baseline X (own mRNA)
    xb = yP.reshape(-1, 1)                 # M direction baseline X (own protein)
    pP0 = oof_ridge(xa, yP, folds)
    pP1 = oof_ridge(np.hstack([xa, Wv]), yP, folds)
    pM0 = oof_ridge(xb, xM, folds)
    pM1 = oof_ridge(np.hstack([xb, Wv]), xM, folds)

    r_p, r_m = yP - pP0, xM - pM0
    d_p, d_m = pP1 - pP0, pM1 - pM0
    sst_p = float(np.sum((yP - yP.mean()) ** 2))
    sst_m = float(np.sum((xM - xM.mean()) ** 2))
    inv = [nv, sst_p, sst_m, r_p.sum(), r_m.sum(), r_p @ r_p, r_m @ r_m, r_p @ r_m]
    obs = np.concatenate([np.array(inv, float), inc_stats(d_p, d_m, r_p, r_m)])

    n = len(prot_col)
    oof = np.full((4, n), np.nan)
    for k, arr in enumerate((pP0, pP1, pM0, pM1)):
        oof[k, v] = arr
    fid = np.full(n, -1, dtype=np.int8)
    for k, (_, te) in enumerate(folds):
        fid[np.flatnonzero(v)[te]] = k

    perm = None
    if _W["B"] > 0 and _W["perm_mask"][gi]:
        PP = _W["PP"]
        perm = np.empty((_W["B"], len(PERM_STATS)))
        for b in range(_W["B"]):
            Wp = PP[b][v]                  # == RA: wsi_pcs[perm][valid]
            dpb = oof_ridge(np.hstack([xa, Wp]), yP, folds) - pP0
            dmb = oof_ridge(np.hstack([xb, Wp]), xM, folds) - pM0
            perm[b] = inc_stats(dpb, dmb, r_p, r_m)
    return gi, oof, fid, obs, perm


# --------------------------------------------------------------------------------------
# data: published loaders + the extra raw-count pass
# --------------------------------------------------------------------------------------
def load_rna_depth(extra_cases=None):
    """Second pass over the SAME files RA.load_rna_matrix reads (same manifest filters,
    same replicate averaging), keeping raw `unstranded` counts and the per-file TPM.
    Returns dict(cnt_case, tpm_case, file_logtpm, lib, nfiles)."""
    manifest = pd.read_csv(RA.RNA_MANIFEST, sep="\t")
    if "sample_type" in manifest.columns:
        manifest = manifest[manifest["sample_type"] == "Primary Tumor"]
    manifest = manifest[~manifest["case_submitter_id"].isin(RA.NON_CCRCC_CASES)]
    want = {"gene_id", "gene_name", "unstranded", "tpm_unstranded"}
    per_cnt, per_tpm, per_logtpm, lib = {}, {}, {}, {}
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
    return dict(cnt_case=cnt_case, tpm_case=tpm_case, file_logtpm=per_logtpm, lib=libdf)


def icc_oneway(groups):
    """One-way random-effects ICC(1) for unbalanced groups; vectorised over genes.
    groups: list of (k_i, G) arrays (k_i >= 2)."""
    a = len(groups)
    ns = np.array([g.shape[0] for g in groups], float)
    N = ns.sum()
    means = np.stack([g.mean(0) for g in groups])
    grand = np.concatenate(groups).mean(0)
    msb = (ns[:, None] * (means - grand) ** 2).sum(0) / (a - 1)
    msw = sum(((g - g.mean(0)) ** 2).sum(0) for g in groups) / (N - a)
    k0 = (N - (ns ** 2).sum() / N) / (a - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        icc = (msb - msw) / (msb + (k0 - 1) * msw)
    return icc, np.sqrt(msw), a, float(N / a)


def depth_table(genes, common, prot_df, rna_log, depth):
    """(a): depth + completeness for every common gene, on each gene's fit patients."""
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
    return out, cnt


# --------------------------------------------------------------------------------------
# storage with checkpoint / resume
# --------------------------------------------------------------------------------------
def _open(path, shape, dtype, resume, fill=np.nan):
    path = Path(path)
    if resume and path.exists():
        return np.lib.format.open_memmap(path, mode="r+")
    mm = np.lib.format.open_memmap(path, mode="w+", dtype=dtype, shape=shape)
    if fill is not None:
        step = max(1, 2_000_000 // max(1, int(np.prod(shape[1:]))))
        for i in range(0, shape[0], step):
            mm[i:i + step] = fill
        mm.flush()
    return mm


def sha256(path, nbytes=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read() if nbytes is None else f.read(nbytes))
    return h.hexdigest()


# --------------------------------------------------------------------------------------
# main pieces
# --------------------------------------------------------------------------------------
def build_perms(n, B, scheme, seed):
    if B == 0:
        return np.zeros((0, n), dtype=np.int16)
    if scheme == "ra":                     # RA: RNG = RandomState(0); sequential draws
        rng = np.random.RandomState(seed)
        perms = [rng.permutation(n) for _ in range(B)]
    elif scheme == "indexed":              # transcriptome_baseline.py: RandomState(i)
        perms = [np.random.RandomState(seed + i).permutation(n) for i in range(B)]
    else:
        raise ValueError(scheme)
    return np.asarray(perms, dtype=np.int16)


def _selftest(fit_idx, prot, mrna, pcs, k=5):
    """oof_ridge must reproduce RA.cv_r2 exactly (protein~mRNA and protein~mRNA+PCs)."""
    worst = 0.0
    for gi in fit_idx[:k]:
        v = ~np.isnan(prot[:, gi])
        y = prot[v, gi]
        xa = mrna[v, gi].reshape(-1, 1)
        for X in (xa, np.hstack([xa, pcs[v]])):
            a = r2_of(y, oof_ridge(X, y, folds_for(int(v.sum()))))
            b = RA.cv_r2(X, y, CV_FOLDS)
            worst = max(worst, abs(a - b))
    log(f"[selftest] oof_ridge vs RA.cv_r2 on {min(k, len(fit_idx))} genes: max |diff| = {worst:.2e}")
    if worst > 1e-10:
        sys.exit("[selftest] FATAL: oof_ridge diverges from the published cv_r2")


def finalize(out, tag, genes_fit, depth_df, B, sig_csv, args, t0, meta_extra):
    G = len(genes_fit)
    obs = np.load(out / f"obs_scalars_{tag}.npy")
    df = pd.DataFrame(obs, columns=OBS_COLS)
    df.insert(0, "gene", genes_fit)
    for d in ("p", "m"):
        sst = df[f"sst_{d}"]
        df[f"r2_rna_{d}"] = 1 - df[f"ss_r_{d}"] / sst
        df[f"incr_{d}"] = (2 * df[f"s_d{d}r{d}"] - df[f"s_d{d}d{d}"]) / sst
        df[f"r2_both_{d}"] = df[f"r2_rna_{d}"] + df[f"incr_{d}"]
    if B > 0 and (out / f"perm_stats_{tag}.npy").exists():
        P = np.load(out / f"perm_stats_{tag}.npy", mmap_mode="r")
        for d in ("p", "m"):
            sst = df[f"sst_{d}"].to_numpy()
            dd, dr = PERM_STATS.index(f"s_d{d}d{d}"), PERM_STATS.index(f"s_d{d}r{d}")
            pv = np.full(G, np.nan)
            for i in range(0, G, 500):
                blk = np.asarray(P[i:i + 500], dtype=float)
                null = (2 * blk[:, :, dr] - blk[:, :, dd]) / sst[i:i + 500, None]
                ok = ~np.isnan(null[:, 0])
                ge = (null >= df[f"incr_{d}"].to_numpy()[i:i + 500, None]).sum(1)
                pv[i:i + 500] = np.where(ok, (ge + 1) / (B + 1), np.nan)
            df[f"pval_{d}"] = pv
    full = depth_df.merge(df, on="gene", how="left")
    full.to_csv(out / f"gene_table_{tag}.csv", index=False)

    chk = {}
    if sig_csv and os.path.exists(sig_csv):
        pin = pd.read_csv(sig_csv).set_index("gene")
        j = df.set_index("gene").join(pin[["r2_rna", "incremental_r2"] +
                                           (["pval"] if "pval" in pin.columns else [])],
                                       how="inner")
        chk["n_compared"] = int(len(j))
        chk["max_abs_diff_incr"] = float(np.nanmax(np.abs(j["incr_p"] - j["incremental_r2"])))
        chk["max_abs_diff_r2_rna"] = float(np.nanmax(np.abs(j["r2_rna_p"] - j["r2_rna"])))
        if "pval_p" in df.columns and "pval" in j.columns and B == 1000 and args.perm_scheme == "ra":
            jj = j.dropna(subset=["pval_p"])
            chk["n_pval_compared"] = int(len(jj))
            chk["max_abs_diff_pval"] = float(np.nanmax(np.abs(jj["pval_p"] - jj["pval"])))
        chk["verdict"] = ("MATCHES published pipeline" if chk["max_abs_diff_incr"] < 1e-4
                          else "DRIFT vs pinned CSV: do not use until explained")
        log(f"[repro] {chk}")
    (out / f"repro_check_{tag}.json").write_text(json.dumps(chk, indent=1))
    man = dict(meta_extra, args=vars(args), wall_seconds=round(time.time() - t0, 1),
               python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__,
               host=platform.node(), finished=time.strftime("%Y-%m-%d %H:%M:%S"),
               sha256_script=sha256(__file__), sha256_residual_analysis=sha256(RA.__file__))
    (out / f"run_manifest_{tag}.json").write_text(json.dumps(man, indent=1, default=str))
    log(f"[done] gene_table_{tag}.csv written ({len(full)} genes)")


def main(args):
    t0 = time.time()
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
    rna, protein = rna.loc[common, common_genes], protein.loc[common, common_genes]
    wsi_raw = wsi.loc[common].values
    pca = PCA(n_components=min(RA.N_PCS, len(common) - 1, wsi_raw.shape[1]), svd_solver="full")
    pcs = pca.fit_transform(wsi_raw)
    n = len(common)
    log(f"n={n} patients, {len(common_genes)} common genes, PCA {pcs.shape[1]} comps "
        f"({pca.explained_variance_ratio_.sum():.1%} var)")

    # ---- raw-count pass; must reproduce RA's RNA matrix exactly ----
    depth = load_rna_depth()
    chk = np.log2(depth["tpm_case"].fillna(0) + 1).reindex(index=common, columns=common_genes)
    dmax = float(np.nanmax(np.abs(chk.to_numpy() - rna.to_numpy())))
    log(f"[depth] max |log2(TPM+1) re-read - RA matrix| = {dmax:.2e}")
    if not dmax < 1e-9:
        sys.exit("[depth] FATAL: raw-count pass selected different files/cases than RA.load_rna_matrix")
    depth_df, cnt = depth_table(common_genes, common, protein, rna, depth)
    libdf = depth["lib"].reindex(common)
    libdf.index.name = "patient"
    libdf.to_csv(out / f"patient_table_{tag}.csv")

    # ---- replicate reliability (cohorts with >1 RNA file per case) ----
    groups = []
    for c in common:
        fl = depth["file_logtpm"].get(c, [])
        if len(fl) >= 2:
            groups.append(pd.concat(fl, axis=1).reindex(common_genes).fillna(0).to_numpy().T)
    if len(groups) >= 8:
        icc, sdw, a, kbar = icc_oneway(groups)
        pd.DataFrame({"gene": common_genes, "icc1_log2tpm1": icc, "within_sd_log2tpm1": sdw,
                      "n_cases_with_replicates": a, "mean_files_per_case": kbar}
                     ).to_csv(out / f"replicate_icc_{tag}.csv", index=False)
        log(f"[icc] {a} cases with replicates, mean {kbar:.2f} files/case")
    else:
        log(f"[icc] only {len(groups)} cases with >=2 RNA files in the analysis set: no ICC file")

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

    sig_set = set()
    if sig_csv and os.path.exists(sig_csv):
        md = pd.read_csv(sig_csv)
        if {"fdr", "incremental_r2"}.issubset(md.columns):
            sig_set = set(md.loc[(md["fdr"] < 0.05) & (md["incremental_r2"] > 0), "gene"].astype(str))
            if tag in EXPECTED and len(sig_set) != EXPECTED[tag][1]:
                sys.exit(f"[sanity] STOP: {sig_csv} has {len(sig_set)} significant genes, "
                         f"expected {EXPECTED[tag][1]} for {tag}")
    else:
        log("[sanity] WARNING: MORPHO_SIG_CSV not set/found: no pinned-reproduction check")

    # ---- which genes get the permutation null ----
    B = args.B
    if args.perm_genes == "all" or B == 0:
        perm_mask = np.ones(G, bool)
    elif args.perm_genes == "sig":
        if not sig_set:
            sys.exit("--perm-genes sig needs MORPHO_SIG_CSV")
        perm_mask = np.array([g in sig_set for g in genes_fit])
        if args.perm_bg > 0:
            bg = np.flatnonzero(~perm_mask)
            pick = np.random.RandomState(args.bg_seed).choice(bg, size=min(args.perm_bg, len(bg)),
                                                              replace=False)
            perm_mask[pick] = True
    else:
        raise ValueError(args.perm_genes)
    log(f"permutation null for {int(perm_mask.sum())}/{G} genes, B={B}, scheme={args.perm_scheme}")

    _selftest(range(G), prot_f, mrna_f, pcs)

    perms = build_perms(n, B, args.perm_scheme, args.seed)
    np.save(out / f"perms_{tag}.npy", perms)
    PP = np.stack([pcs[p.astype(int)] for p in perms]) if B else np.zeros((0, n, pcs.shape[1]))

    np.savez_compressed(out / f"inputs_{tag}.npz", patients=np.array(common), genes=np.array(genes_fit),
                        mrna_log2tpm1=mrna_f.astype(np.float32), protein=prot_f.astype(np.float32),
                        counts=cnt[:, fit_cols].astype(np.float32), pcs=pcs)

    # ---- stores + resume ----
    cfg = dict(cohort=tag, B=B, scheme=args.perm_scheme, seed=args.seed, n=n, G=G,
               perm_genes=args.perm_genes, perm_bg=args.perm_bg, bg_seed=args.bg_seed,
               genes_sha=hashlib.sha256("|".join(genes_fit).encode()).hexdigest(),
               perm_dtype=args.perm_dtype)
    cfg_fp = out / f"config_{tag}.json"
    if args.resume and cfg_fp.exists() and json.loads(cfg_fp.read_text()) != cfg:
        sys.exit("--resume: stored config differs from this run (cohort/B/scheme/seed/genes)")
    cfg_fp.write_text(json.dumps(cfg, indent=1))
    oof_mm = _open(out / f"oof_pred_{tag}.npy", (G, 4, n), np.float64, args.resume)
    fold_mm = _open(out / f"oof_fold_{tag}.npy", (G, n), np.int8, args.resume, fill=-1)
    obs_mm = _open(out / f"obs_scalars_{tag}.npy", (G, len(OBS_COLS)), np.float64, args.resume)
    perm_mm = (_open(out / f"perm_stats_{tag}.npy", (G, B, len(PERM_STATS)),
                     np.dtype(args.perm_dtype), args.resume) if B else None)
    done_fp = out / f"done_{tag}.npy"
    done = np.load(done_fp) if (args.resume and done_fp.exists()) else np.zeros(G, bool)
    todo = [i for i in range(G) if not done[i]]
    log(f"{len(todo)} genes to do ({int(done.sum())} already done)")

    meta_extra = dict(cohort=tag, n_patients=n, n_common_genes=len(common_genes), n_fit_genes=G,
                      n_perm_genes=int(perm_mask.sum()), B=B, n_pcs=int(pcs.shape[1]),
                      RNA_columns_used_by_RA="tpm_unstranded -> log2(x+1)",
                      raw_count_column="unstranded", obs_columns=OBS_COLS, perm_stats=PERM_STATS,
                      oof_axes=OOF_AXES, config=cfg)

    t1 = time.time()
    nd = 0
    ctx = mp.get_context()
    with ctx.Pool(args.workers, initializer=_init_worker,
                  initargs=(mrna_f, prot_f, pcs, PP, perm_mask, B)) as pool:
        for gi, oof, fid, obs, perm in pool.imap_unordered(_gene_job, todo, chunksize=1):
            oof_mm[gi] = oof
            fold_mm[gi] = fid
            obs_mm[gi] = obs
            if perm is not None and perm_mm is not None:
                perm_mm[gi] = perm.astype(perm_mm.dtype)
            done[gi] = True
            nd += 1
            if nd % 200 == 0 or nd == len(todo):
                for mm in (oof_mm, fold_mm, obs_mm, perm_mm):
                    if mm is not None:
                        mm.flush()
                np.save(str(done_fp) + ".tmp.npy", done)
                os.replace(str(done_fp) + ".tmp.npy", done_fp)
                el = time.time() - t1
                log(f"  {nd}/{len(todo)} genes  elapsed {el/60:.1f} min  "
                    f"ETA {el/nd*(len(todo)-nd)/60:.1f} min")
    for mm in (oof_mm, fold_mm, obs_mm, perm_mm):
        if mm is not None:
            mm.flush()
    np.save(done_fp, done)
    finalize(out, tag, genes_fit, depth_df, B, sig_csv, args, t0, meta_extra)


def bench(args):
    """Time the kernels on SYNTHETIC noise (no project data) to extrapolate on the target node."""
    rs = np.random.RandomState(1)
    n, K = args.bench_n, 20
    y, x, W = rs.randn(n), rs.randn(n, 1), rs.randn(n, K)
    folds = folds_for(n)
    X = np.hstack([x, W])
    for _ in range(20):
        oof_ridge(X, y, folds)
    reps = 400
    t = time.perf_counter()
    for _ in range(reps):
        oof_ridge(X, y, folds)
    per_fit = (time.perf_counter() - t) / reps
    per_draw = 2 * per_fit
    print(f"n={n}, 21 features, 5 folds: {per_fit*1e3:.3f} ms per OOF ridge fit")
    print(f"per gene per permutation draw (2 directions): {per_draw*1e3:.3f} ms")
    for G in (2000, 10000):
        for B in (200, 500, 1000):
            cpu_h = G * B * per_draw / 3600
            print(f"  G={G:>6} B={B:>5}: {cpu_h:6.2f} CPU-h  -> {cpu_h/args.workers:6.2f} h wall on {args.workers} workers")


def parse():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cohort", choices=sorted(EXPECTED), help="tag used in output file names")
    ap.add_argument("--out", default=None, help="output dir (default: RA.OUT_DIR = MORPHO_OUT)")
    ap.add_argument("--B", type=int, default=200, help="permutations (0 = observed only)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--perm-scheme", choices=["ra", "indexed"], default="ra",
                    help="ra: RandomState(seed) sequential draws (== RA.RNG); "
                         "indexed: RandomState(seed+i) per draw (== transcriptome_baseline.py)")
    ap.add_argument("--perm-genes", choices=["all", "sig"], default="all",
                    help="genes receiving the permutation null; 'sig' = pinned-significant "
                         "(+ --perm-bg random background genes)")
    ap.add_argument("--perm-bg", type=int, default=0)
    ap.add_argument("--bg-seed", type=int, default=20261008)
    ap.add_argument("--perm-dtype", choices=["float32", "float64"], default="float32")
    ap.add_argument("--genes-file", default=None, help="restrict fit genes (whitespace-separated)")
    ap.add_argument("--max-genes", type=int, default=0, help="smoke test: first N fit genes")
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
