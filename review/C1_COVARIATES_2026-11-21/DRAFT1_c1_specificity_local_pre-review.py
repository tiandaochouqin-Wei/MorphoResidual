#!/usr/bin/env python
"""Post hoc C1 specificity re-reads: the local stages.

Pre-specified in review/POSTHOC_C1_SPECIFICITY_PRESPEC_2026-11-21.md (sections cited as "PS n").
Stage runs are fail-closed: they refuse to start unless the sidecar
review/POSTHOC_C1_SPECIFICITY_PRESPEC_2026-11-21.md.sha256 exists and every file it lists (the
pre-specification, this script, the kernel, every input) still has the listed sha256.

Modes
  --selftest        synthetic data only: batched kernel against the reference kernel, BH against a
                    direct implementation, `assemble` against scipy's Mann-Whitney statistic.
  --validate-only   published quantities only (allowed before freezing, PS 9):
      V1  the code path used for every new discovery set, fed the published block (20 morphology
          components over the gene's own mRNA, the published permutations), reproduces the
          published per-gene increment, p-value, FDR and the published selection set.
      V2  `assemble` on the saved confirmatory matrices reproduces the published C1-UCEC and
          C1-LUAD readings (AUC, null mean, null SD, one-sided p, rho, n selected / background).
  --stage L1        discovery sets that need no file from the cluster (S_R20, S_OP, S_M|R20) and the
                    cells read on the saved confirmatory matrices (PS 6, forward cells for R20, OP).
  --stage L2        after Stage H: the composition-score sets (S_C3, S_M|C3) and every cell that
                    needs a Stage H matrix (PS 6).  Gate V3 first reproduces the published
                    composition controls from the transferred discovery scores.

The estimator is the published one (residual_analysis.py): per gene, on the patients with a
protein value, 5-fold KFold(shuffle, random_state 0), ridge alpha 1 on columns standardised on
the training fold, increment = R2(baseline + block) - R2(baseline), null = the block's rows
permuted with the published permutations (RandomState(0), 1,000 sequential draws), one-sided
p = (1 + #{null >= observed}) / 1001, Benjamini-Hochberg over all fit genes,
selected = FDR < 0.05 and increment > 0.  `bh_fdr` is copied unchanged from
residual_analysis.py and `assemble` from c1_run_test.py.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import collections
import hashlib
import json
import multiprocessing as mp
import sys
import time

import numpy as np
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))                # MorphoResidual_paper/
PRESPEC_REL = "review/POSTHOC_C1_SPECIFICITY_PRESPEC_2026-11-21.md"
SIDECAR = os.path.join(PAPER, "review", "POSTHOC_C1_SPECIFICITY_PRESPEC_2026-11-21.md.sha256")
KERNEL_DIR = os.path.join(PAPER, "review", "MEAS_ERROR_2026-10-08")
KERNEL_SHA = "6c95cd119d28e811e09967946e9f191c69f07f54f653991917e817ef99ec2734"
OUT_DIR = os.path.join(PAPER, "server_export", "results", "c1_specificity")
STAGE_H_DIR = os.path.join(OUT_DIR, "stage_h")                          # files transferred from the cluster

sys.path.insert(0, KERNEL_DIR)
import theta_kernel as tk  # noqa: E402  (reference ridge, folds, published permutations)

B_PERM = 1000
FDR_LEVEL = 0.05
MIN_SET = 50            # PS 5: a set or a background with fewer genes in the universe is untestable
MAX_OP_LEVELS = 12      # PS 2: operator dummies, as residual_analysis_sitepack.py
MIN_OP_GROUP = 3
TOL_REF = 1e-9          # batched kernel against the reference kernel
TOL_PUB = 1e-10         # reference kernel against the published per-gene values
TOL_AUC = 1e-12         # `assemble` against the published confirmatory readings

COHORTS = {
    "ucec": dict(
        inputs="review/MEAS_ERROR_2026-10-08/export/ucec/inputs_ucec.npz",
        inputs_sha="e28b86fadd41ec86363577fffad1db1f2f8a5162e732c7acb0e4133c0982044e",
        conf="server_export/scripts/pinned/c1_ucec_d4_stratified.npz",
        conf_sha="ef75f1c06d4fc33843a79e9a8d1bb4f4224acd201d30784badb97cf3be13384a",
        arms=[("unrestricted", "M_unrestricted"), ("within-operator", "M_within_operator"),
              ("residualised", "M_within_operator_residualized")],
        published="server_export/scripts/pinned/c1_ucec_d4_followup.json",
        compctl="server_export/results/results__composition_controls_ucec.csv",
        n_disc=100, n_fit=10512, n_sel=2566, n_univ=10146, n1=2536, n0=7610,
        n_op_dummies=12, n_op_covered=73),
    "luad": dict(
        inputs="review/MEAS_ERROR_2026-10-08/export/luad/inputs_luad.npz",
        inputs_sha="b64cdf3b31717a32e0d42c42f3f1dde74786d5598008f0d214297a71adaba095",
        conf="server_export/scripts/pinned/luad_c1/c1_luad_primary_primary.npz",
        conf_sha="2034e71143c2e065675d3159832b1a7be17aa38c7feeaacfc82862c411ea2aab",
        arms=[("unrestricted", "M")],
        published="server_export/scripts/pinned/luad_c1/c1_luad_primary_result_primary.json",
        compctl="server_export/results/results__composition_controls_luad.csv",
        n_disc=105, n_fit=10724, n_sel=2310, n_univ=10406, n1=2287, n0=8119,
        n_op_dummies=12, n_op_covered=75),
}
STAGE_H_ARMS = ("A_C3", "A_R20", "A_R20cap", "V_C3", "V_R20")           # PS 4
SCORE_COLS = ("stroma", "immune", "proliferation")                      # PS 2


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def P(rel):
    return os.path.join(PAPER, *rel.split("/"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------
# fail-closed gate
# ------------------------------------------------------------------------------------------
def check_sidecar(required):
    """Every line '<sha256> *<path relative to MorphoResidual_paper/>' is verified; the
    pre-specification and every file in `required` must be listed."""
    if not os.path.exists(SIDECAR):
        sys.exit(f"REFUSED: sidecar {SIDECAR} not found -- the pre-specification is not frozen.")
    listed = set()
    for line in open(SIDECAR, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        want, rel = line.split(None, 1)
        rel = rel.strip().lstrip("*")
        path = os.path.normcase(os.path.abspath(P(rel)))
        if not os.path.exists(path) or sha256(path) != want.lower():
            sys.exit(f"REFUSED: {rel} missing or changed since freezing")
        listed.add(path)
    for f in [P(PRESPEC_REL)] + list(required):
        if os.path.normcase(os.path.abspath(f)) not in listed:
            sys.exit(f"REFUSED: {f} is not listed in the sidecar (not frozen)")
    log("sidecar verified:", len(listed), "files")


def required_files(stage):
    req = [os.path.abspath(__file__), os.path.abspath(tk.__file__)]
    for c in COHORTS.values():
        req += [P(c["inputs"]), P(c["conf"]), P(c["published"])]
        if stage == "L2":
            req.append(P(c["compctl"]))
    return req


# ------------------------------------------------------------------------------------------
# copied unchanged
# ------------------------------------------------------------------------------------------
def bh_fdr(pvals):
    """Identical to residual_analysis.bh_fdr."""
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    fdr = np.empty(n)
    fdr[order] = np.clip(ranked, 0, 1)
    return fdr


def assemble(M, xs, sel_mask):
    """Identical to c1_run_test.assemble."""
    n1, n0 = int(sel_mask.sum()), int((~sel_mask).sum())
    ranks = rankdata(M, axis=0)
    auc = (ranks[sel_mask].sum(axis=0) - n1 * (n1 + 1) / 2.0) / (n1 * n0)
    rs = rankdata(xs)
    rs = (rs - rs.mean()) / rs.std()
    rc = (ranks - ranks.mean(axis=0)) / ranks.std(axis=0)
    rho = (rs[:, None] * rc).mean(axis=0)
    return auc, rho, n1, n0


def emp_p(obs, null):
    return float((1 + int(np.sum(null >= obs))) / (1 + len(null)))


# ------------------------------------------------------------------------------------------
# the estimator: reference (unbatched) and batched over permutations
# ------------------------------------------------------------------------------------------
def ref_increment(y, base, block_v, folds):
    """Reference: R2([base, block]) - R2(base) with theta_kernel.oof_ridge (= RA.cv_r2)."""
    r2_base = tk.r2_of(y, tk.oof_ridge(base, y, folds))
    r2_full = tk.r2_of(y, tk.oof_ridge(np.hstack([base, block_v]), y, folds))
    return r2_base, r2_full - r2_base


def batched_increments(y, base, block_full, v, idx, alpha=tk.RIDGE_ALPHA):
    """Increment of the block over the baseline for every row order in `idx`.

    y           (nv,)   protein of the gene's fit patients
    base        (nv,q)  baseline columns on the fit patients (never permuted)
    block_full  (n,k)   block rows of all cohort patients
    v           (n,)    boolean fit mask
    idx         (R,n)   row orders; row 0 must be the identity.  As in residual_analysis.py the
                        block is permuted over all n patients first and then restricted to the
                        fit patients: block_full[perm][v].
    Returns r2_base (reference kernel) and the R increments."""
    nv = int(v.sum())
    folds = tk.folds_for(nv)
    r2_base = tk.r2_of(y, tk.oof_ridge(base, y, folds))
    blk = block_full[idx][:, v, :]                               # (R, nv, k)
    R = blk.shape[0]
    q = base.shape[1]
    p = q + blk.shape[2]
    eye = alpha * np.eye(p)
    preds = np.empty((R, nv))
    for tr, te in folds:
        Xtr = np.empty((R, len(tr), p))
        Xte = np.empty((R, len(te), p))
        Xtr[:, :, :q] = base[tr]
        Xte[:, :, :q] = base[te]
        Xtr[:, :, q:] = blk[:, tr, :]
        Xte[:, :, q:] = blk[:, te, :]
        mu = Xtr.mean(axis=1, keepdims=True)
        sd = Xtr.std(axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        Xtr = (Xtr - mu) / sd
        Xte = (Xte - mu) / sd
        ym = y[tr].mean()
        Xt = Xtr.transpose(0, 2, 1)
        G = np.matmul(Xt, Xtr) + eye
        rhs = np.matmul(Xt, (y[tr] - ym)[None, :, None])
        w = np.linalg.solve(G, rhs)
        preds[:, te] = np.matmul(Xte, w)[:, :, 0] + ym
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot > 0:
        r2 = 1.0 - np.sum((y[None, :] - preds) ** 2, axis=1) / ss_tot
    else:
        r2 = np.zeros(R)
    return r2_base, r2 - r2_base


_S = {}


def _init(mrna, prot, base_extra, block, idx, n_refcheck):
    _S.update(mrna=mrna, prot=prot, base_extra=base_extra, block=block, idx=idx, n_refcheck=n_refcheck)


def _job(g):
    y_all = _S["prot"][:, g]
    v = ~np.isnan(y_all)
    nv = int(v.sum())
    y = y_all[v]
    base = _S["mrna"][v, g].reshape(-1, 1)
    if _S["base_extra"] is not None:
        base = np.hstack([base, _S["base_extra"][v]])
    r2_base, inc = batched_increments(y, base, _S["block"], v, _S["idx"])
    folds = tk.folds_for(nv)
    # the batched kernel is checked against the reference kernel on every gene: the observed
    # increment, and the first n_refcheck permutations
    worst = 0.0
    for r in range(min(1 + _S["n_refcheck"], len(inc))):
        ref = ref_increment(y, base, _S["block"][_S["idx"][r]][v], folds)[1]
        worst = max(worst, abs(ref - inc[r]))
    null = inc[1:]
    exceed = int(np.sum(null >= inc[0]))
    return g, nv, r2_base, float(inc[0]), exceed, float(null.mean()) if len(null) else np.nan, worst


def run_block(d, block, base_extra, perms, workers, tag, n_refcheck=2):
    """Per-gene increment of `block` over [own mRNA, base_extra], its permutation p and BH FDR
    over all fit genes of the cohort."""
    n, G = d["protein"].shape
    assert block.shape[0] == n and perms.shape[1] == n
    idx = np.vstack([np.arange(n, dtype=np.int64)[None, :], perms])
    t0 = time.time()
    out = np.full((G, 6), np.nan)
    with mp.Pool(workers, initializer=_init,
                 initargs=(d["mrna_log2tpm1"], d["protein"], base_extra, block, idx, n_refcheck)) as pool:
        for k, r in enumerate(pool.imap_unordered(_job, range(G), chunksize=16)):
            out[r[0]] = r[1:]
            if (k + 1) % 2000 == 0:
                log(f"  {tag}: {k + 1}/{G} genes, {time.time() - t0:.0f} s")
    worst = float(np.nanmax(out[:, 5]))
    if not worst <= TOL_REF:
        sys.exit(f"FATAL: {tag}: batched kernel differs from the reference kernel by {worst:.3e}")
    B = perms.shape[0]
    pval = (1 + out[:, 3]) / (1 + B)
    fdr = bh_fdr(pval)
    res = dict(n=out[:, 0].astype(int), r2_base=out[:, 1], inc=out[:, 2], exceed=out[:, 3].astype(int),
               null_mean=out[:, 4], pval=pval, fdr=fdr, selected=(fdr < FDR_LEVEL) & (out[:, 2] > 0),
               worst_ref=worst, k_block=int(block.shape[1]),
               q_base=1 + (0 if base_extra is None else int(base_extra.shape[1])))
    floor = int(np.sum(out[:, 3] == 0))
    kstar = int(np.ceil(G / (FDR_LEVEL * (B + 1))))
    log(f"  {tag}: block {res['k_block']} cols over a {res['q_base']}-col baseline; selected "
        f"{int(res['selected'].sum())}/{G}; BH rejections before the sign filter "
        f"{int((fdr < FDR_LEVEL).sum())}; genes at the permutation floor {floor} (k* = {kstar}); "
        f"batched-vs-reference worst {worst:.2e}; {time.time() - t0:.0f} s")
    res.update(n_floor=floor, kstar=kstar, n_bh=int((fdr < FDR_LEVEL).sum()))
    return res


# ------------------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------------------
def load_inputs(key):
    c = COHORTS[key]
    path = P(c["inputs"])
    if sha256(path) != c["inputs_sha"]:
        sys.exit(f"FATAL: {c['inputs']} is not the frozen file")
    z = np.load(path, allow_pickle=False)
    if bool(z["is_synthetic"]):
        sys.exit("FATAL: synthetic inputs")
    d = {k: z[k] for k in z.files}
    for k in ("patients", "genes", "label_operator"):
        d[k] = d[k].astype(str)
    d["pcs"] = d["pcs"].astype(np.float64)        # exact upcast, as the published hstack does
    n, G = d["protein"].shape
    if (n, G) != (c["n_disc"], c["n_fit"]) or d["pcs"].shape != (n, 20) or d["rna_pcs20"].shape != (n, 20):
        sys.exit(f"FATAL: unexpected shapes in {c['inputs']}")
    if int(d["ra_n_perm"]) != B_PERM or int(d["ra_seed"]) != 0 or int(d["cv_folds"]) != tk.CV_FOLDS \
            or float(d["ridge_alpha"]) != tk.RIDGE_ALPHA:
        sys.exit("FATAL: estimator constants in the inputs file differ from the published ones")
    if len(set(d["genes"])) != G or np.sum(~np.isnan(d["protein"]), axis=0).min() < 30:
        sys.exit("FATAL: duplicate genes or a gene with fewer than 30 protein values")
    pub_sel = (d["pub_fdr"] < FDR_LEVEL) & (d["pub_incremental_r2"] > 0)
    if int(pub_sel.sum()) != c["n_sel"]:
        sys.exit("FATAL: the published selection set is not the frozen size")
    d["pub_sel"] = pub_sel
    return d


def load_conf(key):
    c = COHORTS[key]
    path = P(c["conf"])
    if sha256(path) != c["conf_sha"]:
        sys.exit(f"FATAL: {c['conf']} is not the frozen file")
    z = np.load(path, allow_pickle=False)
    genes = z["genes"].astype(str)
    out = dict(genes=genes, sel_mask=z["sel_mask"].astype(bool), M={})
    for name, k in c["arms"]:
        M = z[k]
        if M.shape != (c["n_univ"], B_PERM + 1) or not np.isfinite(M).all():
            sys.exit(f"FATAL: unexpected confirmatory matrix {k}")
        out["M"][name] = M
    if len(genes) != c["n_univ"] or len(set(genes)) != len(genes) or int(out["sel_mask"].sum()) != c["n1"]:
        sys.exit(f"FATAL: unexpected confirmatory universe in {c['conf']}")
    return out


def published_perms(d):
    return tk.build_perms(d["protein"].shape[0], B_PERM, "ra", 0)


def operator_dummies(labels, key):
    """PS 2: one dummy for each of the (at most 12) largest operator levels with at least 3
    cases; every other case, including an unlabelled one, has an all-zero row."""
    cnt = collections.Counter(x for x in labels if x != "")
    ordered = sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))
    lv = [k for k, m in ordered[:MAX_OP_LEVELS] if m >= MIN_OP_GROUP]
    if len(ordered) > MAX_OP_LEVELS and ordered[MAX_OP_LEVELS][1] >= MIN_OP_GROUP \
            and ordered[MAX_OP_LEVELS][1] == ordered[MAX_OP_LEVELS - 1][1]:
        sys.exit("FATAL: tie at the operator-dummy cut-off")
    if len(lv) == len(cnt) and lv:
        lv = lv[:-1]
    D = np.zeros((len(labels), len(lv)))
    for j, k in enumerate(lv):
        D[:, j] = (labels == k).astype(float)
    c = COHORTS[key]
    if D.shape[1] != c["n_op_dummies"] or int(D.sum()) != c["n_op_covered"]:
        sys.exit(f"FATAL: operator dummies {D.shape[1]} cols / {int(D.sum())} cases differ from PS 2")
    return D


def load_scores(key, d):
    """Discovery composition scores exported in Stage H (PS 2, PS 4)."""
    path = os.path.join(STAGE_H_DIR, f"c1_spec_discovery_scores_{key}.tsv")
    rows = [ln.rstrip("\n").split("\t") for ln in open(path, encoding="utf-8") if ln.strip() and not ln.startswith("#")]
    head = rows[0]
    if head[:4] != ["patient"] + list(SCORE_COLS):
        sys.exit(f"FATAL: unexpected header in {path}")
    tab = {r[0]: [float(x) for x in r[1:4]] for r in rows[1:]}
    if len(tab) != len(rows) - 1 or set(tab) != set(d["patients"]):
        sys.exit(f"FATAL: {path} does not hold exactly the cohort's discovery patients")
    C = np.asarray([tab[p] for p in d["patients"]], dtype=np.float64)
    if not np.isfinite(C).all():
        sys.exit(f"FATAL: non-finite composition score in {path}")
    return C, path


def load_stage_h(key, arm, conf):
    path = os.path.join(STAGE_H_DIR, f"c1_spec_{key}_{arm}.npz")
    z = np.load(path, allow_pickle=False)
    genes = z["genes"].astype(str)
    M = z["M"]
    if list(genes) != list(conf["genes"]) or M.shape != (len(genes), B_PERM + 1) or not np.isfinite(M).all():
        sys.exit(f"FATAL: {path} is not on the published confirmatory universe")
    return M, path


# ------------------------------------------------------------------------------------------
# cells (PS 5)
# ------------------------------------------------------------------------------------------
def cell(M, universe, set_genes, xs_by_gene=None, background=None):
    """AUC of the set against its background on one confirmatory matrix.  background=None: every
    other gene of the universe (the C1-LUAD H2 rule)."""
    sel = np.fromiter((g in set_genes for g in universe), bool, len(universe))
    if background is None:
        keep = np.ones(len(universe), bool)
    else:
        keep = sel | np.fromiter((g in background for g in universe), bool, len(universe))
    n1, n0 = int(sel[keep].sum()), int((~sel[keep]).sum())
    if n1 < MIN_SET or n0 < MIN_SET:
        return dict(status="untestable", n1=n1, n0=n0)
    Mk = M[keep]
    if xs_by_gene is None:
        xs = np.arange(Mk.shape[0], dtype=float)
    else:
        xs = np.asarray([xs_by_gene[g] for g in universe[keep]], dtype=float)
    auc, rho, n1, n0 = assemble(Mk, xs, sel[keep])
    sd = float(auc[1:].std())
    out = dict(status="ok", n1=n1, n0=n0, auc=float(auc[0]), null_mean=float(auc[1:].mean()), null_sd=sd,
               margin=float(auc[0] - auc[1:].mean()), z=float((auc[0] - auc[1:].mean()) / sd) if sd > 0 else None,
               p=emp_p(auc[0], auc[1:]), exceed=int(np.sum(auc[1:] >= auc[0])))
    if xs_by_gene is not None:
        out.update(rho=float(rho[0]), rho_null_mean=float(rho[1:].mean()), rho_null_sd=float(rho[1:].std()),
                   p_rho=emp_p(rho[0], rho[1:]))
    return out


def partition_cells(M, universe, sm, sz):
    """PS 5, descriptive: shared, histology-only and covariate-only genes, each against the
    genes that neither analysis selected."""
    uni = set(universe)
    sm, sz = sm & uni, sz & uni
    neither = uni - sm - sz
    return {"shared": cell(M, universe, sm & sz, background=neither),
            "histology_only": cell(M, universe, sm - sz, background=neither),
            "covariate_only": cell(M, universe, sz - sm, background=neither)}


def set_summary(d, res, universe):
    g = d["genes"]
    sel = set(g[res["selected"]])
    return dict(n_fit=int(len(g)), n_selected=len(sel), n_in_universe=len(sel & set(universe)),
                n_shared_with_published=len(sel & set(g[d["pub_sel"]])),
                n_bh_before_sign_filter=res["n_bh"], n_floor=res["n_floor"], kstar=res["kstar"],
                k_block=res["k_block"], q_base=res["q_base"], batched_vs_reference_worst=res["worst_ref"])


def write_sets(path, d, results):
    names = list(results)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        cols = ["gene", "n", "pub_incremental_r2", "pub_pval", "pub_fdr", "pub_selected"]
        for nm in names:
            cols += [f"{nm}_increment", f"{nm}_pval", f"{nm}_fdr", f"{nm}_selected"]
        fh.write(",".join(cols) + "\n")
        for i, g in enumerate(d["genes"]):
            row = [g, str(int(d["pub_n"][i])), repr(float(d["pub_incremental_r2"][i])),
                   repr(float(d["pub_pval"][i])), repr(float(d["pub_fdr"][i])), str(int(d["pub_sel"][i]))]
            for nm in names:
                r = results[nm]
                row += [repr(float(r["inc"][i])), repr(float(r["pval"][i])), repr(float(r["fdr"][i])),
                        str(int(r["selected"][i]))]
            fh.write(",".join(row) + "\n")


def read_sets(path, names):
    rows = [ln.rstrip("\n").split(",") for ln in open(path, encoding="utf-8")]
    head = rows[0]
    out = {}
    for nm in names:
        j = head.index(f"{nm}_increment")
        k = head.index(f"{nm}_selected")
        out[nm] = dict(genes=set(r[0] for r in rows[1:] if r[k] == "1"),
                       xs={r[0]: float(r[j]) for r in rows[1:]})
    return out


# ------------------------------------------------------------------------------------------
# validation (published quantities only)
# ------------------------------------------------------------------------------------------
def validate_v1(key, d, workers):
    res = run_block(d, d["pcs"], None, published_perms(d), workers, f"V1 {key}")
    dif = dict(increment=float(np.max(np.abs(res["inc"] - d["pub_incremental_r2"]))),
               r2_rna=float(np.max(np.abs(res["r2_base"] - d["pub_r2_rna"]))),
               pval=float(np.max(np.abs(res["pval"] - d["pub_pval"]))),
               fdr=float(np.max(np.abs(res["fdr"] - d["pub_fdr"]))),
               n=int(np.max(np.abs(res["n"] - d["pub_n"]))))
    same_set = bool(np.array_equal(res["selected"], d["pub_sel"]))
    ok = (dif["r2_rna"] <= TOL_PUB and dif["increment"] <= TOL_REF and dif["pval"] <= 1e-12
          and dif["fdr"] <= 1e-12 and dif["n"] == 0 and same_set)
    log(f"V1 {key}: max abs differences {dif}; selection set identical: {same_set}; "
        f"selected {int(res['selected'].sum())} -> {'VALIDATED' if ok else 'FAILED'}")
    return dict(ok=ok, max_abs_diff=dif, selection_identical=same_set, n_selected=int(res["selected"].sum()))


def validate_v2(key, d, conf):
    c = COHORTS[key]
    pub = json.load(open(P(c["published"]), encoding="utf-8"))
    xs = dict(zip(d["genes"], d["pub_incremental_r2"]))
    sm = set(d["genes"][d["pub_sel"]])
    mask = np.fromiter((g in sm for g in conf["genes"]), bool, len(conf["genes"]))
    ok = bool(np.array_equal(mask, conf["sel_mask"]))
    rows = {}
    for name, _ in c["arms"]:
        r = cell(conf["M"][name], conf["genes"], sm, xs_by_gene=xs)
        if key == "ucec":
            w = pub["arms"][name]
            want = dict(auc=w["auc"]["obs"], null_mean=w["auc"]["null_mean"], null_sd=w["auc"]["null_sd"],
                        p=w["auc"]["p_one_sided"], rho=w["rho"]["obs"], p_rho=w["rho"]["p_one_sided"])
        else:
            w = pub["h1"]
            want = dict(auc=w["auc_obs"], null_mean=w["auc_null_mean"], null_sd=w["auc_null_sd"],
                        p=w["p_auc"], rho=w["rho_obs"], p_rho=w["p_rho"])
        worst = max(abs(r[k] - want[k]) for k in want)
        good = worst <= TOL_AUC and (r["n1"], r["n0"]) == (c["n1"], c["n0"])
        ok = ok and good
        rows[name] = dict(got={k: r[k] for k in want}, n1=r["n1"], n0=r["n0"], worst_abs_diff=worst, ok=good)
        log(f"V2 {key} {name}: AUC {r['auc']:.10f} null {r['null_mean']:.10f} SD {r['null_sd']:.10f} "
            f"p {r['p']:.6f} rho {r['rho']:.10f}; n {r['n1']}/{r['n0']}; worst diff {worst:.1e} "
            f"-> {'VALIDATED' if good else 'FAILED'}")
    return dict(ok=ok, mask_identical=bool(np.array_equal(mask, conf["sel_mask"])), arms=rows)


def validate_v3(key, d, C):
    """The transferred discovery scores reproduce the published composition controls
    (results__composition_controls_<cohort>.csv): for every published-set gene, the gain of
    each score alone and of the three together over the gene's own mRNA.  These four
    columns involve no morphology component, so they test the scores and nothing else."""
    c = COHORTS[key]
    rows = [ln.rstrip("\n").split(",") for ln in open(P(c["compctl"]), encoding="utf-8")]
    head = rows[0]
    jg = head.index("gene")
    cols = [(head.index(f"{f}_gain"), [j]) for j, f in enumerate(SCORE_COLS)] + [(head.index("ALL_gain"), [0, 1, 2])]
    gi = {g: i for i, g in enumerate(d["genes"])}
    worst, m = 0.0, 0
    for r in rows[1:]:
        if r[jg] not in gi or any(r[j] == "" for j, _ in cols):
            continue
        g = gi[r[jg]]
        v = ~np.isnan(d["protein"][:, g])
        y = d["protein"][v, g]
        folds = tk.folds_for(int(v.sum()))
        xr = d["mrna_log2tpm1"][v, g].reshape(-1, 1)
        r2_r = tk.r2_of(y, tk.oof_ridge(xr, y, folds))
        for j, which in cols:
            r2_c = tk.r2_of(y, tk.oof_ridge(np.hstack([xr, C[v][:, which]]), y, folds))
            worst = max(worst, abs((r2_c - r2_r) - float(r[j])))
        m += 1
    ok = m == len(rows) - 1 and m >= 0.99 * c["n_sel"] and worst <= TOL_REF
    log(f"V3 {key}: {m} of {len(rows) - 1} published composition-control genes, worst abs difference "
        f"{worst:.2e} -> {'VALIDATED' if ok else 'FAILED'}")
    return dict(ok=ok, n_genes=m, n_rows=len(rows) - 1, worst_abs_diff=worst)


# ------------------------------------------------------------------------------------------
# stages
# ------------------------------------------------------------------------------------------
def tiers(p):
    return "p<0.05" if p < 0.05 else ("0.05<=p<0.15" if p < 0.15 else "p>=0.15")


def add_cell(cells, cohort, family, set_name, matrix, arm, r):
    row = dict(cohort=cohort, family=family, set=set_name, matrix=matrix, null=arm)
    row.update(r)
    if r.get("status") == "ok":
        row["tier"] = tiers(r["p"])
    cells.append(row)


def write_cells(path, cells):
    cols = ["cohort", "family", "set", "matrix", "null", "status", "n1", "n0", "auc", "null_mean", "null_sd",
            "margin", "z", "p", "exceed", "tier", "rho", "rho_null_mean", "rho_null_sd", "p_rho"]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(",".join(cols) + "\n")
        for r in cells:
            fh.write(",".join("" if r.get(k) is None else (repr(r[k]) if isinstance(r[k], float) else str(r[k]))
                              for k in cols) + "\n")


def stage_l1(workers, out_dir):
    cells, summary = [], {}
    for key in COHORTS:
        d, conf = load_inputs(key), load_conf(key)
        perms = published_perms(d)
        D = operator_dummies(d["label_operator"], key)
        res = {"S_R20": run_block(d, d["rna_pcs20"], None, perms, workers, f"{key} S_R20"),
               "S_OP": run_block(d, D, None, perms, workers, f"{key} S_OP"),
               "S_MgR20": run_block(d, d["pcs"], d["rna_pcs20"], perms, workers, f"{key} S_M|R20")}
        write_sets(os.path.join(out_dir, f"sets_L1_{key}.csv"), d, res)
        uni = conf["genes"]
        sm = set(d["genes"][d["pub_sel"]])
        summary[key] = {nm: set_summary(d, r, uni) for nm, r in res.items()}
        for arm, _ in COHORTS[key]["arms"]:
            M = conf["M"][arm]
            add_cell(cells, key, "reference", "S_M", "M_pub", arm,
                     cell(M, uni, sm, xs_by_gene=dict(zip(d["genes"], d["pub_incremental_r2"]))))
            for nm in ("S_R20", "S_OP"):
                sz = set(d["genes"][res[nm]["selected"]])
                add_cell(cells, key, "forward", nm, "M_pub", arm,
                         cell(M, uni, sz, xs_by_gene=dict(zip(d["genes"], res[nm]["inc"]))))
                for part, r in partition_cells(M, uni, sm, sz).items():
                    add_cell(cells, key, "partition", f"{nm}:{part}", "M_pub", arm, r)
    write_cells(os.path.join(out_dir, "cells_L1.csv"), cells)
    return dict(stage="L1", sets=summary, n_cells=len(cells))


def stage_l2(workers, out_dir):
    l1_dir = os.path.join(OUT_DIR, "L1")
    cells, summary, gates, files = [], {}, {}, {}
    for key in COHORTS:
        d, conf = load_inputs(key), load_conf(key)
        perms = published_perms(d)
        C, spath = load_scores(key, d)
        files[os.path.relpath(spath, PAPER).replace(os.sep, "/")] = sha256(spath)
        gates[key] = validate_v3(key, d, C)
        if not gates[key]["ok"]:
            sys.exit(f"FATAL: V3 failed for {key}: the transferred scores are not the published ones")
        res = {"S_C3": run_block(d, C, None, perms, workers, f"{key} S_C3"),
               "S_MgC3": run_block(d, d["pcs"], C, perms, workers, f"{key} S_M|C3")}
        write_sets(os.path.join(out_dir, f"sets_L2_{key}.csv"), d, res)
        uni = conf["genes"]
        summary[key] = {nm: set_summary(d, r, uni) for nm, r in res.items()}
        sets = {nm: dict(genes=set(d["genes"][r["selected"]]), xs=dict(zip(d["genes"], r["inc"])))
                for nm, r in res.items()}
        sets.update(read_sets(os.path.join(l1_dir, f"sets_L1_{key}.csv"), ("S_R20", "S_MgR20")))
        sets["S_M"] = dict(genes=set(d["genes"][d["pub_sel"]]), xs=dict(zip(d["genes"], d["pub_incremental_r2"])))
        H = {}
        for arm in STAGE_H_ARMS:
            H[arm], hp = load_stage_h(key, arm, conf)
            files[os.path.relpath(hp, PAPER).replace(os.sep, "/")] = sha256(hp)

        def rd(family, set_name, matrix, M, null="unrestricted"):
            add_cell(cells, key, family, set_name, matrix, null,
                     cell(M, uni, sets[set_name]["genes"], xs_by_gene=sets[set_name]["xs"]))

        for arm, _ in COHORTS[key]["arms"]:
            rd("forward", "S_C3", "M_pub", conf["M"][arm], arm)
            for part, r in partition_cells(conf["M"][arm], uni, sets["S_M"]["genes"], sets["S_C3"]["genes"]).items():
                add_cell(cells, key, "partition", f"S_C3:{part}", "M_pub", arm, r)
        rd("primary", "S_M", "A_C3", H["A_C3"])
        rd("primary", "S_M", "A_R20", H["A_R20"])
        rd("capacity", "S_M", "A_R20cap", H["A_R20cap"])
        rd("matched", "S_MgC3", "A_C3", H["A_C3"])
        rd("matched", "S_MgR20", "A_R20", H["A_R20"])
        rd("reverse", "S_M", "V_C3", H["V_C3"])
        rd("reverse", "S_M", "V_R20", H["V_R20"])
        rd("covariate-to-covariate", "S_C3", "V_C3", H["V_C3"])
        rd("covariate-to-covariate", "S_R20", "V_R20", H["V_R20"])
    write_cells(os.path.join(out_dir, "cells_L2.csv"), cells)
    return dict(stage="L2", sets=summary, gates=gates, stage_h_files=files, n_cells=len(cells))


# ------------------------------------------------------------------------------------------
# selftest (synthetic)
# ------------------------------------------------------------------------------------------
def selftest():
    from scipy.stats import mannwhitneyu
    rng = np.random.RandomState(1)
    n, G, k = 60, 40, 7
    mrna = rng.randn(n, G)
    block = rng.randn(n, k)
    block[:, 0] = (rng.rand(n) < 0.1).astype(float)            # a sparse dummy-like column
    extra = rng.randn(n, 3)
    prot = 0.5 * mrna + 0.4 * block[:, [1]] + rng.randn(n, G)
    prot[rng.rand(n, G) < 0.1] = np.nan
    perms = tk.build_perms(n, 25, "ra", 0)
    idx = np.vstack([np.arange(n)[None, :], perms])
    worst = 0.0
    for g in range(G):
        v = ~np.isnan(prot[:, g])
        y = prot[v, g]
        folds = tk.folds_for(int(v.sum()))
        for base in (mrna[v, g].reshape(-1, 1), np.hstack([mrna[v, g].reshape(-1, 1), extra[v]])):
            r2b, inc = batched_increments(y, base, block, v, idx)
            for r in range(idx.shape[0]):
                rb, ri = ref_increment(y, base, block[idx[r]][v], folds)
                worst = max(worst, abs(ri - inc[r]), abs(rb - r2b))
        nv, r2r, r2both, pinc = tk.published_increment(mrna[:, g], prot[:, g], block, v)
        worst = max(worst, abs(pinc - batched_increments(y, mrna[v, g].reshape(-1, 1), block, v, idx[:1])[1][0]))
        null = tk.published_null_increments(mrna[:, g], prot[:, g], block, perms, v)
        worst = max(worst, float(np.max(np.abs(null - batched_increments(y, mrna[v, g].reshape(-1, 1), block, v, idx)[1][1:]))))
    assert worst <= 1e-10, worst
    p = rng.rand(500)
    p[:40] *= 1e-4
    m = len(p)
    direct = np.array([min(1.0, min(p[p >= pi] * m / np.array([np.sum(p <= pj) for pj in p[p >= pi]]))) for pi in p])
    assert np.max(np.abs(bh_fdr(p) - direct)) <= 1e-12
    M = rng.randn(300, 11)
    sel = rng.rand(300) < 0.3
    auc, rho, n1, n0 = assemble(M, rng.randn(300), sel)
    for j in range(M.shape[1]):
        u = mannwhitneyu(M[sel, j], M[~sel, j]).statistic / (n1 * n0)
        assert abs(u - auc[j]) <= 1e-12
    uni = np.array([f"g{i}" for i in range(300)])
    c = cell(M, uni, set(uni[sel]))
    assert abs(c["auc"] - auc[0]) <= 1e-12 and c["n1"] == n1
    assert cell(M, uni, set(uni[:10]))["status"] == "untestable"
    labels = np.array(["a"] * 5 + ["b"] * 4 + ["c"] * 2 + [""] * 1)
    cnt = collections.Counter(x for x in labels if x != "")
    assert sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0] == ("a", 5)
    print("selftest passed; batched-vs-reference worst", worst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--stage", choices=["L1", "L2"])
    ap.add_argument("--workers", type=int, default=max(1, min(24, (os.cpu_count() or 2) - 2)))
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    if sha256(tk.__file__) != KERNEL_SHA:
        sys.exit("FATAL: theta_kernel.py is not the frozen kernel")
    if args.validate_only:
        out = {}
        for key in COHORTS:
            d, conf = load_inputs(key), load_conf(key)
            out[key] = dict(V2=validate_v2(key, d, conf), V1=validate_v1(key, d, args.workers))
        ok = all(v["V1"]["ok"] and v["V2"]["ok"] for v in out.values())
        os.makedirs(OUT_DIR, exist_ok=True)
        with open(os.path.join(OUT_DIR, "validate_only.json"), "w", encoding="utf-8", newline="\n") as fh:
            json.dump(dict(ok=ok, script_sha256=sha256(os.path.abspath(__file__)), cohorts=out), fh, indent=1)
        log("VALIDATED" if ok else "VALIDATION FAILED")
        sys.exit(0 if ok else 1)
    if not args.stage:
        ap.error("one of --selftest, --validate-only, --stage is required")
    check_sidecar(required_files(args.stage))
    out_dir = os.path.join(OUT_DIR, args.stage)
    if os.path.exists(out_dir):
        sys.exit(f"REFUSED: {out_dir} exists -- a completed stage is never repeated (PS 9)")
    if args.stage == "L2" and not os.path.exists(os.path.join(OUT_DIR, "L1", "result_L1.json")):
        sys.exit("REFUSED: Stage L1 has not completed")
    os.makedirs(out_dir)
    t0 = time.time()
    res = stage_l1(args.workers, out_dir) if args.stage == "L1" else stage_l2(args.workers, out_dir)
    res.update(script_sha256=sha256(os.path.abspath(__file__)), sidecar_sha256=sha256(SIDECAR),
               minutes=round((time.time() - t0) / 60, 1), workers=args.workers)
    with open(os.path.join(out_dir, f"result_{args.stage}.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(res, fh, indent=1)
    log(f"Stage {args.stage} complete: {res['n_cells']} cells, {res['minutes']} min")


if __name__ == "__main__":
    main()
