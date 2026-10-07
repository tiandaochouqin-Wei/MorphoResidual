#!/usr/bin/env python
"""Post hoc C1-UCEC re-reads with alternative discovery selection sets.

Pre-specified in review/POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md. The full run is
fail-closed: it refuses to run unless the sidecar
review/POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md.sha256 exists and every file listed in it
(the prespec, this script, every input) still has the listed sha256.

Local only. Reads the permutation matrices saved by the D4 stratified re-read
(pinned/c1_ucec_d4_stratified.npz: confirmatory incremental R^2 for 10,146 usable genes;
column 0 observed, columns 1..1000 patient<->slide permutations; three null arms) and
the discovery per-gene tables of the post hoc batch-axis and data-subset arms.

--validate-only (touches only the published 2,566-gene set):
  V1  the three saved arms reproduce c1_ucec_d4_followup.json exactly (|diff| <= 1e-12)
      for AUC and rho: observed, null mean, null SD, one-sided p; n 2,536 / 7,610; npz sha256.
  V2  the code path used for every new set (read_set -> universe -> mask -> xs -> assemble),
      fed the published discovery table, reproduces the published mask, xs and V1 numbers.
Statistic code is copied unchanged from c1_run_test.assemble.
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))                     # server_export/
PAPER = os.path.abspath(os.path.join(ROOT, ".."))                    # MorphoResidual_paper/
NPZ = os.path.join(HERE, "pinned", "c1_ucec_d4_stratified.npz")
FOLLOWUP = os.path.join(HERE, "pinned", "c1_ucec_d4_followup.json")
GENE_LEVEL = os.path.join(ROOT, "pinned", "ucec_c1", "c1_ucec_gene_level.csv")
PUBLISHED_DISC = "pinned/ucec/residual_results_tumoronly.csv"
SIDECAR = os.path.join(PAPER, "review", "POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md.sha256")
OUT_DIR = os.path.join(ROOT, "results", "c1_ucec_posthoc_sets")
NPZ_SHA = "ef75f1c06d4fc33843a79e9a8d1bb4f4224acd201d30784badb97cf3be13384a"
ARMS = [("unrestricted", "M_unrestricted"), ("within-operator", "M_within_operator"),
        ("residualised", "M_within_operator_residualized")]
ARM_PRIMARY = ("unrestricted", "within-operator")                    # both must pass
TOL = 1e-12

PRIMARY = {"S6": ("pinned/posthoc/ucec_ffpe_only.csv", "fdr", "incremental_r2")}
SECONDARY = {
    "S1": ("results/ucec__sitepack_operator.csv", "fdr_batchresid", "incremental_r2_batchresid"),
    "S2": ("results/ucec__sitepack_operator.csv", "fdr_both", "incremental_r2_both"),
    "S3": ("pinned/posthoc/ucec_medium.csv", "fdr_batchresid", "incremental_r2_batchresid"),
    "S4": ("pinned/posthoc/ucec_plex.csv", "fdr_batchresid", "incremental_r2_batchresid"),
    "S5": ("pinned/posthoc/ucec_plex.csv", "fdr_both", "incremental_r2_both"),
    "S7": ("pinned/posthoc/ucec_arm_Q.csv", "fdr", "incremental_r2"),
}
UNTESTABLE = {"S3both": ("pinned/posthoc/ucec_medium.csv", "fdr_both", "incremental_r2_both")}
REFERENCE = {f"D{i:02d}": (f"pinned/posthoc/ucec_arm_D{i:02d}.csv", "fdr", "incremental_r2")
             for i in range(1, 20)}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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
    if not np.isfinite(obs) or not np.isfinite(null).all():
        sys.exit("FATAL: non-finite statistic")
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))


def read_table(path, fdr_col, inc_col):
    df = pd.read_csv(os.path.join(ROOT, path))
    if "cohort" in df.columns:
        df = df[df["cohort"].str.lower() == "ucec"]
    if df["gene"].duplicated().any():
        sys.exit(f"FATAL: duplicated genes in {path}")
    df = df.set_index("gene")
    df = df[df[fdr_col].notna() & df[inc_col].notna()]
    selected = set(df.index[(df[fdr_col] < 0.05) & (df[inc_col] > 0)])
    return df, selected


def read_set(z, path, fdr_col, inc_col, exclude=None):
    """Universe = npz genes tested in this discovery table; mask = selected (minus exclude)."""
    genes = np.asarray(z["genes"])
    df, selected = read_table(path, fdr_col, inc_col)
    keep = np.isin(genes, df.index.values)
    g_u = genes[keep]
    if exclude is not None:
        selected = selected - set(exclude)
    sel = np.isin(g_u, list(selected))
    xs = df.reindex(g_u)[inc_col].values.astype(float)
    if not np.isfinite(xs).all():
        sys.exit(f"FATAL: non-finite discovery increments in {path}")
    return keep, sel, xs, g_u


def stats_for(z, keep, sel, xs):
    out = {}
    for arm, key in ARMS:
        M = z[key][keep]
        if not np.isfinite(M).all():
            sys.exit(f"FATAL: non-finite values in {key}")
        auc, rho, n1, n0 = assemble(M, xs, sel)
        out[arm] = dict(n_selected=n1, n_background=n0,
                        auc_obs=float(auc[0]), auc_null_mean=float(auc[1:].mean()),
                        auc_null_sd=float(auc[1:].std()), p_auc=emp_p(auc[0], auc[1:]),
                        margin=float(auc[0] - auc[1:].mean()),
                        rho_obs=float(rho[0]), rho_null_mean=float(rho[1:].mean()),
                        rho_null_sd=float(rho[1:].std()), p_rho=emp_p(rho[0], rho[1:]))
    return out


def check_equal(label, got, want):
    if abs(got - want) > TOL:
        sys.exit(f"FATAL: validation {label}: {got!r} != {want!r}")


def validate(z):
    if sha256(NPZ) != NPZ_SHA:
        sys.exit("FATAL: npz sha256 changed")
    follow = json.load(open(FOLLOWUP))
    gl = pd.read_csv(GENE_LEVEL)
    if gl["gene"].duplicated().any():
        sys.exit("FATAL: duplicated genes in gene-level table")
    gl = gl.set_index("gene").reindex(z["genes"])
    xs_pub = gl["discovery_incremental_r2"].values.astype(float)
    if not np.isfinite(xs_pub).all():
        sys.exit("FATAL: non-finite published discovery increments")
    # the CSV round-trip differs from the npz in the last ulp (max |diff| 4.4e-16)
    if not np.allclose(z["M_unrestricted"][:, 0], gl["confirmatory_incremental_r2"].values, rtol=0, atol=1e-14):
        sys.exit("FATAL: npz column 0 != published confirmatory increments")
    # V1
    v1 = stats_for(z, np.ones(len(z["genes"]), bool), np.asarray(z["sel_mask"]), xs_pub)
    for arm, _ in ARMS:
        ref = follow["arms"][arm]
        for stat, key in (("auc", "auc"), ("rho", "rho")):
            check_equal(f"{arm} {stat} obs", v1[arm][f"{key}_obs"], ref[stat]["obs"])
            check_equal(f"{arm} {stat} null_mean", v1[arm][f"{key}_null_mean"], ref[stat]["null_mean"])
            check_equal(f"{arm} {stat} null_sd", v1[arm][f"{key}_null_sd"], ref[stat]["null_sd"])
            check_equal(f"{arm} {stat} p", v1[arm][f"p_{key}"], ref[stat]["p_one_sided"])
    if (v1["unrestricted"]["n_selected"], v1["unrestricted"]["n_background"]) != (2536, 7610):
        sys.exit("FATAL: published n_selected/n_background not 2536/7610")
    # V2: the new-set code path on the published discovery table
    keep, sel, xs, _ = read_set(z, PUBLISHED_DISC, "fdr", "incremental_r2")
    if not keep.all() or not np.array_equal(sel, np.asarray(z["sel_mask"])):
        sys.exit("FATAL: read_set does not reproduce the published universe/mask")
    if not np.allclose(xs, xs_pub, rtol=0, atol=TOL):
        sys.exit("FATAL: read_set does not reproduce the published discovery increments")
    v2 = stats_for(z, keep, sel, xs)
    for arm, _ in ARMS:
        for k in ("auc_obs", "auc_null_mean", "auc_null_sd", "p_auc", "rho_obs", "p_rho"):
            check_equal(f"V2 {arm} {k}", v2[arm][k], v1[arm][k])
    print("VALIDATED: V1 (three arms, exact) and V2 (new-set code path) reproduce the published reading")
    for arm, _ in ARMS:
        print(f"  {arm:16s} AUC {v1[arm]['auc_obs']:.10f} null {v1[arm]['auc_null_mean']:.10f}"
              f" sd {v1[arm]['auc_null_sd']:.10f} p {v1[arm]['p_auc']:.6f} rho {v1[arm]['rho_obs']:.10f}")
    return v1


def check_sidecar():
    if not os.path.exists(SIDECAR):
        sys.exit(f"REFUSED: sidecar {SIDECAR} not found -- the prespec is not frozen.")
    for line in open(SIDECAR, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        want, rel = line.split(None, 1)
        path = os.path.join(PAPER, rel.lstrip("*"))
        if not os.path.exists(path) or sha256(path) != want:
            sys.exit(f"REFUSED: {rel} missing or changed since freezing")
    print("sidecar verified")


def holm(pvals):
    order = np.argsort(pvals)
    m, running, adj = len(pvals), 0.0, np.empty(len(pvals))
    for k, i in enumerate(order):
        running = max(running, min(1.0, (m - k) * pvals[i]))
        adj[i] = running
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate-only", action="store_true")
    args = ap.parse_args()
    z = np.load(NPZ, allow_pickle=False)
    v1 = validate(z)
    if args.validate_only:
        return
    check_sidecar()

    published = set(np.asarray(z["genes"])[np.asarray(z["sel_mask"])])
    rows = []

    def record(sid, role, path, fdr_col, inc_col, exclude=None, background_note="universe minus set"):
        keep, sel, xs, g_u = read_set(z, path, fdr_col, inc_col, exclude=exclude)
        if sel.sum() == 0:
            rows.append(dict(set=sid, role=role, arm="", n_selected=0, note="untestable: 0 genes selected",
                             source=path))
            return
        res = stats_for(z, keep, sel, xs)
        overlap = int(np.isin(g_u[sel], list(published)).sum())
        for arm, r in res.items():
            rows.append(dict(set=sid, role=role, arm=arm, universe=int(keep.sum()),
                             overlap_with_published=overlap, background=background_note,
                             source=path, source_sha256=sha256(os.path.join(ROOT, path)), **r))

    def record_new(sid, role, path, fdr_col, inc_col):
        """Sk-new: Sk minus the published set, against background universe minus Sk (AUC only)."""
        keep, sel_full, xs, g_u = read_set(z, path, fdr_col, inc_col)
        sel_new = sel_full & ~np.isin(g_u, list(published))
        drop = sel_full & ~sel_new                    # shared genes removed from both sides
        keep2 = keep.copy()
        keep2[np.where(keep)[0][drop]] = False
        sel2 = sel_new[~drop]
        if sel2.sum() == 0:
            rows.append(dict(set=sid, role=role, arm="", n_selected=0, note="no genes outside the published set",
                             source=path))
            return
        res = stats_for(z, keep2, sel2, xs[~drop])
        for arm, r in res.items():
            r = {k: v for k, v in r.items() if not k.startswith("rho") and k != "p_rho"}
            rows.append(dict(set=sid, role=role, arm=arm, universe=int(keep2.sum()),
                             background="universe minus full set", source=path, **r))

    for sid, spec in PRIMARY.items():
        record(sid, "primary", *spec)
    for sid, spec in SECONDARY.items():
        record(sid, "secondary", *spec)
    for sid, spec in UNTESTABLE.items():
        record(sid, "untestable", *spec)
    for sid, spec in {**PRIMARY, **SECONDARY}.items():
        record_new(sid + "new", "descriptive", *spec)
    for sid, spec in REFERENCE.items():
        record(sid, "reference", *spec)
    out = pd.DataFrame(rows)

    # secondary family: p_set = max(p_unrestricted, p_within-operator); Holm over S1-S5, S7
    sec = out[(out["role"] == "secondary") & (out["arm"].isin(ARM_PRIMARY))]
    pset = sec.groupby("set")["p_auc"].max()
    adj = pd.Series(holm(pset.values), index=pset.index)
    out["p_set_secondary"] = out["set"].map(pset)
    out["p_set_holm_secondary"] = out["set"].map(adj)
    # reference comparison for S6: rank of S6's within-operator margin among S6 + D01..D19
    wo = out[(out["arm"] == "within-operator") & out["role"].isin(["primary", "reference"])]
    s6_margin = float(wo.loc[wo["set"] == "S6", "margin"].iloc[0])
    d_margins = wo.loc[wo["role"] == "reference", "margin"].values
    s6_rank = dict(s6_margin=s6_margin, n_reference=int(len(d_margins)),
                   n_reference_below_s6=int((d_margins < s6_margin).sum()),
                   reference_min=float(d_margins.min()), reference_max=float(d_margins.max()))
    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(os.path.join(OUT_DIR, "result.csv"), index=False)
    json.dump(dict(validation_published=v1, s6_reference_rank=s6_rank, npz_sha256=NPZ_SHA,
                   script_sha256=sha256(__file__), rows=out.to_dict(orient="records")),
              open(os.path.join(OUT_DIR, "result.json"), "w"), indent=1, default=str)
    print(out[["set", "role", "arm", "n_selected", "auc_obs", "auc_null_mean", "margin", "p_auc"]]
          .to_string(index=False))
    print("S6 reference rank:", s6_rank)
    print("wrote", OUT_DIR)


if __name__ == "__main__":
    main()
