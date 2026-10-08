#!/usr/bin/env python
"""Set-level read-out for the BRCA (CPTAC-2, PDC000173) replication arms: AUC, Spearman rho, one-sided permutation p.

Same statistic and same code path as server_export/scripts/c1_ucec_posthoc_sets.py / c1_run_test.assemble (rank-based
Mann-Whitney AUC of the selected set against the rest of the universe; Spearman rho between the discovery score and the BRCA
increment as the mean product of standardised ranks; one-sided empirical p = (#null >= obs + 1) / (B + 1), with the null being the
B = 1000 patient<->slide permutations that cptac2_residual.py (patch v3) saved column-wise, every gene on the same permutation, so
gene-gene correlation is preserved).

FAIL-CLOSED. A real run refuses unless the sidecar
    review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256
exists and (a) every file it lists still has the listed sha256, (b) it lists THIS script, the primary set file, every secondary set
file and at least one .md pre-specification. The matrices need not be in the sidecar (they do not exist when the pre-specification is
frozen); their sha256 is recorded in the output instead. Nothing in this file can be changed after freezing without REFUSED.

The constants below ARE the pre-specification (POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md): PRIMARY_SET, SECONDARY_SETS, PRIMARY_ARM,
EXPECT_B, the DX gate tolerances and the case-count rules. They are frozen by the sidecar; any change after freezing is REFUSED.

Universe (review M1): the BRCA-tested genes (rows of an arm's matrix) among the 8,309 genes tested in all five CPTAC-3
discovery cohorts (= the rows of every frozen set file); "background" = the other genes of that universe.

v3 (2026-10-08, adversarial review of the pre-specification, before freezing):
  B2  Published-file check (cross-review): before the DX gate tiers are applied, the published results/brca_protein_results.csv and
      results/summary_brca.csv must each hash to their line in dx_arm_inputs_outputs.sha256 (written on the server at RUNBOOK step 0,
      brought back, passed with --dx-published-sha256), and the published summary must read n_cases 102, tested 10,031, significant 0
      and negative_frac_all rounding to 0.991 (constants PUB_* below). Otherwise the gate is "fail". The script also refuses to
      overwrite an existing setlevel_result.* (the analysis is run once; an invocation that exits FATAL before writing is not a run).
  M3  DX reproduction gate, three tiers, written into every dx output row and into the json:
        exact : same gene set and max|published - re-run increment| <= 1e-12, and the summary (n_cases, tested, significant,
                negative_frac_all to 4 decimals) identical;
        host  : same gene set, max|delta| <= 1e-4, the same summary identity;
        fail  : anything else -> dx rows labelled "DX re-run (not identical to published)", supplement only.
      The re-run summary must belong to the dx matrix (tested = rows, negative fraction = mean(col 0 < 0)) or the run stops.
  M5  Case sets: dx must hold exactly 102 distinct cases; ts and ts_tile must hold the same cases, a subset of dx; a TS case
      without a usable slide leaves the TS arms only; n_ts >= 95 -> main text, n_ts < 95 -> supplement only.
  N5  A SECONDARY set covered < 50% in the matrix is reported "untestable (coverage)"; the PRIMARY set still stops the run.

Usage (local, after the bring-back; matrices are 3 x ~80 MB):
  python local_brca_setlevel.py --matrix ts=.../brca_incr_matrix_ts.npz --matrix dx=.../brca_incr_matrix.npz \\
         --matrix ts_tile=.../brca_incr_matrix_ts_tilepool.npz \\
         --dx-published-csv .../results/brca_protein_results.csv --dx-published-summary .../results/summary_brca.csv \\
         --dx-repro-summary .../results_dx_repro/summary_brca.csv --dx-published-sha256 .../dx_arm_inputs_outputs.sha256
  python local_brca_setlevel.py --selftest        # synthetic data only; also exercises the fail-closed paths
Arm labels: ts = TS slide-mean (primary), dx = DX re-run, ts_tile = TS tile-weighted. The npz metadata must agree with the label.
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))                      # MorphoResidual_paper/
SIDECAR = os.path.join(PAPER, "review", "POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256")
OUT_DIR = os.path.join(PAPER, "review", "BRCA_TS_ARM_2026-10-08", "setlevel", "results")
SELF_REL = "review/BRCA_TS_ARM_2026-10-08/local_brca_setlevel.py"

# ------------------------------------------------------------------ the pre-specification constants
SETDIR_REL = "review/BRCA_TS_ARM_2026-10-08/setlevel/candidate_sets/"
PRIMARY_SET = SETDIR_REL + "pub_k3.tsv"            # published discovery (UNI): significant in >= 3 of 5 cohorts (777 genes)
SECONDARY_SETS = (SETDIR_REL + "pub_k4.tsv",       # the 170-gene core
                  SETDIR_REL + "pub_k2.tsv",       # >= 2 of 5
                  SETDIR_REL + "op_k2.tsv",        # operator batch-corrected, >= 2 of 5
                  SETDIR_REL + "phi_k3.tsv")       # Phikon-selected (the BRCA arm's encoder), >= 3 of 5
PRIMARY_ARM = "ts"
ARM_META = {"ts": ("ts", "slide"), "dx": ("dx", "slide"), "ts_tile": ("ts", "tile")}   # label -> (npz arm, npz pool)
EXPECT_B = 1000
MIN_GROUP = 10                                     # fewer selected or background genes in the matrix -> untestable, not run
MIN_COVERAGE = 0.5                                 # below this fraction of a set's selected genes in the matrix: PRIMARY -> FATAL,
                                                   # SECONDARY -> untestable (N5)
DISCOVERY_UNIVERSE_N = 8309                        # genes tested in all five discovery cohorts = rows of every frozen set file
EXPECT_DX_CASES = 102                              # the published DX arm
MIN_TS_CASES = 95                                  # TS arms with fewer cases go to the supplement only (prespec section 1)
GATE_EXACT_TOL = 1e-12                             # DX gate tiers (prespec section 1, review M3)
GATE_HOST_TOL = 1e-4
NEGFRAC_DECIMALS = 4                               # negative_frac_all compared to 0.01 percentage points (Methods: 99.11 vs 99.12%)
PUB_N_CASES, PUB_TESTED, PUB_SIGNIFICANT = 102, 10031, 0   # the published summary (prespec section 1, published-file check)
PUB_NEGFRAC_3DP = 0.991                            # negative_frac_all of the published summary must round to this
PUB_FILES = ("results/brca_protein_results.csv", "results/summary_brca.csv")   # path suffixes looked up in dx_arm_inputs_outputs.sha256
TOL = 1e-12
# ---------------------------------------------------------------------------------------------------


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fatal(msg):
    sys.exit(msg)


def check_sidecar():
    """Fail closed. Returns {relpath: sha256} of the frozen files."""
    if not os.path.exists(SIDECAR):
        fatal(f"REFUSED: sidecar {SIDECAR} not found -- the set-level pre-specification is not frozen.")
    listed = {}
    for line in open(SIDECAR, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        want, rel = line.split(None, 1)
        rel = rel.lstrip("*").strip()
        path = os.path.join(PAPER, rel)
        if not os.path.exists(path) or sha256(path) != want:
            fatal(f"REFUSED: {rel} missing or changed since freezing")
        listed[rel] = want
    need = [SELF_REL, PRIMARY_SET, *SECONDARY_SETS]
    for rel in need:
        if rel not in listed:
            fatal(f"REFUSED: {rel} is not listed in the sidecar (script, primary and secondary set files must all be frozen)")
    if not any(r.endswith(".md") for r in listed):
        fatal("REFUSED: the sidecar lists no .md pre-specification")
    print(f"sidecar verified: {len(listed)} files")
    return listed


def assemble(M, xs, sel_mask):
    """Identical to c1_run_test.assemble / c1_ucec_posthoc_sets.assemble."""
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
        fatal("FATAL: non-finite statistic")
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))


def holm(pvals):
    pvals = np.asarray(pvals, float)
    order = np.argsort(pvals)
    m, running, adj = len(pvals), 0.0, np.empty(len(pvals))
    for k, i in enumerate(order):
        running = max(running, min(1.0, (m - k) * pvals[i]))
        adj[i] = running
    return adj


def load_matrix(path, label, expect_b=EXPECT_B):
    z = np.load(path, allow_pickle=False)
    need = ["genes", "M", "B", "perm_seeds", "arm", "pool", "project", "n_pcs"]
    miss = [k for k in need if k not in z.files]
    if miss:
        fatal(f"FATAL: {path} lacks {miss} (not a cptac2_residual.py patch-v3 matrix)")
    if "smoke" in z.files and bool(z["smoke"]):
        fatal(f"FATAL: {path} is a smoke-test matrix")
    arm, pool = ARM_META[label]
    if str(z["arm"]) != arm or str(z["pool"]) != pool:
        fatal(f"FATAL: {path}: label '{label}' expects arm={arm}/pool={pool}, file says {z['arm']}/{z['pool']}")
    if str(z["project"]) != "TCGA-BRCA":
        fatal(f"FATAL: {path} is project {z['project']}")
    M = np.asarray(z["M"], float); genes = np.asarray(z["genes"]).astype(str)
    if int(z["B"]) != expect_b or M.shape[1] != expect_b + 1:
        fatal(f"FATAL: {path}: B={int(z['B'])}, M has {M.shape[1]} columns, expected B={expect_b}")
    if not np.array_equal(np.asarray(z["perm_seeds"]), np.arange(expect_b)):
        fatal(f"FATAL: {path}: permutation seeds are not 0..{expect_b - 1}")
    if int(z["n_pcs"]) != 20:
        fatal(f"FATAL: {path}: n_pcs != 20")
    if len(np.unique(genes)) != len(genes):
        fatal(f"FATAL: duplicated genes in {path}")
    if M.shape[0] != len(genes) or not np.isfinite(M).all():
        fatal(f"FATAL: {path}: shape/non-finite problem")
    cases = np.asarray(z["cases"]).astype(str) if "cases" in z.files else None
    return dict(path=path, sha256=sha256(path), genes=genes, M=M, n_cases=int(len(cases)) if cases is not None else -1,
                cases=cases, n=np.asarray(z["n"]) if "n" in z.files else None,
                fdr=np.asarray(z["fdr"], float) if "fdr" in z.files else None,
                n_slides=int(z["n_slides"]) if "n_slides" in z.files else -1)


def read_set_file(rel_or_abs):
    path = rel_or_abs if os.path.isabs(rel_or_abs) else os.path.join(PAPER, rel_or_abs)
    df = pd.read_csv(path, sep="\t")
    for c in ("gene", "selected", "score"):
        if c not in df.columns:
            fatal(f"FATAL: {path} lacks column {c}")
    if df["gene"].duplicated().any():
        fatal(f"FATAL: duplicated genes in {path}")
    if not np.isfinite(df["score"].values.astype(float)).all():
        fatal(f"FATAL: non-finite discovery score in {path}")
    if not set(df["selected"].unique()) <= {0, 1}:
        fatal(f"FATAL: 'selected' is not 0/1 in {path}")
    return path, df


def run_set(mat, sdf, primary=True):
    """Universe = the BRCA-tested genes (matrix rows) among the set file's genes (the 8,309 genes tested in all five discovery
    cohorts); background = the other genes of that universe (review M1)."""
    pos = pd.Series(np.arange(len(mat["genes"])), index=mat["genes"])
    s = sdf[sdf["gene"].isin(pos.index)]
    n_sel_total = int(sdf["selected"].sum())
    n_sel_in = int(s["selected"].sum())
    out = dict(n_set_file_selected=n_sel_total, n_set_file_universe=int(len(sdf)),
               n_selected=n_sel_in, n_background=int(len(s) - n_sel_in), universe=int(len(s)),
               coverage_selected=round(n_sel_in / max(n_sel_total, 1), 4), coverage_universe=round(len(s) / len(sdf), 4))
    if n_sel_total and n_sel_in / n_sel_total < MIN_COVERAGE:
        if primary:
            fatal(f"FATAL: only {n_sel_in}/{n_sel_total} selected genes of the PRIMARY set are in the matrix -- symbol mapping "
                  f"problem; stop and inspect")
        out["note"] = f"untestable (coverage): only {n_sel_in}/{n_sel_total} selected genes in the matrix (< {MIN_COVERAGE:.0%})"
        return out
    if n_sel_in < MIN_GROUP or (len(s) - n_sel_in) < MIN_GROUP:
        out["note"] = "untestable: too few selected or background genes in the matrix"
        return out
    M = mat["M"][pos[s["gene"].values].values]
    sel = s["selected"].values.astype(bool)
    xs = s["score"].values.astype(float)
    auc, rho, n1, n0 = assemble(M, xs, sel)
    a_obs, a_null, r_obs, r_null = float(auc[0]), auc[1:], float(rho[0]), rho[1:]
    sd = float(a_null.std())
    out.update(auc_obs=a_obs, auc_null_mean=float(a_null.mean()), auc_null_sd=sd, margin=float(a_obs - a_null.mean()),
               z_auc=float((a_obs - a_null.mean()) / sd) if sd > 0 else float("nan"), p_auc=emp_p(a_obs, a_null),
               rho_obs=r_obs, rho_null_mean=float(r_null.mean()), rho_null_sd=float(r_null.std()), p_rho=emp_p(r_obs, r_null),
               median_incr_selected=float(np.median(M[sel, 0])), median_incr_background=float(np.median(M[~sel, 0])),
               frac_negative_incr_universe=float((M[:, 0] < 0).mean()))
    return out


def compare_dx_published(mat, csv):
    pub = pd.read_csv(csv).set_index("gene")
    g = pd.Series(mat["M"][:, 0], index=mat["genes"])
    same_set = set(pub.index) == set(g.index)
    diff = float(np.max(np.abs(g.reindex(pub.index).values - pub["incremental_r2"].values))) if same_set else float("nan")
    return dict(csv=csv, n_published=int(len(pub)), n_matrix=int(len(g)), same_gene_set=bool(same_set), max_abs_diff=diff,
                verdict="MATCH" if same_set and diff <= TOL else "MISMATCH (different node/library or code path -- report, do not hide)")


def _summary_row(path):
    df = pd.read_csv(path)
    if len(df) != 1:
        fatal(f"FATAL: {path} is not a one-row summary_brca.csv")
    r = df.iloc[0]
    for c in ("n_cases", "tested", "significant", "negative_frac_all"):
        if c not in df.columns:
            fatal(f"FATAL: {path} lacks column {c}")
    return dict(n_cases=int(r["n_cases"]), tested=int(r["tested"]), significant=int(r["significant"]),
                negative_frac_all=float(r["negative_frac_all"]))


def _hash_lines(path):
    """{path-as-listed: sha256} from a `sha256sum` file (lines 'hash  path' or 'hash *path')."""
    out = {}
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\r\n")
        if not line.strip():
            continue
        parts = line.split(None, 1)
        if len(parts) != 2 or len(parts[0]) != 64:
            fatal(f"FATAL: {path}: not a sha256sum line: {line[:80]!r}")
        out[parts[1].lstrip("*").strip()] = parts[0].lower()
    return out


def published_file_check(pub_csv, pub_summary, ps, sha_file):
    """Cross-review B2. Returns (ok, reasons): the two published files hash to their dx_arm_inputs_outputs.sha256 lines and the
    published summary reads n_cases/tested/significant/negative_frac_all as published."""
    reasons = []
    listed = _hash_lines(sha_file)
    for suffix, f in zip(PUB_FILES, (pub_csv, pub_summary)):
        hits = [h for pth, h in listed.items() if pth.replace("\\", "/").endswith("/" + suffix) or pth.replace("\\", "/") == suffix]
        if len(hits) != 1:
            reasons.append(f"{suffix}: {len(hits)} line(s) in {os.path.basename(sha_file)}, expected exactly 1")
        elif sha256(f) != hits[0]:
            reasons.append(f"{suffix}: sha256 differs from the step-0 line in {os.path.basename(sha_file)}")
    if (ps["n_cases"], ps["tested"], ps["significant"]) != (PUB_N_CASES, PUB_TESTED, PUB_SIGNIFICANT):
        reasons.append(f"published summary reads n_cases/tested/significant {ps['n_cases']}/{ps['tested']}/{ps['significant']}, "
                       f"expected {PUB_N_CASES}/{PUB_TESTED}/{PUB_SIGNIFICANT}")
    if abs(round(ps["negative_frac_all"], 3) - PUB_NEGFRAC_3DP) > 1e-9:
        reasons.append(f"published negative_frac_all {ps['negative_frac_all']:.6f} does not round to {PUB_NEGFRAC_3DP}")
    return (not reasons), reasons


def dx_gate(mat, pub_csv, pub_summary, repro_summary, pub_sha256):
    """Three-tier DX reproduction gate (review M3): exact / host / fail. Published files are checked first (cross-review B2)."""
    c = compare_dx_published(mat, pub_csv)
    ps, rs = _summary_row(pub_summary), _summary_row(repro_summary)
    pub_ok, pub_reasons = published_file_check(pub_csv, pub_summary, ps, pub_sha256)
    m0 = mat["M"][:, 0]
    # the re-run summary must belong to this matrix (otherwise the wrong file was brought back)
    if rs["tested"] != len(m0) or abs(rs["negative_frac_all"] - float((m0 < 0).mean())) > TOL:
        fatal(f"FATAL: {repro_summary} does not belong to the dx matrix {mat['path']}")
    if mat.get("fdr") is not None and rs["significant"] != int(((mat["fdr"] < 0.05) & (m0 > 0)).sum()):
        fatal(f"FATAL: {repro_summary}: significant count disagrees with the dx matrix")
    if mat["n_cases"] >= 0 and rs["n_cases"] != mat["n_cases"]:
        fatal(f"FATAL: {repro_summary}: n_cases disagrees with the dx matrix")
    same_summary = (ps["n_cases"] == rs["n_cases"] and ps["tested"] == rs["tested"] and ps["significant"] == rs["significant"]
                    and round(ps["negative_frac_all"], NEGFRAC_DECIMALS) == round(rs["negative_frac_all"], NEGFRAC_DECIMALS))
    d = c["max_abs_diff"]
    if not pub_ok:
        tier = "fail"
    elif c["same_gene_set"] and d <= GATE_EXACT_TOL and same_summary:
        tier = "exact"
    elif c["same_gene_set"] and d <= GATE_HOST_TOL and same_summary:
        tier = "host"
    else:
        tier = "fail"
    return dict(tier=tier, max_abs_diff=d, same_gene_set=c["same_gene_set"], same_summary=bool(same_summary),
                published_files_ok=bool(pub_ok), published_files_problems=pub_reasons, published_sha256_file=pub_sha256,
                published_summary=ps, rerun_summary=rs, published_csv=pub_csv, published_summary_csv=pub_summary,
                rerun_summary_csv=repro_summary,
                label="DX re-run (reproduces published)" if tier != "fail" else "DX re-run (not identical to published)")


def check_cases(mats):
    """Review M5/M6: dx = 102 distinct cases; ts and ts_tile identical, a subset of dx; n_ts < 95 -> supplement only."""
    info = {}
    for lab, m in mats.items():
        if m["cases"] is None:
            fatal(f"FATAL: {m['path']} records no case ids")
        if len(set(m["cases"])) != len(m["cases"]):
            fatal(f"FATAL: {m['path']}: duplicated case ids")
    dxc = set(mats["dx"]["cases"])
    if len(dxc) != EXPECT_DX_CASES:
        fatal(f"FATAL: the dx matrix holds {len(dxc)} cases, the published arm {EXPECT_DX_CASES}")
    if "ts" in mats and "ts_tile" in mats and set(mats["ts"]["cases"]) != set(mats["ts_tile"]["cases"]):
        fatal("FATAL: ts and ts_tile matrices hold different case sets (same embeddings, so they must agree)")
    if "ts" in mats and "ts_tile" in mats and mats["ts"]["n_slides"] != mats["ts_tile"]["n_slides"]:
        fatal("FATAL: ts and ts_tile matrices were built from different numbers of slides")
    for lab in ("ts", "ts_tile"):
        if lab not in mats:
            continue
        c = set(mats[lab]["cases"])
        if not c <= dxc:
            fatal(f"FATAL: the {lab} matrix holds cases outside the published DX case set: {sorted(c - dxc)[:5]}")
        info[lab] = dict(n_cases=len(c), missing_cases=sorted(dxc - c), n_slides=mats[lab]["n_slides"],
                         status="main" if len(c) >= MIN_TS_CASES else f"supplement only (n < {MIN_TS_CASES})")
    info["dx"] = dict(n_cases=len(dxc), missing_cases=[], n_slides=mats["dx"]["n_slides"], status="main")
    return info


def real_run(matrices, dx_csv, out_dir, dx_pub_summary="", dx_repro_summary="", dx_pub_sha256=""):
    listed = check_sidecar()
    for fn in ("setlevel_result.csv", "setlevel_result.json"):
        if os.path.exists(os.path.join(out_dir, fn)):
            fatal(f"REFUSED: {os.path.join(out_dir, fn)} exists -- the set-level read-out is run once")
    if PRIMARY_ARM not in matrices:
        fatal(f"FATAL: the primary arm '{PRIMARY_ARM}' matrix was not supplied")
    if "dx" not in matrices:
        fatal("FATAL: the dx matrix was not supplied (needed for the DX segment of the Results sentence and the case check)")
    if not (dx_csv and dx_pub_summary and dx_repro_summary and dx_pub_sha256):
        fatal("FATAL: the DX gate needs --dx-published-csv, --dx-published-summary, --dx-repro-summary and --dx-published-sha256")
    missing_arms = [k for k in ARM_META if k not in matrices]
    if missing_arms:
        print(f"WARNING: pre-declared arm(s) not supplied: {missing_arms} -- report as a deviation")
    mats = {lab: load_matrix(p, lab) for lab, p in matrices.items()}
    cases = check_cases(mats)
    gate = dx_gate(mats["dx"], dx_csv, dx_pub_summary, dx_repro_summary, dx_pub_sha256)
    print(f"DX reproduction gate: {gate['tier']} (max|delta| {gate['max_abs_diff']:.3g}; summary identical: {gate['same_summary']}; "
          f"published files ok: {gate['published_files_ok']} {gate['published_files_problems']})")
    sets = {"PRIMARY": PRIMARY_SET, **{f"S{i + 1}": s for i, s in enumerate(SECONDARY_SETS)}}
    rows = []
    for lab, mat in mats.items():
        for role, rel in sets.items():
            path, sdf = read_set_file(rel)
            if sha256(path) != listed[rel]:
                fatal(f"REFUSED: {rel} changed after the sidecar check")
            if len(sdf) != DISCOVERY_UNIVERSE_N:
                fatal(f"FATAL: {rel} has {len(sdf)} genes, not the {DISCOVERY_UNIVERSE_N}-gene discovery universe")
            r = run_set(mat, sdf, primary=(role == "PRIMARY"))
            rows.append(dict(arm=lab, role=role, set=os.path.basename(rel)[:-4], matrix_sha256=mat["sha256"],
                             n_cases=mat["n_cases"], n_slides=mat["n_slides"], n_matrix_genes=int(len(mat["genes"])),
                             arm_status=cases[lab]["status"],
                             dx_gate=gate["tier"] if lab == "dx" else "",
                             arm_label=gate["label"] if lab == "dx" else lab, **r))
    out = pd.DataFrame(rows)
    # secondary family: Holm within arm over S1..Sm (one-sided AUC p)
    out["p_auc_holm_secondary"] = np.nan
    for lab in mats:
        idx = out.index[(out["arm"] == lab) & (out["role"] != "PRIMARY") & out["p_auc"].notna()] if "p_auc" in out else []
        if len(idx):
            out.loc[idx, "p_auc_holm_secondary"] = holm(out.loc[idx, "p_auc"].values)
    os.makedirs(out_dir, exist_ok=True)
    out.to_csv(os.path.join(out_dir, "setlevel_result.csv"), index=False, lineterminator="\n")
    meta = dict(script_sha256=sha256(__file__), sidecar=listed, primary_arm=PRIMARY_ARM, primary_set=PRIMARY_SET,
                secondary_sets=list(SECONDARY_SETS), matrices={k: dict(path=v["path"], sha256=v["sha256"]) for k, v in mats.items()},
                universe=f"BRCA-tested genes among the {DISCOVERY_UNIVERSE_N} genes tested in all five discovery cohorts",
                arms_not_supplied=missing_arms, cases=cases, dx_gate=gate)
    with open(os.path.join(out_dir, "setlevel_result.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(dict(meta=meta, rows=out.to_dict(orient="records")), fh, indent=1, default=str)
    for lab, ci in cases.items():
        print(f"arm {lab}: {ci['n_cases']} cases ({ci['status']}); missing vs dx: {ci['missing_cases']}")
    cols = [c for c in ["arm", "role", "set", "arm_status", "dx_gate", "n_selected", "n_background", "coverage_selected", "auc_obs",
                        "auc_null_mean", "auc_null_sd", "p_auc", "p_auc_holm_secondary", "rho_obs", "p_rho", "note"]
            if c in out.columns]
    pd.set_option("display.width", 220)
    print(out[cols].to_string(index=False))
    print("wrote", out_dir)
    return out


# =====================================================================================================
# synthetic self-test: statistic vs brute force, calibration, power, one-sidedness, fail-closed paths
# =====================================================================================================
def _synthetic_matrix(rng, G, B, shift_sel=None, sel=None, n_factors=1, load=0.7, ties=False):
    """Columns exchangeable (observed ~ null); gene-gene correlation via a shared factor per column."""
    f = rng.standard_normal((n_factors, B + 1))
    L = rng.standard_normal((G, n_factors)) * load
    M = L @ f + rng.standard_normal((G, B + 1))
    if ties:
        M = np.round(M, 1)
    if shift_sel is not None:
        M[sel, 0] += shift_sel
    return M


def _syn_cases(n, offset=0):
    return np.array([f"TCGA-SY-{i:04d}" for i in range(offset, offset + n)])


def _write_npz(path, genes, M, arm="ts", pool="slide", project="TCGA-BRCA", B=None, seeds=None, cases=None, fdr=None,
               n_slides=150):
    B = M.shape[1] - 1 if B is None else B
    extra = {} if fdr is None else dict(fdr=np.asarray(fdr, float))
    np.savez_compressed(path, genes=np.array(genes, dtype=str), M=M, n=np.full(len(genes), 100),
                        cases=_syn_cases(EXPECT_DX_CASES) if cases is None else np.asarray(cases),
                        perm_seeds=np.arange(B) if seeds is None else seeds, B=np.int64(B), n_pcs=np.int64(20),
                        arm=np.array(arm), pool=np.array(pool), project=np.array(project), smoke=np.bool_(False),
                        n_slides=np.int64(n_slides), **extra)


def selftest():
    from scipy.stats import mannwhitneyu, spearmanr
    rng = np.random.RandomState(20261008)
    # --- 1. statistic vs brute force (continuous and tied) ---
    for ties in (False, True):
        G, B = 240, 60
        sel = np.zeros(G, bool); sel[rng.choice(G, 55, replace=False)] = True
        M = _synthetic_matrix(rng, G, B, ties=ties)
        xs = rng.standard_normal(G) + 0.3 * sel
        if ties:
            xs = np.round(xs, 1)
        auc, rho, n1, n0 = assemble(M, xs, sel)
        for j in (0, 1, 7, B):
            u = mannwhitneyu(M[sel, j], M[~sel, j]).statistic
            assert abs(auc[j] - u / (n1 * n0)) < TOL, (ties, j, auc[j], u / (n1 * n0))
            assert abs(rho[j] - spearmanr(xs, M[:, j]).correlation) < TOL, (ties, j)
        p = emp_p(auc[0], auc[1:])
        assert p == (np.sum(auc[1:] >= auc[0]) + 1) / (len(auc)), p
    print("selftest 1 ok: AUC == Mann-Whitney U/(n1 n0), rho == scipy spearmanr (<=1e-12), with and without ties")
    # --- 2. calibration under the null with gene-gene correlation; null SD wider than the iid formula ---
    G, B, R = 300, 199, 400
    sel = np.zeros(G, bool); sel[:60] = True
    xs = rng.standard_normal(G)
    ps, prs, sds = [], [], []
    Lsel = np.zeros((G, 2)); Lsel[:60] = 0.4              # selected genes co-vary (a shared module), as co-regulated gene sets do
    for _ in range(R):
        M = _synthetic_matrix(rng, G, B, load=0.5, n_factors=2) + Lsel @ rng.standard_normal((2, B + 1))
        auc, rho, n1, n0 = assemble(M, xs, sel)
        ps.append(emp_p(auc[0], auc[1:])); prs.append(emp_p(rho[0], rho[1:])); sds.append(auc[1:].std())
    ps, prs = np.array(ps), np.array(prs)
    iid_sd = np.sqrt((n1 + n0 + 1) / (12.0 * n1 * n0))
    assert 0.02 <= (ps <= 0.05).mean() <= 0.085 and 0.35 < ps.mean() < 0.65, ((ps <= 0.05).mean(), ps.mean())
    assert 0.02 <= (prs <= 0.05).mean() <= 0.085, (prs <= 0.05).mean()
    assert np.mean(sds) > 1.5 * iid_sd, (np.mean(sds), iid_sd)       # correlated genes => iid SD would be anticonservative
    print(f"selftest 2 ok: type-I {np.mean(ps <= 0.05):.3f} (AUC) / {np.mean(prs <= 0.05):.3f} (rho) over {R} reps; "
          f"null SD {np.mean(sds):.4f} vs iid formula {iid_sd:.4f} (x{np.mean(sds) / iid_sd:.1f})")
    # --- 3. power and one-sidedness ---
    G, B = 400, 999
    sel = np.zeros(G, bool); sel[:80] = True
    M = _synthetic_matrix(rng, G, B, shift_sel=2.0, sel=sel)
    auc, rho, n1, n0 = assemble(M, rng.standard_normal(G), sel)
    assert emp_p(auc[0], auc[1:]) == 1.0 / (B + 1) and auc[0] > 0.8
    M = _synthetic_matrix(rng, G, B, shift_sel=-2.0, sel=sel)
    auc, rho, n1, n0 = assemble(M, rng.standard_normal(G), sel)
    assert emp_p(auc[0], auc[1:]) == 1.0 and auc[0] < 0.2           # one-sided upper tail: opposite sign is not significant
    print("selftest 3 ok: planted +shift -> p = 1/(B+1); planted -shift -> p = 1 (one-sided)")
    # --- 4. fail-closed paths + end-to-end on synthetic npz/set files in a temp PAPER ---
    tmp = tempfile.mkdtemp(prefix="brca_setlevel_selftest_")
    saved = {k: globals()[k] for k in ("PAPER", "SIDECAR", "PRIMARY_SET", "SECONDARY_SETS", "PRIMARY_ARM", "DISCOVERY_UNIVERSE_N",
                                       "PUB_TESTED", "PUB_SIGNIFICANT", "PUB_NEGFRAC_3DP")}
    try:
        globals().update(PAPER=tmp, SIDECAR=os.path.join(tmp, "review", "BRCA_TS_ARM_2026-10-08", "S.md.sha256"),
                         PRIMARY_SET="sets/primary.tsv", SECONDARY_SETS=("sets/sec1.tsv", "sets/sec2.tsv"),
                         DISCOVERY_UNIVERSE_N=1500)
        os.makedirs(os.path.join(tmp, "sets")); os.makedirs(os.path.dirname(SIDECAR))
        shutil.copy(__file__, os.path.join(tmp, "local_brca_setlevel.py"))
        globals()["SELF_REL"] = "local_brca_setlevel.py"
        genes = [f"G{i}" for i in range(1200)]
        # discovery universe has 1500 genes; 1200 are 'tested' in the synthetic BRCA matrix (coverage 0.8)
        ug = [f"G{i}" for i in range(1500)]
        zscore = rng.standard_normal(1500)

        def P(f):
            return os.path.join(tmp, f)

        def write_set(nm, a0, a1):
            sv = np.zeros(len(ug), int); sv[a0:a1] = 1
            pd.DataFrame({"gene": ug, "selected": sv, "score": zscore, "k": -1}).to_csv(
                P(f"sets/{nm}.tsv"), sep="\t", index=False, lineterminator="\n")
        for nm, (a0, a1) in (("primary", (0, 150)), ("sec1", (0, 60)), ("sec2", (1000, 1250))):
            write_set(nm, a0, a1)
        sel = np.zeros(1200, bool); sel[:150] = True
        Mts = _synthetic_matrix(rng, 1200, EXPECT_B, shift_sel=0.9, sel=sel)      # planted: primary set shifted in the observed column
        Mdx = _synthetic_matrix(rng, 1200, EXPECT_B)                              # null arm
        fdr_dx = np.ones(1200); fdr_dx[:3] = 0.01                                 # 3 FDR-passing rows (the sign decides 'significant')
        sig_dx = int(((fdr_dx < 0.05) & (Mdx[:, 0] > 0)).sum())
        negf = float((Mdx[:, 0] < 0).mean())
        # the published-summary constants are the synthetic arm's own values for the duration of the selftest
        globals().update(PUB_TESTED=1200, PUB_SIGNIFICANT=sig_dx, PUB_NEGFRAC_3DP=round(negf, 3))
        _write_npz(P("ts.npz"), genes, Mts, "ts", "slide")
        _write_npz(P("ts_tile.npz"), genes, Mts + 0.01 * rng.standard_normal(Mts.shape), "ts", "tile")
        _write_npz(P("dx.npz"), genes, Mdx, "dx", "slide", fdr=fdr_dx, n_slides=107)
        _write_npz(P("bad_b.npz"), genes, Mdx[:, :500], "dx", "slide")
        _write_npz(P("bad_arm.npz"), genes, Mdx, "dx", "tile")
        _write_npz(P("ts94.npz"), genes, Mts, "ts", "slide", cases=_syn_cases(94), n_slides=170)
        _write_npz(P("ts94_tile.npz"), genes, Mts, "ts", "tile", cases=_syn_cases(94), n_slides=170)
        _write_npz(P("ts94_tile_other.npz"), genes, Mts, "ts", "tile", cases=_syn_cases(94, offset=1), n_slides=170)
        _write_npz(P("ts_out.npz"), genes, Mts, "ts", "slide", cases=_syn_cases(102, offset=5))
        _write_npz(P("dx101.npz"), genes, Mdx, "dx", "slide", fdr=fdr_dx, cases=_syn_cases(101))
        matrices = {"ts": P("ts.npz"), "dx": P("dx.npz"), "ts_tile": P("ts_tile.npz")}
        # DX gate inputs: published per-gene csv + published / re-run summaries

        def write_summary(fn, n_cases=EXPECT_DX_CASES, tested=1200, significant=sig_dx, negative_frac_all=negf):
            pd.DataFrame([dict(project="TCGA-BRCA", n_cases=n_cases, tested=tested, significant=significant, pct_sig=0.0,
                               median_incr_sig=np.nan, negative_frac_all=negative_frac_all)]).to_csv(
                P(fn), index=False, lineterminator="\n")
        write_summary("pub_summary.csv"); write_summary("rerun_summary.csv")
        write_summary("pub_summary_sigdiff.csv", significant=sig_dx + 1)
        write_summary("rerun_summary_wrong.csv", tested=1199)
        pd.DataFrame({"gene": genes, "incremental_r2": Mdx[:, 0]}).to_csv(P("pub.csv"), index=False, lineterminator="\n")
        pd.DataFrame({"gene": genes, "incremental_r2": Mdx[:, 0] + rng.uniform(-5e-5, 5e-5, 1200)}).to_csv(
            P("pub_host.csv"), index=False, lineterminator="\n")
        pd.DataFrame({"gene": genes, "incremental_r2": Mdx[:, 0] + rng.uniform(-1e-3, 1e-3, 1200)}).to_csv(
            P("pub_fail.csv"), index=False, lineterminator="\n")
        def write_shafile(fn, csv, summ, lines_for=("results/brca_protein_results.csv", "results/summary_brca.csv")):
            """Mimic RUNBOOK step 0: sha256sum of the published inputs/results (server paths), here for the given files."""
            with open(P(fn), "w", encoding="utf-8", newline="\n") as fh:
                fh.write(f"{sha256(P('pub.csv'))}  /srv/cptac2_brca/protein_brca.csv\n")        # an unrelated line, as on the server
                for suffix, f in zip(("results/brca_protein_results.csv", "results/summary_brca.csv"), (csv, summ)):
                    if suffix in lines_for:
                        fh.write(f"{sha256(f)}  /srv/cptac2_brca/{suffix}\n")
        write_shafile("pub.sha", P("pub.csv"), P("pub_summary.csv"))
        write_shafile("pub_host.sha", P("pub_host.csv"), P("pub_summary.csv"))
        write_shafile("pub_fail.sha", P("pub_fail.csv"), P("pub_summary.csv"))
        write_shafile("pub_sigdiff.sha", P("pub.csv"), P("pub_summary_sigdiff.csv"))
        write_shafile("pub_nosummaryline.sha", P("pub.csv"), P("pub_summary.csv"), lines_for=("results/brca_protein_results.csv",))
        shutil.copy(P("pub.sha"), P("pub_stale.sha"))                                      # hashes of pub.csv, used with a changed file
        G3 = (P("pub.csv"), P("pub_summary.csv"), P("rerun_summary.csv"), P("pub.sha"))
        outn = [0]

        def run(mx=None, g=G3):
            outn[0] += 1                                  # a fresh output folder per call: a real run refuses to overwrite
            return real_run(matrices if mx is None else mx, g[0], P(f"o{outn[0]}"), g[1], g[2], g[3])

        def expect_exit(fn, needle):
            try:
                fn()
            except SystemExit as e:
                assert needle in str(e.code), (needle, e.code)
                return
            raise AssertionError(f"expected exit containing {needle!r}")

        expect_exit(run, "REFUSED: sidecar")                                                            # (a) no sidecar
        files = ["local_brca_setlevel.py", "sets/primary.tsv", "sets/sec1.tsv", "sets/sec2.tsv"]
        with open(P("prespec.md"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("prespec\n")

        def write_sidecar(fl):
            with open(SIDECAR, "w", encoding="utf-8", newline="\n") as fh:
                for f in fl:
                    fh.write(f"{sha256(os.path.join(tmp, f))} *{f}\n")
        write_sidecar(files[:3] + ["prespec.md"])
        expect_exit(run, "sec2.tsv is not listed")                                                      # (b) secondary not frozen
        write_sidecar(files)
        expect_exit(run, "no .md pre-specification")                                                    # (c) no prespec listed
        write_sidecar(files + ["prespec.md"])
        with open(P("sets/sec1.tsv"), "a", encoding="utf-8", newline="\n") as fh:                       # (d) tampered after freezing
            fh.write("GX\t0\t0.0\t-1\n")
        expect_exit(run, "changed since freezing")
        write_sidecar(files + ["prespec.md"])
        expect_exit(run, "not the 1500-gene discovery universe")                                        # set file != universe
        write_set("sec1", 0, 60)
        write_sidecar(files + ["prespec.md"])
        expect_exit(lambda: run({"dx": matrices["dx"]}), "primary arm 'ts'")
        expect_exit(lambda: run({"ts": matrices["ts"]}), "dx matrix was not supplied")
        expect_exit(lambda: run(g=(P("pub.csv"), "", P("rerun_summary.csv"), P("pub.sha"))), "DX gate needs")
        expect_exit(lambda: run(g=(P("pub.csv"), P("pub_summary.csv"), P("rerun_summary.csv"), "")), "DX gate needs")
        expect_exit(lambda: load_matrix(P("bad_b.npz"), "dx"), "B=")
        expect_exit(lambda: load_matrix(P("bad_arm.npz"), "dx"), "expects arm=dx/pool=slide")
        # case-set rules (review M5/M6)
        expect_exit(lambda: run({"ts": P("ts.npz"), "dx": P("dx101.npz")}), "dx matrix holds 101 cases")
        expect_exit(lambda: run({"ts": P("ts94.npz"), "ts_tile": P("ts94_tile_other.npz"), "dx": P("dx.npz")}),
                    "different case sets")
        expect_exit(lambda: run({"ts": P("ts94.npz"), "ts_tile": P("ts_tile.npz"), "dx": P("dx.npz")}), "different case sets")
        expect_exit(lambda: run({"ts": P("ts_out.npz"), "dx": P("dx.npz")}), "outside the published DX case set")
        o94 = run({"ts": P("ts94.npz"), "ts_tile": P("ts94_tile.npz"), "dx": P("dx.npz")}).set_index(["arm", "role"])
        assert o94.loc[("ts", "PRIMARY"), "arm_status"] == f"supplement only (n < {MIN_TS_CASES})"
        assert o94.loc[("ts_tile", "PRIMARY"), "arm_status"].startswith("supplement only")
        assert o94.loc[("dx", "PRIMARY"), "arm_status"] == "main" and o94.loc[("ts", "PRIMARY"), "n_cases"] == 94
        # DX gate tiers (review M3)
        expect_exit(lambda: run(g=(P("pub.csv"), P("pub_summary.csv"), P("rerun_summary_wrong.csv"), P("pub.sha"))), "does not belong")
        # published-file check (cross-review B2): a published file that does not hash to its step-0 line, a missing line, or a
        # published summary that does not read as published, makes the gate fail even if the re-run matches the (wrong) file
        pd.DataFrame({"gene": genes, "incremental_r2": Mdx[:, 0]}).to_csv(P("pub_edited.csv"), index=False, lineterminator="\n",
                                                                         float_format="%.15g")     # same values, other bytes
        assert sha256(P("pub_edited.csv")) != sha256(P("pub.csv"))
        write_summary("pub_summary_ntested.csv", tested=1199)
        write_summary("rerun_summary_ntested.csv", tested=1200)
        for g, tier, why in (((P("pub.csv"), P("pub_summary.csv"), P("rerun_summary.csv"), P("pub.sha")), "exact", ""),
                             ((P("pub_host.csv"), P("pub_summary.csv"), P("rerun_summary.csv"), P("pub_host.sha")), "host", ""),
                             ((P("pub_fail.csv"), P("pub_summary.csv"), P("rerun_summary.csv"), P("pub_fail.sha")), "fail", ""),
                             ((P("pub.csv"), P("pub_summary_sigdiff.csv"), P("rerun_summary.csv"), P("pub_sigdiff.sha")), "fail", "significant"),
                             ((P("pub_edited.csv"), P("pub_summary.csv"), P("rerun_summary.csv"), P("pub_stale.sha")), "fail", "sha256 differs"),
                             ((P("pub.csv"), P("pub_summary.csv"), P("rerun_summary.csv"), P("pub_nosummaryline.sha")), "fail", "0 line(s)"),
                             ((P("pub.csv"), P("pub_summary_ntested.csv"), P("rerun_summary.csv"), P("pub.sha")), "fail", "sha256 differs")):
            og = run(g=g)
            dxr = og[og.arm == "dx"]
            assert (dxr["dx_gate"] == tier).all(), (tier, dxr["dx_gate"].tolist())
            assert (og.loc[og.arm != "dx", "dx_gate"] == "").all()
            want_label = "DX re-run (not identical to published)" if tier == "fail" else "DX re-run (reproduces published)"
            assert (dxr["arm_label"] == want_label).all()
            with open(P(f"o{outn[0]}/setlevel_result.json"), encoding="utf-8") as fh:
                gj = json.load(fh)["meta"]["dx_gate"]
            assert gj["tier"] == tier and gj["published_files_ok"] == (tier != "fail" or why == ""), (tier, why, gj["published_files_problems"])
            if why:
                assert any(why in r for r in gj["published_files_problems"]) or why == "significant", (why, gj["published_files_problems"])
        # a real run is made once: the output folder of a finished run is not overwritten
        expect_exit(lambda: real_run(matrices, *(G3[0], P(f"o{outn[0]}"), G3[1], G3[2], G3[3])), "run once")
        out = run()                                                                                      # (e) valid run
        r = out.set_index(["arm", "role"])
        assert r.loc[("ts", "PRIMARY"), "n_selected"] == 150 and r.loc[("ts", "PRIMARY"), "coverage_selected"] == 1.0
        assert r.loc[("ts", "S2"), "n_selected"] == 200 and r.loc[("ts", "S2"), "coverage_selected"] == 0.8     # symbols 1200-1249 absent
        assert r.loc[("ts", "PRIMARY"), "p_auc"] == 1.0 / (EXPECT_B + 1) and r.loc[("ts", "PRIMARY"), "auc_obs"] > 0.6
        assert r.loc[("dx", "PRIMARY"), "p_auc"] > 0.02                                                       # null arm: not planted
        assert np.isnan(r.loc[("ts", "PRIMARY"), "p_auc_holm_secondary"])
        assert (out.loc[out.arm == "dx", "dx_gate"] == "exact").all() and (out["arm_status"] == "main").all()
        # the full pipeline equals a direct computation on the matrix
        sdf = pd.read_csv(P("sets/primary.tsv"), sep="\t"); sdf = sdf[sdf.gene.isin(genes)]
        sel2 = sdf.selected.values.astype(bool)
        auc, rho, n1, n0 = assemble(Mts, sdf.score.values, sel2)
        assert abs(r.loc[("ts", "PRIMARY"), "auc_obs"] - auc[0]) < TOL and abs(r.loc[("ts", "PRIMARY"), "rho_obs"] - rho[0]) < TOL
        # coverage guard: a SECONDARY set below 50% -> untestable (N5); the PRIMARY set below 50% -> FATAL
        write_set("sec1", 1250, 1500)
        write_sidecar(files + ["prespec.md"])
        oc = run().set_index(["arm", "role"])
        assert str(oc.loc[("ts", "S1"), "note"]).startswith("untestable (coverage)") and np.isnan(oc.loc[("ts", "S1"), "p_auc"])
        assert oc.loc[("ts", "PRIMARY"), "p_auc"] == r.loc[("ts", "PRIMARY"), "p_auc"]               # primary unaffected
        assert np.isnan(oc.loc[("ts", "S1"), "p_auc_holm_secondary"]) and np.isfinite(oc.loc[("ts", "S2"), "p_auc_holm_secondary"])
        assert oc.loc[("ts", "S2"), "p_auc_holm_secondary"] == oc.loc[("ts", "S2"), "p_auc"]          # Holm over the 1 testable set
        write_set("sec1", 0, 60); write_set("primary", 1250, 1500)
        write_sidecar(files + ["prespec.md"])
        expect_exit(run, "PRIMARY set are in the matrix")
        # DX reproduction comparator (also used by test_synthetic_pipeline.py)
        assert compare_dx_published(load_matrix(matrices["dx"], "dx"), P("pub.csv"))["verdict"] == "MATCH"
        print("selftest 4 ok: REFUSED without sidecar / unfrozen secondary / no prespec / tampered file; set file != universe, "
              "missing ts or dx arm, missing gate inputs, wrong B, wrong arm label -> FATAL; case rules: dx != 102, ts vs ts_tile "
              "mismatch, ts case outside dx -> FATAL, n_ts = 94 -> supplement only; DX gate exact / host / fail (incl. summary "
              "mismatch) written into the dx rows and the json; published-file check (sha256 line, missing line, summary constants) -> fail; foreign re-run summary -> FATAL; existing result -> REFUSED (run once); valid run reproduces the direct "
              "computation; planted arm detected, null arm not; secondary < 50% coverage -> untestable, primary -> FATAL")
    finally:
        globals().update(saved)
        globals()["SELF_REL"] = "review/BRCA_TS_ARM_2026-10-08/local_brca_setlevel.py"
        shutil.rmtree(tmp, ignore_errors=True)
    print("SELFTEST PASSED (synthetic data only)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--matrix", action="append", default=[], help="label=path/to/brca_incr_matrix*.npz (labels: ts, dx, ts_tile)")
    ap.add_argument("--dx-published-csv", default="", help="published results/brca_protein_results.csv (DX gate)")
    ap.add_argument("--dx-published-summary", default="", help="published results/summary_brca.csv (DX gate)")
    ap.add_argument("--dx-repro-summary", default="", help="re-run results_dx_repro/summary_brca.csv (DX gate)")
    ap.add_argument("--dx-published-sha256", default="", help="dx_arm_inputs_outputs.sha256 from RUNBOOK step 0 (published-file check)")
    ap.add_argument("--out-dir", default=OUT_DIR)
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    if not a.matrix:
        fatal("no --matrix given (see --help)")
    matrices = {}
    for m in a.matrix:
        lab, _, p = m.partition("=")
        if lab not in ARM_META or not p:
            fatal(f"bad --matrix {m!r}; use label=path with label in {sorted(ARM_META)}")
        matrices[lab] = p
    real_run(matrices, a.dx_published_csv, a.out_dir, a.dx_published_summary, a.dx_repro_summary, a.dx_published_sha256)


if __name__ == "__main__":
    main()
