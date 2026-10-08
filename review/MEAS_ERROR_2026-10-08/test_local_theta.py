#!/usr/bin/env python3
"""
End-to-end test of local_theta.py / theta_kernel.py on SYNTHETIC cohorts only (generated here with
theory/sim_signtest.py's simulate_cohort and small purpose-built generators; no project data).
Revised 2026-10-08 for the FINAL reading rules of theory/DERIVATION.md (section 9).  Checks

  K1  batched/weighted ridge == the reference oof_ridge (permutation draw, bootstrap draw == duplicated rows)
  K2  local statistics == theory/sim_signtest.gene_stats (independently written) on complete genes
  K3  within-stratum permutations == residual_analysis_sitepack.py's construction (re-implemented here)
  K4  one-way ICC estimator recovers a known ICC
  S1  TESTED statistic = sign(rho) corr(rp_hat, M) = kernel Nc/|rho| for the observed value AND for permuted draws
      (independent re-computation with the reference ridge); the set row sums exactly this
  S2  kappa(s, lambda_lo) table spot checks incl. band boundaries
  S3  eligibility / strata: base sets contain only complete genes (n_fit = n); missingness strata have the right
      n_fit range; strata < 15 genes dropped; mnar_exposed flag; outputs written, reading.csv carries numbers and
      no verdict column
  V1  variant path 'var:linear' == kernel statistics; spline theta ~ linear theta when the truth is linear;
      RNA-PC / spline variants run through local_theta.main
  V2  shared RNA artefact read by morphology: primary z_N << 0, the 20-RNA-PC baseline removes it (|z_N| small)
  V3  purity channel: primary z_N >> 0, the purity-residualised variant removes it
  D1  nested-increment draws (phi_op machinery) == published_null_increments; observed == published_increment
  F1  phi_op: ~0 without batch structure, >> 0 when W encodes operator and P carries operator means
  M1  MNAR (collider) toy: truncation of P at the bottom 10% pushes z_N of the missingness stratum below the
      complete-gene stratum (numbers printed)
  P1  population values of DERIVATION.md section 2.1 at large n (theta_pc = f_X lam/(1-lam) | 0 | -1, mixed =
      f_X lam/(1-lam)) through the full local_theta pipeline (complete genes, and the MCAR-missing stratum)
  P2  n = 100, K = 20 set-level tests of DERIVATION.md section 7 under the operator-stratified PRIMARY null and the
      unrestricted sensitivity: N_c,set z >> 0 under H_noise, ~ 0 under H_post, << 0 under H_err; D_set z >> 0;
      diagnostic Wm z large only under H_noise
  G1  fail-closed gate: no sidecar -> REFUSED; wrong hash -> REFUSED; script / inputs / side file not listed ->
      REFUSED; --synthetic on a file not flagged synthetic -> REFUSED; --repro-only works; a valid sidecar runs
      the same numbers as the synthetic mode; spawn-multiprocessing == serial
  R1  repro check: matches on exported-like inputs, refuses when a published value is perturbed by 1e-6

Usage:  python test_local_theta.py <scratch dir> [--quick]
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "theory"))
import local_theta as lt   # noqa: E402
import theta_kernel as tk  # noqa: E402
import sim_signtest as sim  # noqa: E402

tmp = Path(sys.argv[1]).resolve()
QUICK = "--quick" in sys.argv
tmp.mkdir(parents=True, exist_ok=True)
FAIL = []
FAST = ["--variants", "none", "--no-diag", "--no-phi-delta"]


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL.append(name)


def assemble(path, W, M, P, rng, synthetic=True, with_pub=True, op_labels=None, rna_pcs=None):
    """Write an inputs_<c>.npz-like file from arrays (W (n,K), M (n,G), P (n,G) with NaN = missing protein)."""
    n, G = M.shape
    lab_op = op_labels if op_labels is not None else np.array([f"op{x}" for x in rng.integers(0, 5, size=n)])
    lab_plex = np.array([f"px{x}" for x in rng.integers(0, 8, size=n)])
    patients = np.array([f"S{i:04d}" for i in range(n)])
    genes = np.array([f"G{i:04d}" for i in range(G)])
    cnt = rng.poisson(30 * np.exp(rng.normal(size=(1, G)) * 1.2 + 3), size=(n, G)).astype(float)
    z = dict(schema_version=np.int64(2), is_synthetic=np.bool_(synthetic), cohort=np.array("syn"),
             patients=patients, genes=genes, mrna_log2tpm1=M, protein=P, counts=cnt, tpm=cnt / 10.0,
             pcs=W, wsi_pooled=W, pca_evr=np.ones(W.shape[1]) / W.shape[1],
             rna_pcs20=(rna_pcs if rna_pcs is not None else rng.normal(size=(n, 20))),
             n_rna_files=np.ones(n), icc_log2tpm1=np.full(G, np.nan), ridge_alpha=np.float64(1.0),
             cv_folds=np.int64(5), ra_seed=np.int64(0), ra_n_perm=np.int64(1000),
             label_operator=np.asarray(lab_op).astype(str), label_plex=lab_plex,
             label_scanner=np.array(["sc1"] * n), label_status=np.array("{}"))
    if with_pub:
        pub = np.array([tk.published_increment(M[:, g], P[:, g], W) for g in range(G)])
        z.update(pub_n=pub[:, 0], pub_r2_rna=pub[:, 1], pub_r2_both=pub[:, 2], pub_incremental_r2=pub[:, 3],
                 pub_pval=np.full(G, np.nan), pub_fdr=np.where(pub[:, 3] > 0.1, 0.01, 0.5))
    np.savez_compressed(path, **z)


def make_npz(path, hyp, n, G, lam, R2, delta, fX=0.5, seed=0, miss=0.12, synthetic=True, with_pub=True, mnar_q=0.0):
    """simulate_cohort data; the first half of the genes complete, the second half with MCAR missing protein
    (fraction `miss`) or, if mnar_q > 0, with the bottom mnar_q of P truncated (MNAR)."""
    rng = np.random.default_rng(seed)
    W, M, P, pop = sim.simulate_cohort(rng, n, 20, G, lam, R2, -0.5, hyp, delta, fX)
    P = P.copy()
    half = G // 2
    if mnar_q > 0:
        for g in range(half, G):
            P[P[:, g] <= np.quantile(P[:, g], mnar_q), g] = np.nan
    else:
        mask = rng.random((n, G)) < miss
        mask[:, :half] = False
        P[mask] = np.nan
    P = P * 3.0 + 1.0                                     # arbitrary scale/location
    M = M * 1.7 + 4.0
    assemble(path, W, M, P, rng, synthetic, with_pub)
    return pop, half


def scenario(path, kind, n=300, G=60, seed=5, op_labels=None):
    """Purpose-built generators: 'rna' (shared RNA artefact f in M read by morphology), 'purity' (f in both M and
    P), 'post' + operator means (phi_op), 'plain' (no structure).  Returns (purity vector or None)."""
    rng = np.random.default_rng(seed)
    W = rng.normal(size=(n, 20))
    zdir = np.zeros(20)
    zdir[:5] = 1 / np.sqrt(5)
    f = np.sqrt(0.7) * (W @ zdir) + np.sqrt(0.3) * rng.normal(size=n)
    l = rng.uniform(0.6, 1.4, G)
    X = rng.normal(size=(n, G))
    pur = None
    rna = None
    if kind == "rna":
        M = X + 1.0 * l * f[:, None] + 0.5 * rng.normal(size=(n, G))
        P = 0.8 * X + 0.6 * rng.normal(size=(n, G))
        from sklearn.decomposition import PCA
        Z = (M - M.mean(0)) / M.std(0)
        rna = PCA(n_components=20, svd_solver="full").fit_transform(Z)
    elif kind == "purity":
        M = X + 1.0 * l * f[:, None] + 0.5 * rng.normal(size=(n, G))
        P = 0.8 * X + 1.0 * l * f[:, None] + 0.6 * rng.normal(size=(n, G))
        pur = f + 0.1 * rng.normal(size=n)
    elif kind == "opmeans":                                # W encodes operator; P has operator means; no within-op signal
        lev = np.array(sorted(set(op_labels)))
        code = np.array([int(np.where(lev == o)[0][0]) for o in op_labels])
        mu = np.array([-1.5, -0.5, 0.5, 1.5])[code % 4]
        W[:, 0] = mu + 0.3 * rng.normal(size=n)
        M = X + 0.5 * rng.normal(size=(n, G))
        P = 0.8 * X + 0.8 * mu[:, None] * rng.uniform(0.6, 1.4, G)[None, :] + 0.6 * rng.normal(size=(n, G))
    elif kind == "plain":
        M = X + 0.5 * rng.normal(size=(n, G))
        P = 0.8 * X + 0.6 * rng.normal(size=(n, G)) + 0.3 * (W[:, [1]] * rng.normal(size=(1, G)))
    else:
        raise ValueError(kind)
    assemble(path, W, M, P, rng, with_pub=False, op_labels=op_labels, rna_pcs=rna)
    return pur


def row_of(sets, name, arm=None):
    s = sets[sets["set"] == name]
    if arm is not None:
        s = s[s["arm"] == arm]
    return s.iloc[0]


def run():
    # -------------------------------------------------------------------------------- K1
    rs = np.random.RandomState(3)
    n = 97
    Wf = rs.randn(n, 20)
    Mf = rs.randn(n)
    Pf = 0.4 * Mf + 0.5 * Wf[:, 0] + rs.randn(n)
    v = rs.rand(n) > 0.1
    nv = int(v.sum())
    fid = tk.fold_ids_for(nv)
    lit = tk.gene_stats_literal(Wf[v], Mf[v], Pf[v], secondary=True)
    bat = tk.gene_stats_batch(Wf[v][None], Mf[v], Pf[v], fid, np.ones((1, nv)), secondary=True)
    check("K1a batched == reference (observed)", max(abs(lit[k] - bat[k][0]) for k in lit) < 1e-12)
    perms = tk.build_perms(n, 30)
    bb = tk.gene_stats_batch(Wf[perms][:, v], Mf[v], Pf[v], fid, np.ones((1, nv)))
    worst = max(abs(tk.gene_stats_literal(Wf[perms[i]][v], Mf[v], Pf[v])[k] - bb[k][i]) for i in (0, 11, 29) for k in ("N", "Nc", "D", "S_cross"))
    check("K1b batched == reference (permuted W)", worst < 1e-12, f"{worst:.1e}")
    cnt = tk.boot_counts(n, 4)[:, v]
    bb = tk.gene_stats_batch(Wf[v][None], Mf[v], Pf[v], fid, cnt)
    worst = 0
    for i in range(4):
        idx = np.repeat(np.arange(nv), cnt[i].astype(int))
        fd = fid[idx]
        folds = [(np.where(fd != k)[0], np.where(fd == k)[0]) for k in range(5)]
        l = tk.gene_stats_literal(Wf[v][idx], Mf[v][idx], Pf[v][idx], folds=folds)
        worst = max(worst, max(abs(l[k] - bb[k][i]) for k in ("N", "Nc", "D", "S_cross")))
    check("K1c bootstrap weights == duplicated rows (fold inherited)", worst < 1e-12, f"{worst:.1e}")

    # -------------------------------------------------------------------------------- K3
    lab = np.array(["a"] * 30 + ["b"] * 40 + ["c"] * 1 + ["d"] * 26)
    mine, info = tk.build_perms_stratified(lab, 5, seed=0)
    levels = sorted(set(lab))
    RNG = np.random.RandomState(0)
    ref = []
    for _ in range(5):                                      # residual_analysis_sitepack.py main(), verbatim logic
        p = np.arange(len(lab))
        for lv in levels:
            m = np.where(lab == lv)[0]
            if len(m) > 1:
                p[m] = RNG.permutation(m)
        ref.append(p)
    check("K3 stratified permutations == sitepack construction", np.array_equal(mine, np.array(ref)) and info["n_fixed_singletons"] == 1)
    check("K3b unstratified perms == RA sequence (RandomState(0), sequential draws)",
          np.array_equal(tk.build_perms(50, 3)[0], np.random.RandomState(0).permutation(50)) and
          np.array_equal(tk.build_perms(50, 3)[1], (lambda r: (r.permutation(50), r.permutation(50))[1])(np.random.RandomState(0))))
    check("K3c default seed is fresh (not RandomState(0) = the selection draws)", lt.SEED_DEFAULT != 0)

    # -------------------------------------------------------------------------------- K4
    rng = np.random.default_rng(1)
    icc_true = 0.8
    groups = []
    for _ in range(60):
        k = int(rng.integers(2, 4))
        b = rng.normal(scale=np.sqrt(icc_true), size=(1, 30))
        groups.append(b + rng.normal(scale=np.sqrt(1 - icc_true), size=(k, 30)))
    icc, _, a, kbar = tk.icc_oneway(groups)
    check("K4 ICC(1) recovers 0.8", abs(icc.mean() - 0.8) < 0.04, f"mean {icc.mean():.3f}")

    # -------------------------------------------------------------------------------- K2 + S1
    p_syn = tmp / "syn_cross.npz"
    make_npz(p_syn, "mix", 100, 40, 0.7, 0.3, 0.2, fX=0.3, seed=11)
    d = lt.load_inputs(str(p_syn))
    keys = ["N", "Nc", "D", "S_cross"]
    obs, worst = lt.observed(d, keys, False)
    G = d["protein"].shape[1]
    comp = [g for g in range(G) if not np.isnan(d["protein"][:, g]).any()]
    Mz = (d["mrna_log2tpm1"][:, comp] - d["mrna_log2tpm1"][:, comp].mean(0)) / d["mrna_log2tpm1"][:, comp].std(0)
    Pz = (d["protein"][:, comp] - d["protein"][:, comp].mean(0)) / d["protein"][:, comp].std(0)
    simst = sim.gene_stats(d["pcs"], Mz, Pz, tk.folds_for(100), 1.0)
    w = max(np.max(np.abs(simst["N"] - obs["N"][comp])), np.max(np.abs(simst["N_corr"] - obs["Nc"][comp])),
            np.max(np.abs(simst["D"] - obs["D"][comp])), np.max(np.abs(simst["S_cross"] - obs["S_cross"][comp])))
    check("K2 local stats == sim_signtest.gene_stats (complete genes, same folds)", w < 1e-10, f"{w:.1e}")

    # S1: independent re-computation of the tested statistic
    w1 = 0.0
    for g in comp[:12]:
        gg = tk.prep_gene(d["mrna_log2tpm1"][:, g], d["protein"][:, g])
        rph = tk.oof_ridge(d["pcs"], gg["rp"], tk.folds_for(100))
        exp = np.sign(gg["rho"]) * tk._corr(rph, gg["Mz"])
        w1 = max(w1, abs(exp - obs["Ncs"][g]), abs(obs["Ncs"][g] - obs["Nc"][g] / abs(obs["rho"][g])))
    check("S1a observed N_c = sign(rho) corr(rp_hat, M) = Nc/|rho| (independent reference ridge)", w1 < 1e-10, f"{w1:.1e}")
    pm = tk.build_perms(100, 3, "ra", 5)
    nl = lt.run_draws(d, d["pcs"][pm], None, ["N", "Nc", "D"], False, 1, kind="stats", genes=comp[:6])
    nl["rho"] = obs["rho"]
    lt.add_sign_weighted(nl)
    w1 = 0.0
    for g in comp[:6]:
        gg = tk.prep_gene(d["mrna_log2tpm1"][:, g], d["protein"][:, g])
        for b in range(3):
            rph = tk.oof_ridge(d["pcs"][pm[b]], gg["rp"], tk.folds_for(100))
            w1 = max(w1, abs(np.sign(gg["rho"]) * tk._corr(rph, gg["Mz"]) - nl["Ncs"][b, g]))
    check("S1b null draws: N_c = sign(rho) corr(rp_hat(permuted W), M) (rho is a function of (M,P) only)", w1 < 1e-10, f"{w1:.1e}")
    # the set sum is exactly the sum of the sign-weighted gene statistics, and differs from the rho-weighted one
    out1 = lt.main(["--inputs", str(p_syn), "--synthetic", "--strata", "none", "--B", "10", "--boot", "0", "--workers", "1",
                    "--out", str(tmp / "out_S1")] + FAST)
    r_all = row_of(out1, "all")
    okc = np.array(comp)[np.isfinite(obs["N"][comp])]
    check("S1c set row N_c,set = sum_g sign(rho) corr(rp_hat, M) over the complete genes",
          abs(r_all["Nc_set"] - obs["Ncs"][okc].sum()) < 1e-9 and r_all["n_genes"] == len(okc), f"{r_all['Nc_set']:.4f} n={r_all['n_genes']}")
    check("S1d rho-weighted form kept only as diagnostic (Nc_rhow_*)", "Nc_rhow_z" in out1.columns and "Nc_rhow_p2" not in out1.columns)

    # -------------------------------------------------------------------------------- R1
    rep = lt.repro_published(d)
    check("R1a repro matches exported-like inputs", rep["ok"], str(rep.get("worst")))
    d2 = dict(d)
    d2["pub_incremental_r2"] = d["pub_incremental_r2"].copy()
    d2["pub_incremental_r2"][3] += 1e-6
    check("R1b repro refuses a 1e-6 perturbation", not lt.repro_published(d2)["ok"])

    # -------------------------------------------------------------------------------- S2 kappa table
    kt = [((0.60, 0.50), 0.65), ((0.50, 0.50), 0.65), ((0.25, 0.60), 0.45), ((0.25, 0.70), 0.35), ((0.20, 0.55), 0.55),
          ((0.15, 0.50), 0.45), ((0.10, 0.65), 0.40), ((0.07, 0.80), 0.25), ((0.03, 0.70), 0.25), ((0.03, 0.50), 0.35),
          ((0.019, 0.50), 0.30), ((0.0, 0.90), 0.25)]
    bad = [(a, lt.kappa_of(*a), e) for a, e in kt if abs(lt.kappa_of(*a) - e) > 1e-12]
    check("S2 kappa(s, lambda_lo) table spot checks (12 cells incl. band boundaries)", not bad, str(bad))
    check("S2b kappa NaN without lambda_lo / s", np.isnan(lt.kappa_of(0.2, np.nan)) and np.isnan(lt.kappa_of(np.nan, 0.5)))

    # -------------------------------------------------------------------------------- S3 eligibility / strata / outputs
    print("\n--- S3 eligibility, strata, outputs ---")
    p_s3 = tmp / "s3.npz"
    make_npz(p_s3, "post", 100, 120, 0.6, 0.3, 0.2, seed=31)
    d3 = lt.load_inputs(str(p_s3))
    ok3 = np.isfinite(lt.observed(d3, ["N", "Nc", "D", "S_cross"], False)[0]["N"])
    nfit3 = (~np.isnan(d3["protein"])).sum(0)
    depth3 = lt.sets_depth(d3)
    sets3 = lt.build_sets(d3, ok3, nfit3, depth3, None)
    prim = [k for k, v in sets3.items() if v["kind"] == "primary"]
    check("S3a base sets and depth tertiles contain only complete genes (n_fit = n)",
          all((nfit3[v["idx"]] == 100).all() for v in sets3.values() if v["kind"] in ("primary", "depth")), ",".join(prim))
    frac = {k: (1 - nfit3[v["idx"]] / 100.0) for k, v in sets3.items() if v["kind"] == "miss"}
    ranges = {"miss0-5": (1e-9, 0.05), "miss5-20": (0.05, 0.20), "miss20+": (0.20, 1.0), "miss<=5": (0.0, 0.05)}
    check("S3b missingness strata hold the right n_fit range",
          all(((f > ranges[k.split("|")[1]][0] - 1e-9) & (f <= ranges[k.split("|")[1]][1] + 1e-9)).all() for k, f in frac.items()) and len(frac) > 0,
          ",".join(f"{k}:{len(f)}" for k, f in frac.items()))
    check("S3c strata below 15 genes are dropped, base sets need >= 5",
          all(len(v["idx"]) >= (lt.MIN_SET if v["kind"] == "primary" else lt.MIN_STRATUM) for v in sets3.values()))
    pin = sets3.get("pinned_sig")
    check("S3d pinned set: s = n_sel/G_tested of complete genes, mnar flag = (n_complete < 20)",
          pin is not None and abs(pin["s"] - pin["n_sel"] / pin["G_tested"]) < 1e-12 and pin["mnar_exposed"] == (pin["n_complete"] < 20),
          f"n_sel {pin['n_sel']:.0f} G_tested {pin['G_tested']:.0f} mnar {pin['mnar_exposed']}")
    # tiny pinned set -> mnar_exposed and read on miss<=5
    p_s3b = tmp / "s3b.npz"
    make_npz(p_s3b, "post", 100, 40, 0.6, 0.3, 0.2, seed=32, miss=0.03)
    d3b = lt.load_inputs(str(p_s3b))
    ok3b = np.isfinite(lt.observed(d3b, ["N", "Nc", "D", "S_cross"], False)[0]["N"])
    s3b = lt.build_sets(d3b, ok3b, (~np.isnan(d3b["protein"])).sum(0), lt.sets_depth(d3b), None)
    n_c3b = {k: v["n_complete"] for k, v in s3b.items() if v["kind"] == "primary"}
    check("S3e mnar_exposed == (fewer than 20 complete genes in the base set), incl. a pinned set below 20",
          all(v["mnar_exposed"] == (v["n_complete"] < 20) for v in s3b.values()) and any(v["mnar_exposed"] for v in s3b.values())
          and any(not v["mnar_exposed"] for v in s3b.values()), str(n_c3b))
    t0 = time.time()
    out3 = lt.main(["--inputs", str(p_s3), "--synthetic", "--strata", "operator,none", "--B", "20", "--boot", "20", "--workers", "1",
                    "--out", str(tmp / "out_S3")])
    print(f"   S3 full-feature run: {time.time()-t0:.0f}s")
    od = tmp / "out_S3"
    files = ["theta_syn_sets.csv", "theta_syn_reading.csv", "theta_syn_hetero.csv", "theta_syn_variants.csv", "theta_syn_pergene.csv",
             "theta_syn_manifest.json", "theta_syn_boot.npz"]
    check("S3f all output files written", all((od / f).exists() for f in files), ",".join(f for f in files if not (od / f).exists()))
    rd = pd.read_csv(od / "theta_syn_reading.csv")
    need = ["z_D_minus_gate", "D_z", "Nc_z", "Nc_z_unres", "theta_pc", "theta_pc_ci_lo", "theta_pc_ci_hi", "s", "kappa", "lambda_lo",
            "lambda_lo_source", "lambda_hi", "theta_noise_min", "margin_theta_pc_minus_tnm", "margin_ci_hi_minus_tnm",
            "r1_band_limit", "margin_abs_theta_pc_minus_band", "fX_lo", "fX_hi", "post_share_lb", "phi_op_D", "phi_op_dp",
            "Wm_z", "spline_theta_pc", "rnapc_theta_pc", "mnar_exposed", "n_complete"]
    check("S3g reading.csv carries the numbers the rules read", all(c in rd.columns for c in need), ",".join(c for c in need if c not in rd.columns))
    check("S3h no verdict / rule-label column anywhere (numbers only)",
          not any(("verdict" in c.lower() or c.lower() in ("rule", "reading", "call", "r1", "r2", "r3", "r4", "r4w", "r5", "r6", "r7")) for c in rd.columns)
          and not any(rd[c].astype(str).isin(["R1", "R2", "R3", "R4", "R4w", "R5", "R6", "R7"]).any() for c in rd.columns if rd[c].dtype == object))
    allr = rd[rd["set"] == "all"].iloc[0]
    check("S3i lambda_lo falls back to rho2bar without ICC and kappa = 1 for an unselected set",
          allr["lambda_lo_source"] == "rho2bar" and allr["kappa"] == 1.0 and abs(allr["lambda_lo"] - allr["rho2bar"]) < 1e-9)
    pr = rd[rd["set"] == "pinned_sig"]
    if len(pr):
        pr = pr.iloc[0]
        check("S3j selected set: kappa from the table at its own s", abs(pr["kappa"] - lt.kappa_of(pr["s"], pr["lambda_lo"])) < 1e-12,
              f"s={pr['s']:.3f} lam_lo={pr['lambda_lo']:.3f} kappa={pr['kappa']}")
        check("S3k theta_noise_min = kappa lam_lo/(1-lam_lo)", abs(pr["theta_noise_min"] - pr["kappa"] * pr["lambda_lo"] / (1 - pr["lambda_lo"])) < 1e-12)
    mnar = pd.read_csv(od / "theta_syn_sets.csv")
    check("S3l theta_reading_ok flag = (n_genes >= 20) in sets.csv", ((mnar["n_genes"] >= 20) == mnar["theta_reading_ok"]).all())
    check("S3m primary-arm rows in reading.csv are the operator arm", set(pd.read_csv(od / "theta_syn_sets.csv")["arm"]) == {"operator", "none"}
          and (rd["arm"] == "operator").all())

    # -------------------------------------------------------------------------------- V1 variants
    print("\n--- V1 variants ---")
    ob = lt.observed(d, ["N", "Nc", "D", "S_cross"], False)[0]
    vv = lt.run_draws(d, d["pcs"][None], None, ["N", "Nc", "D"], False, 1, kind="var:linear", genes=comp[:10])
    w1 = max(max(abs(vv["N"][0, g] - ob["N"][g]), abs(vv["D"][0, g] - ob["D"][g]), abs(vv["Nc"][0, g] - ob["Ncs"][g])) for g in comp[:10])
    check("V1a variant path 'linear' == kernel N, D and N_c = Nc/|rho|", w1 < 1e-10, f"{w1:.1e}")
    q = np.linspace(-2, 2, 400)
    Bs = lt.ns_basis(q, np.array([-2.0, -0.6, 0.6, 2.0]))
    check("V1b natural cubic spline basis: 3 columns, linear beyond the boundary knots",
          Bs.shape == (400, 3) and np.allclose(np.diff(lt.ns_basis(np.array([2.0, 3.0, 4.0]), np.array([-2.0, -0.6, 0.6, 2.0])), 2, axis=0), 0, atol=1e-9))
    pth = tmp / "pop_0.npz"        # created below in P1, reuse by regenerating here
    pop0, _ = make_npz(pth, "noise", 1500, 24, 0.6, 0.3, 0.2, fX=0.5, seed=100, with_pub=False)
    sv = lt.main(["--inputs", str(pth), "--synthetic", "--strata", "none", "--B", "20", "--boot", "0", "--workers", "1",
                  "--variants", "spline,rnapc", "--no-diag", "--no-phi-delta", "--out", str(tmp / "out_V1")])
    vdf = pd.read_csv(tmp / "out_V1" / "theta_syn_variants.csv")
    rs_ = row_of(sv, "all")
    sp = vdf[(vdf["variant"] == "spline") & (vdf["set"] == "all")].iloc[0]
    rn_ = vdf[(vdf["variant"] == "rnapc") & (vdf["set"] == "all")].iloc[0]
    check("V1c linear truth: theta_pc^spline ~ theta_pc (|diff| < 0.1)", abs(sp["theta_pc"] - rs_["theta_pc"]) < 0.1,
          f"spline {sp['theta_pc']:+.3f} linear {rs_['theta_pc']:+.3f}")
    check("V1d spline nested increment is reported next to the published one", np.isfinite(sp["dp_spline_mean"]) and np.isfinite(sp["dp_linear_published_mean"]) is not None)
    check("V1e irrelevant RNA PCs leave theta ~ unchanged (|diff| < 0.15)", abs(rn_["theta_pc"] - rs_["theta_pc"]) < 0.15,
          f"rnapc {rn_['theta_pc']:+.3f} linear {rs_['theta_pc']:+.3f}")

    # -------------------------------------------------------------------------------- V2 / V3 scenarios
    print("\n--- V2 shared RNA artefact, V3 purity ---")
    p_rna = tmp / "scen_rna.npz"
    scenario(p_rna, "rna")
    o_rna = lt.main(["--inputs", str(p_rna), "--synthetic", "--strata", "none", "--B", "60", "--boot", "0", "--workers", "1",
                     "--variants", "rnapc", "--no-diag", "--no-phi-delta", "--out", str(tmp / "out_V2")])
    v2 = pd.read_csv(tmp / "out_V2" / "theta_syn_variants.csv").iloc[0]
    r2 = row_of(o_rna, "all")
    print(f"   rna artefact: primary Nc_z {r2['Nc_z']:+.1f} theta_pc {r2['theta_pc']:+.2f}; 20-RNA-PC Nc_z {v2['Nc_z']:+.1f} theta_pc {v2['theta_pc']:+.2f}")
    check("V2a primary reads the shared RNA artefact as anti-alignment (N_c z <= -3, theta_pc <= -0.3)", r2["Nc_z"] <= -3 and r2["theta_pc"] <= -0.3)
    check("V2b the 20-RNA-PC baseline removes it (|N_c z| < 3)", abs(v2["Nc_z"]) < 3)
    p_pur = tmp / "scen_pur.npz"
    pur = scenario(p_pur, "purity")
    ppath = tmp / "purity.tsv"
    pd.DataFrame({"patient": [f"S{i:04d}" for i in range(len(pur))], "purity": pur}).to_csv(ppath, sep="\t", index=False)
    o_pur = lt.main(["--inputs", str(p_pur), "--synthetic", "--strata", "none", "--B", "60", "--boot", "0", "--workers", "1",
                     "--variants", "purity", "--purity-file", str(ppath), "--no-diag", "--no-phi-delta", "--out", str(tmp / "out_V3")])
    v3 = pd.read_csv(tmp / "out_V3" / "theta_syn_variants.csv").iloc[0]
    r3 = row_of(o_pur, "all")
    print(f"   purity channel: primary Nc_z {r3['Nc_z']:+.1f} theta_pc {r3['theta_pc']:+.2f}; purity-residualised Nc_z {v3['Nc_z']:+.1f} theta_pc {v3['theta_pc']:+.2f}")
    check("V3a primary reads purity mixing as alignment (N_c z >= 4, theta_pc > 0.3)", r3["Nc_z"] >= 4 and r3["theta_pc"] > 0.3)
    check("V3b the purity-residualised variant removes it (|N_c z| < 3)", abs(v3["Nc_z"]) < 3)

    # -------------------------------------------------------------------------------- D1 nested increment
    print("\n--- D1 nested increment machinery ---")
    pm = tk.build_perms(100, 6, "ra", 9)
    gl = comp[:5]
    dn = lt.run_draws(d, d["pcs"][pm], None, ["dp"], False, 1, kind="delta", genes=gl)["dp"]
    do = lt.run_draws(d, d["pcs"][None], None, ["dp"], False, 1, kind="delta", genes=gl)["dp"][0]
    w1 = 0.0
    for g in gl:
        v_, nv_, M_, P_ = lt.gene_arrays(d, g)
        ref_null = tk.published_null_increments(d["mrna_log2tpm1"][:, g], d["protein"][:, g], d["pcs"], pm, v_)
        w1 = max(w1, np.max(np.abs(ref_null - dn[:, g])), abs(do[g] - tk.published_increment(d["mrna_log2tpm1"][:, g], d["protein"][:, g], d["pcs"], v_)[3]))
    check("D1 batched nested-increment null/observed == published_null_increments / published_increment", w1 < 1e-9, f"{w1:.1e}")

    # -------------------------------------------------------------------------------- F1 phi_op
    print("\n--- F1 phi_op ---")
    nloc = 300
    ops = np.array([f"op{x}" for x in np.random.default_rng(2).integers(0, 4, size=nloc)])
    p_plain = tmp / "scen_plain.npz"
    scenario(p_plain, "plain", n=nloc, op_labels=ops)
    p_opm = tmp / "scen_opm.npz"
    scenario(p_opm, "opmeans", n=nloc, op_labels=ops)
    ph = {}
    for nm, pp in (("plain", p_plain), ("opmeans", p_opm)):
        o = lt.main(["--inputs", str(pp), "--synthetic", "--strata", "operator,none", "--B", "40", "--boot", "0", "--workers", "1",
                     "--variants", "none", "--no-diag", "--out", str(tmp / f"out_F1_{nm}")])
        rdf = pd.read_csv(tmp / f"out_F1_{nm}" / "theta_syn_reading.csv")
        ph[nm] = rdf[rdf["set"] == "all"].iloc[0]
        print(f"   {nm:8s}: phi_op_D {ph[nm]['phi_op_D']:+.2f}  phi_op_dp {ph[nm]['phi_op_dp']:+.2f}  z_N within {ph[nm]['Nc_z']:+.1f} unres {ph[nm]['Nc_z_unres']:+.1f}")
    check("F1a no batch structure: |phi_op| < 0.25", abs(ph["plain"]["phi_op_D"]) < 0.25 and abs(ph["plain"]["phi_op_dp"]) < 0.25)
    check("F1b W encodes operator, P carries operator means: phi_op > 0.5 (D and nested increment)",
          ph["opmeans"]["phi_op_D"] > 0.5 and ph["opmeans"]["phi_op_dp"] > 0.5)

    # -------------------------------------------------------------------------------- M1 MNAR toy
    print("\n--- M1 MNAR toy (H_post, bottom 10% of P truncated in the second half) ---")
    p_mn = tmp / "mnar.npz"
    G_m = 100 if QUICK else 160
    make_npz(p_mn, "post", 100, G_m, 0.7, 0.3, 0.2, seed=41, mnar_q=0.10, with_pub=False)
    om = lt.main(["--inputs", str(p_mn), "--synthetic", "--strata", "none", "--B", "100", "--boot", "0", "--workers", "1"] + FAST
                 + ["--out", str(tmp / "out_M1")])
    rc = row_of(om, "all")
    rm = row_of(om, "all|miss5-20")
    print(f"   complete genes (n={rc['n_genes']}): N_c z {rc['Nc_z']:+.2f} theta_pc {rc['theta_pc']:+.3f};  "
          f"missingness 5-20% (n={rm['n_genes']}): N_c z {rm['Nc_z']:+.2f} theta_pc {rm['theta_pc']:+.3f}")
    check("M1a complete-gene stratum is unaffected (|N_c z| < 3.5)", abs(rc["Nc_z"]) < 3.5)
    check("M1b truncated stratum reads lower than the complete stratum (collider)", rm["Nc_z"] < rc["Nc_z"] - 1.0, f"{rm['Nc_z']:+.2f} vs {rc['Nc_z']:+.2f}")

    # -------------------------------------------------------------------------------- P1 population values
    print("\n--- P1 population values (n large) ---")
    scen = [("noise", 0.6, 0.3, 0.2, 0.5), ("noise", 0.8, 0.3, 0.2, 0.5), ("post", 0.6, 0.3, 0.2, 0.5),
            ("err", 0.6, 0.3, 0.2, 0.5), ("mix", 0.6, 0.3, 0.2, 0.3), ("mix", 0.8, 0.3, 0.2, 0.5)]
    n_big, G_big = (1500, 24) if QUICK else (4000, 40)
    for i, (hyp, lam, R2, dl, fX) in enumerate(scen):
        pth = tmp / f"pop_{i}.npz"
        pop, half = make_npz(pth, hyp, n_big, G_big, lam, R2, dl, fX=fX, seed=100 + i, with_pub=False)
        sets = lt.main(["--inputs", str(pth), "--synthetic", "--strata", "none", "--B", "8", "--boot", "0",
                        "--workers", "1", "--out", str(tmp / f"out_pop_{i}")] + FAST)
        row = row_of(sets, "all")
        th_pop = pop["theta_pop"]
        tol = 0.12 * max(1.0, abs(th_pop)) if abs(th_pop) < 2 else 0.2 * abs(th_pop)
        check(f"P1.{i} {hyp} lam={lam} fX={fX if hyp=='mix' else '-'}: theta_pc {row['theta_pc']:+.3f} vs population {th_pop:+.3f}",
              abs(row["theta_pc"] - th_pop) < tol, f"tol {tol:.2f}; n_genes {row['n_genes']}; lambda_hat {row['lambda_hat_noise']:.3f}")
        mrow = sets[(sets["set"] == "all|miss5-20")]
        if len(mrow):
            check(f"P1.{i} MCAR-missing stratum (5-20%) gives the same theta_pc", abs(mrow.iloc[0]["theta_pc"] - th_pop) < 1.5 * tol,
                  f"{mrow.iloc[0]['theta_pc']:+.3f}")
        # S_cross sign rule: S > 0 iff theta > rho^2/(1-rho^2) ~ R2
        if abs(th_pop - R2) > 0.15:
            check(f"P1.{i} S_cross sign matches (theta > R2_rna)", (row["S_cross_mean"] > 0) == (th_pop > R2),
                  f"S {row['S_cross_mean']:+.3f}, R2 {R2}")

    # -------------------------------------------------------------------------------- P2 set-level tests at n = 100
    print("\n--- P2 set-level tests, n=100, K=20, primary = operator-stratified, sensitivity = unrestricted ---")
    Bp = 60 if QUICK else 200
    Gp = 80 if QUICK else 120
    z_out = {}
    for hyp, lam, fX in (("noise", 0.6, 0.5), ("post", 0.6, 0.5), ("err", 0.6, 0.5)):
        pth = tmp / f"p2_{hyp}.npz"
        pop, half = make_npz(pth, hyp, 100, Gp, lam, 0.3, 0.2, fX=fX, seed=7, with_pub=False, miss=0.0)
        t0 = time.time()
        sets = lt.main(["--inputs", str(pth), "--synthetic", "--strata", "operator,none", "--B", str(Bp), "--boot", "100",
                        "--workers", "1", "--variants", "none", "--no-phi-delta", "--out", str(tmp / f"out_p2_{hyp}")])
        print(f"   {hyp}: {time.time()-t0:.0f}s  (theta_pop {pop['theta_pop']:+.2f}, feasible {pop['feasible']})")
        for arm in ("operator", "none"):
            r = row_of(sets, "all", arm)
            z_out[(hyp, arm)] = r
            print(f"   {hyp:5s} {arm:8s}: Nc_z {r['Nc_z']:+7.2f} (rho-weighted {r['Nc_rhow_z']:+6.2f})  N_z {r['N_z']:+7.2f}  D_z {r['D_z']:+7.2f}  "
                  f"theta_pc {r['theta_pc']:+.2f} [{r['theta_pc_ci_lo']:+.2f},{r['theta_pc_ci_hi']:+.2f}]  S_z {r['S_cross_z']:+.2f}  Wm_z {r['Wm_z']:+.2f}  "
                  f"band below/above {r['frac_Nc_below_band']:.2f}/{r['frac_Nc_above_band']:.2f}")
    for arm in ("operator", "none"):
        check(f"P2 H_noise {arm}: N_c,set z > 10", z_out[("noise", arm)]["Nc_z"] > 10, f"{z_out[('noise', arm)]['Nc_z']:+.1f}")
        check(f"P2 H_post  {arm}: -3 < N_c,set z < 2", -3 < z_out[("post", arm)]["Nc_z"] < 2, f"{z_out[('post', arm)]['Nc_z']:+.1f}")
        check(f"P2 H_err   {arm}: N_c,set z < -10", z_out[("err", arm)]["Nc_z"] < -10, f"{z_out[('err', arm)]['Nc_z']:+.1f}")
        check(f"P2 D_set z > 5 under every hypothesis ({arm})", all(z_out[(h, arm)]["D_z"] > 5 for h in ("noise", "post", "err")))
    check("P2 theta_pc bootstrap CI excludes 0 under H_noise", z_out[("noise", "operator")]["theta_pc_ci_lo"] > 0)
    check("P2 theta_pc bootstrap CI contains 0 under H_post",
          z_out[("post", "operator")]["theta_pc_ci_lo"] <= 0 <= z_out[("post", "operator")]["theta_pc_ci_hi"],
          f"[{z_out[('post', 'operator')]['theta_pc_ci_lo']:+.2f}, {z_out[('post', 'operator')]['theta_pc_ci_hi']:+.2f}]")
    check("P2 theta_pc ~ 0 under H_post and < -0.3 under H_err (operator arm)",
          abs(z_out[("post", "operator")]["theta_pc"]) < 0.15 and z_out[("err", "operator")]["theta_pc"] < -0.3)
    check("P2 diagnostic Wm: large positive z only under H_noise (W reads M)",
          z_out[("noise", "operator")]["Wm_z"] > 5 and abs(z_out[("post", "operator")]["Wm_z"]) < 3.5,
          f"noise {z_out[('noise', 'operator')]['Wm_z']:+.1f} post {z_out[('post', 'operator')]['Wm_z']:+.1f}")

    # -------------------------------------------------------------------------------- G1 gate
    print("\n--- G1 fail-closed gate ---")
    real = tmp / "real_flag.npz"
    make_npz(real, "noise", 100, 40, 0.6, 0.3, 0.2, seed=21, synthetic=False)

    def refused(argv, **kw):
        try:
            lt.main(argv, **kw)
        except SystemExit as e:
            return str(e.code)
        return None

    common = ["--strata", "none", "--B", "6", "--boot", "0", "--workers", "1", "--out", str(tmp / "out_gate")] + FAST
    msg = refused(["--inputs", str(real)] + common)
    check("G1a real file, no sidecar -> REFUSED", msg is not None and msg.startswith("REFUSED") and "sidecar" in msg, msg or "ran!")
    msg = refused(["--inputs", str(real), "--synthetic"] + common)
    check("G1b --synthetic on a file not flagged synthetic -> REFUSED", msg is not None and msg.startswith("REFUSED"), msg or "ran!")
    msg = refused(["--inputs", str(real), "--repro-only", "--out", str(tmp / "out_gate")])
    check("G1c --repro-only runs without sidecar (exit 0, published quantities only)", msg in (None, "0") or msg == "0", str(msg))
    root = tmp / "root"
    (root / "review").mkdir(parents=True, exist_ok=True)
    prespec = root / lt.PRESPEC_REL
    prespec.write_text("stub pre-specification\n")

    def write_sidecar(path, files, hashes=None):
        with open(path, "w") as f:
            for i, fp in enumerate(files):
                f.write(f"{(hashes or {}).get(i) or lt.sha256(fp)} *{fp}\n")

    sc = root / "review" / "sidecar.sha256"
    base_files = [str(prespec), str(HERE / "local_theta.py"), str(HERE / "theta_kernel.py"), str(real)]
    write_sidecar(sc, base_files)
    msg = refused(["--inputs", str(real)] + common, sidecar=str(sc), root=str(root))
    check("G1d valid sidecar -> runs", msg is None, str(msg))
    write_sidecar(sc, base_files, hashes={3: "0" * 64})
    msg = refused(["--inputs", str(real)] + common, sidecar=str(sc), root=str(root))
    check("G1e wrong inputs hash in sidecar -> REFUSED", msg is not None and msg.startswith("REFUSED"), msg or "ran!")
    write_sidecar(sc, [f for f in base_files if f != str(real)])
    msg = refused(["--inputs", str(real)] + common, sidecar=str(sc), root=str(root))
    check("G1f inputs file not listed -> REFUSED", msg is not None and msg.startswith("REFUSED") and "not listed" in msg, msg or "ran!")
    write_sidecar(sc, [f for f in base_files if not f.endswith("local_theta.py")])
    msg = refused(["--inputs", str(real)] + common, sidecar=str(sc), root=str(root))
    check("G1g script not listed -> REFUSED", msg is not None and msg.startswith("REFUSED"), msg or "ran!")
    write_sidecar(sc, [f for f in base_files if f != str(prespec)])
    msg = refused(["--inputs", str(real)] + common, sidecar=str(sc), root=str(root))
    check("G1h pre-spec not listed -> REFUSED", msg is not None and msg.startswith("REFUSED"), msg or "ran!")
    # side files (purity / sets) must be frozen too
    ppath2 = tmp / "purity_real.tsv"
    pd.DataFrame({"patient": [f"S{i:04d}" for i in range(100)], "purity": np.linspace(0, 1, 100)}).to_csv(ppath2, sep="\t", index=False)
    write_sidecar(sc, base_files)
    msg = refused(["--inputs", str(real), "--purity-file", str(ppath2)] + common, sidecar=str(sc), root=str(root))
    check("G1l purity file not listed in the sidecar -> REFUSED", msg is not None and msg.startswith("REFUSED") and "not listed" in msg, msg or "ran!")
    write_sidecar(sc, base_files + [str(ppath2)])
    msg = refused(["--inputs", str(real), "--purity-file", str(ppath2)] + common, sidecar=str(sc), root=str(root))
    check("G1m purity file listed -> runs", msg is None, str(msg))
    # real-mode numbers == synthetic-mode numbers on the same arrays
    write_sidecar(sc, base_files)
    real_syn = tmp / "real_as_syn.npz"
    zz = dict(np.load(real, allow_pickle=False))
    zz["is_synthetic"] = np.bool_(True)
    np.savez_compressed(real_syn, **zz)
    a = lt.main(["--inputs", str(real), "--strata", "none", "--B", "6", "--boot", "0", "--workers", "1", "--out", str(tmp / "o1")] + FAST,
                sidecar=str(sc), root=str(root))
    b = lt.main(["--inputs", str(real_syn), "--synthetic", "--strata", "none", "--B", "6", "--boot", "0", "--workers", "1", "--out", str(tmp / "o2")] + FAST)
    check("G1i real-mode (valid sidecar) == synthetic-mode numbers", np.allclose(a.select_dtypes("number").values, b.select_dtypes("number").values, equal_nan=True))
    full = ["--strata", "operator,none", "--B", "6", "--boot", "5", "--variants", "spline,rnapc"]
    c = lt.main(["--inputs", str(real_syn), "--synthetic", "--workers", "2", "--out", str(tmp / "o3")] + full)
    c1 = lt.main(["--inputs", str(real_syn), "--synthetic", "--workers", "1", "--out", str(tmp / "o4")] + full)
    check("G1k multiprocessing (spawn, 2 workers) == serial, incl. bootstrap, stratified arm, diagnostic, phi_op and variants",
          np.allclose(c.select_dtypes("number").values, c1.select_dtypes("number").values, equal_nan=True) and
          np.allclose(c.attrs["reading"].select_dtypes("number").values, c1.attrs["reading"].select_dtypes("number").values, equal_nan=True))
    # production CLI constant points where the C1 script's convention puts it, and does not exist yet
    check("G1j production SIDECAR constant absent (nothing frozen)", not os.path.exists(lt.SIDECAR) and lt.SIDECAR.endswith(".md.sha256"), lt.SIDECAR)

    print("\nFAILED: " + ", ".join(FAIL) if FAIL else "\nALL CHECKS PASSED")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":   # required: local_theta uses spawn-based multiprocessing
    run()
