#!/usr/bin/env python3
"""
End-to-end test of local_theta.py / theta_kernel.py on SYNTHETIC cohorts only (generated here with
theory/sim_signtest.py's simulate_cohort; no project data).  Checks

  K1  batched/weighted ridge == the reference oof_ridge (permutation draw, bootstrap draw == duplicated rows)
  K2  local statistics == theory/sim_signtest.gene_stats (independently written) on complete genes
  K3  within-stratum permutations == residual_analysis_sitepack.py's construction (re-implemented here)
  K4  one-way ICC estimator recovers a known ICC
  P1  population values of DERIVATION.md section 2.1 at large n (theta_set = f_X lam/(1-lam) | 0 | -1,
      mixed = f_X lam/(1-lam)), through the full local_theta pipeline (missing-protein patterns included)
  P2  n = 100, K = 20 set-level tests of DERIVATION.md section 7: N_c,set z >> 0 under H_noise,
      ~ 0 under H_post, << 0 under H_err; D_set z >> 0; both the unstratified and a within-batch arm
  G1  fail-closed gate: no sidecar -> REFUSED; wrong hash -> REFUSED; script not listed -> REFUSED;
      --synthetic on a file not flagged synthetic -> REFUSED; --repro-only works; a valid sidecar runs
      the same numbers as the synthetic mode
  R1  repro check: matches on exported-like inputs, refuses when a published value is perturbed by 1e-6

Usage:  python test_local_theta.py <scratch dir> [--quick]
"""
import json
import os
import shutil
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


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL.append(name)


def make_npz(path, hyp, n, G, lam, R2, delta, fX=0.5, seed=0, miss=0.12, synthetic=True, with_pub=True):
    rng = np.random.default_rng(seed)
    W, M, P, pop = sim.simulate_cohort(rng, n, 20, G, lam, R2, -0.5, hyp, delta, fX)
    P = P.copy()
    half = G // 2
    mask = rng.random((n, G)) < miss
    mask[:, :half] = False                                # first half complete, second half with missing protein
    P[mask] = np.nan
    P[:, :G] = P * 3.0 + 1.0                              # arbitrary scale/location
    M = M * 1.7 + 4.0
    lab_op = rng.integers(0, 5, size=n)
    lab_plex = rng.integers(0, 8, size=n)
    patients = np.array([f"S{i:04d}" for i in range(n)])
    genes = np.array([f"G{i:04d}" for i in range(G)])
    cnt = rng.poisson(30 * np.exp(rng.normal(size=(1, G)) * 1.2 + 3), size=(n, G)).astype(float)
    z = dict(schema_version=np.int64(2), is_synthetic=np.bool_(synthetic), cohort=np.array("syn"),
             patients=patients, genes=genes, mrna_log2tpm1=M, protein=P, counts=cnt, tpm=cnt / 10.0,
             pcs=W, wsi_pooled=W, pca_evr=np.ones(20) / 20, rna_pcs20=rng.normal(size=(n, 20)),
             n_rna_files=np.ones(n), icc_log2tpm1=np.full(G, np.nan), ridge_alpha=np.float64(1.0),
             cv_folds=np.int64(5), ra_seed=np.int64(0), ra_n_perm=np.int64(1000),
             label_operator=np.array([f"op{x}" for x in lab_op]), label_plex=np.array([f"px{x}" for x in lab_plex]),
             label_status=np.array("{}"))
    if with_pub:
        pub = np.array([tk.published_increment(M[:, g], P[:, g], W) for g in range(G)])
        z.update(pub_n=pub[:, 0], pub_r2_rna=pub[:, 1], pub_r2_both=pub[:, 2], pub_incremental_r2=pub[:, 3],
                 pub_pval=np.full(G, np.nan), pub_fdr=np.where(pub[:, 3] > 0.1, 0.01, 0.5))
    np.savez_compressed(path, **z)
    return pop, half


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
          np.array_equal(tk.build_perms(50, 3)[1], np.random.RandomState(0).permutation(50)[:0].tolist() or
                         (lambda r: (r.permutation(50), r.permutation(50))[1])(np.random.RandomState(0))))

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

    # -------------------------------------------------------------------------------- K2
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

    # -------------------------------------------------------------------------------- R1
    rep = lt.repro_published(d)
    check("R1a repro matches exported-like inputs", rep["ok"], str(rep.get("worst")))
    d2 = dict(d)
    d2["pub_incremental_r2"] = d["pub_incremental_r2"].copy()
    d2["pub_incremental_r2"][3] += 1e-6
    check("R1b repro refuses a 1e-6 perturbation", not lt.repro_published(d2)["ok"])

    # -------------------------------------------------------------------------------- P1 population values
    print("\n--- P1 population values (n large) ---")
    scen = [("noise", 0.6, 0.3, 0.2, 0.5), ("noise", 0.8, 0.3, 0.2, 0.5), ("post", 0.6, 0.3, 0.2, 0.5),
            ("err", 0.6, 0.3, 0.2, 0.5), ("mix", 0.6, 0.3, 0.2, 0.3), ("mix", 0.8, 0.3, 0.2, 0.5)]
    n_big, G_big = (1500, 24) if QUICK else (4000, 40)
    for i, (hyp, lam, R2, dl, fX) in enumerate(scen):
        pth = tmp / f"pop_{i}.npz"
        pop, half = make_npz(pth, hyp, n_big, G_big, lam, R2, dl, fX=fX, seed=100 + i, with_pub=False)
        sets = lt.main(["--inputs", str(pth), "--synthetic", "--strata", "none", "--B", "8", "--boot", "0",
                        "--workers", "1", "--out", str(tmp / f"out_pop_{i}")])
        row = sets[sets["set"] == "all"].iloc[0]
        th_pop = pop["theta_pop"]
        tol = 0.12 * max(1.0, abs(th_pop)) if abs(th_pop) < 2 else 0.2 * abs(th_pop)
        check(f"P1.{i} {hyp} lam={lam} fX={fX if hyp=='mix' else '-'}: theta_set {row['theta_set']:+.3f} vs population {th_pop:+.3f}",
              abs(row["theta_set"] - th_pop) < tol, f"tol {tol:.2f}; lambda_hat {row['lambda_hat_noise']:.3f}")
        # S_cross sign rule: S > 0 iff theta > rho^2 = R2
        if abs(th_pop - R2) > 0.15:
            check(f"P1.{i} S_cross sign matches (theta > R2_rna)", (row["S_cross_mean"] > 0) == (th_pop > R2),
                  f"S {row['S_cross_mean']:+.3f}, R2 {R2}")

    # -------------------------------------------------------------------------------- P2 set-level tests at n = 100
    print("\n--- P2 set-level tests, n=100, K=20 ---")
    Bp = 60 if QUICK else 200
    Gp = 80 if QUICK else 120
    z_out = {}
    for hyp, lam, fX in (("noise", 0.6, 0.5), ("post", 0.6, 0.5), ("err", 0.6, 0.5)):
        pth = tmp / f"p2_{hyp}.npz"
        pop, half = make_npz(pth, hyp, 100, Gp, lam, 0.3, 0.2, fX=fX, seed=7, with_pub=False)
        t0 = time.time()
        sets = lt.main(["--inputs", str(pth), "--synthetic", "--strata", "none,operator", "--B", str(Bp), "--boot", "100",
                        "--workers", "1", "--out", str(tmp / f"out_p2_{hyp}")])
        print(f"   {hyp}: {time.time()-t0:.0f}s  (theta_pop {pop['theta_pop']:+.2f}, feasible {pop['feasible']})")
        for arm in ("none", "operator"):
            r = sets[(sets["set"] == "all") & (sets["arm"] == arm)].iloc[0]
            z_out[(hyp, arm)] = r
            print(f"   {hyp:5s} {arm:8s}: Nc_z {r['Nc_z']:+7.2f}  N_z {r['N_z']:+7.2f}  D_z {r['D_z']:+7.2f}  "
                  f"theta_set {r['theta_set']:+.2f} [{r['theta_ci_lo']:+.2f},{r['theta_ci_hi']:+.2f}]  S_z {r['S_cross_z']:+.2f}  "
                  f"band below/above {r['frac_Nc_below_band']:.2f}/{r['frac_Nc_above_band']:.2f}")
    for arm in ("none", "operator"):
        check(f"P2 H_noise {arm}: N_c,set z > 10", z_out[("noise", arm)]["Nc_z"] > 10, f"{z_out[('noise', arm)]['Nc_z']:+.1f}")
        check(f"P2 H_post  {arm}: |N_c,set z| < 3.5", abs(z_out[("post", arm)]["Nc_z"]) < 3.5, f"{z_out[('post', arm)]['Nc_z']:+.1f}")
        check(f"P2 H_err   {arm}: N_c,set z < -10", z_out[("err", arm)]["Nc_z"] < -10, f"{z_out[('err', arm)]['Nc_z']:+.1f}")
        check(f"P2 D_set z > 5 under every hypothesis ({arm})", all(z_out[(h, arm)]["D_z"] > 5 for h in ("noise", "post", "err")))
    check("P2 theta bootstrap CI excludes 0 under H_noise", z_out[("noise", "none")]["theta_ci_lo"] > 0)
    check("P2 theta bootstrap CI contains 0 under H_post",
          z_out[("post", "none")]["theta_ci_lo"] <= 0 <= z_out[("post", "none")]["theta_ci_hi"],
          f"[{z_out[('post', 'none')]['theta_ci_lo']:+.2f}, {z_out[('post', 'none')]['theta_ci_hi']:+.2f}]")

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


    common = ["--strata", "none", "--B", "6", "--boot", "0", "--workers", "1", "--out", str(tmp / "out_gate")]
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
    # real-mode numbers == synthetic-mode numbers on the same arrays
    write_sidecar(sc, base_files)
    real_syn = tmp / "real_as_syn.npz"
    zz = dict(np.load(real, allow_pickle=False))
    zz["is_synthetic"] = np.bool_(True)
    np.savez_compressed(real_syn, **zz)
    a = lt.main(["--inputs", str(real), "--strata", "none", "--B", "6", "--boot", "0", "--workers", "1", "--out", str(tmp / "o1")],
                sidecar=str(sc), root=str(root))
    b = lt.main(["--inputs", str(real_syn), "--synthetic", "--strata", "none", "--B", "6", "--boot", "0", "--workers", "1", "--out", str(tmp / "o2")])
    check("G1i real-mode (valid sidecar) == synthetic-mode numbers", np.allclose(a.select_dtypes("number").values, b.select_dtypes("number").values, equal_nan=True))
    c = lt.main(["--inputs", str(real_syn), "--synthetic", "--strata", "none,operator", "--B", "6", "--boot", "5", "--workers", "2", "--out", str(tmp / "o3")])
    c1 = lt.main(["--inputs", str(real_syn), "--synthetic", "--strata", "none,operator", "--B", "6", "--boot", "5", "--workers", "1", "--out", str(tmp / "o4")])
    check("G1k multiprocessing (spawn, 2 workers) == serial, incl. bootstrap and stratified arm",
          np.allclose(c.select_dtypes("number").values, c1.select_dtypes("number").values, equal_nan=True))
    # production CLI constant points where the C1 script's convention puts it, and does not exist yet
    check("G1j production SIDECAR constant absent (nothing frozen)", not os.path.exists(lt.SIDECAR) and lt.SIDECAR.endswith(".md.sha256"), lt.SIDECAR)

    print("\nFAILED: " + ", ".join(FAIL) if FAIL else "\nALL CHECKS PASSED")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":   # required: local_theta uses spawn-based multiprocessing
    run()
