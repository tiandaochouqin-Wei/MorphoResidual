#!/usr/bin/env python3
"""
End-to-end plumbing test of the v2 EXPORT job (measurement_error_checks.py) on a SYNTHETIC cohort
(random numbers made up here; no project data).  Builds fake STAR-Counts files (with replicate
files for 20 cases), a TMT table, a crosswalk, a slide map and .pt embeddings in a temp dir, points
the REAL published loaders (server_export/scripts/residual_analysis.py) at them, computes the
"published" per-gene results with RA's own cv_r2 and RA's own permutation sequence, runs the export,
and checks
  E1  the raw-count re-read reproduces RA's RNA matrix (asserted inside the script)
  E2  inputs_<c>.npz: float64 arrays (pcs/wsi keep RA's dtype), patients/genes alignment, labels,
      published arrays, RNA PCs == transcriptome_baseline.py's recipe, n_rna_files
  E3  repro_check says EXACT; the published increments re-derived from the file ON DISK equal RA's
      to 1e-10 for EVERY fit gene (and the job would exit 3 otherwise: tested with a perturbed CSV)
  E4  local_theta.py --repro-only / --check-pca / --check-icc / --repro-pvals on the exported file
  E5  the default job writes NO discriminator file (blind by construction); --B>0 without
      --prespec-sha256 is refused; with it the optional W-only path equals the reference statistics
      from RA's own permutation sequence; --resume reproduces the same numbers
  E6  export_sha256 matches; failed label export is recorded, not fatal; the exported (real-flagged)
      file is REFUSED by local_theta.py without a sidecar (CLI, production sidecar constant)
Usage:  python test_synthetic.py <tmpdir> [path-to-residual_analysis-dir]
"""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

tmp = Path(sys.argv[1]).resolve()
ra_dir = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else \
    Path(__file__).resolve().parents[2] / "server_export" / "scripts"
here = Path(__file__).resolve().parent
sys.path.insert(0, str(here))
tmp.mkdir(parents=True, exist_ok=True)
FAIL = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL.append(name)


rs = np.random.RandomState(7)
N, G = 100, 150
cases = [f"C3L-{i:05d}" for i in range(1, N + 1)]
genes = [f"GENE{i}" for i in range(G)]
rna_dir = tmp / "omics" / "rna"
rna_dir.mkdir(parents=True, exist_ok=True)
(tmp / "omics" / "protein").mkdir(parents=True, exist_ok=True)
emb_dir = tmp / "WSI" / "emb"
emb_dir.mkdir(parents=True, exist_ok=True)

lat = rs.randn(N, 3)
load = rs.randn(3, G)
true = lat @ load
rows = []
for ci, c in enumerate(cases):
    nfile = 2 if ci % 5 == 0 else 1
    for k in range(nfile):
        fn = f"{c}_{k}.rna_seq.augmented_star_gene_counts.tsv"
        expr = np.exp(0.6 * true[ci] + rs.randn(G) * 0.3 + 2.0)
        cnt = rs.poisson(expr * 20).astype(int)
        tpm = cnt / cnt.sum() * 1e6
        body = ["gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded"]
        for nm, v in (("N_unmapped", 5), ("N_multimapping", 7), ("N_noFeature", 9), ("N_ambiguous", 3)):
            body.append(f"{nm}\t\t\t{v}\t{v}\t{v}\t\t\t")
        for j, g in enumerate(genes):
            body.append(f"ENSG{j:08d}.1\t{g}\tprotein_coding\t{cnt[j]}\t0\t0\t{tpm[j]:.4f}\t{tpm[j]:.4f}\t{tpm[j]:.4f}")
        (rna_dir / fn).write_text("# gene-model: GENCODE v36\n" + "\n".join(body) + "\n")
        rows.append((f"fid{ci}_{k}", fn, c, "Primary Tumor"))
man = tmp / "manifest.tsv"
pd.DataFrame(rows, columns=["file_id", "file_name", "case_submitter_id", "sample_type"]).to_csv(man, sep="\t", index=False)

prot = 0.5 * true + 0.8 * lat[:, [0]] * rs.randn(1, G) + rs.randn(N, G) * 0.7
prot[rs.rand(N, G) < 0.12] = np.nan
cols = {f"CPT{ci:07d} Log Ratio": prot[ci] for ci in range(N)}
pd.DataFrame(cols, index=genes).to_csv(tmp / "omics" / "protein" / "fake.tmt10.tsv", sep="\t")
pd.DataFrame({"aliquot_id": [f"CPT{ci:07d}" for ci in range(N)], "case_id": cases,
              "sample_type": "Primary Tumor"}).to_csv(tmp / "xwalk.tsv", sep="\t", index=False)
slides = []
for ci, c in enumerate(cases):
    for s in (21, 22):
        sid = f"{c}-{s}"
        slides.append((sid, "Primary Tumor", c))
        torch.save({"embeddings": torch.tensor(rs.randn(30, 32) + lat[ci, 0] * 0.5 + 0.3 * true[ci, :32 % G].mean(), dtype=torch.float32)},
                   emb_dir / f"{sid}.pt")
pd.DataFrame(slides, columns=["slide_submitter_id", "sample_type", "case_submitter_id"]).to_csv(
    tmp / "slides.tsv", sep="\t", index=False)
lab_df = pd.DataFrame({"case_id": cases, "operator": [f"op{i % 5}" for i in range(N)],
                       "plex": [f"plex{i % 9:02d}" for i in range(N)]})
lab_df.to_csv(tmp / "labels.tsv", sep="\t", index=False)

env = dict(os.environ, MORPHO_ROOT=str(tmp), MORPHO_RNA_MANIFEST=str(man),
           MORPHO_ALIQUOT_XWALK=str(tmp / "xwalk.tsv"), MORPHO_SLIDE_MAP=str(tmp / "slides.tsv"),
           MORPHO_PROTEIN_TSV=str(tmp / "omics" / "protein" / "fake.tmt10.tsv"),
           MORPHO_WSI_EMB_DIR=str(emb_dir), MORPHO_OUT=str(tmp / "out_ra"),
           MORPHO_SCRIPTS_DIR=str(ra_dir), PYTHONIOENCODING="utf-8")
env.pop("MORPHO_SIG_CSV", None)

# ---------------- "published" results: RA's own loaders + RA.cv_r2 + RA's own permutation sequence ----------------
sys.path.insert(0, str(ra_dir))
os.environ.update(env)
import residual_analysis as RA   # noqa: E402
from sklearn.decomposition import PCA   # noqa: E402
import theta_kernel as tk   # noqa: E402

rna, protein, wsi = RA.load_rna_matrix(), RA.load_protein_matrix(), RA.load_wsi_embeddings()
common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
cg = sorted(set(rna.columns) & set(protein.columns))
rna_full = rna.loc[common]
rna, protein = rna.loc[common, cg], protein.loc[common, cg]
pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.loc[common].values)
rng = np.random.RandomState(0)
perms_ra = [rng.permutation(len(common)) for _ in range(RA.N_PERM)]
pub = []
for g in cg:
    y = protein[g].values
    v = ~np.isnan(y)
    if v.sum() < 30:
        continue
    yv, xr = y[v], rna[g].values[v].reshape(-1, 1)
    r2_rna = RA.cv_r2(xr, yv, 5)
    r2_both = RA.cv_r2(np.hstack([xr, pcs[v]]), yv, 5)
    pub.append(dict(gene=g, n=int(v.sum()), r2_rna=r2_rna, r2_mrna_wsi=r2_both, incremental_r2=r2_both - r2_rna))
pub = pd.DataFrame(pub)
pub["pval"] = np.nan
for i in (0, 1, 2):                                    # three genes get RA's genuine 1000-permutation p-value
    g = pub.gene.iloc[i]
    y = protein[g].values
    v = ~np.isnan(y)
    xr = rna[g].values[v].reshape(-1, 1)
    null = np.array([RA.cv_r2(np.hstack([xr, pcs[p][v]]), y[v], 5) for p in perms_ra]) - pub.r2_rna.iloc[i]
    pub.loc[pub.index[i], "pval"] = (np.sum(null >= pub.incremental_r2.iloc[i]) + 1) / (len(null) + 1)
pub["fdr"] = np.where(pub.incremental_r2 > 0.02, 0.01, 0.4)
pub_csv = tmp / "published_snapshot.csv"
pub.to_csv(pub_csv, index=False)

# ---------------- run the export ----------------
out = tmp / "out"
cmd = [sys.executable, "-u", str(here / "measurement_error_checks.py"), "--cohort", "ucec", "--out", str(out),
       "--workers", "2", "--labels-file", str(tmp / "labels.tsv"), "--skip-expected-check"]
env_run = dict(env, MORPHO_SIG_CSV=str(pub_csv))
r = subprocess.run(cmd, env=env_run, capture_output=True, text=True)
print(r.stdout[-2500:])
print(r.stderr[-1500:])
check("export job exit code 0", r.returncode == 0, f"rc={r.returncode}")
assert r.returncode == 0, "export failed"

z = np.load(out / "inputs_ucec.npz", allow_pickle=False)
fit_genes = list(z["genes"])
check("E1 depth re-read line present", "max |log2(TPM+1) re-read - RA matrix| = 0.00e+00" in r.stdout)
check("E2a float64 arrays", all(z[k].dtype == np.float64 for k in ("mrna_log2tpm1", "protein", "counts", "tpm", "rna_pcs20")))
check("E2b pcs/wsi keep RA's dtype", z["pcs"].dtype == pcs.dtype and z["wsi_pooled"].dtype == wsi.loc[common].values.dtype,
      f"pcs {z['pcs'].dtype}")
check("E2c patients == common, genes == fit genes", list(z["patients"]) == common and fit_genes == [g for g in cg if (~np.isnan(protein[g].values)).sum() >= 30])
check("E2d pcs bit-identical to the published pipeline's PCA", np.array_equal(z["pcs"], pcs))
check("E2e mRNA/protein bit-identical", np.array_equal(z["mrna_log2tpm1"], rna[fit_genes].values) and
      np.array_equal(z["protein"], protein[fit_genes].to_numpy(float), equal_nan=True))
check("E2f counts consistent with TPM export (non-negative, finite)", np.isfinite(z["counts"]).all() and (z["counts"] >= 0).all() and np.isfinite(z["tpm"]).all())
check("E2g labels", list(z["label_operator"]) == lab_df.operator.tolist() and list(z["label_plex"]) == lab_df.plex.tolist())
check("E2h is_synthetic flag False (a real export must not be mistakable for test data)", not bool(z["is_synthetic"]))
# RNA PCs per transcriptome_baseline.py
rv = rna_full.values.astype(float)
sdv = rv.std(0)
rz = (rv - rv.mean(0)) / np.where(sdv == 0, 1.0, sdv)
tb = PCA(n_components=min(20, len(common) - 10, rz.shape[1]), svd_solver="full").fit_transform(rz)
check("E2i RNA PCs == transcriptome_baseline.py recipe", np.max(np.abs(tb - z["rna_pcs20"])) < 1e-9, f"{np.max(np.abs(tb - z['rna_pcs20'])):.1e}")
pubz = pub.set_index("gene").loc[fit_genes]
check("E2j published arrays aligned", np.array_equal(z["pub_incremental_r2"], pubz.incremental_r2.values) and
      np.array_equal(z["pub_fdr"], pubz.fdr.values))
nf = pd.Series({c: 2 if i % 5 == 0 else 1 for i, c in enumerate(cases)}).reindex(common).values
check("E2k n_rna_files", np.array_equal(z["n_rna_files"], nf))

rep = json.loads((out / "repro_check_ucec.json").read_text())
check("E3a repro_check verdict EXACT", rep["verdict"].startswith("EXACT"), json.dumps({k: rep[k] for k in rep if k.startswith("max")}))
# every gene, independently, with RA.cv_r2 on the arrays read back from disk
worst = 0.0
for gi, g in enumerate(fit_genes):
    y = z["protein"][:, gi]
    v = ~np.isnan(y)
    xr = z["mrna_log2tpm1"][v, gi].reshape(-1, 1)
    a = RA.cv_r2(xr, y[v], 5)
    b = RA.cv_r2(np.hstack([xr, z["pcs"][v]]), y[v], 5)
    worst = max(worst, abs((b - a) - pubz.incremental_r2.iloc[gi]), abs(a - pubz.r2_rna.iloc[gi]))
check("E3b every fit gene: increment from disk == RA's, <= 1e-10", worst <= 1e-10, f"{worst:.1e}")
# a perturbed published CSV must fail the job with exit code 3 (use the first 30 genes only)
bad = pub.copy()
bad.loc[bad.index[2], "incremental_r2"] += 1e-6
bad_csv = tmp / "published_bad.csv"
bad.to_csv(bad_csv, index=False)
rb = subprocess.run(cmd[:6] + [str(tmp / "out_bad")] + cmd[7:] + ["--max-genes", "30"], env=dict(env, MORPHO_SIG_CSV=str(bad_csv)),
                    capture_output=True, text=True)
check("E3c perturbed published CSV -> job exits 3 and says DRIFT/CLOSE", rb.returncode == 3 and ("DRIFT" in rb.stdout or "CLOSE" in rb.stdout), f"rc={rb.returncode}")

# ---------------- local_theta repro-only on the exported file ----------------
lt_cmd = [sys.executable, str(here / "local_theta.py"), "--inputs", str(out / "inputs_ucec.npz")]
r4 = subprocess.run(lt_cmd + ["--repro-only", "--repro-pvals", "3", "--check-pca", "--check-icc"], capture_output=True, text=True,
                    env=dict(os.environ, PYTHONIOENCODING="utf-8"))
print(r4.stdout[-1800:], r4.stderr[-800:])
check("E4a local_theta --repro-only exit 0", r4.returncode == 0)
ro = json.loads((out / "local_theta" / "repro_only_ucec.json").read_text())
check("E4b published p-values re-derived (RA sequence, 3 genes)", ro["pvals"]["n_genes"] == 3 and ro["pvals"]["n_mismatch"] == 0, str(ro["pvals"]))
check("E4c PCA re-derivable from wsi_pooled", ro["pca_max_abs_diff_abs_scores"] < 1e-4, f"{ro['pca_max_abs_diff_abs_scores']:.1e}")
check("E4d ICC re-derivable from replicate_files", ro["icc"] is not None and ro["icc"]["max_abs_diff_icc"] < 1e-10 and ro["icc"]["n_cases"] == 20, str(ro["icc"]))

# ---------------- blindness / optional path ----------------
nm = sorted(p.name for p in out.iterdir())
check("E5a default job wrote no discriminator file", not any(x.startswith(("wonly_", "perms_", "perm_stats_")) for x in nm), str(nm))
r5 = subprocess.run(cmd[:6] + [str(tmp / "out_refuse")] + cmd[7:] + ["--B", "3"], env=env_run, capture_output=True, text=True)
check("E5b --B>0 without --prespec-sha256 refused", r5.returncode != 0 and "REFUSED" in (r5.stdout + r5.stderr))
out_b = tmp / "out_B"
cmd_b = cmd[:6] + [str(out_b)] + cmd[7:] + ["--B", "3", "--prespec-sha256", "deadbeef" * 8]
r6 = subprocess.run(cmd_b, env=env_run, capture_output=True, text=True)
check("E5c optional path runs", r6.returncode == 0, r6.stderr[-300:])
P_ = np.load(out_b / "wonly_perm_ucec.npy")
O_ = np.load(out_b / "wonly_obs_ucec.npy")
perms3 = tk.build_perms(len(common), 3, "ra", 0)
worst_o = worst_p = 0.0
for gi in range(0, len(fit_genes), 9):
    y = z["protein"][:, gi]
    v = ~np.isnan(y)
    M = z["mrna_log2tpm1"][v, gi]
    lit = tk.gene_stats_literal(z["pcs"][v], M, y[v])
    worst_o = max(worst_o, max(abs(lit[k] - O_[gi, j]) for j, k in enumerate(["rho", "N", "Nc", "D", "S_cross"])))
    for b in range(3):
        l = tk.gene_stats_literal(z["pcs"][perms3[b]][v], M, y[v])
        worst_p = max(worst_p, max(abs(l[k] - P_[gi, b, j]) for j, k in enumerate(["N", "Nc", "D", "S_cross"])))
check("E5d optional on-cluster observed == reference", worst_o < 1e-10, f"{worst_o:.1e}")
check("E5e optional on-cluster null == reference on RA's permutation sequence (float32 storage)", worst_p < 1e-5, f"{worst_p:.1e}")
done = np.load(out_b / "done_ucec.npy")
done[-30:] = False
np.save(out_b / "done_ucec.npy", done)
before = np.load(out_b / "wonly_perm_ucec.npy").copy()
r7 = subprocess.run(cmd_b + ["--resume"], env=env_run, capture_output=True, text=True)
check("E5f --resume reproduces the same numbers", r7.returncode == 0 and np.array_equal(before, np.load(out_b / "wonly_perm_ucec.npy"), equal_nan=True))

# ---------------- checksums, failed labels, fail-closed CLI ----------------
bad_sum = 0
for line in (out / "export_sha256_ucec.txt").read_text().splitlines():
    h, nme = line.split(" *")
    bad_sum += hashlib.sha256((out / nme).read_bytes()).hexdigest() != h
check("E6a export_sha256 matches every file", bad_sum == 0 and len(nm) > 5, f"{len(nm)} files")
out_nl = tmp / "out_nolabels"
r8 = subprocess.run([sys.executable, "-u", str(here / "measurement_error_checks.py"), "--cohort", "ucec", "--out", str(out_nl), "--workers", "2",
                     "--skip-expected-check", "--max-genes", "20"], env=env_run, capture_output=True, text=True)
zn = np.load(out_nl / "inputs_ucec.npz", allow_pickle=False)
status = json.loads(str(zn["label_status"]))
check("E6b failed label export is recorded, not fatal", r8.returncode == 0 and "_error" in status and (zn["label_operator"] == "").all(), str(status)[:120])
r9 = subprocess.run(lt_cmd[:2] + ["--inputs", str(out_nl / "inputs_ucec.npz"), "--repro-only"], capture_output=True, text=True)
check("E6c repro-only on a partial export (20 genes) still works", r9.returncode == 0)
r10 = subprocess.run(lt_cmd + ["--B", "4", "--boot", "0", "--workers", "1", "--out", str(tmp / "should_not_exist")], capture_output=True, text=True,
                     env=dict(os.environ, PYTHONIOENCODING="utf-8"))
check("E6d CLI on the real-flagged export WITHOUT a sidecar is REFUSED and computes nothing",
      r10.returncode != 0 and "REFUSED" in (r10.stdout + r10.stderr) and not (tmp / "should_not_exist").exists(), (r10.stderr or r10.stdout)[-160:])

print("\nFAILED: " + ", ".join(FAIL) if FAIL else "\nALL CHECKS PASSED")
sys.exit(1 if FAIL else 0)
