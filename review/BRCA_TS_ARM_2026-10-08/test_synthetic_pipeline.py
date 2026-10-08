#!/usr/bin/env python
"""Synthetic end-to-end test of the BRCA set-level plumbing (NO project data).

Builds a tiny synthetic cohort (60 cases, 14 genes + 2 genes with < 30 values) in a scratch folder, then:
  T1  runs the PUBLISHED cptac2_residual.py (MorphoResidual_paper/cptac2_residual.py) and the PATCHED one (dx arm, default flags)
      and requires identical results csv / summary csv  (patch v2 does not change a number);
  T2  the patched script's npz: shape genes x 1001, column 0 == csv increments, pval reproduced from the matrix, seeds 0..999,
      row order == csv order, the n<30 genes absent;
  T3  ts arm (sample-code filter drops the 11A slide), slide- and tile-pooling both write their own matrix;
  T4  --perms 7 --max-genes 3 writes *_smoke* files and local_brca_setlevel.load_matrix refuses it;
  T5  load_matrix accepts the three real matrices, label/arm consistency enforced;
  T6  brca_tested_genes.py lists exactly the genes in the results csv;
  T7  set-level statistics computed from the saved matrix equal an independent re-computation from the csv + re-derived nulls.
Usage:  python test_synthetic_pipeline.py [scratch_dir]
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
PUBLISHED = os.path.join(PAPER, "cptac2_residual.py")
PATCHED = os.path.join(HERE, "scripts", "cptac2_residual.py")
TESTED = os.path.join(HERE, "scripts", "brca_tested_genes.py")
sys.path.insert(0, HERE)
import local_brca_setlevel as L  # noqa: E402

scratch = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="brca_syn_")
base = os.path.join(scratch, "ZW"); root = os.path.join(base, "cptac2_brca")
for d in ("rna", "emb_phikon", "emb_phikon_ts"):
    os.makedirs(os.path.join(root, d), exist_ok=True)
rng = np.random.RandomState(1)
N, D = 60, 64
cases = [f"TCGA-AA-{i:04d}" for i in range(1, N + 1)]
genes = [f"GENE{i}" for i in range(14)]
lat = rng.standard_normal((N, 3))
rna_expr = rng.standard_normal((N, 14)) * 2 + 5
prot = rna_expr * 0.5 + lat[:, [0]] * 0.8 * (np.arange(14) % 2) + rng.standard_normal((N, 14)) * 0.7
P = pd.DataFrame(prot, index=cases, columns=genes)
P["SPARSE1"] = np.where(rng.rand(N) < 0.35, rng.standard_normal(N), np.nan)       # < 30 non-missing -> must not be tested
P["SPARSE2"] = np.where(rng.rand(N) < 0.30, rng.standard_normal(N), np.nan)
P["NO_RNA_GENE"] = rng.standard_normal(N)                                         # protein only -> must not be tested
P.index = [c + "-01" for c in cases]                                              # case-indexed csv: tcga_case() truncates
P.to_csv(os.path.join(root, "protein_brca.csv"))
rows = []
for i, c in enumerate(cases):
    fn = f"{c}.rna.tsv"
    d = pd.DataFrame({"gene_id": [f"ENSG{j}" for j in range(14 + 3)], "gene_name": genes + ["SPARSE1", "SPARSE2", "N_unmapped"],
                      "tpm_unstranded": list(2 ** rna_expr[i] - 1) + [5.0, 6.0, 7.0]})
    with open(os.path.join(root, "rna", fn), "w") as fh:
        fh.write("# gene-model: GENCODE v36\n"); d.to_csv(fh, sep="\t", index=False)
    rows.append((c, fn))
pd.DataFrame(rows, columns=["case", "filename"]).to_csv(os.path.join(root, "rna_case_file.tsv"), sep="\t", index=False)
open(os.path.join(root, "dx_case_list.tsv"), "w").write("case\n" + "\n".join(cases) + "\n")


def emb(i, ntile, shift=0.0):
    base_v = lat[i] @ np.random.RandomState(5).standard_normal((3, D))
    return torch.tensor(base_v + rng.standard_normal((ntile, D)) * 0.5 + shift, dtype=torch.float32)


for i, c in enumerate(cases):
    torch.save({"embeddings": emb(i, int(rng.randint(5, 12)))}, os.path.join(root, "emb_phikon", f"{c}-01Z-00-DX1.abc.pt"))
    torch.save({"embeddings": emb(i, int(rng.randint(5, 12)))}, os.path.join(root, "emb_phikon_ts", f"{c}-01A-01-TS1.def.pt"))
    if i % 7 == 0:
        torch.save({"embeddings": emb(i, 6)}, os.path.join(root, "emb_phikon_ts", f"{c}-01A-02-BS2.ghi.pt"))   # 2nd tumour slide
    if i % 5 == 0:
        torch.save({"embeddings": emb(i, 6, shift=3.0)}, os.path.join(root, "emb_phikon_ts", f"{c}-11A-01-TS1.jkl.pt"))  # normal -> dropped


def run(script, extra, out, env_extra=None):
    env = dict(os.environ, MORPHO_OUT=out, PYTHONIOENCODING="utf-8"); env.pop("MORPHO_WSI_EMB_DIR", None)
    r = subprocess.run([sys.executable, "-u", script, "--project", "TCGA-BRCA", "--base", base, *extra], env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:]); raise SystemExit(f"FAILED: {script} {extra}")
    return r.stdout


out_pub, out_dx, out_ts, out_tile, out_smoke = (os.path.join(scratch, x) for x in ("pub", "dx", "ts", "tile", "smoke"))
run(PUBLISHED, [], out_pub)
run(PATCHED, [], out_dx)
a = pd.read_csv(os.path.join(out_pub, "brca_protein_results.csv")); b = pd.read_csv(os.path.join(out_dx, "brca_protein_results.csv"))
pd.testing.assert_frame_equal(a, b, check_exact=True)
pd.testing.assert_frame_equal(pd.read_csv(os.path.join(out_pub, "summary_brca.csv")), pd.read_csv(os.path.join(out_dx, "summary_brca.csv")),
                              check_exact=True)
assert len(a) == 14, len(a)
print(f"T1 ok: published vs patched dx-arm csv/summary identical (check_exact), {len(a)} genes tested (2 sparse + 1 RNA-less excluded)")

z = np.load(os.path.join(out_dx, "brca_incr_matrix.npz"), allow_pickle=False)
M = z["M"]; assert M.shape == (14, 1001) and int(z["B"]) == 1000 and np.array_equal(z["perm_seeds"], np.arange(1000))
assert list(z["genes"]) == list(a["gene"]) and np.allclose(M[:, 0], a["incremental_r2"].values, rtol=0, atol=1e-14)
p = (np.sum(M[:, 1:] >= M[:, [0]], axis=1) + 1) / 1001
assert np.allclose(p, a["pval"].values, rtol=0, atol=1e-15) and not bool(z["smoke"]) and str(z["arm"]) == "dx"
assert len(z["cases"]) == N and not any(g.startswith("SPARSE") for g in z["genes"])
print("T2 ok: matrix 14 x 1001, col 0 == csv increments, pval reproduced exactly, seeds 0..999, row order == csv")

run(PATCHED, ["--slide-type", "ts", "--case-list", os.path.join(root, "dx_case_list.tsv")], out_ts)
run(PATCHED, ["--slide-type", "ts", "--pool", "tile", "--case-list", os.path.join(root, "dx_case_list.tsv")], out_tile)
pm = pd.read_csv(os.path.join(out_ts, "pool_manifest_brca_ts.csv"))
assert not pm["slides"].str.contains("-11A-").any() and pm["n_slides"].max() == 2
for o, tag, pool in ((out_ts, "_ts", "slide"), (out_tile, "_ts_tilepool", "tile")):
    zz = np.load(os.path.join(o, f"brca_incr_matrix{tag}.npz"), allow_pickle=False)
    assert str(zz["arm"]) == "ts" and str(zz["pool"]) == pool and zz["M"].shape == (14, 1001)
assert not np.array_equal(np.load(os.path.join(out_ts, "brca_incr_matrix_ts.npz"))["M"],
                          np.load(os.path.join(out_tile, "brca_incr_matrix_ts_tilepool.npz"))["M"])
print("T3 ok: ts arm drops 11A, writes brca_incr_matrix_ts.npz / _ts_tilepool.npz (slide vs tile pooling differ)")

run(PATCHED, ["--perms", "7", "--max-genes", "3"], out_smoke)
sm = [f for f in os.listdir(out_smoke) if "smoke" in f]
assert len(sm) == 3, sm
try:
    L.load_matrix(os.path.join(out_smoke, "brca_incr_matrix_smoke.npz"), "dx")
    raise AssertionError("smoke matrix was accepted")
except SystemExit as e:
    assert "smoke" in str(e.code) or "B=" in str(e.code), e.code
print("T4 ok: --perms/--max-genes outputs are tagged _smoke and refused by the set-level script")

lab = {"dx": os.path.join(out_dx, "brca_incr_matrix.npz"), "ts": os.path.join(out_ts, "brca_incr_matrix_ts.npz"),
       "ts_tile": os.path.join(out_tile, "brca_incr_matrix_ts_tilepool.npz")}
mats = {k: L.load_matrix(v, k) for k, v in lab.items()}
try:
    L.load_matrix(lab["ts"], "dx"); raise AssertionError("label mismatch accepted")
except SystemExit as e:
    assert "expects arm" in str(e.code)
cmp = L.compare_dx_published(mats["dx"], os.path.join(out_pub, "brca_protein_results.csv"))
assert cmp["verdict"] == "MATCH", cmp
print("T5 ok: three matrices load (B=1000); label/arm mismatch refused; DX comparator vs 'published' csv MATCH")

out = subprocess.run([sys.executable, TESTED, "--root", root, "--case-list", os.path.join(root, "dx_case_list.tsv"),
                      "--out", os.path.join(scratch, "tested.tsv"), "--expect", "14"], capture_output=True, text=True)
assert out.returncode == 0, out.stderr
tg = pd.read_csv(os.path.join(scratch, "tested.tsv"), sep="\t")
assert list(tg["gene"]) == list(a["gene"]) and (tg["n"].values == a["n"].values).all() and "DIFFERS" not in out.stdout
print("T6 ok: brca_tested_genes.py == genes of the results csv (labels/counts only)")

# T7: set-level statistic from the saved matrix vs an independent recomputation (scipy), on a synthetic 'selection set'
from scipy.stats import mannwhitneyu, spearmanr
g = list(mats["ts"]["genes"]); sel = np.array([i % 3 == 0 for i in range(len(g))])
score = rng.standard_normal(len(g))
sdf = pd.DataFrame({"gene": g, "selected": sel.astype(int), "score": score})
MIN_G = L.MIN_GROUP; L.MIN_GROUP = 2
r = L.run_set(mats["ts"], sdf); L.MIN_GROUP = MIN_G
Mt = mats["ts"]["M"]
u = mannwhitneyu(Mt[sel, 0], Mt[~sel, 0]).statistic / (sel.sum() * (~sel).sum())
assert abs(r["auc_obs"] - u) < 1e-12 and abs(r["rho_obs"] - spearmanr(score, Mt[:, 0]).correlation) < 1e-12
nulls = [mannwhitneyu(Mt[sel, j], Mt[~sel, j]).statistic / (sel.sum() * (~sel).sum()) for j in range(1, 1001)]
assert abs(r["p_auc"] - (np.sum(np.array(nulls) >= u) + 1) / 1001) < 1e-12
print(f"T7 ok: set-level AUC/rho/p from the saved matrix == scipy recomputation (AUC {r['auc_obs']:.4f}, p {r['p_auc']:.4f}; synthetic)")
print("ALL SYNTHETIC PIPELINE TESTS PASSED; scratch:", scratch)
