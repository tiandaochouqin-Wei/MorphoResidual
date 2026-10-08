#!/usr/bin/env python3
"""
End-to-end plumbing test of measurement_error_checks.py on a SYNTHETIC cohort (random
numbers made up here; no project data).  Builds fake STAR-Counts files, a TMT table, a
crosswalk, a slide map and .pt embeddings in a temp dir, points the REAL published
loaders (server_export/scripts/residual_analysis.py) at them, runs the script, and checks
  1. the raw-count re-read reproduces RA's RNA matrix,
  2. oof_ridge == RA.cv_r2  (also asserted inside the script),
  3. the observed incremental R^2 == RA._process_gene's,
  4. the first B permutation statistics reproduce RA's null increments (scheme 'ra'),
  5. the sufficient statistics rebuild R^2 from the saved OOF vectors,
  6. --resume works after an interruption.
Usage:  python test_synthetic.py <tmpdir> [path-to-residual_analysis-dir]
"""
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
tmp.mkdir(parents=True, exist_ok=True)

rs = np.random.RandomState(7)
N, G, B = 100, 150, 4
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

env = dict(os.environ, MORPHO_ROOT=str(tmp), MORPHO_RNA_MANIFEST=str(man),
           MORPHO_ALIQUOT_XWALK=str(tmp / "xwalk.tsv"), MORPHO_SLIDE_MAP=str(tmp / "slides.tsv"),
           MORPHO_PROTEIN_TSV=str(tmp / "omics" / "protein" / "fake.tmt10.tsv"),
           MORPHO_WSI_EMB_DIR=str(emb_dir), MORPHO_OUT=str(tmp / "out_ra"),
           MORPHO_SCRIPTS_DIR=str(ra_dir), PYTHONIOENCODING="utf-8")
env.pop("MORPHO_SIG_CSV", None)
out = tmp / "out"
cmd = [sys.executable, "-u", str(here / "measurement_error_checks.py"), "--cohort", "ucec",
       "--out", str(out), "--B", str(B), "--workers", "2"]
r = subprocess.run(cmd, env=env, capture_output=True, text=True)
print(r.stdout[-3500:])
print(r.stderr[-1500:])
assert r.returncode == 0, "script failed"

# ---- independent check against the published code path ----
sys.path.insert(0, str(ra_dir))
os.environ.update(env)
import residual_analysis as RA
from sklearn.decomposition import PCA

rna, protein, wsi = RA.load_rna_matrix(), RA.load_protein_matrix(), RA.load_wsi_embeddings()
common = sorted(set(rna.index) & set(protein.index) & set(wsi.index))
cg = sorted(set(rna.columns) & set(protein.columns))
rna, protein = rna.loc[common, cg], protein.loc[common, cg]
pcs = PCA(n_components=20, svd_solver="full").fit_transform(wsi.loc[common].values)
rng = np.random.RandomState(0)
perms = [rng.permutation(len(common)) for _ in range(B)]

tab = pd.read_csv(out / "gene_table_ucec.csv")
inp = np.load(out / "inputs_ucec.npz", allow_pickle=True)
oof = np.load(out / "oof_pred_ucec.npy")
P = np.load(out / "perm_stats_ucec.npy")
fit_genes = list(inp["genes"])
worst_obs = worst_null = worst_rebuild = 0.0
for gi in range(0, len(fit_genes), 7):
    g = fit_genes[gi]
    y = protein[g].values
    v = ~np.isnan(y)
    yv, xr = y[v], rna[g].values[v].reshape(-1, 1)
    r2_rna = RA.cv_r2(xr, yv, 5)
    r2_both = RA.cv_r2(np.hstack([xr, pcs[v]]), yv, 5)
    row = tab[tab.gene == g].iloc[0]
    worst_obs = max(worst_obs, abs(row["incr_p"] - (r2_both - r2_rna)), abs(row["r2_rna_p"] - r2_rna))
    for b in range(B):
        ra_null = RA.cv_r2(np.hstack([xr, pcs[perms[b]][v]]), yv, 5) - r2_rna
        sst = row["sst_p"]
        mine = (2 * P[gi, b, 5] - P[gi, b, 2]) / sst
        worst_null = max(worst_null, abs(ra_null - mine))
    # rebuild from saved OOF vectors
    p0, p1 = oof[gi, 0][v], oof[gi, 1][v]
    sst = np.sum((yv - yv.mean()) ** 2)
    worst_rebuild = max(worst_rebuild, abs((np.sum((yv - p0) ** 2) - np.sum((yv - p1) ** 2)) / sst - row["incr_p"]))
print(f"max|incr_p - RA|={worst_obs:.2e}   max|null - RA null|={worst_null:.2e} (float32 storage)   "
      f"max|rebuild from OOF|={worst_rebuild:.2e}")
assert worst_obs < 1e-10 and worst_null < 1e-5 and worst_rebuild < 1e-10
print("columns:", list(tab.columns))
print(tab.head(3).T)
print("replicate ICC file:", (out / "replicate_icc_ucec.csv").exists())

# ---- resume: drop the last 30 genes from the done mask, rerun with --resume ----
done = np.load(out / "done_ucec.npy")
done[-30:] = False
np.save(out / "done_ucec.npy", done)
before = np.load(out / "obs_scalars_ucec.npy").copy()
r2 = subprocess.run(cmd + ["--resume"], env=env, capture_output=True, text=True)
assert r2.returncode == 0, r2.stdout[-2000:] + r2.stderr[-2000:]
after = np.load(out / "obs_scalars_ucec.npy")
assert np.allclose(before, after, equal_nan=True)
print("resume OK; files:", sorted(p.name for p in out.iterdir()))
print("ALL CHECKS PASSED")
