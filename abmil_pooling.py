#!/usr/bin/env python3
"""abmil_pooling.py (SERVER, gpu02 or CPU) — aggregation ablation (Tier-2 ⑥): does the
morphology->residual signal depend on mean-pooling? Trains attention-based MIL (ABMIL,
gated attention) on per-tile UNI embeddings to predict the per-patient measured family
residual score, with patient-level 5-fold CV (seed 0), and compares out-of-fold Pearson r /
R^2 against the main-analysis aggregation (mean-pool -> PCA(20) -> ridge lambda=1) on the
SAME folds.  Also saves the attention weights of the best-predicted patients for a figure.

Run on gpu02 via LSF (NOT the login node -- it was OOM-killed there; memory is now bounded but
GPU is much faster):
    bsub -q interactive -gpu "num=1:mode=exclusive_process" -J abmil \
      -o /public/home/fjhui/ZW/results/abmil.out -e /public/home/fjhui/ZW/results/abmil.err \
      bash -c 'cd /public/home/fjhui/ZW/scripts; source morpho_env.sh ccrcc; \
      export MORPHO_COHORT=ccrcc MORPHO_SCORES_CSV=/public/home/fjhui/ZW/results/clinical_link_ccrcc_scores.csv MORPHO_OUT=/public/home/fjhui/ZW/results; \
      /public/home/fjhui/miniconda3/bin/python -u abmil_pooling.py'
Then WinSCP  $MORPHO_OUT/abmil_<cohort>.csv  (+ abmil_attention_<cohort>.npz)  to figures/figdata/ .
"""
import os, glob, numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
import residual_analysis as RA

C = os.environ.get("MORPHO_COHORT", "cohort")
EMB = os.environ.get("MORPHO_WSI_EMB_DIR") or os.environ.get("WSI_EMB_DIR")
SCORES = os.environ["MORPHO_SCORES_CSV"]; OUT = os.environ.get("MORPHO_OUT", ".")
TARGET = os.environ.get("MORPHO_FAMILY", "translation") + "_meas"
MAX_TILES, EPOCHS, LR, WD, SEED = 2000, 40, 1e-4, 1e-4, 0
dev = "cuda" if torch.cuda.is_available() else "cpu"
torch.manual_seed(SEED); np.random.seed(SEED)


def patient_of(stem):                      # CPTAC slide 'C3L-00011-21' -> patient 'C3L-00011'
    return "-".join(stem.split("-")[:2])


# ---- patients (tumour-only set = whatever the main analysis kept) + target ----
wsi_mean = RA.load_wsi_embeddings()          # patients x D (tumour-only, mean-pooled)
sc = pd.read_csv(SCORES, index_col=0)
y_all = pd.to_numeric(sc[TARGET], errors="coerce").reindex(wsi_mean.index)
keep = y_all.notna().values
patients = list(wsi_mean.index[keep]); y = y_all.values[keep].astype(float)
print(f"[abmil] {C}: {len(patients)} patients with {TARGET}; device={dev}")

# ---- per-patient tile bags from raw .pt, MEMORY-BOUNDED: subsample per slide at load time
# (cap MAX_TILES per patient split across its slides) and keep float16; never hold the whole
# cohort's tiles in RAM (that was OOM-killed on the login node). ----
rng = np.random.RandomState(SEED)
pset = set(patients); slides_of = {}
for fp in sorted(glob.glob(f"{EMB}/*.pt")):
    p = patient_of(os.path.splitext(os.path.basename(fp))[0])
    if p in pset:
        slides_of.setdefault(p, []).append(fp)
bags = {}
for p, fps in slides_of.items():
    cap = max(200, MAX_TILES // len(fps)); parts = []
    for fp in fps:
        e = torch.load(fp, map_location="cpu")["embeddings"]
        if e.shape[0] > cap:
            e = e[torch.from_numpy(rng.choice(e.shape[0], cap, replace=False))]
        parts.append(e.to(torch.float16)); del e
    bags[p] = torch.cat(parts, 0)
patients = [p for p in patients if p in bags]; y = np.array([y_all[p] for p in patients], float)
D = bags[patients[0]].shape[1]
print(f"[abmil] bags: {len(patients)} patients, median {int(np.median([b.shape[0] for b in bags.values()]))} tiles, {D}d")


class ABMIL(nn.Module):
    def __init__(self, d, h=256, a=128):
        super().__init__()
        self.f = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(0.1))
        self.att_v = nn.Sequential(nn.Linear(h, a), nn.Tanh())
        self.att_u = nn.Sequential(nn.Linear(h, a), nn.Sigmoid())
        self.att_w = nn.Linear(a, 1)
        self.head = nn.Linear(h, 1)

    def forward(self, x):                   # x: [N, d] one bag
        h = self.f(x)
        a = self.att_w(self.att_v(h) * self.att_u(h)).squeeze(-1)      # [N]
        a = torch.softmax(a, 0)
        return self.head((a[:, None] * h).sum(0)).squeeze(), a


def train_abmil(tr_idx, ymu, ysd):
    m = ABMIL(D).to(dev); opt = torch.optim.Adam(m.parameters(), lr=LR, weight_decay=WD)
    for ep in range(EPOCHS):
        m.train(); order = rng.permutation(tr_idx)
        for i in order:
            x = bags[patients[i]].float().to(dev); t = torch.tensor((y[i] - ymu) / ysd, dtype=torch.float32, device=dev)
            pred, _ = m(x); loss = (pred - t) ** 2
            opt.zero_grad(); loss.backward(); opt.step()
    return m


def r2(a, b):
    return 1 - np.sum((a - b) ** 2) / np.sum((a - a.mean()) ** 2)


kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
oof_ab = np.full(len(patients), np.nan); oof_mp = np.full(len(patients), np.nan); att_store = {}
Xmean = wsi_mean.loc[patients].values.astype(float)   # full-tile mean = EXACTLY the main analysis (not the subsampled bag mean)
for k, (tr, te) in enumerate(kf.split(Xmean)):
    ymu, ysd = y[tr].mean(), y[tr].std() + 1e-8
    # ABMIL
    m = train_abmil(tr, ymu, ysd); m.eval()
    with torch.no_grad():
        for i in te:
            pred, a = m(bags[patients[i]].float().to(dev))
            oof_ab[i] = float(pred) * ysd + ymu; att_store[patients[i]] = a.cpu().numpy()
    # mean-pool -> PCA20 -> ridge (main analysis aggregation), same folds
    pca = PCA(20, svd_solver="full").fit(Xmean[tr]); Ptr, Pte = pca.transform(Xmean[tr]), pca.transform(Xmean[te])
    xm, xs = Ptr.mean(0), Ptr.std(0); xs[xs == 0] = 1; A = (Ptr - xm) / xs
    w = np.linalg.solve(A.T @ A + np.eye(20), A.T @ (y[tr] - ymu)); oof_mp[te] = ((Pte - xm) / xs) @ w + ymu
    print(f"[fold {k+1}] ABMIL r={np.corrcoef(y[te], oof_ab[te])[0,1]:+.3f}  mean-pool r={np.corrcoef(y[te], oof_mp[te])[0,1]:+.3f}")

rows = [dict(cohort=C, method="ABMIL (gated attention)", r=np.corrcoef(y, oof_ab)[0, 1], R2=r2(y, oof_ab), n=len(y)),
        dict(cohort=C, method="mean-pool + PCA20 + ridge (main)", r=np.corrcoef(y, oof_mp)[0, 1], R2=r2(y, oof_mp), n=len(y))]
res = pd.DataFrame(rows); res.to_csv(f"{OUT}/abmil_{C}.csv", index=False)
print("\n" + res.round(3).to_string(index=False))
# attention of the 3 best-predicted patients (for an attention-map figure)
err = np.abs(y - oof_ab); best = [patients[i] for i in np.argsort(err)[:3]]
np.savez(f"{OUT}/abmil_attention_{C}.npz", **{p: att_store[p] for p in best}, patients=np.array(best))
print(f" -> {OUT}/abmil_{C}.csv , abmil_attention_{C}.npz (best-predicted: {best})")
print("READ: if ABMIL ~ mean-pool (or better), the signal does not hinge on the pooling choice.")
