#!/usr/bin/env python3
r"""wsi_spatial_map_oof.py (SERVER, compute-only, NO matplotlib) -- T1-8.

Out-of-fold rescoring of the whole-slide residual maps.

WHAT WAS WRONG WITH wsi_spatial_map.py
Two circularities, and the second is the worse one:

  1. The morphology->residual direction w was fitted on ALL patients, including the
     two whose slides are then painted with it. Every tile of those slides is scored
     by a model that has already seen that patient's own label.

  2. The two slides were then CHOSEN as the argmax and argmin of that same in-sample
     score. Selecting the extremes of an in-sample fit picks precisely the points
     where the fit is most optimistic, so the printed contrast between the two maps
     is an upper bound on what the direction can do, not an estimate of it.

Because of this the manuscript describes the maps at patient level only and the
figure sits in the Supplement, with the roadmap recording an out-of-fold rescoring
as the condition for saying anything more. This script is that rescoring.

WHAT IT DOES INSTEAD
The same 5-fold KFold(shuffle=True, random_state=0) the rest of the study uses.
PCA(20) and the ridge direction are refitted inside each training split, so a
patient's tiles are always scored by a model fitted without that patient. Slides
are then ranked on those out-of-fold scores, and the extremes are chosen from that
ranking. Both fixes matter separately, so the script reports the in-sample answer
alongside and says explicitly whether the same slides survive selection.

It also redoes the variance decomposition out-of-fold. The claim the figure rests
on -- that the signal is largely between patients rather than within a slide,
roughly a 10% within-slide share -- was computed in-sample, and an in-sample
direction inflates between-patient variance by construction. If the within-slide
share moves materially out-of-fold, the figure's caption has to move with it.

USAGE (one cohort at a time; ccRCC shown, the other four differ only in the env)
    cd /public/home/fjhui/ZW/scripts && source morpho_env.sh ccrcc
    export MORPHO_SCORES_CSV=/public/home/fjhui/ZW/results/clinical_link_ccrcc_scores.csv
    export MORPHO_OUT=/public/home/fjhui/ZW/results
    /public/home/fjhui/miniconda3/bin/python -u wsi_spatial_map_oof.py

    Optional: MORPHO_FAMILY (default translation), MORPHO_FOLDS (default 5).
    CPU only, one .pt held in memory at a time. Runtime is dominated by reading
    every tile file once; expect minutes, not hours.

    Then WinSCP  $MORPHO_OUT/wsi_spatial_data_oof.npz  back to figures/figdata/ .

NOTE ON MORPHO_ROOT: morpho_env.sh sets it to the CURRENT cohort's directory. This
script never reads it -- it takes the embedding directory from the same variable
residual_analysis.py uses, so it cannot drift to another cohort the way
export_package.py once did.
"""
import os
import re
import sys
import glob

import numpy as np
import pandas as pd
import torch
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

import residual_analysis as RA

RIDGE_ALPHA = getattr(RA, "RIDGE_ALPHA", 1.0)
N_PC = 20
FOLDS = int(os.environ.get("MORPHO_FOLDS", "5"))
EMB = os.environ.get("MORPHO_WSI_EMB_DIR") or os.environ.get("WSI_EMB_DIR")
SCORES = os.environ.get("MORPHO_SCORES_CSV", "clinical_link_ccrcc_scores.csv")
FAMILY = os.environ.get("MORPHO_FAMILY", "translation")
OUT = os.environ.get("MORPHO_OUT", ".")
if not EMB:
    sys.exit("set MORPHO_WSI_EMB_DIR (or WSI_EMB_DIR) to the raw per-slide .pt directory")

SLIDE_RE = re.compile(r"^(C3[A-Z]-\d{5}-\d+)")
CASE_RE = re.compile(r"^(C3[A-Z]-\d{5})")


def load_slide(fp):
    d = torch.load(fp, map_location="cpu")
    return d["embeddings"].float().numpy(), np.asarray(d["coords"], float)


def fit_one(X, y):
    """PCA(20) + standardise + closed-form ridge, exactly as the main analysis."""
    pca = PCA(n_components=N_PC, svd_solver="full").fit(X)
    P = pca.transform(X)
    xm, xs = P.mean(0), P.std(0)
    xs[xs == 0] = 1.0
    ym = float(y.mean())
    Ps = (P - xm) / xs
    w = np.linalg.solve(Ps.T @ Ps + RIDGE_ALPHA * np.eye(N_PC), Ps.T @ (y - ym))
    return (pca, xm, xs, ym, w)


def score_tiles(emb, model):
    pca, xm, xs, ym, w = model
    pc = (emb - pca.mean_) @ pca.components_.T
    return ((pc - xm) / xs) @ w + ym


def build_models():
    """One model per fold; every patient gets the model that did NOT see them."""
    wsi = RA.load_wsi_embeddings()
    sc = pd.read_csv(SCORES, index_col=0)
    meas = pd.to_numeric(sc[f"{FAMILY}_meas"], errors="coerce").reindex(wsi.index).values
    keep = np.isfinite(meas)
    wsi, meas = wsi[keep], meas[keep]
    X, y = wsi.values, meas
    n = len(wsi)
    if n < FOLDS * 4:
        sys.exit(f"only {n} patients with a {FAMILY} score -- too few for {FOLDS} folds")

    insample = fit_one(X, y)
    per_case, fold_of = {}, {}
    kf = KFold(n_splits=FOLDS, shuffle=True, random_state=0)
    for k, (tr, te) in enumerate(kf.split(X)):
        m = fit_one(X[tr], y[tr])
        for i in te:
            per_case[wsi.index[i]] = m
            fold_of[wsi.index[i]] = k
    print(f"[fit] {n} patients, {FOLDS} folds; a patient is never scored by a model "
          f"that saw them")

    # What is the direction actually worth out of fold? Comparing the coefficient
    # vectors across folds would be meaningless -- each fold refits its own PCA and
    # component signs are arbitrary -- so compare the PREDICTIONS, on the same
    # per-patient mean embeddings the main analysis uses. This is the estimand the
    # maps illustrate, and it belongs in the caption.
    pred_o = np.array([score_tiles(X[i:i + 1], per_case[wsi.index[i]])[0] for i in range(n)])
    pred_i = np.array([score_tiles(X[i:i + 1], insample)[0] for i in range(n)])
    r_o = float(np.corrcoef(pred_o, y)[0, 1])
    r_i = float(np.corrcoef(pred_i, y)[0, 1])
    r2_o = 1.0 - float(((y - pred_o) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    r2_i = 1.0 - float(((y - pred_i) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    print(f"      per-patient direction: in-sample r = {r_i:+.3f} (R2 {r2_i:+.3f}), "
          f"out-of-fold r = {r_o:+.3f} (R2 {r2_o:+.3f})")
    return per_case, fold_of, insample, wsi.index, dict(r_oof=r_o, r_ins=r_i,
                                                        r2_oof=r2_o, r2_ins=r2_i)


def main():
    per_case, fold_of, insample, cases, perf = build_models()
    case_set = set(cases)

    # per-slide streaming statistics under both models
    files = sorted(glob.glob(f"{EMB}/*.pt"))
    rows = []
    acc = {}          # case -> [n, sum, sumsq] of OOF tile scores
    for fp in files:
        stem = os.path.basename(fp)[:-3]
        cm = CASE_RE.match(stem)
        if not cm or cm.group(1) not in case_set:
            continue
        case = cm.group(1)
        try:
            emb, _ = load_slide(fp)
        except Exception as e:                                    # noqa: BLE001
            print(f"  !! {stem}: {str(e)[:80]}")
            continue
        ts_o = score_tiles(emb, per_case[case])
        ts_i = score_tiles(emb, insample)
        a = acc.setdefault(case, [0, 0.0, 0.0])
        a[0] += len(ts_o); a[1] += float(ts_o.sum()); a[2] += float((ts_o ** 2).sum())
        rows.append(dict(file=fp, slide=stem, case=case, fold=fold_of[case],
                         n_tiles=len(ts_o),
                         oof_mean=float(ts_o.mean()), ins_mean=float(ts_i.mean()),
                         oof_sd=float(ts_o.std()), ins_sd=float(ts_i.std())))
        if len(rows) % 100 == 0:
            print(f"  scored {len(rows)} slides", flush=True)
    if not rows:
        sys.exit("no slides matched the analysis cases -- check MORPHO_WSI_EMB_DIR")
    df = pd.DataFrame(rows)
    print(f"[scored] {len(df)} slides, {int(df.n_tiles.sum())} tiles, "
          f"{df.case.nunique()} cases")

    r = np.corrcoef(df.oof_mean, df.ins_mean)[0, 1]
    shrink = df.oof_mean.std() / df.ins_mean.std()
    print(f"\n[in-sample vs out-of-fold, per slide]")
    print(f"  correlation of slide means      r = {r:+.3f}")
    print(f"  spread of slide means           OOF/in-sample sd = {shrink:.3f}")

    # ---- selection: the extremes, chosen honestly and chosen the old way
    hi_o, lo_o = df.loc[df.oof_mean.idxmax()], df.loc[df.oof_mean.idxmin()]
    hi_i, lo_i = df.loc[df.ins_mean.idxmax()], df.loc[df.ins_mean.idxmin()]
    print(f"\n[selection]")
    print(f"  in-sample picks : high {hi_i.slide} ({hi_i.ins_mean:+.2f})  "
          f"low {lo_i.slide} ({lo_i.ins_mean:+.2f})")
    print(f"  out-of-fold picks: high {hi_o.slide} ({hi_o.oof_mean:+.2f})  "
          f"low {lo_o.slide} ({lo_o.oof_mean:+.2f})")
    same = (hi_o.slide == hi_i.slide) and (lo_o.slide == lo_i.slide)
    print(f"  same slides survive selection: {'YES' if same else 'NO'}")
    print(f"  the old picks, rescored out-of-fold: high {hi_i.slide} "
          f"{df.set_index('slide').loc[hi_i.slide, 'oof_mean']:+.2f}, "
          f"low {lo_i.slide} "
          f"{df.set_index('slide').loc[lo_i.slide, 'oof_mean']:+.2f}  "
          f"(this is the honest version of the printed contrast)")

    # ---- variance decomposition, out-of-fold
    ns = np.array([a[0] for a in acc.values()], float)
    ss = np.array([a[1] for a in acc.values()], float)
    sq = np.array([a[2] for a in acc.values()], float)
    mu = ss / ns
    grand = ss.sum() / ns.sum()
    within = float((sq - ss ** 2 / ns).sum())
    between = float((ns * (mu - grand) ** 2).sum())
    tot = within + between
    print(f"\n[variance of tile scores, out-of-fold]")
    print(f"  within patient  {within/tot:6.1%}")
    print(f"  between patient {between/tot:6.1%}")
    print(f"  (the figure's caption currently cites a within-slide share of ~10%, "
          f"computed in-sample)")

    # ---- write the two OOF-selected slides for plotting, same keys as before
    data, allv, files_out = {}, [], []
    for i, fp in enumerate([hi_o.file, lo_o.file]):
        emb, coords = load_slide(fp)
        case = CASE_RE.match(os.path.basename(fp)[:-3]).group(1)
        ts = score_tiles(emb, per_case[case])
        data[f"coords{i}"] = coords.astype(np.float32)
        data[f"scores{i}"] = ts.astype(np.float32)
        files_out.append(os.path.basename(fp))
        allv.append(ts)
    allv = np.concatenate(allv)
    vc = float(np.median(allv))
    data["files"] = np.array(files_out)
    data["family"] = np.array([FAMILY])
    data["vc"] = np.array([vc])
    data["vspan"] = np.array([float(np.percentile(np.abs(allv - vc), 96))])
    data["oof"] = np.array([True])
    for k, v in perf.items():
        data[k] = np.array([v])
    data["same_as_insample"] = np.array([same])
    data["within_share"] = np.array([within / tot])
    data["slide_mean_r"] = np.array([r])
    data["n_slides"] = np.array([len(df)])
    path = f"{OUT}/wsi_spatial_data_oof.npz"
    np.savez(path, **data)
    df.to_csv(f"{OUT}/wsi_spatial_slide_scores_oof.csv", index=False)
    print(f"\n[out] {path}")
    print(f"[out] {OUT}/wsi_spatial_slide_scores_oof.csv  ({len(df)} slides)")
    print("      WinSCP both to figures/figdata/ .")


if __name__ == "__main__":
    main()
