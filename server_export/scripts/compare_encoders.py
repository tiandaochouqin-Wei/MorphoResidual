#!/usr/bin/env python3
"""compare_encoders.py — encoder-robustness table (Tier-1 ③): UNI (pinned) vs Phikon vs
H-optimus-0 (3rd pathology FM) vs ResNet-50 ImageNet (non-pathology weak baseline).
Per cohort and encoder: n significant, median incremental R^2 of the significant set; and
versus UNI: Jaccard / recall of the UNI-significant set / per-gene incremental-R^2 correlation.
Expectation: pathology FMs agree with UNI (high Jaccard/r); ImageNet is markedly weaker ->
the residual signal needs a pathology-specific representation (or is robust across FMs).
Run on the server after each encoder's residual_analysis has produced
results_<enc>/<cohort>/residual_results_tumoronly.csv .  Missing encoders are skipped.
    python compare_encoders.py            -> /public/home/fjhui/ZW/results/encoder_comparison.csv
"""
import os, pandas as pd

COH = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
ENC = {"UNI": "/public/home/fjhui/ZW/scripts/pinned/results_snapshot_20260901/{c}/residual_results_tumoronly.csv",
       "Phikon": "/public/home/fjhui/ZW/results_phikon/{c}/residual_results_tumoronly.csv",
       "Phikon-v2": "/public/home/fjhui/ZW/results_phikonv2/{c}/residual_results_tumoronly.csv",
       "Hibou-B": "/public/home/fjhui/ZW/results_hibou/{c}/residual_results_tumoronly.csv",
       "H-optimus-0": "/public/home/fjhui/ZW/results_hoptimus0/{c}/residual_results_tumoronly.csv",  # gated; skipped if absent
       "ResNet50-ImageNet": "/public/home/fjhui/ZW/results_resnet50/{c}/residual_results_tumoronly.csv"}
OUT = os.environ.get("MORPHO_OUT", "/public/home/fjhui/ZW/results")


def sig_set(df):
    return set(df[(df.fdr < 0.05) & (df.incremental_r2 > 0)].gene.astype(str))


rows = []
for c in COH:
    ref = None
    for e, pat in ENC.items():
        p = pat.format(c=c)
        if not os.path.exists(p):
            print(f"  [{c}] {e}: missing ({p})"); continue
        d = pd.read_csv(p); s = sig_set(d)
        row = dict(cohort=c, encoder=e, n_sig=len(s),
                   med_incr_sig=float(d[d.gene.isin(s)].incremental_r2.median()) if s else float("nan"))
        if e == "UNI":
            ref = (d, s)
        elif ref is not None:
            u, su = ref; inter, union = su & s, su | s
            m = u[["gene", "incremental_r2"]].merge(d[["gene", "incremental_r2"]], on="gene", suffixes=("_uni", "_e"))
            row.update(jaccard_vs_uni=len(inter) / len(union) if union else float("nan"),
                       recall_of_uni_sig=len(inter) / len(su) if su else float("nan"),
                       r_incr_vs_uni=float(m.incremental_r2_uni.corr(m.incremental_r2_e)))
        rows.append(row)
res = pd.DataFrame(rows)
os.makedirs(OUT, exist_ok=True); res.to_csv(f"{OUT}/encoder_comparison.csv", index=False)
pd.set_option("display.width", 160)
print(res.round(3).to_string(index=False))
print(f"\n -> {OUT}/encoder_comparison.csv")
print("READ: pathology FMs should show high jaccard/recall/r vs UNI; ResNet50-ImageNet should be clearly lower.")
