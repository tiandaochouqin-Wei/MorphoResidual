#!/usr/bin/env python3
"""
clinical_cox.py — multivariable survival for MorphoResidual (Batch A-1, Tier 1).

Defends against the "morphology trivially predicts grade -> circular" objection:
does the H&E-only residual pathway score add INDEPENDENT prognostic value beyond
grade / stage / age?  Reads the scores + clinical CSVs already written by
clinical_link.py (no data reload).  Cox PH fitted from scratch in numpy
(Breslow ties, Newton-Raphson, small L2 for stability) -- no lifelines needed.

For each morphology-only score and each cohort:
    uni   : Cox(survival ~ score)
    multi : Cox(survival ~ score + grade + stage + age)   # score HR adjusted
Reports hazard ratio, 95% CI and p for the score.

Run (reads the CSVs already on the server; light, any node):
    for c in ccrcc luad pdac gbm ucec; do
      MORPHO_CLIN_OUT=/public/home/fjhui/ZW/results/clinical_link_$c \
      python -u clinical_cox.py; done
"""
import os, sys, numpy as np, pandas as pd
from scipy import stats


def cox_fit(X, time, event, l2=1e-3, max_iter=100, tol=1e-9):
    """Breslow-tie Cox PH by Newton-Raphson. Returns beta, se, p (ridge-stabilised)."""
    X = np.asarray(X, float); time = np.asarray(time, float); event = np.asarray(event, int)
    n, p = X.shape
    beta = np.zeros(p)
    ev_times = np.unique(time[event == 1])
    risk_masks = [(time >= t) for t in ev_times]
    death_masks = [((time == t) & (event == 1)) for t in ev_times]
    d_counts = [int(m.sum()) for m in death_masks]
    Xdeath_sum = [X[m].sum(0) for m in death_masks]
    for _ in range(max_iter):
        eta = X @ beta
        eta -= eta.max()
        w = np.exp(eta)
        grad = -l2 * beta
        H = l2 * np.eye(p)
        for rm, d, xds in zip(risk_masks, d_counts, Xdeath_sum):
            Xr = X[rm]; wr = w[rm]
            S0 = wr.sum()
            if S0 <= 0:
                continue
            S1 = (Xr * wr[:, None]).sum(0)
            xbar = S1 / S0
            S2 = (Xr[:, :, None] * Xr[:, None, :] * wr[:, None, None]).sum(0)
            grad += xds - d * xbar
            H += d * (S2 / S0 - np.outer(xbar, xbar))
        try:
            step = np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            break
        beta += step
        if np.max(np.abs(step)) < tol:
            break
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)
    z = np.divide(beta, se, out=np.full(p, np.nan), where=se > 0)
    pval = 2 * stats.norm.sf(np.abs(z))
    return beta, se, pval


def zstd(a):
    a = np.asarray(a, float)
    s = np.nanstd(a)
    return (a - np.nanmean(a)) / (s if s > 0 else 1.0)


def main():
    prefix = os.environ.get("MORPHO_CLIN_OUT")
    if not prefix:
        sys.exit("set MORPHO_CLIN_OUT to the clinical_link_<cohort> prefix")
    cohort = os.path.basename(prefix).replace("clinical_link_", "")
    scores = pd.read_csv(f"{prefix}_scores.csv", index_col=0)
    clin = pd.read_csv(f"{prefix}_clinical.csv", index_col=0)
    df = scores.join(clin[[c for c in ["grade", "stage", "age", "time", "event"] if c in clin.columns]])
    df = df[df["time"].notna() & (df["time"] >= 0)]

    morph_cols = [c for c in scores.columns if c.endswith("_morph")]
    n_death = int(df["event"].sum())
    print(f"\n===== {cohort}: n={len(df)}, deaths={n_death} =====")
    if n_death < 8:
        print(" too few deaths for Cox; skipping")
        return

    # covariate availability
    adj = [c for c in ["grade", "stage", "age"] if c in df.columns and df[c].notna().sum() >= 0.6 * len(df)]
    print(f" adjusting for: {adj if adj else '(none available)'}")

    rows = []
    for col in morph_cols:
        base = df[[col, "time", "event"] + adj].copy()
        base = base.dropna()
        if base["event"].sum() < 8 or len(base) < 20:
            rows.append({"score": col, "n": len(base), "note": "insufficient"})
            continue
        t = base["time"].values; e = base["event"].values
        # univariable
        Xu = zstd(base[col].values).reshape(-1, 1)
        bu, seu, pu = cox_fit(Xu, t, e)
        # multivariable (score + adjusters), all z-standardised
        cols = [col] + adj
        Xm = np.column_stack([zstd(base[c].values) for c in cols])
        bm, sem, pm = cox_fit(Xm, t, e)
        rows.append({
            "score": col, "n": len(base),
            "uni_HR": float(np.exp(bu[0])), "uni_p": float(pu[0]),
            "adj_HR": float(np.exp(bm[0])), "adj_p": float(pm[0]),
            "adj_for": "+".join(adj) if adj else "-",
        })
    res = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    out = f"{prefix}_cox.csv"
    res.to_csv(out, index=False)
    print(f" -> {out}")
    print(" READ: adj_HR/adj_p = score's hazard ratio AFTER adjusting for grade/stage/age;")
    print("       adj_p<0.05 => the H&E residual score carries prognostic info beyond grade (not circular).")


if __name__ == "__main__":
    main()
