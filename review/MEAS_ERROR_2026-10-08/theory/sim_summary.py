"""sim_summary.py -- compact tables from sim_grid_*.csv (synthetic results only).
Usage: python sim_summary.py sim_grid_n100_G200.csv [more.csv ...]"""
import sys
import numpy as np
import pandas as pd

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 400)

for path in sys.argv[1:]:
    d = pd.read_csv(path)
    d["fX"] = d["fX"].fillna(-1)
    print(f"\n=== {path}: n={d.n.iloc[0]} K={d.K.iloc[0]} G={d.G.iloc[0]} B={d.B_null.iloc[0]} ===")
    cols = ["hyp", "R2rna", "lam", "delta_target", "fX", "feasible", "n_sel", "incr_p_mean_sel",
            "theta_pop", "theta_set_sel", "theta_set_sel_jkse", "lamhat_noise_sel",
            "N_set_z", "Nc_set_z", "N_corr_pg_rej_sel", "S_cross_agg_z", "S_cross_mean_sel", "S_cross_null_mean_sel",
            "S_cross_pg_power_sel", "N_pg_rej_sel", "S_pp_null_mean_sel", "S_dd_agg_z", "S_drm_mean_sel",
            "ratio_m_over_p_sel"]
    cols = [c for c in cols if c in d.columns]
    t = d[cols].copy()
    t = t.rename(columns={"delta_target": "dP", "incr_p_mean_sel": "oofInc_sel", "theta_set_sel_jkse": "th_se",
                          "lamhat_noise_sel": "lamhat", "S_cross_agg_z": "S_z", "S_cross_mean_sel": "S_sel",
                          "S_cross_null_mean_sel": "S_null", "S_cross_pg_power_sel": "S_pgpow",
                          "N_pg_rej_sel": "N_pgrej", "Nc_set_z": "Nc_z", "N_corr_pg_rej_sel": "Nc_pgrej", "S_pp_null_mean_sel": "Spp_null", "S_dd_agg_z": "Sdd_z",
                          "S_drm_mean_sel": "Sdrm_sel", "ratio_m_over_p_sel": "dM/dP"})
    print(t.round(3).to_string(index=False))
    # headline aggregates over feasible, selected>=10 scenarios
    ok = d[(d.n_sel >= 10)]
    print("\n-- feasible-or-not, n_sel>=10: per-hypothesis summary of set-level statistics --")
    for hyp, g in ok.groupby("hyp"):
        print(f"{hyp:5s}  scenarios={len(g):2d}  N_set_z: min={g.N_set_z.min():+6.1f} med={g.N_set_z.median():+6.1f} max={g.N_set_z.max():+6.1f} | "
              f"S_z: min={g.S_cross_agg_z.min():+6.1f} med={g.S_cross_agg_z.median():+6.1f} max={g.S_cross_agg_z.max():+6.1f} | "
              f"theta_set_sel: min={g.theta_set_sel.min():+5.2f} med={g.theta_set_sel.median():+5.2f} max={g.theta_set_sel.max():+5.2f} | "
              f"per-gene S power med={g.S_cross_pg_power_sel.median():.2f}  per-gene N rej med={g.N_pg_rej_sel.median():.2f}")
    print("\n-- S_cross null mean (selected genes), S_pp null mean, S_drm observed: mechanical offsets --")
    print(f"S_cross null mean: {ok.S_cross_null_mean_sel.mean():+.4f} (sd over scenarios {ok.S_cross_null_mean_sel.std():.4f})")
    print(f"S_pp    null mean: {ok.S_pp_null_mean_sel.mean():+.4f}")
    print(f"S_drm   observed : {ok.S_drm_mean_sel.mean():+.4f} (sd {ok.S_drm_mean_sel.std():.4f}) -- same under every hypothesis")
    print(f"N per-gene rejection rate under H_post (should be ~0.05): cov form {ok[ok.hyp=='post'].N_pg_rej_sel.mean():.3f}"
          + (f" | corr form {ok[ok.hyp=='post'].N_corr_pg_rej_sel.mean():.3f}" if "N_corr_pg_rej_sel" in ok.columns else ""))
    if "Nc_set_z" in ok.columns:
        for hyp, g in ok.groupby("hyp"):
            print(f"{hyp:5s} Nc_set_z: min={g.Nc_set_z.min():+6.1f} med={g.Nc_set_z.median():+6.1f} max={g.Nc_set_z.max():+6.1f}  "
                  f"Nc_pgrej med={g.N_corr_pg_rej_sel.median():.2f}")
