#!/usr/bin/env python3
"""Two follow-ups to D4 item 1, specified in review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md
section 7.4 BEFORE this script was run, and reported whatever they show.

Pure arithmetic on the already-saved c1_ucec_completeness.npz. No permutation is
drawn here: the 1,001 columns of M (column 0 observed, 1..1000 the frozen
crc32(b"c1_ucec|perm") null) are reused exactly as c1_ucec_completeness_rerun.py
wrote them. auc_strat is IMPORTED from that script, not reimplemented, so the
estimator cannot drift.

1. Confirmatory-completeness stratification. Pre-fixed rule, by exact analogy with
   the discovery rule (which splits at the discovery cohort size, 100): "complete"
   iff confirm_n == 138, the confirmatory cohort size. Reported: that split, the
   joint 2x2 with discovery completeness, and a split at every distinct value.

2. A paired null for the level drop. D = AUC_pooled - AUC_strat, computed column by
   column so the null distribution of D is paired. If D observed sits inside the
   null spread of D, the 0.597 -> 0.561 fall is the structural removal of
   cross-completeness pairs and nothing more.

Run (locally, on a downloaded copy of the npz -- no HPC, no cluster job):
  python c1_ucec_completeness_followup.py --npz <path to c1_ucec_completeness.npz>
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# The HPC keeps split_replication_power.py beside these scripts; this repo keeps it at
# the top level. Put both on the path so the same file runs in either place, exactly as
# c1_run_test_luad.py does.
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))
import c1_ucec_completeness_rerun as CR  # noqa: E402  (imported for auc_strat, unchanged)

N_CONFIRMATORY_CASES = 138  # the confirmatory cohort size; the pre-fixed split point
N_DISCOVERY_CASES = 100     # what the frozen discovery split used


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def emp_p_two_sided(obs, null):
    """Two-sided empirical p by distance from the null mean, with the usual +1."""
    c = abs(obs - null.mean())
    return float((np.sum(np.abs(null - null.mean()) >= c) + 1) / (len(null) + 1))


def describe(label, auc_all, log):
    obs, null = auc_all[0], auc_all[1:]
    z = (obs - null.mean()) / null.std()
    p = float((np.sum(null >= obs) + 1) / (len(null) + 1))
    log(f"  {label:<34} AUC {obs:.6f}  null {null.mean():.6f} +/- {null.std():.6f}"
        f"  z {z:5.2f}  one-sided p {p:.6f}")
    return dict(auc_obs=float(obs), null_mean=float(null.mean()),
                null_sd=float(null.std()), z=float(z), p_one_sided=p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out-json", default=str(HERE / "pinned" / "c1_ucec_completeness_followup.json"))
    args = ap.parse_args()

    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    z = np.load(args.npz, allow_pickle=True)
    genes, sel_mask = z["genes"], z["sel_mask"].astype(bool)
    strata_disc, M, confirm_n = z["strata"], z["M"], z["confirm_n"].astype(int)
    log(f"npz {args.npz}")
    log(f"  sha256 {sha256(args.npz)}")
    log(f"  genes {len(genes)}  selected {sel_mask.sum()}  background {(~sel_mask).sum()}"
        f"  M {M.shape}")
    if M.shape[1] != 1001:
        sys.exit(f"FATAL: M has {M.shape[1]} columns, expected 1001 (observed + B=1000).")

    out = dict(npz=args.npz, sha256_npz=sha256(args.npz),
               n_genes=int(len(genes)), n_selected=int(sel_mask.sum()),
               n_background=int((~sel_mask).sum()))

    # ---- reproduce the two reference statistics from the same M --------------
    log("\n=== reference (must match the HPC run: pooled 0.597443 p 0.003996; "
        "strat 0.560924 p 0.006993) ===")
    pooled_all, _ = CR.auc_strat(M, sel_mask, np.array(["all"] * len(genes)))
    out["pooled"] = describe("pooled (single stratum)", pooled_all, log)
    strat_all, strat_info = CR.auc_strat(M, sel_mask, strata_disc)
    out["strat_discovery"] = describe("stratified by discovery n", strat_all, log)

    # ---- follow-up 1: confirmatory completeness ------------------------------
    log(f"\n=== follow-up 1: confirmatory completeness "
        f"(pre-fixed split: complete iff confirm_n == {N_CONFIRMATORY_CASES}) ===")
    log(f"  confirm_n range {confirm_n.min()}-{confirm_n.max()}, "
        f"{len(np.unique(confirm_n))} distinct values, "
        f"{(confirm_n == N_CONFIRMATORY_CASES).sum()} genes at {N_CONFIRMATORY_CASES}")
    strata_conf = np.where(confirm_n == N_CONFIRMATORY_CASES, "complete", "incomplete")
    conf_all, conf_info = CR.auc_strat(M, sel_mask, strata_conf)
    out["strat_confirmatory"] = describe("stratified by confirmatory n", conf_all, log)

    strata_joint = np.array([f"{a}|{b}" for a, b in zip(strata_disc, strata_conf)])
    joint_all, joint_info = CR.auc_strat(M, sel_mask, strata_joint)
    out["strat_joint"] = describe("stratified by both (2x2)", joint_all, log)

    strata_exact = np.array([f"c{v}" for v in confirm_n])
    exact_all, exact_info = CR.auc_strat(M, sel_mask, strata_exact)
    out["strat_confirmatory_exact"] = describe(
        f"stratified by exact confirm_n ({len(np.unique(confirm_n))} levels)", exact_all, log)

    log("\n  per-stratum detail (observed only; nulls are in the JSON):")
    per = {}
    for name, info, lab_strata in (("confirmatory", conf_info, strata_conf),
                                    ("joint", joint_info, strata_joint)):
        per[name] = {}
        for lab in sorted(info):
            d = info[lab]
            m = lab_strata == lab
            entry = dict(n_selected=int(d["n1"]), n_background=int(d["n0"]),
                          pair_weight=float(d["n1"] * d["n0"]),
                          frac_positive_selected=float((M[m & sel_mask, 0] > 0).mean())
                          if (m & sel_mask).sum() else None,
                          frac_positive_background=float((M[m & ~sel_mask, 0] > 0).mean())
                          if (m & ~sel_mask).sum() else None)
            if d["auc"] is not None:
                a = d["auc"]
                entry.update(auc_obs=float(a[0]), null_mean=float(a[1:].mean()),
                              null_sd=float(a[1:].std()),
                              z=float((a[0] - a[1:].mean()) / a[1:].std()))
            per[name][lab] = entry
            tot = sum(v["n1"] * v["n0"] for v in info.values())
            log(f"    [{name}] {lab:<22} sel {d['n1']:>5} bg {d['n0']:>5} "
                f"weight {d['n1'] * d['n0'] / tot:.5f}"
                + (f"  AUC {d['auc'][0]:.6f}  z {entry['z']:5.2f}"
                   if d["auc"] is not None else "  (dropped: empty side)"))
    out["per_stratum"] = per

    # ---- follow-up 2: paired null for the level drop -------------------------
    log("\n=== follow-up 2: paired null for the level drop "
        "D = AUC_pooled - AUC_strat(discovery) ===")
    D = pooled_all - strat_all
    D_obs, D_null = D[0], D[1:]
    zD = (D_obs - D_null.mean()) / D_null.std()
    pD = emp_p_two_sided(D_obs, D_null)
    log(f"  D observed      {D_obs:+.6f}   (0.597443 - 0.560924)")
    log(f"  D null          {D_null.mean():+.6f} +/- {D_null.std():.6f}"
        f"   [{np.percentile(D_null, 2.5):+.6f}, {np.percentile(D_null, 97.5):+.6f}]")
    log(f"  D observed sits at z = {zD:+.2f} within the null of D; "
        f"two-sided empirical p = {pD:.6f}")
    log("  pre-stated reading: |z| small / p large => the fall is the structural "
        "removal of cross-completeness pairs and nothing more;")
    log("                      D observed above the null => the pooled statistic lost "
        "more than the null did, and the fall is not purely structural.")
    out["level_drop"] = dict(D_obs=float(D_obs), D_null_mean=float(D_null.mean()),
                              D_null_sd=float(D_null.std()),
                              D_null_ci95=[float(np.percentile(D_null, 2.5)),
                                           float(np.percentile(D_null, 97.5))],
                              z=float(zD), p_two_sided=pD)

    exc_p = pooled_all[0] - pooled_all[1:].mean()
    exc_s = strat_all[0] - strat_all[1:].mean()
    log(f"\n  excess over own null: pooled {exc_p:+.6f}, stratified {exc_s:+.6f} "
        f"({100 * (exc_s - exc_p) / exc_p:+.1f}%)")
    out["excess_over_null"] = dict(pooled=float(exc_p), stratified=float(exc_s),
                                   pct_change=float(100 * (exc_s - exc_p) / exc_p))

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    log(f"\nwrote {args.out_json}")
    with open(str(args.out_json).replace(".json", ".log"), "w", encoding="utf-8",
              newline="\n") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
