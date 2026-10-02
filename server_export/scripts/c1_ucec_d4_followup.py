#!/usr/bin/env python3
"""Two follow-ups to D4 item 3, specified in review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md
section 7.4 BEFORE this script was run, and reported whatever they show.

Pure arithmetic on the already-saved c1_ucec_d4_stratified.npz, which holds all three
arms as full 1,001-column matrices (column 0 observed, 1..1000 the frozen nulls):
M_unrestricted, M_within_operator, M_within_operator_residualized. No permutation is
drawn here. `assemble` is IMPORTED from c1_run_test.py, not reimplemented, so the AUC
and the Spearman statistic cannot drift from the published ones.

1. The missing rho arm. c1_ucec_stratified_reread.py computes rho_op (L535) and
   rho_resid (L565) and stores neither, so the published two-statistic family
   (AUC and rho) was only half re-read. Reported here for all three arms.

2. Paired nulls for both transitions. All three arms share the same column index, so
   D can be formed column by column and its own permutation distribution read off.
   D observed inside the null spread means that transition moved the statistic no more
   than the permutations themselves move it.

Run (locally, on a downloaded copy of the npz -- no HPC, no cluster job):
  python c1_ucec_d4_followup.py --npz pinned/c1_ucec_d4_stratified.npz \
      --discovery-results ../pinned/ucec/residual_results_tumoronly.csv
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# The HPC keeps split_replication_power.py beside these scripts; this repo keeps it at
# the top level. Put both on the path so the same file runs in either place.
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))
import c1_run_test as C1  # noqa: E402  (imported for assemble, unchanged)

PUBLISHED_AUC = 0.5974425564900907
PUBLISHED_RHO = 0.2935675639399924
TOL = 1e-9

ARMS = [("unrestricted", "M_unrestricted"),
        ("within-operator", "M_within_operator"),
        ("residualised", "M_within_operator_residualized")]


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def emp_p_one_sided(obs, null):
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))


def emp_p_two_sided(obs, null):
    c = abs(obs - null.mean())
    return float((np.sum(np.abs(null - null.mean()) >= c) + 1) / (len(null) + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--discovery-results", required=True,
                    help="the discovery gene table; its incremental R^2 is `xs`, the "
                         "variable the published Spearman is computed against")
    ap.add_argument("--out-json", default=str(HERE / "pinned" / "c1_ucec_d4_followup.json"))
    args = ap.parse_args()

    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    z = np.load(args.npz, allow_pickle=True)
    genes, sel_mask = z["genes"], z["sel_mask"].astype(bool)
    log(f"npz {args.npz}")
    log(f"  sha256 {sha256(args.npz)}")
    log(f"  genes {len(genes)}  selected {int(sel_mask.sum())}  "
        f"background {int((~sel_mask).sum())}")

    disc = pd.read_csv(args.discovery_results).set_index("gene")["incremental_r2"]
    missing = [g for g in genes if g not in disc.index]
    if missing:
        sys.exit(f"FATAL: {len(missing)} genes absent from {args.discovery_results} "
                 f"(e.g. {missing[:5]}) -- cannot rebuild xs.")
    xs = np.array([disc[g] for g in genes])
    log(f"  xs (discovery incremental R^2) rebuilt for all {len(xs)} genes "
        f"from {args.discovery_results}")

    out = {"npz": args.npz, "sha256_npz": sha256(args.npz),
           "n_genes": int(len(genes)), "arms": {}, "transitions": {}}

    stats = {}
    log("\n=== per-arm statistics (AUC reproduced; rho is the half never stored) ===")
    for label, key in ARMS:
        if key not in z.files:
            sys.exit(f"FATAL: {args.npz} has no array {key!r}")
        auc, rho, n1, n0 = C1.assemble(z[key], xs, sel_mask)
        stats[label] = (auc, rho)
        e = {}
        for name, arr in (("auc", auc), ("rho", rho)):
            obs, null = arr[0], arr[1:]
            e[name] = dict(obs=float(obs), null_mean=float(null.mean()),
                           null_sd=float(null.std()),
                           p_one_sided=emp_p_one_sided(obs, null),
                           exceedances=int(np.sum(null >= obs)))
            log(f"  {label:<16} {name.upper():<4} {obs:.6f}  null {null.mean():.6f} "
                f"+/- {null.std():.6f}  p {e[name]['p_one_sided']:.6f} "
                f"({e[name]['exceedances']}/1000 exceed)")
        out["arms"][label] = e

    log("\n=== anchors: the unrestricted arm must reproduce the published values ===")
    ok = True
    for name, got, want in (("AUC", stats["unrestricted"][0][0], PUBLISHED_AUC),
                            ("rho", stats["unrestricted"][1][0], PUBLISHED_RHO)):
        good = abs(got - want) <= TOL
        ok &= good
        log(f"  {name}: got {got!r} want {want!r} -> {'OK' if good else 'MISMATCH'}")
    if not ok:
        sys.exit("FATAL: the unrestricted arm does not reproduce the published result; "
                 "refusing to report the other arms.")
    out["reproduces_published"] = True

    log("\n=== paired nulls for the two transitions ===")
    log("  pre-stated reading: D observed INSIDE the null spread of D means the")
    log("  transition moved the statistic no more than the permutations do; ABOVE it")
    log("  means the observed statistic lost more than the null did.")
    for tname, a, b in (("restriction (unrestricted -> within-operator)",
                         "unrestricted", "within-operator"),
                        ("residualisation (within-operator -> residualised)",
                         "within-operator", "residualised")):
        log(f"\n  -- {tname} --")
        tr = {}
        for i, sname in ((0, "auc"), (1, "rho")):
            D = stats[a][i] - stats[b][i]
            obs, null = D[0], D[1:]
            zz = (obs - null.mean()) / null.std()
            p2 = emp_p_two_sided(obs, null)
            tr[sname] = dict(D_obs=float(obs), D_null_mean=float(null.mean()),
                             D_null_sd=float(null.std()),
                             D_null_ci95=[float(np.percentile(null, 2.5)),
                                          float(np.percentile(null, 97.5))],
                             z=float(zz), p_two_sided=p2,
                             frac_of_drop_reproduced_by_null=(
                                 float(null.mean() / obs) if obs != 0 else None))
            log(f"     {sname.upper():<4} D obs {obs:+.6f}   D null {null.mean():+.6f} "
                f"+/- {null.std():.6f}  [{np.percentile(null, 2.5):+.6f}, "
                f"{np.percentile(null, 97.5):+.6f}]")
            # D observed is exactly 0 for the restriction transition, because that arm
            # changes only the null; print z and p there too rather than nothing, and
            # say so, instead of leaving a silent blank line in the record.
            tail = (f"   null reproduces {100 * null.mean() / obs:.0f}% of the observed"
                    f" drop" if obs != 0 else
                    "   (D observed is 0 by construction: this arm changes only the null)")
            log(f"          z {zz:+.2f}   two-sided p {p2:.6f}{tail}")
        out["transitions"][tname] = tr

    log("\n=== excess over each arm's own null, for the record ===")
    for label, _ in ARMS:
        a = stats[label][0]
        log(f"  {label:<16} AUC excess {a[0] - a[1:].mean():+.6f}")
    log("  like-for-like retention after residualisation "
        f"= {100 * (stats['residualised'][0][0] - stats['residualised'][0][1:].mean()) / (stats['within-operator'][0][0] - stats['within-operator'][0][1:].mean()):.1f}% "
        "of the within-operator excess")

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    log(f"\nwrote {args.out_json}")
    with open(str(args.out_json).replace(".json", ".log"), "w", encoding="utf-8",
              newline="\n") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
