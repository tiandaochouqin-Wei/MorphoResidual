#!/usr/bin/env python3
"""figdata/fig6v7_a_covariates.csv: the four primary cells and the capacity placebos of the post hoc covariate re-reads
(pre-specification 2026-11-21; results record tag c1-covariates-posthoc-results-2026-11-23), for the rows added to Fig. 6a.
Every value is copied from review/C1_COVARIATES_2026-11-21/results_record/L2/cells_L2.csv; nothing is computed except
the count e (checked against the file's own e)."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
rows = list(csv.DictReader(open(ROOT / "review/C1_COVARIATES_2026-11-21/results_record/L2/cells_L2.csv", encoding="utf-8")))
out = []
for coh in ("ucec", "luad"):
    for cov, mat in (("C3", "A_C3"), ("R20", "A_R20")):
        hit = [r for r in rows if r["cohort"] == coh and r["family"] == "primary" and r["set"] == "S_M" and r["matrix"] == mat]
        assert len(hit) == 1
        r = hit[0]
        assert r["null"] == "unrestricted" and r["status"] == "ok" and int(r["n_null"]) == 1000
        rec = dict(cohort=coh, covariate=cov, set="S_M", matrix=mat, n_set=r["n1"], n_background=r["n0"], n_perm=r["n_null"],
                   obs_auc=r["auc"], null_mean=r["null_mean"], null_sd=r["null_sd"], margin=r["margin"],
                   p_one_sided=r["p"], n_exceed=r["exceed"], m_R="", e="", p_abs="", m_j="")
        if cov == "R20":
            mj = [float(x) for x in r["m_j"].split(";")]
            assert len(mj) == 119 and int(r["e"]) == sum(1 for x in mj if x <= float(r["m_R"]))
            rec.update(m_R=r["m_R"], e=r["e"], p_abs=r["p_abs"], m_j=r["m_j"])
        out.append(rec)
dst = Path(__file__).resolve().parent / "figdata" / "fig6v7_a_covariates.csv"
with open(dst, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()), lineterminator="\n")
    w.writeheader()
    w.writerows(out)
print("wrote", dst, len(out), "rows")
