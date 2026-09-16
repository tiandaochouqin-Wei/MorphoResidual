#!/usr/bin/env python3
"""make_score_cv_table.py (LOCAL) -- P2-19 fix (adversarial review, 2026-09-06).
Per-cohort out-of-fold cross-validated R^2 of the H&E-only score predicting the
measured residual-pathway score (translation_morph vs translation_meas, etc.),
computed directly from figdata/scores_<cohort>.csv (translation_morph IS already
the 5-fold OOF ridge prediction, per Methods sec:clinical -- this is 1 - SS_res/SS_tot
of measured vs morph, no re-fitting). Fixes the abstract/Results claim that the
ER-secretion range was 0.17-0.49 (in fact -0.07-0.49; negative in GBM)."""
import glob, os, numpy as np, pandas as pd

DD = "figdata"
FAMS = ["translation", "secretion_ER", "matrisome"]
rows = []
for p in sorted(glob.glob(f"{DD}/scores_*.csv")):
    c = os.path.basename(p).replace("scores_", "").replace(".csv", "").upper()
    d = pd.read_csv(p, index_col=0)
    row = {"cohort": c}
    for fam in FAMS:
        meas, morph = f"{fam}_meas", f"{fam}_morph"
        if meas in d.columns and morph in d.columns:
            y, pr = d[meas].values, d[morph].values
            m = np.isfinite(y) & np.isfinite(pr)
            ss_res = np.sum((y[m] - pr[m]) ** 2); ss_tot = np.sum((y[m] - y[m].mean()) ** 2)
            row[fam] = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rows.append(row)
res = pd.DataFrame(rows).set_index("cohort")
print(res.round(3).to_string())
for fam in FAMS:
    if fam in res:
        print(f"{fam}: range {res[fam].min():+.2f} to {res[fam].max():+.2f}")

BS = chr(92); EOL = BS + BS
NICE = {"translation": "Translation", "secretion_ER": "ER-secretion", "matrisome": "Matrisome"}
present = [f for f in FAMS if f in res.columns]
L = [BS + "begin{table}[htbp]" + BS + "centering" + BS + "small",
     BS + "caption{Out-of-fold cross-validated $R^2$ of the H" + BS + "&E-only residual-pathway score "
          "predicting the measured score (1 minus the residual sum of squares over the total sum of "
          "squares of measured versus H" + BS + "&E-only, per cohort). Matrisome is reported only for "
          "CCRCC and PDAC, where it was one of the enriched residual families ("
          + BS + "S" + BS + "ref{sec:clinical}); a dash means the family was not tested in that cohort.}",
     BS + "label{tab:clinical_r2}", BS + "begin{tabular}{l" + "r" * len(present) + "}" + BS + "toprule",
     "Cohort & " + " & ".join(NICE[f] for f in present) + " " + EOL + " " + BS + "midrule"]
for c, r in res.iterrows():
    L.append(c + " & " + " & ".join((f"{r[f]:+.2f}" if pd.notna(r[f]) else "---") for f in present) + " " + EOL)
L += [BS + "bottomrule", BS + "end{tabular}", BS + "end{table}"]
open("../SuppTable_clinical_r2.tex", "w", encoding="utf-8").write(chr(10).join(L) + chr(10))
print("\nwrote SuppTable_clinical_r2.tex")
