"""Build SuppTable_measerr.tex: the post hoc own-mRNA measurement-error test, one row per
cohort x set or stratum (POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md section 8.5, table
tab:measerr; Amendment 1 for the UCEC rerun folder and the a1 script).

Inputs (read only):
  review/MEAS_ERROR_2026-10-08/export/<c>/<folder>/theta_<c>_reading.csv   every printed number
  review/MEAS_ERROR_2026-10-08/export/<c>/<folder>/theta_<c>_sets.csv      z_N of the plex arm
with <folder> = local_theta (CCRCC frozen script; LUAD, GBM, PDAC a1 script) and
local_theta_rerun1 (UCEC, the one rerun of Amendment 1).

The readings (classes, riders, operator clauses) are applied by hand under section 7 and are
NOT computed here: READINGS below is copied from the final decision record of 2026-10-10
(final_reading.json, two independent readers reconciled), and the script only checks that
it covers exactly the read rows (pinned sets and strata with >= 20 genes and a theta).

Printing (section 7, Formats): theta, interval, lambda, kappa, theta_noise,min and phi_op to
2 decimals, z to 1 decimal; decimals are added while the rounded value lies on the other side
of, or on, a threshold that applies to it (the convention used for every sentence of the
note). theta_noise,min is always printed to 2 decimals.

Usage:  python make_supptable_measerr.py      (writes SuppTable_measerr.tex next to it)
"""
import math
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXPORT = os.path.join(HERE, "review", "MEAS_ERROR_2026-10-08", "export")
OUT = os.path.join(HERE, "SuppTable_measerr.tex")

COHORTS = [  # display name, code, output folder, patients, header note
    ("CCRCC", "ccrcc", "local_theta", 103,
     r"$\lambda_{lo}$ measured (median replicate ICC); primary null within operator; frozen script"),
    ("LUAD", "luad", "local_theta", 105,
     r"$\lambda_{lo}$ transferred from CCRCC and UCEC; primary null within operator; purity for 103 of 105 patients"),
    ("UCEC", "ucec", "local_theta_rerun1", 100,
     r"$\lambda_{lo}$ measured (median replicate ICC); primary null within operator; rerun under Amendment~1; purity for 96 of 100 patients"),
    ("GBM", "gbm", "local_theta", 99,
     r"$\lambda_{lo}$ transferred from CCRCC and UCEC; primary null within operator; purity not available"),
    ("PDAC", "pdac", "local_theta", 137,
     r"$\lambda_{lo}$ transferred from CCRCC and UCEC; primary null within operator $\times$ scanner; purity for 97 of 137 patients"),
]

# (class, rider, operator-clause form, unrestricted reading) per read row; final decision record.
SD, OM = "same direction", "operator-mediated transcript alignment"
READINGS = {
    "CCRCC": {
        "pinned_sig": ("R3", "", SD, "R3"),
        "pinned_sig|depthT1": ("R3", "", SD, "R3"),
        "pinned_sig|depthT2": ("R3", "", SD, "R3"),
        "pinned_sig|depthT3": ("R3", "", SD, "R3"),
        "pinned_sig&ecm": ("R3", "", None, None),
        "pinned_sig&folding": ("R3", "", SD, "R3"),
        "pinned_sig&folding|depthT1": ("R3", "", SD, "R3"),
        "pinned_sig&folding|depthT3": ("R3", "", None, None),
        "pinned_sig&ptm": ("R3", "", SD, "R3"),
        "pinned_sig&ptm|depthT1": ("R3", "", SD, "R3"),
        "pinned_sig&ptm|depthT2": ("R3", "", SD, "R3"),
        "pinned_sig&ptm|depthT3": ("R3", "", SD, "R3"),
        "pinned_sig&secretion": ("R3", "", SD, "R3"),
        "pinned_sig&secretion|depthT1": ("R3", "", SD, "R3"),
        "pinned_sig&secretion|depthT2": ("R3", "", SD, "R3"),
        "pinned_sig&secretion|depthT3": ("R1", "spline-sensitive", None, None),
        "pinned_sig&splicing": ("R3", "", SD, "R3"),
        "pinned_sig&splicing|depthT1": ("R3", "", SD, "R3"),
        "pinned_sig&splicing|depthT2": ("R3", "", SD, "R3"),
        "pinned_sig&splicing|depthT3": ("R3", "", None, None),
        "pinned_sig&translation": ("R3", "", SD, "R3"),
        "pinned_sig&translation|depthT1": ("R3", "", SD, "R3"),
        "pinned_sig&translation|depthT2": ("R3", "", SD, "R3"),
        "pinned_sig&translation|depthT3": ("R3", "", SD, "R3"),
    },
    "LUAD": {
        "pinned_sig": ("R2", "", SD, "R2"),
        "pinned_sig|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig|depthT2": ("R2", "", SD, "R2"),
        "pinned_sig|depthT3": ("R2", "", SD, "R2"),
        "pinned_sig&ecm": ("R2", "carried by purity", None, None),
        "pinned_sig&ecm|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig&ecm|depthT2": ("R1", "spline-stable", None, None),
        "pinned_sig&ecm|depthT3": ("R1", "spline-stable", None, None),
        "pinned_sig&folding": ("R1", "spline-sensitive", None, None),
        "pinned_sig&ptm": ("R2", "", SD, "R2"),
        "pinned_sig&ptm|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig&ptm|depthT2": ("R2", "", None, None),
        "pinned_sig&ptm|depthT3": ("R2", "", SD, "R2"),
        "pinned_sig&secretion": ("R2", "", SD, "R2"),
        "pinned_sig&secretion|depthT1": ("R2", "", SD, "R3"),
        "pinned_sig&secretion|depthT2": ("R2", "", None, None),
        "pinned_sig&secretion|depthT3": ("R2", "", SD, "R2"),
        "pinned_sig&splicing": ("R2", "carried by purity", SD, "R2"),
        "pinned_sig&splicing|depthT1": ("R2", "carried by purity", None, None),
        "pinned_sig&splicing|depthT2": ("R2", "carried by purity", SD, "R2"),
        "pinned_sig&splicing|depthT3": ("R2", "carried by purity", SD, "R2"),
        "pinned_sig&translation": ("R2", "", SD, "R2"),
        "pinned_sig&translation|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig&translation|depthT2": ("R2", "", SD, "R2"),
        "pinned_sig&translation|depthT3": ("R2", "", SD, "R3"),
    },
    "UCEC": {
        "pinned_sig": ("R2", "", None, None),
        "pinned_sig|depthT1": ("R2", "", SD, "R3"),
        "pinned_sig|depthT2": ("R2", "carried by purity", None, None),
        "pinned_sig|depthT3": ("R2", "carried by purity", None, None),
        "pinned_sig&ecm": ("R2", "carried by purity", None, None),
        "pinned_sig&folding": ("R1", "spline-sensitive", OM, "R2"),
        "pinned_sig&folding|depthT1": ("R3", "", SD, "R3"),
        "pinned_sig&folding|depthT2": ("R1", "spline-stable", None, None),
        "pinned_sig&folding|depthT3": ("R1", "spline-stable", OM, "R2"),
        "pinned_sig&ptm": ("R2", "", SD, "R3"),
        "pinned_sig&ptm|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig&ptm|depthT2": ("R2", "carried by purity", None, None),
        "pinned_sig&ptm|depthT3": ("R1", "spline-sensitive", None, None),
        "pinned_sig&secretion": ("R1", "spline-sensitive; monotone fall with missingness", OM, "R3"),
        "pinned_sig&secretion|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig&secretion|depthT2": ("R7", "", None, None),
        "pinned_sig&secretion|depthT3": ("R1", "spline-sensitive", None, None),
        "pinned_sig&splicing": ("R7", "", None, None),
        "pinned_sig&splicing|depthT1": ("R1", "spline-sensitive", OM, "R2"),
        "pinned_sig&splicing|depthT2": ("R1", "spline-sensitive", None, None),
        "pinned_sig&splicing|depthT3": ("R1", "spline-stable", None, None),
        "pinned_sig&translation": ("R2", "carried by purity", SD, "R2"),
        "pinned_sig&translation|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig&translation|depthT2": ("R2", "carried by purity", None, None),
        "pinned_sig&translation|depthT3": ("R1", "spline-sensitive", None, None),
    },
    "GBM": {
        "pinned_sig": ("R2", "", None, None),
        "pinned_sig|depthT1": ("R2", "", None, None),
        "pinned_sig|depthT2": ("R2", "", None, None),
        "pinned_sig|depthT3": ("R2", "", None, None),
        "pinned_sig&ptm": ("R2", "", None, None),
        "pinned_sig&ptm|depthT1": ("R2", "", None, None),
        "pinned_sig&ptm|depthT2": ("R2", "", SD, "R2"),
        "pinned_sig&ptm|depthT3": ("R2", "", None, None),
        "pinned_sig&secretion": ("R2", "", None, None),
        "pinned_sig&secretion|depthT1": ("R2", "", None, None),
        "pinned_sig&secretion|depthT2": ("R1", "spline-sensitive", None, None),
        "pinned_sig&secretion|depthT3": ("R1", "spline-stable", None, None),
        "pinned_sig&splicing": ("R2", "", None, None),
        "pinned_sig&translation": ("R2", "", None, None),
        "pinned_sig&translation|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig&translation|depthT2": ("R2", "", None, None),
        "pinned_sig&translation|depthT3": ("R2", "", None, None),
    },
    "PDAC": {
        "pinned_sig": ("R2", "", SD, "R3"),
        "pinned_sig|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig|depthT2": ("R2", "", SD, "R2"),
        "pinned_sig|depthT3": ("R2", "", SD, "R3"),
        "pinned_sig&ecm": ("R2", "", SD, "R2"),
        "pinned_sig&ecm|depthT1": ("R2", "", SD, "R2"),
        "pinned_sig&ecm|depthT2": ("R2", "carried by purity", None, None),
        "pinned_sig&ecm|depthT3": ("R2", "", SD, "R2"),
        "pinned_sig&folding": ("R2", "carried by purity", SD, "R3"),
        "pinned_sig&ptm": ("R2", "", SD, "R3"),
        "pinned_sig&ptm|depthT1": ("R2", "carried by purity", None, None),
        "pinned_sig&ptm|depthT2": ("R2", "carried by purity", SD, "R2"),
        "pinned_sig&ptm|depthT3": ("R2", "", SD, "R2"),
        "pinned_sig&secretion": ("R2", "", SD, "R3"),
        "pinned_sig&secretion|depthT1": ("R2", "carried by purity", SD, "R3"),
        "pinned_sig&secretion|depthT2": ("R2", "", SD, "R3"),
        "pinned_sig&secretion|depthT3": ("R2", "", SD, "R3"),
        "pinned_sig&splicing": ("R1", "spline-sensitive", None, None),
        "pinned_sig&translation": ("R2", "", SD, "R3"),
        "pinned_sig&translation|depthT1": ("R2", "carried by purity", SD, "R3"),
        "pinned_sig&translation|depthT2": ("R2", "carried by purity", SD, "R2"),
        "pinned_sig&translation|depthT3": ("R2", "", SD, "R3"),
    },
}
RIDER_ABBR = {"carried by purity": "pur.", "spline-sensitive": "spl.-sens.", "spline-stable": "spl.-stab.",
              "monotone fall with missingness": r"miss.\ fall"}
K_EXPECTED = {"CCRCC": 23, "LUAD": 24, "UCEC": 24, "GBM": 16, "PDAC": 21}  # read strata (section 4)

Z_TH = [(-3.0, "le"), (2.0, "lt"), (3.0, "lt")]


def _side(v, t, kind):
    return v < t if kind == "lt" else (v <= t if kind == "le" else v > t)


def num(v, ths=(), d0=2):
    """Round to d0 decimals, adding decimals while the rounded value is on, or on the other
    side of, any threshold in ths; '---' for a missing value; minus sign typeset."""
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "---"
    v = float(v)
    d = d0
    while d < 8:
        p = round(v, d)
        if all(p != t and _side(p, t, k) == _side(v, t, k) for t, k in ths):
            break
        d += 1
    s = f"{round(v, d):.{d}f}"
    if s.startswith("-"):
        s = "$-$" + s[1:]
    return s


def band(n):
    return 0.05 if n >= 100 else 0.10


def theta_th(n, tnm):
    return [(band(n), "le"), (-band(n), "lt"), (0.0, "gt"), (tnm, "lt"), (-0.3, "le")]


def name_tex(s):
    t = s.replace("_", r"\_").replace("&", r"\allowbreak\&").replace("|", r"\allowbreak|")
    t = t.replace("<=", r"$\leq$")
    return r"\texttt{" + t + "}"


def fin(x):
    return x is not None and isinstance(x, (int, float)) and math.isfinite(float(x))


def row_cells(c, r, plex_z, reading, opclause):
    n = int(r["n_genes"])
    tnm = float(r["theta_noise_min"])
    th = theta_th(n, tnm)
    s = num(r["s"], [(x, "lt") for x in (0.02, 0.05, 0.10, 0.20, 0.50)]) if fin(r["s"]) else "---"
    lam_lo = num(r["lambda_lo"], [(0.55, "le"), (0.65, "le")])
    theta = num(r["theta_pc"], th)
    if fin(r["theta_pc_ci_lo"]):
        ci = "[" + num(r["theta_pc_ci_lo"], [(tnm, "le")]) + ", " + num(r["theta_pc_ci_hi"], [(tnm, "lt")]) + "]"
    else:
        ci = "---"
    zD = num(r["D_z"], [(3.0, "lt")], 1)
    zN = num(r["Nc_z"], Z_TH, 1)
    zU = num(r["Nc_z_unres"], Z_TH, 1)
    zP = num(plex_z, Z_TH, 1)
    zS = num(r["S_cross_z"], (), 1)
    phi = num(r["phi_op_D"]) + "/" + num(r["phi_op_dp"])
    # purity variant
    if c == "GBM":
        pur = "not available"
    else:
        npur = r["purity_n_genes"]
        if not fin(npur) or int(npur) < 5:
            pur = "not available ($<$30 matched)"
        elif not fin(r["purity_theta_pc"]):
            pur = "n.r. (" + num(r["purity_Nc_z"], [(2.0, "lt")], 1) + ")"
            if int(npur) < n:
                pur += " [%d]" % int(npur)
        else:
            pur = (num(r["purity_theta_pc"], [(band(n), "le"), (-band(n), "lt")]) + " ("
                   + num(r["purity_Nc_z"], [(2.0, "lt")], 1) + ")")
            if int(npur) < n:
                pur += " [%d]" % int(npur)
    spl = (num(r["spline_theta_pc"], [(band(n), "le"), (-band(n), "lt")]) + " ("
           + num(r["spline_Nc_z"], [(2.0, "lt")], 1) + ")") if fin(r["spline_theta_pc"]) else (
        "--- (" + num(r["spline_Nc_z"], [(2.0, "lt")], 1) + ")")
    rnp = (num(r["rnapc_theta_pc"], [(-0.3, "le")]) + " ("
           + num(r["rnapc_Nc_z"], [(-3.0, "le")], 1) + ")") if fin(r["rnapc_theta_pc"]) else (
        "--- (" + num(r["rnapc_Nc_z"], [(-3.0, "le")], 1) + ")")
    raw = num(r["theta_raw"])
    return [name_tex(r["set"]), str(n), s, lam_lo, num(r["lambda_hi"]), num(r["kappa"]),
            num(tnm), theta, ci, zD, zN, zU, zP, zS, phi, pur, spl, rnp, raw, reading, opclause]


def _stk(*lines):
    return r"\shortstack[r]{" + r"\\".join(lines) + "}"


HEAD = [r"Set or stratum", r"$n$", r"$s$", r"$\lambda_{lo}$", r"$\lambda_{hi}$", r"$\kappa$",
        r"$\theta_{\min}$", r"$\theta_{pc}$", r"95\% interval", r"$z_D$", r"$z_N$",
        r"$z_N^{\mathrm{unr}}$", r"$z_N^{\mathrm{plex}}$", r"$z_S$", _stk(r"$\phi_{op}$", r"$D$\,/\,$\Delta_p$"),
        _stk("Purity", r"$\theta$ ($z$)"), _stk("Spline", r"$\theta$ ($z$)"), _stk("RNA-PC", r"$\theta$ ($z$)"),
        _stk("Raw", r"$\theta$"), r"Reading", _stk("Operator", "clause")]
NCOL = len(HEAD)

CAPTION = (
    r"\textbf{Own-mRNA measurement-error test: every set and stratum.} Post hoc; specified on "
    r"8 October 2026 (UTC) in a dated file released publicly before it ran, with Amendment~1 "
    r"(9 October 2026, UTC) for the UCEC rerun and the LUAD, GBM and PDAC runs (Supplementary "
    r"Note~\ref{snote:measerr}); reported whatever the result. The sets are the "
    r"batch-uncorrected significant sets. One row per set or stratum written by the analysis "
    r"(strata with $\geq15$ genes; $\theta$ read with $\geq20$). Set names: \texttt{pinned\_sig}, "
    r"the published significant set restricted to fully quantified proteins (the primary set); "
    r"\texttt{\&family}, its intersection with a Reactome signature family; \texttt{all} and an "
    r"unprefixed family name, unselected sets (descriptive); \texttt{|depthT1}--\texttt{3}, "
    r"tertiles of raw-count depth (fully quantified proteins); \texttt{|miss0-5}, "
    r"\texttt{|miss5-20}, \texttt{|miss20+}, \texttt{|miss$\leq$5}, proteins with that percentage "
    r"of missing values (with fully quantified ones for $\leq$5; MNAR diagnostic, not read). "
    r"$n$, genes; $s$, selection strength; $\lambda_{lo}$, lower reliability bound (cohort "
    r"header); $\lambda_{hi}$, median Poisson bound (used by no reading); $\kappa$, selection "
    r"floor (1 for unselected sets); $\theta_{\min}=\kappa\lambda_{lo}/(1-\lambda_{lo})$, the "
    r"dilution minimum; $\theta_{pc}$, permutation-centred alignment, with its 95\% "
    r"patient-bootstrap interval (200 resamples); $z_D$, $z$ of the morphology-only signal; "
    r"$z_N$, $z_N^{\mathrm{unr}}$, $z_N^{\mathrm{plex}}$, $z$ of the sign-weighted alignment "
    r"under the primary (within operator; within operator $\times$ scanner in PDAC), "
    r"unrestricted and within-TMT-plex (descriptive) nulls; $z_S$, sign test (reported only); "
    r"$\phi_{op}$, share of the permutation-centred signal $D$ and of the nested increment "
    r"$\Delta_p$ that lies between operators; purity, spline and 20-RNA-PC baselines, "
    r"$\theta_{pc}$ and $z_N$ of each variant (purity: [number of genes with $\geq30$ "
    r"purity-matched patients] where fewer than $n$; n.r., not readable, fewer than 20 such "
    r"genes); raw "
    r"$\theta$, not permutation-centred. Reading (rules in Supplementary "
    r"Note~\ref{snote:measerr}): R1, not aligned with the measured transcript; R2, alignment "
    r"reached what dilution produces; R3, aligned below the dilution minimum; R7, undecidable, "
    r"no reading (\textbf{bold}, the primary set); riders: pur., carried by purity; spl.-sens.\ and "
    r"spl.-stab., spline-sensitive and spline-stable; miss.\ fall, a monotone fall of $z_N$ with "
    r"missing protein. Operator clause, written where the unrestricted $z_N\geq2$ or $\leq-3$ for a "
    r"set or stratum with a reading: same, same direction; op.-med., operator-mediated transcript "
    r"alignment; the unrestricted reading in parentheses. MNAR diag., missingness stratum "
    r"(diagnostic); descr., unselected set; unselected sets, missingness strata and selected sets or "
    r"strata with fewer than 20 genes receive no reading. No set or stratum read R1w, R4, R4w or R6. Values to two "
    r"decimals ($z$ to one), with decimals added where rounding would place a value on, or on "
    r"the other side of, a threshold that applies to it; --- , not computed."
)


def main():
    out = []
    out.append("% Supplementary Table: post hoc own-mRNA measurement-error test (tab:measerr), governed by")
    out.append("% review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md (sidecar written 2026-10-08T01:04:28Z) and its")
    out.append("% Amendment 1 (sidecar written 2026-10-09T19:32:41Z). Generated by make_supptable_measerr.py from")
    out.append("% review/MEAS_ERROR_2026-10-08/export/<c>/{local_theta,local_theta_rerun1}/theta_<c>_{reading,sets}.csv.")
    out.append(r"\begin{landscape}")
    out.append(r"\begingroup")
    out.append(r"\fontsize{5}{6}\selectfont")
    out.append(r"\setlength{\tabcolsep}{1.0pt}")
    out.append(r"\renewcommand{\arraystretch}{0.95}")
    out.append(r"\setlength{\LTcapwidth}{\linewidth}")
    out.append(r"\begin{longtable}{l" + "r" * 18 + "ll}")
    out.append(r"\captionsetup{font={scriptsize,stretch=1.0}}")
    out.append(r"\caption{" + CAPTION + r"}\label{tab:measerr}\\")
    out.append(r"\toprule")
    out.append(" & ".join(HEAD) + r" \\")
    out.append(r"\midrule")
    out.append(r"\endfirsthead")
    out.append(r"\multicolumn{%d}{l}{\emph{Supplementary Table~\ref{tab:measerr}, continued}}\\" % NCOL)
    out.append(r"\toprule")
    out.append(" & ".join(HEAD) + r" \\")
    out.append(r"\midrule")
    out.append(r"\endhead")
    out.append(r"\bottomrule")
    out.append(r"\endlastfoot")

    for disp, code, folder, npat, note in COHORTS:
        base = os.path.join(EXPORT, code, folder)
        rd = pd.read_csv(os.path.join(base, f"theta_{code}_reading.csv"))
        st = pd.read_csv(os.path.join(base, f"theta_{code}_sets.csv"))
        plex = st[st["arm"] == "plex"].set_index("set")["Nc_z"]
        readings = READINGS[disp]
        # coverage check: read rows = pinned sets/strata (not missingness) with >= 20 genes and a theta
        read_rows = rd[rd["set"].str.startswith("pinned_sig") & (rd["kind"] != "miss")
                       & (rd["n_genes"] >= 20) & rd["theta_pc"].notna()]["set"].tolist()
        assert sorted(read_rows) == sorted(readings), (disp, set(read_rows) ^ set(readings))
        assert len(readings) - 1 == K_EXPECTED[disp], disp
        sel = rd["set"].str.startswith("pinned_sig")
        blocks = [
            ("Primary set and selected sets and strata", rd[sel & (rd["kind"] != "miss")]),
            ("Missingness strata of the selected sets (MNAR diagnostic)", rd[sel & (rd["kind"] == "miss")]),
            ("Unselected sets and their strata (descriptive)", rd[~sel]),
        ]
        out.append(r"\midrule")
        out.append(r"\multicolumn{%d}{l}{\textbf{%s} (%d patients; %s)}\\" % (NCOL, disp, npat, note))
        for title, df in blocks:
            out.append(r"\multicolumn{%d}{l}{\emph{%s}}\\" % (NCOL, title))
            for _, r in df.iterrows():
                key = r["set"]
                if key in readings:
                    cls, rider, form, ucls = readings[key]
                    rid = ", ".join(RIDER_ABBR[x] for x in rider.split("; ")) if rider else ""
                    reading = (r"\textbf{" + cls + "}" if key == "pinned_sig" else cls) + (", " + rid if rid else "")
                    if form == SD:
                        opc = "same (%s)" % ucls
                    elif form == OM:
                        opc = "op.-med. (%s)" % ucls
                    else:
                        opc = ""
                elif r["kind"] == "miss":
                    reading, opc = "MNAR diag.", ""
                elif key.startswith("pinned_sig"):
                    reading, opc = r"not read ($n<20$)", ""
                else:
                    reading, opc = "descr.", ""
                pz = plex.get(key, float("nan"))
                out.append(" & ".join(row_cells(disp, r, pz, reading, opc)) + r" \\")
    out.append(r"\end{longtable}")
    out.append(r"\endgroup")
    out.append(r"\end{landscape}")
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
