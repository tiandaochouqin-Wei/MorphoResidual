#!/usr/bin/env python3
"""mtor_mechanism.py (LOCAL) — mechanism anchor: is the morphology-readable TRANSLATION
residual the mTORC1 translational programme?  Uses CPTAC tumour phosphoproteomics
(LinkedOmics site-level matrices in figdata/phospho_<cohort>.txt).

Per cohort:
  mTORC1 score  = mean z of canonical mTORC1 output phosphosites
                  (RPS6 S235/S236/S240/S244, EIF4EBP1 T37/T46/S65/T70, EIF4B S422, RPS6KB1 T389)
  ERK score     = mean z of MAPK1 T185/Y187, MAPK3 T202/Y204   (specificity control)
  then Spearman rho of each score with the MEASURED translation residual, the H&E-ONLY
  translation residual, and the ER-secretion residual; per-site rho heatmap; and a
  mediation-style check: partial correlation of the H&E-only score with the measured
  residual controlling for the mTORC1 score (how much of the morphology link mTORC1 explains).
Outputs: mtor_results.csv, mtor_sites.csv, Fig_mtor.{svg,pdf,png,tiff}, ../SuppTable_mtor.tex.
Loader is defensive about layout (sites x samples or samples x sites) and site-ID styles
('RPS6_S235', 'RPS6-S235s', 'RPS6:S235', 'NP_...:S235 RPS6', 'RPS6_S235S236')."""
import os, re, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as st
import mrstyle as S

DD = "figdata"
COH = S.COH
# Phosphosite parsing, the mTORC1/ERK/ribosomal definitions and partial_r now live
# in mtor_common.py, shared with mtor_protein_normalised.py (T1-4) so the two
# analyses cannot drift apart. Extraction was verified output-identical.
from mtor_common import (MTOR, ERK, RIBO_PAT, norm_id, parse_site,
                         load_phospho, site_score, partial_r)  # noqa: F401

rows, site_rows, panels, ribo_bg_rows, ribo_partial_rows = [], [], {}, [], []
for c in COH:
    p = f"{DD}/phospho_{c.lower()}.txt"
    if not os.path.exists(p):
        print(f"[{c}] missing {p} (skipped)"); continue
    ph = load_phospho(p)
    mt, mtz = site_score(ph, MTOR); er, _ = site_score(ph, ERK)
    s6, _ = site_score(ph, {"RPS6": MTOR["RPS6"]})                 # S6K1 arm (RPS6 phosphorylation)
    bp, _ = site_score(ph, {"EIF4EBP1": MTOR["EIF4EBP1"]})         # 4E-BP1 arm
    if mt is None:
        print(f"[{c}] no mTORC1 sites found; index examples: {list(ph.index[:5])}"); continue
    sc = pd.read_csv(f"{DD}/scores_{c.lower()}.csv", index_col=0); sc.index = [norm_id(i) for i in sc.index]
    d = sc.join(mt.rename("mtor"), how="inner")
    if er is not None:
        d = d.join(er.rename("erk"), how="left")
    if s6 is not None:
        d = d.join(s6.rename("s6"), how="left")
    if bp is not None:
        d = d.join(bp.rename("bp1"), how="left")
    n = len(d)
    def rho(a, b):
        m = d[a].notna() & d[b].notna()
        return st.spearmanr(d.loc[m, a], d.loc[m, b]) if m.sum() >= 10 else (np.nan, np.nan)
    r_meas, p_meas = rho("mtor", "translation_meas")
    r_morph, p_morph = rho("mtor", "translation_morph")
    r_sec, p_sec = rho("mtor", "secretion_ER_meas") if "secretion_ER_meas" in d else (np.nan, np.nan)
    r_erk, p_erk = rho("erk", "translation_meas") if "erk" in d else (np.nan, np.nan)
    r_s6, p_s6 = rho("s6", "translation_meas") if "s6" in d else (np.nan, np.nan)
    r_bp, p_bp = rho("bp1", "translation_meas") if "bp1" in d else (np.nan, np.nan)
    r_full = st.pearsonr(d["translation_meas"], d["translation_morph"])[0]
    r_part, n_part = partial_r(d["translation_morph"].values, d["translation_meas"].values, d["mtor"].values)
    rows.append(dict(cohort=c, n=n, n_mtor_sites=mtz.shape[1],
                     rho_mtor_vs_translation_meas=r_meas, p_meas=p_meas,
                     rho_mtor_vs_translation_HE=r_morph, p_HE=p_morph,
                     rho_mtor_vs_secretion_meas=r_sec, p_sec=p_sec,
                     rho_ERK_vs_translation_meas=r_erk, p_erk=p_erk,
                     rho_S6arm_vs_translation_meas=r_s6, p_s6=p_s6,
                     rho_4EBP1arm_vs_translation_meas=r_bp, p_bp=p_bp,
                     r_HE_vs_meas=r_full, r_HE_vs_meas_given_mtor=r_part,
                     pct_explained_by_mtor=100 * (1 - r_part / r_full) if r_full else np.nan))
    for s in mtz.columns:
        m = d.index.intersection(mtz.index); zz = mtz.loc[m, s]; tm = d.loc[m, "translation_meas"]
        mm = zz.notna() & tm.notna()
        site_rows.append(dict(cohort=c, site=s, rho=st.spearmanr(zz[mm], tm[mm])[0] if mm.sum() >= 10 else np.nan, n=int(mm.sum())))
    panels[c] = d
    # ---- P0-4 control: non-mTORC1 ribosomal phosphosites, same patients/score, same >=50%-complete filter ----
    # skip gene-level-only cohorts (CCRCC): row ids are bare gene symbols, so a trailing digit in the
    # symbol itself (e.g. RPS10 -> "S10") is mistaken for a residue by parse_site -- not a real phosphosite.
    is_gene_level_only = any("(all sites)" in col for col in mtz.columns)
    seen, ribo_cols = set(), {}
    for rid in ([] if is_gene_level_only else ph.index):
        parsed = parse_site(rid)
        if not parsed:
            continue
        g, res_set = parsed
        if not RIBO_PAT.match(g) or g == "RPS6":
            continue
        for rr in sorted(res_set):      # sorted, not set order: Python randomises string
                                        # hashing per process, so an unsorted set made the row
                                        # order of mtor_ribosomal_control.csv differ between
                                        # runs. Values were unaffected (every downstream use is
                                        # a mean or a median) but the file was not reproducible.
            sid = f"{g} {rr}"
            if sid in seen:
                continue
            seen.add(sid)
            zvals = ph.loc[rid]
            zvals = (zvals - zvals.mean()) / zvals.std(ddof=0)
            zj = zvals.reindex(d.index)
            mm = zj.notna() & d["translation_meas"].notna()
            if mm.sum() < max(10, 0.5 * len(d)):
                continue
            rho_i, _ = st.spearmanr(zj[mm], d.loc[mm, "translation_meas"])
            ribo_bg_rows.append(dict(cohort=c, site=sid, rho=rho_i, n=int(mm.sum())))
            ribo_cols[sid] = zj
    if ribo_cols:
        ribo_proxy = pd.DataFrame(ribo_cols).mean(axis=1)               # proxy for ribosomal-protein abundance
        d = d.join(ribo_proxy.rename("ribo_proxy"), how="left")
        r_ribo, _ = rho("ribo_proxy", "translation_meas")
        s6_given_ribo = r_bp2 = np.nan
        if "s6" in d:
            s6_given_ribo, n_p1 = partial_r(d["s6"].values, d["translation_meas"].values, d["ribo_proxy"].values)
            ribo_given_s6, _ = partial_r(d["ribo_proxy"].values, d["translation_meas"].values, d["s6"].values)
        else:
            ribo_given_s6 = np.nan
        bp_given_ribo = np.nan
        if "bp1" in d:
            bp_given_ribo, _ = partial_r(d["bp1"].values, d["translation_meas"].values, d["ribo_proxy"].values)
        ribo_partial_rows.append(dict(cohort=c, n_ribo_sites=len(ribo_cols), rho_ribo_proxy=r_ribo,
                                      rho_S6arm_given_riboproxy=s6_given_ribo,
                                      rho_riboproxy_given_S6arm=ribo_given_s6,
                                      rho_4EBP1arm_given_riboproxy=bp_given_ribo))
    print(f"[{c}] n={n} sites={mtz.shape[1]}  mTORC1~transl_meas rho={r_meas:+.2f} (p={p_meas:.1e})  "
          f"~transl_HE rho={r_morph:+.2f}  ~secretion rho={r_sec:+.2f}  ERK~transl rho={r_erk:+.2f}  "
          f"HE-meas r={r_full:.2f} -> {r_part:.2f} given mTORC1")

res = pd.DataFrame(rows); sites = pd.DataFrame(site_rows); ribo_bg = pd.DataFrame(ribo_bg_rows)
ribo_partial = pd.DataFrame(ribo_partial_rows)
if not len(res):
    raise SystemExit("no cohort processed")
res.to_csv("mtor_results.csv", index=False); sites.to_csv("mtor_sites.csv", index=False)
ribo_bg.to_csv(f"{DD}/mtor_ribosomal_control.csv", index=False)
ribo_partial.to_csv(f"{DD}/mtor_ribosomal_partial.csv", index=False)
pd.set_option("display.width", 200); print(res.round(3).to_string(index=False))
if len(ribo_bg):
    for c_, g in ribo_bg.groupby("cohort"):
        s6r = res.loc[res.cohort == c_, "rho_S6arm_vs_translation_meas"]
        s6r = float(s6r.iloc[0]) if len(s6r) and pd.notna(s6r.iloc[0]) else np.nan
        pct = 100 * (g["rho"] < s6r).mean() if pd.notna(s6r) else np.nan
        print(f"[{c_}] ribosomal-background control: n_sites={len(g)} median={g.rho.median():+.2f} "
              f"p90={g.rho.quantile(.9):+.2f} max={g.rho.max():+.2f} ({g.loc[g.rho.idxmax(),'site']})  "
              f"S6-arm percentile in this background: {pct:.0f}%")
    print(ribo_partial.round(3).to_string(index=False))

# ---------------- figure: a) score x cohort rho heatmap (BH-starred), b) per-site heatmap, --------
# ---------------- c) best-cohort scatter, d) P0-4 ribosomal-phosphosite background control ------
def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (q) for a flat array; NaN entries stay NaN."""
    p = np.asarray(p, float); q = np.full(p.shape, np.nan)
    m = np.isfinite(p); pv = p[m]; n = pv.size
    if n:
        o = np.argsort(pv); r = pv[o] * n / np.arange(1, n + 1)
        r = np.minimum(np.minimum.accumulate(r[::-1])[::-1], 1.0)
        qq = np.empty(n); qq[o] = r; q[m] = qq
    return q


def _compact(res):
    """'S235/S236/S240' -> 'S235/236/240' (drop a residue letter that repeats the previous
    one); 'S65/T70' is left alone. Display only -- mtor_sites.csv keeps the full ids."""
    parts = res.split("/")
    out, prev = [parts[0]], parts[0][:1]
    for p in parts[1:]:
        out.append(p[1:] if p[:1] == prev else p)
        prev = p[:1]
    return "/".join(out)


def short_site(s):
    """Display-only site label: gene aliases already used in panel a (EIF4EBP1 = 4E-BP1,
    RPS6KB1 = S6K1) plus residue compaction, so the labels fit at 7 pt."""
    g, _, rest = s.partition(" ")
    g = {"EIF4EBP1": "4E-BP1", "RPS6KB1": "S6K1", "EIF4B": "eIF4B"}.get(g, g)
    if rest and not rest.startswith("("):
        rest = _compact(rest)
    return f"{g} {rest}" if rest else g


# Canvas: 7-pt text throughout (submission floor) needs ~10 pt of row pitch for panel b's
# 18 sites and ~22 pt of column pitch for the in-cell numbers, so the panel is taller and
# marginally narrower than the 9.2 x 2.9 in version that was tuned for 5-6 pt labels.
# The 0.1-in default white border of bbox_inches="tight" is 0.2 in of dead canvas width,
# which \includegraphics then pays for as a smaller scale factor; 0.02 in is enough.
plt.rcParams["savefig.pad_inches"] = 0.02
fig = plt.figure(figsize=(8.95, 3.5))
# 7 columns = 4 panels + 3 explicit spacer columns (widths in inches, wspace=0), sized at 7 pt from the
# rendered text: the a-b gap holds panel b's site labels (0.92 in), the b-c gap holds panel b's colourbar
# ticks+label AND panel c's y-axis label+ticks (0.38 + 0.45 in), the c-d gap only panel d's y axis (0.45 in).
gs = fig.add_gridspec(1, 7, width_ratios=[2.12, 1.00, 1.50, 0.93, 0.99, 0.54, 1.046], wspace=0, left=0.088, right=0.996, top=0.90, bottom=0.185)
a = fig.add_subplot(gs[0])
A_ROWS = [("rho_mtor_vs_translation_meas", "p_meas", "mTORC1 composite"),
          ("rho_S6arm_vs_translation_meas", "p_s6", "S6K1 arm (pRPS6)"),
          ("rho_4EBP1arm_vs_translation_meas", "p_bp", "4E-BP1 arm"),
          ("rho_mtor_vs_translation_HE", "p_HE", "composite vs H&E-only score"),
          ("rho_ERK_vs_translation_meas", "p_erk", "ERK (specificity control)")]
cohs_a = [c_ for c_ in COH if c_ in res.cohort.values]
ra = res.set_index("cohort").reindex(cohs_a)
R_a = ra[[rc for rc, _, _ in A_ROWS]].T.values.astype(float)        # 5 scores x 5 cohorts
P_a = ra[[pc for _, pc, _ in A_ROWS]].T.values.astype(float)
Q_a = bh_fdr(P_a.ravel()).reshape(P_a.shape)                       # BH across all 25 cells
STAR_a = np.isfinite(Q_a) & (Q_a < 0.05)
# carve a thin colourbar out of the right end of panel a's gridspec slot, so the a-b gap stays free for b's site labels
pa = a.get_position(); CBW, CBGAP, CBROOM = 0.0055, 0.005, 0.049       # figure fractions: bar width, gap, tick+label room
# CBROOM measured on the rendered PDF: at 7 pt the colourbar ticks + rotated label need 0.049 of the
# figure width; at the old 0.042 the label overran the slot and touched panel b's site labels.
a.set_position([pa.x0, pa.y0, pa.width - CBW - CBGAP - CBROOM, pa.height])
cax_a = fig.add_axes([pa.x0 + pa.width - CBW - CBROOM, pa.y0, CBW, pa.height])
im_a = a.imshow(R_a, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
for i in range(R_a.shape[0]):
    for j in range(R_a.shape[1]):
        v = R_a[i, j]
        if np.isfinite(v):
            a.text(j, i, f"{v:+.2f}" + ("*" if STAR_a[i, j] else ""), ha="center", va="center",
                   fontsize=7.0, color="white" if abs(v) > 0.35 else "#333")
a.set_xticks(range(len(cohs_a))); a.set_xticklabels(cohs_a, fontsize=7.0, rotation=30, ha="right", rotation_mode="anchor")
# two of the row labels are wrapped onto two lines at 7 pt so they stay inside the left margin;
# WRAP_A is display-only, the printed table below still uses the one-line names.
WRAP_A = {"mTORC1 composite": "mTORC1\ncomposite",
          "S6K1 arm (pRPS6)": "S6K1 arm\n(pRPS6)",
          "composite vs H&E-only score": "composite vs\nH&E-only score",
          "ERK (specificity control)": "ERK (specificity\ncontrol)"}
a.set_yticks(range(len(A_ROWS))); a.set_yticklabels([WRAP_A.get(lab, lab) for _, _, lab in A_ROWS], fontsize=7.0)
for s_ in a.spines.values():
    s_.set_visible(False)
cb_a = fig.colorbar(im_a, cax=cax_a)
cb_a.set_label(r"Spearman $\rho$ with translation residual", fontsize=7.0, labelpad=2); cb_a.ax.tick_params(labelsize=7.0, length=2, width=0.5)
a.set_title("Phospho-S6 tracks the residual", fontsize=7.6, fontweight="bold", loc="left")
a.text(-0.45, 1.06, "a", transform=a.transAxes, fontsize=10, fontweight="bold")
print("panel a: Spearman rho (rows = scores, cols = cohorts); * = BH q<0.05 across all 25 cells")
print(pd.DataFrame([[f"{R_a[i, j]:+.2f}{'*' if STAR_a[i, j] else ''}" for j in range(R_a.shape[1])] for i in range(R_a.shape[0])],
                   index=[lab for _, _, lab in A_ROWS], columns=cohs_a).to_string())

b = fig.add_subplot(gs[2])
piv = sites.pivot_table(index="site", columns="cohort", values="rho").reindex(columns=[c for c in COH if c in res.cohort.values])
im = b.imshow(piv.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
b.set_xticks(range(piv.shape[1])); b.set_xticklabels(piv.columns, fontsize=7.0, rotation=30, ha="right", rotation_mode="anchor")
b.set_yticks(range(piv.shape[0])); b.set_yticklabels([short_site(s) for s in piv.index], fontsize=7.0)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        if np.isfinite(v):
            b.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7.0, color="white" if abs(v) > 0.35 else "#333")
for s_ in b.spines.values():
    s_.set_visible(False)
cb = fig.colorbar(im, ax=b, fraction=0.05, pad=0.03); cb.set_label(r"$\rho$ with translation residual", fontsize=7.0, labelpad=2); cb.ax.tick_params(labelsize=7.0, length=2, width=0.5)
b.set_title("Per-site view", fontsize=7.6, fontweight="bold", loc="left")
b.text(-0.66, 1.06, "b", transform=b.transAxes, fontsize=10, fontweight="bold")

best = res.sort_values("rho_mtor_vs_translation_meas", ascending=False).iloc[0]["cohort"]
d = panels[best]; cax = fig.add_subplot(gs[4])
cax.scatter(d["mtor"], d["translation_meas"], s=9, color=S.ORG[best], alpha=0.7, edgecolors="none")
m = d["mtor"].notna() & d["translation_meas"].notna()
b1, b0 = np.polyfit(d.loc[m, "mtor"], d.loc[m, "translation_meas"], 1); xr = np.array([d["mtor"].min(), d["mtor"].max()])
cax.plot(xr, b1 * xr + b0, color=S.INK, lw=1.0)
r_, p_ = st.spearmanr(d.loc[m, "mtor"], d.loc[m, "translation_meas"])
cax.text(0.05, 0.93, f"{best}\n$\\rho$={r_:.2f}\nP={p_:.1e}", transform=cax.transAxes, fontsize=7.0, va="top")
cax.set_xlabel("mTORC1 phospho score", fontsize=7); cax.set_ylabel("measured translation residual", fontsize=7)
cax.set_title("Example cohort", fontsize=7.6, fontweight="bold", loc="left")
cax.text(-0.45, 1.06, "c", transform=cax.transAxes, fontsize=10, fontweight="bold")

dax = fig.add_subplot(gs[6])
if len(ribo_bg):
    rng = np.random.default_rng(0)
    cohs_d = [c_ for c_ in COH if c_ in ribo_bg.cohort.values]
    for i, c_ in enumerate(cohs_d):
        vals = ribo_bg.loc[ribo_bg.cohort == c_, "rho"].values
        jitter = rng.uniform(-0.32, 0.32, size=len(vals))
        dax.scatter(i + jitter, vals, s=4, color=S.LGREY, alpha=0.6, zorder=2, linewidths=0)
        dax.scatter(i, np.median(vals), marker="_", s=260, linewidths=1.6, color=S.INK, zorder=4)
        s6r = res.loc[res.cohort == c_, "rho_S6arm_vs_translation_meas"]
        if len(s6r) and pd.notna(s6r.iloc[0]):
            dax.scatter(i, s6r.iloc[0], marker="*", s=110, color="#0B7A6E", zorder=5,
                        edgecolors=S.INK, linewidths=0.4)
    dax.set_xticks(range(len(cohs_d))); dax.set_xticklabels(cohs_d, fontsize=7.0, rotation=30, ha="right", rotation_mode="anchor")
    dax.axhline(0, color="k", lw=0.5)
    dax.set_ylabel(r"$\rho$ with translation residual", fontsize=7)
    dax.set_title("Phospho-S6 vs. other\nribosomal sites", fontsize=7, fontweight="bold", loc="left")
    # legend note wrapped onto three lines: at 7 pt the one-line version was ~3.6 in wide and
    # ran off the canvas on both sides of this 1.06-in panel.
    dax.text(0.5, -0.16, "$\\bigstar$ phospho-S6 (S6K1 arm)\ngrey = other ribosomal sites\ndash = median",
             transform=dax.transAxes, fontsize=7.0, ha="center", va="top", color=S.GREY, linespacing=1.3)
else:
    dax.text(0.5, 0.5, "ribosomal-background\ncontrol pending", ha="center", va="center",
             fontsize=7.0, color=S.LGREY, transform=dax.transAxes, style="italic")
dax.text(-0.43, 1.06, "d", transform=dax.transAxes, fontsize=10, fontweight="bold")

S.save_pub(fig, "Fig_mtor")

BS = chr(92); EOL = BS + BS
L = [BS + "begin{table}[htbp]" + BS + "centering" + BS + "small",
     BS + "caption{Phospho-S6 tracks the translation residual. Per cohort: Spearman correlation of a canonical mTORC1 phospho-output score (RPS6 S235/236/240/244, 4E-BP1 T37/46/S65/T70, EIF4B S422, S6K1 T389; CPTAC phosphoproteomics) with the measured and the H"+BS+"&E-only translation residual scores and with the ER-secretion residual; an ERK phospho score (MAPK1/3 activation sites) as specificity control; and the Pearson correlation between the H"+BS+"&E-only and measured translation residual before and after partialling out the mTORC1 score. This co-variation does not by itself establish mTORC1 specificity -- see Table~"+BS+"ref{tab:mtorribo} and "+BS+"S"+BS+"ref{sec:mtorres}.}",
     BS + "label{tab:mtor}", BS + "resizebox{" + BS + "textwidth}{!}{" + BS + "begin{tabular}{lrrrrrrrrrr}" + BS + "toprule",
     "Cohort & $n$ & sites & composite & S6K1 arm (pRPS6) & 4E-BP1 arm & composite vs H"+BS+"&E score & vs secretion & ERK control & $r$(H"+BS+"&E, meas.) & given mTORC1 " + EOL + " " + BS + "midrule"]
for _, r in res.iterrows():
    L.append(f"{r.cohort} & {int(r.n)} & {int(r.n_mtor_sites)} & {r.rho_mtor_vs_translation_meas:+.2f} & {r.rho_S6arm_vs_translation_meas:+.2f} & {r.rho_4EBP1arm_vs_translation_meas:+.2f} & {r.rho_mtor_vs_translation_HE:+.2f} & {r.rho_mtor_vs_secretion_meas:+.2f} & {r.rho_ERK_vs_translation_meas:+.2f} & {r.r_HE_vs_meas:.2f} & {r.r_HE_vs_meas_given_mtor:.2f} " + EOL)
L += [BS + "bottomrule", BS + "end{tabular}}", BS + "end{table}"]
open("../SuppTable_mtor.tex", "w", encoding="utf-8").write(chr(10).join(L) + chr(10))

if len(ribo_bg) and len(ribo_partial):
    perc = {}
    for c_, g in ribo_bg.groupby("cohort"):
        s6r = res.loc[res.cohort == c_, "rho_S6arm_vs_translation_meas"]
        s6r = float(s6r.iloc[0]) if len(s6r) and pd.notna(s6r.iloc[0]) else np.nan
        perc[c_] = (100 * (g["rho"] < s6r).mean() if pd.notna(s6r) else np.nan, len(g),
                    g.rho.median(), g.rho.max(), g.loc[g.rho.idxmax(), "site"])
    rp = ribo_partial.set_index("cohort")
    L2 = [BS + "begin{table}[htbp]" + BS + "centering" + BS + "small",
          BS + "caption{Ribosomal-phosphosite background control (P0-4, adversarial review): is phospho-S6 "
          "special among phosphosites on ribosomal proteins, or typical of the class? For each site-level "
          "cohort (CCRCC excluded: gene-level phosphoproteome only): the number of non-RPS6 ribosomal "
          "phosphosites tested ($\\geq$50\\% patient coverage), their median Spearman $\\rho$ with the "
          "measured translation residual, the percentile rank of the S6K1 (pRPS6) arm within this "
          "background, the correlation of a ribosomal-proxy score (mean of these sites) with the residual, "
          "and the partial correlations of the S6K1 and 4E-BP1 arms controlling for this ribosomal proxy.}",
          BS + "label{tab:mtorribo}", BS + "begin{tabular}{lrrrrrr}" + BS + "toprule",
          "Cohort & sites & median $\\rho$ & S6-arm pctile & ribo-proxy $\\rho$ & S6-arm$|$ribo & 4E-BP1$|$ribo " + EOL + " " + BS + "midrule"]
    for c_ in [c for c in COH if c in perc]:
        pct, n_s, med, mx, mxsite = perc[c_]
        row = rp.loc[c_] if c_ in rp.index else None
        rr = f"{row.rho_ribo_proxy:+.2f}" if row is not None else "--"
        s6g = f"{row.rho_S6arm_given_riboproxy:+.2f}" if row is not None else "--"
        bpg = f"{row.rho_4EBP1arm_given_riboproxy:+.2f}" if row is not None and pd.notna(row.rho_4EBP1arm_given_riboproxy) else "--"
        L2.append(f"{c_} & {int(n_s)} & {med:+.2f} & {pct:.0f}\\% & {rr} & {s6g} & {bpg} " + EOL)
    L2 += [BS + "bottomrule", BS + "end{tabular}", BS + "end{table}"]
    open("../SuppTable_mtor_ribocontrol.tex", "w", encoding="utf-8").write(chr(10).join(L2) + chr(10))
    print("wrote SuppTable_mtor_ribocontrol.tex")

print("wrote Fig_mtor, mtor_results.csv, mtor_sites.csv, SuppTable_mtor.tex")
