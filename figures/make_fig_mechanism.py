#!/usr/bin/env python3
"""Composite Figure — generalisation and mechanism (Fig_mechanism).  All panels are re-drawn
from the saved outputs of the standalone scripts (no recomputation of anything that is stored):
  row 1  a, b   leave-one-cancer-out dumbbells        <- loco_generalization.py   (loco_results.csv)
  row 2  c      H&E-only score by molecular subtype   <- subtype_classify.py      (figdata/scores_*, figdata/subtype_*,
                                                                                   subtype_results.csv for the title stats)
  row 3  d, e, f mTORC1 phospho anchor                <- mtor_mechanism.py        (mtor_results.csv, mtor_sites.csv,
                                                                                   figdata/mtor_ribosomal_control.csv)
The plotting blocks are copied from those scripts (same data, encoding, colours, annotations); only the layout,
panel letters and sizes differ.  The standalone mtor panel c (LUAD example scatter) is deliberately omitted."""
import os, textwrap, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mrstyle as S

COH, ORGAN, ORG, INK, GREY, LGREY = S.COH, S.ORGAN, S.ORG, S.INK, S.GREY, S.LGREY
DD = "figdata"
FS_TITLE = 7.6


def lab(fig, ax, s, dx=0.0, dy=0.008):
    """bold 10-pt panel letter at the top-left of an axes, in figure coordinates (dx = shift left, fig fraction)."""
    bb = ax.get_position()
    fig.text(bb.x0 - dx, bb.y1 + dy, s, fontsize=10, fontweight="bold", ha="left", va="bottom")


fig = plt.figure(figsize=(7.2, 8.2))

# ============================ row 1: LOCO dumbbells (a, b) ============================
res = pd.read_csv("loco_results.csv")
gs1 = fig.add_gridspec(1, 2, left=0.085, right=0.985, top=0.965, bottom=0.81, wspace=0.12)
axes1 = [fig.add_subplot(gs1[0, 0]), fig.add_subplot(gs1[0, 1], sharey=None)]
tgts = [t for t in ["translation", "secretion_ER"] if t in set(res.target)]
PANEL = {"translation": "Translation residual programme", "secretion_ER": "ER-secretion residual programme"}
ypos = np.arange(len(COH))[::-1]                      # Kidney on top ... Pancreas at bottom
XLIM, XTICKS, OFF = (-0.22, 0.85), [-0.2, 0, 0.2, 0.4, 0.6, 0.8], 0.046
for ax, tgt in zip(axes1, tgts):
    tr = res[res.target == tgt].set_index("held_out").loc[COH]
    ax.axvline(0, color=INK, lw=0.6, ls=(0, (3, 2)), zorder=0)
    for y, c in zip(ypos, COH):
        wi, zs, col = float(tr.within_r[c]), float(tr.zero_shot_r[c]), ORG[c]
        ax.plot([wi, zs], [y, y], color=GREY, lw=1.0, solid_capstyle="round", zorder=1)   # connector = drop
        ax.scatter([zs], [y], s=30, facecolor="white", edgecolor=col, linewidths=1.1, zorder=2)  # open = zero-shot
        ax.scatter([wi], [y], s=30, facecolor=col, edgecolor=col, linewidths=0.8, zorder=3)    # filled = within-organ
        # if the two dots coincide (< 1 marker width apart) add a white halo so the ring stays visible on the disc
        if abs(zs - wi) * ax.get_position().width * fig.get_figwidth() * 72 / (XLIM[1] - XLIM[0]) < 6:
            ax.scatter([zs], [y], s=30, facecolor="none", edgecolor="white", linewidths=3.2, zorder=4)
        ax.scatter([zs], [y], s=30, facecolor="none", edgecolor=col, linewidths=1.1, zorder=5)   # ring drawn on top
        side = -1 if zs <= wi else 1                     # label on the side away from the filled dot
        ax.text(zs + side * OFF, y, f"{zs:.2f}".replace("-", "−"), ha="right" if side < 0 else "left", va="center",
                fontsize=7, color=INK, zorder=4, bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
    ax.set_xlim(*XLIM); ax.set_xticks(XTICKS)
    ax.set_ylim(-0.65, len(COH) - 0.35)
    ax.set_yticks(ypos); ax.set_yticklabels([ORGAN[c] for c in COH])
    ax.tick_params(axis="y", length=0)
    ax.set_title(PANEL[tgt], fontsize=FS_TITLE, fontweight="bold", loc="left", pad=4)
axes1[1].tick_params(axis="y", labelleft=False)
lab(fig, axes1[0], "a", dx=0.075); lab(fig, axes1[1], "b", dx=0.03)
fig.text((axes1[0].get_position().x0 + axes1[1].get_position().x1) / 2, 0.777,
         "Pearson $r$, predicted vs measured residual score", fontsize=7, ha="center", va="center")
handles = [Line2D([], [], marker="o", ls="none", ms=5, mfc=INK, mec=INK, label="filled = within-organ 5-fold CV"),
           Line2D([], [], marker="o", ls="none", ms=5, mfc="white", mec=INK, mew=1.0,
                  label="open = zero-shot (trained on the other 4 organs)")]
fig.legend(handles=handles, fontsize=7, frameon=False, ncol=2, loc="center",
           bbox_to_anchor=((axes1[0].get_position().x0 + axes1[1].get_position().x1) / 2, 0.758),
           handletextpad=0.4, columnspacing=1.8)

# ============================ row 2: subtype boxplot strips (c) ============================
# per-patient data rebuilt exactly as subtype_classify.py does (its CONFIG, norm_id and load_labels copied);
# the title statistics (KW P, AUROC, permutation P) are taken from subtype_results.csv, not recomputed.
CONFIG = {  # cohort: (label file, id column or None=index, subtype column, scheme name)
    "LUAD": (f"{DD}/subtype_luad.tsi", None, "NMF.consensus", "multi-omic NMF C1-C4"),
    "UCEC": (f"{DD}/subtype_ucec.txt", "Proteomics_Participant_ID", "Genomics_subtype", "TCGA genomic subtype"),
    "CCRCC": (f"{DD}/subtype_ccrcc.csv", None, "immune_subtype", "immune subtype (Clark 2019)"),
    "GBM": (f"{DD}/subtype_gbm.csv", None, "nmf_subtype", "nmf1/2/3 (Wang 2021)"),
    "PDAC": (f"{DD}/subtype_pdac.csv", None, "subtype", "Moffitt classical/basal-like (Cao 2021 annotation)"),
}


def norm_id(s):
    s = str(s).strip().replace(".", "-")
    return "-".join(s.split("-")[:2]) if s.startswith("C3") else s


def load_labels(path, id_col, sub_col):
    sep = "\t" if path.endswith((".tsi", ".txt", ".tsv")) else ","
    df = pd.read_csv(path, sep=sep, index_col=0 if id_col is None else None, dtype=str, encoding="latin-1")
    if sub_col not in df.columns and sub_col in df.index:      # LinkedOmics transposed layout
        df = df.T
    if sub_col not in df.columns:
        raise SystemExit(f"{path}: subtype column '{sub_col}' not found. columns={list(df.columns)[:15]} index={list(df.index)[:8]}")
    ids = df.index if id_col is None else df[id_col]
    return pd.Series(df[sub_col].values, index=[norm_id(i) for i in ids]).dropna()


sres = pd.read_csv("subtype_results.csv").set_index("cohort")
panels = []
for c, (path, id_col, sub_col, scheme) in CONFIG.items():
    if not os.path.exists(path) or c not in sres.index:
        continue
    labs_ = load_labels(path, id_col, sub_col)
    sc = pd.read_csv(f"{DD}/scores_{c.lower()}.csv", index_col=0)
    sc.index = [norm_id(i) for i in sc.index]
    df = sc.join(labs_.rename("subtype"), how="inner").dropna(subset=["subtype"])
    df = df[df["subtype"].astype(str).str.lower().isin(["nan", "na", ""]) == False]
    vc = df["subtype"].value_counts(); keep = vc[vc >= 5].index
    df = df[df["subtype"].isin(keep)]
    if len(keep) < 2:
        continue
    r = sres.loc[c]
    best = f"{r.best_score}_morph"
    order = sorted(keep, key=str)
    panels.append(dict(c=c, scheme=scheme, best=best, order=order, kw=float(r.kw_p_best), auc=float(r.cv_auroc),
                       p=float(r.auroc_perm_p), n=len(df),
                       data=[df.loc[df.subtype == s, best].dropna().values for s in order]))
    assert len(df) == int(r.n), f"{c}: rebuilt n={len(df)} != subtype_results n={int(r.n)}"

# wspace 0.62 -> 0.78: at the 7 pt floor the CCRCC immune-subtype labels still touch
gs2 = fig.add_gridspec(1, len(panels), left=0.075, right=0.985, top=0.675, bottom=0.535, wspace=0.78)
axes2 = [fig.add_subplot(gs2[0, i]) for i in range(len(panels))]
for ax, P in zip(axes2, panels):
    c = P["c"]
    bp = ax.boxplot(P["data"], widths=0.55, showfliers=False, patch_artist=True, medianprops=dict(color=INK, lw=1.2))
    for b in bp["boxes"]:
        b.set_facecolor(ORG[c]); b.set_alpha(0.35); b.set_edgecolor(ORG[c])
    for i, d in enumerate(P["data"]):
        ax.scatter(np.random.RandomState(i).normal(i + 1, 0.06, len(d)), d, s=6, color=ORG[c], alpha=0.65, edgecolors="none")
    # Abbreviate aggressively and rotate to 45 deg: at the 6.5 pt floor these labels
    # need for legibility at final print size, the CCRCC immune-subtype names collide
    # at 30 deg in a ~1.1 in wide mini-panel. Full names are given in the caption.
    # At the 7 pt print floor even "CD8- infl." / "Metab. desert" touch at 45 deg in a
    # ~1.1 in mini-panel, so the CCRCC immune subtypes are reduced to their
    # distinguishing token; the caption carries the full names.
    SHORT = {"CD8+ Inflamed": "CD8+", "CD8- Inflamed": "CD8−",
             "Metabolic Immune Desert": "Metab.", "VEGF Immune Desert": "VEGF"}
    labs = [SHORT.get(s, s.replace(" Immune ", " Imm. ").replace("Metabolic", "Metab."))
            for s in P["order"]]
    ax.set_xticks(range(1, len(labs) + 1))
    ax.set_xticklabels(labs, fontsize=7, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel(f"H&E-only {P['best'].replace('_morph', '')} residual", fontsize=7)
    sig = "*" if P["p"] < 0.05 else ""
    # same title text as the standalone; the "KW ... AUROC ..." line is broken at the separator because the
    # composite mini-panels are ~1.1 in wide (the standalone's are 2.2 in)
    scheme = P["scheme"].split(" (")[0]
    # composite-only display shortening: "PDAC: Moffitt classical/basal-like" is 1.36 in at 6.2 pt bold but the rightmost
    # mini-panel has only 0.98 in to the page edge, so it wrapped to a 4th title line; the two class names are already
    # this panel's x tick labels, so nothing is lost
    scheme = {"Moffitt classical/basal-like": "Moffitt subtype"}.get(scheme, scheme)
    head = textwrap.fill(f"{c}: {scheme}", 30, break_on_hyphens=False)  # all five heads now fit on one line
    ax.set_title(f"{head}\nKW P={P['kw']:.1g} ·\nAUROC={P['auc']:.2f} (P={P['p']:.3f}){sig}",
                 fontsize=7, fontweight="bold", loc="left", linespacing=1.15)
    ax.tick_params(axis="y", labelsize=7)
lab(fig, axes2[0], "c", dx=0.065, dy=0.058)
fig.text(0.075, 0.733, "Does the H&E-only residual score track published molecular subtypes?",
         fontsize=FS_TITLE, fontweight="bold", ha="left", va="bottom")

# ============================ row 3: mTORC1 mechanism (d, e, f) ============================
def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (q) for a flat array; NaN entries stay NaN."""
    p = np.asarray(p, float); q = np.full(p.shape, np.nan)
    m = np.isfinite(p); pv = p[m]; n = pv.size
    if n:
        o = np.argsort(pv); r = pv[o] * n / np.arange(1, n + 1)
        r = np.minimum(np.minimum.accumulate(r[::-1])[::-1], 1.0)
        qq = np.empty(n); qq[o] = r; q[m] = qq
    return q


mres = pd.read_csv("mtor_results.csv"); sites = pd.read_csv("mtor_sites.csv")
ribo_bg = pd.read_csv(f"{DD}/mtor_ribosomal_control.csv")
# 5 slots = 3 panels + 2 spacer columns (wspace=0); the d-e gap holds e's long site labels, the e-f gap holds
# e's colourbar plus f's y-axis label.
gs3 = fig.add_gridspec(1, 5, width_ratios=[1.85, 0.85, 1.00, 0.85, 1.15], wspace=0,
                       left=0.13, right=0.985, top=0.425, bottom=0.075)

# ---- d : score x cohort Spearman-rho heatmap, BH-starred (= standalone panel a) ----
a = fig.add_subplot(gs3[0])
A_ROWS = [("rho_mtor_vs_translation_meas", "p_meas", "mTORC1 composite"),
          ("rho_S6arm_vs_translation_meas", "p_s6", "S6K1 arm (pRPS6)"),
          ("rho_4EBP1arm_vs_translation_meas", "p_bp", "4E-BP1 arm"),
          ("rho_mtor_vs_translation_HE", "p_HE", "composite vs\nH&E-only score"),   # same label text as the
          ("rho_ERK_vs_translation_meas", "p_erk", "ERK (specificity\ncontrol)")]     # standalone, wrapped to 2 lines
cohs_a = [c_ for c_ in COH if c_ in mres.cohort.values]
ra = mres.set_index("cohort").reindex(cohs_a)
R_a = ra[[rc for rc, _, _ in A_ROWS]].T.values.astype(float)        # 5 scores x 5 cohorts
P_a = ra[[pc for _, pc, _ in A_ROWS]].T.values.astype(float)
Q_a = bh_fdr(P_a.ravel()).reshape(P_a.shape)                       # BH across all 25 cells
STAR_a = np.isfinite(Q_a) & (Q_a < 0.05)
pa = a.get_position(); CBW, CBGAP, CBROOM = 0.0075, 0.006, 0.050       # figure fractions: bar width, gap, tick+label room
a.set_position([pa.x0, pa.y0, pa.width - CBW - CBGAP - CBROOM, pa.height])
cax_a = fig.add_axes([pa.x0 + pa.width - CBW - CBROOM, pa.y0, CBW, pa.height])
im_a = a.imshow(R_a, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
for i in range(R_a.shape[0]):
    for j in range(R_a.shape[1]):
        v = R_a[i, j]
        if np.isfinite(v):
            a.text(j, i, f"{v:+.2f}" + ("*" if STAR_a[i, j] else ""), ha="center", va="center",
                   fontsize=7, color="white" if abs(v) > 0.35 else "#333")
a.set_xticks(range(len(cohs_a))); a.set_xticklabels(cohs_a, fontsize=7, rotation=30, ha="right")
a.set_yticks(range(len(A_ROWS))); a.set_yticklabels([l for _, _, l in A_ROWS], fontsize=7)
for s_ in a.spines.values():
    s_.set_visible(False)
cb_a = fig.colorbar(im_a, cax=cax_a)
cb_a.set_label(r"Spearman $\rho$ with translation residual", fontsize=7, labelpad=2)
cb_a.ax.tick_params(labelsize=7, length=2, width=0.5)
a.set_title("Phospho-S6 tracks the residual", fontsize=FS_TITLE, fontweight="bold", loc="left")
lab(fig, a, "d", dx=0.12)

# ---- e : per-site heatmap (= standalone panel b) ----
b = fig.add_subplot(gs3[2])
piv = sites.pivot_table(index="site", columns="cohort", values="rho").reindex(columns=[c for c in COH if c in mres.cohort.values])
im = b.imshow(piv.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
b.set_xticks(range(piv.shape[1])); b.set_xticklabels(piv.columns, fontsize=7, rotation=30, ha="right")
b.set_yticks(range(piv.shape[0])); b.set_yticklabels(piv.index, fontsize=7)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        if np.isfinite(v):
            b.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7, color="white" if abs(v) > 0.35 else "#333")
for s_ in b.spines.values():
    s_.set_visible(False)
pb = b.get_position()
cax_b = fig.add_axes([pb.x1 + 0.006, pb.y0 + 0.25 * pb.height, 0.0075, 0.5 * pb.height])
cb = fig.colorbar(im, cax=cax_b); cb.set_label(r"$\rho$ with translation residual", fontsize=7, labelpad=2)
cb.ax.tick_params(labelsize=7, length=2, width=0.5)
b.set_title("Per-site view", fontsize=FS_TITLE, fontweight="bold", loc="left")
lab(fig, b, "e", dx=0.115)

# ---- f : ribosomal-phosphosite background strip (= standalone panel d) ----
dax = fig.add_subplot(gs3[4])
rng = np.random.default_rng(0)
cohs_d = [c_ for c_ in COH if c_ in ribo_bg.cohort.values]
for i, c_ in enumerate(cohs_d):
    vals = ribo_bg.loc[ribo_bg.cohort == c_, "rho"].values
    jitter = rng.uniform(-0.32, 0.32, size=len(vals))
    dax.scatter(i + jitter, vals, s=4, color=LGREY, alpha=0.6, zorder=2, linewidths=0)
    dax.scatter(i, np.median(vals), marker="_", s=260, linewidths=1.6, color=INK, zorder=4)
    s6r = mres.loc[mres.cohort == c_, "rho_S6arm_vs_translation_meas"]
    if len(s6r) and pd.notna(s6r.iloc[0]):
        dax.scatter(i, s6r.iloc[0], marker="*", s=110, color="#0B7A6E", zorder=5, edgecolors=INK, linewidths=0.4)
dax.set_xticks(range(len(cohs_d))); dax.set_xticklabels(cohs_d, fontsize=7)
dax.axhline(0, color="k", lw=0.5)
dax.set_ylabel(r"$\rho$ with translation residual", fontsize=7)
dax.set_title("Phospho-S6 vs. other\nribosomal sites", fontsize=FS_TITLE, fontweight="bold", loc="left")   # 2 lines, same size/weight
# as d/e; "ribosomal sites (non-mTORC1)" is 1.56 in at 7.6 pt bold vs 1.35 in from f's left edge to the page edge, so the
# non-mTORC1 qualifier is carried only by the key below the panel
dax.text(0.5, -0.115, r"$\bigstar$ phospho-S6 (S6K1 arm)" "\n" "grey = non-mTORC1 ribosomal sites" "\n" "dash = median",
         transform=dax.transAxes, fontsize=7, ha="center", va="top", color=GREY, linespacing=1.3)
lab(fig, dax, "f", dx=0.075)

S.save_pub(fig, "Fig_mechanism")
print("wrote Fig_mechanism")
print("panel d (rows = scores, cols = cohorts; * = BH q<0.05 across 25 cells)")
print(pd.DataFrame([[f"{R_a[i, j]:+.2f}{'*' if STAR_a[i, j] else ''}" for j in range(R_a.shape[1])] for i in range(R_a.shape[0])],
                   index=[l for _, _, l in A_ROWS], columns=cohs_a).to_string())
print("panel e"); print(piv.round(2).to_string())
for P in panels:
    print(f"panel c {P['c']}: n={P['n']} order={P['order']} sizes={[len(d) for d in P['data']]} "
          f"KW P={P['kw']:.1g} AUROC={P['auc']:.2f} P={P['p']:.3f}")
for c_ in cohs_d:
    v = ribo_bg.loc[ribo_bg.cohort == c_, "rho"]
    print(f"panel f {c_}: n_sites={len(v)} median={v.median():+.3f} S6-arm={float(mres.set_index('cohort').loc[c_, 'rho_S6arm_vs_translation_meas']):+.3f}")
