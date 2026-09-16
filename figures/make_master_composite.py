#!/usr/bin/env python3
"""Master composite Figure 1 — the finding. Page-filling multi-panel figure, all native,
all OUR real data.  a workflow | b motivation (mRNA explains little) | c ranked counts |
d incr-R2 distribution | e UNI-vs-Phikon robustness | f embedding t-SNE | g morphology
recovers the residual."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy import stats as st
from sklearn.manifold import TSNE
import mrstyle as S

ORG, INK, GREY, LGREY = S.ORG, S.INK, S.GREY, S.LGREY
COH, ORGAN = S.COH, S.ORGAN
DD = "figdata"
RBASE = "../server_export/pinned"
RPATH = {c: f"{RBASE}/{c.lower()}/residual_results_tumoronly.csv" for c in COH}
# (was D:/claude/pajinsen/MorphoResidual/pilot_ccrcc/results_server/... -- a pre-dates-this-repo
# pilot path; the LUAD/UCEC/GBM/PDAC files there no longer exist (script would crash if rerun),
# and the CCRCC file that does still exist there is byte-identical (sha256 verified 2026-09-17)
# to server_export/pinned/ccrcc/residual_results_tumoronly.csv, so this is a path fix only --
# the CCRCC number this panel already reported was correct.)


def lab(ax, s, x=-0.16, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


incr = {c: pd.read_csv(f"{DD}/merged_incr_{c.lower()}.csv") for c in COH}
sig = {c: d[(d.fdr_uni < 0.05) & (d.incremental_r2_uni > 0)] for c, d in incr.items()}
counts = {c: len(sig[c]) for c in COH}
tested = {c: int(incr[c]["incremental_r2_uni"].notna().sum()) for c in COH}
medincr = {c: float(sig[c]["incremental_r2_uni"].median()) for c in COH}
order = sorted(COH, key=lambda c: counts[c], reverse=True)
r2rna = {c: pd.read_csv(p)["r2_rna"].values for c, p in RPATH.items()}
pooled = np.concatenate([v[np.isfinite(v)] for v in r2rna.values()])

fig = plt.figure(figsize=(7.2, 8.2))
gs = fig.add_gridspec(4, 2, height_ratios=[0.66, 1.0, 1.0, 1.0], hspace=0.5, wspace=0.30,
                      left=0.09, right=0.965, top=0.985, bottom=0.058)

# ===== a : workflow =====
axa = fig.add_subplot(gs[0, :]); axa.imshow(mpimg.imread("Fig1_flagship.png"), aspect="auto")
axa.axis("off"); lab(axa, "a", x=-0.02, y=1.02)

# ===== b : motivation (mRNA explains little) =====
axb = fig.add_subplot(gs[1, 0]); lab(axb, "b")
axb.hist(np.clip(pooled, -0.1, 1.0), bins=52, range=(-0.1, 1.0), color="#9ecae1",
         edgecolor="white", linewidth=0.3)
med = np.median(pooled); flow = 100 * np.mean(pooled < 0.25)
axb.axvspan(-0.1, 0.25, color=ORG["PDAC"], alpha=0.10, zorder=0)
axb.axvline(med, color=INK, lw=1.1, ls="--")
axb.text(med + 0.03, axb.get_ylim()[1] * 0.9, f"median\nR²={med:.2f}", fontsize=7, color=INK)
axb.text(0.30, axb.get_ylim()[1] * 0.6, f"{flow:.0f}% of pairs\nR² < 0.25",
         fontsize=7, color=ORG["PDAC"])
axb.set_xlabel("mRNA" + "→" + "protein R² (per gene" + "–" + "cohort pair)", fontsize=7)
axb.set_ylabel("gene" + "–" + "cohort pairs", fontsize=7); axb.set_xlim(-0.1, 1.0)
axb.set_title("mRNA explains little of the proteome", fontsize=7.6, fontweight="bold", loc="left")

# ===== c : bubble scatter — effect size x breadth x count, per organ =====
axc = fig.add_subplot(gs[1, 1]); lab(axc, "c")
pct = {c: 100 * counts[c] / tested[c] for c in COH}
SK = 0.28                                   # pt^2 per protein (bubble area ∝ count)
for c in COH:
    axc.scatter(medincr[c], pct[c], s=counts[c] * SK, color=ORG[c], alpha=0.85,
                edgecolors="white", linewidths=0.5, zorder=3)
# label offsets (points) chosen so LUAD/CCRCC (adjacent) and UCEC never collide
off = {"LUAD": (-16, 0, "right", "center"), "CCRCC": (-4, -19, "left", "top"),
       "UCEC": (0, -19, "center", "top"), "PDAC": (0, -18, "left", "top"),
       "GBM": (11, 0, "left", "center")}
for c in COH:
    dx, dy, ha, va = off[c]
    axc.annotate(f"{ORGAN[c]}  {counts[c]:,}", (medincr[c], pct[c]), xytext=(dx, dy),
                 textcoords="offset points", ha=ha, va=va, fontsize=7, color=INK, zorder=4)
axc.set_xlim(0.04, 0.145); axc.set_ylim(0, 30)
axc.set_xlabel(r"median incremental R² (significant set)", fontsize=7)
axc.set_ylabel("% of tested proteins significant" + chr(10) + "(FDR<0.05)", fontsize=7)
axc.set_title("Breadth and effect size per organ", fontsize=7.6, fontweight="bold", loc="left")
# bubble-size key (bottom right, empty region)
for kx, kn in [(0.121, 500), (0.134, 2500)]:
    axc.scatter(kx, 6.0, s=kn * SK, color="none", edgecolors=GREY, linewidths=0.6, zorder=2)
    axc.text(kx, 1.2, f"{kn:,}", ha="center", va="center", fontsize=7, color=GREY)
axc.text(0.1255, 11.0, "significant proteins", ha="center", va="bottom", fontsize=7, color=GREY)

# ===== d : incremental-R2 distribution =====
axd = fig.add_subplot(gs[2, 0]); lab(axd, "d")
data = [sig[c]["incremental_r2_uni"].values for c in order]
vp = axd.violinplot(data, showextrema=False, widths=0.85)
for b, c in zip(vp["bodies"], order):
    b.set_facecolor(ORG[c]); b.set_alpha(0.75); b.set_edgecolor("none")
for i, c in enumerate(order):
    axd.plot([i + 1 - 0.28, i + 1 + 0.28], [medincr[c]] * 2, color=INK, lw=1.1, zorder=5)
axd.set_xticks(range(1, 6)); axd.set_xticklabels([ORGAN[c] for c in order], fontsize=7)
axd.set_ylabel(r"incremental R² (significant)", fontsize=7)
axd.set_ylim(0, np.percentile(np.concatenate(data), 99))

# ===== e : UNI vs Phikon agreement =====
axe = fig.add_subplot(gs[2, 1]); lab(axe, "e")
allu = np.concatenate([incr[c]["incremental_r2_uni"].values for c in COH])
allp = np.concatenate([incr[c]["incremental_r2_phi"].values for c in COH])
m = np.isfinite(allu) & np.isfinite(allp); allu, allp = allu[m], allp[m]
axe.hexbin(allu, allp, gridsize=45, bins="log", cmap="Greys", mincnt=1, extent=(-0.05, 0.35, -0.05, 0.35))
rr = st.pearsonr(allu, allp)[0]
axe.plot([-0.05, 0.35], [-0.05, 0.35], color="#B22222", lw=0.9, ls="--")
axe.text(0.05, 0.92, f"r = {rr:.2f}", transform=axe.transAxes, fontsize=7, color="#B22222")
axe.set_xlabel(r"incremental R² (UNI)", fontsize=7); axe.set_ylabel(r"incremental R² (Phikon)", fontsize=7)
axe.set_xlim(-0.05, 0.35); axe.set_ylim(-0.05, 0.35)

# ===== f : embedding t-SNE =====
axf = fig.add_subplot(gs[3, 0]); lab(axf, "f")
emb = pd.read_csv(f"{DD}/emb_pca50.csv")
Y = TSNE(n_components=2, random_state=0, init="pca", perplexity=30,
         learning_rate="auto").fit_transform(emb[[f"pc{i}" for i in range(1, 51)]].values)
for c in COH:
    mm = emb["cohort"].values == c.lower()
    axf.scatter(Y[mm, 0], Y[mm, 1], s=3, color=ORG[c], alpha=0.7, edgecolors="none", label=ORGAN[c])
axf.set_xticks([]); axf.set_yticks([]); axf.set_xlabel("t-SNE 1", fontsize=7); axf.set_ylabel("t-SNE 2", fontsize=7)
axf.legend(fontsize=7, markerscale=2.0, loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=5,
           handletextpad=0.2, borderpad=0.2, columnspacing=0.9, frameon=False)

# ===== g : morphology recovers the residual =====
axg = fig.add_subplot(gs[3, 1]); lab(axg, "g")
xs, ys, cs = [], [], []
for c in COH:
    s = pd.read_csv(f"{DD}/scores_{c.lower()}.csv", index_col=0)
    xs.append(s["translation_meas"].values); ys.append(s["translation_morph"].values); cs += [ORG[c]] * len(s)
xs = np.concatenate(xs); ys = np.concatenate(ys)
axg.scatter(xs, ys, s=6, c=cs, alpha=0.65, edgecolors="none")
r2 = st.pearsonr(xs, ys); b1, b0 = np.polyfit(xs, ys, 1); xr = np.array([xs.min(), xs.max()])
axg.plot(xr, b1 * xr + b0, color=INK, lw=1.1)
axg.text(0.05, 0.92, f"r = {r2[0]:.2f}\nP = {r2[1]:.0e}", transform=axg.transAxes, fontsize=7, color=INK, va="top")
axg.set_xlabel("measured translation residual\n(from proteomics)", fontsize=7)
axg.set_ylabel("morphology-only score\n(from H&E)", fontsize=7)

S.save_pub(fig, "Fig_master")
print(f"wrote Fig_master | median r2_rna={med:.3f}, %<0.25={flow:.0f}")
