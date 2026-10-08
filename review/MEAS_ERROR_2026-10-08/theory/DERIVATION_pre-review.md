# Own-mRNA measurement error vs post-transcriptional signal: what the estimator measures, and a signed discriminator

MorphoResidual, Limitations (v). Statistical design note, 2026-10-08. Synthetic-data checks in
`sim_population_check.py` (closed forms), `sim_signtest.py` (n = 100, 20 PCs, joint permutation
null; grids `sim_grid_v2_*.csv`), `sim_mixture.py` (null + signal gene mixtures, replicate seeds;
`sim_mixture_main.csv`, `sim_mixture_size.csv`) and `sim_batchnull.py` (operator-batch structure,
unrestricted vs within-operator null; `sim_batchnull.csv`); tables in `section7_tables.md`
(`sim_section7_tables.py`) and section 7. No project data were used anywhere in this folder; the
only external input is public GDC file metadata (counts of RNA files per case, section 7.5).

中文摘要：在 X（真表达）、M = X + e（测得 mRNA，可靠性 λ）、P = bX + u + ε、W = aX + cu + 噪声 的生成模型下，
(1) 发表的增量 R² 在 H_noise 下收敛到 R²_rna·(1−λ)²ρ_W/(λ(1−ρ_Wλ))（≤ 稀释上界 R²_rna(1−λ)/λ），在 H_post 下收敛到
ρ_U·σ_u²/σ_P²，与 λ 无关；(2) 镜像增量在两种假设下都为正，不能判别；F18 的符号检验 corr(r̂_p, r_m) 其实等价于
把一个可识别量 θ = β·cov(r̂_p, M)/cov(r̂_p, r_p) 与 R²_rna 比较（θ 在 H_noise 下 = λ/(1−λ)，H_post 下 = 0，
"形态读到 mRNA 测量误差"下 = −1，混合下 = f_X·λ/(1−λ)），所以负号是很强的 post 证据，正号却很弱；
应直接报告 θ（及由它推出的 λ̂_noise = θ/(1+θ)）；(3) 用嵌套模型增量 d_p 做符号检验在任何假设下都机械为负，必须用
仅含 W 的岭回归预测 r̂_p；联合置换零分布须保留 (M,P) 配对、折、岭参数，只打乱 W 的行；(4) 模拟（第 7 节）：n=100、20 PC、
B=200 时，集合层面 N_c 统计量在 H_noise 下 z ≥ +10、H_err 下 z ≤ −9，在 H_post 下 z 落在 [−1.5, +1.9]（48 个模拟集合无一 ≥ 2）；
选择（赢家诅咒）把 θ_set 压到总体值的 0.35–0.95 倍，故读数规则按"同样选择下的 H_noise 预测值"校准（§7.6）；
(5) 主零分布定为**操作员（批次）内置换**：模拟中若 RNA 批次效应可被形态经操作员身份读出，无限制置换会把 H_post 误判为
H_err（7/8–8/8 次），而操作员内置换恢复校准（§7.4）；无限制置换只作敏感性分析。

---

## 0. Notation and the published estimator

Per gene, patients i = 1..n (one row per patient; the 20 morphology PCs W_i are shared across genes).

| symbol | meaning |
|---|---|
| X | true transcript abundance (centred; var X = σ_X², set to 1 WLOG) |
| M = X + e | measured mRNA (one RNA-seq library); var e = σ_e²; reliability λ = var X / var M = 1/(1+σ_e²) |
| P = bX + u + ε | protein; u = post-transcriptional component (⊥ X), ε = proteomic noise; σ_P² = b² + σ_u² + σ_ε² |
| W = aX + cu + η | morphology PC block (K = 20); a, c ∈ ℝ^K; η ⊥ (X, u, e, ε); Σ_W = cov W |
| ρ² = R²_rna | population mRNA→protein R² = cov(M,P)²/(var M var P) = b²λ/σ_P² |
| β = cov(P,M)/var M = bλ | OLS slope of P on M; r_p = P − βM |
| δ = cov(M,P)/var P = b/σ_P² | OLS slope of M on P; r_m = M − δP |
| H_noise | c = 0: morphology reads X only; the beyond-mRNA signal exists only because λ < 1 |
| H_post | a = 0: morphology reads u only |
| H_err | morphology reads the mRNA measurement error itself: W = d·e + η (shared RNA-quality artefact) |

Published statistic: Δ̂_p = OOF-R²(P | M, W) − OOF-R²(P | M), 5-fold, ridge α = 1 on in-fold standardised
features, patient↔slide permutation of W for the null. Its population target is the OLS increment
Δ_p = C_p^⊤ Σ_{W⊥M}^{-1} C_p / σ_P², with C_p = cov(W, r_p) and Σ_{W⊥M} = cov of W residualised on M.

Two identities used throughout (exact for in-sample OLS residuals, O(1/n) otherwise):

- (I1) r_m = (1 − ρ²)·M − δ·r_p. Proof: P = r_p + βM ⇒ r_m = M − δr_p − δβM and δβ = ρ².
- (I2) corr(r_p, r_m) = −ρ. (cov = −cov(M,P)(1−ρ²), var r_p = σ_P²(1−ρ²), var r_m = var M (1−ρ²).)

(I1) is the key fact: the mirror residual carries no information beyond (M, r_p). Everything a
"mirror" analysis can see is a linear combination of cov(W, M) and cov(W, r_p).

## 1. What the incremental R² converges to

Residuals under the model (using β = bλ, σ_e² = (1−λ)/λ):

    r_p = b(1−λ)X + u + ε − bλe
    r_m = (1 − ρ²/λ)X + e − δu − δε               [1 − δb = (σ_u²+σ_ε²)/σ_P² = 1 − ρ²/λ]

Covariances with the morphology block (vector form; · is the scalar in the one-direction case):

    C_p := cov(W, r_p) = a·b(1−λ) + c·σ_u²                  (the e-term: cov(aX+cu+η, −bλe) = 0)
    C_M := cov(W, M)   = a
    C_m := cov(W, r_m) = (1−ρ²)·a − δ·C_p                    (from I1)

Write ρ_W := a^⊤Σ_W^{-1}a (R² of the oracle linear prediction of X from W; "morphology as a proxy for
the true transcript") and ρ_U := σ_u² c^⊤Σ_W^{-1}c (R² of u from W). With isotropic η and one signal
direction (the case the simulation uses; the general case only replaces scalars by quadratic forms):

**H_noise (c = 0).**

    Δ_p^noise = ρ² · (1−λ)² ρ_W / ( λ (1 − ρ_W λ) )          ≤ ρ² (1−λ)/λ   (equality at ρ_W = 1)

The right-hand bound is the "dilution bound" of review/recovered_2026-09-27_dilution_bound (λ* =
r²/(r²+Δ)). Two features matter: it is proportional to R²_rna, and it vanishes as (1−λ)² as λ → 1.
Reaching Δ_p = 0.10 with R²_rna = 0.2 needs λ ≤ 0.67 even with a perfect morphological proxy of the
true transcript; Δ_p = 0.20 needs λ ≤ 0.5.

**H_post (a = 0).**

    Δ_p^post = ρ_U · σ_u² / σ_P²                              (independent of λ and of R²_rna)

**Mixed.** With a ⊥ c in W-space (so that the two channels add in quadrature),

    Δ_p^mix = [ b²(1−λ)² ρ_W/(1−ρ_W λ) + σ_u² ρ_U ] / σ_P²,

and we define the X-share f_X := (X-channel term)/(total). If a and c are not orthogonal the cross
term 2ab(1−λ)·cσ_u² (a·c direction cosine) is added; it can have either sign.

**H_err** (W = d·e + η, morphology reads a shared mRNA artefact): C_p = −bλ·d σ_e² ⇒ Δ_p^err =
ρ_E · b² λ(1−λ)/σ_P² = ρ_E · ρ²(1−λ), with ρ_E the R² of e from W. Also vanishes as λ → 1.

**Mirror increment** Δ_m = C_m^⊤ Σ_{W⊥P}^{-1} C_m / var M. One-direction forms (var M = 1/λ):

    H_noise: Δ_m = λ ρ_W (1 − ρ²/λ)² / (1 − ρ_W ρ²/λ)   → large (morphology predicts the mRNA residual well)
    H_post : Δ_m = λ ρ_U σ_u² b²/σ_P⁴ / (1 − ρ_U σ_u²/σ_P²) ≈ ρ² · Δ_p^post   → small but POSITIVE
    H_err  : Δ_m = ρ_E (1−λ)

So the "reverse incremental R² symmetry" fails: Δ_m > 0 under both H_noise and H_post, because r_m
contains −δu (I1). A positive mirror increment is not evidence that morphology proxies the transcript.
The *ratio* Δ_m/Δ_p is (λ−ρ²)²/((1−λ)²ρ²) under H_noise and ≈ ρ² under H_post — unsigned, and the
H_noise value depends on the unknown λ, so it is a weak discriminator (it is reported as a secondary
quantity in the simulation only).

Finite-sample note. With K = 20 PCs, n ≈ 100 and ridge α = 1 on standardised features (effectively
OLS), the OOF increment carries a capacity penalty of about −(K/n_train)(1 − R²_full) ≈ −0.15 to −0.20,
so Δ̂_p ≈ Δ_p − 0.17 and the published permutation null is centred below zero. Genes called significant
therefore have population increments of roughly 0.1–0.3, not 0.03–0.1; the simulation grid uses
Δ_p ∈ {0.10, 0.20}. (This also means H_noise must explain increments of that size, which under the
dilution bound requires λ ≲ 0.5–0.7 together with ρ_W ≈ 1.) A second finite-sample effect works
against H_noise specifically: when ρ_W ≈ 1 the block W is nearly collinear with M, the nested [M, W]
fit has inflated variance and the OOF increment is penalised more (simulation: −0.33 at population
Δ_p = 0.10, λ = 0.6, versus −0.18 under H_post at the same Δ_p), so H_noise genes are *harder* to select
than H_post genes with the same population increment. The published selection is therefore biased, if
anything, against H_noise genes — a point in the paper's favour that is independent of the sign test.

## 2. A signed discriminator; why the mirror sign test is one-sided; what the null must preserve

### 2.1 The identified quantity

Let r̂_p be the OOF ridge prediction of r_p from the **raw** W block (W only, not residualised on M,
not the nested increment). In population r̂_p = w^⊤W with w = Σ_W^{-1}C_p, hence

    cov(r̂_p, M)   = w^⊤a
    cov(r̂_p, r_p) = w^⊤C_p = w^⊤a · b(1−λ) + w^⊤c · σ_u²

Define the X-share of the prediction's covariance, f_X := w^⊤a·b(1−λ) / w^⊤C_p (∈ [0,1] when both
channels pull the same way; equals the X-share of Δ_p in the orthogonal case), and the statistic

    θ := β · cov(r̂_p, M) / cov(r̂_p, r_p)            ("own-transcript alignment")

Using β = bλ:

    θ = f_X · λ/(1−λ)        (H_noise: f_X = 1 ⇒ θ = λ/(1−λ);  H_post: f_X = 0 ⇒ θ = 0)
    θ = −1                   (H_err, exactly: r̂_p ∝ −e, cov(r̂_p,M) = −κσ_e², cov(r̂_p,r_p) = κbλσ_e²)

θ is invariant to the scale of M and P (β·cov(r̂,M) and cov(r̂,r_p) scale identically), to ridge
shrinkage (both covariances are linear in r̂_p) and, by (I1), it is the only functional of the mirror
construction that is not mechanical. Reparametrised as a reliability,

    λ̂_noise := θ/(1+θ)   equals λ under H_noise, 0 under H_post, and f_Xλ/(1−λ+f_Xλ) ≤ λ in the mixed case.

Reading: λ̂_noise is "the RNA-seq reliability the data would need if everything morphology adds were
the true transcript". If it is far below any plausible single-library reliability for that gene's
depth (section 3), H_noise is excluded quantitatively; for an assumed λ the post-transcriptional
share is 1 − f_X = 1 − (1−λ)θ/λ.

### 2.2 The F18 mirror sign test is θ compared with ρ²/(1−ρ²) (≈ R²_rna)

The cross form S := sign(β)·corr(r̂_p, r_m). By (I1),

    cov(r̂_p, r_m) = (1−ρ²)·cov(r̂_p, M) − δ·cov(r̂_p, r_p)
                   ∝ (1−ρ²)·θ/β − δ        (dividing by cov(r̂_p, r_p) > 0 when morphology has signal)

so with βδ = ρ²:   **S < 0 ⇔ θ < ρ²/(1−ρ²)** (≈ R²_rna for small R²_rna; 0.43 at R²_rna = 0.3), i.e.

    S < 0  ⇔  f_X < (1−λ)ρ² / (λ(1−ρ²)).

Numerically (λ = 0.8, ρ² = 0.2): the sign flips to "post" only when the X-channel is below 6% of the
morphology–residual covariance; at λ = 0.6, ρ² = 0.3: 29%. (Simulation check, section 7.2: on selected
genes the sign of S agrees with the sign of θ_set − ρ̂²/(1−ρ̂²) in 152 of 153 grid cells, the exception
having |z_S| < 2; with the threshold ρ̂² alone the agreement is 138/153.) Hence:

- S significantly negative is strong evidence for a non-transcript channel (θ < ρ²/(1−ρ²) ≪ any plausible λ/(1−λ)).
- S positive is weak evidence either way: it only says f_X exceeds a small threshold, and is
  expected under any mixed model with a modest X-channel. The sign test therefore cannot "confirm
  H_noise"; it can only fail to exclude it. θ itself should be reported.
- H_err gives S < 0 as well (θ = −1), so a negative sign test does not separate "post-transcriptional"
  from "morphology reads an RNA-quality artefact". θ does: 0 vs −1, and the sign of N := β·cov(r̂_p, M)
  (0 vs negative). The transcriptome-wide (20 RNA PCs) baseline removes shared RNA artefacts and should
  be the robustness variant.

Verified at n = 20,000 (`sim_population_check.py`): θ matches f_Xλ/(1−λ) − f_E to ±0.02; S changes sign
where θ crosses ρ²/(1−ρ²) = 0.25 (mix f_X = 0.05, λ = 0.8: θ = 0.198, S = −0.02; f_X = 0.1: θ = 0.41,
S = +0.06). Under selection the sign test inherits the winner's curse of θ (section 7.2): it compares
the *deflated* θ_set with ρ̂²/(1−ρ̂²), so it is biased toward "post" on selected sets.

### 2.3 Where the mechanical correlation enters, and two wrong versions

1. Prediction–prediction form corr(r̂_p, r̂_m). Ridge is linear in its target, so with the same folds,
   features and α: r̂_m = (1−ρ²)·M̂ − δ·r̂_p exactly (M̂ := W-only OOF prediction of M). Then
   cov(r̂_p, r̂_m) = (1−ρ²)cov(r̂_p, M̂) − δ·var(r̂_p). The term −δ·var(r̂_p) is strictly negative even when W
   is pure noise (var r̂_p > 0 always), which is the "biased toward post-transcriptional" effect noted
   in F18. At n = 20,000 it reaches −0.96 under H_post and +0.97 under H_noise; at n = 100 with 20 PCs
   its **null mean is −0.46 (R²_rna = 0.3) to −0.31 (R²_rna = 0.1)** (simulation, selected genes): the
   OOF ridge prediction of a 20-PC block at n_train = 80 has a large variance even without signal, so
   −δ·var(r̂_p) dominates. Read without its null this form says "post-transcriptional" for pure noise.
   The cross form uses r_m, not r̂_m, and its null is centred at ≈ +0.01 (OOF mean effects), so prefer
   the cross form, always with the joint null.
2. Nested-increment form corr(d_p, r_m) with d_p = p̂(P | M, W) − p̂(P | M) (the quantity the
   published pipeline and `measurement_error_checks.py` naturally produce). In-sample the full-minus-
   baseline fitted value is the projection on W⊥M, so cov(d_p, M) = 0 identically and, by (I1),
   cov(d_p, r_m) = −δ·cov(d_p, r_p) < 0 under **every** hypothesis with signal. Simulation: −0.10 under
   H_noise, H_post, H_err and all mixes alike. This version is uninformative and would "confirm
   post-transcriptional" for free. **The sign/θ statistics must be computed from the W-only prediction
   of r_p**, which needs one extra OOF ridge per gene per permutation (hat matrix shared across genes;
   computable locally from `inputs_<c>.npz` + `perms_<c>.npy` of the HPC export — see section 4.5).
   corr(d_p, d_m) (both nested increments) does have the right population sign pattern, but inherits the
   mechanical −δ·var term and the same ρ² threshold; it is a valid but inferior variant.

### 2.4 The joint permutation null

H0: the morphology block carries no information on this patient's (M, P). The null must preserve:

- the (M_i, P_i) pair of every patient for every gene — hence r_p, r_m, β̂, δ̂, ρ̂, the −ρ̂ correlation
  (I2), and the gene–gene dependence of the residuals;
- the W block as a whole (all 20 PCs permuted together, same permutation for every gene and for every
  target r_p, r_m, M), its covariance and its PCA basis;
- the folds, the standardisation and the ridge α (the capacity penalty is part of the null);
- the slide↔patient mapping structure (patients with several slides keep their pooled embedding;
  within-batch permutation for the batch-stratified variant, as in the published controls).

Permuting rows of W is exactly the published patient↔slide permutation; the only change is that the
same draw must be applied to both targets (and to the N statistic) so that aggregate statistics have a
correct joint null. Under this null, N = β·cov(r̂_p, M) and D = cov(r̂_p, r_p) are both centred near 0,
θ = N/D is undefined under the null and must not be tested directly: test N (two-sided) and D (upper)
separately, and report θ with a confidence interval, not a p-value.

One subtlety found in simulation: the permutation null is a null of *no signal at all*, whereas the
hypothesis we want to test with N is *no X-alignment given that morphology has a (u-) signal*. The
covariance form β·cov(r̂_p, M) has sampling variance ∝ var(r̂_p)·var(M)/n, and var(r̂_p) is larger under
H_post than under the permutation null, so the per-gene test of the covariance form over-rejects under
H_post (7–15% at the nominal 5% on selected genes, section 7.2). The **correlation form N_c :=
β̂·corr(r̂_p, M)** is scale-free: under H_post, W ⊥ M, so corr(r̂_p, M) has the same distribution as
under permutation of W, whatever var(r̂_p) is (4–5% rejection in the same cells). Test N_c (per gene
and Σ_g N_c,g at set level); keep the covariance ratio for θ.

A second subtlety concerns *which* permutation. "W ⊥ M under H_post" holds for the biological
signal, not for acquisition structure: if the morphology block encodes the scanning operator (it
does: the twelve largest operator dummies remove 25–40% of embedding variance in the paper) and the
RNA libraries carry an operator-aligned batch effect, then the W-only prediction r̂_p reproduces
the operator means of r_p, which contain −ρ̂·e_op, and cov(r̂_p, M) is negative without any gene-level
alignment. Under the unrestricted null this is read as H_err; under a within-operator null (patients
permuted only within operator, singletons fixed) the operator means are preserved and the null is
centred correctly. Section 7.4 quantifies this and section 4.3 makes the within-operator null primary.

## 3. Depth / expression-level stratification

Single-library reliability rises with counts (Poisson floor var[log(count)] ≈ 1/count plus library
overdispersion), so λ increases with mean expression; ρ_W, ρ_U need not. Predictions for strata of
increasing depth (A = abundance quantile):

| quantity | H_noise | H_post | H_err |
|---|---|---|---|
| Δ_p | ∝ R²_rna(1−λ)²/λ… — falls to ~0 at high depth; no signal possible where λ ≈ 1 | flat in depth (if ρ_U σ_u² is) | falls with depth |
| Δ_m/Δ_p | rises steeply with depth ((λ−ρ²)²/((1−λ)²ρ²)) | ≈ ρ², flat | (1−λ)/(ρ²(1−λ)) = 1/ρ², flat |
| θ (λ̂_noise) | λ/(1−λ), rises with depth; λ̂_noise tracks the depth-implied λ | 0, flat | −1, flat |
| S (sign test) | positive, magnitude shrinks with depth | negative, flat | negative |
| fraction of significant genes | concentrated at low depth / low counts | no depth preference beyond protein-side power | low depth |

The sharpest contrast is the top-depth stratum (ribosomal, translocon, spliceosome transcripts in the
pinned sets; these are among the most abundant mRNAs): under H_noise it should contain essentially no
significant genes and θ should be large there; under H_post it contains the full effect and θ ≈ 0.

Confounders that can mimic or mask the H_noise profile, and mitigations:

1. **Protein completeness / missingness.** Low-abundance proteins have more TMT missing values ⇒ fewer
   fit patients ⇒ larger capacity penalty and lower power; missingness is also not at random
   (truncation at the low end attenuates b). Both push significance toward abundant genes, i.e. the
   opposite of H_noise, and so cannot explain an H_noise-like profile but can hide one. Mitigation:
   stratify on (mRNA depth) × (protein completeness) jointly, or restrict to genes with full
   quantification; use the existing gene_covariate_logistic with mRNA depth, SD and fraction TPM<1 added.
2. **Capacity penalty.** The OOF offset ≈ −(K/n)(1−R²_full) is smaller when R²_rna is larger, and R²_rna
   rises with abundance; a raw Δ̂_p profile therefore drifts upward with depth for a reason unrelated to
   λ. Mitigation: use the permutation-centred increment (observed minus its own null mean — the
   "noise-adjusted advantage" already in the paper), and θ, which is invariant to the offset.
3. **Selection on significance (winner's curse).** Conditioning on Δ̂_p passing FDR inflates D =
   cov(r̂_p, r_p) more where power is lower (low depth, high missingness), so (i) selected-set increments
   are biased up unevenly across strata and (ii) θ is biased **toward 0, i.e. toward H_post** (simulation:
   θ_set 1.2 vs population 1.5 under H_noise at Δ = 0.2). Mitigations: compute θ and S on held-out fold
   seeds (the seed_reliability s1–s9 construction), or on the full tested set with permutation-centred
   statistics, and read θ against its simulation-calibrated value under H_noise *with the same
   selection* (section 7 tabulates it), not against the population λ/(1−λ).
4. **Ratio compression in TMT.** Co-isolation compresses log-ratios multiplicatively; θ and the sign
   are invariant to a common scale of P, but compression is worse for low-abundance peptides and adds a
   plex-shared interference component (a shared ε, read by morphology if plex composition correlates with
   morphology). A shared protein-side artefact enters r_p with +1 and r_m with −δ, i.e. it looks exactly
   like u: the sign test and θ cannot separate biological u from a technical protein artefact.
   Mitigation: this is the domain of the batch/plex-stratified permutation and plex residualisation
   already in the paper; run θ under the batch-stratified null as the primary variant.
5. **Depth proxy.** TPM confounds counts with gene length and library size. Use the raw STAR counts
   (median count across patients) as the depth axis, and the Poisson-implied floor
   λ_Poisson = 1 − mean_i(1/count_i)/var_i(log count) as a per-gene *upper* bound on λ (it ignores
   library overdispersion). Under H_noise, λ̂_noise ≤ λ ≤ λ_Poisson must hold gene by gene; λ̂_noise >
   λ_Poisson is impossible under H_noise and λ̂_noise ≈ 0 ≪ λ_Poisson ≈ 1 is the H_post signature for
   abundant genes.
6. **R²_rna rising with depth.** Under H_noise Δ_p ∝ R²_rna(1−λ)², with the two factors moving in
   opposite directions along the depth axis, so a flat Δ̂_p profile does not by itself exclude H_noise;
   the depth profile of Δ̂_p is a weak test. θ, which does not involve R²_rna, is the right stratified
   quantity.
7. **RNA replicates in CCRCC/UCEC.** Where several RNA files exist per case (CCRCC 39 cases, UCEC 67
   cases), the one-way ICC of log2(TPM+1) across files estimates the reliability of one RNA file as a
   measure of the case-level transcript. Public GDC metadata (queried 2026-10-08, CPTAC-3 primary-tumour
   STAR-counts files) show that multi-file cases are almost always *distinct tumour samples* of the
   same case (e.g. 323 of 323 multi-file "adenomas and adenocarcinomas" cases, 41/41 gliomas; only 12
   ductal/lobular cases have two aliquots of one sample; none is a re-run of one aliquot), so the ICC
   counts intra-tumour sampling variation *plus* library noise. Under H_noise the transcript that a
   whole-slide morphology embedding can read is the case-level one, so the relevant reliability is
   var(X_case)/var(M) = this ICC when RNA and protein come from different portions, and higher (library
   noise only) when they are co-extracted from one portion. Either way **ICC ≤ λ: the replicate ICC
   is a lower bound on λ**, the direction needed (a lower bound on λ gives an upper bound on the
   X-share f_X = (1−λ)θ/λ). For averaged cases (k files) the pipeline's reliability is kICC/(1+(k−1)ICC).
   Per-gene ICCs from 39–67 pairs have SE ≈ 0.06 (ICC 0.8) to 0.12 (ICC 0.5); use set- or stratum-level
   means. If the set-level ICC is ≥ 0.8 for the pinned genes while the selection-calibrated λ̂_noise is
   ≤ 0.3, H_noise is excluded (section 7.6, rule R3) without any external value.

## 4. Proposal

### 4.1 Statistics (per gene g, same folds / α / standardisation as the published pipeline)

Standardise M_g and P_g over the fit patients (so β̂ = δ̂ = ρ̂). Compute r_p = P − ρ̂M, r_m = M − ρ̂P
(OLS on the fit patients; the OOF baseline residual P − p̂(P|M) of the HPC script is equivalent to O(1/n)).
Let r̂_p := W-only 5-fold OOF ridge prediction of r_p.

    N_g  = ρ̂ · cov(r̂_p, M)                own-transcript alignment numerator (for θ)
    N_c,g = ρ̂ · corr(r̂_p, M)              the same, scale-free: this is the quantity that is TESTED
    D_g  = cov(r̂_p, r_p)                   ≈ permutation-centred increment of the W-only model
    θ_g  = N_g / D_g,  λ̂_noise,g = θ_g/(1+θ_g)
    S_g  = sign(ρ̂) · corr(r̂_p, r_m)       F18 sign test, cross form (equivalent to θ_g vs ρ̂²/(1−ρ̂²))

Secondary: Δ̂_m (mirror increment) and Δ̂_m/Δ̂_p; corr(d_p, d_m) for continuity with F18.

### 4.2 Per-gene vs aggregated

Per gene, N_g/D_g is a ratio of two O(Δ) quantities estimated from ~100 patients and is heavy-tailed;
per-gene calls are not recommended (F18 already says so). Report:

- **Set level** (primary): θ_set = Σ_g N_g / Σ_g D_g over the pinned significant set of each cohort,
  reported as the **permutation-centred** ratio θ_pc = (Σ N_g − E₀Σ N_g)/(Σ D_g − E₀Σ D_g) with E₀
  the mean over the joint-null draws (raw θ_set as a sensitivity), with N_c,set = Σ_g N_c,g tested
  two-sided against the joint permutation null (B = 200 draws of one patient↔slide permutation applied
  to every gene), and a patient-bootstrap CI for θ_set (resample patients, refit; the gene-jackknife SE
  in the simulation is an approximation that ignores shared W). Do **not** compute θ over all tested
  genes: null genes contribute D ≈ 0 and the ratio is unstable (section 7.3, θ_set,all of −847 and −47
  in two mixture replicates). The selection-free estimate of θ for the pinned sets is the one computed
  in the **confirmatory lung cohort** on the discovery-pinned 2,310 / 726 sets (no selection there);
  in the discovery cohorts θ_set is read against its selection-calibrated H_noise value (section 7.6).
- **Family level**: the same by Reactome family (ribosome, translocon, spliceosome, …) and by depth
  tertile (raw counts) — family × depth is the table that distinguishes the hypotheses (section 3).
- **Per gene**: the fraction of pinned genes with N_c,g above / below the 2.5–97.5% band of their own
  null, reported as two proportions with the null-expected 2.5% each, not as gene calls.

### 4.3 Null — decision: within-operator (batch-stratified) permutation is primary

Joint patient↔slide permutation of the W rows, identical draw for all genes and all targets, folds and
ridge fixed. Two schemes are run; the **primary** null for N_c, D, S and for θ_pc is the
**within-operator** permutation (patients permuted only within scanning operator read from the slide
headers, `aperio.User` of `build_acquisition_meta.py`; singleton or unlabelled patients held fixed, as
in the paper's registered lung test); the **unrestricted** patient↔slide permutation, which is the
published null for the selection itself, is reported as a sensitivity. Observed statistics are computed
on the raw PCs in both cases (as the pinned sets were selected); residualising W on operator dummies is
*not* used for θ (section 7.4: it biases θ upward under H_post by +0.07 to +0.14).

Justification (simulation in section 7.4, `sim_batchnull.py`):

1. Validity. The null for N_c relies on W ⊥ M under H_post. When the embedding encodes operator
   (ν_W = 0.3 of PC variance) and RNA libraries carry an operator-aligned batch effect (30–50% of the
   mRNA error variance), the unrestricted null reads a pure H_post cohort as H_err in 7/8 (ν_e = 0.3)
   and 8/8 (ν_e = 0.5) replicates (z_N = −3.9 ± 1.0 and −5.5 ± 0.7); between-operator transcript
   differences (accrual) read it as "mixed" in 8/8 (z_N = +5.0 ± 1.1). The within-operator null is
   centred in all batch configurations (mean z_N = −0.5 to +0.9 under H_post; 53/55 H_post seeds R1,
   the two exceptions at z = 2.0 and 2.3 with θ_pc = 0.03).
2. Direction of the error. Operator aliasing pushes the unrestricted test *away* from R1 in both
   directions, so the published null is not "the null that favours the paper's reading"; the stratified
   null is the valid one, and it is also the one that removes the artefact pathway (shared RNA batch =
   H_err via operator) that θ cannot otherwise separate from biology.
3. Power. Under H_noise and H_err the stratified z is 1–4 units smaller than the unrestricted z
   (e.g. +15.0 vs +16.0; −23.0 vs −26.1 without batch) against |z| ≥ 9 in every decisive cell: no
   decision changes. Under accrual-X the stratified null removes between-operator transcript signal and
   z_N falls from +31 to +22 — the within-operator statement is the conservative and intended one.
4. Precedent. The paper's confirmatory lung test and Fig. controls already use within-operator nulls on
   raw statistics ("only the reference moves"), and the UCEC site confounder lies inside the operator axis.
5. θ itself is null-free except through centring: θ_pc under the within-operator null recovers the
   no-batch value (H_post with RNA batch: raw −0.09 / −0.13 → θ_pc,strat −0.02 / −0.01; accrual-X raw
   +0.11 → +0.005; H_noise with RNA batch raw 0.50 → 0.65 against 0.62 without batch).

Disagreement rule: a result that is R4 (or R2/R3) under the unrestricted null but R1 under the
within-operator null is reported as "operator-mediated alignment" — the signature of the batch
configurations above — not as H_err or H_noise. Draws for this null must be *fresh* (not the
RandomState(0) draws of `perms_<c>.npy` that produced the selection), generated locally from the
operator labels. B = 200 is enough for set-level z-scores (simulated |z| ≥ 9 for the decisive
outcomes); B = 1000 if per-gene bands are wanted. PDAC: the scanner axis (SS1553 vs SS7559) is added as
a second stratum level there (permute within operator × scanner).

### 4.4 Decision / reading rules (numeric thresholds calibrated in section 7.6) and pre-writable sentences

Inputs per set (pinned significant set of a cohort; or a family × depth stratum), all under the
within-operator null of 4.3 with the unrestricted null as sensitivity:

| symbol | definition | minimum set size |
|---|---|---|
| z_D | set-level z of Σ D_g against the joint null | 5 |
| z_N | set-level z of N_c,set = Σ ρ̂_g·corr(r̂_p, M) | 5 |
| θ_pc | permutation-centred θ_set = (ΣN − E₀ΣN)/(ΣD − E₀ΣD); λ̂ = θ_pc/(1+θ_pc) | 20 |
| ρ̄² | mean ρ̂² over the set (its R²_rna) | — |
| λ_lo | set-level RNA-replicate ICC (CCRCC, UCEC; lower bound on λ, section 3 item 7); for LUAD/GBM/PDAC the CCRCC/UCEC ICC of the matched raw-count depth stratum, flagged "transferred"; if neither, ρ̄² (rigorous but weak) | — |
| λ_hi | set-level λ_Poisson (upper bound on λ) | — |
| κ | winner's-curse floor: the smallest ratio θ_set,sel/θ_pop that an X-channel produced in simulation (section 7.3, 7.5): **0.4** for n_sel ≥ 20 (0.3 for 5–19); central value 0.7 (n_sel 20–100) / 0.9 (n_sel ≥ 200 when more than half of the signal genes are selected). The deflation is set by selection intensity, which is unknown in the real data, so the floor is used for bounds and the central value for point estimates | — |
| θ_noise,min(λ) | κ·λ/(1−λ): the smallest selected-set θ that pure H_noise at reliability λ produced | — |

Rules (numbers from section 7; false-classification rates in 7.6):

- **R0 — gate.** z_D ≥ 3, else **R6** (the W-only model has no signal on this set; should not happen
  for pinned genes — stop and check `repro_check_<c>.json`). Simulation: z_D ≥ +7 in every selected
  set of ≥ 5 genes, ≥ +10 for n_sel ≥ 10.
- **R4 — z_N ≤ −3 under both nulls and θ_pc ≤ −0.3** ⇒ H_err (morphology reads an RNA-quality /
  shared mRNA artefact). Simulation: H_err gives z_N ∈ [−36, −9] and θ ∈ [−1.02, −0.46] in all 23
  cells; no other hypothesis reached z_N ≤ −3 under the within-operator null. If z_N ≤ −3 only under
  the unrestricted null ⇒ "operator-mediated alignment" (section 7.4), reported as such, and the
  20-RNA-PC transcriptome-wide baseline decides. "The morphology-predicted residual is anti-aligned
  with the measured transcript (θ = …), the signature of a shared mRNA measurement artefact rather
  than of post-transcriptional regulation; under the transcriptome-wide baseline θ = …."
- **R1 — −3 < z_N < 2** (< 3 when n_sel ≥ 1,000, section 7.5) ⇒ H_post; consistency check
  |θ_pc| ≤ 0.05 (n_sel ≥ 100) or ≤ 0.10 (20 ≤ n_sel < 100), which is the operative criterion for the
  large baseline sets. Simulation: H_post gives z_N ≤ 1.9 in 32/32 grid cells and 14/15 mixture sets
  with n_sel ≥ 10, ≤ 1.4 at n_sel ≤ 650, +1.1 to +2.8 at n_sel ≈ 1,500 (an O(1/n) coupling through
  the in-sample ρ̂), and 53/55 batch seeds under the within-operator null (two at +2.0 and +2.3); θ_pc
  never exceeded 0.033 in any H_post set. The blind spot is a mixed model with θ_pop < 0.1 (f_X < 0.1 at
  λ ≤ 0.9), i.e. one in which ≥ 90% of the morphology–residual covariance is non-transcript anyway —
  section 7.6 gives the rates.
  "The morphology-predicted protein residual is uncorrelated with the measured transcript (own-
  transcript alignment θ = …, 95% CI …; within-operator permutation z = …), whereas regression
  dilution of the RNA-seq measurement would require θ ≥ θ_noise,min(λ_lo) = … at the reliability
  λ_lo = … that these genes' own RNA replicates show. The beyond-mRNA signal is therefore not an
  artefact of own-mRNA measurement error; it enters protein abundance through a channel independent
  of the transcript."
- **z_N ≥ 2** (2 ≤ z_N < 3 reported as "weak alignment"): the X-channel is detected; split by θ_pc.
  - **R2 — θ_pc ≥ θ_noise,min(λ_lo)** ⇒ H_noise not excluded (or dominant). Report λ̂ and the X-share
    range f_X ∈ [(1−λ_hi)θ_pc/(λ_hi), (1−λ_lo)θ_pc/(κλ_lo)] ∩ [0,1]. "A substantial part of the
    beyond-mRNA signal is aligned with the measured transcript (θ = …, implying a reliability of λ̂ = …
    if the signal were entirely transcript dilution), which is within the range these genes' replicate
    libraries allow; we therefore cannot exclude that morphology recovers true transcript abundance
    lost to RNA-seq noise, and we describe the estimand as 'beyond the measured transcript'." (Keep
    Limitations (v), strengthen it.)
  - **R3 — 0 < θ_pc < θ_noise,min(λ_lo)** ⇒ mixed: pure H_noise is excluded for every λ ≥ λ_lo and
    the post-transcriptional share is at least 1 − (1−λ_lo)θ_pc/(κλ_lo) (conservative, κ = 0.4; central
    with κ = 0.7). Worked numbers: θ_pc = 0.2 with λ_lo = 0.8 ⇒ f_X ≤ 0.125 (central 0.07) ⇒ post
    share ≥ 0.87; with λ_lo = 0.6 ⇒ f_X ≤ 0.33 ⇒ post share ≥ 0.67; with only λ_lo = ρ̄² = 0.3 the bound
    is vacuous (f_X ≤ 1.2) and the ICC is indispensable. "The alignment θ = … implies that at most
    f_X = … of the morphology–residual covariance can be transcript dilution for any reliability
    ≥ …, the lower bound set by these genes' RNA replicates; the remaining ≥ … is independent of the
    transcript."
- **R5 — heterogeneity.** Apply the rules per Reactome family × raw-count depth tertile with the per-
  stratum n_sel; if the strata fall into different rules, or θ_pc differs between strata by more than
  twice its patient-bootstrap SE, report by stratum and do not pool. "Own-transcript alignment is
  family-specific: … ." Under H_noise θ_pc must *rise* with depth (λ rises) and the top-depth stratum
  must contain few selected genes; under H_post neither holds (section 3).

Sign test alone (if only S is computed): S_set with z_S ≤ −3 ⇒ θ_set < ρ̄²/(1−ρ̄²) ⇒ f_X ≤
(1−λ_lo)ρ̄²/(κλ_lo(1−ρ̄²)) (e.g. ρ̄² = 0.3, λ_lo = 0.8: f_X ≤ 0.27, central 0.15); z_S > 0 ⇒ only "the
sign test does not exclude transcript dilution" (never "confirms"). Simulation: S is negative
(z_S ≤ −3) in every H_post set and positive (z_S ≥ +3) in every H_noise set, but straddles zero across
the mixes, so S never replaces θ.

### 4.5 Compute

Per gene and permutation: one W-only OOF ridge of r_p (hat matrix H_k = W_te(W_tr^⊤W_tr + I)^{-1}W_tr^⊤
shared by all genes for a given fold and permutation ⇒ one (n_te × n_tr) @ (n_tr × G) product per fold),
plus the already-planned nested fits. This adds negligible time to `measurement_error_checks.py`
(two extra permutation sums per gene: Σ r̂_p M and Σ r̂_p r_p; the script currently stores the nested
increment d_p, from which θ *cannot* be recovered — section 2.3). The recommended route is local:
from `inputs_<c>.npz` + `oof_fold_<c>.npy` (+ the operator label per patient from
`build_acquisition_meta.py`, joined on the slide → case map), ≈ 2,500 genes × 2 null schemes × 200
fresh permutations × 5 folds of a 20×20 solve and one matrix product — seconds to minutes (the
synthetic run at G = 2,000, n = 100, 400 permutations takes 40 s single-threaded). `perms_<c>.npy`
(the selection draws) is then not reused for the discriminator null.

## 5. Simulation plan (synthetic only)

`sim_signtest.py`: one cohort per scenario (n patients, K = 20 PCs shared across G genes), hypotheses
H_noise / H_post / H_err / H_mix(f_X ∈ {0.1, 0.3, 0.5}), λ ∈ {0.5, …, 0.9}, population increment
Δ_p ∈ {0.10, 0.20}, R²_rna ∈ {0.1, 0.3}, σ_u² = half of the non-mRNA protein variance; published
estimator (5-fold, ridge α = 1, in-fold standardisation); selection by a separate B = 200 permutation
p-value with BH FDR < 0.05 and Δ̂_p > 0; then N, D, θ, S_cross, S_pp, corr(d_p, d_m), corr(d_p, r_m),
Δ̂_m against a joint B = 200 null. Outputs: per-gene power of each statistic, set-level z on the
selected set, θ_set under selection (the calibration needed for reading rule R2/R3). Runs: n = 100 with
G = 200 and G = 40 (small pinned sets), n = 60 (cohorts with heavy missingness). Infeasible cells
(the hypothesis cannot produce Δ_p at that λ even with ρ = 0.95) are kept and flagged, because they
are themselves the dilution-bound statement.

`sim_mixture.py`: the same model with G = 2,000 (or 5,000) genes of which only a fraction π₁ (0.05–0.5,
main 0.2) carries the signal and the rest are pure null, four (size study: three) replicate seeds per
scenario, so that selection runs against a real null background (false discoveries, FDR) and the
*distribution* of the set-level statistics across replicates — hence false-classification rates of the
reading rules — is estimated. Because BH on a permutation p with B = 200 cannot resolve below 0.005 at
G = 2,000 (it passes nothing or everything; the paper notes the same granularity at B = 1,000), the
selection p for BH is the normal-tail p of the permutation z-score; the permutation p is also stored.

`sim_batchnull.py`: the same model with 12 operators of uneven size plus 10 singleton patients, the
operator offset carrying ν_W = 0.3 of the PC variance, and optional operator-aligned components in the
mRNA error (ν_e = 0.3 / 0.5 of var e), in the protein noise (ν_P = 0.3) or in the transcript itself
(ν_X = 0.3, "accrual"); selection under the unrestricted null as published, then the discriminators
under (a) the unrestricted and (b) the within-operator null, eight replicates per cell.

## 6. Limits of the design

- θ separates "enters P through X" from "enters P not through X"; it does not separate biological u
  from a shared technical protein artefact (section 3, item 4) — that remains the batch-stratified
  null's job.
- A gene-specific mRNA artefact that morphology can read (not shared across genes) behaves like H_err
  (θ = −1) and is detectable; a gene-specific artefact that morphology cannot read is, by definition,
  not what the estimator picks up.
- If a and c are strongly anti-aligned in W-space the X-share f_X can be negative and θ < 0 with
  θ > −1; such a value is reported as "between H_post and H_err" and resolved by the 20-RNA-PC variant.
- Nonlinear f(mRNA): (I1) holds for the OLS/linear residual. With a spline baseline the residuals
  are no longer linear combinations and the mechanical terms change; the derivation applies to the
  published linear estimator, which is the one being defended.
- λ̂_noise is a lower bound on the reliability needed under H_noise only when the winner's-curse bias
  (toward 0) has been accounted for; section 7 gives the selected-set values to compare against.
- Simulation limits. The grids have one cohort per cell (cell-level sampling noise ≈ ±2 in z, ±0.1 in
  θ_set for n_sel ≈ 20); replicate rates come from the mixture and batch studies (4–8 seeds per cell).
  The batch model is stylised (one operator offset shared by all PCs, gene-specific operator effects
  in e and ε drawn independently); it demonstrates the mechanism and the size of the error, not the
  real magnitude, which depends on how much of the operator-aligned embedding variance is matched by
  operator-aligned RNA variance in each cohort. The selection in the mixture study uses a normal-tail
  approximation of the permutation p (section 5). A gene-specific nonlinear mRNA→protein map and
  TMT missingness are not simulated.

## 7. Simulation results (synthetic data only)

### 7.1 Design and files

Grids (`sim_signtest.py`, re-run 2026-10-08 as `sim_grid_v2_*` so that the correlation-form columns
are present; the earlier `sim_grid_*.csv` lack them but agree on every shared column): one cohort per
cell, hypotheses H_noise / H_post / H_err / mix f_X ∈ {0.1, 0.3, 0.5} × λ ∈ {0.5, …, 0.9} ×
population increment Δ_p ∈ {0.10, 0.20} × R²_rna ∈ {0.1, 0.3}; n = 100, K = 20, G = 200 genes all
carrying the signal (also n = 60; and G = 40), σ_u² = half the non-mRNA protein variance; selection by
B = 200 permutation p, BH FDR < 0.05 and Δ̂_p > 0; discriminators against a separate B = 200 joint null.
Cells where the hypothesis cannot reach Δ_p even with ρ = 0.95 are run at the cap and marked †.
Mixtures (`sim_mixture.py`): G = 2,000, 20% signal genes, 4 seeds (set-size study: G = 2,000 / 5,000,
5–50% signal, 3 seeds). Batch (`sim_batchnull.py`): 12 operators + 10 singletons, 8 seeds. Full tables:
`section7_tables.md`. Cell-level sampling noise in the grids is about ±2 in z and ±0.1 in θ_set at
n_sel ≈ 20, which the replicate studies confirm.

### 7.2 Grids: feasibility, selection, set-level z, θ under selection, S

**Table 7.1 — n = 100, G = 200, R²_rna = 0.3: selection rate n_sel/G · set-level z of N_c,set**
(– = no gene selected; † = infeasible, run at ρ = 0.95)

| hypothesis | Δ0.1 λ0.5 | Δ0.1 λ0.6 | Δ0.1 λ0.7 | Δ0.1 λ0.8 | Δ0.1 λ0.9 | Δ0.2 λ0.5 | Δ0.2 λ0.6 | Δ0.2 λ0.7 | Δ0.2 λ0.8 | Δ0.2 λ0.9 |
|---|---|---|---|---|---|---|---|---|---|---|
| H_noise | 0.17 · +16 | 0.12 · +16 | 0.10 · +16 | –† | –† | 0.65 · +33 | 0.53 · +34† | 0.21 · +21† | –† | –† |
| mix f_X=0.5 | 0.14 · +10 | 0.18 · +16 | 0.47 · +25† | 0.19 · +14† | –† | 0.89 · +29 | 0.91 · +32† | 0.72 · +25† | 0.46 · +16† | 0.30 · +8† |
| mix f_X=0.3 | 0.17 · +5 | 0.13 · +7 | 0.33 · +20 | 0.27 · +17† | –† | 0.68 · +20 | 0.85 · +27 | 0.80 · +26† | 0.62 · +18† | 0.43 · +8† |
| mix f_X=0.1 | 0.12 · +3 | – | 0.15 · +7 | 0.10 · +7 | 0.17 · +6† | 0.73 · +7 | 0.66 · +11 | 0.64 · +17 | 0.76 · +19† | 0.64 · +10† |
| H_post | – | 0.20 · −1 | 0.10 · +1 | 0.10 · −1 | 0.17 · +1 | 0.62 · −0† | 0.67 · +1 | 0.58 · +1 | 0.60 · +0 | 0.59 · +2 |
| H_err | 0.33 · −23 | 0.39 · −26 | 0.26 · −20† | –† | –† | 0.95 · −36† | 0.55 · −31† | 0.17 · −17† | –† | –† |

Feasibility is the dilution bound as a selection statement: at R²_rna = 0.1 neither H_noise nor H_err
can produce Δ_p = 0.1 at any λ ≥ 0.5 (every cell †, no gene selected), while H_post selects 34–44% of
genes at Δ_p = 0.2 at every λ with z_N ∈ [−0.4, +1.1]; mixes with f_X = 0.1 give z_N = +4.5 to +13 and
θ_set/θ_pop = 0.09/0.10, 0.12/0.15, 0.14/0.15, 0.07/0.10, 0.04/0.05 (λ = 0.5 … 0.9). A pinned set
dominated by genes with R²_rna ≈ 0.1 cannot be H_noise at all. At R²_rna = 0.3 H_noise is feasible
only for λ ≤ 0.7 (Δ_p = 0.1) or λ = 0.5 (Δ_p = 0.2), and its selection rate falls with λ (0.17 → 0.10)
where H_post's is flat — the second finite-sample effect of section 1 (collinearity of W with M).

**Table 7.2 — same grid: θ_set under selection / population θ (winner's curse)**

| hypothesis | Δ0.1 λ0.5 | Δ0.1 λ0.6 | Δ0.1 λ0.7 | Δ0.1 λ0.8 | Δ0.1 λ0.9 | Δ0.2 λ0.5 | Δ0.2 λ0.6 | Δ0.2 λ0.7 | Δ0.2 λ0.8 | Δ0.2 λ0.9 |
|---|---|---|---|---|---|---|---|---|---|---|
| H_noise | +0.56/+1.00 | +0.74/+1.50 | +1.09/+2.33 | –† | –† | +0.96/+1.00 | +1.18/+1.50† | +1.15/+2.33† | –† | –† |
| mix f_X=0.5 | +0.29/+0.50 | +0.54/+0.75 | +0.78/+0.99† | +0.51/+0.89† | –† | +0.48/+0.50 | +0.67/+0.65† | +0.57/+0.63† | +0.36/+0.50† | +0.19/+0.28† |
| mix f_X=0.3 | +0.14/+0.30 | +0.26/+0.45 | +0.56/+0.70 | +0.58/+0.68† | –† | +0.28/+0.30 | +0.45/+0.45 | +0.50/+0.48† | +0.33/+0.37† | +0.15/+0.20† |
| mix f_X=0.1 | +0.10/+0.10 | – | +0.22/+0.23 | +0.27/+0.40 | +0.19/+0.31† | +0.08/+0.10 | +0.13/+0.15 | +0.23/+0.23 | +0.29/+0.29† | +0.14/+0.16† |
| H_post | – | −0.03/0 | +0.02/0 | −0.03/0 | +0.01/0 | −0.01/0† | +0.01/0 | +0.01/0 | +0.00/0 | +0.02/0 |
| H_err | −0.76/−1 | −0.74/−1 | −0.64/−1† | –† | –† | −1.02/−1† | −0.77/−1† | −0.67/−1† | –† | –† |

Selection deflates θ (D = cov(r̂_p, r_p) is inflated by selecting on Δ̂_p; N is not selected on): the
ratio θ_set,sel/θ_pop over the three grids (cells with n_sel ≥ 5, X-channel hypotheses) is 0.45–0.63
(median 0.49) when ≤ 25% of genes are selected and 0.67–1.10 (median 0.93) when > 50% are; by set
size, min/median = 0.42/0.72 (n_sel 10–20), 0.48/0.70 (20–50), 0.41/0.79 (50–100), 0.79/0.93
(100–200). H_post stays at θ_set ∈ [−0.03, +0.05] and H_err at [−1.02, −0.46] regardless. The
oracle θ over *unselected* signal genes is, if anything, above the population value at n = 100
(1.09–1.15 at θ_pop = 1.0; 1.65–1.95 at 1.5; section 7.3), so the deflation is entirely a
selection effect and the calibration must be made against the selected-set values, as in 4.4.

**Table 7.3 — same grid: S_cross mean on selected genes (set-level z)**

| hypothesis | Δ0.1 λ0.5 | Δ0.1 λ0.6 | Δ0.1 λ0.7 | Δ0.1 λ0.8 | Δ0.1 λ0.9 | Δ0.2 λ0.5 | Δ0.2 λ0.6 | Δ0.2 λ0.7 | Δ0.2 λ0.8 | Δ0.2 λ0.9 |
|---|---|---|---|---|---|---|---|---|---|---|
| H_noise | +0.08 (+3) | +0.13 (+3) | +0.19 (+7) | –† | –† | +0.19 (+13) | +0.24 (+17)† | +0.19 (+8)† | –† | –† |
| mix f_X=0.5 | −0.03 (−1) | +0.06 (+2) | +0.11 (+6)† | +0.04 (+1)† | –† | +0.04 (+1) | +0.11 (+9)† | +0.07 (+5)† | −0.02 (−3)† | −0.09 (−6)† |
| mix f_X=0.3 | −0.13 (−6) | −0.06 (−2) | +0.05 (+3) | +0.06 (+3)† | –† | −0.08 (−7) | +0.02 (+2) | +0.04 (+2)† | −0.04 (−5)† | −0.12 (−9)† |
| mix f_X=0.1 | −0.16 (−6) | – | −0.09 (−5) | −0.06 (−2) | −0.10 (−5)† | −0.19 (−15) | −0.16 (−13) | −0.11 (−9) | −0.05 (−6)† | −0.13 (−11)† |
| H_post | – | −0.22 (−11) | −0.21 (−7) | −0.24 (−9) | −0.20 (−9) | −0.24 (−18)† | −0.24 (−19) | −0.24 (−17) | −0.24 (−18) | −0.23 (−17) |
| H_err | −0.44 (−23) | −0.47 (−26) | −0.44 (−20)† | –† | –† | −0.55 (−36)† | −0.49 (−32)† | −0.43 (−17)† | –† | –† |

S is one-sided exactly as derived: negative under H_post (z_S ≤ −3.5 in every cell) and H_err, positive
under H_noise (z_S ≥ +3), and crossing zero across the mixes where θ_set crosses ρ̂²/(1−ρ̂²) (152/153
cells over the three grids). Mechanical offsets: the cross-form null mean is +0.011 (sd 0.003 over
cells); the prediction–prediction form S_pp has null mean −0.46 (−0.31 at R²_rna = 0.1); the
nested-increment form corr(d_p, r_m) is −0.20 ± 0.05 under every hypothesis alike (−0.25 at n = 60).

**Table 7.4 — summary over the three grids (cells with n_sel ≥ 10)**

| grid | hypothesis | cells | z_N min / med / max | z_S min / med / max | θ_set min / med / max | θ_set/θ_pop med |
|---|---|---|---|---|---|---|
| n100 G200 | H_noise | 6 | +15.9 / +18.8 / +33.6 | +3.2 / +7.4 / +17.3 | +0.56 / +1.03 / +1.18 | 0.53 |
| | mix 0.5 / 0.3 / 0.1 | 15 / 15 / 14 | +3.2 / +15.6 / +31.5 · +2.4 / +16.3 / +26.8 · +3.1 / +8.3 / +19.4 | −6.3 / +1.6 / +8.7 · −8.6 / +1.7 / +7.5 · −15.4 / −4.7 / +2.8 | +0.08 / +0.29 / +0.78 · +0.04 / +0.23 / +0.58 · +0.04 / +0.13 / +0.29 | 0.72 · 0.80 · 0.87 |
| | H_post | 15 | −1.4 / +0.5 / +1.9 | −19.3 / −10.1 / −3.5 | −0.03 / +0.01 / +0.02 | — |
| | H_err | 6 | −35.9 / −24.3 / −17.3 | −36.2 / −24.3 / −17.4 | −1.02 / −0.75 / −0.64 | 0.75 |
| n60 G200 | H_noise | 2 | +9.7 / +12.2 / +14.6 | +3.1 / +3.8 / +4.6 | +0.63 / +0.65 / +0.67 | 0.54 |
| | mix 0.5 / 0.3 / 0.1 | 4 / 7 / 9 | min +5.6 · +5.1 · +2.1 | — | med +0.36 · +0.29 · +0.10 | 0.66 · 0.71 · 0.68 |
| | H_post | 7 | −1.5 / +0.8 / +1.5 | −10.4 / −8.7 / −3.1 | −0.02 / +0.01 / +0.04 | — |
| | H_err | 3 | −19.5 / −15.3 / −10.6 | −18.6 / −15.5 / −13.3 | −0.87 / −0.64 / −0.46 | 0.64 |
| n100 G40 | H_noise | 3 | +11.5 / +18.0 / +19.4 | +5.2 / +6.0 / +6.7 | +0.99 / +1.17 / +1.18 | 0.79 |
| | mix 0.5 / 0.3 / 0.1 | 10 / 11 / 10 | min +5.8 · +2.8 · +1.8 | — | med +0.46 · +0.21 · +0.11 | 0.83 · 0.90 · 0.79 |
| | H_post | 10 | −0.7 / +0.6 / +1.8 | −9.8 / −6.5 / −3.6 | −0.02 / +0.01 / +0.05 | — |
| | H_err | 6 | −22.4 / −14.5 / −8.8 | −20.9 / −13.3 / −10.5 | −1.02 / −0.74 / −0.59 | 0.74 |

**Per-gene over-rejection of the covariance form under H_post** (two-sided, nominal 5%, selected
genes): covariance form N = ρ̂·cov(r̂_p, M) rejects 0.105 (range 0.05–0.15; n = 100, G = 200, 15 cells),
0.069 (n = 60) and 0.069 (G = 40); the correlation form N_c rejects 0.049 (0.00–0.14), 0.077 and 0.041
in the same cells. In the mixtures (section 7.3) the gap is larger on strongly selected sets (cov
0.09–0.18 vs corr 0.03–0.07 at n_sel ≈ 230) and both forms are at 0.04–0.06 on pure-null genes.
At set level both forms are calibrated (max |z| under H_post 2.1 vs 1.9), so the choice matters for
per-gene bands only; N_c is used throughout.

### 7.3 Mixture cohorts (null + signal genes, replicate seeds): winner's curse, FDR, classification

**Table 7.5 — `sim_mixture_main.csv`: G = 2,000, 400 signal genes, n = 100, 4 seeds** (mean ± sd over
seeds; "false" = selected null genes; θ oracle = ratio of sums over all 400 signal genes, no selection;
θ_pc = permutation-centred; classes by the rules of 4.4; † = run at the ρ = 0.95 cap)

| hypothesis | R² | λ | Δ_p | n_sel (false) | θ_pop | θ oracle | θ_set,sel | θ_pc | z_N | z_S | class |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H_post | 0.3 | 0.6 | 0.1 | 9 (0) [3–23] | 0 | +0.01 | +0.10±0.16 | +0.01±0.02 | +0.1±0.3 | −8.5 | R6×2, R1×2 |
| H_post | 0.3 | 0.6 | 0.2 | 227 (6) | 0 | 0.00 | 0.00±0.01 | 0.00±0.01 | +0.4±0.9 | −27 | R1×4 |
| H_post | 0.3 | 0.8 | 0.1 | 4 (0) [1–10] | 0 | +0.01 | +0.02±0.02 | — | — | — | R6×3, R1×1 |
| H_post | 0.3 | 0.8 | 0.2 | 231 (5) | 0 | 0.00 | 0.00±0.01 | 0.00±0.01 | +0.2±0.6 | −26 | R1×4 |
| H_post | 0.1 | 0.6 | 0.2 | 108 (4) | 0 | +0.01 | +0.01±0.01 | +0.01±0.01 | +1.3±0.8 (max 2.3) | −14 | R1×3, R3×1 |
| H_noise | 0.3 | 0.5 | 0.1 | 6 (0) | +1.00 | +1.13 | +0.40±0.05 | +0.38 | +7.4±1.0 | −0.5 | R3×3, R6×1 |
| H_noise | 0.3 | 0.6 | 0.1 | 6 (0) | +1.50 | +1.80 | +0.54±0.07 | +0.54 | +9.5±2.4 | +2.3 | R2×3, R6×1 |
| H_noise | 0.3 | 0.7 | 0.1 | 8 (0) | +2.33 | +3.43 | +0.79±0.13 | +0.77 | +11.5±2.6 | +3.8 | R2×3, R6×1 |
| H_noise | 0.3 | 0.5 | 0.2 | 233 (5) | +1.00 | +1.05 | +0.83±0.03 | +0.79 | +36±1 | +20 | R2×4 |
| H_err | 0.3 | 0.5 | 0.1 | 80 (4) | −1 | −1.09 | −0.65±0.03 | −0.63 | −27±2 | −31 | R4×4 |
| H_err | 0.3 | 0.6 | 0.1 | 79 (3) | −1 | −1.10 | −0.64±0.05 | −0.62 | −25±2 | −31 | R4×4 |
| mix f_X=0.1 | 0.3 | 0.6 | 0.2 | 250 (6) | +0.15 | +0.15 | +0.13±0.01 | +0.13 | +14.5 | −22 | R3×4 |
| mix f_X=0.1† | 0.3 | 0.8 | 0.2 | 321 (8) | +0.29 | +0.32 | +0.28±0.01 | +0.27 | +26 | −12 | R3×4 |
| mix f_X=0.3 | 0.3 | 0.6 | 0.2 | 340 (11) | +0.45 | +0.47 | +0.43±0.01 | +0.41 | +33 | −0.7 | R3×4 |
| mix f_X=0.3† | 0.3 | 0.8 | 0.2 | 242 (7) | +0.37 | +0.40 | +0.32±0.01 | +0.31 | +24 | −7.8 | R3×4 |
| mix f_X=0.5† | 0.3 | 0.6 | 0.2 | 373 (10) | +0.65 | +0.69 | +0.64±0.02 | +0.61 | +38 | +14 | R2×4 |
| mix f_X=0.5† | 0.3 | 0.8 | 0.2 | 92 (2) | +0.50 | +0.59 | +0.36±0.05 | +0.35 | +18 | −2.4 | R3×4 |
| mix f_X=0.1 | 0.1 | 0.6 | 0.2 | 144 (4) | +0.15 | +0.17 | +0.12±0.01 | +0.12 | +17 | +0.6 | R3×4 |
| mix f_X=0.3† | 0.1 | 0.6 | 0.2 | 98 (2) | +0.23 | +0.25 | +0.16±0.01 | +0.16 | +18 | +5.8 | R3×4 |

Readings. (i) Selection against an 80% null background realises FDR 0.02–0.04 and, at Δ_p = 0.1,
picks only 1–6% of the signal genes (3–23 of 400), so a pinned gene has a population increment of
roughly 0.2 or more — the capacity-penalty statement of section 1 again. (ii) Winner's curse: with
1–3% of signal genes selected, θ_set,sel/θ_pop is 0.34–0.44 under H_noise (0.28–0.44 pooling the
grids' n_sel 10–20 cells) and θ_pc does not repair it (0.33–0.42); with ~58% selected (Δ_p = 0.2) it is
0.80–0.87 (H_noise) and 0.83–1.02 (mixes); mixes at n_sel 50–100 give 0.61–0.85. The oracle θ over
unselected signal genes is at or above θ_pop (1.05–3.43 vs 1.0–2.33), so the deflation is a pure
selection effect, and its size is governed by selection intensity, not by n_sel as such. (iii) H_post
is calibrated: θ_set,sel = 0.00 ± 0.01 and z_N ∈ [−0.5, +1.8] at n_sel ≈ 230, one seed at z = 2.3
(n_sel = 98, R²_rna = 0.1, θ_pc = 0.016). (iv) The sign test on the selected set reproduces the
threshold ρ̂²/(1−ρ̂²) = 0.43 applied to the *deflated* θ: mix f_X = 0.3 at λ = 0.6 (θ_set 0.43) sits
exactly on it (z_S = −0.7 ± 1.9). (v) Per-gene tests on selected H_post genes: covariance form rejects
0.117, correlation form 0.066 (nominal 0.05); on pure-null genes 0.049 / 0.048; fractions outside the
2.5 / 97.5% null bands on selected H_post genes 0.027 / 0.030 (correlation form) against 0.068 / 0.069
(covariance form) — the per-gene band proportions of 4.2 must use N_c. (vi) θ over *all* tested genes
is unusable in a mixture (−221 ± 418 and +27 ± 165 in two cells; D_all ≈ 0).

### 7.4 Operator-batch structure: unrestricted vs within-operator null

Table 7.6 (`sim_batchnull.csv`; n = 100, G = 200, R²_rna = 0.3, λ = 0.5, Δ_p = 0.10, 12 operators of
uneven size + 10 singletons, 8 seeds per cell; selection under the unrestricted null as published;
values mean ± sd over seeds; classes by the rules of 4.4 with θ_noise,min at λ = 0.5 taken as 0.5).

| config (ν_W, ν_e, ν_P, ν_X) | hyp | n_sel | θ_pop | θ_set raw | θ_pc within-op | z_N unrestricted | z_N within-op | class unres → within-op |
|---|---|---|---|---|---|---|---|---|
| no batch (0,0,0,0) | H_post | 22 | 0 | +0.01±0.03 | +0.01 | +0.4±1.0 | +0.4±1.1 | R1×7 → R1×7 (one set < 5 genes) |
| | H_noise | 32 | +1.00 | +0.62±0.05 | +0.60 | +16.0±2.0 | +15.0±1.4 | R2×8 → R2×8 |
| | H_err | 75 | −1 | −0.75±0.05 | −0.71 | −26.1±2.2 | −23.0±1.5 | R4×8 → R4×8 |
| | mix 0.3 | 32 | +0.30 | +0.20±0.03 | +0.20 | +7.0±1.5 | +6.4±1.4 | R3×8 → R3×8 |
| W only (0.3,0,0,0) | H_post | 40 | 0 | −0.01±0.02 | −0.01 | −0.3±0.9 | −0.3±0.9 | R1×8 → R1×8 |
| | H_noise | 45 | +1.00 | +0.67±0.05 | +0.66 | +19.4±2.1 | +18.6±2.3 | R2×8 → R2×8 |
| RNA batch (0.3,0.3,0,0) | H_post | 51 | 0 | **−0.09±0.02** | −0.02 | **−3.9±1.0** | +0.0±1.1 | **R4×7**, R1×1 → R1×7, R3×1 |
| | H_noise | 46 | +1.00 | +0.50±0.07 | +0.65 | +15.3±1.6 | +18.2±1.8 | R3×5, R2×3 → R3×5, R2×3 |
| | H_err | 67 | −1 | −0.73±0.06 | −0.69 | −23.7±1.0 | −20.0±0.9 | R4×8 → R4×8 |
| | mix 0.3 | 51 | +0.30 | +0.12±0.02 | +0.21 | +4.9±0.9 | +9.5±1.4 | R3×8 → R3×8 |
| RNA batch strong (0.3,0.5,0,0) | H_post | 62 | 0 | **−0.13±0.02** | −0.01 | **−5.5±0.7** | +0.8±0.7 | **R4×8** → R1×8 |
| | H_noise | 49 | +1.00 | +0.41±0.03 | +0.63 | +13.0±1.8 | +18.3±1.2 | R3×8 → R3×8 |
| | mix 0.3 | 58 | +0.30 | +0.08±0.04 | +0.21 | +3.3±1.6 | +9.8±1.4 | R3×6, R1×2 → R3×8 |
| protein batch (0.3,0,0.3,0) | H_post | 58 | 0 | +0.01±0.01 | +0.01 | +0.4±0.7 | +0.5±0.9 | R1×8 → R1×8 |
| | H_noise | 62 | +1.00 | +0.61±0.03 | +0.68 | +21.6±1.4 | +20.7±2.3 | R2×8 → R2×8 |
| both (0.3,0.3,0.3,0) | H_post | 68 | 0 | −0.06±0.02 | +0.01 | −3.0±1.1 | +0.9±1.1 | R4×4, R1×4 → R1×7, R3×1 |
| | H_noise | 59 | +1.00 | +0.45±0.04 | +0.67 | +16.3±1.3 | +21.1±1.5 | R3×7, R2×1 → R3×7, R2×1 |
| accrual X (0.3,0,0,0.3) | H_post | 61 | 0 | **+0.11±0.02** | +0.01 | **+5.0±1.1** | −0.5±1.1 | **R3×8** → R1×8 |
| | H_noise | 119 | +1.00 | +0.86±0.08 | +0.80 | +31.0±1.7 | +22.0±1.6 | R2×8 → R2×8 |
| | H_err | 109 | −1 | −0.58±0.04 | −0.80 | −20.2±1.0 | −26.0±2.1 | R4×8 → R4×8 |
| | mix 0.3 | 98 | +0.30 | +0.38±0.01 | +0.28 | +19.5±1.2 | +12.7±1.3 | R3×8 → R3×8 |

Readings. (i) Without operator-aligned RNA or transcript structure the two nulls agree (z within 1–3
units, identical classes). (ii) An RNA batch effect that morphology can read through operator identity
makes the unrestricted test anti-conservative in the H_err direction: a pure H_post cohort is read as
H_err in 7/8 and 8/8 seeds, and raw θ_set is pulled to −0.09 / −0.13; the within-operator null is
centred (z_N = 0.0 ± 1.1, +0.8 ± 0.7) and θ_pc under it is −0.02 / −0.01. (iii) Between-operator
transcript differences (accrual) push the other way: H_post is read as mixed in 8/8 seeds (z_N = +5.0,
raw θ +0.11) under the unrestricted null and as H_post in 8/8 under the within-operator null
(θ_pc +0.01). (iv) Residualising W on the twelve largest operators does *not* fix θ: under H_post with
RNA batch it gives +0.07 / +0.14 (over-correction with the smaller operators left in), so θ_pc under
the within-operator null is the estimator of record and residualised-W θ is not reported. (v) Power:
under H_noise / H_err the within-operator z is 1–4 units below the unrestricted z against |z| ≥ 13 —
no decision changes; under accrual-X it removes the between-operator transcript signal (+31 → +22),
which is the intended within-operator statement. (vi) The H_noise cells at λ = 0.5 sit on the R2/R3
boundary (θ_set 0.41–0.67 vs θ_noise,min(0.5) = 0.5 with κ = 0.5 used here; 0.2 with κ = 0.4): at a
reliability of 0.5 the mixed/H_noise distinction is weak by construction (θ_pop is only 1), and the
reading rules therefore make θ_noise,min depend on λ_lo rather than on a fixed number.

### 7.5 Set-size dependence

**Table 7.7 — `sim_mixture_size.csv`: R²_rna = 0.3, Δ_p = 0.2, 3 seeds** (SD = sd of θ_set,sel over
seeds; JK = mean gene-jackknife SE; z_N range over seeds)

| hypothesis | G | signal % | n_sel | θ_pop | θ_set,sel (SD) | θ_pc | JK SE | z_N | class |
|---|---|---|---|---|---|---|---|---|---|
| H_post λ=0.6 | 2,000 | 5 | 40 | 0 | +0.012 (0.013) | +0.011 | 0.022 | +0.1 … +1.1 | R1×3 |
| | 2,000 | 50 | 649 | 0 | +0.001 (0.006) | +0.001 | 0.006 | −0.5 … +1.4 | R1×3 |
| | 5,000 | 50 | 1,559 | 0 | +0.005 (0.003) | +0.005 | 0.004 | +1.1 … +2.8 | R1×1, R3×2 |
| H_noise λ=0.5 | 2,000 | 5 | 31 | +1.00 | +0.69 (0.06) | +0.66 | 0.052 | +18.7 … +20.8 | R2×3 |
| | 2,000 | 50 | 630 | +1.00 | +0.87 (0.03) | +0.83 | 0.018 | +43 … +49 | R2×3 |
| | 5,000 | 50 | 1,623 | +1.00 | +0.87 (0.02) | +0.83 | 0.011 | +47 … +47 | R2×3 |
| mix f_X=0.3 λ=0.6 | 2,000 | 5 | 68 | +0.45 | +0.41 (0.05) | +0.39 | 0.028 | +19 … +22 | R3×3 |
| | 2,000 | 50 | 835 | +0.45 | +0.44 (0.02) | +0.42 | 0.008 | +38 … +39 | R3×3 |
| | 5,000 | 50 | 2,120 | +0.45 | +0.44 (0.01) | +0.42 | 0.005 | +43 … +45 | R3×3 |

Readings. (i) Precision: the seed-to-seed SD of θ_set is 0.013 / 0.006 / 0.003 under H_post and
0.06 / 0.03 / 0.02 under H_noise at n_sel ≈ 40 / 650 / 1,600, and the gene-jackknife SE tracks it
(0.022 / 0.006 / 0.004; 0.052 / 0.018 / 0.011), so the jackknife is an acceptable SE and the ±0.05
consistency band of rule R1 is > 3 SD wide for n_sel ≥ 40. (ii) Large sets: under H_post the selected
genes carry a per-gene excess of N_c over its null of +0.001 to +0.004 — an O(1/n) coupling, because
ρ̂ and r_p are fitted in-sample and so contain the test-fold M — which is negligible against the
per-gene null sd (0.05–0.06) but makes the *sum's* z reach +1.1 to +2.8 at n_sel ≈ 1,500 (2 of 3 seeds
≥ 2; θ_pc = 0.002–0.009). For sets of ≥ 1,000 genes (the paper's baseline pinned sets are 702–2,566)
the alignment threshold is therefore z_N ≥ 3 and the operative H_post criterion is |θ_pc| ≤ 0.05;
computing ρ̂ and r_p within training folds (the OOF baseline residual that `measurement_error_checks.py`
already exports) removes the coupling and is the recommended implementation. (iii) Deflation at 50%
selection is 0.86–0.90 (H_noise) and 0.95–1.02 (mix) independent of G; at 5% signal prevalence and
n_sel 30–70 it is 0.63–0.75 and 0.80–1.01.

### 7.6 Calibrated thresholds and false-classification rates

**Thresholds** (as in 4.4; all under the within-operator null; n_sel ≥ 20 for any θ reading):

| quantity | threshold | simulation support |
|---|---|---|
| z_D gate | ≥ 3 | min over all selected sets of ≥ 5 genes: 6.2 (mixture, n_sel 5–9), ≥ 10 for n_sel ≥ 10, ≥ 13 batch study |
| H_err | z_N ≤ −3 and θ_pc ≤ −0.3 | H_err: z_N ∈ [−36, −9], θ ∈ [−1.02, −0.46] in 23 grid/mixture cells and 56 batch seeds; no other hypothesis gave z_N ≤ −3 under the within-operator null (0/395 sets) |
| H_post | −3 < z_N < 2 (n_sel < 1,000) or < 3 (n_sel ≥ 1,000), and |θ_pc| ≤ 0.05 (n_sel ≥ 100) / 0.10 (20–99) | H_post: θ_set ∈ [−0.06, +0.05] in every set with n_sel ≥ 10 (θ_pc ≤ 0.033); z_N ≤ 1.9 in 32/32 grid cells, ≤ 1.8 in 13/14 mixture sets (one 2.3), ≤ 1.4 at n_sel ≤ 650 and 1.1–2.8 at n_sel ≈ 1,500; 53/55 batch seeds < 2 under the within-operator null (max 2.3) |
| alignment | z_N ≥ 2 (≥ 3 for n_sel ≥ 1,000); 2 ≤ z_N < 3 = "weak" | mixes with θ_pop ≥ 0.1: z_N ≥ 2.1 in 98/99 grid cells with n_sel ≥ 10 (one 1.8 at n_sel = 10), ≥ 14 in 32/32 mixture sets, ≥ 3.3 in 64/64 batch seeds (within-operator) |
| θ_noise,min(λ) = κ·λ/(1−λ) | κ = 0.4 (n_sel ≥ 20), 0.3 (5–19); central 0.7 (20–100), 0.9 (≥ 200 with > 50% selection) | θ_set,sel/θ_pop under H_noise: 0.34–0.44 at n_sel 5–9, 0.28–0.63 at 10–50 (grids + mixture), 0.63–0.75 at 27–35 with 5% prevalence, 0.80–0.90 at ≥ 230 with 50–58% selection; mixes 0.41–1.17 (n_sel 10–100), 0.80–1.02 (≥ 200) |
| sign test | z_S ≤ −3 ⇒ θ_set < ρ̄²/(1−ρ̄²) | sign of S agrees with sign(θ_set − ρ̂²/(1−ρ̂²)) in 152/153 grid cells |

**False-classification rates** (rules applied mechanically to every simulated set with n_sel ≥ 10;
batch study under the within-operator null unless stated):

| truth | sets | read as H_post (R1) | read as aligned (R2/R3) | read as H_err (R4) |
|---|---|---|---|---|
| H_post | 32 grid + 14 mixture + 9 size + 55 batch = 110 | 105 | 5 (z = 2.3, 2.0, 2.8, 2.0, 2.3; θ_pc ≤ 0.033; 0 with the joint criterion z ≥ 3 for n_sel ≥ 1,000 **or** θ_pc > 0.05) | 0 |
| H_post, unrestricted null, RNA-batch / both / accrual configurations | 32 batch | 5 | 8 (accrual) | 19 (RNA batch, both) |
| H_err | 15 grid + 8 mixture + 56 batch = 79 | 0 | 0 | 79 |
| H_noise | 11 grid + 13 mixture + 9 size + 56 batch = 89 | 0 | 89 (R3 rather than R2 in 3 mixture seeds at λ = 0.5, n_sel = 6, and in 20 batch seeds at λ = 0.5 with the fixed θ threshold 0.5 used in the scripts; the λ-dependent θ_noise,min(0.5) = 0.4 of 4.4 reads all but one of them R2) | 0 |
| mixes, θ_pop ≥ 0.1 | 99 grid + 32 mixture + 9 size + 56 batch = 196 | 1 (grid, f_X = 0.1, λ = 0.5, R² = 0.1, n_sel = 10, z = 1.8) | 195 | 0 |
| mix f_X = 0.3, unrestricted null, strong RNA batch | 8 batch | 2 | 6 | 0 |

So, under the within-operator null, the H_post reading has a false-positive rate of ≈ 5% (5/110) at
the z ≥ 2 threshold and 0/110 with the joint criterion, mixes with an X-share of ≥ 10% of the
covariance are missed in ≈ 0.5% of sets (1/196), and H_err is never confused with anything (0/79
either way). The blind spot is a mix with
θ_pop < 0.1 (f_X < 0.1 at λ ≤ 0.9), which is read as H_post; this is a model in which ≥ 90% of the
morphology–residual covariance is non-transcript, so the pre-written sentence R1 remains true in
substance ("not an artefact of own-mRNA measurement error") and the quantitative statement that goes
with it is "f_X ≤ 0.1". Sets with fewer than 5 selected genes (9 of 28 mixture replicates at Δ_p =
0.1) return R6 by construction and are reported as "no set".

### 7.7 Pre-registration checklist for the real run

1. Compute per cohort, pinned set, and family × depth strata: n_sel, z_D, z_N (within-operator; and
   unrestricted), θ_pc (within-operator) with patient-bootstrap CI, λ̂, ρ̄², λ_lo (ICC; transferred by
   depth stratum for LUAD/GBM/PDAC), λ_hi (λ_Poisson), S_set z; the C1-LUAD confirmatory cohort on the
   2,310 and 726 sets as the selection-free θ.
2. Apply R0 → R4 → R1 → R2/R3 → R5 of 4.4 in that order; report the disagreement rule of 4.3 when
   the two nulls differ.
3. Write the pre-written sentence of the rule that fires, with the numbers; never per-gene calls.

## 8. Change log (2026-10-08, section 7 author)

Files added in this folder: `sim_mixture.py`, `sim_batchnull.py`, `sim_section7_tables.py`,
`sim_grid_v2_{n100_G200,n100_G40,n60_G200}.csv/.log`, `sim_mixture_main.csv/.log`,
`sim_mixture_size.csv/.log`, `sim_batchnull.csv/.log`, `section7_tables.md`. No file outside this
folder was touched; no project data were read; the one external call was a public GDC metadata
query (file counts per case, section 3 item 7). Nothing was deleted; the original `sim_grid_*.csv`
are kept.

Changes to the text, with reasons:

1. **Sign-test threshold corrected** (2.2, 4.1, 4.4): S < 0 ⇔ θ < ρ²/(1−ρ²), not θ < ρ². The f_X
   formula already carried the (1−ρ²); the prose did not. Verified on the grids: 152/153 cells agree
   with the corrected threshold on the selected-set θ, 138/153 with the old one. `sim_signtest.py`
   now uses the corrected threshold for its expected-sign bookkeeping; this affects only the
   `S_*_pg_power_sel` columns (not used in section 7), and the `sim_grid_v2_*.csv` files were produced
   before that edit (noted in the script).
2. **Primary null decided** (4.3, 2.4, abstract): within-operator (batch-stratified) permutation of the
   W rows is primary for z_N, z_D, z_S and for θ_pc; the unrestricted patient↔slide permutation is the
   sensitivity analysis; a disagreement rule is added. Reason: section 7.4 — with an operator-aligned
   RNA batch effect the unrestricted null reads pure H_post as H_err in 15/16 seeds, with accrual
   structure as "mixed" in 8/8; the within-operator null is centred in all configurations and loses
   1–4 z-units of power against |z| ≥ 13. Draws must be fresh (not `perms_<c>.npy`).
3. **θ of record = permutation-centred θ_pc under the within-operator null** (4.1/4.2, 7.4); the
   residualised-W θ is dropped (biased +0.07 to +0.14 under H_post with RNA batch); θ over all tested
   genes is dropped (unstable in a null-dominated mixture, 7.3). The confirmatory lung cohort on the
   discovery-pinned sets is named as the selection-free θ.
4. **Reading rules made numeric** (4.4 rewritten; 7.6): gate z_D ≥ 3; H_err z_N ≤ −3 and θ_pc ≤ −0.3;
   H_post −3 < z_N < 2 (< 3 for n_sel ≥ 1,000) with |θ_pc| ≤ 0.05/0.10; alignment z_N ≥ 2 split by
   θ_noise,min(λ_lo) = κλ_lo/(1−λ_lo) with the winner's-curse floor κ = 0.4 (n_sel ≥ 20) / 0.3; the
   reliability band [λ_lo, λ_hi] = [replicate ICC (lower bound; transferred by depth stratum where no
   replicates), λ_Poisson (upper bound)], fallback λ_lo = ρ̄²; the post-transcriptional share bound
   1 − (1−λ_lo)θ_pc/(κλ_lo). False-classification rates tabulated in 7.6.
5. **Replicate ICC re-read as a lower bound on λ** (3 item 7), on the public-metadata finding that
   multi-file CPTAC-3 cases are distinct tumour samples (not re-runs), so the ICC includes sampling
   heterogeneity; the earlier wording ("direct estimate … bounds λ from below") is kept in substance
   and made explicit, with the case-level reading under H_noise and the k-file averaging formula.
6. **Large-set calibration caveat** (7.5, 4.4 R1): the in-sample ρ̂ couples test-fold M into r_p
   (O(1/n)); at n_sel ≈ 1,500 the set-level z_N under H_post reaches +1.1 to +2.8. Hence z ≥ 3 for
   sets ≥ 1,000 genes, |θ_pc| ≤ 0.05 as the operative criterion, and the fold-wise (OOF) baseline
   residual as the recommended implementation.
7. **Section 5 and 6 extended** with the mixture and batch designs and the limits of the simulations;
   abstract item (4) updated and item (5) added; 4.5 compute note updated (local route with operator
   labels; `perms_<c>.npy` not reused).
8. **Grid files**: `sim_grid_v2_n100_G200.csv` reproduces `sim_grid_n100_G200.csv` exactly on all
   shared columns (max |diff| = 0); the v2 n60/G40 files differ from the older ones (those were made
   with an earlier script version) and are the ones cited in section 7.
