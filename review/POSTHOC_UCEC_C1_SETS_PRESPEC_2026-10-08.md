# Post hoc C1-UCEC re-reads with alternative discovery selection sets: pre-specification (2026-10-08)

**Status.** Written 2026-10-08 (all times UTC unless stated), after the published C1-UCEC
reading (2026-09-16), the D4 stratified re-read (2026-10-01) and the post hoc batch-axis
analyses (`POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md`). On 2026-10-08 the corresponding
author asked for the UCEC confirmatory checks to be completed and transferred the
discovery tables listed in section 1. Frozen by the sha256 of this file, the script and
every input, listed in the sidecar `POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md.sha256`;
released publicly on GitHub before the single full run (section 5). An adversarial review
of the first draft (same day, before freezing) is reflected throughout.

## 0. Scope, exceptions and prior exposure

- **Post hoc, reported whatever the result.** Nothing below was pre-registered. After
  freezing, no set, arm, background, threshold, reference or reading is added, dropped or
  changed.
- **No effect on the registered or published readings.** The published C1-UCEC reading
  (AUC 0.597; unrestricted null 0.498 +/- 0.036, p = 0.004) and the registered C1-LUAD
  reading stand as reported. Nothing here is pooled with either.
- **Exceptions to the UCEC rule.** Every set below re-derives the selection set, which
  section 2 of `C1_UCEC_FROZEN_RULE_2026-09-15.md` fixed as not re-derived. S6 also
  contradicts Addendum 2, item 2, which excluded re-deriving it on FFPE-only discovery
  cases. The D4 re-read's handling of 8 cases whose operator field is the all-zero UUID (one
  operator stratum) is inherited unchanged.
- **What S6 asks.** Every confirmatory slide that enters C1-UCEC is FFPE (169 slides; all
  465 PathDB slides of the confirmatory cohort are FFPE), so the published set has already
  been read on FFPE confirmatory slides. Of the 100 discovery cases, 22 rest on a single OCT
  section; these 22 are also the Ukrainian accrual and nearly one operator stratum
  (operator nested in OCT, 0.99; supplement, embedding-medium note). In the post hoc
  embedding-medium analysis the FFPE-only discovery arm was indeterminate at the count level
  (rank 9 of 20; alarm not triggered) and its median increment on the published set fell to
  -0.43 times the original, below 16 of 19 random-removal arms. S6 therefore asks whether
  the *discovery selection* depends on that OCT/Ukrainian/operator stratum; it does not
  re-test the confirmatory slides.
- **Prior exposure (disclosed; times UTC).**
  1. S1 and S2 are not blind. 2026-10-07T18:03 (agent `ab2b9ddcf0726ff90`, workflow
     `wf_9d75bc4c-fc4`): unrestricted AUC without a null, 0.649 (198 genes) and 0.658
     (248 genes). 2026-10-07T19:13 (agent `af73eb6ff5f6d29f2`, workflow `wf_93af2dfb-fc0`;
     also read by agent `af2d7ec99953d1d0a`): AUC, null mean and SD and p under all three
     arms, all p = 0.001 (residualised AUC for S1 0.633), for two backgrounds: every other
     usable gene, and the published 7,610 background genes only (AUC 0.646-0.675). The
     background fixed in section 3 was chosen after both had been seen; its justification
     is the C1-LUAD H2 rule, not those values.
  2. On the same confirmatory matrix, agent `ab2b9ddcf0726ff90` also computed unrestricted
     AUCs for discovery top-k sets (k = 100/250/500/1,000: 0.541/0.579/0.565/0.575) and for
     the selected sets of the other four discovery organs (0.564-0.585; 0.539-0.582 after
     removing UCEC-selected genes). The cross-organ values were read by the session that
     wrote this file.
  3. S3, S4, S5, S6, S7, all Sk-new and D01-D19 are blind: no set-level confirmatory
     statistic (AUC, rho, null, p) has been computed for them by anyone. The per-gene
     confirmatory increments themselves have been in the published gene-level table since
     2026-09-16.

## 1. Inputs

All paths relative to `MorphoResidual_paper/`; sha256 of each in the sidecar.

- Confirmatory matrices: `server_export/scripts/pinned/c1_ucec_d4_stratified.npz`
  (sha256 `ef75f1c0...`, = the value in `c1_ucec_d4_followup.json`): 10,146 usable genes x
  1,001 columns (column 0 observed confirmatory incremental R^2, columns 1-1000
  patient<->slide permutations), arms `M_unrestricted`, `M_within_operator`,
  `M_within_operator_residualized`.
- Published reading, for validation: `server_export/scripts/pinned/c1_ucec_d4_followup.json`,
  `server_export/pinned/ucec_c1/c1_ucec_gene_level.csv`,
  `server_export/pinned/ucec/residual_results_tumoronly.csv`.
- Discovery tables (each transferred unchanged from
  `/public/home/fjhui/ZW/ucec/results/...`; for the data-subset arms the sha256 must equal
  the value recorded in `server_export/pinned/posthoc/posthoc_summary_ucec.csv`):
  - operator axis: `server_export/results/ucec__sitepack_operator.csv`;
  - embedding-medium axis: `server_export/pinned/posthoc/ucec_medium.csv`;
  - TMT-plex axis: `server_export/pinned/posthoc/ucec_plex.csv`;
  - arm F (FFPE-only discovery): `server_export/pinned/posthoc/ucec_ffpe_only.csv`
    (sha256 `7801a033...`; 10,405 tested, 1,567 selected, 950 shared with the published
    set, median increment of the published set -0.051, all as in the summary);
  - arm Q (PDC Disqualified aliquots excluded): `server_export/pinned/posthoc/ucec_arm_Q.csv`;
  - arms D01-D19 (22 random discovery cases removed each):
    `server_export/pinned/posthoc/ucec_arm_D01.csv` ... `ucec_arm_D19.csv`.

## 2. Sets

**Inclusion rule.** Every UCEC discovery selection set whose count the manuscript or its
supplementary tables report as an alternative to the published 2,566-gene set: Table 1's
Batch-corr. (S1), Strictest (S2) and Plex-corr. (S4) columns, the embedding-medium count
(S3), the plex `both` count in Supplementary Table posthoc (S5), and the data-subset arms F
(S6) and Q (S7). Not reported, hence not included: the scan-quarter axis. The D arms are the
reference for S6, not sets in their own right.

Selected = `fdr < 0.05` and increment `> 0` on the named column pair (the operator-axis
`batchresid` pair is the C1-LUAD H2 definition). "Universe" = usable confirmatory genes
tested in that discovery analysis; "shared" = genes also in the published set.

| Set | Role | Discovery analysis | Columns | Selected | In universe | Shared | New | Universe |
|---|---|---|---|---|---|---|---|---|
| S6 | primary | arm F, FFPE-only cases (78) | `fdr`, `incremental_r2` | 1,567 | 1,550 | 940 | 610 | 10,095 |
| S1 | secondary | operator axis | `fdr_batchresid`, `incremental_r2_batchresid` | 201 | 198 | 171 | 27 | 10,146 |
| S2 | secondary | operator axis | `fdr_both`, `incremental_r2_both` | 253 | 248 | 200 | 48 | 10,146 |
| S3 | secondary | embedding-medium axis | `fdr_batchresid`, `incremental_r2_batchresid` | 847 | 833 | 738 | 95 | 10,146 |
| S4 | secondary | TMT-plex axis | `fdr_batchresid`, `incremental_r2_batchresid` | 2,124 | 2,112 | 1,725 | 387 | 10,146 |
| S5 | secondary | TMT-plex axis | `fdr_both`, `incremental_r2_both` | 2,358 | 2,328 | 1,749 | 579 | 10,146 |
| S7 | secondary | arm Q | `fdr`, `incremental_r2` | 2,360 | 2,328 | 1,720 | 608 | 10,135 |
| S3both | untestable | embedding-medium axis | `fdr_both`, `incremental_r2_both` | 0 | | | | |

(The S7 cells were filled on 2026-10-08 from the transferred arm-Q table, discovery side
only, before freezing. All 20 transferred tables, Q and D01-D19, match the sha256 and
selected counts recorded in `posthoc_summary_ucec.csv`.)

- **Sk-new (descriptive).** Sk minus the published set, against the background
  universe minus Sk (the shared genes are removed from both sides, so Sk-new and the shared
  part face the same background). AUC only; rho is not reported for Sk-new.
- **D01-D19 (reference for S6).** Each computed exactly as S6.

## 3. Statistic and nulls

- **Background.** Every other gene in the set's universe, including published-set genes
  not in the set (C1-LUAD H2 rule). No other background is computed.
- **Statistic.** AUC of the confirmatory incremental R^2 of the set against the background
  (rank Mann-Whitney) by `assemble`, copied unchanged from `c1_run_test.py`. Secondary:
  Spearman rho between that analysis's discovery increment column and the confirmatory
  increment over the universe. Margin = observed AUC minus that arm's null mean (null means
  differ between sets, so margins, not raw AUCs, are compared).
- **Nulls.** The three saved arms, B = 1,000 each, p = (1 + #{null >= observed})/(1 + B),
  one-sided; floor 1/1,001. Unrestricted: patient<->slide permutation. Within-operator: the
  same permutations kept inside operator strata (its observed AUC equals the unrestricted
  one). Residualised: morphology PCs residualised on dummies for the 12 largest operator
  levels with at least 3 cases (83 of 138 cases covered, including the all-zero-UUID
  stratum), increments recomputed, same within-operator permutations; its observed AUC
  differs. It is reported descriptively, as in D4. This departs from C1-LUAD H2, whose
  decision used residualised PCs with the within-operator null.
- **Not stratified by plex or medium.** No confirmatory null is stratified by TMT plex or
  embedding medium; a set that holds says nothing about the confirmatory ranking's
  robustness to those axes.

## 4. Reading rules and manuscript wording, fixed now

Tiers, as in the registered rule's section 5: p < 0.05 significant; 0.05 <= p < 0.15
suggestive; p >= 0.15 not significant. p_u = unrestricted, p_w = within-operator,
n_below = number of the 19 D arms whose within-operator margin is below S6's.

**Primary: S6.** Results, in the C1-UCEC paragraph (after the D4 sentences), labelled post
hoc, exactly one of:
- (A) p_u < 0.05 and p_w < 0.05: "Post hoc, a selection set re-derived from the 78
  discovery cases whose slides are FFPE-prepared (1,550 genes in the confirmatory universe)
  also ranked higher in the confirmatory cohort (AUC [a]; within-operator null [m] +/- [s],
  p = [p_w])."
- (B) p_u < 0.05 and p_w >= 0.05: "Post hoc, a selection set re-derived from the 78
  discovery cases whose slides are FFPE-prepared ranked higher in the confirmatory cohort
  under the unrestricted null (AUC [a], p = [p_u]) but not under the within-operator null
  (p = [p_w]). The within-stratum null keeps every case inside its stratum, so any
  association carried by differences between strata, including genuine biology aliased
  with accrual, is credited to chance, and the fixed cases lower its power. A
  non-significant result therefore means the data cannot separate the signal from stratum
  structure, not that it is an artefact."
- (C) p_u >= 0.05: "Post hoc, a selection set re-derived from the 78 discovery cases whose
  slides are FFPE-prepared was [suggestive / not significant] in the confirmatory cohort
  (AUC [a], p = [p_u])."

followed in every case by: "Its margin over the within-operator null ([margin]) exceeded
that of [n_below] of 19 re-derivations that each removed 22 random discovery cases", and
exactly one of:
- n_below = 0: "; restricting discovery to the FFPE-prepared cases lowered the
  confirmatory ranking more than any random removal of 22 cases, consistent with the
  published selection drawing on the OCT/Ukrainian/operator stratum."
- 1 <= n_below <= 18: "; the change is within the range of random removals and cannot be
  attributed to the OCT/Ukrainian/operator stratum rather than to the smaller discovery
  sample."
- n_below = 19: "; restricting discovery to the FFPE-prepared cases did not lower the
  confirmatory ranking relative to random removals."

Abstract: unchanged, except in case (C) with n_below = 0, where the UCEC clause gains
"; this did not hold when discovery was restricted to FFPE-prepared cases (post hoc)". No
positive post hoc result is added to the abstract.

**Descriptive: S6new** (supplement only): if S6 is (A) or (B) and S6new p_u >= 0.05: "The
ranking was carried by the 940 genes shared with the published set; the 610 genes new to
the FFPE-only set did not rank higher on their own (AUC [a], p = [p_u])." If S6new
p_u < 0.05: "including the 610 genes not in the published set (AUC [a], p = [p_u])." If S6
is (C), S6new is reported in the table only.

**Secondary family: S1-S5, S7.** p_set = max(p_u, p_w) (both nulls must pass); Holm over
the six; a set "holds" if its Holm-adjusted p_set < 0.05. Raw p_u, p_w, p_set, Holm value,
margin, overlap and Sk-new rows go to the supplement table. Main text, one sentence after
the S6 sentences: "The other batch-corrected and data-subset discovery sets reported in
the paper ([h] of 6) also ranked higher in the confirmatory cohort under both nulls
(Holm-adjusted)[; [names] did not]; they share 74-89% of their genes with the published
set, so these are not independent replications." Supplement adds the sentence on plex and
medium from section 3. No abstract change in any case.

**Untestable:** S3both is reported as 0 genes selected.

## 5. Freeze and run order

1. `python server_export/scripts/c1_ucec_posthoc_sets.py --validate-only` (reads only the
   published set). 2026-10-08, before freezing: V1, all three arms reproduce
   `c1_ucec_d4_followup.json` to <= 1e-12 (AUC 0.5974425565 / null 0.4981945451 / SD
   0.0362962998 / p 0.003996; within-operator null 0.5059362959; residualised AUC
   0.5671994242, p 0.005994; rho 0.2935675639), n 2,536 / 7,610; V2, the new-set code path
   reproduces the published mask and numbers. VALIDATED.
2. The S7 cells of section 2 were filled from the transferred arm-Q table (discovery side
   only); the sidecar is then written, and this file, the script and the sidecar are
   released on GitHub.
3. One full run; it refuses to start unless every sidecar hash matches. Outputs:
   `server_export/results/c1_ucec_posthoc_sets/result.{csv,json}`. The manuscript sentences
   of section 4 are then filled in exactly as written.
