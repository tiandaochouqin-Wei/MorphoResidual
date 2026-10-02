#!/usr/bin/env python3
"""
Step 0b (2026-09-29) branch "protein", part 2/2: discovery-only protein-route
equivalence check. Runs on the HPC, where the discovery CDAP TSV lives.

Compares two protein routes for the LUAD **DISCOVERY** cohort (PDC000153) ONLY --
never PDC000489 (confirmatory):

  route "pdc"  : PDC GraphQL quantDataMatrix(pdc_study_id="PDC000153",
                 data_type="log2_ratio"), pulled locally (not here) by
                 review/C1_LUAD_STEP0B_2026-09-29/pull_discovery_protein_via_pdc.py
                 and uploaded to this host as a case-indexed CSV (--pdc-csv).
                 This is the SAME query method c1_pull_omics.py uses for the
                 confirmatory PDC000489 pull, exercised here on discovery data.
  route "cdap" : residual_analysis.load_protein_matrix() under the LUAD
                 discovery environment (`source morpho_env.sh luad`), i.e. the
                 CDAP '<aliquot> Log Ratio' TSV columns, excluding 'Unshared'.
                 This is what the LUAD discovery pipeline actually used to
                 produce server_export/results/luad__residual_results_tumoronly.csv.

Neither route is treated as ground truth; this reports where and how much they
agree, plus exactly how each route handles duplicate gene symbols -- it does not
edit c1_pull_omics.py or residual_analysis.py (forbidden-edit list), it only
imports residual_analysis by reference.

Usage (real, on HPC, after `source morpho_env.sh luad`):
  python protein_route_check.py --pdc-csv pinned/luad_discovery_protein_via_pdc.csv \
      --tested-genes-csv pinned/luad/residual_results_tumoronly.csv \
      --out-json .../protein_route_check.json --out-csv .../protein_route_check_pergene.csv

Usage (dry run, anywhere, no HPC files needed):
  python protein_route_check.py --self-test
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FORBIDDEN_STUDY = "PDC000489"


def load_pdc_matrix(csv_path):
    df = pd.read_csv(csv_path, index_col=0)
    df.index = df.index.astype(str)
    return df


def load_cdap_matrix():
    """Real route B: residual_analysis.load_protein_matrix() under whatever
    MORPHO_ALIQUOT_XWALK / MORPHO_PROTEIN_TSV morpho_env.sh luad has exported.
    Imported, not copy-pasted, so any future fix to that function is picked up
    automatically -- and so this script cannot silently drift from what the
    discovery pipeline actually ran."""
    sys.path.insert(0, str(HERE))
    import residual_analysis as RA  # noqa: E402  (reused; file not edited)
    mat = RA.load_protein_matrix()
    mat.index = mat.index.astype(str)
    return mat


def load_tested_genes(csv_path):
    df = pd.read_csv(csv_path, usecols=["gene"])
    genes = set(df["gene"].astype(str))
    return genes


def _dedupe_columns_for_comparison(df, label):
    """Average any duplicate-labelled columns so the comparison itself is not
    broken by upstream duplicate-handling differences (see report note). Returns
    (deduped_df, sorted list of duplicate gene labels found)."""
    dup_mask = df.columns.duplicated(keep=False)
    dup_genes = sorted(set(df.columns[dup_mask]))
    if not dup_genes:
        return df, dup_genes
    print(f"  [{label}] {len(dup_genes)} duplicate-labelled gene column(s) found "
          f"(e.g. {dup_genes[:5]}); averaging duplicates for this comparison only "
          f"-- the underlying loader is not changed.")
    deduped = df.T.groupby(level=0).mean().T
    return deduped, dup_genes


def compute_comparison(pdc_df, cdap_df, tested_genes, min_pairs_for_corr=10):
    assert FORBIDDEN_STUDY not in getattr(pdc_df, "_source_study", ""), (
        "FATAL: pdc_df must never originate from the confirmatory study")

    pdc_df, pdc_dup_genes = _dedupe_columns_for_comparison(pdc_df, "pdc")
    cdap_df, cdap_dup_genes = _dedupe_columns_for_comparison(cdap_df, "cdap")

    pdc_genes = set(pdc_df.columns)
    cdap_genes = set(cdap_df.columns)
    union_genes = pdc_genes | cdap_genes
    inter_genes = pdc_genes & cdap_genes
    jaccard = len(inter_genes) / len(union_genes) if union_genes else float("nan")

    tested_missing_from_pdc = sorted(tested_genes - pdc_genes)
    tested_missing_from_cdap = sorted(tested_genes - cdap_genes)

    pdc_cases = set(pdc_df.index)
    cdap_cases = set(cdap_df.index)
    shared_cases = sorted(pdc_cases & cdap_cases)
    shared_genes = sorted(inter_genes)

    pdc_zeros = int((pdc_df.values == 0).sum())
    cdap_numeric = cdap_df.apply(pd.to_numeric, errors="coerce")
    cdap_sentinel_cells = int(cdap_numeric.isna().sum().sum() - cdap_df.isna().sum().sum())
    cdap_zeros = int((cdap_numeric.values == 0).sum())

    result = {
        "n_pdc_cases": len(pdc_cases), "n_cdap_cases": len(cdap_cases),
        "n_shared_cases": len(shared_cases),
        "cases_only_pdc": sorted(pdc_cases - cdap_cases),
        "cases_only_cdap": sorted(cdap_cases - pdc_cases),
        "n_pdc_genes": len(pdc_genes), "n_cdap_genes": len(cdap_genes),
        "gene_jaccard": jaccard,
        "n_tested_genes": len(tested_genes),
        "n_tested_missing_from_pdc": len(tested_missing_from_pdc),
        "tested_missing_from_pdc": tested_missing_from_pdc,
        "n_tested_missing_from_cdap": len(tested_missing_from_cdap),
        "tested_missing_from_cdap": tested_missing_from_cdap,
        "pdc_duplicate_gene_labels": pdc_dup_genes,
        "cdap_duplicate_gene_labels": cdap_dup_genes,
        "cdap_duplicate_handling_note": (
            "residual_analysis.load_protein_matrix() does NOT deduplicate gene-symbol "
            "rows from the CDAP TSV: sub.T.groupby(level=0).mean() groups by the ROW "
            "index (case), not by gene column, so a repeated gene symbol in the TSV "
            "survives as two identically-labelled columns in the returned DataFrame. "
            "c1_pull_omics.py's quantDataMatrix route instead builds a plain dict "
            "keyed by gene name, so a repeated gene symbol there silently keeps only "
            "the LAST row (dict overwrite). The two routes therefore handle duplicate "
            "gene symbols differently by construction, not just by data content -- this "
            "check averages duplicates on both sides before comparing (see above) so "
            "that difference does not masquerade as a data disagreement."
        ),
        "pdc_exact_zeros": pdc_zeros,
        "cdap_exact_zeros": cdap_zeros,
        "cdap_non_numeric_sentinel_cells": cdap_sentinel_cells,
    }

    if not shared_cases or not shared_genes:
        result["warning"] = "no shared cases/genes -- delta and correlation stats skipped"
        return result, pd.DataFrame(columns=["gene", "in_pdc", "in_cdap", "n_pairs",
                                              "pearson_r", "median_abs_delta"])

    pdc_shared = pdc_df.loc[shared_cases, shared_genes]
    cdap_shared = cdap_df.loc[shared_cases, shared_genes]

    pdc_nan = pdc_shared.isna().values
    cdap_nan = cdap_shared.isna().values
    n_cells = pdc_nan.size
    both_nan = int((pdc_nan & cdap_nan).sum())
    neither_nan = int((~pdc_nan & ~cdap_nan).sum())
    pdc_only_nan = int((pdc_nan & ~cdap_nan).sum())
    cdap_only_nan = int((~pdc_nan & cdap_nan).sum())
    result["nan_pattern"] = {
        "n_cells": n_cells, "both_nan": both_nan, "neither_nan": neither_nan,
        "nan_only_in_pdc": pdc_only_nan, "nan_only_in_cdap": cdap_only_nan,
        "agreement_rate": (both_nan + neither_nan) / n_cells if n_cells else float("nan"),
    }

    delta = pdc_shared.values - cdap_shared.values
    valid = ~np.isnan(delta)
    abs_delta = np.abs(delta[valid])
    result["delta_over_shared_cells"] = {
        "n_valid": int(valid.sum()),
        "median_abs_delta": float(np.median(abs_delta)) if abs_delta.size else None,
        "p999_abs_delta": float(np.percentile(abs_delta, 99.9)) if abs_delta.size else None,
    }

    delta_df = pd.DataFrame(delta, index=shared_cases, columns=shared_genes)
    per_case_median_delta = delta_df.median(axis=1, skipna=True)
    result["per_case_median_delta"] = {c: (float(v) if pd.notna(v) else None)
                                        for c, v in per_case_median_delta.items()}

    per_gene_rows = []
    all_pdc_genes_sorted = sorted(union_genes)
    for gene in all_pdc_genes_sorted:
        in_pdc = gene in pdc_genes
        in_cdap = gene in cdap_genes
        row = {"gene": gene, "in_pdc": in_pdc, "in_cdap": in_cdap,
               "n_pairs": 0, "pearson_r": None, "median_abs_delta": None}
        if in_pdc and in_cdap and gene in shared_genes:
            x = pdc_shared[gene].values
            y = cdap_shared[gene].values
            pair_valid = ~np.isnan(x) & ~np.isnan(y)
            n_pairs = int(pair_valid.sum())
            row["n_pairs"] = n_pairs
            if n_pairs >= min_pairs_for_corr and np.std(x[pair_valid]) > 0 and np.std(y[pair_valid]) > 0:
                row["pearson_r"] = float(np.corrcoef(x[pair_valid], y[pair_valid])[0, 1])
            if n_pairs > 0:
                row["median_abs_delta"] = float(np.median(np.abs(x[pair_valid] - y[pair_valid])))
        per_gene_rows.append(row)
    pergene_df = pd.DataFrame(per_gene_rows)

    r_vals = pergene_df["pearson_r"].dropna().values
    result["per_gene_pearson"] = {
        "n_genes_with_r": int(len(r_vals)),
        "median_r": float(np.median(r_vals)) if len(r_vals) else None,
        "p05_r": float(np.percentile(r_vals, 5)) if len(r_vals) else None,
    }

    return result, pergene_df


def synthetic_cdap_from_pdc(pdc_df, rng, drop_cases=2, drop_genes=5,
                             offset_scale=0.05, noise_scale=0.02, extra_nan_frac=0.01):
    """Synthetic stand-in for RA.load_protein_matrix(), built by perturbing a
    copy of the real PDC route matrix with KNOWN, injected differences, so the
    comparison logic above can be checked against ground truth without the
    real CDAP TSV (which lives only on the HPC)."""
    cdap = pdc_df.copy().astype(float)
    dropped_cases = list(cdap.index[:drop_cases])
    cdap = cdap.drop(index=dropped_cases)
    dropped_genes = list(cdap.columns[:drop_genes])
    cdap = cdap.drop(columns=dropped_genes)

    offset = pd.Series(rng.normal(0, offset_scale, size=len(cdap.index)), index=cdap.index)
    noise = pd.DataFrame(rng.normal(0, noise_scale, size=cdap.shape),
                          index=cdap.index, columns=cdap.columns)
    cdap = cdap.add(offset, axis=0) + noise

    # Simulate a duplicate gene symbol surviving load_protein_matrix() (see the
    # duplicate-handling note above): append a second, noisy copy of one column
    # under the SAME label.
    dup_gene = cdap.columns[10]
    dup_col = cdap[[dup_gene]].copy()
    dup_col.iloc[:, 0] = dup_col.iloc[:, 0] + rng.normal(0, noise_scale, size=len(dup_col))
    cdap = pd.concat([cdap, dup_col], axis=1)  # now dup_gene appears twice

    # Simulate a few extra route-specific missing cells.
    mask = rng.random(cdap.shape) < extra_nan_frac
    cdap_vals = cdap.values.copy()
    cdap_vals[mask] = np.nan
    cdap = pd.DataFrame(cdap_vals, index=cdap.index, columns=cdap.columns)

    return cdap, dropped_cases, dropped_genes, offset


def self_test():
    print("=== SELF-TEST: synthetic stand-in for RA.load_protein_matrix() ===")
    pdc_csv = HERE / "pinned" / "luad_discovery_protein_via_pdc.csv"
    if not pdc_csv.exists():
        sys.exit(f"FATAL: self-test needs the real PDC-route CSV at {pdc_csv} "
                  f"(written by pull_discovery_protein_via_pdc.py) as its base matrix -- "
                  f"it perturbs a copy of it, it does not invent gene/case names from "
                  f"scratch, so the comparison exercises real gene/case identifiers.")
    pdc_df = load_pdc_matrix(pdc_csv)
    rng = np.random.RandomState(0)
    cdap_df, dropped_cases, dropped_genes, offset = synthetic_cdap_from_pdc(pdc_df, rng)
    tested_genes = set(pdc_df.columns[:50])  # arbitrary stand-in "tested" set

    result, pergene_df = compute_comparison(pdc_df, cdap_df, tested_genes)

    print(json.dumps({k: v for k, v in result.items()
                       if k not in ("tested_missing_from_pdc", "tested_missing_from_cdap",
                                     "per_case_median_delta")}, indent=2, default=str))

    ok = True
    if result["n_cdap_cases"] != pdc_df.shape[0] - len(dropped_cases):
        print("FAIL: cdap case count does not reflect the injected drop"); ok = False
    flagged_missing = set(pergene_df.loc[~pergene_df["in_cdap"], "gene"])
    if not (set(dropped_genes) <= flagged_missing):
        print("FAIL: injected dropped genes not correctly flagged as in_cdap=False"); ok = False
    if not cdap_df.columns.duplicated().any():
        print("FAIL: injected duplicate gene column not present in synthetic cdap_df"); ok = False
    if not result["cdap_duplicate_gene_labels"]:
        print("FAIL: duplicate-gene detection did not fire on the synthetic duplicate"); ok = False
    r_med = result["per_gene_pearson"]["median_r"]
    if r_med is None or r_med < 0.9:
        print(f"FAIL: median per-gene Pearson r ({r_med}) unexpectedly low for a "
              f"noise-only perturbation (offset_scale=0.05, noise_scale=0.02)"); ok = False
    observed_case_offsets = pd.Series(result["per_case_median_delta"]).dropna()
    # delta = pdc - cdap = pdc - (pdc + offset + noise) = -offset - noise (median~-offset)
    corr_with_injected = np.corrcoef(
        observed_case_offsets.reindex(offset.index).values, (-offset).values)[0, 1]
    if corr_with_injected < 0.8:
        print(f"FAIL: recovered per-case median delta does not track the injected "
              f"per-case offset (corr={corr_with_injected:.3f})"); ok = False
    else:
        print(f"  per-case median delta recovers the injected offset: corr={corr_with_injected:.3f}")

    print("SELF-TEST " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdc-csv")
    ap.add_argument("--tested-genes-csv")
    ap.add_argument("--out-json")
    ap.add_argument("--out-csv")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    for name in ("pdc_csv", "tested_genes_csv", "out_json", "out_csv"):
        if getattr(args, name) is None:
            sys.exit(f"FATAL: --{name.replace('_', '-')} is required outside --self-test")

    print(f"[route pdc] loading {args.pdc_csv}")
    pdc_df = load_pdc_matrix(args.pdc_csv)
    print(f"  {pdc_df.shape[0]} cases x {pdc_df.shape[1]} genes")

    print("[route cdap] loading via residual_analysis.load_protein_matrix() "
          "(requires `source morpho_env.sh luad` to have been run first)")
    cdap_df = load_cdap_matrix()
    print(f"  {cdap_df.shape[0]} cases x {cdap_df.shape[1]} gene columns "
          f"({cdap_df.columns.duplicated().sum()} duplicate-labelled)")

    print(f"[tested genes] loading {args.tested_genes_csv}")
    tested_genes = load_tested_genes(args.tested_genes_csv)
    print(f"  {len(tested_genes)} tested genes")

    result, pergene_df = compute_comparison(pdc_df, cdap_df, tested_genes)

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)
    pergene_df.to_csv(args.out_csv, index=False)
    print(f"\nwrote {args.out_json}")
    print(f"wrote {args.out_csv}")
    print(json.dumps({k: v for k, v in result.items()
                       if not isinstance(v, list) and k != "per_case_median_delta"},
                      indent=2, default=str))


if __name__ == "__main__":
    main()
