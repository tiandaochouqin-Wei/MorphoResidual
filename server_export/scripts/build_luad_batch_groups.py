#!/usr/bin/env python3
"""
Build server_export/scripts/pinned/luad_batch_groups.csv -- the D2 P0 file
for the LUAD discovery split-half pre-check (C1-LUAD Step 0b, 2026-09-29,
metadata-only; see review/C1_LUAD_STEP0B_2026-09-29/ and the frozen-rule
draft's D2 section, review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md).

Local files only. No network access. Does not import or edit
batch_leak_check.py, split_replication_power.py, or any other already-run
script named in the Step 0b task -- it reimplements their documented logic
(see make_batch_groups.py's docstring for exactly what is reproduced and what
is NOT: the emb/*.pt "did this slide actually get embedded" filter is not
available locally and is not applied here; see the printed NOTE).

Schema (8 columns as of the frozen-rule draft v3 rewrite, 2026-09-29; the
first 7 match ucec_batch_groups.csv's schema exactly, the 8th is LUAD-only):
  case_submitter_id,
  q_early_vs_late      -- sel: scan_quarter in {2016Q3,2016Q4,2017Q1,2017Q2}
                           test: scan_quarter in {2017Q3,2017Q4}
  q_2017Q3_vs_rest      -- test: scan_quarter == 2017Q3; sel: everything else
  op_balanced_0..4      -- greedy operator-balanced sel/test, seeds 0-4
                           (make_batch_groups.balance_operators)
  stream_c3n_vs_c3l     -- DESCRIPTIVE ONLY (D2 P0 / D6d), not one of the
                           7 UCEC-schema columns and never a GO/MARGINAL/NO-GO
                           gate (review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md
                           D2 P2). sel: C3N-prefixed case_submitter_id (71 of
                           the 105 analysis cases); test: C3L-prefixed (34).
                           Pure string check on the ID itself
                           (make_batch_groups.stream_side) -- no metadata
                           lookup, so no case is ever blank here.
  This 8th column was appended after the first 7 were already built and
  pinned (2026-09-29 07:27); rebuilding must not perturb columns 1-7's values,
  algorithm, or seeds -- only append column 8.

Analysis case set: the 105 case_submitter_ids in
figures/figdata/scores_luad.csv (verified count == 105 below; this is the
discovery tri-modal set per review/C1_LUAD_STEP0_2026-09-27.md Sec.5 / Sec.1,
n=105 modal in luad__residual_results_tumoronly.csv). All 105 are present in
slide_type_map_luad.tsv. Any case with no operator_uuid / scan_quarter on any
of its Primary Tumor slides is left blank in every affected column and listed
below -- it is not dropped from the case_submitter_id column.
"""
import hashlib
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_batch_groups as mbg  # local, new module written for Step 0b

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PINNED = os.path.join(HERE, "pinned")
OUT_CSV = os.path.join(PINNED, "luad_batch_groups.csv")
OUT_SHA = OUT_CSV + ".sha256"

SCORES_LUAD = os.path.join(ROOT, "figures", "figdata", "scores_luad.csv")

EARLY_Q = {"2016Q3", "2016Q4", "2017Q1", "2017Q2"}
LATE_Q = {"2017Q3", "2017Q4"}

# The 7 columns already pinned (2026-09-29 07:27) before the frozen-rule
# draft v3 rewrite added stream_c3n_vs_c3l as an 8th column. This rebuild must
# not perturb these 7 -- see the byte-identity guard in main().
LEGACY_7_COLS = ["case_submitter_id", "q_early_vs_late", "q_2017Q3_vs_rest",
                  "op_balanced_0", "op_balanced_1", "op_balanced_2",
                  "op_balanced_3", "op_balanced_4"]


def main():
    # --- byte-identity guard for columns 1-7 (fail-closed) -----------------
    # This rebuild's only change is to APPEND an 8th column,
    # stream_c3n_vs_c3l. Columns 1-7 (and row order) must come out identical,
    # byte for byte, to the file already pinned on disk. There is no git
    # history here to diff against (this project directory is not a git
    # repo), so snapshot the pre-existing file now, before it is overwritten
    # below, and compare after rebuilding, before writing anything.
    old_df = None
    if os.path.exists(OUT_CSV):
        old_df = pd.read_csv(OUT_CSV, dtype=str).fillna("")
        print(f"[byte-identity guard] snapshotted pre-existing {OUT_CSV} "
              f"({len(old_df)} rows, {len(old_df.columns)} columns) for "
              f"post-build comparison.")
    else:
        print(f"[byte-identity guard] no pre-existing {OUT_CSV} found -- "
              f"nothing to compare columns 1-7 against.")

    # --- analysis case set ---------------------------------------------
    scores = pd.read_csv(SCORES_LUAD, index_col=0)
    cases = sorted(scores.index.astype(str))
    n_cases = len(cases)
    print(f"[case set] {SCORES_LUAD}: {n_cases} rows/case_submitter_ids")
    if n_cases != 105:
        print(f"  DISCREPANCY: expected 105 (review/C1_LUAD_STEP0_2026-09-27.md "
              f"Sec.5/Sec.1), got {n_cases}. Not adjusting anything -- reporting "
              f"as-is.")
    else:
        print("  matches the task's stated 105. OK.")

    smap = pd.read_csv(mbg.slide_map_path("luad"), sep="\t", dtype=str)
    slidemap_cases = set(smap["case_submitter_id"])
    missing_from_slidemap = [c for c in cases if c not in slidemap_cases]
    if missing_from_slidemap:
        print(f"  NOTE: {len(missing_from_slidemap)} of the {n_cases} analysis "
              f"cases are not in slide_type_map_luad.tsv at all: "
              f"{missing_from_slidemap}")
    else:
        print(f"  all {n_cases} analysis cases are present in "
              f"slide_type_map_luad.tsv. OK.")

    # --- per-case labels (operator, quarter) -----------------------------
    op_label, op_no_label_all = mbg.load_case_operator_labels("luad")
    q_label, q_no_label_all = mbg.load_case_quarter_labels("luad")
    op_no_label = [c for c in cases if c not in op_label]
    q_no_label = [c for c in cases if c not in q_label]
    print(f"[labels] of the {n_cases} analysis cases: "
          f"{len(op_no_label)} have no operator_uuid label, "
          f"{len(q_no_label)} have no scan_quarter label.")
    if op_no_label:
        print(f"  no operator label: {op_no_label}")
    if q_no_label:
        print(f"  no scan_quarter label: {q_no_label}")
    print("  LIMITATION (repeated from make_batch_groups.py): these labels use "
          "ALL Primary Tumor slides per case from slide_acquisition_meta.tsv / "
          "slide_type_map_luad.tsv. batch_leak_check.case_batch_labels further "
          "restricts to slides actually present as embedding files under "
          "<root>/WSI/emb on the HPC filesystem, which this local machine "
          "cannot see. The real gate for any P1/P2 run must intersect this "
          "file's case_submitter_id column against split_replication_power.py's "
          "own rna/protein/wsi-derived `common` set on HPC before use -- this "
          "file's 105-row case list is the local scores_luad.csv proxy for that "
          "set, not a verified reproduction of it.")

    n_op_labelled = n_cases - len(op_no_label)
    distinct_ops = sorted(set(op_label[c] for c in cases if c in op_label))
    print(f"[operators] {len(distinct_ops)} distinct operators among the "
          f"{n_op_labelled} operator-labelled analysis cases.")

    # --- quarter columns (direct rule, no inference needed) --------------
    def q_early_late(c):
        if c not in q_label:
            return ""
        q = q_label[c]
        if q in EARLY_Q:
            return "sel"
        if q in LATE_Q:
            return "test"
        return ""  # should not happen -- every observed LUAD quarter is one of the six above

    def q_2017q3(c):
        if c not in q_label:
            return ""
        return "test" if q_label[c] == "2017Q3" else "sel"

    df = pd.DataFrame({"case_submitter_id": cases})
    df["q_early_vs_late"] = [q_early_late(c) for c in cases]
    df["q_2017Q3_vs_rest"] = [q_2017q3(c) for c in cases]

    # --- op_balanced_0..4 via the generic builder ------------------------
    op_df, _ = mbg.build("luad", cases=cases)
    for s in mbg.SEEDS:
        col = f"op_balanced_{s}"
        df[col] = op_df.set_index("case_submitter_id").loc[cases, col].values

    # --- stream_c3n_vs_c3l (D2 P0 / D6d, descriptive only) -----------------
    # Pure function of the case_submitter_id prefix -- no metadata lookup, so
    # this column cannot be blank for any of the 105 analysis cases (checked
    # below). Appended as column 8; does not touch columns 1-7 above.
    df["stream_c3n_vs_c3l"] = [mbg.stream_side(c) for c in cases]
    n_stream_blank = int((df["stream_c3n_vs_c3l"] == "").sum())
    if n_stream_blank:
        blank_ids = df.loc[df["stream_c3n_vs_c3l"] == "", "case_submitter_id"].tolist()
        print(f"  DISCREPANCY: {n_stream_blank} case_submitter_id(s) match neither "
              f"the C3N- nor C3L- prefix: {blank_ids}. Not adjusting anything -- "
              f"reporting as-is.")
    else:
        print(f"  all {n_cases} analysis cases match C3N- or C3L- exactly. OK.")

    df = df[["case_submitter_id", "q_early_vs_late", "q_2017Q3_vs_rest",
             "op_balanced_0", "op_balanced_1", "op_balanced_2", "op_balanced_3",
             "op_balanced_4", "stream_c3n_vs_c3l"]]

    # --- per-column checks -------------------------------------------------
    print(f"\n[per-column report] N analysis cases = {n_cases}")
    op_sizes_by_op = {}
    for op in distinct_ops:
        op_sizes_by_op[op] = [c for c in cases if op_label.get(c) == op]

    any_floor_fail = False
    for col in df.columns[1:]:
        vc = df[col].value_counts()
        n_sel, n_test = int(vc.get("sel", 0)), int(vc.get("test", 0))
        n_blank = int(vc.get("", 0)) if "" in vc.index else int((df[col] == "").sum())
        floor_ok = n_sel >= 40 and n_test >= 30
        any_floor_fail = any_floor_fail or not floor_ok

        if col.startswith("op_balanced"):
            # operator split-check: does any operator appear on both sides?
            split_ops = []
            for op, members in op_sizes_by_op.items():
                sides = set(df.set_index("case_submitter_id").loc[members, col])
                sides.discard("")
                if len(sides) > 1:
                    split_ops.append(op)
            n_ops_sel = len({op_label[c] for c in cases
                              if df.set_index("case_submitter_id").loc[c, col] == "sel"})
            n_ops_test = len({op_label[c] for c in cases
                               if df.set_index("case_submitter_id").loc[c, col] == "test"})
            print(f"  {col}: n_sel={n_sel} n_test={n_test} n_blank={n_blank} "
                  f"operators_sel={n_ops_sel} operators_test={n_ops_test} "
                  f"any_operator_split={bool(split_ops)} floor_ok={floor_ok}")
            if split_ops:
                print(f"    SPLIT OPERATORS (should be empty): {split_ops}")
        elif col == "stream_c3n_vs_c3l":
            # Descriptive only (D2 P0 / D6d) -- never a GO/MARGINAL/NO-GO gate
            # (frozen-rule draft D2 P2). Rule text pins n_sel=71/n_test=34
            # exactly (105 analysis cases = 71 C3N + 34 C3L); flag, don't
            # silently accept, any deviation from that exact count.
            expect_ok = (n_sel == 71 and n_test == 34)
            print(f"  {col}: n_sel={n_sel} n_test={n_test} n_blank={n_blank} "
                  f"floor_ok={floor_ok} descriptive_only=True "
                  f"matches_rule_71_34={expect_ok}")
            if not expect_ok:
                print(f"    DISCREPANCY: expected n_sel=71/n_test=34 exactly "
                      f"(review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md D2 P0), "
                      f"got n_sel={n_sel}/n_test={n_test}. Not adjusting anything "
                      f"-- reporting as-is.")
        else:
            print(f"  {col}: n_sel={n_sel} n_test={n_test} n_blank={n_blank} "
                  f"floor_ok={floor_ok}")

    if any_floor_fail:
        print("  DISCREPANCY: at least one column fails the split_replication_power.py "
              "L181 floor (n_sel>=40, n_test>=30). Not adjusting the algorithm to force "
              "a pass -- reporting as-is.")

    blank_any = df[(df.iloc[:, 1:] == "").any(axis=1)]["case_submitter_id"].tolist()
    print(f"\n[cases with >=1 blank label] {blank_any}")

    # --- byte-identity check: columns 1-7 unchanged (fail-closed) ----------
    if old_df is not None:
        legacy_cols_on_disk = [c for c in old_df.columns if c != "stream_c3n_vs_c3l"]
        if legacy_cols_on_disk != LEGACY_7_COLS:
            print(f"  DISCREPANCY: pre-existing file's non-stream columns "
                  f"{legacy_cols_on_disk} do not match the expected 7 "
                  f"{LEGACY_7_COLS}. Skipping the byte-identity comparison "
                  f"and proceeding (nothing pinned to compare against).")
        elif list(old_df["case_submitter_id"]) != list(df["case_submitter_id"]):
            sys.exit(f"FATAL: row order (case_submitter_id) differs between the "
                      f"pre-existing {OUT_CSV} and this rebuild -- refusing to "
                      f"overwrite a pinned file with a reordered one.")
        else:
            new_legacy = df[LEGACY_7_COLS].astype(str)
            old_legacy = old_df[LEGACY_7_COLS].astype(str)
            same = new_legacy.equals(old_legacy)
            print(f"[byte-identity guard] columns 1-7 identical to the "
                  f"pre-existing pinned file: {same}")
            if not same:
                diff_mask = new_legacy != old_legacy
                bad_cols = [c for c in LEGACY_7_COLS if diff_mask[c].any()]
                sys.exit(f"FATAL: columns 1-7 changed by this rebuild in "
                          f"{bad_cols} -- refusing to overwrite the pinned "
                          f"luad_batch_groups.csv. The frozen-rule draft "
                          f"(D2 P0, section 0) requires columns 1-7 to be "
                          f"byte-identical across this 8th-column-append "
                          f"rebuild; investigate before rerunning.")

    # --- write outputs -----------------------------------------------------
    # LF only, explicitly. Opening the file handle with newline="\n" is NOT
    # enough on its own: pandas' csv writer composes lines with an explicit
    # "\r\n" terminator by default (the classic "to_csv produces CRLF on
    # Windows" gotcha), which newline= translation does not touch, since it
    # only rewrites a bare "\n" the writer emits, not a "\r\n" it already
    # wrote. lineterminator="\n" makes the writer itself emit LF only, and
    # newline="" on the handle stops Python's text layer from adding a
    # second, redundant translation on top of that (the standard combination
    # per the csv/pandas docs). A CRLF file would break the HPC's
    # `sha256sum` comparison in lsf_precheck_luad.sh (bash reads LF; the two
    # hash differently) -- this exact bug was caught by review before it
    # ever reached signature.
    os.makedirs(PINNED, exist_ok=True)
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
        df.to_csv(f, index=False, lineterminator="\n")
    with open(OUT_CSV, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    with open(OUT_SHA, "w", encoding="utf-8", newline="\n") as f:
        f.write(digest + "\n")
    print(f"\nwrote {OUT_CSV} ({len(df)} rows)")
    print(f"wrote {OUT_SHA}")
    print(f"sha256: {digest}")


if __name__ == "__main__":
    main()
