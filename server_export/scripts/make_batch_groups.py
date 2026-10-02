#!/usr/bin/env python3
"""
Generic, deterministic operator-balanced sel/test group builder.

Written 2026-09-29 for C1-LUAD Step 0b (metadata-only; see
review/C1_LUAD_STEP0B_2026-09-29/). Its trigger: no script anywhere in this
project, and none in any local Claude scratchpad tree under
C:/Users/wei/AppData/Local/Temp/claude/, builds `ucec_batch_groups.csv`'s
`op_balanced_0..4` columns -- only the CSV output survives, plus its consumer
(`split_replication_power.py --group-csv/--group-col`, which only READS a
group CSV, see its L129-137, L155-165). This script infers the schema from
`ucec_batch_groups.csv` (case_submitter_id + N `op_balanced_<seed>` columns,
each valued 'sel'/'test', with every operator kept whole within a column --
verified: no case_submitter_id's operator ever splits across sel/test in that
file) and reimplements a generic builder from scratch.

IMPORTANT: this is NOT a recovered original. There is no way to verify
byte-for-byte reproduction of `ucec_batch_groups.csv` without the original
code, which does not exist locally. What is reproduced is the documented
*schema and properties*: each operator kept whole; roughly balanced sel/test
sizes; the floors `split_replication_power.py` itself enforces at L181
(n_sel >= 40, n_test >= 30). Do not use this script's op_balanced_* output for
UCEC as a replacement for the existing `ucec_batch_groups.csv` -- that file is
untouched and remains the one already used by any prior UCEC run.

Per-case operator label: the modal `operator_uuid` over that case's
`sample_type == "Primary Tumor"` slides (per the cohort's
`slide_type_map[_<cohort>].tsv`) in `slide_acquisition_meta.tsv`, mirroring
`batch_leak_check.case_batch_labels`'s "operator" axis
(server_export/scripts/batch_leak_check.py L120-142) -- EXCEPT the "slide
actually entered the case's mean-pooled embedding" filter there
(`tumor & present`, using an `emb/*.pt` file listing on the HPC filesystem at
ROOT=/public/home/fjhui/ZW, L104-118, L127-128), which this local machine
cannot see. This script uses ALL Primary Tumor slides for the case instead of
that HPC-verified subset. Every caller of this module must state this
limitation; it is not silently assumed away. batch_leak_check.py itself is
untouched -- this is an independent reimplementation, not an import of it
(importing it would require the HPC env vars / residual_analysis module it
sets up in load_wsi(), which are not meaningful here).

CLI usage: python make_batch_groups.py <cohort> <output_csv_path>
  Writes case_submitter_id + op_balanced_0..4 for every case in
  figures/figdata/scores_<cohort>.csv (fallback: every operator-labelled case
  found in the slide map) to <output_csv_path>.

Importable API: load_case_operator_labels, load_case_quarter_labels,
balance_operators, build, stream_side.

Added 2026-09-29 (frozen-rule draft v3, D2 P0 / D6d):
`stream_side` -- the LUAD-only, descriptive `stream_c3n_vs_c3l` column. Unlike
every function above, it reads no metadata file at all: it is a pure function
of the case_submitter_id's own accrual-stream prefix. Every case in the
105-case discovery analysis set used here is C3N-nnnnn or C3L-nnnnn; a third
prefix, 11LU- (4 Taiwanese cases), exists elsewhere in the LUAD cohort but is
discovery-only per the frozen rule (section 6) and never reaches this
function's input. An 11LU- (or any other unrecognised) prefix comes back
blank, not silently misclassified. This column is descriptive only (D2 P2);
it never gates GO/MARGINAL/NO-GO and is not one of the seven columns that
mirror `ucec_batch_groups.csv`'s schema.
"""
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_META = os.path.join(HERE, "slide_acquisition_meta.tsv")
FIGDATA = os.path.normpath(os.path.join(HERE, "..", "..", "figures", "figdata"))

SEEDS = [0, 1, 2, 3, 4]
N_SEL_FLOOR = 40   # split_replication_power.py L181 ("n_sel < 40 or n_test < 30")
N_TEST_FLOOR = 30


def slide_map_path(cohort):
    # Mirrors batch_leak_check.cohort_paths's suffix rule (L60), applied to
    # the *unpinned* server_export/scripts/ copy, since that is what the
    # Step 0b task named explicitly for LUAD.
    suffix = "" if cohort == "ccrcc" else f"_{cohort}"
    return os.path.join(HERE, f"slide_type_map{suffix}.tsv")


def scores_path(cohort):
    return os.path.join(FIGDATA, f"scores_{cohort}.csv")


def _tumor_slide_meta(cohort, meta_path=None, slide_map_path_=None):
    meta_path = meta_path or DEFAULT_META
    slide_map_path_ = slide_map_path_ or slide_map_path(cohort)
    meta = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("")
    meta = meta[meta["cohort"] == cohort]
    smap = pd.read_csv(slide_map_path_, sep="\t", dtype=str)
    tumor_slides = set(smap.loc[smap["sample_type"] == "Primary Tumor", "slide_submitter_id"])
    return meta[meta["slide_submitter_id"].isin(tumor_slides)]


def _modal_label(meta_t, column):
    label, no_label = {}, []
    for case, grp in meta_t.groupby("case_submitter_id"):
        vals = [v for v in grp[column].tolist() if v]
        if not vals:
            no_label.append(case)
            continue
        label[case] = Counter(vals).most_common()[0][0]
    return label, sorted(no_label)


def load_case_operator_labels(cohort, meta_path=None, slide_map_path_=None):
    """Per-case modal operator_uuid. Returns (label_dict, no_label_sorted_list)."""
    meta_t = _tumor_slide_meta(cohort, meta_path, slide_map_path_)
    return _modal_label(meta_t, "operator_uuid")


def load_case_quarter_labels(cohort, meta_path=None, slide_map_path_=None):
    """Per-case modal scan_quarter. Returns (label_dict, no_label_sorted_list)."""
    meta_t = _tumor_slide_meta(cohort, meta_path, slide_map_path_)
    return _modal_label(meta_t, "scan_quarter")


def stream_side(case_submitter_id):
    """C3N/C3L accrual-stream split (D2 P0, D6d, descriptive only).

    sel = C3N-prefixed cases, test = C3L-prefixed cases. Pure string check on
    the case_submitter_id itself -- no slide_acquisition_meta.tsv lookup, no
    slide map, no operator/quarter label, so it cannot be blank for any case
    that has ever been GDC/PDC-catalogued under the CPTAC3 LUAD naming scheme
    (every one of the 105 discovery analysis cases matches one of the two
    prefixes; see build_luad_batch_groups.py's verification print). A
    case_submitter_id matching neither prefix returns "" (blank), the same
    convention every other column here uses for "no label".
    """
    if case_submitter_id.startswith("C3N-"):
        return "sel"
    if case_submitter_id.startswith("C3L-"):
        return "test"
    return ""


def balance_operators(op_label, cases, seed):
    """Greedy longest-operator-first bin balancing; operator kept whole.

    Orders operators by group size descending, ties broken by a per-operator
    RandomState(seed) draw. Assigns each whole operator group to whichever
    side (test/sel) is currently further below its ~50/50 target of the
    labelled case count. Cases in `cases` with no entry in `op_label` are
    absent from the returned dict (no label -> no side).
    """
    labelled = [c for c in cases if c in op_label]
    groups = {}
    for c in labelled:
        groups.setdefault(op_label[c], []).append(c)

    n = len(labelled)
    target_test = round(n / 2.0)
    target_sel = n - target_test

    op_names = sorted(groups)  # fixed order before any randomisation
    rng = np.random.RandomState(seed)
    tiebreak = {op: rng.random() for op in op_names}
    order = sorted(op_names, key=lambda op: (-len(groups[op]), tiebreak[op]))

    side = {}
    n_test = n_sel = 0
    for op in order:
        size = len(groups[op])
        deficit_test = target_test - n_test
        deficit_sel = target_sel - n_sel
        if deficit_test >= deficit_sel:
            side[op] = "test"
            n_test += size
        else:
            side[op] = "sel"
            n_sel += size

    return {c: side[op_label[c]] for c in labelled}


def build(cohort, seeds=SEEDS, cases=None):
    """Returns (DataFrame[case_submitter_id, op_balanced_<seed>...], info_dict)."""
    op_label, op_no_label = load_case_operator_labels(cohort)
    if cases is None:
        scores_f = scores_path(cohort)
        if os.path.exists(scores_f):
            cases = sorted(pd.read_csv(scores_f, index_col=0).index.astype(str))
            case_source = scores_f
        else:
            cases = sorted(op_label)
            case_source = ("no scores_<cohort>.csv found; fell back to every "
                            "operator-labelled case in the slide map")
    else:
        case_source = "caller-supplied case list"

    columns = {}
    for s in seeds:
        assign = balance_operators(op_label, cases, s)
        columns[f"op_balanced_{s}"] = [assign.get(c, "") for c in cases]

    df = pd.DataFrame({"case_submitter_id": cases, **columns})
    info = dict(case_source=case_source, n_cases=len(cases),
                op_no_label=[c for c in cases if c not in op_label],
                op_label=op_label)
    return df, info


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: python make_batch_groups.py <cohort> <output_csv_path>")
    cohort, out_path = sys.argv[1], sys.argv[2]
    df, info = build(cohort)
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"cohort={cohort} n_cases={info['n_cases']} case_source={info['case_source']}")
    print(f"cases with no operator label ({len(info['op_no_label'])}): {info['op_no_label']}")
    for col in df.columns:
        if col == "case_submitter_id":
            continue
        vc = df[col].value_counts()
        n_sel, n_test = int(vc.get("sel", 0)), int(vc.get("test", 0))
        ok = n_sel >= N_SEL_FLOOR and n_test >= N_TEST_FLOOR
        print(f"  {col}: n_sel={n_sel} n_test={n_test} floor_ok={ok}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
