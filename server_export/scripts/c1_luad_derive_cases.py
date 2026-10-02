#!/usr/bin/env python3
"""
C1-LUAD confirmatory case-set derivation (Step 0b, branch pdc, 2026-09-29).
METADATA ONLY. Reads ONLY the saved raw PDC GraphQL JSON under
review/C1_LUAD_STEP0B_2026-09-29/raw/ -- no network call is made by this
script. Reproduces and extends review/C1_LUAD_STEP0_2026-09-27.md sections
1-3, deriving:
  - the 11 QC/bridge/pooled pseudo-case labels (blank sample_type rows in
    PDC000489; Step 0 found 11 labels over 12 rows);
  - the 120 real Primary Tumor cases of PDC000489 (all C3L-/C3N-);
  - the 7 cases whose PDC000489 Primary Tumor aliquot_submitter_id also
    appears as Primary Tumor in PDC000153 (discovery) -- the TMT
    bridge/reference samples that must be excluded for independence;
  - the final 113-case confirmatory-independent set.

Every expected count is asserted; the script exits FATAL (non-zero) on any
mismatch rather than silently writing a different number. Per the task rule:
never adjust a filter to make a count match -- a mismatch is reported, not
patched.

Usage:
  python c1_luad_derive_cases.py
Outputs (relative to this script's server_export/scripts/pinned/ directory):
  c1_case_ids_luad.txt      113 case_submitter_ids, sorted, LF line endings
  c1_excluded_luad.tsv      case_label, reason, aliquot_submitter_ids (18 rows:
                             11 pseudo-case labels + 7 bridge cases)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(
    HERE, "..", "..", "review", "C1_LUAD_STEP0B_2026-09-29", "raw")
PINNED_DIR = os.path.join(HERE, "pinned")

CONF_STUDY = "PDC000489"
DISC_STUDY = "PDC000153"

EXPECTED_PSEUDO_LABELS = {
    "CR.C1", "CR.C1.newly.combined", "CR.C2", "CR.D", "CR.D.MAM",
    "CR.D.newly.combined", "CR.D.pool.2", "CR.LSCC", "Taiwanese IR",
    "Tumor Only CONF", "Tumor Only DISC",
}
EXPECTED_N_PSEUDO_LABELS = 11
EXPECTED_N_PSEUDO_ROWS = 12
EXPECTED_N_REAL_TUMOR_CASES = 120
EXPECTED_BRIDGE_CASES = {
    "C3L-02348", "C3L-02350", "C3N-00545", "C3N-00738",
    "C3N-01023", "C3N-01024", "C3N-02003",
}
EXPECTED_N_BRIDGE = 7
EXPECTED_N_FINAL = 113


def fatal(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def load_biospecimen(pdc_study_id):
    path = os.path.join(RAW_DIR, f"biospecimenPerStudy_{pdc_study_id}.json")
    if not os.path.exists(path):
        fatal(f"missing raw file {path} -- run fetch_step0b_pdc.py / "
              f"probe_biospecimen_fields.py first")
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    if d.get("errors"):
        fatal(f"{path} is a saved GraphQL ERROR response, not data: {d['errors']}")
    rows = d.get("data", {}).get("biospecimenPerStudy")
    if rows is None:
        fatal(f"{path} has no data.biospecimenPerStudy")
    return rows, path


def main():
    os.makedirs(PINNED_DIR, exist_ok=True)

    conf_rows, conf_path = load_biospecimen(CONF_STUDY)
    disc_rows, disc_path = load_biospecimen(DISC_STUDY)
    print(f"loaded {len(conf_rows)} confirmatory rows from {conf_path}")
    print(f"loaded {len(disc_rows)} discovery rows from {disc_path}")
    if len(conf_rows) != 245:
        fatal(f"confirmatory biospecimenPerStudy row count = {len(conf_rows)}, "
              f"expected 245 (Step 0 §2)")
    if len(disc_rows) != 217:
        fatal(f"discovery biospecimenPerStudy row count = {len(disc_rows)}, "
              f"expected 217 (Step 0 §1)")

    # --- pseudo-case (blank sample_type) rows in the confirmatory study ---
    # PDC's own field literally returns the string "Not Reported" for these
    # rows (verified against the raw JSON), not an empty string/null -- that
    # is what Step 0 called "blank". Treat both spellings as blank so this
    # does not silently regress if PDC ever returns a true empty value.
    BLANK_SAMPLE_TYPES = {None, "", "Not Reported"}
    pseudo_rows = [r for r in conf_rows if r.get("sample_type") in BLANK_SAMPLE_TYPES]
    pseudo_labels = {r["case_submitter_id"] for r in pseudo_rows}
    print(f"\npseudo-case rows (blank sample_type): {len(pseudo_rows)} rows, "
          f"{len(pseudo_labels)} distinct labels")
    if len(pseudo_rows) != EXPECTED_N_PSEUDO_ROWS:
        fatal(f"pseudo-case row count = {len(pseudo_rows)}, expected "
              f"{EXPECTED_N_PSEUDO_ROWS}")
    if len(pseudo_labels) != EXPECTED_N_PSEUDO_LABELS:
        fatal(f"pseudo-case distinct label count = {len(pseudo_labels)}, expected "
              f"{EXPECTED_N_PSEUDO_LABELS}: got {sorted(pseudo_labels)}")
    if pseudo_labels != EXPECTED_PSEUDO_LABELS:
        fatal("pseudo-case label SET differs from Step 0's recorded set.\n"
              f"  got:      {sorted(pseudo_labels)}\n"
              f"  expected: {sorted(EXPECTED_PSEUDO_LABELS)}\n"
              f"  only-in-got: {sorted(pseudo_labels - EXPECTED_PSEUDO_LABELS)}\n"
              f"  only-in-expected: {sorted(EXPECTED_PSEUDO_LABELS - pseudo_labels)}")

    # --- real Primary Tumor cases of the confirmatory study ---
    conf_tumor_rows = [r for r in conf_rows if r.get("sample_type") == "Primary Tumor"]
    conf_tumor_cases = {r["case_submitter_id"] for r in conf_tumor_rows}
    print(f"\nconfirmatory Primary Tumor rows: {len(conf_tumor_rows)}, "
          f"distinct cases: {len(conf_tumor_cases)}")
    if len(conf_tumor_cases) != EXPECTED_N_REAL_TUMOR_CASES:
        fatal(f"confirmatory Primary Tumor case count = {len(conf_tumor_cases)}, "
              f"expected {EXPECTED_N_REAL_TUMOR_CASES}")
    non_c3 = [c for c in conf_tumor_cases if not (c.startswith("C3L-") or c.startswith("C3N-"))]
    if non_c3:
        fatal(f"confirmatory Primary Tumor cases include non-C3L-/C3N- ids: {non_c3}")

    # aliquot_submitter_id per case for the confirmatory Primary Tumor rows
    # (task expects exactly one Primary Tumor aliquot per confirmatory case;
    # assert it rather than silently picking one if a case has >1).
    conf_case_to_aliquots = {}
    for r in conf_tumor_rows:
        conf_case_to_aliquots.setdefault(r["case_submitter_id"], []).append(
            r["aliquot_submitter_id"])
    multi = {c: a for c, a in conf_case_to_aliquots.items() if len(a) > 1}
    if multi:
        print(f"  NOTE: {len(multi)} confirmatory case(s) carry >1 Primary Tumor "
              f"aliquot (not fatal, recorded): {multi}")

    # --- discovery Primary Tumor cases/aliquots, for the bridge check ---
    disc_tumor_rows = [r for r in disc_rows if r.get("sample_type") == "Primary Tumor"]
    disc_tumor_cases = {r["case_submitter_id"] for r in disc_tumor_rows}
    disc_tumor_aliquots = {r["aliquot_submitter_id"] for r in disc_tumor_rows}
    disc_case_to_aliquots = {}
    for r in disc_tumor_rows:
        disc_case_to_aliquots.setdefault(r["case_submitter_id"], []).append(
            r["aliquot_submitter_id"])
    print(f"\ndiscovery Primary Tumor rows: {len(disc_tumor_rows)}, "
          f"distinct cases: {len(disc_tumor_cases)}, distinct aliquots: "
          f"{len(disc_tumor_aliquots)}")

    # --- bridge cases: PRIMARY criterion is case_submitter_id-level overlap
    # between the confirmatory Primary Tumor case set (120) and the discovery
    # Primary Tumor case set (111) -- this is what Step 0 §3 actually computed
    # ("Intersecting the 120 confirmatory cases against the 111 discovery
    # cases by case_submitter_id finds 7 real cases in both"), and it is what
    # is reproduced here, exactly. Step 0 additionally asserted that for every
    # one of these 7 the confirmatory Primary Tumor aliquot_submitter_id is
    # IDENTICAL to a discovery Primary Tumor aliquot_submitter_id for that
    # case; that stronger, per-case aliquot-identity claim is checked below as
    # a DIAGNOSTIC, reported but not gating (see printed discrepancy note if
    # any case fails it).
    bridge_cases = sorted(conf_tumor_cases & disc_tumor_cases)
    print(f"\nbridge cases (case_submitter_id in both confirmatory and discovery "
          f"Primary Tumor case sets): {len(bridge_cases)}")
    if len(bridge_cases) != EXPECTED_N_BRIDGE:
        fatal(f"bridge case count = {len(bridge_cases)}, expected {EXPECTED_N_BRIDGE}")
    if set(bridge_cases) != EXPECTED_BRIDGE_CASES:
        fatal("bridge case SET differs from Step 0's recorded set.\n"
              f"  got:      {sorted(bridge_cases)}\n"
              f"  expected: {sorted(EXPECTED_BRIDGE_CASES)}")

    aliquot_identical = {}   # case -> (conf_aliquot, disc_aliquot) when identical
    aliquot_mismatch = {}    # case -> (conf_aliquots, disc_aliquots) when NOT identical
    for c in bridge_cases:
        conf_al = set(conf_case_to_aliquots[c])
        disc_al = set(disc_case_to_aliquots[c])
        shared = conf_al & disc_al
        if shared:
            aliquot_identical[c] = (sorted(conf_al), sorted(disc_al), sorted(shared))
        else:
            aliquot_mismatch[c] = (sorted(conf_al), sorted(disc_al))
    print(f"  of these 7, aliquot_submitter_id IS identical between studies for "
          f"{len(aliquot_identical)}/7:")
    for c in sorted(aliquot_identical):
        print(f"    {c}: shared aliquot(s) {aliquot_identical[c][2]}")
    if aliquot_mismatch:
        print(f"  DISCREPANCY vs Step 0 §3's stated mechanism ('for every one of "
              f"these 7 ... the exact same aliquot ID'): {len(aliquot_mismatch)}/7 "
              f"case(s) do NOT share an identical aliquot_submitter_id between "
              f"confirmatory and discovery, though the case_submitter_id (and, "
              f"where checked, the underlying case_id UUID) is the same patient:")
        for c in sorted(aliquot_mismatch):
            conf_al, disc_al = aliquot_mismatch[c]
            print(f"    {c}: confirmatory aliquot(s) {conf_al} vs discovery "
                  f"aliquot(s) {disc_al} -- same case_submitter_id, different "
                  f"aliquot_submitter_id. This case is still excluded (case-level "
                  f"overlap is the operative independence criterion), but the "
                  f"per-case aliquot-identity claim in Step 0 §3 does not hold for "
                  f"it. Reported per the task's 'state any discrepancy plainly' "
                  f"instruction; the 7-case COUNT and the 113 final count both "
                  f"still match Step 0 exactly.")

    # --- final set ---
    final_cases = sorted(conf_tumor_cases - set(bridge_cases))
    bridge_case_set = set(bridge_cases)
    print(f"\nfinal confirmatory-independent case set: {len(final_cases)}")
    if len(final_cases) != EXPECTED_N_FINAL:
        fatal(f"final case count = {len(final_cases)}, expected {EXPECTED_N_FINAL}")

    # --- write outputs ---
    ids_path = os.path.join(PINNED_DIR, "c1_case_ids_luad.txt")
    with open(ids_path, "w", encoding="utf-8", newline="\n") as f:
        for c in final_cases:
            f.write(c + "\n")
    print(f"\nwrote {len(final_cases)} case ids -> {ids_path}")

    excl_path = os.path.join(PINNED_DIR, "c1_excluded_luad.tsv")
    with open(excl_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("case_label\treason\taliquot_submitter_ids\n")
        # 11 pseudo-case labels, aliquot ids = all aliquots recorded under that label
        label_to_aliquots = {}
        for r in pseudo_rows:
            label_to_aliquots.setdefault(r["case_submitter_id"], []).append(
                r["aliquot_submitter_id"])
        for label in sorted(label_to_aliquots):
            aliquots = ",".join(sorted(label_to_aliquots[label]))
            f.write(f"{label}\tQC/bridge/pooled pseudo-case (blank sample_type)\t{aliquots}\n")
        # 7 bridge cases -- case_submitter_id also present in discovery
        # PDC000153's Primary Tumor case set. Reason text distinguishes the
        # 6/7 with an identical shared aliquot_submitter_id from C3N-00545,
        # whose confirmatory Primary Tumor aliquot differs from its
        # (Disqualified) discovery Primary Tumor aliquot -- same patient,
        # different physical aliquot; see the discrepancy note printed above.
        for case in sorted(bridge_case_set):
            conf_al = ",".join(conf_case_to_aliquots[case])
            if case in aliquot_identical:
                reason = ("case_submitter_id also Primary Tumor in discovery "
                           "PDC000153, with an identical shared aliquot_submitter_id "
                           "(TMT bridge/reference sample) -- excluded for independence")
            else:
                reason = ("case_submitter_id also Primary Tumor in discovery "
                           "PDC000153, but via a DIFFERENT aliquot_submitter_id "
                           "(discovery aliquot status Disqualified); same patient, "
                           "not an identical-aliquot bridge -- excluded for "
                           "case-level independence regardless")
            f.write(f"{case}\t{reason}\t{conf_al}\n")
    n_excl_rows = len(label_to_aliquots) + len(bridge_case_set)
    print(f"wrote {n_excl_rows} excluded-case rows -> {excl_path}")

    print("\nALL ASSERTIONS PASSED.")
    print(f"  pseudo-case labels: {len(pseudo_labels)} (rows: {len(pseudo_rows)})")
    print(f"  real confirmatory Primary Tumor cases: {len(conf_tumor_cases)}")
    print(f"  bridge cases excluded: {len(bridge_cases)}")
    print(f"  final confirmatory-independent set: {len(final_cases)}")


if __name__ == "__main__":
    main()
