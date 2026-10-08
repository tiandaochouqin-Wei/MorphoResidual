#!/usr/bin/env python3
"""ADAPTED from review/BATCH_AXES_SCOPING_2026-10-02/build_plex_maps.py (not edited there).
Case -> TMT plex map for the 138 C1-UCEC confirmatory cases (PDC000439). Labels and
metadata only: no protein, RNA or embedding value is read, no statistic is computed.

Route 1 (per-channel design): raw/r2_sed_full_PDC000439.json (studyExperimentalDesign,
  tmt_126..tmt_131c -> aliquot_submitter_id) + raw/r3_biospecimen_PDC000439.json
  (biospecimenPerStudy: aliquot -> case, sample_type, aliquot_status).
Route 2 (independent): raw/r5_pcsa_arm_PDC000439.json via route2_cases_PDC000439.json
  (paginatedCasesSamplesAliquots: case -> sample -> aliquot -> aliquot_run_metadata{id,
  label, experiment_number}). Its case->aliquot assignment and sample_type do NOT come
  from biospecimenPerStudy, and a (run, channel) is accepted only if the SAME
  aliquot_run_metadata_id carries the same aliquot, label TMT_<channel> and the design
  row's experiment_number.

Resolution rule for a case with more than one measured Primary Tumor aliquot or an
aliquot measured in more than one plex (same rule as the C1-LUAD / UCEC-discovery plex
maps, POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md section P): the plex of the
non-Disqualified aliquot. If that still leaves several plexes: majority by number of
measured tumour channels, ties -> lowest plex (build_plex_maps.py's rule); such a case is
flagged `ambiguous` and listed, so the P3 sensitivity run can drop it.

Cases for which the PDC "Primary Tumor" definition (build_plex_maps.py) and the
definition that c1_pull_omics.py used to build the protein matrix (every aliquot whose
sample_type does not contain "Normal", first row per aliquot) differ are reported.

Outputs (this folder): channel_map_ucec_c1.tsv, plex_map_ucec_c1.tsv,
case_to_plex_ucec_c1.tsv (2 columns case_id, plex; the file to upload),
plex_summary_ucec_c1.json.
"""
import json
import os
from collections import Counter, defaultdict

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
RAW = os.path.join(HERE, "raw")
PID = "PDC000439"
CASE_IDS = os.path.join(PAPER, "server_export", "scripts", "pinned", "c1_case_ids_ucec.txt")
TMT = ["tmt_126", "tmt_127n", "tmt_127c", "tmt_128n", "tmt_128c", "tmt_129n",
       "tmt_129c", "tmt_130n", "tmt_130c", "tmt_131", "tmt_131c"]
MAX_LEVELS, MIN_GROUP = 12, 3   # residual_analysis_sitepack.py / c1_luad_plex_reread MIN_GROUP


def main():
    design = json.load(open(os.path.join(RAW, f"r2_sed_full_{PID}.json"), encoding="utf-8"))
    design = design["data"]["studyExperimentalDesign"]
    bio = json.load(open(os.path.join(RAW, f"r3_biospecimen_{PID}.json"), encoding="utf-8"))
    bio = bio["data"]["biospecimenPerStudy"]
    r2 = json.load(open(os.path.join(HERE, f"route2_cases_{PID}.json"), encoding="utf-8"))
    cases = sorted(l.strip() for l in open(CASE_IDS, encoding="utf-8") if l.strip())
    assert len(cases) == 138 and len(set(cases)) == 138

    # --- biospecimen index (route 1) ---
    a2 = defaultdict(set)
    for r in bio:
        a2[r["aliquot_submitter_id"]].add((r["case_submitter_id"], (r["sample_type"] or "").strip(),
                                           r.get("aliquot_status") or ""))
    tumor_aliquots_of = defaultdict(set)
    for a, s in a2.items():
        for case, st, _ in s:
            if st == "Primary Tumor":
                tumor_aliquots_of[case].add(a)
    status_of = {a: sorted({s[2] for s in v}) for a, v in a2.items()}

    # c1_pull_omics.build_aliquot_xwalk rule: sample_type without "Normal", first row wins
    pull_xwalk = {}
    for r in bio:
        a = r.get("aliquot_submitter_id")
        st = (r.get("sample_type") or "").strip()
        if not a or "Normal" in st or a in pull_xwalk:
            continue
        pull_xwalk[a] = (r.get("case_submitter_id"), st)
    pull_aliquots_of = defaultdict(set)
    for a, (c, _st) in pull_xwalk.items():
        pull_aliquots_of[c].add(a)

    # --- channel map (route 1) ---
    rows = []
    id2run = {}   # aliquot_run_metadata_id -> (plex, channel, aliquot, experiment_number)
    for run in design:
        plex = run["plex_dataset_name"] or run["study_run_metadata_submitter_id"]
        for ch in TMT:
            for e in (run.get(ch) or []):
                a = e["aliquot_submitter_id"]
                info = sorted(a2.get(a, []))
                cs = sorted({i[0] for i in info})
                st = sorted({i[1] for i in info})
                stat = sorted({i[2] for i in info})
                if not info:
                    role = "reference/QC (not a biospecimen of this study)"
                elif cs and all(c.startswith(("C3L-", "C3N-")) for c in cs):
                    role = "patient:" + "/".join(st)
                else:
                    role = "pseudo-case (" + "/".join(st) + ")"
                id2run[e["aliquot_run_metadata_id"]] = (plex, ch.replace("tmt_", ""), a,
                                                        run.get("experiment_number"))
                rows.append(dict(plex=plex, plex_idx=plex[:2], experiment_number=run.get("experiment_number"),
                                 channel=ch.replace("tmt_", ""), aliquot_submitter_id=a,
                                 aliquot_run_metadata_id=e["aliquot_run_metadata_id"],
                                 case_submitter_id=";".join(cs), sample_type=";".join(st),
                                 aliquot_status=";".join(stat), role=role))
    ch_df = pd.DataFrame(rows)
    ch_df.to_csv(os.path.join(HERE, "channel_map_ucec_c1.tsv"), sep="\t", index=False)
    aliq_runs = ch_df.groupby("aliquot_submitter_id")["plex"].apply(lambda s: sorted(set(s))).to_dict()
    nonpat = ch_df[~ch_df.role.str.startswith("patient")]
    ref_by_channel = {ch: dict(Counter(g["aliquot_submitter_id"]).most_common())
                      for ch, g in nonpat.groupby("channel")}

    # --- route 2: case -> aliquot -> run, accepted only on the full (id, aliquot, label, exp) match ---
    r2_hit, r2_miss = [], Counter()
    r2_aliq_case = {}
    for r in r2:
        hit = id2run.get(r["aliquot_run_metadata_id"])
        if hit is None:
            continue   # run metadata belonging to other studies/fractions
        plex, ch, a, expn = hit
        ok = (a == r["aliquot"] and r["label"].lower() == f"tmt_{ch}" and
              str(r["experiment_number"]) == str(expn))
        if not ok:
            r2_miss["id_matches_but_aliquot_label_or_experiment_differs"] += 1
            continue
        r2_hit.append(dict(r, plex=plex))
        r2_aliq_case[a] = (r["case"], r["sample_type"])
    r2_by_case = defaultdict(lambda: defaultdict(set))   # case -> aliquot -> plexes
    for r in r2_hit:
        if r["sample_type"] == "Primary Tumor":
            r2_by_case[r["case"]][r["aliquot"]].add(r["plex"])

    out = []
    for c in cases:
        tal = sorted(tumor_aliquots_of.get(c, []))
        measured = [(a, r) for a in tal for r in aliq_runs.get(a, [])]
        plexes = sorted({r for _, r in measured})
        disq = sorted(a for a in tal if "Disqualified" in status_of.get(a, []))
        meas_aliq = sorted({a for a, _ in measured})
        qual_measured = [(a, r) for a, r in measured if a not in disq]
        qual_plexes = sorted({r for _, r in qual_measured})
        resolution = ""
        if len(plexes) == 1:
            lab, status = plexes[0], "single"
        elif len(plexes) == 0:
            lab, status = "", "no_tumor_channel"
        else:
            if len(qual_plexes) == 1:
                lab, status = qual_plexes[0], "multi_plex"
                resolution = "plex of the non-Disqualified aliquot"
            else:
                pool = qual_measured if qual_measured else measured
                cnt = Counter(r for _, r in pool)
                top = max(cnt.values())
                lab, status = sorted(r for r, k in cnt.items() if k == top)[0], "ambiguous"
                resolution = "majority of measured tumour channels, ties -> lowest plex"
        # build_plex_maps.py's original rule (for comparison only)
        if len(plexes) > 1:
            cnt0 = Counter(r for _, r in measured)
            top0 = max(cnt0.values())
            old = sorted(r for r, k in cnt0.items() if k == top0)[0]
        else:
            old = plexes[0] if plexes else ""
        # route 2
        r2_plexes = sorted({p for a, ps in r2_by_case.get(c, {}).items() for p in ps})
        r2_aliq = sorted(r2_by_case.get(c, {}))
        r2_qual = sorted({p for a, ps in r2_by_case.get(c, {}).items()
                          if a not in disq for p in ps})
        r2_lab = ""
        if len(r2_plexes) == 1:
            r2_lab = r2_plexes[0]
        elif len(r2_qual) == 1:
            r2_lab = r2_qual[0]
        elif r2_plexes:
            cnt2 = Counter(p for a, ps in r2_by_case[c].items() if a not in disq for p in ps) \
                or Counter(p for a, ps in r2_by_case[c].items() for p in ps)
            top2 = max(cnt2.values())
            r2_lab = sorted(p for p, k in cnt2.items() if k == top2)[0]
        pull_al = sorted(pull_aliquots_of.get(c, []))
        out.append(dict(case_id=c, plex=lab, plex_idx=lab[:2] if lab else "", status=status,
                        resolution=resolution,
                        n_tumor_aliquots_pdc=len(tal), tumor_aliquots_pdc=";".join(tal),
                        disqualified_aliquots=";".join(disq),
                        measured_aliquots=";".join(meas_aliq),
                        n_tumor_channels=len(measured), all_plexes=";".join(plexes),
                        plex_original_build_rule=old,
                        differs_from_original_rule=(old != lab),
                        route2_aliquots=";".join(r2_aliq), route2_plexes=";".join(r2_plexes),
                        route2_plex=r2_lab, route2_agrees=(r2_lab == lab),
                        pull_rule_aliquots=";".join(pull_al),
                        pull_rule_differs=(set(pull_al) != set(tal))))
    pm = pd.DataFrame(out)
    pm.to_csv(os.path.join(HERE, "plex_map_ucec_c1.tsv"), sep="\t", index=False)

    labelled = pm[pm.plex != ""]
    unl = pm[pm.plex == ""].case_id.tolist()
    with open(os.path.join(HERE, "case_to_plex_ucec_c1.tsv"), "w", encoding="utf-8", newline="\n") as f:
        f.write("case_id\tplex\n")
        for _, r in pm.sort_values("case_id").iterrows():
            if r.plex:
                f.write(f"{r.case_id}\t{r.plex}\n")

    counts = Counter(labelled.plex)
    sizes = dict(sorted(counts.items()))
    singletons = sorted(p for p, k in counts.items() if k == 1)
    n_fixed = len(unl) + len(singletons)
    dummy = [lv for lv, k in counts.most_common() if k >= MIN_GROUP]
    left = sorted(lv for lv in counts if lv not in set(dummy))
    dropped = None
    if dummy and len(labelled) == len(pm) and not left:
        dropped = dummy[-1]
        dummy = dummy[:-1]
    cases_in_ge3 = sum(k for k in counts.values() if k >= MIN_GROUP)
    summary = dict(
        pdc_study_id=PID, n_runs=len(design), experiment_type=sorted({r["experiment_type"] for r in design}),
        n_channel_rows=len(ch_df), n_patient_channels=int(ch_df.role.str.startswith("patient").sum()),
        non_patient_channels_by_channel=ref_by_channel,
        analysed_cases=len(pm), labelled=len(labelled), unlabelled=unl,
        status_counts=dict(Counter(pm.status)),
        n_plexes=len(counts), cases_per_plex=sizes,
        sizes_sorted=sorted(counts.values(), reverse=True),
        n_strata_ge3=sum(1 for k in counts.values() if k >= 3),
        n_singleton_strata=len(singletons), singleton_strata=singletons,
        cases_held_fixed=n_fixed, frac_cases_fixed=n_fixed / len(pm),
        cases_in_strata_ge3=cases_in_ge3,
        plex_dummies_if_residualised=len(dummy), dummy_dropped_for_full_rank=dropped,
        levels_left_to_intercept=left,
        cases_with_multiple_tumor_aliquots=pm[pm.n_tumor_aliquots_pdc > 1].case_id.tolist(),
        cases_with_disqualified_aliquot=pm[pm.disqualified_aliquots != ""].case_id.tolist(),
        cases_multi_plex=pm[pm.status == "multi_plex"].case_id.tolist(),
        cases_ambiguous=pm[pm.status == "ambiguous"].case_id.tolist(),
        cases_differing_from_original_build_rule=pm[pm.differs_from_original_rule].case_id.tolist(),
        route2_agree=int(pm.route2_agrees.sum()), route2_disagree=pm[~pm.route2_agrees].case_id.tolist(),
        route2_channel_matches=len(r2_hit), route2_channel_mismatch=dict(r2_miss),
        design_channels_total=len(id2run),
        design_channels_found_in_route2=len({r["aliquot_run_metadata_id"] for r in r2_hit}),
        pull_rule_differs=pm[pm.pull_rule_differs].case_id.tolist(),
        patient_aliquots_in_more_than_one_run={
            a: r for a, r in aliq_runs.items() if len(r) > 1 and a in tumor_aliquots_of_all(tumor_aliquots_of)},
    )
    with open(os.path.join(HERE, "plex_summary_ucec_c1.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1, default=str)
    for k, v in summary.items():
        print(f"{k}: {v}")


def tumor_aliquots_of_all(d):
    s = set()
    for v in d.values():
        s |= v
    return s


if __name__ == "__main__":
    main()
