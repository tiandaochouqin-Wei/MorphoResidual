#!/usr/bin/env python
"""Write review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256 (the sidecar local_brca_setlevel.py checks).

Run this ONLY when the pre-specification .md is written and signed off, and BEFORE RUNBOOK step 6 (no matrix exists yet, none is listed).
It lists (v3, adversarial review M4 -- everything section 4 of the pre-specification says is frozen): the pre-specification and
its pre-review copy, local_brca_setlevel.py, the primary and every secondary set file, the inventory script and its json/tsv
output, every server script of the arm (cptac2_replicate.py, cptac2_residual.py, extract_phikon.py, ts_slide_qc.py,
brca_tested_genes.py, both LSF wrappers), RUNBOOK.sh, FREEZE_CHECKLIST.txt, this script, test_synthetic_pipeline.py, the frozen
189-file TS list (dryrun/primary_all/{ts_slide_selection.tsv, manifest_brca_ts_slides.txt, md5_brca_ts.txt}), vial_match.tsv
and its summary (and vial_match.py), the GDC inventory (gdc_inventory_ts.py, inventory_ts_files.tsv, inventory_ts_case_summary.tsv,
inventory_summary.json, inventory_dx_cases.tsv), the discovery-universe files (setlevel/{universe_common.tsv, symbol_qc.json,
hgnc_drift_summary.tsv}, hgnc_drift_precheck.py), and the step-0c output brca_tested_genes.{tsv,txt} as it arrived (hashed only,
not opened), plus any --extra file (e.g. an HGNC crosswalk). The "# written" line of the sidecar is the UTC time (and date) that
the paper text quotes as "specified on [date]". Not listed (they do not exist yet): the server's dx_arm_inputs_outputs.sha256 and
the matrices.
Paths are relative to MorphoResidual_paper/. Refuses to overwrite.
  python make_setlevel_sidecar.py --prespec review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md [--extra REL ...]
"""
import argparse
import datetime
import hashlib
import os
import sys

import local_brca_setlevel as L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prespec", required=True, help="relative path (from MorphoResidual_paper/) of the pre-specification .md")
    ap.add_argument("--extra", action="append", default=[])
    a = ap.parse_args()
    if not a.prespec.endswith(".md"):
        sys.exit("--prespec must be the .md pre-specification")
    d = "review/BRCA_TS_ARM_2026-10-08/"
    files = [a.prespec, a.prespec[:-3] + "_pre-review.md", L.SELF_REL, L.PRIMARY_SET, *L.SECONDARY_SETS,
             d + "brca_selection_set_inventory.py",
             d + "setlevel/candidate_sets_inventory.json", d + "setlevel/candidate_sets_inventory.tsv",
             d + "setlevel/universe_common.tsv", d + "setlevel/symbol_qc.json", d + "setlevel/hgnc_drift_summary.tsv",
             d + "hgnc_drift_precheck.py",
             d + "scripts/cptac2_replicate.py", d + "scripts/cptac2_residual.py", d + "scripts/extract_phikon.py",
             d + "scripts/ts_slide_qc.py", d + "scripts/brca_tested_genes.py",
             d + "scripts/lsf_cptac2_ts_extract.sh", d + "scripts/lsf_cptac2_ts_residual.sh",
             d + "RUNBOOK.sh", d + "FREEZE_CHECKLIST.txt", d + "make_setlevel_sidecar.py", d + "test_synthetic_pipeline.py",
             d + "dryrun/primary_all/ts_slide_selection.tsv", d + "dryrun/primary_all/manifest_brca_ts_slides.txt",
             d + "dryrun/primary_all/md5_brca_ts.txt", d + "vial_match.py", d + "vial_match.tsv", d + "vial_match_summary.json",
             d + "gdc_inventory_ts.py", d + "inventory_ts_files.tsv", d + "inventory_ts_case_summary.tsv",
             d + "inventory_summary.json", d + "inventory_dx_cases.tsv",
             d + "brca_tested_genes.tsv", d + "brca_tested_genes.txt", *a.extra]
    if os.path.exists(L.SIDECAR):
        sys.exit(f"REFUSED: {L.SIDECAR} exists; a frozen sidecar is not overwritten (move it aside deliberately if the freeze is redone)")
    lines = [f"# sha256 of the frozen files for {os.path.basename(a.prespec)} (paths relative to MorphoResidual_paper/)",
             f"# written {datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}"]
    for rel in dict.fromkeys(files):
        p = os.path.join(L.PAPER, rel)
        if not os.path.exists(p):
            sys.exit(f"missing: {rel}")
        lines.append(f"{L.sha256(p)} *{rel}")
    open(L.SIDECAR, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
