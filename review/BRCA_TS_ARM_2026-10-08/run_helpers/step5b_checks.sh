#!/bin/bash
# Step 5b of RUNBOOK.sh (frozen; sidecar review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256), lines 24-26 and 112-143
# copied byte for byte, so the block can be run with: bash step5b_checks.sh  (output appended to $ROOT/prerun_checks.txt).
# This helper is not part of the freeze; it adds nothing to the frozen commands.
ZW=/public/home/fjhui/ZW
ROOT=$ZW/cptac2_brca            # existing DX-arm folder: protein_brca.csv rna_case_file.tsv rna/ slides/ emb_phikon*/ results/
S=$ZW/scripts
{                                                    # everything up to the closing brace is appended to prerun_checks.txt
date -u; hostname
cd $ZW/scripts_ts_arm && chmod +x *.sh && sha256sum * | tee $ROOT/deployed_scripts_ts_arm.sha256
#   M4: every deployed script must equal its frozen copy:
grep '/scripts/' POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256 | grep -v '_v1\.' | sed 's#\*review/BRCA_TS_ARM_2026-10-08/scripts/#*#' > $ROOT/frozen_scripts.sha256
grep 'dryrun/primary_all/manifest_brca_ts_slides.txt' POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256 \
  | sed 's#\*review/BRCA_TS_ARM_2026-10-08/dryrun/primary_all/manifest_brca_ts_slides.txt#*frozen_manifest_brca_ts_slides.txt#' >> $ROOT/frozen_scripts.sha256
sha256sum -c $ROOT/frozen_scripts.sha256                      # expect 8 x OK (the 7 deployed scripts + the frozen list)
#   M7 (a): which DX embedding folders exist, and their tiling density read from the saved tile coordinates (x-step 256 px =
#   stride 0 / dense; 512 = stride 512 / sparse, ~4x fewer tiles):
ls -d $ROOT/emb_phikon*; ls -la --time-style=full-iso -d $ROOT/emb_phikon* $ROOT/results/brca_protein_results.csv
for d in $(ls -d $ROOT/emb_phikon*); do
/public/home/fjhui/miniconda3/bin/python - "$d" <<'PY'
import sys, glob
import numpy as np, torch
d = sys.argv[1]
for f in sorted(glob.glob(d + "/*.pt"))[:5]:
    c = torch.load(f, map_location="cpu").get("coords")
    if c is None:
        print(d, f.split("/")[-1][:28], "no coords"); continue
    xs = np.unique(np.asarray(c)[:, 0]); step = int(np.gcd.reduce(np.diff(xs))) if len(xs) > 1 else -1
    print(d, f.split("/")[-1][:28], "n_tiles", len(c), "x-step", step)
PY
done
#   M7 (b): which folder the PUBLISHED DX job read (the published script prints no EMB path): its log / submission line and
#   any MORPHO_WSI_EMB_DIR in the shell history; and the tile counts of results/dx_slide_qc.tsv (step 5) for the same slides.
ls -la --time-style=full-iso $ROOT/*.log $ZW/*.log 2>/dev/null | head -50
grep -h -n -E "MORPHO_WSI_EMB_DIR|cptac2_residual|extract_phikon.*cptac2_brca" ~/.bash_history 2>/dev/null | tail -20
#   M7 (c): the published copies in $S must be the scripts the paper describes (CR-normalised sha256):
tr -d '\r' < $S/cptac2_residual.py | sha256sum      # expect 1e12432f4a3aedf7b5091cb57c66e224c0ac2e3c85054dd3f8787ff77f2a8fda (= MorphoResidual_paper/cptac2_residual.py)
tr -d '\r' < $S/extract_phikon.py  | sha256sum      # expect 4c616526b6e03a3dbd8ab6f40f85d524af562b118dfb4341594249a7b0ef273c (= server_export/scripts/extract_phikon.py)
} 2>&1 | tee -a $ROOT/prerun_checks.txt
