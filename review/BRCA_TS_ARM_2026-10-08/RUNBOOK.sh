#!/bin/bash
# RUNBOOK -- TCGA-BRCA frozen Tissue-Slide (TS) arm on the CPTAC-2 BRCA (PDC000173) replication.
# NOT meant to be run end to end: paste step by step on the HPC (login node mn02 has internet; GPU nodes may not).
# Governing document: review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md (all times there are UTC). Sizes/counts in comments come from the
# GDC query of 2026-10-07 UTC (local date 2026-10-08); GB = 10^9 bytes.
#
# STATE AT FREEZING (v3, 2026-10-07/08 UTC): steps 0, 0c and 1 have run on the server and step 2 (download) was started: step 1
# printed 189/232 files kept for 102/102 cases, 59.3 GB; step 0c printed 'tested (n>=30, in RNA): 10031'. Steps 3-9 have not run.
# Five of the seven deployed scripts (cptac2_replicate.py, extract_phikon.py, ts_slide_qc.py, brca_tested_genes.py,
# lsf_cptac2_ts_extract.sh) are byte-identical to the sidecar copies. The v2 cptac2_residual.py and lsf_cptac2_ts_residual.sh deployed
# with them have not run, are not kept locally, and are replaced by the v3 copies at step 5b.
#
# Blindness: steps 0-5 compute no statistic relating morphology to protein or mRNA. Step 0 hashes the published inputs/results
# without opening them; step 0c reads per-protein non-missing counts and RNA gene-symbol labels only (no morphology involved);
# steps 1-2 read GDC/PDC metadata and the case column of rna_case_file.tsv; steps 3-5 touch images, slide headers and tensor
# shapes/coordinates only. Do NOT open any *_results*.csv / summary*.csv / *_incr_matrix*.npz (step 6+) before the sidecar
# POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256 is written and released.
# v1 -> RUNBOOK_v1.sh; v2 added step 0c, the increment-matrix output and step 9; v3 (2026-10-07 UTC, adversarial review of the
# pre-specification) adds step 4b (skip list), step 5b (server checks M4/M7), host pinning (N6), the DX gate inputs (M3) and the
# corrected sidecar name; the freeze point moves to before any step-5 output is read. v3 final (2026-10-08 UTC, cross-file review)
# adds step 5c (read-only coverage), the host echo at step 6, the failed-run rules at step 7, the published-file check inputs at
# step 8/9 (dx_arm_inputs_outputs.sha256 + results/brca_protein_results.csv), and the one-run local read-out.

ZW=/public/home/fjhui/ZW
ROOT=$ZW/cptac2_brca            # existing DX-arm folder: protein_brca.csv rna_case_file.tsv rna/ slides/ emb_phikon*/ results/
S=$ZW/scripts

# ---------------------------------------------------------------- 0. deploy patched scripts to a SEPARATE folder + inspect DX state  [RAN]
# (Windows side) scp review/BRCA_TS_ARM_2026-10-08/scripts/{cptac2_replicate.py,cptac2_residual.py,extract_phikon.py,
#   ts_slide_qc.py,brca_tested_genes.py,lsf_cptac2_ts_extract.sh,lsf_cptac2_ts_residual.sh}  ->  $ZW/scripts_ts_arm/     (same file names; the published
#   copies in $S stay untouched; gdc_download.py is NOT patched -- it takes any manifest -- so keep using $S/gdc_download.py;
#   do NOT copy the *_v1.* files)
#   (the v2 cptac2_residual.py and lsf_cptac2_ts_residual.sh deployed here have not run and are replaced by the v3 copies at step 5b)
mkdir -p $ZW/scripts_ts_arm && chmod +x $ZW/scripts_ts_arm/*.sh
ls -la $ROOT; ls $ROOT/emb_phikon* -d; ls $ROOT/emb_phikon | wc -l            # expect 107 DX .pt (97x1 + 5x2 slides)
ls ~/.cache/huggingface/hub | grep -i phikon                                   # must exist (HF_HUB_OFFLINE=1 in the extract job)
diskquota                                                                      # inode count is ~98% of 2,000,000; this run adds ~400 files
cut -f1 $ROOT/rna_case_file.tsv | tail -n +2 | sort -u > $ROOT/dx_case_list.tsv ; wc -l $ROOT/dx_case_list.tsv   # expect 102
sha256sum $ROOT/protein_brca.csv $ROOT/rna_case_file.tsv $ROOT/results/*.csv > $ROOT/dx_arm_inputs_outputs.sha256   # fingerprint of the published arm
#   (the line above only HASHES the published results; do not open them)

# ---------------------------------------------------------------- 0c. BRCA tested-gene list (labels + non-missing counts ONLY; ~1 min, login node ok)  [RAN: 10031]
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
/public/home/fjhui/miniconda3/bin/python $ZW/scripts_ts_arm/brca_tested_genes.py --root $ROOT --case-list $ROOT/dx_case_list.tsv --out $ROOT/brca_tested_genes.tsv
#   expect: 'tested (n>=30, in RNA): 10031  [published: 10031]'. A different number = stop and tell the analyst.
#   (Exception to "no protein/mRNA statistic": this step reads per-protein non-missing counts and RNA symbol labels only.)
#   The selection set was already fixed in the pre-specification BEFORE this output existed; brca_tested_genes.{tsv,txt} are
#   hashed into the sidecar as they arrived and are not opened before freezing. After freezing, the overlap of each set with the
#   BRCA-tested symbols is REPORTED only, by the read-only command of step 5c; it changes nothing except through the MIN_COVERAGE
#   rule (primary < 50% -> stop and amend before step 6; secondary < 50% -> untestable). DO NOT use
#   'brca_selection_set_inventory.py --brca-genes' for this: it rewrites the frozen set and inventory files.

# ---------------------------------------------------------------- 1. TS manifest (metadata only, minutes, login node, internet)  [RAN: 189/232, 102/102, 59.3 GB]
cd $ZW/scripts_ts_arm
python cptac2_replicate.py --project TCGA-BRCA --pdc-study PDC000173 --slide-type ts --ts-rule primary_all
#   expect: [cases] ... tri-modal 102 ; [ts] rule=primary_all: 189/232 files kept for 102/102 cases, 59.3 GB
#   (aborts if the re-derived case set != rna_case_file.tsv). The frozen list is dryrun/primary_all/ts_slide_selection.tsv
#   (189 files; sorted-file-name sha256 9576e733...e049): any difference is reported, never substituted, and the analysis uses
#   the intersection (the residual wrapper enforces this).
sha256sum $ROOT/manifest_brca_ts_slides.txt $ROOT/ts_slide_selection.tsv | tee $ROOT/ts_manifest.sha256
tail -n +2 $ROOT/manifest_brca_ts_slides.txt | cut -f2 | tr -d '\r' | LC_ALL=C sort | sha256sum      # expect 9576e7337466cb20960c281db252b1fa36a80dece76f8199ca1bed988364e049

# ---------------------------------------------------------------- 2. download (login node, background, resumable; 59.3 GB, 189 files)  [RUNNING at freezing]
nohup python $S/gdc_download.py $ROOT/manifest_brca_ts_slides.txt $ROOT/slides_ts 4 > $ROOT/sl_ts.log 2>&1 &
#   progress: ls $ROOT/slides_ts | wc -l ; tail -2 $ROOT/sl_ts.log ;  re-run the same line to retry FAILs (finished files skip)
#   verify when done:
cd $ROOT/slides_ts && md5sum -c $ROOT/md5_brca_ts.txt | grep -v ': OK' ; echo md5-check-done      # expect no lines before md5-check-done
ls $ROOT/slides_ts/*.svs | wc -l                                                                    # expect 189
#   A file that still fails after retries is NOT substituted: record it at step 4b with reason download_failed.

# ======================= FREEZE POINT (v3): the pre-specification is frozen HERE, while step 2 runs =======================
#   Locally: python make_setlevel_sidecar.py --prespec review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md
#   -> writes review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256 (UTC time in its header); release the .md, the sidecar and every
#   listed file on GitHub. This happens before any step-5 output (magnification, tile counts) exists or is read and before any
#   matrix exists. local_brca_setlevel.py refuses to run without the sidecar.

# ---------------------------------------------------------------- 3. smoke test, 3 slides (GPU, ~10-20 min) -> calibrates time + tile counts
bsub -q interactive -gpu "num=1" -o $ROOT/ts_smoke.log -e $ROOT/ts_smoke.err $ZW/scripts_ts_arm/lsf_cptac2_ts_extract.sh 0/1 --limit 3
#   read ts_smoke.log: 'saved N tiles x 768d' and wall time per slide; extrapolate x189.
#   NB --limit writes real .pt files into emb_phikon_ts; they are kept (skip-if-exists) by the full run.

# ---------------------------------------------------------------- 4. full extraction: 2 shards (interactive queue, 24 h, 2 jobs/user)
bsub -q interactive -gpu "num=1" -o $ROOT/ts_ext0.log -e $ROOT/ts_ext0.err $ZW/scripts_ts_arm/lsf_cptac2_ts_extract.sh 0/2
bsub -q interactive -gpu "num=1" -o $ROOT/ts_ext1.log -e $ROOT/ts_ext1.err $ZW/scripts_ts_arm/lsf_cptac2_ts_extract.sh 1/2
#   resume = resubmit the same line.  Done when:  ls $ROOT/emb_phikon_ts/*.pt | wc -l   -> 189 minus the slides of step 4b
#   No minimum tile count per slide: any slide with >= 1 tissue tile is used.

# ---------------------------------------------------------------- 4b. skip list (every slide without an embedding, with its reason; no substitution)
#   One line per frozen slide that has no .pt after all resubmissions: <file stem without .svs><TAB><reason>, reason one of
#   no_tissue ('WARNING: no tissue tiles' in ts_ext*.log), extract_failed ('FAILED on' in ts_ext*.log, still failing after
#   resubmission), download_failed (md5/download failure after retries). An EMPTY file means no slide was skipped.
grep -h -E "WARNING: no tissue tiles|FAILED on" $ROOT/ts_ext*.log          # evidence for the reasons
comm -23 <(tail -n +2 $ROOT/manifest_brca_ts_slides.txt | cut -f2 | tr -d '\r' | sed 's/\.svs$//' | LC_ALL=C sort) \
         <(ls $ROOT/emb_phikon_ts | grep '\.pt$' | sed 's/\.pt$//' | LC_ALL=C sort)       # the stems that need a line
#   then write $ROOT/ts_skipped_slides.txt by hand from the two outputs above (the wrapper refuses to run the TS arm unless the
#   embeddings equal exactly the frozen list minus these lines). Every skipped slide is listed in the supplement.

# ---------------------------------------------------------------- 5. header census (CPU, ~minutes; headers + tile counts only)
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
python $ZW/scripts_ts_arm/ts_slide_qc.py --raw-dir $ROOT/slides_ts --emb-dir $ROOT/emb_phikon_ts --out $ROOT/results_ts/ts_slide_qc.tsv
python $ZW/scripts_ts_arm/ts_slide_qc.py --raw-dir $ROOT/slides --emb-dir $ROOT/emb_phikon --out $ROOT/results/dx_slide_qc.tsv
#   Both mpp distributions go to the supplement as they are; no re-run, restriction or stratification by magnification follows,
#   whatever they show (pre-specification section 1).

# ---------------------------------------------------------------- 5b. SERVER CHECKS before step 6 (review M4, M7; record everything in $ROOT/prerun_checks.txt)
#   (Windows side) re-deploy the two v3 files and the frozen list, and copy the sidecar (none of the four has run on the server; the
#   other five deployed scripts stay as they are):
#     scp review/BRCA_TS_ARM_2026-10-08/scripts/{cptac2_residual.py,lsf_cptac2_ts_residual.sh}  -> $ZW/scripts_ts_arm/
#     scp review/BRCA_TS_ARM_2026-10-08/dryrun/primary_all/manifest_brca_ts_slides.txt          -> $ZW/scripts_ts_arm/frozen_manifest_brca_ts_slides.txt
#     scp review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256                                    -> $ZW/scripts_ts_arm/
#   (do not strip carriage returns from frozen_manifest_brca_ts_slides.txt: it is hashed byte for byte; the .py/.sh are LF already)
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
#   DECISION: if every check matches (deployed = frozen; the published DX job read $ROOT/emb_phikon; its x-step is 256; $S hashes
#   as expected) proceed to step 6. If ANYTHING differs, stop: write a dated amendment to the pre-specification (UTC time,
#   what differed, what is changed -- e.g. DX_EMB_DIR=<folder the published job read> for the dx job, or a Methods correction),
#   release it, and only then submit step 6.

# ---------------------------------------------------------------- 5c. COVERAGE of the frozen sets among the BRCA-tested symbols (LOCAL, after the freeze, before step 6)
#   READ-ONLY: reads the frozen set files and brca_tested_genes.txt, writes only the screen/tee file. Do NOT use
#   brca_selection_set_inventory.py --brca-genes here (it rewrites the frozen candidate_sets/*.tsv and candidate_sets_inventory.*, and
#   local_brca_setlevel.py would then REFUSE). Rules (pre-specification section 2): PRIMARY (pub_k3) < 50% -> stop; a dated amendment
#   is needed before any set-level statistic is read. SECONDARY < 50% -> 'untestable (coverage)'. Save the output next to the other
#   bring-back files (coverage_check.txt) and add it to prerun_checks.txt.
#   cd review/BRCA_TS_ARM_2026-10-08
#   python - <<'PY' | tee coverage_check.txt
#   import pandas as pd
#   brca = {l.strip() for l in open("brca_tested_genes.txt", encoding="utf-8") if l.strip()}
#   for s in ("pub_k3", "pub_k4", "pub_k2", "op_k2", "phi_k3"):
#       d = pd.read_csv("setlevel/candidate_sets/%s.tsv" % s, sep="\t"); sel = d[d["selected"] == 1]
#       c = float(sel["gene"].isin(brca).mean())
#       print(s, "selected", len(sel), "in_BRCA_tested", int(sel["gene"].isin(brca).sum()), "coverage", round(c, 4),
#             "universe_in_BRCA_tested", int(d["gene"].isin(brca).sum()), "of", len(d),
#             "PRIMARY: STOP (<50%)" if s == "pub_k3" and c < 0.5 else "untestable (coverage)" if c < 0.5 else "ok")
#   PY

# ---------------------------------------------------------------- 6. DX reproduction check (same code path as published; same host as step 7)
#   N6: pick ONE smp execution host and use it for all three jobs (the full-SVD PCA agrees across hosts only to ~3e-5):
bqueues -l smp | grep -A2 -i hosts                  # list the queue's hosts; choose one that is 'ok' in bhosts
HOST=__choose_one_smp_host__                        # <-- edit: the chosen host
echo "HOST=$HOST" | tee -a $ROOT/prerun_checks.txt  # record the pinned host (N6)
bsub -q smp -n 8 -m $HOST -o $ROOT/dx_repro.log -e $ROOT/dx_repro.err $ZW/scripts_ts_arm/lsf_cptac2_ts_residual.sh dx
#   FAILED OR INTERRUPTED JOBS (pre-specification section 4, item 6): a job that stops before writing its summary csv (the last file it
#   writes) is reported with its log and repeated ONCE, unchanged, on the same host: first rename its partial outputs with the suffix
#   _failed1 (kept, listed, brought back), then resubmit the same line with -o/-e ...rerun1.log/.err. If $HOST is unavailable BEFORE this
#   dx job has run, choose and record another host. If it is unavailable AFTER it has run, repeat this dx job (the host's gate) first on
#   another recorded host and rerun all three jobs there (earlier outputs kept and listed, not analysed); that is a dated amendment.
#   At most one repeat of any job and one change of host; if a repeat also fails the analysis is reported as not completed and no
#   further run is made. A DX gate of 'fail' is not a failed job. The wrapper accepts no extra argument for dx and only '--pool tile'
#   for ts.
#   (the wrapper clears MORPHO_WSI_EMB_DIR/MORPHO_OUT, reads $ROOT/emb_phikon unless an amendment set DX_EMB_DIR, and writes
#    $ROOT/results_dx_repro/, so the published results/ is never overwritten; it also writes results_dx_repro/brca_incr_matrix.npz
#    = the published arm's genes x 1001 increment matrix (col 0 observed, cols 1..1000 = RandomState(0..999) permutations).
#    The three-tier gate (exact / host / fail) is applied LOCALLY by local_brca_setlevel.py at step 9; nothing is opened here.)

# ---------------------------------------------------------------- 7. TS residual analysis (~5-6 h CPU each; same host as step 6)
bsub -q smp -n 8 -m $HOST -o $ROOT/ts_resid.log -e $ROOT/ts_resid.err $ZW/scripts_ts_arm/lsf_cptac2_ts_residual.sh ts                    # primary
bsub -q smp -n 8 -m $HOST -o $ROOT/ts_resid_tile.log -e $ROOT/ts_resid_tile.err $ZW/scripts_ts_arm/lsf_cptac2_ts_residual.sh ts --pool tile   # pre-declared sensitivity
#   the wrapper stops unless emb_phikon_ts holds exactly the frozen list (intersected with the step-1 manifest) minus
#   ts_skipped_slides.txt, and passes --expect-kept to cptac2_residual.py, which also stops on any non-TS/BS/MS slide label.
#   outputs: $ROOT/results_ts/{brca_protein_results_ts.csv, summary_brca_ts.csv, pool_manifest_brca_ts.csv, brca_incr_matrix_ts.npz,
#                              brca_protein_results_ts_tilepool.csv, summary_brca_ts_tilepool.csv, brca_incr_matrix_ts_tilepool.npz}
#   (the tilepool call rewrites pool_manifest_brca_ts.csv with an identical file -- it does not depend on --pool)
#   the three .npz are ~80 MB each (10,031 genes x 1001 float64, compressed); each job asserts col 0 == csv increments and that the
#   saved nulls reproduce every pval before writing them.

# ---------------------------------------------------------------- 8. what to bring back (small files < 50 MB; the 3 matrices ~240 MB)
#   $ROOT/{manifest_brca_ts_slides.txt, md5_brca_ts.txt, ts_slide_selection.tsv, ts_manifest.sha256, dx_case_list.tsv,
#          dx_arm_inputs_outputs.sha256, sl_ts.log, ts_smoke.log, ts_ext0.log, ts_ext1.log, ts_resid*.log, dx_repro*.log,
#          ts_skipped_slides.txt, prerun_checks.txt, deployed_scripts_ts_arm.sha256, frozen_scripts.sha256}
#          (+ any *_rerun1.log/.err and *_failed1 partial outputs of a repeated job)
#   $ROOT/results_ts/*      $ROOT/results_dx_repro/*
#   PUBLISHED DX FILES, USED ONLY BY THE GATE (do not open): $ROOT/results/brca_protein_results.csv  $ROOT/results/summary_brca.csv
#      and $ROOT/dx_arm_inputs_outputs.sha256 (step 0; local_brca_setlevel.py checks both files against it, and the published summary
#      against 102 cases / 10,031 tested / 0 significant / negative_frac_all 0.991 -- otherwise the gate is 'fail')
#   coverage_check.txt (step 5c, local)
#   SET-LEVEL MATRICES (needed by local_brca_setlevel.py; each also holds genes, n, r2_rna, r2_both, pval, fdr, the case ids, seeds, n_slides):
#      $ROOT/results_ts/brca_incr_matrix_ts.npz            -> --matrix ts=...        (PRIMARY arm)
#      $ROOT/results_dx_repro/brca_incr_matrix.npz         -> --matrix dx=...        (DX re-run = reproduction of the published arm)
#      $ROOT/results_ts/brca_incr_matrix_ts_tilepool.npz   -> --matrix ts_tile=...   (pre-declared tile-pool sensitivity)
#   $ROOT/results_ts/ts_slide_qc.tsv  $ROOT/results/dx_slide_qc.tsv
#   sha256sum of everything above (including the three .npz);  cat /proc/cpuinfo | grep 'model name' | head -1 ;
#   python -c "import numpy,sklearn,torch;print(numpy.__version__,sklearn.__version__,torch.__version__)"

# ---------------------------------------------------------------- 9. LOCAL set-level read-out (Windows; synthetic self-tests first)
#   cd review/BRCA_TS_ARM_2026-10-08
#   python local_brca_setlevel.py --selftest                          # synthetic only; must print SELFTEST PASSED
#   python test_synthetic_pipeline.py <scratch>                       # synthetic end-to-end: patched script -> npz -> statistic
#   python local_brca_setlevel.py --matrix ts=<..>/brca_incr_matrix_ts.npz --matrix dx=<..>/brca_incr_matrix.npz \
#          --matrix ts_tile=<..>/brca_incr_matrix_ts_tilepool.npz \
#          --dx-published-csv <..>/results/brca_protein_results.csv --dx-published-summary <..>/results/summary_brca.csv \
#          --dx-repro-summary <..>/results_dx_repro/summary_brca.csv --dx-published-sha256 <..>/dx_arm_inputs_outputs.sha256
#   (REFUSES unless review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256 exists and verifies, and refuses to overwrite an existing
#    setlevel_result.*; applies the published-file check, the DX gate, the case-set rules and the coverage rules; writes
#    setlevel/results/setlevel_result.{csv,json}; run ONCE. An invocation that exits FATAL before writing setlevel_result.* is not a
#    run: record it, correct the input and repeat; nothing in the frozen files changes.)
