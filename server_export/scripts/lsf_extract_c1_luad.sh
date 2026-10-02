#!/bin/bash
# C1-LUAD confirmatory UNI extraction from the IDC DICOM copy. Models on the EXECUTED
# C1-UCEC precedent lsf_extract_c1.sh, adapted per
# review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 1.3 "Extraction" and the
# blinding order of section 1.6 (this step is 1.6 step 2: it may run only after D3
# registration -- step 1 -- and before any RNA/protein pull -- step 5).
#
# scp upload (this file only -- extract_c1_dicom.py and extract_features.py are already
# on the cluster and are NOT modified by this task):
#   scp server_export/scripts/lsf_extract_c1_luad.sh user@cluster:/public/home/fjhui/ZW/scripts/
#   ssh cluster chmod +x /public/home/fjhui/ZW/scripts/lsf_extract_c1_luad.sh
#
# --dicom-root below is a PLACEHOLDER path, analogous to the UCEC precedent's
# c1_ucec_slides_dicom/cptac_ucec. Nothing under it exists yet: the corresponding
# author must create this directory (or point --dicom-root at wherever the section 1.6
# step 2 download actually lands) and confirm the exact path BEFORE submitting. This
# script does not download anything itself -- it only reads already-downloaded DICOM.
#
# Submit from the login node. interactive queue: 24h limit, 2 jobs per user -> two
# shards (same reasoning as the UCEC precedent). Re-submitting the same shard resumes
# (extract_c1_dicom.py skips slides it already finished):
#   bsub -q interactive -gpu "num=1" -o c1_luad_extract_s0.log -e c1_luad_extract_s0.err \
#        /public/home/fjhui/ZW/scripts/lsf_extract_c1_luad.sh --shard 0 --nshards 2
#   bsub -q interactive -gpu "num=1" -o c1_luad_extract_s1.log -e c1_luad_extract_s1.err \
#        /public/home/fjhui/ZW/scripts/lsf_extract_c1_luad.sh --shard 1 --nshards 2
#
# Section 1.3 requires a --dry-run pass over all downloaded series before extraction
# ("The script resolves every series before its dry-run loop and stops at the first one
# it cannot resolve"):
#   bsub -q interactive -o c1_luad_extract_dryrun.log -e c1_luad_extract_dryrun.err \
#        /public/home/fjhui/ZW/scripts/lsf_extract_c1_luad.sh --dry-run
#
# HF_TOKEN is read from ~/.hf_token, not inherited from the submitting shell: the
# interactive queue was observed re-sourcing shell rc files inconsistently across
# submissions (HF_TOKEN_len 39/37/14 on three otherwise-identical bsub calls of the
# same --only smoke command on 2026-09-16), so an exported value cannot be trusted to
# survive to the job. `chmod 600 ~/.hf_token` beforehand; only this file is authoritative.
export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
if [ ! -s "$HOME/.hf_token" ]; then
  echo "FATAL: $HOME/.hf_token missing or empty" >&2
  exit 1
fi
HF_TOKEN=$(cat "$HOME/.hf_token")
export HF_TOKEN

# --only-file is FROZEN to the 135-ID IDC-availability list (section 1.3: "Run with
# --only-file pinned/c1_primary_slide_ids_luad_idc.txt. That list is the 136-ID list
# minus C3N-02142-23, which has no IDC series; the unchanged script exits on any listed
# ID it cannot find."). This is deliberately the 135-ID list, NOT
# pinned/c1_primary_slide_ids_luad.txt (the 136-ID list including C3N-02142-23, which
# has no IDC series and would make extract_c1_dicom.py exit on an unresolvable ID).
# --only-file, --dicom-root and --out-dir are hard-coded below rather than left to
# "$@", and rejected outright if repeated on the command line, so a stray flag in the
# submitted arguments cannot silently widen, narrow or redirect the frozen confirmatory
# slide set (Appendix item 5: "No argument pass-through that could override frozen
# parameters").
ONLY_FILE=/public/home/fjhui/ZW/scripts/pinned/c1_primary_slide_ids_luad_idc.txt
for a in "$@"; do
  case "${a}" in
    --only*|--dicom-root*|--out-dir*)
      echo "FATAL: '${a}' is frozen by this wrapper (--only-file/--dicom-root/--out-dir) and must not be passed on the command line -- see the comment above ONLY_FILE." >&2
      exit 1
      ;;
  esac
done

echo "=== job start: host=$(hostname) date=$(date) HF_TOKEN_len=${#HF_TOKEN} args=$* ==="
nvidia-smi -L
cd /public/home/fjhui/ZW/scripts || exit 1
/public/home/fjhui/miniconda3/bin/python -u extract_c1_dicom.py \
  --dicom-root /public/home/fjhui/ZW/c1_luad_slides_dicom/cptac_luad \
  --out-dir /public/home/fjhui/ZW/luad_c1/WSI/emb \
  --only-file "${ONLY_FILE}" "$@"
rc=$?
echo "=== job end rc=${rc} date=$(date) ==="
exit ${rc}
