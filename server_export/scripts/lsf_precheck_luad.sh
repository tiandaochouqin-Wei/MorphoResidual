#!/bin/bash
# C1-LUAD split-half pre-check (D2, review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md
# section 0) under LSF. Discovery-cohort data ONLY -- no confirmatory slide, RNA or
# protein file is touched by this script or by split_replication_power.py.
#
# scp upload list (this file only -- split_replication_power.py, residual_analysis.py
# and morpho_env.sh are already on the cluster and are NOT modified by this branch;
# pinned/luad_batch_groups.csv + its .sha256 sidecar are P0 deliverables, built by a
# separate step per the frozen-rule draft D2, and must be uploaded before any `group`
# mode run -- this script FATALs cleanly if they are missing):
#   scp server_export/scripts/lsf_precheck_luad.sh user@cluster:/public/home/fjhui/ZW/scripts/
#
# Submit from the login node. Modes: smoke | random | random34 | group <col>.
#
#   Smoke (minutes; seed 99 means this NEVER computes a real P1 split -- it only
#   proves the pipeline runs end to end):
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_smoke.out \
#       lsf_precheck_luad.sh smoke
#
#   P1 -- 10 random splits, ONE bsub job (see the comment above the `random` case
#   below for why this must not be split across parallel per-seed jobs):
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p1.out \
#       lsf_precheck_luad.sh random
#
#   P1b -- 10 random splits at the stream split's size (n_test=34), the matched-n
#   reference for the descriptive `stream_c3n_vs_c3l` column (frozen-rule draft
#   section 0, D2 P0/D6d). ONE bsub job, same reasoning as `random` above (a
#   different --seed so its permutation draws do not collide with P1's):
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p1b_random34.out \
#       lsf_precheck_luad.sh random34
#
#   P2 -- one bsub job per group column. Seven columns gate GO/MARGINAL/NO-GO
#   (D2 P0: q_early_vs_late, q_2017Q3_vs_rest, op_balanced_0..op_balanced_4):
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p2_q_early_vs_late.out \
#       lsf_precheck_luad.sh group q_early_vs_late
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p2_q_2017Q3_vs_rest.out \
#       lsf_precheck_luad.sh group q_2017Q3_vs_rest
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p2_op_balanced_0.out \
#       lsf_precheck_luad.sh group op_balanced_0
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p2_op_balanced_1.out \
#       lsf_precheck_luad.sh group op_balanced_1
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p2_op_balanced_2.out \
#       lsf_precheck_luad.sh group op_balanced_2
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p2_op_balanced_3.out \
#       lsf_precheck_luad.sh group op_balanced_3
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p2_op_balanced_4.out \
#       lsf_precheck_luad.sh group op_balanced_4
#
#   An 8th group column, `stream_c3n_vs_c3l`, is descriptive only (D2 P0/D6d) and
#   never gates GO/MARGINAL/NO-GO (D2 P2); it is read against the P1b reference
#   above. It works through the same generic `group <col>` mode, no code change
#   needed once pinned/luad_batch_groups.csv carries the column:
#     bsub -q smp -n 8 -R "span[hosts=1]" -o precheck_luad_p2_stream_c3n_vs_c3l.out \
#       lsf_precheck_luad.sh group stream_c3n_vs_c3l
set -u

SCRIPTS_DIR=/public/home/fjhui/ZW/scripts
SPLIT_PY="${SCRIPTS_DIR}/split_replication_power.py"
GROUP_CSV="${SCRIPTS_DIR}/pinned/luad_batch_groups.csv"
GROUP_CSV_SHA_SIDECAR="${GROUP_CSV}.sha256"
# D2 (review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 0) fixes this against
# "the one recorded by the C1-UCEC run" as a LITERAL, not a value read back out of a
# mutable JSON file that a later re-run could silently overwrite. This is the exact
# sha256 c1_ucec_result.json recorded on 2026-09-16 -- kept here so a fixed check
# does not depend on that file staying untouched.
WANT_SPLIT_HASH=79530c075c0f7e75bf06d3c2f1edfe8792cbf4d6c71ae1d32e60bf44f31a0397
REF_UCEC_C1_RESULT=/public/home/fjhui/ZW/ucec_c1/results/c1_ucec_result.json
PY=/public/home/fjhui/miniconda3/bin/python

fatal() { echo "FATAL: $*" >&2; exit 1; }

MODE="${1:-}"
[ -n "${MODE}" ] || fatal "usage: lsf_precheck_luad.sh smoke|random|random34|group <col>"
shift || true
GROUP_COL=""
if [ "${MODE}" = "group" ]; then
  GROUP_COL="${1:-}"
  [ -n "${GROUP_COL}" ] || fatal "group mode needs a column name, e.g. op_balanced_0"
  shift || true
fi
# D2's parameters (--n-test/--splits/--perms/--top-frac/--seed) are frozen and must
# not be overridable from the command line: argparse keeps the LAST value it sees,
# so a stray extra argument here would silently change the pre-check without anyone
# noticing. Refuse any argument beyond mode [col].
[ $# -eq 0 ] || fatal "unexpected extra argument(s) after mode${GROUP_COL:+ <col>}: '$*' -- D2's parameters are frozen and cannot be overridden from the command line."

echo "=== job start: host=$(hostname) date=$(date) cores=$(nproc) mode=${MODE} ${GROUP_COL} ==="

export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:${LD_LIBRARY_PATH:-}
unset MORPHO_ROOT MORPHO_RNA_MANIFEST MORPHO_ALIQUOT_XWALK MORPHO_SLIDE_MAP \
      MORPHO_PROTEIN_TSV MORPHO_WSI_EMB_DIR MORPHO_OUT MORPHO_N_WORKERS MORPHO_SIG_CSV

cd "${SCRIPTS_DIR}" || fatal "cannot cd to ${SCRIPTS_DIR}"
# shellcheck disable=SC1091
source ./morpho_env.sh luad
rc=$?
[ ${rc} -eq 0 ] || fatal "morpho_env.sh luad returned ${rc} -- see its own message above."

echo "MORPHO_ROOT=${MORPHO_ROOT:-<unset>}"
echo "MORPHO_OUT=${MORPHO_OUT:-<unset>}"

# --- environment sanity: the exact bug this guards against -----------------------
# On 2026-09-15 a stale MORPHO_OUT inherited from an earlier `source morpho_env.sh
# <other cohort>` in the same shell silently redirected a whole UCEC run's output
# into pdac/results (see lsf_c1_test.sh's comment on the same incident). Refuse to
# run unless MORPHO_ROOT really is the luad tree and MORPHO_OUT is really inside it.
[ -n "${MORPHO_ROOT:-}" ] || fatal "MORPHO_ROOT is unset after sourcing morpho_env.sh luad"
[ -n "${MORPHO_OUT:-}" ] || fatal "MORPHO_OUT is unset after sourcing morpho_env.sh luad"
root_base=$(basename "${MORPHO_ROOT}")
[ "${root_base}" = "luad" ] || fatal "basename(MORPHO_ROOT)='${root_base}', expected 'luad'."
case "${MORPHO_OUT}" in
  "${MORPHO_ROOT}"|"${MORPHO_ROOT}"/*) : ;;
  *) fatal "MORPHO_OUT='${MORPHO_OUT}' does not lie inside the luad tree (MORPHO_ROOT='${MORPHO_ROOT}')." ;;
esac
echo "environment OK: MORPHO_ROOT and MORPHO_OUT both confirmed under the luad tree."

# --- hash checks (unconditional) --------------------------------------------------
[ -f "${SPLIT_PY}" ] || fatal "missing ${SPLIT_PY}"
got_split_hash=$(sha256sum "${SPLIT_PY}" | awk '{print $1}')
if [ "${got_split_hash}" != "${WANT_SPLIT_HASH}" ]; then
  fatal "split_replication_power.py sha256 mismatch: got=${got_split_hash} want(D2 literal)=${WANT_SPLIT_HASH}. Refusing to run an unverified copy of the shared statistic code."
fi
echo "hash OK: split_replication_power.py matches the D2-fixed literal (${got_split_hash})."
# Secondary, non-fatal cross-check: does the published UCEC result still agree with
# the literal above? A mismatch here would mean c1_ucec_result.json was regenerated
# after 2026-09-16 -- report it, but the literal above is what gates this run, per D2.
if [ -f "${REF_UCEC_C1_RESULT}" ]; then
  ucec_json_hash=$("${PY}" -c "
import json, sys
print(json.load(open(sys.argv[1]))['sha256_split_replication_power'])
" "${REF_UCEC_C1_RESULT}" 2>/dev/null) || ucec_json_hash="<unreadable>"
  if [ "${ucec_json_hash}" != "${WANT_SPLIT_HASH}" ]; then
    echo "WARNING: ${REF_UCEC_C1_RESULT}'s sha256_split_replication_power (${ucec_json_hash}) no longer matches the D2 literal -- the published UCEC result may have been regenerated. This does not block the LUAD pre-check, but flag it before trusting any UCEC comparison." >&2
  fi
fi

if [ "${MODE}" = "group" ]; then
  [ -f "${GROUP_CSV}" ] || fatal "group mode: missing ${GROUP_CSV} (D2 P0 deliverable -- build and upload it first, it is not produced by this script)."
  [ -f "${GROUP_CSV_SHA_SIDECAR}" ] || fatal "group mode: missing sidecar ${GROUP_CSV_SHA_SIDECAR}."
  got_csv_hash=$(sha256sum "${GROUP_CSV}" | awk '{print $1}')
  want_csv_hash=$(awk '{print $1}' "${GROUP_CSV_SHA_SIDECAR}")
  if [ "${got_csv_hash}" != "${want_csv_hash}" ]; then
    fatal "luad_batch_groups.csv sha256 mismatch: got=${got_csv_hash} want(sidecar)=${want_csv_hash}. The P0 file on disk does not match what was recorded when it was built -- refusing to run group='${GROUP_COL}' on a possibly-changed split."
  fi
  echo "hash OK: luad_batch_groups.csv matches its .sha256 sidecar (${got_csv_hash})."
fi

export MORPHO_N_WORKERS="${LSB_DJOB_NUMPROC:-8}"
mkdir -p "${MORPHO_OUT}"

case "${MODE}" in
  random)
    # This MUST be ONE job, not 10 parallel per-seed jobs. Inside
    # split_replication_power.py's main(), split s uses
    # RandomState(a.seed * 1000 + s) to draw the sel/test partition, but the
    # PERMUTATION seed for that same split is
    # RandomState(zlib.crc32(f"perm|{cohort}|{label}".encode())), keyed on the split
    # LABEL (here the integer s, stringified), not on --seed. Running each split as
    # its own job with a different --seed to "parallelise" would therefore give each
    # job the SAME permutation seed per label collision risk aside, it does not
    # actually parallelise anything: --splits still means "how many splits this one
    # invocation loops over with s=0..--splits-1", so 10 separate single-split jobs
    # would just repeat s=0 ten times, not cover s=0..9. One job, --splits 10.
    "${PY}" -u "${SPLIT_PY}" \
      --n-test 50 --splits 10 --perms 200 --top-frac 0.15 --seed 0 \
      --workers "${MORPHO_N_WORKERS}" --tag p1_random
    ;;
  random34)
    # P1b (frozen-rule draft section 0, D2 P0/D6d): the matched-n reference for
    # the descriptive `stream_c3n_vs_c3l` group column, at the same n_test=34
    # the stream split itself uses. Identical one-job reasoning to `random`
    # above -- split s's permutation seed is keyed on the split LABEL (the
    # stringified integer s), not on --seed, so 10 parallel per-seed jobs would
    # not cover s=0..9 and would collide on the permutation draw. --seed is 1
    # here (not 0) only so this job's own sel/test partition draws differ from
    # P1's; it does not change the per-split permutation-seed keying above.
    "${PY}" -u "${SPLIT_PY}" \
      --n-test 34 --splits 10 --perms 200 --top-frac 0.15 --seed 1 \
      --workers "${MORPHO_N_WORKERS}" --tag p1b_random34
    ;;
  group)
    "${PY}" -u "${SPLIT_PY}" \
      --group-csv "${GROUP_CSV}" --group-col "${GROUP_COL}" \
      --perms 200 --top-frac 0.15 \
      --workers "${MORPHO_N_WORKERS}" --tag "p2_${GROUP_COL}"
    ;;
  smoke)
    # Seed 99 means this never computes a real P1 split (the real P1 uses --seed 0);
    # this mode exists only to prove the environment/hash checks and the pipeline run
    # end to end in a couple of minutes, not to produce a usable pre-check number.
    "${PY}" -u "${SPLIT_PY}" \
      --n-test 50 --splits 1 --perms 5 --seed 99 \
      --workers "${MORPHO_N_WORKERS}" --tag smoke
    ;;
  *)
    fatal "unknown mode '${MODE}', expected smoke|random|random34|group <col>"
    ;;
esac
rc=$?
echo "=== job end rc=${rc} date=$(date) ==="
exit ${rc}
