#!/bin/bash
# C1-LUAD signed rule section 1.4 reader QC, run AFTER unblinding.
# Deviation recorded in review/C1_LUAD_ADDENDA.md Entry 4: the rule places this check
# before feature extraction; it was omitted there and is run now, whatever it shows.
#
# This wrapper only orchestrates. The comparison and its thresholds are the frozen
# qc_reader_equivalence.py (sha256 45aea4d8..., rule section 7.2), hash-checked below;
# the base VOLUME instance is chosen by the frozen extractor's own resolve_series().
#   Arm (i):  openslide (DICOM driver, base VOLUME) vs wsidicom (series directory) on
#             identical level-0 regions of the three confirmatory series the rule names.
#   Arm (ii): discovery C3L-00001-21, original SVS vs its IDC Leica-pathway DICOM copy,
#             openslide on both.
#
# Usage, from /public/home/fjhui/ZW/scripts:
#   bsub -q smp -n 1 -o qc_reader_luad.out -e qc_reader_luad.err \
#     ./lsf_qc_reader_c1_luad.sh <discovery C3L-00001-21 .svs path> <dir holding its IDC series>
set -u
cd /public/home/fjhui/ZW/scripts || exit 1
PY=/public/home/fjhui/miniconda3/bin/python
QC_SHA=45aea4d84d768f95870cc937079afee947dac7a6b16ea594422ab12fd1e315dc
got=$(sha256sum qc_reader_equivalence.py | cut -d' ' -f1)
if [ "$got" != "$QC_SHA" ]; then
  echo "FATAL: qc_reader_equivalence.py sha256 ${got} != ${QC_SHA}"
  exit 1
fi
echo "hash OK: qc_reader_equivalence.py (${got:0:12}...)"
if [ $# -ne 2 ]; then
  echo "FATAL: need <discovery svs> <discovery dicom root>"
  exit 1
fi

CONF_ROOT=/public/home/fjhui/ZW/c1_luad_slides_dicom/cptac_luad
DISC_SVS="$1"
DISC_ROOT="$2"
OUT=/public/home/fjhui/ZW/luad_c1/reader_qc
mkdir -p "${OUT}"
echo "=== reader QC start: host=$(hostname) date=$(date) ==="

base_of() {
  "${PY}" -c "import sys; from pathlib import Path; sys.path.insert(0, '.'); import extract_c1_dicom as X; r = X.resolve_series(Path(sys.argv[1])); print(r['base_file']); print(r['slide_id'], r['width'], r['height'], r['transfer_syntax'], file=sys.stderr)" "$1"
}
series_dir() {
  find "$1" -type d -name "SM_$2" | head -1
}

rc_all=0
for pair in "C3L-00444-23:1.3.6.1.4.1.5962.99.1.241014027.1310308374.1640918553867.2.0" \
            "C3L-03717-21:1.3.6.1.4.1.5962.99.1.247772177.1973902576.1640925312017.2.0" \
            "C3L-02513-23:1.3.6.1.4.1.5962.99.1.251556113.1760283742.1640929095953.2.0"; do
  sid=${pair%%:*}
  uid=${pair#*:}
  d=$(series_dir "${CONF_ROOT}" "${uid}")
  if [ -z "${d}" ]; then echo "FATAL: series directory for ${sid} not found"; exit 1; fi
  b=$(base_of "${d}") || { echo "FATAL: resolve_series failed for ${sid}"; exit 1; }
  echo ""
  echo "=== arm (i) ${sid}: openslide ${b} vs wsidicom ${d} ==="
  "${PY}" -u qc_reader_equivalence.py --a "${b}" --b "wsidicom:${d}" --png "${OUT}/qc_i_${sid}.png"
  r=$?
  echo "arm (i) ${sid} exit=${r}"
  [ ${r} -eq 0 ] || rc_all=1
done

d=$(series_dir "${DISC_ROOT}" "1.3.6.1.4.1.5962.99.1.250299211.778750893.1640927839051.2.0")
if [ -z "${d}" ]; then echo "FATAL: discovery C3L-00001-21 DICOM series not found under ${DISC_ROOT}"; exit 1; fi
b=$(base_of "${d}") || { echo "FATAL: resolve_series failed for C3L-00001-21"; exit 1; }
echo ""
echo "=== arm (ii) C3L-00001-21: SVS ${DISC_SVS} vs IDC copy ${b} ==="
"${PY}" -u qc_reader_equivalence.py --a "${DISC_SVS}" --b "${b}" --png "${OUT}/qc_ii_C3L-00001-21.png"
r=$?
echo "arm (ii) exit=${r}"
[ ${r} -eq 0 ] || rc_all=1

echo ""
echo "=== reader QC end: overall rc=${rc_all} (0 = every arm EQUIVALENT) date=$(date) ==="
exit ${rc_all}
