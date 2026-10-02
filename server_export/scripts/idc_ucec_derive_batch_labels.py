#!/usr/bin/env python3
"""idc_ucec_derive_batch_labels.py -- POST HOC C1-UCEC operator-label extraction and
validation, D4 item 3 steps 1-2 (D6b).

Authority: review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 0, decision D4,
item 3 ("Batch-stratified re-read (RUN)"):

    1. Validate, then read. On the IDC copies of UCEC discovery slides, validate
       that the DICOM ImageComments `User =` field equals `aperio.User` in
       slide_acquisition_meta.tsv, as section 4 did for LUAD. Then read the same
       field from the headers of the UCEC confirmatory Primary Tumor slides. Those
       files are already on the HPC and were already analysed, so this is not a
       blinding question.
    2. Labels. A case's label is the mode over its slides in
       scripts/pinned/slide_type_map_ucec_c1.tsv (the 169 analysed Primary Tumor
       slides) whose header carries the field. Ties go to the first slide in
       sorted slide-ID order. Labels are frozen before the re-read. The D1
       qualification criteria are evaluated on the published UCEC analysed set.
    - If the attribute does not validate, or any D1 criterion fails on the UCEC
      analysed set, Methods state that C1-UCEC was tested under the unrestricted
      null only and carries no batch qualifier.

UCEC confirmatory data is ALREADY DOWNLOADED, ALREADY ANALYSED and ALREADY
PUBLISHED (review/C1_UCEC_FROZEN_RULE_2026-09-15.md, run 2026-09-16). Querying its
IDC DICOM headers is therefore explicitly authorised to touch the real network
(collection cptac_ucec) -- this is not a blinding concern, unlike anything on the
LUAD side of this project. Only DICOM header bytes are ever fetched (never
PixelData), following the identical range-GET / stop_before_pixels pattern already
validated for LUAD in review/C1_LUAD_STEP0B_2026-09-29/idc_step3_operator_header_probe.py
and idc_step3c_derive_batch_labels.py:
  - HTTP Range: bytes=0-262143 (256 KiB) first; escalate to 1 MiB then 4 MiB (cap)
    only if the header parse is incomplete.
  - Parse from BytesIO with pydicom stop_before_pixels=True, tolerant of
    truncation; bytes are discarded after parsing, never written to disk.
  - Every fetch is logged (bytes requested, URL, HTTP status, response sha256) to
    review/C1_UCEC_D4_2026-09-29/query_log_idc.tsv (utc_time, endpoint,
    query_summary, http_status, response_sha256, saved_path -- same columns as
    review/C1_LUAD_STEP0B_2026-09-29/query_log_idc.tsv).

Unlike the LUAD Step 0b scripts (which did a blind full-tag sweep to DISCOVER which
DICOM attribute carries the operator/date), this script tests ONLY ImageComments
(0020,4000)'s Aperio "User = " / "Date = " sub-fields directly -- that attribute is
already established (same CPTAC/Aperio acquisition pipeline, same discovery slide
population as LUAD; validated LUAD-side against every checked slide). D4 item 3's
own text asks for exactly this narrower validate-then-read, not a re-discovery
sweep.

What this script is NOT: it does not decide C1-UCEC's confirmatory-test reading
(unchanged, per D4's closing sentence -- "In every case, UCEC's published section 5
reading is unchanged"). It only derives labels plus a qualification verdict that
c1_ucec_stratified_reread.py (the second D4-item-3 script) consumes to decide
whether it may compute a stratified re-read, or must print the fallback sentence
above and stop.

Run: LOCALLY on this machine (this is a metadata/header probe against the public
IDC index, exactly like the LUAD Step 0b scripts -- there is no HPC dependency and
no RNA/protein/WSI data is touched). NOT wrapped in an LSF script for that reason
(compare lsf_c1_ucec_strat.sh, which exists only because c1_ucec_stratified_reread.py
needs the HPC-resident RNA/protein/WSI matrices and multiprocessing workers).

  python idc_ucec_derive_batch_labels.py

idc_index has a known import-time bug on this box: idc_index_data's bundled SQL
text is read without an explicit encoding, and the Windows default codepage (gbk
here) cannot decode it (UnicodeDecodeError at import). Set PYTHONUTF8=1 first:
  (cmd)        set PYTHONUTF8=1 && python idc_ucec_derive_batch_labels.py
  (PowerShell) $env:PYTHONUTF8=1; python idc_ucec_derive_batch_labels.py
This is an environment fix, not a code change to any pinned/frozen file. Confirmed
2026-09-29: `python -c "import idc_index"` fails without it, succeeds with it, on
this exact interpreter.

`from idc_index import IDCClient` is imported lazily, inside main() / _get_idc_client(),
specifically so the pure functions below (parse_aperio_string, scan_quarter_from_date,
modal_case_label, evaluate_d1_criteria) can be imported and unit-tested on synthetic
data without a working idc_index install -- see the tie-break / D1-criteria synthetic
test run alongside this file's delivery.

Outputs:
  server_export/scripts/pinned/c1_ucec_batch_labels.tsv
      case_submitter_id, operator_uuid, operator_purity, scan_quarter,
      scan_quarter_purity, n_primary_slides_used -- IDENTICAL schema to
      pinned/c1_luad_batch_labels.tsv. Written ONLY if the operator attribute
      validates on discovery (mirrors idc_step3c_derive_batch_labels.py: no
      validation -> no labels file).
  review/C1_UCEC_D4_2026-09-29/idc_ucec_discovery_validation.tsv
      per-slide discovery validation rows (slide_id, case, matched, operator_uuid,
      scan_date, truth_operator_uuid, truth_scan_date, agree_operator, agree_date).
  review/C1_UCEC_D4_2026-09-29/idc_ucec_confirmatory_per_slide_labels.tsv
      per-slide confirmatory rows (slide_id, case, operator_uuid, scan_date,
      scan_quarter, note) -- mirrors idc_step3c_per_slide_labels.tsv.
  review/C1_UCEC_D4_2026-09-29/c1_ucec_d4_summary.json
      validated_operator, validated_date, discovery coverage/agreement counts,
      operator_axis / scan_quarter_axis D1 diagnostic reports (coverage, strata,
      largest-stratum share, permutable share, and the four D1 booleans),
      "qualifies" (the single field c1_ucec_stratified_reread.py gates on:
      validated_operator AND all four D1 criteria hold for the operator axis on
      the published 138-case analysed set), sha256_output_labels (sha256 of
      c1_ucec_batch_labels.tsv, when written -- c1_ucec_stratified_reread.py
      asserts this against the labels file it actually reads, so a stale or
      hand-edited labels file cannot silently be substituted), and
      generated_at_utc. Named c1_ucec_d4_summary.json (not idc_ucec_...) to
      match this project's c1_<cohort>_<thing> convention for pinned,
      cohort-facing artifacts and the filename c1_ucec_stratified_reread.py
      actually expects.
  review/C1_UCEC_D4_2026-09-29/query_log_idc.tsv
      one row per HTTP range-read (see above).
"""
import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PAJINSEN = Path(r"D:/claude/pajinsen/MorphoResidual_paper")
SCRIPTS_DIR = PAJINSEN / "server_export/scripts"
PINNED = SCRIPTS_DIR / "pinned"
REVIEW_OUT_DEFAULT = PAJINSEN / "review/C1_UCEC_D4_2026-09-29"

N_CONF_SLIDES_EXPECTED = 169
N_CASES_EXPECTED = 138
N_DISC_SLIDES_EXPECTED = 393  # slide_type_map_ucec.tsv, sample_type == "Primary Tumor"

FIRST_RANGE = 262144           # 256 KiB
STEP_SIZES = [262144, 1024 * 1024, 4 * 1024 * 1024]  # 256 KiB, 1 MiB, 4 MiB cap
MAX_RANGE = STEP_SIZES[-1]

MIN_PERMUTABLE = 0.50  # sitepack / section 4 literal, reused for the D1 report


# =============================================================================
# Pure functions -- no I/O, no network, no globals. Importable and unit-testable
# without idc_index (which is imported lazily, only inside the network-touching
# functions below).
# =============================================================================

def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def parse_aperio_string(s):
    """Parse 'Aperio Image Library ...|AppMag = 20|...|Date = 07/06/16|...|
    User = <uuid>|...' into a dict. Returns {} if not a recognisable Aperio
    ImageComments string. Identical logic to
    review/C1_LUAD_STEP0B_2026-09-29/idc_step3c_derive_batch_labels.py's
    parse_aperio_string, copied rather than imported (script is self-contained;
    review/ is not shipped anywhere this script runs)."""
    if not s or "|" not in s:
        return {}
    out = {}
    for part in s.split("|"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def scan_quarter_from_date(date_str):
    """date_str like '07/06/16' (MM/DD/YY, Aperio convention)."""
    try:
        dt = datetime.strptime(date_str, "%m/%d/%y")
    except (ValueError, TypeError):
        return None
    q = (dt.month - 1) // 3 + 1
    return f"{dt.year}Q{q}"


def pick_base_volume_instance(smi_sub):
    """Mirror extract_c1_dicom.py's resolve_series() / Step 0b's own copy of it:
    VOLUME instances only (image_type[2] == 'VOLUME'), largest
    TotalPixelMatrixColumns*Rows wins. `smi_sub` is a DataFrame of
    sm_instance_index rows already filtered to one SeriesInstanceUID."""
    vols = []
    for _, r in smi_sub.iterrows():
        it = r["ImageType"]
        it_list = list(it) if it is not None else []
        if len(it_list) > 2 and it_list[2] == "VOLUME":
            area = int(r["TotalPixelMatrixColumns"] or 0) * int(r["TotalPixelMatrixRows"] or 0)
            vols.append((area, r))
    if not vols:
        return None
    vols.sort(key=lambda t: t[0], reverse=True)
    return vols[0][1]


def modal_case_label(rows, val_key="value"):
    """The D1 'Stratum' / D4 item 3 tie-break rule, applied literally: "the mode
    over its slides ... with ties going to the first slide in sorted slide-ID
    order." `rows` is a list of {"slide_id": str, val_key: str-or-None} dicts for
    ONE case; the caller must already have restricted it to that case's slides.

    Unlike idc_step3c_derive_batch_labels.py's modal_per_case (which took
    Counter(vals).most_common()[0], breaking ties by first-INSERTION order in
    whatever order the header-fetch loop happened to visit slides), this sorts by
    slide_id FIRST so the tie-break is deterministic and matches the rule's text
    exactly regardless of fetch/iteration order. On the one case this project has
    actually seen a tie (C3N-02281, LUAD, D1), first-insertion order and first-
    sorted-slide-ID order coincided (the header probe visited slides in ID order),
    so this is not a behaviour change there -- it is a correctness fix made
    explicit for this fresh script, where no such coincidence can be assumed.

    Returns (label_or_None, purity_or_None, n_used) where n_used is the count of
    slides with a non-null/non-empty value (the denominator of purity), matching
    idc_step3c's n_primary_slides_used semantics exactly.
    """
    ordered = sorted(rows, key=lambda r: r["slide_id"])
    vals = [r[val_key] for r in ordered if r.get(val_key)]
    n_used = len(vals)
    if n_used == 0:
        return None, None, 0
    counts = Counter(vals)
    max_count = max(counts.values())
    chosen = None
    for v in vals:  # already in sorted slide-id order
        if counts[v] == max_count:
            chosen = v
            break
    purity = counts[chosen] / n_used
    return chosen, purity, n_used


def evaluate_d1_criteria(labels, n_cases_total):
    """The D1 'Stratum' qualification bullet, as a pure function on a plain
    iterable of per-case label strings (empty string / None = unlabelled).
    `n_cases_total` is the denominator for coverage (the published analysed-set
    size, e.g. 138), which may exceed len(labels) if some cases have no row at
    all (defensive; in practice the caller passes exactly n_cases_total labels,
    one per case, using "" for unlabelled).

    Mirrors c1_run_test_luad.py's build_within_operator_perms diagnostics (same
    four booleans, same thresholds, same MIN_PERMUTABLE=0.50 literal) and
    idc_step3c_derive_batch_labels.py's axis_report -- reimplemented here as a
    side-effect-free function so it can be unit-tested directly, per the task's
    own instruction.

    "Permutable" follows section 4's definition used throughout this project:
    a case is permutable if its stratum has >= 2 members (singletons and
    unlabelled cases are necessarily fixed under a within-stratum permutation).
    """
    labels = list(labels)
    if len(labels) != n_cases_total:
        raise ValueError(f"len(labels)={len(labels)} != n_cases_total={n_cases_total}")
    vals = [v for v in labels if v]
    n_labelled = len(vals)
    counts = Counter(vals)
    n_levels = len(counts)
    n_strata_ge3 = sum(1 for v, n in counts.items() if n >= 3)
    largest = max(counts.values(), default=0)
    largest_frac = (largest / n_labelled) if n_labelled else 1.0
    coverage = n_labelled / n_cases_total if n_cases_total else 0.0
    fixed = (n_cases_total - n_labelled) + sum(1 for v, n in counts.items() if n == 1)
    permutable_frac = 1.0 - fixed / n_cases_total if n_cases_total else 0.0

    d1_coverage_ge_80pct = coverage >= 0.80
    d1_ge_2_strata_ge_3_cases = n_strata_ge3 >= 2
    d1_largest_stratum_le_90pct = largest_frac <= 0.90
    d1_ge_50pct_permutable = permutable_frac >= MIN_PERMUTABLE
    qualifies = (d1_coverage_ge_80pct and d1_ge_2_strata_ge_3_cases
                 and d1_largest_stratum_le_90pct and d1_ge_50pct_permutable)

    return {
        "n_cases_total": n_cases_total,
        "n_labelled": n_labelled,
        "coverage_pct": round(100 * coverage, 1),
        "n_levels": n_levels,
        "n_singleton_cases": sum(1 for v, n in counts.items() if n == 1),
        "n_strata_with_ge3_cases": n_strata_ge3,
        "largest_stratum_share_of_labelled_pct": round(100 * largest_frac, 1) if n_labelled else None,
        "pct_cases_permutable": round(100 * permutable_frac, 1),
        "level_sizes": dict(sorted(counts.items(), key=lambda kv: -kv[1])),
        "D1_coverage_ge_80pct": d1_coverage_ge_80pct,
        "D1_ge_2_strata_ge_3_cases": d1_ge_2_strata_ge_3_cases,
        "D1_largest_stratum_le_90pct": d1_largest_stratum_le_90pct,
        "D1_ge_50pct_permutable": d1_ge_50pct_permutable,
        "D1_all_four_hold": qualifies,
    }


def slide_to_case(slide_id):
    m = re.match(r"^(C3[A-Z]-\d{5})-\d+$", str(slide_id))
    return m.group(1) if m else None


# =============================================================================
# I/O + network. idc_index / pydicom / requests are imported inside these
# functions (or lazily at call time), never at module level.
# =============================================================================

class QueryLog:
    """Same six columns / same append-as-you-go convention as
    review/C1_LUAD_STEP0B_2026-09-29/_idc_query_log.py's log_row, reimplemented
    (not imported) so this script has no runtime dependency on anything under
    review/."""
    COLS = ["utc_time", "endpoint", "query_summary", "http_status", "response_sha256", "saved_path"]

    def __init__(self, path: Path):
        self.path = path
        self._wrote_header = path.exists()

    def row(self, endpoint, query_summary, http_status, response_sha256, saved_path=""):
        is_new = not self._wrote_header
        with open(self.path, "a", encoding="utf-8", newline="\n") as f:
            if is_new:
                f.write("\t".join(self.COLS) + "\n")
                self._wrote_header = True
            vals = [utcnow(), str(endpoint), str(query_summary).replace("\t", " ").replace("\n", " "),
                    str(http_status), str(response_sha256), str(saved_path)]
            f.write("\t".join(vals) + "\n")
        print(f"[log] {endpoint} :: {query_summary} :: status={http_status} "
              f"sha256={str(response_sha256)[:12]}...")


def to_https(s3_url: str) -> str:
    assert s3_url.startswith("s3://"), s3_url
    rest = s3_url[len("s3://"):]
    bucket, key = rest.split("/", 1)
    return f"https://{bucket}.s3.amazonaws.com/{key}"


def fetch_header_bytes(https_url, series_uid, label, qlog):
    """Range-GET the DICOM header, enlarging stepwise up to MAX_RANGE. Identical
    protocol to idc_step3_operator_header_probe.py's fetch_header_bytes: 256 KiB
    first, escalate to 1 MiB then 4 MiB only if pydicom's stop_before_pixels parse
    is incomplete (missing SOPClassUID/Modality -- i.e. the file_meta group did
    not fully arrive), cap at 4 MiB, log every attempt. Returns
    (raw_bytes, http_status, bytes_requested, header_complete)."""
    import io
    import pydicom
    import requests
    data, status = b"", None
    for i, size in enumerate(STEP_SIZES):
        headers = {"Range": f"bytes=0-{size - 1}"}
        resp = requests.get(https_url, headers=headers, timeout=30)
        status = resp.status_code
        data = resp.content
        ok = status in (200, 206)
        try:
            ds = pydicom.dcmread(io.BytesIO(data), stop_before_pixels=True, force=True)
            complete = "SOPClassUID" in ds and "Modality" in ds
        except Exception:
            complete = False
        qlog.row(https_url, f"{label} header range-read attempt {i} bytes=0-{size - 1} series={series_uid}",
                  status, sha256_bytes(data), f"(discarded after parse; {len(data)} bytes)")
        if ok and complete:
            return data, status, size, True
        if size >= MAX_RANGE:
            return data, status, size, complete
    return data, status, STEP_SIZES[-1], False


def extract_image_comments(raw_bytes):
    """Parse with pydicom (tolerant of truncation) and return the ImageComments
    (0020,4000) string, or None if absent / unparseable."""
    import io
    import pydicom
    try:
        ds = pydicom.dcmread(io.BytesIO(raw_bytes), stop_before_pixels=True, force=True)
    except Exception:
        return None
    ic = ds.get("ImageComments", None)
    if ic is None:
        return None
    return str(ic)


def build_idc_series_lookup(idc_client, collection_id):
    """idx (c.index) filtered to this collection's SM series, left-joined to
    sm_index for ContainerIdentifier -- the GDC-slide-ID <-> IDC-series join key,
    exactly as review/C1_LUAD_STEP0B_2026-09-29/idc_step2_slide_join.py's
    `sm_luad` (there built for cptac_luad; here for --collection, default
    cptac_ucec, per D4 item 3's explicit authorisation)."""
    idx = idc_client.index
    sm = idc_client.sm_index
    sm_coll = idx[(idx["collection_id"] == collection_id) & (idx["Modality"] == "SM")].copy()
    sm_coll = sm_coll.merge(sm[["SeriesInstanceUID", "ContainerIdentifier"]],
                             on="SeriesInstanceUID", how="left")
    return sm_coll


def resolve_series_for_slide(sm_coll, slide_id, label):
    matches = sm_coll[sm_coll["ContainerIdentifier"] == slide_id]
    if len(matches) == 0:
        return None
    if len(matches) > 1:
        print(f"WARNING: {label} slide {slide_id} matches {len(matches)} IDC series "
              f"(using the first; all others ignored)")
    return matches.iloc[0]["SeriesInstanceUID"]


def read_slide_image_comments(idc_client, sm_coll, slide_id, label, qlog):
    """Resolve slide_id -> IDC series -> base VOLUME instance -> range-read its
    header -> ImageComments string. Returns
    (status, image_comments_or_None, parsed_dict) where status is one of
    "no_idc_series", "no_volume_instance", "ok"."""
    series_uid = resolve_series_for_slide(sm_coll, slide_id, label)
    if series_uid is None:
        return "no_idc_series", None, {}
    smi_sub = idc_client.sm_instance_index
    smi_sub = smi_sub[smi_sub["SeriesInstanceUID"] == series_uid]
    base = pick_base_volume_instance(smi_sub)
    if base is None:
        return "no_volume_instance", None, {}
    idx = idc_client.index
    series_row = idx[idx["SeriesInstanceUID"] == series_uid].iloc[0]
    s3_url = series_row["series_aws_url"][:-1] + base["crdc_instance_uuid"] + ".dcm"
    https_url = to_https(s3_url)
    raw, _status, _size, _complete = fetch_header_bytes(https_url, series_uid, label, qlog)
    ic = extract_image_comments(raw)
    parsed = parse_aperio_string(ic) if ic else {}
    return "ok", ic, parsed


def _get_idc_client():
    """Lazy import + one client, fetching only the two indexes this script needs
    (sm_index for the ContainerIdentifier join, sm_instance_index for base-VOLUME
    resolution) -- same two indexes idc_step3_operator_header_probe.py fetches."""
    from idc_index import IDCClient
    c = IDCClient()
    c.fetch_index("sm_index")
    c.fetch_index("sm_instance_index")
    return c


# =============================================================================
# main
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--collection", default="cptac_ucec",
                     help="IDC collection_id (D4 item 3 explicitly authorises "
                          "cptac_ucec; default and only value this project uses)")
    ap.add_argument("--conf-slide-map", default=str(PINNED / "slide_type_map_ucec_c1.tsv"),
                     help="the 169 analysed Primary Tumor UCEC confirmatory slides")
    ap.add_argument("--disc-slide-map", default=str(SCRIPTS_DIR / "slide_type_map_ucec.tsv"),
                     help="discovery UCEC slide map (filtered to sample_type == "
                          "'Primary Tumor')")
    ap.add_argument("--acq-meta", default=str(SCRIPTS_DIR / "slide_acquisition_meta.tsv"),
                     help="discovery ground truth (aperio.User/aperio.Date), "
                          "filtered to cohort == 'ucec'")
    ap.add_argument("--case-ids", default=str(PINNED / "c1_case_ids_ucec.txt"),
                     help="the published 138-case UCEC confirmatory analysed set")
    ap.add_argument("--out-labels", default=str(PINNED / "c1_ucec_batch_labels.tsv"))
    ap.add_argument("--out-dir", default=str(REVIEW_OUT_DEFAULT))
    ap.add_argument("--limit", type=int, default=None,
                     help="debug only: cap the number of discovery/confirmatory "
                          "slides probed (full run uses all 393 discovery + 169 "
                          "confirmatory slides)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    qlog = QueryLog(out_dir / "query_log_idc.tsv")

    print("=== resolved paths ===")
    for k in ("collection", "conf_slide_map", "disc_slide_map", "acq_meta", "case_ids",
              "out_labels", "out_dir"):
        print(f"  {k:<16}: {getattr(args, k)}")
    print("=" * 74)

    # ---- load inputs ---------------------------------------------------------
    conf_map = pd.read_csv(args.conf_slide_map, sep="\t", dtype=str)
    if not (conf_map["sample_type"] == "Primary Tumor").all():
        sys.exit(f"FATAL: {args.conf_slide_map} has non-'Primary Tumor' rows -- "
                  f"wrong file, this must be the 169 analysed slides only.")
    if len(conf_map) != N_CONF_SLIDES_EXPECTED:
        sys.exit(f"FATAL: {args.conf_slide_map} has {len(conf_map)} rows, expected "
                  f"{N_CONF_SLIDES_EXPECTED}.")
    conf_map = conf_map.rename(columns={"slide_submitter_id": "slide_id"})

    case_ids = sorted({l.strip() for l in open(args.case_ids, encoding="utf-8") if l.strip()})
    if len(case_ids) != N_CASES_EXPECTED:
        sys.exit(f"FATAL: {args.case_ids} has {len(case_ids)} cases, expected "
                  f"{N_CASES_EXPECTED}.")
    conf_cases_from_map = set(conf_map["case_submitter_id"])
    if not conf_cases_from_map.issubset(set(case_ids)):
        extra = sorted(conf_cases_from_map - set(case_ids))
        sys.exit(f"FATAL: {args.conf_slide_map} names {len(extra)} cases not on "
                  f"the published 138-case list, e.g. {extra[:5]}.")

    disc_map = pd.read_csv(args.disc_slide_map, sep="\t", dtype=str)
    disc_pt = disc_map[disc_map["sample_type"] == "Primary Tumor"].rename(
        columns={"slide_submitter_id": "slide_id"})
    if len(disc_pt) != N_DISC_SLIDES_EXPECTED:
        print(f"WARNING: {args.disc_slide_map} has {len(disc_pt)} Primary Tumor "
              f"rows, expected {N_DISC_SLIDES_EXPECTED} -- proceeding on what is "
              f"actually there (this is a discovery-side sanity count, not a "
              f"frozen constant).")

    acq_meta = pd.read_csv(args.acq_meta, sep="\t")
    acq_ucec = acq_meta[acq_meta["cohort"] == "ucec"].set_index("slide_submitter_id")

    disc_ids = disc_pt["slide_id"].tolist()
    conf_ids = conf_map["slide_id"].tolist()
    if args.limit:
        disc_ids = disc_ids[:args.limit]
        conf_ids = conf_ids[:args.limit]
        print(f"--limit {args.limit}: probing {len(disc_ids)} discovery / "
              f"{len(conf_ids)} confirmatory slides only (DEBUG; not a real run)")

    # ---- IDC client + collection join -----------------------------------------
    c = _get_idc_client()
    sm_coll = build_idc_series_lookup(c, args.collection)
    print(f"IDC collection {args.collection!r}: {len(sm_coll)} SM series total")

    # ---- Step 1a: discovery validation ----------------------------------------
    print(f"\n=== Step 1: discovery validation ({len(disc_ids)} slides, "
          f"collection={args.collection}) ===")
    disc_rows = []
    n_checked = n_agree_op = n_agree_date = 0
    for i, sid in enumerate(disc_ids, 1):
        status, ic, parsed = read_slide_image_comments(c, sm_coll, sid, "discovery", qlog)
        op, date = parsed.get("User"), parsed.get("Date")
        truth_op = truth_date = None
        agree_op = agree_date = None
        if sid in acq_ucec.index:
            row = acq_ucec.loc[sid]
            truth_op = None if pd.isna(row["operator_uuid"]) else str(row["operator_uuid"])
            truth_date = None if pd.isna(row["scan_date"]) else str(row["scan_date"])
            if op is not None and truth_op is not None:
                n_checked += 1
                agree_op = (op == truth_op)
                agree_date = (date == truth_date) if (date is not None and truth_date is not None) else False
                n_agree_op += int(agree_op)
                n_agree_date += int(agree_date)
        disc_rows.append(dict(slide_id=sid, case=slide_to_case(sid), idc_status=status,
                               operator_uuid=op, scan_date=date,
                               truth_operator_uuid=truth_op, truth_scan_date=truth_date,
                               agree_operator=agree_op, agree_date=agree_date))
        if i % 25 == 0 or i == len(disc_ids):
            print(f"  [{i}/{len(disc_ids)}] discovery validated so far: "
                  f"{n_agree_op}/{n_checked} operator agree, "
                  f"{n_agree_date}/{n_checked} date agree")

    disc_df = pd.DataFrame(disc_rows)
    disc_out = out_dir / "idc_ucec_discovery_validation.tsv"
    disc_df.to_csv(disc_out, sep="\t", index=False)
    print(f"wrote {disc_out}")

    validated_operator = n_checked > 0 and n_agree_op == n_checked
    validated_date = n_checked > 0 and n_agree_date == n_checked
    print(f"\nDiscovery validation: n_checked={n_checked} "
          f"(slide has both an ImageComments value and a slide_acquisition_meta.tsv "
          f"ground-truth row)")
    print(f"  operator_uuid: {n_agree_op}/{n_checked} exact agreement -> "
          f"{'VALIDATED' if validated_operator else 'NOT VALIDATED'}")
    print(f"  scan_date    : {n_agree_date}/{n_checked} exact agreement -> "
          f"{'VALIDATED' if validated_date else 'NOT VALIDATED'}")

    if not validated_operator:
        print("\nNO VALIDATED operator attribute on UCEC discovery slides -- "
              "stopping per D4 item 3 step 1 / the rule's fallback clause. "
              "c1_ucec_batch_labels.tsv is NOT written.")
        summary = dict(
            collection=args.collection, validated_operator=False, validated_date=validated_date,
            n_discovery_checked=n_checked, n_agree_operator=n_agree_op, n_agree_date=n_agree_date,
            qualifies=False,
            reason="operator attribute (ImageComments User=) did not validate on "
                   "UCEC discovery slides",
            output_labels_path=None, sha256_output_labels=None,
            generated_at_utc=utcnow(),
        )
        with open(out_dir / "c1_ucec_d4_summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
        print(f"wrote {out_dir / 'c1_ucec_d4_summary.json'}")
        return

    # ---- Step 1b/2: read confirmatory headers, derive per-case labels ---------
    print(f"\n=== Step 2: confirmatory label extraction ({len(conf_ids)} slides, "
          f"collection={args.collection}) ===")
    conf_rows = []
    for i, sid in enumerate(conf_ids, 1):
        status, ic, parsed = read_slide_image_comments(c, sm_coll, sid, "confirmatory", qlog)
        op = parsed.get("User")
        date = parsed.get("Date")
        q = scan_quarter_from_date(date) if date else None
        note = "" if status == "ok" else f"idc_status={status}"
        conf_rows.append(dict(slide_id=sid, case=slide_to_case(sid), operator_uuid=op,
                               scan_date=date, scan_quarter=q, note=note))
        if i % 25 == 0 or i == len(conf_ids):
            n_with_op = sum(1 for r in conf_rows if r["operator_uuid"])
            print(f"  [{i}/{len(conf_ids)}] confirmatory headers read so far: "
                  f"{n_with_op}/{i} with an operator_uuid")

    conf_df = pd.DataFrame(conf_rows)
    conf_out = out_dir / "idc_ucec_confirmatory_per_slide_labels.tsv"
    conf_df.to_csv(conf_out, sep="\t", index=False)
    print(f"wrote {conf_out} ({len(conf_df)} rows, "
          f"{conf_df['operator_uuid'].notna().sum()} with operator_uuid, "
          f"{conf_df['scan_quarter'].notna().sum()} with scan_quarter)")

    op_rows_out, q_rows_out = [], []
    for case in sorted(case_ids):
        case_slides = conf_df[conf_df["case"] == case].to_dict("records")
        op_lab, op_purity, op_n = modal_case_label(case_slides, "operator_uuid")
        q_lab, q_purity, q_n = modal_case_label(case_slides, "scan_quarter")
        op_rows_out.append(dict(case=case, label=op_lab or "", purity=op_purity, n_used=op_n))
        q_rows_out.append(dict(case=case, label=q_lab or "", purity=q_purity, n_used=q_n))

    labels_df = pd.DataFrame({
        "case_submitter_id": [r["case"] for r in op_rows_out],
        "operator_uuid": [r["label"] for r in op_rows_out],
        "operator_purity": [round(r["purity"], 3) if r["purity"] is not None else "" for r in op_rows_out],
        "scan_quarter": [r["label"] for r in q_rows_out],
        "scan_quarter_purity": [round(r["purity"], 3) if r["purity"] is not None else "" for r in q_rows_out],
        "n_primary_slides_used": [r["n_used"] for r in op_rows_out],
    })
    out_labels_path = Path(args.out_labels)
    out_labels_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_labels_path, "w", encoding="utf-8", newline="\n") as f:
        labels_df.to_csv(f, sep="\t", index=False, lineterminator="\n")
    sha256_output_labels = hashlib.sha256(out_labels_path.read_bytes()).hexdigest()
    print(f"\nwrote {out_labels_path} ({len(labels_df)} case rows) -- schema "
          f"identical to pinned/c1_luad_batch_labels.tsv; sha256="
          f"{sha256_output_labels[:16]}...")

    # ---- D1 qualification, on the published 138-case analysed set -------------
    op_report = evaluate_d1_criteria(labels_df["operator_uuid"].tolist(), N_CASES_EXPECTED)
    q_report = evaluate_d1_criteria(labels_df["scan_quarter"].tolist(), N_CASES_EXPECTED)
    print("\n=== Operator axis vs D1 criteria (published 138-case analysed set) ===")
    print(json.dumps(op_report, indent=2))
    print("\n=== Scan-quarter axis vs D1 criteria (descriptive only; D4 item 3 "
          "never falls back to it) ===")
    print(json.dumps(q_report, indent=2))

    qualifies = bool(validated_operator and op_report["D1_all_four_hold"])
    print(f"\n{'QUALIFIES' if qualifies else 'DOES NOT QUALIFY'}: operator axis "
          f"validated={validated_operator}, D1_all_four_hold={op_report['D1_all_four_hold']}")
    if not qualifies:
        print("Per D4 item 3's fallback clause: \"Methods state that C1-UCEC was "
              "tested under the unrestricted null only and carries no batch "
              "qualifier.\" c1_ucec_stratified_reread.py must not compute a "
              "stratified AUC when it reads qualifies=false from this summary.")

    summary = dict(
        collection=args.collection,
        validated_operator=validated_operator, validated_date=validated_date,
        n_discovery_checked=n_checked, n_agree_operator=n_agree_op, n_agree_date=n_agree_date,
        operator_axis=op_report, scan_quarter_axis=q_report,
        qualifies=qualifies,
        output_labels_path=str(out_labels_path),
        sha256_output_labels=sha256_output_labels,
        discovery_validation_path=str(disc_out),
        confirmatory_per_slide_path=str(conf_out),
        n_cases=N_CASES_EXPECTED, n_confirmatory_slides=N_CONF_SLIDES_EXPECTED,
        limit=args.limit,
        generated_at_utc=utcnow(),
    )
    summary_path = out_dir / "c1_ucec_d4_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"\nwrote {summary_path}")


if __name__ == "__main__":
    main()
