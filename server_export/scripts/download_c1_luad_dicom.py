#!/usr/bin/env python3
"""C1-LUAD confirmatory cohort: download the 135 confirmatory primary-slide DICOM
series from IDC (collection cptac_luad) to local disk, ahead of extract_c1_dicom.py.

Authority: review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md section 0 D3 is complete
(signed 2026-10-01, OSF https://osf.io/2d6mb/, GitHub release
c1-luad-signed-rule-2026-10-01) -- confirmatory data may now be downloaded. This is
section 1.6 step 2 of the rule ("Extraction ... may run only after D3 registration").

The download list is FROZEN, not re-derived here: pinned/c1_luad_confirmatory_series.tsv
(135 rows, slide_id -> case -> SeriesInstanceUID), built from the already-signed,
hash-pinned pinned/c1_primary_slide_ids_luad_idc.txt and Step 0b's
idc_step2_confirmatory_pt_join.tsv join, and itself hash-verified to reproduce that
135-ID list exactly before being written. This script hash-checks it at start and
refuses to run on any other input.

extract_c1_dicom.py discovers series by globbing --dicom-root recursively for *.dcm
and grouping by parent directory ("one slide = one DICOM series directory"); it does
not care about the exact subdirectory layout beyond that. idc_index's
download_from_selection writes one subdirectory per series under downloadDir, which
satisfies this directly.

Usage:
  # smoke: download the first 2 series only, confirm the mechanism works
  python download_c1_luad_dicom.py --dicom-root <DEST> --limit 2

  # full run (about 53 GB, 135 series; safe to re-run -- each series is verified by
  # file count and byte size against the manifest before being (re-)fetched, so an
  # already-complete series from a prior run is skipped, not re-downloaded)
  python download_c1_luad_dicom.py --dicom-root <DEST>

One series at a time, not one bulk parallel call. 2026-10-01: a bulk
download_from_selection() over all 135 series died twice in a row (non-zero exit
from its internal s5cmd subprocess, each time within ~1-2 minutes, each time with
every series' small instances present but most series missing their large VOLUME
instances) -- consistent with a cross-border long-lived connection to the
idc-open-data S3 bucket (likely us-east-1) being reset partway through a large
transfer, not a quota/host-resource problem (dmesg showed no OOM-kill or signal
matching this process). Looping one series per download_from_selection() call keeps
each individual transfer short-lived, makes a reset cost only that one series'
retry rather than the whole batch, and lets an already-finished series from a prior
invocation be skipped outright without even attempting a network call.

PYTHONUTF8=1 is needed for idc_index's own import on Windows (Step 0b's own note);
harmless and unnecessary on Linux, kept for parity with idc_ucec_derive_batch_labels.py.
"""
import argparse
import hashlib
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PINNED = HERE / "pinned"
MANIFEST = PINNED / "c1_luad_confirmatory_series.tsv"
MANIFEST_SHA256 = "3279e0b57054f5348dfe965209a811e206f312e235105d7a8a437dd57b34a97f"
N_EXPECTED = 135


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log(msg):
    print(msg, flush=True)


MAX_ATTEMPTS = 3
RETRY_DELAY_S = 10


def scan_series(dest, uids, expected):
    """Per-series (n_files, n_bytes) for every uid found under dest, keyed by
    SeriesInstanceUID (parsed from the 'SM_<uid>' directory name idc_index writes)."""
    by_series = {}
    for p in dest.rglob("SM_*"):
        if not p.is_dir():
            continue
        uid = p.name[len("SM_"):]
        if uid not in expected:
            continue
        files = list(p.glob("*.dcm"))
        by_series[uid] = (len(files), sum(f.stat().st_size for f in files))
    return by_series


def is_complete(uid, by_series, expected):
    if uid not in by_series:
        return False
    n_files, n_bytes = by_series[uid]
    want_files = expected[uid]["instanceCount"]
    want_bytes = expected[uid]["series_size_MB"] * 1e6
    return n_files == want_files and n_bytes >= 0.9 * want_bytes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dicom-root", required=True,
                    help="destination directory (matches lsf_extract_c1_luad.sh's "
                         "--dicom-root; nothing under it needs to pre-exist)")
    ap.add_argument("--limit", type=int, default=0,
                    help="smoke-test mode: download only the first N series")
    args = ap.parse_args()

    got_hash = sha256(MANIFEST)
    if got_hash != MANIFEST_SHA256:
        sys.exit(f"FATAL: {MANIFEST} sha256 mismatch: got={got_hash} "
                 f"want={MANIFEST_SHA256} -- refusing to download against an "
                 f"unverified list.")
    log(f"hash OK: {MANIFEST.name} ({got_hash[:12]}...)")

    df = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    df["series_size_MB"] = df["series_size_MB"].astype(float)
    df["instanceCount"] = df["instanceCount"].astype(int)
    if len(df) != N_EXPECTED or not df["SeriesInstanceUID"].is_unique:
        sys.exit(f"FATAL: manifest has {len(df)} rows (want {N_EXPECTED}) or "
                 f"duplicate series -- this should be impossible given the hash "
                 f"check above; stopping rather than guessing.")
    expected = df.set_index("SeriesInstanceUID")[["instanceCount", "series_size_MB"]].to_dict("index")
    series = sorted(df["SeriesInstanceUID"].tolist())
    if args.limit:
        series = series[:args.limit]
        log(f"--limit {args.limit}: smoke-test mode, not a real run")
    log(f"{len(series)} series to fetch, collection cptac_luad, into "
        f"{args.dicom_root}, one series per call, up to {MAX_ATTEMPTS} attempts each "
        f"(2026-10-01: a single bulk call over all 135 died twice running mostly on "
        f"large-file transfers; see the module docstring)")

    from idc_index import IDCClient
    c = IDCClient()

    dest = Path(args.dicom_root)
    dest.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    already_done = failed = fetched = 0
    for i, uid in enumerate(series, 1):
        by_series = scan_series(dest, [uid], expected)
        if is_complete(uid, by_series, expected):
            already_done += 1
            if i % 20 == 0 or i == len(series):
                log(f"[{i}/{len(series)}] (skipped {already_done} already-complete "
                    f"so far)")
            continue

        ok = False
        for attempt in range(1, MAX_ATTEMPTS + 1):
            try:
                c.download_from_selection(downloadDir=str(dest), seriesInstanceUID=[uid],
                                          show_progress_bar=False)
            except Exception as e:
                log(f"[{i}/{len(series)}] {uid} attempt {attempt}/{MAX_ATTEMPTS}: "
                    f"exception {type(e).__name__}: {e}")
            by_series = scan_series(dest, [uid], expected)
            if is_complete(uid, by_series, expected):
                ok = True
                break
            n_files, n_bytes = by_series.get(uid, (0, 0))
            want_files = expected[uid]["instanceCount"]
            want_bytes = expected[uid]["series_size_MB"] * 1e6
            log(f"[{i}/{len(series)}] {uid} attempt {attempt}/{MAX_ATTEMPTS}: "
                f"{n_files}/{want_files} files, {n_bytes/1e6:.1f}/{want_bytes/1e6:.1f} MB "
                f"-- {'retrying' if attempt < MAX_ATTEMPTS else 'giving up'}")
            if attempt < MAX_ATTEMPTS:
                time.sleep(RETRY_DELAY_S)

        if ok:
            fetched += 1
            log(f"[{i}/{len(series)}] {uid} OK ({time.time()-t0:.0f}s elapsed total)")
        else:
            failed += 1
            log(f"[{i}/{len(series)}] {uid} FAILED after {MAX_ATTEMPTS} attempts")

    log(f"\nrun summary: {already_done} already complete, {fetched} fetched this run, "
        f"{failed} failed after {MAX_ATTEMPTS} attempts each, "
        f"{time.time()-t0:.0f}s elapsed")

    # Final whole-set verification, independent of the per-series bookkeeping above.
    by_series = scan_series(dest, series, expected)
    missing = [u for u in series if u not in by_series]
    incomplete = [u for u in series if u not in missing and not is_complete(u, by_series, expected)]
    total_got_bytes = sum(by_series.get(u, (0, 0))[1] for u in series)
    total_want_bytes = sum(expected[u]["series_size_MB"] for u in series) * 1e6
    log(f"final check: {len(series) - len(missing)}/{len(series)} series directories "
        f"found, {total_got_bytes / 1e9:.2f} GB / {total_want_bytes / 1e9:.2f} GB on disk")

    if missing:
        log(f"  {len(missing)} series MISSING entirely: {missing[:10]}"
            + (" ..." if len(missing) > 10 else ""))
    if incomplete:
        log(f"  {len(incomplete)} series INCOMPLETE: {incomplete[:10]}"
            + (" ..." if len(incomplete) > 10 else ""))

    if missing or incomplete:
        log("NOT COMPLETE -- re-run this same command to resume (already-complete "
            "series are skipped without a network call; only the "
            f"{len(missing)+len(incomplete)} remaining series will be attempted). Do "
            "NOT proceed to extract_c1_dicom.py until this check passes clean.")
        sys.exit(1)
    log(f"all {len(series)} series complete and verified (file count and byte size "
        f"both match the manifest). Next: run extract_c1_dicom.py --dry-run over "
        f"the full 135-ID confirmatory set (section 1.3/1.4 of the frozen rule) "
        f"before any real extraction.")


if __name__ == "__main__":
    main()
