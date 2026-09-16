#!/usr/bin/env python3
"""
Census the acquisition metadata of every CPTAC slide we hold, from the SVS
headers only (no image decoding).

Why: the reviewer's fatal objection to C1 is that the WSI embedding might be
reading a scanner/centre fingerprint rather than morphology, and the project
notes recorded the metadata situation as an unsolved blocker -- GDC returns
tissue_source_site = _missing for all 1866 CPTAC-3 cases (re-verified), and
"scanner info: try extracting from SVS/DICOM headers" was marked never tried.

It works. Aperio SVS files carry an ImageDescription of the form
  Aperio Image Library v12.0.11 | AppMag = 20 | StripeWidth = 2032 |
  ScanScope ID = SS1553 | Filename = 102077 | Date = 11/14/16 | Time = ...
so ScanScope ID (a physical scanner serial), the scanning User UUID, the scan
date/time and the scan serial number are all recoverable per slide.

A 30-slide spot check found 29/30 on a single scanner (SS1553) across all five
cancers. If that holds over the full ~2100 slides it does something better than
any correction could: it makes the multi-scanner confound physically impossible,
because there is only one scanner. This script exists to establish that number
honestly over the whole corpus rather than from a spot check -- so it must
account for EVERY slide, and report failures loudly instead of dropping them.

Output: scripts/slide_acquisition_meta.tsv
"""
import csv
import glob
import multiprocessing as mp
import os
import re
import sys

# openslide first: this codebase has been bitten before by timm/torchvision
# grabbing libjpeg/libtiff ahead of it, which surfaces as the misleading
# "Couldn't locate OpenSlide shared library". Nothing else heavy is imported
# here, but keep the ordering habit.
import openslide

ROOT = "/public/home/fjhui/ZW"
OUT = f"{ROOT}/scripts/slide_acquisition_meta.tsv"

COHORT_RAW = {
    "ccrcc": f"{ROOT}/WSI/raw",
    "luad": f"{ROOT}/luad/WSI/raw",
    "ucec": f"{ROOT}/ucec/WSI/raw",
    "gbm": f"{ROOT}/gbm/WSI/raw",
    "pdac": f"{ROOT}/pdac/WSI/raw",
}

SLIDE_RE = re.compile(r"^(C3[A-Z]-\d{5})-(\d+)$")
LIB_RE = re.compile(r"Aperio Image Library (v[\d.]+)")

FIELDS = ["slide_submitter_id", "case_submitter_id", "cohort", "scanner_id",
          "operator_uuid", "scan_date", "scan_time", "mpp", "appmag",
          "scan_serial", "dsr_id", "image_lib_version", "vendor", "status"]


def one(args):
    cohort, path = args
    stem = os.path.basename(path)[:-4]
    m = SLIDE_RE.match(stem)
    case = m.group(1) if m else ""
    row = {f: "" for f in FIELDS}
    row.update(slide_submitter_id=stem, case_submitter_id=case, cohort=cohort)
    try:
        sl = openslide.OpenSlide(path)
        p = dict(sl.properties)
        sl.close()
    except Exception as e:
        row["status"] = f"OPEN_FAILED:{type(e).__name__}"
        return row
    desc = p.get("tiff.ImageDescription", "") or ""
    lib = LIB_RE.search(desc)
    row.update(
        scanner_id=p.get("aperio.ScanScope ID", "") or "",
        operator_uuid=p.get("aperio.User", "") or "",
        scan_date=p.get("aperio.Date", "") or "",
        scan_time=p.get("aperio.Time", "") or "",
        mpp=p.get("aperio.MPP", "") or "",
        appmag=p.get("aperio.AppMag", "") or "",
        scan_serial=p.get("aperio.Filename", "") or "",
        dsr_id=p.get("aperio.DSR ID", "") or "",
        image_lib_version=lib.group(1) if lib else "",
        vendor=p.get("openslide.vendor", "") or "",
        status="ok" if p.get("aperio.ScanScope ID") else "ok_no_scanner_field",
    )
    return row


def quarter(datestr):
    # aperio dates are MM/DD/YY
    try:
        mm, _dd, yy = datestr.split("/")
        return f"20{yy}Q{(int(mm) - 1) // 3 + 1}"
    except Exception:
        return ""


def tally(rows, key, limit=25):
    t = {}
    for r in rows:
        t[r[key]] = t.get(r[key], 0) + 1
    items = sorted(t.items(), key=lambda kv: -kv[1])
    print(f"\n--- {key}: {len(t)} distinct values ---")
    for k, v in items[:limit]:
        print(f"   {v:>6}  {k!r}")
    if len(items) > limit:
        print(f"   ... {len(items) - limit} more")
    return t


def main():
    tasks = []
    expected = {}
    for cohort, raw in COHORT_RAW.items():
        files = sorted(glob.glob(raw + "/**/*.svs", recursive=True))
        expected[cohort] = len(files)
        tasks.extend((cohort, f) for f in files)
    print("SVS files found per cohort:", expected)
    print(f"total {len(tasks)} slides to read\n")

    n_workers = int(os.environ.get("META_WORKERS", "8"))
    rows = []
    with mp.Pool(n_workers) as pool:
        for i, r in enumerate(pool.imap_unordered(one, tasks, chunksize=4), 1):
            rows.append(r)
            if i % 200 == 0:
                print(f"  {i}/{len(tasks)} headers read...", flush=True)

    for r in rows:
        r["scan_quarter"] = quarter(r["scan_date"])

    fields = FIELDS + ["scan_quarter"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["cohort"], r["slide_submitter_id"])):
            w.writerow(r)
    print(f"\nwrote {len(rows)} rows -> {OUT}")

    # ---- accounting: every slide must be accounted for ----
    got = {}
    for r in rows:
        got[r["cohort"]] = got.get(r["cohort"], 0) + 1
    print("\n=== ACCOUNTING ===")
    ok = True
    for c in COHORT_RAW:
        if expected[c] != got.get(c, 0):
            ok = False
            print(f"  MISMATCH {c}: expected {expected[c]} got {got.get(c, 0)}")
        else:
            print(f"  {c}: {got.get(c, 0)}/{expected[c]} OK")
    failed = [r for r in rows if r["status"].startswith("OPEN_FAILED")]
    print(f"  open failures: {len(failed)}")
    for r in failed[:20]:
        print(f"      {r['cohort']} {r['slide_submitter_id']} {r['status']}")
    noscan = [r for r in rows if r["status"] == "ok_no_scanner_field"]
    print(f"  parsed but no ScanScope ID: {len(noscan)}")
    for r in noscan[:20]:
        print(f"      {r['cohort']} {r['slide_submitter_id']}")

    # ---- the decisive numbers ----
    print("\n=== SCANNER CONCENTRATION (the question this script exists to answer) ===")
    st = tally(rows, "scanner_id")
    usable = sum(v for k, v in st.items() if k)
    if usable:
        top, topn = max(((k, v) for k, v in st.items() if k), key=lambda kv: kv[1])
        print(f"\n  dominant scanner {top!r}: {topn}/{usable} = {topn / usable:.2%} "
              f"of slides with a scanner field")
    print("\n  per-cohort scanner breakdown:")
    for c in COHORT_RAW:
        sub = [r for r in rows if r["cohort"] == c]
        t = {}
        for r in sub:
            t[r["scanner_id"]] = t.get(r["scanner_id"], 0) + 1
        print(f"    {c}: {dict(sorted(t.items(), key=lambda kv: -kv[1]))}")

    tally(rows, "mpp")
    tally(rows, "appmag")
    tally(rows, "vendor")
    tally(rows, "image_lib_version")
    tally(rows, "scan_quarter", limit=40)
    tally(rows, "operator_uuid", limit=15)

    print(f"\nOVERALL ACCOUNTING: {'OK' if ok and not failed else 'INCOMPLETE - investigate'}")
    sys.exit(0 if ok and not failed else 1)


if __name__ == "__main__":
    main()
