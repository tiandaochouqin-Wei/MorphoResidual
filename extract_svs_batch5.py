#!/usr/bin/env python3
r"""
extract_svs_batch5.py -- per-slide scanner / site metadata from the raw SVS headers,
for all five cohorts in one run.

WHY THIS MATTERS NOW (it did not before)
The nuisance-block analysis found a site-level confounder in UCEC that the
manuscript's batch suite does not correct. The number of slides banked per patient
predicts UCEC's translation residual at out-of-fold R^2 = 0.30 -- three quarters of
what the 20 morphology components achieve -- and it is not disease burden
(uncorrelated with grade and stage; conditioning on both leaves it untouched) and
not a mean-pooling artefact (morphology predicts worse, not better, with more
slides). Screening the CPTAC UCEC clinical table put it on accrual country
(p = 2.7e-12) and OCT embedding (p = 3.7e-10): patients accrued in one country
contributed a median of six slides against two to three elsewhere and have
systematically higher residuals. Country and embedding medium sit outside the
operator / quarter / scanner axes the batch suite uses.

That makes the physical scanner identity in each SVS header worth extracting:
it is the one site proxy available per slide rather than per cohort, and GDC does
not expose CPTAC tissue source site (all None).

WHAT CHANGED FROM THE PILOT VERSION
  * five cohorts, not CCRCC only. The layout is not uniform: ccRCC predates the
    per-cancer directory scheme and lives at the ZW root with slides under
    WSI/raw, while the other four are at ZW/<cohort>/WSI/raw. Same branch as
    morpho_env.sh, and for the same reason.
  * does NOT read MORPHO_ROOT. morpho_env.sh sets that to the CURRENT COHORT's
    directory, and a shell that had sourced it once sent export_package.py hunting
    for all five cohorts inside ZW/pdac. Pass --root only if the root moves.
  * the hardcoded non-ccRCC exclusion list is now keyed by cohort, so it cannot
    silently drop cases from the other four.
  * per-slide output carries the cohort, and the per-case file is written once for
    all cohorts, which is what a site-adjusted re-run will consume.

USAGE
  cd /public/home/fjhui/ZW/scripts
  export LD_PRELOAD=$CONDA_PREFIX/lib/libopenslide.so.1   # openslide's find_library()
                                                          # ignores LD_LIBRARY_PATH
  bsub -q smp -n 4 -o svs_headers.out "python -u extract_svs_batch5.py"

  Header reads only -- no tile decoding -- so this is minutes, not hours. Add
  --dry-run to list what it would scan and exit.

  Then transfer wsi_scanner_map.tsv and wsi_scanner_by_case.tsv back.

NOTE ON IMPORT ORDER: openslide must be imported before torch/timm (libjpeg symbol
clash). This script imports nothing heavy, so it is safe as written -- do not add
a torch import above it.
"""
import argparse
import re
from collections import Counter
from pathlib import Path

import openslide          # keep first (libjpeg ordering)
import pandas as pd

ZW = "/public/home/fjhui/ZW"
COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
CASE_RE = re.compile(r"(C3[A-Z]-\d{5})")
SLIDE_RE = re.compile(r"(C3[A-Z]-\d{5}-\d+)")
# NOT .dcm. ZW/WSI/raw/ccrcc holds a DICOM mirror of the ccRCC download in which one
# slide is a UID-named directory of ~250 per-instance files; matching .dcm turned 390
# slides into 99098 "slides" and would have sent openslide at every instance in turn.
# The analysis slides are SVS in every cohort.
EXTS = (".svs", ".tif", ".tiff", ".ndpi")

# Cases excluded from the analysis cohort, per cohort. ccRCC's seven are the
# officially non-ccRCC histologies (Clark 2019); the others have no such list, and
# applying ccRCC's to them would silently drop real cases.
EXCLUDE = {
    "ccrcc": {"C3L-00359", "C3N-00313", "C3N-00435", "C3N-00492",
              "C3N-00832", "C3N-01175", "C3N-01180"},
}


def cohort_dirs(root):
    """Raw-slide directory and slide-type map per cohort, honouring ccRCC's legacy layout."""
    out = {}
    for c in COHORTS:
        if c == "ccrcc":
            raw = root / "WSI" / "raw"
            smap = root / "scripts" / "pinned" / "slide_type_map.tsv"
            alt = root / "scripts" / "slide_type_map.tsv"
        else:
            raw = root / c / "WSI" / "raw"
            smap = root / "scripts" / "pinned" / f"slide_type_map_{c}.tsv"
            alt = root / "scripts" / f"slide_type_map_{c}.tsv"
        out[c] = (raw, smap if smap.exists() else alt)
    return out


def props(path):
    try:
        s = openslide.OpenSlide(str(path))
    except Exception as e:                                    # noqa: BLE001
        return {"error": str(e)[:120]}
    p = dict(s.properties)
    s.close()
    return {
        "scanscope_id": p.get("aperio.ScanScope ID", ""),
        "vendor": p.get("openslide.vendor", ""),
        "appmag": p.get("aperio.AppMag", p.get("openslide.objective-power", "")),
        "mpp_x": p.get("openslide.mpp-x", ""),
        "mpp_y": p.get("openslide.mpp-y", ""),
        "scan_date": p.get("aperio.Date", ""),
        "make": p.get("tiff.Make", ""),
        "model": p.get("tiff.Model", ""),
        "user": p.get("aperio.User", ""),
        "error": "",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=ZW,
                    help="data root. Deliberately NOT read from MORPHO_ROOT, which "
                         "morpho_env.sh sets to the current cohort's own directory.")
    ap.add_argument("--out", default=None, help="output dir (default <root>/results)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    root = Path(a.root)
    if not root.is_dir():
        raise SystemExit(f"root not found: {root}")
    out_dir = Path(a.out) if a.out else root / "results"

    dirs = cohort_dirs(root)
    plan = {}
    for c, (raw, smap) in dirs.items():
        files = sorted(f for f in raw.rglob("*") if f.suffix.lower() in EXTS) if raw.is_dir() else []
        plan[c] = (raw, smap, files)
        print(f"  {c.upper():6s} {str(raw):46s} {len(files):5d} slides"
              f"{'' if smap.exists() else '   (no slide-type map: tumour filter skipped)'}")
        if a.dry_run and files:
            # which subtree each slide came from -- ccRCC's raw dir holds more than one
            # copy of the cohort, and a mirror counted twice would silently double rows
            sub = Counter((f.relative_to(raw).parts or ("<top>",))[0] if len(
                f.relative_to(raw).parts) > 1 else "<top>" for f in files)
            for k, v in sorted(sub.items(), key=lambda kv: -kv[1]):
                print(f"         {k:44s} {v:5d}")
            dup = [k for k, v in Counter(f.stem for f in files).items() if v > 1]
            if dup:
                print(f"         !! {len(dup)} slide names occur more than once "
                      f"(e.g. {dup[:3]}) -- two copies of the same cohort?")
    total = sum(len(f) for _, _, f in plan.values())
    print(f"\ntotal {total} raw WSI files")
    if not total:
        raise SystemExit("nothing found -- check --root")
    if a.dry_run:
        print("--dry-run: nothing written")
        return
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for c, (raw, smap, files) in plan.items():
        st = {}
        if smap.exists():
            m = pd.read_csv(smap, sep="\t")
            col = "slide_submitter_id" if "slide_submitter_id" in m.columns else m.columns[0]
            typ = "sample_type" if "sample_type" in m.columns else m.columns[-1]
            st = m.set_index(col)[typ].to_dict()
        for i, f in enumerate(files):
            cm, sm = CASE_RE.search(f.name), SLIDE_RE.search(f.name)
            slide = sm.group(1) if sm else f.stem
            rows.append({"cohort": c.upper(), "file": f.name,
                         "case_id": cm.group(1) if cm else "",
                         "slide_id": slide, "sample_type": st.get(slide, ""),
                         **props(f)})
            if (i + 1) % 200 == 0:
                print(f"    {c.upper()}: {i+1}/{len(files)}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "wsi_scanner_map.tsv", sep="\t", index=False)
    print(f"\nwrote per-slide -> {out_dir/'wsi_scanner_map.tsv'}  ({len(df)} slides)")

    ok = df[df["error"] == ""]
    bad = len(df) - len(ok)
    if bad:
        print(f"  {bad} slides could not be opened; first few: "
              f"{df[df['error']!=''][['cohort','file','error']].head(3).to_dict('records')}")

    print("\nScanScope ID by cohort (blank = header absent):")
    tab = ok.assign(sid=ok["scanscope_id"].replace("", "<blank>")).groupby(
        ["cohort", "sid"]).size().unstack(fill_value=0)
    print(tab.to_string())
    print("\nMPP-x by cohort (median, and how many distinct values):")
    for c, g in ok.groupby("cohort"):
        v = pd.to_numeric(g["mpp_x"], errors="coerce").dropna()
        print(f"  {c:6s} n={len(v):4d} median={v.median() if len(v) else float('nan'):.4f} "
              f"distinct={v.round(4).nunique()}")
    if ok["scanscope_id"].replace("", pd.NA).isna().all():
        print("\n!!! ScanScope ID blank for every slide: these SVS lack Aperio headers. "
              "Fall back to vendor/make/model/appmag as the coarse scanner proxy.")

    # per-case dominant scanner over tumour slides, matching the analysis cohort
    per_case = []
    for (c, case), g in ok[ok["case_id"] != ""].groupby(["cohort", "case_id"]):
        if case in EXCLUDE.get(c.lower(), set()):
            continue
        gg = g[g["sample_type"] == "Primary Tumor"] if (g["sample_type"] != "").any() else g
        if not len(gg):
            continue
        ids = [x for x in gg["scanscope_id"] if x] or [x for x in gg["vendor"] if x]
        dom = Counter(ids).most_common(1)
        mpp = pd.to_numeric(gg["mpp_x"], errors="coerce").dropna()
        per_case.append({"cohort": c, "case_id": case,
                         "scanner_id": dom[0][0] if dom else "",
                         "n_scanners": len(set(ids)), "n_slides": len(gg),
                         "mpp_x": float(mpp.median()) if len(mpp) else float("nan"),
                         "appmag": gg["appmag"].mode().iloc[0] if len(gg["appmag"].mode()) else "",
                         "user": gg["user"].mode().iloc[0] if len(gg["user"].mode()) else ""})
    pc = pd.DataFrame(per_case)
    pc.to_csv(out_dir / "wsi_scanner_by_case.tsv", sep="\t", index=False)
    print(f"\nwrote per-case -> {out_dir/'wsi_scanner_by_case.tsv'}  ({len(pc)} cases)")
    if len(pc):
        print("\nper-case scanner_id counts by cohort:")
        print(pc.assign(s=pc["scanner_id"].replace("", "<blank>")).groupby(
            ["cohort", "s"]).size().unstack(fill_value=0).to_string())
        print(f"\ncases whose tumour slides span >1 scanner: "
              f"{int((pc['n_scanners']>1).sum())}/{len(pc)}")


if __name__ == "__main__":
    main()
