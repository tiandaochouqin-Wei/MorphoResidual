#!/usr/bin/env python3
"""
Build an authoritative per-slide sample_type map for a CPTAC-3 cohort from the
GDC biospecimen slide entities.

Why this exists: build_slide_map.py stamps EVERY slide "Primary Tumor". That is
defensible for GBM (brain carries no adjacent-normal sections) but false for UCEC
and PDAC, which do have Solid Tissue Normal sections already embedded on disk
(e.g. UCEC C3L-00006-26, PDAC C3L-00017-24). Mean-pooling those into the
per-patient WSI vector would recreate exactly the tumor/normal contamination that
residual_analysis.py documents as a fixed bug for ccRCC.

Why the slide entity and nothing else: GDC no longer hosts CPTAC-3 "Slide Image"
FILES (the data_type facet returns 0), so the file endpoint is useless here. The
tempting shortcut -- slide CASE-2N belongs to sample CASE-0N -- was tested and
REJECTED: for LUAD C3N-02000 the registered entity puts slide -22 under sample
-03 (Solid Tissue Normal) and slide -23 under sample -02 (Primary Tumor), i.e.
the +20 offset is violated and would have flipped both labels. C3L-02365-23 fails
the same way. Only samples.portions.slides is authoritative.

Consequence: slides with no registered slide entity get NO row and are dropped by
residual_analysis.py's "unmatched -> exclude to stay conservative" branch. That is
the same treatment LUAD already received (247 of 504 slides mapped), so the
cohorts stay methodologically consistent.

Usage:
  python build_slide_map_gdc.py validate            # must pass before building
  python build_slide_map_gdc.py build ucec pdac
"""
import glob
import json
import os
import re
import sys
import time
import urllib.request

CASES_API = "https://api.gdc.cancer.gov/cases"
SCRIPTS = "/public/home/fjhui/ZW/scripts"
ROOT = "/public/home/fjhui/ZW"
# Where to WRITE the map. The scripts/ directory is shared with other people on
# this account, and the legacy all-Primary-Tumor build_slide_map.py writes the
# same filenames -- on 2026-09-01 a concurrent run of it replaced the UCEC and
# PDAC maps mid-session, silently re-admitting 15 and 79 adjacent-normal slides.
# Point MORPHO_MAP_DIR at a private directory to keep this pipeline's inputs
# immune to that, without overwriting anyone else's files.
OUTDIR = os.environ.get("MORPHO_MAP_DIR", SCRIPTS)

COHORTS = {
    "ccrcc": {"emb": f"{ROOT}/WSI/emb/ccrcc_real", "ref": f"{SCRIPTS}/slide_type_map.tsv"},
    "luad":  {"emb": f"{ROOT}/luad/WSI/emb",       "ref": f"{SCRIPTS}/slide_type_map_luad.tsv"},
    "gbm":   {"emb": f"{ROOT}/gbm/WSI/emb",        "ref": f"{SCRIPTS}/slide_type_map_gbm.tsv"},
    "ucec":  {"emb": f"{ROOT}/ucec/WSI/emb",       "ref": None},
    "pdac":  {"emb": f"{ROOT}/pdac/WSI/emb",       "ref": None},
}

SLIDE_RE = re.compile(r"^(C3[A-Z]-\d{5})-(\d+)$")


def local_slides(cancer):
    out = {}
    bad = 0
    for p in sorted(glob.glob(COHORTS[cancer]["emb"] + "/*.pt")):
        sid = os.path.basename(p)[:-3]
        m = SLIDE_RE.match(sid)
        if not m:
            bad += 1
            continue
        out[sid] = m.group(1)
    if bad:
        print(f"  note: {bad} embedding files have non-CPTAC slide ids, skipped")
    return out


def fetch_slide_entities(case_ids, chunk=60):
    """slide_submitter_id -> (sample_type, case, section_location), authoritative."""
    slide_type = {}
    conflicts = []
    normal_samples = {}   # case -> set of non-tumor sample types present (for reporting)
    case_ids = sorted(case_ids)
    for i in range(0, len(case_ids), chunk):
        batch = case_ids[i:i + chunk]
        payload = {
            "filters": {"op": "in", "content": {"field": "submitter_id", "value": batch}},
            "expand": "samples,samples.portions,samples.portions.slides",
            "fields": "submitter_id",
            "size": str(len(batch) + 10),
        }
        for attempt in range(1, 6):
            try:
                req = urllib.request.Request(
                    CASES_API, data=json.dumps(payload).encode(),
                    headers={"Content-Type": "application/json"})
                with urllib.request.urlopen(req, timeout=180) as r:
                    d = json.load(r)
                break
            except Exception as e:
                print(f"    GDC attempt {attempt}/5 failed: {e}")
                time.sleep(4 * attempt)
        else:
            sys.exit("GDC unreachable -- aborting rather than writing a partial map")

        for h in d["data"]["hits"]:
            case = h["submitter_id"]
            for s in h.get("samples", []):
                stype = s.get("sample_type")
                if stype and "Tumor" not in stype:
                    normal_samples.setdefault(case, set()).add(stype)
                for p in s.get("portions", []):
                    for sl in p.get("slides", []):
                        sid = sl.get("submitter_id")
                        if not sid:
                            continue
                        rec = (stype, case, sl.get("section_location"))
                        if sid in slide_type and slide_type[sid][0] != stype:
                            conflicts.append((sid, slide_type[sid][0], stype))
                        slide_type[sid] = rec
        print(f"  fetched {min(i + chunk, len(case_ids))}/{len(case_ids)} cases")
    if conflicts:
        print(f"  WARNING: {len(conflicts)} slides map to conflicting sample types")
        for c in conflicts[:10]:
            print(f"      {c}")
    return slide_type, normal_samples


def survey(cancer):
    print(f"\n=== {cancer} ===")
    slides = local_slides(cancer)
    cases = set(slides.values())
    print(f"  {len(slides)} slides on disk across {len(cases)} cases")
    entities, normals = fetch_slide_entities(cases)
    mapped = {sid: entities[sid] for sid in slides if sid in entities}
    tally = {}
    for stype, _c, _loc in mapped.values():
        tally[stype] = tally.get(stype, 0) + 1
    print(f"  slide entities registered in GDC: {len(entities)}")
    print(f"  of our {len(slides)} on-disk slides, {len(mapped)} matched, "
          f"{len(slides) - len(mapped)} unmatched (will be dropped)")
    print(f"  tally: {tally}")
    cases_with_normal_tissue = {c for c, t in normals.items() if "Solid Tissue Normal" in t}
    print(f"  cases having a Solid Tissue Normal sample: "
          f"{len(cases_with_normal_tissue)}/{len(cases)}")
    kept_cases = {c for _s, (_t, c, _l) in
                  [(k, v) for k, v in mapped.items() if v[0] == "Primary Tumor"]}
    print(f"  cases retaining >=1 Primary Tumor slide: {len(kept_cases)}")
    return slides, mapped


def validate_against_ref(cancer, mapped):
    ref_path = COHORTS[cancer]["ref"]
    if not ref_path or not os.path.exists(ref_path):
        return None
    ref = {}
    with open(ref_path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                ref[parts[0]] = parts[1]
    common = set(ref) & set(mapped)
    mism = [(s, ref[s], mapped[s][0]) for s in sorted(common) if ref[s] != mapped[s][0]]
    print(f"  ref {os.path.basename(ref_path)}: {len(ref)} rows | overlap {len(common)} "
          f"| mismatches {len(mism)} | only-in-ref {len(set(ref) - set(mapped))} "
          f"| only-in-new {len(set(mapped) - set(ref))}")
    for s, a, b in mism[:15]:
        print(f"      MISMATCH {s}: ref={a!r} new={b!r}")
    return len(mism) == 0 and len(common) > 0


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"

    if mode == "validate":
        ok_all = True
        for cancer in ["ccrcc", "luad"]:
            _slides, mapped = survey(cancer)
            ok = validate_against_ref(cancer, mapped)
            print(f"  [{cancer}] entity map {'REPRODUCES' if ok else 'CONFLICTS WITH'} reference")
            ok_all = ok_all and bool(ok)
        for cancer in ["ucec", "pdac"]:
            survey(cancer)
        print(f"\nOVERALL: {'VALIDATED' if ok_all else 'FAILED -- do not build'}")
        sys.exit(0 if ok_all else 1)

    if mode == "build":
        for cancer in sys.argv[2:]:
            _slides, mapped = survey(cancer)
            os.makedirs(OUTDIR, exist_ok=True)
            out = f"{OUTDIR}/slide_type_map_{cancer}.tsv"
            with open(out, "w") as f:
                f.write("slide_submitter_id\tsample_type\tcase_submitter_id\tsection_location\n")
                for sid in sorted(mapped):
                    stype, case, loc = mapped[sid]
                    f.write(f"{sid}\t{stype}\t{case}\t{loc if loc else ''}\n")
            print(f"  wrote {len(mapped)} rows -> {out}")
        sys.exit(0)

    sys.exit(f"unknown mode {mode!r}")


if __name__ == "__main__":
    main()
