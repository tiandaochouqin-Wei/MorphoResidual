#!/usr/bin/env python3
"""
C1-UCEC confirmatory cohort: per-slide UNI tile embeddings from the IDC DICOM-WSI copy.

Tiling (256 px, level 0, dense), the tissue filter, the encoder and its transform are
imported unchanged from extract_features.py, the script that produced the discovery
embeddings. Only two things differ, both forced by the container format:
  * one slide = one DICOM series directory (LABEL/OVERVIEW/THUMBNAIL + 2-4 VOLUME
    instances), so slides are enumerated per series directory, not per file;
  * openslide is handed the series' base VOLUME instance, and the output is named by
    that instance's ContainerIdentifier (the CPTAC slide barcode, e.g. C3L-00086-21),
    the naming residual_analysis.load_wsi_embeddings() parses.
extract_features.load_slide() would send .dcm to wsidicom; it is overridden so pixels
go through openslide, the same library that read the discovery SVS files.

  python extract_c1_dicom.py --dicom-root R --out-dir O --dry-run          # no GPU
  python extract_c1_dicom.py --dicom-root R --out-dir O --only C3L-00086-21
  python extract_c1_dicom.py --dicom-root R --out-dir O --shard 0 --nshards 2
"""
import openslide  # 必须先于 torch/timm 导入, 避免 libjpeg/libtiff 冲突
import argparse
import hashlib
import os
import re
import sys
import time
from collections import Counter
from pathlib import Path

import pydicom
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import extract_features as EF  # noqa: E402

LEVEL = 0
SLIDE_RE = re.compile(r"^(C3[A-Z]-\d{5})-\d+$")
MANIFEST_COLS = ["slide_id", "case", "status", "n_tiles", "seconds", "width", "height",
                 "mpp_x", "transfer_syntax", "n_volume", "series_uid", "base_file"]


def _openslide_loader(path):
    return ("openslide", openslide.OpenSlide(str(path)))


EF.load_slide = _openslide_loader


def resolve_series(series_dir):
    vols = []
    for f in sorted(series_dir.glob("*.dcm")):
        ds = pydicom.dcmread(str(f), stop_before_pixels=True)
        image_type = list(ds.get("ImageType", []))
        if len(image_type) > 2 and image_type[2] == "VOLUME":
            area = int(ds.TotalPixelMatrixColumns) * int(ds.TotalPixelMatrixRows)
            vols.append((area, f, ds))
    if not vols:
        raise ValueError(f"no VOLUME instance in {series_dir}")
    vols.sort(key=lambda t: t[0], reverse=True)
    if len(vols) > 1 and vols[0][0] == vols[1][0]:
        raise ValueError(f"two base-size VOLUME instances (concatenation?) in {series_dir}")
    _, base, ds = vols[0]
    slide_id = str(ds.get("ContainerIdentifier", "")).strip()
    case = str(ds.get("PatientID", "")).strip()
    m = SLIDE_RE.match(slide_id)
    if not m or m.group(1) != case:
        raise ValueError(f"ContainerIdentifier {slide_id!r} does not match PatientID {case!r} "
                         f"in {series_dir}")
    return {
        "slide_id": slide_id, "case": case, "base_file": base,
        "width": int(ds.TotalPixelMatrixColumns), "height": int(ds.TotalPixelMatrixRows),
        "transfer_syntax": str(ds.file_meta.TransferSyntaxUID), "n_volume": len(vols),
        "series_uid": str(ds.SeriesInstanceUID),
    }


def check_readable(rec):
    s = openslide.OpenSlide(str(rec["base_file"]))
    try:
        dims = s.level_dimensions[LEVEL]
        if dims != (rec["width"], rec["height"]):
            raise ValueError(f"openslide level-0 {dims} != DICOM matrix "
                             f"{(rec['width'], rec['height'])}")
        s.read_region((rec["width"] // 2, rec["height"] // 2), LEVEL, (EF.TILE_SIZE, EF.TILE_SIZE))
        return s.properties.get(openslide.PROPERTY_NAME_MPP_X, "")
    finally:
        s.close()


def append_manifest(path, row):
    new = not path.exists()
    with open(path, "a", encoding="utf-8") as f:
        if new:
            f.write("\t".join(MANIFEST_COLS) + "\n")
        f.write("\t".join(str(row.get(c, "")) for c in MANIFEST_COLS) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dicom-root", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--model", default="uni")
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--only", default="", help="comma-separated slide ids")
    ap.add_argument("--only-file", default="", help="file with one slide id per line")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve every series and test-read one region; no model, no GPU")
    args = ap.parse_args()

    ef_path = Path(EF.__file__).resolve()
    print(f"extract_features.py sha256={hashlib.sha256(ef_path.read_bytes()).hexdigest()} "
          f"TILE_SIZE={EF.TILE_SIZE} SAT={EF.TISSUE_SAT_THRESH} "
          f"MIN_FRAC={EF.MIN_TISSUE_FRACTION} level={LEVEL}", flush=True)

    series_dirs = sorted({p.parent for p in Path(args.dicom_root).rglob("*.dcm")})
    recs = sorted((resolve_series(d) for d in series_dirs), key=lambda r: r["slide_id"])
    dup = sorted(k for k, v in Counter(r["slide_id"] for r in recs).items() if v > 1)
    if dup:
        sys.exit(f"duplicate slide ids across series: {dup}")
    print(f"{len(recs)} series/slides across {len({r['case'] for r in recs})} cases", flush=True)

    wanted = {s.strip() for s in args.only.split(",") if s.strip()}
    if args.only_file:
        wanted |= {line.strip() for line in open(args.only_file, encoding="utf-8") if line.strip()}
    if wanted:
        missing = wanted - {r["slide_id"] for r in recs}
        if missing:
            sys.exit(f"--only ids not found: {sorted(missing)}")
        recs = [r for r in recs if r["slide_id"] in wanted]
    recs = recs[args.shard::args.nshards]
    print(f"this run: {len(recs)} slides (shard {args.shard}/{args.nshards})", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        bad = 0
        for r in recs:
            try:
                mpp = check_readable(r)
                print(f"  ok  {r['slide_id']}  {r['width']}x{r['height']}  mpp={mpp}  "
                      f"{r['base_file'].name}", flush=True)
            except Exception as e:
                bad += 1
                print(f"  BAD {r['slide_id']}: {e}", flush=True)
        print(f"dry-run: {len(recs) - bad} readable, {bad} unreadable")
        sys.exit(1 if bad else 0)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}", flush=True)
    model, tfm = EF.build_encoder(args.model, device)
    manifest = out_dir.parent / f"c1_extract_manifest_shard{args.shard}of{args.nshards}.tsv"

    ok = skip = fail = 0
    for i, r in enumerate(recs, 1):
        tag = f"[{i}/{len(recs)}] {r['slide_id']}"
        out_path = out_dir / f"{r['slide_id']}.pt"
        if out_path.exists():
            print(f"{tag} skip (done)", flush=True)
            skip += 1
            continue
        t0 = time.time()
        row = dict(r, base_file=str(r["base_file"]))
        try:
            row["mpp_x"] = check_readable(r)
            emb, coords = EF.embed_slide(r["base_file"], model, tfm, device,
                                         batch_size=args.batch_size, level=LEVEL)
            if emb is None:
                row.update(status="no_tissue", n_tiles=0)
                fail += 1
                print(f"{tag} no tissue tiles", flush=True)
            else:
                tmp = out_path.with_suffix(".pt.partial")
                torch.save({"embeddings": emb, "coords": coords,
                            "source": str(r["base_file"]), "series_uid": r["series_uid"]}, tmp)
                os.replace(tmp, out_path)
                row.update(status="ok", n_tiles=int(emb.shape[0]))
                ok += 1
                print(f"{tag} {emb.shape[0]} tiles in {time.time() - t0:.0f}s", flush=True)
        except Exception as e:
            row["status"] = f"error:{type(e).__name__}:{e}".replace("\t", " ").replace("\n", " ")
            fail += 1
            print(f"{tag} FAILED: {e}", flush=True)
        row["seconds"] = f"{time.time() - t0:.0f}"
        append_manifest(manifest, row)

    print(f"\ndone. ok={ok} skip={skip} fail={fail} out_dir={out_dir}", flush=True)
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
