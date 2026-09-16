#!/usr/bin/env python3
"""
pathomics.py — interpretable nuclear morphometry per slide (B-9, CPU).

Hand-crafted H&E features (color deconvolution + connected components; no deep
segmentation network, no GPU) to interpret WHAT morphology encodes the residual.
Same tiling and OUTPUT FORMAT as extract_phikon, but each tile's "embedding" is a
6-vector of morphometric features, so residual_analysis.load_wsi_embeddings
mean-pools them into a per-patient morphometry matrix -- just point
MORPHO_WSI_EMB_DIR at the output dir.

Features per tile (Ruifrok-Johnston H&E deconvolution via skimage.rgb2hed):
  0 nuclear_density   fraction of hematoxylin-positive pixels
  1 nucleus_count     connected nuclei per 1000 px
  2 mean_nuc_area     mean nucleus area (px)
  3 nuc_area_std      nucleus-area SD (pleomorphism)
  4 chromatin_od      mean hematoxylin optical density inside nuclei
  5 eosin_fraction    fraction of eosin-positive (cytoplasm) pixels

Run on CPU (smp). openslide imported before torch.
    python pathomics.py --raw-dir <.../WSI/raw> --out-dir <.../WSI/pathomics/...> --stride 1024
"""
import argparse, sys, warnings
from pathlib import Path
import numpy as np
import openslide  # before torch
import torch
from skimage.color import rgb2hed
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops

warnings.filterwarnings("ignore")

TILE_SIZE = 256
TISSUE_SAT_THRESH = 15
MIN_TISSUE_FRACTION = 0.35
FEAT_NAMES = ["nuclear_density", "nucleus_count", "mean_nuc_area",
              "nuc_area_std", "chromatin_od", "eosin_fraction"]


def slide_dims(s, level=0):
    return s.level_dimensions[level]


def read_region(s, x, y, size, level=0):
    return s.read_region((x, y), level, (size, size)).convert("RGB")


def is_tissue(t):
    hsv = np.array(t.convert("HSV"))
    return (hsv[:, :, 1] > TISSUE_SAT_THRESH).mean() > MIN_TISSUE_FRACTION


def iter_tissue_tiles(s, tile=TILE_SIZE, level=0, stride=None):
    stride = stride or tile
    w, h = slide_dims(s, level)
    for y in range(0, h - tile, stride):
        for x in range(0, w - tile, stride):
            t = read_region(s, x, y, tile, level)
            if is_tissue(t):
                yield x, y, t


def _otsu(a):
    a = np.asarray(a, float)
    if a.std() < 1e-6:
        return a.mean()
    try:
        return threshold_otsu(a)
    except Exception:
        return a.mean()


def tile_features(tile_rgb):
    rgb = np.asarray(tile_rgb, dtype=float) / 255.0
    hed = rgb2hed(rgb)
    H, E = hed[:, :, 0], hed[:, :, 1]
    mask = H > _otsu(H)
    npx = H.size
    lbl = label(mask)
    areas = np.array([r.area for r in regionprops(lbl)]) if lbl.max() > 0 else np.array([])
    return [
        float(mask.mean()),                                   # nuclear_density
        float(len(areas)) / (npx / 1000.0),                   # nucleus_count per 1000px
        float(areas.mean()) if len(areas) else 0.0,           # mean_nuc_area
        float(areas.std()) if len(areas) > 1 else 0.0,        # nuc_area_std (pleomorphism)
        float(H[mask].mean()) if mask.any() else 0.0,         # chromatin_od
        float((E > _otsu(E)).mean()),                         # eosin_fraction
    ]


def slide_morphometry(path, level=0, stride=None):
    slide = openslide.OpenSlide(str(path))
    feats, coords = [], []
    for x, y, tile in iter_tissue_tiles(slide, level=level, stride=stride):
        feats.append(tile_features(tile))
        coords.append((x, y))
    if not feats:
        return None, None
    return torch.tensor(feats, dtype=torch.float32), torch.tensor(coords, dtype=torch.long)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--level", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1024,
                    help="tile step in px (coarse is fine for morphometric means)")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    raw_dir, out_dir = Path(args.raw_dir), Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    wsi = sorted(p for p in raw_dir.rglob("*")
                 if p.suffix.lower() in (".svs", ".tif", ".tiff", ".ndpi"))
    if args.limit:
        wsi = wsi[:args.limit]
    print(f"found {len(wsi)} slides under {raw_dir}")

    ok = fail = 0
    for i, path in enumerate(wsi, 1):
        out_path = out_dir / f"{path.stem}.pt"
        if out_path.exists():
            print(f"[{i}/{len(wsi)}] skip (done): {path.stem}")
            ok += 1
            continue
        print(f"[{i}/{len(wsi)}] morphometry {path.stem} ...")
        try:
            feats, coords = slide_morphometry(path, args.level, args.stride or None)
            if feats is None:
                print(f"  WARNING: no tissue tiles in {path.stem}, skipping")
                fail += 1
                continue
            torch.save({"embeddings": feats, "coords": coords, "source": str(path),
                        "feat_names": FEAT_NAMES}, out_path)
            print(f"  saved {feats.shape[0]} tiles x {feats.shape[1]} feats -> {out_path}")
            ok += 1
        except Exception as e:
            print(f"  FAILED on {path.stem}: {e}")
            fail += 1
    print(f"\ndone. ok={ok} fail={fail} out_dir={out_dir}")
    sys.exit(1 if fail and not ok else 0)


if __name__ == "__main__":
    main()
