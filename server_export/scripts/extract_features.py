#!/usr/bin/env python3
"""
Tile CPTAC-CCRCC H&E WSIs and extract per-tile embeddings with a frozen
pathology foundation model, for the beyond-mRNA residual pipeline.

Auto-detects .svs (openslide) vs .dcm (wsidicom) since it's not certain which
format NBIA Data Retriever will hand back for this collection -- both readers
are tried so the script doesn't silently fail on whichever format shows up.

Output per slide: WSI/emb/ccrcc/<slide_id>.pt containing
    {"embeddings": FloatTensor[N_tiles, D], "coords": LongTensor[N_tiles, 2]}
This is the per-slide bag that plan.md's ABMIL aggregation consumes.

Prereq: a HuggingFace token with accepted license for the chosen FM
(default: MahmoodLab/UNI). Set HF_TOKEN in the environment before running.
Mechanical pipeline -- fine to run on Sonnet.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import openslide  # 必须先于 torch/timm 导入, 避免 libjpeg/libtiff 冲突
import torch
from PIL import Image

TILE_SIZE = 256
TISSUE_SAT_THRESH = 15   # HSV saturation threshold; below this = background/glass
MIN_TISSUE_FRACTION = 0.35


def load_slide(path):
    suffix = path.suffix.lower()
    if suffix in (".svs", ".tif", ".tiff", ".ndpi"):
        import openslide
        return ("openslide", openslide.OpenSlide(str(path)))
    if suffix == ".dcm":
        from wsidicom import WsiDicom
        return ("wsidicom", WsiDicom.open(str(path)))
    raise ValueError(f"unrecognized WSI format: {path}")


def slide_dims(kind, slide, level=0):
    if kind == "openslide":
        return slide.level_dimensions[level]
    return (slide.size.width, slide.size.height)


def read_region(kind, slide, x, y, size, level=0):
    if kind == "openslide":
        return slide.read_region((x, y), level, (size, size)).convert("RGB")
    return slide.read_region(location=(x, y), level=level, size=(size, size)).convert("RGB")


def is_tissue(tile_rgb):
    hsv = np.array(tile_rgb.convert("HSV"))
    sat = hsv[:, :, 1]
    return (sat > TISSUE_SAT_THRESH).mean() > MIN_TISSUE_FRACTION


def iter_tissue_tiles(kind, slide, tile_size=TILE_SIZE, level=0, stride=None):
    stride = stride or tile_size
    w, h = slide_dims(kind, slide, level)
    for y in range(0, h - tile_size, stride):
        for x in range(0, w - tile_size, stride):
            tile = read_region(kind, slide, x, y, tile_size, level)
            if is_tissue(tile):
                yield x, y, tile


def build_encoder(model_name, device):
    if model_name == "uni":
        import timm
        token = os.environ.get("HF_TOKEN")
        if not token:
            sys.exit("HF_TOKEN not set. Request access to MahmoodLab/UNI on "
                      "HuggingFace, generate a token, and `export HF_TOKEN=...`.")
        model = timm.create_model(
            "hf-hub:MahmoodLab/UNI", pretrained=True, init_values=1e-5,
            dynamic_img_size=True,
        )
        from torchvision import transforms
        tfm = transforms.Compose([
            transforms.Resize(224),
            transforms.ToTensor(),
            transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
        ])
    else:
        sys.exit(f"unknown --model {model_name}; add its loader here before use")
    model.eval().to(device)
    return model, tfm


@torch.no_grad()
def embed_slide(slide_path, model, tfm, device, batch_size=64, level=0):
    kind, slide = load_slide(slide_path)
    tiles, coords = [], []
    embeddings = []

    def flush():
        if not tiles:
            return
        batch = torch.stack([tfm(t) for t in tiles]).to(device)
        out = model(batch).cpu()
        embeddings.append(out)
        tiles.clear()

    for x, y, tile in iter_tissue_tiles(kind, slide, level=level):
        tiles.append(tile)
        coords.append((x, y))
        if len(tiles) >= batch_size:
            flush()
    flush()

    if kind == "openslide":
        slide.close()
    else:
        slide.close()

    if not embeddings:
        return None, None
    return torch.cat(embeddings, dim=0), torch.tensor(coords, dtype=torch.long)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default=os.environ.get("MORPHO_WSI_DIR", "./WSI") + "/raw/ccrcc_real")
    ap.add_argument("--out-dir", default=os.environ.get("MORPHO_WSI_DIR", "./WSI") + "/emb/ccrcc_real")
    ap.add_argument("--model", default="uni")
    ap.add_argument("--level", type=int, default=0)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--limit", type=int, default=0,
                     help="only process the first N slides (0 = all) -- use for a smoke test "
                          "before committing to the full run")
    ap.add_argument("--delete-raw", action="store_true",
                     help="delete the source WSI after successful embedding "
                          "(off by default -- opt in once you trust the pipeline)")
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}")
    model, tfm = build_encoder(args.model, device)

    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    wsi_paths = sorted(
        p for p in raw_dir.rglob("*")
        if p.suffix.lower() in (".svs", ".tif", ".tiff", ".ndpi", ".dcm")
    )
    if args.limit:
        wsi_paths = wsi_paths[:args.limit]
    print(f"found {len(wsi_paths)} slide files under {raw_dir}"
          f"{f' (limited to first {args.limit})' if args.limit else ''}")

    ok, fail = 0, 0
    for i, path in enumerate(wsi_paths, 1):
        slide_id = path.stem
        out_path = out_dir / f"{slide_id}.pt"
        if out_path.exists():
            print(f"[{i}/{len(wsi_paths)}] skip (done): {slide_id}")
            ok += 1
            continue
        print(f"[{i}/{len(wsi_paths)}] embedding {slide_id} ...")
        try:
            emb, coords = embed_slide(path, model, tfm, device,
                                       batch_size=args.batch_size, level=args.level)
            if emb is None:
                print(f"  WARNING: no tissue tiles found in {slide_id}, skipping")
                fail += 1
                continue
            torch.save({"embeddings": emb, "coords": coords, "source": str(path)}, out_path)
            print(f"  saved {emb.shape[0]} tiles -> {out_path}")
            ok += 1
            if args.delete_raw:
                path.unlink()
        except Exception as e:
            print(f"  FAILED on {slide_id}: {e}")
            fail += 1

    print(f"\ndone. ok={ok} fail={fail} out_dir={out_dir}")
    sys.exit(1 if fail and not ok else 0)


if __name__ == "__main__":
    main()
