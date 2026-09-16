#!/usr/bin/env python3
"""
extract_phikon.py — second foundation-model robustness for MorphoResidual (B-8).

Self-contained clone of extract_features.py with the encoder swapped to
owkin/phikon (ViT-B, 768-d). Same tile size, same HSV tissue filter, same output
format {embeddings: FloatTensor[N,768], coords: LongTensor[N,2], source}, so
residual_analysis.py needs ZERO changes -- just point MORPHO_WSI_EMB_DIR here.

Run on GPU (gpu02 / LSF). openslide is imported BEFORE torch (libjpeg/libtiff clash).
Needs internet the first time to download owkin/phikon (run once on a node with net,
or pre-cache ~/.cache/huggingface).

    export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
    python extract_phikon.py --raw-dir <.../WSI/raw> --out-dir <.../WSI/emb_phikon/...>
"""
import argparse, sys
from pathlib import Path
import numpy as np
import openslide  # MUST precede torch/torchvision (libjpeg/libtiff conflict)
import torch
import torch.nn as nn
from transformers import AutoModel
from torchvision import transforms

TILE_SIZE = 256
TISSUE_SAT_THRESH = 15
MIN_TISSUE_FRACTION = 0.35


def slide_dims(slide, level=0):
    return slide.level_dimensions[level]


def read_region(slide, x, y, size, level=0):
    return slide.read_region((x, y), level, (size, size)).convert("RGB")


def is_tissue(tile_rgb):
    hsv = np.array(tile_rgb.convert("HSV"))
    return (hsv[:, :, 1] > TISSUE_SAT_THRESH).mean() > MIN_TISSUE_FRACTION


def iter_tissue_tiles(slide, tile_size=TILE_SIZE, level=0, stride=None):
    stride = stride or tile_size
    w, h = slide_dims(slide, level)
    for y in range(0, h - tile_size, stride):
        for x in range(0, w - tile_size, stride):
            tile = read_region(slide, x, y, tile_size, level)
            if is_tissue(tile):
                yield x, y, tile


def build_phikon(device):
    backbone = AutoModel.from_pretrained("owkin/phikon")

    class Wrap(nn.Module):
        def __init__(self, m):
            super().__init__()
            self.m = m

        def forward(self, x):
            return self.m(pixel_values=x).last_hidden_state[:, 0]  # CLS token, 768-d

    model = Wrap(backbone).eval().to(device)
    tfm = transforms.Compose([
        transforms.Resize(224),
        transforms.ToTensor(),
        transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
    ])
    return model, tfm


@torch.no_grad()
def embed_slide(path, model, tfm, device, batch_size=64, level=0, stride=None):
    slide = openslide.OpenSlide(str(path))
    embs, coords, batch, bxy = [], [], [], []

    def flush():
        if not batch:
            return
        x = torch.stack(batch).to(device)
        embs.append(model(x).float().cpu())
        coords.extend(bxy)
        batch.clear()
        bxy.clear()

    for x, y, tile in iter_tissue_tiles(slide, level=level, stride=stride):
        batch.append(tfm(tile))
        bxy.append((x, y))
        if len(batch) >= batch_size:
            flush()
    flush()
    if not embs:
        return None, None
    return torch.cat(embs, 0), torch.tensor(coords, dtype=torch.long)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--level", type=int, default=0)
    ap.add_argument("--stride", type=int, default=0,
                    help="tile step in px; 0=tile_size (dense, same as UNI). "
                         "512 = every-other-tile subsample (~4x faster, mean-pool ~unchanged)")
    ap.add_argument("--limit", type=int, default=0, help="0 = all; N = smoke test")
    ap.add_argument("--shard", default="", help="i/n: process only every n-th slide starting at i "
                                                "(multi-GPU: one process per card, e.g. 0/2 and 1/2)")
    ap.add_argument("--reverse", action="store_true", help="walk the slide list backwards "
                                                            "(second job meets a running one in the middle)")
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}")
    model, tfm = build_phikon(device)

    raw_dir, out_dir = Path(args.raw_dir), Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    wsi = sorted(p for p in raw_dir.rglob("*")
                 if p.suffix.lower() in (".svs", ".tif", ".tiff", ".ndpi"))
    if args.limit:
        wsi = wsi[:args.limit]
    if args.shard:
        i, n = (int(v) for v in args.shard.split("/"))
        wsi = wsi[i::n]
    if args.reverse:
        wsi = wsi[::-1]
    print(f"found {len(wsi)} slides under {raw_dir}"
          f"{f' (limit {args.limit})' if args.limit else ''}"
          f"{f' shard {args.shard}' if args.shard else ''}{' reversed' if args.reverse else ''}")

    ok = fail = 0
    for i, path in enumerate(wsi, 1):
        out_path = out_dir / f"{path.stem}.pt"
        if out_path.exists():
            print(f"[{i}/{len(wsi)}] skip (done): {path.stem}")
            ok += 1
            continue
        print(f"[{i}/{len(wsi)}] embedding {path.stem} ...")
        try:
            emb, coords = embed_slide(path, model, tfm, device, args.batch_size, args.level,
                                      stride=(args.stride or None))
            if emb is None:
                print(f"  WARNING: no tissue tiles in {path.stem}, skipping")
                fail += 1
                continue
            torch.save({"embeddings": emb, "coords": coords, "source": str(path)}, out_path)
            print(f"  saved {emb.shape[0]} tiles x {emb.shape[1]}d -> {out_path}")
            ok += 1
        except Exception as e:
            print(f"  FAILED on {path.stem}: {e}")
            fail += 1
    print(f"\ndone. ok={ok} fail={fail} out_dir={out_dir}")
    sys.exit(1 if fail and not ok else 0)


if __name__ == "__main__":
    main()
