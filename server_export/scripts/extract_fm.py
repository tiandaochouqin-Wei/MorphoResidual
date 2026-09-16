#!/usr/bin/env python3
"""
extract_fm.py — third-encoder robustness for MorphoResidual (Tier-1 ③).

Drop-in clone of extract_phikon.py with a --model switch:
  resnet50_imagenet : torchvision ResNet-50 (IMAGENET1K_V2), 2048-d pooled features.
                      A NON-pathology WEAK BASELINE: if the residual signal needs a
                      pathology-specific representation this should underperform.
  hoptimus0         : bioptimus/H-optimus-0 (ViT-g, 1536-d CLS) via timm, open on HF.
                      A THIRD pathology foundation model (after UNI + Phikon).
Same tile size / HSV tissue filter / output format {embeddings, coords, source} as UNI &
Phikon, so residual_analysis.py needs ZERO changes -- point MORPHO_WSI_EMB_DIR here.

GPU. openslide imported BEFORE torch. First run needs internet to fetch weights (do it on
mn02 once, weights cache to ~/.cache/{torch,huggingface}); then on gpu02:
    export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
    export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
    python extract_fm.py --model hoptimus0 --raw-dir <.../WSI/raw> --out-dir <.../emb_hoptimus0> --stride 512
    python extract_fm.py --model resnet50_imagenet --raw-dir <...> --out-dir <.../emb_resnet50> --stride 512
"""
import argparse, sys
from pathlib import Path
import numpy as np
import openslide  # MUST precede torch/torchvision (libjpeg/libtiff conflict)
import torch
import torch.nn as nn
from torchvision import transforms

TILE_SIZE = 256
TISSUE_SAT_THRESH = 15
MIN_TISSUE_FRACTION = 0.35


def is_tissue(tile_rgb):
    hsv = np.array(tile_rgb.convert("HSV"))
    return (hsv[:, :, 1] > TISSUE_SAT_THRESH).mean() > MIN_TISSUE_FRACTION


def iter_tissue_tiles(slide, tile_size=TILE_SIZE, level=0, stride=None):
    stride = stride or tile_size
    w, h = slide.level_dimensions[level]
    for y in range(0, h - tile_size, stride):
        for x in range(0, w - tile_size, stride):
            tile = slide.read_region((x, y), level, (tile_size, tile_size)).convert("RGB")
            if is_tissue(tile):
                yield x, y, tile


def build(model_name, device):
    if model_name == "resnet50_imagenet":
        from torchvision.models import resnet50, ResNet50_Weights
        m = resnet50(weights=ResNet50_Weights.IMAGENET1K_V2)
        m.fc = nn.Identity()                                   # 2048-d pooled features
        model = m.eval().to(device)
        mean, std = (0.485, 0.456, 0.406), (0.229, 0.224, 0.225)
    elif model_name == "hoptimus0":       # GATED on HF (403 without approved access) -- kept for completeness
        import timm
        model = timm.create_model("hf-hub:bioptimus/H-optimus-0", pretrained=True,
                                  init_values=1e-5, dynamic_img_size=False).eval().to(device)
        mean, std = (0.707223, 0.578729, 0.703617), (0.211883, 0.230117, 0.177517)  # model-card values
    elif model_name == "phikon_v2":       # owkin/phikon-v2, ViT-L 1024-d CLS, open (same channel as Phikon)
        from transformers import AutoModel
        backbone = AutoModel.from_pretrained("owkin/phikon-v2")

        class Wrap(nn.Module):
            def __init__(self, m):
                super().__init__(); self.m = m

            def forward(self, x):
                return self.m(pixel_values=x).last_hidden_state[:, 0]
        model = Wrap(backbone).eval().to(device)
        mean, std = (0.485, 0.456, 0.406), (0.229, 0.224, 0.225)
    elif model_name == "hibou_b":         # histai/hibou-b, ViT-B 768-d, open (different lab)
        from transformers import AutoModel
        backbone = AutoModel.from_pretrained("histai/hibou-b", trust_remote_code=True)

        class Wrap(nn.Module):
            def __init__(self, m):
                super().__init__(); self.m = m

            def forward(self, x):
                out = self.m(pixel_values=x)
                return out.pooler_output if getattr(out, "pooler_output", None) is not None else out.last_hidden_state[:, 0]
        model = Wrap(backbone).eval().to(device)
        mean, std = (0.7068, 0.5755, 0.722), (0.195, 0.2316, 0.1816)   # hibou model-card values
    else:
        sys.exit(f"unknown --model {model_name}")
    tfm = transforms.Compose([transforms.Resize(224), transforms.ToTensor(),
                              transforms.Normalize(mean=mean, std=std)])
    return model, tfm


@torch.no_grad()
def embed_slide(path, model, tfm, device, batch_size=64, level=0, stride=None):
    slide = openslide.OpenSlide(str(path))
    embs, coords, batch, bxy = [], [], [], []

    def flush():
        if not batch:
            return
        x = torch.stack(batch).to(device)
        embs.append(model(x).float().cpu()); coords.extend(bxy)
        batch.clear(); bxy.clear()

    for x, y, tile in iter_tissue_tiles(slide, level=level, stride=stride):
        batch.append(tfm(tile)); bxy.append((x, y))
        if len(batch) >= batch_size:
            flush()
    flush()
    if not embs:
        return None, None
    return torch.cat(embs, 0), torch.tensor(coords, dtype=torch.long)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=["resnet50_imagenet", "hoptimus0", "phikon_v2", "hibou_b"])
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--level", type=int, default=0)
    ap.add_argument("--stride", type=int, default=0, help="0=dense; 512=every-other-tile (~4x faster)")
    ap.add_argument("--limit", type=int, default=0, help="0=all; N=smoke test")
    ap.add_argument("--shard", default="", help="i/n: every n-th slide from i (one process per GPU)")
    ap.add_argument("--reverse", action="store_true", help="walk slide list backwards")
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device} model={args.model}")
    model, tfm = build(args.model, device)
    raw_dir, out_dir = Path(args.raw_dir), Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    wsi = sorted(p for p in raw_dir.rglob("*") if p.suffix.lower() in (".svs", ".tif", ".tiff", ".ndpi"))
    if args.limit:
        wsi = wsi[:args.limit]
    if args.shard:
        i, n = (int(v) for v in args.shard.split("/"))
        wsi = wsi[i::n]
    if args.reverse:
        wsi = wsi[::-1]
    print(f"found {len(wsi)} slides under {raw_dir}"
          f"{f' shard {args.shard}' if args.shard else ''}{' reversed' if args.reverse else ''}")

    ok = fail = 0
    for i, path in enumerate(wsi, 1):
        out_path = out_dir / f"{path.stem}.pt"
        if out_path.exists():
            ok += 1; continue
        print(f"[{i}/{len(wsi)}] {path.stem} ...")
        try:
            emb, coords = embed_slide(path, model, tfm, device, args.batch_size, args.level,
                                      stride=(args.stride or None))
            if emb is None:
                print("  WARNING: no tissue tiles, skipping"); fail += 1; continue
            torch.save({"embeddings": emb, "coords": coords, "source": str(path)}, out_path)
            print(f"  saved {emb.shape[0]} x {emb.shape[1]}d"); ok += 1
        except Exception as e:
            print(f"  FAILED: {e}"); fail += 1
    print(f"\ndone. ok={ok} fail={fail} out_dir={out_dir}")
    sys.exit(1 if fail and not ok else 0)


if __name__ == "__main__":
    main()
