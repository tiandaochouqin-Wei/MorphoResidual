#!/usr/bin/env python3
"""
Pixel-level reader-equivalence check for the C1 DICOM route. CPU only, no model.

Each of --a / --b is one of:
  path/to/slide.svs             openslide
  path/to/base_volume.dcm       openslide (DICOM driver)
  wsidicom:path/to/series_dir   wsidicom
The same level-0 256 px regions -- drawn from the 256 px tiling grid and kept only if
extract_features.is_tissue() accepts them on reader A -- are read from both and compared.

Verdict threshold, fixed before any run:
  EQUIVALENT  pooled mean |A-B| < 1.0 grey level AND every channel's mean signed
              difference within +-0.5
  DIFFERENT   otherwise (or level-0 dimensions differ)
A PNG strip is written for a colour check: row 1 = A, row 2 = B, row 3 = |A-B| x10.
"""
import openslide  # 必须先于 torch 导入
import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
import extract_features as EF  # noqa: E402


def open_reader(spec):
    if spec.startswith("wsidicom:"):
        from wsidicom import WsiDicom
        s = WsiDicom.open(spec[len("wsidicom:"):])

        def read(x, y, n):
            return s.read_region(location=(x, y), level=0, size=(n, n)).convert("RGB")
        return (s.size.width, s.size.height), read, s.close
    s = openslide.OpenSlide(spec)

    def read(x, y, n):
        return s.read_region((x, y), 0, (n, n)).convert("RGB")
    return s.level_dimensions[0], read, s.close


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--k", type=int, default=64)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--png", required=True)
    args = ap.parse_args()

    n = EF.TILE_SIZE
    dims_a, read_a, close_a = open_reader(args.a)
    dims_b, read_b, close_b = open_reader(args.b)
    print(f"A {args.a}\n  level-0 {dims_a}\nB {args.b}\n  level-0 {dims_b}")
    if dims_a != dims_b:
        print("VERDICT: DIFFERENT (level-0 dimensions differ; coordinates not comparable)")
        sys.exit(2)

    w, h = dims_a
    rng = np.random.RandomState(args.seed)
    xs, ys = np.arange(0, w - n, n), np.arange(0, h - n, n)
    pairs, tries = [], 0
    while len(pairs) < args.k and tries < args.k * 200:
        tries += 1
        x, y = int(rng.choice(xs)), int(rng.choice(ys))
        tile_a = read_a(x, y, n)
        if EF.is_tissue(tile_a):
            pairs.append((x, y, np.asarray(tile_a, dtype=np.int16),
                          np.asarray(read_b(x, y, n), dtype=np.int16)))
    close_a()
    close_b()
    if not pairs:
        sys.exit("no tissue regions found")

    diffs = np.stack([a - b for _, _, a, b in pairs])
    abs_d = np.abs(diffs)
    pooled_mean = float(abs_d.mean())
    chan_bias = diffs.reshape(-1, 3).mean(axis=0)
    identical = sum(int(not d.any()) for d in diffs)
    print(f"regions compared: {len(pairs)} (grid draws {tries})")
    print(f"pooled mean |diff|: {pooled_mean:.4f}")
    print(f"99.9th pct |diff|: {float(np.percentile(abs_d, 99.9)):.1f}   max: {int(abs_d.max())}")
    print(f"mean signed diff R/G/B (A-B): {chan_bias[0]:+.4f} {chan_bias[1]:+.4f} {chan_bias[2]:+.4f}")
    print(f"bit-identical regions: {identical}/{len(pairs)}")

    m = min(8, len(pairs))
    strip = np.concatenate([
        np.concatenate([p[2] for p in pairs[:m]], axis=1),
        np.concatenate([p[3] for p in pairs[:m]], axis=1),
        np.concatenate([np.clip(np.abs(p[2] - p[3]) * 10, 0, 255) for p in pairs[:m]], axis=1),
    ], axis=0)
    Image.fromarray(strip.astype(np.uint8)).save(args.png)
    print(f"png: {args.png}")

    ok = pooled_mean < 1.0 and bool(np.all(np.abs(chan_bias) < 0.5))
    print(f"VERDICT: {'EQUIVALENT' if ok else 'DIFFERENT'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
