#!/usr/bin/env python3
r"""
export_ccrcc_slides.py -- top-up for the one cohort export_package.py missed.

WHY A SEPARATE SCRIPT
export_package.py looked for each cohort at ZW/<cohort>/. That is right for four
of them and wrong for CCRCC, which -- as morpho_env.sh says in its own comment --
"predates the per-cancer directory scheme and lives at the ZW root, with
embeddings under WSI/emb/ccrcc_real". So the big export returned slide-level
embeddings and pathomics for LUAD, UCEC, GBM and PDAC but nothing for CCRCC, the
primary cohort. Rather than re-shipping and re-testing the 400-line exporter for
a one-off gap, this does just the missing piece.

WHAT IT DOES
Mean-pools every tile .pt under the CCRCC embedding directories to one vector per
slide, and does the same for the CCRCC pathomics directory if present. Output
format is byte-identical to what export_package.py produced for the other four
cohorts, so the files drop straight in alongside them.

  slide_mean_ccrcc_<dirname>.npz   slide_ids, matrix, n_tiles, cohort, source_dir
  slide_features_ccrcc_<dirname>.csv   cohort, slide_id, n_tiles, 6 morphometry cols

Read-only. Nothing on the server is modified except the files it writes into --out.

USAGE
  cd /public/home/fjhui/ZW/scripts
  bsub -q smp -n 4 -o ccrcc_slides.out "python -u export_ccrcc_slides.py"

  It discovers the directories itself; pass --emb-dir / --pathomics-dir only if
  discovery reports nothing. --dry-run lists what it would read and exits, which
  costs seconds and is worth doing first if you want to check the paths.

  Then transfer ccrcc_slides_<stamp>.tar.gz back.
"""
import argparse
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

ZW = "/public/home/fjhui/ZW"
PATHOMICS_COLS = ["nuclear_density", "nucleus_count", "mean_nuc_area",
                  "nuc_area_std", "chromatin_od", "eosin_fraction"]


def discover(root, explicit, patterns, label):
    """Directories containing tile .pt files. Explicit paths win; else glob."""
    if explicit:
        d = Path(explicit)
        if not d.is_dir():
            print(f"  [FATAL] --{label} {d} is not a directory")
            return []
        return [d] if any(d.glob("*.pt")) else []
    out = []
    for pat in patterns:
        for p in sorted(root.glob(pat)):
            if p.is_dir() and any(p.glob("*.pt")) and p not in out:
                out.append(p)
    return out


def pool(d, label):
    """One mean-pooled vector per slide. Reads one file at a time; never aborts."""
    import torch
    files = sorted(d.glob("*.pt"))
    ids, vecs, ntiles, bad = [], [], [], 0
    t0 = time.time()
    for i, f in enumerate(files):
        try:
            obj = torch.load(f, map_location="cpu", weights_only=False)
            e = obj["embeddings"] if isinstance(obj, dict) else obj
            e = np.asarray(e, dtype=np.float32)
            if e.ndim != 2 or e.shape[0] == 0:
                bad += 1
                continue
            ids.append(f.stem)
            vecs.append(e.mean(0))
            ntiles.append(int(e.shape[0]))
        except Exception as exc:                      # noqa: BLE001
            bad += 1
            if bad <= 3:
                print(f"    could not read {f.name} ({type(exc).__name__}); continuing")
        if (i + 1) % 100 == 0:
            print(f"    {label}: {i+1}/{len(files)} ({time.time()-t0:.0f}s)", flush=True)
    if not vecs:
        print(f"  [FATAL] {label}: no readable .pt in {d}")
        return None
    M = np.vstack(vecs)
    print(f"  {label}: {M.shape[0]} slides x {M.shape[1]}d "
          f"({bad} unreadable, {time.time()-t0:.0f}s)")
    return np.array(ids), M, np.array(ntiles)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=ZW)
    ap.add_argument("--emb-dir", default=None)
    ap.add_argument("--pathomics-dir", default=None)
    ap.add_argument("--out", default=".")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    root = Path(a.root)
    if not root.is_dir():
        raise SystemExit(f"root not found: {root}")

    # CCRCC's legacy layout: cohort root IS the ZW root, embeddings one level deeper
    embs = discover(root, a.emb_dir,
                    ["WSI/emb/ccrcc*", "WSI/emb/*", "WSI/emb", "WSI/emb_*", "WSI/emb_*/*"],
                    "emb-dir")
    paths = discover(root, a.pathomics_dir,
                     ["WSI/pathomics/ccrcc*", "WSI/pathomics/*", "WSI/pathomics"],
                     "pathomics-dir")

    print(f"[ccrcc] root: {root}")
    print(f"[ccrcc] embedding dirs found: {len(embs)}")
    for d in embs:
        print(f"    {d}   ({len(list(d.glob('*.pt')))} slides)")
    print(f"[ccrcc] pathomics dirs found: {len(paths)}")
    for d in paths:
        print(f"    {d}   ({len(list(d.glob('*.pt')))} slides)")
    if not embs and not paths:
        raise SystemExit(
            "[ccrcc] nothing found. morpho_env.sh says the embeddings are at\n"
            f"  {root}/WSI/emb/ccrcc_real\n"
            "Pass it explicitly:  --emb-dir /public/home/fjhui/ZW/WSI/emb/ccrcc_real")
    if a.dry_run:
        print("[ccrcc] --dry-run: nothing written")
        return

    stamp = time.strftime("%Y%m%d_%H%M")
    stage = Path(a.out).resolve() / f"ccrcc_slides_{stamp}"
    (stage / "embeddings").mkdir(parents=True, exist_ok=True)
    (stage / "pathomics").mkdir(parents=True, exist_ok=True)

    def tag_of(d):
        """Path-derived tag, unique across encoder directories.

        The first run named outputs by the leaf directory alone. CCRCC nests one
        level and BOTH encoders' directories are called 'ccrcc_real'
        (WSI/emb/ccrcc_real and WSI/emb_phikon/ccrcc_real), so the two collided
        and Phikon silently overwrote the 1024-d UNI embeddings that the main
        analysis uses. Build the tag from every path segment after 'WSI'.
        """
        parts = list(d.parts)
        i = max((k for k, p in enumerate(parts) if p.upper() == "WSI"), default=-1)
        return "_".join(parts[i + 1:]) if i >= 0 else d.name

    for d in embs:
        tag = tag_of(d)
        r = pool(d, f"emb:{tag}")
        if r is None:
            continue
        ids, M, nt = r
        outp = stage / "embeddings" / f"slide_mean_ccrcc_{tag}.npz"
        assert not outp.exists(), f"output name collision on {outp.name}: two directories map to the same tag"
        np.savez_compressed(outp, slide_ids=ids, matrix=M, n_tiles=nt,
                            cohort="ccrcc", source_dir=str(d))
        print(f"    -> {outp.name}  ({M.shape[1]}-d)")

    for d in paths:
        tag = tag_of(d)
        r = pool(d, f"pathomics:{tag}")
        if r is None:
            continue
        ids, M, nt = r
        cols = PATHOMICS_COLS[:M.shape[1]] if M.shape[1] <= 6 else \
            [f"f{i}" for i in range(M.shape[1])]
        df = pd.DataFrame(M, columns=cols)
        df.insert(0, "n_tiles", nt)
        df.insert(0, "slide_id", ids)
        df.insert(0, "cohort", "CCRCC")
        df.to_csv(stage / "pathomics" / f"slide_features_ccrcc_{tag}.csv", index=False)

    tar = Path(a.out).resolve() / f"ccrcc_slides_{stamp}.tar.gz"
    with tarfile.open(tar, "w:gz") as tf:
        tf.add(stage, arcname=stage.name)
    total = sum(f.stat().st_size for f in stage.rglob("*") if f.is_file())
    print(f"\n[ccrcc] wrote {tar}  ({total/1048576:.1f} MB staged)")
    print(f"[ccrcc] staging dir left at {stage} -- delete it once the tarball is transferred")


if __name__ == "__main__":
    main()
