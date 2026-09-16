#!/usr/bin/env python3
"""
export_package.py -- MorphoResidual T1-8 / R12(B): one-shot server export.

WHY THIS EXISTS
Several remaining analyses are blocked on material that only exists on the HPC,
and each of them would otherwise cost a separate login-and-transfer round trip:
  T1-3  batch-corrected per-gene tables (make_fig_controls.py hardcodes counts only)
  T1-5  composition_controls per-gene increments; pathomics per-slide features
  T1-8  per-slide embeddings, for the within- vs between-slide variance decomposition
  T0-4  scanner/operator counts, slide manifest
  R12   scripts + pinned-snapshot README, for Code/Data Availability
This packages all of it in a single run.

DESIGN RULES (learned from previous round trips)
  * READ-ONLY on everything except the output tarball. Nothing on the server is
    modified, moved or deleted.
  * NEVER CRASH on a missing path. Anything absent is recorded in the manifest as
    missing and the run continues. A partial export plus an accurate inventory is
    worth far more than a traceback.
  * NO DEPENDENCE on residual_analysis / morpho_env.sh. This script does not import
    the pipeline and does not need a cohort environment sourced, because the whole
    point is that it runs once for all five cohorts. It needs only numpy, pandas
    and torch (torch solely to read the .pt tile files).
  * INVENTORY FIRST. The manifest lists everything found -- including things not
    packaged -- so that if something is missing we learn it from the manifest
    rather than from a failed analysis a week later.
  * SIZE GUARDED. Per-file and total caps, so a stray directory cannot produce a
    multi-gigabyte tarball.

WHAT IT PRODUCES
  morpho_export_<stamp>.tar.gz  containing
    MANIFEST.json          machine-readable: every item, size, sha256, status
    INVENTORY.txt          human-readable summary, including what was MISSING
    scripts/               *.py, *.sh, *.tsv mapping files from the scripts dir
    pinned/                pinned snapshot README + per-cohort pinned result CSVs
    results/               every per-gene CSV found under the results directories
    embeddings/            per-slide mean-pooled vectors + slide ids + tile counts
    pathomics/             per-slide 6-feature morphometry table per cohort
    wsi_meta/              SVS header table (scanner/operator/MPP) if present

USAGE
  cd /public/home/fjhui/ZW/scripts
  python -u export_package.py 2>&1 | tee export_package.log

  The embedding pass reads every .pt tile file, which is tens of GB of I/O. It is
  reduced to one 1024-d vector per slide immediately, so memory stays small, but
  on a busy login node prefer a small batch job:
      bsub -q smp -n 4 -o export_pkg.out "python -u export_package.py"
  To skip that pass entirely and get everything else in seconds:
      python -u export_package.py --skip-embeddings

  Options:
    --root PATH          data root (default /public/home/fjhui/ZW). NOT taken from
                         MORPHO_ROOT: morpho_env.sh sets that to the current
                         cohort's directory, which sent one run hunting for all
                         five cohorts inside ZW/pdac.
    --out PATH           output directory for the tarball (default: cwd)
    --skip-embeddings    skip the per-slide mean-pooling pass
    --max-file-mb N      per-file cap for copied files (default 200)
    --max-total-mb N     total cap for the staged export (default 2000)
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

COHORTS = ["ccrcc", "luad", "ucec", "gbm", "pdac"]
SCRIPT_SUFFIXES = {".py", ".sh", ".tsv", ".md", ".txt", ".yaml", ".yml", ".json"}


# ------------------------------------------------------------------ utilities
def sha12(path, chunk=1 << 20):
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for b in iter(lambda: fh.read(chunk), b""):
                h.update(b)
        return h.hexdigest()[:12]
    except OSError:
        return "unreadable"


def mb(n):
    return n / (1024 * 1024)


class Log:
    """Collects everything found, missing and skipped, for the manifest."""

    def __init__(self):
        self.items, self.missing, self.skipped, self.notes = [], [], [], []

    def found(self, kind, src, dst, size):
        self.items.append(dict(kind=kind, source=str(src), packaged_as=str(dst),
                               size_bytes=int(size), sha256_12=sha12(src) if Path(src).is_file() else ""))

    def miss(self, kind, where, why=""):
        self.missing.append(dict(kind=kind, looked_in=str(where), why=why))
        print(f"  [MISSING] {kind}: {where} {('-- ' + why) if why else ''}")

    def skip(self, kind, where, why):
        self.skipped.append(dict(kind=kind, path=str(where), why=why))
        print(f"  [skipped] {kind}: {where} -- {why}")

    def note(self, s):
        self.notes.append(s)
        print(f"  [note] {s}")


def copy_in(src, dstdir, log, kind, max_file_mb, rename=None):
    src = Path(src)
    if not src.is_file():
        log.miss(kind, src, "not a file")
        return False
    size = src.stat().st_size
    if mb(size) > max_file_mb:
        log.skip(kind, src, f"{mb(size):.0f} MB exceeds the {max_file_mb} MB per-file cap")
        return False
    dstdir.mkdir(parents=True, exist_ok=True)
    dst = dstdir / (rename or src.name)
    shutil.copy2(src, dst)
    log.found(kind, src, dst.relative_to(dst.parents[1]), size)
    return True


def first_existing(*cands):
    for c in cands:
        if c and Path(c).exists():
            return Path(c)
    return None


# --------------------------------------------------------------- discovery
def find_cohort_dirs(root, log):
    """Cohort data directories, discovered rather than assumed."""
    out = {}
    for c in COHORTS:
        hits = [p for p in root.glob("*") if p.is_dir() and p.name.lower() == c]
        if not hits:
            hits = [p for p in root.glob("*") if p.is_dir() and c in p.name.lower()]
        if hits:
            out[c] = sorted(hits, key=lambda p: len(p.name))[0]
        else:
            log.miss("cohort directory", f"{root}/*{c}*")
    return out


def find_emb_dirs(cdir):
    """Tile-embedding directories under a cohort: WSI/emb, WSI/emb_*, emb*."""
    cands = []
    for pat in ("WSI/emb", "WSI/emb_*", "WSI/embeddings", "emb", "emb_*"):
        cands += [p for p in cdir.glob(pat) if p.is_dir()]
    seen, out = set(), []
    for p in cands:
        if p not in seen and any(p.glob("*.pt")):
            seen.add(p)
            out.append(p)
    return out


def find_pathomics_dirs(cdir):
    out = []
    for pat in ("WSI/pathomics", "WSI/pathomics/*", "pathomics", "pathomics/*"):
        out += [p for p in cdir.glob(pat) if p.is_dir() and any(p.glob("*.pt"))]
    return sorted(set(out))


# ------------------------------------------------------- per-slide reduction
def pool_slides(emb_dir, log, label):
    """Mean-pool every .pt in a directory to one vector per slide.

    Each file is {"embeddings": [n_tiles, d], "coords": [...], "source": str}, the
    format written by extract_phikon / extract_fm / pathomics. Files are read one
    at a time and reduced immediately, so peak memory is one slide.
    """
    try:
        import torch
    except ImportError:
        log.skip("embeddings", emb_dir, "torch unavailable in this environment")
        return None
    files = sorted(emb_dir.glob("*.pt"))
    if not files:
        log.miss("tile .pt files", emb_dir)
        return None
    ids, vecs, ntiles, bad = [], [], [], 0
    t0 = time.time()
    for i, f in enumerate(files):
        try:
            d = torch.load(f, map_location="cpu", weights_only=False)
            e = d["embeddings"] if isinstance(d, dict) else d
            e = np.asarray(e, dtype=np.float32)
            if e.ndim != 2 or e.shape[0] == 0:
                bad += 1
                continue
            ids.append(f.stem)
            vecs.append(e.mean(0))
            ntiles.append(int(e.shape[0]))
        except Exception as exc:                       # noqa: BLE001 - never abort the export
            bad += 1
            if bad <= 3:
                log.note(f"{label}: could not read {f.name} ({type(exc).__name__}); continuing")
        if (i + 1) % 250 == 0:
            print(f"    {label}: {i+1}/{len(files)} slides ({time.time()-t0:.0f}s)", flush=True)
    if not vecs:
        log.miss("readable tile files", emb_dir, f"{bad} unreadable")
        return None
    M = np.vstack(vecs)
    print(f"    {label}: {M.shape[0]} slides x {M.shape[1]}d "
          f"({bad} unreadable, {time.time()-t0:.0f}s)")
    return dict(slide_ids=np.array(ids), matrix=M, n_tiles=np.array(ntiles), n_unreadable=bad)


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/public/home/fjhui/ZW",
                    help=("data root. Deliberately NOT read from MORPHO_ROOT, because "
                          "morpho_env.sh sets that variable to the current cohort's own "
                          "directory: a shell that had sourced it sent one run hunting for "
                          "all five cohorts inside ZW/pdac. Pass --root if the root moves."))
    ap.add_argument("--out", default=".")
    ap.add_argument("--skip-embeddings", action="store_true")
    ap.add_argument("--max-file-mb", type=float, default=200.0)
    ap.add_argument("--max-total-mb", type=float, default=2000.0)
    a = ap.parse_args()

    root = Path(a.root)
    if not root.is_dir():
        sys.exit(f"data root not found: {root}\nPass the right one with --root")
    stamp = time.strftime("%Y%m%d_%H%M")
    stage = Path(a.out).resolve() / f"morpho_export_{stamp}"
    stage.mkdir(parents=True, exist_ok=True)
    log = Log()
    print(f"[export] root   : {root}")
    print(f"[export] staging: {stage}\n")

    # ---------------------------------------------------------- 1. inventory
    print("[1/6] cohort directories")
    cdirs = find_cohort_dirs(root, log)
    for c, d in cdirs.items():
        print(f"  {c.upper():6s} -> {d}")

    # ------------------------------------------------------------ 2. scripts
    print("\n[2/6] scripts, environment files and ID mappings")
    sdir = first_existing(root / "scripts", root / "script", Path.cwd())
    if sdir is None:
        log.miss("scripts directory", root / "scripts")
    else:
        print(f"  scripts dir: {sdir}")
        n = 0
        for f in sorted(sdir.iterdir()):
            if f.is_file() and f.suffix.lower() in SCRIPT_SUFFIXES:
                if copy_in(f, stage / "scripts", log, "script", a.max_file_mb):
                    n += 1
        print(f"  copied {n} script/mapping files")
        for extra in ("morpho_env.sh", "residual_analysis.py", "confound_check.py"):
            if not (stage / "scripts" / extra).exists():
                log.miss("key script", sdir / extra,
                         "expected by Code Availability; may live elsewhere")

    # -------------------------------------------------------------- 3. pinned
    print("\n[3/6] pinned result snapshot")
    pin = None
    for cand in sorted(root.glob("**/results_snapshot_*")):
        if cand.is_dir():
            pin = cand
            break
    if pin is None:
        log.miss("pinned snapshot", f"{root}/**/results_snapshot_*")
    else:
        print(f"  snapshot: {pin}")
        for f in sorted(pin.rglob("*")):
            if f.is_file() and f.suffix.lower() in {".csv", ".md", ".txt", ".json"}:
                rel = f.relative_to(pin)
                copy_in(f, stage / "pinned" / rel.parent, log, "pinned result", a.max_file_mb)

    # ------------------------------------------------------------- 4. results
    print("\n[4/6] per-gene result tables (batch suite, composition controls, all of them)")
    res_dirs = [p for p in [root / "results", root / "scripts", *cdirs.values()] if p.is_dir()]
    seen = set()
    n_csv = 0
    for rd in res_dirs:
        for f in sorted(rd.glob("*.csv")) + sorted(rd.glob("results*/*.csv")):
            if not f.is_file() or f.resolve() in seen:
                continue
            seen.add(f.resolve())
            # take every small CSV; sorting out which is which is cheaper locally
            # than another round trip to fetch the one we guessed wrong
            if mb(f.stat().st_size) > a.max_file_mb:
                log.skip("result csv", f, f"{mb(f.stat().st_size):.0f} MB over cap")
                continue
            if copy_in(f, stage / "results", log, "result csv", a.max_file_mb,
                       rename=f"{rd.name}__{f.name}"):
                n_csv += 1
    print(f"  copied {n_csv} CSVs")

    # ------------------------------------------------------------- 5. WSI meta
    print("\n[5/6] slide metadata (scanner / operator / MPP)")
    got_meta = False
    for pat in ("**/svs_header*.csv", "**/slide_meta*.csv", "**/scanner*.csv",
                "**/*svs*batch*.csv"):
        for f in sorted(root.glob(pat)):
            if f.is_file() and copy_in(f, stage / "wsi_meta", log, "wsi metadata", a.max_file_mb):
                got_meta = True
    if not got_meta:
        log.miss("SVS header table", f"{root}/**/svs_header*.csv",
                 "run extract_svs_batch.py first, or say so and we work without scanner labels")

    # ------------------------------------------- 6. per-slide embeddings/pathomics
    print("\n[6/6] per-slide embeddings and pathomics features")
    if a.skip_embeddings:
        log.note("embedding pass skipped by --skip-embeddings")
    else:
        for c, cdir in cdirs.items():
            for ed in find_emb_dirs(cdir):
                tag = f"{c}_{ed.name}"
                total = sum(f.stat().st_size for f in ed.glob("*.pt"))
                print(f"  {c.upper()} {ed}  ({len(list(ed.glob('*.pt')))} slides, "
                      f"{mb(total)/1024:.1f} GB to read)")
                r = pool_slides(ed, log, tag)
                if r is None:
                    continue
                outp = stage / "embeddings" / f"slide_mean_{tag}.npz"
                outp.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(outp, slide_ids=r["slide_ids"], matrix=r["matrix"],
                                    n_tiles=r["n_tiles"], cohort=c, source_dir=str(ed))
                log.found("slide embeddings", ed, outp.relative_to(stage), outp.stat().st_size)

            for pd_ in find_pathomics_dirs(cdir):
                tag = f"{c}_{pd_.name}"
                r = pool_slides(pd_, log, "pathomics_" + tag)
                if r is None:
                    continue
                cols = ["nuclear_density", "nucleus_count", "mean_nuc_area",
                        "nuc_area_std", "chromatin_od", "eosin_fraction"]
                M = r["matrix"]
                df = pd.DataFrame(M, columns=cols[:M.shape[1]] if M.shape[1] <= 6
                                  else [f"f{i}" for i in range(M.shape[1])])
                df.insert(0, "n_tiles", r["n_tiles"])
                df.insert(0, "slide_id", r["slide_ids"])
                df.insert(0, "cohort", c.upper())
                outp = stage / "pathomics" / f"slide_features_{tag}.csv"
                outp.parent.mkdir(parents=True, exist_ok=True)
                df.to_csv(outp, index=False)
                log.found("pathomics features", pd_, outp.relative_to(stage), outp.stat().st_size)

    # ------------------------------------------------------------- manifest
    total_bytes = sum(f.stat().st_size for f in stage.rglob("*") if f.is_file())
    manifest = dict(created=stamp, root=str(root), cohort_dirs={k: str(v) for k, v in cdirs.items()},
                    total_bytes=int(total_bytes), items=log.items,
                    missing=log.missing, skipped=log.skipped, notes=log.notes)
    (stage / "MANIFEST.json").write_text(json.dumps(manifest, indent=1), encoding="utf8")

    lines = [f"MorphoResidual server export  {stamp}",
             f"root: {root}", f"staged size: {mb(total_bytes):.1f} MB", "",
             f"PACKAGED ({len(log.items)} items)"]
    by_kind = {}
    for it in log.items:
        by_kind.setdefault(it["kind"], []).append(it)
    for k, v in sorted(by_kind.items()):
        lines.append(f"  {k:22s} {len(v):4d} files, {mb(sum(i['size_bytes'] for i in v)):8.1f} MB")
    lines += ["", f"MISSING ({len(log.missing)}) -- tell the assistant about these"]
    lines += [f"  {m['kind']}: {m['looked_in']} {m['why']}" for m in log.missing] or ["  none"]
    lines += ["", f"SKIPPED ({len(log.skipped)})"]
    lines += [f"  {s['kind']}: {s['path']} -- {s['why']}" for s in log.skipped] or ["  none"]
    if log.notes:
        lines += ["", "NOTES"] + [f"  {n}" for n in log.notes]
    (stage / "INVENTORY.txt").write_text("\n".join(lines) + "\n", encoding="utf8")

    if mb(total_bytes) > a.max_total_mb:
        print(f"\n[export] STOP: staged {mb(total_bytes):.0f} MB exceeds the "
              f"{a.max_total_mb:.0f} MB total cap. Nothing was archived; the staging "
              f"directory is left at {stage} for inspection. Re-run with a higher "
              f"--max-total-mb once you have looked at INVENTORY.txt.")
        sys.exit(2)

    tarpath = Path(a.out).resolve() / f"morpho_export_{stamp}.tar.gz"
    with tarfile.open(tarpath, "w:gz") as tf:
        tf.add(stage, arcname=stage.name)
    shutil.rmtree(stage, ignore_errors=True)

    print("\n" + "=" * 74)
    print("\n".join(lines))
    print("=" * 74)
    print(f"\n[export] wrote {tarpath}  ({mb(tarpath.stat().st_size):.1f} MB)")
    print("[export] transfer that one file back, and paste the MISSING section above.")


if __name__ == "__main__":
    main()
