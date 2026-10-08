#!/usr/bin/env python3
"""ts_slide_qc.py -- per-slide acquisition header census + tile counts for the TS (or DX) arm.

Reads ONLY SVS headers (openslide properties) and the SHAPE of each embedding tensor. It never looks at
embedding values, protein or RNA, so it is safe to run before the residual analysis and is the place where
magnification / MPP (REMAINING_WORK C3: "记录 MPP") is recorded. Output: <out>.tsv one row per slide.

    export LD_LIBRARY_PATH=/public/home/fjhui/miniconda3/lib:$LD_LIBRARY_PATH
    python ts_slide_qc.py --raw-dir <ROOT>/slides_ts --emb-dir <ROOT>/emb_phikon_ts --out <ROOT>/results_ts/ts_slide_qc.tsv
"""
import argparse, os, re, sys
from pathlib import Path
import openslide  # MUST precede torch (libjpeg/libtiff clash, see extract_phikon.py)
import torch

PROPS = ["openslide.vendor", "openslide.objective-power", "openslide.mpp-x", "openslide.mpp-y",
         "aperio.AppMag", "aperio.MPP", "aperio.ScanScope ID", "aperio.Date", "aperio.Filename"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--emb-dir", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rx = re.compile(r"^(TCGA-\w\w-\w{4})-(\d\d[A-Z])-(\d\d)-(\w+?)\.")
    rows = []
    for p in sorted(Path(a.raw_dir).rglob("*.svs")):
        m = rx.match(p.name)
        case, code, portion, label = m.groups() if m else (p.name[:12], "", "", "")
        row = {"case": case, "sample_code": code, "portion": portion, "slide_label": label,
               "file": p.name, "bytes": p.stat().st_size}
        try:
            s = openslide.OpenSlide(str(p))
            row.update({"width": s.dimensions[0], "height": s.dimensions[1], "levels": s.level_count})
            for k in PROPS:
                row[k] = s.properties.get(k, "")
            s.close()
        except Exception as e:
            row["error"] = str(e)[:80]
        if a.emb_dir:
            ep = Path(a.emb_dir) / (p.stem + ".pt")
            if ep.exists():
                row["n_tiles"] = int(torch.load(ep, map_location="cpu")["embeddings"].shape[0])
            else:
                row["n_tiles"] = ""
        rows.append(row)
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")
    print(f"{len(rows)} slides -> {a.out}")
    mag = {}
    for r in rows:
        k = (r.get("openslide.objective-power", ""), r.get("openslide.mpp-x", "")[:5])
        mag[k] = mag.get(k, 0) + 1
    print("objective-power / mpp-x (5 chars) counts:", mag)


if __name__ == "__main__":
    main()
