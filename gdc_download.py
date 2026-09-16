#!/usr/bin/env python3
"""gdc_download.py — robust GDC downloader (no gdc-client needed). Reads a gdc-client
manifest (id/filename/md5/size/state), streams each open-access file from the GDC data
endpoint with 4-way parallelism, RESUMES by skipping files already the right size, and
retries. TCGA diagnostic slides are open access (no token).
    nohup python gdc_download.py manifest_kirc_slides.txt slides > gdc_dl.log 2>&1 &
Progress:  ls slides | wc -l"""
import sys, os, csv, time
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

MAN = sys.argv[1] if len(sys.argv) > 1 else "manifest_kirc_slides.txt"
DEST = sys.argv[2] if len(sys.argv) > 2 else "slides"
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 4
DATA = "https://api.gdc.cancer.gov/data/"
os.makedirs(DEST, exist_ok=True)

rows = []
with open(MAN) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        rows.append(row)
print(f"{len(rows)} files -> {DEST} ({WORKERS} workers)", flush=True)


def dl(row):
    fid, fn, size = row["id"], row["filename"], int(row.get("size", 0) or 0)
    out = os.path.join(DEST, fn)
    if os.path.exists(out) and (size == 0 or os.path.getsize(out) == size):
        return (fn, "skip", 0)
    for attempt in range(4):
        try:
            with requests.get(DATA + fid, stream=True, timeout=180) as r:
                r.raise_for_status()
                tmp = out + ".part"
                n = 0
                with open(tmp, "wb") as w:
                    for chunk in r.iter_content(1 << 20):
                        w.write(chunk); n += len(chunk)
                if size and n != size:
                    raise IOError(f"size {n}!={size}")
                os.rename(tmp, out)
                return (fn, "done", n)
        except Exception as e:
            if attempt == 3:
                return (fn, f"FAIL:{e}", 0)
            time.sleep(5 * (attempt + 1))


done = fail = skip = 0
tot = 0
t0 = time.time()
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futs = [ex.submit(dl, r) for r in rows]
    for i, fu in enumerate(as_completed(futs), 1):
        fn, st, n = fu.result()
        tot += n
        if st == "done":
            done += 1
        elif st == "skip":
            skip += 1
        else:
            fail += 1
        if i % 10 == 0 or st.startswith("FAIL"):
            print(f"[{i}/{len(rows)}] done={done} skip={skip} fail={fail} "
                  f"{tot/1e9:.1f}GB {(time.time()-t0)/60:.0f}min  {st}: {fn}", flush=True)
print(f"\nALL DONE. done={done} skip={skip} fail={fail} {tot/1e9:.1f} GB", flush=True)
if fail:
    print("re-run the same command to retry the FAILED files (completed ones skip).")
