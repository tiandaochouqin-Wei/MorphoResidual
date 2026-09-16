#!/usr/bin/env python3
"""
流式抽取: 逐张切片「下载 -> 抽特征 -> 立刻删原始文件」，raw 目录永远只有 1 张片子。

为什么要有这个脚本
------------------
原来的两步走 (download_wsi.py 全下完, 再 extract_features.py 全抽完) 在 inode 上
是最坏的做法: NBIA 每个 series 解压成一个目录装 N 个 .dcm, 全量下载会让 raw/
峰值达到几十万个文件。华农 /public 的文件数配额只有 200 万且已用到 98%,
超线正在跑的作业会被 GPFS 直接杀掉 (管理员 2026-09 提醒)。

本脚本把峰值压到 O(1): 任何时刻磁盘上只有当前这一张片子的原始文件，
抽完存成单个 .pt 就立刻 rm -rf。393 张片子的 inode 占用 = 393 (只有 .pt)。

用法
----
    export HF_TOKEN=...
    python stream_extract.py --cancer gbm

    # 断点续跑: 已有 .pt 的 series 自动跳过，直接重跑即可
    # 空跑看会处理哪些: --dry-run

LSF 提交见 lsf_stream.sh。
"""
import argparse
import io
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

import openslide  # 必须先于 torch/timm 导入, 避免 libjpeg/libtiff 冲突  # noqa: F401
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_features import build_encoder, embed_slide  # noqa: E402

API_URL = "https://services.cancerimagingarchive.net/nbia-api/services/v4/getImage"
WSI_SUFFIXES = (".svs", ".tif", ".tiff", ".ndpi", ".dcm")


def parse_manifest(path):
    uids, in_list = [], False
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line == "ListOfSeriesToDownload=":
                in_list = True
                continue
            if in_list and line:
                uids.append(line)
    return uids


def fetch_series(uid, out_dir, retries=4, timeout=300):
    """下载并解压一个 series 到 out_dir。失败返回 False。"""
    url = f"{API_URL}?SeriesInstanceUID={uid}"
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            out_dir.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(io.BytesIO(data)) as zf:
                zf.extractall(out_dir)
            return True
        except (urllib.error.URLError, zipfile.BadZipFile, TimeoutError, OSError) as e:
            print(f"    下载失败 {attempt}/{retries}: {e}", flush=True)
            shutil.rmtree(out_dir, ignore_errors=True)
            time.sleep(3 * attempt)
    return False


def slide_files(d):
    return sorted(p for p in Path(d).rglob("*") if p.suffix.lower() in WSI_SUFFIXES)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancer", required=True, help="癌种名, 例如 gbm / pdac")
    ap.add_argument("--zw-root", default=os.environ.get("ZW_ROOT", "/public/home/fjhui/ZW"))
    ap.add_argument("--manifest", default=None, help="默认 <zw>/<cancer>/WSI/manifests/<cancer>_pathology.tcia")
    ap.add_argument("--model", default="uni")
    ap.add_argument("--level", type=int, default=0)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--keep-raw", action="store_true",
                    help="抽完后保留原始文件(默认删除)。只在调试单张片子时用。")
    ap.add_argument("--dry-run", action="store_true", help="只列出待处理 series，不下载不抽取")
    args = ap.parse_args()

    base = Path(args.zw_root) / args.cancer / "WSI"
    manifest = Path(args.manifest) if args.manifest else \
        base / "manifests" / f"{args.cancer}_pathology.tcia"
    scratch = base / "raw_stream"      # 临时区，任何时刻只装一张片子
    out_dir = base / "emb"

    if not manifest.exists():
        sys.exit(f"manifest 不存在: {manifest}")
    out_dir.mkdir(parents=True, exist_ok=True)

    uids = parse_manifest(manifest)
    if args.limit:
        uids = uids[:args.limit]
    print(f"manifest={manifest}  共 {len(uids)} 个 series", flush=True)
    print(f"输出 -> {out_dir}   临时区 -> {scratch} (抽完即删)", flush=True)

    # 断点续跑: 已经有 .pt 的跳过
    done_uids = {p.stem for p in out_dir.glob("*.pt")}
    if args.dry_run:
        todo = [u for u in uids if u not in done_uids]
        print(f"[dry-run] 已完成 {len(uids) - len(todo)}，待处理 {len(todo)}")
        for u in todo[:20]:
            print(f"  {u}")
        if len(todo) > 20:
            print(f"  ... 还有 {len(todo) - 20} 个")
        return

    # 上一轮如果被 LSF 杀掉，临时区可能有残留，先清干净
    if scratch.exists():
        print(f"清理上一轮残留: {scratch}", flush=True)
        shutil.rmtree(scratch, ignore_errors=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}", flush=True)
    model, tfm = build_encoder(args.model, device)

    ok = skip = fail = 0
    for i, uid in enumerate(uids, 1):
        tag = f"[{i}/{len(uids)}] {uid}"
        series_dir = scratch / uid
        try:
            # 1) 下载
            print(f"{tag} 下载中 ...", flush=True)
            if not fetch_series(uid, series_dir):
                print(f"{tag} 下载失败，跳过", flush=True)
                fail += 1
                continue

            paths = slide_files(series_dir)
            if not paths:
                print(f"{tag} 解压后没有可识别的 WSI 文件，跳过", flush=True)
                fail += 1
                continue

            # 2) 抽特征 (一个 series 可能含多张片子)
            for path in paths:
                slide_id = path.stem
                out_path = out_dir / f"{slide_id}.pt"
                if out_path.exists():
                    print(f"{tag} 已存在，跳过 {slide_id}", flush=True)
                    skip += 1
                    continue
                emb, coords = embed_slide(path, model, tfm, device,
                                          batch_size=args.batch_size, level=args.level)
                if emb is None:
                    print(f"{tag} {slide_id} 无组织区域，跳过", flush=True)
                    fail += 1
                    continue
                # 先写临时名再 rename，避免作业被杀时留下半截 .pt
                tmp = out_path.with_suffix(".pt.partial")
                torch.save({"embeddings": emb, "coords": coords, "series": uid}, tmp)
                tmp.rename(out_path)
                print(f"{tag} 保存 {emb.shape[0]} tiles -> {out_path.name}", flush=True)
                ok += 1

        except Exception as e:  # noqa: BLE001
            print(f"{tag} 出错: {e}", flush=True)
            fail += 1
        finally:
            # 3) 无论成败都删掉原始文件，这是本脚本存在的全部意义
            if not args.keep_raw:
                shutil.rmtree(series_dir, ignore_errors=True)

    if not args.keep_raw:
        shutil.rmtree(scratch, ignore_errors=True)

    print(f"\n完成。ok={ok} skip={skip} fail={fail}", flush=True)
    print(f"本次新增 inode ≈ {ok} (每张片子一个 .pt，原始文件已全部回收)", flush=True)
    sys.exit(1 if fail and not ok else 0)


if __name__ == "__main__":
    main()
