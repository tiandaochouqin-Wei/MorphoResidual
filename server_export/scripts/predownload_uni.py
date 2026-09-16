#!/usr/bin/env python3
"""
Run this on mn02 (login node, has internet) BEFORE submitting the GPU job.
gpu02 has no internet access ("Network is unreachable" confirmed live
2026-07-15), so timm's on-the-fly HF download inside extract_features.py
fails there. Pre-populating the local HF cache here means the GPU job
finds the weights already on disk and never needs to reach huggingface.co.

Also: huggingface.co itself is unreachable from this network (confirmed via
curl, both IPv4 and IPv6 -- not an IPv6-routing quirk, the domain is genuinely
blocked). Route through the hf-mirror.com mirror instead (confirmed reachable,
HTTP 200) via HF_ENDPOINT.
"""
import os
import sys

os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")

token = os.environ.get("HF_TOKEN")
if not token:
    sys.exit("HF_TOKEN not set. Run this on mn02 after `source ~/.bashrc`.")

import timm

print("downloading MahmoodLab/UNI weights + config into the local HF cache...")
model = timm.create_model(
    "hf-hub:MahmoodLab/UNI", pretrained=True, init_values=1e-5,
    dynamic_img_size=True,
)
print("done. cached at ~/.cache/huggingface/hub -- gpu02 should now load offline.")
