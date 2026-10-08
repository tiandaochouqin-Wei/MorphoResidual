#!/usr/bin/env python3
"""Tiny logged client for the public PDC GraphQL endpoint (metadata only).

Every call appends one row to pdc_query_log.tsv and saves the raw JSON response
under raw/<tag>.json. Nothing user-specific is sent: only public study IDs and
GraphQL field names.
"""
import datetime as _dt
import hashlib
import json
import os
import time

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
LOG = os.path.join(HERE, "pdc_query_log.tsv")
URL = "https://pdc.cancer.gov/graphql"

os.makedirs(RAW, exist_ok=True)
if not os.path.exists(LOG):
    with open(LOG, "w", encoding="utf-8", newline="\n") as f:
        f.write("utc_time\tendpoint\ttag\tquery_summary\thttp_status\tn_errors\t"
                "first_error\tresponse_bytes\tresponse_sha256\tsaved_path\n")


def gql(query, tag, summary, timeout=300, retries=2):
    last = None
    for attempt in range(retries + 1):
        try:
            r = requests.post(URL, json={"query": query}, timeout=timeout)
            break
        except requests.RequestException as e:  # network hiccup -> retry
            last = e
            time.sleep(3)
    else:
        raise last
    body = r.content
    path = os.path.join(RAW, f"{tag}.json")
    with open(path, "wb") as f:
        f.write(body)
    try:
        d = r.json()
    except ValueError:
        d = {"_nonjson": body[:500].decode("utf-8", "replace")}
    errs = d.get("errors") or []
    first = errs[0].get("message", "") if errs else ""
    with open(LOG, "a", encoding="utf-8", newline="\n") as f:
        f.write("\t".join([
            _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            URL, tag, summary.replace("\t", " "), str(r.status_code), str(len(errs)),
            first.replace("\t", " ").replace("\n", " ")[:400], str(len(body)),
            hashlib.sha256(body).hexdigest(), os.path.relpath(path, HERE).replace("\\", "/"),
        ]) + "\n")
    return r.status_code, d
