#!/usr/bin/env python
"""Writes the sidecar review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256 (PS 9 freeze).

It reads SIDECAR_FILES.txt (beside this script; one path per line relative to
MorphoResidual_paper/, '#' lines are comments), adds every file under
review/C1_COVARIATES_2026-11-21/gate_logs/ (recursively, sorted), refuses if any listed file is
missing, if gate_logs/ is missing or empty, or if the sidecar already exists (a frozen sidecar is
never overwritten).  Before writing it also imports c1_covariates_local.py, parses CLOSED_PARTS.txt and
PINNED_HOSTS.txt with the local strict parsers, compares the full list of files the stages require (for the
parsed CLOSED_PARTS.txt) with SIDECAR_FILES.txt and refuses if any required path is missing from the list
or if a cohort that is not fully closed has no pinned host.  It writes

    # sidecar of review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md (written by make_sidecar.py)
    # written <UTC ISO time>
    <sha256> *<path relative to MorphoResidual_paper/>
    ...

with UTF-8, LF line endings and no BOM.  It prints the sidecar's own sha256.

    python make_sidecar.py --dry-run     lists what would be hashed and every missing file; writes nothing
    python make_sidecar.py               writes the sidecar (run only at the freeze)

Afterwards `python c1_covariates_local.py --check-sidecar` verifies the written sidecar against
the files the stages require.
"""
import argparse
import datetime
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
LIST_PATH = os.path.join(HERE, "SIDECAR_FILES.txt")
PRESPEC_REL = "review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md"
SIDECAR_REL = PRESPEC_REL + ".sha256"
GATE_LOGS_REL = "review/C1_COVARIATES_2026-11-21/gate_logs"


def load_local():
    """c1_covariates_local.py, imported without writing any bytecode."""
    sys.dont_write_bytecode = True
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import c1_covariates_local as local
    return local


def required_problems(rels, local=None):
    """Required paths (for the parsed CLOSED_PARTS.txt) that are not in `rels`, plus any parse failure of
    CLOSED_PARTS.txt or PINNED_HOSTS.txt, as a list of strings."""
    L = local or load_local()
    ctx = L.real_ctx()
    try:
        closed = L.parse_closed(ctx.closed_path, ctx.cohorts)
        L.parse_pinned(ctx.pinned_path, ctx.cohorts)
    except L.Fatal as e:
        return ["control file: %s" % e]
    return ["required by c1_covariates_local.py but not listed in SIDECAR_FILES.txt: %s" % r
            for r in ctx.required_rels(closed) if r not in rels]


def unpinned_cohorts(local=None):
    """Cohorts that are not fully closed and have no pinned line in PINNED_HOSTS.txt."""
    L = local or load_local()
    ctx = L.real_ctx()
    closed = L.parse_closed(ctx.closed_path, ctx.cohorts)
    pinned = L.parse_pinned(ctx.pinned_path, ctx.cohorts)
    return [k for k in ctx.cohorts if not L.fully_closed(closed, k) and k not in pinned]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def P(rel):
    return os.path.join(PAPER, *rel.split("/"))


def read_list():
    rels, problems = [], []
    with open(LIST_PATH, "rb") as fh:
        text = fh.read().decode("utf-8")
    for n, ln in enumerate(text.split("\n"), 1):
        line = ln.rstrip("\r").strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("/") or "\\" in line or ".." in line.split("/") or " " in line:
            problems.append("line %d: %r is not a relative path with forward slashes" % (n, line))
        elif line in rels:
            problems.append("line %d: %s is listed twice" % (n, line))
        elif line == SIDECAR_REL:
            problems.append("line %d: the sidecar cannot list itself" % n)
        else:
            rels.append(line)
    return rels, problems


def gate_logs():
    """Every file under gate_logs/, as relative paths; None if the folder does not exist."""
    root = P(GATE_LOGS_REL)
    if not os.path.isdir(root):
        return None
    out = []
    for dp, dn, fn in os.walk(root):
        dn.sort()
        for f in sorted(fn):
            out.append(os.path.relpath(os.path.join(dp, f), PAPER).replace(os.sep, "/"))
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dry-run", action="store_true", help="list what would be hashed and what is missing; write nothing")
    args = ap.parse_args()
    rels, problems = read_list()
    problems += required_problems(rels)
    gl = gate_logs()
    missing = [r for r in rels if not os.path.isfile(P(r))]
    if gl is None:
        missing.append(GATE_LOGS_REL + "/ (folder does not exist)")
        gl = []
    elif not gl:
        missing.append(GATE_LOGS_REL + "/ (folder is empty)")
    for r in gl:
        if r not in rels:
            rels.append(r)
    print("SIDECAR_FILES.txt: %d listed paths; gate_logs/: %d files" % (len(rels) - len(gl), len(gl)))
    for p in problems:
        print("PROBLEM: %s" % p)
    if args.dry_run:
        for r in rels:
            print("  %-8s %s" % ("ok" if os.path.isfile(P(r)) else "MISSING", r))
        for m in missing:
            if m.startswith(GATE_LOGS_REL):
                print("  %-8s %s" % ("MISSING", m))
        print("dry run: %d files would be hashed; %d missing; %d problems; sidecar %s; nothing written"
              % (len(rels) - len([m for m in missing if not m.startswith(GATE_LOGS_REL)]), len(missing),
                 len(problems), "EXISTS (would be refused)" if os.path.exists(P(SIDECAR_REL)) else "not yet written"))
        if missing:
            print("missing:")
            for m in missing:
                print("  " + m)
        unp = unpinned_cohorts()
        if unp:
            print("note: PINNED_HOSTS.txt has no pinned line yet for %s (required before the sidecar is written)"
                  % ", ".join(unp))
        return 0
    if problems or missing:
        print("REFUSED: %d problems, %d missing files: %s" % (len(problems), len(missing), "; ".join(missing[:10])))
        return 1
    unp = unpinned_cohorts()
    if unp:
        print("REFUSED: PINNED_HOSTS.txt has no pinned line for %s, which is not fully closed (the lines are added "
              "from the H-V1 gate record before the sidecar is written)" % ", ".join(unp))
        return 1
    if os.path.exists(P(SIDECAR_REL)):
        print("REFUSED: %s exists; a frozen sidecar is never overwritten" % SIDECAR_REL)
        return 1
    written = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = ["# sidecar of %s (written by make_sidecar.py)" % PRESPEC_REL, "# written %s" % written]
    lines += ["%s *%s" % (sha256(P(r)), r) for r in rels]
    data = ("\n".join(lines) + "\n").encode("utf-8")
    with open(P(SIDECAR_REL), "xb") as fh:          # exclusive: never overwrite
        fh.write(data)
    print("wrote %s: %d files, # written %s" % (SIDECAR_REL, len(rels), written))
    print("sidecar sha256 %s" % sha256(P(SIDECAR_REL)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
