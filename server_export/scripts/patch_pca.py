#!/usr/bin/env python3
"""
Pin the WSI PCA to an exact solver across the whole analysis pipeline.

On a ~100x1024 matrix scikit-learn's auto solver picks randomized SVD with
random_state=None, so every script below produced different numbers on identical
inputs. Measured on GBM: four runs of the unchanged pipeline gave 699/700/704/701
significant proteins (SD 2.2) and only 90.4% of the significant set was common to
all four. Every other random source in these scripts is already seeded
(RandomState(0), KFold(random_state=0)), so this was the single remaining source
of run-to-run drift.

Backs each file up to <name>.pre_pcafix.bak, then byte-compiles to prove the edit
is syntactically sound before reporting success.
"""
import os
import py_compile
import shutil
import sys

SCRIPTS = "/public/home/fjhui/ZW/scripts"

REPLACEMENTS = [
    ("PCA(n_components=min(N_PCS, len(common) - 1, wsi_raw.shape[1]))",
     'PCA(n_components=min(N_PCS, len(common) - 1, wsi_raw.shape[1]),\n'
     '              svd_solver="full")'),
    ("PCA(n_components=min(N_PCS, len(common) - 1, d))",
     'PCA(n_components=min(N_PCS, len(common) - 1, d), svd_solver="full")'),
]

COMMENT = (
    "# svd_solver='full' is load-bearing: at this matrix shape sklearn's auto\n"
    "# solver is randomized SVD with random_state=None, which made this script\n"
    "# non-reproducible (GBM: 699/700/704/701 significant across four identical\n"
    "# runs, 90.4% set overlap). 'full' is exact and deterministic.\n"
)

FILES = ["residual_analysis.py", "residual_analysis_adjusted.py",
         "residual_analysis_wsibatch.py", "randproj_control.py",
         "confound_check.py", "mediation_analysis.py"]


def patch(path):
    with open(path, encoding="utf-8") as f:
        src = f.read()
    if 'svd_solver="full"' in src:
        return "already patched"
    hits = [old for old, _new in REPLACEMENTS if old in src]
    if not hits:
        return "NO MATCH -- inspect by hand"
    shutil.copy2(path, path + ".pre_pcafix.bak")
    out = src
    for old, new in REPLACEMENTS:
        if old in out:
            out = out.replace(old, new)
    # drop the explanatory comment immediately above the (now multi-line) call,
    # matching the indentation of the line that holds it
    lines = out.split("\n")
    for i, ln in enumerate(lines):
        if "PCA(n_components=" in ln:
            indent = " " * (len(ln) - len(ln.lstrip()))
            lines[i:i] = [indent + c for c in COMMENT.rstrip("\n").split("\n")]
            break
    out = "\n".join(lines)
    with open(path, "w", encoding="utf-8") as f:
        f.write(out)
    try:
        py_compile.compile(path, doraise=True)
    except py_compile.PyCompileError as e:
        shutil.copy2(path + ".pre_pcafix.bak", path)
        return f"SYNTAX ERROR, reverted: {e}"
    return "patched + compiles"


ok = True
for name in FILES:
    p = os.path.join(SCRIPTS, name)
    if not os.path.exists(p):
        print(f"  {name:<34} MISSING")
        ok = False
        continue
    status = patch(p)
    print(f"  {name:<34} {status}")
    if "ERROR" in status or "NO MATCH" in status:
        ok = False

print("\n=== verification: every PCA call in the pipeline ===")
import glob
import re
for p in sorted(glob.glob(SCRIPTS + "/*.py")):
    with open(p, encoding="utf-8") as f:
        txt = f.read()
    for m in re.finditer(r"PCA\(n_components=[^\n]*", txt):
        line = m.group(0)
        good = 'svd_solver="full"' in line or 'svd_solver="full"' in txt[m.start():m.start() + 200]
        print(f"  {'OK ' if good else 'BAD'} {os.path.basename(p):<34} {line[:72]}")
        ok = ok and good

sys.exit(0 if ok else 1)
