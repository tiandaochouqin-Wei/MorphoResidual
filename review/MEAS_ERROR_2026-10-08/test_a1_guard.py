#!/usr/bin/env python3
"""
test_a1_guard.py -- SYNTHETIC-ONLY test of Amendment 1 (review/POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md):
local_theta_a1.py (nv == 0 guard in _job, amendment SIDECAR) against the frozen local_theta.py.  No project data
is read; every inputs file is built here with test_local_theta.py's own generator (scenario "purity", flagged
is_synthetic=True) and the frozen code is imported and run unchanged.

  T0  local_theta_a1.py differs from local_theta.py only by the header block, the guard and SIDECAR (line diff)
  T1  partial purity (genes below 30 purity-matched patients):
        frozen  --variants purity               -> crashes with the KFold n_samples=0 error (pool path)
        frozen  --variants spline,rnapc,purity  -> crashes the same way (the production command shape)
        a1      --variants purity / spline,rnapc,purity -> complete; workers 2 == workers 1
      gene level: frozen _job raises on a gene with nv = 0; a1 returns no statistic for exactly the genes with
      fewer than 30 purity-matched patients and finite statistics for all others
  T2  complete purity: frozen and a1 (--variants spline,rnapc,purity) write identical outputs (every theta_* csv
      byte-identical, every npz array-identical, manifests identical after dropping timestamps / script hash / --out)
  T3  partial purity: the a1 purity rows of variants.csv use exactly the set genes with >= 30 matched patients, and
      every other output equals a frozen run with --variants spline,rnapc (csv byte-identical; reading.csv
      identical after dropping its purity_* columns; variants.csv identical after dropping the purity rows)
  T4  a1's fail-closed gate on a file NOT flagged synthetic (synthetic arrays) refuses and names the amendment sidecar
  T5  test_local_theta.py --quick pointed at local_theta_a1.py (a scratch copy with four literal substitutions:
      HERE, the import, the script path in the gate test, the script-name filter of G1g) == the frozen suite run
      on the frozen local_theta.py, check by check and number by number (timing lines and the path-bearing details
      of the G1 gate checks excluded), apart from G1j (the frozen run fails G1j because the original sidecar now
      exists; G1j passes for a1 only while the amendment sidecar does not exist)

Usage:  python -B test_a1_guard.py <scratch dir> [--log <file>]      (exit 0 iff every check passes)
"""
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np      # noqa: E402
import pandas as pd     # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
LOG = Path(sys.argv[sys.argv.index("--log") + 1]).resolve() if "--log" in sys.argv else OUT / "test_a1_guard.log"
FROZEN = HERE / "local_theta.py"
A1 = HERE / "local_theta_a1.py"
TLT = HERE / "test_local_theta.py"
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")

_logf = open(LOG, "w", encoding="ascii", errors="replace", newline="\n")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    _logf.write(s + "\n")
    _logf.flush()


FAIL = []


def check(name, cond, detail=""):
    say(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL.append(name)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# the frozen test's generator; test_local_theta reads sys.argv[1] (its scratch dir) at import time
_argv = sys.argv
sys.argv = [_argv[0], str(OUT / "tlt_import_scratch")]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "theory"))
import test_local_theta as tlt   # noqa: E402
import local_theta as lt         # noqa: E402  (frozen)
import local_theta_a1 as a1      # noqa: E402
sys.argv = _argv

N_PAT, G_GENES, PUR_N, MIN_PUR = 120, 60, 40, 30


def build_inputs():
    base = OUT / "scen_pur_base.npz"
    pur = tlt.scenario(base, "purity", n=N_PAT, G=G_GENES, seed=5)
    z = dict(np.load(base, allow_pickle=False))
    assert bool(z["is_synthetic"])
    P = z["protein"].copy()
    rng = np.random.default_rng(11)
    for g in range(30, 45):            # 11 of the 40 purity patients missing -> 29 matched (< 30): nv = 0 in the purity variant
        P[rng.choice(PUR_N, 11, replace=False), g] = np.nan
    for g in range(45, 60):            # 11 missing among the patients without purity -> 40 matched (>= 30)
        P[PUR_N + rng.choice(N_PAT - PUR_N, 11, replace=False), g] = np.nan
    z["protein"] = P
    syn = OUT / "syn_a1.npz"
    np.savez_compressed(syn, **z)
    pats = z["patients"].astype(str)
    part = OUT / "purity_partial.tsv"
    full = OUT / "purity_complete.tsv"
    pd.DataFrame({"patient": pats[:PUR_N], "purity": pur[:PUR_N]}).to_csv(part, sep="\t", index=False)
    pd.DataFrame({"patient": pats, "purity": pur}).to_csv(full, sep="\t", index=False)
    genes = z["genes"].astype(str)
    fam = OUT / "families_syn.tsv"
    rows = [("famA", g) for g in genes[:40]] + [("famB", g) for g in genes[20:]] \
        + [("famC", g) for g in list(genes[:10]) + list(genes[30:45])]     # famC|miss5-20 = only genes below 30 matched
    pd.DataFrame(rows, columns=["set", "gene"]).to_csv(fam, sep="\t", index=False)
    has_pur = np.zeros(N_PAT, bool)
    has_pur[:PUR_N] = True
    matched = (~np.isnan(P) & has_pur[:, None]).sum(0)
    return syn, part, full, fam, matched


def run(script, tag, variants, purity, fam, syn, workers):
    out = OUT / tag
    argv = [sys.executable, "-B", str(script), "--inputs", str(syn), "--synthetic", "--strata", "operator,none,plex",
            "--primary-arm", "operator", "--B", "30", "--seed", "20261007", "--boot", "20", "--boot-seed", "20261008",
            "--sets-file", str(fam), "--save-null", "--workers", str(workers), "--variants", variants, "--out", str(out)]
    if purity is not None:
        argv += ["--purity-file", str(purity)]
    t0 = time.time()
    p = subprocess.run(argv, capture_output=True, text=True, env=ENV, cwd=str(OUT))
    (OUT / f"{tag}.stdout.txt").write_text(p.stdout, encoding="utf-8")
    (OUT / f"{tag}.stderr.txt").write_text(p.stderr, encoding="utf-8")
    say(f"   run {tag}: {Path(script).name} --variants {variants} purity={'-' if purity is None else Path(purity).name} "
        f"workers={workers} -> exit {p.returncode} ({time.time()-t0:.0f}s)")
    return p, out


def last_error_line(stderr):
    lines = [x for x in stderr.strip().splitlines() if x.strip()]
    return lines[-1] if lines else ""


def outputs(d):
    return sorted(x.name for x in Path(d).iterdir() if x.name.startswith("theta_"))


def norm_manifest(path, drop_args=()):
    m = json.load(open(path))
    for k in ("started", "finished", "sha256_script"):
        m.pop(k, None)
    m["args"].pop("out", None)
    for k in drop_args:
        m["args"].pop(k, None)
    return m


def same_npz(a, b):
    za, zb = np.load(a), np.load(b)
    if sorted(za.files) != sorted(zb.files):
        return False
    return all(np.array_equal(za[k], zb[k], equal_nan=True) for k in za.files)


def compare_dirs(da, db, skip=(), drop_args=()):
    """file -> (identical?, how)."""
    res = {}
    fa, fb = outputs(da), outputs(db)
    for f in sorted(set(fa) | set(fb)):
        if f in skip:
            continue
        if f not in fa or f not in fb:
            res[f] = (False, "missing in one run")
        elif f.endswith(".csv"):
            res[f] = ((Path(da) / f).read_bytes() == (Path(db) / f).read_bytes(), "bytes")
        elif f.endswith(".npz"):
            res[f] = (same_npz(Path(da) / f, Path(db) / f), "arrays")
        elif f.endswith("_manifest.json"):
            res[f] = (norm_manifest(Path(da) / f, drop_args) == norm_manifest(Path(db) / f, drop_args),
                      "json minus started/finished/sha256_script/args.out" + (f"/args.{','.join(drop_args)}" if drop_args else ""))
        else:
            res[f] = (False, "unexpected file")
    return res


def main():
    say("test_a1_guard.py  " + time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + f"  host={platform.node()}")
    import sklearn
    say(f"python {platform.python_version()}  numpy {np.__version__}  pandas {pd.__version__}  sklearn {sklearn.__version__}")
    for p in (FROZEN, A1, HERE / "theta_kernel.py", TLT, Path(__file__).resolve()):
        say(f"sha256 {sha(p)}  {p.name}")
    say(f"scratch: {OUT}")
    check("frozen local_theta.py is the frozen file (2d10449a...)",
          sha(FROZEN) == "2d10449a1d30b95cbebb6a73e13aee3ea88d77c6b8d739523f9c4b71dcae7632")

    # ------------------------------------------------------------------ T0
    say("\n--- T0 the amendment's code diff ---")
    fl = FROZEN.read_bytes().split(b"\n")
    al = A1.read_bytes().split(b"\n")
    import difflib
    sm = difflib.SequenceMatcher(a=fl, b=al, autojunk=False)
    ops = [(t, fl[i1:i2], al[j1:j2]) for t, i1, i2, j1, j2 in sm.get_opcodes() if t != "equal"]
    code_ops = [o for o in ops if not (o[0] == "insert" and o[2][0].startswith(b"local_theta_a1.py -- AMENDMENT 1"))]
    hdr = [o for o in ops if o not in code_ops]
    check("T0a header block inserted at the top of the docstring (one insert, 10 lines, before line 3)",
          len(hdr) == 1 and len(hdr[0][2]) == 10 and al[2].startswith(b"local_theta_a1.py") and al[12] == fl[2])
    want = [("replace", [b'SIDECAR = os.path.join(PAPER, "review", "POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md.sha256")'],
             [b'SIDECAR = os.path.join(PAPER, "review", "POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md.sha256")']),
            ("insert", [], [b"    if nv == 0:", b"        return gi, None"])]
    check("T0b the only code changes are SIDECAR and the two guard lines", code_ops == want, repr(code_ops)[:300])
    gi = al.index(b"        return gi, None", al.index(b"def _job(gi):"))
    check("T0c the guard sits immediately after 'nv = int(v.sum())' in _job", al[gi - 2] == b"    nv = int(v.sum())")
    check("T0d no CR byte in local_theta_a1.py (LF, as local_theta.py)", b"\r" not in A1.read_bytes())
    check("T0e a1 keeps PRESPEC_REL and names the amendment sidecar",
          a1.PRESPEC_REL == lt.PRESPEC_REL and a1.SIDECAR.endswith(os.path.join("review", "POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md.sha256")))

    syn, part, full, fam, matched = build_inputs()
    check("synthetic inputs flagged is_synthetic=True", lt.is_synthetic(str(syn)))
    low = np.flatnonzero(matched < MIN_PUR)
    say(f"   inputs: n={N_PAT}, G={G_GENES}; partial purity for {PUR_N}/{N_PAT} patients; genes with < {MIN_PUR} "
        f"purity-matched patients: {len(low)} ({low.min()}..{low.max()}); matched counts {sorted(set(matched.tolist()))}")

    # ------------------------------------------------------------------ T1
    say("\n--- T1 partial purity: frozen crashes, a1 completes ---")
    KF = "n_splits=5 greater than the number of samples: n_samples=0"
    p1, o1 = run(FROZEN, "frozen_partial_purity", "purity", part, fam, syn, 2)
    check("T1a frozen --variants purity crashes with the KFold n_samples=0 error", p1.returncode != 0 and KF in p1.stderr,
          last_error_line(p1.stderr))
    check("T1a' ... raised from _job at tk.fold_ids_for(nv)", "in _job" in p1.stderr and "fid = tk.fold_ids_for(nv)" in p1.stderr)
    check("T1a'' ... and wrote no reading file", not (o1 / "theta_syn_reading.csv").exists())
    p2, o2 = run(FROZEN, "frozen_partial_all3", "spline,rnapc,purity", part, fam, syn, 2)
    check("T1b frozen --variants spline,rnapc,purity (production shape) crashes the same way",
          p2.returncode != 0 and KF in p2.stderr and not (o2 / "theta_syn_reading.csv").exists(), last_error_line(p2.stderr))
    check("T1b' ... after writing the null files only (the UCEC pattern)",
          outputs(o2) == ["theta_syn_null_none.npz", "theta_syn_null_operator.npz", "theta_syn_null_plex.npz"], str(outputs(o2)))
    p3, o3 = run(A1, "a1_partial_purity", "purity", part, fam, syn, 2)
    check("T1c a1 --variants purity completes (exit 0, reading + variants written)",
          p3.returncode == 0 and (o3 / "theta_syn_reading.csv").exists() and (o3 / "theta_syn_variants.csv").exists(),
          last_error_line(p3.stderr))
    p4, o4 = run(A1, "a1_partial_all3_w2", "spline,rnapc,purity", part, fam, syn, 2)
    check("T1d a1 --variants spline,rnapc,purity completes (workers 2)", p4.returncode == 0, last_error_line(p4.stderr))
    p4s, o4s = run(A1, "a1_partial_all3_w1", "spline,rnapc,purity", part, fam, syn, 1)
    check("T1e a1 workers 1 completes", p4s.returncode == 0, last_error_line(p4s.stderr))
    cw = compare_dirs(o4, o4s, drop_args=("workers",))
    check("T1f a1 workers 2 == workers 1 on every output", all(v[0] for v in cw.values()),
          "; ".join(f"{k}:{'=' if v[0] else 'DIFF'}" for k, v in cw.items()))

    # gene level, in-process (workers 1)
    d = lt.load_inputs(str(syn))
    dd, npur = lt.purity_residualised(d, str(part))
    allg = np.arange(G_GENES)
    try:
        lt.run_draws(dd, dd["pcs"][None], None, ["N", "Nc", "D"], False, 1, kind="var:linear", genes=allg)
        frozen_err = None
    except ValueError as e:
        frozen_err = str(e)
    check("T1g gene level: frozen run_draws(var:linear) on the purity-residualised data raises the KFold error",
          frozen_err is not None and "n_samples=0" in frozen_err, str(frozen_err))
    da1 = a1.load_inputs(str(syn))
    dd1, _ = a1.purity_residualised(da1, str(part))
    vo = a1.run_draws(dd1, dd1["pcs"][None], None, ["N", "Nc", "D"], False, 1, kind="var:linear", genes=allg)
    fin = np.isfinite(vo["N"][0]) & np.isfinite(vo["D"][0]) & np.isfinite(vo["Nc"][0])
    check("T1h gene level: a1 gives no statistic exactly for the genes with < 30 matched patients, finite for all others",
          np.array_equal(~fin, matched < MIN_PUR), f"no statistic for {int((~fin).sum())} genes; < 30 matched: {int((matched < MIN_PUR).sum())}")
    nvz = (~np.isnan(dd1["protein"])).sum(0)
    check("T1i purity_residualised leaves either 0 or >= 30 fit patients per gene (so nv in 1..29 never occurs)",
          bool(np.all((nvz == 0) | (nvz >= MIN_PUR))), str(sorted(set(nvz.tolist()))))

    # ------------------------------------------------------------------ T2
    say("\n--- T2 complete purity: frozen == a1 ---")
    p6, o6 = run(FROZEN, "frozen_complete_all3", "spline,rnapc,purity", full, fam, syn, 1)
    p7, o7 = run(A1, "a1_complete_all3", "spline,rnapc,purity", full, fam, syn, 1)
    check("T2a both complete", p6.returncode == 0 and p7.returncode == 0)
    c2 = compare_dirs(o6, o7)
    for k, v in c2.items():
        say(f"      {k}: {'identical' if v[0] else 'DIFFERENT'} ({v[1]})")
    check("T2b every output file identical", len(c2) >= 10 and all(v[0] for v in c2.values()), f"{len(c2)} files")
    vv = pd.read_csv(o7 / "theta_syn_variants.csv")
    check("T2c the purity variant ran on every set (complete purity: no gene dropped)",
          set(vv[vv["variant"] == "purity"]["set"]) == set(vv[vv["variant"] == "spline"]["set"]))

    # ------------------------------------------------------------------ T3
    say("\n--- T3 partial purity: a1 (spline,rnapc,purity) vs frozen (spline,rnapc) ---")
    p5, o5 = run(FROZEN, "frozen_partial_spline_rnapc", "spline,rnapc", part, fam, syn, 1)
    check("T3a frozen --variants spline,rnapc completes", p5.returncode == 0, last_error_line(p5.stderr))
    c3 = compare_dirs(o4s, o5, skip=("theta_syn_reading.csv", "theta_syn_variants.csv"), drop_args=("variants",))
    for k, v in c3.items():
        say(f"      {k}: {'identical' if v[0] else 'DIFFERENT'} ({v[1]})")
    check("T3b every output other than reading/variants identical", len(c3) >= 8 and all(v[0] for v in c3.values()),
          f"{len(c3)} files")
    ra = pd.read_csv(o4s / "theta_syn_reading.csv", dtype=str, keep_default_na=False)
    rb = pd.read_csv(o5 / "theta_syn_reading.csv", dtype=str, keep_default_na=False)
    pcols = [c for c in ra.columns if c.startswith("purity_")]
    check("T3c reading.csv identical after dropping the purity_* columns",
          ra.drop(columns=pcols).equals(rb) and not [c for c in rb.columns if c.startswith("purity_")], f"dropped {pcols}")
    va = pd.read_csv(o4s / "theta_syn_variants.csv", dtype=str, keep_default_na=False)
    vb = pd.read_csv(o5 / "theta_syn_variants.csv", dtype=str, keep_default_na=False)
    va_np = va[va["variant"] != "purity"].reset_index(drop=True)
    vb_cols = list(vb.columns)
    extra_cols = [c for c in va.columns if c not in vb_cols]
    check("T3d variants.csv identical after dropping the purity rows",
          va_np[vb_cols].equals(vb) and all((va_np[c] == "").all() for c in extra_cols), f"extra empty cols {extra_cols}")
    # expected purity n_genes per set: set genes (as built by local_theta) with >= 30 matched patients
    obs, _ = a1.observed(da1, list(a1.tk.STAT_KEYS), False)
    ok = np.isfinite(obs["N"]) & np.isfinite(obs["D"])
    nfit = (~np.isnan(da1["protein"])).sum(0)
    sets = a1.build_sets(da1, ok, nfit, a1.sets_depth(da1), str(fam))
    exp = {}
    for nm, meta in sets.items():
        k = int((matched[meta["idx"]] >= MIN_PUR).sum())
        if k >= a1.MIN_SET:
            exp[nm] = k
    vpur = pd.read_csv(o4s / "theta_syn_variants.csv")
    vpur = vpur[vpur["variant"] == "purity"]
    got = dict(zip(vpur["set"], vpur["n_genes"].astype(int)))
    say("      purity rows (set: n_genes a1 / set genes with >= 30 matched / set size): "
        + ", ".join(f"{nm}: {got.get(nm, '-')}/{exp.get(nm, '-')}/{len(sets[nm]['idx'])}" for nm in sets))
    check("T3e a1 purity rows = exactly the sets with >= 5 genes having >= 30 matched patients, n_genes = that count",
          got == exp, f"got {got} expected {exp}")
    check("T3f the test exercises the guard: some set loses genes and some purity row is dropped",
          any(exp.get(nm, 0) < len(sets[nm]["idx"]) for nm in sets) and len(exp) < len(sets))

    # ------------------------------------------------------------------ T4
    say("\n--- T4 a1 gate ---")
    zz = dict(np.load(syn, allow_pickle=False))
    zz["is_synthetic"] = np.bool_(False)
    nf = OUT / "syn_flag_false.npz"
    np.savez_compressed(nf, **zz)
    p8 = subprocess.run([sys.executable, "-B", str(A1), "--inputs", str(nf), "--strata", "none", "--B", "4", "--boot", "0",
                         "--variants", "none", "--workers", "1", "--out", str(OUT / "a1_gate")],
                        capture_output=True, text=True, env=ENV, cwd=str(OUT))
    msg = (p8.stdout + p8.stderr).strip().splitlines()[-1] if (p8.stdout + p8.stderr).strip() else ""
    exp_sc = os.path.join(str(HERE.parent.parent), "review", "POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md.sha256")
    check("T4a a1 on a non-synthetic-flagged file refuses (no amendment sidecar yet) and names the amendment sidecar",
          p8.returncode != 0 and msg.startswith("REFUSED") and "AMENDMENT1_2026-10-10.md.sha256" in msg, msg)
    check("T4b the amendment sidecar does not exist yet (nothing frozen)", not os.path.exists(exp_sc))

    # ------------------------------------------------------------------ T5
    say("\n--- T5 test_local_theta.py --quick against local_theta_a1.py (scratch copy) vs frozen suite on local_theta.py ---")
    src = TLT.read_bytes()
    subs = [(b"HERE = Path(__file__).resolve().parent\n", b"HERE = Path(r\"" + str(HERE).encode() + b"\")\n"),
            (b"import local_theta as lt   # noqa: E402\n", b"import local_theta_a1 as lt   # noqa: E402\n"),
            (b'str(HERE / "local_theta.py")', b'str(HERE / "local_theta_a1.py")'),
            (b'f.endswith("local_theta.py")', b'f.endswith("local_theta_a1.py")')]
    for old, new in subs:
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    sdir = OUT / "suite"
    sdir.mkdir(exist_ok=True)
    patched = sdir / "test_local_theta_on_a1.py"
    patched.write_bytes(src)
    say(f"   patched copy {patched.name} sha256 {sha(patched)} (4 substitutions)")
    t0 = time.time()
    sa = subprocess.run([sys.executable, "-B", str(patched), str(sdir / "tmp_a1"), "--quick"], capture_output=True,
                        text=True, env=ENV, cwd=str(sdir))
    say(f"   suite on a1: exit {sa.returncode} ({time.time()-t0:.0f}s)")
    t0 = time.time()
    sf = subprocess.run([sys.executable, "-B", str(TLT), str(sdir / "tmp_frozen"), "--quick"], capture_output=True,
                        text=True, env=ENV, cwd=str(sdir))
    say(f"   frozen suite on frozen local_theta.py: exit {sf.returncode} ({time.time()-t0:.0f}s)")
    (sdir / "suite_a1.stdout.txt").write_text(sa.stdout + "\n--- stderr ---\n" + sa.stderr, encoding="utf-8")
    (sdir / "suite_frozen.stdout.txt").write_text(sf.stdout + "\n--- stderr ---\n" + sf.stderr, encoding="utf-8")

    def checks(out):
        return [ln for ln in out.splitlines() if ln.startswith("PASS ") or ln.startswith("FAIL ")]

    def numbers(out):
        keep = []
        for ln in out.splitlines():
            if re.match(r"^\d\d:\d\d:\d\d ", ln) or re.match(r"^\s+\w+: \d+s  \(theta_pop", ln) \
                    or re.match(r"^\s+S3 full-feature run: \d+s$", ln):
                continue                      # log lines (timings, paths) and the S3 / P2 timing lines
            if ln.startswith("PASS G1j") or ln.startswith("FAIL G1j"):
                continue
            if ln.startswith("PASS G1") or ln.startswith("FAIL G1"):
                ln = ln.split("  [", 1)[0]    # gate details carry scratch paths / sidecar names
            keep.append(ln)
        return keep
    ca, cf = checks(sa.stdout), checks(sf.stdout)
    fails_a = [x for x in ca if x.startswith("FAIL")]
    fails_f = [x for x in cf if x.startswith("FAIL")]
    for x in fails_a:
        say("   a1 suite: " + x)
    for x in fails_f:
        say("   frozen suite: " + x[:200])
    say(f"   a1 suite: {len(ca)} checks, {len(fails_a)} FAIL; last line: {sa.stdout.strip().splitlines()[-1] if sa.stdout.strip() else ''}")
    say(f"   frozen suite: {len(cf)} checks, {len(fails_f)} FAIL; last line: {sf.stdout.strip().splitlines()[-1] if sf.stdout.strip() else ''}")
    check("T5a test_local_theta.py --quick on local_theta_a1.py: ALL CHECKS PASSED",
          sa.returncode == 0 and sa.stdout.strip().endswith("ALL CHECKS PASSED"))
    check("T5b the frozen suite on the frozen script fails only G1j (the original sidecar now exists)",
          len(fails_f) == 1 and fails_f[0].startswith("FAIL G1j"))
    na, nf_ = numbers(sa.stdout), numbers(sf.stdout)
    na = [x for x in na if x not in ("ALL CHECKS PASSED",) and not x.startswith("FAILED:")]
    nf_ = [x for x in nf_ if x not in ("ALL CHECKS PASSED",) and not x.startswith("FAILED:")]
    check("T5c every other check line and every printed number identical between the two suite runs (G1 details "
          "and timing lines excluded)",
          na == nf_, f"{len(na)} lines compared")

    say("\nFAILED: " + ", ".join(FAIL) if FAIL else "\nALL A1 CHECKS PASSED")
    _logf.close()
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
