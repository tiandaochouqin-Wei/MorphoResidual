"""reclass.py -- re-apply the FINAL kappa(strength, lambda) floor table to the F3a/F3b/F5 cohorts (from logs)."""
import re
from collections import Counter

def kappa(s, lam):
    col = 0 if lam <= 0.55 else (1 if lam <= 0.65 else 2)
    rows = [(0.50, (0.65, 0.65, 0.55)), (0.20, (0.55, 0.45, 0.35)), (0.10, (0.45, 0.40, 0.25)),
            (0.05, (0.40, 0.35, 0.25)), (0.02, (0.35, 0.35, 0.25)), (0.0, (0.30, 0.30, 0.25))]
    for lo, k in rows:
        if s >= lo: return k[col]

def classify(zD, zN, th, n, lam, s):
    if n < 5 or zD < 3: return "R6"
    if zN <= -3: return "R4" if th <= -0.3 else "R4w"
    band = 0.05 if n >= 100 else 0.10
    if zN < 2: return "R1" if abs(th) <= band else "R7"
    thmin = kappa(s, lam) * lam / (1 - lam)
    return ("R2" if th >= thmin else "R3") + ("w" if zN < 3 else "")

tally = {}
for fn in ("vt_F3a.log", "vt_F3b.log"):
    for line in open(fn):
        m = re.match(r"(\w+)\s+([\d.]+) ([\d.]+) ([\d.]+) \| (\d) +(\d+) ([\d.]+) \| +([+\-\d.na]+) +([+\-\d.na]+) ([+\-\d.na]+) ([+\-\d.]+) \| (\w+)", line)
        m2 = re.match(r"(\w+)\s+([\d.]+) fX=([\d.]+) \| (\d) n_sel= *(\d+) \(false (\d+)\) strength=([\d.]+) \| +([+\-\d.]+) +([+\-\d.]+) ([+\-\d.]+) ([+\-\d.]+) \| (\w+)", line)
        if m:
            hyp, lam, delta, fX, sd, n, s, zD, zN, th, thp, old = m.groups(); tag = hyp if hyp != "mix" else f"mix{fX}"
        elif m2:
            hyp, lam, fX, sd, n, fal, s, zD, zN, th, thp, old = m2.groups(); tag = "MIX-" + (hyp if hyp != "mix" else f"mix{fX}")
        else:
            continue
        n = int(n); lam = float(lam); s = float(s)
        if "nan" in zN: new = "R6"
        else:
            new = classify(float(zD), float(zN), float(th), n, lam, s)
        thmin = kappa(s, lam) * lam / (1 - lam)
        tally.setdefault(tag, []).append(new)
        print(f"{tag:12s} lam={lam} n_sel={n:4d} s={s:.2f} zN={zN:>7s} th={th:>7s} thmin={thmin:.2f} old={old:4s} new={new}")
print("\n-- tally with the final kappa table --")
for k, v in tally.items(): print(f"{k:12s} n={len(v):2d} {dict(Counter(v))}")
