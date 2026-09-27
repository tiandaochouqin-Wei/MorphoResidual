"""
tile_qc_probe.py

Mechanical QC probe for the 16 H&E tile PNGs used in Fig 3f candidate set
(figures/figdata/tile_{hi,lo}_{i}_{j}.png, i in 0..3, j in 0..1, 256x256 RGB).

Purpose: test whether a STATED, REPRODUCIBLE automated artefact filter
reproduces the authors' by-eye selection of 8 shown / 8 excluded tiles.

Image library: PIL (Pillow) for I/O, numpy for all array math.
No cv2, no scipy — Laplacian and gradient are implemented by hand with
numpy slicing (3x3 kernel convolution done explicitly), as requested.

Run: python tile_qc_probe.py
"""

import os
import numpy as np
from PIL import Image

# ---------------------------------------------------------------------------
# Ground truth from the manuscript's Fig 3f caption / legend text
# ---------------------------------------------------------------------------

SHOWN = {
    ("hi", 1, 0), ("hi", 2, 0), ("hi", 2, 1), ("hi", 3, 1),
    ("lo", 0, 1), ("lo", 1, 1), ("lo", 3, 0), ("lo", 3, 1),
}

REASONS = {
    ("hi", 0, 0): "glass+blurred",
    ("hi", 0, 1): "blood",
    ("hi", 1, 1): "pale stroma",
    ("hi", 3, 0): "sparse stroma",
    ("lo", 0, 0): "pale stroma",
    ("lo", 1, 0): "blood",
    ("lo", 2, 0): "dark+blurred",
    ("lo", 2, 1): "dark+saturated",
    ("hi", 1, 0): "shown",
    ("hi", 2, 0): "shown",
    ("hi", 2, 1): "shown",
    ("hi", 3, 1): "shown",
    ("lo", 0, 1): "shown",
    ("lo", 1, 1): "shown",
    ("lo", 3, 0): "shown",
    ("lo", 3, 1): "shown",
}

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figdata")

GROUPS = ["hi", "lo"]
IS = list(range(4))
JS = list(range(2))


def load_tile(group, i, j):
    path = os.path.join(DATA_DIR, f"tile_{group}_{i}_{j}.png")
    if not os.path.exists(path):
        return None, path
    im = Image.open(path).convert("RGB")
    arr = np.asarray(im).astype(np.float64) / 255.0  # (H, W, 3) in 0-1
    return arr, path


def rgb_to_hsv_np(rgb):
    """Vectorized RGB->HSV, rgb in 0-1, returns H,S,V each in 0-1."""
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    maxc = np.max(rgb, axis=-1)
    minc = np.min(rgb, axis=-1)
    v = maxc
    delta = maxc - minc
    s = np.where(maxc == 0, 0.0, delta / np.where(maxc == 0, 1.0, maxc))

    # Hue (not used downstream but computed for completeness)
    with np.errstate(divide="ignore", invalid="ignore"):
        rc = (maxc - r) / np.where(delta == 0, 1.0, delta)
        gc = (maxc - g) / np.where(delta == 0, 1.0, delta)
        bc = (maxc - b) / np.where(delta == 0, 1.0, delta)
    h = np.zeros_like(maxc)
    h = np.where(maxc == r, bc - gc, h)
    h = np.where(maxc == g, 2.0 + rc - bc, h)
    h = np.where(maxc == b, 4.0 + gc - rc, h)
    h = (h / 6.0) % 1.0
    h = np.where(delta == 0, 0.0, h)
    return h, s, v


def laplacian_3x3(gray):
    """
    Explicit 3x3 Laplacian via numpy slicing (kernel [[0,1,0],[1,-4,1],[0,1,0]]),
    computed on the interior (no padding) to avoid edge artefacts.
    """
    center = gray[1:-1, 1:-1]
    up = gray[0:-2, 1:-1]
    down = gray[2:, 1:-1]
    left = gray[1:-1, 0:-2]
    right = gray[1:-1, 2:]
    lap = up + down + left + right - 4.0 * center
    return lap


def gradient_magnitude(gray):
    """Simple central-difference gradient magnitude on the interior."""
    gy = gray[2:, 1:-1] - gray[0:-2, 1:-1]
    gx = gray[1:-1, 2:] - gray[1:-1, 0:-2]
    mag = np.sqrt(gx ** 2 + gy ** 2)
    return mag


def to_gray(rgb):
    # standard luminance weights
    return 0.299 * rgb[..., 0] + 0.587 * rgb[..., 1] + 0.114 * rgb[..., 2]


def compute_stats(rgb):
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    gray = to_gray(rgb)

    # glass / near-white
    near_white_mask = (r > 0.85) & (g > 0.85) & (b > 0.85)
    frac_near_white = float(near_white_mask.mean())

    # dark pixels + mean luminance
    dark_mask = gray < 0.25
    frac_dark = float(dark_mask.mean())
    mean_luminance = float(gray.mean())

    # HSV saturation
    h, s, v = rgb_to_hsv_np(rgb)
    mean_sat = float(s.mean())
    p10_sat = float(np.percentile(s, 10))
    p90_sat = float(np.percentile(s, 90))

    # focus measures
    lap = laplacian_3x3(gray)
    lap_var = float(lap.var())
    grad = gradient_magnitude(gray)
    grad_var = float(grad.var())
    grad_mean = float(grad.mean())

    # redness / blood proxy
    mean_r_minus_g = float((r - g).mean())
    strong_red_mask = (r - g > 0.15) & (r - b > 0.10)
    frac_strong_red = float(strong_red_mask.mean())

    # extra cheap stats
    std_luminance = float(gray.std())
    frac_saturated_dark = float(((v < 0.35) & (s > 0.45)).mean())  # "dark+saturated" proxy
    mean_value = float(v.mean())

    return dict(
        frac_near_white=frac_near_white,
        frac_dark=frac_dark,
        mean_luminance=mean_luminance,
        std_luminance=std_luminance,
        mean_sat=mean_sat,
        p10_sat=p10_sat,
        p90_sat=p90_sat,
        lap_var=lap_var,
        grad_var=grad_var,
        grad_mean=grad_mean,
        mean_r_minus_g=mean_r_minus_g,
        frac_strong_red=frac_strong_red,
        frac_saturated_dark=frac_saturated_dark,
        mean_value=mean_value,
    )


def main():
    rows = []
    missing = []
    for group in GROUPS:
        for i in IS:
            for j in JS:
                arr, path = load_tile(group, i, j)
                key = (group, i, j)
                shown = key in SHOWN
                reason = "shown" if shown else REASONS.get(key, "?")
                if arr is None:
                    missing.append(path)
                    continue
                stats = compute_stats(arr)
                row = dict(group=group, i=i, j=j, shown=shown, reason=reason, path=path)
                row.update(stats)
                rows.append(row)

    if missing:
        print("MISSING/UNREADABLE FILES:")
        for m in missing:
            print("  ", m)
    else:
        print("All 16 target tile files found and read successfully with PIL.")

    stat_names = [
        "frac_near_white", "frac_dark", "mean_luminance", "std_luminance",
        "mean_sat", "p10_sat", "p90_sat", "lap_var", "grad_var", "grad_mean",
        "mean_r_minus_g", "frac_strong_red", "frac_saturated_dark", "mean_value",
    ]

    # ---------------- 1. Full 16-row table ----------------
    print("\n" + "=" * 140)
    print("1. FULL 16-ROW TABLE")
    print("=" * 140)
    header = f"{'grp':4}{'i':3}{'j':3}{'shown':7}{'reason':16}" + "".join(f"{s:>14}" for s in stat_names)
    print(header)
    for row in rows:
        line = f"{row['group']:4}{row['i']:3}{row['j']:3}{str(row['shown']):7}{row['reason']:16}"
        line += "".join(f"{row[s]:14.4f}" for s in stat_names)
        print(line)

    shown_rows = [r for r in rows if r["shown"]]
    excl_rows = [r for r in rows if not r["shown"]]

    # ---------------- 2. Best single statistic ----------------
    print("\n" + "=" * 140)
    print("2. BEST SINGLE-STATISTIC THRESHOLD SEARCH")
    print("=" * 140)

    def best_threshold_for_stat(stat):
        vals = np.array([r[stat] for r in rows])
        shown_flags = np.array([r["shown"] for r in rows])
        order = np.argsort(vals)
        sorted_vals = vals[order]
        # candidate thresholds: midpoints between consecutive sorted values,
        # plus below-min and above-max
        candidates = []
        candidates.append(sorted_vals[0] - 1e-9)
        for k in range(len(sorted_vals) - 1):
            candidates.append((sorted_vals[k] + sorted_vals[k + 1]) / 2.0)
        candidates.append(sorted_vals[-1] + 1e-9)

        best = None
        for thr in candidates:
            for direction in ("above_shown", "below_shown"):
                if direction == "above_shown":
                    pred_shown = vals > thr
                else:
                    pred_shown = vals < thr
                n_correct = int((pred_shown == shown_flags).sum())
                n_wrong = 16 - n_correct
                if best is None or n_wrong < best[0]:
                    best = (n_wrong, thr, direction, pred_shown.copy())
        return best  # (n_wrong, thr, direction, pred_shown)

    stat_results = []
    for stat in stat_names:
        n_wrong, thr, direction, pred_shown = best_threshold_for_stat(stat)
        stat_results.append((stat, n_wrong, thr, direction, pred_shown))
    stat_results.sort(key=lambda x: x[1])

    for stat, n_wrong, thr, direction, pred_shown in stat_results:
        print(f"  {stat:20} best_thr={thr:10.4f}  rule='{direction}'  misclassified={n_wrong}/16")

    best_stat, best_n_wrong, best_thr, best_direction, best_pred = stat_results[0]
    print(f"\n  BEST SINGLE STATISTIC: {best_stat}")
    print(f"  Rule: predict 'shown' if {best_stat} {'>' if best_direction=='above_shown' else '<'} {best_thr:.4f}")
    print(f"  Misclassified: {best_n_wrong}/16")

    # confusion matrix for best stat
    tp = fp = tn = fn = 0
    misclassified_names = []
    for r, pred in zip(rows, best_pred):
        actual = r["shown"]
        name = f"{r['group']}({r['i']},{r['j']})"
        if actual and pred:
            tp += 1
        elif actual and not pred:
            fn += 1
            misclassified_names.append((name, "actual=shown, predicted=excluded"))
        elif not actual and pred:
            fp += 1
            misclassified_names.append((name, "actual=excluded, predicted=shown"))
        else:
            tn += 1
    print(f"  Confusion matrix: TP(shown->shown)={tp}  FN(shown->excl)={fn}  FP(excl->shown)={fp}  TN(excl->excl)={tn}")
    if misclassified_names:
        print("  Misclassified tiles:")
        for name, desc in misclassified_names:
            print(f"    {name}: {desc}")

    # ---------------- 3. Conjunction of up to 3 thresholds ----------------
    print("\n" + "=" * 140)
    print("3. CONJUNCTION-OF-THRESHOLDS SEARCH (<=3 round-number thresholds)")
    print("=" * 140)

    # Round-number candidate thresholds per statistic (chosen to span plausible
    # cutoffs; kept small so the search is exhaustive and reportable)
    round_candidates = {
        "frac_near_white": [0.01, 0.02, 0.05, 0.10, 0.20, 0.30],
        "frac_dark": [0.01, 0.02, 0.05, 0.10, 0.20, 0.30],
        "mean_luminance": [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70],
        "mean_sat": [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40],
        "lap_var": [0.0005, 0.001, 0.002, 0.003, 0.005, 0.01],
        "grad_var": [0.0005, 0.001, 0.002, 0.003, 0.005, 0.01],
        "mean_r_minus_g": [0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20],
        "frac_strong_red": [0.01, 0.02, 0.05, 0.10, 0.20],
        "frac_saturated_dark": [0.01, 0.02, 0.05, 0.10, 0.20],
    }
    # directions: whether "shown" requires stat to be BELOW or ABOVE the threshold
    directions_map = {
        "frac_near_white": "below",   # too much glass -> excluded
        "frac_dark": "below",         # too dark -> excluded
        "mean_luminance": "above",    # too pale/dark are both bad, but as single dir test "above" a floor
        "mean_sat": "above",          # pale stroma has low saturation -> excluded
        "lap_var": "above",           # blur -> low focus -> excluded
        "grad_var": "above",          # blur -> low grad var -> excluded
        "mean_r_minus_g": "below",    # blood -> high R-G -> excluded
        "frac_strong_red": "below",   # blood -> excluded
        "frac_saturated_dark": "below",  # dark+saturated -> excluded
    }

    shown_flags = np.array([r["shown"] for r in rows])

    def evaluate_conjunction(stat_thr_list):
        pred = np.ones(len(rows), dtype=bool)
        for stat, thr in stat_thr_list:
            vals = np.array([r[stat] for r in rows])
            direction = directions_map[stat]
            if direction == "above":
                cond = vals > thr
            else:
                cond = vals < thr
            pred = pred & cond
        n_wrong = int((pred != shown_flags).sum())
        return n_wrong, pred

    import itertools

    best_combo = None
    stat_list = list(round_candidates.keys())

    # try single, pairs, and triples of statistics (each with its own threshold)
    for r in (1, 2, 3):
        for combo_stats in itertools.combinations(stat_list, r):
            thr_grids = [round_candidates[s] for s in combo_stats]
            for thr_combo in itertools.product(*thr_grids):
                stat_thr_list = list(zip(combo_stats, thr_combo))
                n_wrong, pred = evaluate_conjunction(stat_thr_list)
                if best_combo is None or n_wrong < best_combo[0]:
                    best_combo = (n_wrong, stat_thr_list, pred.copy())

    n_wrong, stat_thr_list, pred = best_combo
    rule_desc = " AND ".join(
        f"{s} {'>' if directions_map[s]=='above' else '<'} {t}" for s, t in stat_thr_list
    )
    print(f"  Best conjunction found (searched 1-3 statistics x round-number grids):")
    print(f"  Rule: predict 'shown' if [{rule_desc}]")
    print(f"  Misclassified: {n_wrong}/16")
    if n_wrong == 0:
        print("  --> PERFECT SEPARATION ACHIEVED")
    else:
        print("  --> Not perfect. Misclassified tiles:")
        for r_, p_ in zip(rows, pred):
            if r_["shown"] != p_:
                name = f"{r_['group']}({r_['i']},{r_['j']})"
                actual = "shown" if r_["shown"] else f"excluded ({r_['reason']})"
                predicted = "shown" if p_ else "excluded"
                print(f"    {name}: actual={actual}, predicted={predicted}")

    # ---------------- 4. Rank-order sanity check ----------------
    print("\n" + "=" * 140)
    print("4. RANK-ORDER SANITY CHECK (per-'patient' best-focus pick after glass/dark filter)")
    print("=" * 140)
    # "patient" = (group, i) pair here since j=0,1 are the two candidate tiles per patient row
    glass_thr = 0.10   # exclude tiles with >10% near-white/glass
    dark_thr = 0.10     # exclude tiles with >10% dark pixels
    picked = []
    for group in GROUPS:
        for i in IS:
            candidates = [r for r in rows if r["group"] == group and r["i"] == i]
            passing = [r for r in candidates if r["frac_near_white"] <= glass_thr and r["frac_dark"] <= dark_thr]
            pool = passing if passing else candidates  # fallback if filter kills both
            best = max(pool, key=lambda r: r["lap_var"])
            picked.append(best)
            filtered_out = [r for r in candidates if r not in pool]
            note = ""
            if filtered_out:
                note = "  (other tile excluded by glass/dark filter: " + \
                    ",".join(f"{r['group']}({r['i']},{r['j']})" for r in filtered_out) + ")"
            print(f"  patient {group}-{i}: candidates j=0,1 -> passing={len(passing)}/2 -> "
                  f"picked j={best['j']} (lap_var={best['lap_var']:.5f}){note}")

    picked_keys = {(r["group"], r["i"], r["j"]) for r in picked}
    match = picked_keys & SHOWN
    print(f"\n  Rank-order-picked 8 tiles: {sorted(picked_keys)}")
    print(f"  Authors' shown 8 tiles:    {sorted(SHOWN)}")
    print(f"  Overlap: {len(match)}/8 match authors' selection")
    print(f"  Matching tiles: {sorted(match)}")
    print(f"  Non-matching (rank-order picked but author excluded): {sorted(picked_keys - SHOWN)}")
    print(f"  Non-matching (author shown but rank-order missed):    {sorted(SHOWN - picked_keys)}")


if __name__ == "__main__":
    main()
