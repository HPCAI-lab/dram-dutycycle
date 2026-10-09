"""Fits and summary statistics for the fixed-work sweep.

For every window and arrival process: the linear model
    energy per completed request = a + b * cycles per completed request
fitted on the training strides by least squares on relative error (the fit used
in the paper; "fits") and by ordinary least squares ("fits_ols"), with R^2,
leave-one-out error, the error on the held-out strides (inside the fitted range)
and on the extrapolated strides (outside it). Operating points at or past
saturation (completed throughput within 1% of the window's highest) are flagged,
and the held-out error is also given without them. Also: the
energy-delay-product minimum under three delay definitions and, at every
stride, the spread of each energy share across the kernel windows.

Writes <results>/all_runs.csv and <results>/summary.json.
Usage: python3 analysis/analyze_sweep.py [results_dir]
"""
import csv, glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_run import load

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results")
RUNS = os.path.join(RES, "runs")
TRAIN, HOLD, EXTRA = [100, 50, 25, 12, 6], [10, 9, 8, 7], [200, 5]
KERNEL = ["g1_nb25", "g2_nb50", "g3_nb75", "g4_spread", "g5_fft", "g6_gather"]
# Bounds on the time-proportional energy per cycle (simulator units), from the
# DDR4_8Gb_x8_3200.ini currents: two ranks of eight devices, VDD 1.2 V,
# IDD2N 37 mA (precharge standby), IDD3N 52 mA (active standby), refresh
# (IDD5AB - IDD3N) over tRFC = 560 cycles once per tREFI = 12480 cycles per rank.
REFRESH = 2 * 1.2 * (250 - 52) * 560 * 8 / 12480
BRACKET = (2 * 1.2 * 37 * 8 + REFRESH, 2 * 1.2 * 52 * 8 + REFRESH)

rows = [load(d) for d in sorted(glob.glob(os.path.join(RUNS, "*", "*", "s*")))
        if os.path.exists(os.path.join(d, "dramsim3.json"))]
if not rows:
    sys.exit(f"no runs under {RUNS}")
with open(os.path.join(RES, "all_runs.csv"), "w", newline="") as f:
    keys = [k for k in rows[0] if not isinstance(rows[0][k], (list, dict))]
    w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
    w.writeheader(); w.writerows(rows)

def fit(pts, relative=True):
    """Least squares on relative error (weights 1/E) or ordinary least squares."""
    x = np.array([p["c_per_req"] for p in pts]); y = np.array([p["E_per_req"] for p in pts])
    b, a = np.polyfit(x, y, 1, w=1 / y if relative else None); pred = a + b * x
    return a, b, 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)

def evaluate(pts, relative):
    tr = [pts[s] for s in TRAIN]; a, b, r2 = fit(tr, relative)
    err = lambda p: 100 * abs(a + b * p["c_per_req"] - p["E_per_req"]) / p["E_per_req"]
    loo = []
    for i in range(len(tr)):
        a2, b2, _ = fit(tr[:i] + tr[i + 1:], relative); p = tr[i]
        loo.append(100 * abs(a2 + b2 * p["c_per_req"] - p["E_per_req"]) / p["E_per_req"])
    peak = max(64 / (p["c_per_req"] * 0.63) for p in pts.values())
    sat = sorted((s for s, p in pts.items() if 64 / (p["c_per_req"] * 0.63) >= 0.99 * peak), reverse=True)
    hold = [err(pts[s]) for s in HOLD if s in pts]
    hold_unsat = [err(pts[s]) for s in HOLD if s in pts and s not in sat]
    cfit = [p["c_per_req"] for p in tr]
    return {"intercept": a, "slope": b, "r2": r2, "insample_mape": float(np.mean([err(p) for p in tr])),
            "loo_mape": float(np.mean(loo)),
            "holdout_mape": float(np.mean(hold)) if hold else None, "holdout_max": float(np.max(hold)) if hold else None,
            "holdout_mape_unsaturated": float(np.mean(hold_unsat)) if hold_unsat else None,
            "holdout_max_unsaturated": float(np.max(hold_unsat)) if hold_unsat else None,
            "extrapolation_err": {str(s): err(pts[s]) for s in EXTRA if s in pts},
            "saturated_strides": sat,
            "c_per_req": {str(s): pts[s]["c_per_req"] for s in sorted(pts, reverse=True)},
            "fitted_c_range": [min(cfit), max(cfit)],
            "slope_within_bracket": bool(BRACKET[0] <= b <= BRACKET[1])}

res = {"bracket_units_per_cycle": BRACKET, "fit_method": "least squares on relative error (weights 1/E)",
       "integrity_failures": [], "fits": {}, "fits_ols": {}, "edp_min_stride": {}, "spread_across_kernel_windows": {}}
res["integrity_failures"] = [(r["window"], r["arrival"], r["stride"]) for r in rows
                             if not r["cpu_ok"] or r["completed"] != r["records"] or r["rw_fallbacks"] != 0]
for g in sorted({(r["window"], r["arrival"]) for r in rows}):
    pts = {r["stride"]: r for r in rows if (r["window"], r["arrival"]) == g}
    if not all(s in pts for s in TRAIN):
        continue
    key = "/".join(g)
    res["fits"][key] = evaluate(pts, relative=True)
    res["fits_ols"][key] = evaluate(pts, relative=False)
    nine = [pts[s] for s in TRAIN + HOLD if s in pts]
    res["edp_min_stride"][key] = {
        "energy_x_completion_time": min(nine, key=lambda p: p["E_per_req"] * p["c_per_req"])["stride"],
        "energy_x_read_latency_accept": min(nine, key=lambda p: p["E_per_req"] * p["lat_acc_mean"])["stride"],
        "energy_x_read_latency_arrival": min(nine, key=lambda p: p["E_per_req"] * p["lat_arr_mean"])["stride"]}
for arrival in sorted({r["arrival"] for r in rows}):
    for s in sorted({r["stride"] for r in rows}):
        sel = [r for r in rows if r["arrival"] == arrival and r["stride"] == s and r["window"] in KERNEL]
        if len(sel) < 2:
            continue
        res["spread_across_kernel_windows"][f"{arrival}/s{s}"] = {
            k: [min(r[k] for r in sel), max(r[k] for r in sel)]
            for k in ("share_residency", "share_read_xfer", "share_activate", "read_row_hit", "rank_active_frac")}
with open(os.path.join(RES, "summary.json"), "w") as f:
    json.dump(res, f, indent=1)
print(f"{len(rows)} runs; integrity failures: {res['integrity_failures'] or 'none'}")
for name, fits in (("least squares on relative error", res["fits"]), ("ordinary least squares", res["fits_ols"])):
    print(f"\n{name}")
    print(f"{'window/arrival':22s}{'a':>10s}{'b':>10s}{'R2':>10s}{'LOO%':>7s}{'held%':>7s}{'max%':>7s}{'unsat%':>7s}{'ext200%':>8s}{'ext5%':>7s}  saturated")
    for k, f in fits.items():
        e = f["extrapolation_err"]
        print(f"{k:22s}{f['intercept']:>10.1f}{f['slope']:>10.2f}{f['r2']:>10.6f}{f['loo_mape']:>7.2f}{f['holdout_mape']:>7.2f}"
              f"{f['holdout_max']:>7.2f}{f['holdout_mape_unsaturated']:>7.2f}{e.get('200', float('nan')):>8.2f}{e.get('5', float('nan')):>7.2f}  {f['saturated_strides']}")
