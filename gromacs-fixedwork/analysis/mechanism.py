"""Confirmatory analyses for the GROMACS sweep (evenly spaced arrivals).
 A. Slope decomposition: energy per request = command part + time-proportional
    part; each regressed on cycles per request with the same relative-error
    weights over the fitted strides. If b exceeds the time-proportional bound,
    which part carries the excess?
 B. Requests without a command: reads completed at latency 1 (served from the
    write buffer, as DRAMsim3 does for a read whose address is queued for
    writing) against reads done minus read commands.
 C. Sensitivity: refit on strides 100, 50, 25, 12 only (points where the stock
    simulator agrees with the patched one) and score strides 10, 9, 8.
 D. Arrival-process sensitivity of the coefficients (Poisson against uniform).
Usage: python3 mechanism.py <gromacs_results> <gemm_results> <analysis_dir>
"""
import glob, json, os, sys
import numpy as np
GMX, GEMM, ANA = sys.argv[1:4]
sys.path.insert(0, ANA)
from read_run import load, ACT_E, RD_E, WR_E
TRAIN = [100, 50, 25, 12, 6]

def tree(res):
    out = {}
    for d in sorted(glob.glob(os.path.join(res, "runs", "*", "*", "s*"))):
        r = load(d); r["dir"] = d; out.setdefault((r["window"], r["arrival"]), {})[r["stride"]] = r
    return out

def wfit(x, y, w):
    b, a = np.polyfit(np.array(x), np.array(y), 1, w=np.array(w)); return a, b

out = {"A": {}, "B": {}, "C": {}, "D": {}}
for label, res in (("GROMACS", GMX), ("GEMM-AI", GEMM)):
    t = tree(res)
    print(f"\n== {label}")
    print("A. slope decomposition (fitted strides, relative-error weights): b_total = b_command + b_time")
    for (w, arr), pts in sorted(t.items()):
        if arr != "uniform" or not all(s in pts for s in TRAIN): continue
        tr = [pts[s] for s in TRAIN]
        x = [p["c_per_req"] for p in tr]; E = [p["E_per_req"] for p in tr]; wt = [1 / e for e in E]
        cmdE = [(p["share_read_xfer"] + p["share_write_xfer"] + p["share_activate"]) * p["E_per_req"] for p in tr]
        timeE = [p["share_residency"] * p["E_per_req"] for p in tr]
        _, bt = wfit(x, E, wt); ac, bc = wfit(x, cmdE, wt); at, bti = wfit(x, timeE, wt)
        cpr = [(p["read_cmds"] + p["write_cmds"]) / p["completed"] for p in tr]
        out["A"][f"{label}/{w}"] = {"b_total": bt, "b_command": bc, "b_time": bti, "a_command": ac,
                                   "cmds_per_req_s100_to_s6": [cpr[0], cpr[-1]]}
        print(f"  {w:12s} b={bt:8.1f}  b_cmd={bc:7.1f}  b_time={bti:8.1f}   commands per request s100 {cpr[0]:.3f} to s6 {cpr[-1]:.3f}")
    if label != "GROMACS": continue
    print("B. reads without a read command, against reads completed at latency 1")
    for (w, arr), pts in sorted(t.items()):
        if arr != "uniform": continue
        for s in (200, 50, 6):
            p = pts[s]; j = json.load(open(os.path.join(p["dir"], "dramsim3.json")))["0"]
            lat1 = int(j["read_latency"].get("1", 0)); gap = p["reads_done"] - p["read_cmds"]
            out["B"].setdefault(w, {})[str(s)] = {"reads_minus_read_cmds": gap, "reads_latency_1": lat1}
        print(f"  {w:12s} " + "  ".join(f"s{s}: gap {out['B'][w][str(s)]['reads_minus_read_cmds']:>6d} lat1 {out['B'][w][str(s)]['reads_latency_1']:>6d}" for s in (200, 50, 6)))
    print("C. refit on strides 100, 50, 25, 12 (stock-equivalent); held out 10, 9, 8; extrapolated 200")
    for (w, arr), pts in sorted(t.items()):
        if arr != "uniform": continue
        tr = [pts[s] for s in (100, 50, 25, 12)]
        x = [p["c_per_req"] for p in tr]; E = [p["E_per_req"] for p in tr]
        a, b = wfit(x, E, [1 / e for e in E])
        err = lambda p: 100 * abs(a + b * p["c_per_req"] - p["E_per_req"]) / p["E_per_req"]
        h = [err(pts[s]) for s in (10, 9, 8)]
        out["C"][w] = {"a": a, "b": b, "holdout_mape": float(np.mean(h)), "holdout_max": float(max(h)), "ext200": err(pts[200])}
        print(f"  {w:12s} a={a:7.1f} b={b:7.1f} held-out MAPE {np.mean(h):.2f}% max {max(h):.2f}%  stride 200 {err(pts[200]):.2f}%")
    print("D. coefficient change, exponential against evenly spaced arrivals (relative-error fit)")
    s = json.load(open(os.path.join(GMX, "summary.json")))["fits"]
    for w in sorted({k.split('/')[0] for k in s}):
        if f"{w}/poisson" not in s: continue
        u, p = s[f"{w}/uniform"], s[f"{w}/poisson"]
        da = 100 * (p["intercept"] / u["intercept"] - 1); db = 100 * (p["slope"] / u["slope"] - 1)
        dh = p["holdout_mape"] - u["holdout_mape"]; de = max(abs(p["extrapolation_err"][k] - u["extrapolation_err"][k]) for k in u["extrapolation_err"])
        out["D"][w] = {"da_pct": da, "db_pct": db, "d_holdout_pts": dh, "d_ext_pts_max": de}
        print(f"  {w:12s} a {da:+.2f}%  b {db:+.2f}%  held-out {dh:+.2f} pts  extrapolation max {de:.2f} pts")
json.dump(out, open(os.path.join(GMX, "mechanism.json"), "w"), indent=1)
