"""Effect of patch change 4 (read scheduled instead of a stall on a write blocked
by an older read to the same address) on the GROMACS sweep.

For every evenly spaced run, the unmodified simulator runs the same converted
trace for a fixed number of cycles: the patched run's final cycle plus a margin.
If the stock run deadlocks, it stops completing requests. If the blocked state
clears by itself, the stock run completes (nearly) the same work, and its
counters, compared over the same cycle span, measure what change 4 alters.
Usage: python3 stock_check.py <repo> <results> <out_dir>
"""
import glob, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
REPO, RES, OUT = sys.argv[1:4]
STOCK = os.path.join(REPO, "build/bin/dramsim3main-stock")
CFG = os.path.join(REPO, "build/DRAMsim3/configs/DDR4_8Gb_x8_3200.ini")
AWK = os.path.join(REPO, "conversion/convert_pin_to_dramsim.awk")

def one(d):
    meta = json.load(open(f"{d}/meta.json")); cpu = json.load(open(f"{d}/cpu_stats.json"))
    pat = json.load(open(f"{d}/dramsim3.json"))["0"]
    w, s = meta["window"], meta["stride"]
    o = os.path.join(OUT, w, f"s{s:03d}"); os.makedirs(o, exist_ok=True)
    tr = os.path.join(o, "trace.txt")
    with open(tr, "w") as f:
        subprocess.run(["awk", "-f", AWK, "-v", f"stride={s}", os.path.join(REPO, "data/windows", f"{w}.pin")],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    cyc = int(pat["num_cycles"])
    subprocess.run([STOCK, CFG, "-c", str(cyc), "-t", tr, "-o", o], capture_output=True)
    os.remove(tr)
    st = json.load(open(f"{o}/dramsim3.json"))["0"]
    e = lambda c: float(c["total_energy"])
    return {"window": w, "stride": s, "fallbacks": cpu["rw_fallbacks"], "cycles": cyc,
            "patched_done": int(pat["num_reads_done"]) + int(pat["num_writes_done"]),
            "stock_done": int(st["num_reads_done"]) + int(st["num_writes_done"]),
            "patched_cmds": int(pat["num_read_cmds"]) + int(pat["num_write_cmds"]),
            "stock_cmds": int(st["num_read_cmds"]) + int(st["num_write_cmds"]),
            "patched_act": int(pat["num_act_cmds"]), "stock_act": int(st["num_act_cmds"]),
            "energy_diff_pct": 100 * (e(st) / e(pat) - 1)}

runs = sorted(glob.glob(os.path.join(RES, "runs", "*", "uniform", "s*")))
with ThreadPoolExecutor(2) as ex:
    rows = list(ex.map(one, runs))
json.dump(rows, open(os.path.join(OUT, "stock_check.json"), "w"), indent=1)
print(f"{'window':12s}{'s':>4s}{'fb':>5s}{'done pat':>10s}{'done stock':>11s}{'cmds pat':>10s}{'cmds stk':>10s}{'act pat':>9s}{'act stk':>9s}{'dE %':>8s}")
for r in rows:
    print(f"{r['window']:12s}{r['stride']:>4d}{r['fallbacks']:>5d}{r['patched_done']:>10d}{r['stock_done']:>11d}{r['patched_cmds']:>10d}"
          f"{r['stock_cmds']:>10d}{r['patched_act']:>9d}{r['stock_act']:>9d}{r['energy_diff_pct']:>8.3f}")
