"""Idle-state-dependent accounting of the fixed-work sweep.

A fixed-work run charges no energy once its work is done. Over a common time
horizon H the device keeps drawing its idle power P after the run ends, so the
total energy of run i is

    E_H(i; P) = E_i + (H - T_i) * P                         (measured form)
    E_H(c; P) = N*a + H*P + N*c*(b - P)                      (model form)

with N requests, T_i = N*c_i cycles, and a, b the fitted coefficients
(summary.json, least squares on relative error). The fixed-work results are the
case P = 0. Raising the duty cycle (lowering c) saves energy if and only if
b > P to first order; exactly, running at a faster operating point and idling
saves energy against a slower one if and only if

    P < P*(slow, fast) = (E_slow - E_fast) / (T_slow - T_fast),

the measured break-even idle power. A constant energy-per-byte model calibrated
at one stride errs at every other stride by exactly the term N*(c - c_cal)*(b - P).

Idle powers are derived from the configuration's currents with DRAMsim3's own
energy formulas (not measured; entry and exit costs of low-power states are
not included). The horizon H is the duration of the stride-200 run of the same
window and arrival process.

Writes <results>/idle_criterion.json.
Usage: python3 analysis/idle_criterion.py [results_dir] [config.ini]
(default configuration: data/DDR4_8Gb_x8_3200.ini, a copy of DRAMsim3's file)
"""
import configparser, glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_run import load

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results")
CFG = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "data", "DDR4_8Gb_x8_3200.ini")
if not os.path.exists(CFG):
    sys.exit(f"missing {CFG}; pass the DRAMsim3 configuration file as the second argument")
SLOW, MATCHED, FAST = 200, 50, 8          # 2%, 8% and 50% of channel peak
BASELINE_STRIDES = [200, 100, 50, 25, 12, 10, 9, 8, 7, 6, 5]


def idle_powers(path):
    c = configparser.ConfigParser(inline_comment_prefixes=(";", "#"))
    c.read(path)
    g = lambda sec, k: float(c[sec][k].split()[0])
    vdd = g("power", "VDD")
    devices = g("system", "bus_width") / g("dram_structure", "device_width")
    # ranks as DRAMsim3 derives them: channel size over the capacity of one rank
    d = lambda k: g("dram_structure", k)
    megs_per_rank = d("rows") * d("columns") * d("bankgroups") * d("banks_per_group") * d("device_width") * devices / 8 / 2**20
    ranks = g("system", "channel_size") / megs_per_rank
    per_cycle = lambda idd: ranks * vdd * idd * devices
    refresh = ranks * vdd * (g("power", "IDD5AB") - g("power", "IDD3N")) * g("timing", "tRFC") * devices / g("timing", "tREFI")
    return {
        "self-refresh": per_cycle(g("power", "IDD6x")),
        "precharge power-down": per_cycle(g("power", "IDD2P")) + refresh,
        "precharge standby": per_cycle(g("power", "IDD2N")) + refresh,
        "active standby": per_cycle(g("power", "IDD3N")) + refresh,
    }, {"ranks": ranks, "devices_per_rank": devices, "VDD": vdd, "refresh_per_cycle": refresh}


P, cfg_used = idle_powers(CFG)
fits = json.load(open(os.path.join(RES, "summary.json")))["fits"]
groups = {}
for d in sorted(glob.glob(os.path.join(RES, "runs", "*", "*", "s*"))):
    r = load(d)
    groups.setdefault(f"{r['window']}/{r['arrival']}", {})[r["stride"]] = r

out = {"idle_power_units_per_cycle": P, "derived_from": os.path.basename(CFG), "config_values": cfg_used,
       "horizon": f"duration of the stride-{SLOW} run", "comparisons": {}, "per_byte_baseline": {}, "groups": {}}
for key, pts in sorted(groups.items()):
    if not all(s in pts for s in (SLOW, MATCHED, FAST)) or key not in fits:
        continue
    a, b = fits[key]["intercept"], fits[key]["slope"]
    N = pts[MATCHED]["units"]
    E = {s: p["E_total"] for s, p in pts.items()}; T = {s: p["cycles"] for s, p in pts.items()}
    H = T[SLOW]
    EH = lambda s, p: E[s] + (H - T[s]) * p
    g = {"N": N, "H_cycles": H, "fit": {"a": a, "b": b},
         "first_order_margin_b_minus_P": {k: b - v for k, v in P.items()},
         "break_even_P": {f"{slow}_vs_{fast}": (E[slow] - E[fast]) / (T[slow] - T[fast])
                          for slow, fast in ((SLOW, FAST), (MATCHED, FAST), (SLOW, 6), (MATCHED, 6), (SLOW, MATCHED), (FAST, 6)) if fast in E},
         "energy_change_pct": {}, "per_byte_error_pct": {}, "model_vs_measured_pct": {}}
    # energy change of the faster point against the slower one, over the slower run's duration
    for slow, fast in ((SLOW, FAST), (MATCHED, FAST)):
        Hs = T[slow]
        g["energy_change_pct"][f"{fast}_vs_{slow}"] = {
            k: 100 * ((E[fast] + (Hs - T[fast]) * v) / E[slow] - 1) for k, v in [("none (fixed work)", 0.0)] + list(P.items())}
    for name, p in [("none (fixed work)", 0.0)] + list(P.items()):
        cal = EH(MATCHED, p)
        g["per_byte_error_pct"][name] = {str(s): 100 * (cal / EH(s, p) - 1) for s in BASELINE_STRIDES if s in E}
        g["model_vs_measured_pct"][name] = {str(s): 100 * ((N * a + H * p + N * pts[s]["c_per_req"] * (b - p)) / EH(s, p) - 1)
                                            for s in BASELINE_STRIDES if s in E}
    g["c_per_req"] = {str(s): pts[s]["c_per_req"] for s in sorted(pts, reverse=True)}
    out["groups"][key] = g

with open(os.path.join(RES, "idle_criterion.json"), "w") as f:
    json.dump(out, f, indent=1)

print("idle power, units per cycle (derived from %s):" % os.path.basename(CFG))
for k, v in P.items():
    print(f"  {k:22s}{v:9.1f}")
for key, g in out["groups"].items():
    print(f"\n{key}: b = {g['fit']['b']:.1f}; break-even idle power " +
          ", ".join(f"{k} {v:.1f}" for k, v in g["break_even_P"].items()))
    for comp, d in g["energy_change_pct"].items():
        print(f"  stride {comp.replace('_vs_', ' against stride ')}: " + ", ".join(f"{k} {v:+.1f}%" for k, v in d.items()))
    for name, d in g["per_byte_error_pct"].items():
        print(f"  per-byte baseline (stride {MATCHED}), idle {name}: stride 200 {d['200']:+.1f}%, stride 6 {d['6']:+.1f}%")
