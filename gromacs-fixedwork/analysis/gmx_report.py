"""GROMACS fixed-work sweep against the GEMM-AI sweep.

Reads the run trees of both sweeps (read_run.load, which checks every energy
sum against the simulator's total) and the summary.json / idle_criterion.json
of the GROMACS tree, and prints:
  1. per GROMACS window: hit rate, coefficients, errors, saturation
  2. requests per command and bus occupancy against 4/c (duty-cycle definition)
  3. standby-and-refresh and read-transfer shares by stride, GROMACS windows
     against the GEMM-AI 85% and 48% patterns, and the cross-workload spread
  4. race-to-idle criterion for the GROMACS windows
Writes gmx_report.json next to the GROMACS summary.
Usage: python3 gmx_report.py <gromacs_results> <gemm_results> <repo_analysis_dir>
"""
import glob, json, os, sys
GMX, GEMM, ANA = sys.argv[1:4]
sys.path.insert(0, ANA)
from read_run import load

def tree(res):
    out = {}
    for d in sorted(glob.glob(os.path.join(res, "runs", "*", "*", "s*"))):
        r = load(d); out.setdefault((r["window"], r["arrival"]), {})[r["stride"]] = r
    return out

g, m = tree(GMX), tree(GEMM)
summ = json.load(open(os.path.join(GMX, "summary.json")))
idle = json.load(open(os.path.join(GMX, "idle_criterion.json")))
STR = [200, 100, 50, 25, 12, 10, 9, 8, 7, 6, 5]
P85 = ["w1_q", "w2_k", "w3_v", "w4_o", "w6_ff2"]
rep = {"windows": {}, "occupancy": {}, "shares": {}, "criterion": {}}

print("1. GROMACS windows, evenly spaced arrivals, relative-error fit")
print(f"{'window':14s}{'hit@50':>8s}{'act/req@50':>11s}{'a':>10s}{'b':>9s}{'R2':>9s}{'LOO':>6s}{'held':>6s}{'hmax':>6s}{'e200':>6s}{'e5':>6s}  saturated")
for (w, arr), pts in sorted(g.items()):
    if arr != "uniform": continue
    f = summ["fits"][f"{w}/{arr}"]; e = f["extrapolation_err"]
    print(f"{w:14s}{100*pts[50]['read_row_hit']:>8.2f}{pts[50]['act_per_req']:>11.3f}{f['intercept']:>10.1f}{f['slope']:>9.1f}{f['r2']:>9.5f}"
          f"{f['loo_mape']:>6.2f}{f['holdout_mape']:>6.2f}{f['holdout_max']:>6.2f}{e.get('200',float('nan')):>6.2f}{e.get('5',float('nan')):>6.2f}  {f['saturated_strides']}")
    rep["windows"][w] = {"read_row_hit_s50": pts[50]["read_row_hit"], "act_per_req_s50": pts[50]["act_per_req"], **f,
                         "poisson": summ["fits"].get(f"{w}/poisson"), "ols": summ["fits_ols"].get(f"{w}/{arr}"),
                         "edp": summ["edp_min_stride"].get(f"{w}/{arr}"), "edp_poisson": summ["edp_min_stride"].get(f"{w}/poisson")}

print("\n2. Requests per command and bus occupancy (4 x commands / cycles) against 4/c, evenly spaced")
for (w, arr), pts in sorted(g.items()):
    if arr != "uniform": continue
    row = []
    for s in STR:
        p = pts[s]; cmds = p["read_cmds"] + p["write_cmds"]
        occ = 4 * cmds / p["cycles"]; dc = 4 / p["c_per_req"]
        row.append((s, p["completed"] / cmds, occ, dc))
        rep["occupancy"].setdefault(w, {})[str(s)] = {"req_per_cmd": p["completed"] / cmds, "req_per_read_cmd": p["req_per_read_cmd"],
                                                      "bus_occupancy": occ, "four_over_c": dc}
    print(f"{w:14s}" + "  ".join(f"s{s}:{rq:.3f}/{o:.3f}/{d:.3f}" for s, rq, o, d in row[::2]))

print("\n3. Standby-and-refresh share (%) by stride, evenly spaced")
print(f"{'stride':>7s}" + "".join(f"{w[:10]:>11s}" for (w, a) in sorted(g) if a == "uniform") + f"{'GEMM85':>14s}{'GEMM48':>9s}{'spread all':>11s}{'spread 85+G':>12s}")
for s in STR:
    gv = {w: 100 * pts[s]["share_residency"] for (w, a), pts in sorted(g.items()) if a == "uniform"}
    m85 = [100 * m[(w, "uniform")][s]["share_residency"] for w in P85]
    m48 = 100 * m[("w5_ff1", "uniform")][s]["share_residency"]
    allv = list(gv.values()) + m85 + [m48]
    sp_all = max(allv) - min(allv); sp85 = max(list(gv.values()) + m85) - min(list(gv.values()) + m85)
    rd = {w: 100 * pts[s]["share_read_xfer"] for (w, a), pts in sorted(g.items()) if a == "uniform"}
    rep["shares"][str(s)] = {"residency_gromacs": gv, "read_xfer_gromacs": rd, "residency_gemm85": [min(m85), max(m85)],
                             "residency_gemm48": m48, "spread_all_windows": sp_all, "spread_gromacs_and_gemm85": sp85}
    print(f"{s:>7d}" + "".join(f"{v:>11.2f}" for v in gv.values()) + f"{min(m85):>8.2f}-{max(m85):<5.2f}{m48:>9.2f}{sp_all:>11.2f}{sp85:>12.2f}")

print("\n4. Race-to-idle criterion (stride 8 against stride 50, horizon = slower run)")
P = idle["idle_power_units_per_cycle"]
for key, gr in sorted(idle["groups"].items()):
    if not key.endswith("/uniform"): continue
    ch = gr["energy_change_pct"]["8_vs_50"]; be = gr["break_even_P"]
    print(f"{key:22s} b={gr['fit']['b']:.1f}  P*(50,8)={be['50_vs_8']:.1f}  P* range={min(be.values()):.1f}..{max(be.values()):.1f}  " +
          ", ".join(f"{k.split()[0]} {v:+.1f}%" for k, v in ch.items()))
    pb = gr["per_byte_error_pct"]
    print(f"{'':22s} per-byte error s200/s6: " + ", ".join(f"{k.split()[0]} {d['200']:+.1f}/{d['6']:+.1f}" for k, d in pb.items()))
    rep["criterion"][key] = {"b": gr["fit"]["b"], "break_even": be, "energy_change_8_vs_50": ch, "per_byte_error": pb}
json.dump(rep, open(os.path.join(GMX, "gmx_report.json"), "w"), indent=1)
