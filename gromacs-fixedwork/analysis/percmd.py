"""Refit the fixed-work model per DRAM command instead of per request.
x = cycles per command (read + write commands), y = energy per command.
Same strides, same relative-error weights as the paper's fit."""
import json, sys, numpy as np
sys.path.insert(0, "repo_gmx/analysis"); from read_run import load
TR, HO, EX = [100, 50, 25, 12, 6], [10, 9, 8, 7], [200, 5]
BR = (880.98, 1168.98)
def run(res, wins):
    for w in wins:
        P = {s: load(f"{res}/runs/{w}/uniform/s{s:03d}") for s in TR + HO + EX}
        for mode in ("request", "command"):
            n = lambda p: p["completed"] if mode == "request" else p["read_cmds"] + p["write_cmds"]
            x = {s: P[s]["cycles"] / n(P[s]) for s in P}; y = {s: P[s]["E_total"] / n(P[s]) for s in P}
            b, a = np.polyfit([x[s] for s in TR], [y[s] for s in TR], 1, w=[1 / y[s] for s in TR])
            e = {s: 100 * abs(a + b * x[s] - y[s]) / y[s] for s in P}
            print(f"{w:12s} per {mode:8s} a={a:7.1f} b={b:7.1f} inside={BR[0] <= b <= BR[1]!s:5s} "
                  f"held MAPE {np.mean([e[s] for s in HO]):.2f}% max {max(e[s] for s in HO):.2f}%  s200 {e[200]:.2f}%  s5 {e[5]:.2f}%")
print("== GROMACS"); run("repo_gmx/results", ["g1_nb25", "g2_nb50", "g3_nb75", "g4_spread", "g5_fft", "g6_gather", "g2_nb50_1M"])
print("== GEMM-AI"); run("repo/results", ["w1_q", "w5_ff1"])
