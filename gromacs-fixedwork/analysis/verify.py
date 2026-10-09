"""Independent recomputation of headline GROMACS numbers from raw dramsim3.json only."""
import json, glob, numpy as np
R = "repo_gmx/results/runs"
def raw(w, a, s):
    c = json.load(open(f"{R}/{w}/{a}/s{s:03d}/dramsim3.json"))["0"]
    n = c["num_reads_done"] + c["num_writes_done"]
    resid = sum(c["act_stb_energy"].values()) + sum(c["pre_stb_energy"].values()) + sum(c["sref_energy"].values()) + c["ref_energy"] + c["refb_energy"]
    return n, c["num_cycles"], c["total_energy"], resid / c["total_energy"]
W = ["g1_nb25", "g2_nb50", "g3_nb75", "g4_spread", "g5_fft", "g6_gather", "g2_nb50_1M"]
PS, SR, AS = 2*1.2*37*8 + 2*1.2*(250-52)*560*8/12480, 2*1.2*30*8, 2*1.2*52*8 + 2*1.2*(250-52)*560*8/12480
hold_all = []; b_all = []; ext = []
for w in W:
    d = {s: raw(w, "uniform", s) for s in [200,100,50,25,12,10,9,8,7,6,5]}
    x = lambda s: d[s][1] / d[s][0]; y = lambda s: d[s][2] / d[s][0]
    tr = [100, 50, 25, 12, 6]
    b, a = np.polyfit([x(s) for s in tr], [y(s) for s in tr], 1, w=[1/y(s) for s in tr])
    e = lambda s: 100 * abs(a + b*x(s) - y(s)) / y(s)
    h = np.mean([e(s) for s in (10, 9, 8, 7)]); hold_all.append(h); b_all.append(b); ext += [e(200), e(5)]
    ch = lambda P: 100 * ((d[8][2] + (d[50][1] - d[8][1]) * P) / d[50][2] - 1)
    print(f"{w:12s} a={a:7.1f} b={b:7.1f} held={h:.2f}% e200={e(200):.2f}% e5={e(5):.2f}%  resid s50={100*d[50][3]:.2f}% s6={100*d[6][3]:.2f}%  "
          f"dE(8 vs 50): none {ch(0):+.1f} SR {ch(SR):+.1f} PS {ch(PS):+.1f} AS {ch(AS):+.1f}")
print(f"held-out MAPE range {min(hold_all):.2f}-{max(hold_all):.2f}; b range {min(b_all):.1f}-{max(b_all):.1f}; max extrapolation {max(ext):.2f}")
