"""Figure: energy per completed request against offered load for one window of
each access pattern, fitted model, held-out and extrapolated points, and the
minimum of energy times mean read latency. Values are read from the run
directories and the fits in summary.json.
Usage: python3 analysis/plot_sweep.py [results_dir]"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_run import load
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results")
FITS = json.load(open(os.path.join(RES, "summary.json")))["fits"]
os.makedirs(os.path.join(RES, "figures"), exist_ok=True)
TRAIN, HOLD, EXTRA = [100, 50, 25, 12, 6], [10, 9, 8, 7], [200, 5]
SER = [("w1_q", "Attention projection window (85% row-buffer hits)", "#2a78d6", "-"),
       ("w5_ff1", "Feed-forward expansion window (48% row-buffer hits)", "#eb6834", "--")]
plt.rcParams.update({"font.family": "serif", "font.size": 7.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})
fig, ax = plt.subplots(figsize=(3.5, 2.75))
ink, muted, grid = "#0b0b0b", "#52514e", "#e4e3df"
for w, label, col, ls in SER:
    f = FITS[f"{w}/uniform"]; a, b = f["intercept"], f["slope"]
    pts = {s: load(os.path.join(RES, "runs", w, "uniform", f"s{s:03d}")) for s in TRAIN + HOLD + EXTRA}
    order = sorted(pts, reverse=True)
    x = lambda s: pts[s]["bw_nominal"]
    ax.plot([x(s) for s in order], [a + b * pts[s]["c_per_req"] for s in order], ls, color=col, lw=1.0, zorder=1)
    ax.scatter([x(s) for s in TRAIN], [pts[s]["E_per_req"] for s in TRAIN], s=22, color=col, edgecolor="white", lw=0.6, zorder=3)
    ax.scatter([x(s) for s in HOLD], [pts[s]["E_per_req"] for s in HOLD], s=22, facecolor="white", edgecolor=col, lw=1.0, zorder=3)
    ax.scatter([x(s) for s in EXTRA], [pts[s]["E_per_req"] for s in EXTRA], s=26, marker="^", color=col, edgecolor="white", lw=0.6, zorder=3)
    edp = min((s for s in TRAIN + HOLD), key=lambda s: pts[s]["E_per_req"] * pts[s]["lat_acc_mean"])
    ax.scatter([x(edp)], [pts[edp]["E_per_req"]], s=95, marker="s", facecolor="none", edgecolor=ink, lw=0.8, zorder=4)
ax.axvline(25.4, color=muted, lw=0.6, ls=":", zorder=0)
ax.text(25.4 * 0.93, 1.25e4, "channel peak", ha="right", va="bottom", rotation=90, color=muted, fontsize=6.3)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Offered load (GB/s)", color=ink)
ax.set_ylabel("Energy per completed request\n(simulator units)", color=ink)
ax.grid(True, which="major", color=grid, lw=0.5); ax.set_axisbelow(True)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
from matplotlib.lines import Line2D
h = [Line2D([], [], color=c, ls=l, lw=1.0, marker="o", ms=4, mec="white", label=lab.split(" (")[0].replace(" window", ""))
     for _, lab, c, l in SER]
h += [Line2D([], [], ls="", marker="o", ms=4, color="#777", mec="white", label="fitted point"),
      Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec="#777", label="held out"),
      Line2D([], [], ls="", marker="^", ms=4.5, color="#777", mec="white", label="extrapolated"),
      Line2D([], [], ls="", marker="s", ms=6, mfc="none", mec=ink, label="min. energy-latency product")]
ax.legend(handles=h, loc="upper right", bbox_to_anchor=(0.90, 1.0), frameon=False, fontsize=6.3, handlelength=1.8, labelspacing=0.3)
fig.tight_layout(pad=0.3)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(RES, "figures", f"fig_fixedwork_sweep.{ext}"), dpi=300, metadata={"Creator": None, "Producer": None} if ext == "pdf" else None)
print("saved", os.path.join(RES, "figures"))
