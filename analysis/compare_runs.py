"""Compares two run trees (for example a fresh build/runs against results/runs):
every counter and energy field of every common run must be identical.
Usage: python3 analysis/compare_runs.py <runs_a> <runs_b>"""
import glob, json, os, sys
a, b = sys.argv[1], sys.argv[2]
common = sorted(os.path.relpath(p, a) for p in glob.glob(os.path.join(a, "*", "*", "s*", "dramsim3.json"))
                if os.path.exists(os.path.join(b, os.path.relpath(p, a))))
if not common:
    sys.exit("no runs in common")
diff = []
for rel in common:
    x = json.load(open(os.path.join(a, rel)))["0"]; y = json.load(open(os.path.join(b, rel)))["0"]
    keys = [k for k in x if k.endswith("energy") or k.startswith("num_") or k in ("rank_active_cycles", "read_latency")]
    if any(x[k] != y.get(k) for k in keys):
        diff.append(os.path.dirname(rel))
print(f"{len(common)} runs compared; identical: {len(common) - len(diff)}; different: {len(diff)}")
for d in diff:
    print("  differs:", d)
sys.exit(1 if diff else 0)
