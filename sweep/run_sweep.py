"""Fixed-work offered-load sweep on the patched DRAMsim3.

For every window, arrival process and stride, the window is converted to a
DRAMsim3 trace and simulated until every request has completed (-f). Each run
directory receives the simulator outputs and a meta.json recording the exact
commands and the SHA-256 of the converted trace and of the simulator binary.
Existing run directories are skipped.

Example:
    python3 sweep/run_sweep.py --out build/runs
"""
import argparse, gzip, hashlib, itertools, json, os, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WINDOWS = ["w1_q", "w2_k", "w3_v", "w4_o", "w5_ff1", "w6_ff2"]
STRIDES = [200, 100, 50, 25, 12, 10, 9, 8, 7, 6, 5]

def rel(p):
    return os.path.relpath(p, ROOT)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def window_file(windows_dir, w):
    for name in (f"{w}.pin", f"{w}.pin.gz"):
        p = os.path.join(windows_dir, name)
        if os.path.exists(p):
            return p
    sys.exit(f"window {w} not found in {windows_dir}")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--windows-dir", default=os.path.join(ROOT, "data", "windows"))
    ap.add_argument("--out", default=os.path.join(ROOT, "build", "runs"))
    ap.add_argument("--sim", default=os.path.join(ROOT, "build", "bin", "dramsim3main"))
    ap.add_argument("--config", default=os.path.join(ROOT, "build", "DRAMsim3", "configs", "DDR4_8Gb_x8_3200.ini"))
    ap.add_argument("--windows", default=",".join(WINDOWS))
    ap.add_argument("--strides", default=",".join(map(str, STRIDES)))
    ap.add_argument("--arrivals", default="uniform,poisson")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args()
    awk = os.path.join(ROOT, "conversion", "convert_pin_to_dramsim.awk")
    pois = os.path.join(ROOT, "conversion", "convert_poisson.py")
    seeds = json.load(open(os.path.join(ROOT, "data", "poisson_seeds.json")))["base"]
    for p in (a.sim, a.config):
        if not os.path.exists(p):
            sys.exit(f"missing {p}; run simulator/build_dramsim3.sh first")
    sim_sha = sha256(a.sim)
    jobs = list(itertools.product(a.windows.split(","), a.arrivals.split(","), [int(s) for s in a.strides.split(",")]))

    def one(job):
        w, arrival, s = job
        d = os.path.join(a.out, w, arrival, f"s{s:03d}")
        if os.path.exists(os.path.join(d, "meta.json")):
            return d + " (exists)"
        os.makedirs(d, exist_ok=True)
        src = window_file(a.windows_dir, w)
        trace = os.path.join(d, "trace.txt")
        if arrival == "uniform":
            if src.endswith(".gz"):
                plain = os.path.join(d, "window.pin")
                with gzip.open(src, "rb") as fi, open(plain, "wb") as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                plain = src
            conv = ["awk", "-f", awk, "-v", f"stride={s}", plain]
        else:
            base = w.replace("_1M", "")
            conv = ["python3", pois, src, str(s), str(seeds[base] + s)]
        with open(trace, "w") as f:
            subprocess.run(conv, stdout=f, stderr=subprocess.DEVNULL, check=True)
        n = sum(1 for _ in open(trace))
        cap = 40 * n * s + 10_000_000
        cmd = [a.sim, a.config, "-c", str(cap), "-f", "-t", trace, "-o", d]
        r = subprocess.run(cmd, capture_output=True, text=True)
        meta = {"window": w, "arrival": arrival, "stride": s, "records": n,
                "trace_sha256": sha256(trace), "sim_sha256": sim_sha,
                "cmd": [rel(x) if os.path.isabs(x) else x for x in cmd],
                "converter": [rel(x) if os.path.isabs(x) else x for x in conv],
                "returncode": r.returncode, "stderr_tail": r.stderr[-500:]}
        for tmp in (trace, os.path.join(d, "window.pin")):
            if os.path.exists(tmp):
                os.remove(tmp)
        with open(os.path.join(d, "meta.json"), "w") as f:
            json.dump(meta, f, indent=1)
        return d

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for d in ex.map(one, jobs):
            print("done", rel(d) if os.path.isabs(d.split(" ")[0]) else d, flush=True)

if __name__ == "__main__":
    main()
