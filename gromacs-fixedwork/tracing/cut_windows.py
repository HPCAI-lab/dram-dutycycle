"""Cut fixed-work sweep windows from one traced GROMACS MD step.

Reads the binary step file and profile written by gmx_filter, ranks functions
by the records their own instructions issue in that step, and for each chosen
function writes the LENGTH consecutive records of that function centred in its
records for the step (pinatrace-style text, as the GEMM-AI windows).

Usage:
  python3 cut_windows.py profile.json step.bin outdir --list
  python3 cut_windows.py profile.json step.bin outdir name=fid[:length] ...
"""
import json, os, sys
import numpy as np

REC = np.dtype([("addr", "<u8"), ("ipoff", "<u4"), ("fid", "<u2"), ("op", "u1"), ("size", "u1")])
prof_path, step_path, outdir = sys.argv[1:4]
prof = json.load(open(prof_path))
base = int(prof["base"], 16)
col = prof["target_step"] + 1
recs = np.memmap(step_path, dtype=REC, mode="r")
n = len(recs)
names = {f["fid"]: f["name"] for f in prof["functions"]}
fids, counts = np.unique(np.asarray(recs["fid"]), return_counts=True)
order = np.argsort(-counts)
if len(sys.argv) > 4 and sys.argv[4] == "--list":
    print(f"step {prof['target_step']}: {n:,} records; step record counts {prof['records_per_step']}")
    cum = 0
    for i in order[:40]:
        cum += counts[i]
        print(f"{fids[i]:6d} {counts[i]:>14,} {100*counts[i]/n:6.2f}% cum {100*cum/n:6.2f}%  {names.get(int(fids[i]), '?')[:150]}")
    sys.exit()
os.makedirs(outdir, exist_ok=True)
meta = {"step": prof["target_step"], "step_records": int(n), "windows": {}}
fidcol = recs["fid"]
cnt_by_fid = dict(zip(fids.tolist(), counts.tolist()))
for spec in sys.argv[4:]:
    # name=fid[:length[:position]]; position is the fraction of the function's
    # records at which the window is centred (default 0.5)
    name, rest = spec.split("=")
    fid, length, pos = (rest.split(":") + ["100000", "0.5"][len(rest.split(":")) - 1:])[:3]
    length, pos = int(length), float(pos)
    if fid == "all":     # contiguous records of every function, centred in the step
        fid, m = -1, n
        start = int(m * pos) - length // 2
        sel = recs[start:start + length]
    else:
        fid = int(fid)
        CH = 50_000_000
        m = int(cnt_by_fid[fid])
        start = int(m * pos) - length // 2
        parts, seen = [], 0          # occurrence ranks [start, start + length), chunk by chunk
        for c0 in range(0, n, CH):
            loc = np.flatnonzero(fidcol[c0:c0 + CH] == fid)
            lo, hi = max(start - seen, 0), min(start + length - seen, len(loc))
            if lo < hi:
                parts.append(np.asarray(recs[c0 + loc[lo:hi]]))
            seen += len(loc)
            if seen >= start + length:
                break
        sel = np.concatenate(parts)
    with open(os.path.join(outdir, f"{name}.pin"), "w") as f:
        for r in sel:
            ip = base + int(r["ipoff"]) if r["ipoff"] != 0xFFFFFFFF else 0
            f.write(f"0x{ip:x}: {chr(r['op'])} 0x{int(r['addr']):x} {int(r['size'])}\n")
    meta["windows"][name] = {"fid": fid, "function": names.get(fid), "function_records_in_step": int(m),
                             "share_of_step": m / n, "position": pos, "offset": int(start), "length": int(len(sel)),
                             "reads": int((sel["op"] == ord("R")).sum()), "writes": int((sel["op"] == ord("W")).sum())}
    print(name, meta["windows"][name])
json.dump(meta, open(os.path.join(outdir, "windows_meta.json"), "w"), indent=1)
