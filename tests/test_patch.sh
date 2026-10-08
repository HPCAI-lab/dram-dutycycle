#!/usr/bin/env bash
# Tests of the fixed-work changes against the unmodified simulator, on small
# constructed traces (written to a temporary directory and removed afterwards).
#  1. end of trace: the last record is submitted once, tail writes reach the
#     device, and the run stops when the work is done
#  2. before end of trace the two builds produce identical statistics
#  3. the write-blocked-by-read deadlock: the stock build stalls, the patched
#     build completes and counts one use of the fix
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CFG="$ROOT/build/DRAMsim3/configs/DDR4_8Gb_x8_3200.ini"
STOCK="$ROOT/build/bin/dramsim3main-stock"; FW="$ROOT/build/bin/dramsim3main"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
random.seed(1)
with open("t1.trace", "w") as f:
    for i in range(1000):
        op = "WRITE" if i >= 995 or random.random() < 0.1 else "READ"
        f.write(f"0x{random.randrange(0, 1 << 33) & ~63:x} {op} {i * 50}\n")
rows = [r << 18 for r in range(1, 41)]          # one bank, 40 different rows
lines = [(a, "READ", c) for c, a in enumerate(rows)]
lines.append((rows[30], "WRITE", len(lines)))  # write to the address of a queued read
for k in range(60):                            # then fill the write buffer
    lines.append((((k % 4) << 13) | ((k // 4 % 4) << 15) | (1 << 17) | ((k + 100) << 18), "WRITE", len(lines)))
with open("t3.trace", "w") as f:
    for a, op, c in lines:
        f.write(f"0x{a:x} {op} {c}\n")
PY
mkdir s1 f1 s2 f2 s3 f3
"$STOCK" "$CFG" -c 200000 -t t1.trace -o s1 > /dev/null
"$FW" "$CFG" -c 200000 -f -t t1.trace -o f1 > /dev/null
"$STOCK" "$CFG" -c 30000 -t t1.trace -o s2 > /dev/null
"$FW" "$CFG" -c 30000 -t t1.trace -o f2 > /dev/null
"$STOCK" "$CFG" -c 300000 -t t3.trace -o s3 > /dev/null
"$FW" "$CFG" -c 300000 -f -t t3.trace -o f3 > /dev/null
python3 - <<'PY'
import json, filecmp, sys
J = lambda d: json.load(open(f"{d}/dramsim3.json"))["0"]
C = lambda d: json.load(open(f"{d}/cpu_stats.json"))
reads1 = sum(1 for l in open("t1.trace") if " READ " in l); writes1 = 1000 - reads1
s1, f1, f3, s3, c1, c3 = J("s1"), J("f1"), J("f3"), J("s3"), C("f1"), C("f3")
checks = [
 ("1a stock submits one extra request at end of trace", s1["num_reads_done"] == reads1 + 1),
 ("1b stock leaves tail writes unissued", s1["num_write_cmds"] < writes1),
 ("1c patched completes exactly the trace", f1["num_reads_done"] == reads1 and f1["num_writes_done"] == writes1),
 ("1d patched issues every write to the device", f1["num_write_cmds"] == writes1),
 ("1e patched stops at completion", c1["memory_idle"] and c1["trace_exhausted"] and f1["num_cycles"] < 60000),
 ("2  identical statistics before end of trace", filecmp.cmp("s2/dramsim3.json", "f2/dramsim3.json", shallow=False)),
 ("3a stock stalls on the constructed deadlock", s3["num_reads_done"] < 40),
 ("3b patched completes it", f3["num_reads_done"] == 40 and f3["num_writes_done"] == 61 and c3["memory_idle"]),
 ("3c patched counts the fix", c3["rw_fallbacks"] >= 1)]
for name, ok in checks:
    print(("PASS " if ok else "FAIL ") + name)
sys.exit(0 if all(ok for _, ok in checks) else 1)
PY
