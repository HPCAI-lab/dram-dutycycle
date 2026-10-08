#!/usr/bin/env bash
# Build the GEMM-AI workload exactly as specified, record every data access with
# valgrind lackey, and keep
#   (a) the first BRIDGE_N records of the whole program (bridge check), and
#   (b) windows of records issued by the gemm kernel, centred in each of its
#       six calls, plus a 1,000,000-record window centred in call 5.
# A full run takes about 45 minutes on one core and writes about 1 GB.
#
# Environment (all optional):
#   TRACE_OUT          output directory (default build/traces)
#   BRIDGE_N           records in bridge.pin (default 40000050; 0 disables it)
#   WINDOWS            which windows to cut (default: all seven, names as below)
#   LACKEY_STOP_EARLY  1 ends tracing as soon as every requested output is
#                      complete (the summary then reports "complete": false)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${TRACE_OUT:-$ROOT/build/traces}"
BRIDGE_N="${BRIDGE_N:-40000050}"
WINDOWS="${WINDOWS-w1_q w2_k w3_v w4_o w5_ff1 w6_ff2 w5_ff1_1M}"
WIN=100000
BIG=1000000
mkdir -p "$OUT" "$ROOT/build/bin"
OUT="$(cd "$OUT" && pwd)"
BIN="$ROOT/build/bin/gemm_transformer"
FILTER="$ROOT/build/bin/lackey_filter"

gcc -O1 -fno-tree-vectorize -o "$BIN" "$ROOT/workload/gemm_transformer.c" -lm
gcc -O2 -Wall -o "$FILTER" "$ROOT/tracing/lackey_filter.c"

# valgrind loads a position-independent amd64 executable at 0x108000; the gemm
# kernel is the address range from its symbol to the next text symbol.
BASE=$((0x108000))
read -r OFF_LO OFF_HI < <(nm -n "$BIN" | awk '$2 ~ /^[tT]$/ {a[n]=$1; s[n]=$3; n++}
    END {for (i = 0; i < n; i++) if (s[i] == "gemm") {print a[i], a[i+1]; exit}}')
[ -n "${OFF_LO:-}" ] || { echo "gemm symbol not found" >&2; exit 1; }
LO=$(printf '%x' $((BASE + 0x$OFF_LO))); HI=$(printf '%x' $((BASE + 0x$OFF_HI)))
echo "gemm kernel at 0x$LO to 0x$HI"

# One call of gemm(m, k, n) issues 2mkn + mn + 3m + 14 records with gcc 13.3 at
# -O1 (loads of both operands, stores of the result, and the stack accesses of
# the prologue, epilogue and row loop; measured counts are in
# data/trace_summary.json). The offset below places each window within one
# record of the centre of its call.
offset() { local m=$1 k=$2 n=$3 len=$4; echo $(( (2*m*k*n + m*n + 3*m + 13) / 2 - len / 2 )); }
ARGS=()
for w in $WINDOWS; do
    case "$w" in
        w1_q)      ARGS+=(1 "$(offset 128 512 512 $WIN)"  $WIN "$OUT/$w.pin") ;;
        w2_k)      ARGS+=(2 "$(offset 128 512 512 $WIN)"  $WIN "$OUT/$w.pin") ;;
        w3_v)      ARGS+=(3 "$(offset 128 512 512 $WIN)"  $WIN "$OUT/$w.pin") ;;
        w4_o)      ARGS+=(4 "$(offset 128 512 512 $WIN)"  $WIN "$OUT/$w.pin") ;;
        w5_ff1)    ARGS+=(5 "$(offset 128 512 2048 $WIN)" $WIN "$OUT/$w.pin") ;;
        w6_ff2)    ARGS+=(6 "$(offset 128 2048 512 $WIN)" $WIN "$OUT/$w.pin") ;;
        w5_ff1_1M) ARGS+=(5 "$(offset 128 512 2048 $BIG)" $BIG "$OUT/$w.pin") ;;
        *) echo "unknown window $w" >&2; exit 2 ;;
    esac
done

# What the program does before main, and so the bridge segment and every stack
# address, depends on the environment and on the path of the executable. The
# windows are counted per gemm call and do not depend on it, apart from their
# few stack-address records; tracing/compare_traces.py compares a regenerated
# window with a stored one on that basis.
set +e
valgrind --tool=lackey --trace-mem=yes --log-fd=9 "$BIN" 9>&1 1>"$OUT/program_stdout.txt" 2>/dev/null | \
  "$FILTER" "$LO" "$HI" "$LO" "$BRIDGE_N" "$OUT/bridge.pin" "$OUT/trace_summary.json" ${ARGS[@]+"${ARGS[@]}"}
st=("${PIPESTATUS[@]}")
set -e
[ "${st[1]}" -eq 0 ] || { echo "lackey_filter failed (status ${st[1]})" >&2; exit 1; }
if [ "${st[0]}" -ne 0 ] && [ "${LACKEY_STOP_EARLY:-0}" != 1 ]; then
    echo "valgrind failed (status ${st[0]})" >&2; exit 1
fi
echo "trace outputs in $OUT"
