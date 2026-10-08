#!/usr/bin/env bash
# Bridge check: the first 40,000,050 records of the regenerated trace, converted
# at stride 50 and run under the original 2-billion-cycle budget (no -f), so that
# exactly 40,000,000 records are submitted, as in the matched-load runs.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${1:-$ROOT/build/traces/bridge.pin}"
OUT="${2:-$ROOT/build/bridge}"
CFG="$ROOT/build/DRAMsim3/configs/DDR4_8Gb_x8_3200.ini"
mkdir -p "$OUT"
awk -f "$ROOT/conversion/convert_pin_to_dramsim.awk" -v stride=50 "$SRC" > "$OUT/bridge.trace"
"$ROOT/build/bin/dramsim3main" "$CFG" -c 2000000000 -t "$OUT/bridge.trace" -o "$OUT" > /dev/null
rm -f "$OUT/bridge.trace" "$OUT/dramsim3epoch.json"
python3 "$ROOT/analysis/compare_bridge.py" "$OUT"
