#!/usr/bin/env bash
# Fetch DRAMsim3 at the pinned commit and build it twice:
#   build/bin/dramsim3main-stock   unmodified, from build/DRAMsim3
#   build/bin/dramsim3main         with the fixed-work patch, from build/DRAMsim3-fixedwork
# DRAMsim3's CMake build writes libdramsim3.so into the source tree and records
# that directory in the executable's run path, so each build needs its own tree
# and both trees must stay in place. build/DRAMsim3 also provides configs/.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
SRC="${DRAMSIM3_DIR:-$ROOT/build/DRAMsim3}"
FW="$ROOT/build/DRAMsim3-fixedwork"
COMMIT="$(cat "$HERE/DRAMSIM3_COMMIT")"
mkdir -p "$ROOT/build/bin"

# CMake 4 no longer accepts DRAMsim3's cmake_minimum_required line on its own.
POLICY=()
[ "$(cmake --version | awk 'NR == 1 {split($3, v, "."); print v[1]}')" -ge 4 ] && POLICY=(-DCMAKE_POLICY_VERSION_MINIMUM=3.5)

build() {   # build <source tree> <destination of the executable>
    cmake -S "$1" -B "$1/build" -DCMAKE_BUILD_TYPE=Release ${POLICY[@]+"${POLICY[@]}"} > /dev/null
    cmake --build "$1/build" --target dramsim3main -j"$(nproc)" > /dev/null
    cp "$1/build/dramsim3main" "$2"
}

[ -d "$SRC/.git" ] || git clone https://github.com/umd-memsys/DRAMsim3.git "$SRC"
git -C "$SRC" checkout -q "$COMMIT"
git -C "$SRC" diff --quiet || { echo "$SRC has local changes; it must be the unmodified tree" >&2; exit 1; }
build "$SRC" "$ROOT/build/bin/dramsim3main-stock"

rm -rf "$FW"
git clone -q "$SRC" "$FW"
git -C "$FW" checkout -q "$COMMIT"
git -C "$FW" apply "$HERE/fixedwork.patch"
build "$FW" "$ROOT/build/bin/dramsim3main"

echo "built build/bin/dramsim3main-stock and build/bin/dramsim3main (DRAMsim3 $COMMIT)"
