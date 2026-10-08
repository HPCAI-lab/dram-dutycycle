# Fixed-work offered-load sweep of device-level DRAM energy

This repository contains the trace generation, simulator modifications, sweep
driver, raw simulator outputs and analysis behind the offered-load sweep in the
paper *Predicting Device-Level DRAM Energy from Operating Duty Cycle: A
Fixed-Work Model for HPC Memory Systems*.

A memory trace of a transformer-layer general matrix multiplication (GEMM)
workload is presented to a simulated double data rate fourth-generation (DDR4)
channel at eleven offered loads, under two arrival processes. Every run
completes a fixed amount of work: the run ends when every request in the trace
has been served and every buffered write has reached the device. Energy per
completed request is then fitted against cycles per completed request, and the
fit is tested on held-out and extrapolated loads.

Everything needed to regenerate the published numbers is here, at three levels
of depth:

| Starting point | Command | Time on two cores |
|---|---|---|
| Stored simulator outputs (`results/runs`) | `make analysis figure` | seconds |
| Stored trace windows (`data/windows`) | `make simulator test sweep compare` | about 4 minutes |
| Nothing but the workload source | `make trace check-trace bridge`, then a sweep over `build/traces` (Section 9) | about one hour |

## 1. Repository layout

```
workload/gemm_transformer.c        the GEMM-AI workload (one transformer layer)
tracing/lackey_filter.c            valgrind lackey output to pinatrace-style records
tracing/trace_gemm.sh              builds the workload, traces it, cuts the windows
tracing/compare_traces.py          regenerated window against a stored one, record by record
simulator/DRAMSIM3_COMMIT          pinned DRAMsim3 commit
simulator/fixedwork.patch          the four simulator changes (Section 4)
simulator/build_dramsim3.sh        builds the unmodified and the patched simulator
conversion/convert_pin_to_dramsim.awk   evenly spaced arrivals
conversion/convert_poisson.py      randomly spaced (exponential-gap) arrivals, seeded
sweep/run_sweep.py                 runs every window, arrival process and stride
sweep/run_bridge.sh                bridge check against the matched-load tables
analysis/read_run.py               reads one run directory, checks its energy sum
analysis/analyze_sweep.py          fits, held-out and extrapolation errors, EDP minima
analysis/plot_sweep.py             the sweep figure
analysis/compare_bridge.py         bridge run against the published table rows
analysis/compare_runs.py           two run trees must agree counter for counter
tests/test_patch.sh                tests of the simulator changes on constructed traces
data/windows/*.pin.gz              the seven trace windows used in the paper
data/SHA256SUMS                    checksums of the uncompressed windows
data/trace_summary.json            record counts of the full trace, per gemm call
data/program_stdout.txt            the workload's checksum line
data/poisson_seeds.json            seeds of the random-arrival converter
results/runs/                      simulator outputs of all 143 sweep runs
results/bridge/                    simulator outputs of the bridge run and its comparison
results/summary.json, all_runs.csv, figures/   analysis outputs
```

Everything the scripts build or download goes to `build/`, which is not tracked.

## 2. Requirements

- Linux on x86-64. The kernel's address range is computed from valgrind's load
  address for an amd64 position-independent executable.
- gcc, CMake, git, make and awk. Tested with CMake 3.28; with CMake 4 the build
  script adds the policy setting that DRAMsim3's minimum-version line needs.
- valgrind 3.22 or later (only for `make trace`).
- Python 3 with numpy and matplotlib (`pip install -r requirements.txt`); tested
  with Python 3.13, numpy 2.5 and matplotlib 3.11.

The stored traces were produced with gcc 13.3.0, glibc 2.39 and valgrind 3.22.0
on Ubuntu 24.04. A different compiler or C library changes the instruction
stream and the heap layout, so a regenerated trace can be compared with the
stored one only on the same toolchain. Even then, stack addresses shift with the
size of the environment and the path of the executable, and so does the number
of records issued before `main`. The windows are counted from the start of each
`gemm` call, so only their stack-address records are affected:
`make check-trace` (`tracing/compare_traces.py`) requires every other field of
every record to agree. A regenerated `w1_q` window differed from the stored one
in 3 of its 100,000 records, all stack addresses, and gave simulator counters
identical to the stored runs at every stride under both arrival processes.

The simulator is deterministic: a fresh sweep over the stored windows reproduces
the stored runs counter for counter (`make compare`).

## 3. Pipeline

1. **Workload.** `workload/gemm_transformer.c` computes one transformer layer
   (model dimension 512, feed-forward dimension 2048, sequence length 128) with a
   single scalar `gemm` routine: four attention projections of shape
   128 by 512 by 512, the feed-forward expansion (128 by 512 by 2048) and the
   contraction (128 by 2048 by 512). It is built with
   `gcc -O1 -fno-tree-vectorize` and prints one checksum line, which is the same
   natively and under valgrind (`data/program_stdout.txt`).
2. **Tracing.** valgrind's lackey tool (`--trace-mem=yes`) records every data
   access with no timing information. `lackey_filter` rewrites each access as a
   pinatrace-style record, `<instruction address>: R|W <data address> <size>`;
   a modify access becomes a read followed by a write. The full program issues
   809,771,252 records.
3. **Windows.** Each window is 100,000 consecutive records issued by the
   instructions of `gemm`, counted from the start of one call and centred in it
   (Section 5). A 1,000,000-record window centred in the feed-forward expansion
   checks that the window length does not matter.
4. **Arrival times.** A window carries no timing, so the converter supplies it.
   The evenly spaced converter issues one request every *stride* cycles. The
   random converter draws exponential gaps whose mean is the stride, so both
   processes offer the same average load at the same stride.
5. **Simulation.** The patched DRAMsim3 runs each converted trace with `-f`, so
   the run stops when the work is done (Section 4).
6. **Analysis.** `analyze_sweep.py` fits the model, scores it, finds the
   energy-delay-product (EDP) minima and measures how much the energy breakdown
   differs between windows at each load.

## 4. Simulator and its four changes

DRAMsim3 at commit `29817593b3389f1337235d63cac515024ab8fd6e`, with its stock
`configs/DDR4_8Gb_x8_3200.ini`: one 64-bit channel of two ranks, eight x8
devices per rank, open-page policy, clock period 0.63 nanoseconds, peak
bandwidth 25.4 gigabytes per second. `simulator/build_dramsim3.sh` builds the
unmodified simulator (`build/bin/dramsim3main-stock`, from `build/DRAMsim3`) and
the patched one (`build/bin/dramsim3main`, from `build/DRAMsim3-fixedwork`). Each
executable loads `libdramsim3.so` from its own source tree, so the two builds
use separate trees, and both trees must stay in place.

The unmodified simulator runs for a fixed number of cycles, so a run can end
with requests still queued, or keep idling after its trace has ended; either
distorts energy per completed request. The patch adds a fixed-work mode
(change 2) and fixes three defects that such a run would otherwise expose:

1. **The last trace record is submitted once.** The stock trace reader
   resubmits the final record after the failed read at end of file.
2. **Fixed-work stop.** With `-f`, the run ends when the trace is exhausted and
   every queue and buffer is empty. `-c` remains as a safety cap.
3. **End-of-trace write drain.** The stock controller drains the write buffer
   only when it holds more than eight writes, so up to eight trailing writes
   never reach the device and are never charged. Once the trace is exhausted,
   the buffer is drained whenever the command queue is empty.
4. **Write blocked by a read to the same address.** When the controller is
   draining writes and the first write depends on an older pending read to the
   same address, the stock scheduler issues nothing that cycle, and the
   condition can persist indefinitely. The patch schedules a read instead, and
   counts every use in `cpu_stats.json` (`rw_fallbacks`).

The patch also writes `cpu_stats.json` for every run: requests added and
completed, whether the trace was exhausted and the memory system idle at the
end, stall cycles, uses of change 4, and read latency measured from the
request's arrival time in the trace.

`tests/test_patch.sh` checks these changes on small constructed traces (not on
any result):

| Check | What it establishes |
|---|---|
| 1a, 1b | The stock build submits one extra read at end of trace and leaves tail writes unissued |
| 1c to 1e | The patched build completes exactly the trace, issues every write and stops at completion |
| 2 | Stopped before the end of the trace, the two builds write byte-identical statistics |
| 3a to 3c | On a constructed write-after-read conflict the stock build stalls, the patched build completes and records the fix |

In every sweep run and in the bridge run, change 4 was used zero times, so up
to the end of each trace the patched simulator behaved exactly as the stock
one. `analyze_sweep.py` treats any nonzero count as an integrity failure.

## 5. Trace windows and the bridge check

`tracing/trace_gemm.sh` locates `gemm` with `nm` and adds valgrind's load
address (0x108000) to obtain its instruction range. It counts records from the
start of each call. One call of `gemm(m, k, n)` issues 2mkn + mn + 3m + 14
records with this build; each window starts at half of 2mkn + mn + 3m + 13 less
half the window length, which places it within one record of the centre of its
call.

| Window | gemm call | Records in call | Offset of first record | Length |
|---|---|---|---|---|
| `w1_q`, `w2_k`, `w3_v`, `w4_o` | 1 to 4 (attention projections) | 67,174,798 | 33,537,398 | 100,000 |
| `w5_ff1` | 5 (feed-forward expansion) | 268,697,998 | 134,298,998 | 100,000 |
| `w6_ff2` | 6 (feed-forward contraction) | 268,501,390 | 134,200,694 | 100,000 |
| `w5_ff1_1M` | 5 | 268,697,998 | 133,848,998 | 1,000,000 |

The windows fall into two access patterns. Five of them have read row-buffer
hit rates between 82% and 87%, depending on load; the feed-forward expansion
has between 47% and 48%.

**Bridge check.** To confirm that the regenerated trace stands in for the trace
behind the paper's matched-load tables, `sweep/run_bridge.sh` runs the first
40,000,050 records of the program at stride 50 under the original budget of
2,000,000,000 cycles with the unmodified stop rule, so that 40,000,000 requests
arrive within the budget. `results/bridge/comparison.txt` holds the result:
every energy share is within 0.03 percentage points of the published value,
and the read row-buffer hit rate is 84.96% against 84.86%. The bridge segment
starts before `main`, so a regenerated one differs from the stored one in its
early records (Section 2). One regenerated in a different environment and
directory gave the same request counts, 12 more read commands out of 38.57
million, and the same energy shares and hit rate to two decimal places.

## 6. Sweep design

| Stride (cycles between requests) | 200 | 100 | 50 | 25 | 12 | 10 | 9 | 8 | 7 | 6 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Offered load (gigabytes per second) | 0.51 | 1.02 | 2.03 | 4.06 | 8.47 | 10.16 | 11.29 | 12.70 | 14.51 | 16.93 | 20.32 |
| Percent of channel peak | 2 | 4 | 8 | 16 | 33 | 40 | 44 | 50 | 57 | 67 | 80 |
| Role | extrapolated | fitted | fitted | fitted | fitted | held out | held out | held out | held out | fitted | extrapolated |

Offered load is 64 bytes per request divided by stride times 0.63 nanoseconds.
The held-out strides lie inside the fitted range, so they test interpolation;
strides 200 and 5 lie outside it.

- **Runs:** six kernel windows, two arrival processes and eleven strides give
  132 runs; the 1,000,000-record window adds eleven runs with evenly spaced
  arrivals, for 143 in all.
- **Seeds:** the random converter's seed for window *w* at stride *s* is
  base(*w*) + *s*, with the bases in `data/poisson_seeds.json`. The 1,000,000-record
  window has no random-arrival runs.
- **Cycle cap:** 40 times records times stride plus 10,000,000 cycles, a safety
  limit only. Every run ended by completing its work.
- **Records:** each run directory holds `dramsim3.json` and `dramsim3.txt`
  (simulator statistics), `cpu_stats.json` (Section 4) and `meta.json` (window,
  arrival process, stride, the converter and simulator commands, and SHA-256
  checksums of the converted trace and of the simulator binary). The binary
  checksum identifies the build used for the stored runs; a rebuild on another
  machine will differ there and nowhere else.

## 7. Energy accounting and units

DRAMsim3 accumulates energy as current times voltage times a cycle count and
never multiplies by the clock period. Its energy unit is therefore the
milliwatt-cycle, which is 0.63 picojoules for this device. All energies in this
repository are in simulator units; `read_run.py` also reports picojoules.
Energy shares do not depend on the unit. Absolute energies scale directly with
the configuration's current and voltage values, so they should be read as
model-dependent.

From the configuration file (derived analytically, not measured):

- **Per-command energy:** activate 6,720.0; read 4,454.4; write 3,763.2.
- **Refresh:** 170.6 per cycle on average, two ranks together.
- **Time-proportional energy per cycle:** between 881.0 (both ranks in
  precharge standby) and 1,169.0 (both in active standby), refresh included.
  These two values bound the slope b of the model in Section 8.

`read_run.py` checks in every run that the components sum to the simulator's own
total energy.

## 8. Provenance

The material falls into three classes, which should not be confused.

### 8.1 Designed artifacts (correct by construction or inspection)

- The workload, the trace filter, the two converters, the patch, the sweep
  driver and the analysis scripts.
- The constructed test traces in `tests/test_patch.sh`.
- The slope bounds and per-command energies of Section 7, which follow from the
  configuration file.

### 8.2 Asserted methodology (how the experiment was run)

- Toolchain, tracing and window placement as in Sections 2, 3 and 5.
- Simulator, configuration and patch as in Section 4.
- Strides, arrival processes, seeds and cycle cap as in Section 6.
- The fitted model, energy per completed request = a + b times cycles per
  completed request, fitted by least squares on strides 100, 50, 25, 12 and 6.

### 8.3 Measured results

All in `results/`, regenerated from the run directories by `make analysis`.
For the five windows with 82% to 87% read row-buffer hits (range across the
five) and for the feed-forward expansion, with evenly spaced arrivals:

| | Five 82% to 87% hit windows | Feed-forward expansion |
|---|---|---|
| Intercept a (simulator units) | 5,859.8 to 5,882.2 | 8,043.2 |
| Slope b (simulator units per cycle) | 1,120.5 to 1,121.1 | 1,148.7 |
| Coefficient of determination over fitted points | 0.999917 to 0.999928 | 0.999995 |
| Leave-one-out error | 1.77% to 1.87% | 0.39% |
| Held-out mean absolute percentage error | 1.66% to 1.75% | 0.36% |
| Held-out maximum error | 2.20% to 2.31% | 0.46% |
| Extrapolation error, stride 200 | 1.60% to 1.65% | 0.89% |
| Extrapolation error, stride 5 | 2.88% to 3.08% | 0.50% |

- **Random arrivals** give held-out errors of 1.66% to 1.78% and 0.34%, and
  extrapolation errors of at most 3.13% and 0.90%. Every slope, under either
  arrival process, lies inside the bounds of Section 7.
- **The errors are systematic.** The model over-predicts at both ends of the
  load range and under-predicts in the middle (strides 50 and 25). For the
  82% to 87% hit windows the error grows steadily toward high load, reaching
  about 3% at stride 5.
- **Window length:** the 1,000,000-record window gives 8,048.9, 1,148.7 and
  0.38% held-out error.
- **Energy-delay product.** With delay taken as mean read latency from
  acceptance by the controller, the minimum falls at stride 10 or 9 for the
  82% to 87% hit windows, with strides 12 to 8 within 1.5% of the minimum for
  evenly spaced arrivals and within 5.9% for random arrivals. For the
  feed-forward expansion it falls at stride 12, with strides 10 and 9 within
  6.1% (evenly spaced) or 8.8% (random). With delay taken as completion time,
  the product keeps falling up to stride 6 in every window and has no interior
  minimum.
- **Energy breakdown across windows.** The spread of the combined
  standby-and-refresh share across the six kernel windows, in percentage
  points, with evenly spaced arrivals:

  | Stride | 200 | 100 | 50 | 25 | 12 | 10 | 9 | 8 | 7 | 6 | 5 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | Spread | 0.91 | 1.77 | 3.36 | 5.75 | 8.57 | 9.10 | 9.36 | 9.59 | 9.77 | 6.15 | 4.64 |

  The two access patterns have nearly the same breakdown at low load; the gap
  grows with load until saturation compresses it.

## 9. Reproducing

```sh
pip install -r requirements.txt

# from the stored simulator outputs
make analysis figure                 # results/summary.json, all_runs.csv, figures/

# from the stored windows
make simulator                       # build/bin/dramsim3main-stock and dramsim3main
make test                            # checks of Section 4
make sweep                           # 143 runs into build/runs
make compare                         # build/runs against results/runs, counter for counter
make analysis figure RESULTS=build   # the same analysis on the fresh runs

# from the workload source
make trace                           # build/traces (about 45 minutes, about 1 GB)
make check-trace                     # regenerated windows against data/windows
make bridge                          # bridge run and comparison (about 10 minutes)
make sweep compare WINDOWS_DIR=build/traces OUT=build/from-trace
make analysis figure RESULTS=build/from-trace
```

`make sweep` skips run directories that already exist, so an interrupted sweep
can be resumed. `make trace` accepts `WINDOWS=w1_q`, `BRIDGE_N=0` and
`LACKEY_STOP_EARLY=1` to cut a single window and stop as soon as it is written;
for `w1_q` this takes about two minutes, and `make check-trace` then checks that
window alone.

## 10. Limitations

- Every processor load and store is presented to the device, with no cache
  between them, as in the paper's matched-load experiments.
- The trace drives the device open-loop; the only back-pressure is a full
  controller queue. Processor stalls are not modelled.
- One transformer layer, one device configuration and one channel.
- At stride 5 the offered load exceeds what the channel completes for every
  window, so latency measured from arrival grows with run length there and is
  not a steady-state quantity.
- The linear form of the model mirrors the simulator's own energy accounting
  (per-command energy plus time-proportional energy). The fits therefore show
  that the effects the model omits are small in this simulator; they do not
  validate the simulator against hardware.

## 11. License and citation

The patch in `simulator/fixedwork.patch` modifies DRAMsim3, which is
distributed under the MIT License; see `NOTICE`. DRAMsim3 itself is fetched at
build time and is not included.

The accompanying paper is under review; citation details will be added here.
