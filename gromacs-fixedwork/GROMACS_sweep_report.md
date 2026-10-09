# GROMACS fixed-work sweep: results and proposed edits to IPDPS_V3_3.tex

8 October 2026. No change has been made to `IPDPS_V3_3.tex`. Every edit below waits for approval.

## 1. What was done

### 1.1 Designed artifacts (correct by construction or inspection)

- **Workload input.** 31,720 SPC water molecules (95,160 atoms) from `gmx solvate` in a 9.865-nanometer cube, OPLS-AA water parameters.
  - Electrostatics: particle-mesh Ewald (PME), 1.0-nanometer cutoffs, Fourier spacing 0.125 nanometers (80 by 80 by 80 grid), Verlet buffer.
  - Two MD steps (`md.mdp`).
  - This system stands in for the matched-load study's GROMACS input, whose system and version were never recorded.
- **Trace filter.** `gmx_filter.c` follows `lackey_filter.c`: one record per data access, and a modify becomes a read followed by a write. It also:
  - attributes every record to the function whose instructions issued it;
  - marks MD steps by entries to `gmx::do_force`;
  - stores every record of step 1 in binary;
  - writes a per-function, per-step record profile.
- **Window cutter** (`cut_windows.py`) and **dumper** (`dump_range.c`).
- **Stock-equivalence check** (`stock_check.py`), **cross-workload report** (`gmx_report.py`), **mechanism analyses** (`mechanism.py`) and **independent recomputation** (`verify.py`).

### 1.2 Asserted methodology (how the experiment was run)

- **Build.** GROMACS 2025.3 (GitHub tag `v2025.3`), built with gcc 13.3.0 at `-O1 -fno-tree-vectorize`.
  - Scalar and single-threaded: `GMX_SIMD=None`, no OpenMP or MPI, fftpack, static.
  - GROMACS adds its own `-funroll-all-loops -fexcess-precision=fast`.
  - The nonbonded interactions therefore run in the plain-C 4x4 reference kernel, so every access reaches the trace, as with GEMM-AI.
- **Tracing.** valgrind 3.22.0 Lackey `--trace-mem=yes`, the same as GEMM-AI. It took 4 hours 25 minutes for setup and steps 0 and 1.
- **Target step.** Step 1, a steady-state step with neither pair search nor energy evaluation. It holds 1,211,292,238 records.
  - Pair search runs on step 0, in `nbnxn_make_pairlist_part` and its helpers, about 0.8 billion records.
  - With `nstlist` raised to 50 by mdrun, that averages about 1.3% of a step's records, so search is not windowed.
- **Window rule.** Each function that issues at least 1% of the step's records gets one 100,000-record window centered in its records.

  | Function | Share of step-1 records |
  |---|---|
  | `nbnxn_kernel_4x4_ElecQSTab_VdwLJ_F_ref` (nonbonded) | 88.28% |
  | `spread_on_grid` (PME spreading) | 2.64% |
  | `fftpack_cfftf1` (3D-FFT complex pass) | 2.49% |
  | `gather_f_bsplines` (PME gathering) | 1.85% |

  - Together they cover 95.26% of the step.
  - The kernel carries 88% of the step, so it gets windows at one quarter, one half and three quarters of its records.
  - Windows: `g1_nb25`, `g2_nb50`, `g3_nb75`, `g4_spread`, `g5_fft`, `g6_gather`.
  - A 1,000,000-record window centered in the kernel (`g2_nb50_1M`) checks window length.
- **Sweep.** Identical to GEMM-AI: same eleven strides, same fit, held-out and extrapolated strides, two arrival processes.
  - New Poisson seed bases, 11000 to 16000.
  - Same patched DRAMsim3 and the same analysis scripts, with only the window list changed.
  - 143 runs in total.
- **Bridge.** 40,000,050 contiguous records centered in step 1, at stride 50, under the 2-billion-cycle budget with the unmodified stop rule.
  - The first 40 million records of the program are setup, not MD, so they cannot serve as the bridge.
  - The step-centered segment lies entirely in the nonbonded kernel.
- **Simulator check.** The rebuilt simulator passes all nine patch tests and reproduces eight stored GEMM-AI runs counter for counter.

### 1.3 Measured results (simulator output)

All values below are recomputed independently from the raw `dramsim3.json` files, and they match the pipeline.

## 2. Where the results land

Evenly spaced arrivals and the relative-error fit throughout, unless stated otherwise.

| Window | Hit rate, stride 50 | a | b | Held-out MAPE | Max held-out | Error, stride 200 | Error, stride 5 |
|---|---|---|---|---|---|---|---|
| g1_nb25 | 91.43% | 3,891.0 | 1,165.5 | 0.42% | 0.55% | 0.87% | 2.18% |
| g2_nb50 | 88.68% | 4,004.4 | 1,164.4 | 0.70% | 0.89% | 1.54% | 3.14% |
| g3_nb75 | 82.76% | 4,083.8 | 1,169.2 | 0.97% | 1.10% | 1.76% | 3.97% |
| g4_spread | 91.91% | 4,144.0 | 1,157.1 | 0.88% | 1.01% | 3.01% | 3.04% |
| g5_fft | 96.85% | 3,572.5 | 1,152.2 | 2.25% | 2.77% | 3.70% | 3.42% |
| g6_gather | 80.09% | 4,256.2 | 1,174.8 | 1.08% | 1.81% | 3.35% | 5.03% |
| g2_nb50_1M | 87.40% | 4,001.7 | 1,165.8 | 0.90% | 1.12% | 1.77% | 3.44% |

### 2.1 Claim C1, the linear model: the form holds, the accuracy degrades

- **Fit quality.** Coefficient of determination 0.99922 to 0.99998.
- **Held-out MAPE is 0.42% to 2.25%, against 0.03% to 0.13% for GEMM-AI.** That is 3 to 17 times larger.
- **Extrapolation errors reach 3.70% at stride 200 and 5.03% at stride 5.** The abstract's "errors below 2.7%" does not hold for GROMACS.
- **Stride 5 is a real test here.** The pipeline flags stride 5 as saturated in every window, but that is an artifact of the definition: each window's maximum throughput is at stride 5 itself. Five windows still complete one request every 5.00 cycles there; only `g5_fft` saturates, at 5.70. So for GROMACS, stride 5 does test extrapolation, unlike GEMM-AI.
- **Ordinary least squares is better here.** Its held-out MAPE is 0.20% to 1.78%, which reverses the GEMM-AI rationale for the relative-error fit.
- **Robustness.**
  - Exponential arrivals change a by at most 1.51% and b by at most 0.08%, and errors by at most 0.52 points. For GEMM-AI the coefficient bound was 0.46%.
  - The 1,000,000-record window gives coefficients within 0.12% of `g2_nb50`.

### 2.2 Why the accuracy degrades (mechanism confirmed in the source)

- **In the DRAMsim3 source.** `Controller::AddTransaction` completes three kinds of request without issuing a command:
  - a read whose address is queued for writing (served from the write buffer, latency 1 cycle);
  - a read to an address that already has a read pending (merged);
  - a write to an address that already has a write pending (merged).
- **In the GROMACS runs.** Commands per request fall with load, from 0.917 at stride 100 to 0.82 to 0.87 at stride 6 in the kernel windows, and to 0.79 in `g5_fft`. In GEMM-AI they are exactly 1.000 at every stride. GROMACS accumulates forces by read-modify-write, so its trace is full of repeated addresses.
- **Effect on the coefficients.** Command energy per request therefore rises with cycles per request.
  - That adds 10.6 to 27.8 units per cycle to the fitted slope (GEMM-AI: 0.6 to 2.8).
  - It bends energy per request away from a straight line.
- **Effect on the slope bound.** The time-proportional slope alone, `b_time` = 1,127.1 to 1,154.9, lies inside the 881.0 to 1,169.0 bound in every window. The total slopes of `g3_nb75` (1,169.2) and `g6_gather` (1,174.8) lie above it.
- **Effect on duty cycle.** The paper's identity "duty cycle equals 4/c" does not hold for GROMACS.
  - At stride 50, bus occupancy (4 times commands over cycles) is 0.073 to 0.079, against 4/c = 0.080.
  - At stride 5 it is 0.53 to 0.70, against 0.80.
- **Which reads skip a command.** At low load, every read without a command is a write-buffer hit (latency 1). At high load, merged duplicate reads add to them.

### 2.3 Simulator integrity: the deadlock fix is exercised (protocol flag)

- **How often.** Change 4 of the patch (schedule a read when the oldest write is blocked by an older read to the same address) fired in 114 of 143 runs.
  - From 0 to 7,841 times per run.
  - At most 1,053 times per run at stride 8 and above.
  - `analyze_sweep.py` counts every such run as an integrity failure, by the protocol's own rule.
- **What it changes.** I ran stock DRAMsim3 on every evenly spaced run, capped at the patched run's duration.
  - **Strides 8 and above (up to 50% of peak), all windows:** stock completes the same work, and total energy differs by at most 0.013%. Across every run that stock completes, the largest difference is 0.034%. The stalls clear on their own, so the fix does not alter these results.
  - **Strides 7 to 5:** stock deadlocks in eight runs. `g3_nb75` at strides 6 and 5 completes 37.5% of its requests; `g4_spread` at stride 5, 5.1%; `g5_fft` at strides 6 and 5, 98.5% and 21.5%; `g2_nb50_1M` at strides 7, 6 and 5, 44.2%, 90.9% and 15.0%.
  - These points exist only because of the patch, and stride 6 is a fitted point.
- **Sensitivity refit** on strides 100, 50, 25 and 12 only, where stock and patched agree, scored on held-out strides 10, 9 and 8:
  - held-out MAPE 0.27% to 2.19% (kernel windows 0.27% to 0.58%);
  - stride-200 error 0.58% to 2.33%.

### 2.4 Claim C2, the race-to-idle criterion: transfers unchanged

Stride 8 (50% of peak) against stride 50 (8%), over the slower run's duration:

| Idle state | P_idle | Energy change, GROMACS (6 windows) | GEMM-AI, for comparison |
|---|---|---|---|
| None (fixed work) | 0 | -78.1% to -78.6% | -73.8% to -76.7% |
| Self-refresh | 576.0 | -38.8% to -39.7% | -37.0% to -38.0% |
| Precharge power-down | 650.6 | -33.7% to -34.7% | -32.2% to -32.9% |
| Precharge standby | 881.0 | -17.8% to -19.3% | -17.4% to -17.5% |
| Active standby | 1,169.0 | -0.0% to +2.0% | +1.0% to +2.0% |

- **The outcome is still set by the idle state.** Accounting without idle energy overstates the saving by a factor of 2.0 to 4.4, the same factor as for GEMM-AI.
- **Break-even power, stride 50 against stride 8:** 1,140.1 to 1,169.4. Over all pairs among strides 200, 50, 8 and 6 it is 1,097.8 to 1,410.7 (GEMM-AI: 1,092.5 to 1,162.2).
- **Active standby lies at break-even for `g6_gather`.** That window's b exceeds 1,169.0, so it shows the active-standby sign sensitivity the Limitations paragraph already describes.

### 2.5 Claim C3, per-byte error: transfers unchanged

- With no idle energy, the per-byte error runs from -72.8% to -73.5% at stride 200 and from +464.0% to +492.5% at stride 6.
- With precharge standby it runs from -14.4% to -17.3% at stride 200 and from +6.6% to +7.1% at stride 6.

### 2.6 Claim C4, composition across access patterns: extends across workloads

Standby-and-refresh share, evenly spaced arrivals:

| Stride (percent of peak) | GROMACS, 6 windows | GEMM-AI 85% pattern | GEMM-AI 48% pattern | Spread, GROMACS and GEMM-AI 85% | Spread, all 12 windows |
|---|---|---|---|---|---|
| 200 (2%) | 97.52 to 97.97 | 97.50 | 96.60 | 0.47 | 1.37 |
| 50 (8%) | 91.47 to 92.96 | 91.23 to 91.25 | 87.89 | 1.73 | 5.07 |
| 12 (33%) | 73.85 to 77.18 | 72.17 to 72.20 | 63.63 | 5.02 | 13.55 |
| 6 (67%) | 62.26 to 67.21 | 56.50 to 56.55 | 50.40 | 10.71 | 16.81 |

- At 8% of peak, the GROMACS windows and the GEMM-AI 85% pattern agree within 1.73 points. That is consistent with Section V's "under two percentage points".
- The spread grows with load, so the invariance at low load and the divergence at high load hold across two workloads.

### 2.7 Energy-delay basin: does not transfer (null result for GROMACS)

- With delay as mean read latency from acceptance, over the nine fitted and held-out strides, EDP is lowest at:
  - stride 6 for the three kernel windows and `g6_gather`, the highest-load point evaluated, so there is no interior minimum;
  - stride 7 for `g4_spread`;
  - stride 12 for `g5_fft`.
- With completion time as delay, EDP falls monotonically, as for GEMM-AI.

### 2.8 Bridge against the GROMACS rows of Tables III and IV

| Quantity | Tables III and IV GROMACS | Step-1 segment (this trace) |
|---|---|---|
| Read row-buffer hit rate | 86.36% | 87.98% |
| Write share of requests | 35.28% | 22.56% |
| Read requests per read command | 1.162 | 1.061 |
| Active standby | 74.79% | 76.06% |
| Precharge standby | 4.53% | 2.93% |
| Refresh | 13.87% | 13.70% |
| Read data transfer | 4.03% | 5.22% |
| Write data transfer | 1.87% | 1.12% |
| Activate | 0.92% | 0.97% |
| Standby and refresh | 93.19% | 92.69% |

- The composition is close: the combined standby-and-refresh share is within 0.50 points, and the largest single gap is 1.60 points (precharge standby).
- It is not a reproduction, the way the GEMM-AI bridge was (0.10 and 0.03 points).
- The paper should present this as a new GROMACS trace, not as the Table III trace.
- Change 4 fired 310 times in the bridge run.

## 3. Proposed edits (for approval; none applied)

**Page budget first.** Page 9 is full in both columns, so these edits add about 0.4 to 0.5 page. If references are excluded from the 10-page limit (confirm against the IPDPS 2027 call), the body can run onto page 10. Otherwise I'll need cuts to match, and I'll propose them.

### E1. Section VI-A: new paragraph after line 629 (the GEMM-AI window description)

> To test the protocol on a second, less regular workload, we traced one steady-state molecular-dynamics step of GROMACS~\cite{abraham2015gromacs} (version 2025.3, built for scalar execution with the same compiler and flags, so that its nonbonded interactions run in the plain-C reference kernel) on a box of 95{,}160 atoms of water with particle-mesh Ewald electrostatics. The step issues 1.21 billion records, 88.3\% of them from the nonbonded kernel and 7.0\% from charge spreading, the complex Fourier-transform pass, and force gathering of the particle-mesh Ewald method. One window is centered in each of these functions, and two more at one quarter and three quarters of the kernel's records, giving six windows with read row-buffer hit rates of 80\% to 97\% at stride 50. This trace is not the one behind Tables~\ref{tab:char} and~\ref{tab:energy}: a 40-million-record segment of the step, simulated under the matched-load protocol, gives a read row-buffer hit rate of 87.98\% against 86.36\% and a standby-and-refresh share of 92.69\% against 93.19\%.

New reference, to verify before use. The authors, title, journal, year and DOI come from mdrun's own citation notice; the volume and pages are my recollection:

```latex
\bibitem{abraham2015gromacs}
M.~J. Abraham, T.~Murtola, R.~Schulz, S.~P\'all, J.~C. Smith, B.~Hess, and
E.~Lindahl, ``GROMACS: High performance molecular simulations through
multi-level parallelism from laptops to supercomputers,'' \emph{SoftwareX},
vol.~1--2, pp.~19--25, 2015, doi:~10.1016/j.softx.2015.06.001.
```

### E2. Section VI-A, line 635 to 636: the sentence becomes false

- Current: "We also fixed a scheduling deadlock, which none of the reported runs encountered."
- Proposed: "We also fixed a scheduling deadlock, which no GEMM-AI run encountered. In the GROMACS runs the fix is used, but up to 50\% of channel peak the unmodified simulator completes the same work with energies within 0.013\%; above that, it deadlocks on some windows (Section~\ref{sec:fw-validation})."

### E3. Section VI-C: new paragraph after line 706 (GROMACS validation)

> For GROMACS the same protocol gives held-out MAPEs of 0.42\% to 2.25\% and extrapolation errors of up to 3.70\% at stride 200 and 5.03\% at stride 5, which here completes one request every 5.00 cycles in five of the six windows and so tests extrapolation. The larger errors have an identifiable cause. The simulator's controller completes a read whose address is queued for writing, and merges a request whose address is already pending, without issuing a command, and the force accumulation of molecular dynamics makes such repeats common. The commands issued per request therefore fall with load, from 0.92 at stride 100 to between 0.79 and 0.87 at stride 6, against exactly 1.00 for GEMM-AI. Command energy per request then varies with $c$, contributes 11 to 28 units per cycle to the fitted slope, and bends the line. The time-proportional part of the slope alone, 1{,}127.1 to 1{,}154.9 units per cycle, stays inside the bound of Section~\ref{sec:fw-model}, although the total slopes of two windows, 1{,}169.2 and 1{,}174.8, exceed it. For the same reason, $4/c$ overstates bus occupancy for GROMACS, by 1\% to 9\% at stride 50. Fitted only on strides 100 to 12, where the unmodified simulator agrees with the modified one, the model predicts strides 10, 9, and 8 with MAPEs of 0.27\% to 2.19\%.

### E4. Table V: add a GROMACS column

- The column spec goes from `lrr` to `lrrr`, the "common to both patterns" header becomes "common to all", and the caption gains "and GROMACS".
- Rows: windows 6; hit rate 80.1 to 96.9%; a 3,572.5 to 4,256.2; b 1,152.2 to 1,174.8; R squared 0.99922 to 0.99998; leave-one-out MAPE 0.59 to 3.71%; held-out MAPE 0.42 to 2.25%; maximum held-out error 0.55 to 2.77%; stride-200 error 0.87 to 3.70%; stride-5 error 2.18 to 5.03%.
- More rows: cycles per request at strides 7, 6 and 5 of 7.0, 6.0 and 5.0 to 5.7; energy per request 116,506 to 120,036 at stride 100 and 10,294 to 11,137 at stride 6.
- More rows: read-transfer share 2.59 to 3.25% at stride 100 and 20.95 to 31.05% at stride 6; standby-and-refresh share 95.38 to 96.18% at stride 100 and 62.26 to 67.21% at stride 6; EDP-minimizing stride 6, 7 or 12.
- **Risk:** a third column may not fit in a single-column table at `\footnotesize`. The fallback is a separate small GROMACS table of six rows.

### E5. Table VI and Section VI-F: add GROMACS

- Add an "Energy change, GROMACS" column with the ranges from Section 2.4 above.
- Add one sentence in Section VI-F: "For the six GROMACS windows the faster operating point saves 38.8\% to 39.7\%, 33.7\% to 34.7\%, and 17.8\% to 19.3\% in the three low-power states and changes energy by $-0.0$\% to $+2.0$\% in active standby, so the outcome is again set by the idle state."
- The Table VI column spec `lrrrrr` grows by one column. The same width risk applies, and a single "GROMACS range" row under the table is the fallback.

### E6. Section VI-E: one sentence after line 826

> Across workloads the pattern repeats: the standby-and-refresh shares of the six GROMACS windows and the five windows of the 85\% pattern differ by 0.47 percentage points at 2\% of peak, 1.73 at 8\%, and 10.71 at 67\%.

### E7. Section VI-D: one sentence after line 801

> For GROMACS the latency-based product has no interior minimum in four of six windows; it falls up to stride 6, the highest load evaluated.

### E8. Consistency edits that the new results force

| Where | Current | Proposed |
|---|---|---|
| Abstract | "on six windows of a GEMM-AI trace forming two access patterns" | "on six windows each of a GEMM-AI trace and a GROMACS molecular-dynamics trace" |
| Abstract | "0.03\% to 0.13\% ... errors below 2.7\%" | "0.03\% to 0.13\% for GEMM-AI and 0.42\% to 2.25\% for GROMACS, and points outside the fitted range with errors below 2.7\% and 5.1\%" |
| Section I, line 145 | "it equals $4/c$" | "it equals $4/c$ when every request issues a command, as for GEMM-AI" |
| Contribution 1, line 168 | "Fitted on five operating points in each of six GEMM-AI trace windows" | Add GROMACS with its MAPE range |
| Section I, line 197 | "The sweep covers one DDR4 configuration and one workload." | "The sweep covers one DDR4 configuration and two workloads." |
| Section III-C, line 453 | "For GEMM-AI it varies" | "For GEMM-AI and GROMACS it varies" |
| Section VII-B, line 1007 | "with at most 0.13\% MAPE" | "with at most 0.13\% MAPE for GEMM-AI and 2.25\% for GROMACS" |
| Section VIII, scope, lines 1047 to 1059 | "covers six windows of one workload ... carry over to STREAM, GROMACS, and BLAST" | Two workloads; the linear form held for GROMACS with lower accuracy, owing to request merging; STREAM, BLAST, DDR5 and HBM remain |
| Section VIII, trace-driven load, line 1108 | "the BLAST and GEMM-AI matched-load runs were bounded" | "the GEMM-AI matched-load run was bounded" (separate BLAST fix; the data confirm BLAST stopped at 1,963,429,800 cycles) |
| Conclusions, line 1129 | "on six trace windows forming two access patterns" | Add GROMACS |
| Section VI-A, line 609 | "Two of the four matched-load runs" | "One of the four matched-load runs" (separate BLAST fix, same evidence) |

## 4. Decisions needed

1. **Include GROMACS in the paper at all?** It broadens the evidence to a second, irregular workload. It also widens the headline accuracy numbers by about 10 times and adds a simulator-behavior caveat. My recommendation is to include it, because a reviewer will ask about generality, and an identified mechanism is a stronger answer than "future work".
2. **Which fit is the headline for GROMACS?** Either the identical protocol (strides 100 to 6, recommended for comparability), or the stock-equivalent refit (strides 100 to 12).
3. **Page budget.** Can the body run onto page 10? If not, which cuts?
4. **Placement.** A third column in Tables V and VI, or a separate small GROMACS table?
