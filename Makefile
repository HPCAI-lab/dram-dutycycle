# Reproduction targets; README.md, Section 9, describes them.
#   OUT          where a sweep writes its runs (OUT/runs); default build
#   RESULTS      directory holding runs/ that analysis and figure read; default results
#   WINDOWS_DIR  trace windows for the sweep; default data/windows
#   TRACE_DIR    where make trace writes, and check-trace and bridge read; default build/traces
PYTHON      ?= python3
OUT         ?= build
RESULTS     ?= results
WINDOWS_DIR ?= data/windows
TRACE_DIR   ?= build/traces
JOBS        ?= $(shell nproc)
SIM         := build/bin/dramsim3main

.PHONY: all simulator test trace check-trace bridge sweep compare analysis figure clean

all: simulator test sweep compare

simulator: $(SIM)

$(SIM): simulator/fixedwork.patch simulator/DRAMSIM3_COMMIT simulator/build_dramsim3.sh
	bash simulator/build_dramsim3.sh

test: $(SIM)
	bash tests/test_patch.sh

trace:
	TRACE_OUT=$(TRACE_DIR) bash tracing/trace_gemm.sh

check-trace:
	@status=0; for f in $(TRACE_DIR)/w*.pin; do \
	    $(PYTHON) tracing/compare_traces.py $$f data/windows/$$(basename $$f).gz || status=1; \
	done; exit $$status

bridge: $(SIM)
	bash sweep/run_bridge.sh $(TRACE_DIR)/bridge.pin build/bridge

sweep: $(SIM)
	$(PYTHON) sweep/run_sweep.py --windows-dir $(WINDOWS_DIR) --out $(OUT)/runs --jobs $(JOBS)
	$(PYTHON) sweep/run_sweep.py --windows-dir $(WINDOWS_DIR) --out $(OUT)/runs --jobs $(JOBS) \
	    --windows w5_ff1_1M --arrivals uniform

compare:
	$(PYTHON) analysis/compare_runs.py $(OUT)/runs results/runs

analysis:
	$(PYTHON) analysis/analyze_sweep.py $(RESULTS)

figure:
	$(PYTHON) analysis/plot_sweep.py $(RESULTS)

clean:
	rm -rf build
