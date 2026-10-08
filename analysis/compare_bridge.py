"""Bridge check: the regenerated GEMM-AI trace at stride 50 under the 2-billion-
cycle budget, against the GEMM-AI rows of Tables III and IV of the paper.
Usage: python3 analysis/compare_bridge.py <bridge_run_dir>"""
import json, sys
d = sys.argv[1]
c = json.load(open(f"{d}/dramsim3.json"))["0"]
vs = lambda k: sum(float(v) for v in c[k].values())
comp = {"Activate standby": vs("act_stb_energy"), "Precharge standby": vs("pre_stb_energy"),
        "Refresh": float(c["ref_energy"]) + float(c["refb_energy"]), "Read data transfer": float(c["read_energy"]),
        "Write data transfer": float(c["write_energy"]), "Activate": float(c["act_energy"])}
tot = sum(comp.values()) + vs("sref_energy")
assert abs(tot - float(c["total_energy"])) < 1e-9 * tot
paper = {"Activate standby": 70.99, "Precharge standby": 6.59, "Refresh": 13.71, "Read data transfer": 6.91,
         "Write data transfer": 0.22, "Activate": 1.58}
print(f"{'quantity':28s}{'paper':>14s}{'this trace':>14s}")
for n, p, v in [("read requests done", 38573371, int(c["num_reads_done"])), ("write requests done", 1426629, int(c["num_writes_done"])),
                ("read commands", 38571501, int(c["num_read_cmds"])), ("cycles", 2000000000, int(c["num_cycles"]))]:
    print(f"{n:28s}{p:>14,d}{v:>14,d}")
print(f"{'read row-buffer hit rate %':28s}{84.86:>14.2f}{100 * int(c['num_read_row_hits']) / int(c['num_read_cmds']):>14.2f}")
print(f"{'total energy (1e12 units)':28s}{2.49:>14.2f}{tot / 1e12:>14.2f}")
print("energy shares, percent of total:")
for k in paper:
    print(f"  {k:26s}{paper[k]:>14.2f}{100 * comp[k] / tot:>14.2f}")
