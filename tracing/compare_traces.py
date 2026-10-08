"""Compares a regenerated trace file with a stored one, record by record.

Records must agree in instruction address, access type and size. Data addresses
must agree too, except stack addresses, which shift with the size of the
environment the program was started in. Exit status 0 when every record agrees
in that sense and the two files have the same number of records.
Usage: python3 tracing/compare_traces.py <new.pin[.gz]> <stored.pin[.gz]>"""
import gzip, sys

STACK = 0x1000000000          # client stack under valgrind on amd64 lies above this

def records(path):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        for line in f:
            p = line.split()
            if len(p) == 4 and not line.startswith("#"):
                yield p[0], p[1], int(p[2], 16), p[3]

new, old = list(records(sys.argv[1])), list(records(sys.argv[2]))
same = stack_only = other = 0
for a, b in zip(new, old):
    if a == b:
        same += 1
    elif a[0] == b[0] and a[1] == b[1] and a[3] == b[3] and a[2] >= STACK and b[2] >= STACK:
        stack_only += 1
    else:
        other += 1
ok = len(new) == len(old) and other == 0
print(f"{sys.argv[1]}: {len(new)} records against {len(old)}; identical {same}, "
      f"stack address only {stack_only}, other differences {other}: {'agree' if ok else 'DIFFER'}")
sys.exit(0 if ok else 1)
