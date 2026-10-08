"""convert_poisson.py -- converts pinatrace-style records to DRAMsim3 trace
format with exponentially distributed gaps between requests (mean = stride
cycles), so the average offered load equals that of the uniform converter at
the same stride, but arrivals are random rather than evenly spaced.
Arrival cycles are the floor of the running sum of the gaps, so two requests
can share a cycle; the simulator's trace reader accepts at most one per
cycle. Seeded: the same window, stride and seed give the same trace.
Usage: python3 convert_poisson.py <in.pin[.gz]> <stride> <seed> > out.trace"""
import sys, random, gzip
src, stride, seed = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed)
t = 0.0
out = sys.stdout
with (gzip.open(src, 'rt') if src.endswith('.gz') else open(src)) as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.split()
        if len(p) < 3:
            continue
        op = 'READ' if p[1] in ('R', 'r') else 'WRITE'
        out.write(f"{p[2]} {op} {int(t)}\n")
        t += rng.expovariate(1.0 / stride)
