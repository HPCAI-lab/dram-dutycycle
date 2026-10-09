"""read_run.py -- reads one run directory and returns every derived quantity,
with each energy total checked against the simulator's own total."""
import json, os
# DUTYCYCLE_UNIT=command normalizes energy and cycles per read or write command
# instead of per completed request (they coincide when every request issues a command).
UNIT = os.environ.get('DUTYCYCLE_UNIT', 'request')
TCK_NS = 0.63
ACT_E, RD_E, WR_E = 6720.0, 4454.4, 3763.2   # per-command increments, from the .ini (simulator units)
def load(d):
    j = json.load(open(f"{d}/dramsim3.json"))
    assert list(j.keys()) == ['0'], "expected exactly one channel"
    c = j['0']; cpu = json.load(open(f"{d}/cpu_stats.json")); meta = json.load(open(f"{d}/meta.json"))
    vs = lambda k: sum(float(v) for v in c[k].values())
    comp = {"act_stb": vs('act_stb_energy'), "pre_stb": vs('pre_stb_energy'), "sref": vs('sref_energy'),
            "refresh": float(c['ref_energy']) + float(c['refb_energy']), "read_xfer": float(c['read_energy']),
            "write_xfer": float(c['write_energy']), "activate": float(c['act_energy'])}
    tot = sum(comp.values())
    assert abs(tot - float(c['total_energy'])) <= 1e-9 * tot, "component sum differs from total_energy"
    n = int(c['num_reads_done']) + int(c['num_writes_done'])
    cyc = int(c['num_cycles']); ranks = len(c['rank_active_cycles'])
    lat = {int(k): int(v) for k, v in c['read_latency'].items()}
    nl = sum(lat.values()); mean_acc = sum(k * v for k, v in lat.items()) / nl if nl else 0.0
    def pct(q):
        tgt = q * (nl - 1); run = 0
        for k in sorted(lat):
            run += lat[k]
            if run - 1 >= tgt: return k
    r = dict(meta)
    r.update({"completed": n, "cycles": cyc, "c_per_req": cyc / n, "E_total": tot, "E_per_req": tot / n,
              "E_per_req_pJ": tot / n * TCK_NS, **{f"share_{k}": v / tot for k, v in comp.items()},
              "share_residency": (comp['act_stb'] + comp['pre_stb'] + comp['sref'] + comp['refresh']) / tot,
              "reads_done": int(c['num_reads_done']), "writes_done": int(c['num_writes_done']),
              "read_cmds": int(c['num_read_cmds']), "write_cmds": int(c['num_write_cmds']),
              "act_cmds": int(c['num_act_cmds']), "ref_cmds": int(c['num_ref_cmds']),
              "read_row_hit": int(c['num_read_row_hits']) / max(1, int(c['num_read_cmds'])),
              "write_row_hit": int(c['num_write_row_hits']) / max(1, int(c['num_write_cmds'])),
              "act_per_req": int(c['num_act_cmds']) / n,
              "req_per_read_cmd": int(c['num_reads_done']) / max(1, int(c['num_read_cmds'])),
              "rank_active_frac": vs('rank_active_cycles') / (ranks * cyc),
              "lat_acc_mean": mean_acc, "lat_acc_p95": pct(0.95), "lat_acc_p99": pct(0.99),
              "lat_arr_mean": cpu['read_lat_from_arrival_mean'], "lat_arr_p99": cpu['read_lat_from_arrival_p99'],
              "bw_reported": float(c['average_bandwidth']), "bw_nominal": 64 / (meta['stride'] * TCK_NS),
              "rw_fallbacks": cpu['rw_fallbacks'], "stall_cycles": cpu['stall_cycles'],
              "cpu_ok": cpu['trace_exhausted'] and cpu['memory_idle'] and cpu['added_reads'] == cpu['done_reads']
                        and cpu['added_writes'] == cpu['done_writes'] and n == meta['records']})
    if UNIT == 'command':
        k = r['read_cmds'] + r['write_cmds']
        r.update({'c_per_req': cyc / k, 'E_per_req': tot / k, 'E_per_req_pJ': tot / k * TCK_NS, 'units': k})
    else:
        r['units'] = n
    return r
