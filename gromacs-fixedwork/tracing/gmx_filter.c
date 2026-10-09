/* gmx_filter.c -- valgrind lackey --trace-mem output of GROMACS mdrun to
 *   (1) a per-function, per-step record profile (JSON),
 *   (2) every record of one target MD step, in compact binary, for cutting
 *       windows centred in each dominant function afterwards,
 *   (3) the first BRIDGE_N records of the whole program as pinatrace-style text
 *       ("<ip>: R|W <addr> <size>"), for the bridge check.
 * Records follow lackey_filter.c: one per data access, a modify (M) becoming a
 * read followed by a write.
 * A step begins at each entry of the step-marker function (do_force); step 0 is
 * everything from the first entry to the second, and so on. Records before the
 * first entry are step -1 (setup). Tracing stops once the target step has ended.
 *
 * Usage: gmx_filter symbols.txt marker_entry_hex target_step bridge_n bridge.pin step.bin profile.json
 *   symbols.txt: lines "start_hex size_dec name", absolute addresses, sorted.
 * step.bin record (16 bytes): u64 addr; u32 ip - base (0xffffffff if outside
 *   the executable); u16 function index (0xffff if outside); u8 'R'/'W'; u8 size.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXSTEP 8
typedef struct { uint64_t lo, hi; char *name; } sym_t;
typedef struct __attribute__((packed)) { uint64_t addr; uint32_t ipoff; uint16_t fid; uint8_t op, size; } rec_t;
static sym_t *S; static int ns;
static int lookup(uint64_t ip) {
    int lo = 0, hi = ns - 1;
    while (lo <= hi) { int m = (lo + hi) / 2; if (ip < S[m].lo) hi = m - 1; else if (ip >= S[m].hi) lo = m + 1; else return m; }
    return -1;
}
int main(int argc, char **argv) {
    if (argc != 8) { fprintf(stderr, "bad args\n"); return 2; }
    FILE *sf = fopen(argv[1], "r"); char nm[4096]; uint64_t a, sz; int cap = 1 << 16;
    S = malloc(cap * sizeof *S);
    while (fscanf(sf, "%lx %lu %4095s", &a, &sz, nm) == 3) {
        if (ns == cap) { cap *= 2; S = realloc(S, cap * sizeof *S); }
        S[ns].lo = a; S[ns].hi = a + (sz ? sz : 1); S[ns].name = strdup(nm); ns++;
    }
    fclose(sf);
    if (ns >= 0xffff) { fprintf(stderr, "too many symbols\n"); return 2; }
    uint64_t base = S[0].lo, top = S[ns - 1].hi;
    uint64_t marker = strtoull(argv[2], 0, 16); int target = atoi(argv[3]);
    uint64_t bridge_n = strtoull(argv[4], 0, 10);
    FILE *bf = bridge_n ? fopen(argv[5], "w") : NULL;
    FILE *out = fopen(argv[6], "wb");
    static char obuf[1 << 22]; setvbuf(out, obuf, _IOFBF, sizeof obuf);
    uint64_t (*cnt)[MAXSTEP + 1] = calloc(ns + 1, sizeof *cnt);   /* index ns = outside */
    uint64_t steprec[MAXSTEP + 1] = {0}, total = 0, instr = 0;
    int step = -1; uint64_t ip = 0; int fid = -1;
    char line[256]; static char ibuf[1 << 22]; setvbuf(stdin, ibuf, _IOFBF, sizeof ibuf);
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == 'I') {
            ip = strtoull(line + 3, 0, 16); instr++;
            if (ip == marker) { step++; if (step > target) break; }
            if (fid < 0 || ip < S[fid].lo || ip >= S[fid].hi) fid = (ip >= base && ip < top) ? lookup(ip) : -1;
            continue;
        }
        if (line[0] != ' ' || (line[1] != 'L' && line[1] != 'S' && line[1] != 'M')) continue;
        char *comma; uint64_t addr = strtoull(line + 3, &comma, 16); int size = atoi(comma + 1);
        char ops[2]; int nops = 0;
        if (line[1] == 'L') ops[nops++] = 'R'; else if (line[1] == 'S') ops[nops++] = 'W'; else { ops[nops++] = 'R'; ops[nops++] = 'W'; }
        int col = step + 1; if (col > MAXSTEP) col = MAXSTEP;
        for (int o = 0; o < nops; o++) {
            if (bf && total < bridge_n) fprintf(bf, "0x%lx: %c 0x%lx %d\n", ip, ops[o], addr, size);
            total++;
            if (bf && total == bridge_n) { fclose(bf); bf = NULL; }
            cnt[fid >= 0 ? fid : ns][col]++; steprec[col]++;
            if (step == target) {
                rec_t r = { addr, (ip >= base && ip < base + 0xffffffffULL) ? (uint32_t)(ip - base) : 0xffffffffu,
                            (uint16_t)(fid >= 0 ? fid : 0xffff), (uint8_t)ops[o], (uint8_t)size };
                fwrite(&r, sizeof r, 1, out);
            }
        }
    }
    if (bf) fclose(bf);
    fclose(out);
    FILE *pf = fopen(argv[7], "w");
    fprintf(pf, "{\n \"base\": \"0x%lx\", \"marker\": \"0x%lx\", \"target_step\": %d, \"steps_seen\": %d,\n \"instructions\": %lu, \"records_total\": %lu,\n \"records_per_step\": {", base, marker, target, step, instr, total);
    for (int c = 0; c <= MAXSTEP; c++) fprintf(pf, "%s\"%d\": %lu", c ? ", " : "", c - 1, steprec[c]);
    fprintf(pf, "},\n \"functions\": [\n");
    int first = 1;
    for (int i = 0; i <= ns; i++) {
        uint64_t s = 0; for (int c = 0; c <= MAXSTEP; c++) s += cnt[i][c];
        if (!s) continue;
        fprintf(pf, "%s  {\"fid\": %d, \"name\": \"%s\", \"per_step\": [", first ? "" : ",\n", i < ns ? i : 0xffff, i < ns ? S[i].name : "<outside executable>");
        for (int c = 0; c <= MAXSTEP; c++) fprintf(pf, "%s%lu", c ? ", " : "", cnt[i][c]);
        fprintf(pf, "]}"); first = 0;
    }
    fprintf(pf, "\n ]\n}\n"); fclose(pf);
    return 0;
}
