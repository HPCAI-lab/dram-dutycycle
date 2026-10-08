/* lackey_filter.c -- turns valgrind lackey --trace-mem output into
 * pinatrace-style records ("<ip>: R|W <addr> <size>"), one record per data
 * access, a modify (M) becoming a read followed by a write, as pinatrace
 * records an access that both reads and writes memory.
 *  - bridge file: the first BRIDGE_N records of the whole program
 *  - window files: records issued by instructions inside [lo,hi) (the gemm
 *    function), counted per call of gemm, from offset S for length L
 * Usage: lackey_filter lo hi entry bridge_n bridge_path summary_path
 *                      [call offset length path]...
 * With LACKEY_STOP_EARLY=1 the filter stops reading (ending the traced
 * program) as soon as every requested output is complete; the summary then
 * reports "complete": false and per-call counts cover only what was read.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXW 16
#define MAXC 16
int main(int argc, char **argv) {
    if (argc < 7 || (argc - 7) % 4) { fprintf(stderr, "bad args\n"); return 2; }
    uint64_t lo = strtoull(argv[1], 0, 16), hi = strtoull(argv[2], 0, 16), entry = strtoull(argv[3], 0, 16);
    uint64_t bridge_n = strtoull(argv[4], 0, 10);
    FILE *bf = bridge_n ? fopen(argv[5], "w") : NULL;
    FILE *sf = fopen(argv[6], "w");
    int nw = (argc - 7) / 4;
    int wcall[MAXW]; uint64_t woff[MAXW], wlen[MAXW], wdone[MAXW]; FILE *wf[MAXW];
    for (int i = 0; i < nw; i++) {
        wcall[i] = atoi(argv[7 + 4*i]); woff[i] = strtoull(argv[8 + 4*i], 0, 10);
        wlen[i] = strtoull(argv[9 + 4*i], 0, 10); wf[i] = fopen(argv[10 + 4*i], "w"); wdone[i] = 0;
    }
    char line[256]; uint64_t ip = 0, total = 0, instr = 0;
    int call = 0; uint64_t kcount[MAXC] = {0}, kreads[MAXC] = {0}, kwrites[MAXC] = {0};
    uint64_t greads = 0, gwrites = 0, nmod = 0;
    static char buf[1 << 20]; setvbuf(stdin, buf, _IOFBF, sizeof buf);
    const char *se = getenv("LACKEY_STOP_EARLY"); int stop_early = se && se[0] == '1'; int complete = 1;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == 'I') { ip = strtoull(line + 3, 0, 16); instr++; if (ip == entry && call < MAXC - 1) call++; continue; }
        if (line[0] != ' ' || (line[1] != 'L' && line[1] != 'S' && line[1] != 'M')) continue;
        char *comma; uint64_t addr = strtoull(line + 3, &comma, 16); int size = atoi(comma + 1);
        char ops[2]; int nops = 0;
        if (line[1] == 'L') ops[nops++] = 'R'; else if (line[1] == 'S') ops[nops++] = 'W'; else { ops[nops++] = 'R'; ops[nops++] = 'W'; nmod++; }
        int inside = (ip >= lo && ip < hi);
        for (int o = 0; o < nops; o++) {
            if (ops[o] == 'R') greads++; else gwrites++;
            if (bf && total < bridge_n) fprintf(bf, "0x%lx: %c 0x%lx %d\n", ip, ops[o], addr, size);
            total++;
            if (bf && total == bridge_n) { fclose(bf); bf = NULL; }
            if (inside && call > 0) {
                uint64_t k = kcount[call]++;
                if (ops[o] == 'R') kreads[call]++; else kwrites[call]++;
                for (int i = 0; i < nw; i++)
                    if (wcall[i] == call && k >= woff[i] && k < woff[i] + wlen[i]) {
                        fprintf(wf[i], "0x%lx: %c 0x%lx %d\n", ip, ops[o], addr, size);
                        if (++wdone[i] == wlen[i]) { fclose(wf[i]); wf[i] = NULL; }
                    }
            }
        }
        if (stop_early && bf == NULL) {
            int all = 1;
            for (int i = 0; i < nw; i++) if (wf[i]) all = 0;
            if (all) { complete = 0; break; }
        }
    }
    if (bf) fclose(bf);
    fprintf(sf, "{\n  \"complete\": %s,\n  \"instructions\": %lu,\n  \"records_total\": %lu,\n  \"reads_total\": %lu,\n  \"writes_total\": %lu,\n  \"modify_ops\": %lu,\n  \"gemm_calls\": %d,\n", complete ? "true" : "false", instr, total, greads, gwrites, nmod, call);
    fprintf(sf, "  \"per_call\": [");
    for (int c = 1; c <= call; c++) fprintf(sf, "%s{\"call\": %d, \"records\": %lu, \"reads\": %lu, \"writes\": %lu}", c > 1 ? ", " : "", c, kcount[c], kreads[c], kwrites[c]);
    fprintf(sf, "],\n  \"windows\": [");
    for (int i = 0; i < nw; i++) fprintf(sf, "%s{\"call\": %d, \"offset\": %lu, \"length\": %lu, \"written\": %lu}", i ? ", " : "", wcall[i], woff[i], wlen[i], wdone[i]);
    fprintf(sf, "]\n}\n"); fclose(sf);
    for (int i = 0; i < nw; i++) if (wf[i]) fclose(wf[i]);
    return 0;
}
