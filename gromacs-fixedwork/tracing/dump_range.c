/* dump_range.c -- records [start, start+count) of a gmx_filter step file as pinatrace-style text */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef struct __attribute__((packed)) { uint64_t addr; uint32_t ipoff; uint16_t fid; uint8_t op, size; } rec_t;
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"); uint64_t base = strtoull(argv[2], 0, 16);
    uint64_t start = strtoull(argv[3], 0, 10), count = strtoull(argv[4], 0, 10);
    fseeko(f, (off_t)(start * sizeof(rec_t)), SEEK_SET); rec_t r; static char ob[1 << 22]; setvbuf(stdout, ob, _IOFBF, sizeof ob);
    for (uint64_t i = 0; i < count && fread(&r, sizeof r, 1, f) == 1; i++)
        printf("0x%lx: %c 0x%lx %d\n", r.ipoff == 0xffffffffu ? 0 : base + r.ipoff, r.op, r.addr, r.size);
    return 0;
}
