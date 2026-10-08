# convert_pin_to_dramsim.awk
#
# Converts a raw Intel Pin "pinatrace" output file into the trace record format
# consumed by DRAMsim3.
#
# The essential fact about this converter, and about the study that depends on
# it, is that the input carries no timing information. A pinatrace record states
# which instruction touched memory, whether the access was a read or a write,
# which address it touched, and how many bytes it moved. It does not state when.
#
# This converter therefore manufactures the timing. It starts a cycle counter at
# zero and advances it by a fixed stride on every record, reads and writes
# alike. The stride is uniform, so the resulting trace presents requests to the
# simulator at a constant rate regardless of the behavior of the application
# that produced it.
#
# The consequence must be understood before any result is interpreted: the
# offered load seen by the simulated device is a property of this script, not of
# the workload. Two workloads with wildly different access patterns will be
# presented at the same rate. That is what makes the comparison in the paper a
# controlled one, and it is also why no statement about application-level
# request rate can be drawn from these runs.
#
# Usage:
#   gawk -f convert_pin_to_dramsim.awk -v stride=50 pinatrace.out > trace.dramsim3
#
# Output format, one record per line:
#   <hexadecimal address> <READ|WRITE> <cycle>

BEGIN {
    if (stride == "") stride = 50
    cycle = 0
    emitted = 0
    skipped = 0
}

# pinatrace writes a header and a trailer as comment lines.
/^#/ { skipped++; next }

# A well formed record has at least an instruction address, an access type, and
# a memory address. Anything shorter is not a record and is counted, not
# guessed at.
NF < 3 { skipped++; next }

{
    type = ($2 == "R" || $2 == "r") ? "READ" : "WRITE"
    print $3, type, cycle
    cycle += stride
    emitted++
}

END {
    printf("converted %d records, skipped %d non-record lines, ", emitted, skipped) > "/dev/stderr"
    printf("final synthetic cycle %d, stride %d\n", cycle, stride) > "/dev/stderr"
}
