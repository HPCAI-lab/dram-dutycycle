#!/bin/bash
cd /home/claude/work/gmx_run
s=$(date +%s)
valgrind --tool=lackey --trace-mem=yes --log-fd=9 /home/claude/gmx/build/bin/gmx -quiet mdrun -s md.tpr -deffnm traced 9>&1 1>program_stdout.txt 2>program_stderr.txt | /home/claude/work/gmx_trace/gmx_filter /home/claude/work/gmx_trace/syms.txt $(cat /home/claude/work/gmx_trace/marker.txt) 1 40000050 bridge.pin step1.bin profile.json
echo "status ${PIPESTATUS[*]} secs $(( $(date +%s)-s ))" > trace_done.txt
