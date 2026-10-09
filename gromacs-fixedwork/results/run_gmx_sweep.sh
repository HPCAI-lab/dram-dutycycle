#!/bin/bash
cd /home/claude/work/repo_gmx
s=$(date +%s)
python3 sweep/run_sweep.py --windows-dir data/windows --out results/runs --jobs 2 --windows g1_nb25,g2_nb50,g3_nb75,g4_spread,g5_fft,g6_gather > sweep.log 2>&1
python3 sweep/run_sweep.py --windows-dir data/windows --out results/runs --jobs 2 --windows g2_nb50_1M --arrivals uniform >> sweep.log 2>&1
echo "DONE $(( $(date +%s)-s )) s" >> sweep.log
