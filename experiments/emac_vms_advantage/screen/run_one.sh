#!/bin/bash
# One screening run: serial, 25-minute cap, one line per run in runs.log.
config=$1
name=$(basename "$config" .toml)
start=$(date +%s)
OMP_NUM_THREADS=1 timeout 1500 python /home/user/repo/scripts/run_simulation.py "$config" \
    > "logs/$name.console" 2>&1 </dev/null
status=$?
echo "$name $status $(( $(date +%s) - start ))" >> runs.log
