#!/bin/bash
start=$(date +%s)
OMP_NUM_THREADS=1 timeout 7000 python driver.py "$@" > /dev/null 2>&1 </dev/null
echo "$(echo "$@" | tr ' ' _) $? $(( $(date +%s) - start ))" >> jobs_smooth.log
