#!/bin/bash
start=$(date +%s)
OMP_NUM_THREADS=1 timeout 7000 python driver_e2.py $1 $2 > /dev/null 2>&1 </dev/null
echo "E2_$1_$2 $? $(( $(date +%s) - start ))" >> jobs_e2.log
