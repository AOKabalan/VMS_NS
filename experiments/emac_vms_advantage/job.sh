#!/bin/bash
start=$(date +%s)
OMP_NUM_THREADS=1 timeout 3600 python driver.py $1 $2 $3 > /dev/null 2>&1 </dev/null
echo "$1_$2_$3 $? $(( $(date +%s) - start ))" >> jobs.log
