#!/bin/bash
# Matched-dissipation runs: τ_m scaled by s with τ_c fixed (README: open question).
cd "$(dirname "$0")"
. /opt/fd/venv312/bin/activate
export OMP_NUM_THREADS=1
run() {
    local start
    start=$(date +%s)
    bash -c "$1" > /dev/null 2>&1 < /dev/null
    echo "$? $(( $(date +%s) - start )) $1" >> jobs_matched.log
}
export -f run
xargs -P 4 -I{} bash -c 'run "$1"' _ {} < jobs_matched.txt
echo MATCHED_DONE > matched.state
