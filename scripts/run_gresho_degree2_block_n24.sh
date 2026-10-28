#!/usr/bin/env bash
# Fresh 24x24 P2/P2 Gresho runs; the stopped 48x48 results stay in place.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source /home/ali/venv-firedrake/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export XDG_CACHE_HOME=/tmp/gresho_exact_cache
export PYOP2_CACHE_DIR=/tmp/gresho_exact_pyop2
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR=/tmp/gresho_exact_tsfc
mkdir -p results
exec 9>results/gresho_degree2_block_n24.lock
flock -n 9 || { echo 'This queue is already running.' >&2; exit 2; }
exec bash scripts/run_validation_queue.sh --configs configs/gresho_degree_2_block_n24 --runner simulation "$@"
