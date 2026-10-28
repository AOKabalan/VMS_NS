#!/usr/bin/env bash
# Run sequentially to avoid doubling the memory required by the direct solver.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export FIREDRAKE_DISABLE_OPTIONS_LEFT=1
case "${LINEAR_SOLVER:-direct}" in
    direct) suffix="" ;;
    reuse_lu) suffix="_reuse_lu" ;;
    reuse_lu_tuned) suffix="_reuse_lu_tuned" ;;
    block_oss) suffix="_block_oss" ;;
    *) echo "LINEAR_SOLVER must be direct, reuse_lu, reuse_lu_tuned or block_oss" >&2; exit 2 ;;
esac
for form in skew1 emac_vms; do
    mpiexec -n "${NP:-4}" python experiments/kelvin_helmholtz_2020/run.py \
        "experiments/kelvin_helmholtz_2020/configs/kh2020_${form}_p1_n196${suffix}.toml" "$@"
done
