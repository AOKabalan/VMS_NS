#!/usr/bin/env bash
# Ordinary EMAC Galerkin first, then the owner's degree-2 VMS configurations.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source /home/ali/venv-firedrake/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export XDG_CACHE_HOME=/tmp/gresho_exact_cache
export PYOP2_CACHE_DIR=/tmp/gresho_exact_pyop2
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR=/tmp/gresho_exact_tsfc
mkdir -p results
exec 9>results/gresho_galerkin_then_degree2.lock
flock -n 9 || { echo 'This queue is already running.' >&2; exit 2; }
configs=(configs/gresho_exact/emac_galerkin.toml configs/gresho_degree_2/*.toml)

# Reject invalid configs before launching any run; never change the owner's run settings.
python - "${configs[@]}" <<'PY'
import sys
from pathlib import Path
from navier_stokes_vms.config import load_config, check_diagnostics_known
from navier_stokes_vms.cases import CASES
from navier_stokes_vms.postprocessing.diagnostics import DIAGNOSTICS
for argument in sys.argv[1:]:
    config = load_config(Path(argument))
    assert config.case in CASES, config.case
    check_diagnostics_known(config, DIAGNOSTICS)
    print(f'Validated {argument}', flush=True)
PY

log=results/gresho_galerkin_then_degree2_queue.log
failed=0
for config in "${configs[@]}"; do
    folder=$(python - "$config" <<'PY'
import sys
from pathlib import Path
from navier_stokes_vms.config import load_config
c = load_config(Path(sys.argv[1]))
print(f'results/{c.case}/{c.run_name}')
PY
)
    if [[ -e $folder ]]; then
        printf '%s\t%s\tskipped: output exists\n' "$(date -Is)" "$config" | tee -a "$log"
        continue
    fi
    printf '%s\t%s\tstarted\n' "$(date -Is)" "$config" | tee -a "$log"
    if mpiexec -n 4 python scripts/run_simulation.py "$config" </dev/null; then
        status=0
    else
        status=$?
        failed=$((failed + 1))
    fi
    printf '%s\t%s\texit=%s\n' "$(date -Is)" "$config" "$status" | tee -a "$log"
done
printf '%s\tqueue finished; failed=%s\n' "$(date -Is)" "$failed" | tee -a "$log"
[[ $failed == 0 ]]
