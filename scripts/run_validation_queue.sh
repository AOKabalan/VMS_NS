#!/usr/bin/env bash
# Sequential, resumable validation queue; the owner launches the full queue.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
mode=()
root=results
if [[ ${1:-} == --smoke && $# == 1 ]]; then
    mode=(--smoke)
    root=results/validation_smoke
elif [[ $# != 0 ]]; then
    echo "Usage: $0 [--smoke]" >&2
    exit 2
fi
mkdir -p "$root"
rm -f "$root/VALIDATION_DONE"
log="$root/validation_queue.log"
failures=0
completed=0
skipped=0
# Within each case group, cheaper runs first: mesh size (cells, or the gmsh file's size) times
# degree² times the number of steps. The order is deterministic.
while IFS=$'\t' read -r config case_name run_name; do
    folder="$root/$case_name/$run_name"
    if [[ -f "$folder/validation.csv" ]]; then
        printf '%s\t0\t0\tskipped\n' "$config" >> "$log"
        skipped=$((skipped+1))
        continue
    fi
    if [[ -e "$folder" ]]; then
        mv -- "$folder" "$folder.incomplete-$(date +%Y%m%dT%H%M%S%N)" || exit 1
    fi
    start=$SECONDS
    # mpiexec reads stdin; without </dev/null it consumes the remaining queue entries.
    OMP_NUM_THREADS=1 mpiexec -n 4 python scripts/run_validation.py "$config" \
        --run-only --results-root "$root" "${mode[@]}" </dev/null
    status=$?
    if [[ $status == 0 ]]; then
        OMP_NUM_THREADS=1 python scripts/run_validation.py "$config" --analyse-only \
            --results-root "$root" </dev/null
        status=$?
    fi
    printf '%s\t%d\t%d\tran\n' "$config" "$status" "$((SECONDS-start))" >> "$log"
    if [[ $status == 0 ]]; then completed=$((completed+1)); else failures=$((failures+1)); fi
done < <(python - <<'PY'
from pathlib import Path
import tomllib
entries = []
groups = {'lid_driven_cavity':0, 'lattice_vortex':1, 'cylinder_dfg_2d':2}
for path in Path('configs/validation').glob('*.toml'):
    config = tomllib.loads(path.read_text())
    mesh = config.get('mesh', {})
    cells = mesh['cells'][0] * mesh['cells'][1] if 'cells' in mesh else Path(mesh['file']).stat().st_size
    degree = config.get('discretization',{}).get('velocity_degree',1)
    steps = config['time']['t_end'] / config['time']['dt']
    size = cells * degree**2 * steps
    entries.append(((groups[config['case']],size,str(path)), path,config['case'],config['run_name']))
for _,path,case,name in sorted(entries):
    print(f'{path}\t{case}\t{name}')
PY
)
summary="Validation queue: $completed completed, $skipped skipped, $failures failed"
printf '%s\n' "$summary" > "$root/VALIDATION_DONE"
if command -v notify-send >/dev/null 2>&1; then notify-send "$summary" || true; fi
printf '%s\n' "$summary"
[[ $failures == 0 ]]
