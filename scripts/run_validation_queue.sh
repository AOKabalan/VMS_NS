#!/usr/bin/env bash
# Sequential, resumable run queue; the owner launches every full queue.
# Default: each config in configs/validation through scripts/run_validation.py.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
usage="Usage: $0 [--smoke] [--configs DIR] [--runner validation|simulation] [--dry-run]"
configs=configs/validation
runner=validation
smoke=0
dry_run=0
while [[ $# -gt 0 ]]; do
    case $1 in
        --smoke) smoke=1 ;;
        --dry-run) dry_run=1 ;;
        --configs|--runner)
            if [[ $# -lt 2 ]]; then echo "$usage" >&2; exit 2; fi
            if [[ $1 == --configs ]]; then configs=${2%/}; else runner=$2; fi
            shift ;;
        *) echo "$usage" >&2; exit 2 ;;
    esac
    shift
done
if [[ $runner != validation && $runner != simulation ]]; then
    echo "$usage: unknown runner '$runner'" >&2
    exit 2
fi
if [[ ! -d $configs ]]; then
    echo "$usage: no config directory '$configs'" >&2
    exit 2
fi
# run_simulation.py always writes to results/; only the validation runner has a smoke mode.
if [[ $smoke == 1 && $runner != validation ]]; then
    echo "$usage: --smoke needs the validation runner" >&2
    exit 2
fi
queue=$(basename "$configs")
mode=()
root=results
if [[ $smoke == 1 ]]; then
    mode=(--smoke)
    root=results/${queue}_smoke
fi
log="$root/${queue}_queue.log"
done_file="$root/${queue^^}_DONE"
if [[ $dry_run == 0 ]]; then
    mkdir -p "$root"
    rm -f "$done_file"
fi
failures=0
completed=0
skipped=0
# A validation run is complete once it has validation.csv. A simulation run is complete once
# run_info.json records an exit status, so a run that stopped without converging is kept as a
# recorded failure and not repeated.
finished() {
    if [[ $runner == validation ]]; then
        [[ -f "$1/validation.csv" ]]
    else
        grep -Eq '"exit_status": [0-9]+' "$1/run_info.json" 2>/dev/null
    fi
}
# Cases in a fixed group order, then any other case alphabetically; within each case group,
# cheaper runs first. Estimated cost: vertices × unknowns per vertex (P2 has about 4 nodes per
# vertex; OSS doubles the fields with ξ, ζ) × outer iterations per step (about 7 with nonlinear
# splitting) × steps. The order is deterministic.
while IFS=$'\t' read -r config case_name run_name; do
    folder="$root/$case_name/$run_name"
    if finished "$folder"; then
        if [[ $dry_run == 1 ]]; then
            printf '%s\tskip\n' "$config"
            continue
        fi
        printf '%s\t0\t0\tskipped\n' "$config" >> "$log"
        skipped=$((skipped+1))
        continue
    fi
    if [[ $dry_run == 1 ]]; then
        if [[ -e "$folder" ]]; then printf '%s\tmove incomplete, run\n' "$config"
        else printf '%s\trun\n' "$config"; fi
        continue
    fi
    if [[ -e "$folder" ]]; then
        mv -- "$folder" "$folder.incomplete-$(date +%Y%m%dT%H%M%S%N)" || exit 1
    fi
    start=$SECONDS
    # mpiexec reads stdin; without </dev/null it consumes the remaining queue entries.
    if [[ $runner == validation ]]; then
        OMP_NUM_THREADS=1 mpiexec -n 4 python scripts/run_validation.py "$config" \
            --run-only --results-root "$root" "${mode[@]}" </dev/null
        status=$?
        if [[ $status == 0 ]]; then
            OMP_NUM_THREADS=1 python scripts/run_validation.py "$config" --analyse-only \
                --results-root "$root" </dev/null
            status=$?
        fi
    else
        OMP_NUM_THREADS=1 mpiexec -n 4 python scripts/run_simulation.py "$config" </dev/null
        status=$?
    fi
    printf '%s\t%d\t%d\tran\n' "$config" "$status" "$((SECONDS-start))" >> "$log"
    if [[ $status == 0 ]]; then completed=$((completed+1)); else failures=$((failures+1)); fi
done < <(python - "$configs" <<'PY'
from pathlib import Path
import sys
import tomllib
groups = {'lid_driven_cavity': 0, 'lattice_vortex': 1, 'cylinder_dfg_2d': 2}
nodes_per_vertex = {1: 1, 2: 4}

def vertices(mesh):
    if 'cells' in mesh:
        return (mesh['cells'][0] + 1) * (mesh['cells'][1] + 1)
    lines = Path(mesh['file']).read_text().splitlines()
    return int(lines[lines.index('$Nodes') + 1])

entries = []
for path in Path(sys.argv[1]).glob('*.toml'):
    config = tomllib.loads(path.read_text())
    degrees, stabilization = config['discretization'], config['stabilization']
    vms = stabilization['method'] == 'vms'
    fields = 2 * nodes_per_vertex[degrees['velocity_degree']] + nodes_per_vertex[degrees['pressure_degree']]
    if vms and stabilization['subscale_space'] == 'oss':
        fields *= 2
    outer = 7 if vms and stabilization['scale_splitting'] == 'nonlinear' else 1
    steps = config['time']['t_end'] / config['time']['dt']
    cost = vertices(config['mesh']) * fields * outer * steps
    case = config['case']
    entries.append(((groups.get(case, len(groups)), case, cost, str(path)), path, case, config['run_name']))
for _, path, case, name in sorted(entries):
    print(f'{path}\t{case}\t{name}')
PY
)
if [[ $dry_run == 1 ]]; then exit 0; fi
summary="${queue^} queue: $completed completed, $skipped skipped, $failures failed"
printf '%s\n' "$summary" > "$done_file"
if command -v notify-send >/dev/null 2>&1; then notify-send "$summary" || true; fi
printf '%s\n' "$summary"
[[ $failures == 0 ]]
