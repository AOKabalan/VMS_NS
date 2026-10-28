#!/usr/bin/env bash
# Run the other exact-Gresho VMS convection forms; EMAC–VMS and SKEW1 are excluded.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

usage="Usage: $0 [--dry-run] [--ranks N]"
dry_run=0
ranks=4
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) dry_run=1 ;;
        --ranks)
            if [[ $# -lt 2 || ! $2 =~ ^[1-9][0-9]*$ ]]; then
                printf '%s\n' "$usage: --ranks requires a positive integer" >&2
                exit 2
            fi
            ranks=$2
            shift
            ;;
        --help|-h) printf '%s\n' "$usage"; exit 0 ;;
        *) printf '%s\n' "$usage: unknown argument '$1'" >&2; exit 2 ;;
    esac
    shift
done

forms=(nonconservative conservative skew_symmetric_1 emac_vms)
root=results/gresho_degree_2
log="$root/remaining_vms_queue.log"
for form in "${forms[@]}"; do
    config="configs/gresho_degree_2/${form}_oss.toml"
    if [[ ! -f $config ]]; then
        printf 'Missing configuration: %s\n' "$config" >&2
        exit 2
    fi
done
if [[ $dry_run == 0 ]]; then
    command -v python >/dev/null
    command -v mpiexec >/dev/null
    mkdir -p "$root"
fi

# Ctrl-C stops the queue rather than starting another case after an interrupted MPI run.
trap 'exit 130' INT
trap 'exit 143' TERM
completed=0
skipped=0
failed=0
for form in "${forms[@]}"; do
    name="${form}_oss"
    config="configs/gresho_degree_2/$name.toml"
    # An existing directory may belong to a running, completed, or interrupted run.
    # Never move it or rerun it automatically.
    if [[ -e $root/$name ]]; then
        printf 'SKIP %s (output already exists)\n' "$name"
        skipped=$((skipped + 1))
        continue
    fi
    if [[ $dry_run == 1 ]]; then
        printf 'OMP_NUM_THREADS=1 mpiexec -n %s python scripts/run_simulation.py %s\n' "$ranks" "$config"
        continue
    fi
    printf 'RUN %s (%s MPI ranks)\n' "$name" "$ranks"
    printf '%s\t%s\tstarted\n' "$(date -Is)" "$name" >> "$log"
    started=$SECONDS
    if OMP_NUM_THREADS=1 mpiexec -n "$ranks" python scripts/run_simulation.py "$config" </dev/null; then
        status=0
        completed=$((completed + 1))
    else
        status=$?
        failed=$((failed + 1))
    fi
    printf '%s\t%s\texit=%s\tseconds=%s\n' "$(date -Is)" "$name" "$status" "$((SECONDS - started))" >> "$log"
    printf 'FINISHED %s: exit %s\n' "$name" "$status"
done
if [[ $dry_run == 1 ]]; then exit 0; fi
printf 'Remaining VMS queue: %s completed, %s skipped, %s failed\n' "$completed" "$skipped" "$failed"
[[ $failed == 0 ]]
