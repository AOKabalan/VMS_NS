# Local results layout

The existing results were grouped on 3 November 2026. The directory name
`bousinesq` follows the owner's requested spelling.

| Location | Contents |
| --- | --- |
| `results/bousinesq/cases/` | Thermal manufactured runs, heated cavities, free-slip checks and stratified KH runs |
| `results/bousinesq/verification/manufactured_solution/2026-11-03/` | Latest MMS/time-refinement figures, CSV/JSON measurements and raw test artifacts |
| `results/bousinesq/clouds/` | KH and stable gravity-wave cloud visualization, with grouped setup checks/previews |
| `results/bousinesq/configs/` | Previous local thermal smoke configuration |
| `results/navier_stokes/` | Isothermal benchmark runs, including unstratified KH |
| `results/analysis/` | Conservation studies, solver profiling and KH preparation experiments |
| `results/presentation/` | Existing presentation figures, data and generation scripts |
| `results/reports/` | Combined results report |
| `results/logs/` | Historical queue logs, completion markers, PID and lock files |

`results/README.md` and `results/bousinesq/README.md` provide local clickable
indexes. `results/relocation_manifest.json` records every old/new path and each
moved file's SHA-256 checksum. All 15,894 moved files were checked unchanged.
Local analysis scripts that depended on folder depth were subsequently adjusted;
the manifest records these edits and the preserved original script copies under
`results/logs/pre_relocation_scripts/`. Numerical results were not edited.
Incomplete and failed runs were retained with their original names and status.
Historical PID/lock files are archived records, not evidence of a running job.

Old case and analysis paths are compatibility symlinks, so existing commands,
provenance references and overwrite checks still find the same data. The solver
continues to use `results/<case>/<run_name>/`; for existing cases the links route
new runs into the grouped directories. New case names still use the original
default layout. No output defaults changed. The warm compiler cache at
`results/.b0-test-cache/` was left in place.

The downloaded external reference dataset remains under
`data/references/kelvin_helmholtz_unstratified/`. It is input reference data,
separate from simulation results. Committed verification reports remain in
`docs/verification/`.

## Regenerating verification figures

```bash
python scripts/plot_boussinesq_verification.py
```

This reads the committed
`docs/verification/boussinesq_mms_recheck_20261103.json` and writes three PNG/PDF
figures plus CSV/JSON data to the manufactured-solution folder above. It does
not run a simulation. `--summary` and `--output` select other paths.

Spatial plots show the coupled manufactured solution at P1/P2. The temporal
plots show the separate buoyancy-driven semidiscrete eigenmode check, which
isolates time discretization. Figure captions and the B5 report retain the
representable-temperature and zero-subscale qualifications. Results remain
gitignored; the plotting script and measured summary are versioned.
