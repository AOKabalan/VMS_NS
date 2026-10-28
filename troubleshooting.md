# Gresho P2 block-solver troubleshooting

Recorded before removing trial outputs at the user’s request on 2026-10-26.

## Production results retained

All four `results/gresho_degree_2/*_n24_block_gauge_oss/` runs completed 1000 steps through t=10 with exit status 0. These are the results used in the presentation plots. Their full fields, diagnostic histories and provenance are retained.

The original 48×48 CONS run in `results/gresho_degree_2/conservative_oss/` was explicitly stopped by the user at step 795, t=7.95. Its saved results remain intact (not a restart checkpoint). The completed 48×48 P2/P1 `results/gresho_exact/emac_galerkin/` run is also retained.

## Diagnosis and correction

- The first 24×24 attempts reused the preconditioner for 20 Newton iterations and used adaptive linear tolerances. CONS failed its linear solve at step 2; EMAC-VMS reached 23 successful steps before the nonlinear solve stalled at step 24. The CONV trial was stopped manually after 10 successful steps while the solver settings were investigated; this is not a recorded convergence failure.
- Fresh preconditioners every Newton iteration, fixed FGMRES rtol=1e-10 / atol=1e-14 and a 300-iteration limit did not suffice: CONS failed at step 1 and EMAC-VMS at step 2.
- The first diagnostic invocation incorrectly supplied PETSc viewer options as Boolean true, causing an unsupported-viewer error. Empty-string viewer options fixed that instrumentation error.
- The subsequent diagnostic CONS run hit DIVERGED_ITS at 300 linear iterations. The recorded true residual stagnated around 2.25e-11 (relative about 2.41e-9), above the requested linear tolerance.
- Enabling `fieldsplit_1_mat_mumps_icntl_24=1` (null-pivot detection) in the velocity–pressure block resolved these observed failures. This supports a pressure-gauge/preconditioner issue; failed linear convergence alone is not evidence of physical divergence.
- With that setting, CONS and EMAC-VMS passed 50-step checks (t=0.5), and CONV and SKEW1 passed 10-step checks (t=0.1). All four subsequent production runs completed t=10.
- The block split is xi (field 2), solved with CG/Jacobi, versus coupled velocity–pressure (fields 0,1), solved with MUMPS LU; the outer solver is FGMRES with full Schur factorization. This is not a fully iterative velocity/pressure split.
- No viscosity, boundary-condition, formulation, conservation-diagnostic or nonlinear-tolerance changes were made to fix the block solver. The smaller runs use 24×24 P2/P2, zero viscosity, free-slip walls, pure Crank–Nicolson, dt=0.01. They start fresh, rather than continuing the stopped 48×48 state.

## Removed trial runs

| Folder under `results/gresho_degree_2/` | Exit status | Successful steps | Last successful time | Requested end time |
|---|---:|---:|---:|---:|
| `conservative_oss_n24_block_diagnose` | 1 | 0 | none | 0.02 |
| `conservative_oss_n24_block_diagnose2` | 1 | 1 | 0.01 | 0.02 |
| `conservative_oss_n24_block_fresh_check` | 1 | 0 | none | 0.5 |
| `conservative_oss_n24_block_nullpivot_check` | 0 | 50 | 0.5 | 0.5 |
| `conservative_oss_n24_block_oss` | 1 | 1 | 0.01 | 10.0 |
| `emac_vms_oss_n24_block_fresh_check` | 1 | 1 | 0.01 | 0.5 |
| `emac_vms_oss_n24_block_nullpivot_check` | 0 | 50 | 0.5 | 0.5 |
| `emac_vms_oss_n24_block_oss` | 1 | 23 | 0.23 | 10.0 |
| `nonconservative_oss_n24_block_nullpivot_check` | 0 | 10 | 0.1 | 0.1 |
| `nonconservative_oss_n24_block_oss` | None | 10 | 0.1 | 10.0 |
| `skew_symmetric_1_oss_n24_block_nullpivot_check` | 0 | 10 | 0.1 | 0.1 |

`None` means the interrupted trial did not record an exit status. Failed-step field dumps, if any, are not counted as successful steps.

Also removed `results/gresho_degree_2_block_checks/`, which contained temporary check configs and console logs. The configuration differences and diagnostic evidence are preserved here.

## Exact production reference configuration (CONS)

```toml
case = "gresho_degree_2"
run_name = "conservative_oss_n24_block_gauge_oss"

[case_parameters]
viscosity = 0.0

[mesh]
cells = [
    24,
    24,
]

[discretization]
velocity_degree = 2
pressure_degree = 2

[time]
scheme = "crank_nicolson"
dt = 0.01
t_end = 10.0
startup = "none"

[stabilization]
method = "vms"
subscale_space = "oss"
subscale_dynamics = "dynamic"
scale_splitting = "nonlinear"
c1 = 4.0
c2 = 2.0
pressure_subscale = false

[convection]
form = "conservative"

[viscosity]
form = "stress"

[solver]
linearization = "newton"
newton_tau = "frozen"
picard_iterations_before_newton = 0
inner_rtol = 1e-12
inner_atol = 1e-12
inner_max_iterations = 200
outer_rtol = 1e-11
outer_max_iterations = 100
subscale_update = "accelerated"
subscale_rtol = 1e-13
subscale_max_iterations = 50

[solver.petsc]
mat_type = "aij"
ksp_type = "fgmres"
ksp_pc_side = "right"
ksp_norm_type = "unpreconditioned"
ksp_rtol = 1e-10
ksp_atol = 1e-14
ksp_max_it = 300
ksp_gmres_restart = 50
snes_lag_jacobian = 1
snes_lag_preconditioner = 1
snes_lag_preconditioner_persists = true
pc_type = "fieldsplit"
pc_fieldsplit_type = "schur"
pc_fieldsplit_schur_fact_type = "full"
pc_fieldsplit_schur_precondition = "a11"
pc_fieldsplit_0_fields = "2"
pc_fieldsplit_1_fields = "0,1"
fieldsplit_0_ksp_type = "cg"
fieldsplit_0_ksp_rtol = 0.01
fieldsplit_0_ksp_atol = 1e-14
fieldsplit_0_ksp_max_it = 50
fieldsplit_0_pc_type = "jacobi"
fieldsplit_1_ksp_type = "preonly"
fieldsplit_1_pc_type = "lu"
fieldsplit_1_pc_factor_mat_solver_type = "mumps"
fieldsplit_1_mat_mumps_icntl_24 = 1

[output]
vtk_every = 100
diagnostics = [
    "kinetic_energy",
    "linear_momentum",
    "angular_momentum",
    "velocity_error_l2",
    "pressure_error_l2",
    "divergence_l2",
    "force",
    "region_velocity_l2",
    "total_kinetic_energy",
    "subscale_kinetic_energy",
    "region_subscale_l2",
    "energy_budget",
    "momentum_budget",
    "linear_momentum_defect",
    "angular_momentum_defect",
    "convective_energy_defect",
]
console_level = "INFO"
```

## Trial configuration differences

Each entry compares the removed trial’s resolved configuration with its own formulation’s retained production configuration. Omitted settings are identical. JSON `<absent>` denotes a missing option, not a literal solver value.

### conservative_oss_n24_block_diagnose

```json
{
  "output.vtk_every": {
    "trial": 0,
    "production": 100
  },
  "run_name": {
    "trial": "conservative_oss_n24_block_diagnose",
    "production": "conservative_oss_n24_block_gauge_oss"
  },
  "solver.petsc.fieldsplit_1_ksp_converged_reason": {
    "trial": true,
    "production": "<absent>"
  },
  "solver.petsc.fieldsplit_1_mat_mumps_icntl_24": {
    "trial": "<absent>",
    "production": 1
  },
  "solver.petsc.ksp_converged_reason": {
    "trial": true,
    "production": "<absent>"
  },
  "solver.petsc.ksp_monitor_true_residual": {
    "trial": true,
    "production": "<absent>"
  },
  "time.t_end": {
    "trial": 0.02,
    "production": 10.0
  }
}
```

### conservative_oss_n24_block_diagnose2

```json
{
  "output.vtk_every": {
    "trial": 0,
    "production": 100
  },
  "run_name": {
    "trial": "conservative_oss_n24_block_diagnose2",
    "production": "conservative_oss_n24_block_gauge_oss"
  },
  "solver.petsc.fieldsplit_1_ksp_converged_reason": {
    "trial": "",
    "production": "<absent>"
  },
  "solver.petsc.fieldsplit_1_mat_mumps_icntl_24": {
    "trial": "<absent>",
    "production": 1
  },
  "solver.petsc.ksp_converged_reason": {
    "trial": "",
    "production": "<absent>"
  },
  "solver.petsc.ksp_monitor_true_residual": {
    "trial": "",
    "production": "<absent>"
  },
  "time.t_end": {
    "trial": 0.02,
    "production": 10.0
  }
}
```

### conservative_oss_n24_block_fresh_check

```json
{
  "output.vtk_every": {
    "trial": 0,
    "production": 100
  },
  "run_name": {
    "trial": "conservative_oss_n24_block_fresh_check",
    "production": "conservative_oss_n24_block_gauge_oss"
  },
  "solver.petsc.fieldsplit_1_mat_mumps_icntl_24": {
    "trial": "<absent>",
    "production": 1
  },
  "time.t_end": {
    "trial": 0.5,
    "production": 10.0
  }
}
```

### conservative_oss_n24_block_nullpivot_check

```json
{
  "output.vtk_every": {
    "trial": 0,
    "production": 100
  },
  "run_name": {
    "trial": "conservative_oss_n24_block_nullpivot_check",
    "production": "conservative_oss_n24_block_gauge_oss"
  },
  "time.t_end": {
    "trial": 0.5,
    "production": 10.0
  }
}
```

### conservative_oss_n24_block_oss

```json
{
  "run_name": {
    "trial": "conservative_oss_n24_block_oss",
    "production": "conservative_oss_n24_block_gauge_oss"
  },
  "solver.petsc.fieldsplit_1_mat_mumps_icntl_24": {
    "trial": "<absent>",
    "production": 1
  },
  "solver.petsc.ksp_max_it": {
    "trial": 200,
    "production": 300
  },
  "solver.petsc.snes_ksp_ew": {
    "trial": true,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_rtol0": {
    "trial": 0.1,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_rtolmax": {
    "trial": 0.5,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_version": {
    "trial": 2,
    "production": "<absent>"
  },
  "solver.petsc.snes_lag_preconditioner": {
    "trial": 20,
    "production": 1
  }
}
```

### emac_vms_oss_n24_block_fresh_check

```json
{
  "output.vtk_every": {
    "trial": 0,
    "production": 100
  },
  "run_name": {
    "trial": "emac_vms_oss_n24_block_fresh_check",
    "production": "emac_vms_oss_n24_block_gauge_oss"
  },
  "solver.petsc.fieldsplit_1_mat_mumps_icntl_24": {
    "trial": "<absent>",
    "production": 1
  },
  "time.t_end": {
    "trial": 0.5,
    "production": 10.0
  }
}
```

### emac_vms_oss_n24_block_nullpivot_check

```json
{
  "output.vtk_every": {
    "trial": 0,
    "production": 100
  },
  "run_name": {
    "trial": "emac_vms_oss_n24_block_nullpivot_check",
    "production": "emac_vms_oss_n24_block_gauge_oss"
  },
  "time.t_end": {
    "trial": 0.5,
    "production": 10.0
  }
}
```

### emac_vms_oss_n24_block_oss

```json
{
  "run_name": {
    "trial": "emac_vms_oss_n24_block_oss",
    "production": "emac_vms_oss_n24_block_gauge_oss"
  },
  "solver.petsc.fieldsplit_1_mat_mumps_icntl_24": {
    "trial": "<absent>",
    "production": 1
  },
  "solver.petsc.ksp_max_it": {
    "trial": 200,
    "production": 300
  },
  "solver.petsc.snes_ksp_ew": {
    "trial": true,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_rtol0": {
    "trial": 0.1,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_rtolmax": {
    "trial": 0.5,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_version": {
    "trial": 2,
    "production": "<absent>"
  },
  "solver.petsc.snes_lag_preconditioner": {
    "trial": 20,
    "production": 1
  }
}
```

### nonconservative_oss_n24_block_nullpivot_check

```json
{
  "output.vtk_every": {
    "trial": 0,
    "production": 100
  },
  "run_name": {
    "trial": "nonconservative_oss_n24_block_nullpivot_check",
    "production": "nonconservative_oss_n24_block_gauge_oss"
  },
  "time.t_end": {
    "trial": 0.1,
    "production": 10.0
  }
}
```

### nonconservative_oss_n24_block_oss

```json
{
  "run_name": {
    "trial": "nonconservative_oss_n24_block_oss",
    "production": "nonconservative_oss_n24_block_gauge_oss"
  },
  "solver.petsc.fieldsplit_1_mat_mumps_icntl_24": {
    "trial": "<absent>",
    "production": 1
  },
  "solver.petsc.ksp_max_it": {
    "trial": 200,
    "production": 300
  },
  "solver.petsc.snes_ksp_ew": {
    "trial": true,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_rtol0": {
    "trial": 0.1,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_rtolmax": {
    "trial": 0.5,
    "production": "<absent>"
  },
  "solver.petsc.snes_ksp_ew_version": {
    "trial": 2,
    "production": "<absent>"
  },
  "solver.petsc.snes_lag_preconditioner": {
    "trial": 20,
    "production": 1
  }
}
```

### skew_symmetric_1_oss_n24_block_nullpivot_check

```json
{
  "output.vtk_every": {
    "trial": 0,
    "production": 100
  },
  "run_name": {
    "trial": "skew_symmetric_1_oss_n24_block_nullpivot_check",
    "production": "skew_symmetric_1_oss_n24_block_gauge_oss"
  },
  "time.t_end": {
    "trial": 0.1,
    "production": 10.0
  }
}
```

## Selected original console evidence

### conservative.log

```text
   DIVERGED_LINEAR_SOLVE. Iteration history:
2026-10-26 05:04:19 ERROR   Run stopped: step 1 did not converge: Nonlinear solve failed to converge after 1 nonlinear iterations.
   DIVERGED_LINEAR_SOLVE; last state written to results/gresho_degree_2/conservative_oss_n24_block_fresh_check/fields at t = 0.005
```

### diagnose.log

```text
[0] Unsupported viewer true
```

### diagnose2.log

```text
  300 KSP unpreconditioned resid norm 4.070041624456e-12 true resid norm 2.252965108664e-11 ||r(i)||/||b|| 2.414152120816e-09
    Linear firedrake_1_ solve did not converge due to DIVERGED_ITS iterations 300
   DIVERGED_LINEAR_SOLVE. Iteration history:
2026-10-26 05:05:56 ERROR   Run stopped: step 2 did not converge: Nonlinear solve failed to converge after 0 nonlinear iterations.
   DIVERGED_LINEAR_SOLVE; last state written to results/gresho_degree_2/conservative_oss_n24_block_diagnose2/fields at t = 0.015
```

### emac_vms.log

```text
   DIVERGED_LINEAR_SOLVE. Iteration history:
2026-10-26 05:04:24 ERROR   Run stopped: step 2 did not converge: Nonlinear solve failed to converge after 1 nonlinear iterations.
   DIVERGED_LINEAR_SOLVE; last state written to results/gresho_degree_2/emac_vms_oss_n24_block_fresh_check/fields at t = 0.015
```

## Diagnostic replay for the EMAC–VMS/SKEW1 dissipation investigation

Production runs were not changed. Auxiliary projection checks initially failed
in a separate SKEW replay. The diagnostic quadrature-buffer interpolation of
`-grad(p_h)` produced infinite values even though direct integration of its
squared norm was finite. Changing the diagnostic mass solver alone did not fix
this. Direct weak-form assembly of each projection RHS, without this intermediate
interpolation, passed finite-value checks and reproduced the actual subscale
squared norm to 3.4e-21. The valid SKEW measurements are in
`results/conservation_analysis/gresho_emac_skew_cause/replay_skew_symmetric_1_direct/`.
The invalid sparse file is named `sampled.invalid.json` in the earlier replay
and is excluded by the analysis script. The separate analysis wrapper also had
an MPI folder-existence race, fixed using a collective broadcast before creating
the output directory. Detailed methods and limitations are in
`results/conservation_analysis/gresho_emac_skew_cause/methods.md`.
