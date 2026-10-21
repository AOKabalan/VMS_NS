# navier-stokes-vms

Variational multiscale (VMS) stabilized finite element solvers for the incompressible
Navier–Stokes equations in [Firedrake](https://www.firedrakeproject.org): a Taylor–Hood Galerkin
baseline and the eight VMS cases (ASGS or OSS × static or dynamic subscales × linear or nonlinear
scale splitting), in 2D and 3D.

- `docs/formulation_spec.md` defines every equation and is the source of truth.
- `docs/implementation_plan.md` fixes the code structure, tests and phases.
- `docs/verification/VERIFICATION.md` records the independent checks of the spec.
- `docs/progress_log.md` records each completed phase.

## Install

Firedrake is installed separately, following its own instructions. Inside the activated
Firedrake virtual environment, from the repository root:

```bash
pip install -e ".[dev]"
```

## Run

```bash
python scripts/run_simulation.py configs/<file>.toml
mpiexec -n 4 python scripts/run_simulation.py configs/<file>.toml
```

A run is described by one TOML file (plan §7). Every field has a default except `case`, and
`run_name` defaults to the file name without extension. The run writes to
`results/<case>/<run_name>/`, relative to the working directory, and refuses to overwrite an
existing folder. It writes `config_resolved.toml`, `run_info.json` and `run.log` (plan §9). An
invalid configuration stops the run with a message naming the offending fields.

## Test and lint

```bash
pytest tests/unit
ruff check . && ruff format --check .
```

## Spec symbols in code

The `formulation/` layer uses the spec's symbols so that it reads side by side with the spec.

| Spec | Code | Meaning |
| --- | --- | --- |
| u_h, p_h | `u`, `p` | FE velocity and pressure (trial functions or current iterate) |
| v_h, q_h | `v`, `q` | Test functions |
| a | `a` | Advection velocity |
| ũ | `u_sub` | Velocity subscale expression (UFL) |
| ũ* | `u_sub_iterate` | Stored subscale iterate inside `a` (nonlinear splitting) |
| ũ_hist | `u_sub_history` | BDF history combination of stored subscales |
| p̃ | `p_sub` | Pressure subscale expression |
| τ_m, τ_c, τ_t | `tau_m`, `tau_c`, `tau_t` | Stabilization parameters |
| ξ_h, ζ_h | `xi`, `zeta` | OSS projections |
| R_u, R̂_u, R_p | `momentum_residual`, `reduced_momentum_residual`, `continuity_residual` | Residuals |
| L_a, L*_a | `momentum_operator`, `adjoint_operator` | Operators |
| Δ_h | `elementwise_laplacian` | Laplacian evaluated cell by cell |
| c(a; u, v) | `convection_form` | Galerkin convective form |
| D_t | `time_derivative` | Discrete time derivative |
| α₀, g_hist | `alpha0`, `history` | Time-scheme coefficients |
| χ_A, χ_D | `is_asgs`, `is_dynamic` | Spec switches |
