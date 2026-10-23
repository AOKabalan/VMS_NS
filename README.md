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

The Firedrake version is pinned: the code is tested with Firedrake 2026.10.0 (PETSc 3.26,
UFL 2026.2). `solver/reusable_assembly.py` relies on Firedrake internals
(`BaseFormAssembler.preprocess_base_form` and the `is_base_form_preprocessed` flag of
`get_assembler`) to preprocess repeatedly assembled forms once. After any Firedrake upgrade,
`pytest tests/unit/test_reusable_assembly.py` must pass, since it checks that these assemblers
give exactly what `fd.assemble` and `Function.interpolate` give. Every run records the versions
it used in `run_info.json`.

## Run

```bash
python scripts/run_simulation.py configs/<file>.toml
mpiexec -n 4 python scripts/run_simulation.py configs/<file>.toml
```

A run is described by one TOML file (plan §7). Every field has a default except `case`, and
`run_name` defaults to the file name without extension. The run writes to
`results/<case>/<run_name>/`, relative to the working directory, and refuses to overwrite an
existing folder. It writes the files of plan §9: `config_resolved.toml`, `run_info.json`,
`run.log`, `inner_iterations.csv`, `outer_iterations.csv`, `time_steps.csv`, `diagnostics.csv`
and, when `vtk_every > 0`, `fields/velocity.pvd` and `fields/pressure.pvd`. The Crank–Nicolson
pressure is written at t^(n+½). An invalid configuration stops the run before anything is
written, with a message naming the offending fields. If a nonlinear solve does not converge, the
run writes the last state to `fields/`, logs the iteration history and exits with status 1.

Available diagnostics: `kinetic_energy`, `subscale_kinetic_energy` (VMS only),
`linear_momentum`, `angular_momentum`, `force`, `divergence_l2`, and, for cases with an exact
solution, `velocity_error_l2`, `velocity_gradient_error_l2` and `pressure_error_l2`.
All integrals are global and use the common quadrature rule; error norms use twice its degree.
Momentum vectors expand into `_x`, `_y`, `_z` columns; angular momentum in 2D is one scalar.
A case declares `Problem.force_boundaries`, a tuple of boundary names; `force` writes
`force_<boundary>_<axis>` for each. It measures the fluid-on-body force
−∫σn ds with σ = −pI + ν(∇u + ∇uᵀ), without coefficient normalization. A force request without
declared boundaries is rejected. Forces on an individual open face depend on the pressure gauge;
closed-boundary forces do not.

In `diagnostics.csv`, `time` labels the endpoint velocity diagnostics; `pressure_time` labels
the pressure error and force (including its velocity gradient). For Crank–Nicolson, these use
the half-level pressure and velocity, except at the backward-Euler startup step.
`subscale_time` labels subscale energy: dynamic subscales use their stored endpoint value;
static subscales use the converged solve-level value, stored for nonlinear splitting and
evaluated from the expression otherwise. This time is blank for Galerkin.

The `ethier_steinman_3d` case uses [−1, 1]³, a = π/4, d = π and ν = 1 by default, with exact
Dirichlet data and all six faces declared for force measurements. Its default mesh is 4³;
`mesh.cells = [8, 8, 8]` selects the finer smoke-test mesh. The three-step BDF2 check uses
δt = 0.01 and T = 0.03, for Galerkin P2/P1 and Dyn-Nl-OSS P1/P1, in serial and on two MPI ranks.
The choice d = π keeps the interpolated Dirichlet trace at zero net flux for P1 as well as P2.
With d = π/2 on the default tetrahedral mesh, the P1 trace has nonzero flux, making the
constant-test continuity equation inconsistent (evidence in the Phase 9 progress entry).
For all-Dirichlet problems, each solve level checks the interpolated boundary flux and rejects
data whose net flux exceeds 1e-10 times the integral of its absolute normal flux.
In 3D, results with VMS depend on the process count at the level of the quadrature error.

## Validation

Three validation cases compare against published data: `lid_driven_cavity` (Ghia, Ghia and
Shin 1982, Re = 100 and 1000), `cylinder_dfg_2d` (DFG 2D-2, FEATFLOW reference) and
`lattice_vortex` (exact kinetic-energy decay, MSc thesis §6.2). Their reference data and
sources are in `data/references/<case>/`; the cylinder meshes are generated from
`meshes/cylinder_dfg_2d.geo` (a coarse test mesh and production levels 1–3, see the `.geo`).

```bash
# One run and its comparison, written to validation.csv and validation.md in the results folder:
OMP_NUM_THREADS=1 mpiexec -n 4 python scripts/run_validation.py configs/validation/<file>.toml
# Re-analyse a finished run without running it again:
python scripts/run_validation.py configs/validation/<file>.toml --analyse-only
# Every config in configs/validation, shortest first, 4 ranks each (hours; the owner launches it):
bash scripts/run_validation_queue.sh
# The same queue with 3 steps per run, into results/validation_smoke:
bash scripts/run_validation_queue.sh --smoke
```

The queue appends one line per run to `results/validation_queue.log`, keeps going after a
failed run, skips runs that already have `validation.csv` and moves an interrupted run's folder
to `<folder>.incomplete-<timestamp>` before running it again. It ends by writing
`results/VALIDATION_DONE` and, where `notify-send` exists, a desktop notification.

## Test and lint

```bash
pytest tests/unit                      # seconds
OMP_NUM_THREADS=1 pytest tests/verification -m slow -n 4  # V1–V4 and the 3D serial/MPI check
ruff check . && ruff format --check .
```

A convergence study runs a base configuration on several meshes or step counts and tabulates the
final value and rate of every configured diagnostic in
`results/<case>/<run_name>_study/convergence.csv`; each level is an ordinary run in
`results/<case>/<run_name>_<level>/`. Run it serially:

```bash
python scripts/run_convergence_study.py configs/<file>.toml --cells 8x8 16x16 32x32
python scripts/run_convergence_study.py configs/<file>.toml --steps 20 40 80 160
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
| L_a, L*_a | `momentum_operator`; L*_a inside `vms_terms` | Operators; the convective part of L*_a is `convective_subscale_coupling` |
| Δ_h | `elementwise_laplacian` | Laplacian evaluated cell by cell |
| c(a; u, v) | `convection_form` | Galerkin convective form |
| D_t | `time_derivative` | Discrete time derivative |
| α₀, g_hist | `alpha0`, `history` | Time-scheme coefficients |
| χ_A, χ_D | `is_asgs`, `is_dynamic` | Spec switches |
