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
`run.log`, `inner_iterations.csv`, `outer_iterations.csv`, `time_steps.csv`, `diagnostics.csv`,
`initial_state.json` (|Ω|, the dimension, and the kinetic energy and linear and angular momentum
of the interpolated initial velocity) and, when `vtk_every > 0`, `fields/velocity.pvd` and
`fields/pressure.pvd`. The Crank–Nicolson
pressure is written at t^(n+½). An invalid configuration stops the run before anything is
written, with a message naming the offending fields. If a nonlinear solve does not converge, the
run writes the last state to `fields/`, logs the iteration history and exits with status 1.

`[convection] form` selects the convective form of the spec section "Convective forms
(Chapter 5)": `nonconservative`, `conservative`, `skew_symmetric_1` (default),
`skew_symmetric_2`, `emac` or `emac_vms`. `[viscosity] form` is `laplacian` (default) or
`stress`. `emac_vms` requires `stabilization.method = "vms"` and `viscosity.form = "stress"`.
`emac` and `emac_vms` require Newton with `picard_iterations_before_newton = 0`. Every form
works with both scale splittings and with Galerkin, except `emac_vms` with Galerkin.

With `emac` and `emac_vms` the pressure unknown is the modified pressure P_h, and the null space
and its zero-mean normalization act on P_h. Every reported pressure is the physical
p_h = P_h + ½|u_h|², with u_h the velocity of the solve. For Crank–Nicolson that is the half-level
velocity, not an average of the endpoint corrections. When the null space is active, p_h gets its
own zero mean. `pressure_error_l2` and `force` use this expression directly.
`fields/pressure.pvd` and `TimeStepper.pressure` hold it exactly in continuous Lagrange of degree
max(k_p, 2k_u), not projected onto the pressure space. With an outflow boundary, the natural
condition of these two forms is written in P_h, so their outflow results differ from the other
four forms.

Available diagnostics: `kinetic_energy`, `subscale_kinetic_energy` (VMS only),
`linear_momentum`, `angular_momentum`, `force`, `divergence_l2`, and, for cases with an exact
solution, `velocity_error_l2`, `velocity_gradient_error_l2` and `pressure_error_l2`. For cases
that declare `Problem.monitor_region`, `region_velocity_l2` and `region_subscale_l2` (VMS only)
are √∫_R|u_h|² and √∫_R|ũ|² over that region, at the same times as the kinetic energies. They
observe the zero-strip assumption and constrain nothing.
All integrals are global and use the common quadrature rule; error norms use twice its degree.
Momentum vectors expand into `_x`, `_y`, `_z` columns; angular momentum in 2D is one scalar.
A case declares `Problem.force_boundaries`, a tuple of boundary names; `force` writes
`force_<boundary>_<axis>` for each. It measures the fluid-on-body force
−∫σn ds with σ = −pI + ν(∇u + ∇uᵀ), without coefficient normalization. A force request without
declared boundaries is rejected. Forces on an individual open face depend on the pressure gauge;
closed-boundary forces do not.

Conservation diagnostics (spec: Conservation checks) are assembled from the forms the run uses,
at the converged solve-level state (`pressure_time`):

- `linear_momentum_defect` (`_x`, `_y`, `_z`) is D_k = c(a; u_h, e_k), and G(e_k) for
  EMAC–VMS.
- `angular_momentum_defect` is D_φ = c(a; u_h, φ) − (ũ, (a·∇)φ) with φ = e_k × x, and G(φ)
  for EMAC–VMS. It is a scalar in 2D, with φ = (−y, x).
- `convective_energy_defect` is 𝒩 = c(a; u_h, u_h) + (ũ, N − (a·∇)u_h), and G(u_h) + (ũ, N)
  for EMAC–VMS.
- `linear_momentum_defect_predicted` and `angular_momentum_defect_predicted` are the interior
  predictions of chapter Tables 5.2–5.3 for the same fields, with
  D_χ(w) = (∇·u_h, w) − χ_N(ũ, ∇w) and, for linear splitting, J = (ũ, (u_h·∇)φ). They are zero
  for EMAC–VMS.

Here c includes the selected outflow term, a = u_h + ũ* is the advector of the residual, and ũ is
the substituted subscale; Galerkin has ũ = 0. The defects include boundary fluxes: they equal the
chapter Tables 5.1–5.2 only with u_h = 0 on the whole boundary. On periodic domains the
linear-momentum and energy fluxes cancel but the angular-momentum fluxes do not.

`energy_budget` is the per-step Crank–Nicolson closure. It is accepted only for Crank–Nicolson
with dynamic VMS subscales. Its columns are:

- `rate`: (E^(n+1) − E^n)/δt, with E = ½‖u_h + ũ‖² for ASGS and ½(‖u_h‖² + ‖ũ‖²) for OSS at the
  endpoints.
- The half-level dissipation terms `viscous` (ν‖∇u_h‖² or 2ν‖∇ˢu_h‖²), `subscale`
  ((τ_m⁻¹ũ, ũ), with τ_m and not τ_t), `viscous_cross` (−2(ũ, L_ν u_h)) and `pressure_subscale`
  ((τ_c(∇·u_h − ζ_h), ∇·u_h)), and their sum `dissipation` 𝒟.
- `power`: (f, u_h + ũ) with the endpoint-average force.
- `defect_predicted`: 𝒩_pred from Tables 5.1–5.2.
- `defect_measured`: 𝒩_meas = power − rate − 𝒟.
- `residual`: r = rate + 𝒟 + 𝒩_pred − power.

With u_h = 0 on the whole boundary, or for the energy on periodic domains, r vanishes to solver
tolerance for every form.

`momentum_budget` is the per-step linear-momentum closure, with the same restrictions and blank
startup row. For each axis k it writes `rate_k` = (M_k^(n+1) − M_k^n)/δt with
M_k = (u_h + ũ, e_k) for ASGS and (u_h, e_k) for OSS at the endpoints, `defect_k` = D_k (as
`linear_momentum_defect`), `forcing_k` = (f, e_k) with the endpoint-average force, and
`residual_k` = rate + D_k − F_k. On periodic domains the residual vanishes to solver tolerance.
With Dirichlet walls it is the consistent wall reaction, the force on the fluid, and is only
reported. The terms are captured after the solve, before the histories advance.
The backward-Euler startup step is excluded and its columns are blank: it would add
‖Δ(u_h + ũ)‖²/(2δt) for ASGS, or (‖Δu_h‖² + ‖Δũ‖²)/(2δt) for OSS.

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

The `gresho` case is the piecewise steady Euler vortex of the spec section "Gresho vortex
(Phase 12b)" on (−0.5, 0.5)², default mesh 32×32, with f = 0 and four no-slip walls, all
declared for force measurements. `[case_parameters] viscosity` defaults to ν = 1e-6,
the authorized fallback for the nominally inviscid OSS families: at ν = 0, τ_c vanishes
where the fluid is at rest, and the frozen-τ OSS Jacobian then has identically
zero ζ rows. The viscous ASGS families set `viscosity = 1e-3`. At either positive viscosity,
velocity and pressure errors are measured against the Euler reference fields, which are not
exact solutions of the unforced viscous problem. Dynamic ASGS at a = 0 and ν ≈ 0 has
τ_m⁻¹ → 0: the subscale absorbs ∂ₜu_h and the resolved mass cancels. This accepted method
limitation is why the inviscid families use OSS only. A case may declare
`Problem.monitor_region`, an indicator (1 inside, 0 outside) evaluated at quadrature points;
Gresho declares the near-wall frame max(|x|, |y|) ≥ 0.45.

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

`--configs DIR` selects another config directory and `--runner simulation` runs each config
with `scripts/run_simulation.py`, on the same 4 ranks, instead of the validation run and
analysis. The log and the end marker are named after the directory: `<dir>_queue.log` and
`<DIR>_DONE` in `results/`. A simulation run counts as finished once its `run_info.json`
records an exit status, so a run that stopped without converging is kept as a recorded failure
and not repeated; a folder without one is treated as interrupted. `--dry-run` prints the order
and the action for each run without running or moving anything. `--smoke` needs the validation
runner.

## Conservation runs (Phase 12b)

The 55 configurations in `configs/chapter5/conservation/` are generated by
`scripts/generate_conservation_configs.py` and committed: 12 inviscid Gresho OSS runs,
12 viscous Gresho ASGS runs (ν = 1e-3), 5 long OSS runs, one refined OSS run, one refined
viscous ASGS run and 24 periodic Taylor–Green runs. They total 10,200 time steps; every
configuration sets `inner_max_iterations = 200`. Regeneration removes the 13 retired
near-inviscid ASGS configs. The owner launches the queue after the validation queue has
finished; it is never started by the agent. From the repository root:

```bash
# Regenerate the configurations (they must match the committed files):
python scripts/generate_conservation_configs.py
# Preview the order:
bash scripts/run_validation_queue.sh --configs configs/chapter5/conservation --runner simulation --dry-run
# Launch it in a detached tmux session that activates the Firedrake environment itself
# (tmux is not installed on the development machine yet: sudo apt install tmux):
tmux new-session -d -c "$PWD" -s conservation "bash -c 'source /home/ali/venv-firedrake/bin/activate && \
    bash scripts/run_validation_queue.sh --configs configs/chapter5/conservation \
    --runner simulation 2>&1 | tee -a results/conservation_queue.console'"
# Follow it, and detach again with Ctrl-b d:
tmux attach -t conservation
```

Runs go to `results/gresho/` and `results/taylor_green_2d/`, one line per run to
`results/conservation_queue.log`, and the summary to `results/CONSERVATION_DONE`.

```bash
# Tables and the long-Gresho plot, into results/conservation_analysis/ (reads the tree only):
python scripts/analyse_conservation.py
```

The analysis writes `<family>.csv` (every column) and `<family>.md` for `gresho`, `gresho_visc`,
`gresho_long`, `gresho_refined`, `gresho_visc_refined` and `taylor_green`, `summary.md` with the
runs that did not complete, and `gresho_long.png`. Each residual, defect and drift is reported in absolute terms
and divided by the fixed scales of its run from `initial_state.json` (spec: Conservation-run
normalization): U_c = √(2E(0)/|Ω|), L_c = |Ω|^(1/d), M_c = |Ω|U_c, A_c = M_c L_c,
F_c = |Ω|U_c²/L_c, torque F_c L_c and P_c = |Ω|U_c³/L_c. Measured defects are compared with
their predictions by |measured − predicted| over the same scale, never by dividing by the
prediction. Closures use steps 2 onward; drifts use the resolved `linear_momentum`,
`angular_momentum` and `kinetic_energy`. A run without `initial_state.json` gets blank scaled
columns. Ten-step prefixes on one core took about 0.5 s per OSS Gresho step,
0.61 s per viscous ASGS Gresho step after startup/compilation, and 1.1 s per Taylor–Green
step. The revised estimate is about 2.1–2.6 h serial for all 55 runs; the progress log gives
the measurements, startup allowance and unmeasured refined-mesh scaling assumption.

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
