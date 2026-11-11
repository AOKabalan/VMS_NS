# Usage and reference

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


## Additional isothermal experiments

The separate [exact Gresho benchmark](gresho_exact.md) uses zero viscosity,
normal-only free-slip walls, and Crank–Nicolson from step one with
`time.startup = "none"`. The `gresho_degree_2` case uses the same exact problem
with separate presets and output paths. These cases retain the original
viscous `gresho` benchmark. `total_kinetic_energy` measures the endpoint norm
of resolved plus subscale velocity, including the cross term, and requires
dynamic VMS.

The [unstratified KH experiment](../experiments/kelvin_helmholtz_2020/README.md)
uses the `kelvin_helmholtz` case, periodic horizontal boundaries and
`FreeSlip(normal_component)` walls. Its experiment drivers include optional
profiling and linear-solver configurations. It is distinct from the coupled
`kelvin_helmholtz_shear_layer` case below.

The optional `emac_vms_full` convection form implements the full-stress
closure described in the [experiment findings](../experiments/emac_vms_advantage/FINDINGS.md)
and formulation specification. Like `emac_vms`, it requires VMS, stress
viscosity and Newton without Picard warm-up. Its reported pressure adds only
the resolved kinetic correction; the additional subscale pressure correction
is omitted. Neither EMAC–VMS variant is supported for thermal runs.

## Thermal Boussinesq baseline (B1)

Run the coupled manufactured problem with the existing entry point:

```bash
python scripts/run_simulation.py configs/boussinesq/manufactured_galerkin.toml
```

B1 uses Taylor–Hood Galerkin velocity/pressure, continuous temperature of the
velocity degree, BDF1 with lower-order startup, Newton without Picard warm-up,
and Laplacian viscosity. Momentum convection may be `skew_symmetric_1` or
`nonconservative`; `[thermal].convection` is `skew_symmetric` (default) or
`nonconservative`. An omitted temperature degree resolves to the velocity degree.
Static and dynamic linear thermal VMS and the remaining time schemes are
described below. Nonlinear thermal splitting, stress viscosity, other momentum
forms and open thermal boundaries remain later phase work.

A thermal case attaches `ThermalProblem` to `Problem`. It supplies positive
thermal diffusivity, nonnegative expansion coefficient, a gravity vector pointing
downwards, constant reference temperature, source, initial temperature and an
independent boundary map. Reuse `Dirichlet` for prescribed temperature and
`DiffusiveFlux` for the outward flux `-kappa grad(T) dot n`; `DiffusiveFlux()`
means insulation. A volumetric heating rate must be divided by `rho_0 c_p`
before supplying the temperature source. Boundary coverage is checked, and
B1 requires zero prescribed normal velocity. `[thermal]` alone does not enable
physics for an isothermal case.

Thermal runs add `fields/temperature.pvd`, `thermal_metadata.json`, and a
`temperature_change` outer-iteration column. Request `temperature_l2`,
`temperature_error_l2` (needs an exact temperature), `thermal_content`, or
`buoyancy_power` in the diagnostic list. BDF1 writes all physical fields at the
endpoint. Failure snapshots include the available temperature without advancing
its successful history. Isothermal runs retain their existing output files and
columns. The complete equations and later gates are in
[boussinesq_formulation.md](boussinesq_formulation.md) and
[boussinesq_implementation_plan.md](boussinesq_implementation_plan.md).

Thermal Newton runs check the momentum, continuity and temperature residuals
separately, in addition to the full mixed-vector norm. Set
`thermal.momentum_residual_scale`, `thermal.continuity_residual_scale` and
`thermal.temperature_residual_scale` to positive characteristic equation scales
when working in different physical units; all default to 1 for the supplied
nondimensional example. Each normalized block must meet its own initial-relative
and absolute tolerance. These are algebraic load-vector norms, not physical
energy norms. Their values appear in `inner_iterations.csv` as
`velocity_residual`, `pressure_residual` and `temperature_residual`, and the
chosen scales are recorded in thermal metadata.


## Static thermal ASGS/OSS and passive scalars (B2)

Run the stabilized manufactured examples through the same entry point:

```bash
python scripts/run_simulation.py configs/boussinesq/manufactured_static_asgs.toml
python scripts/run_simulation.py configs/boussinesq/manufactured_static_oss.toml
```

These use equal-order P1/P1/P1; P2/P2/P2 is also verified. Set
`stabilization.subscale_dynamics = "static"` and `scale_splitting = "linear"`.
B2 introduced these examples with BDF1; B3 extends their time schemes below.
Newton, the supported convection choices and the B1 boundary restrictions apply. Setting `pressure_subscale = false` removes the pressure subscale
and its OSS projection without disabling scalar stabilization. Frozen-τ Newton
freezes the mechanical and scalar parameters consistently. Exact-τ Newton is
available but retains the zero-advector derivative limitation described for the
flow solver.

A case builder can attach a tuple of `PassiveScalarProblem` data to
`Problem.passive_scalars`. Each entry supplies its name, positive diffusivity,
source, initial value, independent boundary map and optional exact solution.
Names must be unique identifiers and must avoid physical/projection field names
and the `eta_` prefix. Temperature supplies the buoyancy; passive fields remain
one-way coupled within the same monolithic system. Their degree, convection
and stabilization switches inherit the thermal choices. B2 also supports these
passive fields in the Galerkin baseline.

For a passive field `q`, OSS unknowns are ordered as
`velocity, pressure, temperature, q, xi, eta_temperature, eta_q, (zeta)`.
Only physical fields receive Dirichlet conditions; the pressure constant is the
sole gauge. Scalar endpoints/history storage, VTK series (`q.pvd`), metadata,
inner residual columns and outer change columns are named independently. The
`passive_scalar_content` and `passive_scalar_l2` diagnostic groups produce
columns such as `passive_scalar_content_q`. They measure resolved scalar values;
conservation statements still require the formulation's boundary, source and
convection qualifications.

Thermal projection residuals use the corresponding thermal/mechanical reference
scales and have independent initial-relative/absolute checks. Empty passive data
preserve temperature-only fields/output, and isothermal configurations retain
their original spaces and output columns. B4 adds nonlinear thermal splitting
below; G4 adds free-slip/channel support below, and B7 supplies the KH application.


## Dynamic thermal subscales and time schemes (B3)

Linear thermal ASGS/OSS now accepts `subscale_dynamics = "dynamic"` and all four
existing schemes: `bdf1`, `bdf2`, `bdf3`, `crank_nicolson`. Galerkin and static
thermal cases accept those time schemes too. Examples:

```bash
python scripts/run_simulation.py configs/boussinesq/manufactured_dynamic_oss_bdf2.toml
python scripts/run_simulation.py configs/boussinesq/manufactured_dynamic_asgs_crank_nicolson.toml
python scripts/run_simulation.py configs/boussinesq/manufactured_dynamic_oss_bdf3.toml
```

Each active scalar has a separate dynamic subscale history on the same
quadrature rule as the velocity subscale. The integrator's actual history
weights, startup sequence and effective parameters apply to both scales.
Crank–Nicolson starts with BDF1, then solves at the half level using endpoint
averages of sources and prescribed boundary data. Resolved scalar and dynamic
subscale endpoints are extrapolated; pressure and projections stay at the solve
level. Temperature/passive VTK and scalar diagnostics use endpoint states.
`buoyancy_power` uses both fields at the solve level, whose time is the existing
`pressure_time` column; the diagnostic is not a kinetic/potential budget.

`startup = "exact_history"` requires exact velocity and exact solutions for
**every** active scalar, including passive tracers. Their negative-time FE
histories are interpolated from those callbacks. Subscale history slots remain
zero. For forced subscale ODEs, zero negative history is generally inconsistent
and can limit startup accuracy; selecting exact resolved histories does not
promise third-order accuracy for arbitrary forced BDF3 runs. Lower-order BDF3
startup retains the existing global second-order limitation. The source-free
order checks isolate the time operator with compatible zero OSS subscales.

Failed steps preserve all resolved and subscale endpoint histories and write
available physical fields at the attempted solve time. The legacy `energy_budget`
and `momentum_budget` groups retain their isothermal meanings. B5 adds a distinct
thermal energy group below; thermal momentum budgets remain B6 work.


## Nonlinear thermal splitting and Picard (B4)

All eight ASGS/OSS × static/dynamic × linear/nonlinear thermal variants now
share the same solver, including passive scalars. For nonlinear scale splitting,
set `stabilization.scale_splitting = "nonlinear"`. The total advector
$u_h+\widetilde u^*$ enters both momentum and every scalar transport equation,
their residuals and stabilization parameters. Scalar subscales remain explicit
algebraic expressions; only the velocity subscale has an outer iterate.

```bash
python scripts/run_simulation.py configs/boussinesq/manufactured_nonlinear_oss.toml
python scripts/run_simulation.py configs/boussinesq/manufactured_nonlinear_asgs_picard.toml
```

Use `solver.linearization = "newton"` or `"picard"`; Newton also accepts
`picard_iterations_before_newton`. Picard freezes the shared total advector and
all stabilization weights while keeping temperature-driven buoyancy coupled
implicitly. Both methods check the true coupled residual. The existing
`newton_tau` choices apply to all transport/projection parameters. Local
velocity updates support `subscale_update = "accelerated"` or `"fixed_point"`.
Static nonlinear updates retain the known coarse-mesh convergence limitation;
a failed local or coupled solve is reported without advancing histories.

For nonlinear thermal runs, `outer_iterations.csv` additionally records
`velocity_subscale_residual`, `<scalar>_subscale_change` (including temperature),
and `<field>_residual` for every physical/projection equation after the local
update. The velocity closure residual is the relative $L^2$ difference between
the algebraic velocity subscale and its stored iterate. Scalar changes compare
their algebraic values before and after the outer solve/update. All these
changes and the velocity closure must meet `outer_rtol`; every scaled FE block
must also pass its inner tolerance after the update. Scalar closures are exact
by elimination and verified independently against the formulation.

The thermal convection choices remain `skew_symmetric` and `nonconservative`;
conservative thermal transport and coupled conservation budgets belong to B5.
These manufactured examples do not start the cavity or Kelvin–Helmholtz studies.


## Thermal budgets and natural convection (B5)

Set `thermal.convection = "conservative"` to use the conservative scalar form.
Momentum convection is selected independently; its thermal support remains
`skew_symmetric_1` and `nonconservative` with Laplacian viscosity. The scalar
default is still `skew_symmetric`.

The conservative Galerkin example includes explicit MUMPS null-pivot detection
for the constant-pressure null space:

```bash
python scripts/run_simulation.py configs/boussinesq/manufactured_galerkin_conservative.toml
```

Its `[solver.petsc]` sets `mat_mumps_icntl_24 = 1`. This handles a singular
factorization without pressure regularization or a formulation change; it does
not change global solver defaults. The pressure null space and spatial-mean
normalization still apply. See the
[PETSc MUMPS option definitions](https://petsc.org/main/src/mat/impls/aij/mpi/mumps/impl/imumps.c.html).

Select these groups in `output.diagnostics` as needed:

| Group | Output and requirements |
| --- | --- |
| `temperature_extrema` | Global FE nodal min/max at the endpoint, without clipping |
| `temperature_gradient_error_l2` | Endpoint gradient error; exact temperature required |
| `thermal_subscale_l2` | Solve-level subscale norm for every scalar; zero for Galerkin |
| `thermal_content_budget` | Resolved/total scalar content, history rate, source/flux and splitting/continuity/iterate terms; no scalar Dirichlet walls |
| `thermal_variance_budget` | Rate, diffusion, subscale damping/cross terms, convection, projection/source/flux and Dirichlet reaction work |
| `boussinesq_energy_budget` | Kinetic/potential storage and rates, mechanical terms, source/boundary work and each exchange remainder |
| `thermal_boundary_flux` | Integrated molecular-gradient flux and variational Dirichlet reaction for each scalar/wall; reaction is blank on natural walls |
| `cavity_nusselt` | Volume, hot/cold gradient and hot/cold reaction estimates; heated-cavity reference scales required |

For example, `thermal_variance_budget_temperature_diffusion_cross` and
`thermal_content_budget_q_total_content` identify their scalar explicitly.
`total_content` includes the dynamic ASGS scalar moment. Budget/flux/subscale
groups use the solve time in `pressure_time`, even when temperature field output
and ordinary temperature diagnostics use the CN endpoint. All budget values are
captured before histories advance.

The rates use the actual BDF/CN history operators. Dynamic ASGS includes cross
storage; static ASGS retains mixed rates. Residuals subtract all named source,
boundary, exchange and finite-solver corrections. Their normalized values divide
by the sum of absolute work/rate terms, with floor `1e-12`, recorded in
`thermal_metadata.json`. A finite difference of energy CSV values is a different
quantity for BDF schemes. The optional G5 stabilization-work ratios are not provided.

Potential energy requires an admissible, single-valued `g·x`: gravity-direction
periodicity and thermal Dirichlet walls are rejected. This budget currently
requires homogeneous full-vector Dirichlet, free-slip or periodic velocity.
Content checks require an
admissible constant scalar test. Matching scalar/pressure degrees are needed to
interpret the measured continuity work as a vanishing FE residual. Nonzero
Dirichlet temperature is supported in the variance budget through reaction work;
it is not silently treated as insulation.

The new closed-cavity case is available through the same entry point:

```bash
python scripts/run_simulation.py configs/boussinesq/cavity_ra1000.toml
python scripts/run_simulation.py configs/boussinesq/cavity_ra10000.toml
```

It uses width `L=1`, configurable height `aspect_ratio`, hot/cold temperatures
`+0.5/-0.5`, insulated horizontal walls, and no-slip velocity. Free-fall units give
`viscosity=sqrt(Pr/Ra)`, `diffusivity=1/sqrt(Ra*Pr)`, unit buoyancy coupling and
`gravity=(0,-1)`. `initial_profile` is `conduction` or `isothermal`; the latter
applies the wall-temperature jump on the first solve. Derived coefficients and
reference scales are recorded. Positive hot-wall Nu is minus outward flux;
positive cold-wall Nu is outward flux. The reaction includes all weak scalar
terms and can differ from a molecular-gradient estimate on a finite mesh.

The three `codina_tall_ra*_smoke.toml` files exercise aspect ratio eight and the
2010 Rayleigh numbers for 20 steps. Their 16×128 triangular meshes, same-scheme
BDF2 subscales and short duration differ from the paper. They are integration
checks, not predictions of oscillation, chaos or long-time benchmark agreement.
The [B5 report](verification/boussinesq_b5.md) contains primary-source conventions,
refinement results, MPI/3D checks and reproducible commands.


## Horizontally periodic free-slip channels (G4)

Case builders can use `periodic_channel(cells, lower, upper)` for a triangular
rectangle periodic only in x. It returns `y_min=3` and `y_max=4`, the bottom and
top walls after translation. `periodic_rectangle()` remains fully periodic.
Construct any coordinate remapping before creating the spaces and BCs.

Attach `FreeSlip(normal_component=1)` to both channel walls. Only the normal
velocity DOFs are constrained; tangent velocity has natural zero tangential
traction. Component indices are zero-based. Flat walls must be normal to the
declared coordinate axis; invalid/tilted wall descriptions are rejected. The
same description works on axis-aligned 3D faces with the corresponding index.

Temperature and passive scalars use `DiffusiveFlux()` on these walls for
insulation. Their OSS projections remain unconstrained. Pressure retains its
constant null space, and there is no zero-mean velocity constraint. Initial
velocity/history fields are clamped only in the slip-normal trace; report any
remaining interpolation divergence when setting up a physical case.

Two runnable preservation examples use the existing entry point:

```bash
python scripts/run_simulation.py configs/boussinesq/free_slip_uniform.toml
python scripts/run_simulation.py configs/boussinesq/free_slip_uniform_stress.toml
```

`free_slip_channel` defaults to a 2×1 channel centered vertically, horizontal
speed 1, T=T_ref=1 and passive tracer q=2, with zero sources/fluxes. Its
`case_parameters` accept finite `speed`, positive `length`/`height`, and boolean
`thermal`. The stress example sets `thermal=false`; thermal stress viscosity
remains B6 work. The uniform examples use explicit MUMPS null-pivot detection
for the pressure null space. Their scalar values are verification data; B7 will
supply stratified shear/moisture initial conditions.

The coupled kinetic–potential budget accepts free-slip wall tests and permits
nonzero tangential wall speed. Its other admissibility checks remain active.
The [G4 report](verification/free_slip_g4.md) records component/gauge tests,
uniform preservation, nonuniform natural-traction decay, failure handling,
serial/MPI checks and compatibility results. This completes the boundary gate
for B7; it does not start the KH growth or cloud studies.

## Kelvin–Helmholtz coupled smoke test (B7.1)

```bash
python scripts/run_simulation.py configs/boussinesq/kh_smoke_static_asgs.toml
python scripts/run_simulation.py configs/boussinesq/kh_smoke_static_oss.toml
```

These run three BDF1 steps on a uniform 24×24 triangular channel mesh with
static linear subscales, temperature and monolithic passive moisture `q`.
They exercise initialization and coupling; this coarse mesh is not a validated
growth-rate or cloud configuration. The next step is B7.2 startup/resolution work.

`kelvin_helmholtz_shear_layer` accepts dimensionless `reynolds`, `prandtl`,
`schmidt`, `richardson`, `n_waves`, `wave_number`, `height`, perturbation
`epsilon`, `epsilon_noise`, `random_seed`, `noise_modes`, wall `wall_buffer`,
`wall_transition_width`, and moisture `moisture` (boolean), `rh_peak`,
`moisture_width`, `saturation_height_rate`, `saturation_temperature_rate`.
Use `[mesh].cells=[nx,nz]` for actual counts. Defaults and the exact profile
definitions are in [the reference conventions](../data/references/kelvin_helmholtz/README.md).
Unknown parameters are rejected; clustering is not exposed at B7.1.

Horizontal velocity and temperature equal tanh in the core and flatten smoothly
near insulated free-slip walls. The compact perturbation and moisture cutoff
also satisfy the wall conditions. Physical pressure starts from the analytic
hydrostatic primitive, interpolated and mean normalized. There is no holding
force. This initializer is optional `Problem.initial_pressure` data; it supplies
physical pressure, including when the solver uses a modified pressure unknown.
Omitting it preserves the previous zero-pressure startup and output schema.

`thermal_metadata.json` includes actual physical/scalar coefficients, case
parameters, cutoff coordinates and generated Fourier phases.
`initial_state.json` additionally records pressure mean, hydrostatic interpolation
defect, divergence, vertical energy, normal-flow and scalar wall-flux norms.
Opt-in `vertical_kinetic_energy` measures the final velocity component;
`advective_cfl` samples `k*dt*|u|/h_min` at the common quadrature points and
takes an MPI-global maximum each step. This is an endpoint sample rather than
a continuous supremum. VTK names the temperature and `q` explicitly; condensate
diagnostics and animations are reserved for B7.5.

The [B7.1 report](verification/kelvin_helmholtz_b71.md) records the coarse startup
motion and verification limits. A nonzero discrete pressure–buoyancy defect can
generate vertical motion much larger than the tiny seed; do not infer instability
growth from this smoke test. Scalar content under linear splitting and skew
transport has the already documented G3 defect; the reported budget includes it.

## Kelvin–Helmholtz mesh and startup checks (B7.2)

```bash
python scripts/run_simulation.py configs/boussinesq/kh_startup_control.toml
python scripts/run_simulation.py configs/boussinesq/kh_startup_seeded.toml
python scripts/analyse_kh_startup.py \
  results/kelvin_helmholtz_shear_layer/b72_startup_control \
  results/kelvin_helmholtz_shear_layer/b72_startup_seeded \
  --output /tmp/kh_startup.json --figure /tmp/kh_startup.png
```

The two base configs compare epsilon=0 and 1e-4 at Re=1000, Pr=1, Ri=0.1 with
dynamic nonlinear OSS, BDF2, dt=0.05 and t_end=0.5. The `uniform_medium`,
`uniform_fine`, `half_dt`, `p2` and `clustered` config pairs under the same
`kh_startup_` prefix reproduce the spatial/time/degree/clustering checks.
All integrations use `run_simulation.py`; the new analysis script reads
completed output only and rejects mismatched or incomplete pairs.

`case_parameters.clustering=0` retains the uniform mesh. A positive value
applies the sinh map to vertex heights, leaving affine triangles. Collapsed
rows and unresolved wall plateaus are rejected. `thermal_metadata.json` records
the actual row heights and global map-Jacobian, cell-measure, edge and aspect
ranges. The shortest mapped edge still controls stabilization and measured CFL.

`output.horizontal_profiles_every` is an optional nonnegative integer (default
0). Positive intervals write `horizontal_profiles.csv` and quadrature metadata.
Each row names `velocity_x`, `velocity_y`, `temperature` and any passive scalars
at a physical height and endpoint time. Output includes initialization, interval
steps and the final/early-stop state. Profiles integrate FE traces on mesh rows
and midrows using collective Gauss-point sampling on horizontal half-cells;
they are horizontal means, not centreline samples. Unsupported cases are
rejected before running. Disabled profile output adds no file and leaves the
old resolved-config output schema unchanged.

The analysis compares each mean to its own initial mean and retains both raw
vertical-energy histories. Their difference is not perturbation energy.
The [B7.2 report](verification/kelvin_helmholtz_b72.md) shows that zero-seed
motion still dominates the tested seeded runs. It establishes no clean fitting
window; further resolution and mode/control checks must precede a B7.3 rate fit.

## Kelvin–Helmholtz Fourier growth measurements (B7.3)

Add `horizontal_fourier_modes = [1, 2]` under `[output]` with a positive
`horizontal_profiles_every` to record sine/cosine coefficients of every named
field alongside its mean. Indices are distinct positive integers, with physical
wavenumber `2*pi*m/domain_length`. The option defaults to an empty tuple and
is omitted from existing resolved configs when unused. Single-wave pilots
measure index 1 and its second harmonic. Horizontal weights are overintegrated;
metadata records their normalization and quadrature.

```bash
python scripts/run_simulation.py configs/boussinesq/kh_growth_pilot_control.toml
python scripts/run_simulation.py configs/boussinesq/kh_growth_pilot_seeded.toml
python scripts/analyse_kh_growth.py \
  results/kelvin_helmholtz_shear_layer/b73_pilot_control \
  results/kelvin_helmholtz_shear_layer/b73_pilot_seeded \
  --window 10 20 --output /tmp/kh_growth.json --figure /tmp/kh_growth.png
```

The analysis reads completed matched runs only. Omitting `--window` exports
energies without fitting. The explicit window must contain at least six positive
energy samples. Results include total-energy and Fourier-energy fits separately,
fit residuals, half-window slopes, shape overlap, control contamination and
recorded base-profile drift. The seeded-minus-control mode energy is computed
from coefficient differences, never by subtracting two energies. Check the
reported second harmonic and perturbation amplitude before interpreting a fit
as linear growth. Its OLS error is not a mesh-convergence certificate.

`kh_growth_coarse_*` and `kh_growth_pilot_*` provide the bounded 32×48/48×72
comparison at Ri=0.1, Re=1000, dt=0.05, to t=20. They are preliminary resolutions.
`kh_growth_reference_timing.toml` measures three steps on 128×192.
The 36 configs under `configs/boussinesq/kh_growth/` prepare the full
Ri=0/0.1/0.3, Re=1000/2000, reference/refined/half-dt matched-pair matrix.
Preparing these inputs does not establish that the corresponding runs passed.

For the bounded t=30 follow-up, `kh_growth_late_*` retains 48×72 and dt=0.05,
`kh_growth_late_half_dt_*` uses dt=0.025, and `kh_growth_late_finer_*` uses
64×96 and dt=0.05. Each has a control/seeded pair and identical physical
sampling intervals. The
[refinement report](verification/kelvin_helmholtz_b73_refinement.md) records
the primary 20–30 and predeclared secondary windows, independent diffusion
reference, costs and reproduction commands. Their local modal rates are
window-conditioned; the full B7.3 matrix remains pending.
