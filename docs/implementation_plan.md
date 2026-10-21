# Implementation plan: VMS-stabilized incompressible Navier–Stokes in Firedrake

This plan turns `docs/formulation_spec.md` into a clean, phased code base. The spec is the source of truth for every equation; this plan fixes the code structure, the interfaces, the tests and the order of work.

## 1. How to use this plan

- Work one phase at a time, in order. A phase is done when its acceptance tests pass, the linter is clean and the progress log is updated. Then stop and report; do not start the next phase until told to.
- Implement exactly what the current phase asks. Do not add stubs, placeholders or options for later phases.
- If the spec is ambiguous, contradicts this plan, or seems wrong, stop and ask. Never resolve it silently in code.
- Equation references in this plan point to sections of the spec (for example "spec: Master formulation").
- Commit once per phase, with the message `Phase N: <title>`.

## 2. Scope

In scope:

- Galerkin baseline (Taylor–Hood, no stabilization).
- The eight VMS cases of the spec: ASGS or OSS × static or dynamic subscales × linear or nonlinear scale splitting.
- Time schemes BDF1, BDF2, BDF3 and Crank–Nicolson, for both scales.
- Inner linearization by Newton or Picard; outer fixed point for subscales and projections.
- Meshes from Firedrake's built-in constructors or from gmsh `.msh` files, with boundaries referred to by name.
- 2D and 3D. Nothing in the core may be 2D-specific.
- Logs of every nonlinear iteration and residual, per-step records, diagnostics and field output.
- The verification tests of section 10 and the validation scripts of phase 10.

Out of scope for this plan: convective forms other than skew-symmetric type 1 and the Chapter 5 conservation study, split OSS, adaptive time stepping, checkpoint and restart, iterative preconditioners and parallel performance tuning, turbulence benchmarks.

## 3. Ground rules for the code

- Python as required by the installed Firedrake. Type hints on all public functions. `ruff` for linting and formatting, line length 100.
- Plain functions and dataclasses. A class only where there is state that persists across time steps (subscale storage, time stepper, output writer). No inheritance hierarchies.
- Registries (plain dicts from name to function) only where variation is real: time schemes, convective forms, cases, diagnostics.
- Naming: the `formulation/` layer uses the spec's symbols (section 6) so it reads side by side with the spec. Everywhere else, names are descriptive.
- Docstrings: one line. Functions that implement a formula cite the spec section, e.g. `"""Velocity subscale ũ (spec: Master formulation)."""`.
- Comments only for a non-obvious reason why. No comments that restate the code, no commented-out code, no TODO or FIXME.
- No unused functions, parameters, imports or options. `ruff` must be clean.
- No `print`. Use the `logging` module (section 9). Only MPI rank 0 writes logs and CSV files.
- Global quantities (norms, energies, forces) are computed with `assemble`, never with NumPy reductions over local arrays, so everything stays correct in parallel.
- Configuration is passed explicitly. No global mutable state.
- Dimension independence: use `mesh.geometric_dimension()`; no hard-coded components, cross products or 2D formulas outside dimension-aware diagnostics.

## 4. Repository layout

```
navier-stokes-vms/
├── README.md
├── pyproject.toml
├── docs/
│   ├── formulation_spec.md            # agreed spec; never edited by the agent
│   ├── implementation_plan.md         # this file
│   ├── progress_log.md                # one entry per phase (section 12)
│   └── verification/VERIFICATION.md
├── src/navier_stokes_vms/
│   ├── __init__.py
│   ├── config.py
│   ├── problem/
│   │   ├── mesh_loading.py
│   │   ├── boundary_conditions.py
│   │   └── problem_definition.py
│   ├── discretization/
│   │   ├── function_spaces.py
│   │   ├── quadrature.py
│   │   └── time_schemes.py
│   ├── formulation/
│   │   ├── operators.py
│   │   ├── convective_forms.py
│   │   ├── stabilization_parameters.py
│   │   ├── galerkin.py
│   │   └── vms.py
│   ├── solver/
│   │   ├── subscale_storage.py
│   │   ├── oss_projections.py
│   │   ├── nonlinear_solver.py
│   │   └── time_stepper.py
│   ├── postprocessing/
│   │   ├── diagnostics.py
│   │   ├── run_logging.py
│   │   └── output.py
│   └── cases/
│       ├── __init__.py                # registry: case name -> build function
│       ├── manufactured_spatial.py
│       ├── manufactured_temporal.py
│       ├── kovasznay.py
│       ├── taylor_green_2d.py
│       ├── ethier_steinman_3d.py
│       ├── lid_driven_cavity.py
│       ├── cylinder_dfg_2d.py
│       └── lattice_vortex.py
├── meshes/
│   ├── cylinder_dfg_2d.geo
│   └── cylinder_dfg_2d.msh
├── configs/                           # one TOML file per experiment
├── scripts/
│   ├── run_simulation.py              # python scripts/run_simulation.py configs/<file>.toml
│   └── run_convergence_study.py
├── tests/
│   ├── data/                          # small meshes for tests
│   ├── unit/                          # seconds
│   └── verification/                  # minutes, marked @pytest.mark.slow
└── results/                           # git-ignored
```

## 5. Module responsibilities

| Module | Owns | Must not |
| --- | --- | --- |
| `config.py` | Option dataclasses, TOML loading, validation of illegal combinations, resolved-config dump | Import Firedrake |
| `problem/mesh_loading.py` | Built-in meshes and `.msh` loading; boundary name → tag tables | Know about equations |
| `problem/boundary_conditions.py` | `Dirichlet(value)`, `Outflow()` descriptions; conversion to Firedrake `DirichletBC` and `ds` measures by name | Know about time schemes |
| `problem/problem_definition.py` | `Problem` dataclass (section 8); pressure null-space decision | Build forms |
| `discretization/function_spaces.py` | Mixed space, unconstrained projection spaces, subscale quadrature spaces | Choose degrees (config does) |
| `discretization/quadrature.py` | The single quadrature rule: `dx`, `ds`, the matching quadrature element, positivity check | Anything else |
| `discretization/time_schemes.py` | α₀, history weights, solve level, CN extrapolation, startup sequence, data level | Import Firedrake |
| `formulation/operators.py` | L_a, L*_a, R_u, R̂_u, R_p, elementwise Laplacian | Store state |
| `formulation/convective_forms.py` | c(a; u, v) including the outflow term | Know about VMS |
| `formulation/stabilization_parameters.py` | τ_m, τ_c, τ_t | Anything else |
| `formulation/galerkin.py` | Galerkin residual form, parameterized by the advection velocity `a` | Know about solvers |
| `formulation/vms.py` | Subscale expression and VMS terms from the three switches | Know about solvers |
| `solver/subscale_storage.py` | ũ* and histories at quadrature points; local subscale update; CN extrapolation of ũ | Assemble global systems |
| `solver/oss_projections.py` | Weighted projections ξ_h, ζ_h | Know about time schemes beyond τ_t |
| `solver/nonlinear_solver.py` | Inner Newton or Picard solve, outer fixed point, convergence tests, iteration records | Know about specific VMS cases |
| `solver/time_stepper.py` | One time step as spec steps 1–6, the time loop, data and boundary updates at the solve level | Contain formulas |
| `postprocessing/diagnostics.py` | Registry of named functionals evaluated per step | Write files |
| `postprocessing/run_logging.py` | `run.log`, iteration and time-step CSV writers | Compute physics |
| `postprocessing/output.py` | Results folder, resolved config, run info, VTK and diagnostics CSV | Compute physics |

## 6. Spec symbols in code

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

The README repeats this table.

## 7. Configuration

One TOML file per run. Every field has a default except `case`. `config.py` converts it to frozen dataclasses, validates it and writes the fully resolved version into the results folder.

```toml
case = "cylinder_dfg_2d"
run_name = "bdf2_dyn_nl_oss"           # results/<case>/<run_name>/

[case_parameters]                      # free-form, read only by the case
reynolds = 100

[mesh]
source = "gmsh"                        # builtin | gmsh
file = "meshes/cylinder_dfg_2d.msh"    # gmsh only
# builtin only: shape = "rectangle" | "periodic_rectangle" | "box" | "periodic_box",
#               cells = [32, 32], lengths = [1.0, 1.0]

[discretization]
velocity_degree = 1
pressure_degree = 1

[time]
scheme = "bdf2"                        # bdf1 | bdf2 | bdf3 | crank_nicolson
dt = 0.01
t_end = 8.0
startup = "lower_order"                # lower_order | exact_history

[stabilization]
method = "vms"                         # none | vms
subscale_space = "oss"                 # asgs | oss
subscale_dynamics = "dynamic"          # static | dynamic
scale_splitting = "nonlinear"          # linear | nonlinear
c1 = 4.0
c2 = 2.0
pressure_subscale = true               # tau_c term on or off

[convection]
form = "skew_symmetric_1"

[solver]
linearization = "newton"               # newton | picard
newton_tau = "frozen"                  # frozen | exact
picard_iterations_before_newton = 0
inner_rtol = 1e-10
inner_atol = 1e-12
inner_max_iterations = 30
outer_rtol = 1e-8
outer_max_iterations = 30
subscale_update = "accelerated"        # fixed_point | accelerated
subscale_rtol = 1e-10
subscale_max_iterations = 50
petsc = { ksp_type = "preonly", pc_type = "lu", pc_factor_mat_solver_type = "mumps" }

[output]
vtk_every = 10                         # 0 disables field output
diagnostics = ["kinetic_energy", "divergence_l2"]
console_level = "INFO"
```

Validation rules (each violation raises `ValueError` with a message naming the fields):

1. `method = "none"` requires `velocity_degree = pressure_degree + 1` and no non-default `[stabilization]` field.
2. `method = "vms"` requires `velocity_degree = pressure_degree`.
3. `startup = "exact_history"` requires a case with an exact solution (checked when the problem is built).
4. `picard_iterations_before_newton > 0` requires `linearization = "newton"`.
5. Only the stabilization combinations implemented so far are accepted; each phase widens this set. An unimplemented combination is rejected with a message naming the phase that adds it.
6. Names in `diagnostics` must exist in the diagnostics registry.

`newton_tau = "frozen"` is the default because ‖a‖ is not differentiable at a = 0, which occurs for flows started from rest. This is the inexact Newton variant of the spec (spec: Inner linearization).

## 8. Problems, meshes and boundaries

A case module exposes one function, `build(config) -> Problem`. `Problem` holds data only:

| Field | Type | Meaning |
| --- | --- | --- |
| `mesh` | Firedrake mesh | From `mesh_loading` |
| `boundary_conditions` | `dict[str, Dirichlet \| Outflow]` | Keyed by boundary name |
| `boundary_tags` | `dict[str, int]` | Name → mesh tag |
| `viscosity` | `float` | ν |
| `forcing` | `(x, t) -> UFL` | f |
| `initial_velocity` | `(x) -> UFL` | u₀ |
| `exact_velocity`, `exact_pressure` | `(x, t) -> UFL` or `None` | For errors and exact-history startup |
| `reference_values` | `dict` or `None` | Validation targets |

- `x` is `SpatialCoordinate(mesh)`; `t` is a Firedrake `Constant` owned by the time stepper. `Dirichlet.value` is also `(x, t) -> UFL`.
- When a case has an exact solution, the forcing is derived from it symbolically with UFL (`f = ∂t u + u·∇u − νΔu + ∇p`, with time handled through `ufl.variable`), never typed by hand. Every case with an exact solution checks that its velocity is divergence-free.
- Boundary names for built-in meshes are `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max`. `mesh_loading.py` holds the name → tag table for each built-in constructor; test U7 confirms it against the installed Firedrake.
- For `.msh` files, names come from the `$PhysicalNames` section. The case refers only to names; a missing name raises an error listing the available ones.
- Periodic directions have no boundary names.
- The pressure null space (constant pressure) is used exactly when no `Outflow` boundary exists.

## 9. Outputs and logs

Every run writes to `results/<case>/<run_name>/`:

| File | Content |
| --- | --- |
| `config_resolved.toml` | Fully resolved configuration |
| `run_info.json` | Git commit and dirty flag, Firedrake and PETSc versions, Python version, MPI size, host, start and end time, wall time, exit status |
| `run.log` | Human-readable log: setup summary, one line per time step, warnings, errors. Also echoed to the console at `console_level` |
| `inner_iterations.csv` | One row per inner (Newton or Picard) iteration |
| `outer_iterations.csv` | One row per outer iteration |
| `time_steps.csv` | One row per time step |
| `diagnostics.csv` | One row per time step, one column per configured diagnostic |
| `fields/` | VTK output of `u_h`, `p_h` every `vtk_every` steps |

Columns:

- `inner_iterations.csv`: `step, time, outer_iteration, inner_iteration, method, residual_norm, relative_residual, update_norm, linear_iterations`
- `outer_iterations.csv`: `step, time, outer_iteration, inner_iterations, velocity_change, pressure_change, subscale_change, xi_change, zeta_change, fe_residual_norm, converged`
- `time_steps.csv`: `step, time, scheme_used, dt, outer_iterations, inner_iterations_total, converged, wall_time`

Definitions:

- `residual_norm` is the ℓ² norm of the assembled nonlinear residual vector of the master FE equation, Dirichlet rows excluded. For Newton it comes from the SNES monitor. For Picard it is assembled explicitly from the full nonlinear form after each linear solve, so Newton and Picard rows are comparable.
- `relative_residual` is `residual_norm` divided by its value at inner iteration 0 of the same inner solve.
- `scheme_used` records the startup sequence (for example `bdf1, bdf2, bdf3, bdf3, …`).
- Changes are relative norms between consecutive outer iterates, computed with `assemble` and the common quadrature rule.

On non-convergence (inner or outer), the run writes the last state to `fields/`, logs an error with the iteration history and exits with a non-zero status. It never continues silently.

## 10. Tests

Unit tests run in seconds and run in every phase. Verification tests are marked `slow` and run in the phases that list them. Rates are measured on the two finest meshes or time steps and must be at least the expected rate minus 0.2.

| ID | File | Checks |
| --- | --- | --- |
| U1 | `unit/test_time_schemes.py` | BDFk differentiates every polynomial of degree ≤ k exactly at t^(n+1); CN differentiates every quadratic exactly at t^(n+½) and the extrapolation reproduces linear functions; startup sequences for `lower_order` and `exact_history`; data level for CN is the endpoint average |
| U2 | `unit/test_quadrature.py` | On the reference triangle and tetrahedron, for q = 3, 6: all weights positive, all monomials of degree ≤ q integrated exactly; a non-polynomial function interpolated into the subscale quadrature space and integrated with `dx` equals the direct `dx` integral to round-off (same points, same weights) |
| U3 | `unit/test_convective_form.py` | c(a; v, v) = 0 to 1e-12 for random FE fields on a periodic mesh (2D and 3D, P1 and P2); with an outflow boundary, c(u, u, u) equals ½∫_{Γ_N}(u·n)\|u\|² |
| U4 | `unit/test_asgs_equals_supg.py` | Static linear ASGS residual assembled from `galerkin.py` + `vms.py` equals an independently written classical form, Galerkin + Σ_K τ_m(D_t u + a·∇u − νΔ_h u + ∇p − f, a·∇v + νΔ_h v + ∇q)_K + Σ_K τ_c(∇·u, ∇·v)_K, to relative 1e-12, for random fields, P1/P1 and P2/P2. The reference is written in the test from this formula, never by calling `vms.py`. This checks the adjoint-sign erratum |
| U5 | `unit/test_dissipation_sign.py` | For random (u, p) with f = 0 and no time terms: the stabilization contribution tested with (u, p) is ≥ 0 for P1; for P2, the viscous Galerkin term plus the stabilization contribution is ≥ 0. Static ASGS from phase 5, static OSS added in phase 8 |
| U6 | `unit/test_oss_projections.py` | Π_τ(w_h) = w_h for random FE w_h and variable τ > 0, including nonzero boundary values; after computing ξ_h with variable τ_t, the assembled vector (ũ, w_h) for all basis functions has relative norm ≤ 1e-10 |
| U7 | `unit/test_meshes_and_boundaries.py` | Built-in rectangle and box: each name maps to the boundary its coordinate says (checked by integrating x over `ds(tag)`); a small `.msh` in `tests/data/` with physical names loads with the right name table; a `DirichletBC` by name changes only that boundary's values |
| V1 | `verification/test_spatial_convergence.py` | Case `manufactured_spatial`: exact solution linear in time, so every time scheme is exact and only the spatial error remains. Unit square, ν = 0.01, meshes 8, 16, 32, two BDF1 steps from the exact initial state. Expected rates (‖u‖_L2, ‖∇u‖_L2, ‖p‖_L2): Taylor–Hood P2/P1 (3, 2, 2); VMS P1/P1 (2, 1, 1); VMS P2/P2 (3, 2, 2). Variants: Galerkin and all static VMS cases |
| V2 | `verification/test_temporal_convergence.py` | Case `manufactured_temporal`: u = g(t)(y, x), p = 0, which lies in every velocity space used, so only the temporal error remains. 4×4 mesh, T = 1, N = 20, 40, 80, 160, error ‖u − u_ex‖_L2 at T. Expected orders: BDF1 1, BDF2 2, BDF3 3 (`exact_history`), CN 2. Variants: Galerkin and all dynamic VMS cases. BDF3 with `lower_order` startup is logged and must show at least 1.8 |
| V3 | `verification/test_steady_equivalence.py` | Case `kovasznay` (Re = 40) marched with BDF1 to steady state (relative step change < 1e-10). Static and dynamic subscales reach the same velocity (relative difference ≤ 1e-8) for each implemented space and splitting; Newton and Picard reach the same velocity (≤ 1e-8); the error against the exact solution drops from mesh 16 to 32 |
| V4 | `verification/test_taylor_green.py` | Case `taylor_green_2d` on a periodic (0, 2π)² mesh, ν = 0.01, BDF2 to T = 1. The velocity error drops by at least a factor 3 from mesh 16 to 32 for P1/P1; the kinetic energy at T is within 1 % of the exact value on mesh 32. Variants: Galerkin and Dyn-Nl-OSS |

V1 covers the static cases only on purpose. Dynamic subscales start at zero and relax on the time scale τ, so a short run would not exercise them; their spatial accuracy follows from V3 (dynamic and static reach the same steady state), and V2 checks their time integration.

## 11. Firedrake notes

Check each of these against the installed Firedrake version; the API changes between releases.

- Quadrature consistency: build the subscale space from a `Quadrature` element of degree q with the same scheme that `dx` uses, and pass the degree explicitly to every measure (`dx(degree=q)` or metadata). τ is not polynomial, so automatic degree estimation must never be used. Test U2 guards this.
- Pointwise subscale updates are interpolations of UFL expressions into the quadrature space. The accelerated update solves (I + τ_t∇u)ũ = rhs pointwise; `ufl.inv` of a d×d matrix is acceptable.
- ‖a‖ = `sqrt(dot(a, a))`. With `newton_tau = "exact"` the derivative is infinite at a = 0; if a NaN appears, raise an error that names this option.
- Time: one `Constant` for the time at the solve level, plus one each for t^n and t^(n+1) where CN averaging needs them. Dirichlet values and forcing are UFL expressions in these constants, updated before each solve.
- Pressure null space: `MixedVectorSpaceBasis` with a constant `VectorSpaceBasis` on the pressure. With MUMPS a null-pivot option may also be needed; V3 and V4 confirm the treatment.
- Newton residual logging: attach a monitor to the solver's SNES through petsc4py; read linear iterations from the SNES. Picard: same solver object with `snes_type = "ksponly"` on the frozen-advector form, followed by an explicit residual assembly.
- `CellDiameter(mesh)` for h. The elementwise Laplacian is `div(grad(u))` inside cell integrals; it vanishes for affine P1.
- Projection spaces carry no `DirichletBC` and no null space.
- Logs and CSV files are written by rank 0 only; use `mesh.comm`.

## 12. Progress log

After each phase, append to `docs/progress_log.md`:

```
## Phase N: <title> (<date>)
Done: <one paragraph>
Tests: <commands run and results; rate tables for verification tests>
Deviations from spec or plan: <none, or each one with the reason and the question asked>
Open questions: <none, or list>
```

## 13. Phases

### Phase 0: Verify the spec addendum

- Goal: independent check of the spec sections added for this plan: Time schemes (BDF3 and Crank–Nicolson for both scales, data averaging, startup, stability notes), Galerkin baseline, and the revised notation and fixed choices.
- Checks: BDF3 coefficients exact for cubics; CN one-leg form equivalent to Codina et al. 2007 eqs. 22–23 with θ = ½ and their definition of f^(n+θ); τ_t for CN; orthogonality of the extrapolated OSS subscale; the order claim for BDF3 with lower-order startup; the Dirichlet averaging; the Galerkin baseline against the master formulation with all stabilization removed.
- Output: a new pass at the top of `docs/verification/VERIFICATION.md` in the same table format as before. No code.

### Phase 1: Skeleton, configuration and run infrastructure

- Files: `pyproject.toml`, `README.md` (install, run, symbol table), `src/navier_stokes_vms/__init__.py`, `config.py`, `postprocessing/run_logging.py`, `postprocessing/output.py`, `scripts/run_simulation.py`, `tests/unit/test_config.py`.
- Details: `run_simulation.py` loads and validates a config, creates the results folder, writes `config_resolved.toml`, `run_info.json` and `run.log`, then exits, because nothing can be solved yet. Validation accepts no stabilization combination yet except `method = "none"`.
- Acceptance: valid configs load; each validation rule of section 7 has a test that triggers it; the run script produces the three files.
- Out of scope: meshes, forms, solvers.

### Phase 2: Problems, meshes and boundaries

- Files: `problem/mesh_loading.py`, `problem/boundary_conditions.py`, `problem/problem_definition.py`, `cases/__init__.py`, `cases/manufactured_spatial.py`, `cases/manufactured_temporal.py`, `cases/kovasznay.py`, `cases/taylor_green_2d.py`, `tests/data/square_named.msh` (generated with gmsh; if gmsh is unavailable, stop and ask), `tests/unit/test_meshes_and_boundaries.py`.
- Details: section 8. `manufactured_spatial`: divergence-free velocity from the stream function ψ = sin²(πx) sin²(πy), times (1 + t); pressure (1 + t) cos(πx) sin(πy); Dirichlet on all four sides. `manufactured_temporal`: u = g(t)(y, x) with g(t) = cos t, p = 0, Dirichlet on all sides. `kovasznay`: Re = 40 on [−0.5, 1] × [−0.5, 1.5], exact Dirichlet data. `taylor_green_2d`: periodic (0, 2π)², u = (−cos x sin y, sin x cos y) e^(−2νt), p = −¼(cos 2x + cos 2y) e^(−4νt).
- Acceptance: U7; each case builds on a coarse mesh; derived forcing for `taylor_green_2d` is zero to round-off (it solves Navier–Stokes exactly).

### Phase 3: Quadrature and time schemes

- Files: `discretization/quadrature.py`, `discretization/time_schemes.py`, `tests/unit/test_quadrature.py`, `tests/unit/test_time_schemes.py`.
- Details: q = 3 × velocity degree. `time_schemes.py` is pure Python: for a scheme name and step index it returns the scheme used at that step (startup), α₀, the history weights, the solve level (1 or ½), the CN extrapolation weights, and the weights that give data at the solve level (spec: Time schemes).
- Acceptance: U1, U2.

### Phase 4: Galerkin baseline

- Spec: Galerkin baseline, Operators and parameters (convective form), Boundary conditions and quadrature, Time schemes, Inner linearization.
- Files: `discretization/function_spaces.py`, `formulation/convective_forms.py`, `formulation/galerkin.py`, `solver/nonlinear_solver.py`, `solver/time_stepper.py`, `postprocessing/diagnostics.py` (velocity and pressure errors, divergence, kinetic energy), `postprocessing/output.py` (VTK, diagnostics CSV), `scripts/run_convergence_study.py`, tests U3, V1, V2, V3, V4 restricted to Galerkin.
- Details: `galerkin.py` builds the residual with `a` as an argument. Newton passes the unknown velocity; Picard passes a frozen `Function` and solves with `ksponly` (section 11). The time stepper implements startup, history shifting, CN half-level solve and extrapolation, data at the solve level, and the logs of section 9. Pressure errors use zero-mean normalization when the null space is active, and compare at t^(n+½) for CN.
- Acceptance: U3; V1, V2, V3 (Newton = Picard part), V4 for Galerkin. `run_simulation.py` runs a Galerkin config end to end and writes every file of section 9.

### Phase 5: Static linear ASGS

- Spec: Operators and parameters, Master formulation, The eight cases (row Sta-Lin-ASGS).
- Files: `formulation/operators.py`, `formulation/stabilization_parameters.py`, `formulation/vms.py`, tests U4, U5 (ASGS), V1 and V3 for this case.
- Details: `vms.py` implements the master formulation as written, driven by the switches. Only the static linear ASGS branch is reachable in this phase; validation rejects the rest.
- Acceptance: U4, U5, V1 (P1/P1 and P2/P2), V3 (Newton = Picard).

### Phase 6: Dynamic subscales

- Spec: Master formulation (χ_D terms), Time schemes, Time-step algorithm step 6.
- Files: `solver/subscale_storage.py`, extensions to `vms.py` and `time_stepper.py`, tests V2 and V3 for Dyn-Lin-ASGS.
- Details: subscale histories at quadrature points for every scheme (three levels for BDF3); stored after convergence; CN extrapolation of ũ; startup sequence identical to the velocity.
- Acceptance: V2 for Dyn-Lin-ASGS (all four schemes, BDF3 with `exact_history`); V3 static = dynamic.

### Phase 7: Nonlinear scale splitting

- Spec: The eight cases (Nl rows), Time-step algorithm steps 3 and 5.
- Files: extensions to `subscale_storage.py` (iterate ũ*, fixed-point and accelerated updates), `nonlinear_solver.py` (outer loop and `outer_iterations.csv`), `vms.py` (a = u + ũ*).
- Acceptance: V1 for Sta-Nl-ASGS; V2 for Dyn-Nl-ASGS; V3 static = dynamic and Newton = Picard; the fixed-point and accelerated subscale updates reach the same converged state (relative difference ≤ 1e-8, checked inside V3).

### Phase 8: Orthogonal subscales

- Spec: OSS projections, Master formulation (χ_A = 0), The eight cases (OSS rows), Time-step algorithm step 4.
- Files: `solver/oss_projections.py`, extensions to `function_spaces.py` (unconstrained spaces), `vms.py` (reduced residual, ξ, ζ) and `nonlinear_solver.py` (projection changes in the outer loop); tests U5 (OSS), U6.
- Acceptance: U5, U6; V1 for Sta-Lin-OSS and Sta-Nl-OSS; V2 for Dyn-Lin-OSS and Dyn-Nl-OSS; V3 for all OSS pairs; V4 for Dyn-Nl-OSS. All eight cases are now accepted by validation.

### Phase 9: Diagnostics and 3D check

- Files: `postprocessing/diagnostics.py` extensions, `cases/ethier_steinman_3d.py`.
- Details: diagnostics registry entries for kinetic energy of u_h and of ũ, linear momentum, angular momentum (a scalar in 2D, a vector in 3D), divergence norm, and drag and lift on a named boundary from the traction σn with σ = −pI + ν(∇u + ∇uᵀ). `ethier_steinman_3d`: the Ethier–Steinman exact solution on the cube [−1, 1]³, forcing derived as in section 8 (it must come out as zero to round-off, which checks the formula).
- Acceptance: 3D smoke run on meshes 4 and 8 for Galerkin and Dyn-Nl-OSS: runs to completion, and the velocity error drops with refinement; the same run on 2 MPI ranks gives the same diagnostics to 1e-10 relative.

### Phase 10: Validation scripts

- Files: `cases/lid_driven_cavity.py`, `cases/cylinder_dfg_2d.py`, `cases/lattice_vortex.py`, `meshes/cylinder_dfg_2d.geo` and `.msh`, configs for each.
- Lid-driven cavity: Re = 100 and 1000 against the centreline data of Ghia, Ghia and Shin (1982), embedded in the case with the citation. If the table is not available to you, stop and ask.
- Cylinder, DFG 2D-2 (Schäfer and Turek, 1996): channel [0, 2.2] × [0, 0.41], cylinder centre (0.2, 0.2), radius 0.05, inflow 4U_m y(0.41 − y)/0.41² with U_m = 1.5, ν = 10⁻³ (Re = 100), do-nothing outflow at x = 2.2, physical names `inlet`, `outlet`, `walls`, `cylinder`. Report c_D,max, c_L,max and the Strouhal number with c_D = 2F_D/(Ū²D), Ū = 1, D = 0.1, against the reference intervals c_D,max ∈ [3.22, 3.24], c_L,max ∈ [0.99, 1.01], St ∈ [0.295, 0.305]. Monitor u·n on the outlet (spec: Boundary conditions and quadrature).
- Lattice vortex: the thesis setup (domain (0, 1)², ν = 7·10⁻⁵, strong Dirichlet data, the thesis' exact solution), reporting kinetic energy against the exact decay.
- Acceptance: each script runs, writes all outputs and a short comparison table into the progress log. These are validation runs, not tests; a miss is reported, not hidden by tuning.

## 14. After this plan

The Chapter 5 work adds convective forms to `convective_forms.py`, the validation rule that rejects forms needing ∇·ũ under nonlinear splitting, and energy and momentum budget tests. It will get its own spec section and plan.
