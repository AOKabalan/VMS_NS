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
- Inner linearization by Newton or Picard, with the OSS projections ξ_h, ζ_h as unknowns of the inner solve; outer fixed point for the nonlinear-splitting subscale ũ*.
- Meshes from Firedrake's built-in constructors or from gmsh `.msh` files, with boundaries referred to by name.
- 2D and 3D. Nothing in the core may be 2D-specific.
- Logs of every nonlinear iteration and residual, per-step records, diagnostics and field output.
- The verification tests of section 10 and the validation scripts of phase 10.

Phase 12 adds the Chapter 5 convective forms and conservation study described below. Otherwise out of scope: split OSS, adaptive time stepping, checkpoint and restart, iterative preconditioners and parallel performance tuning, turbulence benchmarks beyond the optional Phase 12c comparison.

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
- Dimension independence: use `mesh.geometric_dimension` (a property in the installed Firedrake); no hard-coded components, cross products or 2D formulas outside dimension-aware diagnostics.

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
│   │   ├── oss_projections.py
│   │   └── vms.py
│   ├── solver/
│   │   ├── subscale_storage.py
│   │   ├── reusable_assembly.py
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
│   ├── cylinder_dfg_2d.geo            # coarse test mesh and production levels 1–3
│   └── cylinder_dfg_2d.msh            # level 1; _level2.msh, _level3.msh, _coarse.msh
├── data/references/<case>/            # validation reference data, each with a README.md
├── configs/                           # one TOML file per experiment
│   └── validation/                    # one TOML file per validation run (phase 10)
├── scripts/
│   ├── run_simulation.py              # python scripts/run_simulation.py configs/<file>.toml
│   ├── run_convergence_study.py
│   ├── run_validation.py              # a validation run and its comparison
│   └── run_validation_queue.sh        # every validation run in sequence, resumable
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
| `discretization/function_spaces.py` | Mixed space, the OSS space V × Q × V̂ (× Q̂) whose projection components stay unconstrained, subscale quadrature spaces | Choose degrees (config does) |
| `discretization/quadrature.py` | The single quadrature rule: `dx`, `ds`, the matching quadrature element, positivity check | Anything else |
| `discretization/time_schemes.py` | α₀, history weights, solve level, CN extrapolation, startup sequence, data level | Import Firedrake |
| `formulation/operators.py` | L_a, L*_a, R_u, R̂_u, R_p, elementwise Laplacian | Store state |
| `formulation/convective_forms.py` | Resolved convection, convective subscale coupling, residual convection and outflow term for all six Chapter 5 choices | Store state or know about solvers |
| `formulation/stabilization_parameters.py` | τ_m, τ_c, τ_t | Anything else |
| `formulation/galerkin.py` | Galerkin residual form, parameterized by the advection velocity `a` | Know about solvers |
| `formulation/vms.py` | Subscale expression and VMS terms from the three switches | Know about solvers |
| `formulation/oss_projections.py` | The τ-weighted projection equations for ξ_h, ζ_h | Know about solvers |
| `solver/subscale_storage.py` | ũ* and histories at quadrature points; local subscale update; CN extrapolation of ũ | Assemble global systems |
| `solver/reusable_assembly.py` | Assemblers preprocessed once for forms evaluated many times | Know about equations |
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
| L_a, L*_a | `momentum_operator`; L*_a inside `vms_terms` | Operators; the convective part of L*_a is `convective_subscale_coupling` |
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

[mesh]                                 # every field optional; omitted means "use the case's default"
source = "gmsh"                        # "builtin" | "gmsh"; the case may support one or both
file = "meshes/cylinder_dfg_2d.msh"    # "<path>.msh"; gmsh only
# cells = [nx, ny(, nz)]               # builtin only

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
| `force_boundaries` | `tuple[str, ...]`, default empty | Named boundaries measured separately by the force diagnostic |
| `monitor_region` | `(x) -> UFL` or `None` | Indicator of an optional monitor region, 1 inside and 0 outside (Phase 12b) |

- `x` is `SpatialCoordinate(mesh)`; `t` is a Firedrake `Constant` owned by the time stepper. `Dirichlet.value` is also `(x, t) -> UFL`.
- When a case has an exact solution, the forcing is derived from it symbolically with UFL (`f = ∂t u + u·∇u − νΔu + ∇p`, with time handled through `ufl.variable`), never typed by hand. Every case with an exact solution checks that its velocity is divergence-free.
- Boundary names for built-in meshes are `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max`. `mesh_loading.py` holds the name → tag table for each built-in constructor; test U7 confirms it against the installed Firedrake.
- For `.msh` files, names come from the `$PhysicalNames` section. The case refers only to names; a missing name raises an error listing the available ones.
- Periodic directions have no boundary names.
- The case owns the geometry (shape, origin, lengths, periodicity) and its default mesh; the config only overrides resolution or file. A case raises a clear error for a mesh source it does not support.
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
| `diagnostics.csv` | One row per time step; scalar diagnostics and expanded vector/boundary components, with their field times |
| `initial_state.json` | \|Ω\|, d, and the kinetic energy and linear and angular momentum of the interpolated initial velocity, for the Phase 12b run scales |
| `fields/` | VTK output of `u_h`, `p_h` every `vtk_every` steps |

Columns:

- `diagnostics.csv`: `step, time, pressure_time, subscale_time`, then the configured scalar/component columns. `time` is the velocity endpoint; pressure error and force use `pressure_time` and the solve-level velocity. Dynamic subscale energy uses its stored endpoint; static energy uses the solve level (`subscale_time`, blank for Galerkin). Vectors use `_x`, `_y`, `_z`; force uses `force_<boundary>_<axis>`; angular momentum in 2D is scalar.
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
| U2 | `unit/test_quadrature.py` | On the reference triangle and tetrahedron, for q = 3, 6: all weights positive, all monomials of degree ≤ q integrated exactly; a non-polynomial function interpolated into the subscale quadrature space and integrated with `dx` equals the direct `dx` integral to round-off (same points, same weights). The same holds for the production vector factory `subscale_space(mesh, q)` with q = 3k, k = 1, 2, in 2D and 3D, sampling a non-polynomial vector expression and integrated with the production `dx` |
| U3 | `unit/test_convective_form.py` | c(a; v, v) = 0 to 1e-12 for random FE fields on a periodic mesh (2D and 3D, P1 and P2); with an outflow boundary, c(u, u, u) equals ½∫_{Γ_N}(u·n)\|u\|² |
| U4 | `unit/test_asgs_equals_supg.py` | Static linear ASGS residual assembled from `galerkin.py` + `vms.py` equals an independently written classical form, Galerkin + Σ_K τ_m(D_t u + a·∇u − νΔ_h u + ∇p − f, a·∇v + νΔ_h v + ∇q)_K + Σ_K τ_c(∇·u, ∇·v)_K, to relative 1e-12, for random fields, P1/P1 and P2/P2. The reference is written in the test from this formula, never by calling `vms.py`. This checks the adjoint-sign erratum. Dynamic ASGS: the time stepper's own residual at a BDF2 step and at a Crank–Nicolson step (P1/P1 and P2/P2, linear and nonlinear splitting, random velocity and subscale levels) equals a classical reference whose u_hist, ũ_hist, τ_t and (D_t ũ, v) use coefficients written as literals in the test (BDF2: 3/2; 2, −1/2. Crank–Nicolson: 2; 2gⁿ, with averaged data), never taken from `time_schemes.py` |
| U5 | `unit/test_dissipation_sign.py` | Freeze a random advection field a and assemble E(u,p; v,q) = ν(∇u,∇v) + Σ_K (−τ_m L_a(u,p), L*_a(v,q))_K + Σ_K (τ_c ∇·u, ∇·v)_K, with f = 0 and no time terms. Restrict to unconstrained dofs after homogeneous velocity Dirichlet conditions; the symmetric part must have λ_min ≥ −1e-12 λ_max. Check P1/P1 and P2/P2 on 4×4 squares with diagonals `left` and `crossed` and a 2×2×2 box, at viscosities in the viscous and convective regimes. Retain the P2 cell-diameter counterexample as a regression that must detect negative energy. Static ASGS from phase 5, static OSS added in phase 8, with ξ_h and ζ_h eliminated by the Schur complement of the coupled bilinear form |
| U6 | `unit/test_oss_projections.py` | Π_τ(w_h) = w_h for random FE w_h and variable τ > 0, including nonzero boundary values; after the coupled solve, the assembled vector (ũ, w_h) for all basis functions w_h of V̂_h has norm ≤ 1e-10 times that of (τ_t X, w_h), for the subscale expression and the stored ũ, static and dynamic subscales, linear and nonlinear splitting, P1/P1 and P2/P2, at a BDF2 step with non-zero ũ_hist |
| U7 | `unit/test_meshes_and_boundaries.py` | Built-in rectangle and box: each name maps to the boundary its coordinate says (checked by integrating x over `ds(tag)`); a small `.msh` in `tests/data/` with physical names loads with the right name table; a `DirichletBC` by name changes only that boundary's values |
| U8 | `unit/test_subscale_time_integration.py` | The time stepper's dynamic subscale path (step coefficients, histories, τ_t, storage, Crank–Nicolson extrapolation) with prescribed FE fields u_h = U and p_h, so that τ is constant and R = g(t)F(x): ũ at T against the exact y(T)F of ∂t ũ + ũ/τ = gF, ũ(0) = 0, for N = 10, 20, 40, 80, rates within 0.2 of: BDF1 1, BDF2 2, BDF3 with `lower_order` 2, Crank–Nicolson 2. BDF3 with `exact_history` is first order for generic data, because the stored subscale is zero at negative times (spec: Initial subscale) while the smooth extension of the solution is not, an O(δt) defect in the first step's history; it is third order when y vanishes to third order at t = 0 |
| U9 | `unit/test_reusable_assembly.py` | `reusable_assembler` and the classes built on it give results identical to `fd.assemble` and `Function.interpolate`, for Dyn-Nl-ASGS and Dyn-Nl-OSS P2/P2: the time stepper's coupled residual, every outer L² change (ξ_h, ζ_h included), the local subscale update and the stored subscale, with all inputs redrawn between checks. It guards the Firedrake internals the helper uses |
| V1 | `verification/test_spatial_convergence.py` | Case `manufactured_spatial`: exact solution linear in time, so every time scheme is exact and only the spatial error remains. Unit square, ν = 0.01, meshes 8, 16, 32, two BDF1 steps from the exact initial state. Expected rates (‖u‖_L2, ‖∇u‖_L2, ‖p‖_L2): Taylor–Hood P2/P1 (3, 2, 2); VMS P1/P1 (2, 1, 1); VMS P2/P2 (3, 2, 2). Variants: Galerkin and all static VMS cases. Exception: Sta-Nl-ASGS and Sta-Nl-OSS P1/P1 use meshes 16, 32, 64, with rates on 32 → 64, because static nonlinear subscales can be ill-posed on coarse meshes and their local update does not converge on mesh 8 (spec: The eight cases, the paragraph after the table) |
| V2 | `verification/test_temporal_convergence.py` | Case `manufactured_temporal`: u = g(t)(y, −x), p = 0, which lies in every velocity space used, so only the temporal error remains. The field must have nonzero vorticity; a gradient field's time error is absorbed by the pressure and leaves the velocity exact. 4×4 mesh, T = 1, N = 20, 40, 80, 160, error ‖u − u_ex‖_L2 at T, and the pressure error against the exact pressure at the time it is computed (t^(n+1) for BDF, t^(n+½) for Crank–Nicolson), with zero-mean normalization; both must show the expected orders. Expected orders: BDF1 1, BDF2 2, BDF3 3 (`exact_history`), CN 2. Variants: Galerkin and all dynamic VMS cases. BDF3 with `lower_order` startup is logged and must show at least 1.8 |
| V3 | `verification/test_steady_equivalence.py` | Case `kovasznay` (Re = 40) marched with BDF1 to steady state: the relative change between the last two time levels of the velocity, and of the stored subscale for dynamic cases, is below 1e-10. That every step's outer iteration converged is checked separately. Static and dynamic subscales reach the same velocity (relative difference ≤ 1e-8) for each implemented space and splitting; Newton and Picard reach the same velocity (≤ 1e-8); the error against the exact solution drops from mesh 16 to 32 |
| V4 | `verification/test_taylor_green.py` | Case `taylor_green_2d` on a periodic (0, 2π)² mesh, ν = 0.01, BDF2 with δt = 0.025 to T = 1. Galerkin P2/P1: the velocity error drops by at least a factor 3 from mesh 16 to 32, and the kinetic energy at T is within 1 % of the exact value on mesh 32. Dyn-Nl-OSS P1/P1: the velocity error drops by at least a factor 3 from mesh 16 to 32; the kinetic energy at T is within 1 % of the exact value on mesh 64, and its relative error drops by at least a factor 3 from mesh 32 to 64. P1/P1 is not accurate enough for the 1 % bound on mesh 32: every P1/P1 VMS case is 1.7–2.2 % low there (phase 8), with an error that decreases at second order |

V1 covers the static cases only on purpose. Dynamic subscales start at zero and relax on the time scale τ, so a short run would not exercise them; their spatial accuracy follows from V3 (dynamic and static reach the same steady state), and V2 checks their time integration.

## 11. Firedrake notes

Check each of these against the installed Firedrake version; the API changes between releases.

- Quadrature consistency: build the subscale space from a `Quadrature` element of degree q with the same scheme that `dx` uses, and pass the degree explicitly to every measure (`dx(degree=q)` or metadata). τ is not polynomial, so automatic degree estimation must never be used. Test U2 guards this.
- Pointwise subscale updates are interpolations of UFL expressions into the quadrature space. The accelerated update solves (I + τ_t∇u)ũ = rhs pointwise; `ufl.inv` of a d×d matrix is acceptable.
- `fd.assemble` and `Function.interpolate` preprocess the UFL form on every call, which dominated small runs with many outer iterations. Forms evaluated repeatedly (the residual norm, the outer changes, the local update and the stored subscale) use `reusable_assembler`, the same preprocessing done once (phase 8; results unchanged to the last digit).
- ‖a‖ = `sqrt(dot(a, a))`. With `newton_tau = "exact"` the derivative is infinite at a = 0; if a NaN appears, raise an error that names this option.
- Time: one `Constant` for the time at the solve level, plus one each for t^n and t^(n+1) where CN averaging needs them. Dirichlet values and forcing are UFL expressions in these constants, updated before each solve.
- Pressure null space: `MixedVectorSpaceBasis` with a constant `VectorSpaceBasis` on the pressure. With MUMPS a null-pivot option may also be needed; V3 and V4 confirm the treatment. The solver fixes a zero coefficient mean, not the zero spatial mean of the spec, so after each step's solve the time stepper subtracts ∫p dx/|Ω| (phase 7b, R4).
- Newton residual logging: attach a monitor to the solver's SNES through petsc4py; read linear iterations from the SNES. Picard: same solver object with `snes_type = "ksponly"` on the frozen-advector form, followed by an explicit residual assembly.
- `ufl.MinCellEdgeLength(mesh)` for h, the shortest edge of the cell. The elementwise Laplacian is `div(grad(u))` inside cell integrals; it vanishes for affine P1.
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
- Details: section 8. `manufactured_spatial`: divergence-free velocity from the stream function ψ = sin²(πx) sin²(πy), times (1 + t); pressure (1 + t) cos(πx) sin(πy); Dirichlet on all four sides. `manufactured_temporal`: u = g(t)(y, −x) with g(t) = cos t, a rotation with nonzero vorticity, p = 0, Dirichlet on all sides. `kovasznay`: Re = 40 on [−0.5, 1] × [−0.5, 1.5], exact Dirichlet data. `taylor_green_2d`: periodic (0, 2π)², u = (−cos x sin y, sin x cos y) e^(−2νt), p = −¼(cos 2x + cos 2y) e^(−4νt).
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

### Phase 7b: Review fixes

- Source: the independent review of phases 0–7, `docs/reviews/review_phases_0-7.md`, findings R1–R5.
- Files: `nonlinear_solver.py` (R1), `time_stepper.py` (R4), tests U2, U4, U8 (new), V3 and the inner-solver tests in `test_time_stepper.py`.
- R1: a failed Picard linear solve fails the inner solve at once, through the existing failure path; the Picard warm-up counts toward `inner_max_iterations`. R2: U4 for dynamic ASGS and the new U8; the review's two subscale mutations (backward Euler for the subscale under BDF2, no Crank–Nicolson extrapolation of ũ) must each fail at least one test. R3: V3 declares steady state between time levels. R4: with the pressure null space, the pressure has zero spatial mean. R5: U2 through the production `subscale_space`.
- Acceptance: all unit and verification tests pass, and each review mutation fails at least one test.

### Phase 8: Orthogonal subscales

- Spec: OSS projections, Master formulation (χ_A = 0), The eight cases (OSS rows), Time-step algorithm steps 2 and 4, Inner linearization.
- ξ_h and ζ_h (ζ_h only when τ_c is active) are unknowns of the inner solve, with their projection equations; the pressure null space acts on p_h only and no Dirichlet condition applies to ξ_h, ζ_h. Linear OSS cases have no outer loop; nonlinear OSS cases iterate on ũ* only, and the ξ_h, ζ_h changes enter their convergence test.
- Files: `formulation/oss_projections.py` (projection equations), extensions to `function_spaces.py` (the OSS space), `vms.py` (reduced residual, ξ, ζ), `time_stepper.py` (the coupled forms) and `nonlinear_solver.py` (ξ_h, ζ_h changes in the outer records); tests U5 (OSS), U6.
- Acceptance: U5, U6; V1 for Sta-Lin-OSS and Sta-Nl-OSS; V2 for Dyn-Lin-OSS and Dyn-Nl-OSS; V3 for all OSS pairs; V4 for Dyn-Nl-OSS. All eight cases are now accepted by validation.

### Phase 9: Diagnostics and 3D check

- Files: `postprocessing/diagnostics.py` extensions, `cases/ethier_steinman_3d.py`.
- Details: diagnostics registry entries for kinetic energy of u_h and of ũ, linear momentum, angular momentum (a scalar in 2D, a vector in 3D), divergence norm, and raw fluid-on-body force −∫σn ds on each boundary in `Problem.force_boundaries`, with σ = −pI + ν(∇u + ∇uᵀ). Drag/lift coefficients belong to Phase 10 (owner's Phase 9 clarification). Pressure-dependent diagnostics use the solve-level pressure and velocity, including the Crank–Nicolson half level. `ethier_steinman_3d`: the Ethier–Steinman exact solution on the cube [−1, 1]³, forcing derived as in section 8 (it must come out as zero to round-off, which checks the formula).
- Acceptance: 3D smoke run on meshes 4 and 8 for Galerkin and Dyn-Nl-OSS: runs to completion, and the velocity error drops with refinement. Compare every diagnostic at every step with the same run on 2 MPI ranks. For Galerkin, polynomial-integrand diagnostics (energies, momentum, angular momentum, divergence, forces) agree to 1e-10 relative; error norms against exact solutions have non-polynomial integrands and use 1e-8 relative (absolute floor 1e-13 for both). For VMS, each force component uses absolute tolerance 1e-3 times the largest absolute force component in the serial reference run, over all steps and measured faces. Every other VMS diagnostic uses 1e-4 relative with absolute floor 1e-4 in the nondimensional benchmark units. Keep default quadrature; non-polynomial 3D integrands have partition-dependent quadrature error (spec: Boundary conditions and quadrature). Report the actual largest differences and, without asserting it, whether the largest force difference decreases from mesh 4 to mesh 8. Galerkin forces verify the force diagnostic's MPI correctness; VMS pressure carries quadrature-orientation noise through τ, which coarse-mesh boundary tractions amplify. These bounds distinguish that error from O(1) MPI errors.
- For problems with the pressure null space, check the interpolated Dirichlet trace at every solve level: reject |∫u_D·n ds| > 1e-10 ∫|u_D·n| ds with a discrete-flux compatibility error. Do not apply the check with an outflow boundary.

### Phase 10: Validation scripts

Phase 10a delivers everything below except the long runs: reference data in `data/references/<case>/` (Ghia Tables I–II at Re = 100 and 1000, the FEATFLOW DFG 2D-2 results and archives, the thesis' exact energy decay), the three cases, the cylinder meshes (a coarse test mesh and production levels 1–3), `scripts/run_validation.py` (run or `--analyse-only`, writing `validation.csv` and `validation.md`), one config per run in `configs/validation/`, the resumable `scripts/run_validation_queue.sh` with a `--smoke` mode, the U5 check on the coarse cylinder mesh and a runtime estimate. The owner launches the full queue; its results and the comparison tables belong to phase 10b.

- Files: `cases/lid_driven_cavity.py`, `cases/cylinder_dfg_2d.py`, `cases/lattice_vortex.py`, `meshes/cylinder_dfg_2d.geo` and `.msh`, configs for each.
- Lid-driven cavity: Re = 100 and 1000 against the centreline data of Ghia, Ghia and Shin (1982), embedded in the case with the citation. If the table is not available to you, stop and ask.
- Cylinder, DFG 2D-2 (Schäfer and Turek, 1996): channel [0, 2.2] × [0, 0.41], cylinder centre (0.2, 0.2), radius 0.05, inflow 4U_m y(0.41 − y)/0.41² with U_m = 1.5, ν = 10⁻³ (Re = 100), do-nothing outflow at x = 2.2, physical names `inlet`, `outlet`, `walls`, `cylinder`. Report c_D,max, c_L,max and the Strouhal number with c_D = 2F_D/(Ū²D), Ū = 1, D = 0.1, against the reference intervals c_D,max ∈ [3.22, 3.24], c_L,max ∈ [0.99, 1.01], St ∈ [0.295, 0.305]. Monitor u·n on the outlet (spec: Boundary conditions and quadrature). Run the U5 eigenvalue check on a coarse version of the gmsh cylinder mesh, since the P2 stability condition depends on cell shape.
- Lattice vortex: the thesis setup (domain (0, 1)², ν = 7·10⁻⁵, strong Dirichlet data, the thesis' exact solution), reporting kinetic energy against the exact decay.
- Acceptance: each script runs, writes all outputs and a short comparison table into the progress log. These are validation runs, not tests; a miss is reported, not hidden by tuning.

## Phase 12 — Chapter 5

The formulation spec governs this phase. Phase 12a used the Phase 11b spec at commit 1fa820d; Phase 12b updates it against the 23-page chapter PDF with SHA-256 c4995958c5e410591261f676aae661f5854def033feb534dc930e40668491784. docs/chapter5_renumbering.md records the reference changes. No chapter edits belong here. The divergence pairing makes skew-symmetric 2 and EMAC available with nonlinear splitting. The former restriction on those forms is superseded.

### 12a: implementation and algebraic tests

Four stages, each tested and committed before the next:

1. **Forms/operators.** One convection module supplies resolved convection, convective subscale coupling, residual convection and the outflow term. Keep the existing default entry point. Implement the five forms and minimal EMAC–VMS, using the U10 index convention and the weak divergence pairing. Add the Laplacian/stress viscous operators and their adjoints. Test both defect tables, both splittings, arbitrary quadrature subscales, T1–T4/T6, boundary conditions and stress self-adjointness.
2. **Solver wiring.** Pass the selected forms through one assembly path, including the residual, OSS projection and master local update. Add the six choices to [convection].form and [viscosity].form = "laplacian" or "stress", preserving defaults. EMAC–VMS requires VMS and stress viscosity. EMAC and EMAC–VMS require Newton with zero Picard warm-up, per the explicit owner decision. Assemble every stress slot from the substituted expression and differentiate through it; only the residual advector and τ retain the stored subscale. Check Newton Jacobians against finite differences, including projection rows, and default regression against commit 1fa820d on three unit-scale runs in temporary directories.
3. **Physical pressure.** Keep the unknown/null space on P for EMAC variants. Recover p = P + |u|²/2 in a space containing the correction, normalize it separately, and use it for every output, error and force. Midpoint output uses the solve velocity, not the endpoint velocity or averaged endpoint squares. Test both time levels, output space and forces.
4. **Conservation budgets.** Preserve existing diagnostic meanings. Add the actual assembled convective defects and the dynamic midpoint closure r, with explicit previous/current endpoint and solve-level state. Test T5 for all forms from step 2 onward, the weighted pressure and velocity terms, the viscous cross term, and the startup exclusion.

Run ruff check, ruff format --check and the unit suite at each stage. Use at most two pytest workers while the validation queue runs (two is also the default for this work). Tests and regression outputs go under /tmp; do not read or write the repository's results/. Run no validation simulations in 12a. The slow verification suite is eligible only once the queue is finished, and its execution must still respect the unit-scale restriction for this phase. A budget stop is allowed only after a tested, committed stage; report the remaining stages.

### 12b: conservation runs

**12b-0, documents.** Add the Gresho reference, momentum budget, predicted momentum defects and region norms to the spec; independently verify the equations, boundaries, zero-viscosity limits and chapter references. Commit only after ruff check, ruff format check and the full unit suite pass. The Phase 12b-0 verification is a gate: resolve any equation, iteration or diagnostic-definition disagreement before implementation. Its two open items were resolved by owner decisions K-D1 (finite frozen-τ acceptance at rest) and K-D2 (analysis normalization), recorded in the spec and below.

**12b-1, implementation and unit tests.** Commit each numbered step separately with ruff check, ruff format check and the full unit suite passing. Assume the owner's validation queue is active even when ps cannot see it. Use OMP_NUM_THREADS=1 and at most two pytest workers; every test/probe output goes under /tmp. Do not read or write the repository's results/.

1. Add cases/gresho.py with piecewise UFL Euler reference fields, f=0, four no-slip/force walls and the frame max(|x|,|y|)≥0.45. Check divergence, Euler residual, pressure continuity and derived forcing. **K-D1, updated by owner Q1:** require finite residual, frozen-τ Jacobian, solve and every diagnostic through a complete tiny-mesh step for both OSS splittings at ν=10⁻⁶ and both ASGS splittings at ν=10⁻³. Viscosity is a Gresho case parameter, default 10⁻⁶ (the authorized fallback for the singular OSS system at ν=0), not a second module constant for the viscous family. At positive viscosity these Euler fields are reference fields, not exact viscous solutions; say so in the docstring and README. At a=0 and ν≈0, τ_m⁻¹→0: the dynamic ASGS subscale absorbs ∂ₜu_h and the resolved mass cancels. Accept this as a method limitation, not a solver problem; retain the strict nonlinear ASGS xfail at ν≤10⁻⁴. Keep the exact-τ failure check at rest and exact-Jacobian tests only at nonzero advection speed. No regularization of |a|.
2. Add momentum_budget with the energy_budget restrictions and explicit endpoint/solve-level state; add predicted D_k and D_φ; add both optional-region norms and their validation. Check periodic midpoint momentum closure from step 2 for all 24 form/space/splitting combinations, all predicted defects, a finite Gresho reaction and known region norms; extend the step-1 frozen-τ Gresho check to the new diagnostics. Gresho does not enforce a zero strip.
3. Generate and commit the exact 55-file matrix below with a small script; test file membership, configuration validity, unique run names, family viscosities and step totals. Regeneration removes the 13 retired near-inviscid ASGS configurations. Set inner_max_iterations=200 for every Phase 12b configuration (owner Q2). The ordinary OSS Gresho, long-Gresho and Taylor–Green ten-step timings are recorded in the progress log; owner Q1 also authorizes one ten-step prefix of gresho_visc_emac_vms_asgs_nonlinear if its timing is unknown. Use one core and /tmp. Record measured timings and the revised extrapolation, including the refined mesh assumption. These prefixes are the only matrix-related simulations authorized for the agent.
4. Generalize scripts/run_validation_queue.sh to accept a config directory and runner while preserving its default behavior, completion detection, interrupted-folder handling and one-line run logging. Document a tmux command in the README. Test with a temporary dry run.
5. Add scripts/analyse_conservation.py and test it on a synthetic tree in /tmp. Write per-family CSV/Markdown tables for maximum normalized energy/momentum closure from step 2, measured versus predicted defects, M/L/E drifts, final region norms and wall forces, iterations and failures. Plot long-Gresho energy, both momenta and L² velocity error, one curve per form. **K-D2 (owner decision):** normalize by the fixed characteristic scales of each run in spec: Conservation-run normalization, never by a value that can be zero; report every closure residual, defect and drift both absolute and scaled, and measured against predicted as |measured − predicted|/scale next to the scaled prediction. The scales need E(0), M(0), A(0) and |Ω| of the interpolated initial state, which the run must record.

Common settings for all 55 configurations:

| Setting | Value |
| --- | --- |
| Elements | P1/P1 |
| Viscous form | stress |
| Pressure subscale | on; τ_c-off cases excluded |
| Velocity subscales | dynamic |
| Time scheme | Crank–Nicolson with BDF1 startup |
| Linearization | Newton, no Picard warm-up |
| newton_tau | frozen |
| outer_rtol / outer_max_iterations | 1e-11 / 100 |
| inner_rtol | 1e-12 |
| inner_max_iterations | 200 (owner Q2) |
| subscale_rtol | 1e-13 |

In the matrix, {form} ranges in order over nonconservative, conservative, skew_symmetric_1, skew_symmetric_2, emac, emac_vms; {space} over asgs, oss; and {split} over linear, nonlinear. All files are in configs/chapter5/conservation/, named after their unique run names with the .toml suffix.

| Configuration family | Matrix / count | Mesh and viscosity | Time / steps per run | Diagnostics |
| --- | --- | --- | --- | --- |
| gresho_{form}_oss_{split} | 6 × 2 = 12 | 32×32; ν=10⁻⁶, the documented inviscid fallback | dt=0.005, T=0.5 / 100 | Gresho set below |
| gresho_visc_{form}_asgs_{split} | 6 × 2 = 12 | 32×32; ν=10⁻³ | dt=0.005, T=0.5 / 100 | Same Gresho set |
| gresho_long_{form} | nonconservative, conservative, skew_symmetric_1, emac, emac_vms; OSS nonlinear / 5 | 32×32; ν=10⁻⁶ | dt=0.01, T=10 / 1000 | Gresho set; Charnyi et al. 2017 Fig. 2 comparison |
| gresho_refined_emac_vms_oss | emac_vms, OSS, nonlinear / 1 | 64×64; ν=10⁻⁶ | dt=0.0025, T=0.5 / 200 | Gresho set |
| gresho_visc_refined_emac_vms_asgs | emac_vms, ASGS, nonlinear / 1 | 64×64; ν=10⁻³ | dt=0.0025, T=0.5 / 200 | Gresho set |
| taylor_green_{form}_{space}_{split} | 6 × 2 × 2 = 24 | Existing fully periodic taylor_green_2d, 32×32, ν=0.01 | dt=0.005, T=0.5 / 100 | Taylor–Green set below |

The Gresho diagnostic set includes resolved/subscale kinetic energy, resolved linear/angular momentum, energy_budget, momentum_budget, measured/predicted linear/angular/energy defects, force on all walls, both region norms, velocity/physical-pressure error and divergence. The Taylor–Green set includes resolved/subscale kinetic energy, linear momentum, both budgets, measured/predicted linear/energy defects and velocity error. No angular quantity is requested on the periodic domain.

The revised matrix still contains 55 runs and 10,200 time steps: 1200 ordinary OSS Gresho, 1200 viscous ASGS Gresho, 5000 long OSS Gresho, 200 refined OSS Gresho, 200 refined viscous ASGS Gresho and 2400 Taylor–Green. The 13 viscous ASGS runs replace the 13 near-inviscid ASGS runs. The authorized viscous ASGS prefix completed ten steps in 29.4 s including startup/compilation; steps 3–10 averaged 0.61 s. Combined with the existing OSS/Taylor–Green timings, 4–5× per-step cost on refined meshes and 3–25 s startup per run, the revised estimate is 2.1–2.6 h serial on one core. The progress log records the breakdown; refined costs and other forms are extrapolations, and there is no measured four-rank estimate.

Expected outcomes are reported, never tuned, with the K-D2 scales: |r|≤1e-9 P_c for every converged run from step 2; |r_M|≤1e-9 F_c on Taylor–Green; EMAC–VMS measured D_k, D_φ and 𝒩 ≤1e-12 of their scales; for the other forms |measured − predicted| ≤1e-9 of their scales. These criteria apply to the analysis script only; the 12a unit tests keep their own. Preserve nonconverging configurations and iteration counts. For Gresho report M, L, E and velocity-error evolution with wall forces and region norms as context. In ASGS the new momentum budget uses u_h+ũ, while the existing linear_momentum diagnostic continues to report u_h.

**12b-2, owner launch only.** Place this queue after the current validation queue. The owner launches it in tmux using the README command completed in step 4; the agent never launches it. Include every generated configuration, retaining failures. A later authorized phase reads the resulting tree.

**Consistent splitting: no code in this phase.** Chapter §5.7 distinguishes its own natural residual (5.25) from the implemented fixed-coupling residual. For the conservative form the resolved convective expression coincides with fixed coupling, but its natural residual adds (∇·U)u_h; equality of that resolved expression does not equate the full discrete schemes. Consistent non-conservative splitting generally conserves none of the three quantities. Consistent skew splitting needs ∇·ũ in the pointwise residual. Full consistent EMAC convection is the unimplemented full-stress variant in §5.8.5, whose energy-matched subscale model also differs from the natural-residual choice in §5.7. Minimal EMAC–VMS remains Definition 5.1 in §5.8.

### 12c: stabilization comparison (planned, not launched in 12a)

Use the same mesh, velocity/pressure pair, time scheme, stabilization constants and stress viscosity for standard skew-1 OSS, EMAC with standard OSS, and minimal EMAC–VMS. Configuration suffix {model} denotes those three choices. Existing Chapter 10 timings below are estimates for the standard solver, not guarantees for EMAC Newton convergence. Keep the physical-pressure outflow difference explicit in cylinder comparisons.

| Configuration family in configs/chapter5/stabilization/ | Settings and measurements | Estimated cost on 4 ranks | Queue placement |
| --- | --- | --- | --- |
| cylinder_re100_{model}_p1_l1.toml | Existing L1 mesh, Dyn-Nl-OSS P1/P1, dt=0.005, T=30; drag, lift, Strouhal and solver iterations | Existing baseline ≈1.16 h; allow 1–4 h/model, 3–12 h total | After 12b and the existing queue |
| cavity_re100_{model}_p1_64.toml | Existing 64×64 cavity, Dyn-Nl-OSS; same steady-stop rule as Chapter 10 | Minutes/model; allow 0.5 h for the trio | First in 12c by estimated cost |
| cavity_re1000_{model}_p1_64.toml | Dyn-Lin-OSS for all three, matching the feasible existing Re=1000 case; EMAC has no warm-up | Allow 0.5 h for the trio if Newton converges; report failure without changing the model | After Re=100; no nonlinear Re=1000 run is silently substituted |
| manufactured_{model}_p{k}_n{n}.toml | Manufactured spatial solution, k=1,2; n=8,16,32; short small-dt runs, physical-pressure and velocity convergence rates | Allow 0.5–2 h for 18 runs; verify temporal error is negligible | Before cylinder, after the cavity comparison |
| taylor_green_3d_re1600_{model}.toml (optional) | Requires a new 3D case and a separate resolution/time-horizon decision; compare resolved and subscale energy/dissipation | Not responsibly bounded before a unit-scale timing probe; potentially tens of hours | Last, separate optional batch; not part of the current queue |

Do not add these future runs to the active queue or launch them during 12a. Refine costs before each run phase; the Gresho frame is monitored as specified in 12b, and the optional 3D configuration still requires a decision.
