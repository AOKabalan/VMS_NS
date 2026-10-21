# Progress log

## Phase 0: Verify the spec addendum (2026-10-21)
Done: Initialised the git repository with a `.gitignore` for Python and `results/`. Copied the previous `VERIFICATION.md` unchanged to `docs/verification/VERIFICATION.md` and added an "Addendum pass" at the top, with rows A01–A15 in the earlier table format, derivations A-D1 to A-D6, and the disagreements and open questions. The pass checks the following against `formulation_spec.md`:
- BDF3 coefficients, exact for cubics.
- The Crank–Nicolson one-leg form against Codina et al. 2007 eqs. (22)–(23) with θ = ½ and their f^(n+θ).
- τ_t for Crank–Nicolson against their (25)–(26).
- Orthogonality of the extrapolated OSS subscale.
- The second-order limit of BDF3 with lower-order startup.
- The Dirichlet averaging.
- The Galerkin baseline as the master equation with all stabilization removed.
- The revised notation and the fixed-choice rows Variants, Interpolation, Time scheme and Quadrature.

Verdict: every addendum equation agrees and none needs to change. There are two wording errors in the stability paragraph, one missing Crank–Nicolson condition, and several stale references.

Tests: No code, as Phase 0 requires. Hand derivations were cross-checked with scratch calculations outside the repository:
- Exact rational sums of the BDF2 and BDF3 weights on monomials up to degree 5: BDF3 is exact through degree 3, with leading error −δt³g''''/4.
- The BDF3 boundary locus: A(α) with α ≈ 86.03°; unstable left-half-plane set inside −1/12 < Re z < 0, |Im z| < √15/2.

Source equations were read from `pdftotext` output. Pages 2417, 2421 and 2422 of the 2007 paper were checked on rendered images.

Deviations from spec or plan: None.

Open questions:
1. Crank–Nicolson initial compatibility (VERIFICATION.md A06, A08). The extrapolation u^(n+1) = 2u^(n+½) − u^n never damps a defect in u_h⁰. Dirichlet values are exact only if u_h⁰ matches the data. For Galerkin, a u_h⁰ that is not discretely divergence-free leaves a pressure component of size 2φ⁰/δt that flips sign every step. Either (i) require and construct a compatible u_h⁰, or (ii) start Crank–Nicolson with one BDF1 step, which keeps second order. Recommendation: (ii). It changes the spec's "Crank–Nicolson needs no startup". V2 is unaffected; impulsive starts and the thesis lattice-vortex setup are affected.
2. Exact-history startup levels (A10). Does "velocity histories initialised from the exact solution" mean u⁻¹, u⁻² at t = −δt, −2δt (BDF3 from step 1), or u¹, u² (BDF3 from step 3)? Recommendation: negative times. Both give third order in V2.
3. Spec wording (A11, A14). Two stability sentences need replacing (proposed text in VERIFICATION.md), and the references to t^(n+1) and to two history levels at spec lines 5, 100, 117, 162–168 and 239 are stale. The spec is not edited by the agent; please apply these or confirm they can wait.

## Phase 1: Skeleton, configuration and run infrastructure (2026-10-21)
Done: Added `pyproject.toml` (setuptools, Python ≥ 3.11 as Firedrake 2026.10.0 requires, dependency `tomli-w`, extras `dev` = pytest + ruff, ruff line length 100), a README with install, run and the symbol table, and the package `navier_stokes_vms`.
- `config.py` holds frozen option dataclasses for every field of plan §7. Allowed values are `Literal` types that the loader also checks.
- TOML loading fills defaults and rejects unknown sections or keys, wrong types, values outside the allowed set, and out-of-range numbers. It also checks that the mesh source matches the file and the builtin fields, that cells and lengths have one entry per dimension, and that case and run_name are safe folder names.
- Rules 1, 2, 4 and 5 of plan §7 are enforced on load. Rule 5 accepts only `method = "none"` for now and names the phase that adds each VMS case.
- `config_to_toml` writes the resolved configuration, which reloads to an equal Config.
- `postprocessing/output.py` creates `results/<case>/<run_name>/` on rank 0 and refuses an existing folder. It writes `config_resolved.toml`, and writes `run_info.json` at start and again at the end with end time, wall time and exit status.
- `postprocessing/run_logging.py` sends the package log to the console at `console_level` and to `run.log` at INFO, on rank 0 only, and writes the setup summary.
- `scripts/run_simulation.py` loads and validates, prepares the folder, writes the three files, logs "No solver is implemented yet; stopping after setup." and exits with status 0. An invalid configuration or an existing folder gives status 1 and writes nothing.

Tests:
- `ruff check .`: all checks passed.
- `ruff format --check .`: 11 files already formatted.
- `pytest tests/unit`: 39 passed in 1.6 s. The 39 tests cover defaults, a full config with every field, the resolved dump round trip, one test per rule 1–6 (rule 5 for all eight VMS cases with their phase), 18 invalid-input cases, and the run script run as a subprocess: three files written, existing folder refused, invalid config rejected before writing.
- Manual: a serial run and `mpiexec -n 2` (one log copy, files from rank 0 only, `mpi_size = 2` in `run_info.json`, consistent exit status 1 on an existing folder).

Deviations from spec or plan:
1. Rules 3 and 6 are implemented and tested as `check_startup_supported(config, has_exact_solution)` and `check_diagnostics_known(config, available)` in `config.py`, but nothing calls them yet. Rule 3 needs the case's exact solution (problem building, phase 2) and rule 6 needs the diagnostics registry (phase 4), so the run script does not enforce them yet. They are not stubs. Question asked in the phase report: wire them in phases 2 and 4?
2. Plan §7 gives an example rather than a list of defaults. The example's values are taken as the defaults, except the case-specific ones: `run_name` defaults to the config file stem, `case_parameters` to empty, and `[mesh]` to builtin `rectangle`, cells [32, 32], lengths [1.0, 1.0]. As a result, a config that sets only `case` is VMS Dyn-Nl-OSS on P1/P1 and is rejected by rule 5 until phase 8.
3. The plan does not say what to do when a results folder exists; the run refuses it rather than overwriting.
4. The package was installed into the Firedrake venv with `pip install -e ".[dev]"`, which added `tomli-w` and `ruff`.

Open questions:
1. Wire rules 3 and 6 in phases 2 and 4 as described above?
2. Refuse or replace an existing results folder? Refusal is the current, safer choice.
3. Heads-up for phase 2: Kovasznay is posed on [−0.5, 1] × [−0.5, 1.5], but `[mesh]` has only `lengths`. Either phase 2 adds an `origin` field to `[mesh]` or the case shifts coordinates.

## Phase 2: Problems, meshes and boundaries (2026-10-21)
Done: Before this phase the plan was amended (commit "Plan: mesh geometry belongs to the case"). `[mesh]` now has only optional `source`, `cells` (builtin only) and `file` (gmsh only), each defaulting to the case's choice; `config.py` and its tests were updated to match.
- `problem/mesh_loading.py` provides `rectangle`, `periodic_rectangle` (periodic in both directions, no boundary names) and `box`, each on [lower, upper]. It holds the name → tag tables `RECTANGLE_TAGS` (x_min 1, x_max 2, y_min 3, y_max 4) and `BOX_TAGS` (adds z_min 5, z_max 6). `gmsh_mesh` reads the facet-dimension `$PhysicalNames` entries. `builtin_cells` returns the config override or the case default and rejects a gmsh source or a wrong number of cells.
- `problem/boundary_conditions.py` has `Dirichlet(value)` and `Outflow()`, plus `dirichlet_bc(space, name, value, boundary_tags)`. An unknown name raises an error listing the available ones.
- `problem/problem_definition.py` has the frozen `Problem` dataclass of plan §8 and `exact_solution_forcing`, which derives f = ∂t u + u·∇u − νΔu + ∇p in UFL with ∂t through `ufl.variable`. It also has `case_parameters` (case defaults overridden by `[case_parameters]`, unknown names rejected), `check_problem` and `uses_pressure_null_space` (true exactly when no boundary is `Outflow`). `check_problem` requires every named boundary to have exactly one known condition and checks that an exact velocity is divergence-free: ‖div u‖ ≤ 1e-10 ‖∇u‖ at t = 0 and 1.
- `cases/__init__.py` holds the registry `CASES` and `build_problem(config)`. That function rejects an unknown case (listing the available ones), runs `check_problem`, and applies plan §7 rule 3 (exact-history startup needs an exact solution), as decided after phase 1.
- Cases, all built-in only, with Dirichlet exact data on every side except Taylor–Green:

  | Case | Domain | Default cells | Parameter (default) | Solution |
  | --- | --- | --- | --- | --- |
  | `manufactured_spatial` | unit square | 16×16 | viscosity 0.01 | u = (1 + t) curl ψ with ψ = sin²(πx) sin²(πy), built with UFL derivatives; p = (1 + t) cos(πx) sin(πy) |
  | `manufactured_temporal` | unit square | 4×4 | viscosity 0.01 | u = cos t (y, x), p = 0 |
  | `kovasznay` | [−0.5, 1] × [−0.5, 1.5] | 24×32 (h = 1/16) | reynolds 40 | ν = 1/Re, λ = Re/2 − √(Re²/4 + 4π²) |
  | `taylor_green_2d` | periodic (0, 2π)² | 32×32 | viscosity 0.01 | no boundaries |

- `tests/data/square_named.msh` is the unit square with h = 0.25, physical curves left 11, right 12, bottom 13, top 14 and physical surface fluid 20. The tags are deliberately not 1–4. It was generated with gmsh 4.14 (`gmsh -2 -format msh41`).

Tests:
- `ruff check .` and `ruff format --check .` are clean.
- `pytest tests/unit`: 68 passed in 4.9 s (44 config, 24 in `test_meshes_and_boundaries.py`).
- U7:
  - Rectangle and box with a non-zero lower corner: each name's facets lie on its coordinate plane (mean of x_i over ds(tag) equals the plane value, spread 0), and the tagged lengths sum to the whole boundary.
  - The periodic rectangle has no ds and the right area and centre.
  - The gmsh mesh loads {left: 11, right: 12, bottom: 13, top: 14}, each on its side, with the surface name excluded.
  - A `DirichletBC` by name sets exactly the nodes of that boundary and leaves every other node zero (builtin and gmsh).
  - An unknown name lists the available ones.
- Phase 2 acceptance:
  - Every case builds on a 4×4 mesh: 32 cells, initial velocity equal to the exact one at t = 0, finite forcing, boundary conditions covering all named boundaries, pressure null space on.
  - The derived forcing vanishes to round-off for `taylor_green_2d` (and also for `kovasznay`): ‖f‖ < 1e-12 ‖∇p‖, about 1e-15 measured.
  - Extra checks: the derived forcing equals a hand derivation for `manufactured_temporal`; rule 3 is triggered through `build_problem`; invalid case setups are rejected; `check_problem` rejects unknown, unassigned and divergent cases; the null-space decision is correct.
- Manual: all four cases and the gmsh mesh built on 1 and 2 MPI ranks give the same global cell counts and forcing norms.

Deviations from spec or plan:
1. `check_problem` also requires every named boundary to have a condition. The plan does not say this, but otherwise an unassigned boundary would silently act as an outflow without the ½(u·n)u term, and the null-space decision would be wrong.
2. For p = 0 (`manufactured_temporal`), UFL cannot take the gradient of a constant without a mesh, so `exact_solution_forcing` uses a zero ∇p for any pressure without spatial dependence.
3. All meshes are built at the origin and then translated. In Firedrake 2026.10.0, `RectangleMesh(nx, ny, Lx, Ly, originX, originY)` treats Lx and Ly as the upper corner rather than the lengths, and `BoxMesh` and `PeriodicRectangleMesh` take no origin.
4. The `.geo` source of the test mesh is not committed, since the plan lists only the `.msh`; its content is described above.
5. The run script does not build the problem yet; it is not in the phase 2 file list. Rule 3 and an unknown case name are therefore caught only when the time stepper (phase 4) calls `build_problem`.

Open questions:
1. Defaults not fixed by the plan: ν = 0.01 for `manufactured_temporal` (V2 gives no ν; the solution has Δu = 0, so ν only enters τ) and the default cells listed above. Please confirm or change.
2. Firedrake 2026.10.0 makes `topological_dimension` and `geometric_dimension` properties, not methods, so plan §3's `mesh.geometric_dimension()` must be read as `mesh.geometric_dimension`.
3. Noise: when pytest is called with flags such as `-q` or `-x`, PETSc reports them as unused options at exit. This is harmless; a `tests/conftest.py` that initialises petsc4py with no arguments would remove it, if you want that file added.

## Phase 3: Quadrature and time schemes (2026-10-21)
Done: Before this phase, commit "Phase 2 follow-ups" made three changes. Plan §3 now uses the property `mesh.geometric_dimension`. `tests/conftest.py` sets `FIREDRAKE_DISABLE_OPTIONS_LEFT=1`, Firedrake's own switch, so PETSc no longer reports pytest's flags as unused options. `tests/data/square_named.geo` was added and regenerates the committed `.msh` byte for byte.
- `discretization/quadrature.py`:
  - `quadrature_degree(k) = 3k`.
  - `volume_measure(mesh, q)` returns `dx(degree=q, scheme="default")` after `check_positive_weights`; `boundary_measure(mesh, q)` returns the matching `ds`.
  - `quadrature_element(cell, q)` returns a scalar `Quadrature` element with the same scheme.
  - `reference_rule(cell, q)` gives the reference points and weights. They come from the same `finat.quadrature.make_quadrature(cell, q, scheme)` that TSFC uses for `dx` and FInAT uses for the Quadrature element, so the check inspects the rule actually used.
  - The default rules have positive weights for q = 3, 6, 9 on triangles and tetrahedra; the smallest are 8.3e-2, 2.0e-2 and 1.3e-2 on triangles, and 2.0e-2, 1.2e-3 and 7.2e-4 on tetrahedra.
- `discretization/time_schemes.py` is pure Python; importing it loads neither Firedrake nor petsc4py.
  - `scheme_at_step(scheme, startup, step)` gives the startup sequence: lower-order is BDF1 → BDF2 → BDF3; Crank–Nicolson takes one BDF1 first step in both modes; exact-history runs BDF2 and BDF3 at full order from step 1.
  - `step_coefficients(...)` returns a frozen `StepCoefficients` with the scheme used, α₀, the history weights of gⁿ, gⁿ⁻¹, gⁿ⁻², the solve level (1 or ½), the data weights and the extrapolation weights.

  | Scheme | α₀ | History weights | Solve level | Data weights (new, old) | Extrapolation (solve, old) |
  | --- | --- | --- | --- | --- | --- |
  | BDF1 | 1 | (1) | 1 | (1, 0) | (1, 0) |
  | BDF2 | 3/2 | (2, −½) | 1 | (1, 0) | (1, 0) |
  | BDF3 | 11/6 | (3, −3/2, 1/3) | 1 | (1, 0) | (1, 0) |
  | Crank–Nicolson | 2 | (2) | ½ | (½, ½), the endpoint average | (2, −1) |

  With exact-history startup, the number of history weights at step 1 is the number of velocity levels to fill: t⁰, t⁻¹ and t⁻² for BDF3.

Tests:
- `ruff check .` and `ruff format --check .` are clean.
- `pytest tests/unit`: 107 passed in 3.6 s (39 new).
- U1:
  - BDFk is exact for polynomials of degree ≤ k at t^(n+1), and inexact at degree k+1.
  - Crank–Nicolson is exact for quadratics at t^(n+½) when the half-level state is the endpoint average, and inexact for cubics.
  - The extrapolation of midpoint values reproduces linear functions, and not quadratics.
  - Both startup sequences are checked for all four schemes; exact-history needs back to t⁻² for BDF3.
  - The Crank–Nicolson data weights are (½, ½), consistent with the extrapolation. Every step has history weights summing to α₀ (D_t of a constant is zero), and BDF steps have identity data and extrapolation weights.
  - Step 0 is rejected.
- U2:
  - On the reference triangle and tetrahedron, for q = 3 and 6, all weights are positive and every monomial of degree ≤ q integrates exactly (relative 1e-13).
  - A non-polynomial function interpolated into the Quadrature space and integrated with `dx` equals the direct `dx` integral to relative 1e-14, in 2D and 3D. A rule of another degree differs by more than 1e-9, so this agreement is meaningful.
- Extra checks:
  - Both measures carry the explicit degree and scheme in their metadata, and `ds` of degree 6 is exact for x⁶.
  - q = 3k.
  - A rule with a negative weight is rejected. The test uses the degree-3 Strang–Fix rule, the counterexample from VERIFICATION.md R-D2.

Deviations from spec or plan:
1. Plan §3 asks for one-line docstrings, so `StepCoefficients` explains itself through field names (`history_weights`, `data_weight_new`, `extrapolation_weight_solve`, …) rather than a longer docstring.
2. The positivity check covers volume rules only. The spec requires positive weights for the τ-weighted volume mass matrices; the facet rule only needs an explicit degree.

Open questions: None.

## Phase 4a: Galerkin solver, time stepper and outputs (2026-10-21)
Done:
- `discretization/function_spaces.py` builds the Taylor–Hood mixed space and the constant-pressure null space.
- `formulation/convective_forms.py` has `convection_form(a, u_h, u, v, dx, outflow_ds)` = ½(a·∇u, v) − ½(a·∇v, u) + ½⟨(u_h·n)u, v⟩ on Γ_N.
- `formulation/galerkin.py` has `galerkin_residual` (the master equation with every stabilization term removed), `time_derivative` D_t g = (α₀g − g_hist)/δt, and `SolveLevel`.
- `solver/nonlinear_solver.py` has `InnerSolver`:
  - Newton is SNES `newtonls` with the `basic` line search (full steps, u += δ). A petsc4py monitor records the residual, update norm and KSP iterations.
  - Picard uses the same residual builder with `a` a frozen Function, solved with `ksponly`; `picard_iterations_before_newton` Picard steps can precede Newton.
  - The Picard residual_norm is the full nonlinear residual assembled after each linear solve, with the Dirichlet rows zeroed. The test suite confirms that SNES reports the same norm for Newton.
  - One tolerance max(inner_atol, inner_rtol·r₀) applies to every method, where r₀ is the explicit residual at the guess after the BCs are applied. Newton uses a petsc4py convergence test for it.
- `solver/time_stepper.py` has `TimeStepper`, which takes everything time-related from `time_schemes.py`:
  - It fills the exact history at t⁻¹, t⁻² by interpolation, combines the history, sets α₀, and puts both the forcing and the Dirichlet values at the solve level with the data weights (the endpoint average for Crank–Nicolson).
  - Crank–Nicolson solves at the half level and extrapolates; its pressure is reported at t^(n+½).
  - It writes the plan §9 logs and VTK. On non-convergence it writes the last state to `fields/`, logs the iteration history and raises `ConvergenceError`.
- `postprocessing/diagnostics.py` holds the registry `kinetic_energy`, `divergence_l2`, `velocity_error_l2`, `velocity_gradient_error_l2` and `pressure_error_l2`.
  - The error norms use degree 2q; everything else uses the common rule.
  - The pressure error removes the mean in two passes when the null space is active, and is evaluated at the solve-level time.
  - `check_diagnostics_supported` rejects error norms for a case without an exact solution.
- `postprocessing/output.py` has `CsvWriter` (rank 0, flushed per row), `FieldWriter` (fields/velocity.pvd and fields/pressure.pvd, each with its own time) and `open_diagnostics_csv`. `postprocessing/run_logging.py` has the three iteration and time-step CSVs with the plan §9 columns.
- `problem/boundary_conditions.py` adds `named_boundary_measure` (ds by name, keeping the degree).
- `scripts/run_simulation.py` checks rule 6 against the registry and builds the problem (rule 3) before creating the folder, then runs end to end. It exits 1 on non-convergence.
- `config.py` now requires `t_end` to be a whole number of steps. The README run section was updated.

Tests:
- `ruff check .` and `ruff format --check .` are clean.
- `pytest tests/unit`: 124 passed in 18.8 s.
- U3 (`test_convective_form.py`): c(a; w, w) = 0 to 1e-12 relative for random a and w on periodic meshes in 2D and 3D, P1 and P2, computed through the assembled matrix so that nothing cancels symbolically. With outflow on two sides, c(u; u, u) equals ½∫(u·n)|u|² to 1e-12 relative.
- `test_time_stepper.py` (new): short runs of all four cases.
  - The runs cover BDF1 with Newton; Crank–Nicolson with exact history and Picard; BDF2 with 2 Picard steps then Newton; BDF3 with exact history.
  - Every plan §9 file is written. The scheme sequence follows `scheme_at_step`. Inner iteration 0 has relative residual 1 and no update; the last iterate meets the tolerance; the method columns are as configured; there is one outer row per step; the diagnostics columns are right; VTK times are correct.
  - The Crank–Nicolson pressure is written at t^(n+½) (0.375, 0.625) and the velocity at the endpoints.
  - A run with one Picard iteration exits 1, writes the last state to `fields/` with `vtk_every = 0`, logs the iteration history, and records `converged = False` and `exit_status = 1`.
  - Newton, Picard and 3 Picard + Newton reach the same steady Kovasznay solution (≤ 1e-8 relative). The logged final residual equals `InnerSolver.residual_norm()`. Newton logs update norms and KSP iterations.
- `test_config.py`: the run-script tests now run one step and check all files. They also cover rejection before writing for rule 5, rule 6 and an unknown case, and the new whole-steps check.
- Manual checks:
  - On 2 MPI ranks, the Kovasznay run gives the same diagnostics as 1 rank to about 1e-15, CSVs from rank 0 only and parallel VTK (`.pvtu`).
  - Taylor–Green, P2/P1, BDF2, dt = 0.1, T = 0.3: velocity error 2.4e-1, 4.8e-2 and 4.9e-3 on meshes 8, 16 and 32 (tending to O(h³)). The kinetic energy on mesh 32 matches π²e^(−4νT) to 1e-4.
  - MUMPS LU with the null space needs no null-pivot option here.

Deviations from spec or plan:
1. `outer_iterations.csv` is written now, with one row per time step, because phase 4a must write every file of plan §9 and the Galerkin step is a single pass. The row records the velocity and pressure changes from the guess, the final residual and `converged`. The subscale, ξ and ζ changes are empty. Phase 7 adds further outer iterations.
2. VTK: velocity and pressure are separate series with their own times, so the Crank–Nicolson pressure can carry t^(n+½). Frames are written at multiples of `vtk_every` and at the last step, but not at step 0, where the pressure is only the initial guess.
3. `diagnostics.csv` has rows for steps 1 to N, not step 0.
4. Times are rounded to 15 significant digits (step·δt), so they print as 0.3 rather than 0.30000000000000004.
5. Not in the plan: the whole-steps check in `config.py`, the README update, and the new test file `tests/unit/test_time_stepper.py`.

Open questions:
1. V2 cannot measure Galerkin temporal orders with the specified `manufactured_temporal`. Its velocity g(t)(y, x) = g(t)∇(xy) is a gradient, so the temporal truncation error ε(t)(y, x) is a gradient too. The pressure absorbs it, and the velocity comes out exact to round-off for every scheme. A scratch probe (4×4, P2/P1, T = 1, N = 20, 40, 80) gave:

   | Field | Scheme | ‖u − u_ex‖ at N = 20, 40, 80 | Velocity rates | Pressure rates |
   | --- | --- | --- | --- | --- |
   | g(t)(y, x), as specified | BDF1 | 1.4e-16, 2.7e-15, 2.2e-12 | meaningless | 1.02, 1.01 |
   | g(t)(y, x), as specified | BDF2 | 7.3e-15, 1.9e-14, 1.0e-12 | meaningless | 1.98, 1.99 |
   | g(t)(y, x), as specified | BDF3, exact history | 5.5e-16, 2.5e-13, 7.1e-13 | meaningless | 3.06, 3.03 |
   | g(t)(y, x), as specified | Crank–Nicolson | 1.6e-15, 8.5e-14, 6.0e-13 | meaningless | 1.98, 1.99 |
   | g(t)(y, −x), rotation | BDF1 | 3.2e-3, 1.6e-3, 8.2e-4 | 0.99, 0.99 | 1.01, 1.01 |
   | g(t)(y, −x), rotation | BDF2 | 1.1e-4, 2.7e-5, 6.6e-6 | 2.03, 2.01 | 1.98, 1.99 |
   | g(t)(y, −x), rotation | BDF3, exact history | 4.3e-6, 5.2e-7, 6.5e-8 | 3.02, 3.01 | 3.06, 3.03 |
   | g(t)(y, −x), rotation | Crank–Nicolson | 6.5e-5, 1.6e-5, 4.1e-6 | 2.00, 2.00 | 1.99, 2.00 |

   Recommendation: change `manufactured_temporal` to u = g(t)(y, −x), p = 0. It is still linear, divergence-free and in every FE space, and its rotational truncation error cannot be absorbed by the pressure. This changes the plan's Phase 2 case definition and the meaning of V2, so it needs your decision before phase 4b. The alternative is to keep the field and measure the pressure error in V2.

## Phase 4b: Galerkin verification (2026-10-21)
Done: Before this phase, commit "Plan: rotational manufactured_temporal" changed `manufactured_temporal` to u = cos t (y, −x), p = 0, together with its hand-derived forcing test. It also updated plan §13 Phase 2 and §10 V2: the field must have nonzero vorticity, and V2 also asserts pressure orders, comparing against the exact pressure at the solve level (t^(n+1) for BDF, t^(n+½) for Crank–Nicolson) with zero-mean normalization.
- `scripts/run_convergence_study.py` takes a base TOML and either `--cells NxM ...` or `--steps N ...` (dt = t_end/N). For each level it writes the level's resolved configuration into `results/<case>/<run_name>_study/` and runs it with `run_simulation.py` as a subprocess, so every level is an ordinary run with all plan §9 files. It then reads each run's final diagnostics and writes `convergence.csv` (level, refinement, value and rate per diagnostic; rate = log(e_{i−1}/e_i)/log(N_i/N_{i−1})). The table is also logged. The script refuses an existing study folder and stops at the first failed level.
- `TimeStepper` gained read-only `velocity` (endpoint) and `pressure` (solve level) properties, used by V3.
- Verification tests in `tests/verification/`, all marked `slow`, with the helper `study.py`:
  - V1 `test_spatial_convergence.py`: P2/P1, BDF1, two steps of dt = 0.05, meshes 8, 16, 32.
  - V2 `test_temporal_convergence.py`: 4×4, P2/P1, T = 1, N = 20, 40, 80, 160 for BDF1, BDF2, BDF3 (exact history) and Crank–Nicolson; BDF3 with lower-order startup is logged and checked at ≥ 1.8.
  - V3 `test_steady_equivalence.py`: Kovasznay Re = 40, BDF1 with dt = 5 to T = 100, meshes "16" and "32" taken as h = 1/16 and h = 1/32 (cells 24×32 and 48×64).
  - V4 `test_taylor_green.py`: BDF2 with dt = 0.025 to T = 1, meshes 16 and 32.
- The README documents the verification and study commands.

Tests:
- `ruff check .` and `ruff format --check .` are clean.
- `pytest tests/unit`: 124 passed.
- `pytest tests/verification -m slow`: 9 passed in 119 s.

The rates are between the two finest levels; the threshold is the expected rate − 0.2.

V1, spatial (manufactured_spatial, Taylor–Hood P2/P1):

| Mesh | ‖u − u_ex‖ | Rate | ‖∇(u − u_ex)‖ | Rate | ‖p − p_ex‖ | Rate |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 1.436e-2 | | 9.035e-1 | | 1.142e-2 | |
| 16 | 1.884e-3 | 2.93 | 2.327e-1 | 1.96 | 1.910e-3 | 2.58 |
| 32 | 2.093e-4 | **3.17** (≥ 2.8) | 5.091e-2 | **2.19** (≥ 1.8) | 4.442e-4 | **2.10** (≥ 1.8) |

V2, temporal (manufactured_temporal, P2/P1, 4×4, T = 1). Velocity at T; pressure at the solve level, zero mean:

| Scheme | N | ‖u − u_ex‖ | Rate | ‖p − p_ex‖ | Rate |
| --- | --- | --- | --- | --- | --- |
| BDF1 | 20 / 40 / 80 / 160 | 3.223e-3 / 1.627e-3 / 8.181e-4 / 4.102e-4 | 0.99, 0.99, **1.00** (≥ 0.8) | 3.006e-3 / 1.489e-3 / 7.409e-4 / 3.696e-4 | 1.01, 1.01, **1.00** (≥ 0.8) |
| BDF2 | 20 / 40 / 80 / 160 | 1.080e-4 / 2.645e-5 / 6.555e-6 / 1.633e-6 | 2.03, 2.01, **2.01** (≥ 1.8) | 1.342e-4 / 3.394e-5 / 8.534e-6 / 2.140e-6 | 1.98, 1.99, **2.00** (≥ 1.8) |
| BDF3, exact history | 20 / 40 / 80 / 160 | 4.262e-6 / 5.237e-7 / 6.488e-8 / 8.073e-9 | 3.02, 3.01, **3.01** (≥ 2.8) | 3.994e-6 / 4.803e-7 / 5.884e-8 / 7.278e-9 | 3.06, 3.03, **3.02** (≥ 2.8) |
| Crank–Nicolson | 20 / 40 / 80 / 160 | 6.520e-5 / 1.624e-5 / 4.058e-6 / 1.014e-6 | 2.00, 2.00, **2.00** (≥ 1.8) | 1.079e-4 / 2.712e-5 / 6.801e-6 / 1.703e-6 | 1.99, 2.00, **2.00** (≥ 1.8) |
| BDF3, lower-order startup (logged) | 20 / 40 / 80 / 160 | 9.039e-5 / 2.249e-5 / 5.627e-6 / 1.408e-6 | 2.01, 2.00, **2.00** (≥ 1.8) | 2.940e-5 / 7.418e-6 / 1.870e-6 / 4.698e-7 | 1.99, 1.99, **1.99** (≥ 1.8) |

V3, steady Kovasznay (Re = 40, P2/P1, BDF1 dt = 5 to T = 100; final relative step change 0 < 1e-10 in every run):
- Newton and Picard steady velocities differ by 4.4e-13 relative (≤ 1e-8).
- Errors from mesh 16 (24×32) to mesh 32 (48×64): velocity 4.096e-4 → 5.112e-5 (×8.0, rate 3.0); pressure 5.152e-4 → 1.276e-4 (×4.0, rate 2.0). Both drop, as required.

V4, Taylor–Green (ν = 0.01, P2/P1, BDF2 dt = 0.025 to T = 1):

| Mesh | ‖u − u_ex‖ at T | Kinetic energy at T |
| --- | --- | --- |
| 16 | 7.826e-2 | 9.47098 |
| 32 | 4.720e-3 (drop ×16.6, ≥ 3) | 9.48217 (exact π²e^(−0.04) = 9.48261; 0.005 %, ≤ 1 %) |

Deviations from spec or plan:
1. Choices the plan leaves open:
   - V1 uses dt = 0.05; with a solution linear in time, any dt gives the same spatial error.
   - V3 reads "mesh 16 and 32" as h = 1/16 and 1/32 on the 1.5 × 2 Kovasznay domain. It reaches steady state with dt = 5 and T = 100 and checks the final relative step change from `outer_iterations.csv`.
   - V4 uses dt = 0.025.
2. Not in the plan: the verification helper `tests/verification/study.py`, the `TimeStepper.velocity` and `.pressure` properties, and the README section. The study script reuses `run_simulation.py` through a subprocess per level, so the runs are serial even if the study itself is started under MPI; the README says to run it serially.
3. V2's pressure assertion follows the instruction given with the V2 decision. The plan's V2 text was updated with it in the "Plan: rotational manufactured_temporal" commit.

Open questions: None.
