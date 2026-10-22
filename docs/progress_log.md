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

## Phase 5: Static linear ASGS (2026-10-22)
Done:
- `formulation/operators.py` has `elementwise_laplacian` (Δ_h = div grad inside cell integrals), `momentum_operator` L_a(u, p) = a·∇u − νΔ_h u + ∇p, `adjoint_operator` L*_a(v, q) = −a·∇v − νΔ_h v − ∇q, and `momentum_residual` R_u = f − D_t u − L_a(u, p).
- `formulation/stabilization_parameters.py` evaluates τ_m = (c1·k⁴·ν/h² + c2·k·|a|/h)⁻¹ and τ_c = h²/(c1·k⁴·τ_m) pointwise, with |a| = √(a·a), h = `ufl.MinCellEdgeLength` (the shortest cell edge) and no 1/δt.
- `formulation/vms.py` has `vms_terms` = Σ_K (ũ, L*_a(v, q))_K + Σ_K (τ_c ∇·u, ∇·v)_K with ũ = τ_m R_u, the master formulation for Sta-Lin-ASGS (a = u_h, ζ = 0). `pressure_subscale = false` sets τ_c to zero. Expanded, this is +τ_m(D_t u + L_a u − f, a·∇v + νΔ_h v + ∇q), the corrected sign of the spec's Errata row 1.
- `config.py`: rule 5 now accepts exactly (asgs, static, linear).
- `solver/time_stepper.py` adds `vms_terms` to both residuals. Newton evaluates a and τ at the current iterate; Picard evaluates them at the frozen advector. The equal-order spaces Pk/Pk come from the existing `mixed_space`.
- Newton τ options (`solver/nonlinear_solver.py`, `time_stepper.py`):
  - The residual, and therefore every logged residual, always uses τ at the current iterate.
  - `newton_tau = "frozen"`: the Jacobian is `derivative` of the residual with τ built on the frozen Function, then `ufl.replace(frozen → u)`. UFL's `replace` expands the derivative first, so the Jacobian has τ(u) but no τ′ (inexact Newton).
  - `newton_tau = "exact"`: Firedrake's derivative of the full form, including τ′.
  - With `exact`, a post-Jacobian callback checks that the Jacobian is finite. A non-finite residual or Jacobian produces a failure reason naming `solver.newton_tau = 'exact'` ("differentiates |a|, which is undefined at a = 0").
- PETSc error-wrapping fix: petsc4py turns an exception raised inside the Jacobian callback into a `PETSc.Error`. `InnerSolver` records the failure reason before raising. It catches `PETSc.Error` and reports non-convergence only when a reason was recorded, and re-raises otherwise. The named message thus reaches the normal failure path: last state written to `fields/`, iteration history logged, the reason appended to the run's `ConvergenceError`, exit status 1. `InnerSolveResult` gained `failure_reason` for this.

First stop and decision:
- With h = cell diameter, random-sample U5 passed, but a targeted P2 field broke dissipation. The field was on a 4×4 mesh with zero boundary values, p = 0, f = 0, no time terms, ν = 1 and pressure subscales on. Its viscous term was 7.50114e-8, and viscous + stabilization was −1.34586e-8.
- With the cell diameter, V1 P2/P2 lost convergence on mesh 32:

  | Mesh | ‖u − u_ex‖ | ‖∇(u − u_ex)‖ | ‖p − p_ex‖ |
  | --- | --- | --- | --- |
  | 8 | 1.251e-2 | 7.307e-1 | 1.058e-2 |
  | 16 | 1.596e-3 | 1.865e-1 | 6.391e-4 |
  | 32 | 1.941e-3 | 3.443e-1 | 2.020e-3 |
  | Rate 16 → 32 | −0.28 | −0.88 | −1.66 |

- V3 P2/P2 Picard needed more than 100 iterations on mesh 32. Work stopped there.
- Decision (spec owner): for k ≥ 2, the adjoint term contributes −Σ_K τ_m ν² ‖Δ_h u‖²_K. Stability needs c1·k⁴ (64 for P2) to exceed the inverse constant C_inv with ‖Δ_h v‖²_K ≤ C_inv h⁻² ‖∇v‖²_K. On right isosceles triangles C_inv = 48 with h the shortest edge, but 96 with h the cell diameter.
- So h = shortest cell edge (Colomés §7.1.2, minimum element length). c1 = 4 and c2 = 2 are unchanged.
- The spec now has the paragraph "Stability condition for k ≥ 2" and the changed τ-evaluation row. Plan §11 says `ufl.MinCellEdgeLength`, and plan §10 U5 describes the eigenvalue test. `ufl.MinCellEdgeLength` works directly in Firedrake 2026.10.0 in 2D and 3D.

Tests:
- `ruff check .` and `ruff format --check .` are clean (45 files).
- `pytest tests/unit`: 159 passed in 27 s.
- `pytest tests/verification -m slow`: 13 passed in 186 s.
- U4 (`test_asgs_equals_supg.py`): Galerkin + `vms_terms` equals a classical form written in the test, Galerkin + Σ τ_m(D_t u + a·∇u − νΔ_h u + ∇p − f, a·∇v + νΔ_h v + ∇q) + Σ τ_c(∇·u, ∇·v), with an independent τ, to relative 1e-12. The fields u, p, history and f are random, with α₀ = 1.5 and δt = 0.07 so that D_t contributes. It covers P1/P1 and P2/P2, 2D and 3D, τ_c on and off. For P2 it checks that ‖Δ_h u‖ is non-trivial.
- U5 (`test_dissipation_sign.py`), eigenvalue test:
  - E(u, p; v, q) = ν(∇u, ∇v) + Σ_K (−τ_m L_a(u, p), L*_a(v, q))_K + Σ_K (τ_c ∇·u, ∇·v)_K, with a frozen random advector (components in [0.5, 1.5]) and f = 0, no time terms.
  - The symmetric part is restricted to the dofs free of the homogeneous velocity Dirichlet conditions.
  - It covers P1/P1 and P2/P2 on 4×4 squares with diagonals `left` and `crossed` and on a 2×2×2 box, at ν = 10 (viscous regime) and ν = 1e-4 (convective regime).
  - In all 12 cases λ_min is at round-off, |λ_min| ≤ 1.5e-15 (the constant-pressure mode). λ_max ranges from 0.26 to 174, so λ_min ≥ −1e-12 λ_max holds.
- U5 regression: with the cell diameter, the at-rest P2 energy matrix has a negative eigenvalue. Its eigenvector, scaled by 1e-4, gives energy −1.34586e-8 with the diameter τ and +3.39009e-8 with the shortest-edge τ. The test detects the original failure.
- `test_asgs_solver.py` (new):
  - Newton Jacobians for frozen and exact τ, P1 and P2, match central finite differences of the corresponding residual to 1e-8 relative.
  - For Newton (frozen and exact), Picard, and 2 Picard + Newton, every logged residual equals the assembled ASGS residual with τ at the current iterate (1e-7), and the last meets the tolerance.
  - Started from rest, frozen τ converges. Exact τ fails with a `ConvergenceError` naming `solver.newton_tau = 'exact'`, writes `fields/` and records `converged = False`.
- `test_config.py`: rule 5 accepts static linear ASGS at P1/P1 and P2/P2 with the pressure subscale on and off. The other seven cases are still rejected with their phase. The ASGS config round-trips.

V1, Sta-Lin-ASGS (manufactured_spatial, BDF1, two steps of dt = 0.05; rate on 16 → 32, threshold expected − 0.2):

| Pair | Mesh | ‖u − u_ex‖ | Rate | ‖∇(u − u_ex)‖ | Rate | ‖p − p_ex‖ | Rate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P1/P1 | 8 | 2.556e-1 | | 4.927e+0 | | 6.262e-1 | |
| P1/P1 | 16 | 5.988e-2 | 2.09 | 2.467e+0 | 1.00 | 1.562e-1 | 2.00 |
| P1/P1 | 32 | 1.327e-2 | **2.17** (≥ 1.8) | 1.230e+0 | **1.00** (≥ 0.8) | 3.014e-2 | **2.37** (≥ 0.8) |
| P2/P2 | 8 | 1.197e-2 | | 7.291e-1 | | 9.943e-3 | |
| P2/P2 | 16 | 1.519e-3 | 2.98 | 1.843e-1 | 1.98 | 9.236e-4 | 3.43 |
| P2/P2 | 32 | 1.882e-4 | **3.01** (≥ 2.8) | 4.512e-2 | **2.03** (≥ 1.8) | 1.654e-4 | **2.48** (≥ 1.8) |

V1 Galerkin P2/P1 is unchanged: 3.17, 2.19, 2.10.

V3, Sta-Lin-ASGS on steady Kovasznay (Re = 40, BDF1 dt = 5 to T = 100, mesh 16 = 24×32 and mesh 32 = 48×64; final relative step change < 1e-10 in every run). The Newton/Picard difference is ‖u_N − u_P‖/‖u_N‖ in L² with the common rule:

| Pair | Mesh | Newton/Picard | Newton inner (max per step, total) | Picard inner (max per step, total) | ‖u − u_ex‖ | ‖p − p_ex‖ |
| --- | --- | --- | --- | --- | --- | --- |
| P1/P1 | 16 | 6.0e-13 | 6, 27 | 17, 76 | 3.056e-2 | 2.586e-2 |
| P1/P1 | 32 | 7.5e-13 | 4, 16 | 16, 68 | 9.092e-3 | 8.210e-3 |
| P2/P2 | 16 | 1.0e-12 | 2, 8 | 14, 41 | 4.108e-4 | 1.133e-4 |
| P2/P2 | 32 | 4.4e-12 | 2, 5 | 11, 26 | 5.117e-5 | 2.501e-5 |

All Newton/Picard differences are ≤ 1e-8, and both errors drop from mesh 16 to 32 for both pairs. V3 Galerkin is unchanged (Newton/Picard 4.4e-13).

Deviations from spec or plan:
1. h = shortest cell edge instead of the cell diameter, and the U5 eigenvalue test instead of random sampling. These were decided by the spec owner after the first stop and recorded in the spec and plan, which are committed with this phase.
2. New test file `tests/unit/test_asgs_solver.py` (Jacobian finite differences, residual logging, exact-τ failure). It is white-box: it reads private `TimeStepper` and `InnerSolver` members and wraps `InnerSolver._append`.
3. `InnerSolveResult.failure_reason` and the post-Jacobian finiteness callback, needed to name `newton_tau = 'exact'` in the error, and the `PETSc.Error` handling described above.
4. The run log names the case "Sta-Lin-ASGS" for `method = "vms"`, which is correct while it is the only accepted VMS case. Later phases name the case from the switches.
5. V1 for ASGS uses the Galerkin V1 settings; V3 for ASGS uses the Galerkin V3 settings (`inner_max_iterations = 100`).

Open questions: None. For later phases: the stability condition depends on cell shape, and U5 checks it only on the structured test meshes. The phase 10 gmsh cylinder mesh should be checked the same way before P2/P2 runs on it.

## Phase 6: Dynamic subscales (2026-10-22)
Done: Before this phase, commit "Plan: U5 check on the cylinder mesh" added to plan §13 Phase 10 the U5 eigenvalue check on a coarse version of the gmsh cylinder mesh.
- `discretization/function_spaces.py` has `subscale_space(mesh, q)`, a vector space on `quadrature.quadrature_element`, so its values sit at the points of the dx rule.
- `solver/subscale_storage.py` has `SubscaleStorage`:
  - It keeps as many levels as the velocity: 1 for BDF1, 2 for BDF2, 3 for BDF3, 1 for Crank–Nicolson. All start at zero, which gives ũ⁰ = 0 and zero exact-history levels at negative times. Every step shifts them from step 1, including the lower-order startup steps.
  - `form_history(coefficients)` sets ũ_hist = Σ_j history_weights[j] ũ^(n−j) from the step's own `StepCoefficients`.
  - `store(u_sub, coefficients)` runs after convergence. It interpolates the subscale expression at the converged state (the half-level value for Crank–Nicolson), extrapolates it with the step's weights (ũ^(n+1) = 2ũ^(n+½) − ũⁿ for Crank–Nicolson, identity for BDF) and shifts.
- `formulation/stabilization_parameters.py` has `dynamic_parameter` τ_t = (α₀/δt + 1/τ_m)⁻¹, with α₀ the `SolveLevel` Constant that is reset every step.
- `formulation/vms.py` has `velocity_subscale`, ũ = τ_t(R_u + χ_D ũ_hist/δt). `vms_terms` adds χ_D (D_t ũ, v) = (α₀ũ − ũ_hist, v)/δt with ũ substituted. Static subscales pass no history, so τ_t = τ_m and the Phase 5 forms are unchanged.
- `solver/time_stepper.py`: for dynamic subscales the residual uses τ_t(u), and Picard uses τ_t(â).
  - The frozen-τ Newton Jacobian builds τ_t on the frozen advector and then replaces it with u, so it has τ_t(u) without τ_t′. The exact option differentiates τ_t.
  - ũ_hist is formed in `_set_level` before each solve, and ũ^(n+1) is stored in `_advance` after convergence only.
  - The run log names the case from the switches through the new public `config.method_name` ("Galerkin", "Sta-Lin-ASGS", "Dyn-Lin-ASGS", …).
- `config.py`: rule 5 accepts Sta-Lin-ASGS and Dyn-Lin-ASGS.

Tests:
- `ruff check .` and `ruff format --check .` are clean (47 files).
- `pytest tests/unit`: 174 passed in 25 s.
- `pytest tests/verification -m slow`: 23 passed in 476 s.
- `test_subscale_storage.py` (new, the bookkeeping test): all four schemes, both startups, six steps, on the known expression g(t)(sin x, y) with g = t + t² (zero for t ≤ 0).
  - At each step the formed ũ_hist equals Σ_j w_j ũ(t^(n−j)), to 1e-14.
  - The stored solve-level value is the endpoint-average data for Crank–Nicolson and the endpoint value for BDF steps.
  - This checks the shifting, the startup sequence, zero exact-history levels and the Crank–Nicolson extrapolation.
  - Mutation check in scratch: skipping the extrapolation fails Crank–Nicolson, and skipping the shift fails BDF2 and BDF3.
- `test_asgs_solver.py`: the finite-difference Jacobian check now also covers dynamic subscales with a random non-zero ũ_hist, for frozen and exact τ_t, P1 and P2.
- `test_config.py`: rule 5 accepts Dyn-Lin-ASGS (P1/P1, P2/P2, τ_c on and off), and `method_name` gives the spec names.

V2, Dyn-Lin-ASGS (manufactured_temporal, 4×4, T = 1, N = 20 / 40 / 80 / 160). Velocity at T; pressure at the solve level, zero mean. Rates for each refinement; thresholds as for Galerkin:

| Pair | Scheme | ‖u − u_ex‖ | Rates | ‖p − p_ex‖ | Rates |
| --- | --- | --- | --- | --- | --- |
| P1/P1 | BDF1 | 1.593e-3 / 7.986e-4 / 3.998e-4 / 2.001e-4 | 1.00, 1.00, **1.00** | 2.751e-3 / 1.359e-3 / 6.752e-4 / 3.366e-4 | 1.02, 1.01, **1.00** |
| P1/P1 | BDF2 | 2.743e-5 / 7.123e-6 / 1.819e-6 / 4.599e-7 | 1.95, 1.97, **1.98** | 1.339e-4 / 3.391e-5 / 8.530e-6 / 2.139e-6 | 1.98, 1.99, **2.00** |
| P1/P1 | BDF3, exact history | 2.086e-6 / 2.556e-7 / 3.162e-8 / 3.932e-9 | 3.03, 3.02, **3.01** | 3.656e-6 / 4.384e-7 / 5.362e-8 / 6.628e-9 | 3.06, 3.03, **3.02** |
| P1/P1 | Crank–Nicolson | 1.357e-5 / 3.365e-6 / 8.398e-7 / 2.098e-7 | 2.01, 2.00, **2.00** | 1.014e-4 / 2.563e-5 / 6.442e-6 / 1.615e-6 | 1.98, 1.99, **2.00** |
| P1/P1 | BDF3, lower order (logged) | 2.330e-5 / 5.965e-6 / 1.508e-6 / 3.791e-7 | 1.97, 1.98, **1.99** (≥ 1.8) | 1.200e-5 / 2.973e-6 / 7.467e-7 / 1.875e-7 | 2.01, 1.99, **1.99** (≥ 1.8) |
| P2/P2 | BDF1 | 2.890e-3 / 1.457e-3 / 7.320e-4 / 3.669e-4 | 0.99, 0.99, **1.00** | 2.985e-3 / 1.478e-3 / 7.353e-4 / 3.668e-4 | 1.01, 1.01, **1.00** |
| P2/P2 | BDF2 | 9.363e-5 / 2.315e-5 / 5.770e-6 / 1.442e-6 | 2.02, 2.00, **2.00** | 1.332e-4 / 3.365e-5 / 8.461e-6 / 2.122e-6 | 1.98, 1.99, **2.00** |
| P2/P2 | BDF3, exact history | 3.823e-6 / 4.690e-7 / 5.805e-8 / 7.220e-9 | 3.03, 3.01, **3.01** | 3.965e-6 / 4.768e-7 / 5.839e-8 / 7.223e-9 | 3.06, 3.03, **3.02** |
| P2/P2 | Crank–Nicolson | 5.064e-5 / 1.259e-5 / 3.143e-6 / 7.855e-7 | 2.01, 2.00, **2.00** | 1.060e-4 / 2.665e-5 / 6.684e-6 / 1.674e-6 | 1.99, 2.00, **2.00** |
| P2/P2 | BDF3, lower order (logged) | 6.749e-5 / 1.671e-5 / 4.174e-6 / 1.044e-6 | 2.01, 2.00, **2.00** (≥ 1.8) | 2.683e-5 / 6.796e-6 / 1.716e-6 / 4.314e-7 | 1.98, 1.99, **1.99** (≥ 1.8) |

V2 Galerkin is unchanged.

V3, linear ASGS on steady Kovasznay (Re = 40, BDF1 dt = 5 to T = 100, meshes 16 = 24×32 and 32 = 48×64; final relative step change < 1e-10 in every run). The differences are ‖u − u_static,Newton‖/‖u_static,Newton‖ in L² with the common rule:

| Pair | Mesh | Static Picard | Dynamic Newton | Dynamic Picard | Dynamic Newton/Picard |
| --- | --- | --- | --- | --- | --- |
| P1/P1 | 16 | 6.0e-13 | 3.0e-13 | 4.6e-13 | 3.9e-13 |
| P1/P1 | 32 | 7.5e-13 | 3.0e-12 | 3.0e-12 | 5.0e-13 |
| P2/P2 | 16 | 1.0e-12 | 3.7e-12 | 3.6e-12 | 7.9e-13 |
| P2/P2 | 32 | 4.4e-12 | 4.4e-12 | 4.4e-12 | 1.6e-12 |

All are ≤ 1e-8. Static and dynamic steady errors agree to about 10 digits (P1/P1 velocity 3.0560e-2 → 9.0922e-3; P2/P2 velocity 4.1082e-4 → 5.1174e-5). Both drop from mesh 16 to 32.

Deviations from spec or plan:
1. `config.method_name` was made public so the run log can name the case. The internal `_vms_case_name` is kept for the rule 5 message.
2. `vms_terms` takes τ_t (τ_m for static subscales) and an optional `u_sub_history`; None means static (χ_D = 0).
3. New unit test file `tests/unit/test_subscale_storage.py` (the bookkeeping test the acceptance asked for); `test_asgs_solver.py` was extended to dynamic subscales.
4. V2 and V3 tests were restructured to run over variants. V2 runs Galerkin, Dyn-Lin-ASGS P1/P1 and P2/P2 (P2/P2 is cheap on the 4×4 mesh). V3 runs a single ASGS test per degree covering static and dynamic, Newton and Picard, on both meshes. The test names changed accordingly.

Open questions: None.

## Phase 7: Nonlinear scale splitting (2026-10-22)
Done:
- Design change, made first: `formulation/vms.py` has `Switches(is_asgs, is_dynamic, is_nonlinear)` and `switches(stabilization)`.
  - `velocity_subscale` and `vms_terms` receive the switches explicitly, instead of inferring "static" from a missing history.
  - A history given without χ_D = 1, or χ_D = 1 without a history, raises `ValueError`.
- a = u_h + ũ*. ũ* is a Function in the quadrature space `subscale_space(mesh, q)`, frozen during the inner solve.
  - It enters c, L_a in R_u, L*_a and τ_m, τ_c, τ_t. The outflow term keeps u_h.
  - ũ* is a coefficient, never part of the unknown, so no derivative is taken with respect to it.
  - Picard uses â = u^(j−1) + ũ*. The frozen-τ Newton Jacobian builds τ on frozen + ũ* and then replaces frozen with u, as in Phase 6.
- `formulation/vms.py` has `accelerated_subscale`. It gives (I + τ_t∇u_h)⁻¹ τ_t(f − D_t u_h − u_h·∇u_h + νΔ_h u_h − ∇p_h + χ_D ũ_hist/δt), with τ_t from the previous ũ* and `ufl.inv` of the d×d matrix.
- `solver/subscale_storage.py`:
  - `SubscaleIterate` holds ũ* and its update expression, either fixed point ũ* ← τ_t(a)(R_u(a) + χ_D ũ_hist/δt) or accelerated. `update()` interpolates the update into the quadrature space until max|Δũ| ≤ subscale_rtol · max|ũ| over all points.
  - Both maxima are PETSc infinity norms, so they are global over MPI ranks.
  - `SubscaleStorage.latest` returns ũⁿ.
- `solver/nonlinear_solver.py` has `StepSolver`, the outer loop of spec step 5:
  - Each outer iteration runs the inner solve, then the local update. It converges when the relative L² changes of u_h, p_h and ũ* (common rule) are ≤ outer_rtol and the FE residual is ≤ max(inner_atol, outer_rtol · r0), with r0 the initial residual of the step's first inner solve. The limit is outer_max_iterations.
  - Galerkin and linear splitting take one pass, as before.
  - Failures take the existing path: last state written, iteration history logged, `ConvergenceError`, exit status 1. The reasons are the inner solve's own reason, "local subscale update did not converge in N sweeps" and "outer iteration did not converge in N iterations".
  - For static subscales the local-update reason adds "static nonlinear subscales on a coarse mesh can make (I + τ∇u) nearly singular; use dynamic subscales or a finer mesh", as the NaN message names its likely cause.
- `solver/time_stepper.py`:
  - Writes the inner rows of every outer iteration and one real outer_iterations.csv row per outer iteration, with subscale_change filled. time_steps.csv has the outer count and the summed inner count.
  - The step log line gives the outer, inner and local counts.
  - ũ* ← ũⁿ at step 1 of the spec algorithm. For dynamic subscales ũⁿ is the stored level. For static subscales nothing else is stored, so ũⁿ is the previous step's converged iterate, zero at step 1.
  - ũ^(n+1) is stored as in Phase 6: the subscale expression at the converged state, at the half level and then extrapolated for Crank–Nicolson.
- `config.py`: rule 5 accepts Sta-Nl-ASGS and Dyn-Nl-ASGS (all four ASGS cases).
- Spec paragraph "Static nonlinear subscales can be ill-posed on coarse meshes" (after the eight-cases table) is your edit; it is included in this commit, as the Phase 5 spec edit was.
- Plan §10 V1 notes the exception for Sta-Nl-ASGS P1/P1 and cites that paragraph.
- `pyproject.toml`: pytest-xdist added to the dev extras. The README gives the parallel command for the slow suite.

Coarse-mesh finding (V1 Sta-Nl-ASGS P1/P1 on mesh 8, the reason for the plan exception):
- At step 1 the inner Newton solve converges, but the local update does not converge in 50 sweeps, accelerated or fixed point.
- At points where u_h ≈ 0, τ_m is at its viscous limit h²/(c1·k⁴·ν) = (1/8)²/(4 · 0.01) ≈ 0.39. With an eigenvalue of ∇u_h near −2.6, 1 + τ·λ ≈ 0, so (I + τ∇u_h) is nearly singular.
- The update then cycles with period 2: at those points |ũ| jumps between about 0 and about 84, against a converged |ũ| of about 1 elsewhere. After 60 sweeps 77 of 768 points were still changing (scratch probe).
- Under-relaxation by ½, tried in scratch only and not kept, still left a change of 2.4 after 60 sweeps.
- Mesh 16 passes but is near the edge: step 1 takes 25 outer iterations and 264 local sweeps, against 7 and 24 on mesh 64.

Tests:
- `ruff check .` and `ruff format --check .` are clean (47 files).
- `pytest tests/unit`: 197 passed in 40 s.
- `pytest tests/verification -m slow -n 4`: 37 passed, in 10 min 48 s and 11 min 0 s on two runs, with identical numbers (sum of per-test durations about 41 min). The serial run before the V3 trim took 38 min 11 s.
- `test_subscale_storage.py`, fixed point against accelerated (P1, P2; static, dynamic):
  - Data: a velocity bounded away from zero, (1 + sin(3x)y/4, 1 + xy/2), and pressure xy − sin y.
  - History 1.5·velocity, so D_t u_h = 0 with α₀ = 1.5. ũ_hist = 0.01·(sin x, y); ν = 0.01; local rtol 1e-13.
  - Results: both reach the same ũ (≤ 1e-10 relative), and ũ − τ_t(a)(R_u(a) + χ_D ũ_hist/δt) ≤ 1e-13 relative at every point.
  - Sweeps, fixed point / accelerated: P1 static 13 / 13, P1 dynamic 9 / 9, P2 static 10 / 9, P2 dynamic 8 / 8.
- `test_asgs_solver.py`: the finite-difference Jacobian check now also runs with nonlinear splitting, with a random ũ* in [−0.5, 0.5], a = u + ũ* and τ on (frozen or u) + ũ*. It covers frozen and exact τ, static and dynamic, P1 and P2.
- `test_time_stepper.py`:
  - A static and a dynamic Kovasznay run (12×16, P1/P1, BDF2, 2 steps) write at least 2 outer rows per step, with only the last one converged, subscale_change filled, final changes ≤ 1e-8 and consistent totals.
  - `outer_max_iterations = 1` and `subscale_max_iterations = 1` take the failure path with the expected reasons, including the static hint and its absence for dynamic subscales.
- `test_config.py`: rule 5 accepts all four ASGS cases (P1/P1, P2/P2, τ_c on and off), and `method_name` gives their names.
- MPI (scratch, not a test): static and dynamic nonlinear Kovasznay runs (12×16, P1/P1, BDF2, 4 steps) on 2 ranks have the same outer and local counts as serial at every step, and final errors equal to about 16 digits. Inner counts differ by up to 3 per step through rounding of the parallel direct solve near the tolerance.

V1, Sta-Nl-ASGS (manufactured_spatial, two BDF1 steps). Errors on three meshes, then rates for each refinement:

| Pair | Meshes | ‖u − u_ex‖ | Rates | ‖∇(u − u_ex)‖ | Rates | ‖p − p_ex‖ | Rates |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P1/P1 | 16, 32, 64 | 6.056e-2 / 1.327e-2 / 3.173e-3 | 2.19, **2.06** | 2.473 / 1.231 / 6.155e-1 | 1.01, **1.00** | 1.344e-1 / 2.111e-2 / 2.141e-3 | 2.67, **3.30** |
| P2/P2 | 8, 16, 32 | 1.187e-2 / 1.487e-3 / 1.857e-4 | 3.00, **3.00** | 7.318e-1 / 1.819e-1 / 4.474e-2 | 2.01, **2.02** | 8.273e-3 / 9.470e-4 / 1.767e-4 | 3.13, **2.42** |

V2, Dyn-Nl-ASGS (manufactured_temporal, 4×4, T = 1, N = 20 / 40 / 80 / 160). Velocity at T; pressure at the solve level, zero mean:

| Pair | Scheme | ‖u − u_ex‖ | Rates | ‖p − p_ex‖ | Rates |
| --- | --- | --- | --- | --- | --- |
| P1/P1 | BDF1 | 1.574e-3 / 7.895e-4 / 3.953e-4 / 1.978e-4 | 1.00, 1.00, **1.00** | 2.787e-3 / 1.377e-3 / 6.843e-4 / 3.411e-4 | 1.02, 1.01, **1.00** |
| P1/P1 | BDF2 | 2.739e-5 / 7.118e-6 / 1.818e-6 / 4.597e-7 | 1.94, 1.97, **1.98** | 1.343e-4 / 3.403e-5 / 8.561e-6 / 2.147e-6 | 1.98, 1.99, **2.00** |
| P1/P1 | BDF3, exact history | 2.063e-6 / 2.527e-7 / 3.126e-8 / 3.887e-9 | 3.03, 3.02, **3.01** | 3.703e-6 / 4.442e-7 / 5.434e-8 / 6.717e-9 | 3.06, 3.03, **3.02** |
| P1/P1 | Crank–Nicolson | 1.340e-5 / 3.326e-6 / 8.299e-7 / 2.074e-7 | 2.01, 2.00, **2.00** | 1.016e-4 / 2.566e-5 / 6.450e-6 / 1.617e-6 | 1.98, 1.99, **2.00** |
| P1/P1 | BDF3, lower order (logged) | 2.288e-5 / 5.855e-6 / 1.480e-6 / 3.720e-7 | 1.97, 1.98, **1.99** (≥ 1.8) | 1.245e-5 / 3.104e-6 / 7.814e-7 / 1.964e-7 | 2.00, 1.99, **1.99** (≥ 1.8) |
| P2/P2 | BDF1 | 2.902e-3 / 1.464e-3 / 7.352e-4 / 3.685e-4 | 0.99, 0.99, **1.00** | 2.985e-3 / 1.478e-3 / 7.354e-4 / 3.668e-4 | 1.01, 1.01, **1.00** |
| P2/P2 | BDF2 | 9.417e-5 / 2.327e-5 / 5.797e-6 / 1.448e-6 | 2.02, 2.01, **2.00** | 1.332e-4 / 3.367e-5 / 8.465e-6 / 2.123e-6 | 1.98, 1.99, **2.00** |
| P2/P2 | BDF3, exact history | 3.840e-6 / 4.710e-7 / 5.831e-8 / 7.252e-9 | 3.03, 3.01, **3.01** | 3.966e-6 / 4.768e-7 / 5.839e-8 / 7.223e-9 | 3.06, 3.03, **3.02** |
| P2/P2 | Crank–Nicolson | 5.116e-5 / 1.272e-5 / 3.175e-6 / 7.935e-7 | 2.01, 2.00, **2.00** | 1.062e-4 / 2.671e-5 / 6.701e-6 / 1.678e-6 | 1.99, 2.00, **2.00** |
| P2/P2 | BDF3, lower order (logged) | 6.831e-5 / 1.692e-5 / 4.226e-6 / 1.056e-6 | 2.01, 2.00, **2.00** (≥ 1.8) | 2.726e-5 / 6.904e-6 / 1.743e-6 / 4.381e-7 | 1.98, 1.99, **1.99** (≥ 1.8) |

V3, nonlinear ASGS on steady Kovasznay (Re = 40, BDF1 dt = 5 to T = 100; final relative step change < 1e-10 in every run). The equivalences are checked on mesh 16 (24×32). Differences are ‖u − u_static,Newton,accelerated‖/‖u_static,Newton,accelerated‖ in L² with the common rule:

| Pair | Static Picard | Dynamic Newton | Dynamic Picard | Static fixed point | Dynamic Newton/Picard |
| --- | --- | --- | --- | --- | --- |
| P1/P1 | 1.2e-13 | 2.9e-13 | 3.9e-13 | 8.7e-17 | 2.7e-13 |
| P2/P2 | 6.3e-13 | 3.9e-13 | 8.2e-13 | 1.5e-16 | 7.9e-13 |

All are ≤ 1e-8. The linear cases on mesh 16 give the Phase 6 values again: at most 3.7e-12.

Iteration counts over the 20 steps:

| Run | Outer total (max per step) | Inner total | Local sweeps total (max per step) |
| --- | --- | --- | --- |
| P1/P1, mesh 16, static Newton | 54 (9) | 93 | 152 (45) |
| P1/P1, mesh 16, static Picard | 54 (9) | 261 | 152 (45) |
| P1/P1, mesh 16, dynamic Newton | 52 (9) | 89 | 142 (45) |
| P1/P1, mesh 16, dynamic Picard | 52 (9) | 247 | 142 (45) |
| P1/P1, mesh 16, static Newton, fixed point | 54 (9) | 92 | 182 (56) |
| P1/P1, mesh 32, static Newton | 42 (8) | 50 | 96 (30) |
| P1/P1, mesh 32, dynamic Newton | 42 (8) | 50 | 93 (30) |
| P2/P2, mesh 16, static Newton | 32 (7) | 24 | 57 (20) |
| P2/P2, mesh 16, static Picard | 32 (7) | 100 | 57 (20) |
| P2/P2, mesh 16, dynamic Newton | 34 (7) | 25 | 60 (20) |
| P2/P2, mesh 16, dynamic Picard | 34 (7) | 107 | 60 (20) |
| P2/P2, mesh 16, static Newton, fixed point | 32 (7) | 24 | 71 (27) |
| P2/P2, mesh 32, static Newton | 29 (6) | 13 | 44 (15) |
| P2/P2, mesh 32, dynamic Newton | 29 (6) | 15 | 46 (15) |

Runs are accelerated unless marked fixed point. Steady errors, static and dynamic agreeing to about 10 digits, drop from mesh 16 to 32 in every case:

| Pair | Velocity | Pressure |
| --- | --- | --- |
| P1/P1 | 3.2942e-2 → 9.7145e-3 | 1.8403e-2 → 5.8394e-3 |
| P2/P2 | 4.0945e-4 → 5.1126e-5 | 1.0847e-4 → 2.5020e-5 |

Deviations from spec or plan:
1. Explicit switches, as you asked: `vms_terms` and `velocity_subscale` take a `Switches` value and the history, instead of Phase 6's "None means static".
   - χ_A gates the (D_t ũ, v) term. Its other role, dropping D_t u_h from R_u for OSS, comes with R̂_u and ξ_h in Phase 8.
2. FE-residual part of the outer test: ‖r‖ ≤ max(inner_atol, outer_rtol · r0), with r0 the initial residual of the step's first inner solve. The spec says only "below tolerance".
3. Local-update test: max-norm change ≤ subscale_rtol · max|ũ| over all quadrature points. The subscale_change reported in outer_iterations.csv is the relative L² change, like the u and p changes.
4. Every inner solve iterates to its tolerance. The single Picard iteration per outer iteration that spec step 2 allows is not offered (your decision 4).
5. The coarse-mesh hint in the local-update failure appears for static subscales only, because its advice is to use dynamic subscales.
6. Test meshes, accepted by you:
   - The nonlinear run tests use Kovasznay 12×16, because 6×8 hits the same local-update cycling.
   - The local-update unit test uses a velocity bounded away from zero and consistent history data, for the same reason.
7. V1 Sta-Nl-ASGS P1/P1 uses meshes 16, 32, 64, with rates on 32 → 64 (plan §10 V1 updated, your decision 1).
8. V3 runtime (your decision 3):
   - static = dynamic, Newton = Picard and fixed point = accelerated are checked on mesh 16 only, for every ASGS case.
   - Mesh 32 runs only static and dynamic Newton for the error-drop check, so each case keeps its 16 → 32 drop.
   - fixed point = accelerated is checked for static Newton runs.
   - No variant or time scheme was dropped.
9. The slow suite is meant to run in parallel with pytest-xdist: `pytest tests/verification -m slow -n 4`. Four workers were chosen for the 15 GB of memory on this machine. The ten V2 Dyn-Nl-ASGS studies take 100–235 s each and set the wall time.

Open questions: None.

## Phase 7b: Review fixes (2026-10-22)
Done: Fixes for findings R1–R5 of the independent review `docs/reviews/review_phases_0-7.md`, committed with this phase; R6 is Phase 8. The uncommitted lagged Phase 8 work was stashed first, so this phase starts from dad7b55.
- R1, `solver/nonlinear_solver.py`:
  - A Picard linear solve that raises `ConvergenceError` sets the reason "Picard iteration N: <error>", and `InnerSolver.solve` returns failure at once, through the existing failure path. Before, the loop retried without recording a row and never ended.
  - A non-finite residual after a Picard iteration also fails at once, with the existing non-finite message.
  - The Picard warm-up runs at most min(picard_iterations_before_newton, inner_max_iterations) iterations, and Newton gets what is left of inner_max_iterations.
- R3, V3 (`tests/verification/test_steady_equivalence.py`): `march` records the endpoint velocity, and the stored subscale for dynamic cases, before the last step. Steady means both changed by less than 1e-10 (relative L², common rule) in that step. That every step converged is checked separately from time_steps.csv.
- R4, `solver/time_stepper.py`: with the pressure null space active, p ← p − ∫p dx/|Ω| (common rule) right after each step's solve. This comes before the diagnostics, the VTK output, the failure path's last state and the next step's initial guess. Outflow problems are untouched.
- R2 and R5: new and extended tests, below.
- Plan: §10 U2, U4, U8 (new) and V3; §11 pressure null space; §13 Phase 7b.
- README: OMP_NUM_THREADS=1 in the parallel slow-suite command, as you asked after Phase 7.

Tests:
- `ruff check .` and `ruff format --check .` are clean (49 files).
- `pytest tests/unit`: 222 passed, in 70 s serially (while the slow suite ran) and 23 s with `-n 4`.
- `OMP_NUM_THREADS=1 pytest tests/verification -m slow -n 4`: 37 passed in 11 min 56 s, while the unit and mutation runs shared the machine. The V1 and V2 tables are identical to Phase 7.
- R1, `test_time_stepper.py`:
  - `test_a_failed_picard_linear_solve_fails_the_step_at_once`, for Picard and for Newton with a Picard warm-up: a patched Picard linear solve raises `ConvergenceError`.
    - There is exactly one attempt, inner_iterations.csv holds only row 0, the last state is written, and the log reads "Step 1: Picard iteration 1: synthetic linear solver failure".
    - A guard raises after three attempts, so the old endless retry fails the test instead of hanging it.
  - `test_the_picard_warm_up_counts_toward_the_inner_iteration_limit`: a warm-up of 5 with inner_max_iterations = 2 gives exactly 2 Picard rows. The old code gave 5.
- R4, `test_pressure_has_zero_spatial_mean_exactly_with_the_null_space` (Galerkin P2/P1 and static linear ASGS P1/P1, two BDF1 steps on manufactured_spatial 4×4):
  - With the null space, |∫p dx|/|Ω| ≤ 1e-13 ‖p‖.
  - With x_max made an outflow, the mean stays above 1e-6 ‖p‖, so it was not touched. The review measured a mean of 0.0368 before the fix.
- R2(a), U4 (`test_dynamic_asgs_step_equals_independent_classical_form`, 8 cases: BDF2 and Crank–Nicolson; P1/P1 and P2/P2; linear and nonlinear splitting):
  - What is compared: the time stepper's own Newton residual at step 2, after `_set_level`, with random velocity levels, stored subscale levels and state, against a classical reference written in the test. The nonlinear cases use a = u + ũⁿ.
  - The reference's coefficients are literals: BDF2 α₀ = 3/2 with u_hist = 2uⁿ − ½uⁿ⁻¹ and ũ_hist = 2ũⁿ − ½ũⁿ⁻¹; Crank–Nicolson α₀ = 2 with 2uⁿ, 2ũⁿ and averaged forcing. The reference also has τ_t = (α₀/δt + 1/τ_m)⁻¹, the closure ũ = τ_t(R_u + ũ_hist/δt) and (D_t ũ, v) = ((α₀ũ − ũ_hist)/δt, v).
  - The two agree to 1e-12 relative.
- R2(b), U8 (`tests/unit/test_subscale_time_integration.py`, new):
  - What it runs: the time stepper's `_set_level` and `_advance`, without solving, with u_h = (0.2, 0) and p_h = xy (P1) prescribed at every step on manufactured_spatial 2×2.
  - Then τ_m = (4ν/h² + 2|U|/h)⁻¹ = 1/0.96 is the same at every point, D_t u_h = 0, and the forcing g(t)F(x) + ∇p_h gives R = g(t)F(x).
  - Exact solution: ũ = y(t)F, with y′ + y/τ = g and y(0) = 0. The error is max |ũ^N − y(1)F| / max |y(1)F| at the quadrature points, for N = 10, 20, 40, 80.
  - Data "generic" is g = cos 2t. Data "smooth start" is y = sin³ t, which vanishes with y′ and y″ at t = 0.

    | Scheme | Startup | Data | Errors | Rates |
    | --- | --- | --- | --- | --- |
    | BDF1 | lower order | generic | 2.388e-1 / 1.212e-1 / 6.105e-2 / 3.064e-2 | 0.98, 0.99, **0.99** |
    | BDF2 | lower order | generic | 1.579e-2 / 3.290e-3 / 7.391e-4 / 1.743e-4 | 2.26, 2.15, **2.08** |
    | BDF3 | lower order | generic | 1.164e-2 / 2.799e-3 / 6.916e-4 / 1.722e-4 | 2.06, 2.02, **2.01** |
    | BDF3 | exact history | generic | 8.732e-2 / 4.468e-2 / 2.243e-2 / 1.122e-2 | 0.97, 0.99, **1.00** |
    | BDF3 | exact history | smooth start | 2.551e-3 / 3.267e-4 / 4.072e-5 / 5.063e-6 | 2.96, 3.00, **3.01** |
    | Crank–Nicolson | lower order | generic | 9.355e-3 / 2.170e-3 / 5.213e-4 / 1.277e-4 | 2.11, 2.06, **2.03** |

  - Each rate is asserted within 0.2 of the expected order.
  - Expected BDF3 order with exact history (asked for): 1 for generic data. The stored subscale is zero at negative times, while the smooth extension of y has y(−jδt) ≈ −j y′(0) δt. The first step's history is therefore wrong by O(δt), so ũ¹ carries an O(δt) error, about (5/11) y′(0) δt, and the stable recursion carries it to T.
  - The order is 3 when y vanishes to third order at t = 0. That case is included, so BDF3's subscale weights are also checked at third order.
- R2(c), the review's mutations, in scratch git worktrees of HEAD plus the Phase 7b changes, unit suite with `-n 4`:

  | Mutation | Unit result | Failing tests |
  | --- | --- | --- |
  | Backward Euler for the subscale under BDF2 (the review's `bdf2_be.patch`: α₀ = 1 in τ_t, ũ_hist = ũⁿ and α₀ = 1 in (D_t ũ, v)) | 8 failed, 214 passed | U4 dynamic BDF2, all 4 cases (relative residual difference 0.22–0.32 against 1e-12); U8 BDF2 (rate 0.99 instead of 2); the storage-history tests bdf2-lower_order, bdf2-exact_history and bdf3-lower_order, which already failed in the review |
  | No Crank–Nicolson extrapolation of ũ (`cn_subscale.patch`) | 3 failed, 219 passed | U8 Crank–Nicolson (rate −0.09); the storage-history tests for Crank–Nicolson with both startups (already failing in the review) |

- R5, U2 (`test_production_subscale_space_samples_the_production_dx`, 2D and 3D, k = 1, 2):
  - What it checks: a non-polynomial vector field sampled into `subscale_space(mesh, 3k)` and tested against a non-polynomial weight with `volume_measure(mesh, 3k)` integrates to the direct value to 1e-14 relative. A rule of degree 3k + 4 differs.
  - The review's q + 2 mutation now fails this test as well, because the mismatched quadrature points raise an error.
- R3, V3: the change over the last step is 0.0 for the velocity in every run, because the inner solve is converged at its initial residual. For the stored subscale it is ≤ 4.6e-15.
  - A one-step march (Dyn-Nl-ASGS P1/P1, 24×32, δt = 5) is now rejected, with a velocity change of 9.374e-3, the review's number, and a subscale change of 1.0.
  - Under the old check it passed with a last-outer change of 6.9e-11.

Deviations from spec or plan:
1. R1: the failure reason carries a "Picard iteration N:" prefix. A non-finite Picard residual is now a failure too, matching the Newton path.
2. U8 drives the time stepper's private `_set_level` and `_advance` with prescribed fields, so it covers the production subscale path apart from the FE solve. The smooth-start case is an addition, so that BDF3's subscale coefficients are checked at third order.
3. U4 dynamic uses step 2 with `lower_order` startup, where BDF2 and Crank–Nicolson are the step's own schemes. The velocity and subscale levels are set directly instead of being marched.
4. R4 uses the common rule for ∫p dx and |Ω|. The normalization also applies to the state written by the failure path.
5. V3's 1e-10 tolerance applies to the relative change of the stored subscale as well as the velocity.

Open questions:
1. U8 shows that `exact_history` startup gives first-order dynamic subscales in general, because the subscale history at negative times is zero (spec: Initial subscale). V2's BDF3 exact-history rates are third order only because the subscale there comes from the discretization error alone. Should exact-history startup also fill the subscale levels, where a case defines an exact subscale? Nothing in Phase 8 depends on this.

## Phase 8: Orthogonal subscales (2026-10-22)
Done:
- First attempt, lagged projections (spec step 4 as first written). I implemented it, measured it, stopped and reported; it was replaced by coupled projections after your decision and the independent review (finding R6). The evidence, from a frozen-advector power-iteration probe of the map (ξ_h, ζ_h) → solve → project on Kovasznay:

  | Pair | Spectral radius of the lagged map |
  | --- | --- |
  | P1/P1 | 0.83–0.87 for 12×16 and 24×32, δt from 0.005 to 5, static or dynamic, with or without ζ_h |
  | P2/P2, static | 1.09–1.60 (divergent) |
  | P2/P2, dynamic | 0.985 |

  - Projecting the full residual, as Colomés does, gave the same radii (P1 0.852, P2 static 1.165).
  - The review confirmed the radii independently, including −1.186 for static P2/P2 on 12×16, and found that half relaxation leaves 0.993.
  - The real V-test configurations matched: V1 P1/P1 mesh 8 took 238 and 210 outer iterations; V2 Dyn-Lin-OSS P1/P1 took 75–87 per step; V1 and V2 P2/P2 did not converge in 300; V3 static P2/P2 diverged (FE residual about 1e87).
- `formulation/oss_projections.py` has `projection_equations`: (τ_t(ξ_h − X), w_h) + (τ_c(ζ_h − ∇·u_h), r_h), with X = R̂_u + χ_D ũ_hist/δt. It uses consistent mass with the common rule. ζ_h and its equation are omitted when τ_c is off.
- `formulation/vms.py`:
  - `subscale_residual` gives X, with R̂_u = f − L_a(u_h, p_h) for OSS (χ_A = 0), so D_t u_h is dropped from both ũ and ξ_h.
  - `velocity_subscale` and `accelerated_subscale` take ξ_h. `vms_terms` takes ξ_h and ζ_h and gives τ_c(∇·u − ζ_h, ∇·v); OSS has no (D_t ũ, v) term.
  - ξ_h must be given exactly for OSS, and ζ_h only for OSS.
- `discretization/function_spaces.py`:
  - `oss_space` builds V × Q × V̂ (× Q̂), where V̂ and Q̂ have the elements of V and Q. They stay unconstrained because no Dirichlet condition or null space is ever applied to components 2 and 3.
  - `constant_pressure_null_space` puts the constant on p_h only, for any number of components.
- `solver/time_stepper.py`: for OSS the unknown is (u_h, p_h, ξ_h, ζ_h), and the residual adds the projection equations.
  - Newton uses τ(a) in the projection weights; frozen-τ Newton builds the weights on the frozen advector before the replacement, as elsewhere.
  - Picard uses τ̂ and X(â).
  - ξ_h and ζ_h start at zero at step 1 (the spec allows it) and then carry over from the previous step's solution. Under Crank–Nicolson they are part of the half-level solution.
- `solver/nonlinear_solver.py`:
  - The outer loop exists only when a ũ* iterate exists, so linear OSS takes one pass per step and nonlinear OSS iterates on ũ* alone.
  - xi_change and zeta_change come from the solution's components. For a linear OSS row they are the change over the step; for nonlinear OSS they enter the convergence test.
- `config.py`: rule 5 accepts all eight cases. Its check and the phase table are removed, because nothing is left to reject.
- Profiling finding (V2 Dyn-Nl-OSS P1/P1, N = 160, 641 outer iterations, 81 s):
  - Of 155 s under the profiler, about 100 s went to UFL preprocessing on every `fd.assemble` and `Function.interpolate` call: `residual_norm` 73 s (1,282 calls), the local-update and storage interpolations 43 s (1,443 calls), `relative_change` 16 s (3,205 calls).
  - `solver/reusable_assembly.py` has `reusable_assembler`: `fd.assemble`'s own preprocessing (`BaseFormAssembler.preprocess_base_form`), done once, then `get_assembler(..., is_base_form_preprocessed=True)`. It is used by `InnerSolver.residual_norm`, the new `RelativeChange` (the outer changes), `SubscaleIterate.update` and `SubscaleStorage.store`.
  - A plain `get_assembler` or `get_interpolator` on the interpolation failed inside pyop2 on the OSS subscale expression; only fd.assemble's own preprocessing sequence works.
  - Result: that run fell from 81 s to 17.6 s, and Dyn-Nl-ASGS from 84 s to 20.9 s. diagnostics.csv and outer_iterations.csv are byte-identical to before.
- README: the pinned Firedrake version (2026.10.0, PETSc 3.26, UFL 2026.2), why it is pinned, and the slow-suite time.
- Plan: §1 scope, §3 module tree, §5 module table, §10 U5, U6, U9 (new), V1 and V4, §11 the reusable-assembler note, §13 Phase 8.
- Spec: your edits for coupled projections (fixed-choices row "OSS projection", the eight-cases paragraph, OSS projections, Inner linearization, Time-step algorithm and its step table) are committed with this phase.

Tests:
- `ruff check .` and `ruff format --check .` are clean (53 files).
- `pytest tests/unit -n 4`: 279 passed in 27 s.
- `OMP_NUM_THREADS=1 pytest tests/verification -m slow -n 4`: 61 passed in 9 min 11 s (33 min 36 s of CPU time). The previous run, identical apart from V4, took 8 min 57 s.
- Coupled against fully converged lagged (scratch; Kovasznay 12×16, P1/P1, BDF2, δt = 0.5, 2 steps; lagged code restored from the stash into a scratch worktree, with outer_rtol 1e-13 and inner_rtol 1e-12). Relative ℓ² differences of the final fields:

  | Case | Velocity | Pressure | ξ_h | ζ_h | Lagged outer iterations per step | Coupled outer iterations per step |
  | --- | --- | --- | --- | --- | --- | --- |
  | Sta-Lin-OSS | 3.6e-13 | 2.6e-11 | 1.9e-10 | 1.5e-10 | 134, 116 | 1, 1 |
  | Sta-Nl-OSS | 3.1e-13 | 2.4e-11 | 1.8e-10 | 1.2e-10 | 141, 121 | 17, 15 |
  | Dyn-Lin-OSS | 4.2e-13 | 3.4e-11 | 2.5e-10 | 1.7e-10 | 133, 119 | 1, 1 |
  | Dyn-Nl-OSS | 3.0e-13 | 2.7e-11 | 2.0e-10 | 1.2e-10 | 139, 125 | 15, 14 |

- An independent steady monolithic probe (scratch) gave static linear OSS errors on 24×32 → 48×64 of P1 1.80e-2 → 4.51e-3 and P2 4.1226e-4 → 5.1235e-5. V3's P2/P2 steady errors agree to four digits.
- System size and memory (V3 Dyn-Nl-OSS P2/P2, Kovasznay 48×64, the largest coupled case):

  | Case | Unknowns | Max RSS | Wall time |
  | --- | --- | --- | --- |
  | Coupled OSS | 75,078 | 754 MB | 15.2 s |
  | ASGS | 37,539 | 466 MB | 7.2 s |

  On V2's 4×4 mesh the P1/P1 coupled system has 150 unknowns, against 75.
- U6, `tests/unit/test_oss_projections.py` (new):
  - Identity: random FE w_h and r_h with values on every node, boundary included, and pointwise-variable τ_t and τ_c from a random advector; the projection returns them to 1e-12, P1/P2, 2D/3D.
  - Orthogonality after the coupled solve, at step 2 of BDF2 (dynamic ũ_hist non-zero, checked), Kovasznay 12×16: |(ũ, w_h)|/|(τ_t X, w_h)| over all basis functions of V̂_h, for the subscale expression and the stored ũ.
    - Linear splitting: ≤ 5.5e-16 for P1, ≤ 2.1e-14 for P2.
    - Nonlinear splitting: 1.7e-11 to 2.5e-11 for P1 and 1.1e-13 for P2. These are higher because the final local update moves ũ*, inside a and τ_t, slightly after the solve.
- U5, static OSS (`test_static_oss_energy_has_no_negative_eigenvalues`, 12 cases): ξ_h and ζ_h are eliminated by the Schur complement of the coupled bilinear form, which is the column-by-column energy matrix because the projections are linear in (u, p). The worst λ_min/λ_max is −1.4e-16 (P2/P2 crossed, ν = 10: −2.2e-14 against 156), against the bound −1e-12.
- Finite-difference Jacobian (`test_asgs_solver.py`): now 32 cases, ASGS and coupled OSS × degree × dynamics × splitting × τ option, all within 1e-8. This covers frozen τ in the projection weights.
- U9, `tests/unit/test_reusable_assembly.py` (new), Dyn-Nl-ASGS and Dyn-Nl-OSS P2/P2:
  - Compared: the coupled residual, every outer L² change (ξ_h and ζ_h included), one fixed-point sweep of `SubscaleIterate.update` and `SubscaleStorage.store`.
  - All are byte-identical to `fd.assemble` and `Function.interpolate`, with every input redrawn between two checks.
- `test_time_stepper.py`:
  - Nonlinear OSS records xi_change and zeta_change in every outer row, both ≤ 1e-8 at convergence; for ASGS they are empty.
  - `test_linear_oss_solves_each_step_in_one_pass` (static and dynamic): one converged outer row per step.
  - Test ids no longer embed config text.
- `test_config.py`: rule 5 accepts all eight cases × degree × pressure_subscale; the rejection tests are removed.
- MPI (scratch): Dyn-Lin-OSS and Dyn-Nl-OSS Kovasznay 12×16, P1/P1, 4 BDF2 steps, on 2 ranks give the same outer and local counts as serial. One inner count differs by 1, and the final errors agree to 15–16 digits.

V1, static OSS (manufactured_spatial, two BDF1 steps). Errors on three meshes, then rates:

| Case | Pair | Meshes | ‖u − u_ex‖ | Rates | ‖∇(u − u_ex)‖ | Rates | ‖p − p_ex‖ | Rates |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Sta-Lin-OSS | P1/P1 | 8, 16, 32 | 2.525e-1 / 5.811e-2 / 1.250e-2 | 2.12, **2.22** | 5.031 / 2.486 / 1.233 | 1.02, **1.01** | 7.075e-1 / 1.727e-1 / 3.716e-2 | 2.03, **2.22** |
| Sta-Lin-OSS | P2/P2 | 8, 16, 32 | 1.290e-2 / 1.551e-3 / 1.880e-4 | 3.06, **3.04** | 7.673e-1 / 1.858e-1 / 4.500e-2 | 2.05, **2.05** | 1.663e-2 / 1.655e-3 / 2.336e-4 | 3.33, **2.82** |
| Sta-Nl-OSS | P1/P1 | 16, 32, 64 | 5.977e-2 / 1.267e-2 / 2.875e-3 | 2.24, **2.14** | 2.493 / 1.233 / 6.158e-1 | 1.02, **1.00** | 1.770e-1 / 3.755e-2 / 8.428e-3 | 2.24, **2.16** |
| Sta-Nl-OSS | P2/P2 | 8, 16, 32 | 1.245e-2 / 1.510e-3 / 1.864e-4 | 3.04, **3.02** | 7.503e-1 / 1.815e-1 / 4.461e-2 | 2.05, **2.02** | 1.748e-2 / 1.794e-3 / 2.475e-4 | 3.28, **2.86** |

Sta-Nl-OSS P1/P1 fails on mesh 8 with "local subscale update did not converge in 50 sweeps" at step 1, so it takes the same exception as Sta-Nl-ASGS (spec: The eight cases, the paragraph on static nonlinear subscales). On mesh 16, step 1 takes 22 outer iterations and 137 local sweeps.

V2, dynamic OSS (manufactured_temporal, 4×4, T = 1, N = 20 / 40 / 80 / 160). Rates for each refinement, velocity at T and pressure at the solve level:

| Case | Pair | Scheme | Velocity rates | Pressure rates | Finest ‖u − u_ex‖ | Finest ‖p − p_ex‖ |
| --- | --- | --- | --- | --- | --- | --- |
| Dyn-Lin-OSS | P1/P1 | BDF1 | 0.99, 1.00, **1.00** | 1.02, 1.01, **1.00** | 2.701e-4 | 3.530e-4 |
| Dyn-Lin-OSS | P1/P1 | BDF2 | 1.87, 1.91, **1.95** | 1.98, 1.99, **2.00** | 2.574e-7 | 2.215e-6 |
| Dyn-Lin-OSS | P1/P1 | BDF3, exact history | 3.03, 3.01, **3.01** | 3.06, 3.03, **3.02** | 5.308e-9 | 6.951e-9 |
| Dyn-Lin-OSS | P1/P1 | Crank–Nicolson | 2.01, 2.00, **2.00** | 1.98, 1.99, **2.00** | 2.648e-7 | 1.610e-6 |
| Dyn-Lin-OSS | P1/P1 | BDF3, lower order | 1.96, 1.98, **1.99** | 2.01, 1.99, **1.99** | 6.052e-7 | 1.815e-7 |
| Dyn-Nl-OSS | P1/P1 | BDF1 | 0.99, 1.00, **1.00** | 1.02, 1.01, **1.00** | 2.718e-4 | 3.543e-4 |
| Dyn-Nl-OSS | P1/P1 | BDF2 | 1.89, 1.92, **1.95** | 1.98, 1.99, **2.00** | 2.529e-7 | 2.207e-6 |
| Dyn-Nl-OSS | P1/P1 | BDF3, exact history | 3.03, 3.01, **3.01** | 3.06, 3.03, **3.02** | 5.342e-9 | 6.977e-9 |
| Dyn-Nl-OSS | P1/P1 | Crank–Nicolson | 2.01, 2.00, **2.00** | 1.98, 1.99, **2.00** | 2.758e-7 | 1.613e-6 |
| Dyn-Nl-OSS | P1/P1 | BDF3, lower order | 1.96, 1.98, **1.99** | 2.00, 1.99, **1.99** | 6.174e-7 | 1.921e-7 |
| Dyn-Nl-OSS | P2/P2 | BDF1 | 0.99, 0.99, **1.00** | 1.01, 1.01, **1.00** | 3.965e-4 | 3.698e-4 |
| Dyn-Nl-OSS | P2/P2 | BDF2 | 2.04, 2.01, **2.01** | 1.98, 1.99, **2.00** | 1.519e-6 | 2.139e-6 |
| Dyn-Nl-OSS | P2/P2 | BDF3, exact history | 3.02, 3.01, **3.01** | 3.06, 3.03, **3.02** | 7.803e-9 | 7.283e-9 |
| Dyn-Nl-OSS | P2/P2 | Crank–Nicolson | 2.01, 2.00, **2.00** | 1.99, 2.00, **2.00** | 9.374e-7 | 1.693e-6 |
| Dyn-Nl-OSS | P2/P2 | BDF3, lower order | 2.00, 2.00, **2.00** | 1.98, 1.99, **1.99** | 1.310e-6 | 4.730e-7 |

V3, OSS on steady Kovasznay (Re = 40, BDF1 δt = 5 to T = 100; every run steady between its last two time levels and every step converged). Differences on mesh 16 (24×32) relative to static Newton, accelerated:

| Splitting | Pair | Static Picard | Dynamic Newton | Dynamic Picard | Static fixed point | Dynamic Newton/Picard |
| --- | --- | --- | --- | --- | --- | --- |
| Linear | P1/P1 | 1.5e-12 | 2.3e-14 | 1.5e-12 | – | 1.5e-12 |
| Linear | P2/P2 | 9.9e-13 | 8.4e-16 | 9.8e-13 | – | 9.8e-13 |
| Nonlinear | P1/P1 | 2.0e-13 | 2.3e-14 | 2.0e-13 | 1.0e-16 | 2.0e-13 |
| Nonlinear | P2/P2 | 1.1e-12 | 8.7e-16 | 1.1e-12 | 1.8e-16 | 1.1e-12 |

All are ≤ 1e-8. Steady errors, static and dynamic agreeing to at least 12 digits, drop from mesh 16 to 32 in every case:

| Splitting | Pair | Velocity | Pressure |
| --- | --- | --- | --- |
| Linear | P1/P1 | 1.8569e-2 → 4.6640e-3 | 1.4897e-2 → 3.6908e-3 |
| Linear | P2/P2 | 4.1227e-4 → 5.1236e-5 | 1.2571e-3 → 1.4438e-4 |
| Nonlinear | P1/P1 | 1.8516e-2 → 4.6565e-3 | 1.4902e-2 → 3.6889e-3 |
| Nonlinear | P2/P2 | 4.1039e-4 → 5.1174e-5 | 1.2472e-3 → 1.4365e-4 |

Outer iteration counts over the 20 steps (total, with the most in one step):

| Run | Outer | Inner | Local sweeps |
| --- | --- | --- | --- |
| Linear OSS, any pair, mesh, dynamics or linearization | 20 (1) | P1 mesh 16: Newton 16, Picard 75–77; P2 mesh 16: Newton 9–14, Picard 44; mesh 32 Newton: P1 12, P2 6 | – |
| Nonlinear OSS P1/P1, mesh 16, static or dynamic, Newton | 49 (9) | 59 | 113–116 (34) |
| Nonlinear OSS P1/P1, mesh 16, Picard | 49 (9) | 217–223 | 114–117 (34) |
| Nonlinear OSS P1/P1, mesh 16, static Newton, fixed point | 49 (9) | 59 | 151 (52) |
| Nonlinear OSS P1/P1, mesh 32, Newton | 40 (8) | 35 | 81–82 (25) |
| Nonlinear OSS P2/P2, mesh 16, Newton | 33 (7) | 24 | 58 (19) |
| Nonlinear OSS P2/P2, mesh 16, Picard | 33 (7) | 92 | 58 (19) |
| Nonlinear OSS P2/P2, mesh 16, static Newton, fixed point | 33 (7) | 24 | 73 (27) |
| Nonlinear OSS P2/P2, mesh 32, Newton | 29 (6) | 14–15 | 44–46 (15) |

V4, Taylor–Green (periodic (0, 2π)², ν = 0.01, BDF2 δt = 0.025 to T = 1; exact kinetic energy 9.4826):

| Case | Mesh | ‖u − u_ex‖ | Kinetic energy | Relative energy error |
| --- | --- | --- | --- | --- |
| Galerkin P2/P1 | 16 | 7.826e-2 | 9.4710 | 1.2e-3 |
| Galerkin P2/P1 | 32 | 4.720e-3 | 9.4822 | 4.7e-5 |
| Dyn-Nl-OSS P1/P1 | 16 | 2.097e-1 | 8.7596 | 7.62e-2 |
| Dyn-Nl-OSS P1/P1 | 32 | 4.506e-2 | 9.3227 | 1.69e-2 |
| Dyn-Nl-OSS P1/P1 | 64 | 1.000e-2 | 9.4467 | 3.79e-3 |

- The Dyn-Nl-OSS P1/P1 velocity error drops 4.65× from 16 to 32. Its energy is within 1% on mesh 64, and the energy error drops 4.46× from 32 to 64.
- The first run, which still had the 1% bound on mesh 32, failed at 1.69%. Other P1/P1 cases on mesh 32 (scratch): Dyn-Lin-OSS −1.68%, Sta-Lin-OSS −1.72%, Dyn-Nl-ASGS and Dyn-Lin-ASGS −2.16%. Dyn-Nl-OSS P2/P2 on mesh 32: −0.002%.

Deviations from spec or plan:
1. Lagged projections were implemented first, as the spec and your Phase 8 point 3 then said, and replaced by coupled projections on your decision; the spec now prescribes this. The lagging evidence is above.
2. The projection equations live in `formulation/oss_projections.py`, not `solver/oss_projections.py`, because nothing is solved separately any more. Plan §3, §5 and §13 are updated.
3. With linear splitting an OSS row's xi_change and zeta_change are the change over the step, like velocity_change for the cases without an outer loop.
4. Rule 5's check and phase table are removed, because all eight cases are implemented. The old rejection tests are replaced by a test that all eight load.
5. U5 for OSS uses the Schur complement of the coupled bilinear form instead of assembling column by column; the two give the same matrix.
6. U6 measures orthogonality relative to |(τ_t X, w_h)| and runs with inner_rtol 1e-12 and inner_atol 1e-14, so the solver tolerance does not dominate the measured value.
7. V1 adds Sta-Nl-OSS P1/P1 to the meshes 16, 32, 64 exception (plan §10 V1).
8. V4 (your decision 1): for Dyn-Nl-OSS P1/P1 the energy bound is checked on mesh 64, together with a ≥ 3× drop of the energy error from 32 to 64; Galerkin P2/P1 keeps the 1% bound on mesh 32 (plan §10 V4, with the reason).
9. `reusable_assembler` and `RelativeChange` are performance changes that use Firedrake internals. The results are unchanged byte for byte; U9 guards them and the README records the pinned version (your decision 3).
10. The finite-difference Jacobian test and the nonlinear run tests now cover OSS, and U9 is a new row in plan §10.

Open questions: None. Resolved: exact-history startup keeps zero subscale values at negative times (your decision 2). The subscale has no exact history in general, the spec states the first-step defect, and U8 covers BDF3 with a subscale that vanishes smoothly at t = 0.

## Phase 9: Diagnostics and 3D verification (2026-10-22)

Done:
- `postprocessing/diagnostics.py` adds `subscale_kinetic_energy`, `linear_momentum`, `angular_momentum` and `force`. Every integral uses global `assemble`; the common quadrature rule is used except for error norms, which retain degree 2q.
- Subscale energy is ½∫|ũ|². Dynamic cases read the stored endpoint subscale; static nonlinear cases read the converged stored iterate; static linear cases evaluate the converged closure expression. All eight VMS variants are covered. Configuration rejects subscale energy with `method = "none"`.
- Momentum vectors expand into `_x`, `_y`, `_z` columns; 2D angular momentum is the scalar ∫(x u_y − y u_x). `Problem.force_boundaries` is an optional tuple of boundary names. `force` expands into `force_<boundary>_<axis>` and measures −∫σn ds separately on each declared boundary, with σ = −pI + ν(∇u + ∇uᵀ). Requests without declared force boundaries fail setup validation; unknown boundary names are rejected. No drag/lift coefficient normalization is added.
- `diagnostics.csv` now includes `pressure_time` and `subscale_time`. `time` labels endpoint velocity diagnostics. Pressure error and force use `pressure_time` and the solve-level velocity, including the Crank–Nicolson half level. Dynamic subscale energy uses the endpoint time; static energy uses the solve time; Galerkin's subscale timestamp is blank. The convergence-study driver reads expanded diagnostic columns and excludes time metadata. The earlier steady-test CSV reader now skips blank optional metadata.
- `cases/ethier_steinman_3d.py` implements the cyclic exponential–trigonometric exact velocity on [−1,1]³, with a = π/4, **d = π**, ν = 1, exact Dirichlet data on all six faces and all six faces declared for force measurement. For this Beltrami flow p = −|u|²/2; forcing is derived through the existing UFL helper. Reference: [Ethier and Steinman (1994)](https://doi.org/10.1002/fld.1650190502).
- The time stepper checks the interpolated Dirichlet trace at every solve level when the pressure null space is active and physical boundaries exist. The check globally assembles the net flux and the integral of the absolute normal flux, rejecting |∫u_D·n ds| > 1e-10 ∫|u_D·n| ds. The error names the step, solve time, measured flux and scale, and explains discrete flux compatibility. Crank–Nicolson uses the interpolated endpoint average. Outflow problems skip this check; fully periodic problems have no physical boundary flux.
- The default quadrature rule is retained, as explicitly requested. The owner's updated spec paragraph explains partition-dependent tetrahedral quadrature; README records the resulting VMS process-count dependence. No quadrature degree, point set, solver tolerance or refinement threshold was changed to hide discrepancies.

Tests:
- `OMP_NUM_THREADS=1 pytest tests/unit -n 4 -q`: **327 passed in 51.91 s**, including all earlier unit tests.
- `ruff check .`, `ruff format --check .`, `git diff --check`: clean.
- Analytic diagnostics tests, without solves: uniform flow and rigid rotation momentum/angular momentum in 2D/3D; plane Poiseuille wall shear and pressure force signs/magnitudes, including inlet/outlet and spanwise faces; subscale energy evaluated directly and through stored quadrature fields in 2D/3D, k = 1,2. The energy test also distinguishes q from 2q.
- Time-stepper wiring tests cover all eight VMS cases at BDF2 and Crank–Nicolson levels. A separate prescribed-field test distinguishes half-level force velocity from the extrapolated endpoint, and checks pressure error at the half time. The existing output test checks startup and half-level CSV/VTK pressure times.
- Ethier exact forcing and divergence norms are below 1e-13 at t = 0, 0.17 and 1. Interpolated P1/P2 boundary flux is below 1e-13 on meshes 4 and 8. Guard tests reject d = π/2 under both BDF2 and Crank–Nicolson, accept d = π, detect boundary data that become incompatible at the second step, verify endpoint averaging, and skip the check with outflow.
- `OMP_NUM_THREADS=1 pytest tests/verification -m slow -n 4 -q`: **63 passed in 615.35 s**, including all 61 earlier slow tests. The 3D tests complete **all eight runs** before comparing diagnostics: both pairs, meshes 4/8, serial/two ranks. All steps converge and both 3D tests pass the final rules. The preceding run passed all 61 earlier slow tests, with only the two 3D assertions failing under the superseded MPI rules (61 passed, 2 failed in 643.07 s).

3D refinement: three steps with BDF2, δt = 0.01, T = 0.03, lower-order startup (BDF1, BDF2, BDF2), default solver tolerances, all diagnostics enabled and VTK disabled:

| Formulation | Pair | Mesh-4 velocity L² error | Mesh-8 velocity L² error | Drop factor | Observed rate |
| --- | --- | --- | --- | --- | --- |
| Galerkin | P2/P1 | 0.1155107144313311 | 0.016579676135947765 | 6.96700668 | 2.80053895 |
| Dyn-Nl-OSS | P1/P1 | 1.0057415840094566 | 0.2672772792788728 | 3.76291463 | 1.91185056 |

These are coarse-mesh smoke checks, not asymptotic convergence certification. Rates approach the expected 3 and 2; the Galerkin rate is lower on these two meshes. The tests require only a decrease, as requested, and were not tuned to achieve a target rate.

Largest serial/two-rank differences across all three steps in the final acceptance run (raw relative differences, before applying an absolute floor):

| Formulation | Mesh | Largest absolute difference | Largest relative difference | Where |
| --- | --- | --- | --- | --- |
| Galerkin | 4 | 4.28033e-10 | 5.58377e-10 | pressure_error_l2, step 2 |
| Galerkin | 8 | 2.81997e-13 | 3.01627e-11 | Absolute: pressure_error_l2, step 2; relative: linear_momentum_x, step 1 (value 1.18547e-4) |
| Dyn-Nl-OSS | 4 | 2.68987e-3 | 1.58168e-2 | Absolute: force_y_min_y, step 1; relative: force_y_min_y, step 3 |
| Dyn-Nl-OSS | 8 | 1.99661e-4 | 6.21292e-4 | Absolute: force_y_min_y, step 1; relative: force_y_min_y, step 2 |

- Galerkin's polynomial diagnostics agree within 3.01627e-11 relative (the corresponding absolute discrepancy is 3.576e-15). Its pressure error contains the non-polynomial exact solution; mesh 4 gives 0.7665658362827 versus 0.7665658367107329. The other Galerkin errors meet 1e-10 relative.
- VMS diagnostics other than force agree within 5.93452e-5 relative (pressure_error_l2 on mesh 4, step 2), so they meet the revised 1e-4 bound. The absolute floor used in the test is 1e-4 in the nondimensional benchmark units; Galerkin retains 1e-13.
- The maximum VMS force discrepancy is −14.94078829658439 versus −14.938098426301291: 2.68987e-3 absolute and 1.80035e-4 relative. Near cancellation at step 3, the same component is −0.11953153045748921 versus −0.11764092947071975: 1.89060e-3 absolute and 1.58168e-2 relative. Kinetic energy alone (about 6.17e-7 relative difference at step 1) understates the pressure/force sensitivity. The maximum absolute force discrepancy decreases from 2.68987e-3 on mesh 4 to 1.99661e-4 on mesh 8, about 13.5×. This decrease is reported, not asserted. An earlier complete 3D run gave 9.05368e-5 on mesh 8 (29.7× decrease); the table uses the final run. The serial force scales are 14.9407883 and 8.45833573, giving absolute tolerances 0.0149407883 and 0.00845833573, respectively.

Boundary-flux evidence and parameter decision:
- With the initially requested d = π/2, at t = 0.01 the P1 trace has net flux **0.23533385276573374** on mesh 4 and **0.058833463191433504** on mesh 8. P2 flux is round-off zero. The P1 constant-test continuity equation is inconsistent; Newton stalls at residual 0.02104890 despite updates near 1e-13. This is a data-compatibility defect, not failed nonlinear convergence.
- Choosing d = π makes both nodal traces discretely flux-free, while retaining the same exact-solution family and exact nodal boundary values. The owner explicitly accepted this parameter choice. The new guard diagnoses incompatible data before any solve instead of silently projecting away the residual.

Quadrature-orientation evidence (independent analytic integral, no solve): ∫exp(0.7x + 1.1y + 1.7z) dx on the default 4³ tetrahedral mesh of [−1,1]³:

| Degree | Serial | Two ranks |
| --- | --- | --- |
| q = 3 | 16.377583504693483 | 16.377554366139236 |
| q = 6 | 16.382309127302463 | 16.382309128852057 |
| q = 12 | 16.382309179778755 | 16.382309179778770 |

The difference decreases to round-off at q = 12. The installed FIAT default tetrahedral rule is not permutation symmetric, and Firedrake's simplex vertex ordering uses distributed global numbering. Disabling mesh reordering does not remove the effect. This probe supports the spec's quadrature-error explanation; production still uses q = 3k, with errors at 2q, and unchanged default point sets.

Memory, measured on separate serial mesh-8 runs with `/usr/bin/time -v`, including the flux guard, all diagnostics, Python/Firedrake and direct LU:

| Formulation | Unknowns | Maximum RSS | Wall time |
| --- | --- | --- | --- |
| Galerkin P2/P1 | 15,468 | 544,644 KiB = 531.88 MiB | 9.75 s |
| Dyn-Nl-OSS P1/P1 | 5,832 | 432,636 KiB = 422.50 MiB | 12.62 s |

These are process RSS measurements, not isolated factorization memory or performance benchmarks; the slow suite shared the machine.
For the largest case on two ranks, `resource.getrusage(RUSAGE_SELF).ru_maxrss` measured 469,528 and 463,184 KiB (458.52 and 452.33 MiB). The sum of the rank peaks is 910.85 MiB, an upper bound on simultaneous summed process RSS, not a measurement of unique physical RAM or simultaneous peak allocation.

Deviations from spec or plan:
1. d = π replaces π/2 for the documented P1 trace-flux compatibility reason, as accepted by the owner. a = π/4 and ν = 1 are unchanged.
2. Raw boundary forces are implemented; drag/lift coefficients are deferred to Phase 10, per the Phase 9 clarification.
3. Pressure and subscale time columns make half-level versus endpoint diagnostics explicit. Static subscale energy uses the solve level; dynamic energy uses the stored endpoint.
4. The owner accepted MPI rules based on the integrands and force scale, not thresholds fitted to the measured discrepancies. Galerkin polynomial-integrand diagnostics retain 1e-10 relative; error norms against exact solutions use 1e-8 relative because their integrands are non-polynomial. Both use absolute floor 1e-13. VMS forces use absolute tolerance 1e-3 times the largest absolute force component in each serial reference run, over all faces, components and steps. Every other VMS diagnostic retains 1e-4 relative with absolute floor 1e-4 in the nondimensional benchmark units. Plan §13 and the test record these rules. Default quadrature remains unchanged.
5. Galerkin forces agree within 5.7351e-13 relative (maximum absolute difference 2.4869e-14), confirming that the force diagnostic is MPI-correct. VMS force differences come from the pressure, which carries the quadrature-orientation noise through τ; boundary tractions on the coarse mesh amplify it. These bounds allow that quadrature error while detecting O(1) differences from MPI bugs.

Open questions: None. The d = π choice, flux guard, retained default quadrature and diagnostic-specific MPI acceptance rules were explicitly resolved by the owner.

## Phase 10a: Validation cases (2026-10-22)
Done: Taken over mid-phase from another agent, whose work was uncommitted. All of it was read and treated as unverified; the fixes are listed below.
- Reference data, each folder with a README.md (source, table or page, method, date):
  - `data/references/lid_driven_cavity/`: Ghia, Ghia and Shin (1982), Tables I–II, u at x = 0.5 and v at y = 0.5 for Re = 100 and 1000, at the paper's 17 stations.
    - I re-verified every value from 400 dpi renderings of PDF pages 12–13: all 68 values and 34 stations agree.
    - The commonly quoted extremes agree with the PDF, where they are underscored as local extrema: Re 100 u min −0.21090 at y = 0.4531, v max 0.17527 at x = 0.2344, v min −0.24533 at x = 0.8047; Re 1000 u min −0.38289 at y = 0.1719, v max 0.37095 at x = 0.1563, v min −0.51550 at x = 0.9063.
  - `data/references/cylinder_dfg_2d/`: FEATFLOW's finest level (6, k = 1/200) in `featflow.json`. I re-checked all nine values against the page: C_D min/max/mean/amplitude 3.1569/3.2200/3.1884/0.0631; C_L −1.0206/0.9859/−0.0173/2.008; St 0.30188.
    - The 1996 intervals (Schäfer and Turek, Table 4) are kept as context.
    - All five FEATFLOW archives are kept unchanged, with SHA-256 sums; they pass ZIP CRC checks.
    - `featflow_curve.csv` is the level-6, k = 1/1600 force history from the newer 20141124 archive.
  - `data/references/lattice_vortex/`: the thesis §6.2 solution and its exact energy ¼ exp(−16νπ²t). Re-checked against the thesis text; the thesis used 3 fixed Newton iterations per step, while we iterate to the tolerances.
- Cases (each cites its source in the module docstring):
  - `lid_driven_cavity`: unit square, lid u = (1, 0), Re through case_parameters. Both top corner nodes belong to the stationary walls, with the same value on both incident boundaries.
  - `cylinder_dfg_2d`: DFG 2D-2, inflow 6y(0.41 − y)/0.41², ν = 1e-3, do-nothing outlet, `force_boundaries = ("cylinder",)`. The docstring notes that FEATFLOW's ν∇u − pI and our symmetric σ give the same force on a no-slip wall.
  - `lattice_vortex`: the thesis setup, with exact Dirichlet data.
  - Comparison functions live in the case modules: Ghia maximum deviations; FEATFLOW cycle statistics, from the last lift-minimum-to-minimum cycle, with mean = (max + min)/2, amplitude = max − min, St = D/(Ū × cycle length) and the change in max C_L over the last three cycles; curve deviations aligned at the C_L minimum; the outlet u·n minimum; the energy decay error.
- Cylinder meshes (`meshes/cylinder_dfg_2d.geo`, `-setnumber refine r` divides every production size by r):
  - coarse test mesh, 74 vertices;
  - level 1 (r = 1, unchanged): 1,983 vertices;
  - level 2 (r = 2.2): 9,075 vertices;
  - level 3 (r = 4): 29,344 vertices.

  All levels refine at the cylinder and in the wake box (0.25 ≤ x ≤ 1.5, 0.1 ≤ y ≤ 0.3); the cylinder gets about 78, 172 and 314 boundary edges. The coarse mesh and level 1 are reproduced exactly from the `.geo`.
- `scripts/run_validation.py <config> [--analyse-only | --run-only] [--smoke | --steps N]`:
  - Runs a case and writes `validation.csv`, written last as the completion marker, plus `validation.md` and per-case tables.
  - Cavity runs stop at steady state: relative L² change between time levels < 1e-8 on three consecutive steps, through a new `after_step` callback of `TimeStepper.run`.
  - Cylinder runs log the outlet u·n minimum every step.
- `configs/validation/`: 13 configs.
  - Cavity Re 100: Dyn-Nl-OSS P1/P1 64 and 128, P2/P2 64, Taylor–Hood 64.
  - Cavity Re 1000: Dyn-Lin-OSS P1/P1 64 and 128, P2/P2 64 (see deviation 1).
  - Lattice vortex: Dyn-Nl-OSS and Dyn-Nl-ASGS P1/P1, 50×50, Crank–Nicolson δt = 0.01, T = 10.
  - Cylinder level 1: Taylor–Hood, Dyn-Nl-OSS P1/P1 and P2/P2 with BDF2 δt = 0.005, and Dyn-Nl-OSS P1/P1 with δt = 0.0025, all to t = 30. The level is in the run name and the mesh file is explicit.
- `scripts/run_validation_queue.sh [--smoke]`:
  - Runs cheapest first within cavity, lattice vortex and cylinder (mesh size × degree² × steps), each with `OMP_NUM_THREADS=1 mpiexec -n 4`, then the serial analysis.
  - Logs config, status and wall time to `validation_queue.log`, continues after failures, skips runs with `validation.csv`, and moves an interrupted folder to `.incomplete-<timestamp>`.
  - Writes `VALIDATION_DONE` and calls `notify-send` if present. `--smoke` runs 3 steps per config into `results/validation_smoke/`.
- Fixes to the handed-over work:
  - **Point evaluation:** the smoke queue failed with "Point evaluation gave different results across processes". `Function.at` evaluates a point on a partition boundary on every rank that touches it and rejects round-off differences (it is also deprecated). Centreline and outlet sampling now use `fd.PointEvaluator`, a VertexOnlyMesh that gives each point one owner. The MPI "Sorry!" text was only Open MPI's missing help file for a non-zero exit.
  - **Queue stdin:** the queue then stopped after its first run ("0 completed, 0 skipped, 1 failed"). `mpiexec` read the loop's stdin and consumed the remaining config lines; both commands now read `/dev/null`.
  - **Ordering:** the queue order includes the mesh size of gmsh meshes and the number of steps.
  - Ruff was already clean at the handover; the unit run left behind had passed (347).
- README: a validation section (commands, queue behaviour). Plan: §4 layout and §13 Phase 10a.

Tests:
- `ruff check .` and `ruff format --check .` are clean (65 files).
- `pytest tests/unit -n 4`: 348 passed in 35 s.
  - New: `test_validation_cases.py` (9 tests: transcription and stations, comparison functions on synthetic data including a pure-sine lift of known frequency and amplitude, cylinder names and outlet minimum, reference archives and curve, configs) and `test_after_step_callback_can_stop_the_run_early`.
  - U5 and the case-build test now include the cylinder meshes.
- `OMP_NUM_THREADS=1 pytest tests/verification -m slow -n 4`: 63 passed in 9 min 50 s, all earlier verification tests.
- U5 on the coarse cylinder mesh (74 vertices; Dirichlet on every boundary for the eigenvalue check) passes for ASGS and OSS, P1/P1 and P2/P2, ν = 10 and 1e-4. The worst λ_min/λ_max is −1.0e-16 (OSS P2/P2, ν = 10: −1.55e-14 against 153).
- Smoke queue on 4 ranks: 13 completed, 0 failed, in 4 min.
  - Rerun: 12 skipped and 1 rerun after its `validation.csv` was deleted, with the old folder moved to `.incomplete-20261022T175321567086162`.
  - Smoke outputs: the cavity profile comparison runs on 4 ranks; the cylinder reports "insufficient lift minima" instead of inventing cycle statistics; the lattice energy error at t = 0.03 is 0.45%, the interpolation error of the initial state.

Runtime estimate. Probes of each config on 4 ranks with the real mesh, δt and solver settings: 30 steps for the cavities (they stop at steady state), 50 for the lattice vortex and cylinder levels 1–2, 20 for level 3. A projection is the probe's setup plus first two steps, plus the mean time of the later steps times the remaining steps. The cylinder probes cover only the start from rest (t ≤ 0.25), so the shedding phase may cost somewhat differently. RSS is the sum over the 4 ranks.

| Run | Probe s/step | Outer per step | Steps | Projected | Peak RSS |
| --- | --- | --- | --- | --- | --- |
| cavity Re 100 Dyn-Nl-OSS P1/P1 64 | 0.74 | 3.7 | 17 (steady) | 0.01 h | 1.5 GB |
| cavity Re 100 Taylor–Hood 64 | 0.18 | 1 | 17 (steady) | < 0.01 h | 1.4 GB |
| cavity Re 100 Dyn-Nl-OSS P1/P1 128 | 1.94 | 2.4 | 17 (steady) | 0.02 h | 2.0 GB |
| cavity Re 100 Dyn-Nl-OSS P2/P2 64 | 2.08 | 2.4 | 17 (steady) | 0.02 h | 2.0 GB |
| cavity Re 1000 Dyn-Lin-OSS P1/P1 64 | 1.04 | 1 | ~59 | 0.02 h | 1.6 GB |
| cavity Re 1000 Dyn-Lin-OSS P1/P1 128 | 4.44 | 1 | ~59 | 0.10 h | 2.5 GB |
| cavity Re 1000 Dyn-Lin-OSS P2/P2 64 | 3.75 | 1 | ~59 | 0.06 h | 2.7 GB |
| lattice Dyn-Nl-ASGS P1/P1 | 0.88 | 9.0 | 1,000 | 0.24 h | 1.3 GB |
| lattice Dyn-Nl-OSS P1/P1 | 0.88 | 5.7 | 1,000 | 0.25 h | 1.5 GB |
| cylinder L1 Taylor–Hood δt 0.005 | 0.11 | 1 | 6,000 | 0.19 h | 1.4 GB |
| cylinder L1 Dyn-Nl-OSS P1/P1 δt 0.005 | 0.70 | 7.2 | 6,000 | 1.16 h | 1.5 GB |
| cylinder L1 Dyn-Nl-OSS P1/P1 δt 0.0025 | 0.64 | 6.5 | 12,000 | 2.13 h | 1.5 GB |
| cylinder L1 Dyn-Nl-OSS P2/P2 δt 0.005 | 2.71 | 8.4 | 6,000 | 4.53 h | 1.7 GB |
| cylinder L2 Taylor–Hood δt 0.005 | 0.65 | 1 | 6,000 | 1.08 h | 1.7 GB |
| cylinder L2 Dyn-Nl-OSS P1/P1 δt 0.005 | 2.89 | 7.0 | 6,000 | 4.81 h | 1.7 GB |
| cylinder L2 Dyn-Nl-OSS P1/P1 δt 0.0025 | 2.49 | 6.1 | 12,000 | 8.30 h | 1.7 GB |
| cylinder L2 Dyn-Nl-OSS P2/P2 δt 0.005 | 13.80 | 8.1 | 6,000 | 23.0 h | 2.9 GB |
| cylinder L3 Taylor–Hood δt 0.005 | 1.85 | 1 | 6,000 | 3.09 h | 2.4 GB |
| cylinder L3 Dyn-Nl-OSS P1/P1 δt 0.005 | 12.07 | 7.2 | 6,000 | 20.1 h | 2.5 GB |
| cylinder L3 Dyn-Nl-OSS P1/P1 δt 0.0025 | 12.22 | 6.8 | 12,000 | 40.8 h | 2.5 GB |
| cylinder L3 Dyn-Nl-OSS P2/P2 δt 0.005 | 60.07 | 8.2 | 6,000 | 100.1 h | 6.6 GB |

The committed queue (cavity, lattice vortex, cylinder level 1) totals about 9 h; every probed run together is about 210 h.

Deviations from spec or plan:
1. **Re 1000 cavity uses linear splitting (Dyn-Lin-OSS), with 5 Picard iterations before Newton.** With Dyn-Nl-OSS at δt = 5:
   - The P1/P1 128 and P2/P2 64 runs diverge in Newton at step 1, from rest. With the Picard warm-up the inner solve converges, but every variant then fails the local subscale update at step 1: δt = 5 with warm-up, δt = 1, 0.5 and 0.05 fail; only δt = 0.01 works.
   - The update fails even with 1,000 sweeps, accelerated or fixed point.
   - A probe of P1/P1 64 at δt = 0.5 shows 6 of 49,152 quadrature points still changing after 60 sweeps, all in the cells at the two top corners (x ≈ 0 or 1, y ≈ 0.99–1).
   - There the lid jump gives |∇u_h| ≈ 1/h, while τ reaches h²/(4ν) ≈ 0.06 at Re 1000. So τ|∇u_h| ≈ 4, the near-singular local problem of the spec paragraph on static nonlinear subscales. Dynamic subscales do not help unless δt|∇u_h| is moderate, which needs δt ≈ 0.01: about 20,000 steps to steady state, roughly 17 h for P1/P1 64 alone.
   - Re 100 is unaffected (τ|∇u_h| ≈ 0.4).
   - Per the spec, the failure is reported, not repaired. The configs say so in a comment; this is open question 1.
2. Cylinder configs name the mesh level (`_l1_`) and give the mesh file explicitly; they run level 1 until you choose the levels (open question 2).
3. `run_validation.py` has `--run-only` and `--steps N` besides the requested flags: the queue runs the simulation under MPI and the analysis serially, and `--steps` serves the runtime probes.
4. Cavity steadiness is judged between time levels with a 1e-8 relative tolerance on three consecutive steps, at δt = 5 as in V3.
5. Cycle extrema are sampled values. Lift minima, which set the cycle and St, are refined by a local quadratic; FEATFLOW's 1/f values are quantised by their sampling.

Open questions:
1. Re 1000 cavity: keep Dyn-Lin-OSS (feasible, about 10 min for all three); run Dyn-Nl-OSS at δt = 0.01 (about 17 h for P1/P1 64, more for 128 and P2/P2); or drop the Re 1000 VMS runs.
2. Cylinder levels and cuts, against the 24 h budget: the proposal is in the report to the owner and is decided there.

## Phase 10a: final validation run set (2026-10-22)
Done (your decisions on the Phase 10a open questions):
1. **Re 1000 cavity:** it keeps Dyn-Lin-OSS with 5 Picard iterations before Newton.
   - Your spec update records the condition in the paragraph after the eight-cases table: at a boundary-velocity jump τ_m|∇u_h| ≈ Re·h/(c1·k⁴), with the Phase 10a measurements (about 3.9 for Re 1000 on 64×64 P1, where 6 of 49,152 quadrature points failed, all in the two corner cells; about 0.39 at Re 100).
   - The headers of the three `cavity_re1000_linoss_*` configs now cite that paragraph and the value.
2. **Cylinder run set (your choice of the proposed set, t_end = 30 as in FEATFLOW's procedure):** Taylor–Hood on levels 1, 2 and 3; Dyn-Nl-OSS P1/P1 at δt = 0.005 on levels 1 and 2; Dyn-Nl-OSS P2/P2 on level 1; Dyn-Nl-OSS P1/P1 at δt = 0.0025 on level 1.
   - New configs: `cylinder_th_l2_dt0.005`, `cylinder_th_l3_dt0.005`, `cylinder_oss_p1_l2_dt0.005`. They differ from their level-1 versions only in run name and mesh file.
   - With all 7 cavity and both lattice-vortex runs, the projected queue total is about 18 h: cavities 0.2 h, lattice vortex 0.5 h, cylinder 17.0 h. The cylinder projections cover only the start from rest, so allow about ±30 %.
3. **Queue order:** the cost key mesh size × degree² × steps put Taylor–Hood level 3 (3.1 h) after OSS P1/P1 level 2 (4.8 h).
   - The key is now vertices × unknowns per vertex × outer iterations × steps. P2 has about 4 nodes per vertex, OSS doubles the fields for ξ_h and ζ_h, and nonlinear splitting is weighted by 7 outer iterations per step, the mean measured in the probes.
   - The cylinder order is now the measured order: Taylor–Hood L1 (0.19 h), Taylor–Hood L2 (1.08 h), OSS P1/P1 L1 (1.16 h), OSS P1/P1 L1 δt 0.0025 (2.13 h), Taylor–Hood L3 (3.09 h), OSS P2/P2 L1 (4.53 h), OSS P1/P1 L2 (4.81 h).
   - Cavities and lattice vortex still come first; all cavity runs take minutes.

Possible second batch (cut from this queue; probe projections on 4 ranks, t_end = 30):

| Run | Projected | Peak RSS |
| --- | --- | --- |
| Dyn-Nl-OSS P1/P1 L3, δt 0.005 | 20.1 h | 2.5 GB |
| Dyn-Nl-OSS P1/P1 L3, δt 0.0025 | 40.8 h | 2.5 GB |
| Dyn-Nl-OSS P2/P2 L3, δt 0.005 | 100.1 h | 6.6 GB |
| Dyn-Nl-OSS P2/P2 L2, δt 0.005 | 23.0 h | 2.9 GB |
| Dyn-Nl-OSS P1/P1 L2, δt 0.0025 | 8.3 h | 1.7 GB |

If the first batch's periodicity check shows established shedding well before t = 15, these could run to t = 15 at about half the cost.

Tests:
- `ruff check .` and `ruff format --check .` are clean. `pytest tests/unit -n 4`: 348 passed. The config test now expects 16 configs, 7 of them cylinder.
- Smoke queue on 4 ranks from a clean folder: 16 completed, 0 failed, in 3 min 44 s, in the new cost order. The new runs load their meshes: Taylor–Hood L2 80,875 unknowns, Taylor–Hood L3 262,656, OSS P1/P1 L2 54,450.
- The slow suite was not rerun: nothing in `src/` changed since the Phase 10a run (63 passed), only configs, the queue's ordering, a test count and documents.

Deviations from spec or plan: none beyond the queue cost key, item 3.

Open questions: None. Both Phase 10a questions are resolved by your decisions above. The full queue is yours to launch: `bash scripts/run_validation_queue.sh`.
