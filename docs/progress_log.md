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
