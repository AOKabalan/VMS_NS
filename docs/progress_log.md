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
