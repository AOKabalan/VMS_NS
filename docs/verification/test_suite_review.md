# Test suite review (2026-11-03)

Reviewed all 41 test modules: 35 in `tests/unit/` and six in
`tests/verification/`, along with their shared helpers. Most tests protect
distinct numerical or user-facing behavior. The cleanup removes nine test
functions and redundant parameter combinations, reducing collection from
1,141 to 1,083 cases (1,072 to 1,014 unit cases; all 69 verification cases
retained). These counts include the existing strict Gresho xfail.

Production code, numerical tolerances, configurations and reference data are
unchanged. Historical phase reports retain their original test counts.

## Removal decisions

| File in `tests/unit/` | Removed or consolidated | Coverage retained |
| --- | --- | --- |
| `test_config.py` | Direct startup/diagnostic-validator tests and 24 repeated VMS parser combinations. | Problem construction checks exact-history eligibility; CLI tests reject unknown diagnostics before creating results. All eight VMS variants still load, with both degrees and pressure-subscale settings represented in each ASGS/OSS family. Degree and boolean parsing do not branch on dynamics/splitting in the current implementation. Numerical variant matrices remain exhaustive. |
| `test_quadrature.py` | Four generic scalar quadrature-element sampling cases. | The production subscale-space test exercises the same quadrature-element constructor and integration rule on triangles/tetrahedra at P1/P2. Positive weights, independent monomial integrals, boundary integration and negative-weight rejection remain. |
| `test_boussinesq_vms.py` | Two manually constructed weighted-mass projection cases. These built the projection equations inside the test rather than calling the production projection builders. | Production OSS projection identity tests, solved coupled orthogonality at P1/P2 with pressure projection on/off, thermal feedback-sign checks, and solved energy-exchange budgets. The coupled orthogonality test also checks the scalar subscale paired with gravity times velocity tests. |
| `test_chapter5_solver.py` | Twelve periodic runs asserting only convergence and finite velocity/pressure. | The periodic momentum-budget matrix solves all the same convection-form/ASGS/OSS/nonlinear combinations and additionally checks the discrete balance. Open-boundary and Picard runs, independent EMAC residuals and Jacobian tests remain. |
| `test_meshes_and_boundaries.py` | Standalone pressure-null-space boolean check. | Galerkin and ASGS runs check actual pressure normalization with and without outflow; mixed-system gauge tests and free-slip pressure constraints remain. |
| `test_time_schemes.py` | Eight repeated coefficient-sum consistency cases. | Independent polynomial differentiation includes the constant-function identity. Startup sequences and CN extrapolation/data averaging remain; BDF data and endpoint-weight assertions moved into the polynomial tests. |
| `test_time_stepper.py` | Two extra solves solely for CN output timestamps and Newton iteration logging. | Their assertions moved into existing output and Newton/Picard-equivalence tests using the same simulation setups. |
| `test_conservation_configs.py` | Repeated snapshots of mesh sizes, time steps, solver tolerances and aggregate step counts; three duplicate case builds for long/refined presets that were replaced by the same coarse mesh inside the test. | Complete generated-versus-committed manifest/byte comparison, retired-file cleanup, valid configurations, unique run names, known diagnostics, and case support for Gresho OSS, viscous Gresho ASGS and Taylor–Green. Intended preset tuning no longer requires editing a second copy of its numerical values in tests. |
| `test_validation_cases.py` | Hard-coded preset counts and benchmark mesh/time settings. | Every preset still loads, all three validation cases are represented, and run destinations must be unique. Reference comparisons, sampling, geometry and physical checks remain. |

The removed parameter combinations are not claimed to have identical future
coverage: adding a new interaction between options should bring a regression
test for that interaction. The reductions above follow the current code paths.

## Retained coverage by module

Names below omit the `test_` prefix and `.py` suffix. This inventory includes
the unchanged files as well as the surviving tests in edited files.

| Unit modules | Regression protected |
| --- | --- |
| `asgs_equals_supg`, `chapter5_forms`, `convective_form` | Independent residual expansions, tensor indices, closure/defect identities, periodic skew symmetry and outflow traction/flux. |
| `asgs_solver`, `chapter5_solver`, `boussinesq_vms` | Jacobians, frozen/exact stabilization derivatives, true residual logging, projection spaces, passive-scalar feedback and field layout. |
| `dissipation_sign`, `oss_projections`, `quadrature` | Dissipation eigenvalues, the cell-diameter counterexample, weighted orthogonality and matching quadrature/storage rules. |
| `boussinesq_galerkin`, `boussinesq_nonlinear`, `boussinesq_time` | Buoyancy and flux signs, hydrostatics, thermal/isothermal equivalence, nonlinear scalar cross terms, linearization, scalar histories and temporal order. |
| `conservation_budget`, `momentum_budget`, `thermal_budgets` | Solved balances with independently checked terms, actual time histories, closure defects, boundary reactions and admissible tests. These matrices are retained because their options change the identities being checked. |
| `time_schemes`, `subscale_storage`, `subscale_time_integration` | Time-integration order/startup, history rotation, CN endpoint recovery and local-update convergence. |
| `field_layout`, `time_stepper`, `reusable_assembly` | Mixed-field convergence, solver limits and failures, output lifecycle, pressure gauge and assemblers observing changed state. |
| `meshes_and_boundaries`, `free_slip`, `ethier_exact_solution` | Mesh tags, component constraints, natural traction, discrete boundary-flux compatibility and 3D exact fields. |
| `diagnostics`, `physical_pressure` | Analytic moment/force/norm values, output schemas, solve-level versus endpoint data and EMAC physical-pressure conversion. |
| `gresho`, `validation_cases`, `conservation_configs` | Reference fields/data, documented solver limitations, benchmark sampling/comparison and runnable preset generation. |
| `kelvin_helmholtz`, `shear_layer_startup` | Wall-compatible profiles, scaling/noise, affine clustering, independent horizontal integration, dry/moist initialization and startup controls. |
| `config`, `run_queue` | Parsing and invalid inputs, CLI output/refusal behavior, resumability and preservation of failed/interrupted run records. |
| `analyse_conservation`, `analyse_kh_startup` | Analysis against known synthetic data, scaling, plots and rejection of incomplete or mismatched runs. |

All six verification modules remain: `test_spatial_convergence.py`,
`test_temporal_convergence.py`, `test_steady_equivalence.py`,
`test_taylor_green.py`, `test_ethier_steinman.py`, and
`test_boussinesq_convergence.py`. Refinement rates and serial/MPI comparisons
provide evidence that kernel identities and coarse smoke runs cannot replace.
The thermal MPI probe and convergence-study helper also remain.

For future changes, add a test for a concrete failure mode or behavior. Prefer
an independent numerical expectation or observable result; reuse an existing
run when only additional assertions are needed. A formula written independently
from the specification is a useful oracle, even when it resembles the code.
Do not shrink numerical matrices solely to reduce the reported test count.

## Validation

Using the installed Firedrake environment, `python -m pytest tests/unit -q -n 2`
passed **1,013 cases with one expected failure**, in 537.30 seconds. The xfail
is the existing near-inviscid nonlinear dynamic ASGS Gresho limitation.
The run log is `/tmp/vms-test-review-all.log`.

`python -m pytest tests --collect-only -q` collected all 1,083 remaining cases.
The 69 slow verification cases were reviewed and collected, but not rerun for
this test-only change. Ruff checks over `src`, `tests` and `scripts`, formatting
checks over `tests`, and `git diff --check` passed.
