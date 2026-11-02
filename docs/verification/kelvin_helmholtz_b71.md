# B7.1 Kelvin–Helmholtz coupled smoke verification

Date: 3 November 2026. Scope: the B7.0 case definitions needed for the
[B7.1 checkpoint](../boussinesq_implementation_plan.md#b71--early-coupled-channel-smoke-test).
Measured CLI data: [JSON](kelvin_helmholtz_b71_results.json).
Actual profiles, scaling and noise conventions:
[case definitions](../../data/references/kelvin_helmholtz/README.md).

The case builds and evolves a stratified periodic/free-slip channel with
temperature and passive moisture in one monolithic system. These are short
integration checks, not a growth-rate, billow or cloud validation. B7.2 mesh
and startup work is the next checkpoint; B6 compatibility remains independent.

Unit regression coverage is **1,028 passing checks and one existing xfail**
across the full run and a focused correction rerun. The full run recorded
1,027 passed, one failed and one xfailed in 464.98 s: its sole failure expected
`angular_momentum` to be first in the diagnostic list printed by an invalid CLI
request. The new `advective_cfl` sorts before it. After updating that test
expectation, all five focused CLI checks passed. No simulation code was changed
to address that failure. The final KH-only run passed all 44 tests in 19.01 s.
Ruff lint/format, TOML/JSON parsing, local documentation links, display-math
fences and recorded-measurement consistency also passed.

## Profile and initialization checks

The tanh core is unchanged for `|z| <= 7` with default parameters. A nonnegative
degree-five gradient taper joins it to constant wall plateaus beginning at
`|z|=9`; the profile is C² and monotone. A C² compact envelope suppresses the
Gaussian perturbation and moisture near both walls. Tests check monotonicity,
core agreement, matching across transition endpoints, and zero wall derivatives
for transition widths 0.5, 2 and 2.5. Differentiating the analytic pressure
primitive recovers the adjusted base temperature to an L² tolerance of 3e-14.
The endpoint interpolation check uses DG6 to represent the degree-six integrated
gradient taper; a preliminary DG3 check introduced an interpolation jump and
was corrected without changing the profile.

Noise phases use explicitly selected PCG64 with a recorded seed and are
independent of mesh partition. Distinct periodic modes with amplitude
`1/sqrt(number_of_modes)` give the same spatial RMS as the base cosine;
tests include multiple waves, nonzero noise and a changed seed. Parameters are
validated before mesh construction, including physically meaningful zero
Richardson/perturbation/noise controls. No body force holds the profiles steady.

The optional `Problem.initial_pressure` supplies **physical** pressure.
Initialization converts it to the solver pressure convention if needed and
normalizes the gauge. Isothermal tests also cover existing EMAC/EMAC–VMS
modified-pressure reporting. Old cases omit this data and retain the original
startup/output behavior. Dynamic histories and OSS projections retain zero
initial values; static subscales are determined algebraically by the residual,
not an independently prescribed initial condition.

P1/P2 coupled tests check analytic seed solenoidality, measured interpolation
divergence, normal and tangential wall compatibility, scalar insulation,
pressure gauge, finite residuals/Jacobians and the constant-pressure null mode.
OSS tests check velocity, temperature and moisture orthogonality before history
shifting with dual residuals below 1e-9. A duplicate temperature tracer agrees
with temperature below 1e-11 under both ASGS and OSS, including active thermal
buoyancy. Dry, zero-seed and zero-Richardson controls are also exercised.

## Canonical CLI runs

```bash
python scripts/run_simulation.py configs/boussinesq/kh_smoke_static_asgs.toml
python scripts/run_simulation.py configs/boussinesq/kh_smoke_static_oss.toml
```

Both use a uniform 24×24 triangular channel mesh, P1/P1, static linear
subscales, skew-symmetric momentum/scalar transport and three BDF1 steps of
0.02, ending at 0.06. Coefficients are Re=1000, Pr=Sc=1, Ri=0.1. The explicit
MUMPS null-pivot option handles the constant pressure null space.
Results are under `results/kelvin_helmholtz_shear_layer/b71_static_{asgs,oss}/`.

Both start with zero normal flow, zero temperature/moisture wall flux,
initial divergence L² 8.684e-5 and vertical kinetic energy 7.818e-9.
The interpolated hydrostatic pressure–buoyancy mismatch has L² norm 0.10437;
physical pressure mean is -1.338e-14. These numbers use the common quadrature
of the solved equations. They do not assert discrete hydrostatic equilibrium.

| Quantity at t=0.06 | ASGS | OSS |
| --- | ---: | ---: |
| Vertical kinetic energy | 6.24071e-5 | 6.76793e-5 |
| Divergence L² | 1.57243e-2 | 1.66556e-2 |
| Sampled advective CFL | 0.03396493 | 0.03396493 |
| Moisture content | 36.4650871 | 36.4651107 |
| Maximum absolute kinetic–potential budget residual over three steps | 4.369e-14 | 5.969e-14 |

All reported scalar-content/variance and kinetic–potential budget residuals
are below 6.1e-14 in absolute value. The existing identities include the linear
splitting/skew transport content and exchange defects; this is not a claim of
exact resolved moisture-content conservation. Default convection is unchanged.
Normalized residuals of near-zero thermal content use the existing 1e-12 floor;
absolute residuals are the meaningful companion measure here.

The final vertical energies greatly exceed the tiny initial seed, in the
presence of a sizable initial discrete hydrostatic defect. These coarse data
cannot establish a KH growth rate or distinguish startup modes from the seeded
mode. No growth fit is performed. B7.2 must examine zero-seed startup controls,
profile drift and refinement before B7.3 uses a linear-growth window.

VTK writes `velocity`, `pressure`, `temperature` and `q` separately. JSON
metadata records dimensions, coefficients, all parameters, mesh spacing,
cutoff and generated phases. Initial-state JSON adds the compatibility norms;
CSV contains named tracer columns, vertical energy and actual sampled CFL.
The CFL is the global maximum of `k*dt*|u_h|/MinCellEdgeLength` at common
quadrature points, not a continuous supremum. A separate test verifies shortest
edge, velocity-degree weighting and updates after velocity changes.

## Parallel and compatibility checks

A two-rank MPI run of the OSS config (with a distinct run name) completed all
three steps. Thermal/case metadata, including Fourier phases, is identical.
The maximum initial-state absolute difference is 2.217e-12 (moment reductions);
the maximum unnormalized per-step diagnostic difference is 2.274e-13.
Vertical-energy differences are below 3.94e-19 and sampled CFL agrees exactly.
The largest normalized-residual difference, 0.02618 for near-zero thermal
content, corresponds to an absolute residual difference of only 2.618e-14
divided by the existing floor. MPI required local sockets outside the sandbox;
automatic approval allowed the local test.

Four isothermal and three existing thermal baselines were run against an
archive of G4 commit `a0d4b83`. All field arrays and diagnostic/inner/outer CSV,
initial-state JSON and thermal metadata agree bitwise; runtime wall times are
excluded. Artifacts are `/tmp/vms-b71-before-*` and `/tmp/vms-b71-after-*`.
Logs: `/tmp/vms-b71-focus.log`, `/tmp/vms-b71-profile3.log`,
`/tmp/vms-b71-all.log`, `/tmp/vms-b71-{asgs,oss,oss-mpi}.log`.
Final focused logs are `/tmp/vms-b71-final-focus.log` and `/tmp/vms-b71-cli.log`.

This checkpoint does not enable thermal stress/EMAC compatibility, clustering,
production resolution, growth fitting, long integrations or cloud diagnostics.
