# Integrating thermal Boussinesq physics into the existing framework

Status: B0–B5, boundary prerequisite G4 and B7.1–B7.2 completed, 3 November 2026.
B7 now has its boundary prerequisite; B6 remains the separate compatibility
extension. Equations refer to the companion
[thermal formulation](boussinesq_formulation.md). The objective is one shared
solver for isothermal and thermally coupled problems, with small, explicit
thermal components and unchanged isothermal defaults.
The [review disposition](reviews/boussinesq_addendum_review.md) explains why some
budget hypotheses and KH acceptance criteria were corrected before adoption.

## 1. Findings from the current code

The numerical kernels already provide most of the needed mechanics:

- `formulation/galerkin.py` accepts an explicit `SolveLevel.forcing`, so
  $f+b(T_h)$ can be shared by resolved momentum and its residual.
- `formulation/vms.py` implements the three existing switches, local elimination
  and accelerated nonlinear velocity updates. Thermal subscales can mirror this
  small algebra without duplicating a solver.
- `formulation/oss_projections.py` already solves consistent weighted projections
  in the same mixed system. A scalar temperature equation is the missing block.
- `discretization/time_schemes.py` is independent of field type;
  `SubscaleStorage` also operates on arbitrary supplied function spaces despite
  its velocity-specific names and docstrings.
- `InnerSolver` largely accepts arbitrary mixed residuals. The main structural
  assumptions are in `TimeStepper`, mixed-space construction, `StepSolver`
  reporting, histories, diagnostics and output.

Before B0, `TimeStepper.__init__` interpreted every field after $(u,p)$ as
$(\xi,\zeta)$, and `StepSolver._row` did the same with changes. Simply appending
temperature would silently mislabel convergence data and projection fields.
B0 removes these positional assumptions before wiring in thermal equations.

## 2. Architecture and ownership

Keep plain functions and dataclasses, as required by the existing implementation
plan. Do not introduce a generic multiphysics inheritance hierarchy or a second
copy of `TimeStepper`.

### Physical data and numerical configuration

Add an optional `thermal: ThermalProblem | None = None` field to `Problem`.
`ThermalProblem` contains thermal diffusivity, expansion coefficient, gravity
vector, reference temperature, source, initial temperature, independent thermal
boundary descriptions and optional exact temperature. Put this dataclass in
`problem/thermal_problem.py`; it is physical data, not a form builder. Existing
case builders continue to return `Problem` with no thermal data.

For G1, add an optional ordered tuple of passive-scalar data to `Problem`, default
empty. Each named scalar has diffusivity, source, initial value, independent
boundary data and optional exact solution. Use a small `PassiveScalarProblem`
dataclass in `problem/scalar_problem.py`, not an inheritance hierarchy. Reject
duplicate/reserved field names. Introduce this path with the paired-scalar test
in B2 so that it has an immediate consumer. Initially passive scalars require
thermal problem data and inherit its degree, convection form and all three
stabilization switches. Their equations have no feedback into buoyancy; their
unknowns and OSS projections remain in the shared monolithic solve.

Add `ThermalOptions` in `config.py` only for numerical choices:

```toml
[thermal]
# Presence of Problem.thermal selects the physical model.
temperature_degree = 1
convection = "skew_symmetric"  # nonconservative; conservative added after its budget test
```

Resolve an omitted temperature degree to the velocity degree. For VMS initially
require equality of the velocity and temperature degrees; this supports the
OSS buoyancy cancellation in (B16). The Galerkin thermal baseline uses the
velocity degree for temperature, too. Thermal dynamics, subscale space and
splitting inherit `[stabilization]`; do not add independent combinations before
there is a verified use case. If a thermal section is explicitly supplied for
a case with no thermal data, reject it as ineffective configuration rather than
silently ignoring it. Preserve omission information in parsing as needed.
Omit the inactive thermal section when serializing an isothermal configuration,
so a resolved configuration can be loaded again without tripping this check.

Physical values belong in the case's parameters and resolved `ThermalProblem`,
following the existing viscosity convention. Persist the derived physical
parameters and any nondimensional scaling in run metadata. Extend validation
without reusing `case_parameters()` unchanged: it currently requires every value
to be positive, which excludes valid zero expansion, signed temperatures,
reference temperatures and gravity components. Use field-specific validation:
positive finite $\nu,\kappa$, nonnegative finite $\beta$, finite reference/data
values, a gravity vector of the mesh dimension, and complete thermal boundary
coverage. Do not add an independent `enabled` flag that can disagree with the
presence of thermal problem data.

### Mixed-space layout

Add a small immutable `FieldLayout` alongside the space builders. It maps the
named fields to component indices and distinguishes physical and projection
fields. Centralize splitting of functions, UFL expressions and test functions
through this layout. Keep pressure at component 1 to preserve existing null-space
behavior.

| Run | Ordered fields |
| --- | --- |
| Isothermal Galerkin/ASGS | $u,p$ |
| Isothermal OSS | $u,p,\xi,(\zeta)$ |
| Thermal Galerkin/ASGS | $u,p,T$ |
| Thermal OSS | $u,p,T,\xi,\eta,(\zeta)$ |
| Thermal OSS with passive moisture | $u,p,T,q_m,\xi,\eta_T,\eta_q,(\zeta)$ |

Parentheses mean optional. Thermal Dirichlet conditions act only on $T$;
projections remain unconstrained. The pressure null space acts only on $p$.
Transient insulated temperature has a time mass term, so it must not be given
a constant null space. A future genuinely steady pure-Neumann thermal solve
would need a separate compatibility/gauge treatment.

Add a scalar quadrature-space builder using the existing `quadrature_element`.
Both scalar and vector subscales must use the same nodes and measure as the
resolved equations. `FieldLayout` should be the only place that decides optional
component indices.
Additional passive fields occupy an ordered block after temperature, followed
by velocity and scalar projections and optional pressure projection. Name all
scalar histories, convergence records, CSV columns and VTK series explicitly;
empty passive data must preserve the temperature-only layout and output.

### Form composition

Add `formulation/thermal.py` for scalar convection, thermal residual, temperature
subscale and thermal weak terms (B10), (B14), (B15). Add
`formulation/boussinesq.py` for the resolved force (B1), temperature-subscale
momentum feedback (B13), and composition of the thermal and existing momentum
forms. Keep both modules stateless and independent of nonlinear solvers.
Implement the kernels in `thermal.py` with generic scalar arguments (field,
test, diffusivity, source, histories, subscale and projection). Keep temperature
names in the composition and physical data, not inside transport algebra. The
weighted scalar projection term is equally generic. Do not duplicate these
kernels for moisture.

The composition must construct a momentum `SolveLevel` with
`forcing = external_forcing + buoyancy(T)` and pass the **same expression** into
the Galerkin, VMS, OSS and accelerated-update kernels. Add
$+\beta(g\widetilde T,v)$ separately as the resolved feedback. Do not hide
temperature-subscale buoyancy in `SolveLevel.forcing`, because that would change
the local closure in (B6)–(B7).

Reuse `time_derivative` for temperature. Give temperature its own history/source
record, while sharing the same `alpha0` and `dt` objects. Extract the common
inverse advection–diffusion scale from `stabilization_parameters.py` only as
needed; keep its existing public velocity/pressure wrapper unchanged. Extend
`oss_projections.py` with a scalar weighted projection term. Avoid making every
momentum helper accept optional temperature arguments.

Extract form construction from the growing `TimeStepper.__init__` into a small
`solver/system_builder.py` returning a dataclass with residuals/Jacobian,
subscale expressions and any local-update expression. Dispatch once based on
thermal data, then compose the shared momentum kernels. Keep quadrature, spaces,
boundaries and persistent state owned by their existing layers. Use no second
physics-specific time loop.

### Histories, iteration and output

Extract the resolved history bookkeeping into one small stateful helper when
temperature makes it necessary. Instantiate it for velocity and temperature;
reuse `SubscaleStorage` for scalar and vector dynamic histories. No histories
are allocated for pressure or projections, and no thermal state is allocated
for isothermal cases.
Apply this bookkeeping to each named passive scalar as well. Exact-history
startup requires exact histories for every active transported scalar.

The diagonal closure requires one nonlinear outer iterate, for velocity. Store
dynamic temperature history, but do not invent a temperature outer fixed point:
its current value is already an algebraic expression at fixed advector. Extend
convergence checks to $T,\eta$ and the recomputed thermal subscale as well as
the existing fields. Report named field changes; preserve existing isothermal
CSV columns and append thermal columns for thermal runs. Mixed dimensional
residuals must not be treated as a physical norm: add per-equation scaled residual
checks using recorded characteristic scales, alongside PETSc's algebraic norm,
so a converged momentum block cannot hide an unconverged heat equation.

Newton differentiates both directions of coupling. The existing `newton_tau =
"frozen"` option must freeze thermal and mechanical parameters consistently;
`"exact"` differentiates both, retaining the existing zero-speed derivative
limitation. Later Picard freezes the shared transport velocity while retaining
the unknown temperature in resolved buoyancy. No previous-step buoyancy lag is
introduced.

Extend diagnostic state with endpoint and solve-level temperature and the
temperature subscale. Use a separate `postprocessing/thermal_diagnostics.py`
for the thermal functionals and register them in the existing registry. Extend
the existing conservation context with buoyancy work; do not leave its forcing
term equal to external forcing alone. Temperature VTK output uses the endpoint
time, like velocity; pressure retains its solve-level time. Failure snapshots
must include all three physical fields at the actual available time.
Include every active passive scalar in those snapshots as well, without adding
empty scalar columns or files to isothermal runs.

## 3. File-level changes

All source paths below are relative to `src/navier_stokes_vms/`.

| Path | Change |
| --- | --- |
| `config.py` | Thermal options, parsing, validation and resolved output |
| `problem/problem_definition.py` | Optional thermal data; joint problem checks |
| `problem/thermal_problem.py` (new) | Thermal physical data and exact-solution source helpers |
| `problem/scalar_problem.py` (new) | Named passive-scalar data; no buoyancy parameters |
| `problem/boundary_conditions.py` | Scalar temperature/flux descriptions and construction; independent thermal map |
| `problem/mesh_loading.py` | Explicit x-periodic channel constructor with verified wall tags |
| `discretization/function_spaces.py` | Named layout, temperature/projection blocks, scalar quadrature space |
| `discretization/quadrature.py` | Common rule based on maximum transported-field degree |
| `formulation/thermal.py` (new) | Thermal operators, weak terms and scalar subscale |
| `formulation/boussinesq.py` (new) | Buoyancy and composition using existing momentum forms |
| `formulation/stabilization_parameters.py` | Reuse the spatial scale for diffusivity $\kappa$ |
| `formulation/oss_projections.py` | Scalar projection equation with thermal effective weight |
| `solver/system_builder.py` (new) | Shared construction of coupled solve forms |
| `solver/field_history.py` (new when needed) | Resolved field histories, initialization and endpoint recovery |
| `solver/subscale_storage.py` | Reuse for scalar fields; adjust names/docstrings only where needed |
| `solver/nonlinear_solver.py` | Named changes, thermal convergence and coupling-aware Picard |
| `solver/time_stepper.py` | Coordinate additional fields, data levels, history commits and output |
| `postprocessing/thermal_diagnostics.py` (new) | Thermal errors, content, variance, flux and Nusselt diagnostics |
| `postprocessing/diagnostics.py`, `output.py`, `run_logging.py` | Thermal state, buoyancy budgets, field output and logs |
| `cases/__init__.py` and new thermal cases | Manufactured solutions, diffusion, passive transport and cavity |
| `cases/kelvin_helmholtz_shear_layer.py` (new in B7) | Stratified shear layer and optional moisture data |

Keep `scripts/run_simulation.py` as the single entry point. Add thermal TOML
examples under `configs/boussinesq/`, reference data under
`data/references/differentially_heated_cavity/`, and thermal verification under
the existing test directories. Do not change the current agreed flow formulation
to imply that the paper's different closures were already implemented.

## 4. Phased implementation and acceptance gates

### B0 — Prepare shared field handling without changing equations

Introduce the named mixed-field layout and extract form construction/history
handling only where needed by the next phase. Preserve existing signatures
through narrow wrappers if useful. Replace positional projection assumptions in
`TimeStepper` and `StepSolver`.

Acceptance: existing configuration, time-stepper, OSS, subscale, Chapter 5 and
physical-pressure tests pass; representative Galerkin, ASGS and OSS runs retain
their results to solver tolerance. No thermal options or classes should exist
without a consumer in the implemented phase.

B0 implementation: `FieldLayout` binds component functions, UFL expressions and
mixed tests by name, with the pressure index fixed at one. `flow_layout` selects
only the current Galerkin/ASGS/OSS fields. `StepSolver` checks all named changes
and exposes the existing log properties through the named iteration record.
No additional form/history helper is needed for this refactor; extract those
with their first thermal consumer in B1. Numerical equivalence and regression
results are recorded in `progress_log.md`.

### B1 — Monolithic Galerkin thermal baseline

Add physical data, scalar space, independent thermal boundary conditions,
buoyancy force, thermal Galerkin form, temperature histories and output. Start
with Newton, BDF1, Laplacian viscosity and the nonconservative/skew-symmetric
forms on periodic or closed domains. Use Taylor–Hood velocity/pressure.

Acceptance:

1. A diffusion-only analytic mode decays at the correct rate; a linear conduction
   state verifies source and outward-flux signs.
2. A vertically varying conduction state has $u=0$ and $\nabla p=b(T)$ when the
   temperature and hydrostatic pressure are representable. A hot perturbation
   accelerates against gravity.
3. With $\beta=0$, flow matches the isothermal solution and temperature behaves
   as a passive scalar. With no thermal data, legacy behavior is unchanged.
4. A smooth coupled manufactured solution uses
   $f=\partial_tu+u\cdot\nabla u-\nu\Delta u+\nabla p-b(T)$ and
   $s=\partial_tT+u\cdot\nabla T-\kappa\Delta T$; a directional Jacobian check
   detects missing $u\leftrightarrow T$ derivatives.

B1 implementation uses `ThermalProblem`, optional `ThermalOptions`, independent
`Dirichlet`/`DiffusiveFlux` scalar boundaries, and the shared `FieldHistory`
helper. `galerkin_thermal_system` composes the existing momentum kernel with
scalar-generic transport and resolved buoyancy. `manufactured_boussinesq` and
`configs/boussinesq/manufactured_galerkin.toml` exercise the single entry point.
Physical coefficients and resolved thermal choices are recorded in
`thermal_metadata.json`; temperature VTK, error/content/work diagnostics, and
named temperature changes are available. Each equation also has its own
initial-relative/absolute Newton check using configurable positive reference
scales; the normalized algebraic block norms and chosen scales are recorded. Unsupported thermal combinations fail
explicitly. The progress log records acceptance and regression results.

### B2 — Static linear ASGS and OSS

Add temperature residual/subscale, scalar stabilization and projection,
quasi-static pressure closure and temperature-subscale buoyancy feedback.
Use equal-order $P_1$ and $P_2$ spaces. Initially reject unimplemented dynamic
or nonlinear thermal settings explicitly.

Acceptance:

1. Static linear ASGS on affine $P_1$ gives the expected scalar SUPG term with
   the correct sign. A $P_2$ check activates both thermal Laplacian terms.
2. OSS thermal and velocity orthogonality hold with spatially varying effective
   parameters, unconstrained projection spaces and the common quadrature.
3. Temperature-subscale buoyancy vanishes in the matching-space OSS case, but
   is present with the correct sign for ASGS. The direct resolved buoyancy
   contribution cancels in the OSS velocity subscale under (B16)'s assumptions.
4. Turning the pressure subscale off removes only $\widetilde p,\zeta$; the
   thermal projection remains and pressure is still the sole gauge field.
5. With $\beta=0$, a passive scalar with data and diffusivity identical to
   temperature reproduces it to solver tolerance. With no passive data,
   temperature-only results and output are unchanged. All scalar fields have
   distinct projection and history storage.

B2 implementation: `thermal_system` composes the existing mechanical VMS kernels
with scalar-generic static residual/subscale/adjoint/projection kernels. OSS
fields are selected by `FieldLayout`, with `eta_temperature` and one
`eta_<name>` per passive scalar. The optional `Problem.passive_scalars` tuple
is normalized through `transported_scalars`; each scalar has independent data,
history, projection, output and named residual/change records. Passive content
and L² diagnostics expand to scalar-named columns. Both Newton parameter modes
are checked against directional derivatives; pressure projection can be removed
independently. Dynamic/nonlinear thermal choices remain rejected. Phase G4 slip
and partial-periodic mesh support are a separate prerequisite and are not
implemented merely by reaching B2.

### Boundary prerequisite G4 — Before the early B7.1 test

Implemented after B5 at the owner's request. `periodic_channel()` creates an
x-periodic affine triangle rectangle with bottom/top tags `y_min=3`, `y_max=4`;
translated boundary integrals confirm those tags and empty periodic tags 1/2 in
installed Firedrake 2026.10.0. The fully periodic helper keeps its existing
behavior. Any vertex remapping for B7 must precede space/BC construction.

`FreeSlip(normal_component)` imposes zero only on that axis component and leaves
tangential traction natural. Validation rejects invalid component indices and
walls whose normals are not aligned with the declared axis. The same BC builder
serves mixed solves and boundary-normal flux checks. Initialization enforces
only the slip-normal trace on resolved initial/exact-history states. Solver
residual norms and failed-state output use the component BCs consistently.

On flat horizontal walls, $w=0$ implies $\partial_x w=0$, so both Laplacian and
stress viscosity reduce to natural $\partial_z u=0$ for tangential motion.
Curved-wall slip remains outside this boundary type. Pressure retains its
constant null space when no outflow is present; only component 1 is gauged.
No mean-velocity constraint or OSS projection BC is introduced. Temperature
and passive insulation use their existing natural zero-flux descriptions.

The kinetic–potential budget now admits slip test functions: it checks zero
normal flow everywhere and zero velocity only on full-vector Dirichlet walls.
Tangential slip velocity is admissible without adding mechanical boundary work.
Both viscous forms pass isothermal preservation/decay checks; thermal Laplacian
runs verify T/q and budgets. The thermal stress-viscosity gate remains B6.

Acceptance passed: uniform preservation across Galerkin/eight VMS variants,
P1/P2 and all time schemes; translated tags, 3D component constraints, pressure
gauge/projection freedom, mixed-boundary flux rejection, failure output, and
nonuniform shear decay under both viscous forms. Serial/two-rank CLI examples
are in `configs/boussinesq/free_slip_uniform*.toml`.
See [G4 verification](verification/free_slip_g4.md). G4 does not implement the
stratified KH initial state or the B7 growth/cloud studies.

### B3 — Dynamic thermal subscales and all existing time schemes

Reuse scalar subscale histories, same-scheme integration, CN data averaging and
endpoint extrapolation. Extend startup validation to both exact physical fields.
Expose BDF1–3 and CN only after their thermal paths pass.

Acceptance: scalar subscale ODE tests verify effective parameters, history signs
and CN endpoint recovery; source-free temporal tests check the intended order
without contamination from inconsistent subscale startup; a forced test measures
the documented startup limitation. Check dynamic/static steady-state agreement,
temperature Dirichlet endpoint values and half-level buoyancy. Failed coupled
steps must not advance either field's histories.

B3 implementation uses scalar `SubscaleStorage` on the common quadrature rule,
matched dynamic effective parameters for both projection families, and the
existing resolved `FieldHistory` for every scalar. CN endpoint recovery includes
all dynamic scalar subscales. Exact-history validation covers velocity and all
scalars; zero subscale histories are explicitly retained with their documented
forced-startup limitation. Buoyancy work uses matching solve-level fields.
Thermal energy/momentum budget requests are rejected until B5, since the
isothermal budget formulas alone omit coupled work. The progress log records
source-free temporal orders, steady-state checks, failed-step preservation,
matched Jacobians, example runs and regressions.

### B4 — Nonlinear splitting and complete core variant matrix

Use the same $a=u_h+\widetilde u^*$ in both transports and both parameters.
Evaluate temperature subscales at the updated advector; extend convergence
records and residual checks. Adapt fixed-point and accelerated velocity updates
with the resolved buoyancy force included. Add Picard after Newton verification.

Acceptance: all eight thermal variants pass a small smooth coupled problem;
nonlinear thermal cross/subscale flux terms are present once; projections and
both subscale closure residuals converge together; Newton and Picard agree when
converged. Retain honest local-failure reporting for static nonlinear subscales
on coarse meshes, instead of masking the existing singular-update limitation.

B4 implementation passes an explicit shared advector into `thermal_system`,
reuses the existing pointwise velocity update with resolved buoyancy, and leaves
temperature/passive subscales algebraically eliminated. The outer solver samples
each scalar subscale before its solve/update, checks its change afterward, and
checks the relative velocity-closure residual. It reassembles the true coupled
FE residual after every local update and requires each scaled equation block to
pass its inner tolerance in addition to the combined outer residual criterion.
Nonlinear thermal outer CSVs name scalar/projection changes, scalar-subscale
changes, the velocity closure residual, and post-update FE block residuals;
isothermal and linear thermal columns retain their existing schemas.

Newton supports both parameter differentiation modes. Picard freezes the shared
total advector and parameters, preserving implicit buoyancy and every projection
unknown, including during optional warm-up before Newton. Tests independently
expand the nonlinear scalar weak terms, check diagonal scalar/velocity closure
and OSS orthogonality with dynamic history, compare both local update choices
and Newton/Picard solutions, and preserve explicit local-failure reporting.
Small examples live in `manufactured_nonlinear_oss.toml` and
`manufactured_nonlinear_asgs_picard.toml`; they are integration checks rather
than the B5 production validation or B7 shear-layer study.

### B5 — Conservation, verification and natural-convection validation

Implement the thermal and buoyancy diagnostics, including scheme-correct energy
budgets, and the conservative thermal form. Verify (B17) and (B18) only under
their stated assumptions. Include the pressure-subscale terms in mechanical
budgets and the higher-degree thermal diffusion cross term.

Add the verified G3 scalar-content identity (B19) for every splitting/form
combination, including the dynamic ASGS subscale moment. Implement potential
energy and the full tested kinetic–potential budget (B20)–(B24): preserve dynamic
ASGS cross energy, static ASGS mixed time rates, momentum-convection defects,
source and diffusion/boundary work. Separately output the linear-splitting
buoyancy remainder, nonlinear subscale buoyancy remainder and scalar-convection
defect. Only conservative **linear OSS**, with the stated matching spaces, has
zero buoyancy exchange remainder in general; do not implement the rejected
linear-ASGS exact-exchange hypothesis as a test expectation.

Budget fixtures need horizontal periodicity or a closed domain where
$\phi=g\cdot x$ is admissible. Do not use fully periodic gravity directions or
Dirichlet thermal walls for this simple test. Retain actual BDF/CN histories and
report a residual normalized by the sum of absolute work/rate terms, with an
absolute floor, so near-cancelling work does not create a misleading relative
error. These new diagnostics do not silently widen the existing isothermal
`energy_budget`, which currently permits dynamic CN only. G5 is optional: add
signed global/cellwise stabilization-work ratios with denominator validity masks,
using the definitions and qualifications in formulation §7.3.

Run spatial refinement for coupled manufactured $P_1$ and $P_2$ cases and temporal
refinement for supported schemes. Under smoothness and stability assumptions,
target the existing flow convergence rates and temperature $L^2$ order $k_T+1$,
$H^1$ order $k_T$; estimate rates over a refinement sequence instead of accepting
one low error. Include a 3D smoke test and serial/two-rank agreement for global
thermal diagnostics and orthogonality, allowing documented quadrature effects.

Add a differentially heated closed cavity: hot/cold vertical walls, insulated
remaining walls, no-slip velocity. Start with moderate Rayleigh numbers and
check mesh convergence, circulation orientation and consistent hot/cold heat
fluxes. Before asserting agreement with published Nusselt numbers, obtain and
record a primary benchmark's exact geometry, nondimensionalization, wall data,
Prandtl number and reported quantities.

Only then attempt the 2010 §4.1 aspect-ratio-eight unsteady cavity
($Ra=3.45\times10^5,10^6,10^7$). That paper uses bilinear quadrilateral elements,
BDF2 resolved fields, BE subscales and a particular nondimensional time step.
Our different mesh/parameter/time choices preclude calling a direct figure
comparison an exact reproduction. Consult its cited Christon et al. benchmark
for missing physical details before creating that case. Turbulent cavity/channel
comparisons are subsequent studies, not the initial correctness gate.

B5 implementation keeps stateless field diagnostics in `thermal_diagnostics.py`
and the captured, history-dependent budgets in `thermal_budgets.py`. No legacy
isothermal budget is repurposed. New CSV groups expose resolved/total scalar
moments, variance with ASGS cross storage and higher-order diffusion, full
kinetic–potential exchange, finite-iterate/projection terms, and molecular versus
variational boundary fluxes. Invalid constant/potential tests are rejected.
Signed rate/work sums use a recorded absolute normalization floor.

Solved tests cover the 48 core scalar/momentum-form combinations, all schemes
and P1/P2 histories, nonzero reference temperature, prescribed flux, Dirichlet
reactions, periodic geometry, 3D, and serial/two-rank agreement. Spatial and
buoyancy-driven temporal refinement passed. The square-cavity builder uses
explicit width-based free-fall scaling; moderate-Rayleigh studies distinguish
volume, wall-gradient and wall-reaction Nusselt estimates against the primary
reference. Only after those checks were the 2010 tall-cavity smoke cases run.
Their short, coarser triangular calculations are not long-time benchmark
reproductions. Details and measured results are in
[the B5 report](verification/boussinesq_b5.md). G5 remains optional and unimplemented.

### B6 — Compatibility with the remaining framework options

Extend and test stress viscosity, the remaining momentum convection forms,
EMAC/EMAC–VMS, and open thermal boundaries. For EMAC preserve
$p=P+|u_h|^2/2$ and its existing residual convention. There is no automatic scalar
EMAC counterpart: choose the thermal convection form independently and derive
its heat/variance budget. Add external buoyancy force and torque to momentum
budgets; convection conservation is not conservation of a forced flow.

Acceptance: every newly enabled combination has form-level consistency and
budget checks, plus a coupled integration case. Open boundaries include signed
advective heat flux, prescribed diffusive flux, temperature inflow data and an
explicit backflow policy. Until then, validation rejects unsupported thermal
combinations while leaving existing isothermal options available.

### B7 — Stratified Kelvin–Helmholtz shear layer and moisture visualization

This application follows B5 and the G4 boundary gate; it does not require all of
B6. Run only momentum forms already verified for thermal problems. B7.1 alone
may run after B2/G4 as a small integration smoke test. The addendum also labels
B7.3 as early work; this conflicts with its explicit “only B7.1” exception, so
the full B7.3 study stays after B5. No production run is authorized by this
documentation change.

#### B7.0 — Case definition and reference data

Implement `kelvin_helmholtz_shear_layer` using formulation §8, including recorded
dimensionless scaling and boundary-compatible profiles. Store reference
conventions and actual profile/cutoff definitions under
`data/references/kelvin_helmholtz/README.md`. Use the single existing simulation
entry point and configs under `configs/boussinesq/`.

Use the following proposed starting values, not a claim of validated resolution:

| Parameter | Starting value | Study variation |
| --- | --- | --- |
| `reynolds`, `prandtl`, `schmidt` | 1000, 1, 1 | Re = 2000 |
| `richardson` | 0.1 | 0, 0.3 |
| `n_waves`, `wave_number`, `height` | 1, 0.4446, 20 | Four waves for moisture |
| `nx_per_wave`, `nz` | 128, 192 for P1 | 1.5× refinement; half counts for a P2 comparison only |
| `clustering` | 0 | Positive sinh clustering after uniform-mesh validation |
| `dt`, `t_end` | 0.05, 150 | Half dt check; t = 200 for moisture |
| `epsilon`, `epsilon_noise` | 1e-4, 0 | Noise amplitude 1e-2 relative to epsilon for four waves |
| `rh_peak`, `moisture_width` | 0.97, 1.5 | Record initialization near the walls |
| `saturation_height_rate`, `saturation_temperature_rate` | 0.03, 0 | Height rates 0.02, 0.05 |
| Stabilization | Dynamic nonlinear OSS | Budget-study variants below |

Use descriptive case-parameter names and validate their ranges individually;
zero Richardson number, zero clustering and zero noise are valid. Random seed,
Fourier mode list/normalization and wall cutoff must be fixed and recorded.
The stated mesh gives advective CFL of order 0.5 initially; compute actual
maximum CFL throughout the run, especially with clustering and P2. Halving the
P2 cell count is a comparison choice, not an accuracy-equivalence guarantee.

#### B7.1 — Early coupled channel smoke test

After B2 and G4, build a coarse periodic/free-slip case with temperature and
optional duplicate passive scalar. Check profile/boundary compatibility,
hydrostatic initialization, finite forms/Jacobians, pressure gauge, projection
orthogonality, scalar names and initial divergence/net flux. Evolve a few static
linear ASGS/OSS steps. This exercises the coupling but establishes no instability
growth rate, cloud pattern or long-time conservation result.

Implementation: `cases/kelvin_helmholtz_shear_layer.py`, the two
`kh_smoke_static_*.toml` configs and the case reference definitions above.
The B7.0 definitions needed by this test are included. At B7.1 the builder
accepted uniform mesh counts through the existing `[mesh].cells`; optional
clustering is added at B7.2 below. Initial physical pressure is optional problem data; old cases
retain their zero-pressure initialization and original output schema.
Hydrostatic interpolation defects are recorded explicitly. Per-step vertical
energy and quadrature-sampled degree-weighted advective CFL are opt-in diagnostics.
See [B7.1 verification](verification/kelvin_helmholtz_b71.md).
Stop here before B7.2 startup/refinement work or the later application studies.

#### B7.2 — Mesh, startup and reference preparation

Check the uniform mesh first; test positive Jacobians and aspect ratios before
enabling vertex clustering. Verify exact normal-wall compatibility after FE
interpolation, zero scalar wall flux and the recorded base-profile adjustments.
Keep subscales initially zero. Compare initialization transients and monitor
the horizontal-mean profiles rather than treating tanh as a steady viscous
solution. A no-perturbation control distinguishes numerical/hydrostatic startup
motion from the intended tiny perturbation.

Implemented B7.2: nonnegative `clustering` moves vertices before space/measure
construction. Ordered rows guarantee positive vertical-map Jacobians; MPI-global
cell-measure, map-Jacobian, edge and aspect-ratio ranges are recorded. Clustered
meshes must retain a complete cell inside each wall plateau. The original uniform
path and zero-clustering limit remain available.

Opt-in `output.horizontal_profiles_every` writes endpoint velocity, temperature
and each tracer mean on all horizontal rows/midrows, including initialization
and the final/early-stop state. Gauss integration splits each horizontal cell
at its midpoint so it respects triangle intersections at these heights.
Disabled output preserves the prior output/config schema. An analysis-only
script compares completed pairs and exports the measured drift and startup
energy; it does not integrate or fit a growth rate.

Six zero-seed/seeded pairs check uniform refinement, half dt, P2 and clustering.
The zero-seed vertical energy falls with refinement but remains 98–99.98% of
the seeded energy at t=0.5 on the tested meshes. These runs establish **no clean
growth-fitting window**. Initial pressure is only a solve guess; changing that
guess alone does not remove the converged discrete startup motion. Keep the
formulation, zero subscale histories and seed amplitude unchanged. Before B7.3
fits a rate, resolve this contamination through further reference-resolution
checks and seeded-mode/control diagnostics; do not interpret a slope of the raw
vertical energy as KH growth. See [B7.2 verification](verification/kelvin_helmholtz_b72.md).

#### B7.3 — Linear growth study

Run one wave for Ri = 0, 0.1, 0.3 and Re = 1000, 2000, with the reference mesh,
1.5× refinement and a timestep refinement. Fit
$\sigma=\frac12 d\log E_w/dt$, $E_w=\frac12\|w_h\|^2$, only over a demonstrated
linear window after initialization and before nonlinear roll-up/base drift.
Record window endpoints, fit residual/uncertainty and base-profile change;
also track the seeded Fourier component to distinguish it from transient modes.

For Ri = 0 compare against the unstratified inviscid reference near 0.1897–0.1898
in the specified units, accounting for finite walls and viscosity. Re = 1000
need not satisfy a prescribed distance below that value. For Ri = 0.1 measure
whether growth is positive and reduced, comparing with a stability calculation
for the **actual** stratified profile if a quantitative reference is needed.
For Ri = 0.3 test for absence of a sustained growing KH eigenmode in the
appropriate small-amplitude regime; do not require monotonic $E_w$ decay from
Miles–Howard. Oscillations and transient growth are admissible observations.

Use 2% as a proposed discretization-convergence target for clearly positive
growth rates, with an absolute/uncertainty criterion near zero. Report scheme
spread separately: agreement within 2% is a study target, not a theorem about
underresolved VMS variants. A discrepancy triggers resolution/reference checks
and reporting, not tuning the formulation to the expected rate.

#### B7.4 — Dry nonlinear billows and budget matrix

Use one wave at Ri = 0.1, Re = 1000 to t = 150. Report roll-up, saturation and
decay if observed, plus vorticity/temperature snapshots at maxima of $E_w$ and
$E_p$ (with the chosen potential-energy datum). Study

```text
thermal convection {conservative, skew-symmetric, nonconservative}
× momentum convection {already verified thermal forms}
× {ASGS, OSS} × {linear, nonlinear} × {static, dynamic}
```

Estimate cost before the full matrix; begin on a documented coarse mesh and
retain refinement checks on representative combinations. Record local/outer
solver failures, including known static-nonlinear limitations. They are findings,
not missing runs to conceal.

For each successful combination report the instantaneous and accumulated
normalized (B24) residual; all three (B22) remainder terms separately; scalar
content against (B19); variance against (B18) only where applicable; and actual
time-rate, viscous, pressure-subscale and momentum-convection contributions.
Do not attribute every budget discrepancy to G2.2/G2.3: the corrected linear-ASGS
remainder and other momentum/time terms are also required. Expect the full
tested-equation residual, rather than every individual defect, to approach
solver/quadrature accuracy. Report any unexplained remainder and combination.

#### B7.5 — Four-wave moisture visualization

After B5 and the dry study, run four waves with fixed-seed noise, Ri = 0.1,
Re = 1000, to t = 200. Use the project defaults and add conservative thermal
transport as a comparison if the budget study shows a material effect. Repeat
the saturation-height rates 0.02, 0.03, 0.05; these are diagnostic sensitivity
values, not guarantees that one run saturates the whole layer.

Write VTK for velocity, pressure, temperature, moisture, cloud proxy and vorticity;
optional G5 fields retain their solve-level time labels. Produce an animation of
the cloud proxy over temperature contours. CSV records $E_w$, appropriate $KE$,
$E_p$, resolved scalar integrals, budget scalar integrals including dynamic-ASGS
subscales when relevant, maximum cloud proxy, tracer extrema and cloud fraction
with raw threshold $q_c>10^{-6}$. Include optional valid G5 ratios.

Acceptance is scalar-budget closure per (B19), resolved spatial/time sensitivity,
and reproducible diagnostics. Report whether lifted moist regions form, whether
troughs remain cloud-free, the number/spacing of regions, roll-up time, and any
vortex pairing. Do not force four bands, pairing, or whole-layer saturation by
adjusting data or clipping moisture. Document the passive-tracer and saturation
proxy limitations. The 2D result is a cross-section, not a resolved spanwise cloud
field; absence of the proposed morphology is a legitimate result.

#### B7.6 — Optional extensions, separately approved after B7.5

Latent-heating feedback requires a new physical formulation, not just an edit to
`boussinesq.py`: define total water versus vapor, liquid-water potential
temperature, dimensionless latent coefficient, saturation adjustment and its
sources/energy budget. A nonlinear buoyancy function also invalidates the simple
matching-space OSS cancellation (B16) in general and requires new Jacobian and
projection checks. A smooth positive-part law has a nonzero tail and needs a
documented smoothing-sensitivity study. No such feedback is part of B7.5.

A 3D spanwise-periodic extension needs agreed resolution, run-cost/core-count
estimates, boundary construction and a distinct secondary-instability study.

## 5. Completion criteria

The core thermal extension is ready after B0–B5 when the eight variants, four
time schemes, Newton/Picard, manufactured convergence and moderate cavity
validation pass; B6 completes compatibility with the other framework options.
B7 adds the separately gated shear-layer and passive-moisture application.
The default is dynamic nonlinear OSS, using shared switches and the current
additive stabilization scales. Each implemented phase updates documentation and
the progress log, passes its targeted checks and the relevant existing regression
suite, and reports its supported combinations explicitly.

No new time integrator, turbulence closure, full off-diagonal local buoyancy
model, variable-property physics, or free-surface wave system is needed for
this extension. Those would be separate formulations with their own validation.
