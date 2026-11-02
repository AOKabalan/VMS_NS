# Review of the budget and Kelvin–Helmholtz addendum

2 November 2026. Scope: verify the review against the actual flow code and
[proposed thermal equations](../boussinesq_formulation.md), then update the
[integration plan](../boussinesq_implementation_plan.md). No production solver,
case, configuration or diagnostic has been implemented in this change.

The architectural suggestions are adopted. G3 is correct under its stated
space/boundary assumptions, with explicit ASGS storage qualification. G2.1 is
false for general linear ASGS; G2.2 needs a sign convention and the nonlinear
subscale contribution; G2.3's nonlinear conservative remainder is correct.
Several KH outcome requirements overstate what theory or the proposed moisture
proxy guarantees. The equations and global defaults have not been altered to
make these hypotheses pass.

## 1. Item-by-item disposition and code placement

Source paths in this table are relative to `src/navier_stokes_vms/`.

| Item | Decision | Where it lands and current limitation |
| --- | --- | --- |
| G1: generic scalar transport | Adopt | Generic kernels in planned `formulation/thermal.py` and scalar OSS kernel; buoyancy stays in `boussinesq.py`. `Problem` gains an optional named passive tuple, backed by `problem/scalar_problem.py`. Named layouts, histories, solver records and output handle each scalar. None of these thermal/passive paths exists yet. Introduce a second-scalar equality test in B2. |
| G2: potential energy | Adopt the corrected derivation | Formulation §7.2, (B20)–(B24); `postprocessing/thermal_diagnostics.py`, existing conservation context and B5 budget tests. Current `diagnostics.py` already uses cross energy for dynamic ASGS; its total energy closure is restricted to dynamic CN. |
| G2.1: exact linear exchange for both ASGS and OSS | Reject as stated | Conservative linear exchange leaves $-\beta(T_h,g\cdot\widetilde u)$. Matching-space OSS cancels it; ASGS generally does not. Dynamic ASGS potential subscale storage is necessary but insufficient to remove this remainder. |
| G2.2: form defect | Adopt with corrections | Use $D_\phi=(\nabla\cdot u_h,T_h\phi)-\chi_N(\widetilde u,\nabla(T_h\phi))$, $\phi=g\cdot x$. Nonconservative minus conservative convection is $-D_\phi$ on the equation's left; the energy-budget right gains $+\beta D_\phi$. Skew gives half. The nonlinear term cannot be omitted. |
| G2.3: nonlinear closure defect | Adopt with qualifications | The conservative nonlinear exchange remainder is $+\beta(\widetilde T,g\cdot\widetilde u)$ on the combined budget's right. Other momentum, time, diffusion and source terms still belong in the full budget. Restoring off-diagonal buoyancy would change the local model and projection; it is not implemented or planned as a hidden repair. |
| G3: scalar content | Adopt | Formulation §7.1, (B19), generic diagnostics and B5. Nonlinear splitting gives zero convection contribution for all three scalar forms under the admissibility assumptions. Dynamic ASGS conserves the resolved-plus-subscale moment, not necessarily the resolved integral alone. |
| G4: channel and slip walls | Adopt as an explicit prerequisite | `problem/mesh_loading.py`, `boundary_conditions.py`, problem checks and `TimeStepper` BC handling. Installed Firedrake supports x-periodicity/component BCs; repository wrappers support fully periodic rectangles and full-vector Dirichlet/outflow only. Thermal insulation awaits B1. |
| G5: effective dissipation ratios | Optional, with precise interpretation | Formulation §7.3 and optional B5 diagnostics. Use signed stabilization-work proxies with validity masks; include reversible storage/exchange and do not equate every negative value to physical backscatter. |
| B7: dry KH and moisture visualization | Adopt with revised gates | Formulation §8, planned `cases/kelvin_helmholtz_shear_layer.py`, registry entry, reference metadata, configs under `configs/boussinesq/`, existing entry point, and new B7 in the plan. Cloud morphology, pairing and transient stability are observations to measure, not outcomes to force. |

An empty passive-scalar tuple preserves temperature-only layout/output; no
thermal data preserves existing isothermal defaults and output. The expanded
OSS layout is $u,p,T,q_m,\xi,\eta_T,\eta_q,(\zeta)$ with pressure at component 1.
The current positional interpretation of fields after $(u,p)$ is therefore a
known B0 refactor requirement, not a capability the existing code already has.

## 2. Algebraic corrections that affect acceptance criteria

Adding the resolved momentum test and velocity-subscale test yields

$$
M_K=(D_tu_h,u_h)+\chi_A\chi_D(D_t\widetilde u,u_h)
 +\chi_D(D_t\widetilde u,\widetilde u)
 +\chi_A(D_tu_h,\widetilde u).
$$

Consequently, dynamic ASGS uses $\frac12\|u_h+\widetilde u\|^2$; dynamic OSS
uses the sum of the two scale energies; static ASGS retains a mixed rate that
is not the derivative of either energy. This agrees with the current flow
specification and `Conservation.__init__` in `postprocessing/diagnostics.py`.
The addendum's resolved-plus-subscale-square prescription loses the ASGS cross
term. Equations (B21)–(B24) retain all BDF/CN history pairings explicitly.

At zero reference temperature, kinetic buoyancy work is
$-\beta[(T_h,g\cdot u_h)+(T_h,g\cdot\widetilde u)
+(\widetilde T,g\cdot u_h)]$. Conservative potential advection is
$+\beta(T_h+\widetilde T,g\cdot(u_h+\chi_N\widetilde u))$.
Their sum is exactly the first two terms in (B22), making the missing
linear-ASGS term visible without changing any equation.

For G3, use continuity with pressure test $T_h$, integrate only the continuous
resolved velocity by parts, and retain subscale pairings without differentiating
stored subscales. This gives the proposed convection expressions. With
nonconverged continuity residual
$C(T_h)=(T_h,\nabla\cdot u_h)-(\widetilde u,\nabla T_h)$, the nonconservative
expression additionally contains $-C(T_h)$; skew has half. This is a useful
way to separate algebraic solver error from the predicted defect. Pressure
subscale stabilization does not enter that continuity test.

Potential energy requires more precise geometry than “closed or periodic”:
$g\cdot x$ is not a periodic FE function along gravity. An x-periodic channel
with vertical gravity is valid. Periodicity also permits nonzero mean velocity,
so $\int u=0$ cannot be assumed there. With admissible pressure test $g\cdot x$
and impermeable walls, continuity instead gives $(g,u_h+\widetilde u)=0$.
Thus nonzero-reference cancellation in the total budget concerns this **total**
moment, not just $(g,u_h)$. Include source work $\beta(s,g\cdot x)$, absent
from the proposed continuous formula unless $s=0$ is imposed.

No general scalar convection choice simultaneously removes the variance,
potential-exchange, diagonal-closure and momentum-convection contributions.
Keep the global skew default and compare scalar forms in stratified studies.

## 3. KH case corrections

- **Stability criterion:** $Ri_g\geq1/4$ is a sufficient modal-stability
  condition for the relevant inviscid steady parallel-flow problem. It does
  not imply monotonic vertical-energy decay, or instability whenever the
  inequality is violated. Finite viscosity, base-state diffusion and transient
  amplification require separate checks. See
  [Miles (1961)](https://doi.org/10.1017/S0022112061000305) and the primary study
  of [stratified transient growth](https://doi.org/10.1017/jfm.2017.396).
- **Reference growth:** $k_0=0.4446$ and growth near 0.1897–0.1898 are appropriate
  unstratified tanh-layer references in the stated half-jump velocity units.
  Michalke's approximately 0.0949 uses half that velocity amplitude. An
  independent published Rayleigh calculation documents this convention in
  [Callaham et al., §6.2](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/multiscale-model-reduction-for-incompressible-flows/FC550C3A53CEB7174759D0C45C3DFDBB).
  Neither the optimal wavenumber nor the precise growth rate carries unchanged
  to stratified finite-Re flow. A 2% mesh/time target is useful; unconditional
  2% agreement among underresolved stabilization variants is not established.
- **Walls and initialization:** tanh has $\operatorname{sech}^2(10)\simeq
  8.24\times10^{-9}$ wall derivative, not exact zero; a Gaussian streamfunction
  also has nonzero wall trace. Require explicitly documented boundary-compatible
  cutoffs/flattening for exact budget tests, and measure the resulting base
  profile. Hydrostatic/startup motion can overwhelm an $\epsilon=10^{-4}$ seed;
  add an unperturbed control and fit-window/base-drift checks.
- **Moisture:** G3 applies to the total dynamic-ASGS scalar moment. At the
  proposed default nonlinear OSS it applies to the resolved integral directly.
  The passive tracer and diagnostic positive part model no condensation sink
  or latent heating. For $\alpha=0$, saturation is independent of temperature.
  The displacement thresholds in formulation §8 explain possible crest
  activation but prove no particular cloud pattern or whole-layer saturation.
  In fact, for the proposed $\alpha=0$ data the continuum tracer maximum is
  below 0.972 whereas saturation is at least 1 for $z\leq0$; saturation of the
  entire fixed channel cannot follow from source-free maximum-preserving
  transport. Define the region meant by “whole layer” before measuring it.
  Report negative/overshooting numerical moisture rather than clipping it.
- **Optional feedback:** replacing temperature by a nonlinear function of two
  prognostic scalars changes the physics, budgets and generally the simple OSS
  buoyancy cancellation. It is not confined to changing one force routine.
  Retain a separate formulation/approval gate after B7.5, as for the 3D study.
- **Phase ordering:** the addendum says only B7.1 can precede B5, but later
  labels B7.3 as early too. The plan explicitly retains the narrower B7.1 smoke
  exception; the full growth/budget/cloud studies follow B5 and free-slip support.

## 4. Checks performed during this review

1. Read the current momentum/residual/convection kernels, time-scheme
   coefficients, conservation diagnostics, pressure-null-space logic, mesh and
   boundary builders, and the existing flow conservation specification.
2. Re-derived (B19)–(B24) from (B10), (B13)–(B15), including static ASGS
   resolved time derivatives and all dynamic history pairings. Compared the
   interpretation of optional ratios with the local 2010 PDF §3.3, (49)–(51).
3. Assembled an independent identity probe on an x-periodic affine triangular
   mesh with $P_2$ resolved fields and degree-six scalar/vector quadrature
   subscales. It used the existing momentum convection kernels and actual
   `step_coefficients`/`time_derivative`, with distinct nontrivial field histories.
   Across 96 combinations (four schemes, three scalar forms, two spaces,
   two dynamics choices, two splittings), the maximum combined-budget difference
   was $4.0\times10^{-15}$; the content identity including continuity residual
   differed by at most $2.7\times10^{-15}$. This checks assembled algebra on
   arbitrary fields; it does **not** constitute a solved thermal problem or a
   test of OSS projection convergence. A sample linear-ASGS remainder was
   nonzero, about $-0.237$, confirming it is not an algebraic zero.
4. Ran a boundary probe using installed Firedrake 2026.10.0 on a $2\times3$
   rectangle: tags 1/2 have zero measure; tag 3 has length 2 at height 0;
   tag 4 has length 2 at height 3. Applying the vertical-component BC on tags
   3/4 preserved uniform horizontal flow with squared error below $10^{-30}$.
   This verifies boundary construction, not a time-dependent free-slip solve.
5. Independently solved the unstratified inviscid Rayleigh eigenproblem for
   $U=\tanh z$, $z\in[-10,10]$, $k=0.4446$ using Chebyshev collocation and
   homogeneous perturbation wall data. Growth rates for 96, 144 and 192
   intervals were 0.1895231, 0.1896382 and 0.1896404. This supports the scaling
   and shows why finite-wall/reference differences must be distinguished from
   a numerical method's error. It does not validate the planned viscous KH case.

The runtime probes were temporary review calculations, not additions to the
production test suite. Full thermal implementation tests remain the gates in
the plan. Documentation links, equation references and Markdown delimiters are
checked separately before handoff.
