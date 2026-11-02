# Thermal Boussinesq VMS formulation

Status: B5 budgets and G4 periodic/free-slip boundary support implemented,
2 November 2026. This document specifies the thermal
buoyancy model requested for the existing incompressible Navier–Stokes framework;
The implemented subset includes Taylor–Hood Galerkin and equal-order static
or dynamic, linear or nonlinear ASGS/OSS with temperature and optional passive scalars.
Newton and Picard, BDF1–3 and Crank–Nicolson, Laplacian viscosity and the nonconservative/
skew choices, plus conservative scalar transport, are supported on periodic or
impermeable domains. Coupled budget diagnostics and moderate-Rayleigh cavity
validation are implemented; remaining options and applications follow later.
The companion
[integration plan](boussinesq_implementation_plan.md) identifies the code changes
and verification gates. Existing flow conventions remain defined by
[formulation_spec.md](formulation_spec.md).
The [addendum review](reviews/boussinesq_addendum_review.md) records the accepted
suggestions, corrected hypotheses and checks behind the Kelvin–Helmholtz extension.

## 1. Sources and modelling decisions

The primary source is Codina, Principe and Ávila (2010), *Finite element
approximation of turbulent thermally coupled incompressible flows with numerical
sub-grid scale modelling*, IJNMHFF 20(5), 492–516:
[local PDF](papers/codina2010.pdf), [DOI](https://doi.org/10.1108/09615531011048213).
Equation numbers below refer to this published PDF, whose first page is a publisher
cover sheet.

The other requested source, Codina, González-Ondina, Díaz-Hernández and Principe
(2008), *Finite element approximation of the modified Boussinesq equations using
a stabilized formulation*, IJNMF 57, 1249–1268, concerns **dispersive water waves**:
[local PDF](papers/Codina_et_al-2008-IJNMEinFluids.pdf),
[publisher record](https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.1718).
Its unknowns are water elevation and horizontal velocity, not temperature and
incompressible pressure. Its §4 supports the general residual/subscale construction
and dimensional scaling of stabilization; §4.3.2, Remark 3 explains the removal of
resolved time derivatives under orthogonal projection. Its wave stabilization
parameter, dispersive auxiliary fields and Adams–Moulton predictor–corrector are
not part of this thermal extension. This paper is also distinct from the
single-author “Codina 2008” cited in the existing flow specification.

| Item | Source | Proposed choice here |
| --- | --- | --- |
| Incompressible momentum, continuity and temperature | 2010 (1)–(10) | Constant material properties, density variation only in buoyancy |
| Velocity, pressure and temperature splitting | 2010 (12)–(17) | Add a scalar temperature subscale to the existing split |
| Local closure | 2010 (26), (28) | Diagonal velocity/temperature relaxation, dynamic or static |
| Pressure subscale | 2010 (27) contains a time derivative of the continuity residual | Retain the repository's quasi-static pressure closure; this is a deliberate difference |
| Spatial parameters | 2010 (23)–(25) use a square root of sums of squares | Retain the repository's additive inverse time scales, extended to thermal diffusion |
| Convection | 2010 (29)–(39) | Default skew-symmetric form as in the repository; also specify nonconservative and conservative scalar forms |
| OSS | 2010 (32)–(34) | Consistent, effective-parameter-weighted projections to preserve orthogonality with pointwise parameters |
| Time integration | 2010 §2.5 allows subsequent time discretization; §4.1 uses BDF2 for resolved fields and BE for subscales | Reuse the repository's same-scheme integration of all dynamic fields |
| Variant switches | 2010 emphasizes dynamic nonlinear OSS | Extend all eight existing ASGS/OSS × static/dynamic × linear/nonlinear variants; this complete matrix is a project choice |

Thus this is a thermal extension of the current method, informed by the paper,
rather than a claim of exact reproduction of its numerical experiments.

## 2. Continuous problem and sign convention

On a fixed domain $\Omega\subset\mathbb R^d$, $d=2,3$, let $T$ be temperature,
$T_{\rm ref}$ a spatially and temporally constant reference temperature, $\nu>0$
kinematic viscosity, $\kappa>0$ thermal diffusivity, $\beta\geq0$ thermal expansion
coefficient, and $\boldsymbol g$ the gravitational acceleration vector pointing
downwards. Set

$$
\boldsymbol b(T)=-\beta(T-T_{\rm ref})\boldsymbol g.
\tag{B1}
$$

The equations are

$$
\begin{aligned}
\partial_t\boldsymbol u+(\boldsymbol u\cdot\nabla)\boldsymbol u
-\nu\Delta\boldsymbol u+\nabla p
&=\boldsymbol f+\boldsymbol b(T),\\
\nabla\cdot\boldsymbol u&=0,\\
\partial_tT+\boldsymbol u\cdot\nabla T-\kappa\Delta T&=s.
\end{aligned}
\tag{B2}
$$

Here $p$ is kinematic pressure after subtracting the constant-density hydrostatic
part, $\boldsymbol f$ is acceleration excluding thermal buoyancy, and $s$ has units
of temperature per time. A volumetric heating rate must be divided by $\rho_0c_p$
before being supplied as $s$; thermal conductivity must likewise be divided by
$\rho_0c_p$ to obtain $\kappa$. With $\boldsymbol g=(0,-g_0)$ and $T>T_{\rm ref}$,
the buoyancy force points upwards. The approximation assumes
$|\beta(T-T_{\rm ref})|\ll1$; variable density in inertia or continuity is outside
this model.

Initial data are $\boldsymbol u_0,T_0$. Retain the existing velocity boundary
conditions. Independently partition the thermal boundary into prescribed
temperature $\Gamma_T$ and prescribed outward diffusive flux $\Gamma_j$:

$$
T=T_D\quad\text{on }\Gamma_T,\qquad
-\kappa\nabla T\cdot\boldsymbol n=j_{\rm out}\quad\text{on }\Gamma_j.
\tag{B3}
$$

Insulation means $j_{\rm out}=0$; $j_{\rm out}$ has units of temperature times
velocity, and multiplication by $\rho_0c_p$ gives physical heat flux. Thermal
inflow requires temperature data. Initially verify periodic domains and closed
cavities; open thermal boundaries need explicit inflow/backflow treatment.

The existing stress viscosity option replaces $\nu\Delta\boldsymbol u$ by
$\nabla\cdot(2\nu\nabla^s\boldsymbol u)$ consistently in weak and strong
operators. The first thermal implementation should use the existing
skew-symmetric or nonconservative momentum form; EMAC pressure conventions and
conservation budgets require the later compatibility gate in the plan.

### Dimensionless interpretation

For scales $L,U,\Delta T$, set $\theta=(T-T_{\rm ref})/\Delta T$. The dimensionless
diffusion coefficients are $1/Re$ and $1/Pe$, and the force is
$-Ri\,\theta\widehat{\boldsymbol g}$, where

$$
Re=\frac{UL}{\nu},\quad Pe=\frac{UL}{\kappa}=Re\,Pr,\quad
Pr=\frac{\nu}{\kappa},\quad
Ra=\frac{\beta|\boldsymbol g|\Delta T L^3}{\nu\kappa},\quad
Ri=\frac{\beta|\boldsymbol g|\Delta T L}{U^2}
=\frac{Ra}{Re^2Pr}.
\tag{B4}
$$

With the thermal diffusion velocity $U=\kappa/L$, viscosity is $Pr$, thermal
diffusivity is $1$, and buoyancy magnitude is $Ra\,Pr$. Case builders must select
and record one scaling; the formulation receives $\nu,\kappa,\beta,\boldsymbol g$
directly and never guesses whether a run is dimensional.

## 3. Spaces, splitting and residuals

Use $\boldsymbol V_h\times Q_h\times S_h$ for velocity, pressure and temperature.
For VMS start with continuous equal-degree $P_k/P_k/P_k$, $k=1,2$. The Galerkin
baseline uses Taylor–Hood velocity/pressure and continuous temperature of the
velocity degree. Put

$$
\boldsymbol u=\boldsymbol u_h+\widetilde{\boldsymbol u},\quad
p=p_h+\widetilde p,\quad T=T_h+\widetilde T,\qquad
\boldsymbol a=\boldsymbol u_h+\chi_N\widetilde{\boldsymbol u}.
\tag{B5}
$$

$\chi_N$ is one for nonlinear splitting and zero for linear splitting. Also let
$\chi_D$ select dynamic subscales, and $\chi_A$ select ASGS (zero for OSS).
Temperature and velocity share these switches initially. The subscales vanish
on element boundaries in the modelling derivation; their quadrature storage is
not a differentiable finite element field. No implementation may take spatial
derivatives of stored subscales.

Additional passive scalars $C_i$ use the same scalar spaces, transport forms,
stabilization switches and time integration as temperature, with their own
diffusivity, source, initial/boundary data, projection and subscale histories.
Equations (B6), (B8)–(B12), (B14)–(B15) apply with
$(T,\kappa,s,\widetilde T,\eta)$ replaced by
$(C_i,\kappa_i,s_i,\widetilde C_i,\eta_i)$. Their kernels contain no buoyancy.
Only $T$ enters (B1), (B13); solve passive scalars in the same monolithic system.
For moisture, write $q_m$ in the derivations to distinguish it from the pressure
test $q_h$; its output name can be `q`.

With the repository's negative-operator residual convention,

$$
\begin{aligned}
\boldsymbol R_u&=\boldsymbol f+\boldsymbol b(T_h)-D_t\boldsymbol u_h
 -(\boldsymbol a\cdot\nabla)\boldsymbol u_h
 +\nu\Delta_h\boldsymbol u_h-\nabla p_h,\\
R_p&=-\nabla\cdot\boldsymbol u_h,\\
R_T&=s-D_tT_h-\boldsymbol a\cdot\nabla T_h+\kappa\Delta_hT_h.
\end{aligned}
\tag{B6}
$$

The Laplacians are elementwise, zero for affine $P_1$ but retained for higher
degree. For other supported momentum forms, use the existing
`momentum_residual` and its pressure/viscous conventions with the resolved force
$\boldsymbol f+\boldsymbol b(T_h)$.

### What is approximated in the buoyancy coupling

The exact split subscale momentum equation contains
$+\beta\boldsymbol g\widetilde T$ on its left-hand side [2010 (15)]. The paper's
final **diagonal** closure [2010 (26)–(28)] approximates the spatial block operator
by relaxation parameters and does not retain this off-diagonal term explicitly:

$$
\chi_D\partial_t\widetilde{\boldsymbol u}+\tau_u^{-1}\widetilde{\boldsymbol u}
 =\mathcal P_u(\boldsymbol R_u),\qquad
\chi_D\partial_t\widetilde T+\tau_T^{-1}\widetilde T
 =\mathcal P_T(R_T).
\tag{B7}
$$

Use that closure here. **Do not put $\boldsymbol b(T_h+\widetilde T)$ into (B6)**:
that would define a different, explicitly coupled local closure, requiring a new
projection derivation and local block solve. The resolved momentum equation
still includes the temperature-subscale force specified in §5; these are distinct
parts of the approximation. Nonlinear transport still couples the subscales
through $\boldsymbol a$.

## 4. Parameters, time discretization and OSS

Following existing repository choices, with $h$ the shortest cell edge,

$$
\begin{aligned}
\tau_u&=\left(c_1k_u^4\frac{\nu}{h^2}
                  +c_2k_u\frac{|\boldsymbol a|}{h}\right)^{-1},\\
\tau_T&=\left(c_1k_T^4\frac{\kappa}{h^2}
                  +c_2k_T\frac{|\boldsymbol a|}{h}\right)^{-1},\\
\tau_c&=\frac{h^2}{c_1k_u^4\tau_u}.
\end{aligned}
\tag{B8}
$$

Use the existing $c_1=4,c_2=2$ defaults and pressure-subscale switch. For comparison,
2010 (23), (25) use
$[(c_1\mu/h^2)^2+(c_2|\boldsymbol a|/h)^2]^{-1/2}$ with $\mu=\nu,\kappa$.
Do not silently replace (B8) by those expressions or introduce a separately
chosen turbulent Prandtl number. The thermal and mechanical subscale responses
follow from their respective transport equations. Higher-degree inverse-estimate
restrictions from the existing specification also apply to thermal diffusion.

The other closure difference is explicit: 2010 (27) gives
$\widetilde p=\tau_2\mathcal P_p(R_p+\tau_1\partial_tR_p)$.
Here the pressure subscale remains quasi-static even when velocity and
temperature subscales are dynamic; no continuity-residual history is added.

For any resolved or dynamic subscale field $y$, use

$$
D_ty=\frac{\alpha_0y-y_{\rm hist}}{\Delta t},\qquad
\tau_{u,t}=\left(\chi_D\frac{\alpha_0}{\Delta t}+\tau_u^{-1}\right)^{-1},\quad
\tau_{T,t}=\left(\chi_D\frac{\alpha_0}{\Delta t}+\tau_T^{-1}\right)^{-1}.
\tag{B9}
$$

Define $\widehat{\boldsymbol R}_u=\boldsymbol R_u+D_t\boldsymbol u_h$ and
$\widehat R_T=R_T+D_tT_h$. The algebraically eliminated subscales are

$$
\begin{aligned}
\boldsymbol X_u&=\boldsymbol f+\boldsymbol b(T_h)
 -(\boldsymbol a\cdot\nabla)\boldsymbol u_h
 +\nu\Delta_h\boldsymbol u_h-\nabla p_h-\chi_A D_t\boldsymbol u_h
 +\chi_D\widetilde{\boldsymbol u}_{\rm hist}/\Delta t,\\
X_T&=s-\boldsymbol a\cdot\nabla T_h+\kappa\Delta_hT_h
 -\chi_A D_tT_h+\chi_D\widetilde T_{\rm hist}/\Delta t,\\
\widetilde{\boldsymbol u}&=\tau_{u,t}(\boldsymbol X_u-\boldsymbol\xi_h),\\
\widetilde T&=\tau_{T,t}(X_T-\eta_h),\\
\widetilde p&=-\tau_c(\nabla\cdot\boldsymbol u_h-\zeta_h).
\end{aligned}
\tag{B10}
$$

For ASGS set all three projections to zero and allocate no projection unknowns.
For OSS solve simultaneously with the physical fields:

$$
\begin{aligned}
(\tau_{u,t}(\boldsymbol\xi_h-\boldsymbol X_u),\boldsymbol w_h)&=0
 &&\forall\boldsymbol w_h\in\widehat{\boldsymbol V}_h,\\
(\tau_{T,t}(\eta_h-X_T),r_h)&=0&&\forall r_h\in\widehat S_h,\\
(\tau_c(\zeta_h-\nabla\cdot\boldsymbol u_h),z_h)&=0
 &&\forall z_h\in\widehat Q_h.
\end{aligned}
\tag{B11}
$$

The last equation and $\zeta_h$ are absent when $\tau_c=0$. Projection spaces
carry neither Dirichlet constraints nor a pressure gauge. Consistent mass
matrices, current residuals including sources and histories, and identical
quadrature in every equation imply

$$
(\widetilde{\boldsymbol u},\boldsymbol w_h)=0,\qquad
(\widetilde T,r_h)=0.
\tag{B12}
$$

This is why resolved time derivatives can be removed from the OSS residuals,
and why $(D_t\widetilde{\boldsymbol u},\boldsymbol v_h)$ and
$(D_t\widetilde T,s_h)$ vanish with orthogonal initial histories. Weighting by
$\tau_{u,t},\tau_{T,t}$ makes (B12) hold even when these parameters vary within
an element. Lagging the projections or using independent quadrature rules
would break this argument.

Use a common positive-weight quadrature rule with degree
$q=3\max(k_u,k_T)$, including both subscale spaces and all projection terms.
Pointwise parameters make this an integration approximation, but not an
approximation of the algebraic orthogonality identity. Preserve the existing
checks for quadrature positivity and MPI sensitivity of tetrahedral rules.

## 5. Coupled resolved weak equations

Write $(\cdot,\cdot)_h=\sum_K(\cdot,\cdot)_K$ with the common quadrature rule.
For the default Laplacian viscosity and skew-symmetric momentum convection,
find $(\boldsymbol u_h,p_h,T_h)$ and any OSS projections such that

$$
\begin{aligned}
0={}&(D_t\boldsymbol u_h,\boldsymbol v_h)
 +\nu(\nabla\boldsymbol u_h,\nabla\boldsymbol v_h)
 +c_u(\boldsymbol a;\boldsymbol u_h,\boldsymbol v_h)
 -(p_h,\nabla\cdot\boldsymbol v_h)+(q_h,\nabla\cdot\boldsymbol u_h)\\
&-(\boldsymbol f+\boldsymbol b(T_h),\boldsymbol v_h)
 +\chi_A\chi_D(D_t\widetilde{\boldsymbol u},\boldsymbol v_h)_h\\
&-(\widetilde{\boldsymbol u},
     \boldsymbol a\cdot\nabla\boldsymbol v_h
       +\nu\Delta_h\boldsymbol v_h+\nabla q_h)_h
 -(\widetilde p,\nabla\cdot\boldsymbol v_h)_h
 +\beta(\boldsymbol g\widetilde T,\boldsymbol v_h)_h,
\end{aligned}
\tag{B13}
$$

$$
\begin{aligned}
0={}&(D_tT_h,s_h)+\kappa(\nabla T_h,\nabla s_h)
 +c_T(\boldsymbol a;T_h,s_h)-(s,s_h)
 +\langle j_{\rm out},s_h\rangle_{\Gamma_j}\\
&+\chi_A\chi_D(D_t\widetilde T,s_h)_h
 -(\widetilde T,\boldsymbol a\cdot\nabla s_h+\kappa\Delta_hs_h)_h.
\end{aligned}
\tag{B14}
$$

$s_h$ is the temperature test function; $s$ without a subscript is the source.
The plus sign on the outward diffusive flux follows directly from (B3).
$c_u$ is the existing flow convection form, with its existing natural boundary
terms. Define thermal convection by

$$
\begin{aligned}
c_T^{\rm nc}(\boldsymbol a;T,s_h)&=(\boldsymbol a\cdot\nabla T,s_h),\\
c_T^{\rm c}(\boldsymbol a;T,s_h)&=-(T,\boldsymbol a\cdot\nabla s_h)
 +\langle(\boldsymbol u_h\cdot\boldsymbol n)T,s_h\rangle_{\Gamma_j},\\
c_T^{\rm ss}(\boldsymbol a;T,s_h)&=\tfrac12(\boldsymbol a\cdot\nabla T,s_h)
 -\tfrac12(T,\boldsymbol a\cdot\nabla s_h)
 +\tfrac12\langle(\boldsymbol u_h\cdot\boldsymbol n)T,s_h\rangle_{\Gamma_j}.
\end{aligned}
\tag{B15}
$$

Tests vanish on $\Gamma_T$; periodic facets cancel. On physical facets the
advector is $\boldsymbol u_h$ by the bubble assumption. Start with
$c_T=c_T^{\rm ss}$, allowing the nonconservative form for comparison with
2010 (31). Conservative scalar transport is useful when exact discrete heat
balance is the priority. These choices have different discrete conservation
properties when divergence is only weakly enforced.

For nonlinear splitting, (B14) includes both the resolved-temperature cross term
due to $\widetilde{\boldsymbol u}$ inside $c_T$, and the subscale heat-flux term
$-(\widetilde T,\widetilde{\boldsymbol u}\cdot\nabla s_h)_h$.
With nonconservative convection, expansion recovers the transport terms of
2010 (31). Do not add those cross terms again outside (B14). The skew-symmetric
choice follows the alternative formulation in 2010 (36), (39).

### Buoyancy and orthogonality

The last term in (B13) is the resolved feedback of temperature subscales in
2010 (29). For OSS with constant $\beta\boldsymbol g$ and matching scalar
temperature and velocity-component spaces,

$$
\beta(\boldsymbol g\widetilde T,\boldsymbol v_h)_h
 =(\widetilde T,\beta\boldsymbol g\cdot\boldsymbol v_h)_h=0.
\tag{B16}
$$

Retain it in the general definition, or omit it only under the stated
conditions. It does not generally vanish for ASGS, variable coefficients or
incompatible spaces. Also, $\boldsymbol b(T_h)$ lies in the velocity projection
space under these same constant-coefficient/equal-space assumptions, so its
direct contribution to the OSS velocity subscale is projected out. The
resolved buoyancy force remains active. Keep this force in the residual and
projection code to preserve consistency and avoid special-case algebra.

### Implemented B1 baseline

The monolithic Galerkin residual uses the physical force (B1) in the existing
momentum kernel and adds only the resolved terms of (B14). Newton differentiates
both transport and buoyancy. Independent scalar Dirichlet and outward-flux data
use the same solve-level data weights as velocity. B1 exposes only BDF1 with
lower-order startup and the two stated convection choices; dynamic/nonlinear subscale and later time
choices below remain phase-gated. The existing solver's combined algebraic
residual criterion is supplemented by a criterion for each equation, normalized
by a recorded positive reference scale and compared with its own initial norm
and absolute tolerance. These algebraic criteria are distinct from energy norms.

### Implemented B2 extension

Static linear ASGS/OSS uses the existing momentum kernels with exactly the
resolved force $f+b(T_h)$ in their residual and projection. Scalar kernels use
(B10), (B11), (B14) with no dynamic history or dynamic feedback. The separate
$+\beta(g\widetilde T,v_h)$ resolved momentum term remains present for ASGS
and OSS; no off-diagonal force is added to the local closure. Effective static
weights equal their spatial parameters, and pointwise parameters and all
projection/adjoint products use the common positive quadrature rule.

In code, scalar projection names are `eta_temperature` and `eta_<passive_name>`;
physical scalar names occur before the projection block. Disabling the pressure
subscale removes only `zeta` and its equation. Passive scalar data live in
`Problem.passive_scalars`, with independent diffusivity/source/initial/boundary
values and inherited degree, convection and stabilization choices. All active
fields enter the monolithic Newton solve and named convergence checks. Frozen-τ
Newton freezes both scalar and mechanical parameters, including projection
weights; exact-τ Newton retains the existing zero-speed derivative limitation.

### Implemented B3 time integration

Dynamic velocity, temperature and passive subscales use their own quadrature
histories with the same time-scheme coefficients as the resolved fields. Both
residual-history terms and effective-parameter projection weights follow
(B9)–(B11). ASGS retains the scalar subscale time term in (B14); OSS removes it
through matched orthogonality. The quasi-static pressure closure is unchanged.
Sources and prescribed fluxes/Dirichlet data use solve-level endpoint weights.

Crank–Nicolson solves all fields at the endpoint-average half level and recovers
resolved velocity/scalars and dynamic subscales at the endpoint. Pressure and
projections remain at the solve level. Buoyancy and its reported power use
solve-level temperature/velocity together. Zero expansion or gravity gives
exactly zero buoyancy work. All histories are committed only after the coupled
step succeeds.

Exact-history startup requires exact velocity and every scalar, and fills their
negative-time FE histories. Subscale histories remain zero, as specified in §6.
Source-free temporal verification uses discrete diffusion eigenmodes for which
OSS subscales are exactly zero; it isolates temporal order without spatial
error or inconsistent subscale startup. A forced local-ODE check separately
measures the first-order first-step defect of zero negative BDF3 subscale
histories. Full thermal kinetic/potential and momentum budgets stay disabled
until B5 supplies their coupled terms; §7 remains their derivation/specification.

### Implemented B4 nonlinear splitting and Picard

The same $a=u_h+\widetilde u^*$ enters resolved momentum/scalar convection,
their strong residuals, subscale adjoint feedback, and every spatial/effective
parameter and OSS weight. No spatial derivative of a stored quadrature subscale
is taken. For example, nonconservative thermal transport adds
$(\widetilde u^*\cdot\nabla T_h,s_h)$, and scalar subscale feedback adds
$-(\widetilde T,\widetilde u^*\cdot\nabla s_h)_h$, each once. Skew transport uses
the corresponding skew pairing in (B15). These are the scale-interaction terms
in Codina (2010), §2.4, equations (31) and (34), printed pp. 499–500
(PDF pages 9–10); the paper uses nonconservative resolved thermal convection.
This nonlinearity comes from transport by the coupled velocity/scales, with
constant $\kappa$ and linear buoyancy, rather than a nonlinear material law.

Fixed-point and accelerated local velocity updates both include $f+b(T_h)$.
The off-diagonal $-\beta g\widetilde T$ force remains excluded from the local
velocity closure as required by (B7)/(B10); the separate coarse momentum
feedback is unchanged. Each scalar subscale is evaluated algebraically at the
updated advector. It needs no additional local nonlinear iteration.

For thermal nonlinear splitting, convergence requires every named FE/projection
change, the velocity-iterate change, each eliminated scalar-subscale change,
and $\|\widetilde u-\widetilde u^*\|_h/\|\widetilde u\|_h$ to meet `outer_rtol`.
Changes use the existing relative $L^2$ convention (absolute if the new norm is
zero). The true FE residual is reassembled after the local update: both its
combined criterion and each scaled equation's inner criterion must pass.
The scalar closure is satisfied by algebraic elimination and is independently
verified against (B10); scalar-subscale changes check its outer consistency.
CSV records expose these changes, the velocity closure residual, and the
post-update scaled FE residuals. Orthogonality is checked at convergence.

Newton differentiates all coupled unknowns at fixed $\widetilde u^*$; frozen-τ
Newton freezes all parameter/projection weights, while exact-τ Newton retains
their derivatives. Picard freezes $a=u_h^{(k)}+\widetilde u^*$ and the same
parameters in both equations, keeping temperature, buoyancy and projections
unknown in an affine monolithic system. Its convergence check uses the true
coupled residual. Picard warm-up before Newton is also supported.

### Implemented B5 diagnostics and verification

Conservative scalar convection now uses the third form of (B15), independently
of momentum convection. The global scalar default remains skew-symmetric.
`postprocessing/thermal_budgets.py` composes (B19)–(B24) from matched solve-level
fields, the same spatial parameters and the actual resolved/subscale histories.
It captures all work before history shifts or CN endpoint extrapolation.
`thermal_diagnostics.py` provides endpoint norms, gradient error, and global FE
nodal extrema; these are not extrema of an entire higher-order polynomial.

The named groups are `thermal_content_budget`, `thermal_variance_budget`,
`boussinesq_energy_budget`, `thermal_boundary_flux`, `thermal_subscale_l2`, and
`cavity_nusselt`. Scalar groups name temperature and every passive scalar.
Content reports resolved and total moments; variance includes ASGS cross-mass,
the higher-degree diffusion cross term, projection work and boundary reactions.
The kinetic–potential group separately reports every exchange remainder in
(B22), reference-temperature work, momentum convection and pressure stabilization.
At finite outer tolerance it also records the exchange correction
$\beta(t+r,g\cdot(a-U-\chi_N w))$; the convective defect uses the actual $a$.
The content identity includes the corresponding
$-\lambda_T(a-U-\chi_Nw,\nabla t)$ and continuity work. Both disappear under
their stated convergence/space assumptions; they are not silently discarded.

Residuals use signed work/rate terms from the tested discrete equations.
Normalization divides by the sum of their absolute values, with absolute
floor $10^{-12}$, recorded in thermal metadata. BDF rates are history pairings,
not finite differences of plotted energies. Legacy isothermal `energy_budget`
and `momentum_budget` keep their previous meanings; thermal energy uses the new
`boussinesq_energy_budget` name. Thermal momentum budgets remain B6 work.

Content budgets reject scalar Dirichlet walls. Potential budgets additionally
require representable, single-valued $g\cdot x$ and admissible velocity tests
(homogeneous full-vector Dirichlet, free-slip or periodic); gravity-direction
periodicity is rejected. Variance budgets
permit nonzero Dirichlet temperature by retaining variational reaction work.
Gradient wall flux and variational reaction are different diagnostics. Boundary
lifts partition shared Dirichlet nodes so corner reactions are counted once.

The [B5 verification report](verification/boussinesq_b5.md) records solved variant
identities, P1/P2 spatial and coupled temporal rates, 3D and MPI checks, cavity
refinement, and the limited purpose of the tall-cavity smoke runs. The optional
G5 proxy ratios are not implemented in B5.

### Implemented G4 channel boundaries

`periodic_channel()` creates an x-periodic affine triangle mesh with horizontal
walls tagged `y_min=3`, `y_max=4`. This uses coordinate index 1 for the vertical
direction; it is the $z$ direction in the shear-layer notation. Translation
preserves these wall tags. The fully periodic mesh helper is unchanged.

`FreeSlip(normal_component=1)` imposes only the vertical resolved velocity
component on those walls. The tangential component remains a free test/unknown;
zero tangential traction follows naturally from the viscous weak form. On a
flat wall, the tangential derivative of the zero normal trace vanishes, so
Laplacian and stress viscosity have the same tangential Neumann condition.
The geometry check rejects curved/tilted walls and incorrect axis indices.

No pressure datum or mean-horizontal-velocity constraint is added. Pressure is
still component 1 with its existing constant null space. OSS projection fields
remain unconstrained; temperature and passive insulation are natural zero flux.
Initialization clamps only the slip-normal trace, including exact-history
states; it does not project the initial velocity onto a solenoidal space.
Residual-row elimination and failed-state output use the same component BCs.

For (B24), the resolved velocity is an admissible test at a slip wall despite
nonzero tangential speed. The budget therefore checks zero normal flow on all
walls and zero full velocity only on full-vector Dirichlet walls. No additional
mechanical boundary work is required for zero natural tangential traction.
G4 verifies both viscous forms isothermally and thermal Laplacian transport;
thermal stress viscosity remains the separate B6 extension. Measured checks:
[G4 verification](verification/free_slip_g4.md).

### Galerkin and passive-temperature limits

For `method = "none"`, remove every subscale and projection contribution from
(B13)–(B14); retain resolved thermal transport and buoyancy. For $\beta=0$,
momentum must reproduce the existing isothermal method, while temperature is
passively transported. Disabling thermal physics must reproduce the old mixed
space, histories, forms and output. Setting $T_h=T_{\rm ref}$ alone does not
remove arbitrary nonzero temperature subscale histories; the reference-state
limit also requires zero thermal source and compatible initial/boundary data.

## 6. Time stepping and nonlinear solution

Use BDF1–3 or Crank–Nicolson for $\boldsymbol u_h,T_h$ and their dynamic
subscales, with coefficients supplied by the existing `time_schemes.py`.
Pressure and OSS projections are algebraic fields with no time histories.

For Crank–Nicolson, solve every spatial and coupling term at the endpoint-average
half-level, use endpoint averages of prescribed sources and boundary data, and
extrapolate both resolved fields and both dynamic subscales to the endpoint.
Pressure and projections remain at the solve level. In particular,
$\boldsymbol b(T_h)$ uses the unknown half-level temperature, not an old or
endpoint temperature. Retain the initial BDF1 step.

Initialize both subscales to zero. Lower-order BDF startup has the accuracy
limitations already documented for the flow solver. Exact-history startup
requires exact temperature as well as velocity histories, and zero subscale
histories are not exact histories of a forced subscale ODE. Temporal verification
must isolate this issue instead of promising BDF3 order for arbitrary startup.

At a fixed nonlinear-splitting iterate $\widetilde{\boldsymbol u}^{*}$, solve
the physical fields and all OSS projections **monolithically**. The temperature
unknown must enter the momentum Jacobian through buoyancy, and velocity must
enter the thermal Newton Jacobian through transport. Picard freezes the shared advector
and parameters consistently while keeping resolved buoyancy coupled to the
current temperature.

Only the velocity subscale needs an outer iterate for this diagonal closure:
it defines $\boldsymbol a$. At fixed $\boldsymbol a$, $\widetilde T$ is the
explicit expression (B10), even though its OSS projection is a global unknown.
The existing accelerated velocity update
$(I+\tau_{u,t}\nabla\boldsymbol u_h)\widetilde{\boldsymbol u}=
\tau_{u,t}(\boldsymbol X_u|_{\boldsymbol a=\boldsymbol u_h}-\boldsymbol\xi_h)$
continues to apply for the default residual with $\boldsymbol b(T_h)$ included.
After every local update, recompute the coupled FE residual, thermal subscale,
and field-change criteria. Orthogonality is required at convergence, not during
an update with temporarily fixed projections. Store histories only after the
entire coupled step converges.

## 7. Diagnostics and identities to verify

Record temperature extrema, $\|T_h\|_{L^2}$, thermal subscale norm, exact-solution
errors, thermal content, diffusive boundary fluxes, and resolved buoyancy power
$W_{b,h}=(\boldsymbol b(T_h),\boldsymbol u_h)$. For ASGS include the additional
resolved power $(-\beta\boldsymbol g\widetilde T,\boldsymbol u_h)_h$.
An OSS kinetic budget that includes subscale energy must also include
$(\boldsymbol b(T_h),\widetilde{\boldsymbol u})_h$, which vanishes under the
matching-space assumptions. Attribute each work term to the equation actually
being tested; do not reuse an unforced conservation assertion for buoyant flow.

For a smooth continuous solution, thermal content satisfies

$$
\frac{d}{dt}\int_\Omega T
=\int_\Omega s-\int_{\partial\Omega}(\boldsymbol u\cdot\boldsymbol n)T
 -\int_{\partial\Omega}j_{\rm out}.
\tag{B17}
$$

At prescribed-temperature walls, the heat flux is a reaction, not input data.
The precise discrete thermal-content balance is given in §7.1 below, including
its cancellation for nonlinear splitting. It is distinct from temperature
variance conservation and from kinetic–potential exchange.

For homogeneous thermal Dirichlet or periodic conditions, skew-symmetric
convection, dynamic OSS, and affine $P_1$ temperature, adding resolved and
subscale thermal energy equations gives the semidiscrete identity

$$
\frac12\frac{d}{dt}
 \left(\|T_h\|^2+\|\widetilde T\|_h^2\right)
+\kappa\|\nabla T_h\|^2
+(\tau_T^{-1}\widetilde T,\widetilde T)_h
=(s,T_h+\widetilde T)_h.
\tag{B18}
$$

For higher degree, the left side also contains
$-2(\widetilde T,\kappa\Delta_hT_h)_h$, requiring the diffusion/inverse-estimate
check. Boundary work must be restored for nonhomogeneous or open boundaries.
This identity follows the cancellation in 2010 (44), (46), (48), adapted to the
repository parameters and weighted orthogonality. It is a temperature-variance
balance, not physical internal energy, and should not be added to kinetic energy
without dimensional scaling. Fully discrete BDF budgets must use the actual
history terms and time-scheme dissipation; a finite difference of squared norms
alone is not the BDF energy identity. Static subscales have no subscale storage
term. ASGS needs its own cross-mass accounting.

The implemented general scalar variance budget uses

$$
\begin{aligned}
M_T={}&(D_tt,t)+\chi_A\chi_D(D_tr,t)
       +\chi_D(D_tr,r)+\chi_A(D_tt,r),\\
0={}&M_T+\kappa\|\nabla t\|^2+(\tau_T^{-1}r,r)
 -2\kappa(r,\Delta_ht)+c_T(a;t,t)+(\eta,r)\\
 &-(s,t+r)+\langle j_{\rm out},t\rangle-W_{\rm Dirichlet}.
\end{aligned}
$$

Here $\eta=0$ for ASGS. For dynamic ASGS the storage is
$\frac12\|t+r\|_h^2$; static ASGS retains the mixed rate $(D_tt,t+r)$.
$W_{\rm Dirichlet}$ is the full weak scalar residual tested with the FE
boundary trace of $t$, zero at interior nodes. It is measured, not inferred
from the temperature gradient. Orthogonality removes $(\eta,r)$ for converged
matching OSS; its measured value exposes projection error.

For cavity comparisons define the normal convention explicitly, for example
$Nu_{\rm out}=j_{\rm out}L/(\kappa\Delta T)$. At the hot wall, heat enters the
fluid, so the conventional positive hot-wall Nusselt number is
$Nu_{\rm hot}=-Nu_{\rm out}$. Distinguish a gradient-based wall estimate from a
variational reaction-flux measurement when checking discrete balance.

### 7.1 Exact thermal-content identity (review G3)

Assume no normal flow at insulated walls or periodic boundaries, matching
temperature/pressure spaces, a constant admissible scalar test, and converged
subscales. The continuity equation tested with $T_h$ gives
$(T_h,\nabla\cdot u_h)=(\widetilde u,\nabla T_h)_h$.
The pressure gauge does not remove this test: if working in a zero-mean quotient,
subtract the mean and use zero net boundary flux. Define
$J_T=(\widetilde u,\nabla T_h)_h$ and let
$\lambda_T=1,1/2,0$ for nonconservative, skew-symmetric and conservative scalar
convection respectively. Integration by parts of the **resolved** velocity gives

$$
c_T(a;T_h,1)=-\lambda_T(1-\chi_N)J_T,\qquad
D_t\int_\Omega(T_h+\chi_A\chi_D\widetilde T)
=\int_\Omega s-\int_{\Gamma_j}j_{\rm out}
 +\lambda_T(1-\chi_N)J_T.
\tag{B19}
$$

Thus all three forms preserve scalar content for nonlinear splitting with zero
source/flux; conservative convection also does so for linear splitting.
Dynamic ASGS preserves the moment of $T_h+\widetilde T$, not generally that of
$T_h$ alone. OSS has zero subscale constant moment. The pressure subscale enters
only the momentum test and does not modify this result. OSS does not generally
remove $J_T$: $\nabla T_h$ is discontinuous. The same identity applies to every
passive scalar of matching degree.

Before convergence the right side also includes the continuity residual tested
with $T_h$, multiplied by $\lambda_T$. These cancellations use exact integration
of the resolved polynomial product (satisfied by the common rule on affine
cells), and identical quadrature for subscale terms. With thermal Dirichlet data
the constant test is inadmissible; open boundaries require their explicit fluxes.

### 7.2 Potential energy and admissible geometry (review G2)

Let $\phi=\boldsymbol g\cdot\boldsymbol x$ and
$E_p=\beta(T-T_{\rm ref},\phi)$. Use this potential only when $\phi$ is a
single-valued admissible scalar test. Horizontal periodicity with vertical
gravity satisfies that condition; periodicity **in the gravity direction does
not**. A closed domain gives $\int u=0$ for divergence-free impermeable flow;
a horizontally periodic channel can have a nonzero mean horizontal velocity.
The identity needed here is $(\boldsymbol g,u)=0$, obtained with the admissible
test $\phi$, not a blanket zero-mean-velocity assertion.

For the continuous equations and impermeable walls,

$$
\frac{dE_p}{dt}=-(b(T),u)-\beta\kappa(\boldsymbol g,\nabla T)
 -\beta\langle j_{\rm out},\phi\rangle+\beta(s,\phi).
\tag{B20}
$$

Thermal source work must be retained. On Dirichlet thermal boundaries the flux
is a reaction; the simple discrete test below is not admissible there. This is
gravitational potential energy of the buoyancy anomaly per reference density,
not available potential energy obtained by sorting a stratified field.

For the discrete identities set $T_{\rm ref}=0$ initially and use the notation
$U=u_h$, $w=\widetilde u$, $t=T_h$, $r=\widetilde T$,
$a=U+\chi_Nw$. All inner products below use the common quadrature. Define

$$
\begin{aligned}
M_K={}&(D_tU,U)+\chi_A\chi_D(D_tw,U)
              +\chi_D(D_tw,w)+\chi_A(D_tU,w),\\
M_P={}&\beta(D_t(t+\chi_A\chi_Dr),\phi),\\
D_\phi={}&(\nabla\cdot U,t\phi)-\chi_N(w,\nabla(t\phi)),\\
W_b={}&-\beta\big[(t,\boldsymbol g\cdot U)
                    +(t,\boldsymbol g\cdot w)+(r,\boldsymbol g\cdot U)\big],\\
P_{\rm adv}={}&\beta(t+r,\boldsymbol g\cdot a)+\beta\lambda_TD_\phi.
\end{aligned}
\tag{B21}
$$

These follow by testing the resolved momentum equation with $(U,p_h)$, adding
the velocity closure tested with $w$, and testing the scalar equation with
$\beta\phi$. Since $\Delta_h\phi=0$, there is no scalar diffusion-adjoint
contribution. Crucially, the ASGS velocity residual includes $D_tU$ even for
**static** subscales; this produces the last term of $M_K$.

For dynamic ASGS, the kinetic energy is $\frac12\|U+w\|^2$, including the cross
term. For dynamic OSS it is $\frac12(\|U\|^2+\|w\|^2)$, with cross pairings
zero. These are already the conventions of the flow specification and existing
CN diagnostics. Static ASGS has $M_K=(D_tU,U+w)$; it cannot be replaced by the
derivative of resolved energy or by a dynamic subscale energy. The potential
storage is $E_{p,h}=\beta(t+\chi_A\chi_Dr,\phi)$.

The net exchange contribution on the **right-hand side** of the combined budget is

$$
R_{\rm exchange}=W_b+P_{\rm adv}
=\beta\big[(\chi_N-1)(t,\boldsymbol g\cdot w)
       +\chi_N(r,\boldsymbol g\cdot w)+\lambda_TD_\phi\big].
\tag{B22}
$$

Consequences, correcting hypotheses G2.1–G2.3:

- For conservative scalar transport and linear splitting, the remainder is
  $-\beta(t,\boldsymbol g\cdot w)$. It vanishes for matching-space OSS, but
  **not generally for ASGS**, whether static or dynamic. Adding ASGS potential
  subscale storage does not remove this missing exchange term.
- Changing conservative convection to nonconservative adds $-D_\phi$ to the
  convection operator on the left, or $+\beta D_\phi$ to the potential budget
  on the right; skew-symmetric adds half. For linear splitting this reduces to
  $(\nabla\cdot U,t\phi)$ in $D_\phi$. Nonlinear splitting also requires
  $-(w,\nabla(t\phi))$. The pressure test $t\phi$ generally lies outside $Q_h$,
  so continuity cannot cancel the whole defect.
- For conservative convection and nonlinear splitting, the remainder is indeed
  $+\beta(r,\boldsymbol g\cdot w)$, for ASGS and OSS. Restoring the omitted local
  off-diagonal buoyancy term would add its negative to the kinetic work, but
  would require a newly derived closure and projection. No such change is made.

For the displayed momentum forms with Laplacian viscosity, define

$$
\begin{aligned}
\mathcal D_u&=\nu\|\nabla U\|^2+(\tau_u^{-1}w,w)
 -2\nu(w,\Delta_hU)+(\tau_c(\nabla\cdot U-\zeta),\nabla\cdot U),\\
\mathcal N_u&=c_u(a;U,U)-(w,a\cdot\nabla U)+(w,a\cdot\nabla U).
\end{aligned}
\tag{B23}
$$

The deliberately uncombined terms in $\mathcal N_u$ identify resolved convection,
subscale feedback and residual convection. For other momentum forms use the
existing `energy_defect` with its matching residual convention, and use the
matching stress viscous operator when enabled. The full identity is

$$
M_K+M_P+\mathcal D_u+\mathcal N_u
 +\beta\kappa(\boldsymbol g,\nabla t)
 +\beta\langle j_{\rm out},\phi\rangle
 -(f,U+w)-\beta(s,\phi)-R_{\rm exchange}=0.
\tag{B24}
$$

This is a tested discrete equation, not an assertion of positive dissipation
for every term or of vanishing momentum-convection defects. Projection work
vanishes at converged OSS states. Record unresolved solver/projection residuals
separately when checking finite tolerances.

For nonzero $T_{\rm ref}$ the extra combined work is
$\beta T_{\rm ref}(\boldsymbol g,U+w)$, rather than only the resolved moment.
The admissible pressure test $\phi$ and impermeable walls imply
$(\boldsymbol g,U+w)=0$ at continuity convergence. If these conditions fail,
retain the term and the relevant boundary contributions. Choosing zero reference
temperature is convenient, not a replacement for this argument.

Use actual $D_t$ histories in (B19), (B21), (B24). At CN half-levels, dynamic
kinetic rates equal endpoint energy differences. BDF1 has its jump-norm
dissipation, BDF2 its two-level modified-energy identity, and BDF3 no general
positive energy-dissipation decomposition. Potential energy is linear in the
stored fields, so its corresponding history rate is direct. Static ASGS still
requires its mixed rate $M_K$. Evaluate all work before history shifting and CN
extrapolation. Keep skew-symmetric thermal convection as the global default;
study the three scalar forms for stratified problems. No choice alone removes
all closure, convection and time-discretization terms in (B24).

### 7.3 Optional stabilization-work ratios (review G5)

Define $\epsilon^u_{\rm stab}$ as the actual momentum stabilization form tested
with $(U,p_h)$, including dynamic ASGS mass coupling, pressure-gradient coupling
and pressure subscale, but excluding temperature-subscale buoyancy work.
Define $\epsilon^T_{\rm stab}$ from the scalar stabilization form tested with
$T_h$. The requested ratios are
$\nu_{\rm eff}=\epsilon^u_{\rm stab}/\|\nabla U\|^2$,
$\kappa_{\rm eff}=\epsilon^T_{\rm stab}/\|\nabla T_h\|^2$ and
$Pr_{\rm eff}=\nu_{\rm eff}/\kappa_{\rm eff}$.

Adopt these as **signed resolved stabilization-work proxies**, not positive
physical eddy coefficients or a claim of identity with 2010 (49)–(51).
They contain reversible exchange and storage; a negative value alone does not
establish physical backscatter. Pressure work and viscous cross terms can also
produce negative values in static variants. For stress viscosity, the requested
gradient denominator remains a defined proxy; an energy-equivalent viscosity
would instead use $2\|\nabla^sU\|^2$ and must have a distinct name.

Report global ratios and optional DG0 cell averages of the signed work density
and cellwise ratios, using the same original volume integrands. Never obtain the
global ratio by averaging local ratios. Near-zero denominators produce an
undefined value with a validity mask, not a clipped ratio; record the threshold
and normalization. Preserve negative work values. Use solve-level states for
these fields even when endpoint solution fields are written at the same output
event, and label the times.

## 8. Stratified Kelvin–Helmholtz application with a moisture proxy

Use a two-dimensional $(x,z)$ channel, periodic in $x$, with free-slip,
impermeable horizontal walls and insulated temperature/moisture. With scales
$L=\delta$, $U=U_0$, set
$\nu=Re^{-1}$, $\kappa_T=(Re\,Pr)^{-1}$,
$\kappa_q=(Re\,Sc)^{-1}$, $\beta=Ri$, $g=(0,-1)$ and $T_{\rm ref}=0$.
Here the argument named `expansion_coefficient` receives the dimensionless
coupling $Ri$, not a dimensional material coefficient. Then $b(\theta)=(0,Ri\theta)$.

Take $U_0(z)=\tanh z$, $\theta_0(z)=\tanh z$, and zero vertical base velocity.
The initial gradient Richardson number is $Ri\cosh^2z$. Its being below $1/4$
somewhere permits inviscid modal instability but does not prove instability.
The Miles–Howard result does not imply monotonic decay of arbitrary perturbation
energy or directly establish a result for a finite-Reynolds-number diffusing
base state; [research on stratified transient growth](https://doi.org/10.1017/jfm.2017.396)
explicitly distinguishes these mechanisms.

Use $L_z=20$ and $L_x=n_{\rm waves}2\pi/k_0$, with $k_0=0.4446$ an
**unstratified** reference. Michalke's $\sigma\simeq0.0949$ for
$U=(1+\tanh z)/2$ corresponds to about $0.1898$ for $U=\tanh z$ in these units;
see the independently computed Rayleigh comparison in
[Callaham et al., §6.2](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/multiscale-model-reduction-for-incompressible-flows/FC550C3A53CEB7174759D0C45C3DFDBB).
Do not call $k_0$ the most unstable wavenumber of the stratified, viscous,
bounded problem without a matching stability calculation.

Generate a small divergence-free continuum perturbation using
$\psi'=\epsilon A(z)[\cos(k_0x)+\epsilon_{\rm noise}
\sum_m a_m\cos(2\pi mx/L_x+\varphi_m)]$,
$u'=\partial_z\psi'$, $w'=-\partial_x\psi'$. Record modes, normalization,
amplitudes and seed; generate the same phases independently of MPI partition.
Interpolate to the FE space, impose the normal-wall condition and report both
divergence and net boundary flux before the first BDF1 step. The interpolation
need not remain exactly solenoidal. The hydrostatic base pressure satisfies
$p_0'=Ri\theta_0$ (formally $p_0=Ri\log\cosh z$ before wall adjustments).

The exact tanh profile has a small but nonzero wall derivative, and a Gaussian
perturbation has nonzero wall trace. For an exactly compatible channel case,
flatten the base profiles smoothly near the walls, and choose $A(z)$ equal to
$e^{-z^2}$ in the core with $A=A''=0$ at both walls. Record the cutoff and compare
against this actual profile; the exact $Ri\cosh^2z$ formula then applies only in
the unmodified core. Initialize moisture with zero wall-normal derivative too,
or explicitly measure the initialization adjustment. Do not maintain the base
with a body force: it diffuses on the $O(Re)$ time scale, so report its drift
during each growth fit and longer run.

Optional clustering moves triangle vertices by
$z=(L_z/2)\sinh(\lambda\zeta)/\sinh\lambda$ from uniform
$\zeta\in[-1,1]$, with the uniform limit at $\lambda=0$ and $\lambda\geq0$.
Move vertices first and retain piecewise affine cell geometry; do not represent
the sinh map as curved geometry. Check positive cell Jacobians and aspect ratios.
The common shortest-edge $h$ makes stabilization sensitive to clustered spacing;
record both minimum edge and maximum advective CFL.

For the moisture **visualization proxy**, advect/diffuse total tracer $q_m$ and set

$$
q_{\rm sat}(z,\theta)=e^{\alpha\theta-\gamma z},\quad
q_{m,0}=RH_{\rm peak}e^{-(z/\sigma_q)^2}q_{\rm sat}(z,\theta_0),\quad
q_c=\max(q_m-q_{\rm sat},0).
\tag{B25}
$$

This diagnostic removes no tracer and releases no latent heat. For $\alpha=0$
it is independent of temperature; the height dependence is an imposed saturation
proxy, not a simulated adiabatic-cooling law. Along an ideal nondiffusing parcel,
$RH(t)=RH(0)e^{\gamma\Delta z}$, so a parcel starting at the interface with
$RH=0.97$ saturates after upward displacement
$-\log(0.97)/\gamma$, approximately $1.523,1.015,0.609$ for
$\gamma=0.02,0.03,0.05$. These values motivate a sensitivity study; they do not
guarantee four cloud bands, cloud-free troughs or whole-layer saturation.

For these $\alpha=0$ initial data, the continuum tracer maximum is
$q_{m,\max}(0)=RH_{\rm peak}\exp(\gamma^2\sigma_q^2/4)<0.972$ throughout
the proposed parameter range. Source-free advection–diffusion with insulated
walls cannot increase this maximum, while $q_{\rm sat}\geq1$ at $z\leq0$.
Consequently saturation of the entire fixed channel, or of a fixed layer that
includes its lower half, is not an expected physical outcome of this proxy.
An advected moist sheet is a different region and needs an explicit definition
before discussing its saturation fraction.

Store nodal interpolation of $q_c$ for visualization, but compute cloud fraction
from the raw quadrature expression with a documented threshold and quadrature
refinement check. High-order nodal interpolation of a positive-part function
need not remain nonnegative between nodes. The transport method is not
positivity-preserving; report tracer undershoot/overshoot rather than clipping
the conserved scalar to force a cloud pattern. Two-dimensional sections can
show periodic lifted moist regions; spanwise cloud bands and pairing remain
observations to test, not prescribed outcomes.

The equations and limits above, rather than agreement with a single cavity or
cloud picture, are the acceptance basis for the implementation plan.
