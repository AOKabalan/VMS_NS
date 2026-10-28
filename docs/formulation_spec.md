# VMS formulations for incompressible Navier–Stokes: implementation spec

Oct 21, 2026 · @Ali Kabalan

Every VMS variant is one master system at the solve level with three switches: subscale space (ASGS or OSS), subscale dynamics (static or dynamic) and scale splitting (linear or nonlinear). A fourth switch selects the inner linearization (Newton or Picard); it changes the iterations, not the converged solution. All eight cases follow from these switches with no further choices; the Galerkin baseline is the same system with every stabilization term removed. The master equations below display the default convective and viscous forms; the Chapter 5 sections specify the selectable replacements, their pressure conventions, and their linearization restrictions. Revised after independent re-derivations in VERIFICATION.md; split OSS was removed because the thesis does not use it.

## Fixed choices

Every modelling decision is fixed once; "Ours" marks a choice that goes beyond what the cited source states.

| Decision | Choice | Basis |
| --- | --- | --- |
| Variants | Galerkin baseline without stabilization, and ASGS or OSS × static or dynamic × linear or nonlinear splitting: 8 VMS cases | Ours. Colomés §§2.2.1–2.2.3 define the three choices; their Table 1 benchmarks only six (no static-nonlinear cases) |
| Inner linearization | Newton or Picard, selectable; same converged solution. EMAC and EMAC–VMS: Newton only, without Picard warm-up in Phase 12a | Ours. Colomés §4 uses Picard; thesis §3.2 compares both; Chapter 5 owner decision |
| Interpolation | VMS: equal order Pk/Pk, k = 1 first, k = 2 when the viscous residual must be active. Galerkin baseline: Taylor–Hood P(k+1)/Pk | Ours. Admissible in all sources (Codina 2008 also allows unequal orders); Δ\_h vanishes for affine P1 |
| Galerkin convective form | Default skew-symmetric type 1, with ½⟨(u\_h·n)u\_h, v\_h⟩ on outflow; the Chapter 5 forms and EMAC–VMS are selectable (see Convective forms) | Volume form: Colomés eq. 10 and Remark 2.1 (needs no derivative of the discontinuous ũ; volume terms energy-neutral). Outflow term: ours |
| Viscous form | Default Laplacian ν(∇u\_h, ∇v\_h) with ν Δ\_h u\_h in the residual; stress form 2ν(∇ˢu\_h, ∇ˢv\_h) selectable, required for angular momentum and EMAC–VMS | All sources (volume terms) |
| Outflow condition | Default gradient-form do-nothing, ν ∂\_n u − p n = 0, no traction load; stress-form and modified-pressure conditions in Convective forms | Ours; see Boundary conditions and quadrature |
| Subscales on element boundaries | Zero (bubble assumption); interelement terms neglected | Codina et al. 2007 §2; VMS chapter §5.4 |
| Adjoint approximation | (ũ, (∇·a) v\_h) and element-boundary terms neglected | Colomés, text after eq. 14 |
| Pressure subscale | Quasi-static, p̃ = τ\_c P(R\_p) | Colomés eq. 16; Principe et al. eq. 11. Differs from Codina et al. 2007 eq. 12 (extra τ ∂t r\_p term) and eq. 17 (p̃ = 0) |
| τ formulas and constants | c1 = 4, c2 = 2 for k = 1; τ\_c kept (laminar) | Principe et al. 2010 eq. 12 |
| τ evaluation | Pointwise \|a\| at quadrature points, h = shortest edge of the cell, no 1/δt inside τ\_m | Pointwise: ours (VMS chapter eq. 60 and Codina 2008 use the element maximum speed). h: Colomés §7.1.2 (minimum element length); the cell diameter violates the stability condition for k ≥ 2 in Operators and parameters. No δt: Colomés §9 |
| Time scheme | BDF1, BDF2, BDF3 or Crank–Nicolson on a uniform δt, the same scheme for u\_h and ũ; BDF2 and BDF3 start with lower orders, or from exact histories when the case has an exact solution; Crank–Nicolson always starts with one BDF1 step | Same-scheme principle: Codina et al. 2007 §3.1. Crank–Nicolson: 2007 eqs. 22–23 with θ = ½. BDF2 and BDF3 for both scales and the startup: ours |
| Quadrature | One rule of degree q = 3k, k the velocity degree, with all weights positive, for every volume integral, projection and stored subscale | Ours; makes orthogonality exact in the inner product actually computed |
| OSS projection | τ-weighted (τ\_t when dynamic), onto unconstrained spaces, consistent mass; ξ\_h and ζ\_h are unknowns of the inner solve, not lagged | Weights: Colomés eqs. 33, 53. Unconstrained spaces: Codina 2008 eq. 18 and Remark 1. Coupled solve: Codina 2008 §2.3, Colomés §5.2.2; lagging rejected, see OSS projections |
| Force | Kept in the residual and the projection; requires f in L² | Colomés eq. 17 and §2.2.3 |
| Initial subscale | ũ⁰ = 0 | Ours; admissible per Castillo & Codina 2019 §3.1 |

For turbulence benchmarks Colomés et al. use c1 = 12 and τ\_c = 0, with c2 = 2 for isotropic turbulence and Taylor–Green but c2 = 8 for channel flow (§§5.1, 7.1, 8). These are not the defaults here.

## Operators and parameters

The sign convention follows Colomés eqs. 13–14 and Principe et al. eq. 9: the stabilization enters as +Σ\_K (ũ, L\*\_a(v\_h, q\_h))\_K with the adjoint carrying the minus signs.

Notation: u\_h, p\_h are the FE unknowns; ũ, p̃ the subscales; a the advection velocity; Σ\_K a sum over elements; Δ\_h the Laplacian evaluated element by element (zero for P1); V̂\_h and Q̂\_h the FE spaces without Dirichlet or zero-mean constraints. Everything is evaluated at the solve level, t^(n+1) or t^(n+½) for Crank–Nicolson, unless a superscript says otherwise. All formulas hold in d = 2 and d = 3.

```latex
\begin{aligned}
L_a(u,p) &= a\cdot\nabla u - \nu\,\Delta_h u + \nabla p \\
L_a^*(v,q) &= -\,a\cdot\nabla v - \nu\,\Delta_h v - \nabla q \\
R_u &= f - D_t u_h - L_a(u_h, p_h), \qquad R_p = -\nabla\cdot u_h
\end{aligned}
```

Stabilization parameters, evaluated pointwise at quadrature points:

```latex
\tau_m = \left(\frac{c_1\,\nu}{h^2} + \frac{c_2\,|a|}{h}\right)^{-1}, \qquad \tau_c = \frac{h^2}{c_1\,\tau_m} = \nu + \frac{c_2}{c_1}\,|a|\,h
```

With c1 = 4 and c2 = 2 for linear elements. For degree k, replace c1 by c1·k⁴ and c2 by c2·k in both parameters, so τ\_c = h²/(c1·k⁴·τ\_m) (Colomés §2.2; consistent with Castillo & Codina eq. 23). Under nonlinear splitting |a| = |u\_h + ũ|, so τ\_m and τ\_c depend on the subscale.

Stability condition for k ≥ 2. The viscous part of the adjoint term contributes −Σ\_K τ\_m ν² ‖Δ\_h u\_h‖²\_K to the energy, which the Galerkin viscous term must dominate on every element (Principe et al. 2010 §4: "provided c1 is big enough"). Since τ\_m ≤ h²/(c1·k⁴·ν), this holds when c1·k⁴ exceeds the inverse-estimate constant C\_inv of the mesh, defined by ‖Δ\_h v‖²\_K ≤ C\_inv h⁻² ‖∇v‖²\_K for all v in Pk on every cell. For P2 on right isosceles triangles, C\_inv = 48 with h the shortest edge, below c1·k⁴ = 64; with h the cell diameter, C\_inv = 96 and the condition fails, which makes the energy indefinite on fine meshes. Badly shaped cells can still violate the condition, so test U5 checks it on the meshes used. For P1, Δ\_h vanishes and no condition arises.

Time discretization, one formula for both scales on a uniform time step (coefficients for all four schemes in Time schemes):

```latex
D_t g^{n+1} = \frac{\alpha_0\, g^{n+1} - g_{\mathrm{hist}}}{\delta t}
```

The subscale ODE ∂t ũ + τ\_m⁻¹ ũ = P(R\_u) discretized with the same formula gives the effective parameter

```latex
\tau_t = \left(\chi_D\,\frac{\alpha_0}{\delta t} + \frac{1}{\tau_m}\right)^{-1}
```

where χ\_D = 1 for dynamic and 0 for static subscales, so τ\_t = τ\_m in the static case. δt enters the stabilization through τ\_t and the history term ũ\_hist/δt, both from integrating the subscale ODE; τ\_m itself contains no δt (Colomés §9 and footnote 4).

Galerkin convective form (skew-symmetric type 1). On outflow boundaries the advection velocity is u\_h because ũ vanishes on element boundaries:

```latex
c(a;\,u,v) = \tfrac12\,(a\cdot\nabla u,\, v) - \tfrac12\,(a\cdot\nabla v,\, u) + \tfrac12\,\big\langle (u_h\cdot n)\, u,\, v \big\rangle_{\Gamma_N}
```

## Time schemes

Every scheme is written as D\_t g = (α0 g − g\_hist)/δt at a solve level, so the master formulation, τ\_t and the projections need no change; Crank–Nicolson additionally solves at the half level and extrapolates.

| Scheme | Solve level | α0 | g\_hist | After convergence | Order | Basis |
| --- | --- | --- | --- | --- | --- | --- |
| BDF1 | t^(n+1) | 1 | gⁿ | Nothing | 1 | Codina et al. 2007; Colomés §4 |
| BDF2 | t^(n+1) | 3/2 | 2gⁿ − ½gⁿ⁻¹ | Nothing | 2 | Ours for both scales |
| BDF3 | t^(n+1) | 11/6 | 3gⁿ − (3/2)gⁿ⁻¹ + (1/3)gⁿ⁻² | Nothing | 3 | Ours |
| Crank–Nicolson | t^(n+½) | 2 | 2gⁿ | g^(n+1) = 2g^(n+½) − gⁿ | 2 | Codina et al. 2007 eqs. 22–23, θ = ½ |

The same row applies to u\_h and, for dynamic subscales, to ũ with its stored histories. τ\_t = (χ\_D α0/δt + 1/τ\_m)⁻¹ uses the scheme's α0, which gives τ\_t = (2/δt + 1/τ\_m)⁻¹ for Crank–Nicolson.

Crank–Nicolson. The unknowns are the half-level states u\_h^(n+½), p\_h^(n+½) and ũ^(n+½), and every spatial term, τ, a and the OSS projections are evaluated with them. This is the one-leg form of Codina et al. 2007 eqs. 22–23, rewritten with w = g^(n+½), since (g^(n+1) − gⁿ)/δt = (w − gⁿ)/(δt/2); here g^(n+½) denotes the endpoint average ½(gⁿ + g^(n+1)), not g evaluated at t^(n+½). Data at the half level are endpoint averages, as in that paper: f^(n+½) = ½(f^(n+1) + fⁿ), and on Γ\_D the half-level velocity is ½(g^(n+1) + gⁿ), so the extrapolated u\_h^(n+1) matches the Dirichlet data exactly. After convergence u\_h and ũ are extrapolated; the pressure is not, and is reported at t^(n+½). For this scheme, read every superscript n+1 in the master formulation and the algorithm as n+½. Extrapolating two orthogonal subscales keeps ũ^(n+1) orthogonal for OSS.

Startup. With lower-order startup, step 1 uses BDF1 for every scheme; BDF2 and Crank–Nicolson run from step 2, and BDF3 uses BDF2 at step 2 and runs from step 3. The subscale histories follow the same sequence. This keeps BDF2 and Crank–Nicolson second order but limits a BDF3 run to second order globally, because of the BDF1 first step. For Crank–Nicolson the BDF1 first step also matters beyond accuracy: the extrapolation u^(n+1) = 2u^(n+½) − uⁿ never damps a defect in the starting state, so an initial velocity that misses the Dirichlet data or is not discretely divergence-free would leave a sign-alternating error, with a pressure component growing like 1/δt. The backward Euler step removes that defect. Exact-history startup, available when a case provides an exact solution, fills the velocity histories at t⁰, t⁻¹ and t⁻² (negative times) from it, so BDF2 and BDF3 run at full order from step 1; Crank–Nicolson still takes its BDF1 first step. With exact-history startup, BDF2 and BDF3 use zero subscale values at negative times, which are not a history of the subscale equation, so in general problems the first subscale value has an error of order δt·|R\_u|; a BDF1 first step from ũ⁰ = 0 has no such defect. The temporal convergence test is not affected.

Stability. BDF1 and BDF2 are A-stable. BDF3 is only A(α)-stable, with α ≈ 86°: its stability region excludes a small bounded region of the left half-plane next to the imaginary axis (roughly −1/12 < Re z < 0 and |Im z| < 1.94, with z = λδt), so weakly damped, convection-dominated runs may need a smaller δt. Crank–Nicolson is A-stable but does not damp high frequencies. The stability analysis of Codina et al. 2007 covers backward Euler; a remark in that paper notes that one bound their Theorem 3 needs (their eq. 49) is not expected to hold for Crank–Nicolson, so no unconditional stability result is claimed for it here.

## Boundary conditions and quadrature

Verify every variant on periodic and homogeneous all-Dirichlet problems first; periodic cancellation applies to linear momentum and energy, not angular momentum. The table and equations here describe the default skew-symmetric type 1 and Laplacian forms; see Convective forms for the other outflow terms. The outflow condition is used only for the cylinder.

| Boundary | Velocity | Pressure | Terms added |
| --- | --- | --- | --- |
| Periodic | Periodic FE space | Periodic; mean fixed to zero | None |
| Dirichlet Γ\_D | u\_h = g strongly at the solve level (endpoint average for Crank–Nicolson) | Mean fixed to zero if Γ\_D is the whole boundary | None; ũ is not constrained (bubble assumption) |
| Axis-aligned impermeable free-slip Γ\_S | u\_h·n = 0 strongly; tangential velocity free | Mean fixed to zero when no outflow exists | Natural zero tangential traction; ũ retains the bubble assumption |
| Outflow Γ\_N | Natural | Natural | ½⟨(u\_h·n) u\_h, v\_h⟩ in c; no traction load |

For the Kelvin–Helmholtz experiment, `FreeSlip(normal_axis)` restricts the velocity test
space to zero normal trace on a flat coordinate-aligned wall, leaving tangential traction
zero naturally. This means (ν ∂\_n u)\_t = 0 for gradient viscosity and
(2ν ∇ˢu n)\_t = 0 for stress viscosity. On a flat wall with u\_n = 0 these
tangential conditions agree; the discrete volume viscous forms still differ. The wall
alignment is checked, the normal constraint is also applied to the initial velocity,
and no outflow convective correction is added on this impermeable wall.

The natural condition implied on Γ\_N by the gradient-form weak equation is the do-nothing condition:

```latex
\nu\,\partial_n u - p\,n = 0 \quad \text{on } \Gamma_N
```

The outflow term makes the skew-symmetric form behave on Γ\_N like the standard convective form (a·∇u\_h, v\_h), so the natural condition is unchanged. The volume terms of c cancel exactly, but the convective energy contribution is the outflow flux, not zero:

```latex
c(a;\, u_h, u_h) = \tfrac12 \int_{\Gamma_N} (u_h\cdot n)\, |u_h|^2
```

Where backflow occurs (u\_h·n < 0 on Γ\_N) this flux injects energy. Whether it occurs depends on the domain and Reynolds number, so monitor u\_h·n on Γ\_N in the cylinder runs; backflow stabilization would be a separate addition. The projection space Q̂\_h never carries the zero-mean constraint, even when the pressure does.

Quadrature: every volume integral (the FE equation, both projections and the energy diagnostics) uses one rule of degree q = 3k, and ũ and the subscale histories the scheme needs are stored at its points. The degree must be set explicitly because τ is not a polynomial. With one rule, the orthogonality (ũ, w\_h) = 0 holds exactly in the discrete inner product, up to solver tolerance; with mismatched rules it holds only approximately. The rule must have all weights positive so that the τ-weighted mass matrices stay positive definite; the degree alone does not guarantee this when τ varies. Exactness is limited: for linear splitting on affine Pk elements, k = 1 or 2, with τ constant within each element, q = 3k is sufficient for exact volume integration provided the force and any stored subscale histories are polynomials of degree at most 2k − 1 on each element. With pointwise τ, nonlinear splitting or a general force the integrands are not polynomial and the rule is an approximation; the orthogonality identity does not depend on exactness, only on using the same rule everywhere. In 3D the default tetrahedral rules are not invariant under permutations of a cell's vertices, and Firedrake's vertex order depends on how the mesh is partitioned. Non-polynomial integrands (τ, stored subscales, error norms) are therefore evaluated at slightly different physical points in serial and parallel runs, and results differ at the level of the quadrature error; polynomial integrands, and the 2D runs checked so far, agree to round-off. Both results are equally valid discretizations; a permutation-symmetric rule would remove the difference at up to 24 times the points per cell, which is not worth it.

## Master formulation

All eight variants solve the same FE equation and the same pointwise subscale relation at the solve level; only χ\_A, χ\_D, the choice of a and the projections ξ\_h, ζ\_h change.

FE equation, for all v\_h ∈ V\_h,0 and q\_h ∈ Q\_h:

```latex
\begin{aligned}
&(D_t u_h,\, v_h) + \nu\,(\nabla u_h,\, \nabla v_h) + c(a;\, u_h, v_h) - (p_h,\, \nabla\cdot v_h) + (q_h,\, \nabla\cdot u_h) \\
&\quad + \chi_A\,\chi_D\,(D_t \tilde u,\, v_h) + \sum_K \big(\tilde u,\; L_a^*(v_h, q_h)\big)_K + \sum_K \big(\tau_c\,(\nabla\cdot u_h - \zeta_h),\; \nabla\cdot v_h\big)_K = \langle f,\, v_h\rangle
\end{aligned}
```

Velocity subscale, pointwise at every quadrature point:

```latex
\tilde u^{n+1} = \tau_t\left(R_u + \chi_D\,\frac{\tilde u_{\mathrm{hist}}}{\delta t} - \xi_h\right)
```

Pressure subscale, already substituted in the last term on the left: p̃ = −τ\_c (∇·u\_h − ζ\_h).

| Switch | Value | Meaning |
| --- | --- | --- |
| χ\_A | 1 for ASGS, 0 for OSS | Keeps (D\_t ũ, v\_h); it vanishes for OSS because ũ is orthogonal to the FE space |
| χ\_D | 1 for dynamic, 0 for static | Activates ∂t ũ through τ\_t, ũ\_hist and the term above |
| a | u\_h (linear) or u\_h + ũ (nonlinear) | Used in c, L\_a, L\*\_a, τ\_m and τ\_c |
| ξ\_h, ζ\_h | 0 for ASGS; weighted projections for OSS (next sections) | Removes the FE-space component of the residuals |

Rules that keep the variants consistent:

- Wherever ũ appears in the FE equation it is the expression above in terms of u\_h, p\_h, not a stored value. Only ũ\_hist (dynamic) and the ũ inside a (nonlinear) are stored.
- (D\_t ũ, v\_h) = (α0 ũ − ũ\_hist, v\_h)/δt, with ũ substituted.
- For OSS, drop D\_t u\_h from R\_u both in ũ and in ξ\_h. This is exact, because the weighted projection returns any FE function unchanged, so D\_t u\_h − Π(D\_t u\_h) = 0.

Sign check: substituting ũ gives +Σ\_K (τ\_t (D\_t u\_h + L\_a(u\_h, p\_h) − f − …), a·∇v\_h + ∇q\_h + νΔ\_h v\_h)\_K, the familiar SUPG/PSPG sign. With linear splitting and static ASGS on affine P1 the system reduces exactly to SUPG + PSPG + grad-div with weight τ\_c, a useful first test; for k ≥ 2 it also keeps the viscous adjoint products with νΔ\_h v\_h. With dynamic subscales the steady state coincides with the static one, for OSS as well, once all stored states have converged (Codina et al. 2007 §3.2).

## The eight cases

Each case is the master formulation with the switch values in its row; R̂\_u = f − L\_a(u\_h, p\_h) is the OSS residual with D\_t u\_h removed.

| Case | a | Velocity subscale ũ | (D\_t ũ, v\_h) term | ξ\_h, ζ\_h | Stored per quadrature point |
| --- | --- | --- | --- | --- | --- |
| Sta-Lin-ASGS | u\_h | τ\_m R\_u | No | 0, 0 | Nothing |
| Dyn-Lin-ASGS | u\_h | τ\_t (R\_u + ũ\_hist/δt) | Yes | 0, 0 | Subscale history: ũⁿ, plus ũⁿ⁻¹ for BDF2, or ũⁿ⁻¹ and ũⁿ⁻² for BDF3 |
| Sta-Nl-ASGS | u\_h + ũ | τ\_m(a) R\_u(a) | No | 0, 0 | Current iterate of ũ |
| Dyn-Nl-ASGS | u\_h + ũ | τ\_t(a) (R\_u(a) + ũ\_hist/δt) | Yes | 0, 0 | Current iterate and the subscale history |
| Sta-Lin-OSS | u\_h | τ\_m (R̂\_u − ξ\_h) | No | Weight τ\_m, τ\_c | Nothing (ξ\_h, ζ\_h are FE functions) |
| Dyn-Lin-OSS | u\_h | τ\_t (R̂\_u + ũ\_hist/δt − ξ\_h) | No | Weight τ\_t, τ\_c | Subscale history |
| Sta-Nl-OSS | u\_h + ũ | τ\_m(a) (R̂\_u(a) − ξ\_h) | No | Weight τ\_m, τ\_c | Current iterate of ũ |
| Dyn-Nl-OSS | u\_h + ũ | τ\_t(a) (R̂\_u(a) + ũ\_hist/δt − ξ\_h) | No | Weight τ\_t, τ\_c | Current iterate and the subscale history |

The four linear-splitting cases complete a time step with the inner solve alone: their subscale is an explicit expression in u\_h and p\_h, and for OSS the projections ξ\_h, ζ\_h are unknowns of the same solve. The four nonlinear-splitting cases need the outer loop, because ũ sits inside a, τ\_m and τ\_c. Colomés et al. report the linear ASGS cases diverging on coarse turbulent meshes and dynamic OSS as robust and cheapest, with no significant cost difference between its linear and nonlinear forms (§5.2.2). Dyn-Nl-OSS is the target because it is the most complete model, not because the evidence ranks it above Dyn-Lin-OSS.

Static nonlinear subscales can be ill-posed on coarse meshes. With a = u\_h + ũ the pointwise problem is (I + τ ∇u\_h) ũ = τ R, with τ depending on |u\_h + ũ|. Where the flow is nearly at rest, τ\_m approaches its viscous limit h²/(c1·k⁴·ν); when τ\_m times an eigenvalue of −∇u\_h approaches 1, the matrix is nearly singular, the problem may have no well-behaved solution, and no local iteration converges. Dynamic subscales avoid this while δt|∇u\_h| is moderate, because τ\_t ≤ δt/α0, and finer meshes or higher degrees shrink τ\_m. This matches the sources: Colomés et al. leave static nonlinear cases out of their comparison (Table 1) and report static nonlinear ASGS diverging in their small-time-step study (Table 4). A failed local update is therefore reported as a property of the case and mesh, not repaired by relaxation or a different local solver. The same failure reaches dynamic subscales when δt is large, because the bound τ\_t ≤ δt/α0 is then inactive. At a jump in the boundary velocity, such as the lid corners of a driven cavity, |∇u\_h| ≈ 1/h in slow cells where τ\_m ≈ h²/(c1·k⁴·ν), so τ\_m|∇u\_h| ≈ Re·h/(c1·k⁴), with Re = UL/ν of the jump. Measured in Phase 10a: about 3.9 for the Re = 1000 cavity on a 64×64 P1 mesh, where Dyn-Nl-OSS at δt = 5 failed at 6 of 49,152 quadrature points, all in the two corner cells, and about 0.39 at Re = 100, which runs. Linear splitting has no pointwise nonlinear problem and is unaffected.

## Galerkin baseline

The unstabilized reference solves the master FE equation with every stabilization term removed, on an inf-sup stable Taylor–Hood pair.

```latex
(D_t u_h,\, v_h) + \nu\,(\nabla u_h,\, \nabla v_h) + c(a;\, u_h, v_h) - (p_h,\, \nabla\cdot v_h) + (q_h,\, \nabla\cdot u_h) = \langle f,\, v_h\rangle, \qquad a = u_h
```

- Spaces: Taylor–Hood P(k+1)/Pk with k ≥ 1; equal-order pairs are rejected.
- Shared with the VMS cases: the convective form and its outflow term, the time schemes and data treatment, the boundary conditions, and the Newton and Picard linearizations (Picard freezes a = û\_h).
- No subscales, projections or τ. The quadrature degree is 3(k + 1).
- Without convection stabilization it is reliable only at moderate mesh Péclet numbers. It is the reference in the verification tests, not a production method.

## OSS projections

The projections are τ-weighted L² projections onto the unconstrained FE spaces, with weight τ\_t for the velocity so that the matched subscale expression is exactly orthogonal after a converged inner solve. In nonlinear splitting this does not assert orthogonality during a local sweep with a fixed projection; see the qualification following the accelerated update.

```latex
\begin{aligned}
(\tau_t\, \xi_h,\; w_h) &= \Big(\tau_t \big(\hat R_u + \chi_D\, \tilde u_{\mathrm{hist}}/\delta t\big),\; w_h\Big) && \forall\, w_h \in \hat V_h \\
(\tau_c\, \zeta_h,\; r_h) &= (\tau_c\, \nabla\cdot u_h,\; r_h) && \forall\, r_h \in \hat Q_h
\end{aligned}
```

Why this weight. With X = R̂\_u + χ\_D ũ\_hist/δt, the subscale is ũ = τ\_t (X − ξ\_h), so for every FE function w\_h:

```latex
(\tilde u^{n+1},\, w_h) = (\tau_t X,\, w_h) - (\tau_t\, \xi_h,\, w_h) = 0
```

This holds pointwise-variable τ\_t included, which is why the history term sits inside the projected quantity (Colomés eq. 53) and why (D\_t ũ, v\_h) can be dropped for OSS. In the static case τ\_t = τ\_m and the history term is absent.

Implementation rules:

- Project onto V̂\_h and Q̂\_h, the spaces without Dirichlet or zero-mean constraints. Projecting onto the constrained spaces creates spurious numerical boundary layers (Codina 2008, Remark 1). The unconstrained space also reproduces time-dependent Dirichlet values in D\_t u\_h.
- Use the consistent mass matrix and the common quadrature rule. Lumping, or a different rule for the projection, breaks exact orthogonality.
- ξ\_h and ζ\_h are unknowns of the inner solve, together with u\_h and the pressure unknown (Codina 2008 §2.3; Colomés §5.2.2). The pressure null space applies only to p\_h, or to P\_h for EMAC and EMAC–VMS; the projections stay unconstrained, and ζ\_h with its equation is absent when τ\_c is off. Lagging them by one outer iteration, as in Colomés §4.2, has the same fixed point but was measured to fail here: with a frozen Kovasznay advector, the error contraction per outer iteration was 0.85–0.97 for P1/P1, and for static P2/P2 reached 1.19 in magnitude (divergent); half relaxation and lagging the full residual did not fix it (implementation probe and independent review, 22 October 2026). Shrinking δt does not help, because the pressure has no time-derivative term.
- This ξ\_h is not the same field as Colomés' ξ\_m, which projects the full residual: ξ\_m = ξ\_h − D\_t u\_h. The converged solutions agree; iteration counts are not comparable with the paper's lagged algorithm.
- R̂\_u keeps f and the elementwise viscous term νΔ\_h u\_h. Both cancel or vanish in common cases (f in the FE space; P1 elements), and dropping νΔ\_h does not reduce the convergence order, but keeping them is what makes the formulation exactly residual based.
- ζ\_h only matters when τ\_c is active. With τ\_c = 0 (turbulence settings) the pressure subscale, its projection and the grad-div term disappear.

## Inner linearization: Newton and Picard

For the first four convective forms, both solve the same discrete equation F(u\_h, p\_h) = 0 at fixed outer data, so when both converge to the same solution the time step is identical; they differ only in how the advection velocity and τ are treated during the inner iterations. EMAC and EMAC–VMS use Newton only, without Picard warm-up in Phase 12a.

F is the left side minus the right side of the master FE equation, with ũ substituted. Frozen during the inner solve: u\_hist, ũ\_hist, and the stored ũ\* inside a (nonlinear splitting). For OSS, ξ\_h and ζ\_h are unknowns of F together with (u\_h, p\_h), through their projection equations. The unknowns (u\_h, p\_h) enter nonlinearly in two ways: through a = u\_h (+ ũ\*) in c, L\_a and L\*\_a, and through τ\_m, τ\_c and τ\_t, which depend on |a|.

Newton: iteration j solves J δ = −F(u\_h^(j−1), p\_h^(j−1)) and updates (u\_h, p\_h) += δ, where J is the exact Gateaux derivative of F with respect to (u\_h, p\_h), including every occurrence of a and of τ. Under the usual assumptions (smooth F, nonsingular Jacobian) convergence is quadratic near the solution; it can fail from a poor initial guess at high Reynolds number, and the exact derivative of |a| is undefined at a = 0. A permitted variant freezes τ at the previous iterate, dropping its derivative. This is an inexact Newton method: it targets the same equations and avoids the non-smoothness of |a|, but does not guarantee quadratic convergence.

Picard: iteration j evaluates every advection velocity and every τ at â = u\_h^(j−1) (+ ũ\*), with τ̂ = τ(â). The unknown then appears only as the transported field, giving a linear Oseen-type system:

```latex
\begin{aligned}
&(D_t u_h^{j},\, v_h) + \nu\,(\nabla u_h^{j},\, \nabla v_h) + c(\hat a;\, u_h^{j}, v_h) - (p_h^{j},\, \nabla\cdot v_h) + (q_h,\, \nabla\cdot u_h^{j}) + \chi_A\chi_D\,(D_t \tilde u^{j},\, v_h) \\
&\quad + \sum_K \big(\tilde u^{j},\; L^*_{\hat a}(v_h, q_h)\big)_K + \sum_K \big(\hat\tau_c\,(\nabla\cdot u_h^{j} - \zeta_h),\; \nabla\cdot v_h\big)_K = \langle f,\, v_h\rangle \\
&\tilde u^{j} = \hat\tau_t \Big( f - \chi_A\, D_t u_h^{j} - L_{\hat a}(u_h^{j}, p_h^{j}) + \chi_D\, \frac{\tilde u_{\mathrm{hist}}}{\delta t} - \xi_h \Big)
\end{aligned}
```

The factor χ\_A removes D\_t u\_h for OSS (reduced residual). In c, the outflow term uses û\_h = u\_h^(j−1). For OSS the projection equations join the system with weights τ̂\_t and τ̂\_c, and are linear in (u\_h, p\_h, ξ\_h, ζ\_h). Convergence is typically linear and more robust from poor initial guesses; neither method is guaranteed to converge.

For the first four forms, recommended use is a few Picard iterations, then Newton once the residual has dropped (thesis §3.2.2). Picard with one inner iteration per outer iteration has the same structure as the Colomés §4 algorithm (one linearized solve, then the local subscale update and the projection update), but for OSS the iterates differ, because Colomés lag the projections while this spec solves them with the FE unknowns. The fixed points agree. The EMAC and EMAC–VMS Newton residuals additionally depend on the kinetic-pressure gradient and, for EMAC–VMS, every occurrence of the substituted subscale in the symmetric stress; differentiate through these dependencies. Freezing the advector alone does not give an affine Picard equation for them.

## Time-step algorithm

One time step is an inner solve for (u\_h, p\_h), Newton or Picard, inside an outer fixed point on the stored subscale; the four linear-splitting cases need only the inner solve.

Given the velocity and subscale histories the scheme needs (back to level n−2 for BDF3), and ξ\_hⁿ, ζ\_hⁿ:

1. Form u\_hist and ũ\_hist from the Time schemes table. Initialise (u\_h, p\_h) from the previous step, the stored subscale ũ\* ← ũⁿ, and ξ\_h, ζ\_h ← their previous values.
2. Inner solve: Newton or Picard on the master FE equation with ũ substituted, as in the linearization section. Frozen: ũ\* inside a (nonlinear splitting) and ũ\_hist. For OSS, ξ\_h and ζ\_h are unknowns of this solve, with their projection equations. In Picard mode, a single inner iteration per outer iteration is allowed only in cases that run step 5; the two linear-splitting ASGS cases have no outer loop, so step 2 itself must iterate until the FE residual is below tolerance.
3. Subscale update (nonlinear splitting only): at each quadrature point iterate ũ\* ← τ\_t(a) (R + χ\_D ũ\_hist/δt − ξ\_h) with a = u\_h + ũ\*, where R = R\_u for ASGS and R = R̂\_u for OSS, until the pointwise change is below tolerance. The accelerated form below is equivalent.
4. Projections (OSS only): nothing to compute, since ξ\_h and ζ\_h come out of step 2; in nonlinear-splitting OSS cases their changes still enter the convergence test of step 5.
5. Convergence: stop when the relative changes of u\_h, p\_h, ũ\*, ξ\_h and ζ\_h (when τ\_c is active) and the FE residual are all below tolerance; otherwise return to step 2.
6. Store: for dynamic subscales, ũ^(n+1) is the subscale expression evaluated at the converged state; shift the histories of u\_h, ũ, ξ\_h and ζ\_h.

Accelerated subscale update for step 3, a d×d linear solve per quadrature point, with τ\_t evaluated at the previous ũ\* and G\_ij = ∂\_j u\_h,i (UFL grad(u\_h), the transpose of the chapter gradient). Here p\_h is the physical pressure and L\_ν is the selected viscous operator:

```latex
\big(I + \tau_t G\big)\, \tilde u = \tau_t \left(f - \chi_A\, D_t u_h - u_h\cdot\nabla u_h + L_\nu u_h - \nabla p_h + \chi_D\,\frac{\tilde u_{\mathrm{hist}}}{\delta t} - \xi_h\right)
```

This is Colomés eqs. 50 and 52 with the τ\_t factor restored on the left, ξ\_h multiplied by a single τ\_t (not τ\_t² as in thesis eq. 4.39), and χ\_A dropping D\_t u\_h for OSS.

Local-update qualification (all nonlinear cases, including EMAC–VMS). The matrix can be singular for static or dynamic subscales; a time step only protects it while τ\_t‖G‖ is sufficiently small. In two dimensions, with ∇·u\_h = 0 and ω\_h = ∂₁u\_h,2 − ∂₂u\_h,1,

```latex
\det(I+\tau_t G)
=\det(I+\tau_t\nabla^s u_h)+\frac{\tau_t^2\omega_h^2}{4}.
```

Thus pure strain has exactly the same singularity as the former symmetric-matrix closure. Removing that closure does not remove the cavity-corner failure mechanism described in The eight cases. The determinant correction alone gives no general ordering of condition numbers.

For OSS, let X denote the right-hand quantity before −ξ\_h, with the subscale advection Gũ omitted. At a converged inner solve, the matched expression ũ\_expr = τ\_t(X − Gũ\* − ξ\_h) is orthogonal. If u\_h, X, τ\_t and ξ\_h are held fixed while the local system is solved, subtraction gives

```latex
(\tilde u_{\rm new},w_h)
=-(\tau_t G(\tilde u_{\rm new}-\tilde u^*),w_h).
```

This formula assumes the same τ\_t as in the matched projection. If a local sweep changes the weight from τ₀ to τ₁, the right-hand side instead is

```latex
\big((\tau_1-\tau_0)(X-G\tilde u^*-\xi_h),w_h\big)
-\big(\tau_1G(\tilde u_{\rm new}-\tilde u^*),w_h\big).
```

Neither expression is generally zero during unmatched local sweeps. Orthogonality holds for the matched expression after each converged inner solve and at convergence of the outer iteration; store that final expression as in step 6. Do not invoke OSS cancellations for the unmatched local values.

An optional alternative, outside Phase 12, is W = (τ\_t⁻¹I + G)⁻¹ and (Wξ\_h, w\_h) = (WX, w\_h), with ũ = W(X − ξ\_h). This enforces orthogonality when the same W and X are used in both equations. W is generally nonsymmetric, so the projection is oblique, not an orthogonal projection in a W norm; a positive-definite symmetric part of W is sufficient for well-posedness with the positive common quadrature and a full-rank projection space. It does not remove the requirement that the pointwise inverse exist.

Which steps each case runs:

| Case | Steps |
| --- | --- |
| Sta-Lin-ASGS | 2 |
| Dyn-Lin-ASGS | 2, 6 |
| Sta-Nl-ASGS | 2, 3, 5 |
| Dyn-Nl-ASGS | 2, 3, 5, 6 |
| Sta-Lin-OSS | 2 |
| Dyn-Lin-OSS | 2, 6 |
| Sta-Nl-OSS | 2, 3, 5 |
| Dyn-Nl-OSS | 2, 3, 5, 6 |

Initial data: ũ⁰ = 0, which also satisfies the initial orthogonality OSS needs, and ξ\_h⁰, ζ\_h⁰ from projecting the initial residuals (or zero).

## Convective forms (Chapter 5)

These sections specify Chapter 5 of the thesis (revised chapter5\_preview.pdf, 23 pages, SHA-256 c4995958c5e410591261f676aae661f5854def033feb534dc930e40668491784): five convective forms, EMAC–VMS, and the conservation checks; numbers (5.x) refer to that chapter.

The master formulation uses the fixed coupling discussed in chapter §5.7: its common subscale convective coupling is exactly the conservative form's fine-scale part (Lemma 5.3). Minimal EMAC–VMS is the separate replacement in §5.8.

Index convention. The chapter defines (∇v)ᵢⱼ = ∂ᵢvⱼ; UFL's grad uses the transpose, grad(v)\[i, j\] = ∂vᵢ/∂xⱼ. Implement every term from the table below, never by transcribing a chapter formula into UFL letter by letter; test U10 catches a swapped convention.

| Chapter | Index form | UFL |
| --- | --- | --- |
| (a·∇)u | aᵢ ∂ᵢuⱼ | dot(grad(u), a) |
| (∇u) a | ∂ᵢuⱼ aⱼ (free index i) | dot(a, grad(u)) |
| ∇ˢu | ½(∂ᵢuⱼ + ∂ⱼuᵢ) | sym(grad(u)) |
| (a ⊗ u) : ∇v = u·(a·∇)v | aᵢ uⱼ ∂ᵢvⱼ | inner(outer(u, a), grad(v)) |
| A : ∇ˢv, A symmetric | Aᵢⱼ ∂ᵢvⱼ | inner(A, sym(grad(v))) |

The forms, for advection velocity a, transported u and test function v (chapter (5.3)–(5.7)):

```latex
\begin{aligned}
c_{nc}(a;u,v) &= ((a\cdot\nabla)u,\, v) \\
c_{c}(a;u,v) &= -(u,\, (a\cdot\nabla)v) \\
c_{ss1}(a;u,v) &= \tfrac12((a\cdot\nabla)u,\, v) - \tfrac12(u,\, (a\cdot\nabla)v) \\
c_{ss2}(a;u,v) &= ((a\cdot\nabla)u,\, v) + \tfrac12((\nabla\cdot a)\, u,\, v) \\
c_{E}(a;u,v) &= (2(\nabla^s u)\, a,\, v) + ((\nabla\cdot a)\, u,\, v)
\end{aligned}
```

| Form | Config value | Outflow term on Γ\_N | Uses ∇·a | Residual operator N | Pressure unknown |
| --- | --- | --- | --- | --- | --- |
| Non-conservative | nonconservative | None | No | (a·∇)u\_h | p\_h |
| Conservative | conservative | +⟨(u\_h·n) u, v⟩ | No | (a·∇)u\_h | p\_h |
| Skew-symmetric 1 (default) | skew\_symmetric\_1 | +½⟨(u\_h·n) u, v⟩ | No | (a·∇)u\_h | p\_h |
| Skew-symmetric 2 | skew\_symmetric\_2 | None | Yes | (a·∇)u\_h | p\_h |
| EMAC | emac | None | Yes | (a·∇)u\_h + (∇u\_h)u\_h | P\_h |
| EMAC–VMS | emac\_vms | None | No | (a·∇)u\_h + (∇u\_h)u\_h | P\_h |

Outflow terms (ours): they give the first four forms the same natural condition, ν∂\_n u − p n = 0, or (2ν∇ˢu − p I)n = 0 with stress viscosity. The EMAC resolved operator is a volume form, so its natural condition comes from the viscous and pressure terms: ν∂\_n u − P n = 0, or (2ν∇ˢu − P I)n = 0 with stress viscosity. For EMAC–VMS, −(T̃, ∇ˢv) would add −T̃n to the traction if the subscale had a trace; under the zero-subscale-boundary convention T̃n = 0, and the condition is that of stress-form EMAC (chapter Remark 5.14). In physical pressure it is (2ν∇ˢu − p I)n + ½|u|²n = 0. Results on outflow problems carry this difference.

Divergence of the subscale. With a = u\_h + ũ, skew-symmetric 2 and EMAC need ∇·ũ, which stored point values cannot provide. Following chapter Remark 5.3, this spec takes as the definition

```latex
\big((\nabla\cdot a)\, u,\; v\big) := \big((\nabla\cdot u_h)\, u,\; v\big) - \big(\tilde u,\; \nabla(u\cdot v)\big)_h, \qquad a = u_h + \tilde u,
```

and more generally (∇·a, w) := (∇·u\_h, w) − (ũ, ∇w)\_h for any computable w. The continuity equation (q\_h, ∇·u\_h) − (ũ, ∇q\_h) is then exactly (∇·a, q\_h), so the chapter's Lemma 5.2 holds discretely: after a converged step, (∇·a, q\_h) = 0 for every q\_h. This replaces the earlier rule that skew-symmetric 2 and EMAC cannot be used with nonlinear splitting. The neglected term (ũ, (∇·a)v\_h) of the adjoint approximation stays neglected for every form (chapter Remark 5.2).

For linear splitting a = u\_h, use (∇·a, w) = (∇·u\_h, w) in the convective form. Its continuity equation still constrains the pairing for u\_h + ũ; it does not impose (∇·u\_h, q\_h) = 0 separately.

EMAC pressure. With EMAC and EMAC–VMS the pressure unknown is the modified pressure P\_h ≈ p − ½|u|² (chapter (5.9)); the null space and the zero-mean normalization act on P\_h. Every reported pressure (pressure errors, forces, output fields) is the physical p\_h = P\_h + ½|u\_h|², normalized separately where the null space is active. Evaluate the kinetic term and its gradient pointwise, without projecting onto Q\_h. Output the physical field in a space containing both P\_h and ½|u\_h|² (on affine polynomial cells, degree at least max(k\_p, 2k\_u)), so reporting does not discard the correction.

Under Crank–Nicolson, p\_h^(n+½) = P\_h^(n+½) + ½|ū\_h|² with ū\_h = (u\_hⁿ⁺¹ + u\_hⁿ)/2, the velocity of the converged solve. Use this pressure for half-level forces and output; pressure is not extrapolated. Averaging the two endpoint kinetic corrections would add ⅛|u\_hⁿ⁺¹ − u\_hⁿ|² and is incorrect. Manufactured volume forcing is unchanged, since the continuous equations are equivalent; prescribed traction data must respect the chosen pressure convention.

Viscous form. Two options: "laplacian" (default), ν(∇u\_h, ∇v\_h) with L\_ν v = νΔ\_h v; and "stress", 2ν(∇ˢu\_h, ∇ˢv\_h) with L\_ν v = 2ν∇·(∇ˢv) element by element. L\_ν enters L\_a and L\*\_a in place of νΔ\_h; it is formally self-adjoint, with the same element-boundary convention. Thus the residual contains +L\_ν u\_h and the adjoint coupling contains −(ũ, L\_ν v\_h). The stress form is required for angular-momentum checks (chapter §5.4: the Laplacian form adds ν∫(∇×u\_h)\_k) and for EMAC–VMS; its natural outflow condition uses p for the first four forms and P for EMAC and EMAC–VMS, as above. Forces always use the physical σ = −p I + ν(∇u + ∇uᵀ). The U5 eigenvalue check must also pass for the stress form on the test meshes.

Linearization. EMAC and EMAC–VMS use Newton only, with no Picard warm-up in Phase 12a. The standard and symmetric Picard variants of chapter §5.8.6 are deferred: freezing an advector alone does not linearize the kinetic-pressure term or the substituted stress. Any future warm-up must specify its symmetric affine equation, for example (5.42) for c\_E and a fully lagged subscale in T̃. The other forms keep Newton and Picard as before; Picard freezes a in c and in the ∇·a pairing.

## EMAC–VMS (chapter Definition 5.1)

EMAC–VMS replaces the resolved convective term and the convective part of the subscale coupling by a symmetric stress, so the subscales cannot produce a net force or torque; it requires method = "vms" and the stress form.

Resolved equation, for all v\_h, q\_h, at the solve level:

```latex
\begin{aligned}
&(D_t u_h, v_h) + \chi_A\chi_D (D_t \tilde u, v_h) + c_E(u_h;\, u_h, v_h) - (\tilde T,\, \nabla^s v_h) + 2\nu(\nabla^s u_h, \nabla^s v_h) - (P_h, \nabla\cdot v_h) + (q_h, \nabla\cdot u_h) \\
&\quad - \sum_K (\tilde u,\; L_\nu v_h + \nabla q_h)_K + \sum_K \big(\tau_c(\nabla\cdot u_h - \zeta_h),\, \nabla\cdot v_h\big)_K = \langle f, v_h\rangle
\end{aligned}
```

Compared with the master formulation, c(a; u\_h, v\_h) becomes c\_E(u\_h; u\_h, v\_h), with advection velocity u\_h, and the convective subscale coupling −(ũ, (a·∇)v\_h) becomes −(T̃, ∇ˢv\_h). The standard convective coupling must not be added as well.

Minimal subgrid stress and residual, with χ\_N = 1 for nonlinear and 0 for linear splitting (chapter (5.27), (5.30), (5.36)):

```latex
\begin{aligned}
\tilde T &= \tilde u\otimes u_h + u_h\otimes\tilde u + \chi_N\,\tilde u\otimes\tilde u, \qquad a=u_h+\chi_N\tilde u, \\
N &= (a\cdot\nabla)u_h+(\nabla u_h)u_h, \\
R_u &= f-\chi_A D_tu_h-N+L_\nu u_h-\nabla P_h.
\end{aligned}
```

There is no isotropic subgrid stress. The resolved equation uses P\_h, but the residual is exactly the master residual with the physical pressure p\_h = P\_h + ½|u\_h|², since ∇p\_h = ∇P\_h + (∇u\_h)u\_h. Here the gradient has the chapter convention; UFL implements N as `dot(grad(uh), a) + dot(uh, grad(uh))`. The χ\_A factor gives the reduced residual for OSS; it is equivalent after projection to including D\_t u\_h. No kinetic-pressure projection into Q\_h is allowed. The pressure null space, separate physical normalization, output space and midpoint level are specified in Convective forms.

The velocity subscale is the master relation ũ = τ\_t (R\_u + χ\_D ũ\_hist/δt − ξ\_h), with τ evaluated at a = u\_h + χ\_Nũ; ξ\_h, ζ\_h and the switches χ\_A, χ\_D are as for the other cases. There is no new local solve. Step 3 is the master accelerated update with G = grad(u\_h) in the UFL convention and physical p\_h as above, including its singularity and OSS-orthogonality qualifications.

Define G(v) = c\_E(u\_h; u\_h, v) − (T̃, ∇ˢv). For one common subscale field its equivalent expansion (5.32) is

```latex
G(v)=c_E(u_h;u_h,v)-\big(\tilde u,(a\cdot\nabla)v+(\nabla v)u_h\big).
```

The second term is `inner(us, dot(grad(v), a) + dot(uh, grad(v)))` in UFL. The energy cancellation is pointwise for both splittings:

```latex
\begin{aligned}
\tilde T:\nabla^s u_h
&=2\tilde u\cdot(\nabla^s u_h)u_h
  +\chi_N\,\tilde u\cdot(\nabla^s u_h)\tilde u\\
&=\tilde u\cdot\big[(a\cdot\nabla)u_h+(\nabla u_h)u_h\big]
=\tilde u\cdot N.
\end{aligned}
```

It uses 2(∇ˢu\_h)u\_h = (u\_h·∇)u\_h + (∇u\_h)u\_h and ũ·(∇ˢu\_h)ũ = ũ·(ũ·∇)u\_h. With c\_E(u\_h; u\_h, u\_h) = 0 under the energy-flux conditions below, G(u\_h) + (ũ, N) = 0, giving (5.35). Also ∇ˢe\_k = ∇ˢφ = 0 for φ = e\_k × x, so the subgrid stress gives zero force and torque for arbitrary point values of ũ, without a bubble assumption. The resolved EMAC contribution is boundary fluxes only; with u\_h zero on the boundary, G(e\_k) = G(φ) = 0. These algebraic statements do not by themselves supply admissible tests for the solved momentum equations.

Assembly and iteration (owner decision). Build T̃ directly as `outer(us_expr, uh) + outer(uh, us_expr) + chiN*outer(us_expr, us_expr)`, with the substituted subscale expression in every slot, and differentiate through it for Newton. In the residual and τ only, keep the master's stored advector a\* = u\_h + χ\_Nũ\*. The resolved c\_E uses u\_h alone. Do not assemble the expanded coupling with a\*: it mixes ũ\_expr and ũ\* and need not give a symmetric stress during iteration. The expansion is an algebraic test with one common ũ.

T̃ and its Newton derivative are symmetric at every iterate. The substituted stress is nonlinear in u\_h; it would be affine if the subscale were fully lagged (chapter Remark 5.11). Symmetry removes the convective force/torque contribution, up to resolved EMAC fluxes, at every Newton iterate and in its linearization (Proposition 5.12). Full momentum balances, including the ASGS subscale time term, and the energy balance require convergence of the coupled iteration. In particular, the energy cancellation uses the same ũ in T̃ and a, as obtained at outer convergence. The optional matrix-weighted projection described in the master algorithm is outside Phase 12.

Linear splitting (χ\_N = 0): T̃ = ũ⊗u\_h + u\_h⊗ũ and a = u\_h in the residual. The same pointwise energy identity holds without the quadratic term. The subscale is an expression throughout, and no outer loop is needed.

Consistency (Proposition 5.8). At ũ = 0 this is exactly resolved EMAC. Under the chapter's zero-trace model, comparison with the standard nonlinear convective coupling gives a difference ((∇·U)u\_h, v) + (½∇|u\_h|², v), up to boundary fluxes, with U = u\_h + ũ. Pointwise ∇·U = 0 removes the first term and a continuous pressure change absorbs the second. This compares convective forms, not discrete schemes: the continuity equation imposes only orthogonality of ∇·U to Q\_h, and ½|u\_h|² generally does not belong to Q\_h.

## Conservation checks

Distinguish algebraic defect identities, evaluable for any matching fields, from balances of the solved equations. The latter require admissible velocity tests, the stated boundary conditions, and convergence of the coupled iteration. Pointwise tensor identities and matched weighted projections are exact in the common quadrature inner product. Identities that also use integration by parts require exact integration of the remaining FE polynomial products: on affine Pk meshes degree 3k suffices here, with sufficiently exact facet quadrature. Common quadrature alone is not an integration-by-parts theorem on curved or non-polynomial geometry.

Defect functionals, assembled from the forms actually used by the code:

- Linear momentum: D\_k = c(a; u\_h, e\_k); for EMAC–VMS, G(e\_k).
- Angular momentum: D\_φ = c(a; u\_h, φ) − (ũ, (a·∇)φ), with φ = e\_k × x (2D: φ = (−y, x)); for EMAC–VMS, G(φ).
- Total convective energy defect: 𝒩 = c(a; u\_h, u\_h) + (ũ, N − (a·∇)u\_h) for the first five forms; for EMAC–VMS, 𝒩 = G(u\_h) + (ũ, N). Include the selected outflow term in c.

Their interior predicted values are chapter Tables 5.2 and 5.3 (zero in every column for minimal EMAC–VMS), with the divergence pairing for the selected splitting. In particular, writing D\_χ(w) = (∇·u\_h, w) − χ\_N(ũ, ∇w), the energy column is:

| Form | Nonlinear splitting | Linear splitting |
| --- | --- | --- |
| Non-conservative | −½D\_1(u\_h·u\_h) | −½D\_0(u\_h·u\_h) |
| Conservative | +½D\_1(u\_h·u\_h) | +½D\_0(u\_h·u\_h) |
| Skew-symmetric 1 or 2 | 0 | 0 |
| EMAC | (ũ, (u\_h·∇)u\_h) | ½(ũ, ∇(u\_h·u\_h)) |
| EMAC–VMS | 0 | 0 |

Here the argument of D\_χ is the scalar u\_h·u\_h. The resolved-scale EMAC production (ũ, ω\_h × u\_h) is a different quantity (Proposition 5.3), not 𝒩.

Predicted momentum defects. Report these alongside the measured `linear_momentum_defect` and `angular_momentum_defect`, using the same solve-level fields and common quadrature. Write w = u\_h·φ and J = (ũ, (u\_h·∇)φ); in 2D φ = (−y, x) and J = ∫(u\_h,1 ũ₂ − u\_h,2 ũ₁). The chapter's Tables 5.2 and 5.3 give:

| Form | D\_k, nonlinear | D\_φ, nonlinear | D\_k, linear | D\_φ, linear |
| --- | --- | --- | --- | --- |
| Non-conservative | −D\_1(u\_h,k) | −D\_1(w) | −D\_0(u\_h,k) | −D\_0(w) − J |
| Conservative | 0 | 0 | 0 | −J |
| Skew-symmetric 1 or 2 | −½D\_1(u\_h,k) | −½D\_1(w) | −½D\_0(u\_h,k) | −½D\_0(w) − J |
| EMAC | (ũ, ∂\_k u\_h) | (ũ, (φ·∇)u\_h) | 0 | −J |
| EMAC–VMS | 0 | 0 | 0 | 0 |

The boundary qualifications below also apply to these predictions. Use matching subscales after convergence; during nonlinear iteration the stored advector and substituted subscale can differ.

Boundary scope. With u\_h = 0 on all of ∂Ω, all convective fluxes vanish, including for arbitrary quadrature-point ũ. On periodic domains the linear-momentum and energy fluxes cancel, but the angular-momentum fluxes generally do not: φ is not periodic. Use the explicit boundary corrections in the Phase 11a verification (H-D1) otherwise; the corrections are unchanged by the minimal stress. Run angular-defect checks only with zero boundary u\_h or compactly supported flows. For full momentum balances, use a strip of elements in which the solution and subscales vanish, as in the compact Gresho setting, so admissible cutoff tests reproduce e\_k and φ on the active region; account for imposed forces, wall reactions or fluxes if they remain. Homogeneous velocity data alone do not make rigid tests admissible.

OSS moments. The identities ∫ũ = 0 and ∫x × ũ = 0 apply on the whole domain only when the corresponding rigid test belongs to the unconstrained projection space, not on arbitrary patches. In particular φ is not in a periodic projection space. These cancellations require matched converged states; do not use them during unmatched local sweeps.

Variable τ. Chapter assumption (A2), constant τ, is unnecessary for this weighted model (Remark 5.6). Matched τ-weighted velocity projection gives (ũ, ξ\_h) = 0, and the pressure projection gives (τ\_c(d − ζ\_h), ζ\_h) = 0, where d = ∇·u\_h. Therefore (τ\_c(d − ζ\_h), d) = (τ\_c(d − ζ\_h), d − ζ\_h) ≥ 0 with positive quadrature. ASGS uses ξ\_h = ζ\_h = 0 and requires no projection or constant-τ assumption. This variable-weight formulation is not obtained by multiplying an unchanged unweighted OSS projector by a variable τ.

Crank–Nicolson closure. For dynamic subscales, from step 2 onward, use the one-leg midpoint scheme for both scales, fixed FE spaces and quadrature, and the energy-flux/test conditions above. At a converged step let ū, s̄, P̄ (or p̄), ζ̄ and τ̄ be the solve-level fields; let uⁿ, sⁿ and uⁿ⁺¹, sⁿ⁺¹ be the endpoint states. The exact discrete counterpart of (5.20), or (5.35) for EMAC–VMS, is

```latex
\begin{aligned}
E^j_{\rm ASGS}&=\tfrac12\|u^j+s^j\|_h^2,\qquad
E^j_{\rm OSS}=\tfrac12(\|u^j\|_h^2+\|s^j\|_h^2),\\
\mathrm{rate}&=(E^{n+1}-E^n)/\delta t,\qquad
\mathrm{power}=(\bar f,\bar u+\bar s)_h,\\
\mathcal D&=a_\nu(\bar u,\bar u)
  +(\bar\tau_m^{-1}\bar s,\bar s)_h
  -2(\bar s,L_\nu\bar u)_h
  +(\bar\tau_c(\nabla\cdot\bar u-\bar\zeta),\nabla\cdot\bar u)_h,\\
r&=\mathrm{rate}+\mathcal D+\mathcal N_{\rm pred}-\mathrm{power}=0,\\
\mathcal N_{\rm meas}&=\mathrm{power}-\mathrm{rate}-\mathcal D
=\mathcal N_{\rm pred}.
\end{aligned}
```

Here a\_ν(u,u) is ν‖∇u‖² or 2ν‖∇ˢu‖². Every term of 𝒟 and 𝒩 is evaluated at the midpoint, and f̄ is the prescribed endpoint-average force. The residual r must vanish for every convective form, including those with nonzero 𝒩. Equivalently Eⁿ⁺¹ − Eⁿ + δt(𝒟 − power) = −δt𝒩. Do not replace τ\_m by τ\_t in the dissipation. The viscous cross term has coefficient −2. The weighted pressure term above is the implemented meaning of the chapter's (τ\_c P̃(d), d).

The BDF1 startup is excluded from this midpoint check: it adds ‖Δ(u+s)‖²/(2δt) for ASGS or (‖Δu‖² + ‖Δs‖²)/(2δt) for OSS. BDF2 has a two-level modified-energy identity; no general analogous dissipation claim is made for BDF3 (Remark 5.12). A per-step total linear-momentum balance uses the same integrator on u\_h and ũ, or the whole-domain OSS moment cancellation when admissible. With mixed integrators and nonzero subscale moments the tested discrete quantity mixes the two time operators.

Diagnostics must retain the existing names and meanings; add these defect and budget quantities alongside them. Capture the old endpoint before shifting history and keep solve-level and endpoint state explicit, particularly for midpoint pressure, dissipation and power. Algebraic defects may be computed at arbitrary matching fields; the per-step closure is a solved-equation check.

Momentum budget (`momentum_budget`). As for `energy_budget`, accept it only for Crank–Nicolson with dynamic VMS subscales and leave the BDF1 startup columns blank. From step 2, for each component k, define

```latex
\begin{aligned}
M_k^j &= (u_h^j+\chi_A\tilde u^j,e_k)_h,\\
\mathrm{rate}_{M,k} &= (M_k^{n+1}-M_k^n)/\delta t,\\
D_k &= c(\bar a;\bar u_h,e_k)
       \quad\text{or }G(e_k)\text{ for EMAC--VMS},\\
F_k &= (\bar f,e_k)_h,\qquad
r_{M,k}=\mathrm{rate}_{M,k}+D_k-F_k.
\end{aligned}
```

M uses endpoint fields; D and F use the converged solve-level fields, with the same endpoint-average forcing as the equation. D is the actual assembled defect, including the selected outflow term, rather than its interior prediction. The ASGS subscale moment must be included; OSS uses only the resolved moment, as in its discrete equation. Retain the existing `linear_momentum` diagnostic's resolved-velocity meaning.

Since derivatives of e\_k vanish, the viscous, pressure, pressure-subscale and viscous-subscale terms contribute zero. On periodic domains e\_k is admissible, so r\_M = 0 to solve tolerance. With Dirichlet walls it is not admissible: r\_M is the consistent boundary reaction, defined by evaluating the unmodified velocity residual on the constant extension of the prescribed boundary test. Equivalently it is the sum of the residual reactions at constrained velocity degrees of freedom, up to the remaining free-equation residual. Its sign is the force on the fluid. It is reported, not asserted to vanish. It need not equal the negative of the existing `force` diagnostic, which evaluates the physical stress traction directly and reports fluid-on-body force. Homogeneous boundary data alone do not remove reactions.

Monitor-region diagnostics (`region_velocity_l2`, `region_subscale_l2`). An optional case-supplied region in `Problem` defines R\_monitor; requesting either norm without that region is invalid. Report √∫\_Rmonitor|u\_h|² and √∫\_Rmonitor|ũ|², using the common quadrature and the same velocity/subscale time levels as the existing kinetic-energy diagnostics. The subscale norm requires VMS. These are observations of the zero-strip assumption, not constraints on the fields or conservation assertions.

Required Chapter 5 tests (implemented in Phase 12a):

| Test | Requirement |
| --- | --- |
| T1 | Stress and expansion (5.32) give identical G(v) for random u\_h and one common arbitrary quadrature-point ũ, for both splittings. |
| T2 | With zero boundary u\_h, G(e\_k), G(φ) and G(u\_h) + (ũ,N) vanish, including nonzero ũ at boundary-adjacent quadrature points. No bubble assumption is needed for these algebraic tests. |
| T3 | T̃ : ∇ˢu\_h = ũ·N pointwise, with χ\_N on the quadratic contraction. |
| T4 | Setting ũ = 0 gives G(v) = c\_E(u\_h;u\_h,v). |
| T5 | The midpoint closure residual r vanishes for every form on unit-scale dynamic Crank–Nicolson steps from step 2 onward; measured and predicted defects agree. |
| T6 | Reproduce the periodic angular-momentum counterexample in H-D1; periodicity must not suppress this flux. |

Relation to the thesis chapter. The index-convention table and divergence pairing specify the discrete interpretation of the chapter. For degree k the bound is τ\_m ≤ h²/(c1·k⁴·ν); the chapter's c1 > 2C²\_inv (Laplacian) and c1 > 4C²\_inv (stress) are sufficient conditions with c1 replaced by c1·k⁴, while U5 checks positivity directly on the test meshes. The master's local solvability and unmatched-sweep orthogonality limitations apply also to minimal EMAC–VMS.

## Gresho vortex (Phase 12b)

Domain Ω = (−0.5, 0.5)², built-in triangular mesh, no-slip on all four walls, zero forcing, and the constant-pressure null space active. All four walls are force boundaries. Set the initial velocity to the interpolated exact Euler velocity below; the reference fields are steady for every t. With r = √(x²+y²),

```latex
u_\theta(r)=
\begin{cases}
5r,& r<0.2,\\
2-5r,&0.2\le r<0.4,\\
0,&r\ge0.4,
\end{cases}
\qquad
u(x,y)=
\begin{cases}
(-5y,5x),&r<0.2,\\
u_\theta(r)(-y/r,x/r),&r\ge0.2,
\end{cases}
```

```latex
p(r)=
\begin{cases}
12.5r^2+C_1,&r<0.2,\\
12.5r^2-20r+4\ln r+C_2,&0.2\le r<0.4,\\
0,&r\ge0.4,
\end{cases}
\quad
C_2=-12.5(0.4)^2+20(0.4)-4\ln0.4,\quad
C_1=C_2-20(0.2)+4\ln0.2.
```

The vorticity ∂₁u₂ − ∂₂u₁ is 10, 2/r − 10, and 0 in the three regions. The pressure is continuous at both interfaces, with p′ = u\_θ²/r; hence ∇·u = 0 and (u·∇)u + ∇p = 0 region by region. Evaluate the central branch without dividing by r. Pressure errors use the physical pressure and its existing separate mean normalization.

Viscosity is the Gresho case parameter `case_parameters.viscosity`, default ν = 10⁻⁶. This is the authorized fallback for the nominally inviscid OSS families: at ν = 0 the frozen-τ OSS Jacobian has identically zero ζ rows in the resting exterior and the first linear solve fails. The viscous ASGS families use ν = 10⁻³. At either positive viscosity, with f = 0 unchanged, the Euler fields above are reference fields, not exact solutions of the unforced viscous problem.

Let β = τ\_m⁻¹ = c1·k⁴ν/h² + c2·k|a|/h. In the resting exterior, ν = 0 and a = 0 give β = 0, τ\_m = ∞, τ\_t = (α0/δt + β)⁻¹ = δt/α0 for dynamic subscales, and τ\_c = h²β/(c1·k⁴) = 0. The implemented residual and diagnostics use finite reciprocal/dynamic combinations. τ\_c is generally nonzero inside the moving vortex even at ν = 0.

**Dynamic ASGS limitation (owner decision Q1).** Where a = 0 and ν ≈ 0, τ\_m⁻¹ approaches zero (exactly zero in the inviscid limit). The dynamic ASGS subscale absorbs ∂ₜu_h, so its time-derivative coupling cancels the resolved mass. After elimination the effective resolved mass factor is (α0/δt)(1 − α0τ\_t/δt), which vanishes when τ\_m⁻¹ = 0. This is an accepted limitation of the method, not a solver problem; raising iteration limits does not remove the cancellation. The inviscid Gresho families therefore use OSS only: `gresho_{form}_oss_{split}`, `gresho_long_{form}` and `gresho_refined_emac_vms_oss`. ASGS is exercised separately by `gresho_visc_{form}_asgs_{split}` and `gresho_visc_refined_emac_vms_asgs`, at ν = 10⁻³ with the same mesh/time/diagnostic settings as their corresponding short OSS runs.

Acceptance (owner decision K-D1, updated by Q1):

- With frozen τ (`newton_tau = "frozen"`), a complete step on a tiny mesh must be finite: residual, Jacobian, solve and every diagnostic, for both OSS splittings at ν = 10⁻⁶ and both ASGS splittings at ν = 10⁻³.
- Keep the strict expected convergence failure for nonlinear dynamic ASGS at ν ≤ 10⁻⁴ as documentation of the method limitation; the existing test exercises the ν = 10⁻⁶ fallback. Linear ASGS at the fallback remains a regression check, not a matrix family.
- With exact τ the Jacobian at rest is non-finite, because the derivative of |a| is undefined at a = 0 for every ν (Inner linearization). The unit tests keep an explicit expected-failure check for it.
- Exact-τ Jacobians are tested only at states with nonzero advection speed, as in Phase 12a.
- No regularization of |a| and no generalized derivative.
- All Phase 12b configurations use `inner_max_iterations = 200` (owner decision Q2); the general solver default is unchanged.

The monitor region is the near-wall frame max(|x|, |y|) ≥ 0.45. Report both region norms and all four wall forces. The exact initial vortex vanishes outside r = 0.4, but the discrete velocity and subscales are not forced to remain zero in the frame. Momentum and angular-momentum drifts must be interpreted with those monitors and reactions, rather than asserted to vanish from the initial support alone.

## Conservation-run normalization (Phase 12b)

Owner decision K-D2. The Phase 12b analysis normalizes by fixed characteristic scales of each run, never by a value that can be zero, such as an initial momentum or a predicted defect. These scales apply to the analysis script only; the Phase 12a unit tests keep their own criteria.

With E(0) = ½‖u\_h(0)‖², the kinetic energy of the interpolated initial velocity (ũ⁰ = 0, so it is also the initial ASGS and OSS budget energy), and d the dimension:

```latex
U_c=\sqrt{2E(0)/|\Omega|},\qquad L_c=|\Omega|^{1/d}.
```

| Quantity | Scale |
| --- | --- |
| Linear momentum | M\_c = \|Ω\|U\_c |
| Angular momentum | A\_c = \|Ω\|U\_c L\_c |
| Energy | E(0) |
| Force | F\_c = \|Ω\|U\_c²/L\_c |
| Torque | F\_c L\_c |
| Power | P\_c = \|Ω\|U\_c³/L\_c |

Each closure residual, defect and drift is reported both in absolute terms and divided by its scale:

- closures: r/P\_c and r\_M,k/F\_c;
- measured defects: D\_k/F\_c, D\_φ/(F\_c L\_c), 𝒩/P\_c;
- drifts: (M(t) − M(0))/M\_c, (A(t) − A(0))/A\_c, (E(t) − E(0))/E(0).

A measured defect is compared with its prediction by |measured − predicted| divided by the same scale, reported next to the predicted value in the same units. Never divide by the prediction.

The Phase 12b expectations read with these scales:

- |r| ≤ 10⁻⁹ P\_c for every converged run from step 2;
- |r\_M,k| ≤ 10⁻⁹ F\_c on Taylor–Green;
- EMAC–VMS measured D\_k, D\_φ and 𝒩 at most 10⁻¹² of their scales;
- the other forms: |measured − predicted| ≤ 10⁻⁹ of their scales.

## Chapter 5 references

Active references use the revised 23-page PDF identified above. The previous column refers to the 19-page Phase 11b PDF (SHA-256 720be75921b3956325668b85cac1641428b0fb626e9e4ca4bf67e5d248b014bc). The full-stress variant is retained in the chapter for comparison and is not implemented.

| Item | Previous chapter | Revised chapter |
| --- | --- | --- |
| Fixed coupling versus consistent splitting | Remark 5.2 | Remark 5.2, §5.7, Lemma 5.3 |
| Consistent-splitting natural residual and defects | — | (5.25)–(5.26), Propositions 5.6–5.7, Table 5.1 |
| Fixed-coupling defect tables / properties | Tables 5.1–5.3 | Tables 5.2–5.4 |
| EMAC symmetric stress and identity | S(w): (5.21), Lemma 5.3, (5.22) | S(w): (5.22), Lemma 5.4, (5.23); new identity (5.24) |
| EMAC–VMS definition | §5.7, Definition 5.1 | §5.8, Definition 5.1, minimal closure |
| Minimal subgrid stress | (5.23) | (5.27), no isotropic part |
| Resolved equation and G | (5.24)–(5.25) | (5.28)–(5.29) |
| Physical-pressure residual | (5.26) | (5.30), N = (a·∇)u\_h + (∇u\_h)u\_h |
| G as stress and expanded coupling | (5.27)–(5.28) | (5.31)–(5.32) |
| Consistency of forms | Proposition 5.6, (5.29)–(5.30) | Proposition 5.8, (5.33)–(5.34) |
| EMAC–VMS energy balance | (5.31) | (5.35) |
| Linear splitting | (5.32) | (5.36) |
| Momenta, energy, linear splitting, iteration | Propositions 5.7–5.10 | Propositions 5.9–5.12 |
| Momentum boundary scope | §5.7.3 | §5.8.2 |
| Full-stress M\_h variant | §5.7.6, (5.33)–(5.35) | §5.8.5, (5.37)–(5.39), not implemented |
| Interpretation / eliminated subscale | Remarks 5.8–5.9 | Remarks 5.10–5.11 |
| Linearization, DS and symmetric Picard | §5.7.7, (5.36)–(5.38) | §5.8.6, (5.40)–(5.42) |
| Time integration / stabilization | Remarks 5.10–5.11 | Remarks 5.12–5.13 |
| Boundary conditions and pressure | Remark 5.12 | Remark 5.14 |
| Summary / implementation | §§5.8–5.9 | §§5.9–5.10 |
| Implementation Picard / Newton | (5.39)–(5.40) | (5.43)–(5.44) |
| Implementation stress matrix | §5.9.3, (5.41) | §5.10.3, (5.45) |
| Forms, continuity and total energy | (5.1)–(5.20), Lemma 5.2, Corollary 5.1 | Unchanged |
| Projection assumption / (A2) | Remark 5.6 | Unchanged |

## Errata

Eight inconsistencies were found by checking printed equations against their own semi-discrete sources, and an independent re-derivation confirmed the corrections; the formulation above uses the corrected forms.

| Where | As printed | Correct form | Why it matters |
| --- | --- | --- | --- |
| Colomés eqs. 48, 51 (fully discrete) and eq. 38 (semi-discrete energy balance); thesis eqs. 4.36, 4.38 | Adjoint products with a plus sign: +(τ(…), L\*\_a(v\_h, q\_h)) on both sides | Minus sign on both sides, i.e. +(τ(…), a·∇v\_h + ∇q\_h + νΔ\_h v\_h) | Follows from substituting ũ into eq. 13 with L\* from eq. 14. As printed, convection and pressure stabilization would be anti-diffusive. Eq. 39 has the correct sign for its non-force terms but omits the force-transfer term, which is harmless only when that term vanishes |
| Colomés eqs. 50, 52; thesis eq. 4.39 | ũ + ũ·∇u\_h = … | ũ + τ\_t ũ·∇u\_h = … | Dimensionally inconsistent as printed |
| Thesis eq. 4.39 | −τ\_t ξ\_h placed inside a bracket already multiplied by τ\_t, i.e. −τ\_t² ξ\_h | A single −τ\_t ξ\_h, as in Colomés eq. 52 | Thesis only; independent of the previous row |
| Colomés eq. 51; thesis eq. 4.38 | Current u\_h/δt kept in the stabilization on the left, u\_hⁿ⁻¹/δt missing on the right, while eq. 53 projects the full residual | Either restore the history term and keep the full-residual projection, or (this spec) drop D\_t u\_h from R\_u and project the reduced residual | Deleting only the current-time term while keeping the full-residual projection is wrong |
| VMS chapter eqs. 63, 64 (semi-discrete) | Eq. 63 omits (∂t ũ, v\_h); eq. 64 omits ∂t u\_h from the subscale residual, although eqs. 57–58 contain both | Keep both terms for dynamic ASGS | Both vanish for exact OSS, which may be the intended reading, but as printed the pair cannot be used for ASGS |
| Thesis eq. 4.9 | −(∂t u\_h + a·∇u\_h, ṽ) + (−ν∇²u\_h + ∇p\_h, ṽ) | Second bracket with a minus sign | Convection sign is already correct; viscous and pressure signs are wrong |
| Thesis eqs. 4.11, 4.15 | f − ∂t u\_h + a·∇u\_h − ν∇²u\_h + ∇p\_h | f − ∂t u\_h − L\_a(u\_h, p\_h) | All spatial signs reversed; eq. 4.30 later uses the correct sign without reconciling the two |
| Thesis eqs. 4.28–4.29 and 4.36–4.39 | Backward Euler for ũ with BDF2 declared for u\_h; the displayed fully discrete equations use backward Euler for u\_h too | The same BDF for both scales (our requirement) | Mixed schemes are admissible in principle (Codina et al. 2007 §3.1) but leave a first-order defect; the displayed equations do not implement the declared BDF2 at all, so the scheme actually used is unknown |

The first row has the largest effect on any code that implemented the printed equations: it flips the stabilization in the main fully discrete equations of both the paper and the thesis, while the semi-discrete equations they were derived from are correct. The τ\_t, residual-sign and time-history rows also matter for code that implemented those equations literally.

## Further thesis Chapter 4 errors

These derivation errors do not change the formulation above; they need fixing when Chapter 4 is rewritten. Their effect on the original thesis code is unknown.

| Eq. | Problem | Correct form |
| --- | --- | --- |
| 4.1, 4.2 | Wrong test functions: q\_h in the continuous form; v instead of v\_h in the subscale term | q in 4.1; v\_h in 4.2 |
| 4.3 | +(a·∇v\_h, ũ); viscous boundary term written with ν∇²ũ·n; continuity boundary term with a minus sign | −(ũ, a·∇v\_h) − (ũ, (∇·a) v\_h), the second neglected by assumption; ν⟨ũ, ∂\_n v\_h⟩ on ∂K; +⟨q\_h, ũ·n⟩ on ∂K |
| 4.8 | +p n in the viscous traction; an extra convective boundary term; ũ used as the test function | Boundary term ⟨ν ∂\_n(u\_h + ũ) − (p\_h + p̃) n, ṽ⟩ on ∂K only |
| 4.10 | (∇·ũ, q\_s) = −(∇·ũ, q̃) = 0 | (∇·ũ, q̃) = −(∇·u\_h, q̃); neither side is zero in general |
| 4.11–4.18 | Presented as exact: ξ\_h = Π(R\_m) | Exact only if the projection of ∂t ũ + a·∇ũ − νΔũ + ∇p̃ onto the FE space vanishes; state it as a closure assumption |
| Text after 4.34 | Nonlinear subscales called explicit | Implicit in u\_h and p\_h; solved iteratively |
| 4.40–4.41 | Projection equations quantified over V0, Q0; symbols change (ξ\_h to ξ\_m, R^m to R\_u) without explanation | Quantify over the unconstrained FE spaces; one notation throughout |
| 4.43 | +(u⊗u, ∇v) set equal to a negative four-term expansion | All four terms of the expansion have the same sign as the left side; the conservative convective form is −(u⊗u, ∇v\_h) |

## Sources

All equation numbers refer to these documents as uploaded.

- R. Codina, S. Badia, J. Baiges, J. Principe. Variational multiscale methods in computational fluid dynamics. Encyclopedia of Computational Mechanics, 2nd ed., 2017. Section 5.
- R. Codina, J. Principe, O. Guasch, S. Badia. Time dependent subscales in the stabilized finite element approximation of incompressible flow problems. CMAME 196 (2007) 2413–2430.
- R. Codina. Analysis of a stabilized finite element approximation of the Oseen equations using orthogonal subscales. Applied Numerical Mathematics 58 (2008) 264–283.
- J. Principe, R. Codina, F. Henke. The dissipative structure of variational multiscale methods for incompressible flows. CMAME 199 (2010) 791–801.
- O. Colomés, S. Badia, R. Codina, J. Principe. Assessment of variational multiscale models for the large eddy simulation of turbulent incompressible flows. CMAME 285 (2015) 32–63.
- E. Castillo, R. Codina. Dynamic term-by-term stabilized finite element formulation using orthogonal subgrid-scales for the incompressible Navier–Stokes problem. CMAME 349 (2019) 701–721.
- A. Kabalan Ammar. Implementing and testing different formulations of the convective term in HPC finite element code using subgrid scales. MSc thesis, UPC, 2023.
