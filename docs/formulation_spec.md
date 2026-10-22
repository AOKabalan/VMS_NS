# VMS formulations for incompressible Navier–Stokes: implementation spec

Oct 21, 2026 · @Ali Kabalan

Every VMS variant is one master system at the solve level with three switches: subscale space (ASGS or OSS), subscale dynamics (static or dynamic) and scale splitting (linear or nonlinear). A fourth switch selects the inner linearization (Newton or Picard); it changes the iterations, not the converged solution. All eight cases follow from these switches with no further choices; the Galerkin baseline is the same system with every stabilization term removed. Revised after an independent re-derivation (VERIFICATION.md, 21 October 2026); split OSS was removed because the thesis does not use it.

## Fixed choices

Every modelling decision is fixed once; "Ours" marks a choice that goes beyond what the cited source states.

| Decision | Choice | Basis |
| --- | --- | --- |
| Variants | Galerkin baseline without stabilization, and ASGS or OSS × static or dynamic × linear or nonlinear splitting: 8 VMS cases | Ours. Colomés §§2.2.1–2.2.3 define the three choices; their Table 1 benchmarks only six (no static-nonlinear cases) |
| Inner linearization | Newton or Picard, selectable; same converged solution | Ours. Colomés §4 uses Picard; thesis §3.2 compares both |
| Interpolation | VMS: equal order Pk/Pk, k = 1 first, k = 2 when the viscous residual must be active. Galerkin baseline: Taylor–Hood P(k+1)/Pk | Ours. Admissible in all sources (Codina 2008 also allows unequal orders); Δ\_h vanishes for affine P1 |
| Galerkin convective form | Skew-symmetric type 1; on outflow add ½⟨(u\_h·n)u\_h, v\_h⟩ | Volume form: Colomés eq. 10 and Remark 2.1 (needs no derivative of the discontinuous ũ; volume terms energy-neutral). Outflow term: ours |
| Viscous form | ν(∇u\_h, ∇v\_h); residual uses ν Δ\_h u\_h | All sources (volume terms) |
| Outflow condition | Gradient-form do-nothing, ν ∂\_n u − p n = 0, no traction load | Ours; see Boundary conditions and quadrature |
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

Verify every variant on periodic and all-Dirichlet problems first, where no boundary terms arise; the outflow condition is used only for the cylinder.

| Boundary | Velocity | Pressure | Terms added |
| --- | --- | --- | --- |
| Periodic | Periodic FE space | Periodic; mean fixed to zero | None |
| Dirichlet Γ\_D | u\_h = g strongly at the solve level (endpoint average for Crank–Nicolson) | Mean fixed to zero if Γ\_D is the whole boundary | None; ũ is not constrained (bubble assumption) |
| Outflow Γ\_N | Natural | Natural | ½⟨(u\_h·n) u\_h, v\_h⟩ in c; no traction load |

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

The projections are τ-weighted L² projections onto the unconstrained FE spaces, with weight τ\_t for the velocity so that the dynamic subscale stays exactly orthogonal.

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
- ξ\_h and ζ\_h are unknowns of the inner solve, together with u\_h and p\_h (Codina 2008 §2.3; Colomés §5.2.2). The pressure null space applies only to p\_h; the projections stay unconstrained, and ζ\_h with its equation is absent when τ\_c is off. Lagging them by one outer iteration, as in Colomés §4.2, has the same fixed point but was measured to fail here: with a frozen Kovasznay advector, the error contraction per outer iteration was 0.85–0.97 for P1/P1, and for static P2/P2 reached 1.19 in magnitude (divergent); half relaxation and lagging the full residual did not fix it (implementation probe and independent review, 22 October 2026). Shrinking δt does not help, because the pressure has no time-derivative term.
- This ξ\_h is not the same field as Colomés' ξ\_m, which projects the full residual: ξ\_m = ξ\_h − D\_t u\_h. The converged solutions agree; iteration counts are not comparable with the paper's lagged algorithm.
- R̂\_u keeps f and the elementwise viscous term νΔ\_h u\_h. Both cancel or vanish in common cases (f in the FE space; P1 elements), and dropping νΔ\_h does not reduce the convergence order, but keeping them is what makes the formulation exactly residual based.
- ζ\_h only matters when τ\_c is active. With τ\_c = 0 (turbulence settings) the pressure subscale, its projection and the grad-div term disappear.

## Inner linearization: Newton and Picard

Both solve the same discrete equation F(u\_h, p\_h) = 0 at fixed outer data, so when both converge to the same solution the time step is identical; they differ only in how the advection velocity and τ are treated during the inner iterations.

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

Recommended use: a few Picard iterations, then Newton once the residual has dropped (thesis §3.2.2). Picard with one inner iteration per outer iteration has the same structure as the Colomés §4 algorithm (one linearized solve, then the local subscale update and the projection update), but for OSS the iterates differ, because Colomés lag the projections while this spec solves them with the FE unknowns. The fixed points agree. The thesis found Picard converging poorly with EMAC; that belongs to the later convective-term study (Charnyi et al. 2019 analyse the EMAC linearizations), not to this baseline.

## Time-step algorithm

One time step is an inner solve for (u\_h, p\_h), Newton or Picard, inside an outer fixed point on the stored subscale; the four linear-splitting cases need only the inner solve.

Given the velocity and subscale histories the scheme needs (back to level n−2 for BDF3), and ξ\_hⁿ, ζ\_hⁿ:

1. Form u\_hist and ũ\_hist from the Time schemes table. Initialise (u\_h, p\_h) from the previous step, the stored subscale ũ\* ← ũⁿ, and ξ\_h, ζ\_h ← their previous values.
2. Inner solve: Newton or Picard on the master FE equation with ũ substituted, as in the linearization section. Frozen: ũ\* inside a (nonlinear splitting) and ũ\_hist. For OSS, ξ\_h and ζ\_h are unknowns of this solve, with their projection equations. In Picard mode, a single inner iteration per outer iteration is allowed only in cases that run step 5; the two linear-splitting ASGS cases have no outer loop, so step 2 itself must iterate until the FE residual is below tolerance.
3. Subscale update (nonlinear splitting only): at each quadrature point iterate ũ\* ← τ\_t(a) (R + χ\_D ũ\_hist/δt − ξ\_h) with a = u\_h + ũ\*, where R = R\_u for ASGS and R = R̂\_u for OSS, until the pointwise change is below tolerance. The accelerated form below is equivalent.
4. Projections (OSS only): nothing to compute, since ξ\_h and ζ\_h come out of step 2; in nonlinear-splitting OSS cases their changes still enter the convergence test of step 5.
5. Convergence: stop when the relative changes of u\_h, p\_h, ũ\*, ξ\_h and ζ\_h (when τ\_c is active) and the FE residual are all below tolerance; otherwise return to step 2.
6. Store: for dynamic subscales, ũ^(n+1) is the subscale expression evaluated at the converged state; shift the histories of u\_h, ũ, ξ\_h and ζ\_h.

Accelerated subscale update for step 3, a d×d linear solve per quadrature point, with τ\_t evaluated at the previous ũ\* and (∇u\_h)\_ij = ∂u\_h,i/∂x\_j:

```latex
\big(I + \tau_t\, \nabla u_h\big)\, \tilde u = \tau_t \left(f - \chi_A\, D_t u_h - u_h\cdot\nabla u_h + \nu\,\Delta_h u_h - \nabla p_h + \chi_D\,\frac{\tilde u_{\mathrm{hist}}}{\delta t} - \xi_h\right)
```

This is Colomés eqs. 50 and 52 with the τ\_t factor restored on the left, ξ\_h multiplied by a single τ\_t (not τ\_t² as in thesis eq. 4.39), and χ\_A dropping D\_t u\_h for OSS.

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
