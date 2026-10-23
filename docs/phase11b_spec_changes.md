# Phase 11b — spec changes after the Chapter 5 revision

Source of truth: the revised chapter, `chapter5_preview.pdf`
(SHA-256 `720be75921b3956325668b85cac1641428b0fb626e9e4ca4bf67e5d248b014bc`, 19 pages),
built from `chapter5.tex`
(SHA-256 `409d6689f02e8ad6c7b8ba085f1325e86a5e677921e018cea7a02dce3bbc5de6`).
Copy it over the old file in `~/04_applications/FEops_prep/VMS_masters_thesis/`.

The revision adopts the Phase 11a findings and replaces the EMAC–VMS closure with a
simpler one. Reason for the switch: the subscale model is then exactly the standard
nonlinear one, so EMAC–VMS changes only the resolved equation and reuses the subscale
solver. It does **not** remove the singularity or OSS-orthogonality properties of the
local update; those are shared with the standard nonlinear update (item 1b), as your
Phase 11a follow-up pointed out. Edit `docs/formulation_spec.md` as below, then verify the edits
independently as in Phase 11a (re-derive; do not trust this list or the chapter).
No code, no heavy computations, do not touch `results/`.

## 1. EMAC–VMS: replace the closure (chapter Definition 5.1, (5.23)–(5.26))

The M_h variant is now chapter §5.7.6 and is **not** implemented. Replace the section
"EMAC–VMS (chapter Definition 5.1)" as follows.

- **Subgrid stress** (5.23), no isotropic part:
  T̃ = ũ⊗u_h + u_h⊗ũ + χ_N ũ⊗ũ.
  UFL: `outer(us, uh) + outer(uh, us) + chiN*outer(us, us)`, contracted with
  `inner(T, sym(grad(v)))`. T̃ is symmetric, so the order inside `outer` is immaterial here.
- **Resolved equation**: unchanged from the current spec except for T̃. Keep
  c_E(u_h; u_h, v_h) with advection velocity u_h, the stress-form viscous term, and the rule
  that the standard convective coupling −(ũ, (a·∇)v_h) must not be added as well.
  Equivalent form (chapter (5.28)), useful as a test:
  G(v_h) = c_E(u_h; u_h, v_h) − (ũ, (a·∇)v_h + (∇v_h)u_h), with a = u_h + χ_N ũ.
  UFL: `inner(us, dot(grad(v), a) + dot(uh, grad(v)))`.
- **Residual** (5.26): the master formulation's residual written with the physical
  pressure,
  R_u = f − χ_A D_t u_h − (a·∇)u_h + L_ν u_h − ∇P_h − (∇u_h)u_h.
  UFL for the convective part: `dot(grad(uh), a) + dot(uh, grad(uh))`. This is the same
  residual operator as the existing EMAC row of the convective-forms table. Evaluate
  (∇u_h)u_h pointwise; do not project ½|u_h|² onto Q_h.
- **Linear splitting** (χ_N = 0, chapter (5.32)): T̃ = ũ⊗u_h + u_h⊗ũ, a = u_h in the residual.
- **Delete**: N_E, M_h, the local update (I + τ_t M_h)ũ = …, its singularity paragraph,
  and the sentence that the OSS projection keeps ũ orthogonal at every iterate. The
  subscale update is the master formulation's, unchanged, with the residual above.
- **Write into the definition** (pressure convention and energy proof, so the section is
  self-contained):
  - the resolved equation uses the modified pressure P_h;
  - the residual uses the physical pressure p_h = P_h + ½|u_h|², implemented as
    ∇P_h + (∇u_h)u_h;
  - the energy cancellation is the pointwise identity
    T̃ : ∇ˢu_h = 2 ũ·(∇ˢu_h)u_h + χ_N ũ·(∇ˢu_h)ũ = ũ·[(a·∇)u_h + (∇u_h)u_h],
    with a = u_h + χ_N ũ, using 2(∇ˢu_h)u_h = (u_h·∇)u_h + (∇u_h)u_h and
    ũ·(∇ˢu_h)ũ = ũ·(ũ·∇)u_h. This covers both splittings; the ũ in T̃ and in a must
    be the same field for the identity to hold pointwise, which is the case at
    convergence;
    together with c_E(u_h; u_h, u_h) = 0 this gives G(u_h) + (ũ, N) = 0 (chapter
    Proposition 5.8).
- **Assembly of T̃ (owner decision, from your pre-11b code review).** Build T̃
  directly from the substituted subscale expression in all slots:
  T̃ = ũ_expr⊗u_h + u_h⊗ũ_expr + χ_N ũ_expr⊗ũ_expr.
  Differentiate through the expression for Newton. Do not assemble the expansion (5.28)
  with the stored advector a = u_h + ũ*: that mixes ũ* and ũ_expr, giving ũ_expr⊗ũ*,
  which is not symmetric during the iteration. The residual keeps the master
  convention, a* = u_h + ũ* in (a·∇)u_h and τ. The resolved EMAC term uses u_h alone.
  The expansion (5.28) remains valid as an algebraic test with one ũ field (T1).
- **Linearization (owner decision).** emac and emac_vms use Newton only in Phase 12a,
  with no Picard warm-up. Freezing the advector does not linearize the substituted
  stress or the kinetic-pressure term. If a later case needs warm-up, specify a
  symmetric affine form then, e.g. chapter (5.38) for c_E with ũ* in every slot of T̃.
- **Iteration rule**: apart from the assembly of T̃ above, the subscale is handled exactly
  as in the master formulation.
  Replace "T̃ is then linear in u_h and symmetric at every iterate, so the momentum
  statements hold at every iterate" by: T̃ is symmetric at every iterate (affine in u_h if
  ũ is lagged; nonlinear but still symmetric, with symmetric Newton derivative, if ũ is
  substituted as an expression, chapter Remark 5.9). Hence the convective part of every
  linearized resolved equation tested with e_k or φ is fluxes only (chapter
  Proposition 5.10). The full momentum balances (ASGS subscale time term) and the energy
  balance hold at convergence of the coupled iteration.
- Convective-forms table, EMAC–VMS row: residual operator (a·∇)u_h + (∇u_h)u_h, pressure
  unknown P_h, uses ∇·a: no.

## 1b. Master formulation: local-update qualification (pre-existing)

This applies to the standard nonlinear update the spec already has, and hence to every
nonlinear case, EMAC–VMS included. The residual term (ũ·∇)u_h = G ũ, with
G_ij = ∂_j u_{h,i} (UFL `grad(uh)`), makes the implicit local update
(I + τ_t G)ũ = τ_t(…).

- **Solvability.** It can become singular. In 2D with ∇·u_h = 0,
  det(I + τ_t G) = det(I + τ_t ∇ˢu_h) + τ_t²ω_h²/4. So in pure strain it fails exactly
  like the deleted M_h system; this is the mechanism already measured at the cavity lid
  corners. State it for the dynamic update too, not only for static nonlinear subscales.
- **OSS orthogonality.** With the projection held fixed during local sweeps,
  (ũ_new, w_h) = −(τ_t G(ũ_new − ũ*), w_h), as in H-D3 with G in place of M_h.
  Promise orthogonality only for the matched expression after a converged inner solve
  (and at outer convergence). Store that final matched expression. Do not invoke the OSS
  cancellations during unmatched sweeps.
- **Optional exact fix**, not needed for conservation at convergence and not part of
  Phase 12: weight the projection with W = (τ_t⁻¹I + G)⁻¹ in place of τ_t.
  - W is not symmetric, so this is an oblique projection.
  - It is well posed if the symmetric part of W is positive definite.

Chapter §5.7.6 now states this comparison (same system with B = M_h or B = (∇u_h)ᵀ),
and Remark 5.8(3) says EMAC–VMS shares these properties with the standard formulation.

## 1c. Consistency wording

Do not describe EMAC–VMS as "exactly Codina's scheme for divergence-free fields".
Chapter Proposition 5.6(ii) compares convective **forms** under (A1) and pointwise
∇·(u_h + ũ) = 0: they differ by ((∇·U)u_h, v) and the gradient (½∇|u_h|², v), up to
fluxes. The discrete schemes are different, because:
- ∇·U is only orthogonal to Q_h;
- ½|u_h|² is generally not in Q_h, so the gradient is not absorbed by the discrete pressure.

## 2. Outflow paragraph

Replace "EMAC and EMAC–VMS integrate nothing by parts" by: the EMAC resolved operator is a
volume form, so its natural condition comes from the viscous and pressure terms,
(2ν∇ˢu − P I)n = 0. For EMAC–VMS, −(T̃, ∇ˢv) would add −T̃n to the traction if ũ had a
trace; under the zero-subscale-boundary convention, T̃n = 0 and the condition is that of
EMAC (chapter Remark 5.12). Physical-pressure form:
(2ν∇ˢu − p I)n + ½|u|²n = 0.

## 3. Conservation checks

- **Closure residual (H24).** Define r = rate + D + 𝒩_pred − power and require r = 0 for
  **every** form. Report 𝒩_meas = power − rate − D and compare it with 𝒩_pred. D includes
  2ν‖∇ˢu_h‖² or ν‖∇u_h‖², (τ_m⁻¹ũ, ũ), −2(ũ, L_ν u_h) and (τ_c P̃(∇·u_h), ∇·u_h), all at
  the half level. Predicted values:
  - nonlinear splitting: Table 5.1, energy column; EMAC (ũ, (u_h·∇)u_h); EMAC–VMS 0;
  - linear splitting: Table 5.2; EMAC ½(ũ, ∇|u_h|²); EMAC–VMS 0.

  The resolved-scale EMAC production (ũ, ω_h × u_h) is a different quantity
  (Proposition 5.3) and is not 𝒩.
- **Boundary terms (H12, H-D1).** Replace "vanish when u_h = 0 on ∂Ω or the domain is
  periodic" by:
  - with u_h = 0 on all of ∂Ω, all convective fluxes vanish;
  - on periodic domains, the linear-momentum and energy fluxes cancel, but the
    angular-momentum fluxes do not, because φ = e_k × x is not periodic
    (counterexample H-D1).

  Run angular-momentum checks only with u_h = 0 on ∂Ω, or on compactly supported flows
  (Gresho vortex, solution and subscales zero in a strip of elements).
- **Defects vs balances.** Distinguish the algebraic defect functionals, which can be
  evaluated for any state, from the momentum balances of the solved equations. The
  balances need admissible tests, i.e. the strip setting (chapter §5.7.3, "Scope").
- **OSS.** ∫ũ = 0 and ∫x × ũ = 0 hold on the whole domain only, not on patches (C05).
- **Time schemes.** "The linear-momentum balance holds per step for every time scheme"
  becomes: it holds per step if u_h and ũ use the same integrator, or with OSS on the
  whole domain. Otherwise the conserved discrete quantity mixes the two schemes.
- **(A2).** Keep the τ-weighted-projection statement; it now matches chapter Remark 5.6.

## 4. Tests to add or adjust

Algebraic, cheap; no results/ access.

1. **T1.** G(v) assembled from the stress form equals G(v) from the expansion (5.28), to
   round-off, for random u_h and random quadrature-point ũ.
2. **T2.** EMAC–VMS defect functionals vanish for random u_h (zero on ∂Ω) and arbitrary
   ũ, including ũ nonzero at boundary-adjacent quadrature points:
   - G(e_k) = 0;
   - G(φ) = 0;
   - G(u_h) + (ũ, N) = 0 with N = (a·∇)u_h + (∇u_h)u_h.
3. **T3.** The pointwise identity T̃ : ∇ˢu_h = ũ · N holds at every quadrature point.
4. **T4.** With ũ = 0, G(v) = c_E(u_h; u_h, v).
5. **T5.** The closure residual r of §3 is zero to solver tolerance for all forms on a
   Crank–Nicolson run, from step 2 onwards.
6. **T6.** The periodic angular-momentum counterexample of H-D1 is reproduced, as a guard
   against reintroducing the periodic claim.

## 5. Chapter references used by the spec (old → new)

| Item | Old | New |
| --- | --- | --- |
| EMAC–VMS definition | Definition 5.1 | Definition 5.1 (new content) |
| Subgrid stress | (5.23) | (5.23) (new content) |
| Residual | (5.26) with N_E | (5.26) with N = (a·∇)u_h + (∇u_h)u_h |
| Linear splitting | (5.29) | (5.32) |
| EMAC–VMS energy balance | (5.28) | (5.31) |
| Total balance, Corollary | (5.20), Corollary 5.1 | unchanged |
| Projection assumption / (A2) | — | Remark 5.6 |
| Consistency | — | Proposition 5.6 |
| Momenta, energy, linear, iteration | Propositions 5.6, 5.7, 5.8, 5.9 | Propositions 5.7, 5.8, 5.9, 5.10 |
| M_h variant | Definition 5.1, Remark 5.9 | §5.7.6, (5.33)–(5.35), not implemented |
| Time integration | Remark 5.8 | Remark 5.10 |
| Stabilization | Remark 5.10 | Remark 5.11 |
| Boundary conditions and pressure | — | Remark 5.12 |
| Implementation | §5.9.3, (5.37) | §5.9.3, (5.41) |
| Forms, Remarks 5.1–5.3, Lemma 5.2, Tables 5.1–5.2 | unchanged | unchanged |

## 6. Disposition of the Phase 11a chapter issues

All are fixed in the revised chapter; verify each. The follow-up corrections (local
update shared with the standard formulation; forms vs discrete schemes in Proposition 5.6)
are in §5.7.2, §5.7.6 and §5.9.3.

- **C01:** introduction now separates the resolved production (ũ, ω_h × u_h) from the
  total defect (ũ, (u_h·∇)u_h).
- **C02:** Remark 5.12 states the zero-trace convention and the −T̃n term without it.
- **C03:** §5.9.1 gives the skew-1 and conservative outflow conditions and the corrective
  facet terms.
- **C04:** only the resolved term is linearized (Newton or the symmetric Picard (5.38)).
  T̃ uses the lagged subscale. Other symmetric linearizations are noted as different
  iterations with the same fixed point.
- **C05:** §5.3, §5.4 and §5.7.3 restrict the OSS statements to R = Ω.
- **C06:** Remark 5.10 distinguishes BDF1, BDF2 (modified energy) and BDF3.
- **C07:** Remark 5.8(3) and Remark 5.11. With the new closure the residual is literally
  the standard one with the physical pressure. The τ⁻¹‖ũ‖² ≈ τ‖P̃R‖² link is marked as a
  quasi-static heuristic.
- **C08:** moot for the implemented formulation; §5.7.6 records the W-weighted projection
  for the M_h variant.

Also changed:
- patch momentum statements are labelled formal, with global conservation under the strip
  conditions (§5.7.3);
- Proposition 5.4 is restated with the projection property and weighted dissipation terms.

Output as usual: a verification pass at the top of `docs/verification/VERIFICATION.md`, a
progress-log entry, and a commit "Phase 11b: EMAC–VMS minimal closure and conservation-check
corrections". Stop and report.
