**Phase 11b — minimal EMAC–VMS closure and conservation corrections — 23 October 2026**

Reviewed the revised 19-page `~/04_applications/FEops_prep/VMS_masters_thesis/chapter5_preview.pdf` after first checking its SHA-256: `720be75921b3956325668b85cac1641428b0fb626e9e4ca4bf67e5d248b014bc`, exactly the supplied value. Applied items 1–5 of `docs/phase11b_spec_changes.md`, SHA-256 `4bb90f5d7784b6d6253bff03ba53ef71e10a5a092e3e886d2cec07b1308729e2`, including the explicit owner decisions. The resulting `docs/formulation_spec.md` has SHA-256 `f7912ddf478fc276284a061b7d66123b8102ce7dd1d7269e4fcab059c1b339b0`. The supplied change list is retained unchanged. The chapter's equations were re-derived, not accepted on agreement with the change list. PDF pages 12–13 were also inspected as rendered images.

**Verdict: the updated spec passes. No remaining discrepancy requires changing its equations, signs or iteration.** The minimal closure has the claimed total convective energy cancellation for both splittings, reuses the physical-pressure master residual, and removes convective force and torque up to the stated boundary fluxes. It retains the master's possible local singularity and unmatched-sweep OSS orthogonality defect. The gate has no equation-level blocker; work stops after Phase 11b as explicitly requested. Two remaining chapter prose errors are recorded separately below.

No implementation or plan changes, simulations, repository test suites or access to `results/`. A process check found no running validation queue, validation runner, MPI or pytest processes. Independent verification used component algebra, integration by parts and the discrete scale equations, supplemented by a small standalone exact-rational check in `/tmp`: 80 tensor/pressure/midpoint cases across 2D/3D and both splittings, 40 determinant cases, and the two local-update counterexamples below all passed. These are verification probes, not the future implementation tests T1–T6.

Notation: \(u=u_h\), \(s=\tilde u\), \(\chi=\chi_N\), \(A_{ij}=\partial_j u_i\) (UFL gradient), \(D=(A+A^T)/2\), \(d=\operatorname{tr}A\), and \(a=u+\chi s\). Thus the chapter gradient is \(A^T\). \(G(v)\) denotes the minimal resolved convective functional, not the matrix called G in the spec's local update. All discrete products use the common quadrature.

| Item | Spec says | Source says / independent derivation | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| J01 — index convention | Retain all five UFL translations. | Advection is Aa; chapter-gradient times a is Aᵀa; the nonsymmetric contraction is uⱼaᵢ∂ᵢvⱼ. Symmetrization removes the transpose difference. | Agree | The Phase 11a U10 conclusions remain valid; the new residual uses both distinct contractions. |
| J02 — minimal stress and expansion | \(T=s\otimes u+u\otimes s+\chi s\otimes s\); replace both convective terms. | \(T:\nabla^s v=s\cdot[(a\cdot\nabla)v+(\nabla v)u]\) by indices, for one common s. | Agree | Exactly (5.23), (5.25), (5.28), with χ=0 giving (5.32). No isotropic part or second standard coupling. |
| J03 — pointwise energy | \(T:D=s\cdot N\), \(N=Aa+A^Tu\). | \(T:D=2s\cdot Du+\chi s\cdot Ds=s\cdot[(A+A^T)u+\chi As]\). | Agree | χ multiplies the quadratic contraction. No divergence, projection, bubble or integration-by-parts assumption. |
| J04 — integrated energy and momenta | \(G(u)+(s,N)=0\), \(G(e_k)=G(\phi)=0\) under the stated flux conditions. | The energy sum is \(B(u)\); the two rigid tests give \(B(v)+H_b(v)\). Definitions and derivation in J-D1. | Agree with explicit boundary scope | These zeros are not unconditional on open or periodic angular domains. Arbitrary quadrature-point s is allowed, including at boundary-adjacent points. |
| J05 — linear splitting and zero subscale | Drop only the Reynolds stress and set a=u; no outer subscale loop. | \(T_{\rm lin}:D=s\cdot(A+A^T)u\), and s=0 gives exactly \(c_E(u;u,v)\). | Agree | Both conservation identities survive; (5.32) is not the former full-stress linear model. |
| J06 — physical residual | Use modified P in the resolved equation and physical \(p=P+\lvert u\rvert^2/2\) in the residual. | \(\nabla p=\nabla P+A^Tu\), so \(N+\nabla P=Aa+\nabla p\). | Agree | Exactly the master subscale equation with the selected stress viscosity. For OSS, dropping \(D_tu\) is an equivalent projected representative, not equality of the unreduced residual fields. |
| J07 — pressure reporting and midpoint | Null space on P; recover and normalize physical p separately in a rich enough space; use averaged velocity for CN. | Endpoint-average kinetic pressure exceeds the correct half-level correction by \(\lvert u^{n+1}-u^n\rvert^2/8\). | Agree | Projection of the correction into Qh changes the physical field and generally its gradient. Pressure itself is not extrapolated. |
| J08 — stress assembly and Newton | Substitute \(s_{\rm expr}\) in every T slot; residual and τ use stored s*. | T and its full derivative are symmetric; the substituted tensor is nonlinear. The unmatched energy remainder is \(\chi(s_{\rm expr},A(s^*-s_{\rm expr}))\). | Agree | J-D2. At the coupled fixed point this remainder is zero. This is the eliminated-subscale option of Remark 5.9, not the fully lagged matrix (5.41) during iterations. |
| J09 — master local update | \((I+\tau_tA)s=\tau_t(X-\xi)\). | Move \((s\cdot\nabla)u=As\) from the residual to the left; \(\tau_t^{-1}=\chi_D\alpha_0/\delta t+\tau_m^{-1}\). | Agree | Same pointwise system for all nonlinear forms, with physical pressure. No new M-based solve. Dynamic subscales can also fail. |
| J10 — determinant and solvability | In 2D trace-free flow, determinant differs from the symmetric-strain determinant by \(\tau_t^2\omega^2/4\). | Set \(A=\left[\begin{smallmatrix}r&b\\c&-r\end{smallmatrix}\right]\); subtraction gives \(\tau_t^2(c-b)^2/4\). | Agree | Pure strain gives identical failure. The chapter's further conditioning claim does not follow; C09 below. |
| J11 — OSS orthogonality | Matched expression is orthogonal after a converged inner solve, not during unmatched local sweeps. | At fixed τ, \((s_{\rm new},w)=-(\tau_t A(s_{\rm new}-s^*),w)\); variable τ adds the weight-change term. | Agree with the fixed-coefficient qualification added to item 1b | J-D3. Store the final matched expression; CN extrapolation preserves endpoint orthogonality. No claim for raw Newton iterates. |
| J12 — optional matrix projection | W is optional, oblique in general, outside Phase 12; positive symmetric part suffices. | \((W\xi,w)=(WX,w)\) makes \(s=W(X-\xi)\) orthogonal. Positive \(\operatorname{sym}W\) makes the finite projection matrix coercive under the stated quadrature/rank conditions. | Agree | Nonsymmetric W supplies no inner-product norm. The inverse and matched coefficients are prerequisites. The revised PDF now states this correctly. |
| J13 — consistency | Continuous convective forms can agree after a pressure change; discrete schemes need not. | The difference from the standard nonlinear coupling is \(((\nabla\cdot U)u,v)+(\nabla(\lvert u\rvert^2/2),v)\), up to fluxes, U=u+s. | Agree | Weak continuity is not pointwise solenoidality; the kinetic correction generally lies outside Qh. J-D4 re-derives Proposition 5.6. |
| J14 — outflow and viscosity | The first four forms have physical-pressure traction; EMAC uses P; subgrid traction vanishes under the zero-trace convention. | The standard corrective facets cancel their weak-convection surface terms. Minimal stress would contribute \(-Tn\); \(L_\nu u=\nu(\Delta u+\nabla d)\) is formally self-adjoint. | Agree | Stress residual has \(+L_\nu u\); adjoint coupling has \(-(s,L_\nu v)\). EMAC traction in p has \(+\lvert u\rvert^2n/2\). |
| J15 — variable weights and (A2) | Constant τ unnecessary for ASGS and matched weighted OSS. | \((s,\xi)=0\), \((\tau_c(d-\zeta),\zeta)=0\); hence the pressure product is \((\tau_c(d-\zeta),d-\zeta)\). | Agree | No coefficient is pulled through a projector or integral. This is the weighted extension, not the unweighted projector multiplied by variable τ. |
| J16 — closure sign and defect | \(r=\mathrm{rate}+\mathcal D+\mathcal N_{\rm pred}-\mathrm{power}=0\) for every form. | Summing the two scale equations gives exactly this sign; \(\mathcal N_{\rm meas}=\mathrm{power}-\mathrm{rate}-\mathcal D\). | Agree; H24 corrected | EMAC nonlinear defect is \((s,(u\cdot\nabla)u)\); linear defect is \((s,A^Tu)\). Resolved Lamb production is different. |
| J17 — exact CN budget | Endpoint total energy and midpoint dissipation/power, dynamic both scales, step 2 onward. | Midpoint polarization gives the exact rate; subscale ODE contributes \(\tau_m^{-1}\), and the two viscous products sum to \(-2(s,L_\nu u)\). | Agree at coupled convergence | J-D5. BDF1 startup has an additional positive increment term. Existing diagnostics retain their meanings. |
| J18 — defect tables and admissible tests | Tables 5.1–5.2 with boundary corrections; periodic linear momentum/energy only. | The unchanged H-D1 product identity cancels arbitrary s pointwise; J-D1 supplies the minimal closure's zero interior row. | Agree | The continuity pairing controls u+s even for linear splitting. Rigid tests on zero-Dirichlet spaces need the strip/cutoff argument; global OSS does not imply patch moments vanish. |
| J19 — chapter implementation and references | Updated equation numbers and eliminated-subscale assembly. | Differentiating the lagged cross stress gives \(-\int N^B(s_l^*\partial_lN^A\delta_{ij}+s_i^*\partial_jN^A)\); the Reynolds RHS is \(+\int s_i^*s_l^*\partial_lN^A\). | Agree | Exactly (5.41); the former isotropic matrix term is gone. K and D for the five forms retain the correct index formulas. Newton through the eliminated stress adds its derivative, explicitly permitted by Remark 5.9 and §5.9.3. |
| J20 — allowed linearization | EMAC and EMAC–VMS use Newton without Picard warm-up in Phase 12a. | The explicit owner decision and change list select this restriction. Frozen advection alone leaves the kinetic-pressure term and substituted stress nonlinear. | Agree | The later task bullet allowing warm-up is superseded by that explicit decision. No new iteration is invented. |

**Independent derivations and qualifications.**

**J-D1 — minimal tensor, fluxes and arbitrary subscales.** For \(H_{ij}=\partial_jv_i\), symmetry gives

\[
T:\operatorname{sym}H=T:H
=s\cdot H(u+\chi s)+s\cdot H^Tu.
\]

Setting H=A yields \(T:D=s\cdot N\). Setting v to either rigid test gives zero pointwise, since \(\operatorname{sym}H=0\). This uses neither a trace nor a derivative of s. Define

\[
B(v)=\int_{\partial\Omega}(u\cdot n)(u\cdot v),\qquad
H_b(v)=\frac12\int_{\partial\Omega}|u|^2(v\cdot n).
\]

Component differentiation shows the resolved EMAC operator is the divergence of \(u\otimes u+\tfrac12|u|^2I\). Therefore

\[
G(e_k)=B(e_k)+H_b(e_k),\qquad
G(\phi)=B(\phi)+H_b(\phi),\qquad
G(u)+(s,N)=B(u).
\]

For the last identity, \(c_E(u;u,u)=\int_\Omega\nabla\cdot(|u|^2u)\). Thus all three vanish with zero boundary u, whatever the quadrature-point s. Periodicity cancels the first and third, but not generally the second. These arguments hold for χ=0 and χ=1. They prove algebraic defects, not availability of eₖ or φ as tests of a Dirichlet solve. The strip assumption supplies admissible cutoff tests on the active region for the full balances.

The all-form table derivation remains \(((a\cdot\nabla)u,v)+(u,(a\cdot\nabla)v)+D_\chi(u\cdot v)=B(v)\), with \(D_\chi(w)=(d,w)-\chi(s,\nabla w)\). Only the polynomial resolved part requires integration by parts. The 12 interior rows and the six boundary-correction rows in H-D1 below remain valid with the new minimal EMAC–VMS zero row; none of the first five forms changed.

**J-D2 — substituted stress versus stored residual advector.** For an increment v in u and the induced increment ρ in \(s_{\rm expr}\),

\[
\delta T=\rho\otimes u+u\otimes\rho+s\otimes v+v\otimes s
+\chi(\rho\otimes s+s\otimes\rho).
\]

Each pair is symmetric, without assuming ρ independent of u, P or ξ. In contrast, using a*=u+χs* in the expanded coupling replaces the quadratic symmetric tensor by a mixed product. Directly, with s denoting the expression used in T,

\[
G(u)+(s,N(a^*))=B(u)+\chi(s,A(s^*-s)).
\]

This is why the energy budget is a check at the coupled fixed point. Newton symmetry preserves the convective rigid-test cancellation during iteration, but does not establish the full time-dependent balance before convergence. Linear splitting has no such stored-advector remainder.

**J-D3 — shared local update, determinant, and OSS.** Write \(X=f-\chi_A D_tu-Au+L_\nu u-\nabla p+\chi_Ds_{\rm hist}/\delta t\). The matched expression is \(s_{\rm expr}=\tau_0(X-As^*-\xi)\), with its projection satisfying \((s_{\rm expr},w)=0\). A local sweep with τ₁ solves \(s_{\rm new}=\tau_1(X-As_{\rm new}-\xi)\), so subtraction gives

\[
(s_{\rm new},w)
=((\tau_1-\tau_0)(X-As^*-\xi),w)
-(\tau_1A(s_{\rm new}-s^*),w).
\]

The shorter formula in item 1b assumes τ₁=τ₀. For two equally weighted points and projection onto constants, take τ=1, ξ=0, s*=0, X=(1,−1), and one active diagonal component of A equal to (0,1). The matched expression has zero mean; the implicit update is (1,−½), with mean ¼. This is already the master update, without M. The example concerns the local projection implication; it does not assume a solved FE momentum equation. At outer convergence the projection and expression match again.

For \(A=\left[\begin{smallmatrix}r&b\\c&-r\end{smallmatrix}\right]\), \(\omega=c-b\),

\[
\det(I+tA)=1-t^2(r^2+bc),\quad
\det(I+tD)=1-t^2\left(r^2+\frac{(b+c)^2}{4}\right).
\]

Their difference is \(t^2\omega^2/4\). This proves the stated determinant identity and pure-strain comparison, not a global conditioning bound. The dynamic \(t=(\alpha_0/\delta t+\tau_m^{-1})^{-1}\) can still reach singular values.

For the optional W rule, the projected system is generally nonsymmetric. Nevertheless \(z^TMz=(Ww_h,w_h)=(\operatorname{sym}Ww_h,w_h)>0\) for nonzero coefficients z under the stated positivity/rank assumptions, proving uniqueness. It is a sufficient condition, not an assertion that W is symmetric or that the condition is necessary.

**J-D4 — residual, pressure and consistency.** Component differentiation gives \(\partial_i(\tfrac12u_ju_j)=u_j\partial_i u_j=(A^Tu)_i\), which proves the physical-pressure residual equality without a discrete pressure projection. Under the zero-trace model and with \(U=u+s\),

\[
c_{nc}(U;u,v)-(s,(U\cdot\nabla)v)
\doteq ((U\cdot\nabla)U,v)+((\nabla\cdot U)s,v),
\]
\[
G(v)\doteq((U\cdot\nabla)U,v)+((\nabla\cdot U)U,v)
+(\nabla(|u|^2/2),v).
\]

Here \(\doteq\) suppresses boundary fluxes. Subtracting proves Proposition 5.6(ii) and the spec's consistency wording. The divergence pairing only gives \(D_1(q_h)=0\), not \(D_1(u\cdot v)=0\), and it cannot put the kinetic-pressure correction into Qh. This does not identify the two discrete schemes.

**J-D5 — total-energy signs and time level.** Test the resolved equation with u and its pressure equation with the pressure unknown; dot the subscale equation with s and add. The pressure-gradient products cancel, the viscous cross product appears twice with a minus sign, and the forcing is \((f,u+s)\). ASGS gives the time product \((D_t(u+s),u+s)\); weighted OSS gives \((D_tu,u)+(D_ts,s)\), with the cross terms zero at matched states in the fixed FE space. Hence

\[
\mathrm{rate}
+a_\nu(u,u)+(\tau_m^{-1}s,s)-2(s,L_\nu u)
+(\tau_c(d-\zeta),d)+\mathcal N=(f,u+s).
\]

No derivative of τ is introduced. For midpoint states \((z^{n+1}-z^n,\bar z)=(\|z^{n+1}\|_h^2-\|z^n\|_h^2)/2\), proving the per-step formula, its signs, and the need for separate endpoint and solve-level state. For EMAC–VMS \(\mathcal N=0\) by J-D1; for standard EMAC the residual leaves the additional \((s,A^Tu)\) beyond the resolved convective energy. Combining this with the H-D1 product rule gives \((s,Au)\) for nonlinear and \((s,A^Tu)\) for linear splitting. Both skew forms have zero interior defect. The residual including the predicted defect is zero for every form.

**Disagreements and dispositions for the spec/change list.**

- The former spec's isotropic stress and M-based residual/update are replaced, not patched. H19/H-D2 (linearity), H21/H-D3 (orthogonality), H12 (periodic angular flux) and H24 (closure sign) are corrected. The master limitation remains explicitly documented; no solver change or W projection is silently selected.
- Item 1b's orthogonality subtraction needs fixed τ. The spec now states that hypothesis and includes the extra term when local sweeps recompute τ. This qualifies the diagnostic identity without changing the requested iteration.
- The change list's C08 disposition “moot” is too strong: the same mechanism survives in the standard update. The revised chapter itself acknowledges it, and the updated spec follows that qualification.
- Zero defects require the stated flux conditions. No bubble condition on arbitrary s is needed for the minimal tensor tests, whereas a zero-trace model is needed to suppress its natural traction when interpreting the weak equation.
- The explicit owner decision selects Newton without warm-up despite the later generic Part 2 bullet allowing warm-up. The spec records the explicit decision. No implementation choice remains open from this conflict.

**Chapter issues and dispositions, separate from the spec.**

| Issue | Revised chapter location | Independent disposition |
| --- | --- | --- |
| C01 — resolved versus total EMAC energy defect | Introduction, p. 2; Proposition 5.4 | Fixed. The Lamb-vector production and total defect are distinguished; J-D5 confirms the latter. |
| C02 — subgrid boundary term | Remark 5.12, p. 14 | Fixed. Zero trace suppresses Tn; without it the traction has −Tn. The physical-pressure form and midpoint correction have the correct signs. |
| C03 — skew/conservative outflow | §5.9.1, p. 16 | Fixed. The half/full advective natural-traction changes and compensating facet terms are now explicit. |
| C04 — different symmetric iterations | §5.7.7, p. 13; §5.9.3, p. 17 | Fixed. Only the resolved term is linearized in the displayed fully lagged-subscale recipe; (5.38) and other symmetric iterations are distinguished. The eliminated-subscale option is separately acknowledged. |
| C05 — OSS moments on patches | §§5.3–5.4; scope paragraph §5.7.3, p. 11 | Fixed with scope qualifications. Orthogonality is global and depends on projection-space membership; patch statements are formal. Periodic φ is explicitly excluded. |
| C06 — BDF dissipation | Remark 5.10, p. 14 | Fixed. BDF1 adds physical-energy dissipation; BDF2 controls a two-level modified energy; BDF3 has no general identity of that kind. |
| C07 — physical residual and residual bound | Remark 5.8(3), p. 10; Remark 5.11, p. 14 | Fixed for the minimal closure. The residual really is standard with physical p; the τ-residual norm relation is marked quasi-static and heuristic. |
| C08 — local update versus OSS projection | §5.7.6, pp. 12–13 | Addressed as a shared limitation, not eliminated. Both matrices can fail; fixed-projection implicit sweeps need not be orthogonal. The optional W rule is now correctly described as oblique, with a positive symmetric-part sufficient condition. |
| Phase 11a follow-up — shared master properties | Remark 5.8(3), §5.7.6, §5.9.3 | Fixed. Minimal closure reuses the standard local system, including its failure and projection properties. |
| Phase 11a follow-up — forms versus schemes | Proposition 5.6 and following paragraph, p. 10 | Fixed. J-D4 reproduces both extra terms; the text explicitly rejects equality of discrete schemes. |
| Per-iterate full momentum claim | Proposition 5.10 and Remark 5.9, p. 13 | Fixed. The claim is now for the convective part; the full balance requires coupled convergence, including the ASGS time term. |
| C09 — rotation implies better conditioning | §5.7.6, p. 12, immediately after the determinant identity | Still too broad. A positive determinant correction need not improve conditioning and can make an invertible strain matrix singular. Counterexample below. The displayed determinant equation is correct. |
| C10 — resolved energy versus conservative exchange | §5.7.4, p. 12, paragraph after the dissipation bound | Wording needs narrowing. Zero total convective defect does not imply zero resolved-scale production or zero exchange: \(G(u)=-(s,N)\) generally is nonzero. Say that convection transfers energy between scales with zero net contribution to the total balance. The displayed equations already express this correctly. |

For C09, take t=1 and

\[
D=\begin{pmatrix}5/4&0\\0&-5/4\end{pmatrix},\qquad
A=\begin{pmatrix}5/4&-3/4\\3/4&-5/4\end{pmatrix}.
\]

The trace is zero and ω=3/2. The strain matrix I+D is invertible, with determinant −9/16 and 2-norm condition number 9. But I+A is singular: its determinant is zero. The positive correction \(+\omega^2/4=9/16\) moves the determinant to zero. Thus retain the determinant identity and pure-strain comparison, and remove or restrict the general conditioning claim.

C09–C10 require prose corrections only, not a change to the minimal equations, signs or chosen iteration. No chapter file was edited. Phase 12a has not started.

---


**Phase 11a — independent verification of the Chapter 5 addendum — 23 October 2026**

Reviewed `docs/formulation_spec.md`, SHA-256 `2f2022d86c520725a90d8911abeab134d039f7a65abc98b6d9378db1f00f5abf`, against `~/04_applications/FEops_prep/VMS_masters_thesis/chapter5_preview.pdf`, SHA-256 `ee4aed24de3a624f40a2059b1bb2ef0c99cd1ef04aef63735a7af575df765c8c` (17 pages). Spec line numbers below refer to that revision. The check re-derives the identities from indices, integration by parts and the discrete equations; agreement with the chapter alone is not treated as verification. PDF pages 14–15 were also inspected as rendered images to check the implementation matrices and their parentheses. UFL's operator definitions were checked against its [official form-language reference](https://docs.fenicsproject.org/ufl/main/manual/form_language.html#gradient). No implementation code, simulations or test suites were run, and `results/` was not accessed. The spec and chapter were not edited.

**Verdict: the converged volume equations and the entries of Tables 5.1–5.2 pass, with the boundary corrections below. The addendum does not pass as written.** The partial freezing of the subscale does not make the stress linear in the resolved unknown; the local implicit update does not preserve OSS orthogonality with its projection held fixed; periodicity does not remove general angular-momentum fluxes; and the proposed closure diagnostic counts the convective defect incorrectly. The variable-τ extension and the midpoint energy identity themselves are valid. Chapter issues are recorded separately after the spec disagreements.

Notation in this pass: \(u=u_h\), \(s=\tilde u\), \(D u=\nabla^s u\), \(d=\nabla\cdot u\), and all unlabelled volume inner products use the common quadrature. Write \(N_{\rm op}\) for the strong convective residual operator and \(\mathcal N\) for its energy defect, to distinguish the two uses of N in the chapter.

| Item | Spec says | Source says / independent derivation | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| H01 — advection translation, 283 | `dot(grad(u), a)` | Its component j is \((\partial_i u_j)a_i\). | Agree | Exactly \(((a\cdot\nabla)u)_j\). |
| H02 — chapter gradient times vector, 284 | `dot(a, grad(u))` | Its component i is \(a_j\partial_i u_j\). | Agree | Exactly \(((\nabla u)a)_i\) in the chapter convention. |
| H03 — symmetric gradient, 285 | `sym(grad(u))` | Components are \(\tfrac12(\partial_j u_i+\partial_i u_j)\). | Agree | Symmetrization removes the transpose difference. |
| H04 — nonsymmetric outer product, 286 | `inner(outer(u, a), grad(v))` | The contraction is \(u_j a_i\partial_i v_j\). | Agree | Reversing the two factors of `outer` is essential. |
| H05 — symmetric tensor contraction, 287 | `inner(A, sym(grad(v)))` | For symmetric A, relabelling i,j gives \(A_{ij}\partial_i v_j\). | Agree | This row is correct precisely with the stated symmetry. |
| H06 — five volume forms, 293–307 | Forms (5.3)–(5.7) and the specified residual operators. | Direct expansion gives \(c_{ss1}=(c_{nc}+c_c)/2\), \(c_{ss2}=c_{nc}+D_a(u\cdot v)/2\), \(c_E=(2Du\,a,v)+D_a(u\cdot v)\). | Agree | The first four deliberately use the same non-conservative residual; the residual need not be the strong operator of their individual resolved form. H-D1 derives the defect tables without assuming solenoidality. |
| H07 — first four outflow terms, 303–310 | Add \(0,B_N,B_N/2,0\), respectively, to obtain gradient-form do-nothing. | For admissible v, \(c_c=c_{nc}+D_a(u\cdot v)-B_N(v)\), \(c_{ss1}=c_{nc}+D_a(u\cdot v)/2-B_N(v)/2\). | Agree, with the stated zero-subscale-trace model | The additions cancel those surface terms. Thus the resolved Laplacian traction is \(\nu\partial_n u-pn=0\); with stress viscosity it is \((2\nu Du-pI)n=0\). The interior divergence corrections remain. |
| H08 — EMAC outflow, 307–310 | No convective surface addition; traction uses P. | The EMAC resolved operator is already in strong volume form, so the viscous/pressure integration gives \(\nu\partial_n u-Pn=0\), or \((2\nu Du-PI)n=0\). | Agree; EMAC–VMS needs a qualification | In physical pressure the latter is \((2\nu Du-pI)n+\tfrac12\lvert u\rvert ^2n=0\). The subgrid stress in EMAC–VMS is already a weak divergence; its trace matters if the bubble convention is removed (H-D5). |
| H09 — divergence pairing and Lemma 5.2, 315–318 | \(D_a(w)=(d,w)-(s,\nabla w)_h\) for \(a=u+s\). | The assembled q-equation is identically \(D_a(q)=0\); hence \(D_a(w)=D_a(w-\Pi_Qw)\) for any suitable projection into the actual pressure test space. | Agree | This defines a functional, not a pointwise divergence of stored quadrature values. No derivative of s is needed for the table weights. For linear splitting use \(D_0(w)=(d,w)\); continuity still controls \(D_a\), not \(D_0\). |
| H10 — nonlinear defect table, 368–372 | Chapter Table 5.1. | H-D1 derives all eighteen entries, including EMAC's energy defect \((s,(u\cdot\nabla)u)\). | Agree up to the explicit boundary terms | The EMAC resolved-energy Lamb-vector defect and its total-energy defect are different. Conservative momentum cancellations are pointwise even before integration. |
| H11 — linear defect table, 368–372 | Chapter Table 5.2. | H-D1 derives all eighteen entries. The angular defect has the additional \(-J\), \(J=(s,e_k\times u)\), for the five standard forms. | Agree up to the explicit boundary terms | J vanishes for global OSS with componentwise identical velocity spaces. EMAC–VMS removes it also for ASGS. |
| H12 — vanishing boundary terms, 372 | Homogeneous Dirichlet or periodic domain removes the fluxes. | Homogeneous zero u removes the convective fluxes. Periodic u removes fluxes tested with periodic functions, but \(\phi=e_k\times x\) is generally not periodic. | Disagree for angular momentum | H-D1 gives a periodic P2 counterexample. Convective-defect identities must also be distinguished from full momentum laws with wall reactions. |
| H13 — EMAC residual and pressure, 307, 320 | \(N_{\rm op}=(a\cdot\nabla)u+(\nabla u)u\), \(p_h=P_h+\lvert u\rvert ^2/2\). | \((\nabla u)u=\nabla(\lvert u\rvert ^2/2)\), so \(N_{\rm op}+\nabla P_h=(a\cdot\nabla)u+\nabla p_h\). | Agree | Use the physical expression, not merely P, for pressure errors and stresses. A projection into Qh would generally change this identity; \(\lvert u\rvert ^2/2\) need not belong to Qh. At a pressure null space, normalize the recovered physical field separately. |
| H14 — half-level pressure, 320 | Crank–Nicolson uses half-level u. | \(p_h^{n+1/2}=P_h^{n+1/2}+\tfrac12\lvert \bar u\rvert ^2\), \(\bar u=(u^{n+1}+u^n)/2\). | Agree | The average of endpoint kinetic pressures instead exceeds this correction by \(\tfrac18\lvert u^{n+1}-u^n\rvert ^2\). Neither the endpoint u nor an extrapolated pressure is appropriate. Exact manufactured volume forcing stays unchanged; outflow traction data must respect the changed pressure convention. |
| H15 — viscous operator and adjoint, 322 | \(a_\nu=2\nu(Du,Dv)\), \(L_\nu v=2\nu\nabla\cdot Dv\). | For constant ν, \((L_\nu v)_j=\nu(\Delta v_j+\partial_j\nabla\cdot v)\). Two integrations by parts give formal elementwise self-adjointness. | Agree | The momentum operator is \(N_{\rm op}-L_\nu u+\nabla p\); the standard subscale coupling uses \(L_a^*(v,q)=-(a\cdot\nabla)v-L_\nu v-\nabla q\) under the same omitted-divergence/boundary approximation. The strong residual contains \(+L_\nu u\) and the resolved coupling contains \(-(s,L_\nu v)\); do not retain just νΔ in the stress case. |
| H16 — angular viscosity and stability, 322, 376 | Stress form for angular momentum; inverse-estimate restrictions. | \(D\phi=0\), whereas \(\nu(\nabla u,\nabla\phi)=\nu\int(\nabla\times u)_k\). Weighted Young inequalities give the chapter's sufficient constants, with \(c_1k^4\) replacing \(c_1\). | Agree, with scope qualification | Stress is needed for the general patch argument; the Laplacian torque can vanish on the whole zero-boundary or periodic domain. The sufficient bounds do not prove positivity on an arbitrary production mesh; no U5 computation was run here. |
| H17 — EMAC–VMS resolved equation, 334–346 | Replace both resolved convection and the convective subscale coupling. | Adding (5.11) to (5.24), inserting \(\tilde p=-\tau_c(d-\zeta_h)\), gives exactly 334–335, including \(-(s,L_\nu v+\nabla q)\). | Agree | The conventional \(-(s,(a\cdot\nabla)v)\) must not be added a second time. Removing OSS time terms uses orthogonality; static cases are the stated extension of the chapter's dynamic formulation. |
| H18 — stress, residual, \(M_h\), 345–346 | \(\tilde T=DS(s)[u]+\chi_NS(s)\), \(N_E=2Du\,u+du+\chi_NM_hs\), \(M_h=Du+dI/2\). | \(S(u+s)-S(u)=DS(s)[u]+S(s)\); contraction with Du gives \(\tilde T:Du=s\cdot N_E\) pointwise. | Agree | These are (5.23), (5.26); setting \(\chi_N=0\) is exactly (5.29). The OSS residual is the reduced representative after projecting away \(D_tu\), rather than literally the full residual (5.26). |
| H19 — iteration, 352 | Stored s* in τ, \(M_hs\) and quadratic terms makes T linear in u and symmetric. | The actual stress is \(DS(s_{\rm expr}(u,P,\xi))[u]+S(s^*)\). | Agree on symmetry and fixed point; Disagree on linearity | H-D2: the expression depends nonlinearly on u, even with τ frozen. Its exact Newton derivative stays symmetric. At \(s_{\rm expr}=s^*\) the correct energy identity returns. This is not the fully lagged-subscale iteration in chapter §5.9. |
| H20 — local EMAC–VMS update, 355–358 | \((I+\tau_tM_h)s=\tau_t(R_{u,0}+\chi_Ds_{\rm hist}/\delta t-\xi_h)\). | Rearranging \((\chi_D\alpha_0/\delta t+\tau_m^{-1})s+M_hs=R_{u,0}+\chi_Ds_{\rm hist}/\delta t-\xi_h\) gives it exactly. | Agree | It is symmetric. Nonsingularity requires \(1+\tau_t\lambda_i(M_h)\ne0\) for every eigenvalue. The minimum-eigenvalue condition \(1+\tau_t\lambda_{\min}>0\) is sufficient for positive definiteness, not necessary for invertibility. CN has \(\alpha_0=2\). |
| H21 — OSS at every iterate, 358 | Projecting with \(-M_hs^*\) keeps s exactly orthogonal. | It enforces \((s_{\rm expr},w_h)=0\) for the matching converged inner solve. The implicit local update changes that expression while keeping ξ fixed. | Disagree as a claim about every iterate | H-D3 gives the defect and a counterexample. Converged, matched states are orthogonal. Raw Newton iterates need not satisfy a nonlinear projection constraint exactly either. |
| H22 — variable τ and (A2), 374 | Weighted projections remove the constant-τ assumption for energy. | \((s,\xi_h)=0\), and \((\tau_c(d-\zeta_h),\zeta_h)=0\), without pulling either τ through an integral. | Agree | ASGS uses ξ=ζ=0 and requires no constant τ. Variable-τ weighted OSS is an extension of the chapter's unweighted projector, not the identical subscale model away from constant τ (H-D4). |
| H23 — exact midpoint energy identity, 374 | Rate of total endpoint energy plus half-level terms equals power. | The midpoint product rule and the two scale equations give H-D4 exactly, including \(-2(s,L_\nu u)\). | Agree at a converged coupled solve with the stated boundary conditions | Applies to dynamic subscales and fixed spaces/quadrature. The nonlinear coefficient τ does not need a time derivative. The BDF1 start has a separate nonnegative increment term. |
| H24 — closure residual, 374 | Rate + dissipation + N − power is zero for skew forms/EMAC–VMS, but equals N for the other forms. | The same summed equations give rate + dissipation + N − power = 0 for **every** form. | Disagree | Without N, rate + dissipation − power = **−N**. Equivalently power − rate − dissipation = N. H-D4 specifies all terms and their signs. |
| H25 — discrete exactness and momentum in time, 364–376 | Exact identities, up to solver tolerance; momentum balances for every scheme. | Pointwise cancellations and matched projections are exact in the computed inner product; the remaining integration-by-parts identities concern only polynomial FE fields. | Agree with qualifications | On affine Pk meshes, degree 3k suffices for the volume identities here, with sufficiently exact facet flux integration. Identical quadrature alone is not an integration-by-parts theorem on curved/non-polynomial geometry. Time schemes preserve the linear tested balance, including fluxes/reactions and their actual discrete time operator. |

**Spec disagreements and qualifications.**

**H-D1 — exact table entries, boundary corrections, and admissible tests (318, 368–372).** Define, for \(\chi=0\) (linear splitting) or \(\chi=1\) (nonlinear splitting),

\[
a_\chi=u+\chi s,\qquad D_\chi(w)=(d,w)-\chi(s,\nabla w)_h,\qquad
B(v)=\int_{\partial\Omega}(u\cdot n)(u\cdot v),\qquad
H(v)=\tfrac12\int_{\partial\Omega}|u|^2(v\cdot n).
\]

Here v is single-valued, piecewise differentiable, and \(B_N(v)\) restricts B to the outflow. Directly expanding the product derivative gives

\[
((a_\chi\cdot\nabla)u,v)+(u,(a_\chi\cdot\nabla)v)+D_\chi(u\cdot v)
=(u,\nabla(u\cdot v))+(d,u\cdot v)=B(v).
\tag{H-D1a}
\]

Every term involving s cancels at the same quadrature points. Only the last equality uses integration by parts, on the continuous FE polynomial product; internal facet fluxes cancel. This is why variable/non-polynomial τ and stored point values do not invalidate these particular identities on the specified affine meshes. For \(v=u\), (H-D1a) gives \(c_{nc}=-D_\chi(|u|^2)/2+B(u)/2\). With \(\phi=e_k\times x\), write \(w_k=u\cdot\phi\), \(J=(s,e_k\times u)\). Then

\[
(u,(a_\chi\cdot\nabla)\phi)=-\chi J,\qquad
S_\phi=-(s,(a_\chi\cdot\nabla)\phi)=-J,
\]

because \(s\cdot(e_k\times s)=0\) pointwise. Finally,

\[
((\nabla u)a_\chi,e_k)=H(e_k)+\chi(s,\partial_k u),\qquad
((\nabla u)a_\chi,\phi)=H(\phi)+\chi(s,(\phi\cdot\nabla)u).
\tag{H-D1b}
\]

These three identities give the following **interior entries**, before adding the boundary corrections in the next table. They independently reproduce every entry of both chapter tables; set s=0 to recover the Galerkin cases.

| Splitting | Form | Linear momentum | Angular momentum | Total convective energy defect \(\mathcal N\) |
| --- | --- | --- | --- | --- |
| Nonlinear | Non-conservative | \(-D_1(u_k)\) | \(-D_1(w_k)\) | \(-\tfrac12D_1(\lvert u\rvert ^2)\) |
| Nonlinear | Conservative | 0 | 0 | \(+\tfrac12D_1(\lvert u\rvert ^2)\) |
| Nonlinear | Skew 1 | \(-\tfrac12D_1(u_k)\) | \(-\tfrac12D_1(w_k)\) | 0 |
| Nonlinear | Skew 2 | \(-\tfrac12D_1(u_k)\) | \(-\tfrac12D_1(w_k)\) | 0 |
| Nonlinear | EMAC | \((s,\partial_k u)\) | \((s,(\phi\cdot\nabla)u)\) | \((s,(u\cdot\nabla)u)\) |
| Nonlinear | EMAC–VMS | 0 | 0 | 0 |
| Linear | Non-conservative | \(-D_0(u_k)\) | \(-D_0(w_k)-J\) | \(-\tfrac12D_0(\lvert u\rvert ^2)\) |
| Linear | Conservative | 0 | \(-J\) | \(+\tfrac12D_0(\lvert u\rvert ^2)\) |
| Linear | Skew 1 | \(-\tfrac12D_0(u_k)\) | \(-\tfrac12D_0(w_k)-J\) | 0 |
| Linear | Skew 2 | \(-\tfrac12D_0(u_k)\) | \(-\tfrac12D_0(w_k)-J\) | 0 |
| Linear | EMAC | 0 | \(-J\) | \((s,(\nabla u)u)=\tfrac12(s,\nabla\lvert u\rvert ^2)\) |
| Linear | EMAC–VMS | 0 | 0 | 0 |

| Form | Add to the linear-momentum entry | Add to the angular-momentum entry | Add to the energy entry | Additional term from the spec's outflow addition |
| --- | --- | --- | --- | --- |
| Non-conservative | \(B(e_k)\) | \(B(\phi)\) | \(B(u)/2\) | None |
| Conservative | 0 | 0 | \(-B(u)/2\) | \(B_N(v)\) for v = \(e_k,\phi,u\) |
| Skew 1 | \(B(e_k)/2\) | \(B(\phi)/2\) | 0 | \(B_N(v)/2\) for v = \(e_k,\phi,u\) |
| Skew 2 | \(B(e_k)\) | \(B(\phi)\) | \(B(u)/2\) | None |
| EMAC | \(B(e_k)+H(e_k)\) | \(B(\phi)+H(\phi)\) | \(B(u)\) | None |
| EMAC–VMS | \(B(e_k)+H(e_k)\) | \(B(\phi)+H(\phi)\) | \(B(u)\) | None |

Thus the conservative volume momentum entries, skew-1 volume resolved-energy cancellation, and symmetric-subgrid-stress force/torque cancellations are exact even on an open domain. Other replacements by the table entries require the indicated surface terms. In particular, with zero u on the Dirichlet part and only an outflow elsewhere, each of the first four forms has the usual \(B_N(u)/2\) energy flux after the additions. EMAC and EMAC–VMS have \(B(u)\) in this modified-pressure convention. The EMAC energy calculation uses

\[
c_E(a_\chi;u,u)=\chi(s,(u\cdot\nabla)u-(\nabla u)u)+B(u),\qquad
\mathcal N=c_E(a_\chi;u,u)+(s,(\nabla u)u).
\]

For EMAC–VMS the separate identity \(\tilde T:Du=s\cdot N_E\) reduces the full energy defect to \(c_E(u;u,u)=B(u)\), for both splittings.

All these convective boundary terms vanish if u=0 on the entire boundary. On periodic domains the linear-momentum and energy fluxes cancel, but the angular ones need not. An exact polynomial counterexample on the periodic unit square is \(u=(y(1-y),1)\), s=0, \(\phi=(-y,x)\). This is a continuous periodic P2 velocity, with d=0, yet

\[
c_{nc}(u;u,\phi)=-\int_0^1y(1-2y)\,dy=\tfrac16=B(\phi),
\]

whereas the interior entry \(-D_0(w_k)\) is zero. A general “periodic implies no angular flux” assertion is therefore false, including in the actual FE setting.

For zero-mean pressure tests, the nonlinear equal-order momentum inference also uses \(D_1(1)=\int_{\partial\Omega}u\cdot n=0\) to extend the test to \(u_k\) minus its mean; it holds on the closed and periodic domains under discussion. In linear splitting the same argument gives \(D_0(u_k)=(s,\nabla u_k)\), not zero. Likewise \(J=0\) for global OSS because \(e_k\times u\in\hat V_h\). These are global projection statements, not arbitrary-patch ones.

**Required qualification:** replace the blanket periodic-boundary sentence by the explicit flux conditions above. Also distinguish these algebraic defect evaluations from testing the solved momentum equation: \(e_k\) and \(\phi\) are not admissible homogeneous-Dirichlet tests, and \(\phi\) is not generally an admissible periodic test. Full momentum laws on those domains need the corresponding boundary reactions/fluxes (the chapter's \(B_R\)); zero boundary velocity does not imply zero pressure or viscous force on the wall.

**H-D2 — the proposed stress is symmetric but is not linear in u (352).** At fixed stored s*, put

\[
s_{\rm expr}=\tau_t\bigl(R_{u,0}(u,P)+\chi_Ds_{\rm hist}/\delta t-M_h(u)s^*-\xi_h\bigr).
\]

Then the spec's actual stress is

\[
T^*=s_{\rm expr}\otimes u+u\otimes s_{\rm expr}+(u\cdot s_{\rm expr})I+S(s^*)
=DS(s_{\rm expr})[u]+S(s^*).
\]

The cross stress is bilinear in its two arguments, and the first argument here is a nonlinear expression in the solve unknowns. Even with constant τ, \(R_{u,0}\) contains \(2Du\,u+du\), so T* is generically at least cubic in u. The extra dependence of τ on \(u+s^*\) reinforces the failure of linearity. Chapter §5.9 instead lags **every** occurrence of s in this stress, obtaining \(DS(s^*)[u]+S(s^*)\), affine in u with a constant right-hand-side term.

The spec stress and its exact Newton variation are nevertheless symmetric. Consequently their contractions with \(D e_k\) and \(D\phi\) vanish pointwise, which is the actual force/torque argument; linearity is unnecessary for that cancellation. At the common fixed point \(s_{\rm expr}=s^*=s\), all stresses and residual terms recover (5.23)–(5.26), and the energy cancellation holds. With \(\chi_N=0\), (5.29) is recovered directly. For variable τ, this statement refers to the spec's weighted-projection extension of the subscale closure (H-D4).

The energy coefficient check, for either splitting, is explicitly

\[
\begin{aligned}
\tilde T:Du
&=2s\cdot(Du)u+d(u\cdot s)+\chi_N\bigl(s\cdot(Du)s+\tfrac12d|s|^2\bigr)\\
&=s\cdot\bigl(2(Du)u+du+\chi_N(Du+\tfrac12dI)s\bigr)=s\cdot N_E.
\end{aligned}
\]

**Required correction:** remove the claim of linearity, and describe this as Newton differentiation of the substituted, symmetric stress. If the intent is to implement the chapter's linear stress matrix (5.37) alone, freeze s also in the cross terms; that would be a different iteration rule. Do not infer exact nonlinear momentum balances at every raw Newton iterate from symmetry alone: in ASGS the substituted time term \(D_ts_{\rm expr}(u,P,\xi)\) is nonlinear too, so its Newton equation enforces the linearized balance. The full balance follows at convergence of the inner solve, with its evaluated subscale and appropriate boundary terms.

**H-D3 — distinguish the projected expression from the implicit local update (355–358).** At fixed u and τ, set \(X=R_{u,0}+\chi_Ds_{\rm hist}/\delta t\). The inner OSS projection satisfies

\[
(\tau_t\xi_h,w_h)=(\tau_t(X-M_hs^*),w_h),\qquad
s_{\rm expr}=\tau_t(X-M_hs^*-\xi_h),\qquad (s_{\rm expr},w_h)=0.
\]

But step 3 solves \(s_{\rm new}=\tau_t(X-M_hs_{\rm new}-\xi_h)\). Subtracting gives

\[
(s_{\rm new},w_h)=-\bigl(\tau_t M_h(s_{\rm new}-s^*),w_h\bigr),
\tag{H-D3}
\]

which is generally nonzero. A minimal scalar component example uses two equal quadrature weights, projection onto constants, τt=1, s*=0, \(X=(1,-1)\), \(M=(0,1)\). Then ξ=0 and \(s_{\rm expr}=(1,-1)\) is orthogonal, but the local solve gives \(s_{\rm new}=(1,-1/2)\), whose mean is 1/4. Re-evaluating τ during local sweeps introduces another mismatch; it cannot restore the missing identity in general.

The local system itself is the correct rearrangement of the intended subscale equation and has the correct fixed point. The claim that its iterates remain orthogonal is the error. The current algorithm can retain that update, provided it promises orthogonality only for the matched expression after a converged inner solve and at convergence of the outer loop, and stores the final matched expression. Exact orthogonality of every implicit local update would require solving its projection consistently too: for fixed coefficients let \(W=(\tau_t^{-1}I+M_h)^{-1}\), impose \((W\xi_h,w_h)=(WX,w_h)\), and use \(s=W(X-\xi_h)\). In the positive definite regime this is a matrix-weighted projection, rather than the proposed lagged scalar-weighted projection. Equivalently solve the original projection with the current \(M_hs\) coupled. A raw Newton iterate of either nonlinear projection equation satisfies it only to its current residual, not identically.

**H-D4 — energy identity, variable τ, and the diagnostic sign (374).** With ξ=ζ=0 in ASGS, or matched weighted projections in OSS, the discrete pressure-subscale contribution is

\[
Q_c=(\tau_c(d-\zeta_h),d)=(\tau_c(d-\zeta_h),d-\zeta_h)\ge0.
\]

The equality uses the pressure projection tested with ζ; it is not merely a consequence of velocity orthogonality. In OSS, velocity orthogonality gives \((s,\xi_h)=0\), removes \((D_tu,s)\), and, for orthogonal histories, removes \((D_ts,u)\). For ASGS the two cross time terms are retained. These cancellations work for arbitrary positive pointwise τ without differentiating τ or taking it outside an integral. Thus (A2) is unnecessary for the spec. The chapter's original unweighted OSS equation with variable τ would not preserve the range of its projector: it must actually be replaced by these weighted constraints to obtain this extension.

For dynamic subscales, let bars denote the endpoint averages, and define

\[
\begin{aligned}
E^n_{\rm ASGS}&=\tfrac12\|u^n+s^n\|_Q^2,\\
E^n_{\rm OSS}&=\tfrac12(\|u^n\|_Q^2+\|s^n\|_Q^2),\\
\mathcal D(\bar u,\bar s)&=a_\nu(\bar u,\bar u)+(\tau_m^{-1}\bar s,\bar s)_Q
-2(\bar s,L_\nu\bar u)_Q+(\tau_c(\nabla\cdot\bar u-\bar\zeta_h),\nabla\cdot\bar u)_Q.
\end{aligned}
\]

Testing the resolved equation with \((\bar u,\bar\pi_h)\), the subscale equation with \(\bar s\), and adding gives, on the zero-boundary/periodic energy domain,

\[
\boxed{\frac{E^{n+1}-E^n}{\delta t}+\mathcal D(\bar u,\bar s)+\mathcal N(\bar u,\bar s)
=(\bar f,\bar u+\bar s)_Q.}
\tag{H-D4}
\]

For the standard coupling, \(\mathcal N=c(a;\bar u,\bar u)+(\bar s,N_{\rm op}-(a\cdot\nabla)\bar u)_Q\). For EMAC–VMS it is \(G(\bar u)+(\bar s,N_E)_Q=0\) under these energy boundary conditions. The viscous cross term has coefficient **−2**: one copy comes from each scale. This signed cross term belongs in the balance even though its positivity is only a property of the combined dissipation under an inverse-estimate condition. The midpoint time identity is simply

\[
(\delta_t(u+s),\overline{u+s})_Q=(\|u^{n+1}+s^{n+1}\|_Q^2-\|u^n+s^n\|_Q^2)/(2\delta t).
\]

For OSS, orthogonal endpoints make this equal to the stated sum of energies. Endpoint extrapolation preserves orthogonality by linearity, starting from the projected BDF1 result. The half-level projection, residual, stored s at convergence, and τ must all describe the same state. Static subscales have no such two-scale dynamic midpoint identity. Nonzero Dirichlet data require boundary work/lifting terms; on outflow domains the fluxes in H-D1 and the specified boundary work must be retained.

**Required correction:** define the closure residual as the entire left side of (H-D4) minus its right side, and require zero for **all** forms. To measure the predicted defect, use power − rate − \(\mathcal D=\mathcal N\), or label rate + \(\mathcal D\) − power as \(-\mathcal N\). The BDF1 startup instead adds \(\|\Delta(u+s)\|_Q^2/(2\delta t)\) in ASGS, and \((\|\Delta u\|_Q^2+\|\Delta s\|_Q^2)/(2\delta t)\) in OSS; excluding it from the midpoint-only check is correct.

**H-D5 — natural conditions and the scope of the bubble assumption (310, 322).** The first four outflow additions and the resolved EMAC modified-pressure condition pass under the spec's zero-subscale-boundary convention. “EMAC–VMS integrates nothing by parts” is not an adequate explanation: \(-(\tilde T,Dv)\) is already a weak stress divergence. If a sufficiently regular s is allowed a nonzero trace, then

\[
-(\tilde T,Dv)=(\nabla\cdot\tilde T,v)-\langle\tilde T n,v\rangle,
\]

so this term changes the resolved natural traction to \((2\nu Du-PI-\tilde T)n=0\), before any further non-bubble subscale boundary contributions are considered. Under the prescribed bubble/no-subscale-flux model \(\tilde T n=0\), giving the condition printed in the spec. For a mere L² quadrature field the normal trace is not defined, and the omitted flux must be understood as part of that model, not proved from “volume term.” The stress-form sentence at 322 uses physical p for the first four forms; EMAC variants use P in the imposed condition, while reported physical stress still uses the post-processed p.

**Chapter issues, separate from spec issues.** No sign or coefficient error was found in (5.1)–(5.29), the entries of Tables 5.1–5.2 under their flux convention, or the printed element matrices. The following statements or implementation qualifications do need correction.

| Item | Chapter location | Finding and correction |
| --- | --- | --- |
| C01 — two different EMAC energy defects | §5.1, fourth bullet, p. 2 | The introduction identifies “the energy defect” as the Lamb-vector work \((s,\omega_h\times u)\). That is the **resolved convective production** of Proposition 5.3. With the consistent residual, the total-energy defect is \((s,(u\cdot\nabla)u)\), as correctly stated in Corollary 5.1 and Table 5.1. Label the introductory assertion accordingly. |
| C02 — a volume stress can change natural traction | §5.9.3, final sentence, p. 15 | “The subgrid stress, being a volume term, adds none” is false without a zero normal stress trace or an explicitly omitted flux. H-D5 gives the missing \(-\tilde Tn\). The chapter advertises EMAC–VMS without (A1); that is valid for its volume conservation algebra, but cannot by itself imply the unchanged outflow condition. |
| C03 — skew form and outflow | §5.9.1, paragraph after the element matrices, p. 14 | Using skew 1 as written does not add an explicit facet integral, but its natural condition differs from convective-form do-nothing. Without the spec's added half flux, the gradient-form condition is \(\nu\partial_nu-\pi n-\tfrac12(u\cdot n)u=0\). Conservative form similarly has the full flux. State this when using those matrices on outflow problems. This is a boundary qualification, not a wrong volume matrix. |
| C04 — the symmetric Picard recipe is a different iteration | (5.32), §5.7.6 and §5.9.3 | With b the old resolved velocity, v the new one and fixed s, (5.32) uses \(\tfrac12DS(b+s)[v+s]\). The §5.9.3 recipe instead uses \(\tfrac12DS(b)[v]+DS(s)[v]+S(s)\). Their difference is \(\tfrac12DS(s)[v-b]\), generally nonzero. Both are symmetric and have the same fixed point; distinguish them rather than presenting the assembly recipe as (5.32). The analogous Newton recipe does reproduce (5.31), up to the explicitly suppressed fluxes. |
| C05 — local OSS orthogonality is not inherited from global OSS | §5.3 after (5.16); §5.4 after (5.17); Proposition 5.6 | Orthogonality to global continuous FE functions does not imply \(\int_R s=0\) or \(\int_R x\times s=0\) for every element patch. Restriction to R would require testing against \(1_R e_k\) or \(1_R\phi\), which generally are not global FE functions. These deductions are valid on the whole domain when the tests belong to the projection space, or for an explicitly separate patchwise projector, not for the restriction of the globally projected solution. The global OSS space on a periodic domain generally excludes φ as well. The convective moment identities up to boundary flux remain valid independently of this issue. |
| C06 — BDF is not universally extra dissipation of the physical energy | Remark 5.8, p. 12 | The midpoint assertion is correct; “Backward differentiation formulas add numerical dissipation” is too broad. BDF1 has the positive physical-energy increment term in H-D4. BDF2 has a positive remainder with a **two-level modified energy**, rather than the physical energy alone. BDF3 has no general nonnegative dissipation identity of this kind. Restrict the claim to the actual scheme and energy. |
| C07 — residual equivalence is consistency, not equality at a discrete state | Remark 5.10 versus Remark 5.7(2), pp. 10, 12 | In physical pressure, \(R_u=f-\partial_tu-(u\cdot\nabla)u+L_\nu u-\nabla p_h-du-M_hs\). Thus calling (5.26) simply “the physical momentum residual” omits \(-du-M_hs\). Remark 5.7(2) already makes the correct distinction. Also \(\tau_m^{-1}\Vert s\Vert ^2\approx\tau_m\Vert \widetilde P R_u\Vert ^2\) is a quasi-static scale estimate, not a bound furnished by the dynamic energy identity without control of \(\partial_ts\). |
| C08 — lagged OSS implementation needs a matched-projection qualification | Remark 5.9 and §5.9.3 | Lagging the projected residual or \(M_hs\) retains the intended coupled fixed point, but does not guarantee exact orthogonality of each subsequent implicit local update. H-D3 gives the obstruction even at constant τ. The chapter does not explicitly promise orthogonality at each local sweep, but implementations must not invoke the converged OSS cancellations during those unmatched sweeps. |

For C06, the exact BDF2 identity with \(a=g^{n+1},b=g^n,c=g^{n-1}\) is

\[
(\tfrac32a-2b+\tfrac12c,a)
=\tfrac14\bigl(\|a\|^2+\|2a-b\|^2-\|b\|^2-\|2b-c\|^2+\|a-2b+c\|^2\bigr).
\]

This identifies the modified energy precisely; it is not \((\|a\|^2-\|b\|^2)/2\) plus a universally nonnegative remainder.

For completeness, the implementation section was checked directly, without trusting its transcription. In the following table each entry is the integrand, \(A=N^A\), \(B=N^B\), and b is the old resolved velocity. \(K_{ij}^{AB}=c(a;Be_j,Ae_i)\), \(D_{ij}^{AB}=c(Be_j;b,Ae_i)\).

| Form | K integrand re-derived | D integrand re-derived | Printed §5.9 |
| --- | --- | --- | --- |
| Non-conservative | \(\delta_{ij}A a_l\partial_lB\) | \(AB\partial_jb_i\) | Both agree |
| Conservative | \(-\delta_{ij}B a_l\partial_lA\) | \(-B b_i\partial_jA\) | Both agree |
| Skew 1 | \(\tfrac12\delta_{ij}(A a_l\partial_lB-B a_l\partial_lA)\) | \(\tfrac12(AB\partial_jb_i-Bb_i\partial_jA)\) | Both agree |
| Skew 2 | \(\delta_{ij}A(a_l\partial_lB+\tfrac12(\nabla\cdot a)B)\) | \(A(B\partial_jb_i+\tfrac12b_i\partial_jB)\) | Both agree |
| EMAC | \(A\partial_iB\,a_j+\delta_{ij}A(a_l\partial_lB+(\nabla\cdot a)B)\) | \(A(B\partial_i b_j+B\partial_j b_i+b_i\partial_jB)\) | Both agree |

The additional subscale-divergence contribution is \(-\delta_{ij}(s,\nabla(AB))\), with factor 1/2 for skew 2, as printed. Expanding \(-DS(W)[Be_j]:D(Ae_i)\) gives

\[
-B\bigl(W_l\partial_lA\,\delta_{ij}+W_i\partial_jA+W_j\partial_iA\bigr),
\]

which verifies (5.37). Moving the constant stress to the right gives \(s_i s_l\partial_lA+|s|^2\partial_iA/2\), also with the printed sign. The Newton RHS in (5.36) is positive \(c(b;b,v)\); the equation itself has the correct negative subtraction on the left. The local ASGS backward-Euler matrix in §5.9.3 agrees with H20. None of these verified matrices supplies the additional derivative of \(s_{\rm expr}\) required by the spec's different substituted-subscale Newton iteration.

Remaining chapter improvements are qualifications rather than errors in the verified equations: define the quadrature divergence pairing explicitly if exact discrete Lemma 5.2 is intended; state the weighted-projection extension before dropping (A2); and retain fluxes in the energy assertions of Remark 5.12 outside the homogeneous/periodic energy setting. Sufficient inverse bounds are \(c_1k^4>2C_{\rm inv}^2\) (Laplacian) and \(c_1k^4>4C_{\rm inv}^2\) (stress), with the chapter's **unsquared** inverse-estimate constant. This differs in notation from the spec's earlier constant defined directly in a squared norm inequality.

---

**Check of the Phase 0 answers — 21 October 2026**

Reviewed `docs/formulation_spec.md`, SHA-256 `bce98772729e2c9e799f3317fffb71716a55ae86ca53b4d2815612b9ea5e390f`, against the addendum pass below. Every change from the previous revision is wording; no equation changed. I checked the Time schemes section, the Fixed choices time-scheme row, the eight-cases table and the algorithm, plus every remaining history-level reference.

**Result: two items do not match yet; the rest matches.**

| Item | Spec says | Phase 0 finding or answer | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| B01 — Crank–Nicolson BDF1 first step, lines 24, 89 | Step 1 is BDF1 for every scheme in both startup modes; this removes the start-up defect. | Answer 1; A06, A08, A-D3. | Agree | With \(u_h^1\) from BDF1, the endpoint Dirichlet match (line 87) and the orthogonality of \(\tilde u^1\) hold. Second order is kept. |
| B02 — exact-history levels, line 89 | Velocity histories at \(t^0,t^{-1},t^{-2}\); BDF2 and BDF3 at full order from step 1. | Answer 2; A10, A-D4. | Agree | |
| B03 — stability wording, line 91 | BDF3 is A(α)-stable with α ≈ 86° and a bounded excluded region; the 2007 remark concerns eq. 49 for Theorem 3. | A11, A-D5. | Agree | |
| B04 — endpoint-average meaning, line 87 | \(g^{n+1/2}=\tfrac12(g^n+g^{n+1})\), not \(g(t^{n+1/2})\). | A04, A-D1. | Agree | |
| B05 — stale references, lines 5, 100, 162–168, 239 | "solve level"; histories back to n−2 for BDF3. | A14. | Agree | |
| B06 — stale reference, line 117 | "ũ, ũⁿ, ũⁿ⁻¹ are stored at its points." | A14 listed line 117; answer 3 says the n−1 references are fixed. | Disagree | BDF3 also stores \(\tilde u^{n-2}\). Suggested: "ũ and the subscale histories the scheme needs are stored at its points". |
| B07 — subscale accuracy sentence, line 89 | "Subscale histories start at zero in both modes, so in general problems the first subscale value is only first-order accurate." | A10, A-D4: the defect comes from zero values at \(t^{-1},t^{-2}\) under exact-history BDF2/BDF3. A BDF1 first step from \(\tilde u^0=0\) is consistent. | Disagree (wording) | In lower-order mode, and for Crank–Nicolson in both modes, step 1 is BDF1 from the true initial value \(\tilde u^0=0\), so there is no defect. Suggested: "With exact-history startup, BDF2 and BDF3 use zero subscale values at negative times, which are not a history of the subscale equation, so in general problems the first subscale value has an error of order δt·\|R\_u\|; a BDF1 first step from ũ⁰ = 0 has no such defect. The temporal convergence test is not affected." |

Two minor points do not block anything. Line 24 says BDF2 and BDF3 "start with lower orders" without mentioning the exact-history option. Line 5 still says "Every variant ... with three switches", although the Galerkin baseline is not a switch setting (A14). B06 and B07 are wording only; they change no equation and nothing in Phase 1.

**Follow-up:** in the spec revision with SHA-256 `5629e68cbd26a3984023eff16f09f70a50f508e0557bedf511f128dca324b339`, B06, B07, the line-24 time-scheme row and the line-5 Galerkin-baseline sentence all match the suggestions above, and no other text changed, so every item now agrees.

---

**Addendum pass — 21 October 2026**

This pass covers the spec material added for the implementation plan (plan §13, Phase 0): the new sections Time schemes and Galerkin baseline, the revised notation, and the revised fixed-choice rows Variants, Interpolation, Time scheme and Quadrature. The rest of the spec is unchanged from the revision checked in the final pass below, so P01–P08 remain the disposition for it.

**Verdict: the addendum's equations pass; no equation needs to change.** BDF3 differentiates every cubic exactly. The Crank–Nicolson row is exactly Codina et al. 2007 eqs. (22)–(23) at \(\theta=\tfrac12\), including their definition of \(f^{n+\theta}\), and it gives their effective parameter (26). I re-derived τ\_t for Crank–Nicolson, the orthogonality of the extrapolated OSS subscale, the second-order limit of BDF3 with lower-order startup, the Dirichlet averaging and the Galerkin baseline, and all of them agree with the spec.

**Disagreements:** the stability paragraph contains two wording errors. The unstable set of BDF3 is bounded, not a sector. The 2007 remark is about the increment bound (49), not about stability in general. The Crank–Nicolson text also omits a condition: the extrapolation never damps a defect in the initial state. So u\_h⁰ must satisfy the Dirichlet data, and for exact discrete incompressibility it must also be discretely divergence-free; otherwise Crank–Nicolson needs a startup step. Several references were not updated to the solve level or to BDF3's third history level. It is also unclear which time levels exact-history startup fills. Details follow the table.

Reviewed `docs/formulation_spec.md`, SHA-256 `93903adfa709d3be485d9cc17cb2300ed4b9f5f00c2d63158b777d742111ce98`, and line numbers refer to this revision. The spec copy in the thesis folder (SHA-256 `fc581b2ed68377b15ee0a04e260f3acaf56da8df02eb0055c79fa6b75fb53956`) differs from it only by the addendum; I checked the diff line by line. That copy already contains the three corrections the final pass asked for (A15). Source equations come from non-layout `pdftotext` output. Pages 2417, 2421 and 2422 of `time_dependent_2007.pdf` were read as rendered images because the extraction garbles them. Two numerical cross-checks ran as scratch calculations outside the repository: exact rational sums for the BDF weights and the BDF3 boundary locus. No project code was written. The spec was not edited.

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| A01 — one time-scheme template, lines 54–58, 76, 85 | Every scheme is \(D_tg=(\alpha_0g-g_{\rm hist})/\delta t\) at a solve level, so the master formulation, τ\_t and the projections need no change. | Independent derivation A-D1; the earlier D2 shows the master system uses the scheme only through \(\alpha_0\) and \(g_{\rm hist}\). | Agree | In every row the history weights sum to \(\alpha_0\): \(1\); \(2-\tfrac12=\tfrac32\); \(3-\tfrac32+\tfrac13=\tfrac{11}6\); \(2\). So \(D_t\) of a constant state is zero, and the static/dynamic steady-state equivalence (line 153) also holds for BDF3 and Crank–Nicolson. |
| A02 — BDF1 and BDF2 rows, lines 80–81 | Unchanged values: \(\alpha_0=1,\ g^n\) and \(\alpha_0=\tfrac32,\ 2g^n-\tfrac12g^{n-1}\). | `time_dependent_2007.pdf`, (22) with \(\theta=1\) is backward Euler; `colomes_VMS_paper.pdf`, §4 opening, writes its fully discrete problem with backward Euler. BDF2 for both scales is ours. | Agree | Re-checked in A-D1: BDF2 is exact for quadratics and has error \(-\tfrac13\delta t^2g'''\). The attributions are correct. |
| A03 — BDF3 row, line 82 | \(\alpha_0=\tfrac{11}6\), \(g_{\rm hist}=3g^n-\tfrac32g^{n-1}+\tfrac13g^{n-2}\), order 3, ours. | No source uses BDF3; this is checked by derivation A-D1. | Agree | Exact for every polynomial of degree at most 3 at \(t^{n+1}\). Leading error \(-\tfrac14\delta t^3g''''\). |
| A04 — Crank–Nicolson one-leg form and data, lines 83, 87 | Half-level unknowns \(u_h^{n+1/2},p_h^{n+1/2},\tilde u^{n+1/2}\); \(\alpha_0=2\), \(g_{\rm hist}=2g^n\); every spatial term, τ, a and the projections at the half level; \(f^{n+1/2}=\tfrac12(f^{n+1}+f^n)\); u\_h and ũ extrapolated, pressure reported at \(t^{n+1/2}\). | `time_dependent_2007.pdf`, definitions before (22): \(\delta_tf^n=(f^{n+1}-f^n)/\delta t\), \(f^{n+\theta}=\theta f^{n+1}+(1-\theta)f^n\), \(\tfrac12\le\theta\le1\). (22) evaluates every spatial term, including \(u_h^{n+\theta}\cdot\nabla u_h^{n+\theta}\) and the subscale terms, on \(\theta\)-states; its load is \(f^{n+\theta}\) (a time average if f is not continuous in time). (23): \(\delta_t\tilde u^n+\tilde u^{n+\theta}/\tau^{n+\theta}=r_{u,h}^{n+\theta}\). | Agree | A-D1 gives the exact rewriting. The identity at line 87 holds because \(g^{n+1/2}\) means the paper's endpoint average, not \(g(t^{n+1/2})\); the two differ by \(\tfrac18\delta t^2g''\). The paper labels the pressure \(p_h^{n+1}\). A-D2 shows that this unknown approximates \(p(t^{n+1/2})\) to second order but \(p(t^{n+1})\) only to first order, so reporting it at \(t^{n+1/2}\) is correct; `Kabalan_Ali_Thesis.pdf`, (3.30), writes \(\nabla p^{n+\theta}\). The paper never defines \(\tau^{n+\theta}\) beyond its superscript. Evaluating \(\tau(a^{n+1/2})\) is the one-leg reading consistent with (22); averaging the endpoint values of τ would differ by \(O(\delta t^2)\). `colomes_VMS_paper.pdf`, §4 opening, gives further support: its examples in §§5–7 use Crank–Nicolson. |
| A05 — τ\_t for Crank–Nicolson, line 85 | \(\tau_t=(2/\delta t+1/\tau_m)^{-1}\), with \(\tilde u_{\rm hist}=2\tilde u^n\). | `time_dependent_2007.pdf`, (25)–(26): \(\tau_t=(1/(\theta\delta t)+1/\tau^{n+\theta})^{-1}\) and \(\tilde u^{n+\theta}=\tau_t(r^{n+\theta}_{u,h}+\tilde u^n/(\theta\delta t))\). | Agree | At \(\theta=\tfrac12\) these are identical to the spec's master subscale with \(\alpha_0=2\) and \(\tilde u_{\rm hist}=2\tilde u^n\) (A-D1). The static case gives \(\tau_t=\tau_m\). For OSS, \(D_tu_h=(2u_h^{n+1/2}-2u_h^n)/\delta t\) is still an FE function, so the reduced residual is exact (earlier D3b). |
| A06 — Dirichlet averaging, line 87 | On Γ\_D the half-level velocity is \(\tfrac12(g^{n+1}+g^n)\), "so the extrapolated \(u_h^{n+1}\) matches the Dirichlet data exactly." | The 2007 paper uses homogeneous \(V_{0,h}\), so this is the paper's \(f^{n+\theta}\) applied to boundary data (ours). Derivation A-D3. | Agree — condition missing | This holds by induction **only if** \(u_h^0=g_h^0\) on Γ\_D, with the same linear boundary interpolation at every level. Otherwise the boundary defect changes sign every step and never decays (\(e^{n+1}=-e^n\)). Through the mass term it then puts an alternating load of size \(2\lvert e^0\rvert/\delta t\) on the interior equations next to Γ\_D. Using \(g(t^{n+1/2})\) instead would break the exact endpoint match, so the averaging itself is the right choice. |
| A07 — orthogonality of the extrapolated OSS subscale, line 87 | Extrapolating two orthogonal subscales keeps \(\tilde u^{n+1}\) orthogonal. | `colomes_VMS_paper.pdf`, (53), weight and history inside the projection; earlier D3; derivation A-D3. | Agree | The extrapolation is linear and acts pointwise at the common quadrature points. Starting from \(\tilde u^0=0\), \((\tilde u^n,w_h)_Q=0\) for all n, so \((D_t\tilde u,v_h)=0\) and dropping that term for OSS (\(\chi_A=0\)) stays valid with Crank–Nicolson. One qualification: with a finite outer tolerance, the defect follows \(e^{n+1}=2\varepsilon^{n+1/2}-e^n\). It is carried forward, whereas BDF schemes evaluate \(\tilde u^{n+1}\) directly and reset it every step. Tight outer tolerances are enough; no formula changes. |
| A08 — other constraints under extrapolation, lines 87, 89 | "Crank–Nicolson needs no startup." | Derivation A-D3. `time_dependent_2007.pdf`, §4, analyses backward Euler only. | Qualified — missing condition (open question) | The sign flip of A06 applies to every constraint that the half-level equations impose. For Galerkin, \((q_h,\nabla\cdot u_h^{n+1})=-(q_h,\nabla\cdot u_h^n)\). If \(u_h^0\) is not discretely divergence-free, the endpoint velocities keep a gradient component of fixed size and alternating sign. The half-level pressure then carries \(\pm2\varphi^0/\delta t\), which grows as \(\delta t\to0\). The half-level velocity is unaffected. Common causes are interpolating a divergence-free field into Taylor–Hood and an impulsive start. A BDF step imposes the constraint at \(t^{n+1}\) and clears the defect after the startup steps. V2 is not affected: its \(u^0\) lies in every FE space and is exactly divergence-free. For VMS the continuity equation contains ũ, so I have not quantified the effect there. |
| A09 — BDF3 order with lower-order startup, line 89 | BDF1 then BDF2 keeps BDF2 second order but limits BDF3 to second order globally. | No source; standard linear multistep convergence theory; derivation A-D4. | Agree | The BDF1 step has local error \(O(\delta t^2)\). A zero-stable k-step method of order p has global error \(O(\delta t^{p})+O(\max_{j<k}\lvert e^j\rvert)\), and the starting errors are \(O(\delta t^2)\). That gives order 2 for BDF2 and generically 2 for BDF3. The plan's V2 threshold of at least 1.8 for BDF3 with `lower_order` is consistent with this. |
| A10 — exact-history startup, line 89 | Velocity histories may be set from the exact solution; subscale histories start at zero; V2 uses this to observe third order. | No source; derivation A-D4. | Agree for V2 — qualified in general; ambiguous levels | In V2 the exact velocity lies in the FE space, so the subscale is driven only by the temporal truncation error, \(O(\delta t^3)\). Zero subscale histories then cost nothing. In general, zero values at \(t^{-1},t^{-2}\) are not a history of the subscale ODE started from \(\tilde u^0=0\): BDF3 gives \(\tilde u^1\approx\tfrac6{11}\delta t\,R_u\) instead of \(\approx\delta t\,R_u\), a first-order starting error. The spec claims third order only for V2, so this is not a disagreement. The spec does not say whether "histories" means \(t^{-1},t^{-2}\) or the first computed levels; this is open question 2 below. |
| A11 — stability notes, line 91 | BDF1 and BDF2 are A-stable. BDF3 is not: "its stability region misses a sector near the imaginary axis." Crank–Nicolson is A-stable without high-frequency damping, and "the unconditional stability result of Codina et al. 2007 does not extend to it (their remark after the proof of Theorem 2)." | `time_dependent_2007.pdf`, §4 first paragraph: backward Euler only (\(\theta=1\)). Theorem 1's estimate contains \(\sum\lVert\delta u_h^n\rVert^2+\sum\lVert\delta\tilde u^n\rVert^2\). The text after Remark 9 states the bound (49) on those sums, "a by-product of the proof of Theorem 1 (unfortunately, this property is not to be expected in other schemes such as the Crank–Nicolson time integration)", and uses it for Theorem 3. Derivation A-D5. | Agree on the A-stability facts; Disagree with two wordings | (a) BDF3 is A(α)-stable with \(\alpha\approx86.03^\circ\). Its unstable set in the left half-plane is **bounded**: it lies in \(-\tfrac1{12}<\operatorname{Re}z<0\), \(\lvert\operatorname{Im}z\rvert<1.94\), and includes that segment of the imaginary axis. That is not a sector. The spec's practical conclusion (weakly damped convection may need a smaller δt) is right. (b) The remark concerns the increment bound (49), not stability in general. The paper proves nothing for Crank–Nicolson either way. A-D5 shows why (49) is lost: the \(\theta=\tfrac12\) energy identity has no increment term. Replacement wording is given below. |
| A12 — Galerkin baseline, lines 172–183 | The master FE equation with every stabilization term removed and \(a=u_h\). Taylor–Hood P(k+1)/Pk, no equal order. Same convective form, time schemes, boundary conditions, and Newton/Picard (Picard freezes \(a=\hat u_h\)). No ũ, projections or τ; quadrature degree \(3(k+1)\). | No source states this baseline; it is ours. Checked against the master formulation (lines 126–129) and the Picard system (lines 224–228) in derivation A-D6. | Agree | Removing the three stabilization terms is the same as \(\tilde u\equiv0,\ \tau_c\equiv0\). Then \(a=u_h\) under either splitting, and spec line 177 follows term by term. Picard and Newton reduce consistently, and the two Newton τ options coincide because the baseline has no τ. On affine elements with velocity degree \(m=k+1\), every Galerkin integrand has degree at most \(3m\), so the rule is exact for these forms. Taylor–Hood inf-sup stability is standard background with mild mesh conditions. No supplied source treats it, so it is accepted here, not verified. |
| A13 — revised fixed-choice rows, lines 13, 15, 24, 25 | Variants add the Galerkin baseline (ours). Interpolation: Taylor–Hood for Galerkin (ours). Time scheme: four schemes, same for both scales, lower-order starts, with the attributions shown. Quadrature: \(q=3k\), k the velocity degree. | `time_dependent_2007.pdf`, §3.1, text after (23): same scheme for both scales is the "natural" choice, and different integrators are also allowed; (22)–(23) give the θ-method. `colomes_VMS_paper.pdf`, §4 opening. `term_by_term.pdf`, (8)–(9), BDF2 for the FE velocity only. `onseen_2008.pdf`, §1, allows unequal orders. | Agree | Every attribution checks out, and choices beyond the sources are marked as ours. The Galerkin degree \(3(k+1)\) matches "k the velocity degree" and plan Phase 3. The symbol k is the velocity degree in the Quadrature and VMS rows but the pressure degree in "P(k+1)/Pk"; the explicit \(3(k+1)\) removes the ambiguity. |
| A14 — revised notation, lines 36, 54, 121, 87 | Everything is at the solve level (\(t^{n+1}\), or \(t^{n+1/2}\) for Crank–Nicolson). For Crank–Nicolson, read n+1 as n+½ in the master formulation and the algorithm. All formulas hold for d = 2 and 3. | Derived; consistency read of the whole spec. | Agree — stale references remain | No formula is dimension-specific: \(L_a\), \(L_a^*\), τ with cell diameter h, c, the d×d update at line 251, and the projections. Several places were not updated. Line 5 still says "at time \(t^{n+1}\)" and "Every variant ... with three switches", although the baseline is not a switch setting. Line 100 still says "at \(t^{n+1}\)"; the paragraph at line 87 gives the correct half-level value. BDF3 needs \(\tilde u^{n-2}\) and \(u_h^{n-2}\), but lines 162–168, 117 and 239 stop at \(n-1\). Plan Phase 6 already says three levels. |
| A15 — corrections carried over from the final pass | Line 117: force and stored histories of degree at most 2k−1. Line 219: "does not guarantee quadratic convergence". Line 233: "the reduced minus the full subscale", with matched inputs. | Final pass below: P03 replacement and the two wording refinements. | Agree — resolved | All three appear as recommended. |

**A-D1. Time-level identities.**

A formula \(Dg^{n+1}=\delta t^{-1}\sum_{j\ge0}\beta_jg^{n+1-j}\) is exact at \(t^{n+1}\) for every polynomial of degree at most m if and only if

\[
\sum_j\beta_j(-j)^p=\delta_{p1},\qquad p=0,\dots,m.
\]

For BDF3, \(\beta=(\tfrac{11}6,-3,\tfrac32,-\tfrac13)\). The sums for \(p=0,1,2,3\) are

\[
\tfrac{11}6-3+\tfrac32-\tfrac13=0,\quad 3-3+1=1,\quad -3+6-3=0,\quad 3-12+9=0,
\]

and for \(p=4\) the sum is \(-3+24-27=-6\). Hence BDF3 is exact for cubics and \(Dg-g'=-\tfrac6{4!}\delta t^3g''''+O(\delta t^4)=-\tfrac14\delta t^3g''''+O(\delta t^4)\). For BDF2, \(\beta=(\tfrac32,-2,\tfrac12)\) gives \(0,1,0\) and then \(-2\), so the error is \(-\tfrac13\delta t^2g'''\). The \(p=0\) sum says \(\alpha_0\) equals the sum of the history weights, so \(D_t\) of a constant state is zero. At a steady state \(H_s=(\alpha_0/\delta t)\tilde u\), and the dynamic relation reduces to the static one exactly as in the earlier D3.

For Crank–Nicolson, the 2007 paper defines \(g^{n+\theta}=\theta g^{n+1}+(1-\theta)g^n\) and \(\delta_tg^n=(g^{n+1}-g^n)/\delta t\) for every time-dependent quantity. With \(\theta=\tfrac12\) and \(w:=g^{n+1/2}\), we have \(g^{n+1}=2w-g^n\), so

\[
\delta_tg^n=\frac{2w-2g^n}{\delta t}:\qquad \alpha_0=2,\quad g_{\rm hist}=2g^n,\quad g^{n+1}=2w-g^n.
\tag{A-D1a}
\]

This uses w as the endpoint average. Then \(\delta_tg^n=g'(t^{n+1/2})+\tfrac1{24}\delta t^2g'''\), and \(w=g(t^{n+1/2})+\tfrac18\delta t^2g''\). If w were instead read as \(g(t^{n+1/2})\), (A-D1a) would become a one-sided difference over \(\delta t/2\), exact only for linear g, and the extrapolation would not reproduce \(g^{n+1}\).

Substituting (A-D1a) in 2007 (23) with \(\theta=\tfrac12\):

\[
\frac{2\tilde u^{n+1/2}-2\tilde u^n}{\delta t}+\frac{\tilde u^{n+1/2}}{\tau^{n+1/2}}=r_{u,h}^{n+1/2}
\;\Longrightarrow\;
\tilde u^{n+1/2}=\Bigl(\frac2{\delta t}+\frac1{\tau^{n+1/2}}\Bigr)^{-1}\Bigl(r_{u,h}^{n+1/2}+\frac{2\tilde u^n}{\delta t}\Bigr).
\tag{A-D1b}
\]

This is 2007 (25)–(26) at \(\theta=\tfrac12\). It is also the spec's \(\tilde u=\tau_t(R_u+\chi_D\tilde u_{\rm hist}/\delta t-\xi_h)\) with \(\alpha_0=2\), \(\tilde u_{\rm hist}=2\tilde u^n\) and \(\xi_h=0\). Both the ASGS coupling \((\delta_t\tilde u^n,v_h)\) of 2007 (22) and the spec's \((D_t\tilde u,v_h)\) equal \(((2\tilde u^{n+1/2}-2\tilde u^n)/\delta t,v_h)\). The spec's half-level system is therefore 2007 (22)–(23) at \(\theta=\tfrac12\), with the spec's own convective form and subscale closure in place of the paper's. Those spatial differences were verified in the earlier passes and do not depend on the time level.

**A-D2. Which time the Crank–Nicolson pressure approximates.**

For the Stokes part, take \((u^{n+1}-u^n)/\delta t+\nabla P=\bar f:=\tfrac12(f^{n+1}+f^n)\) with \(\nabla\cdot u^{n+1/2}=0\) and \(\nabla\cdot u^n=0\). Then \(\nabla\cdot u^{n+1}=0\), and \(\nabla P\) is the gradient part of \(\bar f\):

\[
\nabla P=\tfrac12\bigl(\nabla p(t^{n+1})+\nabla p(t^n)\bigr)=\nabla p(t^{n+1/2})+O(\delta t^2),\qquad \nabla p(t^{n+1})-\nabla P=O(\delta t).
\]

For the full equations, each of \(\delta_tu^n\), \(u^{n+1/2}\) and \(\bar f\) is within \(O(\delta t^2)\) of its value at \(t^{n+1/2}\), so consistency at \(t^{n+1/2}\) gives the same conclusion. The step's pressure unknown is a second-order approximation of \(p(t^{n+1/2})\) and only a first-order one of \(p(t^{n+1})\). The superscript n+1 in 2007 (22) only names the step's Lagrange multiplier. The spec's choice not to extrapolate the pressure and to report it at \(t^{n+1/2}\) is the accurate one.

**A-D3. What the Crank–Nicolson extrapolation preserves.**

Let \(\ell\) be a linear functional of the state that the half-level equations set to a target \(c^{n+1/2}\), where the target depends linearly on the data and \(c^{n+1/2}=\tfrac12(c^{n+1}+c^n)\). Define the endpoint defect \(e^n:=\ell(u^n)-c^n\). From \(u^{n+1}=2u^{n+1/2}-u^n\),

\[
e^{n+1}=2\bigl(\ell(u^{n+1/2})-c^{n+1/2}\bigr)-e^n .
\tag{A-D3a}
\]

If the half-level equations hold exactly, \(e^{n+1}=-e^n\) and \(\lvert e^n\rvert=\lvert e^0\rvert\) for every n. The extrapolation preserves the constraint if and only if the initial state satisfies it, and it never damps a violation.

- *Dirichlet values* (\(\ell\) = boundary DOFs, \(c^n=g_h^n\)). Exact endpoint data follow if and only if \(u_h^0=g_h^0\) on Γ\_D. Otherwise the boundary DOFs of \(u_h^n\) carry \((-1)^ne^0\). Interior test functions see \(D_tu_h=(2u_h^{n+1/2}-2u_h^n)/\delta t\) through the mass matrix, so the equations next to Γ\_D get an alternating load of size \(2\lvert e^0\rvert/\delta t\). A BDF step imposes \(u_h^{n+1}=g_h^{n+1}\) directly, so \(e^0\) enters only the first k steps.
- *OSS orthogonality* (\(\ell(\tilde u)=(\tilde u,w_h)_Q\), \(c=0\)). Starting from \(\tilde u^0=0\), \((\tilde u^n,w_h)_Q=0\) for all n, because the extrapolation is pointwise at the common quadrature points. With a finite outer tolerance the half-level defect is \(\varepsilon^{n+1/2}\), and (A-D3a) carries it forward, \(e^{n+1}=2\varepsilon^{n+1/2}-e^n\), bounded by \(2\sum\lvert\varepsilon\rvert\). BDF schemes evaluate \(\tilde u^{n+1}\) from the converged state, so their defect is only the current \(\varepsilon\).
- *Galerkin discrete divergence* (\(\ell(u)=(q_h,\nabla\cdot u)\), \(c=0\)). In DOF vectors, split \(u_h^n=w^n+M^{-1}B^{\mathsf T}\varphi^n\) with \(Bw^n=0\). Here M is the velocity mass matrix and B is defined by \(q^{\mathsf T}Bv=(q_h,\nabla\cdot v_h)\); assume compatible Dirichlet data. The history part of \((D_tu_h,v_h)\) then contains \(-(2/\delta t)\,v^{\mathsf T}B^{\mathsf T}\varphi^n\), a discrete gradient. It combines with \(-(p_h,\nabla\cdot v_h)\), so the pressure absorbs it: \(p_h^{n+1/2}=\hat p^{\,n+1/2}-2\varphi^n/\delta t\), where \(\hat p\) solves the system with \(w^n\) in place of \(u_h^n\). Meanwhile \(u_h^{n+1/2}\) depends only on \(w^n\), and \(\varphi^{n+1}=-\varphi^n\). The endpoint velocities therefore keep the gradient component \((-1)^nM^{-1}B^{\mathsf T}\varphi^0\), and the reported pressure carries \((-1)^{n+1}2\varphi^0/\delta t\), which grows as \(\delta t\to0\) on a fixed mesh. Typical sources of \(\varphi^0\ne0\) are interpolating a smooth divergence-free field into Taylor–Hood and an impulsive start (fluid at rest, inflow data imposed at t = 0). V2 is not affected: \(u^0=(y,x)\) lies in every FE space and is exactly divergence-free. For the VMS cases the continuity equation contains ũ, so it does not constrain u\_h alone. The sign-flip mechanism is the same, but I have not quantified its size.

There are two ways to restore the spec's claim. (i) Construct a compatible \(u_h^0\): Dirichlet values applied and, for Galerkin, discretely divergence-free, for example by a Stokes-type projection of \(u_0\). (ii) Start Crank–Nicolson with one BDF1 step. That step imposes every constraint at \(t^1\) and keeps second order, because its local error is \(O(\delta t^2)\). Which one to use is a decision for the spec owner (open question 1).

**A-D4. Startup orders.**

*Lower-order startup.* For a zero-stable linear k-step method of order p, the global error at fixed T is bounded by \(C\bigl(\max_{j<k}\lvert e^j\rvert+\delta t^{p}\bigr)\). From an exact \(u^0\), one BDF1 step gives \(e^1=O(\delta t^2)\). A BDF2 step from \(u^1,u^0\) gives \(e^2=O(\delta t^3)+O(e^1)=O(\delta t^2)\). For BDF2 (k = 2) this means order 2. For BDF3 (k = 3), the starting values \(u^0,u^1,u^2\) carry \(O(\delta t^2)\) errors. These excite the principal root and propagate like a perturbation of the initial data, so the global order is generically 2, as the spec says.

*Exact-history startup with zero subscale histories.* Exact velocity histories make the velocity starting error zero. The subscale has no exact solution. Its model ODE for dynamic ASGS, \(\partial_t\tilde u+\tilde u/\tau_m=R_u\), starts from \(\tilde u(0)=0\) with slope \(R_u(0)\). Unless \(R_u(0)=0\), zero values at \(t^{-1}\) and \(t^{-2}\) are not samples of a smooth continuation. With \(\tau_m\gg\delta t\), BDF3 at step 1 gives

\[
\tilde u^1=\tau_tR_u^1\approx\tfrac6{11}\,\delta t\,R_u
\qquad\text{instead of}\qquad
\tilde u(\delta t)\approx\delta t\,R_u ,
\]

a first-order starting error. A BDF1 first step from \(\tilde u^0=0\) gives \(\tilde u^1\approx\delta t\,R_u\), which is consistent. In V2 the exact velocity lies in the FE space, and the forcing \(f=g'(y,x)+g^2(x,y)\) is in the P1 space. The converged residual \(R_u\) is then only the temporal truncation error, \(O(\delta t^3)\), and the subscale defect is \(O(\delta t^4)\). Third order with exact-history startup is therefore a property of V2's construction, which is all the spec claims.

**A-D5. Stability facts.**

The boundary locus of BDF-k is \(z(\vartheta)=\sum_{j=1}^k(1-e^{-i\vartheta})^j/j\). For BDF2, \(\operatorname{Re}z=(1-\cos\vartheta)^2\ge0\), so the locus stays in the closed right half-plane and BDF2 is A-stable. For BDF3, the locus has \(\operatorname{Re}z<0\) for \(0<\vartheta<1.318\). The arc ends where \(\cos\vartheta=\tfrac14\), at \(\lvert\operatorname{Im}z\rvert=\tfrac{\sqrt{15}}2\approx1.94\), and \(\operatorname{Re}z\) reaches its minimum \(-\tfrac1{12}\) at \(\vartheta=\pi/3\). The scratch computation confirms both values to round-off. So BDF3 is unstable only in a bounded set next to the imaginary axis, including the segment \(0<\lvert\operatorname{Im}z\rvert<1.94\). It is stable for all \(\operatorname{Re}z<-\tfrac1{12}\), and A(α)-stable with \(\alpha\approx86.03^\circ\). On the imaginary axis the growth factor per step is about \(1+y^4/4\) for small \(y=\omega\delta t\): 1.011 at y = 0.5 and 1.044 at y = 1. An undamped mode therefore grows at a rate of about \(\omega^4\delta t^3/4\) per unit time. Reducing δt lowers the rate but never removes it, and a physical damping rate above that value suppresses it.

Both Crank–Nicolson forms, one-leg and trapezoidal, share the stability function \(R(z)=(1+z/2)/(1-z/2)\). It has \(\lvert R\rvert\le1\) exactly for \(\operatorname{Re}z\le0\), and \(R\to-1\) as \(z\to-\infty\): A-stable, with no damping of stiff modes.

The 2007 analysis tests backward Euler with \(g^{n+1}\):

\[
(g^{n+1}-g^n,g^{n+1})=\tfrac12\lVert g^{n+1}\rVert^2-\tfrac12\lVert g^n\rVert^2+\tfrac12\lVert g^{n+1}-g^n\rVert^2 .
\]

This identity gives the increment sums in Theorem 1 and the bound (49). Testing the \(\theta=\tfrac12\) one-leg form with \(g^{n+1/2}\) gives \((g^{n+1}-g^n,g^{n+1/2})=\tfrac12\lVert g^{n+1}\rVert^2-\tfrac12\lVert g^n\rVert^2\), with no increment term. This is why (49), and with it the paper's route to Theorem 3, is unavailable for Crank–Nicolson. The paper neither proves nor disproves a Crank–Nicolson analogue of Theorem 1, and this pass does not prove one.

**A-D6. Galerkin baseline as the master equation with stabilization removed.**

Remove three terms from the master FE equation at lines 126–129: \(\chi_A\chi_D(D_t\tilde u,v_h)\), \(\sum_K(\tilde u,L_a^*(v_h,q_h))_K\) and \(\sum_K(\tau_c(\nabla\cdot u_h-\zeta_h),\nabla\cdot v_h)_K\). This is the same as setting \(\tilde u\equiv0\) and \(\tau_c\equiv0\). Then \(a=u_h+\tilde u=u_h\) under either splitting, so the linear/nonlinear switch has no effect. \(\xi_h\), \(\zeta_h\), \(\tau_m\) and \(\tau_t\) no longer appear. What remains is

\[
(D_tu_h,v_h)+\nu(\nabla u_h,\nabla v_h)+c(u_h;u_h,v_h)-(p_h,\nabla\cdot v_h)+(q_h,\nabla\cdot u_h)=\langle f,v_h\rangle,
\]

which is line 177. The same removal in the Picard system at lines 224–228 gives \(c(\hat a;u_h^j,v_h)\) with \(\hat a=u_h^{j-1}\), including the outflow term, as the spec states. F contains no τ, so the frozen-τ and exact-τ Newton variants coincide. The time schemes carry over unchanged (A-D1), and so does the Crank–Nicolson caveat of A-D3, which is cleanest in exactly this case.

On affine elements with velocity degree \(m=k+1\), the integrand degrees are: convection \(3m-1\), mass \(2m\), viscous and pressure–divergence \(2m-2\), and the outflow facet term \(3m\). A rule of degree \(3m=3(k+1)\) integrates them exactly, and it integrates \(\langle f,v_h\rangle\) exactly for polynomial f of degree at most 2m.

**Disagreements and open questions from this pass**

1. *BDF3 stability wording (A11a), line 91.* Replace "its stability region misses a sector near the imaginary axis" with: "it is A(α)-stable with α ≈ 86°; with z = λδt it is unstable in a bounded region next to the imaginary axis (−1/12 < Re z < 0, |Im z| < 1.94), which includes part of the axis itself".
2. *Crank–Nicolson attribution (A11b), line 91.* Replace "the unconditional stability result of Codina et al. 2007 does not extend to it (their remark after the proof of Theorem 2)" with: "Codina et al. 2007 analyse backward Euler only; their remark after the proof of Theorem 2 notes that the increment bound (49), which backward Euler provides and their Theorem 3 uses, is not expected for Crank–Nicolson".
3. *Missing Crank–Nicolson condition (A06, A08), lines 87 and 89.* "So the extrapolated \(u_h^{n+1}\) matches the Dirichlet data exactly" holds only if \(u_h^0\) does. For Galerkin, a \(u_h^0\) that is not discretely divergence-free leaves a non-decaying pressure component of order \(1/\delta t\). This needs a decision, not just rewording: (i) require and construct a compatible \(u_h^0\), or (ii) start Crank–Nicolson with one BDF1 step. Option (ii) is simpler, keeps second order, and also covers impulsive starts. It changes "Crank–Nicolson needs no startup". V2 is unaffected either way. The thesis lattice vortex (thesis §6) used Crank–Nicolson and is affected.
4. *Exact-history levels (A10), line 89.* The spec does not say which levels exact-history startup fills: \(u^{-1},u^{-2}\) at \(t=-\delta t,-2\delta t\) with BDF3 from step 1, or \(u^1,u^2\) with BDF3 from step 3. Both give third order in V2. The first matches the word "histories" and the "subscale histories start at zero" sentence, and it works because every planned exact solution is defined for t < 0.
5. *Stale references (A14), wording only.* Line 5 says "at time \(t^{n+1}\)" and "Every variant" (the baseline is not a switch setting). Line 100 says "at \(t^{n+1}\)" (should be "at the solve level"). Lines 162–168, 117 and 239 list histories only to \(n-1\) (BDF3 also needs \(n-2\)).
6. *Recommendation, not a disagreement (A07).* With Crank–Nicolson the OSS orthogonality defect is carried forward rather than reset each step, so use tight outer tolerances.

No other disagreement was found. The addendum introduces no change to the master formulation, the eight cases, the OSS projections or the algorithm beyond the time level at which they are evaluated.

---

**Final pass — 21 October 2026**

The master formulation, all eight cases, BDF coefficients and history terms, weighted OSS projections, nonlinear subscale update, and Picard equations pass the algebraic verification under the declared modelling assumptions. The eight errata rows and the additional Chapter 4 corrections also agree with the independent derivation. No additional formulation error was found in this pass.

**One qualification remains in the quadrature-exactness sentence:** “a polynomial force” is insufficient without a degree bound; dynamic histories also need an appropriate polynomial bound. This affects the exact-integration claim, not the validity of the common-quadrature OSS identity or the master equations. Precise replacement wording is given below. The other outstanding findings from the preceding review are resolved.

Reviewed [spec revision](</home/ali/04_applications/FEops_prep/VMS_masters_thesis/VMS formulations for incompressible Navier–Stokes implementation spec.md>) SHA-256: `91f6d5936052f58c551f7e20bfd3e7ac864001427b0a31a3938073bd4e61e94e`. Line references identify this revision. I rechecked the source equations using the non-layout PDF text extractions and the retained derivations. Chapter 5 and the deferred EMAC material remain outside this review. The spec itself was not edited; this pass updates the verification record only.

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| P01 — master and eight cases | Corrected adjoint sign; common BDF; ASGS dynamic mass coupling; reduced-residual OSS; static/dynamic and linear/nonlinear switches. | `colomes_VMS_paper.pdf`, (13)–(19), (22), (24)–(25); `time_dependent_2007.pdf`, (9)–(17); `chapweb09.pdf`, (57)–(61); independent D1–D4 below. | Agree | All eight rows remain consistent. The global instruction to remove \(D_tu_h\) for OSS is essential when reading the compact master subscale expression. |
| P02 — positive common quadrature | Degree \(3k\), positive weights, same rule for coupling, projections, diagnostics and histories. | `colomes_VMS_paper.pdf`, (33), (53), requires weighted projection; the positivity argument below checks the chosen discrete version. | Agree — previous R02 resolved | For affine Pk spaces, degree \(3k\ge2k\) plus positive weights gives nondegenerate sampling; positive \(\tau_t\) then gives a positive-definite weighted mass matrix. Omit the scalar projection when \(\tau_c=0\), as specified. |
| P03 — quadrature exactness, line 103 | With linear splitting, affine Pk, elementwise constant \(\tau\), and a polynomial force, degree \(3k\) integrates the products exactly for k≤2. | `colomes_VMS_paper.pdf`, (17), includes force in the residual; (53) includes subscale history. Neither source supplies an arbitrary-polynomial exactness theorem. Degree count in R-D1 below. | Partially resolved — remaining qualification | Bound the degree of the force and, for dynamic cases, the subscale histories. Merely being polynomial is not enough. See the replacement sentence below. |
| P04 — OSS iteration comparison | Same iteration structure as Colomés, but different lagged full/reduced-residual iterates; same fixed-point equations. | `colomes_VMS_paper.pdf`, (52)–(53); D3b and R-D3 below. | Agree — previous R03 resolved | The displayed positive difference is correct when read as \(s_{\rm red}-s_{\rm full}\), with otherwise matched inputs. Defining that direction explicitly would make the sentence unambiguous. It is not an equality between arbitrary distinct solver trajectories or different time schemes. |
| P05 — Newton/Picard and stopping | Conditions Newton convergence on smoothness and a nonsingular Jacobian; no guaranteed convergence; linear-splitting ASGS must converge its inner solve. | `Kabalan_Ali_Thesis.pdf`, §3.2, (3.27)–(3.29); master equations; R-D4 below. | Agree — substantive R05–R06 findings resolved | Both linearizations target the same equations. “Loses the quadratic rate” for lagged \(\tau\) is more precisely “does not guarantee quadratic convergence”; see wording note below. |
| P06 — boundaries | Explicit zero gradient-form traction; correct boundary energy flux; monitor backflow instead of assuming it cannot occur. | `colomes_VMS_paper.pdf`, (8), (10), with the declared boundary extension; D6 and R-D4 below. | Agree — previous R07 resolved | No claim of a verified cylinder flow or backflow-free solution remains. |
| P07 — all eight errata rows and further thesis errors | Includes both chapter time omissions, both errors in thesis (4.39), distinct thesis residual-sign patterns, and the time-scheme inconsistency. | `colomes_VMS_paper.pdf`, (38)–(39), (48)–(54); `chapweb09.pdf`, (57)–(64); `Kabalan_Ali_Thesis.pdf`, (4.1)–(4.43); D1–D5 and D8 below. | Agree | The retained ambiguity about a possible OSS specialization in chapter (63)–(64) is appropriate. The force-transfer qualification for Colomés (39) is correct. |
| P08 — provenance, scope and code impact | Choices marked as ours; split OSS excluded; impact on original thesis software explicitly unknown. | Sources and distinctions checked in F01–F14, U01–U08 and D1–D8 below. | Agree — previous R08 resolved | Manuscript corrections do not establish what existing software implements. This audit makes no such claim. |

The remaining quadrature statement can be replaced with:

> For linear splitting on affine Pk elements, k=1 or 2, with stabilization parameters constant within each element, degree q=3k is sufficient for exact volume integration provided the force and any stored subscale histories are represented on each element by polynomials of degree at most 2k−1. Otherwise the rule defines a quadrature approximation. The discrete OSS orthogonality identity still follows from using the same rule for the projection and coupling.

This is a sufficient degree restriction, not a claim that it is necessary. If the maximum degree of the force/history contribution is r, the subscale degree is bounded by \(\max(2k-1,r)\). The largest convective-adjoint stabilization product is therefore bounded by

\[
\max(2k-1,r)+(2k-1).
\]

Taking \(r\le2k-1\) gives \(4k-2\le3k\) for k=1,2. Without a restriction on r, exactness does not follow. For example, at k=1, a force component \(x^4\) multiplied by the constant pressure-test gradient produces a degree-four stabilization load; a degree-three rule does not guarantee its exact integration. Histories stored only as quadrature values have no assumed continuous polynomial representation, so an exact continuous-integral claim about them also needs qualification. This does not impair the algebraic projection cancellation at those points.

The new positivity requirement resolves the previous mass-matrix concern. For any nonzero FE polynomial z on an affine element, positive weights and degree \(3k\) give

\[
\sum_j\omega_j z(x_j)^2=\int_K z^2>0.
\]

If \(\tau_j>0\), then

\[
\sum_j\omega_j\tau_j z(x_j)^2
\ge\min_j\tau_j\int_K z^2>0.
\]

This proves positive definiteness of the computed weighted mass matrix on the stated FE space, rather than merely assuming that positivity of the coefficient survives an arbitrary quadrature rule.

Two wording refinements would remove unnecessary absolute claims without changing any equation:

- At line 192, use “does not guarantee quadratic convergence” instead of “loses the quadratic rate.” A lagged-coefficient Jacobian usually gives only linear convergence, but special cases in which the omitted derivative vanishes can retain a faster rate.
- At line 206, write \(s_{\rm red}^i-s_{\rm full}^i=\tau_t^i(\alpha_0/\delta t)(u_h^i-u_h^{i-1})\) and specify otherwise matched inputs. This fixes the direction of the comparison already derived in R-D3; no sign error is being asserted in the current wording.

**Final disposition:** the formulation and source-errata checks pass. The one remaining correctness qualification is the bounded-degree/history condition for the exact-quadrature claim. The two wording refinements above are not formulation defects. No code or implementation plan is part of this verification.

---

**Second review — superseded spec, retained for evidence**

The U/R findings and current-status language in the following review refer to the previous spec hash, not the final-pass revision. Use P01–P08 above for the latest disposition.

**Re-verification of the updated VMS formulation spec — 21 October 2026**

The revised master equations, eight case equations, reduced-residual OSS projections, and new Picard equations pass the algebraic check under the stated assumptions. Most findings from the first review are addressed. Split OSS has been removed, so its previous findings are now outside the current spec's scope, not outstanding defects in it.

Two new statements need correction: constant stabilization parameters do not alone make degree-\(3k\) quadrature exact for all the listed formulations, and one Picard solve per outer iteration does not exactly reproduce the paper's lagged full-residual OSS iteration. The quadrature rule also needs an admissibility condition, and a few convergence and implementation-impact statements need qualification.

This review refers to [the updated spec](</home/ali/04_applications/FEops_prep/VMS_masters_thesis/VMS formulations for incompressible Navier–Stokes implementation spec.md>), SHA-256 `c1cde02b095baabbbaa509f105b0c5bd9053697135b4ed1987f5be702d92bb43`. Line numbers below identify that revision. Source PDFs and the independent derivation D1–D8 remain as in the original review preserved below. I rechecked the relevant Colomés §4 equations and the newly cited thesis §3.2.2 using the existing non-layout `pdftotext` extractions. Chapter 5 and the deferred EMAC discussion remain outside this review. No implementation code was examined or changed.

The status of the original findings is now as follows.

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| U01 — fixed-choice provenance | Explicitly marks our selections, distinguishes eight admissible cases from six benchmarked cases, allows source-dependent interpolation and parameter conventions, and corrects the channel constants. | `colomes_VMS_paper.pdf`, §§2.2, 4, 7.1, Table 1; `onseen_2008.pdf`, §1; `dissipative_2010.pdf`, (12). | Agree — resolved | Original F01–F02, F08–F09, F12–F14 are addressed. New Newton/Picard and quadrature selections are explicitly ours. |
| U02 — boundary convention | Specifies zero gradient-form traction, pressure normalization where appropriate, and the open-boundary energy flux. | `colomes_VMS_paper.pdf`, (8), (10); independent identities D6a–D6b. | Agree — resolved at the formulation level | The old F03–F04 ambiguity is resolved. The claimed guarantee against backflow is a separate issue in R07 below. |
| U03 — master and eight cases | Retains the corrected residual/adjoint signs and same-scheme BDF histories; limits the SUPG/PSPG-only identification to affine P1. | `colomes_VMS_paper.pdf`, (13)–(17); `time_dependent_2007.pdf`, (9)–(17); derivation D1–D4. | Agree | C01–C08 still pass. M06 is corrected. The compact master formula's OSS residual is still interpreted using its explicit instruction to remove \(D_tu_h\). |
| U04 — OSS and stopping conditions | Uses common quadrature, the reduced projection, correct residual in step 3, outer iterations for every OSS case, and checks pressure, both projections, subscale, and FE residual. | `colomes_VMS_paper.pdf`, (33), (53)–(54); D3; `onseen_2008.pdf`, (18), Remark 1. | Agree — resolved, with R02 qualification | Original M08–M11 are addressed. The full/reduced projection relation is now explicit. Exact orthogonality means the common computed bilinear form, up to solve tolerance. |
| U05 — evidence for target variant | Reports no significant cost difference between dynamic linear and nonlinear OSS; identifies Dyn-Nl-OSS as our target choice. | `colomes_VMS_paper.pdf`, §5.2.2 and Fig. 4. | Agree — resolved | Original M12 is addressed. |
| U06 — all eight revised errata rows | Distinguishes semi-discrete (38), includes both chapter omissions, separates thesis (4.9) from (4.11)/(4.15), includes the extra RHS \(\tau_t\), and qualifies the mixed time schemes. | `colomes_VMS_paper.pdf`, (13)–(17), (38)–(39), (48)–(54); `chapweb09.pdf`, (57)–(64); `Kabalan_Ali_Thesis.pdf`, (4.9), (4.11), (4.15), (4.28)–(4.39); D1–D5, D8. | Agree — resolved | The original E01–E08 findings are now represented correctly, including the general force-transfer qualification for (39). |
| U07 — additional Chapter 4 errors | Adds the test-index, integration-by-parts, continuity, projection-closure, implicit-subscale, and tensor-expansion corrections. | `Kabalan_Ali_Thesis.pdf`, (4.1)–(4.18), text after (4.34), (4.40)–(4.43); D8. | Agree — resolved as manuscript findings | Original T01–T12 are covered between the errata and the new table. “Unconstrained” is the selected projection space; the correction forced by (4.40)–(4.41) alone is to identify finite-dimensional FE test spaces. Existing code impact remains unverified; see R08. |
| U08 — split OSS | Removed because it is not used in the thesis. | The original comparison was with `term_by_term.pdf`, (19)–(25). | Removed from scope | Original S01–S09 and D7 are retained below only as historical evidence. They are not current blockers. The remaining reference to (23) for degree scaling is correct. |

These are the findings on the added or revised material.

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| R01 — quadrature exactness, line 103 | For constant \(\tau\), degree \(3k\) integrates Galerkin and stabilization products exactly for k≤2. | `colomes_VMS_paper.pdf`, (15), (17), (25), gives an implicit nonlinear subscale; it does not assert this quadrature exactness. Independent degree count and counterexample R-D1 below. | Disagree as a statement covering all eight cases | It holds under additional polynomial, affine-geometry, and linear-splitting assumptions. Arbitrary force/history or nonlinear subscales can produce nonpolynomial integrands even with constant \(\tau\). |
| R02 — quadrature admissibility, lines 25 and 103 | Choosing a common rule of degree \(3k\) supplies a discrete inner product for weighted projections. | `colomes_VMS_paper.pdf`, (33), (53), and `onseen_2008.pdf`, (13), (18), use positive weighted inner products. A quadrature degree specifies exactness, not positivity with variable weights. | Incomplete / rule-dependent | The rule must make the computed weighted mass matrices positive definite. Positive quadrature weights and nondegenerate FE sampling suffice. Degree alone is not enough; R-D2 gives a concrete counterexample. The algebraic orthogonality cancellation itself remains correct if the projection equations are solved. |
| R03 — exact reproduction of Colomés, line 206 | One Picard inner solve per outer iteration reproduces Colomés §4 “exactly.” | `colomes_VMS_paper.pdf`, (52)–(53), lags a full-residual projection; the revised spec, correctly, lags its reduced-residual projection. Its own line 182 says these iterates differ. | Disagree | It has the same broad iteration structure, not identical OSS iterates. The difference includes \(-\tau_t(\alpha_0/\delta t)(u_h^i-u_h^{i-1})\); see R-D3. BDF2 and the different projection spaces also prevent a general exact-reproduction claim. |
| R04 — new Picard equations, lines 194–204 | Freeze every advecting velocity, stabilization coefficient, and the outflow advector; retain the transported velocity and pressure as unknowns. | `colomes_VMS_paper.pdf`, §4 and (13)–(17); `Kabalan_Ali_Thesis.pdf`, (3.27), for the basic Picard principle. | Agree | The displayed system is affine linear in the new unknowns and reduces to the master system at a common fixed point. The \(\chi_A\) factor correctly removes the coarse time derivative only for OSS. |
| R05 — Newton/Picard solution and rates, lines 188, 192, 204 | The converged time step is identical; Newton is quadratic near the solution and Picard is linear but robust. | Thesis §3.2.2 discusses the usual local Newton/Picard behavior, not uniqueness or unconditional convergence of this VMS system. `colomes_VMS_paper.pdf`, §4, does not prove those guarantees. | Qualified | The methods target the same equations. Identical computed roots require convergence to the same solution branch; rates require the usual smoothness, nonsingular-Jacobian or contraction hypotheses. The exact Newton derivative may be undefined at zero advecting speed; freezing \(\tau\) is an inexact-Newton option, not a quadratic-convergence guarantee. |
| R06 — one Picard solve and ASGS stopping, lines 210–239 | Permits one inner iteration per outer iteration, while linear-splitting ASGS runs only step 2 (and storage for dynamic ASGS). | The master equation still contains \(u_h\cdot\nabla u_h\) and velocity-dependent \(\tau\); `colomes_VMS_paper.pdf`, (24), is linear scale splitting, not a linear Navier–Stokes problem. | Ambiguous | These two ASGS cases need step 2 to converge internally. The one-solve option is valid only if an enclosing iteration still checks and converges the nonlinear residual; it cannot mean stopping after one Picard linear solve. |
| R07 — cylinder backflow, line 101 | A sufficiently long channel avoids outlet backflow. | The cited equations supply the boundary-energy identity, not a guarantee for an unspecified cylinder domain, Reynolds number, or solution. | Unverified assumption | The energy sign is correct. Absence of outlet backflow is a condition on the eventual solution; it has not been established by this source/formulation audit. |
| R08 — effect on existing code, lines 259 and 263 | The first errata row is the one that matters for code; the further Chapter 4 errors “do not change the code.” | The papers and thesis show printed equations. They do not establish which expressions any existing implementation uses. | Unsupported implementation claim | The additional manuscript findings do not alter the already corrected master equations. Their impact on existing software is unknown. The \(\tau_t\), residual, and history errors in other errata rows would also matter if implemented literally. |

The remaining mathematical disagreements can be demonstrated without introducing an implementation plan.

**R-D1. What the degree count actually proves.**

For affine Pk elements with **linear scale splitting**, \(a=u_h\) has degree k and \(a\cdot\nabla u_h\) has degree at most \(2k-1\). If the stabilization parameters are constant on each element and the force/history terms have degree no greater than \(2k-1\), the subscale expression has degree at most \(2k-1\). Consequently the largest adjoint stabilization product has degree

\[
(2k-1)+(2k-1)=4k-2,
\]

whereas the Galerkin convection products have degree \(3k-1\). Degree \(3k\) covers both for k=1,2. These are sufficient conditions for the asserted polynomial exactness; constant \(\tau\) alone is not. For example, static linear ASGS evaluated at zero velocity and pressure with \(f=(\sin x,0)\) has constant \(\tau_m=h^2/(c_1\nu)\) on each element, yet its pressure-stabilization load contains \(\tau_m\sin x\,\partial_xq_h\), which is not polynomial even for P1.

For a counterexample to the nonlinear-splitting claim, even grant a constant scalar \(\tau>0\), take a quadratic velocity \(u_h=(x^2,0)\), and set pressure, force, and time derivative to zero. The static nonlinear subscale closure gives

\[
s_2=0,\qquad s_1=\tau\bigl[-2x(x^2+s_1)+2\nu\bigr],
\]

so

\[
s_1=\frac{2\tau(\nu-x^3)}{1+2\tau x}.
\tag{R-D1}
\]

On an element where the denominator is nonzero this is generically a rational, nonpolynomial function. Its advection and stabilization products therefore have no polynomial degree bound of six. This is a check of integrands for admissible FE arguments, not a claim that this trial velocity solves the global incompressibility equation. Arbitrary L2 force or quadrature-point histories are independent reasons why polynomial exactness need not hold.

The common quadrature identity still holds:

\[
(s,w_h)_Q=(\tau_tX,w_h)_Q-(\tau_t\xi_h,w_h)_Q=0.
\]

It is an algebraic consequence of the projection equations with the same rule. It does not require exact integration of a nonpolynomial continuous integrand, and should not be confused with that stronger claim.

**R-D2. Exactness degree does not ensure a positive weighted mass matrix.**

Consider a triangle and a rule with points given in barycentric coordinates by the centroid and the three permutations of \((3/5,1/5,1/5)\). Give the centroid normalized weight \(-27/48\), and each remaining point weight \(25/48\). Multiplying weights by the element area gives a rule exact for polynomials of degree three; this follows directly by checking the barycentric monomials through degree three.

For a positive coefficient \(\tau\) equal to 4 at the centroid and 1 at the other three points, the computed squared norm of the FE constant function is

\[
(\tau\,1,1)_Q
=|K|\frac{-27\cdot4+3\cdot25}{48}
=-\frac{11}{16}|K|<0.
\tag{R-D2}
\]

Such a coefficient pattern is compatible with the spatial parameter formula: on an equilateral element, an affine velocity that vanishes at the centroid can have equal nonzero speeds at the other three points, giving larger \(\tau_m\) at the centroid. Thus this is not excluded by the pointwise-parameter choice. The example shows why “degree three” does not by itself specify an admissible weighted inner product. It does not assert that the rule eventually chosen will have this defect; that choice has not been given.

**R-D3. Reduced and full lagged projections have different intermediate subscales.**

Let i denote the outer iteration, with time histories held fixed. Write \(\lambda^{i-1}\) for the full-residual projection in the paper and \(\xi^{i-1}\) for the spec's reduced projection. Even under otherwise identical weights and spaces,

\[
\lambda^{i-1}=\xi^{i-1}-D_tu_h^{i-1}.
\]

The correctly assembled full-residual update with that lagged projection is

\[
s_{\mathrm{full}}^i
=\tau_t^i\bigl[f-D_tu_h^i-L_{\hat a}(u_h^i,p_h^i)
+H_s-\lambda^{i-1}\bigr],
\]

whereas the spec's update is

\[
s_{\mathrm{red}}^i
=\tau_t^i\bigl[f-L_{\hat a}(u_h^i,p_h^i)+H_s-\xi^{i-1}\bigr].
\]

Subtracting gives

\[
\boxed{
s_{\mathrm{full}}^i-s_{\mathrm{red}}^i
=-\tau_t^i\bigl(D_tu_h^i-D_tu_h^{i-1}\bigr)
=-\tau_t^i\frac{\alpha_0}{\delta t}(u_h^i-u_h^{i-1}).
}
\tag{R-D3}
\]

This generally survives before convergence even for BE. It vanishes at the fixed point. Therefore the new “exact reproduction” sentence conflicts with both the derivation and the spec's own correct explanation of its projection convention.

**R-D4. Check of the new Picard system and boundary convention.**

At fixed outer data and frozen \(\hat a,\hat\tau\), the subscale in the new displayed Picard system is affine in \((u_h^j,p_h^j)\). Its test operator is fixed, the divergence stabilization is linear, and freezing the normal advector in the outflow term makes that term linear too. Setting consecutive inner and outer iterates equal recovers the master equation. Newton and Picard thus share the target equations; this argument does not establish uniqueness, convergence from every initial guess, or a common basin of attraction.

For smooth a, integration by parts in the convection form shows

\[
\tfrac12(a\cdot\nabla u,v)-\tfrac12(a\cdot\nabla v,u)
+\tfrac12\langle(a\cdot n)u,v\rangle
=(a\cdot\nabla u,v)+\tfrac12((\nabla\cdot a)u,v).
\]

Hence the added half-flux removes the convective boundary contribution associated with writing the skew form in strong form. The remaining natural traction from viscosity and pressure is the declared gradient-form traction. This does not identify the volume form with ordinary convection for a non-solenoidal discrete advector. For discontinuous nonlinear a the identity must be read elementwise with the already declared subscale boundary approximations. Testing the specified convection form with u on both slots directly gives the stated boundary energy flux, without any integration-by-parts assumption.

The current disagreement list is R01 (overbroad quadrature exactness) and R03 (exact algorithm reproduction). R02, R05–R08 are missing conditions or unsupported claims, not errors in the eight case equations. No new sign, BDF-history, or OSS-weight defect was found. The original derivation is retained next so the source checks remain reviewable.

---

**Original review — superseded spec, retained for derivations and source evidence**

Everything from this point onward describes the **pre-update** spec. Its “Spec says,” unresolved findings, and split-OSS discussion are historical. Use U01–U08 and R01–R08 above for the current status.

**Independent verification of the original VMS formulation spec — 21 October 2026**

The master equations and all eight case equations are algebraically consistent with the stated choices, provided the OSS projections are solved to convergence in fixed spaces. The spec is not yet fully verified: its split-OSS parameter formula is wrong for degree greater than one; its split-OSS energy claim does not follow for the variable parameters it permits; several source attributions and errata descriptions need correction; and Chapter 4 contains further errors.

I extracted all seven supplied PDFs with `pdftotext` in its default, non-layout mode. I also inspected rendered original pages where signs, parentheses, overbars, or exponents mattered. The derivation below starts from the semi-discrete equations, not the printed fully discrete equations. This review covers thesis Chapter 4; Chapter 5 is deferred. It does not establish what any existing software actually implements.

Equation references below use the supplied filenames. “Agree” in an errata row means agreement with the spec's diagnosis, even though the printed source equation is wrong. “Qualified” and “Ambiguous” distinguish missing conditions or modelling decisions from algebraic errors. All derivations assume constant viscosity, a fixed spatial mesh and projection spaces, a uniform time step for the stated BDF2 coefficients, and sufficiently regular data for the indicated inner products. The spec does not explicitly state all these conditions.

The following table checks every fixed-choice row, including whether its cited source actually supports the choice.

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| F01 — variants | Eight combinations, citing Colomés Table 1, plus split OSS. | `colomes_VMS_paper.pdf`, §§2.2.1–2.2.3, eqs. (22), (24)–(25), (31)–(32), §4.1: the three independent choices exist. Table 1 actually contains **six** cases, omitting both static-nonlinear cases. `term_by_term.pdf`, (19)–(22), supplies a different term-by-term method. | Qualified | Eight combinations are valid, but Table 1 is not an eight-case validation. Including all eight is our selection from the framework. |
| F02 — interpolation | Equal-order Pk/Pk; P1 first, P2 when the Laplacian is to be active; “all sources assume equal order.” | `time_dependent_2007.pdf`, §2, assumes continuous equal order. `onseen_2008.pdf`, §1 before (6), explicitly allows different polynomial orders and pressure spaces that need not be continuous. | Disagree with attribution | Equal order and the P1/P2 sequence are our choices, supported as admissible choices, not universal assumptions of the sources. The elementwise Laplacian vanishes for affine simplicial P1. |
| F03 — convection | Type-1 skew form with an outflow half-flux; “exactly energy-neutral,” citing Remark 2.1. | `colomes_VMS_paper.pdf`, (10), Remark 2.1 and (29): the volume type-1 form is skew and needs no derivative of the discontinuous advector. The exposition assumes homogeneous velocity on the whole boundary after (5). | Qualified / Disagree as stated | The volume choice is supported. The outflow extension is our decision, not established by the cited remark. With that extension, the energy contribution is the boundary flux, not zero; see D6. “Only common form” should be limited to the alternatives compared in the paper. |
| F04 — viscosity | Gradient form and elementwise Laplacian “match all sources.” | `colomes_VMS_paper.pdf`, (8), (14), (17), (19); `dissipative_2010.pdf`, (9)–(10); `time_dependent_2007.pdf`, (9), (15), (17): these volume operators agree. Colomés (4) initially specifies symmetric-stress traction, then excludes Neumann boundaries in the derivation. | Agree for volume; boundary ambiguous | The spec's open-boundary extension must identify its traction convention and whether traction is zero. A gradient-form weak equation alone has the natural traction \(\nu\partial_n u-pn\). The master RHS supplies no explicit nonzero traction load. |
| F05 — element boundaries | Bubble assumption; subscales zero on element boundaries; neglect interelement terms. | `time_dependent_2007.pdf`, §2 following (9)–(10); `chapweb09.pdf`, §5.4 after (58): this is an explicit approximation. | Agree | Supported modelling assumption. A pointwise algebraic subscale formula is not itself a constructed bubble with a verified trace; the boundary neglect is part of the model. |
| F06 — adjoint approximation | Discard the divergence-of-advector and boundary contributions. | `colomes_VMS_paper.pdf`, explanation following (13)–(14), explicitly introduces these approximations. `dissipative_2010.pdf`, (9), assumes a solenoidal advector. | Agree | The displayed \(L_a^*\) is the adopted approximate adjoint when \(\nabla\cdot a\ne0\), not the complete formal adjoint. See D1. |
| F07 — pressure subscale | Quasi-static pressure; drop the 2007 pressure-residual derivative. | `colomes_VMS_paper.pdf`, (16), and `dissipative_2010.pdf`, (11), use this closure. `time_dependent_2007.pdf`, (12), instead has \(\tilde p=\tau_2(r_p+\tau_1\partial_t r_p)\); (17) sets pressure subscales to zero for its analysis. | Agree | The spec explicitly declares the departure from the 2007 closure. It is a model choice, not an algebraic consequence of the full 2007 system. |
| F08 — stabilization parameters | Pointwise speed, cell diameter, c1=4, c2=2, spatial \(\tau_m\), nonzero \(\tau_c\) for laminar flow. | `dissipative_2010.pdf`, (12) and following text, supports the formula and 4/2 for linear elements. `time_dependent_2007.pdf`, §2 and (13)–(14), uses mesh diameters. `chapweb09.pdf`, (60), and `onseen_2008.pdf`, (11), instead use an element maximum speed. `colomes_VMS_paper.pdf`, §9 and footnote 4, distinguishes spatial \(\tau_m\) from time-integrated \(\tau_{m,t}\). | Qualified | The formula and constants are supported. Pointwise quadrature evaluation and the precise common mesh-length convention must be labelled our choices; the sources do not use one universal convention. Keeping \(\tau_c\) is explicitly a selection. |
| F09 — time scheme | Same BDF1/BDF2 on both scales, BDF1 startup, citing 2007 §3.1. | `time_dependent_2007.pdf`, §3.1 and (22)–(24), motivates the same scheme for commutation, uses a theta method, and explicitly allows different scale integrators. `term_by_term.pdf`, (8)–(9), gives BDF1/BDF2 for the FE velocity, with BE subscales in (24)–(25). | Qualified | Matching schemes is source-supported. Selecting BDF1/BDF2 and the startup is our decision; the citation does not mandate those particular methods. |
| F10 — OSS projection | Consistent weighted projection, dynamic weight \(\tau_t\), no Dirichlet or zero-mean restriction. | `colomes_VMS_paper.pdf`, (33)–(34), (53)–(54), and introduction: weighted consistent projections; its displayed test spaces are constrained \(V_{0,h},Q_{0,h}\). `onseen_2008.pdf`, (18)–(20), projects onto the larger spaces \(V_h,Q_h\); Remark 1 discusses boundary layers from constrained projection. | Agree as a combination of sources | The unconstrained choice is supported by the 2008 paper, not literally by the spaces printed in Colomés. It also permits reproduction of time-dependent Dirichlet lifts in \(D_tu_h\). |
| F11 — force | Keep the force inside the residual and OSS projection; it cancels if it is an FE function. | `colomes_VMS_paper.pdf`, (17), §2.2.3 after (32): projection of the force requires \(f\in L^2\); omitting it is discussed for weaker data. | Agree, with regularity condition | Cancellation is exact in the OSS **subscale forcing** when \(f\in\hat V_h\), not in the Galerkin load. Keeping a pointwise projected residual requires more than an arbitrary \(H^{-1}\) force. |
| F12 — nonlinear solver | Inner Newton in velocity/pressure, outer fixed point in subscale/projections, citing Colomés §4. | `colomes_VMS_paper.pdf`, opening of §4: Picard linearization of the convective term. (49)–(50) use local subscale iteration; text after (50) mentions Newton as an alternative for that local problem. §4.2 lags projections. | Disagree with attribution | Projection lagging and local fixed points are supported. The particular nested Newton architecture is our decision, not the algorithm stated by the citation. |
| F13 — initial subscale | \(\tilde u^0=0\). | `term_by_term.pdf`, §3.1 after (17): zero is admissible; (19)–(20) specify zero initial conditions for both split subscales. | Agree | It is an admissible chosen initial condition, not a universal requirement. It satisfies the initial orthogonality needed in D3. |
| F14 — turbulence parameters | Colomés uses c1=12, c2=2, \(\tau_c=0\) for turbulence. | `colomes_VMS_paper.pdf`, §§5.1, 7.1, 8: 12/2/0 for homogeneous turbulence, but 12/8/0 for the reported turbulent-channel calculations. | Disagree with generalization | The proposed 12/2/0 default can be our choice, but the source does not use it for every turbulent benchmark. |

The next table checks the master system, its eight cases, and statements surrounding them. Here \(R=f-D_tu_h-L_a(u_h,p_h)\), \(\hat R=f-L_a(u_h,p_h)\), and \(H_s=\tilde u_{\rm hist}/\delta t\).

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| M01 — FE master signs | \(+(\tilde u,L_a^*)\), ASGS dynamic mass coupling, and positive pressure stabilization in divergence-minus-projection convention. | `colomes_VMS_paper.pdf`, (13)–(19); `time_dependent_2007.pdf`, (9)–(17); `dissipative_2010.pdf`, (9)–(11). | Agree | Independent substitution gives exactly these signs; see D1–D3. |
| M02 — effective parameter/history | \(\tau_t=(\chi_D\alpha_0/\delta t+\tau_m^{-1})^{-1}\), with history only for dynamic subscales. | Derive from `colomes_VMS_paper.pdf`, (15), or `time_dependent_2007.pdf`, (11), or `chapweb09.pdf`, (58)–(61). | Agree | BDF1 and uniform-step BDF2 histories in the spec are correct. History terms, as well as \(\tau_t\), depend on the time step; the “only place” wording is accurate only if it means the definition of a stabilization **parameter**. |
| M03 — dynamic OSS | Project \(\hat R+H_s\) with weight \(\tau_t\); remove \(D_tu_h\) in both subscale and projection. | `colomes_VMS_paper.pdf`, (53), contains the full residual plus history, weighted by \(\tau_{m,t}\). | Agree | The spec is an exactly equivalent redefinition of the projection unknown at convergence. History cannot generally be removed from the weighted projection. See D3. |
| M04 — pressure projection sign | \(\zeta=\Pi_{\tau_c}(\nabla\cdot u_h)\), \(\tilde p=-\tau_c(\nabla\cdot u_h-\zeta)\). | `colomes_VMS_paper.pdf`, (18), (54), projects \(R_p=-\nabla\cdot u_h\). | Agree | The paper's \(\xi_c=-\zeta\). Its negative RHS pressure-projection term in (51) is not a further sign error. |
| M05 — steady limit | Dynamic and static equations have the same steady state. | `time_dependent_2007.pdf`, §3.2 after (28). | Agree | All stored states must have reached the same value. Even the dynamic weighted OSS projection then reduces to the static closure; see D3. |
| M06 — SUPG/PSPG reduction | Static linear ASGS is SUPG + PSPG + grad-div. | `onseen_2008.pdf`, (15), and `colomes_VMS_paper.pdf`, (14), retain the viscous test-function Laplacian. | Qualified | Exact as stated for affine P1, where \(\Delta_hv_h=0\). At higher order there is also the \(\nu\Delta_hv_h\) adjoint contribution; the blanket equivalence is incomplete. |
| M07 — polynomial scaling | Replace c1 by c1 k^4 and c2 by c2 k. | `colomes_VMS_paper.pdf`, §2.2 following (21); `term_by_term.pdf`, (23), (30), corroborates this scaling when length conventions are matched. | Agree | Applying the c1 replacement in the pressure-parameter denominator too gives \(\tau_c=h^2/(c_1k^4\tau_m)\). This agrees with (23) when \(h_1=h_2=h\). |
| M08 — single Newton solve in linear cases | One Newton solve per time step is exact in all linear cases. | `colomes_VMS_paper.pdf`, §4.2, lags and updates projections; `onseen_2008.pdf`, §2.3, exposes them as additional globally coupled unknowns. | Disagree with blanket claim | Linear **scale splitting** does not remove the OSS projection coupling or the Navier–Stokes nonlinearity. With the spec's frozen projections, linear OSS still needs the listed outer iteration. A monolithic or exactly eliminated projection system would be a different meaning of “one solve.” |
| M09 — nonlinear OSS update | Algorithm step 3 prints \(R_u\); the later accelerated formula explicitly says to omit \(D_tu_h\) for OSS. | `colomes_VMS_paper.pdf`, (52)–(53), consistently uses full-residual projections; the spec instead defines reduced-residual projections. | Ambiguous locally | The earlier global rule makes the intended algorithm correct. Read literally, step 3 with the defined full \(R_u\) and reduced \(\xi\) is inconsistent. The same reduced-residual convention must apply there. |
| M10 — convergence variables | Stop on changes of \(u_h,\tilde u^*,\xi\), without mentioning \(\zeta\). | `colomes_VMS_paper.pdf`, (53)–(54): both projection equations are part of the converged OSS problem. | Qualified | Exact fixed points are unaffected. The stated finite-tolerance criterion alone does not establish that the pressure projection or pressure unknown has converged. This is a limitation of the spec's convergence claim, not a new continuous equation. |
| M11 — exact orthogonality | A consistent weighted mass solve enforces exact orthogonality. | `colomes_VMS_paper.pdf`, introduction, (33), (53). | Agree with integration qualification | Exact in the integral inner product if evaluated exactly; otherwise exact in the common quadrature inner product used for both projection and FE coupling, up to solve tolerance. Consistent mass alone does not eliminate quadrature error. |
| M12 — comparative performance | Dyn-Nl-OSS is the most robust and cheapest per accuracy; linear ASGS variants diverge on coarse turbulent meshes. | `colomes_VMS_paper.pdf`, §5.2.2, Fig. 4, and §10: the two linear ASGS cases fail on the coarse DHIT example, and dynamic OSS performs well, but the dynamic linear and nonlinear OSS cases show no significant cost difference in that comparison. | Qualified / overstated ranking | The reported ASGS failures and advantage of dynamic OSS are supported. A unique superiority of Dyn-Nl-OSS over Dyn-Lin-OSS, or a universal cost-per-accuracy ranking, is not established there. Selecting Dyn-Nl-OSS as the target is our decision. |
| C01 — Sta-Lin-ASGS | \(a=u_h\), \(\tilde u=\tau_mR\), no subscale time term. | `colomes_VMS_paper.pdf`, (13), (17), (22), (24), (31). | Agree | D2 with \(\chi_D=0\), \(\xi=0\). |
| C02 — Dyn-Lin-ASGS | \(a=u_h\), \(\tilde u=\tau_t(R+H_s)\), keep \((D_t\tilde u,v_h)\). | `colomes_VMS_paper.pdf`, (13), (15), (17), (24), (31). | Agree | D2 with \(\chi_D=1\), \(\xi=0\). |
| C03 — Sta-Nl-ASGS | \(a=u_h+\tilde u\), \(\tilde u=\tau_m(a)R(a)\), no time term. | `colomes_VMS_paper.pdf`, (22), (25), (31), and end of §4.1. | Agree | Implicit local nonlinear relation. Not one of Table 1's six benchmark cases. |
| C04 — Dyn-Nl-ASGS | \(a=u_h+\tilde u\), \(\tilde u=\tau_t(a)(R(a)+H_s)\), keep time term. | `colomes_VMS_paper.pdf`, (13), (15), (17), (25), (31). | Agree | D2; this agrees with the semi-discrete sources, not the erroneous adjoint signs in (48). |
| C05 — Sta-Lin-OSS | \(a=u_h\), \(\tilde u=\tau_m(\hat R-\xi)\), no time coupling. | `colomes_VMS_paper.pdf`, (22), (24), (32)–(34). | Agree | Requires the converged weighted projection in D3. |
| C06 — Dyn-Lin-OSS | \(a=u_h\), \(\tilde u=\tau_t(\hat R+H_s-\xi)\), no time coupling. | `colomes_VMS_paper.pdf`, (15), (24), (53). | Agree | Weight \(\tau_t\), projected history, and orthogonal past subscales are essential. |
| C07 — Sta-Nl-OSS | \(a=u_h+\tilde u\), \(\tilde u=\tau_m(a)(\hat R(a)-\xi)\). | `colomes_VMS_paper.pdf`, (22), (25), (32)–(34). | Agree | Valid combination of the choices; not independently benchmarked in Table 1. |
| C08 — Dyn-Nl-OSS | \(a=u_h+\tilde u\), \(\tilde u=\tau_t(a)(\hat R(a)+H_s-\xi)\). | `colomes_VMS_paper.pdf`, (13)–(17), (25), (53). | Agree | At a converged common state for subscale, velocity, pressure, and projections. |

Every original Errata row is checked below. Rows E01–E03 cover its first row, E04 its second, E03 also its third, E05 its fourth, E06–E07 its fifth, and E08 its sixth.

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| E01 — Colomés (38) | Wrong positive adjoint terms; describes it together with fully discrete equations using \(\tau_t u_h/\delta t\). | `colomes_VMS_paper.pdf`, (38), is a **semi-discrete static energy balance** with \(\tau_m P(\partial_tu_h)\) and \(\tau_m P(L_au_h)\), not that fully discrete expression. These terms and its RHS adjoint-force term are printed with plus signs. | Agree on signs; Disagree on transcription/classification | All three adjoint products require minus signs when rearranged as printed. See D4. |
| E02 — Colomés (48), thesis (4.36), and comparison with (39) | Both left and right adjoint products must have minus signs. | `colomes_VMS_paper.pdf`, (48), and `Kabalan_Ali_Thesis.pdf`, (4.36), print plus signs with the negative adjoint. Colomés (39) does use minus signs for the time/residual adjoint contributions. | Agree | D2 independently confirms the correction. The other dynamic mass-coupling terms in (48)/(4.36) have the correct signs. Equation (39) is only a sign cross-check: a general forced version of that energy-transfer identity also needs its force contribution unless it vanishes. |
| E03 — Colomés (51), thesis (4.38) | Wrong adjoint signs and an unbalanced FE time term; remove \(D_tu_h\) from the OSS residual and projection. | `colomes_VMS_paper.pdf`, (51), has \(u_h^n/\delta t\) on the stabilization LHS, no matching \(u_h^{n-1}/\delta t\) on its RHS, and defines the projection from the full residual in (53). `Kabalan_Ali_Thesis.pdf`, (4.38), copies this. | Agree, with an alternative correction | With the original full-residual projection, restore the missing FE history and correct both adjoint signs. Equivalently use the spec's reduced projection and remove the entire FE time derivative. Merely deleting the current-time term while retaining the original projection is wrong. See D3–D4. |
| E04 — Colomés (50), (52); thesis (4.39) | Missing \(\tau_t\) multiplying \(\tilde u\cdot\nabla u_h\) on the left. | `colomes_VMS_paper.pdf`, (50), (52), and `Kabalan_Ali_Thesis.pdf`, (4.39), indeed omit it. | Agree, but thesis correction incomplete | The left side is \((I+\tau_t\nabla u_h)\tilde u\). Thesis (4.39) additionally nests the projection factor inside another \(\tau_t\); see T08 and D8. This additional error is not in Colomés (52). |
| E05 — chapter (64) | Missing \(\partial_tu_h\) in the residual. | `chapweb09.pdf`, (58), contains \((\partial_tu_h,v')\); (64) omits it. Also (57) contains \((\partial_tu',v_h)\), while (63) omits it. | Agree for general dynamic ASGS; incomplete audit | **Both** cross-scale time terms are needed for the general dynamic ASGS model. Both vanish for exact OSS in fixed spaces. The chapter's intent at this transition is not explicit enough to decide whether it silently specialized to OSS or accidentally omitted them. It cannot be used as a generic ASGS pair as printed. |
| E06 — thesis (4.9) | Groups (4.9), (4.11), (4.15) as printing \(f-\partial_tu_h+L_a(u_h,p_h)\). | `Kabalan_Ali_Thesis.pdf`, (4.9), actually prints \(f-(\partial_tu_h+a\cdot\nabla u_h)+(-\nu\Delta_hu_h+\nabla p_h)\) in weak form. | Disagree with exact diagnosis; agree that it is wrong | Its convection sign is already correct; its viscous and pressure signs are wrong. The spec's common replacement residual is correct. See D4. |
| E07 — thesis (4.11), (4.15) | The whole spatial residual sign is reversed. | `Kabalan_Ali_Thesis.pdf`, (4.11), (4.15), print \(f-\partial_tu_h+a\cdot\nabla u_h-\nu\Delta_hu_h+\nabla p_h\). | Agree | These do reverse all spatial terms. Equation (4.30) later uses the corrected sign without reconciling it with (4.15). |
| E08 — thesis (4.28)–(4.29) | BE subscales plus BDF2 velocity is an error; the correct form must use the same BDF. | `Kabalan_Ali_Thesis.pdf`, (4.28), is BE; (4.29) is BDF2. `time_dependent_2007.pdf`, §3.1, explicitly allows different integrators after explaining why matching schemes commute. | Qualified / Disagree with “must” | Matching BDF is the spec's sensible declared choice, not the only mathematically admissible method. BE subscales give a generic first-order temporal defect in the coupled system, so a general second-order claim is unsupported. The thesis also later prints BE for the FE velocity itself in (4.36)–(4.39), conflicting with its stated BDF2 choice. See D5. |

The split method needs its own audit because it is not merely a ninth switch setting of the residual-based master system.

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| S01 — split momentum/pressure couplings | Couple \(\tilde u_1\) only to \(u_h\cdot\nabla v_h\), and \(\tilde u_2\) only to \(\nabla q_h\), both with minus signs; use \(a=u_h\). | `term_by_term.pdf`, (19)–(22), has those two negative FE couplings and the negative projected drivers. Its final FE equation (21) uses \(u_h\). | Agree | Independent BDF discretization gives the two subscale expressions in the spec, for the chosen projection. See D7. |
| S02 — split parameter formula | \(\tau_1^{-1}=c_1\nu/(h_1/k)^2+c_2\lvert u_h\rvert/(h_2/k)\), \(\tau_2=(h_1/k)^2/(c_1\tau_1)\). | `term_by_term.pdf`, (23), actually uses **\(h_1/k^2\)** in both viscous and pressure length factors, but \(h_2/k\) in convection. Equation (30) confirms \(k^4\) viscous scaling. | Disagree | Verified in the original page image, printed p. 707. The error is hidden at k=1 and changes both parameters for k>1. |
| S03 — split pressure projection | Reuses the master \(\zeta\), which is \(\tau_c\)-weighted. | `term_by_term.pdf`, (22), uses the same plain orthogonal projector notation for \(\nabla\cdot u_h\); no separate weighted pressure projection is specified. | Disagree with faithful-reproduction claim | Weighted pressure projection is a defensible modification, but must be marked as ours. For spatially variable \(\tau_c\), it is not identical to the printed split method. |
| S04 — split projection spaces | Plain L2 projection onto unconstrained \(\hat V_h\), “as published.” | `term_by_term.pdf`, §§2.3–2.4, defines \(V_h\subset H_0^1\) and \(Q_h\subset L^2/\mathbb R\); (18)–(22) uses the projector onto the FE space. It does not explicitly introduce unconstrained projection spaces. | Unsupported / Ambiguous | Unconstrained projection is supported by `onseen_2008.pdf`, not explicitly specified in the 2019 paper. The intended boundary treatment in its numerical implementation cannot be inferred from these equations. Label this as our choice rather than “as published.” |
| S05 — variable-weight split orthogonality | Plain projection of each driving term, not history; subscales need not remain orthogonal when \(\tau_t\) varies. | `term_by_term.pdf`, (19)–(20), (24)–(25), prints this structure; the text after (14) acknowledges that omitting a left-side projection is an approximation for variable coefficients. | Agree | The spec correctly identifies this limitation of the printed pointwise closure. Elementwise constant but element-dependent weights can also break global orthogonality. |
| S06 — split energy | The energy balance has no cross terms between FE and subscales. | `term_by_term.pdf`, §4.1, (27)–(29), inserts projections into the FE test operators using orthogonality and takes parameters constant in the patch. Its §4.3 stability analysis uses constant parameters. | Disagree as an unconditional claim | For the spec's unprojected FE test operators and plain projections with variable \(\tau_t\), nonzero projection-transfer terms remain. D7 displays them. The spec cannot simultaneously assert loss of orthogonality and use the orthogonality-based cancellation without qualification. |
| S07 — split time scheme | Uses a common BDF, while the publication gives BE subscales. | `term_by_term.pdf`, (24)–(25), uses BE; (8)–(9) permits FE BDF1/BDF2. | Agree | Clearly marked as our decision in the spec. A common BDF extension follows directly from (19)–(20). |
| S08 — consistency | Not residual based, with consistency error \(O(h^{k+1})\). | `term_by_term.pdf`, §3.2 following (18): this estimate is for sufficiently smooth projected quantities; it discusses using powers of \(\tau_1\) if that smoothness fails. | Qualified | The source supports optimal consistency order under regularity assumptions, not an unconditional \(h^{k+1}\) estimate for arbitrary data. |
| S09 — parameters fixed for comparison | Either use the split parameters everywhere or the master parameters in split OSS. | `term_by_term.pdf`, (23), and `colomes_VMS_paper.pdf`, (20)–(21), define different geometric conventions. | Ambiguous / unresolved decision | Both comparison conventions can be specified, but the spec currently leaves this choice open despite its opening assertion that no further choices remain. No option is selected in this verification. |

These are additional Chapter 4 findings not adequately captured by the original six-row errata table. “Not mentioned” in the second column is intentional.

| Item | Spec says | Source says (file, eq. no.) | Agree / Disagree | Note |
| --- | --- | --- | --- | --- |
| T01 — test-function indices | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, (4.1), uses \(q_h\) in the continuous weak form; (4.2) uses \(v\) in the subscale contribution while the rest is tested by \(v_h\). | Additional notation errors | Use \(q\) in (4.1), \(v_h\) in the final argument of (4.2). Compare `colomes_VMS_paper.pdf`, (7), (13). |
| T02 — bulk integration by parts | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, (4.3), has \(+(a\cdot\nabla v_h,\tilde u)\), inconsistent with its negative adjoint in (4.5). | Additional error | The sign must be negative. A \(-(\tilde u,(\nabla\cdot a)v_h)\) term also occurs unless solenoidality or its omission is explicitly assumed. See D8. |
| T03 — coarse-equation boundary terms | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, (4.3), prints a viscous boundary expression involving \(\nu\nabla^2\tilde u\cdot n\) and a negative continuity boundary term. | Additional errors | Integration of the weak viscous term onto the test function gives \(\nu\langle\tilde u,\partial_n v_h\rangle\), not the printed expression; continuity gives \(+\langle q_h,\tilde u\cdot n\rangle\). See D8. |
| T04 — fine-equation boundary terms | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, (4.8), adds \(+pn\) to the viscous traction, includes a convective boundary term although its volume convection was not integrated by parts, and uses \(\tilde u\) rather than the test function in that boundary term. | Additional errors | The boundary term from (4.6) is \(\langle\nu\partial_n(u_h+\tilde u)-(p_h+\tilde p)n,\tilde v\rangle\). There is no extra convective boundary term in that representation. Bubble test functions can remove the correct boundary term, but do not make the printed identity correct. |
| T05 — continuity split | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, (4.10), has \((\nabla\cdot\tilde u,q_s)=-(\nabla\cdot\tilde u,\tilde q)=0\). | Additional error | It must be \((\nabla\cdot\tilde u,\tilde q)=-(\nabla\cdot u_h,\tilde q)\); neither side is generally zero. Its own (4.7) and (4.12) establish this. |
| T06 — exactness of the projected strong subscale equation | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, (4.11)–(4.18), moves from projected differential equations to an unprojected operator with \(\xi_h=\Pi R_m\), claiming the RHS must belong to the subscale space. | Additional modelling gap / Ambiguous | The differential operator need not preserve that space. Projected equality alone does not establish this particular multiplier; D8 gives the missing term. It can be a closure approximation, but is not an exact rearrangement. |
| T07 — mutually inconsistent time specifications | Mentions the BE/BDF2 mixture only. | `Kabalan_Ali_Thesis.pdf`, (4.29), declares FE BDF2, whereas (4.36), (4.37), (4.38), and (4.39) explicitly use FE BE current/history coefficients. | Additional inconsistency | These displayed equations do not implement even the stated mixed scheme. (4.31) and (4.35) are valid BE subscale formulas, but would also change for common BDF2. The actual scheme used in software remains unknown. |
| T08 — extra RHS factor in (4.39) | Mentions only the missing factor on the left. | `Kabalan_Ali_Thesis.pdf`, (4.39), places \(-\tau_{m,t}^{k-1}\xi_h^{n,i-1}\) **inside** a parenthesis already multiplied by \(\tau_{m,t}^{k-1}\). `colomes_VMS_paper.pdf`, (52), places it outside that bracket. | Additional error | The thesis has an unintended \(-\tau_t^2\xi_h\). Correct it to one factor, independently of restoring the left-side factor. Confirmed visually on printed thesis p. 29. |
| T09 — projection spaces and notation | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, text before (4.40)–(4.41), quantifies over \(V_0,Q_0\), not FE spaces, and changes \(\xi_h,R^m,R^c\) to \(\xi_m,R_u,R_p\) without explaining the aliases. | Additional notation / definition errors | Literally enforcing the projection for every continuous test function is not the stated finite-dimensional projection. FE spaces must be identified; the constrained/unconstrained choice must then be stated. See `colomes_VMS_paper.pdf`, (53)–(54), and `onseen_2008.pdf`, (18). |
| T10 — conservative tensor expansion | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, (4.43), has positive \((u\otimes u,\nabla v)\) on the left and the negative four-term expansion on the right. | Additional error | With \(v=v_h\), expanding the tensor product has all positive signs. If the RHS is retained, the LHS must instead be the negative tensor contraction or the appropriate conservative convective form. Compare `colomes_VMS_paper.pdf`, (30), and `time_dependent_2007.pdf`, (31)–(35). |
| T11 — explicit subscale after model selection | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, text after (4.34), calls the selected subscale an explicit expression, while (4.25), (4.21), (4.37) make it nonlinear and implicit in nonlinear splitting. | Additional overstatement | Explicit elimination in terms of only \(u_h,p_h\) is not supplied for nonlinear splitting. Later iterative treatment acknowledges this. |
| T12 — residual sign changes later | Mentions (4.15), not its contradiction with (4.30). | `Kabalan_Ali_Thesis.pdf`, (4.30), equates \(R^m\) with \(f-\partial_tu_h-L\), contrary to its definition in (4.15). | Additional manifestation of E07 | The second line of (4.30) has the correct spatial signs. This is not a separate new residual-sign mechanism, but the two definitions cannot both be used. |
| T13 — stress definitions around (4.42) | Not mentioned. | `Kabalan_Ali_Thesis.pdf`, §4.6 and (4.42), uses the usual filtered tensor difference with nested overbars in the Leonard/cross/Reynolds decomposition. | No further error found here | Checked the original image because text extraction loses overbars. Do not infer an error from the flattened `pdftotext` output. |

The independent derivation follows. The labels D1–D8 provide the algebra behind the table and identify exactly where the sources or spec disagree.

**D1. Start with the semi-discrete system and distinguish the source models.**

Write \(s=\tilde u\), \(d=\nabla\cdot u_h\), and let \((\cdot,\cdot)_h\) denote the sum of element integrals. Define

\[
L_a(u,p)=a\cdot\nabla u-\nu\Delta_hu+\nabla p,
\qquad
A_a(v,q)=-a\cdot\nabla v-\nu\Delta_hv-\nabla q.
\]

The letter \(A_a\) here is the spec's \(L_a^*\). The convection part of the exact formal adjoint would also contain \(-(\nabla\cdot a)v\). Discarding that contribution and the element-boundary terms gives the adopted operator above.

Let

\[
B_a(u_h,p_h;v_h,q_h)
=\nu(\nabla u_h,\nabla v_h)+c(a;u_h,v_h)
 -(p_h,\nabla\cdot v_h)+(q_h,d).
\]

Colomés (13)–(17) then gives, before any time discretization,

\[
(\partial_tu_h,v_h)+B_a
+\chi_D(\partial_ts,v_h)+(s,A_a(v_h,q_h))_h
-(\tilde p,\nabla\cdot v_h)=(f,v_h),
\tag{D1a}
\]

\[
\chi_D\partial_ts+\tau_m^{-1}s
=f-\partial_tu_h-L_a(u_h,p_h)-\lambda_h,
\tag{D1b}
\]

where \(\lambda_h=0\) for ASGS; for exact OSS it is the FE correction determined together with the orthogonality constraint. The static choice sets \(\chi_D=0\) in **both** equations, even though the resulting algebraic subscale may still change with time. Taking

\[
\tilde p=-\tau_c(d-\zeta_h)
\tag{D1c}
\]

gives the spec's positive divergence-stabilization term. This establishes the basic signs without looking at Colomés (38), (48), or (51).

The 2007 equations give the same velocity structure by a direct regrouping. Their (17) contains

\[
(u_h\cdot\nabla u_h,v_h)+(s\cdot\nabla u_h,v_h)
-(s,s\cdot\nabla v_h)_h
-(s,u_h\cdot\nabla v_h+\nu\Delta_hv_h+\nabla q_h)_h.
\]

The first two terms combine to \((a\cdot\nabla u_h,v_h)\), and the last two combine to \((s,A_a(v_h,q_h))_h\), with \(a=u_h+s\). Its (15) is exactly \(P[f-\partial_tu_h-L_a(u_h,p_h)]\). Replacing its nonconservative Galerkin convection by the spec's type-1 form is the declared Colomés-based modelling choice, not a time-discretization identity.

For pressure, retaining the full 2007 equation (12) and applying a common discrete derivative would give

\[
\tilde p^{n+1}=\tau_c^{n+1}
\left(r_p^{n+1}+\tau_m^{n+1}D_tr_p^{n+1}\right),
\qquad r_p=-P(\nabla\cdot u_h).
\tag{D1d}
\]

That is different from D1c. The spec explicitly discards this derivative and adopts the later pressure closure. The 2007 equation (17) itself instead sets \(\tilde p=0\). These are three distinct pressure choices; none should be silently conflated.

Starting independently from chapter (57)–(58), neglecting the stated boundary terms and applying (59) to the fine-scale spatial operator gives

\[
(\partial_tu_h,v_h)+B_a+(\partial_ts,v_h)
+(s,A_a(v_h,q_h))_h=(f,v_h),
\]

\[
(\partial_ts,w')+(\tau_m^{-1}s,w')_h
=(f-\partial_tu_h-L_a(u_h,p_h),w')_h.
\tag{D1e}
\]

Here the chapter has set pressure subscales to zero. Both time cross terms survive this derivation. For OSS, \((\partial_ts,v_h)=0\) and \((\partial_tu_h,w')=0\) in fixed mutually orthogonal spaces, which explains a possible OSS reading of chapter (63)–(64). For general ASGS, (63) needs \((\partial_ts,v_h)\) and (64) needs \((\partial_tu_h,w')\). The spec identifies only the second omission.

**D2. Apply BDF before eliminating the subscale.**

At the new time level define

\[
\gamma=\frac{\alpha_0}{\delta t},\qquad
H_u=\frac{u_{h,\mathrm{hist}}}{\delta t},\qquad
H_s=\frac{s_{\mathrm{hist}}}{\delta t}.
\]

Thus \(D_tu_h=\gamma u_h-H_u\) and \(D_ts=\gamma s-H_s\). The correct uniform-step coefficients are

\[
\begin{array}{c|c|c}
 &\alpha_0&g_{\mathrm{hist}}\\\hline
\mathrm{BDF1}&1&g^n\\
\mathrm{BDF2}&3/2&2g^n-\tfrac12g^{n-1}.
\end{array}
\]

Directly substituting these definitions into D1b gives

\[
(\chi_D\gamma+\tau_m^{-1})s
=f-(\gamma u_h-H_u)-L_a(u_h,p_h)
+\chi_DH_s-\lambda_h.
\]

Consequently,

\[
\boxed{\tau_t=(\chi_D\gamma+\tau_m^{-1})^{-1},\qquad
s=\tau_t\bigl(f-D_tu_h-L_a(u_h,p_h)+\chi_DH_s-\lambda_h\bigr).}
\tag{D2a}
\]

The corresponding FE equation is

\[
\boxed{
(D_tu_h,v_h)+B_a+\chi_D(D_ts,v_h)
+(s,A_a(v_h,q_h))_h
+(\tau_c(d-\zeta_h),\nabla\cdot v_h)_h=(f,v_h).
}
\tag{D2b}
\]

No additional temporal test operator or unexplained change of adjoint sign occurs. For OSS, D3 justifies removing the subscale mass coupling; for ASGS it remains. Together D2a–D2b reproduce the spec's master system and the eight equations in its table.

For an explicit check of the printed fully discrete ASGS equation, take dynamic subscales, set

\[
Z=\gamma u_h+L_a(u_h,p_h),\qquad
Y=f+H_u+H_s,\qquad s=\tau_t(Y-Z),
\]

and substitute in D2b. With \(G_c=(\tau_cd,\nabla\cdot v_h)_h\), the result is

\[
\begin{aligned}
\gamma(u_h,v_h)+B_a+G_c
&-(\tau_tZ,A_a(v_h,q_h))_h
-\gamma(\tau_tZ,v_h)\\
={}&(f,v_h)+(H_u,v_h)+(H_s,v_h)\\
&-(\tau_tY,A_a(v_h,q_h))_h
-\gamma(\tau_tY,v_h).
\end{aligned}
\tag{D2c}
\]

For BE, \(\gamma=1/\delta t\), \(H_u=u_h^n/\delta t\), \(H_s=s^n/\delta t\). Comparing D2c with Colomés (48) and thesis (4.36) proves that **both adjoint products have the wrong printed sign**, while the additional subscale mass-coupling signs are correct. BDF2 changes all of \(\gamma,H_u,H_s,\tau_t\); changing just the Galerkin mass coefficient would not give D2c.

As a separate check against the 2007 semi-discrete equations, applying a theta rule to D1b gives, with \(s^{n+\theta}=\theta s^{n+1}+(1-\theta)s^n\),

\[
\left(\frac1{\theta\delta t}+\tau_m^{-1}\right)s^{n+\theta}
=r_u^{n+\theta}+\frac{s^n}{\theta\delta t}.
\tag{D2d}
\]

This is obtained from \((s^{n+1}-s^n)/\delta t\), not by copying (25)–(26). It checks the effective-parameter mechanism independently of BDF. The FE residual's time derivative must use the corresponding coarse-scale time discretization.

**D3. Derive the OSS weight, history projection, and reduced residual.**

Let \(\Pi_w\) denote the projection onto the chosen unconstrained FE space with positive weight \(w\):

\[
(w\Pi_wX,z_h)=(wX,z_h)\quad\forall z_h\in\hat V_h.
\]

For the full-residual variable \(\lambda_h\), impose ordinary L2 orthogonality on the subscale obtained from D2a:

\[
0=(s,z_h)=
\bigl(\tau_t[R_u+\chi_DH_s-\lambda_h],z_h\bigr).
\]

Therefore the required equation is

\[
\lambda_h=\Pi_{\tau_t}(R_u+\chi_DH_s).
\tag{D3a}
\]

This proves why the weight must be \(\tau_t\), and why the old subscale must be **inside** the projected quantity. Even when \((H_s,z_h)=0\), it need not follow that \((\tau_tH_s,z_h)=0\). A spatially constant weight is a special case where the history projection vanishes.

On the fixed mesh, \(D_tu_h\in\hat V_h\), including a representable nonhomogeneous boundary lift. A projection reproduces every function in its range, so

\[
\lambda_h
=\Pi_{\tau_t}(\hat R_u+\chi_DH_s-D_tu_h)
=\xi_h-D_tu_h,
\quad
\xi_h:=\Pi_{\tau_t}(\hat R_u+\chi_DH_s).
\tag{D3b}
\]

Substituting this identity into D2a yields

\[
\boxed{s=\tau_t(\hat R_u+\chi_DH_s-\xi_h),\qquad
(\tau_t\xi_h,z_h)=(\tau_t[\hat R_u+\chi_DH_s],z_h).}
\tag{D3c}
\]

This is precisely the spec's reduced formulation. Its \(\xi_h\) is **not numerically the same field** as Colomés' full-residual \(\xi_m\). They differ by \(D_tu_h\). Lagged full- and reduced-projection iterations may consequently differ before convergence, even when their correctly assembled fixed points agree.

For each stored subscale that is orthogonal to the same FE space,

\[
(D_ts^{n+1},v_h)
=\gamma(s^{n+1},v_h)-(H_s,v_h)=0.
\tag{D3d}
\]

Zero initialization establishes the initial case. This proves the spec's \(\chi_A\chi_D\) switch for the residual-based OSS variants. The cancellation is not guaranteed for a subscale computed with a lagged, unconverged projection or with a history transferred to a different mesh without preserving orthogonality.

For pressure, independently impose \((\tilde p,r_h)=0\). With D1c this gives

\[
(\tau_c\zeta_h,r_h)=(\tau_cd,r_h)\qquad\forall r_h\in\hat Q_h.
\tag{D3e}
\]

The paper projects \(R_p=-d\), hence \(\xi_c=-\zeta_h\); there is no pressure-sign contradiction. If \(\tau_c=0\) everywhere, the pressure subscale and this projection equation are simply omitted, as the spec says. It would not be a positive-weight projection problem in that case.

At a steady discrete state, BDF history gives \(H_s=\gamma s\), \(D_tu_h=0\). D2a then reduces to \(\tau_m^{-1}s=\hat R_u-\xi_h\). Combined with \((s,z_h)=0\), this implies \(\xi_h=\Pi_{\tau_m}\hat R_u\), even though the dynamic solve originally used \(\tau_t\). Thus the steady-limit statement is correct for the weighted OSS model too.

One qualification to the notation \(P(R_u)\) matters: for variable coefficients, a fixed plain projector on the right of a pointwise ODE does not by itself preserve the subscale space. D3a is the constraint-compatible correction obtained after the time discretization. It is the weighted modification explicitly motivated by Colomés, not a consequence of simply keeping an unchanged plain projector in (15).

**D4. Resolve the adjoint signs, the OSS history omission, and the residual errata.**

For static subscales, set \(s=\tau_mP(f-\partial_tu_h-L_au_h)\) in D1a and test with \((v_h,q_h)=(u_h,p_h)\). Writing the pressure contribution as \(D_c\), direct expansion gives

\[
\begin{aligned}
\frac12\partial_t\|u_h\|^2+\nu\|\nabla u_h\|^2+c(a;u_h,u_h)
&-(\tau_mP\partial_tu_h,A_a(u_h,p_h))_h\\
&-(\tau_mPL_a(u_h,p_h),A_a(u_h,p_h))_h+D_c\\
={}&(f,u_h)-(\tau_mPf,A_a(u_h,p_h))_h.
\end{aligned}
\tag{D4a}
\]

This corrects Colomés (38). The negative adjoint has not changed definition: moving the residual's spatial and temporal terms to the left produces the minus signs. With \(\nu\Delta_hv_h=0\), \(A_a=-(a\cdot\nabla v_h+\nabla q_h)\), so these minus signs give the usual positive streamline/pressure stabilization. If the printed plus signs were retained, the corresponding quadratic stabilization would have the opposite sign. At higher order the viscous test-Laplacian contribution remains, explaining qualification M06.

For dynamic subscales, substitution of \(s=\tau_m(Pf-P\partial_tu_h-PL_au_h-\partial_ts)\) into the FE-to-subscale transfer gives

\[
\begin{aligned}
\varepsilon_h={}&(\partial_ts,u_h)
-(\tau_m\partial_ts,A_a(u_h,p_h))_h
-(\tau_mP\partial_tu_h,A_a(u_h,p_h))_h\\
&-(\tau_mPL_a(u_h,p_h),A_a(u_h,p_h))_h
+(\tau_mPf,A_a(u_h,p_h))_h+D_c.
\end{aligned}
\tag{D4b}
\]

The signs of the nonforce terms agree with Colomés (39). Its displayed formula omits the force product, which is harmless only when that product vanishes, for example for zero force or OSS with force in the projection space. It should not be used as a complete general forced identity.

For dynamic OSS, substituting D3c in D2b and eliminating the time coupling gives the reduced equation

\[
\begin{aligned}
\gamma(u_h,v_h)+B_a
&-(\tau_tL_a(u_h,p_h),A_a(v_h,q_h))_h
+(\tau_c(d-\zeta_h),\nabla\cdot v_h)_h\\
={}&(f,v_h)+(H_u,v_h)
-(\tau_t[f+H_s-\xi_h],A_a(v_h,q_h))_h.
\end{aligned}
\tag{D4c}
\]

Alternatively, retaining the original full-residual projection \(\lambda_h\) gives the same equation with \(L_a\) on the left replaced by \(\gamma u_h+L_a\), and the right stabilization bracket replaced by

\[
f+H_u+H_s-\lambda_h.
\tag{D4d}
\]

Colomés (51) and thesis (4.38) keep the current \(\gamma u_h\) contribution but omit \(H_u\) from this bracket. Their own projection definition does not justify that omission. Both valid representations D4c and D4d are available; the spec chooses D4c. This choice requires the projection redefinition D3b as well as the corrected signs.

For the thesis residual, moving the coarse terms in (4.6) to the other side gives, without any projection or discretization argument,

\[
f-\partial_tu_h-a\cdot\nabla u_h+\nu\Delta_hu_h-\nabla p_h.
\tag{D4e}
\]

In (4.9) the \(-a\cdot\nabla u_h\) term already has this sign; only the last two spatial terms are wrong. In (4.11) and (4.15) all three spatial terms have the wrong signs. The spec's final replacement is right, but its common “as printed” description is not.

For nonlinear splitting,

\[
L_{u_h+s}(u_h,p_h)=L_{u_h}(u_h,p_h)+(\nabla u_h)s,
\]

where \((\nabla u_h)_{ij}=\partial_j u_{h,i}\). Moving this last term across D2a gives, with \(\tau_t\) frozen during the local update,

\[
\boxed{
(I+\tau_t\nabla u_h)s
=\tau_t\left(f-D_tu_h-L_{u_h}(u_h,p_h)+\chi_DH_s-\lambda_h\right).
}
\tag{D4f}
\]

For reduced OSS replace \(-D_tu_h-\lambda_h\) by \(-\xi_h\). This derives the missing factor in Colomés (50), (52) and thesis (4.39). Without it the left side adds a velocity to an acceleration, so the printed equation also fails a dimensional check.

**D5. What can and cannot be concluded about the thesis time schemes.**

For a smooth subscale, Taylor expansion at the new time gives

\[
\frac{s^{n+1}-s^n}{\delta t}
=s'(t^{n+1})-\frac{\delta t}{2}s''(t^{n+1})+O(\delta t^2),
\]

whereas

\[
\frac{3u_h^{n+1}-4u_h^n+u_h^{n-1}}{2\delta t}
=u_h'(t^{n+1})+O(\delta t^2).
\tag{D5a}
\]

Thus BE in the subscale equation leaves a first-order temporal defect in the coupled system. BDF2 in the coarse equation does not remove it. Under the usual stability assumptions, first-order coupled accuracy is the generic expectation at fixed mesh, although special solutions, vanishing subscales, or other cancellations can make the observed error smaller. The unqualified claim that every solution necessarily loses second order is too strong.

There is also no single common discrete derivative satisfying

\[
D(u_h+s)=D_{\mathrm{BDF2}}u_h+D_{\mathrm{BE}}s
\]

for arbitrary independent histories if \(D\) is chosen as either BDF1 or BDF2 on the unsplit velocity. This is why the simple common-scheme space–time commutation argument does not apply to the mixed method. It does not make every mixed-scale integrator inadmissible; 2007 §3.1 explicitly discusses using different schemes.

The separate and definite thesis inconsistency is that (4.36)–(4.39) do not use its declared BDF2 derivative for the coarse velocity at all. They contain only the current and one old value with coefficients \(1,-1\), while BDF2 requires \(3/2,-2,1/2\). The manuscript alone cannot resolve which choice was actually used.

**D6. Boundary energy and limits of the fixed-choice attributions.**

For the spec's type-1 convection including its proposed boundary term,

\[
c(a;u_h,u_h)
=\frac12\int_{\Gamma_N}(u_h\cdot n)|u_h|^2.
\tag{D6a}
\]

The two volume terms cancel identically, including with discontinuous \(a\). The whole expression is zero only when the boundary flux is zero or cancels. It represents advective energy transport through an open boundary. This is not a defect in the chosen expression, but it contradicts calling the entire extended form “exactly energy-neutral” without the boundary qualification.

Likewise, integration of the chosen viscous and pressure volume terms produces

\[
\nu(\nabla u,\nabla v)-(p,\nabla\cdot v)
=(-\nu\Delta u+\nabla p,v)
+\langle\nu\partial_nu-pn,v\rangle_{\partial\Omega}.
\tag{D6b}
\]

This identifies the gradient-form natural traction. A prescribed symmetric-stress traction \([\nu(\nabla u+\nabla u^T)-pI]n\) is not interchangeable with it without the corresponding weak-form/boundary treatment. Colomés avoids this issue in its derivation by taking the whole boundary Dirichlet. The spec adds \(\Gamma_N\) but does not say whether the missing traction load is intentionally zero, so this remains an ambiguity, not a selected answer.

The other attribution issues have no algebraic ambiguity: Table 1 lists six cases, the Oseen paper expressly allows unequal order, Colomés §4 says Picard for the FE solve, and its channel calculations use c2=8. Selecting eight cases, equal-order P elements, an inner Newton solve, and c2=2 as a common default is possible; those selections must be identified as ours where they go beyond the cited statements. Pointwise parameter evaluation and common BDF/startup choices likewise need that distinction.

**D7. Re-derive split OSS and check its energy claim rather than importing it.**

Let \(\Pi_0\) be the plain L2 projector onto the selected FE velocity space, \(C_1=u_h\cdot\nabla u_h\), and \(C_2=\nabla p_h\). Starting with the semi-discrete equations (19)–(20) in `term_by_term.pdf`, applying the chosen BDF gives

\[
(\gamma+\tau_1^{-1})s_i=H_i-C_i+\Pi_0C_i,
\quad H_i=\frac{s_{i,\mathrm{hist}}}{\delta t},
\]

and hence

\[
\boxed{s_i=\tau_{1,t}(H_i-C_i+\Pi_0C_i),
\qquad \tau_{1,t}=(\gamma+\tau_1^{-1})^{-1}.}
\tag{D7a}
\]

This confirms the spec's signs and history coefficients for its declared plain-projector version. It does not imply orthogonality:

\[
(s_i,z_h)=(\tau_{1,t}H_i,z_h)
-(\tau_{1,t}(I-\Pi_0)C_i,z_h),
\tag{D7b}
\]

and neither term generally vanishes with a spatially varying weight. Even a weight that is constant on each element need not commute with a global FE projection.

To test the energy assertion, retain exactly the spec's FE equation with its **unprojected** test operators. Test it by \((u_h,p_h)\) and add each discrete subscale equation tested by \(s_i\). Let \(p_s=-\tau_c(d-\zeta)\) be the spec's weighted pressure closure, so that \((p_s,\zeta)=0\). Then direct cancellation gives

\[
\begin{aligned}
(D_tu_h,u_h)+\sum_{i=1}^2(D_ts_i,s_i)
&+\nu\|\nabla u_h\|^2+c(u_h;u_h,u_h)\\
&+\sum_{i=1}^2(\tau_1^{-1}s_i,s_i)
+(\tau_c^{-1}p_s,p_s)\\
={}&(f,u_h)+\sum_{i=1}^2(\Pi_0C_i,s_i).
\end{aligned}
\tag{D7c}
\]

Omit the pressure term if \(\tau_c=0\). This identity is already fully discrete; the \((D_tg,g)\) terms can further be written as BDF energy differences and numerical dissipation. Their particular representation cannot cancel the extra spatial projection products on the right.

If the velocity subscales are truly orthogonal, the two extra terms vanish. For the variable-weight plain-projection closure, that premise is unavailable. Section 4.1 of the paper rewrites the FE test operators with \(P_h^\perp\), using orthogonality, and assumes constant parameters in its displayed patch calculation. The unprojected and projected test-operator equations are not automatically equivalent under the variable-weight approximation. This explains the conflict within the spec; no choice between those distinct formulations is made here. Nor does global orthogonality alone justify orthogonality on every arbitrary element patch.

The pressure part also differs from a literal reading of the paper: (22) gives \(p_s=-\tau_2(I-\Pi_0)d\), while the spec has a weighted scalar projection. They coincide for a globally constant weight, but not in general. The spec's pressure modification does make its pressure contribution in D7c dissipative; that does not eliminate the two velocity-projection terms.

Finally, the actual parameter formula printed in (23) is

\[
\boxed{
\tau_1=\left(\frac{c_1\nu}{(h_1/k^2)^2}
+\frac{c_2|u_h|}{h_2/k}\right)^{-1},\qquad
\tau_2=\frac{(h_1/k^2)^2}{c_1\tau_1}.
}
\tag{D7d}
\]

Equivalently, the viscous contribution to \(\tau_1^{-1}\) is \(c_1k^4\nu/h_1^2\). The spec instead gives \(c_1k^2\nu/h_1^2\). At k=2 these coefficients differ by a factor of four. The pressure numerator also differs by a factor of four before accounting for the changed \(\tau_1\). With \(h_1=h_2=h\), D7d matches the polynomial scaling of the master parameters when that scaling is applied throughout.

The consistency statement requires the same qualifications as its source: for a sufficiently smooth quantity \(g\), \(\|(I-\Pi_0)g\|=O(h^{k+1})\). It is not an unconditional statement about every convection, pressure-gradient, or force field regardless of smoothness. The source explicitly discusses stabilization-weight compensation when fewer derivatives are available.

**D8. Derivations of the additional thesis Chapter 4 errors.**

The integration-by-parts identities needed in thesis (4.3), on an element K, are

\[
\begin{aligned}
(a\cdot\nabla s,v_h)_K
&=-(s,a\cdot\nabla v_h)_K
 -(s,(\nabla\cdot a)v_h)_K
 +\langle(a\cdot n)s,v_h\rangle_{\partial K},\\
\nu(\nabla s,\nabla v_h)_K
&=-\nu(s,\Delta v_h)_K
 +\nu\langle s,\partial_n v_h\rangle_{\partial K},\\
(q_h,\nabla\cdot s)_K
&=-(\nabla q_h,s)_K+\langle q_h,s\cdot n\rangle_{\partial K}.
\end{aligned}
\tag{D8a}
\]

These establish the wrong convective volume sign, the wrong viscous boundary expression, and the wrong continuity boundary sign in (4.3). The corrected volume expression agrees with (4.4)–(4.5) after the declared divergence/boundary approximations. Thus (4.4) is not obtained from the printed (4.3) by merely deleting boundary terms.

For (4.8), let \(U=u_h+s\), \(P=p_h+\tilde p\). Starting from (4.6) and integrating just viscosity and pressure gives

\[
(\partial_tU+a\cdot\nabla U-\nu\Delta U+\nabla P,\tilde v)_K
+\langle\nu\partial_n U-Pn,\tilde v\rangle_{\partial K}
=(f,\tilde v)_K.
\tag{D8b}
\]

The pressure boundary sign is negative. There is no added convection boundary integral while \((a\cdot\nabla U,\tilde v)\) remains unchanged in the volume, and any retained term must be linear in \(\tilde v\), not replace the test function by \(s\).

For (4.10), splitting the original continuity equation (4.7) gives immediately

\[
(\nabla\cdot s,\tilde q)=-(\nabla\cdot u_h,\tilde q).
\tag{D8c}
\]

The printed replacement of \(u_h\) by \(s\) on the right, the new undefined \(q_s\), and the final equality to zero are unjustified. The thesis's (4.12) has the correct projected continuity relation.

The projected-operator issue in (4.13)–(4.18) can be isolated without relying on any residual sign. Set

\[
S=\partial_ts+a\cdot\nabla s-\nu\Delta s+\nabla\tilde p,
\quad R=f-\partial_tu_h-L_a(u_h,p_h).
\]

For an L2-orthogonal splitting, the fine equation says \((I-\Pi_0)(S-R)=0\). It implies

\[
S=R-\lambda,\qquad \lambda=\Pi_0(R-S),
\tag{D8d}
\]

not automatically \(\lambda=\Pi_0R\). The latter requires \(\Pi_0S=0\). Orthogonality of \(s\) does not imply orthogonality of its convective or viscous derivatives, or of the subscale pressure gradient. After an operator model is introduced, its compatibility with the subscale space is a separate condition, addressed for the fully discrete weighted closure in D3. The thesis presents the earlier step as exact without stating such a condition. This is a modelling/justification gap rather than a uniquely repairable sign typo.

For (4.39), inspection of the actual parentheses is essential. Its RHS has the form

\[
\tau_tH_s+\tau_t\bigl(f-D_tu_h-L_{u_h}(u_h,p_h)-\tau_t\xi_h\bigr).
\tag{D8e}
\]

D4f requires a single \(-\tau_t\xi_h\), whereas D8e has \(-\tau_t^2\xi_h\). Move the projection term outside the outer bracket or remove its inner factor. This correction is independent of the missing left-side \(\tau_t\). The paper's (52) does not have the extra RHS factor. With the original full-residual projection retained, the corrected equation is

\[
(I+\tau_t\nabla u_h)s
=\tau_t\bigl(H_s+f-D_tu_h-L_{u_h}(u_h,p_h)-\lambda_h\bigr).
\tag{D8f}
\]

For (4.40)–(4.41), the functions called projections are finite-dimensional FE fields. Quantifying their defining equations over the continuous spaces as printed would, by density, demand equality with the generally non-FE projected quantity in the weighted L2 sense, rather than a finite-dimensional best approximation. The FE test spaces and aliases for the symbols therefore need to be specified. Choosing constrained or unconstrained FE spaces is then a separate decision.

For (4.43), plain tensor expansion gives

\[
\begin{aligned}
((u_h+s)\otimes(u_h+s),\nabla v_h)
={}&(u_h\otimes u_h,\nabla v_h)
+(u_h\otimes s,\nabla v_h)\\
&+(s\otimes u_h,\nabla v_h)
+(s\otimes s,\nabla v_h).
\end{aligned}
\tag{D8g}
\]

The conservative convective weak form has the negative of both sides, with the appropriate boundary flux if present. The thesis prints the positive left side and negative right side. Its cross-stress and Reynolds-stress grouping is otherwise consistent with this expansion.

The remaining Chapter 4 issues are direct internal contradictions rather than new differential identities: the wrong test indices in (4.1)–(4.2), the BE formulas after the announced FE BDF2 choice, the corrected residual suddenly used in (4.30), and the claim that nonlinear subscales are explicit. The thesis's printed stress filters around (4.42) do not supply an additional error.

For completeness, the following register lists the disagreements and unresolved qualifications without treating a source-supported modelling choice as an algebraic error.

1. **Source adjoint signs:** Colomés (38), (48), (51) and thesis (4.36), (4.38) disagree with D2c/D4a/D4c. The spec's sign corrections agree with the derivation. Colomés (39) has the correct nonforce signs but lacks the general force-transfer term in D4b (E01–E03).
2. **Source nonlinear update:** Colomés (50), (52) and thesis (4.39) lack the left-side factor in D4f. The spec corrects this, but misses the additional squared projection factor unique to thesis (4.39), derived in D8e–D8f (E04, T08).
3. **Source OSS time imbalance:** Colomés (51) and thesis (4.38) omit coarse history in the full-residual representation. Either D4c or D4d is consistent; the printed hybrid is not. The spec's reduced formulation is valid only with its redefined projection (E03).
4. **Chapter time couplings:** Chapter (63) as well as (64) omits a term required for general dynamic ASGS. The spec's errata is incomplete. An implicit OSS specialization would explain both omissions, but the text does not establish that scope clearly (E05, D1).
5. **Thesis residuals:** (4.9), (4.11), (4.15) disagree with D4e. The spec wrongly assigns the convection error to (4.9), whose convection sign is already right. Equation (4.30) silently contradicts (4.15) (E06–E07, T12).
6. **Time-scheme claims:** The mixed BE/BDF2 method is not forbidden by the 2007 source. It has a generic first-order subscale defect and lacks the simple common-scheme commutation argument. The thesis also fails to put its announced FE BDF2 into (4.36)–(4.39). The actual implemented method is undetermined (E08, T07, D5).
7. **Split parameter exponent:** The spec's \(h_1/k\) disagrees with the printed and independently cross-checked \(h_1/k^2\) in both length factors of D7d (S02).
8. **Split energy:** The unconditional absence of FE/subscale transfer cross terms disagrees with D7c for the variable-weight plain-projection method. Constant-parameter/orthogonality assumptions of the published calculation cannot be silently removed (S05–S06).
9. **Split definition:** Weighted pressure projection and explicitly unconstrained projection spaces are not specified that way in the cited 2019 equations. They need to be declared choices; the source's actual boundary treatment remains ambiguous. The choice between two parameter conventions also remains unresolved in the spec (S03–S04, S09).
10. **Convection and boundaries:** “Exactly energy-neutral” applies to the volume skew form, not the added open-boundary flux. The chosen open-boundary traction/load convention is not specified by the cited homogeneous-boundary derivation (F03–F04, D6).
11. **Fixed-choice provenance:** Colomés Table 1 has six cases, not eight; not every source assumes equal order; the pointwise/common-length convention, BDF/startup details, and inner FE Newton architecture are selections beyond the cited statements. The channel benchmark uses c2=8, contrary to the blanket 12/2/0 attribution (F01–F02, F08–F09, F12, F14).
12. **Claims surrounding the master:** SUPG+PSPG+grad-div alone omits the higher-order viscous adjoint term; a single solve with frozen projections does not solve linear-splitting OSS to its fixed point; step 3 must obey the reduced-residual rule; the stopping variables do not by themselves verify all coupled unknowns. Exact orthogonality is relative to the actual common integration and solve accuracy. The reported performance evidence supports dynamic OSS, but does not establish the claimed unique ranking of Dyn-Nl-OSS over Dyn-Lin-OSS (M06, M08–M12).
13. **Regularity qualifications:** A retained force projection needs an L2-representable force; the split \(h^{k+1}\) estimate needs sufficient regularity. These conditions should accompany the spec's claims (F11, S08).
14. **Thesis integration by parts:** (4.3) disagrees with D8a in its convection sign and both cited boundary expressions; (4.8) disagrees with D8b in pressure traction, extra convection boundary term, and its test variable (T02–T04).
15. **Thesis continuity and projected operator:** (4.10) disagrees with D8c. The asserted exact transition in (4.11)–(4.18) does not follow without the missing condition in D8d (T05–T06).
16. **Thesis tensor algebra and definitions:** (4.43) disagrees with D8g; (4.1)–(4.2) have inconsistent test indices; (4.40)–(4.41) do not state the FE projection test spaces correctly; the blanket explicit-subscale assertion conflicts with nonlinear splitting (T01, T09–T11).
17. **Errata scope:** Colomés (38) and chapter (63)–(64) are semi-discrete statements, not printed fully discrete equations. The spec's introduction to the errata and its common transcription of the first row need that distinction (E01, E05).

No algebraic disagreement was found with the spec's master subscale formula, its corrected adjoint convention, its residual-based weighted OSS equations, or any of the eight case equations under the stated assumptions. That conclusion does not extend to the split-OSS claims or the unresolved boundary and solver statements above.
