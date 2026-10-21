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
