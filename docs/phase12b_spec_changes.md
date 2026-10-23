# Phase 12b — conservation runs: spec and plan additions

Chapter reference: the revised `chapter5_preview.pdf`, 23 pages, SHA-256
`c4995958c5e410591261f676aae661f5854def033feb534dc930e40668491784`. It adds §5.7 "Consistent
splitting"; no implemented equation changed. `chapter5_renumbering.md` maps old to new numbers.

Phase 12b shows numerically what Chapter 5 proves. For each of the six convective forms:
- which convective defects appear;
- whether the energy balance closes;
- how momentum, angular momentum and energy evolve in a flow that should conserve them.

It runs in two parts:
- **12b-1** (agent): the Gresho case, two diagnostics, the configs, the queue and the analysis
  script, all with unit tests.
- **12b-2** (owner): launch the queue. A later analysis phase reads the results.

## A. Spec edits (`formulation_spec.md`)

1. **Chapter references.** Point the Chapter 5 sections to the new hash and update every
   chapter number using `chapter5_renumbering.md`:
   - the fixed-coupling defect tables are now Tables 5.2 and 5.3;
   - EMAC–VMS is §5.8;
   - its propositions moved by two.

   Add one sentence to the "Convective forms" section: the master formulation is the *fixed
   coupling* of chapter §5.7, in which the common subscale coupling equals the conservative
   form's fine-scale part (Lemma 5.3).
2. **Gresho case.** Add a new section, "Gresho vortex (Phase 12b)".
   - **Domain.** (−0.5, 0.5)², built-in mesh, no-slip on all four walls, f = 0. The pressure
     null space is active.
   - **Exact solution.** It is a steady Euler solution, the same for all t; r = |x|.
     - **Swirl velocity u_θ:** 5r for r < 0.2; 2 − 5r for 0.2 ≤ r < 0.4; 0 for r ≥ 0.4.
       The velocity is u = u_θ (−y/r, x/r), and u = (−5y, 5x) for r < 0.2.
     - **Pressure p:**
       - 12.5r² + C₁ for r < 0.2;
       - 12.5r² − 20r + 4 ln r + C₂ for 0.2 ≤ r < 0.4;
       - 0 for r ≥ 0.4;

       with C₂ = −12.5(0.4)² + 20(0.4) − 4 ln 0.4 and C₁ = C₂ − 20(0.2) + 4 ln 0.2.
     - **Vorticity:** 10, 2/r − 10, and 0 in the same three regions.
     - The initial velocity is the exact one, interpolated.
   - **Viscosity (owner Q1).** A Gresho case parameter, default ν = 1e-6, the authorized
     fallback for the singular frozen-τ OSS system at ν = 0. The inviscid families use OSS
     only. The viscous ASGS families use ν = 1e-3; at positive viscosity the Euler fields
     are reference fields, not exact viscous solutions (also stated in the docstring and README).
     Where a = 0 and ν ≈ 0, τ_m⁻¹ → 0 (zero in the inviscid limit), the dynamic ASGS subscale
     absorbs ∂ₜu_h, and the resolved mass cancels. This is an accepted method limitation,
     not a solver problem. Keep the strict nonlinear ASGS xfail at ν ≤ 1e-4. Finite frozen-τ
     steps and diagnostics are required for both OSS splittings at ν = 1e-6 and both ASGS
     splittings at ν = 1e-3; the exact-τ failure at rest remains documented under K-D1.
   - **Diagnostic data.**
     - Force boundaries: all four walls.
     - Monitor region: the near-wall frame max(|x|, |y|) ≥ 0.45.
3. **Momentum budget diagnostic** `momentum_budget`. Same restrictions as `energy_budget`:
   Crank–Nicolson, dynamic VMS subscales, blank at the startup step. At step n ≥ 2:
   - M_k = ∫(u_h + χ_A ũ)_k at the endpoints. These are the quantities the discrete
     equation contains; with OSS, ũ does not enter.
   - rate_k = (M_kⁿ⁺¹ − M_kⁿ)/δt.
   - D_k is the assembled convective term tested with e_k at the solve level; this is the
     12a defect.
   - F_k = (f^{n+½}, e_k).
   - r_M,k = rate_k + D_k − F_k.

   The viscous, pressure, pressure-subscale and viscous-subscale terms vanish for e_k. So:
   - **without boundaries** (periodic), r_M,k = 0 to solver tolerance;
   - **with Dirichlet boundaries**, r_M,k equals the consistent boundary reaction (Hughes &
     Wells). It is reported, never asserted.
4. **Predicted momentum defects.** If 12a provides only the measured D_k and D_φ, add the
   predicted values from chapter Tables 5.2 and 5.3, assembled with the ∇·a pairing. That
   allows measured and predicted to be compared like 𝒩. EMAC–VMS predicts 0.
5. **Monitor-region diagnostics** `region_velocity_l2` and `region_subscale_l2`:
   - the L² norms of u_h and of ũ over the case's monitor region;
   - they require the case to declare the region, via an optional `Problem` field;
   - they show how far Charnyi's "zero strip" assumption holds. They are reported only.

## B. Run matrix and settings (owner decisions)

Common settings for every 12b run:

| Setting | Value |
| --- | --- |
| Elements | P1/P1 |
| Viscous form | `stress` (required for `emac_vms`; also used for the other forms, for comparability) |
| Pressure subscale | on |
| Subscales | dynamic only |
| Time scheme | Crank–Nicolson |
| Linearization | Newton for all forms, no Picard warm-up |
| `newton_tau` | `frozen` (default) |
| `outer_rtol` | 1e-11 |
| `outer_max_iterations` | 100 |
| `inner_rtol` | 1e-12 |
| `inner_max_iterations` | 200 for every Phase 12b configuration (owner Q2) |
| `subscale_rtol` | 1e-13 |

The conservation identities hold at convergence, so the tolerances bound the measurable
residuals. τ_c-off runs are excluded.

The six forms, in this order: nonconservative, conservative, skew_symmetric_1,
skew_symmetric_2, emac, emac_vms. Each "×" below means every combination.

| Family | Matrix | Mesh | Time | Diagnostics |
| --- | --- | --- | --- | --- |
| `gresho_{form}_oss_{split}` | 6 forms × {linear, nonlinear} = 12; OSS | 32×32, ν = 1e-6 | dt 0.005, T 0.5 | kinetic and subscale energy, linear and angular momentum, energy and momentum budgets, measured and predicted defects, force (walls), region monitors, velocity and pressure error, divergence |
| `gresho_visc_{form}_asgs_{split}` | 6 forms × {linear, nonlinear} = 12; ASGS | 32×32, ν = 1e-3 | dt 0.005, T 0.5 | same Gresho diagnostics |
| `gresho_long_{form}` | nonconservative, conservative, skew_symmetric_1, emac, emac_vms; OSS, nonlinear | 32×32, ν = 1e-6 | dt 0.01, T 10 | as above (for comparison with Charnyi et al. 2017, Fig. 2) |
| `gresho_refined_emac_vms_oss` | emac_vms, OSS, nonlinear / 1 | 64×64, ν = 1e-6 | dt 0.0025, T 0.5 | as above |
| `gresho_visc_refined_emac_vms_asgs` | emac_vms, ASGS, nonlinear / 1 | 64×64, ν = 1e-3 | dt 0.0025, T 0.5 | as above |
| `taylor_green_{form}_{space}_{split}` | 6 forms × {asgs, oss} × {linear, nonlinear} = 24 | periodic 32×32, ν = 0.01 | dt 0.005, T 0.5 | kinetic and subscale energy, linear momentum, energy and momentum budgets, measured and predicted defects, velocity error |

The revised total is 55 runs and 10,200 steps: 1200 OSS Gresho + 1200 viscous ASGS Gresho +
5000 long OSS + 200 refined OSS + 200 refined viscous ASGS + 2400 Taylor–Green. The viscous
ASGS runs replace the 13 near-inviscid ASGS runs. The authorized ten-step
`gresho_visc_emac_vms_asgs_nonlinear` prefix completed in 29.4 s, with steps 3–10 averaging
0.61 s. The revised estimate is 2.1–2.6 h serial on one core; the progress log records the
breakdown and refined-mesh/startup assumptions.

- **Taylor–Green** uses the existing `taylor_green_2d` case on a periodic mesh. If the case
  cannot run fully periodic, stop and ask.
- **No angular momentum on periodic domains:** e_k × x is not periodic.

Expected outcomes, for the later analysis phase. These are reported, never tuned:
- **Energy closure:** |r| ≤ 1e-9 × the budget scale from step 2, for every converged run.
- **Momentum closure** (Taylor–Green only): |r_M| ≤ 1e-9 × scale.
- **EMAC–VMS:** measured D_k, D_φ and 𝒩 ≤ 1e-12 relative.
- **The other forms:** measured defect = predicted to 1e-9 relative.
- **Non-converging configurations** (e.g. emac with nonlinear ASGS) are recorded with their
  iteration counts, not removed.
- **Gresho drifts** are reported:
  - over time: M, L, E and the velocity error;
  - for context: the wall force and the region monitors.

## C. Plan edits (`implementation_plan.md`)

1. Replace the provisional Phase 12b brief with the matrix and settings of §B.
2. Split it into:
   - **12b-1:** implementation and unit tests; may run unit-scale prefixes only.
   - **12b-2:** owner launch; the agent never launches it.
3. Note that the consistent splitting needs no code. For the conservative form it coincides
   with the fixed coupling; for non-conservative it conserves nothing; for skew it needs ∇·ũ
   in the pointwise residual; for EMAC it is the not-implemented §5.8.5 variant.
