# EMAC–VMS against the classical convective forms: experiments and findings

This document records every experiment run after Phase 12b to find out whether, and where,
EMAC–VMS (chapter Definition 5.1) is better than the other convective forms under VMS
stabilization. It also covers the full-stress variant of chapter §5.8.5. `README.md` in this
folder explains how to rerun everything; the run data are in `data/`.

## 1. Common setup

Unless a section says otherwise:

| Setting | Value |
| --- | --- |
| Elements | P1/P1 (equal order) |
| Stabilization | VMS, dynamic orthogonal subscales (OSS), nonlinear scale splitting, pressure subscale on |
| Constants | c1 = 4, c2 = 2, so τ_m = (c1 ν/h² + c2 \|a\|/h)⁻¹ and τ_c = h²/(c1 τ_m) |
| Viscous form | stress |
| Time scheme | Crank–Nicolson with one BDF1 startup step |
| Linearization | Newton, frozen τ, no Picard warm-up |
| Tolerances | inner 1e-10, outer 1e-9, subscale 1e-11; up to 200 inner, 100 outer and 500 local subscale iterations |
| Dimension | 2D |
| Software | `main` at c8cf378 plus `emac_vms_full`; Firedrake 2026.10.0, PETSc 3.26, Python 3.12 |

The forms compared, with the short names used below:

| Short name | `convection.form` | Note |
| --- | --- | --- |
| CONV | `nonconservative` | |
| CONS | `conservative` | |
| SKEW | `skew_symmetric_1` | `skew_symmetric_2` gave identical results wherever both were run |
| EMAC | `emac` | resolved EMAC with the standard (fixed) subscale coupling |
| EMAC–VMS | `emac_vms` | chapter Definition 5.1 |
| EMAC–VMS full | `emac_vms_full` | chapter §5.8.5, on branch `experiment/emac-vms-advantage` only |

Quantities. E is the resolved kinetic energy ½‖u_h‖². M and A are the resolved linear and
angular momentum, A about the origin. 𝒩, D_k and D_φ are the measured convective defects of
energy, linear momentum and angular momentum (spec: Conservation checks). The error is the L²
velocity error against the stated reference field. Z is the resolved enstrophy ½‖ω_h‖², with
ω_h the curl of u_h only, because subscales are stored at quadrature points.

What theory predicts on equal-order P1/P1 (chapter Tables 5.2–5.3 and Lemma 5.2). D_k contains
∇·u_h only tested against u_h,k, which lies in the pressure space, where the stabilized
continuity equation makes it vanish. So CONV, CONS and SKEW all conserve linear momentum;
standard EMAC does not, because its D_k contains (ũ, ∂_k u_h). The energy and angular-momentum
arguments (u_h·u_h and u_h·φ) have higher degree than the pressure space, so those defects
remain. The expected pattern is:

| Form | Energy defect 𝒩 | Linear momentum D_k | Angular momentum D_φ |
| --- | --- | --- | --- |
| CONV | nonzero | 0 | nonzero |
| CONS | nonzero | 0 | 0 (nonlinear splitting) |
| SKEW | 0 | 0 | nonzero |
| EMAC | nonzero | nonzero | nonzero |
| EMAC–VMS | 0 | 0 | 0 |

## 2. Experiments

### 2.1 Phase 12b conservation matrix (owner's runs, for context)

The 55 configurations of plan Phase 12b, run by the owner and analysed with
`scripts/analyse_conservation.py`. All 55 completed.

- EMAC–VMS: worst energy closure 1.28e-10 of P_c (limit 1e-9), periodic momentum closure at most
  7.1e-13, all three convective defects below 8e-15.
- Long Gresho, T = 10 (32², δt = 0.01, OSS, ν = 1e-6). Resolved energy lost: CONS 16.8 %,
  SKEW 19.0 %, CONV 20.2 %, EMAC 21.9 %, EMAC–VMS 25.9 %. For EMAC–VMS, 21.8 % is
  velocity-subscale dissipation, 3.9 % pressure-subscale and 0.15 % viscous. Standard EMAC's
  convective defect adds back 5.9 % of E(0), so its gross dissipation is about 27.7 %.
- Long Gresho, L² error at T = 10, read from the plots: EMAC–VMS ≈ EMAC ≈ 0.125, CONS ≈ 0.16,
  SKEW ≈ 0.17, CONV ≈ 0.23. The errors of EMAC and EMAC–VMS grow more slowly after t ≈ 5.
- Refinement (64², δt/2) cuts the short-run EMAC–VMS energy loss from 2.26 % to 0.40 %.
- With ASGS (viscous Gresho, ν = 1e-3) the resolved angular momentum jumps by 1.0–1.2 % in the
  first step for every form: the subscale takes up part of the angular momentum. Resolved-only
  angular momentum is therefore not meaningful for ASGS (Principe et al. 2010 explain why only
  OSS separates the scales cleanly in the energy budget).

### 2.2 Blow-up screening

Question: is there a setting in which the other forms blow up and EMAC–VMS does not, as reported
for Galerkin methods by Charnyi et al. (2017) and Li & Rui (2023)?

Setup: Gresho vortex, ν = 1e-6, δt = 0.01, T = 3. A Galerkin Taylor–Hood P2/P1 anchor on 48²,
the Li & Rui §4.2 resolution (but with no-slip instead of no-penetration walls). Then 72 VMS runs:
six forms × meshes 16², 32² × τ_c on or off × c2 ∈ {2, 32, 512}.

| Galerkin Taylor–Hood 48² | Outcome | Error at the end | ∫𝒩 dt / E(0) |
| --- | --- | ---: | ---: |
| CONV | blows up at t = 2.54 (Li & Rui: ≈ 2.3) | 0.646 | −2.6 |
| CONS | blows up at t = 0.28 | 0.622 | −2.3 |
| SKEW | bounded but fills with grid-scale noise; ‖∇·u_h‖ grows from 0.2 to about 90, and ν = 1e-6 dissipates the noise | 0.496 | 0 |
| EMAC | stable, clean | 0.063 | 0 |

VMS results:

- No energy blow-up in any of the 72 runs: the largest E/E(0) was 1.0004.
- τ_c on with c2 = 512: every form fails at t ≈ 0.3. This is an artefact of the sweep:
  τ_c = ν + (c2/c1) h|a|, so raising c2 strengthens the grad-div term (58 % energy loss already at
  c2 = 32). Those rows are not a "weaker stabilization" test.
- τ_c off with c2 = 512 on 32², the only genuinely weak setting: CONS fails at t = 0.47 (outer
  iteration diverges even with 500 subscale sweeps, energy slightly above E(0)). SKEW1 and SKEW2
  stopped only at the 50-sweep subscale cap and complete with 500 (error 0.031 against 0.029 for
  EMAC–VMS). CONV, EMAC and EMAC–VMS complete.
- Accuracy follows the stabilization strength: weakly stabilized (16², τ_c off, c2 = 512), EMAC and
  EMAC–VMS reach 0.036 against 0.072–0.086; at the default (32², c2 = 2) EMAC–VMS has 0.057
  against 0.036 for CONV and SKEW.

Finding: the Galerkin blow-up of the literature is reproduced, but VMS dissipation removes the
mechanism (energy injection by 𝒩). No survival advantage of EMAC–VMS under VMS.

### 2.3 E1: translating Gresho vortex (linear momentum and energy)

Problem: periodic unit box (−½, ½)², ν = 1e-6, f = 0, initial velocity
u₀ = U₀ + u_G(x) with U₀ = (1, ½) and u_G the Gresho vortex (radius 0.4). By Galilean invariance
the exact Euler solution is u(x, t) = U₀ + u_G(x − U₀t), evaluated with the periodic images;
linear momentum M = U₀|Ω| and the energy are physical invariants (no walls, no forcing). The
vortex centre is tracked as the periodic centroid of ω_h². 32², δt = 0.01, T = 4. Variants:
default OSS, weak OSS (τ_c off, c2 = 32), ASGS.

| Default OSS | E(T)/E(0) | Error | Centre error | ∫𝒩 dt / E(0) | max \|M − M(0)\|/M(0) |
| --- | ---: | ---: | ---: | ---: | ---: |
| CONV | 0.9838 | 0.167 | 0.015 | +2.0e-4 | 4e-13 |
| CONS | 0.9845 | 0.176 | 0.010 | −2.5e-4 | 2e-15 |
| SKEW | 0.9841 | 0.172 | 0.012 | 0 | 4e-15 |
| EMAC | 0.9920 | 0.178 | 0.015 | −9.0e-3 | **4.3e-3** |
| EMAC–VMS | 0.9832 | 0.178 | 0.011 | 4e-16 | 2e-15 |

- Standard EMAC drifts linear momentum by 0.43 %, growing linearly, and injects 0.9 % of E(0);
  it keeps the most energy only because of that injection. EMAC–VMS removes both.
- CONV, CONS and SKEW conserve linear momentum to round-off, as Lemma 5.2 predicts for equal
  order.
- Accuracy is about the same for all forms; EMAC–VMS has the second-best vortex trajectory over
  most of the run.
- Weak OSS: errors 0.250–0.260, EMAC and EMAC–VMS marginally best. ASGS: the same pattern; the
  5e-5 resolved momentum drift of every form is the subscale share, and EMAC drifts 4.1e-3.

### 2.4 E3: off-centre Gresho vortex in a large no-slip box (angular momentum)

Problem: box (−1, 1)² with no-slip walls, ν = 1e-6, the Gresho vortex centred at (0.1, 0.05), so
the gap to the walls is at least 0.5. The exact Euler field is steady, and angular momentum is a
physical invariant as long as the flow stays away from the walls. Monitor region: the wall frame
max(|x|, |y|) ≥ 0.9.

48², T = 10:

| Form | ∫D_φ dt / A(0) | A(10)/A(0) − 1 | E(10)/E(0) | Error | Frame ‖u_h‖ |
| --- | ---: | ---: | ---: | ---: | ---: |
| CONV | +2.9e-2 | **+1.04** | 0.725 | 0.521 | 0.125 |
| CONS | 0 | −1.6e-2 | 0.759 | **0.293** | 0.030 |
| SKEW | +1.6e-2 | −9.5e-2 | 0.745 | 0.343 | 0.037 |
| EMAC | −3.7e-2 | +7.0e-2 | 0.664 | 0.396 | 0.040 |
| EMAC–VMS | 0 | −6.1e-2 | 0.601 | 0.391 | 0.036 |

96², T = 3:

| Form | ∫D_φ dt / A(0) | A(3)/A(0) − 1 | E(3)/E(0) | Error | Frame ‖u_h‖ |
| --- | ---: | ---: | ---: | ---: | ---: |
| CONV | +4.8e-3 | −3.4e-3 | 0.974 | 0.034 | 1.5e-3 |
| CONS | 0 | +1.5e-3 | 0.975 | **0.032** | 1.2e-3 |
| SKEW | +2.4e-3 | −0.9e-3 | 0.975 | 0.033 | 1.3e-3 |
| EMAC | −0.5e-3 | +3.1e-3 | 0.967 | 0.059 | 3.5e-3 |
| EMAC–VMS | 0 | +2.4e-3 | 0.959 | 0.054 | 3.2e-3 |

- The measured spurious torque follows the theory exactly: zero for CONS and EMAC–VMS,
  nonzero for CONV, SKEW and EMAC.
- The angular-momentum drift itself does not follow D_φ: it is dominated by wall torque once
  velocity reaches the walls. On 48² the vortex is under-resolved (h = 1/24) and smears into the
  walls; CONV degrades badly after t ≈ 6. On 96² every form keeps A within 0.4 %.
- On 48² at t = 3 the errors are about even (EMAC 0.136, EMAC–VMS 0.147, CONV 0.151, SKEW 0.163,
  CONS 0.170). On 96², and at T = 10, EMAC and EMAC–VMS are about 1.7× worse than the classical
  forms and leak more velocity towards the walls.
- The weak-OSS variant (48², T = 10) destroys the vortex for every form (error ≈ 0.56); it tells
  nothing.

### 2.5 E2: thin double shear layer (under-resolved robustness)

Problem: doubly periodic unit square, ν = 1e-4, f = 0,
u = tanh(ρ(y − ¼)) for y ≤ ½ and tanh(ρ(¾ − y)) for y > ½, v = δ sin(2πx), with ρ = 80 and
δ = 0.05. This follows the Bell–Colella–Glaz / Minion–Brown thin shear layer, with constants
taken from memory and not checked against those papers; our own 128² runs (SKEW and EMAC–VMS, on
4 MPI ranks) are the reference. δt = 0.005, T = 1.5 on 32² and 64², T = 1.2 on 128². Vorticity
snapshots at t = 0.8, 1.0, 1.2, 1.5.

| Form | E(1.2)/E(0), 32² | 64² | 128² | Z(1.2)/Z(0), 64² | Braid vortices at 64² |
| --- | ---: | ---: | ---: | ---: | --- |
| CONV | 0.930 | 0.949 | | 0.395 | yes |
| CONS | 0.936 | 0.952 | | 0.396 | yes |
| SKEW | 0.932 | 0.950 | 0.961 | 0.393 | yes |
| EMAC | 0.833 | 0.908 | | 0.315 | no |
| EMAC–VMS | 0.828 | 0.905 | 0.946 | 0.313 | no |

- At 128² SKEW and EMAC–VMS agree on the structure: two rolled-up cores and clean, thin braids.
  At 64² CONV, CONS and SKEW grow secondary vortices on the braids, the known under-resolution
  artefact; EMAC and EMAC–VMS keep the correct topology (`presentation/7_shear_layer_resolution.png`).
- EMAC and EMAC–VMS are more dissipative at every resolution (energy gap to SKEW 10 %, 4.5 %,
  1.5 % on 32², 64², 128²) and their cores are weaker. This flow is smooth, so the dissipation
  penalty is not only a Gresho effect (section 3).
- Convective energy injection by T = 1.5: CONS −1.3e-2 and −4.6e-3 of E(0) on 32² and 64², EMAC
  −6.9e-3 and −4.4e-3, CONV +5.6e-3 and +2.6e-3; SKEW and EMAC–VMS at round-off.

### 2.6 E4 and E5: smooth vortex

The "smooth vortex" is not a published benchmark. It was constructed for this study to test
whether the EMAC penalty on Gresho comes from Gresho's non-smooth profile. It is a compactly
supported swirl with the same support radius as Gresho:

```
u(x) = 12 (1 − |x − x_c|²/0.16)₊⁴ (−(y − y_c), x − x_c),   that is
u_θ(r) = 12 r (1 − r²/0.16)⁴ for r < 0.4,   u_θ = 0 for r ≥ 0.4.
```

- Any purely radial swirl u_θ(r) is divergence-free and is a steady solution of the Euler
  equations with the pressure p′(r) = u_θ²/r, so the velocity above is an exact steady Euler field
  (and, as for Gresho, a reference field at ν = 1e-6). The pressure was not needed, since only
  velocity errors are measured.
- The core is a rigid rotation (u_θ ≈ 12r), the peak speed is ≈ 1 at r = 0.4/3, and the profile
  is C³ at r = 0.4. Unlike Gresho, whose u_θ has kinks at r = 0.2 and 0.4, |u|² and hence the
  EMAC pressure P = p − ½|u|² are smooth.
- Being compact, it has zero net circulation; E(0) ≈ 0.045, against 0.084 for Gresho.
- It is steeper than Gresho near its peak, so 48² in the (−1, 1)² box is quite under-resolved
  (energy losses of 18–62 % by t = 3).

E4 is E1 with this vortex (periodic box, U₀ = (1, ½), 32², T = 4); E5 is E3 with it (no-slip box,
centre (0.1, 0.05), 48², T = 3).

| Form | E4 error | E4 centre error | E4 E(4)/E(0) | E5 error | E5 A(3)/A(0) − 1 | E5 E(3)/E(0) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CONV | 0.177 | 0.022 | 0.972 | 0.242 | −6.0e-2 | 0.688 |
| CONS | 0.169 | 0.018 | 0.975 | 0.191 | −3.0e-3 | 0.738 |
| SKEW | 0.174 | 0.020 | 0.974 | 0.219 | −3.3e-2 | 0.707 |
| EMAC | 0.172 | 0.026 | 0.986 | **0.126** | +2.4e-2 | 0.704 |
| EMAC–VMS | **0.165** | **0.014** | 0.973 | 0.136 | **−3.3e-3** | 0.643 |

- E4: EMAC–VMS has the lowest error and the best trajectory at the same energy loss as the others.
- E5: EMAC and EMAC–VMS have 30–45 % lower error than the classical forms although EMAC–VMS loses
  the most energy; EMAC–VMS has the smallest angular-momentum drift, tied with CONS.
- Compared with Gresho on the same 48² box, where the EMAC forms were even at t = 3 and worse
  later and on 96², the ranking flips (`presentation/4_smooth_vs_gresho_accuracy.png`).

### 2.7 Time-step sweep (is the time integrator the source of the extra dissipation?)

Gresho with the long-run settings (32², OSS, ν = 1e-6), T = 2, δt ∈ {0.02, 0.01, 0.005, 0.0025},
EMAC–VMS and SKEW, run with `scripts/run_simulation.py`.

| Form | E(2)/E(0) over the four δt | ∫ velocity-subscale dissipation | ∫ pressure subscale | ∫ viscous | Error |
| --- | --- | ---: | ---: | ---: | ---: |
| SKEW | 0.9502–0.9507 | 3.75 % | 1.15 % | 0.04 % | 0.0345 |
| EMAC–VMS | 0.9203–0.9208 | 6.67 % | 1.17 % | 0.04 % | 0.0411 |

- The energy loss changes by less than 0.05 % over an 8× range of δt; with Crank–Nicolson the
  discrete energy identity is exact (chapter Remark 5.12) and the budget closes to 1e-10. The time
  integrator is ruled out.
- Pressure-subscale and viscous dissipation are identical for the two forms. All of the extra loss
  of EMAC–VMS is velocity-subscale dissipation τ_m⁻¹‖ũ‖², 1.8× that of SKEW
  (`presentation/5_dissipation_budget_dt_sweep.png`).

### 2.8 Full-stress variant (chapter §5.8.5)

`emac_vms_full` uses T̃_full = T̃ + (u_h·ũ + ½|ũ|²) I (Eq. 5.37), the residual
N_E = 2(∇ˢu_h)u_h + (∇·u_h)u_h + M_h ũ with π_h = P_h (Eq. 5.38), and the local-update matrix
M_h = ∇ˢu_h + ½(∇·u_h) I. Verification: 𝒩, D_k and D_φ at 1e-17–1e-16, energy budget closed to
4e-14, and a unit test of the pointwise identity T̃_full : ∇ˢu_h = ũ·N_E for both splittings.

| Case | EMAC–VMS: energy lost / error | Full variant: energy lost / error |
| --- | --- | --- |
| E1 translating Gresho (32², T = 4) | 0.017 / 0.178 (centre 0.011) | 0.025 / 0.194 (centre 0.022) |
| E3 box Gresho (96², T = 3) | 0.041 / 0.054 | 0.049 / 0.055 |
| E3 box Gresho (48², T = 10) | 0.399 / 0.391 | 0.479 / 0.398 |
| E2 shear layer (64², t = 1.2) | 0.095, Z = 0.313 | 0.098, Z = 0.302 |

Finding: same exact conservation, more dissipation, no accuracy gain; on E1 the vortex-centre
error doubles. The full variant is not worth it; Definition 5.1, which keeps the standard subscale
model, is the better choice (`presentation/6_full_stress_variant.png`).

### 2.9 Matched-dissipation comparison (EMAC–VMS against SKEW)

Question: are the advantages of EMAC–VMS (smooth-vortex accuracy, shear-layer topology) structural,
or a side effect of its extra damping? Method: multiply τ_m by s through c1/s and c2/s, which
keeps τ_c = h²/(c1 τ_m) unchanged (verified numerically); EMAC–VMS with s ∈ {1, 0.5, 0.25} and SKEW
with s ∈ {0.25, 0.5, 1, 2, 4}. Compare at equal energy loss and at each form's best τ.

| Case | EMAC–VMS (s = 1, 0.5, 0.25): energy lost / error | SKEW (s = 0.25, 0.5, 1, 2, 4): energy lost / error |
| --- | --- | --- |
| E4 translating smooth vortex | 0.027/0.165, 0.021/0.158, 0.017/0.150 | 0.016/0.153, 0.021/0.164, 0.027/0.174, 0.033/0.181, 0.041/0.197 |
| E5 smooth vortex in the box | 0.357/0.136, 0.247/0.137, 0.178/0.150 | 0.155/0.209, 0.208/0.216, 0.293/0.219, 0.422/0.211, 0.625/0.218 |
| E3 Gresho, 96² | 0.041/0.054, 0.025/0.041, 0.017/0.030 | 0.013/0.023, 0.017/0.027, 0.026/0.033, 0.041/0.040, 0.072/0.050 |
| E2 shear layer, 64² (Z(1.2)/Z(0) in place of error; 128²: 0.396 SKEW, 0.376 EMAC–VMS) | 0.095/0.313, 0.074/0.349, 0.061/0.376 | 0.043/0.421, 0.046/0.410, 0.050/0.393, 0.056/0.366, 0.066/0.334 |

- Both forms improve as τ_m decreases: the default constants over-stabilize these nearly inviscid
  flows, EMAC–VMS more so.
- Gresho, 96²: SKEW is better at every matched loss and at its best τ (0.023 against 0.030).
- Translating smooth vortex: a tie at the best τ (0.153 against 0.150).
- Smooth vortex in the box: EMAC–VMS is about 35 % more accurate at every τ.
- Shear layer: SKEW grows spurious braid vortices at every τ, including τ_m × 4, where it
  dissipates more than EMAC–VMS at τ_m × 0.25; EMAC–VMS never does. At small τ the SKEW enstrophy
  overshoots the 128² value, the signature of those spurious vortices
  (`presentation/8_matched_dissipation.png`, `presentation/9_shear_layer_matched.png`).

Finding: the differences are structural, not damping, in both directions.

### 2.10 E6 and E7: rigid-rotation and shear robustness tests

Motivated by the angular-momentum robustness tests of Brubeck, Parker and Zerbinati (FETC 2026
slides, "On the symmetry constraint and angular momentum conservation in mixed stress formulation")
and by pressure robustness (Linke et al. 2017). Two exact steady Navier–Stokes solutions on
(−½, ½)² with their own Dirichlet data, ν = 1e-6, f = 0, started from the exact field, 10 steps of
δt = 0.01, C ∈ {1, 10, 100}, meshes 16² and 32² (`driver_rotation.py`):

- E6, rigid rotation u = C(−y, x), p = C²r²/2. The convective term (u·∇)u = −∇p is a pure
  gradient that the classical forms' P1 pressure must balance, but p is quadratic. For EMAC,
  ∇ˢu = 0 makes 2(∇ˢu)u + (∇·u)u vanish and P = p − ½|u|² = 0.
- E7, parallel shear u = C(y, 0), p constant. Now (u·∇)u = 0, so the classical forms need no
  pressure, while EMAC must balance ∇½|u|² = C²y e_y with P = −½C²y², which is quadratic.

Relative L² velocity error after 10 steps (16² / 32²):

| Form | E6, C = 1 | E6, C = 10 | E6, C = 100 | E7, C = 1 | E7, C = 10 | E7, C = 100 |
| --- | --- | --- | --- | --- | --- | --- |
| CONV | 4.7e-4 / 9.4e-5 | 1.2e-3 / 2.5e-4 | 6.1e-3 / 1.5e-3 | exact | exact | exact |
| CONS | 4.7e-4 / 9.4e-5 | 1.1e-3 / 2.3e-4 | 2.7e-3 / 6.6e-4 | exact | exact | exact |
| SKEW | 4.7e-4 / 9.4e-5 | 1.0e-3 / 2.2e-4 | 2.6e-3 / 6.4e-4 | exact | exact | exact |
| EMAC | exact | exact | exact | 5.1e-4 / 1.1e-4 | 2.8e-3 / 9.7e-4 | no convergence / 9.7e-3 |
| EMAC–VMS | exact | exact | exact | 5.1e-4 / 1.1e-4 | 2.8e-3 / 9.6e-4 | 2.2e-2 / 9.8e-3 |

Figure: `presentation/10_rotation_shear_robustness.png` (also C = 3 and 30; `robustness_figure.py`).
"Exact" means 1e-17 to 1e-16. The nonzero errors converge like h² (P1 interpolation of a quadratic
pressure), and their relative size grows with C (5–13× for the classical forms from C = 1 to 100,
40× for EMAC–VMS), so their absolute size grows like C² or faster.

Finding: each form leaves a different part of the convective term to the P1 pressure. The classical
forms must balance the centripetal gradient, which dominates in vortices; the EMAC forms must balance
the Bernoulli term ∇½|u|², which dominates in shear. EMAC and EMAC–VMS are exact for rigid rotation
at any speed, the classical forms for parallel shear. This one mechanism explains the earlier
results: EMAC–VMS wins on vortex cores and the smooth vortex (E4, E5, long Gresho), and loses on
Gresho's strained ring and in shear layers through extra dissipation (E2, E3, section 3). The
robustness is a property of the EMAC structure, shared by standard EMAC; EMAC–VMS adds exact
conservation. E6 and E7 together are also a sharp design target: a scheme that is better than both
must be exact in both tests.

## 3. Why EMAC–VMS dissipates more

Principe, Codina and Henke (2010) give the numerical dissipation of OSS as
τ_m‖P⊥(N + ∇π_h)‖² + τ_c‖P⊥(∇·u_h)‖² per element (for linear elements). Only the convective operator
N inside the orthogonal projection depends on the form. For the classical forms it is (a·∇)u_h + ∇p_h;
for EMAC–VMS it is (a·∇)u_h + ∇½|u_h|² + ∇P_h, with P_h approximating p − ½|u|². Continuously
the two residuals are identical, but the P1 pressure must now also represent ½|u_h|². In a shear layer
(u·∇)u = 0 and p is nearly flat, while ½|u|² varies strongly across the layer, so the EMAC residual
leaves a large orthogonal part and hence large τ_m dissipation; in rigid rotation (a vortex core)
the EMAC residual vanishes. Gresho's kinked |u|² makes P harder still to represent, which is
consistent with the Gresho penalty and the smooth-vortex reversal. The δt sweep confirms that all
of the extra loss is this velocity-subscale term. Under refinement τ_m⁻¹ → c1ν/h² and the dissipation
becomes viscous, which is why the gap shrinks with h. Changing N to avoid this would break the exact
energy identity T̃ : ∇ˢu_h = ũ·N, so the remedies are resolution, the polynomial degree or recalibrated
constants (Guasch and Codina 2013 show that the OSS dissipation constant depends on the
stabilization parameters).

## 4. Conclusions

1. **Conservation.** EMAC–VMS is the only form with energy, linear-momentum and angular-momentum
   convective defects all at round-off (`presentation/1_defect_matrix.png`), confirmed in flows where
   each invariant is nonzero.
2. **EMAC–VMS against standard EMAC** is the clearest result: standard EMAC under VMS drifts momentum,
   injects energy and adds spurious torque; EMAC–VMS removes all three at the same accuracy
   (`presentation/2_emac_vs_emac_vms.png`).
3. **No survival advantage under VMS.** The Galerkin blow-ups are reproduced, but stabilization removes
   them for every form (`presentation/3_galerkin_blowup_anchor.png`).
4. **Against SKEW, the result depends on the flow.** SKEW is less dissipative and more accurate on the
   Gresho vortex; EMAC–VMS is more accurate on a smooth vortex and keeps the correct topology of an
   under-resolved shear layer. These differences survive matched dissipation, so they are structural.
   SKEW is the less dissipative, generally safe choice; EMAC–VMS adds exact angular momentum, robust
   shear-layer topology and higher accuracy on smooth vortices.
5. **The extra dissipation of EMAC–VMS** is velocity-subscale dissipation caused by representing ½|u|²
   in the P1 pressure, not the time integrator.
6. **The default stabilization constants are too large** for these nearly inviscid flows, for both forms
   and more so for EMAC–VMS.
7. **The full-stress variant is not worth it.**
8. **Robustness is flow-type dependent and explains everything above** (section 2.10): EMAC and
   EMAC–VMS reproduce rigid rotation exactly at any speed, where the classical forms' error grows
   with the speed; for parallel shear it is the reverse. In a pressure-robustness sense, EMAC is
   "rotation-robust" and the classical forms are "shear-robust".

## 5. Limitations

- 2D only, P1/P1 only, mostly dynamic OSS with nonlinear splitting; one mesh per case in most comparisons.
- Short and moderate times (T ≤ 10); no turbulent flow.
- The smooth vortex is a constructed test, run on 48² and 32² only; whether its advantage persists on
  finer meshes was not tested.
- The E2 constants are not checked against the original papers; the reference is our own 128² runs,
  which differ slightly between SKEW and EMAC–VMS.
- Braid topology was judged visually.
- In this environment the Firedrake unit suite has 115 failures that also occur on `main`: the
  all-Dirichlet flux check compares round-off with round-off for `manufactured_spatial`, and
  `_parse_value` handles `X | None` unions only from Python 3.14. They do not affect these runs.

## 6. Open questions

- A scheme exact in both E6 and E7, for example EMAC–VMS with a projection space enriched by the
  gradients of quadratic edge bubbles (so that the P1-unrepresentable part of ∇½|u_h|² is not
  dissipated), or a pressure-robust (divergence-free reconstruction) variant in the sense of Linke.
- Recalibrate c1 and c2 per form (for example by matching dissipation or by the Guasch–Codina
  argument) and repeat the comparison at calibrated constants.
- The smooth vortex on finer meshes and over longer times.
- A higher-order pair (P2/P2), to test whether the EMAC pressure penalty shrinks.
- The oblique dipole–wall collision (Clercx and Bruneau 2006) with an angular-momentum budget based on
  the consistent wall reaction, and the Kelvin–Helmholtz benchmark of Schroeder et al. (2019) with
  free-slip walls; their exact definitions could not be retrieved in this environment.
- A turbulent case (3D Taylor–Green, Re = 1600), which needs an iterative solver and a cluster.

## 7. References consulted

- Charnyi, Heister, Olshanskii, Rebholz (2017): EMAC and its Galerkin conservation properties
  (cited through the chapter and secondary sources).
- Li & Rui (2023), J. Sci. Comput. 96:69: modified convective formulation; Gresho and lattice-vortex
  setups used for the Galerkin anchor.
- Principe, Codina & Henke (2010), CMAME 199:791–801: the dissipative structure of OSS.
- Guasch & Codina (2013), CMAME 261–262:154–166: statistical behaviour of OSS stabilization terms.
- Brubeck, Parker & Zerbinati (2026), slides, Workshop on Finite Element Tensor Calculus, Sanya:
  symmetry constraint and angular-momentum robustness in mixed stress formulations (motivated E6).
- Linke et al. (2017), SIAM Review: pressure robustness (cited in those slides; not read here).

## 8. Files

| Path | Content |
| --- | --- |
| `presentation/1_…9_*.png` | presentation figures (sections 2.2–2.9) |
| `figures/` | diagnostic panels for E1, E3 and E2 |
| `screen/screening.md`, `matched.md`, `summary.md` | generated tables |
| `driver.py`, `driver_e2.py`, `screen/generate.py` | run drivers |
| `report.py`, `presentation.py`, `matched.py`, `screen/analyse.py` | analysis scripts |
| `data/runs.tar.xz`, `data/screening.tar.xz` | per-run diagnostics, step logs, configs, E2 snapshots |
