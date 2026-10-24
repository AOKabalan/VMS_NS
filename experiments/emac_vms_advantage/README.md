# EMAC–VMS advantage experiments (exploratory, after Phase 12b)

Where does EMAC–VMS differ from the other convective forms, and is the full-stress variant of
chapter §5.8.5 worth it? Everything here is P1/P1, Crank–Nicolson, 2D, mostly dynamic OSS with
nonlinear splitting, the default stabilization (c1 = 4, c2 = 2, pressure subscale on) and
ν = 1e-6 unless stated. These are exploratory runs, not part of the Phase 12b matrix.

## Experiments

| Name | Flow | Invariant tested | Reference |
| --- | --- | --- | --- |
| screening | Gresho, T = 3: Galerkin Taylor–Hood 48² (Li & Rui 2023 §4.2 setup) and 72 VMS runs over mesh, τ_c and c2 | survival | none (blow-up) |
| E1 | Gresho vortex translating with U₀ = (1, ½) in a periodic unit box, 32², T = 4 (OSS, weak OSS, ASGS) | linear momentum, energy | exact u = U₀ + u_G(x − U₀t) |
| E2 | thin double shear layer (ρ = 80, δ = 0.05, ν = 1e-4), periodic, 32², 64², 128² | under-resolved robustness | own 128² runs |
| E3 | off-centre Gresho vortex in a no-slip box (−1, 1)², 48² (T = 10) and 96² (T = 3) | angular momentum | exact steady field |
| E4, E5 | as E1 and E3 with a smooth compact swirl u_θ = 12 r (1 − r²/0.16)⁴ | accuracy without kinks in \|u\|² | exact |
| dt sweep | Gresho long-run settings, 32², T = 2, δt = 0.02 … 0.0025 | time-integrator dependence | — |
| full variant | `emac_vms_full` (chapter Eq. 5.37–5.38) on E1, E2 (64²), E3 (48², 96²) | dissipation, accuracy | — |

The E2 constants follow the Bell–Colella–Glaz / Minion–Brown thin layer from memory and were not
checked against those papers; the 128² runs are the judge.

## Findings

1. **Conservation.** EMAC–VMS is the only form with the energy, linear-momentum and
   angular-momentum convective defects all at round-off (about 1e-16), in flows where each
   invariant is nonzero (`presentation/1_defect_matrix.png`). On equal-order P1/P1 the
   nonconservative, conservative and skew-symmetric forms also conserve linear momentum
   (Proposition 5.1); nonconservative and skew add spurious torque, nonconservative and
   conservative spurious energy.
2. **Against standard EMAC** (fixed coupling), the clearest result: it drifts momentum by 0.4 %,
   injects about 1 % of E(0) on E1 and about 10 % on E3, and adds 3.7 % spurious torque. EMAC–VMS
   removes all three at the same accuracy (`presentation/2_emac_vs_emac_vms.png`).
3. **No blow-up advantage under VMS.** Galerkin Taylor–Hood reproduces the literature
   (nonconservative blows up at t ≈ 2.5, conservative at t ≈ 0.3, skew fills with grid-scale
   noise, EMAC stays clean; `presentation/3_galerkin_blowup_anchor.png`). With VMS no form blows
   up; only with stabilization weakened almost to nothing does the conservative form fail.
4. **Accuracy depends on the flow** (`presentation/4_smooth_vs_gresho_accuracy.png`). On the
   smooth vortex EMAC and EMAC–VMS have 30–45 % lower error than the classical forms (E5), and on
   its translating version EMAC–VMS has the lowest error and best trajectory at equal energy loss
   (E4). On the Gresho vortex they are even on 48² at t = 3 but about 1.7× worse on 96²
   and at T = 10. On the under-resolved shear layer (64²) the classical forms grow spurious braid
   vortices that the 128² runs do not have; EMAC and EMAC–VMS keep the correct topology
   (`presentation/7_shear_layer_resolution.png`), at the cost of weaker cores.
5. **Cost: more dissipation.** EMAC–VMS loses more energy than skew everywhere (shear layer: 10 %,
   4.5 %, 1.5 % more on 32², 64², 128²). The loss is independent of δt over an 8× range and is
   entirely velocity-subscale dissipation, 1.8× that of skew; pressure-subscale and viscous
   dissipation are identical (`presentation/5_dissipation_budget_dt_sweep.png`). With Principe,
   Codina and Henke (2010), the dissipation is τ_m‖P⊥(N + ∇π)‖²: the EMAC pressure must also carry
   ½|u_h|², which P1 resolves poorly wherever there is shear, and Gresho's kinked |u|² makes it
   worse.
6. **Full-stress variant: not worth it.** Its conservation is exact as claimed, and the pointwise
   identity T̃_full : ∇ˢu_h = ũ·N_E holds (unit test), but it dissipates more than Definition 5.1
   in every case and is not more accurate; on E1 its vortex-centre error is twice as large
   (`presentation/6_full_stress_variant.png`). Definition 5.1, which keeps the standard subscale
   model, is the better choice.

7. **Matched dissipation** (`matched.py`, `presentation/8_matched_dissipation.png`,
   `presentation/9_shear_layer_matched.png`, `matched.md`). τ_m is multiplied by s through c1/s
   and c2/s, which keeps τ_c = h²/(c1 τ_m): EMAC–VMS with s = 1, 0.5, 0.25 and skew with
   s = 0.25 … 4, so each form can be compared at equal energy loss and at its own best τ.
   - Both forms improve as τ_m decreases: the default constants over-stabilize these nearly
     inviscid flows, EMAC–VMS more so.
   - Gresho, 96²: skew is better at every matched loss and at its best τ (0.023 against 0.030).
   - Translating smooth vortex: a tie at the best τ (0.153 against 0.150).
   - Smooth vortex in the box: EMAC–VMS is about 35 % more accurate at every τ
     (0.136–0.150 against 0.209–0.219).
   - Shear layer, 64²: skew grows spurious braid vortices at every τ, even when it dissipates
     more than EMAC–VMS; EMAC–VMS never does. At small τ skew's enstrophy overshoots the 128²
     value, the signature of those vortices.

   So the shear-layer topology and the smooth-vortex accuracy are structural, not damping, and
   the Gresho penalty is structural too. Against skew-symmetric 1 the comparison is
   flow-dependent: skew is the less dissipative, generally safe choice, and EMAC–VMS adds exact
   angular momentum, robustness of shear-layer topology and higher accuracy on smooth vortices.

Open: the smooth vortex on a finer mesh; stabilization constants recalibrated per form; a turbulent
case (3D Taylor–Green, Re = 1600), for which Guasch and Codina (2013) give the implicit-LES
framework.

## Code

`emac_vms_full` is a convective form on this branch only: `convection.form = "emac_vms_full"`,
with the same requirements as `emac_vms` (VMS, stress viscosity, Newton). It changes the subgrid
stress to T̃ + (u_h·ũ + χ_N ½|ũ|²) I, the residual to N_E of Eq. (5.38) with π_h = P_h, and the
implicit local update matrix to M_h = ∇ˢu_h + ½(∇·u_h) I. The reported pressure still adds only
½|u_h|², so pressure outputs omit u_h·ũ + ½|ũ|² for this form.

## Reproducing

Run from this folder inside the Firedrake environment, with `OMP_NUM_THREADS=1`:

```bash
# Unpack the recorded per-run data (diagnostics, step logs, resolved configs, E2 snapshots):
tar -xJf data/runs.tar.xz && tar -xJf data/screening.tar.xz
python presentation.py           # presentation/*.png
python report.py                 # summary.md and figures/
python screen/analyse.py         # screen/screening.md and screen/screening.csv
python matched.py                # matched.md and presentation/8_*.png, 9_*.png

# Rerun one experiment, e.g. E1 with EMAC–VMS and default OSS (results/E1_emac_vms_oss):
python driver.py E1 emac_vms oss            # E1, E3, E4, E5: experiment form variant [T] [cells]
python driver.py E5 emac_vms oss_s0.5 3     # variant oss_s<s>: τ_m × s, τ_c unchanged
python driver_e2.py emac_vms 64             # E2: form cells [variant] [T]; MPI-safe
python screen/generate.py                   # screening configs, then run each with run_one.sh
```

The run data were produced at commit c8cf378 of `main` plus the `emac_vms_full` change, with
Firedrake 2026.10.0, PETSc 3.26 and Python 3.12. The dt-sweep configurations are in
`dtsweep/configs/` and run with `scripts/run_simulation.py`.
