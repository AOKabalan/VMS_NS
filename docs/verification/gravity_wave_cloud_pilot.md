# Stable gravity-wave cloud visualization

Completed the owner-requested companion to the KH case on 3 November 2026.
The new case produces stable buoyancy oscillations with a passive moist layer,
illustrating a mechanism associated with undulatus ripple patterns. It remains
an idealized two-dimensional saturation proxy, without ice microphysics or
latent heating; it does not establish a cirrocumulus cloud classification.
The [formulation](../gravity_wave_formulation.md) derives the background,
boundary conditions, standing-wave reference and moisture profile.

## Accepted configuration and results

`configs/boussinesq/gravity_wave_cloud_pilot.toml` uses four wavelengths,
P2 velocity/pressure/temperature/q on 64×16 affine cells, dynamic linear OSS
and skew-symmetric scalar convection. There is no mean shear or forcing to
sustain the wave. Γ=β=k=1, H=2π, A=0.2, momentum/thermal diffusion 0.01,
moisture diffusion 0.0002, initial peak RH 0.98, moisture width 1.2 and
saturation height rate 0.3. Fixed-temperature walls maintain the stable
background, while moisture is insulated. Pressure retains its null space.

BDF2 uses dt=0.08781018413800909 to t=21.0744441931222: 240 steps covering
three reference periods. All steps converged on two local MPI processes in
263.40 s (4.39 minutes). Maximum sampled advective CFL was 0.079392. Field
output every two steps gives 120 frames, with no invented t=0 frame.

| Check | Small-kick reference run | Four-wave moisture pilot |
| --- | ---: | ---: |
| Kick A | 0.01 | 0.2 |
| Mesh, P2 | 16×12, one wavelength | 64×16, four wavelengths |
| Linear reference period | 7.02481473 | 7.02481473 |
| Measured signed-mode period | 7.03926421 | 7.04154542 |
| Relative period error | 0.20569% | 0.23817% |
| Normalized velocity-mode RMS error | 1.33376% | 1.53318% |
| Normalized displacement-mode RMS error | 1.32447% | 1.52870% |
| Minimum dT/dz across saved fields | 0.995056 | 0.904591 |

Both runs pass the specified 2% period, 5% modal RMS and positive-gradient
checks. Periods use same-direction zero crossings of the signed velocity
mode; kinetic energy has twice the oscillation frequency and is not used for
this measurement. The displacement is inferred from the temperature anomaly
relative to Γz. P2 gradient extrema are checked at every parent cell vertex.
These are representative checks at one resolution/time step, not a full
convergence study. The full nonlinear equations are solved; the damped linear
reference is not declared an exact nonlinear solution.

Cloud bands appear and fade over six alternating half-cycle pulses. All
sampled cloud area lies in upward-displaced air, using an independent 10×10
triangle quadrature check. The first detected cloud frame is t=0.5268611.
Maximum nodal q_c is 0.03626568 at t=1.7562037; maximum whole-domain cloud
fraction is 0.01917760 (1.91776%). The four periodic repeats were prescribed
in the initial kick; the calculation is not an independent prediction of
their number or wavelength. Portions at the two plot edges join periodically.

## Bounds, quadrature and data checks

The initial continuum q upper bound is 1.00362073. Across saved frames,
nodal q lies in [-2.61369e-5, 1.00475228] and quadrature samples in
[-2.79298e-5, 1.00470396]. Maximum bound excess is 0.00113155. This transport
scheme is not bound-preserving, and q was not clipped. At the strongest cloud
point, q=0.97906214 and q_sat=0.94279646, below the initial upper bound: that
peak does not require violation of the continuum tracer maximum. The remaining
small transport errors still require refinement before quantitative cloud
morphology is trusted.

Cloud fraction uses the raw original quadratic field, not the displayed nodal
q_c interpolation. Its maximum absolute change between 8×8 and 12×12
Duffy–Gauss rules is 0.000324894 of the whole domain, or 0.03249 percentage
points. Report this sensitivity alongside coverage. The scalar threshold is
q_c>1e-6. Display images subdivide each P2 triangle into four P1 triangles;
native VTK nodal q_c can overshoot between P2 nodes. These display fields do
not change q or the raw quadrature diagnostics.

The maximum relative moisture-content drift over the recorded solver steps,
relative to their first sample, is 1.95659e-7. Linear splitting with skew convection does not guarantee exact
content conservation per B19. The existing `thermal_content_budget` group
rejects fixed-temperature walls, so it was not requested and no unsupported
budget closure is claimed. Saved-field KE, vertical KE, T and q integrals
agree with assembled FE diagnostics within 1.71e-13 absolute. The offline
integrals gather the complete MPI mesh and verify its area.

The five targeted tests passed: linearized PDE/wall compatibility and solved
hydrostatic rest, plus three cloud-kernel tests covering P1 raw/nodal
saturation, analytic coverage, and P2 polynomial derivatives/integrals on
both cell orientations. With the existing mesh/boundary tests, 38 passed.
The old KH visualization command also regenerates its timing frame through
the shared implementation, with agreement against assembled diagnostics.
The final GIF/MP4/VTK frame counts, fields and time sequences are checked
separately; mode, cloud and late-cycle plots were visually inspected.

## Preliminary findings retained

The following setup trials were excluded from the accepted results:

- Coarse P1 hydrostatic controls: quadratic pressure was not representable;
  static stabilization failed its first-step solve, and dynamic stabilization
  produced startup velocity instead of preserving rest. P2 passes the rest
  check. No force was subtracted to disguise this issue.
- μ=0.0002 P2 wave runs, with conservative and skew scalar convection:
  early oscillations follow theory, but later short-scale temperature
  oscillations give negative gradients. Their completed solves do not pass
  stable-wave verification. μ=0.01 defines the more damped physical pilot;
  no general low-diffusion stability claim follows from it.
- An intermediate μ=0.01 run had the thermal and moisture diffusivities
  swapped in the case builder. Saved physical metadata exposed the mistake.
  It was corrected before the accepted verification/pilot. The case test now
  checks coefficient wiring, and the period analyzer rejects unequal
  momentum/thermal diffusion for this reference.
- An early timing request included the unsupported content-budget group at
  fixed-temperature walls. It exited before stepping; the diagnostic request
  was corrected without changing the wall data.

Raw runs remain under `results/bousinesq/cases/internal_gravity_wave/`.
Accepted small-wave plots and rejected preliminary plots are indexed under
`results/bousinesq/clouds/gravity_wave_checks/`; their labels distinguish their
status. The [measured JSON](gravity_wave_cloud_pilot_results.json) contains
the accepted run provenance, configuration and diagnostics.

## Reproduction and artifacts

From the repository root with the Firedrake environment active:

```bash
OMP_NUM_THREADS=1 mpiexec -n 2 python scripts/run_simulation.py \
  configs/boussinesq/gravity_wave_cloud_pilot.toml
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 python scripts/analyse_gravity_wave.py \
  results/bousinesq/cases/internal_gravity_wave/gravity_wave_cloud_pilot \
  --output results/bousinesq/clouds/gravity_wave_pilot
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 python scripts/visualize_clouds.py \
  results/bousinesq/cases/internal_gravity_wave/gravity_wave_cloud_pilot \
  --output results/bousinesq/clouds/gravity_wave_pilot
```

The simulation refuses to overwrite a run; select a new run_name for repeats.
The solver equations, isothermal defaults and existing FE outputs were not
changed. The new case is registered, with shared P1/P2 post-processing and
the original `visualize_kh_clouds.py` command retained for compatibility.

Final outputs are under `results/bousinesq/clouds/gravity_wave_pilot/`:

- [GIF](../../results/bousinesq/clouds/gravity_wave_pilot/clouds.gif) and
  [MP4](../../results/bousinesq/clouds/gravity_wave_pilot/clouds.mp4), 120 frames,
  12-second playback with actual dimensionless time labels.
- [Cloud snapshot](../../results/bousinesq/clouds/gravity_wave_pilot/maximum_cloud.png)
  and [wave/reference comparison](../../results/bousinesq/clouds/gravity_wave_pilot/wave_modes.png).
- Combined `clouds.pvd`/VTU series, PNG frames, CSV diagnostics, cloud-history
  plots, wave-mode CSV and both measured summaries.
