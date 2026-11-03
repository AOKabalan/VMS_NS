# Four-wave cloud visualization pilot

This is the bounded visualization run requested on 3 November 2026 after the
coupled Boussinesq manufactured-solution recheck. It exercises the existing
monolithic velocity/pressure/temperature/moisture solver and applies the B25
saturation proxy to saved fields. The full B7.3 stability matrix, B7.4 dry-budget
matrix and B7.5 resolution/parameter studies remain separate.

## Recorded configuration

The runnable configuration is
`configs/boussinesq/kh_cloud_visualization_pilot.toml`.

| Setting | Value |
| --- | --- |
| Domain | Four wavelengths, length 56.5288826557, height 20 |
| Re, Ri, Pr, Sc | 1000, 0.1, 1, 1 |
| Mesh | 96×64 P1 cells; sinh clustering 1.5 |
| Initial perturbation | epsilon 0.02, relative noise 0.01, PCG64 seed 1729 |
| Stabilization | Static linear OSS |
| Scalar convection | Conservative for temperature and moisture |
| Time | BDF2, dt 0.2, final time 80 |
| Moisture | RH peak 0.97, width 1.5, alpha 0, gamma 0.03 |
| Field output | Every five steps: t=1,2,...,80 |
| Execution | Two local MPI processes, direct MUMPS solves |

The mesh is deliberately coarse and the perturbation is 200 times the proposed
small-amplitude growth-test seed. This shortens the wait for visible roll-up.
The run uses a cheaper verified stabilization variant and a shorter final time
than the full application proposal. These choices are explicit pilot settings;
they do not replace the planned defaults or establish quantitative convergence.
No force holds the initial tanh profile steady.

## Visualization and diagnostics

The conserved tracer q is unchanged. Post-processing forms
`q_c = max(q - exp(alpha*theta - gamma*z), 0)` and cell vorticity
`d_x w - d_z u`. With alpha=0 this is a height-dependent diagnostic saturation
proxy; it does not remove moisture, release latent heat or simulate cloud
microphysics. The images show a two-dimensional cross-section.

The serial post-processor gathers all MPI VTK pieces and verifies the cell area
against the complete domain. It exports a combined VTK series with velocity,
temperature, q, nodal q_c and cell vorticity. Original pressure remains in the
solver's separate field series. Animation frames use the same cloud color scale
and show the central |z|<=5 portion of the height-20 channel. Diagnostics use
the whole channel and actual saved times; there is no fabricated initial frame.

CSV records resolved KE, vertical KE, potential energy, temperature/moisture
content, unmodified tracer extrema, nodal and quadrature cloud maxima, cloud
fraction, and its quadrature sensitivity. The cloud fraction uses raw
`q-exp(alpha*theta-gamma*z)>1e-6` at triangle Duffy–Gauss points, comparing 8
and 12 nodes per reference direction. Nodal q_c is only a display field.
The offline path explicitly supports P1 output. Its snapshot integrals do not
replace the solver's actual BDF history terms in budget identities.

For alpha=0, the continuum initial tracer is bounded above by
`0.97*exp(0.03^2*1.5^2/4) = 0.970491187...`. The source-free continuum
advection–diffusion problem preserves that bound. The discretization is not
bound-preserving, so the history plot shows this bound beside the unmodified
tracer extrema. Cloud morphology must be read alongside any reported overshoot.

## Checks and reproduction

Two analytic tests check saturation that exists inside a P1 triangle while its
nodal cloud values vanish, exact solid-body vorticity, P1 scalar/energy integrals,
and cloud fraction against an analytic height threshold. Both passed. The
saved-field KE, vertical KE, temperature and moisture integrals also match the
assembled FE diagnostics within 1e-11 tolerance in serial and two-rank timing
runs. These checks exercise the global gathered output rather than summing a
single rank's local arrays. Ruff passes for the new code.

From the repository root with the Firedrake environment active:

```bash
OMP_NUM_THREADS=1 mpiexec -n 2 python scripts/run_simulation.py \
  configs/boussinesq/kh_cloud_visualization_pilot.toml
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 python scripts/visualize_kh_clouds.py \
  results/bousinesq/cases/kelvin_helmholtz_shear_layer/b75_cloud_visualization_pilot \
  --output results/bousinesq/clouds/four_wave_pilot
```

The simulation refuses to overwrite an existing run. For a repeat, copy the
TOML and select a new run_name. `--allow-incomplete` on the post-processor is
only for an explicitly labelled in-progress preview. The solver, isothermal
defaults and previous output behavior were not changed for this work.

## Completed pilot results

All 400 steps converged on two MPI processes, reaching t=80 in
1659.88 s (27.66 minutes). Maximum sampled advective CFL was 1.00709.
The complete [measured summary](kelvin_helmholtz_cloud_pilot_results.json)
records parameters, solver provenance, extrema and output checks.

| Observation/check | Measured value |
| --- | ---: |
| Frames | 80, t=1,...,80 |
| Peak vertical KE, evaluated every step | 13.32430 at t=45.8 |
| Saved frames with nonzero cloud fraction | t=28,...,35 |
| Peak nodal cloud proxy | 0.0186754 at t=31 |
| Peak whole-domain cloud fraction | 0.000453162 (0.0453162%) at t=32 |
| Maximum 8-vs-12 quadrature fraction difference | 5.76498e-6 absolute |
| Moisture nodal range over saved frames | [-0.00132221, 1.00215336] |
| Initial continuum tracer upper bound | 0.97049119 |
| Maximum relative moisture-content drift, all steps | 6.99441e-15 |
| Maximum absolute moisture-content budget residual | 5.50948e-15 |
| Maximum absolute tested coupled energy residual | 2.83576e-14 |
| Largest saved-integral vs assembled-FE difference | 6.08935e-12 |

Four repeating billows are visible, with clear overturning by approximately
t=31 and four rolled temperature/vorticity structures retained at t=80. No
completed pairing is apparent in the saved frames. These are observations of
this coarse trajectory, not converged roll-up or pairing times. Clouds are
absent from the saved frames t=36,...,80.

The peak cloud point is x=20.60949, z=0.8857483 at t=31. There q=0.9924529
and q_sat=0.9737775. **The peak cloud value relies on numerical tracer
overshoot:** even the initial continuum upper bound 0.9704912 lies below
q_sat at that point. Exact source-free incompressible transport could not
produce that peak there with the specified data. Small negative tracer values
also develop later. Neither was clipped. Conservation/budget closure does not
imply preservation of tracer bounds or physical accuracy of the cloud pattern.
The cloud visualization workflow is exercised; physical cloud morphology needs
spatial/time refinement and tracer-bound checks. The planned gamma sensitivity
and full dry-budget studies were not launched to force a different picture.

Final artifacts are under `results/bousinesq/clouds/four_wave_pilot/`:

- [Compact GIF preview](../media/kelvin-helmholtz-clouds.gif) and
  [MP4 video](../media/kelvin-helmholtz-clouds.mp4). The MP4 preserves all 80 frames
  and eight-second playback with labelled simulation times.
- [Peak-cloud snapshot](../media/kelvin-helmholtz-clouds.png)
  and [history plots](kelvin_helmholtz_cloud_history.png).
- `clouds.pvd` with 80 combined VTU frames, `cloud_diagnostics.csv`,
  `summary.json` and all PNG frames. Raw velocity/pressure/T/q series remain
  under the case's `fields/` folder.

The GIF and MP4 frame counts, output time sequences, VTK arrays, global
integrals and scalar/energy residuals were checked. Peak, late-time and history
figures were visually inspected. Raw results are gitignored; selected media,
the plotting code, configuration, report and measured summary are versioned.
[Media provenance](../media/README.md) records source runs and preview encoding.
