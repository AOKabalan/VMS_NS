# Exact Gresho benchmark — prepared, not run

`gresho_exact` is separate from the existing `gresho` case and its Phase 12b
configurations. It reuses the existing analytical velocity, pressure and frame
indicator by importing them directly; their expressions and constants are unchanged.

The domain is (−0.5, 0.5)², the default mesh is 48×48 squares split into triangles,
ν = 0 exactly, and f = 0. The Euler fields are therefore an exact steady reference.
The case rejects positive viscosity rather than silently reverting to the old fallback.
`FreeSlip(0)` fixes only u_x on the vertical walls; `FreeSlip(1)` fixes only u_y
on the horizontal walls. Tangential velocity is free, including in the solver test
space. At corners both incident normal conditions apply. No outflow term is added.

These physical and time parameters follow §4.1 of
[Charnyi et al. (2017)](https://www.math.uh.edu/~molshan/ftp/pub/CHOR2017.pdf) and
§4.1 of [Charnyi et al. (2019), preprint](https://arxiv.org/pdf/1712.00857).
The papers use Taylor–Hood P2/P1; the OSS configurations below are additional VMS
comparisons, not reproductions of the papers' discretizations or linearization study.

All files in `configs/gresho_exact/` specify Δt = 0.01, T = 10 (1,000 steps),
Newton with frozen τ derivatives, stress viscosity (zero here), and direct LU/MUMPS.

| Configuration | Elements | Stabilization |
| --- | --- | --- |
| `emac_galerkin.toml` | P2/P1 | None; paper-style EMAC baseline |
| `skew_symmetric_1_galerkin.toml` | P2/P1 | None; matched SKEW1 baseline |
| `emac_vms_oss.toml` | P1/P1 | Dynamic nonlinear OSS, pressure subscale off |
| `skew_symmetric_1_oss.toml` | P1/P1 | Dynamic nonlinear OSS, pressure subscale off |
| `nonconservative_oss.toml` | P1/P1 | Same OSS settings |
| `conservative_oss.toml` | P1/P1 | Same OSS settings |
| `skew_symmetric_2_oss.toml` | P1/P1 | Same OSS settings |
| `emac_oss.toml` | P1/P1 | Same OSS settings; ordinary EMAC convection with VMS |
| `emac_vms_full_oss.toml` | P1/P1 | Same OSS settings; existing full-stress EMAC–VMS form |

`time.startup = "none"` explicitly selects Crank–Nicolson from step 1 for both
velocity scales. The first pressure time is 0.005; the first endpoint time is 0.01.
This new option is restricted to Crank–Nicolson. Existing `lower_order` and
`exact_history` startup behavior is unchanged, including their initial BDF1 step.
The new option is an explicit exception to the legacy startup rule in
`formulation_spec.md`; it uses the same midpoint coefficients and equations.
Initial velocity remains the nodal interpolant, with zero initial subscale.
The interpolant need not satisfy discrete continuity: without BDF1 damping, an
initial divergence defect can alternate between endpoints. No projection or damping
step has been inserted. Record this when interpreting early pressure/divergence data.

## Inviscid stabilization

The OSS files select the already implemented `pressure_subscale = false` variant.
This removes ζ and its weighted equation everywhere, as well as the associated
pressure-subscale term. With the option enabled, τ_c = 0 in the resting exterior
would leave ζ basis functions with zero weighted mass and a singular Jacobian.
This is an explicit numerical stabilization choice, not a viscosity floor, speed
regularization, or alteration of the Euler problem. It differs from Phase 12b,
which kept the pressure subscale and used ν = 1e-6. Velocity-residual OSS stabilization
and its consistent, unconstrained ξ projection remain active.

Dynamic subscales use τ_t = (α₀/Δt + τ_m⁻¹)⁻¹, finite at rest, where τ_m⁻¹ = 0.
Frozen-τ Newton avoids differentiating |a| at a = 0. The existing algebra is used
unchanged, including τ_m⁻¹ in the dissipation. The constant-pressure null space and
zero-mean normalization remain active. No artificial pressure pin or physical
viscosity is introduced.

Static VMS has unbounded τ_m at rest; dynamic ASGS can cancel the resolved mass
there. These variants are not prepared for this benchmark. Do not turn the pressure
subscale back on or change to those variants expecting the same singularity handling.
Initial-matrix checks do not establish convergence of the nonlinear evolution.

## Diagnostics and interpretation

- `kinetic_energy` retains its resolved meaning, ½‖u_h‖².
- The new `total_kinetic_energy` evaluates ½‖u_h + ũ‖² directly, including the
  cross term, with the common quadrature rule and both fields at the endpoint.
  It requires dynamic VMS to avoid mixing endpoint and midpoint fields.
- `subscale_kinetic_energy` remains ½‖ũ‖². For converged OSS with matching
  quadrature, orthogonality makes total energy equal to resolved plus subscale
  energy up to iteration error; the new diagnostic does not assume that equality.
- `linear_momentum` and `angular_momentum` retain their resolved-velocity meanings.
  These are the appropriate OSS balance quantities: its time term uses u_h,
  and ∫ũ = ∫x×ũ = 0 because constants and rigid rotations belong to its unconstrained
  projection space on this nonperiodic box. Endpoint extrapolation preserves those
  moments. For ASGS the theoretical moments would instead use u_h + ũ; no ASGS
  configuration is supplied here.
- EMAC/EMAC–VMS `pressure_error_l2` already uses physical p_h = P_h + ½|ū_h|²,
  where ū_h is the converged midpoint velocity. The error subtracts the spatial
  mean of p_h − p_exact using the error quadrature, so the analytical pressure's
  arbitrary constant does not affect it. Neither the pressure reconstruction nor
  this diagnostic was changed. VTK pressure also contains the physical correction.
- Existing `energy_budget`, `momentum_budget`, and measured convective defects
  retain their definitions. With this startup, midpoint budgets are populated
  from step 1. The OSS budget's energy is ½(‖u_h‖² + ‖ũ‖²), as before.
- Free-slip walls can exert normal forces and torques. Rigid translations and
  rotations are not generally admissible velocity tests on this square. Momentum
  drift and the momentum-budget residual therefore need not vanish; the latter
  reports the consistent wall reaction. The wall `force` diagnostic instead measures
  physical stress traction, with the existing fluid-on-body sign convention.
  Interior-only predicted momentum defects are deliberately not requested, because
  their zero-boundary qualifications need not hold with free tangential traces.
  The retained frame norms show whether numerical motion reaches the near-wall strip.

Each run writes to `results/gresho_exact/<configuration-stem>/`. It reports every
step and exports velocity/physical pressure every 100 steps. `initial_state.json`
contains initial resolved energy and moments; initial total energy is identical
because ũ⁰ = 0. The old Gresho analysis/queue remains separate.

## Manual launch

From the repository root, for one VMS run:

```bash
source /home/ali/venv-firedrake/bin/activate
OMP_NUM_THREADS=1 mpiexec -n 4 python scripts/run_simulation.py configs/gresho_exact/emac_vms_oss.toml
```

For the paper-style EMAC baseline, replace the file with
`configs/gresho_exact/emac_galerkin.toml`. To run every configuration sequentially:

```bash
source /home/ali/venv-firedrake/bin/activate
bash scripts/run_validation_queue.sh --configs configs/gresho_exact --runner simulation
```

Add `--dry-run` to that queue command to preview only. The queue uses four MPI
ranks per run. Individual runs refuse to overwrite an existing output directory.
No benchmark, prefix, Newton update, or time step was executed during preparation.
Validation is limited to configuration/coefficient checks, boundary traces, initial
residual/Jacobian assembly and prescribed-field diagnostics; runtime and evolution
remain unmeasured.

## Remaining VMS convection forms

After the existing EMAC–VMS and SKEW1 runs, launch the other five supported forms:

```bash
source /home/ali/venv-firedrake/bin/activate
bash scripts/run_gresho_exact_remaining_vms.sh --dry-run
bash scripts/run_gresho_exact_remaining_vms.sh
```

The script runs nonconservative, conservative, SKEW2, ordinary EMAC with VMS,
and full-stress EMAC–VMS sequentially, using four MPI ranks each (`--ranks N`
overrides this). It excludes `emac_vms_oss`, `skew_symmetric_1_oss`, and both
Galerkin configurations. The five new files differ from `emac_vms_oss.toml`
only in run name and convection form.

Any existing output directory is skipped, including an active or interrupted run;
no results are moved or overwritten. A failed simulation is recorded and the queue
continues to the next form, returning a nonzero status if any new run failed.
Ctrl-C stops the queue. Start/exit records go to
`results/gresho_exact/remaining_vms_queue.log`; each simulation retains its usual
logs and diagnostics. This script does not wait for the separately running SKEW1
process; launch it after SKEW1 finishes if you want to avoid competing for cores.
Only syntax, configuration equality/validation and dry-run checks were performed
for this addition; no additional simulation was launched during preparation.
