# G4 periodic/free-slip boundary verification

Date: 2 November 2026. Implementation:
[boundary prerequisite G4](../boussinesq_implementation_plan.md#boundary-prerequisite-g4--before-the-early-b71-test).
Measured CLI data: [JSON](free_slip_g4_results.json).

The full unit suite passed **983 tests with 1 known xfail in 448.70 s**.
The xfail remains the existing near-inviscid nonlinear ASGS limitation.
The 48 new G4 tests and the added case-registry check account for the increase
from B5. No KH growth/cloud simulation is included in this phase.

## Mesh and essential constraints

The installed Firedrake 2026.10.0 constructor with `direction="x"` exposes
bottom/top tags 3/4; periodic tags 1/2 have no facets. The test translates a
channel to `[-0.75,1.75] × [-1.5,2.5]` and checks both wall lengths, heights and
horizontal normals. `periodic_channel()` returns `y_min=3`, `y_max=4`; the
existing fully periodic helper is unchanged.

`FreeSlip(normal_component)` is a homogeneous, axis-aligned boundary type.
Components must be integer indices within the geometric dimension. Validation
checks nonempty wall measure and that the mean squared off-axis normal is no
larger than `1e-20`. Incorrect components and tilted walls are rejected. A 3D
box test constrains only the normal component on its top/bottom faces.

A mixed-space test fills every field with nonzero values, applies the component
BCs, and checks that tangential velocity, pressure, temperature, tracer and OSS
projections stay unconstrained. Applying the homogeneous BCs to a residual
removes normal rows only. Jacobian actions confirm that constant pressure is
the null mode, while constant horizontal velocity and velocity projection are
not null modes. Pressure remains mixed component 1.

Initialization clamps only slip-normal boundary DOFs in resolved initial and
exact-history states. It leaves tangent values and interior normal values
alone; this is not a divergence-free projection. The mixed Dirichlet/slip flux
test rejects incompatible prescribed normal flow. A failed coupled solve keeps
the tangential state, writes velocity/pressure/T/q, and preserves endpoint
histories. Slip walls are never classified as outflow.

## Uniform preservation and natural traction

The uniform fixture is a 2×1 channel with `u=(1,0)`, pressure zero, no forcing,
and `nu=0.1`. Thermal runs add `T=T_ref=1`, tracer `q=2`, `kappa=0.1`, zero
sources and insulated walls. The scalar values are analytic preservation data.

The primary matrix checks Galerkin and all eight VMS variants under:

- isothermal Laplacian viscosity;
- isothermal stress viscosity;
- thermal Laplacian viscosity with passive q.

These 27 cases use three steps with BDF2 startup. Eight additional thermal OSS
checks cover P1/P2 and BDF1/BDF2/BDF3/CN, including scalar histories, VTK and
coupled budgets. Preservation assertions require velocity/scalar L2 errors
below `1e-11`, zero normal trace to tolerance, unchanged horizontal mean, and
the pressure gauge. They do not impose zero mean horizontal velocity.

A nonuniform, source-free shear test starts from
`u_x=1+0.2 cos(pi(y+0.5))`, `u_y=0`. The continuum amplitude decays as
`exp(-nu pi^2 t)` and has zero tangential traction at both walls. On the 8×8
P1 static-linear OSS mesh to `t=0.1`, both viscous forms give velocity error
below `0.005`, preserve the horizontal mean and normal velocity, and agree in
L2 below `1e-10`. This exercises an actual diffusion solve as well as constant
state preservation.

The thermal kinetic–potential budget is separately checked on decaying shear
with nonzero tangential wall speed. The velocity is an admissible wall test:
only its normal component vanishes, and tangential traction is natural zero.
The budget now checks full wall speed only on full-vector Dirichlet walls.
Its residual stays below `1e-10` in this test. Thermal stress viscosity retains
the B6 gate; G4 changes boundary support, not that formulation's compatibility.

## CLI, MPI and previous results

Both examples completed four steps, including BDF1 startup and CN half-levels:

| Configuration | Velocity L2 error | Pressure L2 error | T L2 error |
| --- | ---: | ---: | ---: |
| Thermal Laplacian + passive q | 0 | 0 | 1.835e-16 |
| Isothermal stress | 0 | 0 | — |

Kinetic energy stays 1, horizontal mean stays 1, and the thermal/tracer
contents stay 2/4. Maximum unnormalized thermal CLI budget residuals are below
`1e-17` for kinetic–potential energy and `1e-30` for scalar variance; content
residuals are zero at recorded precision.

The same configurations completed on two MPI ranks. Stress diagnostics are
identical to serial; the largest unnormalized thermal diagnostic difference is
`3.669e-18`, with zero reported velocity and pressure error on both partitions.
Near-zero budget terms use the documented `1e-12` normalization floor: their
normalized residuals can differ by `3.669e-6` while their absolute difference is
`3.669e-18`. This is roundoff relative to the floor, not a comparable error in
flow energy or scalar content. The JSON retains both comparisons.

Four isothermal baseline runs and three previous thermal runs have bitwise
identical arrays and diagnostic/inner/outer CSV records before/after G4.
Artifacts are under `/tmp/vms-g4-before-*` and `/tmp/vms-g4-after-*`.
The full regression log is `/tmp/vms-g4-all.log`.

## Reproduction

Use the project Firedrake environment:

```bash
python -m pytest tests/unit/test_free_slip.py -q
python -m pytest tests/unit -n 2 -q
python scripts/run_simulation.py configs/boussinesq/free_slip_uniform.toml
python scripts/run_simulation.py configs/boussinesq/free_slip_uniform_stress.toml
```

For MPI, copy each config with a fresh `run_name`, then invoke the same
`scripts/run_simulation.py` through `mpiexec -n 2`. The measured serial run names
are `g4_uniform_thermal` and `g4_uniform_stress`; MPI adds `_mpi`. Resolved configs,
run information, diagnostics and physical-field output are in
`results/free_slip_channel/<run_name>/`. The runner refuses existing folders.

G4 completes the boundary gate for the later stratified shear-layer application.
Its initial-condition construction, stability study and cloud proxy remain B7.
