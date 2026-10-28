# Kelvin–Helmholtz: SKEW1 versus EMAC–VMS

Prepared to produce a Figure 4-style comparison from Olshanskii & Rebholz,
*Longer time accuracy for incompressible Navier–Stokes simulations with the
EMAC formulation* (2020), section 4.3:
[paper, pp. 15–17](https://arxiv.org/pdf/2002.01416#page=15).
The source was checked against `/home/ali/Downloads/longer_time_accuracy.pdf`.

The experiment asks whether SKEW1 still develops a large enstrophy increase
with dynamic OSS, and how EMAC–VMS compares. VMS may suppress the instability
seen with unstabilized SKEW in the paper; it does not guarantee its absence.
A smoother solution can also be more dissipative. Inspect vorticity, enstrophy
and energy together. These two VMS runs alone do not isolate the causal effect
of adding VMS: that would require an unstabilized control on the same setup.

## Setup

| Setting | Both prepared configurations |
|---|---|
| Domain | Unit square, periodic in x |
| Walls y=0,1 | Strong u_y=0; natural zero tangential traction; u_x is free |
| Initial shear thickness | delta_0 = 1/28 |
| Velocity scale | u_infinity = 1 |
| Perturbation | c_n = 0.001, streamfunction below |
| Reynolds number | Re = delta_0 u_infinity / nu = 1000 |
| Viscosity | nu = 1/28000, approximately 3.5714286e-5 |
| Mesh | 196 x 196 rectangles split into triangles |
| Velocity / pressure | Continuous P1/P1 (owner's revised choice) |
| VMS | Dynamic OSS, nonlinear splitting, pressure subscale on |
| Stabilization constants | c1=4, c2=2 |
| Viscous form | Stress form, for both formulations |
| Time | BDF2, dt=0.001, T=14; one BDF1 startup step |
| Nonlinear solve | Newton with frozen tau; maximum 200 inner iterations |
| Snapshots | t=0, 1, 4.5, 5, 5.5, 6, 6.5, 10 and 14 |

```text
psi = u_infinity exp(-(y-0.5)^2/delta_0^2) [cos(8 pi x) + cos(20 pi x)]
u_0 = (u_infinity tanh((2y-1)/delta_0), 0) + c_n (d_y psi, -d_x psi)
```

The paper uses unstabilized P2/P1 Taylor–Hood and ordinary EMAC. Here both
methods use equal-order P1/P1, selected after the P2/P2 timing prefix proved
too expensive on the owner's machine;
EMAC–VMS is a different stabilized formulation. The physical parameters,
boundary conditions, and main time-stepping parameters match the paper's stated
setup. The default 196x196 P1/P1 experiment uses the paper's mesh count but
lower polynomial degree than its P2/P1 calculation. It is an exploratory comparison, and any
late-time structures need a refinement check before an accuracy conclusion.
We use the repository's uniform triangulation;
identical triangle orientation and initial projection to the authors' code
are not established. We interpolate the initial velocity, enforce zero normal
wall data, and start subscales at zero.

The paper writes Laplacian viscosity, whereas EMAC–VMS requires stress
viscosity here. On flat impermeable walls their tangential free-slip conditions
agree, but the discrete volume forms differ by a divergence term. Both of our
runs use the same stress form. This is a VMS adaptation, not a claim to
reproduce the paper's exact curves. There is no exact time-dependent solution
or independently converged reference supplied for this experiment.

## Run from the repository root

```bash
cd /home/ali/02_projects/VMS_NS
source /home/ali/venv-firedrake/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export FIREDRAKE_CACHE_DIR=/tmp/vms_kh_cache
export MPLCONFIGDIR=/tmp/vms_kh_mpl

# Recommended first: ten steps of each full-size configuration.
# These go into separate *_prefix10 folders and do not overwrite full runs.
NP=4 bash experiments/kelvin_helmholtz_2020/run_pair.sh --max-steps 10

# Full experiment: 14,000 steps per formulation, 28,000 total.
NP=4 bash experiments/kelvin_helmholtz_2020/run_pair.sh

# After completion (also supports clearly marked partial output):
python experiments/kelvin_helmholtz_2020/plot.py
```

The two runs execute sequentially. A failure stops the shell queue and gives
a nonzero exit status. Existing output directories are never overwritten.
Change `--output-root` to repeat a run without removing results. To launch only
one formulation:

```bash
mpiexec -n 4 python experiments/kelvin_helmholtz_2020/run.py \
  experiments/kelvin_helmholtz_2020/configs/kh2020_emac_vms_p1_n196.toml
```

At the restored 196x196 P1/P1 resolution there are 77,224 velocity unknowns
and 231,672 mixed unknowns including OSS projections. The P2/P2 pair has
308,112 and 924,336. The optional 64x64 P1/P1 pair has 8,320 and 24,960.
A direct MUMPS solver is used. Use the prefix's
`time_steps.csv` and `run_info.json` to judge runtime and memory on your machine.
The first steps include compilation; later nonlinear roll-up may cost more.
No full-resolution runtime estimate is asserted from the small validation runs.

## LU reuse assessment

The default `LINEAR_SOLVER=direct` retains the original LU/MUMPS solver.
The separately named `*_reuse_lu.toml` configurations use right-preconditioned
GMRES on the current Jacobian, with an LU/MUMPS preconditioner rebuilt every
five Newton updates. The rebuild counter persists across inner nonlinear
solves. Jacobians are still refreshed every update; nonlinear and subscale
tolerances, spatial discretization, physical parameters and time steps are
unchanged. Linear solves use the unpreconditioned residual norm, rtol=1e-10,
atol=1e-14 and a maximum of 100 GMRES iterations.

This targets repeated numeric factorization, not the number of coupled
unknowns or the LU memory requirement. As the flow evolves, an older
factorization can require more GMRES iterations or fail to converge. A short
timing comparison cannot establish robustness through the vortex mergers.
See PETSc's documentation on
[preconditioner lagging](https://petsc.org/release/manualpages/SNES/SNESSetLagPreconditioner/)
and [persistence across solves](https://petsc.org/release/manualpages/SNES/SNESSetLagPreconditionerPersists/).

Generate both variants and run explicitly bounded comparisons in fresh output
directories (existing directories are protected):

```bash
python experiments/kelvin_helmholtz_2020/generate.py
python experiments/kelvin_helmholtz_2020/generate.py --linear-solver reuse_lu
NP=4 LINEAR_SOLVER=direct bash experiments/kelvin_helmholtz_2020/run_pair.sh \
  --max-steps 4 --output-root results/kh_lu_baseline \
  --petsc-profile --save-final-state
NP=4 LINEAR_SOLVER=reuse_lu bash experiments/kelvin_helmholtz_2020/run_pair.sh \
  --max-steps 4 --output-root results/kh_lu_reuse \
  --petsc-profile --save-final-state
```

`petsc_profile.txt` records event times and counts. `--save-final-state` writes
owned arrays of all coupled fields and the local subscale, one NPZ per MPI
rank. New exports also include FE/subscale coordinates and subscale cell
centres, so arrays can be matched by physical location across different MPI
partitions. Comparisons must identify periodic x=0 and x=1 nodes. These
arrays are comparison data, not restart checkpoints. The run name records
the solver variant, so plotting a reuse run requires explicit `--emac` and
`--skew` paths. Neither option changes the normal default output.

The owner-approved four-step assessment on 2026-10-25 completed successfully
for both forms and both solvers on four MPI ranks. Steps 3–4 averaged:

| Form | Direct LU | Reused LU + GMRES | Speedup |
|---|---:|---:|---:|
| SKEW1 | 9.92 s/step | 3.55 s/step | 2.79x |
| EMAC–VMS | 11.85 s/step | 5.00 s/step | 2.37x |

Across all four steps, numeric LU factorizations fell from 25 to 5 for SKEW1
and from 28 to 6 for EMAC–VMS. Outer and Newton iteration counts were unchanged.
Final velocity nodal values agree within 1.62e-15 absolute, and the largest
relative coefficient-vector difference over all coupled fields and subscales
is 3.56e-11 (the divergence projection). Energy differences are at most
5.56e-17 absolute and enstrophy differences at most 2.85e-14 over the prefix.
The reuse solves required at most 13 GMRES iterations per Newton update.
The comparison covers only t<=0.004, not the late-time vortex dynamics.
Both variants used the same snapshot and profiling options; the first two
steps are excluded from speedup estimates to reduce startup/cache effects.
Saved evidence: `results/kh_lu_assessment/README.md`, `comparison.json`,
`compare.py`, and the four run folders with PETSc profiles and final arrays.

To launch the faster variant later, explicitly select it and a fresh output
root (this does not resume the interrupted direct run):

```bash
NP=4 LINEAR_SOLVER=reuse_lu bash experiments/kelvin_helmholtz_2020/run_pair.sh \
  --output-root results/kh_lu_production
python experiments/kelvin_helmholtz_2020/plot.py \
  --emac results/kh_lu_production/kelvin_helmholtz/kh2020_emac_vms_p1_n196_reuse_lu \
  --skew results/kh_lu_production/kelvin_helmholtz/kh2020_skew1_p1_n196_reuse_lu \
  --output results/kh_lu_production/comparison
```

## Detailed profiling and a further solver option

`--profile-phases` enables PETSc logging with separate `KH_initial_steps`
(steps 0–2) and `KH_measured_steps` (steps 3 onward) stages. It also writes
`phase_timings.json`: per-rank, per-step inclusive and exclusive wall times
for the global solve, local subscale updates, change/residual checks,
diagnostics and history updates. Inclusive times nest; only exclusive
times can be added within a rank. There are no profiling barriers inside
the iteration loop. Setup and final-state export stay outside measured
time-step stages. Instrumentation is confined to this experiment driver.

Ten-step profiling on four ranks at 196x196 P1/P1 found that the reuse-5
baseline spends about 91% of measured step time inside SNES. Numeric LU
factorization takes about 42%, factor applications 29%, Jacobian assembly
17%, local subscales 4–5%, and diagnostics plus snapshots under 1%.
These percentages describe steps 3–10; PETSc event fractions use
maximum-rank event times and are approximate.

The optional `reuse_lu_tuned` configuration combines LU reuse for 20 Newton
updates with Eisenstat–Walker version 2 adaptive linear tolerances
(`rtol0=0.1`, `rtolmax=0.5`). This changes intermediate linear stopping
tolerances. The original final inner/outer nonlinear and local-subscale
convergence tests, current Jacobian, full residual, mesh, physical data,
stabilization and time scheme are retained.

| Solver | SKEW1 seconds/step | EMAC–VMS seconds/step |
|---|---:|---:|
| Reuse 5 (existing faster baseline) | 4.18 | 5.97 |
| Reuse 20 only | 3.57 | 5.01 |
| Reuse 5 + adaptive linear tolerance | 3.98 | 5.83 |
| Reuse 20 + adaptive linear tolerance (`reuse_lu_tuned`) | 2.81 | 4.33 |

All averages exclude steps 1–2. A repeated baseline averaged 4.08 and
6.26 s/step, respectively. Against that baseline, coordinate-matched final
velocity values differ by at most 8.39e-13 absolute; all coupled/subscale
coefficient vectors differ by at most 3.13e-10 relative. Energy histories
differ by at most 1.78e-15 and enstrophy by at most 1.24e-12 absolute.
All final coupled residuals in the combined candidate are below 1e-12.
The option passed these early-time checks through t=0.01; robustness at
the later vortex mergers remains unmeasured. No full simulation was started.

The profile, comparisons and scripts are under `results/kh_phase_profile/`.
To generate and select the further tuned option explicitly:

```bash
python experiments/kelvin_helmholtz_2020/generate.py --linear-solver reuse_lu_tuned
NP=4 LINEAR_SOLVER=reuse_lu_tuned bash experiments/kelvin_helmholtz_2020/run_pair.sh \
  --max-steps 10 --profile-phases --save-final-state \
  --output-root results/kh_tuned_check
```

For a later full run, omit `--max-steps 10` and choose a fresh output root.
Run names end in `_reuse_lu_tuned`, so use explicit run paths when plotting.
The default remains the direct reference solver. Larger speedups would
require further solver work; the profile identifies factor applications
and Jacobian assembly as the leading costs after this tuning.

## Experimental OSS block preconditioning

`generate.py --linear-solver block_oss` generates a further opt-in solver for
this dynamic nonlinear OSS experiment. The full operator retains all 231,672
unknowns at 196x196 P1/P1 and the original consistent weighted projections.
Only the preconditioner is split: fields 2,3 (xi,zeta) form a 115,836-unknown
weighted mass block, and fields 0,1 (u,p) form a 115,836-unknown fluid block.
Flexible GMRES acts on the full current Jacobian with a right Schur block
preconditioner. Projection-block solves use diagonally preconditioned CG
(relative tolerance 1e-2, absolute tolerance 1e-14, maximum 50 iterations).
The fluid Schur preconditioner uses the fluid diagonal block A11 with
LU/MUMPS. The Schur factorization is `full`; all off-diagonal couplings remain
in the actual operator and in the appropriate block correction operations.
See [PETSc field splitting](https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/).

As with `reuse_lu_tuned`, the preconditioner refreshes every 20 Newton
updates and uses Eisenstat–Walker version 2 with initial/max relative linear
tolerances 0.1/0.5. The outer linear iteration cap is 200. Final nonlinear,
subscale and coupled-residual acceptance criteria are unchanged. Approximate
inner mass solves affect preconditioning, not the discrete equations.

Six four-step screens compared different Schur/mass-block approximations.
The selected option and the tuned whole-system LU reference were then run
twice for each formulation, ten steps per run, reversing execution order.
Mean seconds/step over steps 3–10:

| Form | Tuned LU, rerun | Block | Time reduction |
|---|---:|---:|---:|
| SKEW1 | 4.17 | 3.50 | 16% |
| EMAC–VMS | 5.89 | 4.88 | 17% |

Summed per-rank peak RSS averages drop from 3.12 to 2.19 GiB for SKEW1 and
3.10 to 2.22 GiB for EMAC–VMS (rank peaks need not occur simultaneously).
These are paired measurements from the block assessment, not comparisons
against older runs under different machine conditions. Individual reference
timings varied; two repetitions do not establish a statistical confidence
interval or late-time performance.

Coordinate-matched final velocity differences are below 4.69e-12 absolute.
All final coupled/projection/subscale coefficient-vector differences are
below 4.08e-9 relative, energy-history differences below 2.23e-15 absolute,
and enstrophy-history differences below 1.15e-12. All terminal coupled
residuals are below 1e-12. Validation covers t<=0.01 on this P1/P1 case;
later vortex dynamics and other discretizations remain untested.

To select this experimental option explicitly for another short prefix:

```bash
python experiments/kelvin_helmholtz_2020/generate.py --linear-solver block_oss
NP=4 LINEAR_SOLVER=block_oss bash experiments/kelvin_helmholtz_2020/run_pair.sh \
  --max-steps 10 --profile-phases --save-final-state \
  --output-root results/kh_block_check
```

Full run names end in `_block_oss`; plotting requires explicit run paths.
The default solver is unchanged. The benchmark data, screening results,
comparison scripts and figures are in `results/kh_block_profile/`.
Profiling now also writes `solver_view.txt` and `memory_peaks.json`.
The SNES view displays PETSc's built-in tolerances; this application installs
a custom convergence test, so the resolved configuration and recorded
nonlinear/coupled residual checks determine acceptance.

## Outputs

Each run is under `results/kelvin_helmholtz/kh2020_<form>_p1_n196/`:

- `config_resolved.toml`, `experiment.json`, `run_info.json`: configuration,
  actual viscosity, degrees of freedom, software/code provenance and status.
- `diagnostics.csv`: resolved and subscale energy, linear and angular momentum,
  divergence and convective energy defect. No CN energy-closure residual is
  reported for BDF2.
- `enstrophy.csv`: Z=1/2 integral curl(u_h)^2, maximum absolute vorticity,
  and wall-normal velocity norm, including the actual t=0 state.
- `mesh.npz`, `snapshot_*.npz`: triangle-local geometry and signed resolved
  vorticity. P1 velocity has constant curl on each triangle, stored at the three
  DG1 vertex nodes. P2 remains supported with linear curl. Neither is averaged
  across elements. Files are gathered on rank 0.
- Standard nonlinear iteration and step logs. Diagnostics are flushed and
  each snapshot is written when reached, so completed data survive a later failure.

The plotter writes to `results/kelvin_helmholtz/comparison_p1_n196/`:

- `figure4_vorticity.png` and `.pdf`: EMAC–VMS left, SKEW1 right, at the seven
  Figure 4 times; absolute resolved vorticity with one common color scale.
- `histories.png` and `.pdf`: energy, enstrophy, maximum vorticity, momentum,
  angular momentum, and divergence histories.
- `summary.json`: completion status, final energy/enstrophy fractions, the
  enstrophy peak in 5<=t<=7, and unavailable snapshots.

The default color scale saturates at |omega|=80; use `--vmax` to change it.
The maximum-vorticity history remains unclipped. Missing snapshots are labelled,
never replaced with another time. Angular momentum is a diagnostic in this
periodic channel, not an asserted global invariant. Field errors are not reported
without an exact or independently validated reference.

Before plotting, the saved triangles must cover the unit domain and each
snapshot must reproduce its independently assembled enstrophy using the exact
DG1 triangle mass integral. The plotter rejects inconsistent exports. Snapshot
files are published by atomic rename so an ongoing write is not read as a frame.

## Smaller pilot or stabilization comparison

To make an optional still smaller 32x32 pilot (this changes spatial resolution):

```bash
python experiments/kelvin_helmholtz_2020/generate.py --cells 32 \
  --output results/kh_configs_n32
mpiexec -n 4 python experiments/kelvin_helmholtz_2020/run.py \
  results/kh_configs_n32/kh2020_skew1_p1_n32.toml
mpiexec -n 4 python experiments/kelvin_helmholtz_2020/run.py \
  results/kh_configs_n32/kh2020_emac_vms_p1_n32.toml
python experiments/kelvin_helmholtz_2020/plot.py \
  --emac results/kelvin_helmholtz/kh2020_emac_vms_p1_n32 \
  --skew results/kelvin_helmholtz/kh2020_skew1_p1_n32 \
  --output results/kelvin_helmholtz/comparison_p1_n32
```

`generate.py --tau-scale 0.5 --output results/kh_configs_tau05` creates an
optional pair with c1 and c2 both divided by 0.5. At fixed advector this scales
tau_m by 0.5 while preserving the functional formula for tau_c. Since the
computed advector can change, tau_c's realized field need not be identical.
Only compare both forms at matched settings; do not tune just one until a
desired instability appears. Refinement and time-step sensitivity are needed
before interpreting a new enstrophy increase as a robust formulation effect.

## Preparation checks

- 190 targeted tests passed: new KH geometry/initial-condition/free-slip tests,
  mesh and boundary tests, time-stepper regressions, validation cases and configuration tests.
- Both P2/P2 formulations completed four steps on a 16x16 mesh, in serial and
  on two MPI ranks. BDF1 startup followed by BDF2 was confirmed in the logs.
- Serial/MPI energy and enstrophy agree within 2.82e-15 relative in those runs.
- Wall-normal velocity is zero to the checked tolerance (1e-12).
- Saved triangles cover the domain and the snapshot mass integral reproduces
  the independently assembled enstrophy within 1e-10 relative.
- PNG/PDF rendering, short-prefix output names, and missing-frame labels were
  exercised. These are preparation checks, not validation of the late-time KH
  dynamics.

The owner selected P1/P1 and temporarily reduced the default to 64x64 after
stopping the expensive 196x196 run. The original 196x196 mesh is now restored,
retaining P1/P1. Old configurations and results are preserved. The owner
approved only four-step comparisons of the direct and LU-reuse solvers for
both formulations; full simulations remain stopped.
