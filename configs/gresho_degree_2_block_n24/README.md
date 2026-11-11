# Smaller degree-2 exact Gresho runs

Run from the repository root:

```bash
bash scripts/run_gresho_degree2_block_n24.sh
```

Four MPI ranks; sequential order: CONS, EMAC-VMS, CONV, SKEW1.
The configs use 24×24 cells instead of 48×48 (one quarter of the elements).
All original physical, time-integration, discretization, diagnostic, and nonlinear
convergence settings from `configs/gresho_degree_2/` are retained: P2/P2,
zero viscosity, free-slip walls, pure Crank–Nicolson, dt=0.01, T=10.

The assembled block preconditioner splits the OSS projection field xi (field 2)
from velocity and pressure (fields 0,1), with a full Schur factorization.
The projection block uses CG/Jacobi; the velocity–pressure block uses MUMPS LU.
The outer linear solver is FGMRES. Preconditioners are refreshed each Newton
iteration. MUMPS null-pivot detection is enabled on the velocity–pressure block
(`fieldsplit_1_mat_mumps_icntl_24=1`) to handle its pressure gauge singularity.
This changes the preconditioner, not the physical operator or boundary conditions.
This is a block solver with a direct fluid-block solve, not a fully iterative
velocity/pressure split.

New output folders are `results/gresho_degree_2/*_n24_block_gauge_oss/`.
The stopped 48×48 CONS run remains at
`results/gresho_degree_2/conservative_oss/`: diagnostics through step 795,
t=7.95, velocity snapshots through t=7. These are saved results, not a restart
checkpoint. The completed exact EMAC Galerkin run is also preserved.
These smaller simulations start from their own initial conditions.

Earlier unsuccessful block trials and short solver-check outputs were removed
after recording their outcomes and configuration differences in
[`troubleshooting.md`](../../troubleshooting.md).
The queue records statuses in `results/gresho_degree_2_block_n24_queue.log`;
console output is in `results/gresho_degree2_block_n24.console.log`.
The generic queue skips already finished runs, including recorded failures.

Validation before the full queue: CONS and EMAC-VMS each completed 50 steps
(t=0.5); CONV and SKEW1 each completed 10 steps (t=0.1), all with exit status 0.
Check settings and relevant console evidence are recorded in
[`troubleshooting.md`](../../troubleshooting.md).
