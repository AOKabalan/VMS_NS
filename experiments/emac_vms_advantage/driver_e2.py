"""E2 thin double shear layer: python driver_e2.py <form> <cells> [variant] [T].

Doubly periodic unit square, u = tanh(ρ(y − ¼)) for y ≤ ½ and tanh(ρ(¾ − y)) above,
v = δ sin(2πx), ρ = 80, δ = 0.05, ν = 1e-4 (Bell–Colella–Glaz / Minion–Brown type; the
constants are not checked against those papers, so the judge is our own fine-mesh run).
Writes the solver outputs, extra.csv (energy, enstrophy, momentum) and vorticity snapshots.
"""

import csv
import logging
import math
import sys
from dataclasses import replace

import firedrake as fd
import numpy as np
import ufl
from driver import HERE, make_config

from navier_stokes_vms.config import MeshOptions, OutputOptions, TimeOptions, validate_config
from navier_stokes_vms.postprocessing.diagnostics import check_diagnostics_supported
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.problem_definition import Problem, at_initial_time
from navier_stokes_vms.solver.time_stepper import TimeStepper

RHO, DELTA, VISCOSITY = 80.0, 0.05, 1e-4
SNAPSHOTS = (0.8, 1.0, 1.2, 1.5)
DIAGNOSTICS = (
    "kinetic_energy",
    "subscale_kinetic_energy",
    "linear_momentum",
    "energy_budget",
    "momentum_budget",
    "linear_momentum_defect",
    "convective_energy_defect",
    "divergence_l2",
)


def shear_layer(cells: int) -> Problem:
    mesh, tags = mesh_loading.periodic_rectangle((cells, cells), (0.0, 0.0), (1.0, 1.0))

    def velocity(x: ufl.core.expr.Expr, _t: object) -> ufl.core.expr.Expr:
        lower = ufl.tanh(RHO * (x[1] - 0.25))
        upper = ufl.tanh(RHO * (0.75 - x[1]))
        return ufl.as_vector(
            (ufl.conditional(x[1] <= 0.5, lower, upper), DELTA * ufl.sin(2 * math.pi * x[0]))
        )

    return Problem(
        mesh=mesh,
        boundary_conditions={},
        boundary_tags=tags,
        viscosity=VISCOSITY,
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=at_initial_time(velocity),
    )


def main() -> int:
    form, cells = sys.argv[1], int(sys.argv[2])
    variant = sys.argv[3] if len(sys.argv) > 3 else "oss"
    t_end = float(sys.argv[4]) if len(sys.argv) > 4 else 1.5
    dt = 0.005
    config = make_config("E2", form, variant, t_end, cells)
    config = replace(
        config,
        run_name=f"E2_{form}_{variant}_{cells}",
        mesh=MeshOptions(cells=(cells, cells)),
        time=TimeOptions(scheme="crank_nicolson", dt=dt, t_end=t_end, startup="lower_order"),
        output=OutputOptions(vtk_every=0, diagnostics=DIAGNOSTICS),
    )
    validate_config(config)
    problem = replace(shear_layer(cells), viscosity=VISCOSITY)
    check_diagnostics_supported(config.output.diagnostics, problem)
    folder = HERE / "results" / config.run_name
    folder.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=folder / "driver.log",
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
    )
    logger = logging.getLogger(config.run_name)
    stepper = TimeStepper(config, problem, folder, logger)
    mesh = problem.mesh
    dx = fd.dx(domain=mesh, degree=4)
    u = stepper.velocity
    vorticity = u[1].dx(0) - u[0].dx(1)
    dg0 = fd.FunctionSpace(mesh, "DG", 0)
    centres = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)
    )
    rows, snapshots = [], {}

    def gather(local: np.ndarray) -> np.ndarray:
        """Owned DG0 values of every rank, concatenated in rank order (rank 0 keeps them)."""
        owned = local[: dg0.dof_dset.size] if local.ndim == 1 else local[: dg0.dof_dset.size]
        parts = mesh.comm.gather(np.array(owned), root=0)
        return np.concatenate(parts) if mesh.comm.rank == 0 else np.empty(0)

    pending = list(SNAPSHOTS)

    def record(step: int, time: float) -> bool:
        rows.append(
            {
                "step": step,
                "time": time,
                "energy": fd.assemble(ufl.inner(u, u) / 2 * dx),
                "enstrophy": fd.assemble(vorticity**2 / 2 * dx),
                "momentum_x": fd.assemble(u[0] * dx),
                "momentum_y": fd.assemble(u[1] * dx),
            }
        )
        if pending and time >= pending[0] - 1e-9:
            omega = fd.Function(dg0).interpolate(vorticity)
            snapshots[f"t{pending.pop(0):.2f}"] = gather(omega.dat.data_ro)
        return False

    status = 0
    try:
        stepper.run(record)
    except Exception as error:  # noqa: BLE001 — a failed run is a result, recorded below
        logger.error("Run stopped: %s", error)
        status = 1
    with (folder / "extra.csv").open("w", newline="") as stream:
        if rows and mesh.comm.rank == 0:
            writer = csv.DictWriter(stream, list(rows[-1]))
            writer.writeheader()
            writer.writerows(rows)
    all_centres = gather(centres.dat.data_ro)
    if mesh.comm.rank == 0:
        np.savez(folder / "vorticity.npz", centres=all_centres, **snapshots)
        (folder / "status.txt").write_text(f"{status}\n")
    return status


if __name__ == "__main__":
    sys.exit(main())
