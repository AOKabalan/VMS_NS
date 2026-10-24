"""EMAC–VMS advantage experiments: python driver.py <experiment> <form> <variant> [T].

E1 translating Gresho vortex in a periodic unit box (linear momentum and energy invariants);
E3 off-centre Gresho vortex in a large no-slip box (angular momentum far from the walls).
Each run writes the solver's own outputs plus extra.csv, computed after every step.
"""

import csv
import logging
import math
import sys
from dataclasses import replace
from pathlib import Path
from types import MappingProxyType

import firedrake as fd
import ufl

from navier_stokes_vms.cases import gresho
from navier_stokes_vms.config import (
    Config,
    ConvectionOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
    ViscosityOptions,
    validate_config,
)
from navier_stokes_vms.postprocessing.diagnostics import check_diagnostics_supported
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.boundary_conditions import Dirichlet
from navier_stokes_vms.problem.problem_definition import Problem, at_initial_time
from navier_stokes_vms.solver.time_stepper import TimeStepper

HERE = Path(__file__).parent
VISCOSITY = 1e-6
BACKGROUND = (1.0, 0.5)
BOX_CENTRE = (0.1, 0.05)
VARIANTS = {
    # name: (subscale space, pressure subscale, c2)
    "oss": ("oss", True, 2.0),
    "oss_weak": ("oss", False, 32.0),
    "asgs": ("asgs", True, 2.0),
}
SOLVER = SolverOptions(
    linearization="newton",
    newton_tau="frozen",
    inner_rtol=1e-10,
    inner_max_iterations=200,
    outer_rtol=1e-9,
    outer_max_iterations=100,
    subscale_rtol=1e-11,
    subscale_max_iterations=500,
)
COMMON_DIAGNOSTICS = (
    "kinetic_energy",
    "subscale_kinetic_energy",
    "linear_momentum",
    "energy_budget",
    "momentum_budget",
    "linear_momentum_defect",
    "convective_energy_defect",
    "divergence_l2",
)


def smooth_velocity(x: ufl.core.expr.Expr, _t: object) -> ufl.core.expr.Expr:
    """Compact C³ swirl u_θ = 12 r (1 − r²/0.16)⁴ for r < 0.4: a steady Euler field, peak ≈ 1."""
    r2 = ufl.dot(x, x)
    bump = ufl.max_value(1 - r2 / 0.16, 0.0)
    # u = (u_θ/r)(−y, x) with u_θ/r = 12 bump⁴: no division by r.
    return 12 * bump**4 * ufl.as_vector((-x[1], x[0]))


def swirl(smooth: bool):  # noqa: ANN201
    return smooth_velocity if smooth else gresho.velocity


def wrap(value: float) -> float:
    """Periodic image of a coordinate in [−0.5, 0.5)."""
    return value - math.floor(value + 0.5)


def translating_problem(cells: int, smooth: bool = False) -> tuple[Problem, fd.Constant]:
    """Unit periodic box; u = U₀ + Σ images of the Gresho vortex centred at c(t) = U₀t."""
    mesh, tags = mesh_loading.periodic_rectangle((cells, cells), (-0.5, -0.5), (0.5, 0.5))
    centre = fd.Constant((0.0, 0.0))
    background = ufl.as_vector(BACKGROUND)

    def vortex(x: ufl.core.expr.Expr, c: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        # Support radius 0.4 < 0.5, so the images never overlap inside the unit cell.
        return sum(
            swirl(smooth)(x - c - ufl.as_vector((i, j)), 0.0)
            for i in (-1, 0, 1)
            for j in (-1, 0, 1)
        )

    problem = Problem(
        mesh=mesh,
        boundary_conditions={},
        boundary_tags=tags,
        viscosity=VISCOSITY,
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=at_initial_time(
            lambda x, _t: background + vortex(x, ufl.as_vector((0.0, 0.0)))
        ),
    )
    exact = lambda x: background + vortex(x, centre)  # noqa: E731
    return problem, (centre, exact)


def box_problem(cells: int, smooth: bool = False) -> Problem:
    """(−1, 1)², no-slip walls, steady Gresho vortex centred off the origin."""
    mesh, tags = mesh_loading.rectangle((cells, cells), (-1.0, -1.0), (1.0, 1.0))
    shift = ufl.as_vector(BOX_CENTRE)
    return Problem(
        mesh=mesh,
        boundary_conditions={name: Dirichlet(lambda _x, _t: ufl.zero(2)) for name in tags},
        boundary_tags=tags,
        viscosity=VISCOSITY,
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=at_initial_time(lambda x, t: swirl(smooth)(x - shift, t)),
        exact_velocity=lambda x, t: swirl(smooth)(x - shift, t),
        force_boundaries=tuple(tags),
        monitor_region=lambda x: ufl.conditional(
            ufl.ge(ufl.max_value(abs(x[0]), abs(x[1])), 0.9), 1.0, 0.0
        ),
    )


def make_config(experiment: str, form: str, variant: str, t_end: float, cells: int) -> Config:
    space, tau_c, c2 = VARIANTS[variant]
    diagnostics = COMMON_DIAGNOSTICS
    if experiment in ("E3", "E5"):
        diagnostics += (
            "angular_momentum",
            "angular_momentum_defect",
            "force",
            "region_velocity_l2",
            "velocity_error_l2",
        )
    return Config(
        case="gresho",
        run_name=f"{experiment}_{form}_{variant}",
        case_parameters=MappingProxyType({}),
        mesh=MeshOptions(cells=(cells, cells)),
        time=TimeOptions(scheme="crank_nicolson", dt=0.01, t_end=t_end, startup="lower_order"),
        stabilization=StabilizationOptions(
            method="vms",
            subscale_space=space,
            subscale_dynamics="dynamic",
            scale_splitting="nonlinear",
            c2=c2,
            pressure_subscale=tau_c,
        ),
        convection=ConvectionOptions(form=form),
        viscosity=ViscosityOptions(form="stress"),
        solver=SOLVER,
        output=OutputOptions(vtk_every=0, diagnostics=diagnostics),
    )


def main() -> int:
    experiment, form, variant = sys.argv[1:4]
    periodic = experiment in ("E1", "E4")
    t_end = float(sys.argv[4]) if len(sys.argv) > 4 else (4.0 if periodic else 10.0)
    cells = int(sys.argv[5]) if len(sys.argv) > 5 else (32 if periodic else 48)
    config = make_config(experiment, form, variant, t_end, cells)
    if len(sys.argv) > 5:
        config = replace(config, run_name=f"{config.run_name}_{cells}")
    validate_config(config)
    exact = None
    smooth = experiment in ("E4", "E5")
    if periodic:
        problem, exact = translating_problem(cells, smooth)
    else:
        problem = box_problem(cells, smooth)
    check_diagnostics_supported(config.output.diagnostics, problem)
    folder = HERE / "results" / config.run_name
    folder.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(
        filename=folder / "driver.log",
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
    )
    logger = logging.getLogger(config.run_name)
    stepper = TimeStepper(config, problem, folder, logger)

    mesh = problem.mesh
    x = fd.SpatialCoordinate(mesh)
    dx = fd.dx(domain=mesh, degree=8)
    u = stepper.velocity
    vorticity = u[1].dx(0) - u[0].dx(1)
    weight = vorticity**2
    phase = 2 * math.pi
    rows = []

    def record(step: int, time: float) -> bool:
        row = {
            "step": step,
            "time": time,
            "energy": fd.assemble(ufl.inner(u, u) / 2 * dx),
            "momentum_x": fd.assemble(u[0] * dx),
            "momentum_y": fd.assemble(u[1] * dx),
            "enstrophy": fd.assemble(weight / 2 * dx),
        }
        if exact is not None:
            centre, field = exact
            c = (wrap(BACKGROUND[0] * time), wrap(BACKGROUND[1] * time))
            centre.assign(c)
            difference = u - field(x)
            row["error_l2"] = math.sqrt(fd.assemble(ufl.inner(difference, difference) * dx))
            # Periodic centroid of the enstrophy density, from its first Fourier moments.
            moments = [
                complex(
                    fd.assemble(weight * ufl.cos(phase * x[k]) * dx),
                    fd.assemble(weight * ufl.sin(phase * x[k]) * dx),
                )
                for k in (0, 1)
            ]
            found = [math.atan2(m.imag, m.real) / phase for m in moments]
            row["centre_x"], row["centre_y"] = found
            row["centre_error"] = math.hypot(wrap(found[0] - c[0]), wrap(found[1] - c[1]))
        rows.append(row)
        return False

    status = 0
    try:
        stepper.run(record)
    except Exception as error:  # noqa: BLE001 — a failed run is a result, recorded below
        logger.error("Run stopped: %s", error)
        status = 1
    with (folder / "extra.csv").open("w", newline="") as stream:
        if rows:
            writer = csv.DictWriter(stream, list(rows[-1]))
            writer.writeheader()
            writer.writerows(rows)
    (folder / "status.txt").write_text(f"{status}\n")
    return status


if __name__ == "__main__":
    sys.exit(main())
