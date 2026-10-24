"""E6/E7 robustness tests: python driver_rotation.py <form> <C> <cells> [steps] [rotation|shear].

u = C (−y, x) on (−½, ½)² with its own Dirichlet data, ν = 1e-6, f = 0, is an exact steady
Navier–Stokes solution with p = C² r²/2. The velocity is P1, but p is not, so a scheme whose P1
pressure must balance the centripetal gradient (u·∇)u = −∇p has an error that grows with C. For
EMAC the convective term 2(∇ˢu)u + (∇·u)u vanishes (∇ˢu = 0) and P = p − ½|u|² = 0, so EMAC–VMS
should reproduce the rotation exactly. The analogue of a pressure-robustness (no-flow) test.

E7 (shear) is the dual: Couette flow u = C (y, 0) with p constant has (u·∇)u = 0, so the classical
forms need no pressure at all, while EMAC must represent P = −½C²y² in P1.
"""

import logging
import sys
from dataclasses import replace

import ufl
from driver import HERE, VISCOSITY, make_config

from navier_stokes_vms.config import MeshOptions, OutputOptions, TimeOptions, validate_config
from navier_stokes_vms.postprocessing.diagnostics import check_diagnostics_supported
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.boundary_conditions import Dirichlet
from navier_stokes_vms.problem.problem_definition import Problem, at_initial_time
from navier_stokes_vms.solver.time_stepper import TimeStepper

DIAGNOSTICS = ("kinetic_energy", "velocity_error_l2", "divergence_l2")


def rotation_problem(speed: float, cells: int, flow: str = "rotation") -> Problem:
    mesh, tags = mesh_loading.rectangle((cells, cells), (-0.5, -0.5), (0.5, 0.5))

    def velocity(x: ufl.core.expr.Expr, _t: object) -> ufl.core.expr.Expr:
        if flow == "shear":
            return speed * ufl.as_vector((x[1], 0.0))
        return speed * ufl.as_vector((-x[1], x[0]))

    return Problem(
        mesh=mesh,
        boundary_conditions={name: Dirichlet(velocity) for name in tags},
        boundary_tags=tags,
        viscosity=VISCOSITY,
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
    )


def main() -> int:
    form, speed, cells = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    steps = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    flow = sys.argv[5] if len(sys.argv) > 5 else "rotation"
    experiment = "E7" if flow == "shear" else "E6"
    dt = 0.01
    config = replace(
        make_config("E6", form, "oss", steps * dt, cells),
        run_name=f"{experiment}_{form}_C{speed:g}_{cells}",
        mesh=MeshOptions(cells=(cells, cells)),
        time=TimeOptions(scheme="crank_nicolson", dt=dt, t_end=steps * dt, startup="lower_order"),
        output=OutputOptions(vtk_every=0, diagnostics=DIAGNOSTICS),
    )
    validate_config(config)
    problem = rotation_problem(speed, cells, flow)
    check_diagnostics_supported(config.output.diagnostics, problem)
    folder = HERE / "results" / config.run_name
    folder.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(filename=folder / "driver.log", level=logging.INFO)
    stepper = TimeStepper(config, problem, folder, logging.getLogger(config.run_name))
    status = 0
    try:
        stepper.run()
    except Exception as error:  # noqa: BLE001 — a failed run is a result
        logging.getLogger(config.run_name).error("Run stopped: %s", error)
        status = 1
    (folder / "status.txt").write_text(f"{status}\n")
    return status


if __name__ == "__main__":
    sys.exit(main())
