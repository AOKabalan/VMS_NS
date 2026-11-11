"""Single shear layer of Olshanskii & Rebholz (2020), section 4.3, Figure 4.

Unit square, periodic x, impermeable free-slip horizontal walls. The Reynolds
number is based on delta_0, not the box length: nu = delta_0 * u_infinity / Re.
There is no known exact evolving solution supplied by this case.
"""

import math

import firedrake as fd
import ufl

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem.boundary_conditions import FreeSlip
from navier_stokes_vms.problem.mesh_loading import builtin_cells
from navier_stokes_vms.problem.problem_definition import Problem, case_parameters


def build(config: Config) -> Problem:
    """Build the paper's perturbed tanh layer with Re=1000 and delta_0=1/28."""
    values = case_parameters(
        config, {"reynolds": 1000.0, "delta_0": 1 / 28, "u_infinity": 1.0, "noise": 1e-3}
    )
    nx, ny = builtin_cells(config, (196, 196))
    mesh = fd.PeriodicRectangleMesh(nx, ny, 1.0, 1.0, direction="x")
    delta, speed, noise = (values[k] for k in ("delta_0", "u_infinity", "noise"))
    x = fd.SpatialCoordinate(mesh)
    psi = (
        speed
        * ufl.exp(-(((x[1] - 0.5) / delta) ** 2))
        * (ufl.cos(8 * math.pi * x[0]) + ufl.cos(20 * math.pi * x[0]))
    )
    initial = ufl.as_vector(
        (speed * ufl.tanh((2 * x[1] - 1) / delta) + noise * psi.dx(1), -noise * psi.dx(0))
    )
    return Problem(
        mesh=mesh,
        boundary_conditions={"y_min": FreeSlip(1), "y_max": FreeSlip(1)},
        boundary_tags={"y_min": 3, "y_max": 4},
        viscosity=delta * speed / values["reynolds"],
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=lambda _x: initial,
        reference_values=values,
    )
