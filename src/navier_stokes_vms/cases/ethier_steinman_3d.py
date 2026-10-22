"""Ethier–Steinman exact solution (1994, doi:10.1002/fld.1650190502), with zero forcing."""

import math

import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.boundary_conditions import Dirichlet
from navier_stokes_vms.problem.problem_definition import (
    Problem,
    at_initial_time,
    case_parameters,
    exact_solution_forcing,
)

DEFAULT_CELLS = (4, 4, 4)
DEFAULT_PARAMETERS = {"viscosity": 1.0}
LOWER = (-1.0, -1.0, -1.0)
UPPER = (1.0, 1.0, 1.0)


def build(config: Config) -> Problem:
    """Cube [−1, 1]³ with a = π/4, d = π and exact Dirichlet data on every face."""
    viscosity = case_parameters(config, DEFAULT_PARAMETERS)["viscosity"]
    # d = π makes the nodal Dirichlet trace flux-compatible also for P1 on this mesh.
    a, d = math.pi / 4, math.pi

    def velocity(x: Expr, t: Expr) -> Expr:
        """The cyclic exponential–trigonometric velocity, decaying as exp(−νd²t)."""
        return (
            -a
            * ufl.exp(-viscosity * d**2 * t)
            * ufl.as_vector(
                tuple(
                    ufl.exp(a * x[i]) * ufl.sin(a * x[j] + d * x[k])
                    + ufl.exp(a * x[k]) * ufl.cos(a * x[i] + d * x[j])
                    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1))
                )
            )
        )

    def pressure(x: Expr, t: Expr) -> Expr:
        """p = −|u|²/2: this Beltrami velocity has u × curl u = 0."""
        u = velocity(x, t)
        return -ufl.inner(u, u) / 2

    cells = mesh_loading.builtin_cells(config, DEFAULT_CELLS)
    mesh, tags = mesh_loading.box(cells, LOWER, UPPER)
    return Problem(
        mesh=mesh,
        boundary_conditions={name: Dirichlet(velocity) for name in tags},
        boundary_tags=tags,
        viscosity=viscosity,
        forcing=exact_solution_forcing(velocity, pressure, viscosity),
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
        exact_pressure=pressure,
        force_boundaries=tuple(tags),
    )
