"""Decaying Taylor–Green vortex, an exact Navier–Stokes solution with f = 0 (plan test V4)."""

import math

import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.problem_definition import (
    Problem,
    at_initial_time,
    case_parameters,
    exact_solution_forcing,
)

DEFAULT_CELLS = (32, 32)
DEFAULT_PARAMETERS = {"viscosity": 0.01}
LOWER = (0.0, 0.0)
UPPER = (2 * math.pi, 2 * math.pi)


def build(config: Config) -> Problem:
    """Periodic square (0, 2π)², so there are no boundaries."""
    viscosity = case_parameters(config, DEFAULT_PARAMETERS)["viscosity"]

    def velocity(x: Expr, t: Expr) -> Expr:
        """u = (−cos x sin y, sin x cos y) e^(−2νt)."""
        return ufl.exp(-2 * viscosity * t) * ufl.as_vector(
            (-ufl.cos(x[0]) * ufl.sin(x[1]), ufl.sin(x[0]) * ufl.cos(x[1]))
        )

    def pressure(x: Expr, t: Expr) -> Expr:
        """p = −¼(cos 2x + cos 2y) e^(−4νt)."""
        return -ufl.exp(-4 * viscosity * t) * (ufl.cos(2 * x[0]) + ufl.cos(2 * x[1])) / 4

    cells = mesh_loading.builtin_cells(config, DEFAULT_CELLS)
    mesh, tags = mesh_loading.periodic_rectangle(cells, LOWER, UPPER)
    return Problem(
        mesh=mesh,
        boundary_conditions={},
        boundary_tags=tags,
        viscosity=viscosity,
        forcing=exact_solution_forcing(velocity, pressure, viscosity),
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
        exact_pressure=pressure,
    )
