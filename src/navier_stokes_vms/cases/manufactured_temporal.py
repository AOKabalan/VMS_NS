"""Rotation g(t)(y, −x) lies in every FE space, so only the temporal error remains (plan V2)."""

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

DEFAULT_CELLS = (4, 4)
DEFAULT_PARAMETERS = {"viscosity": 0.01}


def velocity(x: Expr, t: Expr) -> Expr:
    """u = g(t)(y, −x) with g(t) = cos t: nonzero vorticity, so the pressure cannot absorb it."""
    return ufl.cos(t) * ufl.as_vector((x[1], -x[0]))


def pressure(_x: Expr, _t: Expr) -> Expr:
    """p = 0."""
    return ufl.zero()


def build(config: Config) -> Problem:
    """Unit square with the exact velocity imposed on all four sides."""
    viscosity = case_parameters(config, DEFAULT_PARAMETERS)["viscosity"]
    cells = mesh_loading.builtin_cells(config, DEFAULT_CELLS)
    mesh, tags = mesh_loading.rectangle(cells, (0.0, 0.0), (1.0, 1.0))
    return Problem(
        mesh=mesh,
        boundary_conditions={name: Dirichlet(velocity) for name in tags},
        boundary_tags=tags,
        viscosity=viscosity,
        forcing=exact_solution_forcing(velocity, pressure, viscosity),
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
        exact_pressure=pressure,
    )
