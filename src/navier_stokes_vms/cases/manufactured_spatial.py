"""Manufactured solution linear in time, so only the spatial error remains (plan test V1)."""

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

DEFAULT_CELLS = (16, 16)
DEFAULT_PARAMETERS = {"viscosity": 0.01}


def velocity(x: Expr, t: Expr) -> Expr:
    """u = (1 + t)(∂ψ/∂y, −∂ψ/∂x) with ψ = sin²(πx) sin²(πy): solenoidal, zero on the boundary."""
    stream_function = ufl.sin(ufl.pi * x[0]) ** 2 * ufl.sin(ufl.pi * x[1]) ** 2
    return (1 + t) * ufl.as_vector((stream_function.dx(1), -stream_function.dx(0)))


def pressure(x: Expr, t: Expr) -> Expr:
    """p = (1 + t) cos(πx) sin(πy), with zero mean on the unit square."""
    return (1 + t) * ufl.cos(ufl.pi * x[0]) * ufl.sin(ufl.pi * x[1])


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
