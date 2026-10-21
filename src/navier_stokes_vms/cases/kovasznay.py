"""Kovasznay flow, a steady exact Navier–Stokes solution with ν = 1/Re (plan test V3)."""

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

DEFAULT_CELLS = (24, 32)
DEFAULT_PARAMETERS = {"reynolds": 40.0}
LOWER = (-0.5, -0.5)
UPPER = (1.0, 1.5)


def build(config: Config) -> Problem:
    """[−0.5, 1] × [−0.5, 1.5] with the exact velocity imposed on all four sides."""
    reynolds = case_parameters(config, DEFAULT_PARAMETERS)["reynolds"]
    decay = reynolds / 2 - math.sqrt(reynolds**2 / 4 + 4 * math.pi**2)

    def velocity(x: Expr, _t: Expr) -> Expr:
        """u = (1 − e^(λx) cos 2πy, λ/(2π) e^(λx) sin 2πy)."""
        growth = ufl.exp(decay * x[0])
        return ufl.as_vector(
            (
                1 - growth * ufl.cos(2 * ufl.pi * x[1]),
                decay / (2 * ufl.pi) * growth * ufl.sin(2 * ufl.pi * x[1]),
            )
        )

    def pressure(x: Expr, _t: Expr) -> Expr:
        """p = ½(1 − e^(2λx))."""
        return (1 - ufl.exp(2 * decay * x[0])) / 2

    viscosity = 1 / reynolds
    cells = mesh_loading.builtin_cells(config, DEFAULT_CELLS)
    mesh, tags = mesh_loading.rectangle(cells, LOWER, UPPER)
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
