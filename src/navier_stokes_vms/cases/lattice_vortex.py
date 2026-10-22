"""Lattice vortex: Kabalan MSc thesis, §6.2, equations (6.3)–(6.5), pp.44–45."""

import math

import numpy as np
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


def build(config: Config) -> Problem:
    """Exact Dirichlet data on the unit square, ν=7e-5 and h=0.02 by default."""
    viscosity = case_parameters(config, {"viscosity": 7e-5})["viscosity"]

    def velocity(x: Expr, t: Expr) -> Expr:
        return ufl.exp(-8 * viscosity * math.pi**2 * t) * ufl.as_vector(
            (
                ufl.sin(2 * math.pi * x[0]) * ufl.sin(2 * math.pi * x[1]),
                ufl.cos(2 * math.pi * x[0]) * ufl.cos(2 * math.pi * x[1]),
            )
        )

    def pressure(x: Expr, t: Expr) -> Expr:
        return (
            -0.5
            * ufl.exp(-16 * viscosity * math.pi**2 * t)
            * (ufl.sin(2 * math.pi * x[0]) ** 2 + ufl.cos(2 * math.pi * x[1]) ** 2)
        )

    mesh, tags = mesh_loading.rectangle(
        mesh_loading.builtin_cells(config, (50, 50)), (0.0, 0.0), (1.0, 1.0)
    )
    return Problem(
        mesh=mesh,
        boundary_conditions={name: Dirichlet(velocity) for name in tags},
        boundary_tags=tags,
        viscosity=viscosity,
        forcing=exact_solution_forcing(velocity, pressure, viscosity),
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
        exact_pressure=pressure,
        reference_values={"initial_energy": 0.25, "energy_decay_rate": 16 * viscosity * math.pi**2},
    )


def compare(times: np.ndarray, energy: np.ndarray, viscosity: float = 7e-5) -> dict[str, float]:
    """Compare global kinetic energy with ¼ exp(−16νπ²t), including the endpoint error."""
    exact = 0.25 * np.exp(-16 * viscosity * math.pi**2 * times)
    relative = np.abs(energy - exact) / exact
    return {
        "energy_relative_error_at_T": float(relative[-1]),
        "energy_max_relative_error": float(np.max(relative)),
    }
