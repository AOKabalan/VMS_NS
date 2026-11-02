"""Closed natural-convection cavity in width-based free-fall units (B5)."""

import math

import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, Dirichlet
from navier_stokes_vms.problem.mesh_loading import builtin_cells, rectangle
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.problem.thermal_problem import ThermalProblem


def build(config: Config) -> Problem:
    """Build the square or tall cavity with explicit Ra, Pr, wall data and initial profile."""
    values = {
        "rayleigh": 1000.0,
        "prandtl": 0.71,
        "aspect_ratio": 1.0,
        "initial_profile": "conduction",
        **config.case_parameters,
    }
    if set(values) - {"rayleigh", "prandtl", "aspect_ratio", "initial_profile"}:
        raise ValueError("unknown differentially_heated_cavity case_parameters")
    for key in ("rayleigh", "prandtl", "aspect_ratio"):
        value = values[key]
        if (
            isinstance(value, bool)
            or not isinstance(value, int | float)
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError(f"case_parameters.{key} must be finite and positive")
    if values["initial_profile"] not in ("conduction", "isothermal"):
        raise ValueError("initial_profile must be conduction or isothermal")
    ra, pr, height = (float(values[key]) for key in ("rayleigh", "prandtl", "aspect_ratio"))
    mesh, tags = rectangle(builtin_cells(config, (16, 16)), (0.0, 0.0), (1.0, height))

    def zero(x: Expr, _t: Expr) -> Expr:
        return 0.0 * x[0]

    def velocity(x: Expr, _t: Expr) -> Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    def hot(x: Expr, _t: Expr) -> Expr:
        return 0.5 + 0.0 * x[0]

    def cold(x: Expr, _t: Expr) -> Expr:
        return -0.5 + 0.0 * x[0]

    def initial(x: Expr) -> Expr:
        return 0.5 - x[0] if values["initial_profile"] == "conduction" else 0.0 * x[0]

    return Problem(
        mesh=mesh,
        boundary_tags=tags,
        boundary_conditions={name: Dirichlet(velocity) for name in tags},
        viscosity=math.sqrt(pr / ra),
        forcing=velocity,
        initial_velocity=lambda x: velocity(x, ufl.as_ufl(0.0)),
        reference_values={
            "rayleigh": ra,
            "prandtl": pr,
            "height": height,
            "length_scale": 1.0,
            "temperature_difference": 1.0,
        },
        thermal=ThermalProblem(
            1 / math.sqrt(ra * pr),
            1.0,
            (0.0, -1.0),
            0.0,
            zero,
            initial,
            {
                "x_min": Dirichlet(hot),
                "x_max": Dirichlet(cold),
                "y_min": DiffusiveFlux(),
                "y_max": DiffusiveFlux(),
            },
            scaling="nondimensional",
        ),
    )
