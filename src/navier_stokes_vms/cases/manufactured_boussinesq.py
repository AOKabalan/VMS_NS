"""Nondimensional coupled manufactured solution linear in time for B1 verification."""

import math

import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config
from navier_stokes_vms.formulation.boussinesq import buoyancy
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.boundary_conditions import Dirichlet
from navier_stokes_vms.problem.problem_definition import (
    Problem,
    at_initial_time,
    exact_solution_forcing,
)
from navier_stokes_vms.problem.thermal_problem import ThermalProblem


def velocity(x: Expr, t: Expr) -> Expr:
    """Solenoidal polynomial velocity vanishing at every boundary."""
    psi = x[0] ** 2 * (1 - x[0]) ** 2 * x[1] ** 2 * (1 - x[1]) ** 2
    return (1 + t) * ufl.as_vector((psi.dx(1), -psi.dx(0)))


def pressure(x: Expr, t: Expr) -> Expr:
    """Linear pressure with zero mean."""
    return (1 + t) * (x[0] - 0.5)


def temperature(x: Expr, t: Expr) -> Expr:
    """Bilinear temperature with time-dependent boundary data."""
    return (1 + t) * (1 + x[0] + x[1] + x[0] * x[1])


def build(config: Config) -> Problem:
    """Build the coupled manufactured case with explicit, finite physical parameters."""
    defaults = dict(
        viscosity=0.1, diffusivity=0.2, expansion_coefficient=1.0, reference_temperature=0.0
    )
    if set(config.case_parameters) - set(defaults):
        raise ValueError("unknown manufactured_boussinesq case_parameters")
    values = {**defaults, **config.case_parameters}
    if any(
        isinstance(v, bool) or not isinstance(v, int | float) or not math.isfinite(v)
        for v in values.values()
    ):
        raise ValueError("thermal case parameters must be finite numbers")
    nu, kappa, beta, reference = (float(values[key]) for key in defaults)
    mesh, tags = mesh_loading.rectangle(
        mesh_loading.builtin_cells(config, (8, 8)), (0.0, 0.0), (1.0, 1.0)
    )
    gravity = (0.0, -1.0)
    mechanical = exact_solution_forcing(velocity, pressure, nu)

    def forcing(x: Expr, t: Expr) -> Expr:
        return mechanical(x, t) - buoyancy(temperature(x, t), beta, gravity, reference)

    def source(x: Expr, t: Expr) -> Expr:
        time = ufl.variable(t)
        field = temperature(x, time)
        return (
            ufl.diff(field, time)
            + ufl.dot(velocity(x, time), ufl.grad(field))
            - kappa * ufl.div(ufl.grad(field))
        )

    return Problem(
        mesh=mesh,
        boundary_tags=tags,
        boundary_conditions={n: Dirichlet(velocity) for n in tags},
        viscosity=nu,
        forcing=forcing,
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
        exact_pressure=pressure,
        thermal=ThermalProblem(
            kappa,
            beta,
            gravity,
            reference,
            source,
            at_initial_time(temperature),
            {n: Dirichlet(temperature) for n in tags},
            temperature,
            scaling="nondimensional",
        ),
    )
