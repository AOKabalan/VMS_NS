"""Uniform x-periodic free-slip channel with optional constant thermal/passive fields (G4)."""

import math

import firedrake as fd
import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, FreeSlip
from navier_stokes_vms.problem.mesh_loading import builtin_cells, periodic_channel
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.problem.scalar_problem import PassiveScalarProblem
from navier_stokes_vms.problem.thermal_problem import ThermalProblem


def build(config: Config) -> Problem:
    """Build the G4 preservation case; thermal stress support remains gated to B6."""
    defaults = {"length": 2.0, "height": 1.0, "speed": 1.0, "thermal": True}
    if set(config.case_parameters) - defaults.keys():
        raise ValueError("unknown free_slip_channel case_parameters")
    values = {**defaults, **config.case_parameters}
    for name in ("length", "height", "speed"):
        value = values[name]
        if (
            isinstance(value, bool)
            or not isinstance(value, int | float)
            or not math.isfinite(value)
        ):
            raise ValueError(f"case_parameters.{name} must be finite")
        if name != "speed" and value <= 0:
            raise ValueError(f"case_parameters.{name} must be positive")
    if not isinstance(values["thermal"], bool):
        raise ValueError("case_parameters.thermal must be boolean")
    length, height, speed = (float(values[key]) for key in ("length", "height", "speed"))
    mesh, tags = periodic_channel(
        builtin_cells(config, (8, 8)), (0.0, -height / 2), (length, height / 2)
    )
    initial = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1)).assign(fd.Constant((speed, 0.0)))

    def velocity(_x: Expr, _t: Expr) -> Expr:
        return initial

    def zero(x: Expr, _t: Expr) -> Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    def temperature(x: Expr, _t: Expr) -> Expr:
        return 1.0 + 0.0 * x[0]

    def moisture(x: Expr, _t: Expr) -> Expr:
        return 2.0 + 0.0 * x[0]

    def scalar_source(x: Expr, _t: Expr) -> Expr:
        return 0.0 * x[0]

    scalar_boundaries = {name: DiffusiveFlux() for name in tags}
    thermal = (
        ThermalProblem(
            0.1,
            0.7,
            (0.0, -1.0),
            1.0,
            scalar_source,
            lambda x: temperature(x, ufl.as_ufl(0.0)),
            scalar_boundaries,
            temperature,
            scaling="nondimensional",
        )
        if values["thermal"]
        else None
    )
    passive = (
        (
            PassiveScalarProblem(
                "q",
                0.1,
                scalar_source,
                lambda x: moisture(x, ufl.as_ufl(0.0)),
                scalar_boundaries,
                moisture,
            ),
        )
        if thermal is not None
        else ()
    )
    return Problem(
        mesh=mesh,
        boundary_tags=tags,
        boundary_conditions={name: FreeSlip(1) for name in tags},
        viscosity=0.1,
        forcing=zero,
        initial_velocity=lambda _x: initial,
        exact_velocity=velocity,
        exact_pressure=scalar_source,
        reference_values={"length_scale": length, "height": height, "speed": speed},
        thermal=thermal,
        passive_scalars=passive,
    )
