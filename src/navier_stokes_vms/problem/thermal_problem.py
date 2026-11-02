"""Thermal problem data and validation for the Boussinesq extension."""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

from ufl.core.expr import Expr

from navier_stokes_vms.config import Config, validate_config
from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, Dirichlet, Outflow
from navier_stokes_vms.problem.scalar_problem import transported_scalars

if TYPE_CHECKING:
    from navier_stokes_vms.problem.problem_definition import Problem


@dataclass(frozen=True)
class ThermalProblem:
    """Constant thermal properties, independent boundary data and optional exact temperature."""

    diffusivity: float
    expansion_coefficient: float
    gravity: tuple[float, ...]
    reference_temperature: float
    source: Callable[[Expr, Expr], Expr]
    initial_temperature: Callable[[Expr], Expr]
    boundary_conditions: Mapping[str, Dirichlet | DiffusiveFlux]
    exact_temperature: Callable[[Expr, Expr], Expr] | None = None
    scaling: Literal["dimensional", "nondimensional"] = "dimensional"


def check_thermal_problem(problem: Problem) -> None:
    """Check material values and complete, independent thermal boundary coverage."""
    thermal = problem.thermal
    if thermal is None:
        if problem.passive_scalars:
            raise ValueError("passive scalars require thermal problem data")
        return
    if thermal.scaling not in ("dimensional", "nondimensional"):
        raise ValueError("thermal scaling must be dimensional or nondimensional")
    for name, value, positive in (
        ("viscosity", problem.viscosity, True),
        ("diffusivity", thermal.diffusivity, True),
        ("expansion_coefficient", thermal.expansion_coefficient, False),
    ):
        if (
            not isinstance(value, int | float)
            or isinstance(value, bool)
            or not math.isfinite(value)
            or (value <= 0 if positive else value < 0)
        ):
            raise ValueError(
                f"thermal {name} must be finite and {'positive' if positive else 'nonnegative'}"
            )
    if (
        not isinstance(thermal.reference_temperature, int | float)
        or isinstance(thermal.reference_temperature, bool)
        or not math.isfinite(thermal.reference_temperature)
    ):
        raise ValueError("thermal reference_temperature must be finite")
    if len(thermal.gravity) != problem.mesh.geometric_dimension or not all(
        isinstance(g, int | float) and not isinstance(g, bool) and math.isfinite(g)
        for g in thermal.gravity
    ):
        raise ValueError("thermal gravity must be finite and match the mesh dimension")
    if set(thermal.boundary_conditions) != set(problem.boundary_tags):
        raise ValueError("thermal boundaries must cover exactly the named mesh boundaries")
    if any(
        not isinstance(bc, Dirichlet | DiffusiveFlux) for bc in thermal.boundary_conditions.values()
    ):
        raise ValueError("thermal boundaries require Dirichlet or DiffusiveFlux conditions")

    names = {"velocity", "pressure", "temperature", "xi", "zeta"}
    for scalar in problem.passive_scalars:
        if (
            not isinstance(scalar.name, str)
            or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", scalar.name)
            or scalar.name in names
            or scalar.name.startswith("eta_")
        ):
            raise ValueError("passive scalar names must be unique, valid and nonreserved")
        names.add(scalar.name)
        if (
            isinstance(scalar.diffusivity, bool)
            or not isinstance(scalar.diffusivity, int | float)
            or not math.isfinite(scalar.diffusivity)
            or scalar.diffusivity <= 0
        ):
            raise ValueError("passive scalar diffusivity must be finite and positive")
        if set(scalar.boundary_conditions) != set(problem.boundary_tags) or any(
            not isinstance(bc, Dirichlet | DiffusiveFlux)
            for bc in scalar.boundary_conditions.values()
        ):
            raise ValueError("passive scalar boundaries must cover the named mesh boundaries")


def check_thermal_supported(config: Config, problem: Problem) -> None:
    """Reject unsupported nonlinear thermal variants and missing startup data."""
    check_thermal_problem(problem)
    if problem.thermal is None:
        if config.thermal is not None:
            raise ValueError("[thermal] requires a case with thermal problem data")
        return
    validate_config(config)
    if config.stabilization.method == "vms" and config.stabilization.scale_splitting != "linear":
        raise ValueError("B3 thermal VMS requires linear splitting")
    if config.time.startup == "exact_history":
        missing = [] if problem.exact_velocity is not None else ["velocity"]
        missing += [s.name for s in transported_scalars(problem) if s.exact_solution is None]
        if missing:
            raise ValueError(
                "exact_history requires exact data for every transported field; "
                f"missing: {', '.join(missing)}"
            )
    if {"energy_budget", "momentum_budget"}.intersection(config.output.diagnostics):
        raise ValueError("thermal energy/momentum budgets require the B5 coupling identities")
    if config.solver.linearization != "newton" or config.solver.picard_iterations_before_newton:
        raise ValueError("B3 thermal runs require Newton without Picard warm-up")
    if config.viscosity.form != "laplacian" or config.convection.form not in (
        "skew_symmetric_1",
        "nonconservative",
    ):
        raise ValueError(
            "B3 thermal runs require Laplacian viscosity and "
            "skew_symmetric_1 or nonconservative momentum convection"
        )
    if any(isinstance(bc, Outflow) for bc in problem.boundary_conditions.values()):
        raise ValueError("B3 thermal runs require periodic or closed velocity boundaries")
    if (
        config.stabilization.method == "none"
        and config.discretization.velocity_degree != config.discretization.pressure_degree + 1
    ):
        raise ValueError("B3 thermal Galerkin runs require Taylor–Hood velocity/pressure degrees")
    degree = config.thermal.temperature_degree if config.thermal is not None else None
    if degree is not None and degree != config.discretization.velocity_degree:
        raise ValueError("B1 temperature degree must equal the velocity degree")
