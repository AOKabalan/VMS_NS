"""Named passive scalar data and a common transport view of the active scalars."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from ufl.core.expr import Expr

from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, Dirichlet

if TYPE_CHECKING:
    from navier_stokes_vms.problem.problem_definition import Problem


@dataclass(frozen=True)
class PassiveScalarProblem:
    """Independent scalar diffusion, source, initial/boundary data and optional exact solution."""

    name: str
    diffusivity: float
    source: Callable[[Expr, Expr], Expr]
    initial_value: Callable[[Expr], Expr]
    boundary_conditions: Mapping[str, Dirichlet | DiffusiveFlux]
    exact_solution: Callable[[Expr, Expr], Expr] | None = None


def transported_scalars(problem: Problem) -> tuple[PassiveScalarProblem, ...]:
    """Temperature followed by passive scalars as transport data without buoyancy coefficients."""
    thermal = problem.thermal
    if thermal is None:
        return ()
    return (
        PassiveScalarProblem(
            "temperature",
            thermal.diffusivity,
            thermal.source,
            thermal.initial_temperature,
            thermal.boundary_conditions,
            thermal.exact_temperature,
        ),
        *problem.passive_scalars,
    )
