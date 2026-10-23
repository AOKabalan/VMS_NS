"""Problem data of a case, forcing from exact solutions and the null-space decision (plan §8)."""

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass

import firedrake as fd
import ufl
from firedrake.mesh import MeshGeometry
from ufl.core.expr import Expr
from ufl.domain import extract_unique_domain

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem.boundary_conditions import Dirichlet, Outflow, boundary_tag

SpaceTimeExpression = Callable[[Expr, Expr], Expr]

_DIVERGENCE_TOLERANCE = 1e-10
_CHECK_QUADRATURE_DEGREE = 6
_CHECK_TIMES = (0.0, 1.0)


@dataclass(frozen=True)
class Problem:
    """Mesh, boundaries, viscosity, data and optional exact solution of one case."""

    mesh: MeshGeometry
    boundary_conditions: Mapping[str, Dirichlet | Outflow]
    boundary_tags: Mapping[str, int]
    viscosity: float
    forcing: SpaceTimeExpression
    initial_velocity: Callable[[Expr], Expr]
    exact_velocity: SpaceTimeExpression | None = None
    exact_pressure: SpaceTimeExpression | None = None
    reference_values: Mapping[str, float] | None = None
    force_boundaries: tuple[str, ...] = ()
    monitor_region: Callable[[Expr], Expr] | None = None


def exact_solution_forcing(
    velocity: SpaceTimeExpression, pressure: SpaceTimeExpression, viscosity: float
) -> SpaceTimeExpression:
    """f = ∂t u + u·∇u − νΔu + ∇p of an exact solution, ∂t taken through ufl.variable."""

    def forcing(x: Expr, t: Expr) -> Expr:
        time = ufl.variable(t)
        u = velocity(x, time)
        p = pressure(x, time)
        if extract_unique_domain(p) is None:
            pressure_gradient = ufl.zero(*u.ufl_shape)
        else:
            pressure_gradient = ufl.grad(p)
        return (
            ufl.diff(u, time)
            + ufl.dot(ufl.grad(u), u)
            - viscosity * ufl.div(ufl.grad(u))
            + pressure_gradient
        )

    return forcing


def at_initial_time(velocity: SpaceTimeExpression) -> Callable[[Expr], Expr]:
    """The initial velocity u₀(x) = u(x, 0) of an exact solution."""
    return lambda x: velocity(x, ufl.as_ufl(0.0))


def case_parameters(config: Config, defaults: Mapping[str, float]) -> dict[str, float]:
    """[case_parameters] over the case's defaults; unknown names or non-positive values fail."""
    unknown = sorted(set(config.case_parameters) - set(defaults))
    if unknown:
        raise ValueError(
            f"case_parameters: {', '.join(unknown)} not used by case {config.case!r}; "
            f"it accepts {', '.join(sorted(defaults))}"
        )
    values = {**defaults, **config.case_parameters}
    for name, value in values.items():
        number = isinstance(value, int | float) and not isinstance(value, bool)
        if not (number and math.isfinite(value) and value > 0):
            raise ValueError(f"case_parameters.{name} = {value!r}: must be a positive number")
    return {name: float(value) for name, value in values.items()}


def check_problem(problem: Problem) -> None:
    """Every named boundary has exactly one known condition; an exact velocity is solenoidal."""
    for name in problem.boundary_conditions:
        boundary_tag(name, problem.boundary_tags)
    for name in problem.force_boundaries:
        boundary_tag(name, problem.boundary_tags)
    missing = sorted(set(problem.boundary_tags) - set(problem.boundary_conditions))
    if missing:
        raise ValueError(f"boundaries without a condition: {', '.join(missing)}")
    if problem.exact_velocity is not None:
        _check_divergence_free(problem.mesh, problem.exact_velocity)


def uses_pressure_null_space(problem: Problem) -> bool:
    """The constant-pressure null space applies exactly when no boundary is an outflow."""
    return not any(isinstance(bc, Outflow) for bc in problem.boundary_conditions.values())


def _check_divergence_free(mesh: MeshGeometry, velocity: SpaceTimeExpression) -> None:
    x = fd.SpatialCoordinate(mesh)
    dx = fd.dx(domain=mesh, degree=_CHECK_QUADRATURE_DEGREE)
    for time in _CHECK_TIMES:
        u = velocity(x, fd.Constant(time))
        divergence = math.sqrt(fd.assemble(ufl.div(u) ** 2 * dx))
        gradient = math.sqrt(fd.assemble(ufl.inner(ufl.grad(u), ufl.grad(u)) * dx))
        if divergence > _DIVERGENCE_TOLERANCE * gradient:
            raise ValueError(
                f"exact velocity is not divergence-free at t = {time}: "
                f"‖div u‖ = {divergence:.3e}, ‖grad u‖ = {gradient:.3e}"
            )
