"""Registry of named per-step diagnostics; error norms use a rule of twice the common degree."""

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import firedrake as fd
import ufl
from ufl.core.expr import Expr
from ufl.measure import Measure

from navier_stokes_vms.problem.boundary_conditions import named_boundary_measure
from navier_stokes_vms.problem.problem_definition import Problem


@dataclass(frozen=True)
class DiagnosticState:
    """What a diagnostic reads after a step: endpoint velocity and solve-level pressure.

    The pressure is the physical p_h, P_h + ½|u_h|² for EMAC variants (spec: EMAC pressure).
    """

    problem: Problem
    velocity: fd.Function
    pressure: Expr
    velocity_time: fd.Constant
    pressure_time: fd.Constant
    dx: Measure
    error_dx: Measure
    pressure_null_space: bool
    solve_velocity: fd.Function
    ds: Measure
    subscale_velocity: Expr | None
    subscale_time: fd.Constant | None


@dataclass(frozen=True)
class Diagnostic:
    """A functional, possibly with named components, and its exact-solution requirement."""

    evaluate: Callable[[DiagnosticState], float | dict[str, float]]
    needs_exact_solution: bool


def _kinetic_energy(state: DiagnosticState) -> float:
    return fd.assemble(ufl.inner(state.velocity, state.velocity) / 2 * state.dx)


def _subscale_kinetic_energy(state: DiagnosticState) -> float:
    return fd.assemble(ufl.inner(state.subscale_velocity, state.subscale_velocity) / 2 * state.dx)


def _components(expression: Expr, measure: Measure) -> dict[str, float]:
    return {
        axis: fd.assemble(expression[i] * measure)
        for i, axis in enumerate("xyz"[: expression.ufl_shape[0]])
    }


def _linear_momentum(state: DiagnosticState) -> dict[str, float]:
    return _components(state.velocity, state.dx)


def _angular_momentum(state: DiagnosticState) -> float | dict[str, float]:
    x = fd.SpatialCoordinate(state.problem.mesh)
    if state.problem.mesh.geometric_dimension == 2:
        return fd.assemble((x[0] * state.velocity[1] - x[1] * state.velocity[0]) * state.dx)
    return _components(ufl.cross(x, state.velocity), state.dx)


def _force(state: DiagnosticState) -> dict[str, float]:
    mesh = state.problem.mesh
    gradient = ufl.grad(state.solve_velocity)
    stress = -state.pressure * ufl.Identity(mesh.geometric_dimension)
    stress += state.problem.viscosity * (gradient + ufl.transpose(gradient))
    traction = -ufl.dot(stress, fd.FacetNormal(mesh))
    values = {}
    for name in state.problem.force_boundaries:
        ds = named_boundary_measure(state.ds, (name,), state.problem.boundary_tags)
        values.update(
            {f"{name}_{axis}": value for axis, value in _components(traction, ds).items()}
        )
    return values


def _divergence_l2(state: DiagnosticState) -> float:
    return math.sqrt(fd.assemble(ufl.div(state.velocity) ** 2 * state.dx))


def _velocity_error(state: DiagnosticState) -> Expr:
    x = fd.SpatialCoordinate(state.velocity.function_space().mesh())
    return state.velocity - state.problem.exact_velocity(x, state.velocity_time)


def _velocity_error_l2(state: DiagnosticState) -> float:
    error = _velocity_error(state)
    return math.sqrt(fd.assemble(ufl.inner(error, error) * state.error_dx))


def _velocity_gradient_error_l2(state: DiagnosticState) -> float:
    error = ufl.grad(_velocity_error(state))
    return math.sqrt(fd.assemble(ufl.inner(error, error) * state.error_dx))


def _pressure_error_l2(state: DiagnosticState) -> float:
    x = fd.SpatialCoordinate(state.problem.mesh)
    error = state.pressure - state.problem.exact_pressure(x, state.pressure_time)
    if state.pressure_null_space:
        area = fd.assemble(fd.Constant(1.0) * state.error_dx)
        error = error - fd.Constant(fd.assemble(error * state.error_dx) / area)
    return math.sqrt(fd.assemble(error**2 * state.error_dx))


DIAGNOSTICS: dict[str, Diagnostic] = {
    "kinetic_energy": Diagnostic(_kinetic_energy, needs_exact_solution=False),
    "subscale_kinetic_energy": Diagnostic(_subscale_kinetic_energy, needs_exact_solution=False),
    "linear_momentum": Diagnostic(_linear_momentum, needs_exact_solution=False),
    "angular_momentum": Diagnostic(_angular_momentum, needs_exact_solution=False),
    "force": Diagnostic(_force, needs_exact_solution=False),
    "divergence_l2": Diagnostic(_divergence_l2, needs_exact_solution=False),
    "velocity_error_l2": Diagnostic(_velocity_error_l2, needs_exact_solution=True),
    "velocity_gradient_error_l2": Diagnostic(
        _velocity_gradient_error_l2, needs_exact_solution=True
    ),
    "pressure_error_l2": Diagnostic(_pressure_error_l2, needs_exact_solution=True),
}


def check_diagnostics_supported(names: Sequence[str], problem: Problem) -> None:
    """Reject diagnostics whose required exact solution or force boundaries are absent."""
    if "force" in names and not problem.force_boundaries:
        raise ValueError("output.diagnostics: force requires the case to declare force_boundaries")
    needing = [name for name in names if DIAGNOSTICS[name].needs_exact_solution]
    if needing and problem.exact_velocity is None:
        raise ValueError(
            f"output.diagnostics: {', '.join(needing)} need an exact solution, "
            "which this case does not provide"
        )


def evaluate_diagnostics(names: Sequence[str], state: DiagnosticState) -> dict[str, float]:
    """The configured diagnostics of the current state, in configuration order."""
    values = {}
    for name in names:
        value = DIAGNOSTICS[name].evaluate(state)
        if isinstance(value, dict):
            values.update({f"{name}_{component}": entry for component, entry in value.items()})
        else:
            values[name] = value
    return values


def diagnostic_columns(names: Sequence[str], problem: Problem) -> tuple[str, ...]:
    """Scalar column names after expanding vectors and the case's measured boundaries."""
    axes = "xyz"[: problem.mesh.geometric_dimension]
    columns = []
    for name in names:
        if name == "linear_momentum" or (name == "angular_momentum" and len(axes) == 3):
            columns.extend(f"{name}_{axis}" for axis in axes)
        elif name == "force":
            columns.extend(
                f"force_{boundary}_{axis}" for boundary in problem.force_boundaries for axis in axes
            )
        else:
            columns.append(name)
    return tuple(columns)
