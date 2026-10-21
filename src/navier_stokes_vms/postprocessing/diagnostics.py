"""Registry of named per-step diagnostics; error norms use a rule of twice the common degree."""

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import firedrake as fd
import ufl
from ufl.core.expr import Expr
from ufl.measure import Measure

from navier_stokes_vms.problem.problem_definition import Problem


@dataclass(frozen=True)
class DiagnosticState:
    """What a diagnostic reads after a step: endpoint velocity and solve-level pressure."""

    problem: Problem
    velocity: fd.Function
    pressure: fd.Function
    velocity_time: fd.Constant
    pressure_time: fd.Constant
    dx: Measure
    error_dx: Measure
    pressure_null_space: bool


@dataclass(frozen=True)
class Diagnostic:
    """A scalar functional of the state and whether it needs the case's exact solution."""

    evaluate: Callable[[DiagnosticState], float]
    needs_exact_solution: bool


def _kinetic_energy(state: DiagnosticState) -> float:
    return fd.assemble(ufl.inner(state.velocity, state.velocity) / 2 * state.dx)


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
    x = fd.SpatialCoordinate(state.pressure.function_space().mesh())
    error = state.pressure - state.problem.exact_pressure(x, state.pressure_time)
    if state.pressure_null_space:
        area = fd.assemble(fd.Constant(1.0) * state.error_dx)
        error = error - fd.Constant(fd.assemble(error * state.error_dx) / area)
    return math.sqrt(fd.assemble(error**2 * state.error_dx))


DIAGNOSTICS: dict[str, Diagnostic] = {
    "kinetic_energy": Diagnostic(_kinetic_energy, needs_exact_solution=False),
    "divergence_l2": Diagnostic(_divergence_l2, needs_exact_solution=False),
    "velocity_error_l2": Diagnostic(_velocity_error_l2, needs_exact_solution=True),
    "velocity_gradient_error_l2": Diagnostic(
        _velocity_gradient_error_l2, needs_exact_solution=True
    ),
    "pressure_error_l2": Diagnostic(_pressure_error_l2, needs_exact_solution=True),
}


def check_diagnostics_supported(names: Sequence[str], problem: Problem) -> None:
    """Reject error diagnostics for a case without an exact solution."""
    if problem.exact_velocity is not None:
        return
    needing = [name for name in names if DIAGNOSTICS[name].needs_exact_solution]
    if needing:
        raise ValueError(
            f"output.diagnostics: {', '.join(needing)} need an exact solution, "
            "which this case does not provide"
        )


def evaluate_diagnostics(names: Sequence[str], state: DiagnosticState) -> dict[str, float]:
    """The configured diagnostics of the current state, in configuration order."""
    return {name: DIAGNOSTICS[name].evaluate(state) for name in names}
