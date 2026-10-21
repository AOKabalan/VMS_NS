"""Time-scheme coefficients, startup sequences and data level (plan U1, spec: Time schemes)."""

from collections.abc import Callable

import pytest

from navier_stokes_vms.config import Startup, TimeScheme
from navier_stokes_vms.discretization.time_schemes import (
    StepCoefficients,
    scheme_at_step,
    step_coefficients,
)

DT = 0.1
T_NEW = 0.7
SHIFT = 0.3


def monomial(power: int) -> tuple[Callable[[float], float], Callable[[float], float]]:
    def value(t: float) -> float:
        return (t - SHIFT) ** power

    def derivative(t: float) -> float:
        return power * (t - SHIFT) ** (power - 1) if power else 0.0

    return value, derivative


def discrete_derivative(
    coefficients: StepCoefficients, solve_value: float, g: Callable[[float], float]
) -> float:
    history = sum(
        weight * g(T_NEW - (j + 1) * DT) for j, weight in enumerate(coefficients.history_weights)
    )
    return (coefficients.alpha0 * solve_value - history) / DT


@pytest.mark.parametrize(("scheme", "order"), [("bdf1", 1), ("bdf2", 2), ("bdf3", 3)])
def test_bdf_differentiates_polynomials_up_to_its_order_exactly(
    scheme: TimeScheme, order: int
) -> None:
    coefficients = step_coefficients(scheme, "lower_order", step=order)
    assert coefficients.scheme == scheme
    assert coefficients.solve_level == 1.0
    for power in range(order + 2):
        g, dg = monomial(power)
        error = discrete_derivative(coefficients, g(T_NEW), g) - dg(T_NEW)
        if power <= order:
            assert abs(error) < 1e-12
        else:
            assert abs(error) > 1e-6


def test_crank_nicolson_differentiates_quadratics_at_the_half_level() -> None:
    coefficients = step_coefficients("crank_nicolson", "lower_order", step=2)
    assert coefficients.scheme == "crank_nicolson"
    assert coefficients.solve_level == 0.5
    t_old = T_NEW - DT
    t_half = t_old + coefficients.solve_level * DT
    for power in range(4):
        g, dg = monomial(power)
        half_state = (
            g(T_NEW) - coefficients.extrapolation_weight_old * g(t_old)
        ) / coefficients.extrapolation_weight_solve
        error = discrete_derivative(coefficients, half_state, g) - dg(t_half)
        if power <= 2:
            assert abs(error) < 1e-12
        else:
            assert abs(error) > 1e-6


def test_crank_nicolson_extrapolation_reproduces_linear_functions() -> None:
    coefficients = step_coefficients("crank_nicolson", "exact_history", step=3)
    t_old = T_NEW - DT
    t_half = t_old + coefficients.solve_level * DT
    for power in range(3):
        g, _ = monomial(power)
        endpoint = coefficients.extrapolation_weight_solve * g(
            t_half
        ) + coefficients.extrapolation_weight_old * g(t_old)
        if power <= 1:
            assert endpoint == pytest.approx(g(T_NEW), abs=1e-14)
        else:
            assert abs(endpoint - g(T_NEW)) > 1e-6


@pytest.mark.parametrize("startup", ["lower_order", "exact_history"])
def test_crank_nicolson_data_level_is_the_endpoint_average(startup: Startup) -> None:
    coefficients = step_coefficients("crank_nicolson", startup, step=2)
    assert (coefficients.data_weight_new, coefficients.data_weight_old) == (0.5, 0.5)
    t_old = T_NEW - DT
    for power in range(4):
        g, _ = monomial(power)
        data = coefficients.data_weight_new * g(T_NEW) + coefficients.data_weight_old * g(t_old)
        endpoint = (
            coefficients.extrapolation_weight_solve * data
            + coefficients.extrapolation_weight_old * g(t_old)
        )
        assert endpoint == pytest.approx(g(T_NEW), abs=1e-14)


@pytest.mark.parametrize(
    ("scheme", "startup", "expected"),
    [
        ("bdf1", "lower_order", ["bdf1"] * 5),
        ("bdf2", "lower_order", ["bdf1"] + ["bdf2"] * 4),
        ("bdf3", "lower_order", ["bdf1", "bdf2"] + ["bdf3"] * 3),
        ("crank_nicolson", "lower_order", ["bdf1"] + ["crank_nicolson"] * 4),
        ("bdf1", "exact_history", ["bdf1"] * 5),
        ("bdf2", "exact_history", ["bdf2"] * 5),
        ("bdf3", "exact_history", ["bdf3"] * 5),
        ("crank_nicolson", "exact_history", ["bdf1"] + ["crank_nicolson"] * 4),
    ],
)
def test_startup_sequence(scheme: TimeScheme, startup: Startup, expected: list[str]) -> None:
    assert [scheme_at_step(scheme, startup, step) for step in range(1, 6)] == expected


@pytest.mark.parametrize(
    ("scheme", "levels"), [("bdf1", 1), ("bdf2", 2), ("bdf3", 3), ("crank_nicolson", 1)]
)
def test_exact_history_first_step_needs_levels_back_to_t_minus_two(
    scheme: TimeScheme, levels: int
) -> None:
    assert len(step_coefficients(scheme, "exact_history", step=1).history_weights) == levels


@pytest.mark.parametrize("scheme", ["bdf1", "bdf2", "bdf3", "crank_nicolson"])
@pytest.mark.parametrize("startup", ["lower_order", "exact_history"])
def test_every_step_is_consistent(scheme: TimeScheme, startup: Startup) -> None:
    for step in range(1, 5):
        coefficients = step_coefficients(scheme, startup, step)
        assert sum(coefficients.history_weights) == pytest.approx(coefficients.alpha0)
        assert coefficients.data_weight_new + coefficients.data_weight_old == 1.0
        assert (
            coefficients.extrapolation_weight_solve + coefficients.extrapolation_weight_old == 1.0
        )
        if coefficients.scheme != "crank_nicolson":
            assert (coefficients.data_weight_new, coefficients.extrapolation_weight_solve) == (
                1.0,
                1.0,
            )


def test_steps_are_numbered_from_one() -> None:
    with pytest.raises(ValueError, match="step = 0: steps are numbered from 1"):
        step_coefficients("bdf2", "lower_order", step=0)
