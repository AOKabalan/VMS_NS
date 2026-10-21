"""Temporal convergence on a rotation that lies in every FE space (plan V2)."""

import logging
from pathlib import Path

import pytest
from study import finest_rate, run_study

pytestmark = pytest.mark.slow

LOGGER = logging.getLogger(__name__)

STEPS = ("--steps", "20", "40", "80", "160")
ERRORS = ("velocity_error_l2", "pressure_error_l2")


def galerkin(scheme: str, startup: str) -> str:
    return f"""case = "manufactured_temporal"

[discretization]
velocity_degree = 2
pressure_degree = 1

[time]
scheme = "{scheme}"
dt = 0.05
t_end = 1.0
startup = "{startup}"

[stabilization]
method = "none"

[output]
vtk_every = 0
diagnostics = ["velocity_error_l2", "pressure_error_l2"]
"""


@pytest.mark.parametrize(
    ("scheme", "startup", "order"),
    [
        ("bdf1", "lower_order", 1),
        ("bdf2", "lower_order", 2),
        ("bdf3", "exact_history", 3),
        ("crank_nicolson", "lower_order", 2),
    ],
)
def test_galerkin_velocity_and_pressure_reach_the_scheme_order(
    tmp_path: Path, scheme: str, startup: str, order: int
) -> None:
    table = run_study(tmp_path, galerkin(scheme, startup), f"galerkin_{scheme}", *STEPS)
    LOGGER.info("V2 Galerkin %s (%s): %s", scheme, startup, table)
    for name in ERRORS:
        assert finest_rate(table, name) >= order - 0.2, name


def test_bdf3_with_lower_order_startup_is_at_least_second_order(tmp_path: Path) -> None:
    table = run_study(
        tmp_path, galerkin("bdf3", "lower_order"), "galerkin_bdf3_lower_order", *STEPS
    )
    LOGGER.info("V2 Galerkin bdf3 (lower_order): %s", table)
    for name in ERRORS:
        assert finest_rate(table, name) >= 1.8, name
