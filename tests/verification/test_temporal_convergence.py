"""Temporal convergence on a rotation that lies in every FE space (plan V2)."""

import logging
from pathlib import Path

import pytest
from study import finest_rate, run_study

pytestmark = pytest.mark.slow

LOGGER = logging.getLogger(__name__)

STEPS = ("--steps", "20", "40", "80", "160")
ERRORS = ("velocity_error_l2", "pressure_error_l2")

DYNAMIC_ASGS = (
    'method = "vms"\nsubscale_space = "asgs"\nsubscale_dynamics = "dynamic"\n'
    'scale_splitting = "linear"'
)
VARIANTS = {
    "galerkin": (2, 1, 'method = "none"'),
    "dyn_lin_asgs_p1": (1, 1, DYNAMIC_ASGS),
    "dyn_lin_asgs_p2": (2, 2, DYNAMIC_ASGS),
}


def base(variant: str, scheme: str, startup: str) -> str:
    velocity_degree, pressure_degree, stabilization = VARIANTS[variant]
    return f"""case = "manufactured_temporal"

[discretization]
velocity_degree = {velocity_degree}
pressure_degree = {pressure_degree}

[time]
scheme = "{scheme}"
dt = 0.05
t_end = 1.0
startup = "{startup}"

[stabilization]
{stabilization}

[output]
vtk_every = 0
diagnostics = ["velocity_error_l2", "pressure_error_l2"]
"""


@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize(
    ("scheme", "startup", "order"),
    [
        ("bdf1", "lower_order", 1),
        ("bdf2", "lower_order", 2),
        ("bdf3", "exact_history", 3),
        ("crank_nicolson", "lower_order", 2),
    ],
)
def test_velocity_and_pressure_reach_the_scheme_order(
    tmp_path: Path, variant: str, scheme: str, startup: str, order: int
) -> None:
    table = run_study(tmp_path, base(variant, scheme, startup), f"{variant}_{scheme}", *STEPS)
    LOGGER.info("V2 %s %s (%s): %s", variant, scheme, startup, table)
    for name in ERRORS:
        assert finest_rate(table, name) >= order - 0.2, name


@pytest.mark.parametrize("variant", VARIANTS)
def test_bdf3_with_lower_order_startup_is_at_least_second_order(
    tmp_path: Path, variant: str
) -> None:
    table = run_study(
        tmp_path, base(variant, "bdf3", "lower_order"), f"{variant}_bdf3_lower_order", *STEPS
    )
    LOGGER.info("V2 %s bdf3 (lower_order): %s", variant, table)
    for name in ERRORS:
        assert finest_rate(table, name) >= 1.8, name
