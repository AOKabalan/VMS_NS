"""Spatial convergence on a solution linear in time, so time stepping is exact (plan V1)."""

import logging
from pathlib import Path

import pytest
from study import finest_rate, run_study

pytestmark = pytest.mark.slow

LOGGER = logging.getLogger(__name__)

GALERKIN = """case = "manufactured_spatial"

[discretization]
velocity_degree = 2
pressure_degree = 1

[time]
scheme = "bdf1"
dt = 0.05
t_end = 0.1

[stabilization]
method = "none"

[output]
vtk_every = 0
diagnostics = ["velocity_error_l2", "velocity_gradient_error_l2", "pressure_error_l2"]
"""

EXPECTED = {"velocity_error_l2": 3, "velocity_gradient_error_l2": 2, "pressure_error_l2": 2}
ASGS_EXPECTED = {
    1: {"velocity_error_l2": 2, "velocity_gradient_error_l2": 1, "pressure_error_l2": 1},
    2: {"velocity_error_l2": 3, "velocity_gradient_error_l2": 2, "pressure_error_l2": 2},
}


def test_galerkin_taylor_hood_reaches_the_expected_spatial_rates(tmp_path: Path) -> None:
    table = run_study(tmp_path, GALERKIN, "galerkin_p2p1", "--cells", "8x8", "16x16", "32x32")
    LOGGER.info("V1 Galerkin P2/P1: %s", table)
    for name, expected in EXPECTED.items():
        assert finest_rate(table, name) >= expected - 0.2, name


@pytest.mark.parametrize("degree", [1, 2])
def test_static_linear_asgs_reaches_the_expected_spatial_rates(tmp_path: Path, degree: int) -> None:
    base = GALERKIN.replace("velocity_degree = 2", f"velocity_degree = {degree}")
    base = base.replace("pressure_degree = 1", f"pressure_degree = {degree}")
    base = base.replace(
        'method = "none"',
        'method = "vms"\nsubscale_space = "asgs"\nsubscale_dynamics = "static"\n'
        'scale_splitting = "linear"',
    )
    table = run_study(
        tmp_path, base, f"asgs_p{degree}p{degree}", "--cells", "8x8", "16x16", "32x32"
    )
    LOGGER.info("V1 Sta-Lin-ASGS P%d/P%d: %s", degree, degree, table)
    for name, expected in ASGS_EXPECTED[degree].items():
        assert finest_rate(table, name) >= expected - 0.2, name
