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


def test_galerkin_taylor_hood_reaches_the_expected_spatial_rates(tmp_path: Path) -> None:
    table = run_study(tmp_path, GALERKIN, "galerkin_p2p1", "--cells", "8x8", "16x16", "32x32")
    LOGGER.info("V1 Galerkin P2/P1: %s", table)
    for name, expected in EXPECTED.items():
        assert finest_rate(table, name) >= expected - 0.2, name
