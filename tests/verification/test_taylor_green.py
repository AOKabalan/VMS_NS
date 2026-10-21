"""Decaying Taylor–Green vortex with BDF2 to T = 1 (plan V4)."""

import logging
import math
from pathlib import Path

import pytest
from study import run_study

pytestmark = pytest.mark.slow

LOGGER = logging.getLogger(__name__)

VISCOSITY = 0.01
T_END = 1.0

GALERKIN = """case = "taylor_green_2d"

[discretization]
velocity_degree = 2
pressure_degree = 1

[time]
scheme = "bdf2"
dt = 0.025
t_end = 1.0

[stabilization]
method = "none"

[output]
vtk_every = 0
diagnostics = ["velocity_error_l2", "kinetic_energy"]
"""


def test_galerkin_error_drops_and_energy_decays_correctly(tmp_path: Path) -> None:
    table = run_study(tmp_path, GALERKIN, "galerkin_p2p1", "--cells", "16x16", "32x32")
    LOGGER.info("V4 Galerkin P2/P1: %s", table)
    coarse, fine = (float(row["velocity_error_l2"]) for row in table)
    assert coarse / fine >= 3
    exact_energy = math.pi**2 * math.exp(-4 * VISCOSITY * T_END)
    assert float(table[-1]["kinetic_energy"]) == pytest.approx(exact_energy, rel=0.01)
