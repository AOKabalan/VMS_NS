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


DYN_NL_OSS = GALERKIN.replace("velocity_degree = 2", "velocity_degree = 1").replace(
    'method = "none"',
    'method = "vms"\nsubscale_space = "oss"\nsubscale_dynamics = "dynamic"\n'
    'scale_splitting = "nonlinear"',
)


EXACT_ENERGY = math.pi**2 * math.exp(-4 * VISCOSITY * T_END)


def velocity_errors(table: list[dict[str, str]]) -> list[float]:
    return [float(row["velocity_error_l2"]) for row in table]


def energy_errors(table: list[dict[str, str]]) -> list[float]:
    return [abs(float(row["kinetic_energy"]) - EXACT_ENERGY) / EXACT_ENERGY for row in table]


def test_galerkin_error_drops_and_energy_decays_correctly(tmp_path: Path) -> None:
    table = run_study(tmp_path, GALERKIN, "galerkin_p2p1", "--cells", "16x16", "32x32")
    LOGGER.info("V4 galerkin_p2p1: %s", table)
    coarse, fine = velocity_errors(table)
    assert coarse / fine >= 3
    assert energy_errors(table)[-1] <= 0.01


def test_dyn_nl_oss_error_drops_and_energy_converges(tmp_path: Path) -> None:
    """P1/P1 meets the 1 % energy bound on mesh 64, not 32 (plan V4)."""
    table = run_study(tmp_path, DYN_NL_OSS, "dyn_nl_oss_p1p1", "--cells", "16x16", "32x32", "64x64")
    velocity, energy = velocity_errors(table), energy_errors(table)
    LOGGER.info("V4 dyn_nl_oss_p1p1: %s; relative energy errors %s", table, energy)
    assert velocity[0] / velocity[1] >= 3
    assert energy[2] <= 0.01
    assert energy[1] / energy[2] >= 3
