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
MESHES = ("8x8", "16x16", "32x32")
# The static nonlinear local update has no convergent iteration on mesh 8 for P1/P1 (spec: The
# eight cases, the paragraph on static nonlinear subscales); plan V1 moves this variant up a level.
STATIC_NONLINEAR_P1_MESHES = ("16x16", "32x32", "64x64")


def test_galerkin_taylor_hood_reaches_the_expected_spatial_rates(tmp_path: Path) -> None:
    table = run_study(tmp_path, GALERKIN, "galerkin_p2p1", "--cells", *MESHES)
    LOGGER.info("V1 Galerkin P2/P1: %s", table)
    for name, expected in EXPECTED.items():
        assert finest_rate(table, name) >= expected - 0.2, name


@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("degree", [1, 2])
def test_static_asgs_reaches_the_expected_spatial_rates(
    tmp_path: Path, degree: int, splitting: str
) -> None:
    base = GALERKIN.replace("velocity_degree = 2", f"velocity_degree = {degree}")
    base = base.replace("pressure_degree = 1", f"pressure_degree = {degree}")
    base = base.replace(
        'method = "none"',
        'method = "vms"\nsubscale_space = "asgs"\nsubscale_dynamics = "static"\n'
        f'scale_splitting = "{splitting}"',
    )
    name = f"asgs_{splitting}_p{degree}p{degree}"
    meshes = STATIC_NONLINEAR_P1_MESHES if (degree, splitting) == (1, "nonlinear") else MESHES
    table = run_study(tmp_path, base, name, "--cells", *meshes)
    LOGGER.info("V1 static %s ASGS P%d/P%d: %s", splitting, degree, degree, table)
    for name, expected in ASGS_EXPECTED[degree].items():
        assert finest_rate(table, name) >= expected - 0.2, name
