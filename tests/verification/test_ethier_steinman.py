"""Coarse 3D BDF2 refinement and equality of serial and two-rank diagnostics (phase 9)."""

import csv
import logging
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    StabilizationOptions,
    TimeOptions,
    config_to_toml,
)

LOGGER = logging.getLogger(__name__)
RUN = Path(__file__).resolve().parents[2] / "scripts" / "run_simulation.py"
ERROR_DIAGNOSTICS = (
    "velocity_error_l2",
    "velocity_gradient_error_l2",
    "pressure_error_l2",
)
DIAGNOSTICS = (
    "kinetic_energy",
    "divergence_l2",
    *ERROR_DIAGNOSTICS,
    "linear_momentum",
    "angular_momentum",
    "force",
)


def run(folder: Path, cells: int, method: str, ranks: int) -> list[dict[str, str]]:
    folder.mkdir()
    names = DIAGNOSTICS + (("subscale_kinetic_energy",) if method == "vms" else ())
    config = Config(
        case="ethier_steinman_3d",
        run_name="smoke",
        mesh=MeshOptions(cells=(cells,) * 3),
        discretization=DiscretizationOptions(velocity_degree=2 if method == "none" else 1),
        time=TimeOptions(scheme="bdf2", dt=0.01, t_end=0.03),
        stabilization=StabilizationOptions(method=method),
        output=OutputOptions(vtk_every=0, diagnostics=names),
    )
    path = folder / "run.toml"
    path.write_text(config_to_toml(config), encoding="utf-8")
    command = [sys.executable, str(RUN), str(path)]
    if ranks > 1:
        command = ["mpiexec", "-n", str(ranks), *command]
    result = subprocess.run(
        command,
        cwd=folder,
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
    )
    assert result.returncode == 0, result.stderr[-6000:]
    output = folder / "results" / config.case / config.run_name
    with (output / "time_steps.csv").open(encoding="utf-8") as stream:
        steps = list(csv.DictReader(stream))
    assert [step["scheme_used"] for step in steps] == ["bdf1", "bdf2", "bdf2"]
    assert all(step["converged"] == "True" for step in steps)
    with (output / "diagnostics.csv").open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


@pytest.mark.slow
@pytest.mark.parametrize("method", ["none", "vms"], ids=["galerkin_p2p1", "dyn_nl_oss_p1p1"])
def test_ethier_refinement_and_mpi_diagnostics(tmp_path: Path, method: str) -> None:
    errors = []
    largest_relative_difference = 0.0
    largest_absolute_difference = 0.0
    differences = []
    absolute_tolerance = 1e-13 if method == "none" else 1e-4
    force_scales = {}
    force_differences = {}
    for cells in (4, 8):
        serial = run(tmp_path / f"cells{cells}_serial", cells, method, 1)
        parallel = run(tmp_path / f"cells{cells}_mpi", cells, method, 2)
        force_scales[cells] = max(
            abs(float(value))
            for row in serial
            for name, value in row.items()
            if name.startswith("force_")
        )
        force_differences[cells] = 0.0
        for expected, actual in zip(serial, parallel, strict=True):
            assert actual.keys() == expected.keys()
            for name, value in expected.items():
                if value == "":
                    assert actual[name] == ""
                else:
                    first, second = float(value), float(actual[name])
                    differences.append((cells, expected["step"], name, first, second))
                    largest_absolute_difference = max(
                        largest_absolute_difference, abs(second - first)
                    )
                    if name.startswith("force_"):
                        force_differences[cells] = max(
                            force_differences[cells], abs(second - first)
                        )
                    if abs(first) > 1e-12:
                        largest_relative_difference = max(
                            largest_relative_difference, abs(second - first) / abs(first)
                        )
        errors.append(float(serial[-1]["velocity_error_l2"]))
    drop = errors[0] / errors[1]
    LOGGER.info(
        "Ethier %s: errors %s, drop %.8f, rate %.8f; largest serial/MPI "
        "relative difference %.3e, absolute %.3e",
        method,
        errors,
        drop,
        math.log2(drop),
        largest_relative_difference,
        largest_absolute_difference,
    )
    LOGGER.info(
        "Ethier %s: largest force differences %s, serial force scales %s; "
        "force difference decreases with refinement: %s",
        method,
        force_differences,
        force_scales,
        force_differences[8] < force_differences[4],
    )
    assert drop > 1
    failures = []
    for cells, step, name, first, second in differences:
        if method == "vms" and name.startswith("force_"):
            tolerance = 1e-3 * force_scales[cells]
        else:
            relative_tolerance = (
                1e-4 if method == "vms" else 1e-8 if name in ERROR_DIAGNOSTICS else 1e-10
            )
            tolerance = max(absolute_tolerance, relative_tolerance * abs(first))
        if abs(second - first) > tolerance:
            failures.append((cells, step, name, first, second, tolerance))
    assert not failures, failures
