"""Kovasznay marched to steady state with BDF1: Newton = Picard, mesh convergence (plan V3)."""

import csv
import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import pytest

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
    validate_config,
)
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.solver.time_stepper import TimeStepper

pytestmark = pytest.mark.slow

LOGGER = logging.getLogger(__name__)

STEADY_TOLERANCE = 1e-10
MESH_16 = (24, 32)
MESH_32 = (48, 64)


def galerkin(cells: tuple[int, int], linearization: str) -> Config:
    config = Config(
        case="kovasznay",
        run_name=f"galerkin_{linearization}_{cells[0]}x{cells[1]}",
        mesh=MeshOptions(cells=cells),
        discretization=DiscretizationOptions(velocity_degree=2, pressure_degree=1),
        time=TimeOptions(scheme="bdf1", dt=5.0, t_end=100.0),
        stabilization=StabilizationOptions(method="none"),
        solver=replace(SolverOptions(), linearization=linearization, inner_max_iterations=100),
        output=OutputOptions(vtk_every=0, diagnostics=("velocity_error_l2", "pressure_error_l2")),
    )
    validate_config(config)
    return config


def march(config: Config, problem: Problem, folder: Path) -> TimeStepper:
    folder.mkdir()
    stepper = TimeStepper(config, problem, folder, LOGGER)
    stepper.run()
    with (folder / "outer_iterations.csv").open(encoding="utf-8") as stream:
        last_change = float(list(csv.DictReader(stream))[-1]["velocity_change"])
    assert last_change < STEADY_TOLERANCE, config.run_name
    return stepper


def final_diagnostics(folder: Path) -> dict[str, float]:
    with (folder / "diagnostics.csv").open(encoding="utf-8") as stream:
        return {key: float(value) for key, value in list(csv.DictReader(stream))[-1].items()}


def test_galerkin_newton_and_picard_reach_the_same_steady_velocity(tmp_path: Path) -> None:
    newton_config = galerkin(MESH_16, "newton")
    problem = build_problem(newton_config)
    newton = march(newton_config, problem, tmp_path / "newton")
    picard = march(galerkin(MESH_16, "picard"), problem, tmp_path / "picard")
    difference = fd.norm(picard.velocity - newton.velocity) / fd.norm(newton.velocity)
    LOGGER.info("V3 Galerkin: relative Newton/Picard velocity difference %.3e", difference)
    assert difference <= 1e-8


def test_galerkin_steady_error_drops_from_mesh_16_to_32(tmp_path: Path) -> None:
    errors = []
    for cells in (MESH_16, MESH_32):
        config = galerkin(cells, "newton")
        folder = tmp_path / config.run_name
        march(config, build_problem(config), folder)
        errors.append(final_diagnostics(folder))
    LOGGER.info("V3 Galerkin steady errors, mesh 16 then 32: %s", errors)
    for name in ("velocity_error_l2", "pressure_error_l2"):
        assert errors[1][name] < errors[0][name], name
