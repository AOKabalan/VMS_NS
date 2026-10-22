"""Kovasznay marched to steady state with BDF1: equivalences and mesh convergence (plan V3)."""

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
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.discretization.time_schemes import StepCoefficients
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.solver.nonlinear_solver import relative_change
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
    """Run to t_end and require a steady state between the last two time levels.

    Steady means the endpoint velocity, and for dynamic subscales the stored ũ, changed by less
    than STEADY_TOLERANCE in the last step. Every step's outer convergence is a separate check.
    """
    folder.mkdir()
    stepper = TimeStepper(config, problem, folder, LOGGER)
    levels = {"velocity": lambda: stepper.velocity}
    if stepper._subscales is not None:
        levels["subscale"] = lambda: stepper._subscales.latest
    previous = {name: fd.Function(level().function_space()) for name, level in levels.items()}
    advance = stepper._advance

    def recording_advance(coefficients: StepCoefficients) -> None:
        for name, level in levels.items():
            previous[name].assign(level())
        advance(coefficients)

    stepper._advance = recording_advance
    stepper.run()
    with (folder / "time_steps.csv").open(encoding="utf-8") as stream:
        assert all(row["converged"] == "True" for row in csv.DictReader(stream)), config.run_name
    changes = {
        name: relative_change(level(), previous[name], stepper._dx)
        for name, level in levels.items()
    }
    LOGGER.info("V3 %s: change over the last step %s", config.run_name, changes)
    assert max(changes.values()) < STEADY_TOLERANCE, (config.run_name, changes)
    return stepper


def final_diagnostics(folder: Path) -> dict[str, float]:
    with (folder / "diagnostics.csv").open(encoding="utf-8") as stream:
        return {
            key: float(value)
            for key, value in list(csv.DictReader(stream))[-1].items()
            if value != ""
        }


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


def vms(
    cells: tuple[int, int],
    degree: int,
    dynamics: str,
    splitting: str,
    linearization: str,
    update: str = "accelerated",
    space: str = "asgs",
) -> Config:
    base = galerkin(cells, linearization)
    config = replace(
        base,
        run_name=(f"{space}_p{degree}_{dynamics}_{splitting}_{linearization}_{update}_{cells[0]}"),
        discretization=DiscretizationOptions(velocity_degree=degree, pressure_degree=degree),
        stabilization=StabilizationOptions(
            subscale_space=space, subscale_dynamics=dynamics, scale_splitting=splitting
        ),
        solver=replace(base.solver, subscale_update=update),
    )
    validate_config(config)
    return config


def iteration_totals(folder: Path) -> tuple[int, int, int]:
    with (folder / "time_steps.csv").open(encoding="utf-8") as stream:
        steps = list(csv.DictReader(stream))
    outer = [int(row["outer_iterations"]) for row in steps]
    return sum(outer), max(outer), sum(int(row["inner_iterations_total"]) for row in steps)


@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("degree", [1, 2])
def test_vms_steady_state_is_independent_of_dynamics_linearization_and_update(
    tmp_path: Path, degree: int, splitting: str, space: str
) -> None:
    label = f"{splitting} {space.upper()} P{degree}/P{degree}"
    runs = [
        (dynamics, linearization, "accelerated")
        for dynamics in ("static", "dynamic")
        for linearization in ("newton", "picard")
    ]
    if splitting == "nonlinear":
        runs.append(("static", "newton", "fixed_point"))
    errors: dict[str, list[dict[str, float]]] = {"static": [], "dynamic": []}
    problem = build_problem(vms(MESH_16, degree, "static", splitting, "newton", space=space))
    dx = volume_measure(problem.mesh, 3 * degree)
    velocities = {}
    for run in runs:
        dynamics, linearization, update = run
        config = vms(MESH_16, degree, dynamics, splitting, linearization, update, space)
        folder = tmp_path / config.run_name
        velocities[run] = march(config, problem, folder).velocity
        LOGGER.info(
            "V3 iterations %s: outer total %d (max per step %d), inner total %d",
            config.run_name,
            *iteration_totals(folder),
        )
        if linearization == "newton" and update == "accelerated":
            errors[dynamics].append(final_diagnostics(folder))
    reference = velocities["static", "newton", "accelerated"]
    for run, velocity in velocities.items():
        difference = relative_change(velocity, reference, dx)
        LOGGER.info("V3 %s mesh 16: %s vs static newton accelerated %.3e", label, run, difference)
        assert difference <= 1e-8, run
    newton_picard = relative_change(
        velocities["dynamic", "picard", "accelerated"],
        velocities["dynamic", "newton", "accelerated"],
        dx,
    )
    LOGGER.info("V3 dynamic %s mesh 16: Newton/Picard %.3e", label, newton_picard)
    assert newton_picard <= 1e-8
    problem = build_problem(vms(MESH_32, degree, "static", splitting, "newton", space=space))
    for dynamics in errors:
        config = vms(MESH_32, degree, dynamics, splitting, "newton", space=space)
        folder = tmp_path / config.run_name
        march(config, problem, folder)
        LOGGER.info(
            "V3 iterations %s: outer total %d (max per step %d), inner total %d",
            config.run_name,
            *iteration_totals(folder),
        )
        errors[dynamics].append(final_diagnostics(folder))
    for dynamics, pair in errors.items():
        LOGGER.info("V3 %s %s steady errors: %s", dynamics, label, pair)
        for name in ("velocity_error_l2", "pressure_error_l2"):
            assert pair[1][name] < pair[0][name], (dynamics, name)
