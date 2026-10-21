"""Galerkin runs end to end, iteration records, field output and failure handling (phase 4a)."""

import csv
import importlib.util
import json
import math
import re
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import Config, MeshOptions, SolverOptions, load_config
from navier_stokes_vms.discretization.function_spaces import (
    constant_pressure_null_space,
    mixed_space,
)
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.discretization.time_schemes import scheme_at_step
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.problem.boundary_conditions import dirichlet_bc
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.solver.nonlinear_solver import InnerSolver, InnerSolveResult

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "run_simulation.py"
_SPEC = importlib.util.spec_from_file_location("run_simulation", SCRIPT)
RUN_SIMULATION = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(RUN_SIMULATION)

ALL_DIAGNOSTICS = (
    "kinetic_energy, divergence_l2, velocity_error_l2, velocity_gradient_error_l2, "
    "pressure_error_l2"
)
RUNS = {
    "manufactured_spatial": (
        "cells = [4, 4]",
        'scheme = "bdf1"\ndt = 0.1\nt_end = 0.2',
        "",
        1,
    ),
    "manufactured_temporal": (
        "",
        'scheme = "crank_nicolson"\ndt = 0.25\nt_end = 0.75\nstartup = "exact_history"',
        'linearization = "picard"',
        2,
    ),
    "kovasznay": (
        "cells = [6, 8]",
        'scheme = "bdf2"\ndt = 0.5\nt_end = 1.0',
        "picard_iterations_before_newton = 2",
        0,
    ),
    "taylor_green_2d": (
        "cells = [8, 8]",
        'scheme = "bdf3"\ndt = 0.1\nt_end = 0.3\nstartup = "exact_history"',
        "",
        3,
    ),
}
FILES = [
    "config_resolved.toml",
    "diagnostics.csv",
    "inner_iterations.csv",
    "outer_iterations.csv",
    "run.log",
    "run_info.json",
    "time_steps.csv",
]


def config_text(case: str) -> str:
    mesh, time, solver, vtk_every = RUNS[case]
    diagnostics = ", ".join(f'"{name}"' for name in ALL_DIAGNOSTICS.split(", "))
    return (
        f'case = "{case}"\n\n[mesh]\n{mesh}\n\n[discretization]\nvelocity_degree = 2\n'
        f'pressure_degree = 1\n\n[time]\n{time}\n\n[stabilization]\nmethod = "none"\n\n'
        f"[solver]\n{solver}\n\n[output]\nvtk_every = {vtk_every}\n"
        f"diagnostics = [{diagnostics}]\n"
    )


def run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, text: str) -> tuple[int, Path]:
    path = tmp_path / "run.toml"
    path.write_text(text, encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    status = RUN_SIMULATION.main([str(path)])
    case = re.search(r'case = "(\w+)"', text).group(1)
    return status, tmp_path / "results" / case / "run"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def pvd_times(path: Path) -> list[float]:
    return [float(value) for value in re.findall(r'timestep="([^"]+)"', path.read_text())]


@pytest.mark.parametrize("case", sorted(RUNS))
def test_short_run_writes_every_output_with_sensible_contents(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, case: str
) -> None:
    status, folder = run(tmp_path, monkeypatch, config_text(case))
    assert status == 0
    config = load_config(folder / "config_resolved.toml")
    vtk_every = config.output.vtk_every
    assert sorted(path.name for path in folder.iterdir()) == sorted(
        FILES + (["fields"] if vtk_every else [])
    )
    steps = round(config.time.t_end / config.time.dt)
    time_steps = rows(folder / "time_steps.csv")
    expected_schemes = [
        scheme_at_step(config.time.scheme, config.time.startup, step)
        for step in range(1, steps + 1)
    ]
    assert [row["scheme_used"] for row in time_steps] == expected_schemes
    assert all(row["converged"] == "True" and row["outer_iterations"] == "1" for row in time_steps)
    inner = rows(folder / "inner_iterations.csv")
    for step in range(1, steps + 1):
        iterations = [row for row in inner if row["step"] == str(step)]
        first, last = iterations[0], iterations[-1]
        assert (first["inner_iteration"], first["relative_residual"]) == ("0", "1.0")
        assert first["update_norm"] == "" and first["linear_iterations"] == ""
        tolerance = max(1e-12, 1e-10 * float(first["residual_norm"]))
        assert float(last["residual_norm"]) <= tolerance
        methods = [row["method"] for row in iterations[1:]]
        if config.solver.linearization == "picard":
            assert set(methods) == {"picard"}
        elif config.solver.picard_iterations_before_newton:
            assert methods[:2] == ["picard", "picard"] and set(methods[2:]) == {"newton"}
        else:
            assert set(methods) == {"newton"}
        assert time_steps[step - 1]["inner_iterations_total"] == str(len(iterations) - 1)
    outer = rows(folder / "outer_iterations.csv")
    assert [row["converged"] for row in outer] == ["True"] * steps
    diagnostics = rows(folder / "diagnostics.csv")
    assert list(diagnostics[0]) == ["step", "time", *config.output.diagnostics]
    assert len(diagnostics) == steps
    assert all(math.isfinite(float(row["velocity_error_l2"])) for row in diagnostics)
    assert json.loads((folder / "run_info.json").read_text())["exit_status"] == 0
    log = (folder / "run.log").read_text()
    assert f"step {steps}, t = {config.time.t_end:g}" in log
    if vtk_every:
        written = [step for step in range(1, steps + 1) if step % vtk_every == 0 or step == steps]
        assert pvd_times(folder / "fields" / "velocity.pvd") == pytest.approx(
            [step * config.time.dt for step in written]
        )


def test_crank_nicolson_pressure_is_written_at_the_half_level(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    status, folder = run(tmp_path, monkeypatch, config_text("manufactured_temporal"))
    assert status == 0
    assert pvd_times(folder / "fields" / "velocity.pvd") == pytest.approx([0.5, 0.75])
    assert pvd_times(folder / "fields" / "pressure.pvd") == pytest.approx([0.375, 0.625])


def test_inner_solve_failure_exits_nonzero_and_writes_the_last_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    text = config_text("kovasznay").replace(
        "picard_iterations_before_newton = 2", 'linearization = "picard"\ninner_max_iterations = 1'
    )
    status, folder = run(tmp_path, monkeypatch, text)
    assert status == 1
    assert pvd_times(folder / "fields" / "velocity.pvd") == [0.5]
    assert pvd_times(folder / "fields" / "pressure.pvd") == [0.5]
    assert rows(folder / "time_steps.csv")[-1]["converged"] == "False"
    assert rows(folder / "outer_iterations.csv")[-1]["converged"] == "False"
    assert json.loads((folder / "run_info.json").read_text())["exit_status"] == 1
    log = (folder / "run.log").read_text()
    assert "Step 1: the inner solve did not converge" in log
    assert "iteration 1 (picard)" in log


def steady_kovasznay_solve(
    problem: Problem, options: SolverOptions
) -> tuple[InnerSolver, fd.Function, InnerSolveResult]:
    space = mixed_space(problem.mesh, 2, 1)
    velocity_space = space.sub(0).collapse()
    dx = volume_measure(problem.mesh, 6)
    x = fd.SpatialCoordinate(problem.mesh)
    level = SolveLevel(
        fd.Constant(0.0),
        fd.Constant(1.0),
        fd.Function(velocity_space),
        problem.forcing(x, fd.Constant(0.0)),
    )
    solution, frozen = fd.Function(space), fd.Function(velocity_space)
    u, p = ufl.split(solution)
    v, q = fd.TestFunctions(space)
    nu = problem.viscosity
    newton = galerkin_residual(u, p, v, q, u, u, nu, level, dx, None)
    picard = galerkin_residual(u, p, v, q, frozen, frozen, nu, level, dx, None)
    bcs = [
        dirichlet_bc(space.sub(0), name, bc.value(x, fd.Constant(0.0)), problem.boundary_tags)
        for name, bc in problem.boundary_conditions.items()
    ]
    null_space = constant_pressure_null_space(space)
    solver = InnerSolver(solution, newton, picard, frozen, bcs, null_space, options)
    return solver, solution, solver.solve()


def kovasznay() -> Problem:
    return build_problem(Config(case="kovasznay", run_name="test", mesh=MeshOptions(cells=(6, 8))))


def test_newton_picard_and_warm_started_newton_reach_the_same_solution() -> None:
    problem = kovasznay()
    options = replace(SolverOptions(), inner_max_iterations=60)
    newton_solver, newton, newton_result = steady_kovasznay_solve(problem, options)
    picard_solver, picard, picard_result = steady_kovasznay_solve(
        problem, replace(options, linearization="picard")
    )
    _, warm, warm_result = steady_kovasznay_solve(
        problem, replace(options, picard_iterations_before_newton=3)
    )
    assert newton_result.converged and picard_result.converged and warm_result.converged
    assert {row.method for row in newton_result.iterations} == {"newton"}
    assert {row.method for row in picard_result.iterations} == {"picard"}
    assert [row.method for row in warm_result.iterations[1:4]] == ["picard"] * 3
    assert {row.method for row in warm_result.iterations[4:]} == {"newton"}
    assert len(newton_result.iterations) < len(picard_result.iterations)
    for solver, result in ((newton_solver, newton_result), (picard_solver, picard_result)):
        assert solver.residual_norm() == pytest.approx(
            result.iterations[-1].residual_norm, rel=1e-6, abs=1e-14
        )
    reference = newton.subfunctions[0]
    for other in (picard, warm):
        difference = fd.norm(other.subfunctions[0] - reference) / fd.norm(reference)
        assert difference < 1e-8


def test_newton_logs_updates_and_linear_iterations() -> None:
    _, _, result = steady_kovasznay_solve(kovasznay(), SolverOptions())
    assert result.converged
    later = result.iterations[1:]
    assert later and all(row.update_norm > 0 and row.linear_iterations == 1 for row in later)
    norms = [row.residual_norm for row in result.iterations]
    assert norms[-1] < 1e-10 * norms[0]
