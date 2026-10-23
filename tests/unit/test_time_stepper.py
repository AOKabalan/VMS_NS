"""Galerkin runs end to end, iteration records, field output and failure handling (phase 4a)."""

import csv
import importlib.util
import json
import logging
import math
import re
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import pytest
import ufl
from firedrake.exceptions import ConvergenceError

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
    load_config,
)
from navier_stokes_vms.discretization.function_spaces import (
    constant_pressure_null_space,
    mixed_space,
)
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.discretization.time_schemes import scheme_at_step
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.problem.boundary_conditions import Outflow, dirichlet_bc
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.solver.nonlinear_solver import InnerSolver, InnerSolveResult
from navier_stokes_vms.solver.time_stepper import TimeStepper

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
    "initial_state.json",
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
    assert list(diagnostics[0]) == [
        "step",
        "time",
        "pressure_time",
        "subscale_time",
        *config.output.diagnostics,
    ]
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


def test_after_step_callback_can_stop_the_run_early(tmp_path: Path) -> None:
    config = Config(
        case="manufactured_spatial",
        run_name="early",
        mesh=MeshOptions(cells=(4, 4)),
        discretization=DiscretizationOptions(velocity_degree=2, pressure_degree=1),
        time=TimeOptions(scheme="bdf1", dt=0.1, t_end=1.0),
        stabilization=StabilizationOptions(method="none"),
        output=OutputOptions(vtk_every=5),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger(__name__))
    seen = []

    def stop_at_step_three(step: int, time: float) -> bool:
        seen.append((step, time))
        return step == 3

    stepper.run(stop_at_step_three)
    assert seen == [(1, 0.1), (2, 0.2), (3, 0.3)]
    assert [row["step"] for row in rows(tmp_path / "time_steps.csv")] == ["1", "2", "3"]
    # The state at the stop is written even though step 3 is not a VTK step.
    assert pvd_times(tmp_path / "fields" / "velocity.pvd") == pytest.approx([0.3])


def test_crank_nicolson_pressure_is_written_at_the_half_level(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    status, folder = run(tmp_path, monkeypatch, config_text("manufactured_temporal"))
    assert status == 0
    assert pvd_times(folder / "fields" / "velocity.pvd") == pytest.approx([0.5, 0.75])
    assert pvd_times(folder / "fields" / "pressure.pvd") == pytest.approx([0.375, 0.625])
    diagnostics = rows(folder / "diagnostics.csv")
    assert [float(row["time"]) for row in diagnostics] == [0.25, 0.5, 0.75]
    assert [float(row["pressure_time"]) for row in diagnostics] == [0.25, 0.375, 0.625]
    assert all(row["subscale_time"] == "" for row in diagnostics)


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


@pytest.mark.parametrize(
    "solver",
    ['linearization = "picard"', "picard_iterations_before_newton = 2"],
    ids=["picard", "newton_after_picard_warm_up"],
)
def test_a_failed_picard_linear_solve_fails_the_step_at_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, solver: str
) -> None:
    text = config_text("kovasznay").replace("picard_iterations_before_newton = 2", solver)
    solve = fd.NonlinearVariationalSolver.solve
    attempts = []

    def failing_linear_solve(self: fd.NonlinearVariationalSolver) -> None:
        if self.snes.getType() != "ksponly":
            return solve(self)
        attempts.append(1)
        if len(attempts) > 3:
            raise RuntimeError("the failed linear solve was retried")
        raise ConvergenceError("synthetic linear solver failure")

    monkeypatch.setattr(fd.NonlinearVariationalSolver, "solve", failing_linear_solve)
    status, folder = run(tmp_path, monkeypatch, text)
    assert status == 1
    assert attempts == [1]
    assert pvd_times(folder / "fields" / "velocity.pvd") == [0.5]
    assert rows(folder / "time_steps.csv")[-1]["converged"] == "False"
    assert [row["inner_iteration"] for row in rows(folder / "inner_iterations.csv")] == ["0"]
    log = (folder / "run.log").read_text(encoding="utf-8")
    assert "Step 1: Picard iteration 1: synthetic linear solver failure" in log


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


def test_the_picard_warm_up_counts_toward_the_inner_iteration_limit() -> None:
    options = replace(SolverOptions(), inner_max_iterations=2, picard_iterations_before_newton=5)
    _, _, result = steady_kovasznay_solve(kovasznay(), options)
    assert not result.converged
    assert [row.method for row in result.iterations[1:]] == ["picard", "picard"]


@pytest.mark.parametrize(
    ("degrees", "stabilization"),
    [
        ((2, 1), StabilizationOptions(method="none")),
        (
            (1, 1),
            StabilizationOptions(
                subscale_space="asgs", subscale_dynamics="static", scale_splitting="linear"
            ),
        ),
    ],
    ids=["galerkin", "asgs"],
)
@pytest.mark.parametrize("outflow", [False, True])
def test_pressure_has_zero_spatial_mean_exactly_with_the_null_space(
    tmp_path: Path,
    degrees: tuple[int, int],
    stabilization: StabilizationOptions,
    outflow: bool,
) -> None:
    config = Config(
        case="manufactured_spatial",
        run_name="mean",
        mesh=MeshOptions(cells=(4, 4)),
        discretization=DiscretizationOptions(*degrees),
        time=TimeOptions(scheme="bdf1", dt=0.05, t_end=0.1),
        stabilization=stabilization,
        output=OutputOptions(vtk_every=0),
    )
    problem = build_problem(config)
    if outflow:
        conditions = {**problem.boundary_conditions, "x_max": Outflow()}
        problem = replace(problem, boundary_conditions=conditions)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger(__name__))
    stepper.run()
    pressure = stepper.pressure
    mean = fd.assemble(pressure * stepper._dx) / fd.assemble(fd.Constant(1.0) * stepper._dx)
    size = fd.norm(pressure)
    if outflow:
        assert abs(mean) > 1e-6 * size
    else:
        assert abs(mean) <= 1e-13 * size


def test_newton_logs_updates_and_linear_iterations() -> None:
    _, _, result = steady_kovasznay_solve(kovasznay(), SolverOptions())
    assert result.converged
    later = result.iterations[1:]
    assert later and all(row.update_norm > 0 and row.linear_iterations == 1 for row in later)
    norms = [row.residual_norm for row in result.iterations]
    assert norms[-1] < 1e-10 * norms[0]


def nonlinear_asgs_text(
    dynamics: str, solver: str = "", space: str = "asgs", splitting: str = "nonlinear"
) -> str:
    return (
        'case = "kovasznay"\n\n[mesh]\ncells = [12, 16]\n\n[discretization]\n'
        'velocity_degree = 1\npressure_degree = 1\n\n[time]\nscheme = "bdf2"\ndt = 0.5\n'
        f't_end = 1.0\n\n[stabilization]\nmethod = "vms"\nsubscale_space = "{space}"\n'
        f'subscale_dynamics = "{dynamics}"\nscale_splitting = "{splitting}"\n\n'
        f'[solver]\n{solver}\n\n[output]\nvtk_every = 0\ndiagnostics = ["velocity_error_l2"]\n'
    )


@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
def test_nonlinear_vms_records_every_outer_iteration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dynamics: str, space: str
) -> None:
    status, folder = run(tmp_path, monkeypatch, nonlinear_asgs_text(dynamics, space=space))
    assert status == 0
    outer, inner = rows(folder / "outer_iterations.csv"), rows(folder / "inner_iterations.csv")
    checked = ["velocity_change", "pressure_change", "subscale_change"]
    if space == "oss":
        checked += ["xi_change", "zeta_change"]
    for step_row in rows(folder / "time_steps.csv"):
        step = step_row["step"]
        iterations = [row for row in outer if row["step"] == step]
        assert len(iterations) >= 2
        assert step_row["outer_iterations"] == str(len(iterations))
        assert [row["converged"] for row in iterations] == ["False"] * (len(iterations) - 1) + [
            "True"
        ]
        assert all(row["subscale_change"] != "" for row in iterations)
        assert all((row["xi_change"] != "") == (space == "oss") for row in iterations)
        last = iterations[-1]
        for name in checked:
            assert float(last[name]) <= 1e-8
        steps_inner = [row for row in inner if row["step"] == step]
        assert {row["outer_iteration"] for row in steps_inner} == {
            row["outer_iteration"] for row in iterations
        }
        assert int(step_row["inner_iterations_total"]) == sum(
            int(row["inner_iterations"]) for row in iterations
        )


COARSE_MESH_HINT = (
    "static nonlinear subscales on a coarse mesh can make (I + τ∇u) nearly singular; "
    "use dynamic subscales or a finer mesh"
)
LOCAL_FAILURE = "local subscale update did not converge in 1 sweeps"


@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
def test_linear_oss_solves_each_step_in_one_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dynamics: str
) -> None:
    text = nonlinear_asgs_text(dynamics, space="oss", splitting="linear")
    status, folder = run(tmp_path, monkeypatch, text)
    assert status == 0
    outer = rows(folder / "outer_iterations.csv")
    assert [row["outer_iteration"] for row in outer] == ["1", "1"]
    assert all(row["converged"] == "True" for row in outer)
    assert all(row["subscale_change"] == "" for row in outer)
    assert all(float(row["xi_change"]) > 0 and float(row["zeta_change"]) > 0 for row in outer)


@pytest.mark.parametrize(
    ("dynamics", "setting", "reason"),
    [
        ("static", "outer_max_iterations = 1", "outer iteration did not converge in 1 iterations"),
        ("static", "subscale_max_iterations = 1", f"{LOCAL_FAILURE}; {COARSE_MESH_HINT}"),
        ("dynamic", "subscale_max_iterations = 1", f"{LOCAL_FAILURE}. Iteration history"),
    ],
    ids=["outer_limit", "static_local_limit", "dynamic_local_limit"],
)
def test_nonlinear_asgs_failures_write_the_last_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dynamics: str, setting: str, reason: str
) -> None:
    status, folder = run(tmp_path, monkeypatch, nonlinear_asgs_text(dynamics, setting))
    assert status == 1
    assert pvd_times(folder / "fields" / "velocity.pvd") == [0.5]
    assert rows(folder / "time_steps.csv")[-1]["converged"] == "False"
    log = (folder / "run.log").read_text(encoding="utf-8")
    assert f"Step 1: {reason}" in log
    assert "outer iteration 1:" in log
