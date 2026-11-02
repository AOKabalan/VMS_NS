"""B1 coupled physics, boundary signs, Jacobian coupling and thermal output."""

import csv
import json
import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl
from firedrake.exceptions import ConvergenceError
from ufl.core.expr import Expr

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    ThermalOptions,
    TimeOptions,
    config_as_dict,
    config_to_toml,
    load_config,
)
from navier_stokes_vms.discretization.function_spaces import FieldLayout, mixed_space
from navier_stokes_vms.formulation.galerkin import SolveLevel
from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, Dirichlet
from navier_stokes_vms.problem.mesh_loading import periodic_rectangle, rectangle
from navier_stokes_vms.problem.problem_definition import Problem, check_problem
from navier_stokes_vms.problem.thermal_problem import ThermalProblem, check_thermal_supported
from navier_stokes_vms.solver.nonlinear_solver import StepSolveResult
from navier_stokes_vms.solver.system_builder import galerkin_thermal_system
from navier_stokes_vms.solver.time_stepper import TimeStepper


def configuration(**changes: object) -> Config:
    return replace(
        Config(
            case="manufactured_boussinesq",
            run_name="b1",
            mesh=MeshOptions(cells=(4, 4)),
            discretization=DiscretizationOptions(2, 1),
            stabilization=StabilizationOptions(method="none"),
            time=TimeOptions(scheme="bdf1", dt=0.01, t_end=0.02),
            output=OutputOptions(vtk_every=0),
        ),
        **changes,
    )


def simple_problem(*, periodic: bool = False, beta: float = 0.0) -> Problem:
    mesh, tags = (periodic_rectangle if periodic else rectangle)((4, 4), (0.0, 0.0), (1.0, 1.0))

    def zero(x: Expr, _t: Expr) -> Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    thermal = ThermalProblem(
        0.1,
        beta,
        (0.0, -1.0),
        0.0,
        lambda x, _t: 0.0 * x[0],
        lambda x: 1.0 + 0.0 * x[0],
        {n: DiffusiveFlux() for n in tags},
    )
    return Problem(
        mesh,
        {n: Dirichlet(zero) for n in tags},
        tags,
        0.1,
        zero,
        lambda x: zero(x, 0.0),
        thermal=thermal,
    )


def run(config: Config, problem: Problem, folder: Path) -> TimeStepper:
    folder.mkdir()
    stepper = TimeStepper(config, problem, folder, logging.getLogger("b1"))
    stepper.run()
    return stepper


def norm(expression: ufl.core.expr.Expr) -> float:
    return fd.assemble(ufl.inner(expression, expression) * fd.dx) ** 0.5


def test_diffusion_mode_decay(tmp_path: Path) -> None:
    problem = simple_problem()

    def mode(x: Expr) -> Expr:
        return ufl.cos(ufl.pi * x[0]) * ufl.cos(ufl.pi * x[1])

    problem = replace(problem, thermal=replace(problem.thermal, initial_temperature=mode))
    stepper = run(configuration(), problem, tmp_path / "diffusion")
    x = fd.SpatialCoordinate(problem.mesh)
    initial = fd.Function(stepper.temperature.function_space()).interpolate(mode(x))
    amplitude = fd.assemble(stepper.temperature * mode(x) * fd.dx) / fd.assemble(
        initial * mode(x) * fd.dx
    )
    expected = (1 + 0.01 * 0.1 * 2 * np.pi**2) ** -2
    assert amplitude == pytest.approx(expected, abs=2e-4)
    assert norm(stepper.velocity) < 1e-12


def test_nonzero_flux_and_source_signs(tmp_path: Path) -> None:
    problem = simple_problem()

    def field(x: Expr, t: Expr) -> Expr:
        return x[0] ** 2 + t

    thermal = replace(
        problem.thermal,
        initial_temperature=lambda x: field(x, 0.0),
        source=lambda _x, _t: ufl.as_ufl(0.8),
        exact_temperature=field,
        boundary_conditions={
            n: DiffusiveFlux((lambda _x, _t: ufl.as_ufl(-0.2)) if n == "x_max" else None)
            for n in problem.boundary_tags
        },
    )
    stepper = run(configuration(), replace(problem, thermal=thermal), tmp_path / "flux")
    x = fd.SpatialCoordinate(problem.mesh)
    assert norm(stepper.temperature - field(x, 0.02)) < 1e-11
    assert fd.assemble(stepper.temperature * fd.dx) == pytest.approx(1 / 3 + 0.02, abs=1e-11)


def test_hydrostatic_vertical_conduction(tmp_path: Path) -> None:
    problem = simple_problem(beta=0.7)
    thermal = replace(
        problem.thermal,
        initial_temperature=lambda x: x[1],
        boundary_conditions={n: Dirichlet(lambda x, _t: x[1]) for n in problem.boundary_tags},
    )
    config = configuration(discretization=DiscretizationOptions(3, 2))
    stepper = run(config, replace(problem, thermal=thermal), tmp_path / "hydrostatic")
    x = fd.SpatialCoordinate(problem.mesh)
    assert norm(stepper.velocity) < 1e-10
    assert norm(stepper.pressure - 0.7 * (x[1] ** 2 / 2 - 1 / 6)) < 1e-10
    assert norm(stepper.temperature - x[1]) < 1e-11


def test_hot_periodic_fluid_accelerates_upward(tmp_path: Path) -> None:
    problem = simple_problem(periodic=True, beta=0.7)
    stepper = run(configuration(), problem, tmp_path / "hot")
    assert norm(stepper.velocity - ufl.as_vector((0.0, 0.014))) < 1e-11
    assert norm(stepper.temperature - fd.Constant(1.0)) < 1e-11


def test_zero_expansion_matches_isothermal_flow(tmp_path: Path) -> None:
    problem = simple_problem(periodic=True)

    def initial(x: Expr) -> Expr:
        return ufl.as_vector((1.0 + 0.0 * x[0], 0.3 + 0.0 * x[1]))

    problem = replace(
        problem,
        initial_velocity=initial,
        thermal=replace(problem.thermal, initial_temperature=lambda x: ufl.sin(2 * ufl.pi * x[0])),
    )
    config = configuration()
    thermal = run(config, problem, tmp_path / "passive")
    flow = run(config, replace(problem, thermal=None), tmp_path / "flow")
    assert norm(thermal.velocity - flow.velocity) < 1e-11
    assert norm(thermal.pressure - flow.pressure) < 1e-11
    assert thermal.temperature is not None and flow.temperature is None
    assert not (tmp_path / "flow" / "thermal_metadata.json").exists()


@pytest.mark.parametrize("convection", ["skew_symmetric", "nonconservative"])
def test_manufactured_solution_and_output(tmp_path: Path, convection: str) -> None:
    config = configuration(
        thermal=ThermalOptions(
            convection=convection,
            momentum_residual_scale=100.0,
            continuity_residual_scale=0.1,
            temperature_residual_scale=0.01,
        ),
        output=OutputOptions(vtk_every=1, diagnostics=("temperature_error_l2", "thermal_content")),
    )
    problem = build_problem(config)
    stepper = run(config, problem, tmp_path / "coupled")
    x = fd.SpatialCoordinate(problem.mesh)
    assert (
        norm(stepper.temperature - problem.thermal.exact_temperature(x, fd.Constant(0.02))) < 3e-5
    )
    with (tmp_path / "coupled" / "outer_iterations.csv").open() as f:
        assert all("temperature_change" in row for row in csv.DictReader(f))
    with (tmp_path / "coupled" / "inner_iterations.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    last_step = rows[-1]["step"]
    step_rows = [row for row in rows if row["step"] == last_step]
    for name in ("velocity", "pressure", "temperature"):
        initial = float(step_rows[0][f"{name}_residual"])
        final = float(step_rows[-1][f"{name}_residual"])
        assert final <= max(config.solver.inner_atol, config.solver.inner_rtol * initial)
    metadata = json.loads((tmp_path / "coupled" / "thermal_metadata.json").read_text())
    assert metadata["scaling"] == "nondimensional"
    assert metadata["gravity"] == [0.0, -1.0] and metadata["temperature_degree"] == 2
    assert (tmp_path / "coupled" / "fields" / "temperature.pvd").exists()


def test_jacobian_contains_both_coupling_directions() -> None:
    problem = simple_problem(beta=0.7)
    space = mixed_space(problem.mesh, 2, 1, 2)
    layout = FieldLayout(("velocity", "pressure", "temperature"))
    solution = fd.Function(space)
    x = fd.SpatialCoordinate(problem.mesh)
    solution.subfunctions[0].interpolate(ufl.as_vector((x[0], x[1])))
    solution.subfunctions[2].interpolate(1 + x[0] + 2 * x[1])
    level = SolveLevel(
        fd.Constant(1.0),
        fd.Constant(0.1),
        fd.Function(space.sub(0).collapse()),
        ufl.as_vector((0.0, 0.0)),
    )
    scalar = replace(level, history=fd.Function(space.sub(2).collapse()), forcing=ufl.as_ufl(0.0))
    residual = galerkin_thermal_system(
        solution,
        layout,
        0.1,
        level,
        scalar,
        problem.thermal,
        fd.dx(domain=problem.mesh, degree=6),
        "skew_symmetric_1",
        "skew_symmetric",
    )
    jacobian = fd.derivative(residual, solution)
    for source, target in ((2, 0), (0, 2)):
        direction = fd.Function(space)
        direction.subfunctions[source].assign(1.0)
        action = fd.assemble(ufl.action(jacobian, direction))
        with action.subfunctions[target].dat.vec_ro as v:
            assert v.norm() > 1e-3
        before = fd.assemble(residual)
        epsilon = 1e-6
        solution.assign(solution + epsilon * direction)
        after = fd.assemble(residual)
        solution.assign(solution - epsilon * direction)
        with before.dat.vec_ro as b, after.dat.vec as a, action.dat.vec_ro as j:
            a.axpy(-1.0, b)
            a.axpy(-epsilon, j)
            assert a.norm() < 1e-10


def test_failed_step_keeps_temperature_history_and_writes_snapshot(tmp_path: Path) -> None:
    problem = simple_problem()
    folder = tmp_path / "failure"
    folder.mkdir()
    stepper = TimeStepper(configuration(), problem, folder, logging.getLogger("b1"))
    stepper._step_solver.solve = lambda: StepSolveResult(False, (), "forced failure")
    with pytest.raises(ConvergenceError):
        stepper.run()
    assert norm(stepper.temperature - fd.Constant(1.0)) < 1e-12
    assert (folder / "fields" / "temperature.pvd").exists()


def test_thermal_config_roundtrip_and_inactive_omission(tmp_path: Path) -> None:
    assert "thermal" not in config_as_dict(configuration())
    config = configuration(thermal=ThermalOptions())
    path = tmp_path / "config.toml"
    path.write_text(config_to_toml(config))
    loaded = load_config(path)
    assert loaded.thermal.temperature_degree == 2
    build_problem(loaded)
    with pytest.raises(ValueError, match="requires a case"):
        check_thermal_supported(config, replace(simple_problem(), thermal=None))


@pytest.mark.parametrize(
    "changes",
    [
        dict(stabilization=StabilizationOptions()),
        dict(time=TimeOptions(scheme="bdf2")),
        dict(solver=SolverOptions(linearization="picard")),
        dict(thermal=ThermalOptions(temperature_degree=1)),
    ],
)
def test_unsupported_thermal_combinations_fail(changes: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        check_thermal_supported(configuration(**changes), simple_problem())


@pytest.mark.parametrize(
    "changes",
    [
        dict(diffusivity=0.0),
        dict(expansion_coefficient=-1.0),
        dict(gravity=(0.0,)),
        dict(reference_temperature=float("inf")),
        dict(boundary_conditions={}),
    ],
)
def test_invalid_physical_data_fail(changes: dict[str, object]) -> None:
    problem = simple_problem()
    with pytest.raises(ValueError):
        check_problem(replace(problem, thermal=replace(problem.thermal, **changes)))


def test_open_normal_velocity_is_rejected(tmp_path: Path) -> None:
    problem = simple_problem()
    problem = replace(
        problem,
        boundary_conditions={
            n: Dirichlet(lambda x, _t: ufl.as_vector((1.0 + 0.0 * x[0], 0.0)))
            for n in problem.boundary_tags
        },
    )
    with pytest.raises(ValueError, match="zero normal velocity"):
        run(configuration(), problem, tmp_path / "open")
