"""Dynamic scalar closure, matched histories, solve-level coupling and temporal accuracy."""

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
from scipy.linalg import eigh

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.discretization.function_spaces import scalar_subscale_space
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.formulation.galerkin import SolveLevel
from navier_stokes_vms.formulation.stabilization_parameters import dynamic_parameter
from navier_stokes_vms.formulation.thermal import scalar_residual, scalar_subscale
from navier_stokes_vms.problem.boundary_conditions import Dirichlet
from navier_stokes_vms.problem.mesh_loading import periodic_rectangle, rectangle
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.problem.scalar_problem import PassiveScalarProblem
from navier_stokes_vms.problem.thermal_problem import ThermalProblem, check_thermal_supported
from navier_stokes_vms.solver.nonlinear_solver import StepSolveResult
from navier_stokes_vms.solver.subscale_storage import SubscaleStorage
from navier_stokes_vms.solver.time_stepper import TimeStepper

SCHEMES = ("bdf1", "bdf2", "bdf3", "crank_nicolson")


def configuration(
    mode: str = "oss", scheme: str = "bdf2", dt: float = 0.05, startup: str = "lower_order"
) -> Config:
    return Config(
        case="manufactured_boussinesq",
        run_name="b3",
        mesh=MeshOptions(cells=(3, 3)),
        discretization=DiscretizationOptions(2, 1)
        if mode == "none"
        else DiscretizationOptions(1, 1),
        time=TimeOptions(scheme=scheme, dt=dt, t_end=4 * dt, startup=startup),
        stabilization=(
            StabilizationOptions(method="none")
            if mode == "none"
            else StabilizationOptions(
                subscale_space=mode, subscale_dynamics="dynamic", scale_splitting="linear"
            )
        ),
        output=OutputOptions(vtk_every=0),
    )


def norm(expr: ufl.core.expr.Expr, dx: ufl.measure.Measure) -> float:
    return fd.assemble(ufl.inner(expr, expr) * dx) ** 0.5


@pytest.mark.parametrize("scheme", SCHEMES)
def test_scalar_subscale_ode_history_and_endpoint_recovery(scheme: str) -> None:
    mesh = fd.UnitSquareMesh(2, 2)
    Q = scalar_subscale_space(mesh, 3)
    storage = SubscaleStorage(Q, 3)
    field = fd.Function(fd.FunctionSpace(mesh, "CG", 1))
    alpha, dt = fd.Constant(1.0), fd.Constant(0.03)
    source = fd.Constant(1.0)
    level = SolveLevel(alpha, dt, fd.Function(field.function_space()), source)
    tau = fd.Constant(0.37)
    rate = dynamic_parameter(tau, alpha, dt)
    residual = scalar_residual(field, ufl.as_vector((0.0, 0.0)), 0.1, level, True, storage.history)
    subscale = scalar_subscale(residual, rate)
    previous = [0.0, 0.0, 0.0]
    for step in range(1, 6):
        coeff = step_coefficients(scheme, "lower_order", step)
        alpha.assign(coeff.alpha0)
        storage.form_history(coeff)
        history = sum(w * y for w, y in zip(coeff.history_weights, previous, strict=False))
        source.assign(0.8 + 0.1 * step)
        solve = (float(source) + history / 0.03) / (coeff.alpha0 / 0.03 + 1 / 0.37)
        expected = (
            coeff.extrapolation_weight_solve * solve + coeff.extrapolation_weight_old * previous[0]
        )
        assert norm(subscale - solve, fd.dx(domain=mesh, degree=3)) < 1e-12
        storage.store(subscale, coeff)
        assert norm(storage.latest - expected, fd.dx(domain=mesh, degree=3)) < 1e-12
        previous = [expected, *previous[:2]]


def test_forced_zero_subscale_history_limits_bdf3_startup() -> None:
    # For y' + y/τ = 1 and y(0)=0, a zero negative history is not an exact history.
    tau = 0.37
    mesh = fd.UnitSquareMesh(2, 2)
    field = fd.Function(fd.FunctionSpace(mesh, "CG", 1))
    history = fd.Function(scalar_subscale_space(mesh, 3))
    coeff = step_coefficients("bdf3", "exact_history", 1)
    errors = []
    for dt in (0.01, 0.005, 0.0025):
        level = SolveLevel(fd.Constant(coeff.alpha0), fd.Constant(dt), field, fd.Constant(1.0))
        residual = scalar_residual(field, ufl.as_vector((0.0, 0.0)), 0.1, level, True, history)
        subscale = scalar_subscale(
            residual, dynamic_parameter(fd.Constant(tau), level.alpha0, level.dt)
        )
        numerical = fd.assemble(subscale * fd.dx(domain=mesh, degree=3))
        exact = tau * (1 - np.exp(-dt / tau))
        errors.append(abs(numerical - exact))
    assert np.log2(errors[-2] / errors[-1]) == pytest.approx(1.0, abs=0.04)


@pytest.mark.parametrize("scheme", SCHEMES)
def test_dynamic_oss_histories_remain_orthogonal(scheme: str, tmp_path: Path) -> None:
    config = configuration(scheme=scheme)
    problem = build_problem(config)
    t = problem.thermal
    scalar = PassiveScalarProblem(
        "q",
        2 * t.diffusivity,
        t.source,
        t.initial_temperature,
        t.boundary_conditions,
        t.exact_temperature,
    )
    problem = replace(problem, passive_scalars=(scalar,))
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b3"))
    stepper._initialise()
    for step in range(1, 4):
        coefficients = step_coefficients(scheme, "lower_order", step)
        stepper._set_level(step, coefficients)
        assert stepper._step_solver.solve().converged
        for name in ("temperature", "q"):
            S = stepper._scalar_states[name].levels[0].function_space()
            weak = fd.assemble(stepper._scalar_subscales[name] * fd.TestFunction(S) * stepper._dx)
            with weak.dat.vec_ro as vector:
                assert vector.norm() < 1e-9
        stepper._advance(coefficients)
        for name, storage in stepper._scalar_subscale_storage.items():
            S = stepper.scalar(name).function_space()
            weak = fd.assemble(storage.latest * fd.TestFunction(S) * stepper._dx)
            with weak.dat.vec_ro as vector:
                assert vector.norm() < 1e-9
    assert norm(stepper._scalar_subscale_storage["q"].latest, stepper._dx) > 1e-8


@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize("scheme", SCHEMES)
def test_dynamic_scalar_copy_and_cn_boundaries(mode: str, scheme: str, tmp_path: Path) -> None:
    config = replace(configuration(mode, scheme), case_parameters={"expansion_coefficient": 0.0})
    problem = build_problem(config)
    t = problem.thermal
    q = PassiveScalarProblem(
        "q",
        t.diffusivity,
        t.source,
        t.initial_temperature,
        t.boundary_conditions,
        t.exact_temperature,
    )
    stepper = TimeStepper(
        config, replace(problem, passive_scalars=(q,)), tmp_path, logging.getLogger("b3")
    )
    stepper.run()
    assert norm(stepper.scalar("q") - stepper.temperature, stepper._dx) < 1e-9
    x = fd.SpatialCoordinate(problem.mesh)
    exact = t.exact_temperature(x, fd.Constant(config.time.t_end))
    for tag in problem.boundary_tags.values():
        assert fd.assemble((stepper.temperature - exact) ** 2 * fd.ds(tag)) < 1e-18
    if scheme == "crank_nicolson":
        assert float(stepper._pressure_time) == pytest.approx(
            config.time.t_end - config.time.dt / 2
        )
        assert norm(stepper._scalar_subscale_storage["q"].latest, stepper._dx) > 1e-8


def test_exact_history_initializes_all_scalars(tmp_path: Path) -> None:
    config = configuration(scheme="bdf3", startup="exact_history")
    problem = build_problem(config)
    t = problem.thermal
    q = PassiveScalarProblem(
        "q",
        t.diffusivity,
        t.source,
        t.initial_temperature,
        t.boundary_conditions,
        t.exact_temperature,
    )
    stepper = TimeStepper(
        config, replace(problem, passive_scalars=(q,)), tmp_path, logging.getLogger("b3")
    )
    stepper._initialise()
    x = fd.SpatialCoordinate(problem.mesh)
    for name, state in stepper._scalar_states.items():
        for index in (1, 2):
            expected = fd.Function(state.levels[index].function_space()).interpolate(
                t.exact_temperature(x, fd.Constant(-index * config.time.dt))
            )
            assert norm(state.levels[index] - expected, stepper._dx) < 1e-12
        for stored in stepper._scalar_subscale_storage[name]._levels:
            assert norm(stored, stepper._dx) == 0.0


def test_exact_history_requires_temperature_and_passive_data() -> None:
    config = configuration(scheme="bdf3", startup="exact_history")
    problem = build_problem(config)
    with pytest.raises(ValueError, match="temperature"):
        check_thermal_supported(
            config, replace(problem, thermal=replace(problem.thermal, exact_temperature=None))
        )
    with pytest.raises(ValueError, match="velocity"):
        check_thermal_supported(config, replace(problem, exact_velocity=None))
    t = problem.thermal
    q = PassiveScalarProblem(
        "q", t.diffusivity, t.source, t.initial_temperature, t.boundary_conditions
    )
    with pytest.raises(ValueError, match="q"):
        check_thermal_supported(config, replace(problem, passive_scalars=(q,)))


@pytest.mark.parametrize(
    "beta,gravity", [(0.0, (0.0, -1.0)), (0.7, (0.0, -1.0)), (0.7, (0.0, 0.0))]
)
def test_cn_buoyancy_power_uses_matching_half_level(
    tmp_path: Path, beta: float, gravity: tuple[float, float]
) -> None:
    config = replace(
        configuration(mode="none", scheme="crank_nicolson", dt=0.1),
        output=OutputOptions(vtk_every=1, diagnostics=("buoyancy_power",)),
    )
    mesh, tags = periodic_rectangle((4, 4), (0.0, 0.0), (1.0, 1.0))

    def zero(x: ufl.core.expr.Expr, _t: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    t = ThermalProblem(
        0.1, beta, gravity, 0.0, lambda x, _t: 1.0 + 0.0 * x[0], lambda x: 1.0 + 0.0 * x[0], {}
    )
    problem = Problem(mesh, {}, tags, 0.1, zero, lambda x: zero(x, 0.0), thermal=t)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b3"))
    stepper.run()
    with (tmp_path / "diagnostics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    beta = -beta * gravity[1]
    # BE first step: T1=1.1, w1=0.1*0.7*T1. Midpoint steps use the average temperature.
    old_velocity = 0.1 * beta * 1.1
    for step, row in enumerate(rows, 1):
        if step == 1:
            expected = beta * 1.1 * old_velocity
        else:
            midpoint_T = 1.0 + 0.1 * (step - 0.5)
            new_velocity = old_velocity + 0.1 * beta * midpoint_T
            expected = beta * midpoint_T * (old_velocity + new_velocity) / 2
            old_velocity = new_velocity
        assert float(row["buoyancy_power"]) == pytest.approx(expected, abs=1e-10)
    assert norm(stepper.temperature - fd.Constant(1.4), stepper._dx) < 1e-10


def diffusion_problem(degree: int) -> tuple[Problem, fd.Function, float]:
    mesh, tags = periodic_rectangle((4, 4), (0.0, 0.0), (1.0, 1.0))
    S = fd.FunctionSpace(mesh, "CG", degree)
    trial, test = fd.TrialFunction(S), fd.TestFunction(S)
    mass = fd.assemble(trial * test * fd.dx).petscmat.convert("dense").getDenseArray().copy()
    stiffness = (
        fd.assemble(0.1 * ufl.dot(ufl.grad(trial), ufl.grad(test)) * fd.dx)
        .petscmat.convert("dense")
        .getDenseArray()
        .copy()
    )
    values, vectors = eigh(stiffness, mass)
    index = np.flatnonzero(values > 1e-8)[0]
    rate = float(values[index])
    mode = fd.Function(S)
    with mode.dat.vec as vector:
        vector.array[:] = vectors[:, index]

    def zero(x: ufl.core.expr.Expr, _t: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    thermal = ThermalProblem(
        0.1,
        0.0,
        (0.0, -1.0),
        0.0,
        lambda x, _t: 0.0 * x[0],
        lambda _x: mode,
        {},
        lambda _x, t: ufl.exp(-rate * t) * mode,
    )
    problem = Problem(
        mesh, {}, tags, 0.1, zero, lambda x: zero(x, 0.0), exact_velocity=zero, thermal=thermal
    )
    return problem, mode, rate


@pytest.mark.parametrize("family", ["none", "static_oss", "dynamic_oss"])
@pytest.mark.parametrize(
    "scheme,order", [("bdf1", 1), ("bdf2", 2), ("bdf3", 3), ("crank_nicolson", 2)]
)
def test_source_free_discrete_eigenmode_temporal_order(
    family: str, scheme: str, order: int, tmp_path: Path
) -> None:
    # A semidiscrete diffusion eigenmode removes spatial error from this time-operator check.
    degree = 2 if family == "none" else 1
    problem, mode, rate = diffusion_problem(degree)
    errors = []
    for dt in (0.05, 0.025, 0.0125):
        config = configuration(
            "none" if family == "none" else "oss",
            scheme,
            dt,
            "exact_history" if scheme in ("bdf2", "bdf3") else "lower_order",
        )
        config = replace(config, time=replace(config.time, t_end=0.2))
        if family == "static_oss":
            config = replace(
                config, stabilization=replace(config.stabilization, subscale_dynamics="static")
            )
        folder = tmp_path / str(dt)
        folder.mkdir()
        stepper = TimeStepper(config, problem, folder, logging.getLogger("b3"))
        stepper.run()
        errors.append(norm(stepper.temperature - np.exp(-rate * 0.2) * mode, stepper._dx))
    observed = np.log2(errors[-2] / errors[-1])
    (tmp_path / "temporal_order.json").write_text(
        json.dumps(
            {
                "family": family,
                "scheme": scheme,
                "diffusion_eigenvalue": rate,
                "time_steps": [0.05, 0.025, 0.0125],
                "errors": errors,
                "observed_order": float(observed),
            },
            indent=2,
        )
        + "\n"
    )
    assert observed == pytest.approx(order, abs=0.25), (errors, observed)


@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize("scheme", SCHEMES)
def test_dynamic_static_steady_state_agreement(mode: str, scheme: str, tmp_path: Path) -> None:
    mesh, tags = rectangle((3, 3), (0.0, 0.0), (1.0, 1.0))

    def zero(x: ufl.core.expr.Expr, _t: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    thermal = ThermalProblem(
        0.1,
        0.0,
        (0.0, -1.0),
        0.0,
        lambda x, _t: 1.0 + 0.0 * x[0],
        lambda x: 0.0 * x[0],
        {n: Dirichlet(lambda x, _t: 0.0 * x[0]) for n in tags},
    )
    problem = Problem(
        mesh,
        {n: Dirichlet(zero) for n in tags},
        tags,
        0.1,
        zero,
        lambda x: zero(x, 0.0),
        thermal=thermal,
    )
    dynamic_config = configuration(mode, scheme)
    static_config = replace(
        dynamic_config,
        stabilization=replace(dynamic_config.stabilization, subscale_dynamics="static"),
    )
    static = TimeStepper(static_config, problem, tmp_path, logging.getLogger("b3"))
    static._initialise()
    static._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    static._alpha0.assign(0.0)
    assert static._step_solver.solve().converged
    reference = static._temperature.copy(deepcopy=True)
    storage_value = fd.Function(scalar_subscale_space(mesh, 3)).interpolate(
        static._scalar_subscales["temperature"]
    )
    if mode == "asgs":
        assert norm(storage_value, static._dx) > 1e-3
    dynamic = TimeStepper(dynamic_config, problem, tmp_path, logging.getLogger("b3"))
    dynamic._initialise()
    dynamic._solution.assign(static._solution)
    for level in dynamic._scalar_states["temperature"].levels:
        level.assign(reference)
    for level in dynamic._scalar_subscale_storage["temperature"]._levels:
        level.assign(storage_value)
    coeff = step_coefficients(scheme, "lower_order", 3)
    dynamic._set_level(3, coeff)
    assert dynamic._step_solver.solve().converged
    assert norm(dynamic._temperature - reference, dynamic._dx) < 1e-9
    assert norm(dynamic._scalar_subscales["temperature"] - storage_value, dynamic._dx) < 1e-9
    dynamic._advance(coeff)
    assert norm(dynamic.temperature - reference, dynamic._dx) < 1e-9
    assert (
        norm(dynamic._scalar_subscale_storage["temperature"].latest - storage_value, dynamic._dx)
        < 1e-9
    )


@pytest.mark.parametrize("mode", ["asgs", "oss"])
def test_failed_dynamic_step_preserves_all_endpoint_histories(mode: str, tmp_path: Path) -> None:
    config = configuration(mode, scheme="crank_nicolson")
    problem = build_problem(config)
    t = problem.thermal
    q = PassiveScalarProblem(
        "q",
        t.diffusivity,
        t.source,
        t.initial_temperature,
        t.boundary_conditions,
        t.exact_temperature,
    )
    problem = replace(problem, passive_scalars=(q,))
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b3"))
    original = stepper._step_solver.solve
    snapshots = []
    fields = []
    calls = 0

    def solve() -> StepSolveResult:
        nonlocal calls
        calls += 1
        if calls == 1:
            return original()
        fields.extend(stepper._velocity_states.levels)
        fields.extend(stepper._subscales._levels)
        for state in stepper._scalar_states.values():
            fields.extend(state.levels)
        for storage in stepper._scalar_subscale_storage.values():
            fields.extend(storage._levels)
        snapshots.extend(field.dat.data_ro.copy() for field in fields)
        stepper._temperature.assign(10.0)
        return StepSolveResult(False, (), "forced failure")

    stepper._step_solver.solve = solve
    with pytest.raises(ConvergenceError):
        stepper.run()
    assert calls == 2
    for field, expected in zip(fields, snapshots, strict=True):
        np.testing.assert_array_equal(field.dat.data_ro, expected)
    assert (tmp_path / "fields" / "q.pvd").exists()


def test_incomplete_thermal_energy_budgets_remain_disabled() -> None:
    config = replace(
        configuration(scheme="crank_nicolson"), output=OutputOptions(diagnostics=("energy_budget",))
    )
    with pytest.raises(ValueError, match="boussinesq_energy_budget"):
        build_problem(config)


def test_cn_scalar_source_and_boundary_use_endpoint_averages(tmp_path: Path) -> None:
    config = configuration(mode="none", scheme="crank_nicolson", dt=0.1)
    mesh, tags = rectangle((3, 3), (0.0, 0.0), (1.0, 1.0))

    def zero(x: ufl.core.expr.Expr, _t: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    thermal = ThermalProblem(
        0.1,
        0.0,
        (0.0, -1.0),
        0.0,
        lambda _x, t: t * t,
        lambda x: 0.0 * x[0],
        {name: Dirichlet(lambda _x, t: t * t) for name in tags},
    )
    problem = Problem(
        mesh,
        {name: Dirichlet(zero) for name in tags},
        tags,
        0.1,
        zero,
        lambda x: zero(x, 0.0),
        thermal=thermal,
    )
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b3"))
    stepper._initialise()
    stepper._set_level(2, step_coefficients("crank_nicolson", "lower_order", 2))
    average = (0.1**2 + 0.2**2) / 2
    residual = fd.assemble(stepper._step_solver._inner._newton_residual)
    S = stepper._temperature.function_space()
    expected = fd.assemble(-average * fd.TestFunction(S) * stepper._dx)
    with residual.subfunctions[2].dat.vec_ro as actual, expected.dat.vec_ro as target:
        difference = actual.copy()
        difference.axpy(-1.0, target)
        assert difference.norm() < 1e-12
    for bc in stepper._step_solver._inner._bcs:
        bc.apply(stepper._solution)
    assert fd.assemble((stepper._temperature - average) ** 2 * fd.ds(domain=mesh)) < 1e-20
