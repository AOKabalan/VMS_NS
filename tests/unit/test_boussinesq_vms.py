"""Static linear thermal subscales, coupled OSS projections and passive-scalar transport."""

import csv
import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.formulation.boussinesq import buoyancy
from navier_stokes_vms.formulation.galerkin import SolveLevel
from navier_stokes_vms.formulation.stabilization_parameters import transport_parameter
from navier_stokes_vms.formulation.thermal import scalar_residual, scalar_subscale, scalar_vms_terms
from navier_stokes_vms.problem.problem_definition import check_problem
from navier_stokes_vms.problem.scalar_problem import PassiveScalarProblem, transported_scalars
from navier_stokes_vms.solver.system_builder import thermal_system
from navier_stokes_vms.solver.time_stepper import TimeStepper


def configuration(mode: str = "oss", degree: int = 1, pressure_projection: bool = True) -> Config:
    return Config(
        case="manufactured_boussinesq",
        run_name="b2",
        mesh=MeshOptions(cells=(3, 3)),
        discretization=DiscretizationOptions(degree, degree),
        time=TimeOptions(scheme="bdf1", dt=0.01, t_end=0.02),
        stabilization=StabilizationOptions(
            subscale_space=mode,
            subscale_dynamics="static",
            scale_splitting="linear",
            pressure_subscale=pressure_projection,
        ),
        output=OutputOptions(vtk_every=0),
    )


def dual_norm(form: ufl.form.Form) -> float:
    with fd.assemble(form).dat.vec_ro as v:
        return v.norm()


@pytest.mark.parametrize("degree", [1, 2])
def test_scalar_stabilization_sign_and_viscous_terms(degree: int) -> None:
    mesh = fd.UnitSquareMesh(2, 2)
    S = fd.FunctionSpace(mesh, "CG", degree)
    x, y = fd.SpatialCoordinate(mesh)
    field = fd.Function(S).interpolate(x + 2 * y if degree == 1 else x * x + 2 * y * y)
    history = fd.Function(S).interpolate(0.1 * x - 0.3 * y)
    test = fd.TestFunction(S)
    a = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        ufl.as_vector((1 + y, 1 + x))
    )
    dx = volume_measure(mesh, 3 * degree)
    kappa = 0.2
    level = SolveLevel(fd.Constant(1.0), fd.Constant(0.1), history, 3 + x)
    tau = transport_parameter(a, kappa, degree, 4.0, 2.0)
    residual = scalar_residual(field, a, kappa, level, True)
    convection = ufl.dot(a, ufl.grad(field))
    laplacian = 0.0 if degree == 1 else 6.0
    independent = 3 + x - (field - history) / 0.1 - convection + kappa * laplacian
    assert fd.assemble((residual - independent) ** 2 * dx) < 1e-20
    sub = scalar_subscale(residual, tau)
    expected = (
        tau
        * ((field - history) / 0.1 + convection - kappa * laplacian - 3 - x)
        * (ufl.dot(a, ufl.grad(test)) + kappa * ufl.div(ufl.grad(test)))
        * dx
    )
    actual = scalar_vms_terms(sub, test, a, kappa, dx)
    assert dual_norm(actual - expected) < 1e-11
    if degree == 2:
        polynomial = fd.Function(S).interpolate(x * x + y * y)
        value = fd.assemble(scalar_vms_terms(sub, polynomial, a, kappa, dx))
        target = fd.assemble(-sub * (ufl.dot(a, ufl.as_vector((2 * x, 2 * y))) + 4 * kappa) * dx)
        assert value == pytest.approx(target, abs=1e-11)


@pytest.mark.parametrize("degree", [1, 2])
def test_pointwise_weighted_projection_cancels_fe_buoyancy(degree: int) -> None:
    mesh = fd.UnitSquareMesh(3, 3)
    S = fd.FunctionSpace(mesh, "CG", degree)
    V = fd.VectorFunctionSpace(mesh, "CG", degree)
    x, y = fd.SpatialCoordinate(mesh)
    a = fd.Function(V).interpolate(ufl.as_vector((1 + x, 1 + y)))
    T = fd.Function(S).interpolate(x + y if degree == 1 else x * y)
    tau = transport_parameter(a, 0.1, degree, 4.0, 2.0)
    dx = volume_measure(mesh, 3 * degree)
    force = buoyancy(T, 0.7, (0.0, -1.0), 0.2)
    trial, test = fd.TrialFunction(V), fd.TestFunction(V)
    xi = fd.Function(V)
    fd.solve(
        tau * ufl.inner(trial, test) * dx == tau * ufl.inner(force, test) * dx,
        xi,
        solver_parameters={"ksp_type": "preonly", "pc_type": "lu"},
    )
    assert fd.assemble(ufl.inner(tau * (force - xi), tau * (force - xi)) * dx) < 1e-22
    source = 1 + x * x + y**3
    eta = fd.Function(S)
    t, s = fd.TrialFunction(S), fd.TestFunction(S)
    fd.solve(
        tau * t * s * dx == tau * source * s * dx,
        eta,
        solver_parameters={"ksp_type": "preonly", "pc_type": "lu"},
    )
    sub = tau * (source - eta)
    assert fd.assemble(sub**2 * dx) > 1e-10
    assert dual_norm(0.7 * sub * ufl.dot(ufl.as_vector((0.0, -1.0)), test) * dx) < 1e-12
    assert dual_norm(0.7 * tau * source * ufl.dot(ufl.as_vector((0.0, -1.0)), test) * dx) > 1e-5


@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("pressure_projection", [False, True])
def test_coupled_oss_subscales_are_orthogonal(
    tmp_path: Path, degree: int, pressure_projection: bool
) -> None:
    config = configuration(degree=degree, pressure_projection=pressure_projection)
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger("b2"))
    # Capture orthogonality before the successful history shift changes the residual expression.
    stepper._initialise()
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    assert stepper._step_solver.solve().converged
    mesh = stepper._problem.mesh
    V = fd.VectorFunctionSpace(mesh, "CG", degree)
    S = fd.FunctionSpace(mesh, "CG", degree)
    assert dual_norm(ufl.inner(stepper._u_sub, fd.TestFunction(V)) * stepper._dx) < 1e-9
    sub = stepper._scalar_subscales["temperature"]
    assert dual_norm(sub * fd.TestFunction(S) * stepper._dx) < 1e-9
    assert (
        dual_norm(sub * ufl.dot(ufl.as_vector((0.0, -1.0)), fd.TestFunction(V)) * stepper._dx)
        < 1e-9
    )
    assert ("zeta" in stepper._layout.names) == pressure_projection
    assert "eta_temperature" in stepper._layout.names
    bcs = stepper._step_solver._inner._bcs
    physical = {stepper._layout.index("velocity"), stepper._layout.index("temperature")}
    assert all(bc.function_space().index in physical for bc in bcs)


def test_asgs_temperature_subscale_buoyancy_sign(tmp_path: Path) -> None:
    config = configuration(mode="asgs")
    problem = build_problem(config)

    def zero(x: ufl.core.expr.Expr, _t: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    problem = replace(
        problem,
        forcing=zero,
        initial_velocity=lambda x: zero(x, 0.0),
        thermal=replace(
            problem.thermal,
            source=lambda x, _t: 2.0 + 0.0 * x[0],
            reference_temperature=0.0,
            initial_temperature=lambda x: 0.0 * x[0],
        ),
    )
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b2"))
    # Zero state leaves no resolved force, momentum residual or mechanical stabilization.
    F = stepper._step_solver._inner._newton_residual
    residual = fd.assemble(F)
    V = stepper._solution.subfunctions[0].function_space()
    expected = fd.assemble(
        stepper._scalar_subscales["temperature"]
        * ufl.dot(ufl.as_vector((0.0, -1.0)), fd.TestFunction(V))
        * stepper._dx
    )
    with residual.subfunctions[0].dat.vec_ro as actual, expected.dat.vec_ro as target:
        difference = actual.copy()
        difference.axpy(-1.0, target)
        assert target.norm() > 1e-4
        assert difference.norm() < 1e-12


@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("newton_tau", ["frozen", "exact"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
def test_thermal_jacobian_with_passive_scalar(
    tmp_path: Path, mode: str, degree: int, newton_tau: str, dynamics: str
) -> None:
    config = replace(configuration(mode, degree), solver=SolverOptions(newton_tau=newton_tau))
    config = replace(
        config, stabilization=replace(config.stabilization, subscale_dynamics=dynamics)
    )
    problem = build_problem(config)
    thermal = problem.thermal
    passive = PassiveScalarProblem(
        "q",
        thermal.diffusivity,
        thermal.source,
        thermal.initial_temperature,
        thermal.boundary_conditions,
    )
    problem = replace(problem, passive_scalars=(passive,))
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b2"))
    stepper._initialise()
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    solution = stepper._solution
    direction = fd.Function(solution.function_space())
    rng = np.random.default_rng(37)
    for field in (*solution.subfunctions, *direction.subfunctions):
        field.dat.data[:] = rng.uniform(0.2, 1.0, field.dat.data.shape)
    if dynamics == "dynamic":
        histories = [
            stepper._subscales.history,
            *(storage.history for storage in stepper._scalar_subscale_storage.values()),
        ]
        for history in histories:
            history.dat.data[:] = rng.uniform(-0.5, 0.5, history.dat.data.shape)
    base = solution.copy(deepcopy=True)
    x = fd.SpatialCoordinate(problem.mesh)
    level = SolveLevel(
        stepper._alpha0,
        fd.Constant(config.time.dt),
        stepper._velocity_states.history,
        problem.forcing(x, stepper._t_new),
    )
    levels = {
        s.name: replace(
            level,
            history=stepper._scalar_states[s.name].history,
            forcing=s.source(x, stepper._t_new),
        )
        for s in transported_scalars(problem)
    }
    reference = thermal_system(
        solution,
        stepper._layout,
        problem,
        level,
        levels,
        stepper._dx,
        config,
        base.subfunctions[0] if newton_tau == "frozen" else None,
        velocity_history=None if dynamics == "static" else stepper._subscales.history,
        scalar_histories={
            name: storage.history for name, storage in stepper._scalar_subscale_storage.items()
        },
    ).residual
    J = stepper._step_solver._inner._newton._problem.J
    action = fd.assemble(ufl.action(J, direction))
    epsilon = 1e-6
    solution.assign(base + epsilon * direction)
    plus = fd.assemble(reference)
    solution.assign(base - epsilon * direction)
    minus = fd.assemble(reference)
    with plus.dat.vec_ro as vp, minus.dat.vec_ro as vm, action.dat.vec_ro as derivative:
        difference = vp.copy()
        difference.axpy(-1.0, vm)
        difference.scale(1 / (2 * epsilon))
        difference.axpy(-1.0, derivative)
        assert difference.norm() < 1e-7 * derivative.norm()


@pytest.mark.parametrize("mode", ["none", "asgs", "oss"])
def test_passive_copy_matches_temperature_and_has_named_output(tmp_path: Path, mode: str) -> None:
    config = configuration("oss" if mode == "none" else mode)
    if mode == "none":
        config = replace(
            config,
            discretization=DiscretizationOptions(2, 1),
            stabilization=StabilizationOptions(method="none"),
        )
    config = replace(
        config,
        case_parameters={"expansion_coefficient": 0.0},
        output=OutputOptions(
            vtk_every=1, diagnostics=("passive_scalar_content", "passive_scalar_l2")
        ),
    )
    problem = build_problem(config)
    thermal = problem.thermal
    scalar = PassiveScalarProblem(
        "q",
        thermal.diffusivity,
        thermal.source,
        thermal.initial_temperature,
        thermal.boundary_conditions,
    )
    problem = replace(problem, passive_scalars=(scalar,))
    check_problem(problem)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b2"))
    stepper.run()
    assert fd.assemble((stepper.scalar("q") - stepper.temperature) ** 2 * stepper._dx) < 1e-18
    assert (tmp_path / "fields" / "q.pvd").exists()
    with (tmp_path / "inner_iterations.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert "q_residual" in rows[-1]
    if mode == "oss":
        assert "eta_q_residual" in rows[-1]
        assert stepper._scalar_states["q"] is not stepper._scalar_states["temperature"]
        assert stepper._layout.names == (
            "velocity",
            "pressure",
            "temperature",
            "q",
            "xi",
            "eta_temperature",
            "eta_q",
            "zeta",
        )
    with (tmp_path / "diagnostics.csv").open() as stream:
        assert "passive_scalar_content_q" in next(csv.DictReader(stream))


@pytest.mark.parametrize("name", ["temperature", "pressure", "eta_q", "bad name"])
def test_passive_names_are_checked(name: str) -> None:
    problem = build_problem(configuration())
    t = problem.thermal
    scalar = PassiveScalarProblem(
        name, t.diffusivity, t.source, t.initial_temperature, t.boundary_conditions
    )
    with pytest.raises(ValueError, match="names"):
        check_problem(replace(problem, passive_scalars=(scalar,)))


def test_passive_scalar_has_no_feedback_into_buoyant_flow(tmp_path: Path) -> None:
    config = configuration()
    problem = build_problem(config)
    thermal = problem.thermal
    scalar = PassiveScalarProblem(
        "q",
        2 * thermal.diffusivity,
        lambda x, _t: 10 + x[0],
        lambda x: x[0],
        thermal.boundary_conditions,
    )
    base_folder = tmp_path / "base"
    passive_folder = tmp_path / "passive"
    base_folder.mkdir()
    passive_folder.mkdir()
    base = TimeStepper(config, problem, base_folder, logging.getLogger("b2"))
    passive = TimeStepper(
        config, replace(problem, passive_scalars=(scalar,)), passive_folder, logging.getLogger("b2")
    )
    base.run()
    passive.run()
    assert (
        fd.assemble(
            ufl.inner(base.velocity - passive.velocity, base.velocity - passive.velocity) * base._dx
        )
        < 1e-18
    )
    assert fd.assemble((base.temperature - passive.temperature) ** 2 * base._dx) < 1e-18


def test_duplicate_passive_names_are_rejected() -> None:
    problem = build_problem(configuration())
    thermal = problem.thermal
    scalar = PassiveScalarProblem(
        "q",
        thermal.diffusivity,
        thermal.source,
        thermal.initial_temperature,
        thermal.boundary_conditions,
    )
    with pytest.raises(ValueError, match="names"):
        check_problem(replace(problem, passive_scalars=(scalar, scalar)))
    with pytest.raises(ValueError, match="require thermal"):
        check_problem(replace(problem, thermal=None, passive_scalars=(scalar,)))


def test_constant_gauge_is_confined_to_pressure(tmp_path: Path) -> None:
    config = configuration()
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger("b2"))
    stepper._initialise()
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    space = stepper._solution.function_space()
    # The pressure constant is a null mode; a temperature/projection constant has its own equation.
    jacobian = stepper._step_solver._inner._newton._problem.J
    for name in ("pressure", "temperature", "xi", "eta_temperature", "zeta"):
        direction = fd.Function(space)
        direction.subfunctions[stepper._layout.index(name)].assign(1.0)
        action = fd.assemble(ufl.action(jacobian, direction))
        for bc in stepper._step_solver._inner._bcs:
            bc.zero(action)
        with action.dat.vec_ro as vector:
            size = vector.norm()
        assert size < 1e-12 if name == "pressure" else size > 1e-8


@pytest.mark.parametrize("dynamics,splitting", [("static", "nonlinear"), ("dynamic", "nonlinear")])
def test_later_thermal_variants_are_rejected(dynamics: str, splitting: str) -> None:
    config = configuration()
    config = replace(
        config,
        stabilization=replace(
            config.stabilization, subscale_dynamics=dynamics, scale_splitting=splitting
        ),
    )
    with pytest.raises(ValueError, match="linear splitting"):
        build_problem(config)
