"""B4 nonlinear coupled splitting, closure, linearization and convergence records."""

import csv
import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import pytest
import ufl
from firedrake.exceptions import ConvergenceError

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    ConvectionOptions,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    ThermalOptions,
    TimeOptions,
)
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.formulation.galerkin import time_derivative
from navier_stokes_vms.formulation.stabilization_parameters import (
    dynamic_parameter,
    transport_parameter,
)
from navier_stokes_vms.problem.scalar_problem import PassiveScalarProblem
from navier_stokes_vms.solver.time_stepper import TimeStepper


def configuration(
    mode: str = "oss", dynamics: str = "dynamic", splitting: str = "nonlinear"
) -> Config:
    return Config(
        case="manufactured_boussinesq",
        run_name="b4",
        mesh=MeshOptions(cells=(3, 3)),
        discretization=DiscretizationOptions(1, 1),
        stabilization=StabilizationOptions(
            subscale_space=mode, subscale_dynamics=dynamics, scale_splitting=splitting
        ),
        time=TimeOptions(scheme="bdf2", dt=0.01, t_end=0.02),
        output=OutputOptions(vtk_every=0),
    )


def stepper_with_passive(config: Config, folder: Path) -> TimeStepper:
    problem = build_problem(config)
    thermal = problem.thermal
    scalar = PassiveScalarProblem(
        "q",
        thermal.diffusivity,
        thermal.source,
        thermal.initial_temperature,
        thermal.boundary_conditions,
        thermal.exact_temperature,
    )
    return TimeStepper(
        config, replace(problem, passive_scalars=(scalar,)), folder, logging.getLogger("b4")
    )


@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize(
    "dynamics,splitting,update",
    [
        ("static", "linear", "accelerated"),
        ("dynamic", "linear", "accelerated"),
        ("static", "nonlinear", "accelerated"),
        ("dynamic", "nonlinear", "accelerated"),
        ("static", "nonlinear", "fixed_point"),
        ("dynamic", "nonlinear", "fixed_point"),
    ],
)
def test_variant_matrix_newton_and_picard_agree(
    tmp_path: Path, mode: str, dynamics: str, splitting: str, update: str
) -> None:
    config = configuration(mode, dynamics, splitting)
    steppers = []
    for method in ("newton", "picard"):
        folder = tmp_path / method
        folder.mkdir()
        options = SolverOptions(linearization=method, subscale_update=update)
        current_config = replace(config, solver=options)
        stepper = (
            TimeStepper(current_config, steppers[0]._problem, folder, logging.getLogger("b4"))
            if steppers
            else stepper_with_passive(current_config, folder)
        )
        stepper.run()
        steppers.append(stepper)
        assert fd.assemble((stepper.temperature - stepper.scalar("q")) ** 2 * stepper._dx) < 1e-20
        with (folder / "inner_iterations.csv").open() as stream:
            inner = list(csv.DictReader(stream))
        assert {row["method"] for row in inner} == {method}
        with (folder / "outer_iterations.csv").open() as stream:
            outer = list(csv.DictReader(stream))
        final = outer[-1]
        assert final["converged"] == "True"
        if splitting == "nonlinear":
            assert float(final["velocity_subscale_residual"]) < options.outer_rtol
            assert float(final["temperature_subscale_change"]) < options.outer_rtol
            assert float(final["q_subscale_change"]) < options.outer_rtol
            for name in stepper._layout.names:
                assert float(final[f"{name}_residual"]) < 1e-10
                assert float(final[f"{name}_change"]) < options.outer_rtol
    for name in steppers[0]._layout.names:
        i = steppers[0]._layout.index(name)
        difference = steppers[0]._solution.subfunctions[i] - steppers[1]._solution.subfunctions[i]
        assert fd.assemble(ufl.inner(difference, difference) * steppers[0]._dx) < 1e-17


@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
def test_nonlinear_subscale_closures_and_oss_orthogonality(
    tmp_path: Path, mode: str, dynamics: str
) -> None:
    config = configuration(mode, dynamics)
    config = replace(config, time=replace(config.time, scheme="crank_nicolson"))
    stepper = stepper_with_passive(config, tmp_path)
    stepper._initialise()
    startup = step_coefficients("crank_nicolson", "lower_order", 1)
    stepper._set_level(1, startup)
    assert stepper._step_solver.solve().converged
    stepper._advance(startup)
    stepper._set_level(2, step_coefficients("crank_nicolson", "lower_order", 2))
    result = stepper._step_solver.solve()
    assert result.converged, result.failure_reason
    dx = stepper._dx
    u, _, temperature, moisture, *_ = ufl.split(stepper._solution)
    a = u + stepper._iterate.value
    assert (
        fd.assemble(
            ufl.inner(
                stepper._u_sub - stepper._iterate.value, stepper._u_sub - stepper._iterate.value
            )
            * dx
        )
        < 1e-20
    )
    for name, field in (("temperature", temperature), ("q", moisture)):
        scalar = stepper._scalar_data[name]
        tau = transport_parameter(
            a, scalar.diffusivity, 1, config.stabilization.c1, config.stabilization.c2
        )
        residual = scalar.source(fd.SpatialCoordinate(stepper._problem.mesh), stepper._t_new)
        # CN uses endpoint averages even when the manufactured source is nonlinear in time.
        residual = (
            residual + scalar.source(fd.SpatialCoordinate(stepper._problem.mesh), stepper._t_old)
        ) / 2
        residual -= ufl.dot(a, ufl.grad(field))
        residual += scalar.diffusivity * ufl.div(ufl.grad(field))
        if mode == "asgs":
            residual -= time_derivative(
                field,
                stepper._scalar_states[name].history,
                stepper._alpha0,
                fd.Constant(config.time.dt),
            )
        if dynamics == "dynamic":
            tau = dynamic_parameter(tau, stepper._alpha0, fd.Constant(config.time.dt))
            residual += stepper._scalar_subscale_storage[name].history / config.time.dt
        if mode == "oss":
            residual -= ufl.split(stepper._solution)[stepper._layout.index(f"eta_{name}")]
        subscale = stepper._scalar_subscales[name]
        assert fd.assemble((subscale / tau - residual) ** 2 * dx) < 1e-20
        if mode == "oss":
            test = fd.TestFunction(stepper._scalar_fields[name].function_space())
            with fd.assemble(subscale * test * dx).dat.vec_ro as vector:
                assert vector.norm() < 1e-11
    if mode == "oss":
        test = fd.TestFunction(stepper._velocity.function_space())
        with fd.assemble(ufl.inner(stepper._u_sub, test) * dx).dat.vec_ro as vector:
            assert vector.norm() < 1e-11


@pytest.mark.parametrize("mode", ["none", "asgs", "oss"])
def test_picard_warmup_and_bdf3_nonlinear_history(tmp_path: Path, mode: str) -> None:
    config = configuration("oss" if mode == "none" else mode)
    if mode == "none":
        config = replace(
            config,
            discretization=DiscretizationOptions(2, 1),
            stabilization=StabilizationOptions(method="none"),
        )
    config = replace(
        config,
        time=replace(config.time, scheme="bdf3", t_end=0.04),
        solver=SolverOptions(picard_iterations_before_newton=1),
    )
    stepper = stepper_with_passive(config, tmp_path)
    stepper.run()
    with (tmp_path / "inner_iterations.csv").open() as stream:
        assert {r["method"] for r in csv.DictReader(stream)} == {"picard", "newton"}
    assert fd.assemble((stepper.temperature - stepper.scalar("q")) ** 2 * stepper._dx) < 1e-20


@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
@pytest.mark.parametrize("convection", ["nonconservative", "skew_symmetric"])
def test_nonlinear_scalar_cross_terms_are_present_once(
    tmp_path: Path, mode: str, dynamics: str, convection: str
) -> None:
    config = replace(configuration(mode, dynamics), thermal=ThermalOptions(convection=convection))
    stepper = stepper_with_passive(config, tmp_path)
    stepper._initialise()
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    fields = stepper._layout.bind(stepper._solution.subfunctions)
    x, y = fd.SpatialCoordinate(stepper._problem.mesh)
    u, temperature = fields["velocity"], fields["temperature"]
    u.interpolate(ufl.as_vector((1 + y, 1 + x)))
    temperature.interpolate(1 + 2 * x + 3 * y)
    stored = stepper._iterate.value
    stored.interpolate(ufl.as_vector((0.2 + 0.1 * y, -0.3 + 0.1 * x)))
    if mode == "oss":
        fields["eta_temperature"].interpolate(0.3 + y)
    scalar = stepper._scalar_data["temperature"]
    kappa = scalar.diffusivity
    dt = fd.Constant(config.time.dt)
    tau = transport_parameter(
        u + stored, kappa, 1, config.stabilization.c1, config.stabilization.c2
    )
    rate = time_derivative(
        temperature, stepper._scalar_states["temperature"].history, stepper._alpha0, dt
    )
    source = scalar.source(fd.SpatialCoordinate(stepper._problem.mesh), stepper._t_new)
    residual = source - ufl.dot(u, ufl.grad(temperature)) - ufl.dot(stored, ufl.grad(temperature))
    residual += kappa * ufl.div(ufl.grad(temperature))
    if mode == "asgs":
        residual -= rate
    if dynamics == "dynamic":
        history = stepper._scalar_subscale_storage["temperature"].history
        history.interpolate(0.1 + 0.2 * x)
        residual += history / dt
        tau = dynamic_parameter(tau, stepper._alpha0, dt)
    if mode == "oss":
        residual -= fields["eta_temperature"]
    subscale = tau * residual
    test = fd.TestFunction(temperature.function_space())
    dx = stepper._dx
    expected = (
        rate * test + kappa * ufl.dot(ufl.grad(temperature), ufl.grad(test)) - source * test
    ) * dx
    forward = (ufl.dot(u, ufl.grad(temperature)) + ufl.dot(stored, ufl.grad(temperature))) * test
    adjoint = ufl.dot(u, ufl.grad(test)) + ufl.dot(stored, ufl.grad(test))
    expected += (
        forward if convection == "nonconservative" else (forward - temperature * adjoint) / 2
    ) * dx
    expected -= subscale * (adjoint + kappa * ufl.div(ufl.grad(test))) * dx
    if dynamics == "dynamic" and mode == "asgs":
        expected += time_derivative(subscale, history, stepper._alpha0, dt) * test * dx
    actual = fd.assemble(stepper._step_solver._inner._newton_residual)
    with actual.subfunctions[2].dat.vec_ro as value, fd.assemble(expected).dat.vec_ro as target:
        difference = value.copy()
        difference.axpy(-1.0, target)
        assert difference.norm() < 1e-12 * target.norm()


@pytest.mark.parametrize("mode", ["asgs", "oss"])
def test_picard_freezes_transport_and_preserves_buoyancy_coupling(
    tmp_path: Path, mode: str
) -> None:
    config = replace(configuration(mode), solver=SolverOptions(linearization="picard"))
    stepper = stepper_with_passive(config, tmp_path)
    stepper._initialise()
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    solution = stepper._solution
    x, y = fd.SpatialCoordinate(stepper._problem.mesh)
    stepper._frozen.interpolate(ufl.as_vector((1 + y, 1 + x)))
    stepper._iterate.value.assign(0.1)
    direction = fd.Function(solution.function_space())
    direction.subfunctions[0].interpolate(ufl.as_vector((1 + x, 1 + y)))
    inner = stepper._step_solver._inner
    jacobian = inner._picard._problem.J
    velocity_action = fd.assemble(ufl.action(jacobian, direction))
    with velocity_action.subfunctions[2].dat.vec_ro as vector:
        assert vector.norm() < 1e-13
    true_action = fd.assemble(
        ufl.action(fd.derivative(inner._newton_residual, solution), direction)
    )
    with true_action.subfunctions[2].dat.vec_ro as vector:
        assert vector.norm() > 1e-4
    direction.assign(0.0)
    direction.subfunctions[2].assign(1.0)
    temperature_action = fd.assemble(ufl.action(jacobian, direction))
    with temperature_action.subfunctions[0].dat.vec_ro as vector:
        assert vector.norm() > 1e-4
    # A frozen-advector Picard residual is affine in all monolithic unknowns.
    base = solution.copy(deepcopy=True)
    direction.assign(0.2)
    center = fd.assemble(inner._picard._problem.F)
    solution.assign(base + direction)
    plus = fd.assemble(inner._picard._problem.F)
    solution.assign(base - direction)
    minus = fd.assemble(inner._picard._problem.F)
    with center.dat.vec_ro as vc, plus.dat.vec_ro as vp, minus.dat.vec_ro as vm:
        difference = vp.copy()
        difference.axpy(1.0, vm)
        difference.axpy(-2.0, vc)
        assert difference.norm() < 1e-12


@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
def test_local_failure_reports_and_preserves_coupled_histories(
    tmp_path: Path, dynamics: str
) -> None:
    config = replace(
        configuration("asgs", dynamics),
        solver=SolverOptions(subscale_max_iterations=1, subscale_rtol=1e-15),
    )
    stepper = stepper_with_passive(config, tmp_path)
    with pytest.raises(ConvergenceError, match="local subscale update did not converge") as error:
        stepper.run()
    if dynamics == "static":
        assert "static nonlinear subscales on a coarse mesh" in str(error.value)
    else:
        for storage in (stepper._subscales, *stepper._scalar_subscale_storage.values()):
            assert fd.assemble(ufl.inner(storage.latest, storage.latest) * stepper._dx) == 0
    x = fd.SpatialCoordinate(stepper._problem.mesh)
    for name, scalar in stepper._scalar_data.items():
        initial = fd.Function(stepper._scalar_fields[name].function_space()).interpolate(
            scalar.initial_value(x)
        )
        assert fd.assemble((stepper.scalar(name) - initial) ** 2 * stepper._dx) < 1e-20
        assert (tmp_path / "fields" / f"{name}.pvd").exists()
    with (tmp_path / "time_steps.csv").open() as stream:
        assert list(csv.DictReader(stream))[-1]["converged"] == "False"


@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
@pytest.mark.parametrize("momentum_form", ["nonconservative", "skew_symmetric_1"])
def test_nonlinear_zero_expansion_matches_isothermal(
    tmp_path: Path, mode: str, dynamics: str, momentum_form: str
) -> None:
    config = replace(
        configuration(mode, dynamics),
        case_parameters={"expansion_coefficient": 0.0},
        convection=ConvectionOptions(form=momentum_form),
    )
    problem = build_problem(config)
    coupled_folder, flow_folder = tmp_path / "coupled", tmp_path / "flow"
    coupled_folder.mkdir()
    flow_folder.mkdir()
    coupled = TimeStepper(config, problem, coupled_folder, logging.getLogger("b4"))
    flow = TimeStepper(config, replace(problem, thermal=None), flow_folder, logging.getLogger("b4"))
    coupled.run()
    flow.run()
    for coupled_field, flow_field in (
        (coupled.velocity, flow.velocity),
        (coupled.pressure, flow.pressure),
    ):
        difference = coupled_field - flow_field
        assert fd.assemble(ufl.inner(difference, difference) * coupled._dx) < 1e-17
