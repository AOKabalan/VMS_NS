"""Global diagnostics on interpolated analytic fields and their time-stepper wiring."""

import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import Config, MeshOptions, StabilizationOptions, TimeOptions
from navier_stokes_vms.discretization.function_spaces import subscale_space
from navier_stokes_vms.discretization.quadrature import boundary_measure, volume_measure
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.postprocessing.diagnostics import (
    DiagnosticState,
    check_diagnostics_supported,
    diagnostic_columns,
    evaluate_diagnostics,
)
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.problem_definition import Problem, check_problem
from navier_stokes_vms.solver.time_stepper import TimeStepper


def analytic_state(dimension: int, degree: int = 2) -> DiagnosticState:
    constructor = mesh_loading.rectangle if dimension == 2 else mesh_loading.box
    mesh, tags = constructor((2,) * dimension, (0.0,) * dimension, (1.0,) * dimension)
    velocity = fd.Function(fd.VectorFunctionSpace(mesh, "CG", degree))
    pressure = fd.Function(fd.FunctionSpace(mesh, "CG", 1))
    problem = Problem(
        mesh=mesh,
        boundary_conditions={},
        boundary_tags=tags,
        viscosity=0.7,
        forcing=lambda _x, _t: ufl.zero(dimension),
        initial_velocity=lambda _x: ufl.zero(dimension),
        force_boundaries=tuple(tags),
    )
    return DiagnosticState(
        problem=problem,
        velocity=velocity,
        pressure=pressure,
        velocity_time=fd.Constant(0.0),
        pressure_time=fd.Constant(0.0),
        dx=volume_measure(mesh, 3 * degree),
        error_dx=fd.dx(domain=mesh, degree=6 * degree),
        pressure_null_space=False,
        solve_velocity=velocity,
        ds=boundary_measure(mesh, 3 * degree),
        subscale_velocity=None,
        subscale_time=None,
    )


@pytest.mark.parametrize("dimension", [2, 3])
@pytest.mark.parametrize("flow", ["uniform", "rotation"])
def test_momenta_of_analytic_fields(dimension: int, flow: str) -> None:
    state = analytic_state(dimension)
    x = fd.SpatialCoordinate(state.problem.mesh)
    if flow == "uniform":
        state.velocity.interpolate(ufl.as_vector((2.0, -3.0, 5.0)[:dimension]))
        momentum = (2.0, -3.0, 5.0)[:dimension]
        angular = (-2.5,) if dimension == 2 else (4.0, -1.5, -2.5)
    elif dimension == 2:
        state.velocity.interpolate(ufl.as_vector((-3 * x[1], 3 * x[0])))
        momentum, angular = (-1.5, 1.5), (2.0,)
    else:
        state.velocity.interpolate(ufl.cross(ufl.as_vector((1.0, 2.0, 3.0)), x))
        momentum, angular = (-0.5, 1.0, -0.5), (-7 / 12, 1 / 3, 5 / 4)
    names = ("linear_momentum", "angular_momentum")
    values = evaluate_diagnostics(names, state)
    assert tuple(values) == diagnostic_columns(names, state.problem)
    assert tuple(values.values()) == pytest.approx((*momentum, *angular), abs=1e-13)


@pytest.mark.parametrize("dimension", [2, 3])
def test_poiseuille_wall_force_sign_shear_and_pressure(dimension: int) -> None:
    state = analytic_state(dimension)
    x = fd.SpatialCoordinate(state.problem.mesh)
    amplitude, pressure_offset = 3.0, 5.0
    viscosity = state.problem.viscosity
    state.velocity.interpolate(
        ufl.as_vector((amplitude * x[1] * (1 - x[1]),) + (0,) * (dimension - 1))
    )
    state.pressure.interpolate(pressure_offset - 2 * viscosity * amplitude * x[0])
    values = evaluate_diagnostics(("force",), state)
    assert tuple(values) == diagnostic_columns(("force",), state.problem)
    mean_pressure = pressure_offset - viscosity * amplitude
    for wall, sign in (("y_min", -1), ("y_max", 1)):
        assert values[f"force_{wall}_x"] == pytest.approx(viscosity * amplitude)
        assert values[f"force_{wall}_y"] == pytest.approx(sign * mean_pressure)
        if dimension == 3:
            assert values[f"force_{wall}_z"] == pytest.approx(0, abs=1e-13)
    assert values["force_x_min_x"] == pytest.approx(-pressure_offset)
    assert values["force_x_max_x"] == pytest.approx(pressure_offset - 2 * viscosity * amplitude)
    assert values["force_x_min_y"] == pytest.approx(0, abs=1e-13)
    assert values["force_x_max_y"] == pytest.approx(0, abs=1e-13)
    if dimension == 3:
        assert values["force_z_min_z"] == pytest.approx(-mean_pressure)
        assert values["force_z_max_z"] == pytest.approx(mean_pressure)


@pytest.mark.parametrize("dimension", [2, 3])
@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("stored", [False, True])
def test_subscale_energy_uses_the_common_rule(dimension: int, degree: int, stored: bool) -> None:
    state = analytic_state(dimension, degree)
    x = fd.SpatialCoordinate(state.problem.mesh)
    expression = ufl.as_vector(tuple(ufl.exp((i + 1) * x[i]) for i in range(dimension)))
    subscale = expression
    if stored:
        subscale = fd.Function(subscale_space(state.problem.mesh, 3 * degree)).interpolate(
            expression
        )
    state = replace(state, subscale_velocity=subscale)
    expected = fd.assemble(
        sum(ufl.exp(2 * (i + 1) * x[i]) for i in range(dimension)) / 2 * state.dx
    )
    value = evaluate_diagnostics(("subscale_kinetic_energy",), state)["subscale_kinetic_energy"]
    assert value == pytest.approx(expected, rel=1e-14)
    other_rule = fd.assemble(ufl.inner(expression, expression) / 2 * state.error_dx)
    assert abs(expected - other_rule) > 1e-9 * expected


@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("scheme", ["bdf2", "crank_nicolson"])
def test_subscale_energy_reads_the_correct_field_and_time(
    tmp_path: Path, space: str, dynamics: str, splitting: str, scheme: str
) -> None:
    config = Config(
        case="manufactured_temporal",
        run_name="diagnostics",
        mesh=MeshOptions(cells=(2, 2)),
        time=TimeOptions(scheme=scheme, dt=0.1, t_end=0.2),
        stabilization=StabilizationOptions(
            subscale_space=space, subscale_dynamics=dynamics, scale_splitting=splitting
        ),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger(__name__))
    stepper._initialise()
    x = fd.SpatialCoordinate(stepper._problem.mesh)
    if stepper._subscales is not None:
        stepper._subscales.latest.interpolate(ufl.as_vector((0.3 + x[0], -0.2)))
    coefficients = step_coefficients(scheme, "lower_order", 2)
    stepper._set_level(2, coefficients)
    stepper._velocity.interpolate(ufl.as_vector((0.2 * x[1], -0.1 * x[0])))
    stepper._pressure.interpolate(x[0] - x[1])
    if stepper._iterate is not None:
        stepper._iterate.value.interpolate(ufl.as_vector((0.1, 0.2 + x[1])))
    reference = stepper._u_sub
    expected_time = 0.15 if scheme == "crank_nicolson" else 0.2
    if dynamics == "dynamic":
        expected_time = 0.2
        if scheme == "crank_nicolson":
            reference = 2 * reference - stepper._subscales.latest
    elif splitting == "nonlinear":
        reference = stepper._iterate.value
    expected = fd.assemble(ufl.inner(reference, reference) / 2 * stepper._dx)
    stepper._advance(coefficients)
    state = stepper._diagnostic_state
    value = evaluate_diagnostics(("subscale_kinetic_energy",), state)["subscale_kinetic_energy"]
    assert value == pytest.approx(expected, rel=1e-13)
    assert float(state.subscale_time) == pytest.approx(expected_time)


def test_crank_nicolson_force_uses_half_velocity_and_pressure(tmp_path: Path) -> None:
    config = Config(
        case="manufactured_temporal",
        run_name="half_level",
        mesh=MeshOptions(cells=(2, 2)),
        time=TimeOptions(scheme="crank_nicolson", dt=0.1, t_end=0.2),
    )
    problem = build_problem(config)
    problem = replace(
        problem,
        force_boundaries=("y_max",),
        exact_pressure=lambda x, t: t * x[0],
    )
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger(__name__))
    stepper._initialise()
    coefficients = step_coefficients("crank_nicolson", "lower_order", 2)
    stepper._set_level(2, coefficients)
    x = fd.SpatialCoordinate(problem.mesh)
    stepper._history[0].interpolate(ufl.as_vector((0.1 * x[1], 0)))
    stepper._velocity.interpolate(ufl.as_vector((0.15 * x[1], 0)))
    stepper._pressure.interpolate(0.15 * x[0])
    stepper._advance(coefficients)
    values = evaluate_diagnostics(("force", "pressure_error_l2"), stepper._diagnostic_state)
    assert values["force_y_max_x"] == pytest.approx(-0.15 * problem.viscosity)
    assert values["force_y_max_y"] == pytest.approx(0.075)
    assert values["pressure_error_l2"] < 1e-14
    assert fd.assemble(stepper.velocity[0] * stepper._dx) == pytest.approx(0.1)


def test_force_boundaries_are_checked_even_for_cases_with_exact_solutions() -> None:
    state = analytic_state(2)
    problem = replace(state.problem, exact_velocity=lambda x, _t: x, force_boundaries=())
    with pytest.raises(ValueError, match="force_boundaries"):
        check_diagnostics_supported(("force",), problem)
    problem = replace(problem, boundary_conditions={}, force_boundaries=("missing",))
    with pytest.raises(ValueError, match="boundary 'missing' does not exist"):
        check_problem(problem)
