"""Preparation checks only: case data, boundary traces, assembly and diagnostics; no solves."""

import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.cases import build_problem, gresho, gresho_exact
from navier_stokes_vms.config import (
    Config,
    MeshOptions,
    OutputOptions,
    StabilizationOptions,
    TimeOptions,
    check_diagnostics_known,
    load_config,
    validate_config,
)
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.formulation.convective_forms import physical_pressure
from navier_stokes_vms.postprocessing.diagnostics import (
    DIAGNOSTICS,
    check_diagnostics_supported,
    diagnostic_columns,
    evaluate_diagnostics,
)
from navier_stokes_vms.problem.boundary_conditions import FreeSlip, free_slip_bc
from navier_stokes_vms.problem.problem_definition import uses_pressure_null_space
from navier_stokes_vms.solver.time_stepper import TimeStepper

CONFIGS = Path(__file__).resolve().parents[2] / "configs" / "gresho_exact"
CONFIG_NAMES = (
    "emac_galerkin",
    "skew_symmetric_1_galerkin",
    "emac_vms_oss",
    "skew_symmetric_1_oss",
)


@pytest.fixture(autouse=True)
def forbid_time_integration(monkeypatch: pytest.MonkeyPatch) -> None:
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("This preparation test must not advance the benchmark")

    monkeypatch.setattr(TimeStepper, "run", forbidden)
    monkeypatch.setattr(TimeStepper, "_step", forbidden)


def prepared_stepper(name: str, tmp_path: Path, cells: int = 4) -> TimeStepper:
    config = replace(load_config(CONFIGS / f"{name}.toml"), mesh=MeshOptions(cells=(cells, cells)))
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger(__name__))
    stepper._initialise()
    stepper._set_level(1, step_coefficients(config.time.scheme, config.time.startup, 1))
    return stepper


@pytest.mark.parametrize("name", CONFIG_NAMES)
def test_prepared_configs_and_first_step(name: str) -> None:
    config = load_config(CONFIGS / f"{name}.toml")
    assert config.case == "gresho_exact"
    assert config.mesh.cells == (48, 48)
    assert config.case_parameters["viscosity"] == 0.0
    assert config.time == TimeOptions(scheme="crank_nicolson", dt=0.01, t_end=10, startup="none")
    first = step_coefficients(config.time.scheme, config.time.startup, 1)
    assert first.scheme == "crank_nicolson"
    assert first.solve_level == 0.5
    assert first.alpha0 == 2.0 and first.history_weights == (2.0,)
    assert first.extrapolation_weight_solve == 2.0 and first.extrapolation_weight_old == -1.0
    check_diagnostics_known(config, DIAGNOSTICS)
    problem = build_problem(replace(config, mesh=MeshOptions(cells=(4, 4))))
    check_diagnostics_supported(config.output.diagnostics, problem)
    if config.stabilization.method == "vms":
        assert config.stabilization.subscale_space == "oss"
        assert config.stabilization.subscale_dynamics == "dynamic"
        assert not config.stabilization.pressure_subscale
        assert config.solver.newton_tau == "frozen"
        assert "total_kinetic_energy" in config.output.diagnostics
    else:
        assert config.discretization.velocity_degree == 2
        assert config.discretization.pressure_degree == 1


def test_case_defaults_keep_exact_fields_and_impermeable_slip_walls() -> None:
    problem = build_problem(Config(case="gresho_exact", run_name="defaults"))
    assert problem.mesh.num_cells() == 2 * 48 * 48
    assert problem.viscosity == 0.0
    assert problem.exact_velocity is gresho.velocity is gresho_exact.velocity
    assert problem.exact_pressure is gresho.pressure is gresho_exact.pressure
    assert uses_pressure_null_space(problem)
    assert all(isinstance(bc, FreeSlip) for bc in problem.boundary_conditions.values())
    x = fd.SpatialCoordinate(problem.mesh)
    assert problem.forcing(x, 0.0) == ufl.zero(2)
    assert gresho.DEFAULT_CELLS == (32, 32)
    assert gresho.DEFAULT_PARAMETERS == {"viscosity": 1e-6}


@pytest.mark.parametrize("viscosity", [1e-6, -1.0, True, "0"])
def test_exact_case_rejects_changed_physics(viscosity: object) -> None:
    with pytest.raises(ValueError, match="exact Euler"):
        build_problem(
            Config(
                case="gresho_exact", run_name="invalid", case_parameters={"viscosity": viscosity}
            )
        )


def test_wall_conditions_leave_tangential_velocity_free() -> None:
    problem = build_problem(
        Config(case="gresho_exact", run_name="walls", mesh=MeshOptions(cells=(4, 4)))
    )
    space = fd.VectorFunctionSpace(problem.mesh, "CG", 2)
    velocity = fd.Function(space).assign(fd.Constant((1.0, 1.0)))
    for name, bc in problem.boundary_conditions.items():
        free_slip_bc(space, name, bc.normal_axis, problem.boundary_tags).apply(velocity)
    for name, bc in problem.boundary_conditions.items():
        ds = fd.ds(problem.boundary_tags[name], domain=problem.mesh, degree=4)
        assert fd.assemble(velocity[bc.normal_axis] ** 2 * ds) == 0.0
        assert fd.assemble(velocity[1 - bc.normal_axis] ** 2 * ds) > 0.5


@pytest.mark.parametrize("name", CONFIG_NAMES)
def test_inviscid_initial_jacobian_has_only_the_pressure_gauge(name: str, tmp_path: Path) -> None:
    # 8×8 includes basis functions supported wholly in the resting exterior. No factorization,
    # linear solve, Newton update, or time step is performed.
    stepper = prepared_stepper(name, tmp_path, cells=8 if name.endswith("oss") else 4)
    problem = stepper._step_solver._inner._newton._problem
    matrix = fd.assemble(problem.J, bcs=problem.bcs, mat_type="aij").petscmat
    dense = matrix.convert("dense").getDenseArray()
    assert np.isfinite(dense).all()
    assert np.all(np.any(dense != 0, axis=1))
    singular = np.linalg.svd(dense, compute_uv=False)
    assert np.count_nonzero(singular < 1e-11 * singular[0]) == 1
    assert np.isfinite(stepper._step_solver._inner.residual_norm())
    if name.endswith("oss"):
        assert len(stepper._solution.subfunctions) == 3  # u, p, ξ; no zero-weight ζ block
        stepper._conservation.capture(midpoint=True)
        budgets = evaluate_diagnostics(
            ("energy_budget", "momentum_budget"), stepper._diagnostic_state
        )
        assert all(value is not None and np.isfinite(value) for value in budgets.values())


def test_total_energy_uses_both_endpoint_fields_and_the_cross_term(tmp_path: Path) -> None:
    stepper = prepared_stepper("emac_vms_oss", tmp_path)
    state = stepper._diagnostic_state
    # Deliberately nonorthogonal data distinguish the requested norm from E_h + E_sub.
    state.velocity.assign(fd.Constant((2.0, 0.0)))
    state.subscale_velocity.assign(fd.Constant((1.0, 0.0)))
    state.solve_velocity.assign(fd.Constant((7.0, 0.0)))
    names = ("kinetic_energy", "subscale_kinetic_energy", "total_kinetic_energy")
    values = evaluate_diagnostics(names, state)
    assert tuple(values) == diagnostic_columns(names, state.problem)
    assert tuple(values.values()) == pytest.approx((2.0, 0.5, 4.5))
    assert state.subscale_time is state.velocity_time


@pytest.mark.parametrize("name", ["emac_galerkin", "emac_vms_oss"])
def test_emac_pressure_error_uses_physical_half_level_pressure_and_removes_gauge(
    name: str, tmp_path: Path
) -> None:
    stepper = prepared_stepper(name, tmp_path)
    state = stepper._diagnostic_state
    x = fd.SpatialCoordinate(state.problem.mesh)
    stepper._pressure.interpolate(x[0] - 2 * x[1] + 7)
    state.solve_velocity.interpolate(ufl.as_vector((1 + x[1], 2 - x[0])))
    state.velocity.interpolate(ufl.as_vector((4 + x[1], -3 * x[0])))
    stepper._normalize_pressure()
    p = stepper._pressure + ufl.inner(state.solve_velocity, state.solve_velocity) / 2
    error = p - gresho.pressure(x, state.pressure_time)
    area = fd.assemble(fd.Constant(1.0) * state.error_dx)
    mean = fd.assemble(error * state.error_dx) / area
    expected = fd.assemble((error - mean) ** 2 * state.error_dx) ** 0.5
    result = evaluate_diagnostics(("pressure_error_l2",), state)["pressure_error_l2"]
    assert result == pytest.approx(expected, rel=1e-13)
    assert float(state.pressure_time) == 0.005
    # An exact physical pressure with an arbitrary gauge must have zero error, even though
    # the modified EMAC pressure differs by the spatially varying kinetic correction.
    exact = gresho.pressure(x, state.pressure_time)
    modified = exact - ufl.inner(state.solve_velocity, state.solve_velocity) / 2 + 23
    reconstructed = physical_pressure(
        modified, state.solve_velocity, stepper._config.convection.form
    )
    shifted = replace(state, pressure=reconstructed)
    assert evaluate_diagnostics(("pressure_error_l2",), shifted)["pressure_error_l2"] < 1e-12


def test_total_energy_rejects_mismatched_static_time_levels() -> None:
    config = Config(
        case="gresho_exact",
        run_name="invalid",
        stabilization=StabilizationOptions(subscale_dynamics="static"),
        output=OutputOptions(diagnostics=("total_kinetic_energy",)),
    )
    with pytest.raises(ValueError, match="both velocities are at the endpoint"):
        validate_config(config)


@pytest.mark.parametrize("scheme", ["bdf1", "bdf2", "bdf3"])
def test_no_startup_is_only_supported_for_crank_nicolson(scheme: str) -> None:
    config = Config(
        case="gresho_exact", run_name="invalid", time=TimeOptions(scheme=scheme, startup="none")
    )
    with pytest.raises(ValueError, match="requires time.scheme"):
        validate_config(config)
