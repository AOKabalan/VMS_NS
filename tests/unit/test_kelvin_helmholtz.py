"""B7.1 profiles, initialization and coarse coupled transport, without growth-rate claims."""

import csv
import json
import logging
import math
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.cases.kelvin_helmholtz_shear_layer import (
    checked_parameters,
    wall_profiles,
)
from navier_stokes_vms.config import (
    Config,
    ConvectionOptions,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.postprocessing.diagnostics import advective_cfl_sampler
from navier_stokes_vms.problem.scalar_problem import PassiveScalarProblem
from navier_stokes_vms.solver.reusable_assembly import scalar_functional
from navier_stokes_vms.solver.time_stepper import TimeStepper


def configuration(mode: str = "oss", degree: int = 1) -> Config:
    options = SolverOptions()
    return Config(
        case="kelvin_helmholtz_shear_layer",
        run_name="b71",
        mesh=MeshOptions(cells=(8, 24)),
        discretization=DiscretizationOptions(degree, degree),
        stabilization=StabilizationOptions(
            subscale_space=mode, subscale_dynamics="static", scale_splitting="linear"
        ),
        time=TimeOptions(scheme="bdf1", dt=0.02, t_end=0.04),
        solver=replace(options, petsc={**options.petsc, "mat_mumps_icntl_24": 1}),
        output=OutputOptions(
            vtk_every=0,
            diagnostics=(
                "kinetic_energy",
                "vertical_kinetic_energy",
                "advective_cfl",
                "thermal_content",
                "passive_scalar_content",
            ),
        ),
    )


def norm(expression: Expr, measure: ufl.Measure) -> float:
    return scalar_functional(ufl.inner(expression, expression) * measure)() ** 0.5


def dual_norm(form: ufl.Form) -> float:
    value = fd.assemble(form)
    with value.dat.vec_ro as vector:
        return vector.norm()


@pytest.mark.parametrize(
    "parameters",
    [
        {"unknown": 1},
        {"reynolds": 0},
        {"prandtl": -1},
        {"schmidt": False},
        {"richardson": -0.1},
        {"epsilon": -1},
        {"epsilon_noise": -1},
        {"height": float("inf")},
        {"wave_number": 0},
        {"n_waves": True},
        {"n_waves": 1.5},
        {"random_seed": -1},
        {"random_seed": False},
        {"moisture": 1},
        {"noise_modes": [1, 1]},
        {"noise_modes": [True]},
        {"noise_modes": [0]},
        {"noise_modes": 1},
        {"noise_modes": [], "epsilon_noise": 1},
        {"rh_peak": 1.1},
        {"moisture_width": 0},
        {"saturation_height_rate": -1},
        {"saturation_temperature_rate": float("nan")},
        {"wall_buffer": 0},
        {"wall_transition_width": 2.6},
        {"height": 6},
    ],
)
def test_invalid_case_parameters_are_rejected(parameters: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        checked_parameters(parameters)


def test_zero_control_parameters_are_valid() -> None:
    parameters = checked_parameters(
        {
            "richardson": 0,
            "epsilon": 0,
            "epsilon_noise": 0,
            "noise_modes": [],
            "random_seed": 0,
            "saturation_height_rate": 0,
            "saturation_temperature_rate": -0.1,
        }
    )
    assert parameters.noise_modes == ()


@pytest.mark.parametrize("width", [0.5, 2.0, 2.5])
def test_profile_primitive_monotonicity_and_wall_compatibility(width: float) -> None:
    parameters = checked_parameters({"wall_transition_width": width})
    mesh = fd.RectangleMesh(2, 160, 1.0, 20.0)
    mesh.coordinates.dat.data[:, 1] -= 10
    x = fd.SpatialCoordinate(mesh)
    base, primitive, envelope = wall_profiles(x[1], parameters)
    dx = fd.dx(domain=mesh, degree=8)
    assert norm(primitive.dx(1) - base, dx) < 3e-14
    core = 10 - parameters.wall_buffer - width
    assert norm(ufl.conditional(abs(x[1]) < core, base - ufl.tanh(x[1]), 0), dx) < 1e-14
    sample_space = fd.FunctionSpace(mesh, "DG", 6)
    slope = fd.Function(sample_space).interpolate(base.dx(1))
    assert np.min(slope.dat.data_ro) >= -1e-14
    cutoff = fd.Function(sample_space).interpolate(envelope)
    assert np.min(cutoff.dat.data_ro) >= -1e-14
    assert np.max(cutoff.dat.data_ro) <= 1 + 1e-14
    walls = fd.ds((3, 4), domain=mesh, degree=8)
    for expression in (base.dx(1), envelope, envelope.dx(1), envelope.dx(1).dx(1)):
        assert norm(expression, walls) < 1e-14
    epsilon = 1e-9
    points = [(0.3, z + offset) for z in (-9, -core, core, 9) for offset in (-epsilon, epsilon)]
    for expression in (
        base,
        base.dx(1),
        base.dx(1).dx(1),
        primitive,
        envelope,
        envelope.dx(1),
        envelope.dx(1).dx(1),
    ):
        values = fd.Function(sample_space).interpolate(expression).at(points)
        assert np.max(np.abs(np.array(values)[::2] - np.array(values)[1::2])) < 1000 * epsilon


def test_physical_scaling_and_reproducible_noise() -> None:
    config = replace(
        configuration(),
        case_parameters={
            "n_waves": 2,
            "epsilon_noise": 0.01,
            "random_seed": 42,
            "noise_modes": [1, 3, 4],
            "reynolds": 2000,
            "prandtl": 2,
            "schmidt": 3,
        },
    )
    first, second = build_problem(config), build_problem(config)
    assert first.case_metadata == second.case_metadata
    assert first.viscosity == 1 / 2000
    assert first.thermal.diffusivity == 1 / 4000
    assert first.passive_scalars[0].diffusivity == 1 / 6000
    assert first.thermal.gravity == (0.0, -1.0)
    assert first.thermal.reference_temperature == 0
    assert first.thermal.scaling == "nondimensional"
    assert first.exact_velocity is None and first.exact_pressure is None
    assert first.thermal.exact_temperature is None
    length = first.case_metadata["domain_length"]
    assert length == pytest.approx(4 * math.pi / 0.4446)
    phases = first.case_metadata["noise_phases"]
    amplitude = first.case_metadata["noise_amplitude_per_mode"]
    grid = np.arange(1024) * length / 1024
    noise = sum(
        amplitude * np.cos(2 * math.pi * m * grid / length + p)
        for m, p in zip((1, 3, 4), phases, strict=True)
    )
    assert np.mean(noise**2) == pytest.approx(0.5, abs=1e-14)
    changed = build_problem(
        replace(config, case_parameters={**config.case_parameters, "random_seed": 43})
    )
    assert phases != changed.case_metadata["noise_phases"]


@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("mode", ["asgs", "oss"])
def test_initial_state_and_finite_coupled_forms(tmp_path: Path, mode: str, degree: int) -> None:
    config = configuration(mode, degree)
    problem = build_problem(config)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b71"))
    stepper._initialise()
    initial = stepper._initial_state()
    assert initial["boundary_normal_l2"] == 0
    assert abs(initial["boundary_flux"]) < 1e-15
    assert initial["temperature_boundary_flux_l2"] < 1e-14
    assert initial["q_boundary_flux_l2"] < 1e-14
    assert norm(stepper.velocity[0].dx(1), stepper._ds) < 1e-13
    assert abs(initial["pressure_mean"]) < 1e-14
    assert initial["vertical_kinetic_energy"] > 0
    assert initial["divergence_l2"] > 0
    assert initial["pressure_buoyancy_defect_l2"] > 0
    assert initial["advective_cfl"] > 0
    x = fd.SpatialCoordinate(problem.mesh)
    assert norm(ufl.div(problem.initial_velocity(x)), stepper._dx) < 1e-13
    target = fd.Function(stepper._pressure.function_space()).interpolate(
        problem.initial_pressure(x)
    )
    target.assign(target - fd.Constant(fd.assemble(target * stepper._dx) / stepper._area))
    assert norm(stepper.pressure - target, stepper._dx) < 1e-14
    assert norm(problem.forcing(x, fd.Constant(0.0)), stepper._dx) == 0
    assert scalar_functional(problem.thermal.source(x, fd.Constant(0.0)) * stepper._dx)() == 0
    assert stepper._layout.names[:4] == ("velocity", "pressure", "temperature", "q")
    if mode == "oss":
        assert norm(stepper._solution.subfunctions[stepper._layout.index("xi")], stepper._dx) == 0
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    inner = stepper._step_solver._inner
    residual = fd.assemble(inner._newton_residual)
    with residual.dat.vec_ro as vector:
        assert math.isfinite(vector.norm())
    matrix = fd.assemble(inner._newton._problem.J, mat_type="aij")
    assert math.isfinite(matrix.petscmat.norm())
    direction = fd.Function(stepper._solution.function_space())
    direction.subfunctions[1].assign(1.0)
    action = fd.assemble(ufl.action(inner._newton._problem.J, direction))
    for bc in inner._bcs:
        bc.zero(action)
    with action.dat.vec_ro as vector:
        assert vector.norm() < 1e-12


@pytest.mark.parametrize("mode", ["asgs", "oss"])
def test_coupled_smoke_preserves_duplicate_transport(tmp_path: Path, mode: str) -> None:
    config = configuration(mode)
    problem = build_problem(config)
    thermal = problem.thermal
    duplicate = PassiveScalarProblem(
        "temperature_copy",
        thermal.diffusivity,
        thermal.source,
        thermal.initial_temperature,
        thermal.boundary_conditions,
    )
    problem = replace(problem, passive_scalars=(*problem.passive_scalars, duplicate))
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b71"))
    stepper.run()
    assert norm(stepper.temperature - stepper.scalar("temperature_copy"), stepper._dx) < 1e-11
    assert abs(fd.assemble(stepper.pressure * stepper._dx)) < 1e-12
    with (tmp_path / "diagnostics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 2
    assert all(math.isfinite(float(value)) for row in rows for value in row.values())
    assert "passive_scalar_content_q" in rows[-1]
    assert "passive_scalar_content_temperature_copy" in rows[-1]
    metadata = json.loads((tmp_path / "thermal_metadata.json").read_text())
    assert metadata["case_metadata"]["parameters"]["random_seed"] == 1729


@pytest.mark.parametrize("degree", [1, 2])
def test_oss_orthogonality_with_moisture(tmp_path: Path, degree: int) -> None:
    config = configuration(degree=degree)
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger("b71"))
    stepper._initialise()
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    assert stepper._step_solver.solve().converged
    mesh = stepper._problem.mesh
    velocity_space = fd.VectorFunctionSpace(mesh, "CG", degree)
    scalar_space = fd.FunctionSpace(mesh, "CG", degree)
    assert (
        dual_norm(ufl.inner(stepper._u_sub, fd.TestFunction(velocity_space)) * stepper._dx) < 1e-9
    )
    for name in ("temperature", "q"):
        assert (
            dual_norm(stepper._scalar_subscales[name] * fd.TestFunction(scalar_space) * stepper._dx)
            < 1e-9
        )


def test_no_perturbation_and_dry_controls(tmp_path: Path) -> None:
    config = replace(
        configuration(), case_parameters={"epsilon": 0, "moisture": False, "richardson": 0}
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger("b71"))
    stepper._initialise()
    initial = stepper._initial_state()
    assert initial["vertical_kinetic_energy"] == 0
    assert initial["divergence_l2"] < 1e-14
    assert initial["pressure_buoyancy_defect_l2"] == 0
    assert stepper.passive_scalars == {}
    stepper.run()
    assert norm(stepper.velocity[1], stepper._dx) < 1e-11


@pytest.mark.parametrize("form", ["skew_symmetric_1", "emac", "emac_vms"])
def test_optional_initial_pressure_is_physical(tmp_path: Path, form: str) -> None:
    config = replace(
        configuration(),
        case="free_slip_channel",
        case_parameters={"thermal": False},
        convection=ConvectionOptions(form=form),
    )
    problem = replace(build_problem(config), initial_pressure=lambda x: x[1])
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b71"))
    stepper._initialise()
    x = fd.SpatialCoordinate(problem.mesh)
    target = x[1] - fd.Constant(fd.assemble(x[1] * stepper._dx) / stepper._area)
    assert norm(stepper.pressure - target, stepper._dx) < 1e-13


def test_cfl_uses_shortest_edge_degree_and_current_velocity() -> None:
    mesh = fd.RectangleMesh(2, 4, 2.0, 1.0)
    velocity = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 2)).assign(fd.Constant((3.0, 4.0)))
    sample = advective_cfl_sampler(velocity, 0.01, 2, 6)
    assert sample() == pytest.approx(2 * 0.01 * 5 / 0.25)
    velocity.assign(fd.Constant((0.0, 0.0)))
    assert sample() == 0
