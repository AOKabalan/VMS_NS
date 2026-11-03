"""B7.2 affine mesh checks, collective horizontal means and startup output lifecycle."""

import csv
import json
import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.cases.kelvin_helmholtz_shear_layer import (
    checked_parameters,
    vertical_mesh_vertices,
)
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    TimeOptions,
    config_to_toml,
    load_config,
    validate_config,
)
from navier_stokes_vms.postprocessing.horizontal_profiles import (
    HorizontalProfiles,
    profile_quadrature,
)
from navier_stokes_vms.problem.problem_definition import HorizontalProfileGeometry
from navier_stokes_vms.solver.time_stepper import TimeStepper


def configuration(clustering: float = 0.0, degree: int = 1) -> Config:
    solver = SolverOptions()
    return Config(
        case="kelvin_helmholtz_shear_layer",
        run_name="b72_test",
        case_parameters={"clustering": clustering},
        mesh=MeshOptions(cells=(4, 48 if clustering else 24)),
        discretization=DiscretizationOptions(degree, degree),
        time=TimeOptions(scheme="bdf2", dt=0.02, t_end=0.06),
        solver=replace(solver, petsc={**solver.petsc, "mat_mumps_icntl_24": 1}),
        output=OutputOptions(vtk_every=0, horizontal_profiles_every=2),
    )


@pytest.mark.parametrize("strength", [-0.1, True, float("nan"), float("inf")])
def test_invalid_clustering_is_rejected(strength: object) -> None:
    with pytest.raises(ValueError, match="clustering"):
        checked_parameters({"clustering": strength})


def test_vertical_map_uniform_limit_and_plateau_resolution() -> None:
    uniform = vertical_mesh_vertices(checked_parameters({}), 48)
    tiny = vertical_mesh_vertices(checked_parameters({"clustering": 1e-12}), 48)
    assert np.max(abs(tiny - uniform)) < 1e-12
    clustered = vertical_mesh_vertices(checked_parameters({"clustering": 2}), 48)
    assert clustered[0] == -10 and clustered[-1] == 10
    assert np.max(abs(clustered + clustered[::-1])) < 1e-14
    spacing = np.diff(clustered)
    assert spacing.min() > 0
    assert spacing[23] < np.diff(uniform).min()
    assert spacing.max() <= 1
    with pytest.raises(ValueError, match="plateau"):
        vertical_mesh_vertices(checked_parameters({"clustering": 2}), 24)
    with pytest.raises(ValueError, match="collapses"):
        vertical_mesh_vertices(checked_parameters({"clustering": 1000}), 48)


@pytest.mark.parametrize("strength", [0.0, 2.0])
@pytest.mark.parametrize("degree", [1, 2])
def test_affine_mesh_metrics_and_interpolated_wall_data(
    tmp_path: Path, strength: float, degree: int
) -> None:
    config = configuration(strength, degree)
    problem = build_problem(config)
    quality = problem.case_metadata["mesh_quality"]
    rows = np.asarray(problem.case_metadata["vertical_vertices"])
    dx = problem.case_metadata["domain_length"] / config.mesh.cells[0]
    dz = np.diff(rows)
    assert quality["minimum_absolute_jacobian"] == pytest.approx(dx * dz.min())
    assert quality["maximum_absolute_jacobian"] == pytest.approx(dx * dz.max())
    assert quality["minimum_map_jacobian"] == pytest.approx(dz.min() / (20 / len(dz)))
    assert quality["minimum_edge"] == pytest.approx(min(dx, dz.min()))
    assert quality["maximum_aspect"] == pytest.approx(max((dx * dx + dz * dz) / (dx * dz)))
    assert problem.mesh.coordinates.ufl_element().degree() == 1
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b72"))
    stepper._initialise()
    state = stepper._initial_state()
    assert state["boundary_normal_l2"] == 0
    assert state["temperature_boundary_flux_l2"] < 1e-13
    assert state["q_boundary_flux_l2"] < 1e-13
    assert fd.assemble(stepper.velocity[0].dx(1) ** 2 * stepper._ds) < 1e-26
    stores = (stepper._subscales, *stepper._scalar_subscale_storage.values())
    for store in stores:
        for field in (*store._levels, store.history):
            with field.dat.vec_ro as vector:
                assert vector.norm() == 0
    for name in ("xi", "eta_temperature", "eta_q", "zeta"):
        field = stepper._solution.subfunctions[stepper._layout.index(name)]
        with field.dat.vec_ro as vector:
            assert vector.norm() == 0
    measured_area = fd.assemble(fd.Constant(1.0) * stepper._dx)
    assert measured_area == pytest.approx(20 * problem.case_metadata["domain_length"])


@pytest.mark.parametrize("degree", [1, 2, 4])
@pytest.mark.parametrize("clustering", [0.0, 2.0])
def test_horizontal_means_match_independent_overintegration(
    tmp_path: Path, degree: int, clustering: float
) -> None:
    config = configuration(clustering)
    problem = build_problem(config)
    mesh = problem.mesh
    scalar_space = fd.FunctionSpace(mesh, "CG", degree)
    vector_space = fd.VectorFunctionSpace(mesh, "CG", degree)
    random = np.random.default_rng(123)
    velocity = fd.Function(vector_space)
    tracer = fd.Function(scalar_space)
    velocity.dat.data[:] = random.normal(size=velocity.dat.data.shape)
    tracer.dat.data[:] = random.normal(size=tracer.dat.data.shape)
    writer = HorizontalProfiles(
        tmp_path, problem.horizontal_profile_geometry, velocity, {"q": tracer}
    )
    writer.write(0, 0.0)
    writer.close()
    geometry = problem.horizontal_profile_geometry
    heights, _points = profile_quadrature(geometry, degree)
    nodes, weights = np.polynomial.legendre.leggauss(5)
    half_width = geometry.length / (2 * geometry.horizontal_cells)
    horizontal = (
        half_width * (np.arange(2 * geometry.horizontal_cells)[:, None] + (nodes + 1) / 2)
    ).ravel()
    points = np.column_stack(
        (np.tile(horizontal, len(heights)), np.repeat(heights, len(horizontal)))
    )
    evaluator = fd.PointEvaluator(mesh, points)
    line_weights = np.tile(weights / 2, 2 * geometry.horizontal_cells) / (
        2 * geometry.horizontal_cells
    )
    expected_q = evaluator.evaluate(tracer).reshape(len(heights), -1) @ line_weights
    expected_u = evaluator.evaluate(velocity).reshape(len(heights), -1, 2)
    with (tmp_path / "horizontal_profiles.csv").open() as stream:
        records = list(csv.DictReader(stream))
    for i, record in enumerate(records):
        assert float(record["q"]) == pytest.approx(expected_q[i], abs=2e-13)
        assert float(record["velocity_x"]) == pytest.approx(
            expected_u[i, :, 0] @ line_weights, abs=2e-13
        )
        assert float(record["velocity_y"]) == pytest.approx(
            expected_u[i, :, 1] @ line_weights, abs=2e-13
        )


@pytest.mark.parametrize(
    "geometry",
    [
        HorizontalProfileGeometry(0, 2, (-1, 0, 1)),
        HorizontalProfileGeometry(1, 0, (-1, 1)),
        HorizontalProfileGeometry(1, 2, (0, 0)),
        HorizontalProfileGeometry(1, 2, (1, 0)),
        HorizontalProfileGeometry(1, 2, (0, float("nan"))),
    ],
)
def test_invalid_profile_geometry_is_rejected(geometry: HorizontalProfileGeometry) -> None:
    with pytest.raises(ValueError, match="horizontal profiles"):
        profile_quadrature(geometry)


@pytest.mark.parametrize("name", ["time", "velocity_x"])
def test_profile_column_conflicts_are_rejected(tmp_path: Path, name: str) -> None:
    problem = build_problem(configuration())
    velocity = fd.Function(fd.VectorFunctionSpace(problem.mesh, "CG", 1))
    scalar = fd.Function(fd.FunctionSpace(problem.mesh, "CG", 1))
    with pytest.raises(ValueError, match="reserved CSV columns"):
        HorizontalProfiles(tmp_path, problem.horizontal_profile_geometry, velocity, {name: scalar})
    assert not (tmp_path / "horizontal_profiles.csv").exists()


@pytest.mark.parametrize("value", [-1, True, 0.5])
def test_profile_output_interval_validation(value: object) -> None:
    config = replace(configuration(), output=OutputOptions(horizontal_profiles_every=value))
    with pytest.raises(ValueError, match="horizontal_profiles_every"):
        validate_config(config)


def test_profile_output_rejects_unsupported_cases() -> None:
    config = replace(configuration(), case="manufactured_spatial", case_parameters={})
    with pytest.raises(ValueError, match="horizontal profile geometry"):
        build_problem(config)


def test_profile_option_roundtrip_and_default_schema(tmp_path: Path) -> None:
    config = configuration()
    path = tmp_path / "profiles.toml"
    path.write_text(config_to_toml(config))
    assert load_config(path) == config
    assert "horizontal_profiles_every" not in config_to_toml(
        replace(config, output=OutputOptions())
    )


@pytest.mark.parametrize("stop_early", [False, True])
def test_profile_lifecycle_and_named_scalar_output(tmp_path: Path, stop_early: bool) -> None:
    config = configuration()
    problem = build_problem(config)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("b72"))
    stepper.run((lambda step, _time: step == 1) if stop_early else None)
    with (tmp_path / "horizontal_profiles.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    steps = {int(row["step"]) for row in rows}
    assert steps == ({0, 1} if stop_early else {0, 2, 3})
    assert {"velocity_x", "velocity_y", "temperature", "q"} <= rows[0].keys()
    assert len(rows) == len(steps) * (2 * config.mesh.cells[1] + 1)
    assert all(np.isfinite(float(value)) for row in rows for value in row.values())
    metadata = json.loads((tmp_path / "horizontal_profiles.json").read_text())
    assert metadata["time_level"].startswith("endpoint")


def test_profiles_are_absent_when_disabled(tmp_path: Path) -> None:
    config = replace(configuration(), output=OutputOptions(vtk_every=0))
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger("b72"))
    stepper.run()
    assert not (tmp_path / "horizontal_profiles.csv").exists()
    assert not (tmp_path / "horizontal_profiles.json").exists()


def test_converged_startup_does_not_depend_on_pressure_guess(tmp_path: Path) -> None:
    config = replace(
        configuration(),
        case_parameters={"epsilon": 0.0, "moisture": False},
        time=TimeOptions(scheme="bdf2", dt=0.02, t_end=0.04),
        output=OutputOptions(vtk_every=0),
    )
    problem = build_problem(config)
    steppers = []
    for name, guess in (("hydrostatic", problem.initial_pressure), ("zero", lambda x: 0 * x[1])):
        folder = tmp_path / name
        folder.mkdir()
        stepper = TimeStepper(
            config, replace(problem, initial_pressure=guess), folder, logging.getLogger("b72")
        )
        stepper.run()
        steppers.append(stepper)
    first, second = steppers
    for difference in (
        first.velocity - second.velocity,
        first.temperature - second.temperature,
        first.pressure - second.pressure,
    ):
        assert fd.assemble(ufl.inner(difference, difference) * first._dx) ** 0.5 < 1e-10
    assert fd.assemble(first.velocity[1] ** 2 * first._dx) > 1e-10


def test_fourier_profiles_match_piecewise_linear_sine_and_cosine(tmp_path: Path) -> None:
    config = replace(
        configuration(),
        mesh=MeshOptions(cells=(8, 24)),
        output=OutputOptions(
            vtk_every=0, horizontal_profiles_every=1, horizontal_fourier_modes=(1, 2)
        ),
    )
    path = tmp_path / "fourier.toml"
    path.write_text(config_to_toml(config))
    assert load_config(path) == config
    assert "horizontal_fourier_modes" not in config_to_toml(configuration())
    problem = build_problem(config)
    geometry = problem.horizontal_profile_geometry
    x, z = fd.SpatialCoordinate(problem.mesh)
    angle = 2 * ufl.pi * x / geometry.length
    velocity = fd.Function(fd.VectorFunctionSpace(problem.mesh, "CG", 1)).interpolate(
        ufl.as_vector((z, 2 * ufl.sin(angle)))
    )
    tracer = fd.Function(fd.FunctionSpace(problem.mesh, "CG", 1)).interpolate(
        3 + ufl.cos(2 * angle)
    )
    writer = HorizontalProfiles(tmp_path, geometry, velocity, {"q": tracer}, (1, 2))
    writer.write(0, 0)
    writer.close()
    with (tmp_path / "horizontal_profiles.csv").open() as stream:
        records = list(csv.DictReader(stream))
    for record in records:
        assert float(record["velocity_y_mode_1_sin"]) == pytest.approx(
            np.sinc(1 / geometry.horizontal_cells) ** 2, abs=2e-13
        )
        assert float(record["q_mode_2_cos"]) == pytest.approx(
            0.5 * np.sinc(2 / geometry.horizontal_cells) ** 2, abs=2e-13
        )
        for name in ("velocity_x_mode_1_sin", "velocity_x_mode_1_cos", "q_mode_1_cos"):
            assert abs(float(record[name])) < 1e-13
    for modes in ((0,), (True,), (1, 1)):
        with pytest.raises(ValueError, match="horizontal_fourier_modes"):
            validate_config(
                replace(config, output=replace(config.output, horizontal_fourier_modes=modes))
            )
    with pytest.raises(ValueError, match="horizontal_profiles_every"):
        validate_config(replace(config, output=OutputOptions(horizontal_fourier_modes=(1,))))
