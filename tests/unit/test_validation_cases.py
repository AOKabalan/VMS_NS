"""Reference transcription, analytic validation comparisons and MPI-compatible point sampling."""

import json
import math
from pathlib import Path
from zipfile import ZipFile

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.cases import (
    build_problem,
    cylinder_dfg_2d,
    lattice_vortex,
    lid_driven_cavity,
)
from navier_stokes_vms.config import Config, MeshOptions, load_config
from navier_stokes_vms.problem.boundary_conditions import check_dirichlet_flux, dirichlet_bc

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("reynolds", [100, 1000])
def test_ghia_transcription_and_comparison(reynolds: int) -> None:
    profiles = lid_driven_cavity.reference_profiles(reynolds)
    expected_y = [
        1,
        0.9766,
        0.9688,
        0.9609,
        0.9531,
        0.8516,
        0.7344,
        0.6172,
        0.5,
        0.4531,
        0.2813,
        0.1719,
        0.1016,
        0.0703,
        0.0625,
        0.0547,
        0,
    ]
    expected_x = [
        1,
        0.9688,
        0.9609,
        0.9531,
        0.9453,
        0.9063,
        0.8594,
        0.8047,
        0.5,
        0.2344,
        0.2266,
        0.1563,
        0.0938,
        0.0781,
        0.0703,
        0.0625,
        0,
    ]
    np.testing.assert_array_equal(profiles["u"][:, 0], expected_y)
    np.testing.assert_array_equal(profiles["v"][:, 0], expected_x)
    np.testing.assert_array_equal(profiles["u"][[0, -1], 1], [1, 0])
    np.testing.assert_array_equal(profiles["v"][[0, -1], 1], [0, 0])
    samples = {name: values[:, 1].copy() for name, values in profiles.items()}
    assert all(value == 0 for value in lid_driven_cavity.compare(samples, reynolds).values())
    samples["u"][5] += 0.123
    samples["v"][10] -= 0.456
    assert lid_driven_cavity.compare(samples, reynolds) == pytest.approx(
        {"u_max_absolute_deviation": 0.123, "v_max_absolute_deviation": 0.456}
    )


def test_cavity_corners_and_global_point_values() -> None:
    problem = build_problem(
        Config(case="lid_driven_cavity", run_name="test", mesh=MeshOptions(cells=(4, 4)))
    )
    space = fd.VectorFunctionSpace(problem.mesh, "CG", 2)
    velocity = fd.Function(space)
    x = fd.SpatialCoordinate(problem.mesh)
    for name, bc in problem.boundary_conditions.items():
        dirichlet_bc(space, name, bc.value(x, fd.Constant(0)), problem.boundary_tags).apply(
            velocity
        )
    corners = fd.PointEvaluator(problem.mesh, [(0, 1), (1, 1), (0.5, 1)]).evaluate(velocity)
    np.testing.assert_allclose(corners, [[0, 0], [0, 0], [1, 0]], atol=1e-14)
    check_dirichlet_flux(velocity, fd.ds(domain=problem.mesh, degree=6))
    velocity.interpolate(ufl.as_vector((x[1], 2 * x[0])))
    samples = lid_driven_cavity.sample_profiles(velocity, 100)
    refs = lid_driven_cavity.reference_profiles(100)
    np.testing.assert_allclose(samples["u"], refs["u"][:, 0], atol=1e-14)
    np.testing.assert_allclose(samples["v"], 2 * refs["v"][:, 0], atol=1e-14)


def test_cylinder_names_geometry_and_interior_outlet_minimum() -> None:
    config = Config(
        case="cylinder_dfg_2d",
        run_name="test",
        mesh=MeshOptions(file=str(ROOT / "meshes/cylinder_dfg_2d_coarse.msh")),
    )
    problem = build_problem(config)
    assert problem.boundary_tags == {"inlet": 11, "outlet": 12, "walls": 13, "cylinder": 14}
    assert problem.force_boundaries == ("cylinder",)
    x = fd.SpatialCoordinate(problem.mesh)
    ds = fd.ds(problem.boundary_tags["outlet"], domain=problem.mesh, degree=4)
    assert fd.assemble(fd.Constant(1) * ds) == pytest.approx(0.41)
    assert fd.assemble((x[0] - 2.2) ** 2 * ds) < 1e-26
    velocity = fd.Function(fd.VectorFunctionSpace(problem.mesh, "CG", 2))
    velocity.interpolate(ufl.as_vector(((x[1] - 0.183) ** 2 - 0.01, 0)))
    evaluator = cylinder_dfg_2d.outlet_evaluator(problem)
    assert cylinder_dfg_2d.outlet_minimum(velocity, evaluator) == pytest.approx(-0.01, abs=1e-13)
    assert problem.reference_values["st"] == 0.30188


def test_sine_lift_cycle_statistics_and_phase_alignment() -> None:
    times = np.arange(0, 4.001, 0.001)
    frequency = 3.0
    lift = 0.2 + 1.5 * np.sin(2 * np.pi * frequency * times + 0.31)
    drag = 3 + 0.1 * np.cos(4 * np.pi * frequency * times + 0.62)
    metrics = cylinder_dfg_2d.compare(times, drag, lift)
    assert metrics["st"] == pytest.approx(0.3, rel=1e-6)
    for name, expected in [
        ("cl_min", -1.3),
        ("cl_max", 1.7),
        ("cl_mean", 0.2),
        ("cl_amplitude", 3),
        ("cd_min", 2.9),
        ("cd_max", 3.1),
        ("cd_mean", 3),
        ("cd_amplitude", 0.2),
    ]:
        assert metrics[name] == pytest.approx(expected, abs=1e-4)
    assert metrics["cl_max_change_last_three_cycles"] < 1e-4
    cycle = cylinder_dfg_2d.last_cycle(times, drag, lift)
    shifted = cylinder_dfg_2d.last_cycle(times + 7, drag, lift)
    assert max(cylinder_dfg_2d.compare_curves(cycle, shifted).values()) < 1e-12
    cycle[:, 1] += 0.25
    assert cylinder_dfg_2d.compare_curves(cycle, shifted)[
        "cd_curve_max_absolute_deviation"
    ] == pytest.approx(0.25)


def test_cylinder_insufficient_cycles_are_not_invented() -> None:
    times = np.linspace(0, 1, 10)
    with pytest.raises(ValueError, match="no complete"):
        cylinder_dfg_2d.compare(times, times, times)


def test_reference_archives_and_new_curve_load() -> None:
    for path in cylinder_dfg_2d.REFERENCES.glob("*.zip"):
        with ZipFile(path) as archive:
            assert archive.testzip() is None
    curve = cylinder_dfg_2d.reference_curve()
    assert 0.32 < curve[-1, 0] < 0.34
    assert np.all(np.diff(curve[:, 0]) > 0)
    assert 3.2 < np.max(curve[:, 1]) < 3.3
    assert 0.9 < np.max(curve[:, 2]) < 1.1
    intervals = json.loads((cylinder_dfg_2d.REFERENCES / "dfg_1996.json").read_text())
    assert intervals == {"cd_max": [3.22, 3.24], "cl_max": [0.99, 1.01], "st": [0.295, 0.305]}
    assert (
        json.loads((ROOT / "data/references/lattice_vortex/decay.json").read_text())[
            "initial_energy"
        ]
        == 0.25
    )


def test_lattice_formula_forcing_flux_and_energy() -> None:
    problem = build_problem(
        Config(case="lattice_vortex", run_name="test", mesh=MeshOptions(cells=(8, 8)))
    )
    x, t = fd.SpatialCoordinate(problem.mesh), fd.Constant(0.7)
    force = problem.forcing(x, t)
    assert fd.assemble(fd.inner(force, force) * fd.dx(degree=8)) ** 0.5 < 1e-13
    velocity = fd.Function(fd.VectorFunctionSpace(problem.mesh, "CG", 1)).interpolate(
        problem.exact_velocity(x, t)
    )
    check_dirichlet_flux(velocity, fd.ds(domain=problem.mesh, degree=3))
    times = np.linspace(0, 10, 21)
    exact = 0.25 * np.exp(-16 * 7e-5 * math.pi**2 * times)
    assert lattice_vortex.compare(times, exact)["energy_relative_error_at_T"] == 0
    assert lattice_vortex.compare(times, 1.1 * exact)["energy_max_relative_error"] == pytest.approx(
        0.1
    )


def test_every_validation_config_loads_with_requested_pairs_and_time_steps() -> None:
    configs = [load_config(path) for path in (ROOT / "configs/validation").glob("*.toml")]
    assert len(configs) == 13
    assert sum(c.case == "lid_driven_cavity" for c in configs) == 7
    assert sum(c.case == "lattice_vortex" for c in configs) == 2
    assert sum(c.case == "cylinder_dfg_2d" for c in configs) == 4
    for config in configs:
        if config.case == "cylinder_dfg_2d":
            assert config.time.t_end == 30 and config.time.dt in (0.005, 0.0025)
        elif config.case == "lattice_vortex":
            assert (
                config.mesh.cells == (50, 50) and config.time.dt == 0.01 and config.time.t_end == 10
            )
