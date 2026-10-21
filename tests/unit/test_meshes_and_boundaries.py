"""Meshes, boundary names, boundary conditions and case problems (plan U7 and phase 2)."""

import math
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl
from firedrake.mesh import MeshGeometry
from ufl.core.expr import Expr

from navier_stokes_vms.cases import CASES, build_problem
from navier_stokes_vms.config import Config, MeshOptions, TimeOptions
from navier_stokes_vms.problem.boundary_conditions import Dirichlet, Outflow, dirichlet_bc
from navier_stokes_vms.problem.mesh_loading import box, gmsh_mesh, periodic_rectangle, rectangle
from navier_stokes_vms.problem.problem_definition import (
    Problem,
    check_problem,
    uses_pressure_null_space,
)

SQUARE_NAMED = Path(__file__).resolve().parents[1] / "data" / "square_named.msh"
LOWER_2D, UPPER_2D = (-0.5, -0.25), (1.0, 1.5)
LOWER_3D, UPPER_3D = (-1.0, 0.5, -0.5), (1.0, 1.5, 0.25)


def side_of(name: str, lower: tuple[float, ...], upper: tuple[float, ...]) -> tuple[int, float]:
    axis = "xyz".index(name[0])
    return axis, lower[axis] if name.endswith("min") else upper[axis]


def mean_and_spread(mesh: MeshGeometry, tag: int, axis: int) -> tuple[float, float, float]:
    x = fd.SpatialCoordinate(mesh)
    ds = fd.ds(tag, domain=mesh, degree=2)
    length = fd.assemble(fd.Constant(1.0) * ds)
    mean = fd.assemble(x[axis] * ds) / length
    return length, mean, fd.assemble((x[axis] - mean) ** 2 * ds)


def zero_velocity(_x: Expr, _t: Expr) -> Expr:
    return ufl.as_vector((0.0, 0.0))


def dirichlet_square() -> Problem:
    mesh, tags = rectangle((2, 2), (0.0, 0.0), (1.0, 1.0))
    return Problem(
        mesh=mesh,
        boundary_conditions={name: Dirichlet(zero_velocity) for name in tags},
        boundary_tags=tags,
        viscosity=1.0,
        forcing=zero_velocity,
        initial_velocity=lambda _x: ufl.as_vector((0.0, 0.0)),
    )


def l2_norm(expression: Expr, mesh: MeshGeometry) -> float:
    return math.sqrt(fd.assemble(ufl.inner(expression, expression) * fd.dx(domain=mesh, degree=6)))


@pytest.mark.parametrize(
    ("build", "cells", "lower", "upper", "names"),
    [
        (rectangle, (3, 4), LOWER_2D, UPPER_2D, ["x_min", "x_max", "y_min", "y_max"]),
        (
            box,
            (2, 3, 2),
            LOWER_3D,
            UPPER_3D,
            ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"],
        ),
    ],
    ids=["rectangle", "box"],
)
def test_builtin_boundary_names_lie_on_their_coordinate_planes(
    build: Callable[..., tuple[MeshGeometry, dict[str, int]]],
    cells: tuple[int, ...],
    lower: tuple[float, ...],
    upper: tuple[float, ...],
    names: list[str],
) -> None:
    mesh, tags = build(cells, lower, upper)
    assert sorted(tags) == sorted(names)
    total = 0.0
    for name, tag in tags.items():
        axis, value = side_of(name, lower, upper)
        length, mean, spread = mean_and_spread(mesh, tag, axis)
        assert mean == pytest.approx(value, abs=1e-12)
        assert spread == pytest.approx(0.0, abs=1e-24)
        total += length
    assert total == pytest.approx(fd.assemble(fd.Constant(1.0) * fd.ds(domain=mesh)))


def test_periodic_rectangle_has_no_boundary_and_the_requested_extent() -> None:
    mesh, tags = periodic_rectangle((4, 4), LOWER_2D, UPPER_2D)
    assert tags == {}
    x = fd.SpatialCoordinate(mesh)
    area = fd.assemble(fd.Constant(1.0) * fd.dx(domain=mesh))
    assert area == pytest.approx(1.5 * 1.75)
    assert fd.assemble(fd.Constant(1.0) * fd.ds(domain=mesh)) == 0.0
    centre = [fd.assemble(x[i] * fd.dx) / area for i in range(2)]
    assert centre == pytest.approx([0.25, 0.625])


def test_gmsh_mesh_loads_facet_physical_names() -> None:
    mesh, tags = gmsh_mesh(SQUARE_NAMED)
    assert tags == {"left": 11, "right": 12, "bottom": 13, "top": 14}
    sides = {"left": (0, 0.0), "right": (0, 1.0), "bottom": (1, 0.0), "top": (1, 1.0)}
    for name, (axis, value) in sides.items():
        length, mean, spread = mean_and_spread(mesh, tags[name], axis)
        assert (length, mean) == pytest.approx((1.0, value))
        assert spread == pytest.approx(0.0, abs=1e-24)


@pytest.mark.parametrize("source", ["builtin", "gmsh"])
def test_dirichlet_bc_by_name_changes_only_that_boundary(source: str) -> None:
    if source == "builtin":
        mesh, tags = rectangle((3, 4), LOWER_2D, UPPER_2D)
        name, axis, value = "x_max", 0, UPPER_2D[0]
    else:
        mesh, tags = gmsh_mesh(SQUARE_NAMED)
        name, axis, value = "top", 1, 1.0
    space = fd.VectorFunctionSpace(mesh, "CG", 2)
    x = fd.SpatialCoordinate(mesh)
    u = fd.Function(space)
    dirichlet_bc(space, name, ufl.as_vector((1 + x[1], x[0] - 2)), tags).apply(u)
    nodes = fd.Function(space).interpolate(x).dat.data_ro
    on = np.isclose(nodes[:, axis], value)
    assert on.any() and (~on).any()
    expected = np.column_stack((1 + nodes[on, 1], nodes[on, 0] - 2))
    np.testing.assert_allclose(u.dat.data_ro[on], expected, atol=1e-13)
    assert not u.dat.data_ro[~on].any()


def test_unknown_boundary_name_lists_the_available_ones() -> None:
    mesh, tags = rectangle((2, 2), (0.0, 0.0), (1.0, 1.0))
    space = fd.VectorFunctionSpace(mesh, "CG", 1)
    with pytest.raises(ValueError, match="'outlet' does not exist; available: x_max, x_min, y_max"):
        dirichlet_bc(space, "outlet", ufl.as_vector((0.0, 0.0)), tags)


@pytest.mark.parametrize("case", sorted(CASES))
def test_each_case_builds_on_a_coarse_mesh(case: str) -> None:
    problem = build_problem(Config(case=case, run_name="test", mesh=MeshOptions(cells=(4, 4))))
    assert problem.mesh.num_cells() == 32
    assert problem.exact_velocity is not None and problem.exact_pressure is not None
    assert set(problem.boundary_conditions) == set(problem.boundary_tags)
    assert uses_pressure_null_space(problem)
    x = fd.SpatialCoordinate(problem.mesh)
    initial_gap = problem.initial_velocity(x) - problem.exact_velocity(x, fd.Constant(0.0))
    assert l2_norm(initial_gap, problem.mesh) == 0.0
    assert math.isfinite(l2_norm(problem.forcing(x, fd.Constant(0.4)), problem.mesh))


@pytest.mark.parametrize("case", ["taylor_green_2d", "kovasznay"])
def test_derived_forcing_vanishes_for_exact_navier_stokes_solutions(case: str) -> None:
    problem = build_problem(Config(case=case, run_name="test", mesh=MeshOptions(cells=(8, 8))))
    x, t = fd.SpatialCoordinate(problem.mesh), fd.Constant(0.7)
    scale = l2_norm(ufl.grad(problem.exact_pressure(x, t)), problem.mesh)
    assert l2_norm(problem.forcing(x, t), problem.mesh) < 1e-12 * scale


def test_derived_forcing_matches_a_hand_derivation() -> None:
    problem = build_problem(Config(case="manufactured_temporal", run_name="test"))
    x, t = fd.SpatialCoordinate(problem.mesh), fd.Constant(0.7)
    by_hand = -ufl.sin(t) * ufl.as_vector((x[1], x[0])) + ufl.cos(t) ** 2 * ufl.as_vector(x)
    assert l2_norm(problem.forcing(x, t) - by_hand, problem.mesh) < 1e-14


def test_case_defaults_and_overrides() -> None:
    default = build_problem(Config(case="manufactured_temporal", run_name="test"))
    assert default.mesh.num_cells() == 2 * 4 * 4
    assert default.viscosity == 0.01
    config = Config(case="kovasznay", run_name="test", case_parameters={"reynolds": 20})
    assert build_problem(config).viscosity == pytest.approx(0.05)


@pytest.mark.parametrize(
    ("config", "message"),
    [
        (Config(case="cavity", run_name="test"), "case = 'cavity': unknown; available: kovasznay"),
        (
            Config(case="kovasznay", run_name="test", mesh=MeshOptions(source="gmsh")),
            "case 'kovasznay' supports only builtin meshes",
        ),
        (
            Config(case="kovasznay", run_name="test", mesh=MeshOptions(cells=(4, 4, 4))),
            r"mesh.cells = \[4, 4, 4\]: case 'kovasznay' needs 2 entries",
        ),
        (
            Config(case="kovasznay", run_name="test", case_parameters={"viscosity": 0.1}),
            "case_parameters: viscosity not used by case 'kovasznay'; it accepts reynolds",
        ),
        (
            Config(case="kovasznay", run_name="test", case_parameters={"reynolds": -5}),
            "case_parameters.reynolds = -5: must be a positive number",
        ),
    ],
    ids=["unknown-case", "gmsh-source", "cell-count", "unknown-parameter", "negative-parameter"],
)
def test_invalid_case_setup_is_rejected(config: Config, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        build_problem(config)


def test_rule3_is_checked_when_the_problem_is_built(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(CASES, "no_exact_solution", lambda _config: dirichlet_square())
    exact_history = TimeOptions(startup="exact_history")
    with pytest.raises(ValueError, match="time.startup = 'exact_history' requires an exact"):
        build_problem(Config(case="no_exact_solution", run_name="test", time=exact_history))
    build_problem(Config(case="no_exact_solution", run_name="test"))
    build_problem(Config(case="manufactured_temporal", run_name="test", time=exact_history))


def test_problem_check_rejects_unknown_and_unassigned_boundaries() -> None:
    problem = dirichlet_square()
    check_problem(problem)
    conditions = dict(problem.boundary_conditions)
    with pytest.raises(ValueError, match="'outlet' does not exist; available: x_max, x_min"):
        check_problem(replace(problem, boundary_conditions={**conditions, "outlet": Outflow()}))
    del conditions["y_max"]
    with pytest.raises(ValueError, match="boundaries without a condition: y_max"):
        check_problem(replace(problem, boundary_conditions=conditions))


def test_problem_check_rejects_a_divergent_exact_velocity() -> None:
    divergent = replace(dirichlet_square(), exact_velocity=lambda x, _t: ufl.as_vector((x[0], 0)))
    with pytest.raises(ValueError, match="exact velocity is not divergence-free at t = 0.0"):
        check_problem(divergent)


def test_pressure_null_space_exactly_without_outflow() -> None:
    problem = dirichlet_square()
    assert uses_pressure_null_space(problem)
    with_outflow = {**problem.boundary_conditions, "x_max": Outflow()}
    assert not uses_pressure_null_space(replace(problem, boundary_conditions=with_outflow))
