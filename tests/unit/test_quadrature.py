"""The single quadrature rule, its quadrature element and its weights (plan U2)."""

import itertools
import math
from types import SimpleNamespace

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.discretization import quadrature
from navier_stokes_vms.discretization.function_spaces import subscale_space
from navier_stokes_vms.discretization.quadrature import (
    boundary_measure,
    check_positive_weights,
    quadrature_degree,
    quadrature_element,
    reference_rule,
    volume_measure,
)

CELLS = {"triangle": ufl.triangle, "tetrahedron": ufl.tetrahedron}


def exact_simplex_integral(exponents: tuple[int, ...]) -> float:
    return math.prod(math.factorial(a) for a in exponents) / math.factorial(
        sum(exponents) + len(exponents)
    )


@pytest.mark.parametrize("degree", [3, 6])
@pytest.mark.parametrize("cell", ["triangle", "tetrahedron"])
def test_rule_has_positive_weights_and_integrates_monomials_exactly(cell: str, degree: int) -> None:
    points, weights = reference_rule(CELLS[cell], degree)
    assert (weights > 0).all()
    dimension = points.shape[1]
    for exponents in itertools.product(range(degree + 1), repeat=dimension):
        if sum(exponents) > degree:
            continue
        computed = weights @ np.prod(points**exponents, axis=1)
        assert computed == pytest.approx(exact_simplex_integral(exponents), rel=1e-13, abs=1e-16)


@pytest.mark.parametrize("degree", [3, 6])
@pytest.mark.parametrize("dimension", [2, 3])
def test_quadrature_space_samples_the_dx_rule(dimension: int, degree: int) -> None:
    mesh = fd.UnitSquareMesh(2, 2) if dimension == 2 else fd.UnitCubeMesh(1, 1, 1)
    x = fd.SpatialCoordinate(mesh)
    f = ufl.exp(x[0]) * ufl.sin(5 * x[1]) / (1 + x[0] * x[dimension - 1])
    space = fd.FunctionSpace(mesh, quadrature_element(mesh.ufl_cell(), degree))
    sampled = fd.Function(space).interpolate(f)
    dx = volume_measure(mesh, degree)
    direct = fd.assemble(f * dx)
    assert fd.assemble(sampled * dx) == pytest.approx(direct, rel=1e-14, abs=1e-15)
    other_rule = fd.assemble(f * fd.dx(domain=mesh, degree=degree + 4))
    assert abs(other_rule - direct) > 1e-9


@pytest.mark.parametrize("velocity_degree", [1, 2])
@pytest.mark.parametrize("dimension", [2, 3])
def test_production_subscale_space_samples_the_production_dx(
    dimension: int, velocity_degree: int
) -> None:
    mesh = fd.UnitSquareMesh(2, 2) if dimension == 2 else fd.UnitCubeMesh(1, 1, 1)
    x = fd.SpatialCoordinate(mesh)
    q = quadrature_degree(velocity_degree)
    stored = ufl.as_vector(
        [ufl.exp(x[i]) * ufl.sin(3 * x[(i + 1) % dimension] + i) for i in range(dimension)]
    )
    weight = ufl.as_vector([1 / (1 + x[i] ** 2 + x[(i + 1) % dimension]) for i in range(dimension)])
    sampled = fd.Function(subscale_space(mesh, q)).interpolate(stored)
    dx = volume_measure(mesh, q)
    direct = fd.assemble(ufl.inner(stored, weight) * dx)
    assert fd.assemble(ufl.inner(sampled, weight) * dx) == pytest.approx(direct, rel=1e-14)
    other_rule = fd.assemble(ufl.inner(stored, weight) * fd.dx(domain=mesh, degree=q + 4))
    assert abs(other_rule - direct) > 1e-9 * abs(direct)


def test_measures_carry_an_explicit_degree_and_scheme() -> None:
    mesh = fd.UnitSquareMesh(2, 2)
    expected = {"quadrature_degree": 6, "quadrature_rule": quadrature.SCHEME}
    assert volume_measure(mesh, 6).metadata() == expected
    assert boundary_measure(mesh, 6).metadata() == expected
    x = fd.SpatialCoordinate(mesh)
    assert fd.assemble(x[0] ** 6 * boundary_measure(mesh, 6)) == pytest.approx(1 + 2 / 7)


def test_quadrature_degree_is_three_times_the_velocity_degree() -> None:
    assert [quadrature_degree(k) for k in (1, 2, 3)] == [3, 6, 9]


def test_a_rule_with_a_negative_weight_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    centroid, third, fifth = [1 / 3, 1 / 3], 3 / 5, 1 / 5
    points = [centroid, [fifth, fifth], [third, fifth], [fifth, third]]
    weights = np.array([-27, 25, 25, 25]) / 96
    strang_fix = SimpleNamespace(point_set=SimpleNamespace(points=points), weights=weights)
    monkeypatch.setattr(quadrature, "make_quadrature", lambda *_args, **_kwargs: strang_fix)
    with pytest.raises(ValueError, match="non-positive weight -2.812e-01"):
        check_positive_weights(ufl.triangle, 3)
