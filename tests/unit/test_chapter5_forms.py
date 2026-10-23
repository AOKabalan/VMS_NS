"""Chapter 5 algebra: index conventions, both defect tables, minimal closure and traction."""

from collections.abc import Iterator

import firedrake as fd
import numpy as np
import pytest
import ufl
from firedrake.mesh import MeshGeometry
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.discretization.function_spaces import subscale_space
from navier_stokes_vms.discretization.quadrature import boundary_measure, volume_measure
from navier_stokes_vms.formulation.convective_forms import (
    PRESSURE_UNKNOWNS,
    convection_form,
    convective_subscale_coupling,
    divergence_pairing,
    minimal_stress,
    outflow_term,
    residual_convection,
    resolved_convection,
)
from navier_stokes_vms.formulation.operators import (
    momentum_residual,
    viscous_bilinear,
    viscous_operator,
)

FORMS = tuple(PRESSURE_UNKNOWNS)


@pytest.fixture(scope="module", params=[(2, 1), (2, 2), (3, 1), (3, 2)])
def fields(request: pytest.FixtureRequest) -> tuple[fd.Function, fd.Function, fd.Function, Measure]:
    dimension, degree = request.param
    mesh = fd.UnitSquareMesh(3, 3) if dimension == 2 else fd.UnitCubeMesh(2, 2, 2)
    space = fd.VectorFunctionSpace(mesh, "CG", degree)
    u, v = fd.Function(space), fd.Function(space)
    s = fd.Function(subscale_space(mesh, 3 * degree))
    rng = np.random.default_rng(12000 + 10 * dimension + degree)
    for field in (u, v, s):
        field.dat.data[:] = rng.standard_normal(field.dat.data.shape)
    fd.DirichletBC(space, 0, "on_boundary").apply(u)
    return u, v, s, volume_measure(mesh, 3 * degree)


def value(form: Form) -> float:
    expanded = ufl.algorithms.expand_derivatives(form)
    return fd.assemble(expanded) if any(i.integrand() != 0 for i in expanded.integrals()) else 0.0


def same(left: Form, right: Form) -> None:
    actual, expected = value(left), value(right)
    assert actual == pytest.approx(expected, rel=3e-11, abs=3e-11)


def rigid_tests(mesh: MeshGeometry) -> Iterator[tuple[int, Expr, Expr]]:
    dimension = mesh.geometric_dimension
    x = fd.SpatialCoordinate(mesh)
    for k in range(dimension):
        e = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1)).assign(
            fd.Constant(tuple(float(i == k) for i in range(dimension)))
        )
        phi = ufl.as_vector((-x[1], x[0])) if dimension == 2 else ufl.cross(e, x)
        yield k, e, phi


@pytest.mark.parametrize("nonlinear", [False, True])
@pytest.mark.parametrize("form", FORMS)
def test_forms_against_component_indices(
    fields: tuple[fd.Function, fd.Function, fd.Function, Measure], form: str, nonlinear: bool
) -> None:
    u, v, s, dx = fields
    dimension = u.ufl_shape[0]
    a = u + s if nonlinear else u
    sub = s if nonlinear else ufl.zero(u.ufl_shape)
    nc = sum(a[i] * u[j].dx(i) * v[j] for i in range(dimension) for j in range(dimension))
    conservative = -sum(
        u[j] * a[i] * v[j].dx(i) for i in range(dimension) for j in range(dimension)
    )
    uv = sum(u[j] * v[j] for j in range(dimension))
    divergence = sum(u[i].dx(i) * uv - sub[i] * uv.dx(i) for i in range(dimension))
    expected = {
        "nonconservative": nc,
        "conservative": conservative,
        "skew_symmetric_1": (nc + conservative) / 2,
        "skew_symmetric_2": nc + divergence / 2,
        "emac": sum(
            (u[j].dx(i) + u[i].dx(j)) * a[i] * v[j]
            for i in range(dimension)
            for j in range(dimension)
        )
        + divergence,
        "emac_vms": sum(
            (u[j].dx(i) + u[i].dx(j)) * u[i] * v[j]
            for i in range(dimension)
            for j in range(dimension)
        )
        + sum(u[i].dx(i) * uv for i in range(dimension)),
    }
    same(resolved_convection(a, u, u, v, dx, form), expected[form] * dx)
    expected_residual = ufl.as_vector(
        [
            sum(a[i] * u[j].dx(i) for i in range(dimension))
            + (
                sum(u[i] * u[i].dx(j) for i in range(dimension))
                if form in ("emac", "emac_vms")
                else 0
            )
            for j in range(dimension)
        ]
    )
    same(
        ufl.inner(residual_convection(a, u, form), v) * dx,
        ufl.inner(expected_residual, v) * dx,
    )
    same(ufl.inner(ufl.outer(u, a), ufl.grad(v)) * dx, -conservative * dx)


@pytest.mark.parametrize("nonlinear", [False, True])
def test_minimal_closure_T1_T2_T3_T4(
    fields: tuple[fd.Function, fd.Function, fd.Function, Measure], nonlinear: bool
) -> None:
    u, v, s, dx = fields
    a = u + s if nonlinear else u
    stress = minimal_stress(u, s, nonlinear)
    coupling = convective_subscale_coupling(a, u, v, s, dx, "emac_vms", nonlinear)
    expansion = -ufl.inner(s, ufl.dot(ufl.grad(v), a) + ufl.dot(u, ufl.grad(v))) * dx
    same(coupling, expansion)
    # Check T3 pointwise, not merely after cancellation of quadrature sums.
    contraction = ufl.inner(stress, ufl.sym(ufl.grad(u)))
    transfer = ufl.inner(s, residual_convection(a, u, "emac_vms"))
    scalar_space = fd.FunctionSpace(
        u.function_space().mesh(), s.function_space().ufl_element().sub_elements[0]
    )
    error = fd.Function(scalar_space).interpolate(contraction - transfer)
    scale = fd.Function(scalar_space).interpolate(abs(contraction) + abs(transfer))
    assert np.max(np.abs(error.dat.data_ro)) < 1e-12 * max(1, np.max(scale.dat.data_ro))
    zero = ufl.zero(u.ufl_shape)
    same(
        resolved_convection(u, u, u, v, dx, "emac_vms")
        + convective_subscale_coupling(u, u, v, zero, dx, "emac_vms", nonlinear),
        resolved_convection(u, u, u, v, dx, "emac"),
    )
    # T2 is checked with all forms in the next test, with arbitrary non-bubble s.


@pytest.mark.parametrize("nonlinear", [False, True])
@pytest.mark.parametrize("form", FORMS)
def test_every_entry_of_both_defect_tables(
    fields: tuple[fd.Function, fd.Function, fd.Function, Measure], form: str, nonlinear: bool
) -> None:
    u, _, s, dx = fields
    mesh = u.function_space().mesh()
    a = u + s if nonlinear else u
    sub = s if nonlinear else ufl.zero(u.ufl_shape)
    factor = {
        "nonconservative": -1,
        "conservative": 0,
        "skew_symmetric_1": -0.5,
        "skew_symmetric_2": -0.5,
    }.get(form, 0)
    for k, e, phi in rigid_tests(mesh):
        momentum = resolved_convection(a, u, u, e, dx, form)
        expected = factor * divergence_pairing(u, sub, u[k], dx)
        if form == "emac" and nonlinear:
            expected = ufl.inner(s, u.dx(k)) * dx
        same(momentum, expected)
        angular = resolved_convection(a, u, u, phi, dx, form)
        angular += convective_subscale_coupling(a, u, phi, s, dx, form, nonlinear)
        expected = factor * divergence_pairing(u, sub, ufl.inner(u, phi), dx)
        if form == "emac" and nonlinear:
            expected = ufl.inner(s, ufl.dot(ufl.grad(u), phi)) * dx
        if not nonlinear and form != "emac_vms":
            expected -= ufl.inner(s, ufl.dot(ufl.grad(phi), u)) * dx
        same(angular, expected)
    energy = resolved_convection(a, u, u, u, dx, form)
    energy += convective_subscale_coupling(a, u, u, s, dx, form, nonlinear)
    energy += ufl.inner(s, residual_convection(a, u, form)) * dx
    energy_factor = {"nonconservative": -0.5, "conservative": 0.5}.get(form, 0)
    expected = energy_factor * divergence_pairing(u, sub, ufl.inner(u, u), dx)
    if form == "emac":
        transfer = ufl.dot(ufl.grad(u), u) if nonlinear else ufl.dot(u, ufl.grad(u))
        expected = ufl.inner(s, transfer) * dx
    same(energy, expected)


def test_periodic_angular_counterexample_T6() -> None:
    mesh = fd.PeriodicUnitSquareMesh(3, 3)
    x, y = fd.SpatialCoordinate(mesh)
    u = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 2)).interpolate(
        ufl.as_vector((y * (1 - y), 1))
    )
    phi = ufl.as_vector((-y, x))
    dx = volume_measure(mesh, 6)
    assert value(resolved_convection(u, u, u, phi, dx, "nonconservative")) == pytest.approx(1 / 6)
    assert value(ufl.div(u) ** 2 * dx) < 1e-26


@pytest.mark.parametrize("form", FORMS)
def test_outflow_correction_and_natural_pressure(form: str) -> None:
    mesh = fd.UnitSquareMesh(2, 2)
    x, y = fd.SpatialCoordinate(mesh)
    u = ufl.as_vector((1 + x + y**2, x * y))
    # v vanishes on every face except the outlet x=1.
    v = x * y * (1 - y) * ufl.as_vector((1 + y, x))
    dx, ds = volume_measure(mesh, 9), boundary_measure(mesh, 9)(2)
    nc = ufl.inner(ufl.dot(ufl.grad(u), u), v) * dx
    correction = ufl.div(u) * ufl.inner(u, v) * dx
    alpha = {
        "nonconservative": 0,
        "conservative": 1,
        "skew_symmetric_1": 0.5,
        "skew_symmetric_2": 0.5,
    }.get(form)
    if alpha is not None:
        same(convection_form(u, u, u, v, dx, ds, form), nc + alpha * correction)
        assert PRESSURE_UNKNOWNS[form] == "p"
    else:
        # No facet correction: the viscous/pressure traction is in P, not physical p.
        assert not outflow_term(u, u, v, ds, form).integrals()
        assert PRESSURE_UNKNOWNS[form] == "P"
        same(
            convection_form(u, u, u, v, dx, ds, form),
            nc + correction + ufl.inner(ufl.grad(ufl.inner(u, u) / 2), v) * dx,
        )
    nu = 0.3
    p = x + y
    n = fd.FacetNormal(mesh)
    for viscous in ("laplacian", "stress"):
        stress = nu * ufl.grad(u) if viscous == "laplacian" else 2 * nu * ufl.sym(ufl.grad(u))
        weak = viscous_bilinear(u, v, nu, dx, viscous) - p * ufl.div(v) * dx
        strong = ufl.inner(-viscous_operator(u, nu, viscous) + ufl.grad(p), v) * dx
        traction = ufl.inner(ufl.dot(stress, n) - p * n, v) * ds
        same(weak, strong + traction)


def test_stress_operator_residual_and_formal_self_adjointness() -> None:
    mesh = fd.UnitSquareMesh(1, 1)
    x, y = fd.SpatialCoordinate(mesh)
    u = ufl.as_vector((x**3 + x * y, y**3 - x * y))
    v = ufl.as_vector((y**2 + x * y, x**2 - y))
    p = x**2 - y
    nu = 0.7
    dx, ds = volume_measure(mesh, 9), boundary_measure(mesh, 9)
    n = fd.FacetNormal(mesh)
    lu, lv = viscous_operator(u, nu, "stress"), viscous_operator(v, nu, "stress")
    same(
        ufl.inner(lu, v) * dx, ufl.inner(nu * (ufl.div(ufl.grad(u)) + ufl.grad(ufl.div(u))), v) * dx
    )
    surface = (
        2
        * nu
        * (
            ufl.inner(ufl.dot(ufl.sym(ufl.grad(u)), n), v)
            - ufl.inner(u, ufl.dot(ufl.sym(ufl.grad(v)), n))
        )
        * ds
    )
    same((ufl.inner(lu, v) - ufl.inner(u, lv)) * dx, surface)
    zero = ufl.zero(u.ufl_shape)
    residual = momentum_residual(u, u, p, nu, zero, zero, "skew_symmetric_1", "stress")
    same(
        ufl.inner(residual, v) * dx, ufl.inner(-ufl.dot(ufl.grad(u), u) + lu - ufl.grad(p), v) * dx
    )
