"""Static linear ASGS equals the independently expanded classical residual (plan U4)."""

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.discretization.function_spaces import mixed_space
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.formulation.stabilization_parameters import stabilization_parameters
from navier_stokes_vms.formulation.vms import Switches, vms_terms

STATIC_LINEAR = Switches(is_asgs=True, is_dynamic=False, is_nonlinear=False)


@pytest.mark.parametrize("k", [1, 2])
@pytest.mark.parametrize("dimension", [2, 3])
@pytest.mark.parametrize("pressure_subscale", [False, True])
def test_asgs_equals_independent_classical_form(
    k: int, dimension: int, pressure_subscale: bool
) -> None:
    mesh = fd.UnitSquareMesh(3, 3) if dimension == 2 else fd.UnitCubeMesh(2, 2, 2)
    space = mixed_space(mesh, k, k)
    state = fd.Function(space)
    u, p = state.subfunctions
    history = fd.Function(u.function_space())
    force = fd.Function(u.function_space())
    rng = np.random.default_rng(20261022)
    for field in (u, p, history, force):
        field.dat.data[:] = rng.standard_normal(field.dat.data.shape)
    v, q = fd.TestFunctions(space)
    dx = volume_measure(mesh, 3 * k)
    nu, alpha0, dt = 0.13, 1.5, 0.07
    level = SolveLevel(fd.Constant(alpha0), fd.Constant(dt), history, force)
    tau_m, tau_c = stabilization_parameters(u, nu, k, 4.0, 2.0)
    actual = galerkin_residual(u, p, v, q, u, u, nu, level, dx, None) + vms_terms(
        u,
        p,
        v,
        q,
        u,
        nu,
        level,
        tau_m,
        tau_c if pressure_subscale else ufl.zero(),
        dx,
        STATIC_LINEAR,
        None,
    )

    h = ufl.MinCellEdgeLength(mesh)
    tau = 1 / (4 * k**4 * nu / h**2 + 2 * k * ufl.sqrt(ufl.dot(u, u)) / h)
    dt_u = (alpha0 * u - history) / dt
    convection = ufl.dot(ufl.grad(u), u)
    test_convection = ufl.dot(ufl.grad(v), u)
    residual = dt_u + convection - nu * ufl.div(ufl.grad(u)) + ufl.grad(p) - force
    test_operator = test_convection + nu * ufl.div(ufl.grad(v)) + ufl.grad(q)
    reference = (
        ufl.inner(dt_u - force, v)
        + nu * ufl.inner(ufl.grad(u), ufl.grad(v))
        + (ufl.inner(convection, v) - ufl.inner(test_convection, u)) / 2
        - p * ufl.div(v)
        + q * ufl.div(u)
        + tau * ufl.inner(residual, test_operator)
    ) * dx
    if pressure_subscale:
        reference += h**2 / (4 * k**4 * tau) * ufl.div(u) * ufl.div(v) * dx
    assembled, expected = fd.assemble(actual), fd.assemble(reference)
    with assembled.dat.vec_ro as vector, expected.dat.vec_ro as ref:
        difference = vector.copy()
        difference.axpy(-1, ref)
        assert difference.norm() <= 1e-12 * ref.norm()
    if k == 2:
        assert fd.assemble(ufl.inner(ufl.div(ufl.grad(u)), ufl.div(ufl.grad(u))) * dx) > 1
