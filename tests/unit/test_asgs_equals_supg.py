"""ASGS residuals equal independently expanded classical forms (plan U4)."""

import logging
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.discretization.function_spaces import mixed_space
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.formulation.stabilization_parameters import stabilization_parameters
from navier_stokes_vms.formulation.vms import Switches, vms_terms
from navier_stokes_vms.solver.time_stepper import TimeStepper

STATIC_LINEAR = Switches(is_asgs=True, is_dynamic=False, is_nonlinear=False)


def assert_same_vector(actual: ufl.Form, reference: ufl.Form) -> None:
    assembled, expected = fd.assemble(actual), fd.assemble(reference)
    with assembled.dat.vec_ro as vector, expected.dat.vec_ro as ref:
        difference = vector.copy()
        difference.axpy(-1, ref)
        assert difference.norm() <= 1e-12 * ref.norm()


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
    assert_same_vector(actual, reference)
    if k == 2:
        assert fd.assemble(ufl.inner(ufl.div(ufl.grad(u)), ufl.div(ufl.grad(u))) * dx) > 1


@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("k", [1, 2])
@pytest.mark.parametrize("scheme", ["bdf2", "crank_nicolson"])
def test_dynamic_asgs_step_equals_independent_classical_form(
    tmp_path: Path, scheme: str, k: int, splitting: str
) -> None:
    """The time stepper's residual at step 2, against coefficients written out here."""
    dt = 0.07
    config = Config(
        case="manufactured_spatial",
        run_name="u4",
        mesh=MeshOptions(cells=(3, 3)),
        discretization=DiscretizationOptions(velocity_degree=k, pressure_degree=k),
        time=TimeOptions(scheme=scheme, dt=dt, t_end=2 * dt),
        stabilization=StabilizationOptions(
            subscale_space="asgs", subscale_dynamics="dynamic", scale_splitting=splitting
        ),
        output=OutputOptions(vtk_every=0),
    )
    problem = build_problem(config)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger(__name__))
    rng = np.random.default_rng(7)
    velocities, subscales = stepper._history, stepper._subscales._levels
    for field in (*velocities, *subscales):
        field.dat.data[:] = rng.standard_normal(field.dat.data.shape)
    stepper._set_level(2, step_coefficients(scheme, "lower_order", 2))
    for field in stepper._solution.subfunctions:
        field.dat.data[:] = rng.standard_normal(field.dat.data.shape)
    actual = stepper._step_solver._inner._newton_residual

    mesh = problem.mesh
    x = fd.SpatialCoordinate(mesh)
    u, p = ufl.split(stepper._solution)
    v, q = fd.TestFunctions(stepper._solution.function_space())
    nu, dx = problem.viscosity, stepper._dx
    if scheme == "bdf2":
        alpha0 = 3 / 2
        u_hist = 2 * velocities[0] - velocities[1] / 2
        u_sub_hist = 2 * subscales[0] - subscales[1] / 2
        force = problem.forcing(x, fd.Constant(2 * dt))
    else:
        # Crank–Nicolson solves at t^(n+½): D_t g = (2 g^(n+½) − 2 gⁿ)/δt, data averaged.
        alpha0 = 2.0
        u_hist = 2 * velocities[0]
        u_sub_hist = 2 * subscales[0]
        force = (problem.forcing(x, fd.Constant(dt)) + problem.forcing(x, fd.Constant(2 * dt))) / 2
    # Nonlinear splitting starts the step from ũ* = ũⁿ (spec: Time-step algorithm step 1).
    a = u + subscales[0] if splitting == "nonlinear" else u
    h = ufl.MinCellEdgeLength(mesh)
    tau_m = 1 / (4 * k**4 * nu / h**2 + 2 * k * ufl.sqrt(ufl.dot(a, a)) / h)
    tau_t = 1 / (alpha0 / dt + 1 / tau_m)
    dt_u = (alpha0 * u - u_hist) / dt
    residual = force - dt_u - ufl.dot(ufl.grad(u), a) + nu * ufl.div(ufl.grad(u)) - ufl.grad(p)
    u_sub = tau_t * (residual + u_sub_hist / dt)
    test_operator = ufl.dot(ufl.grad(v), a) + nu * ufl.div(ufl.grad(v)) + ufl.grad(q)
    reference = (
        ufl.inner(dt_u - force, v)
        + nu * ufl.inner(ufl.grad(u), ufl.grad(v))
        + (ufl.inner(ufl.dot(ufl.grad(u), a), v) - ufl.inner(ufl.dot(ufl.grad(v), a), u)) / 2
        - p * ufl.div(v)
        + q * ufl.div(u)
        + ufl.inner((alpha0 * u_sub - u_sub_hist) / dt, v)
        - ufl.inner(u_sub, test_operator)
        + h**2 / (4 * k**4 * tau_m) * ufl.div(u) * ufl.div(v)
    ) * dx
    assert_same_vector(actual, reference)
