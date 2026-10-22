"""Subscale storage: history bookkeeping and the pointwise update of the iterate ũ*."""

import logging

import firedrake as fd
import numpy as np
import pytest
import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Startup, TimeScheme
from navier_stokes_vms.discretization.function_spaces import mixed_space, subscale_space
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.formulation.galerkin import SolveLevel
from navier_stokes_vms.formulation.stabilization_parameters import (
    dynamic_parameter,
    stabilization_parameters,
)
from navier_stokes_vms.formulation.vms import Switches, accelerated_subscale, velocity_subscale
from navier_stokes_vms.solver.subscale_storage import SubscaleIterate, SubscaleStorage

DT = 0.1
LOGGER = logging.getLogger(__name__)


def known(x: Expr, t: float) -> Expr:
    """g(t)(sin x, y) with g = t + t² for t > 0; zero for t ≤ 0, as ũ⁰ and exact histories."""
    return (t + t**2 if t > 1e-12 else 0.0) * ufl.as_vector((ufl.sin(x[0]), x[1]))


@pytest.mark.parametrize("startup", ["lower_order", "exact_history"])
@pytest.mark.parametrize("scheme", ["bdf1", "bdf2", "bdf3", "crank_nicolson"])
def test_histories_follow_the_startup_and_crank_nicolson_extrapolates(
    scheme: TimeScheme, startup: Startup
) -> None:
    mesh = fd.UnitSquareMesh(2, 2)
    space = subscale_space(mesh, 3)
    x = fd.SpatialCoordinate(mesh)
    levels = max(len(step_coefficients(scheme, startup, s).history_weights) for s in (1, 2, 3))
    storage = SubscaleStorage(space, levels)
    expected = fd.Function(space)
    for step in range(1, 7):
        coefficients = step_coefficients(scheme, startup, step)
        storage.form_history(coefficients)
        expected.interpolate(
            sum(
                weight * known(x, (step - 1 - j) * DT)
                for j, weight in enumerate(coefficients.history_weights)
            )
        )
        np.testing.assert_allclose(storage.history.dat.data_ro, expected.dat.data_ro, atol=1e-14)
        solve_level_value = coefficients.data_weight_new * known(
            x, step * DT
        ) + coefficients.data_weight_old * known(x, (step - 1) * DT)
        storage.store(solve_level_value, coefficients)


@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
@pytest.mark.parametrize("degree", [1, 2])
def test_fixed_point_and_accelerated_updates_reach_the_defining_subscale(
    degree: int, dynamics: str
) -> None:
    mesh = fd.UnitSquareMesh(4, 4)
    x = fd.SpatialCoordinate(mesh)
    state = fd.Function(mixed_space(mesh, degree, degree))
    velocity, pressure = state.subfunctions
    velocity.interpolate(ufl.as_vector((1 + ufl.sin(3 * x[0]) * x[1] / 4, 1 + x[0] * x[1] / 2)))
    pressure.interpolate(x[0] * x[1] - ufl.sin(x[1]))
    space = subscale_space(mesh, 3 * degree)
    history = fd.Function(velocity.function_space()).interpolate(1.5 * velocity)
    level = SolveLevel(fd.Constant(1.5), fd.Constant(0.1), history, ufl.as_vector((x[0], 1.0)))
    chi = Switches(is_asgs=True, is_dynamic=dynamics == "dynamic", is_nonlinear=True)
    u_sub_history = None
    if chi.is_dynamic:
        u_sub_history = fd.Function(space).interpolate(0.01 * ufl.as_vector((ufl.sin(x[0]), x[1])))
    u, p = ufl.split(state)
    nu = 0.01

    def tau_t(u_sub: fd.Function) -> Expr:
        tau_m, _ = stabilization_parameters(u + u_sub, nu, degree, 4.0, 2.0)
        return dynamic_parameter(tau_m, level.alpha0, level.dt) if chi.is_dynamic else tau_m

    def converged_subscale(update: str) -> tuple[fd.Function, int]:
        value = fd.Function(space)
        expression = (
            velocity_subscale(u, p, u + value, nu, level, tau_t(value), chi, u_sub_history)
            if update == "fixed_point"
            else accelerated_subscale(u, p, nu, level, tau_t(value), chi, u_sub_history)
        )
        result = SubscaleIterate(value, expression, 1e-13, 200).update()
        assert result.converged, update
        return value, result.iterations

    fixed_point, fixed_point_sweeps = converged_subscale("fixed_point")
    accelerated, accelerated_sweeps = converged_subscale("accelerated")
    LOGGER.info(
        "Local sweeps: fixed point %d, accelerated %d", fixed_point_sweeps, accelerated_sweeps
    )
    scale = np.abs(accelerated.dat.data_ro).max()
    assert scale > 1e-3
    assert np.abs(fixed_point.dat.data_ro - accelerated.dat.data_ro).max() <= 1e-10 * scale
    defining = fd.Function(space).interpolate(
        accelerated
        - velocity_subscale(
            u, p, u + accelerated, nu, level, tau_t(accelerated), chi, u_sub_history
        )
    )
    assert np.abs(defining.dat.data_ro).max() <= 1e-13 * scale
