"""Subscale history bookkeeping: shifting, startup sequence and Crank–Nicolson extrapolation."""

import firedrake as fd
import numpy as np
import pytest
import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Startup, TimeScheme
from navier_stokes_vms.discretization.function_spaces import subscale_space
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.solver.subscale_storage import SubscaleStorage

DT = 0.1


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
