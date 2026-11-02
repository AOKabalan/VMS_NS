"""Scalar Galerkin transport kernels (Boussinesq formulation B14–B15)."""

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.formulation.galerkin import SolveLevel, time_derivative


def scalar_convection(a: Expr, field: Expr, test: Expr, dx: Measure, form: str) -> Form:
    """Closed-domain nonconservative or skew-symmetric scalar transport (B15)."""
    forward = ufl.dot(a, ufl.grad(field)) * test
    if form == "nonconservative":
        return forward * dx
    if form == "skew_symmetric":
        return (forward - field * ufl.dot(a, ufl.grad(test))) / 2 * dx
    raise ValueError(f"unsupported scalar convection: {form!r}")


def scalar_galerkin_residual(
    field: Expr,
    test: Expr,
    a: Expr,
    diffusivity: float,
    level: SolveLevel,
    dx: Measure,
    convection: str,
) -> Form:
    """Scalar transient advection–diffusion residual before prescribed outward flux (B14)."""
    return (
        time_derivative(field, level.history, level.alpha0, level.dt) * test
        + diffusivity * ufl.dot(ufl.grad(field), ufl.grad(test))
        - level.forcing * test
    ) * dx + scalar_convection(a, field, test, dx, convection)
