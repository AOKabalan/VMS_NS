"""Scalar Galerkin transport kernels (Boussinesq formulation B14–B15)."""

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.formulation.galerkin import SolveLevel, time_derivative


def scalar_convection(a: Expr, field: Expr, test: Expr, dx: Measure, form: str) -> Form:
    """Closed-domain scalar transport in any of the three forms (B15)."""
    forward = ufl.dot(a, ufl.grad(field)) * test
    if form == "conservative":
        return -field * ufl.dot(a, ufl.grad(test)) * dx
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


def scalar_residual(
    field: Expr,
    a: Expr,
    diffusivity: float,
    level: SolveLevel,
    is_asgs: bool,
    subscale_history: Expr | None = None,
) -> Expr:
    """Scalar residual with OSS time cancellation and optional dynamic history (B10)."""
    residual = level.forcing - ufl.dot(a, ufl.grad(field)) + diffusivity * ufl.div(ufl.grad(field))
    if is_asgs:
        residual -= time_derivative(field, level.history, level.alpha0, level.dt)
    if subscale_history is not None:
        residual += subscale_history / level.dt
    return residual


def scalar_subscale(residual: Expr, tau: Expr, projection: Expr | None = None) -> Expr:
    """Eliminated scalar subscale τ_t(X − η), with no spatial differentiation (B10)."""
    return tau * (residual if projection is None else residual - projection)


def scalar_vms_terms(
    subscale: Expr,
    test: Expr,
    a: Expr,
    diffusivity: float,
    dx: Measure,
    subscale_rate: Expr | None = None,
) -> Form:
    """Scalar subscale feedback through the convective and diffusive adjoint (B14)."""
    form = -subscale * (ufl.dot(a, ufl.grad(test)) + diffusivity * ufl.div(ufl.grad(test))) * dx
    if subscale_rate is not None:
        form += subscale_rate * test * dx
    return form
