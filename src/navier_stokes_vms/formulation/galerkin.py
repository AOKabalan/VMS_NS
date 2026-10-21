"""Galerkin residual: the master FE equation with every stabilization term removed."""

from dataclasses import dataclass

import firedrake as fd
import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.formulation.convective_forms import convection_form


@dataclass(frozen=True)
class SolveLevel:
    """Data at the solve level: α₀, δt, the velocity history g_hist and the forcing f."""

    alpha0: fd.Constant
    dt: fd.Constant
    history: fd.Function
    forcing: Expr


def time_derivative(g: Expr, history: Expr, alpha0: Expr, dt: Expr) -> Expr:
    """D_t g = (α₀ g − g_hist)/δt (spec: Operators and parameters)."""
    return (alpha0 * g - history) / dt


def galerkin_residual(
    u: Expr,
    p: Expr,
    v: Expr,
    q: Expr,
    a: Expr,
    u_h: Expr,
    nu: float,
    level: SolveLevel,
    dx: Measure,
    outflow_ds: Measure | None,
) -> Form:
    """(D_t u, v) + ν(∇u, ∇v) + c(a; u, v) − (p, ∇·v) + (q, ∇·u) − (f, v) (spec: Galerkin)."""
    return (
        ufl.inner(time_derivative(u, level.history, level.alpha0, level.dt), v) * dx
        + nu * ufl.inner(ufl.grad(u), ufl.grad(v)) * dx
        + convection_form(a, u_h, u, v, dx, outflow_ds)
        - p * ufl.div(v) * dx
        + q * ufl.div(u) * dx
        - ufl.inner(level.forcing, v) * dx
    )
