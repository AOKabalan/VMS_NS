"""Linear-splitting ASGS terms of the master FE equation (spec: Master formulation)."""

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.formulation.galerkin import SolveLevel, time_derivative
from navier_stokes_vms.formulation.operators import adjoint_operator, momentum_residual


def velocity_subscale(
    u: Expr, p: Expr, a: Expr, nu: float, level: SolveLevel, tau_t: Expr, u_sub_history: Expr | None
) -> Expr:
    """ũ = τ_t (R_u + χ_D ũ_hist/δt); u_sub_history is None for static subscales (χ_D = 0)."""
    dt_u = time_derivative(u, level.history, level.alpha0, level.dt)
    residual = momentum_residual(a, u, p, nu, dt_u, level.forcing)
    if u_sub_history is None:
        return tau_t * residual
    return tau_t * (residual + u_sub_history / level.dt)


def vms_terms(
    u: Expr,
    p: Expr,
    v: Expr,
    q: Expr,
    a: Expr,
    nu: float,
    level: SolveLevel,
    tau_t: Expr,
    tau_c: Expr,
    dx: Measure,
    u_sub_history: Expr | None = None,
) -> Form:
    """χ_D (D_t ũ, v) + (ũ, L*_a(v, q)) + (τ_c ∇·u, ∇·v) with ũ substituted (spec: Master)."""
    u_sub = velocity_subscale(u, p, a, nu, level, tau_t, u_sub_history)
    form = (
        ufl.inner(u_sub, adjoint_operator(a, v, q, nu)) * dx + tau_c * ufl.div(u) * ufl.div(v) * dx
    )
    if u_sub_history is None:
        return form
    return form + ufl.inner(time_derivative(u_sub, u_sub_history, level.alpha0, level.dt), v) * dx
