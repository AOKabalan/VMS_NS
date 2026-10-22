"""VMS terms of the master FE equation and the subscale expressions (spec: Master formulation)."""

from dataclasses import dataclass

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.config import StabilizationOptions
from navier_stokes_vms.formulation.galerkin import SolveLevel, time_derivative
from navier_stokes_vms.formulation.operators import adjoint_operator, momentum_residual


@dataclass(frozen=True)
class Switches:
    """The spec's switches: χ_A (ASGS, not OSS), χ_D (dynamic) and nonlinear scale splitting."""

    is_asgs: bool
    is_dynamic: bool
    is_nonlinear: bool


def switches(stabilization: StabilizationOptions) -> Switches:
    """The switches of a configured VMS case."""
    return Switches(
        is_asgs=stabilization.subscale_space == "asgs",
        is_dynamic=stabilization.subscale_dynamics == "dynamic",
        is_nonlinear=stabilization.scale_splitting == "nonlinear",
    )


def velocity_subscale(
    u: Expr,
    p: Expr,
    a: Expr,
    nu: float,
    level: SolveLevel,
    tau_t: Expr,
    chi: Switches,
    u_sub_history: Expr | None,
) -> Expr:
    """ũ = τ_t (R_u + χ_D ũ_hist/δt) for ASGS (spec: Master formulation)."""
    if chi.is_dynamic != (u_sub_history is not None):
        raise ValueError("ũ_hist must be given exactly for dynamic subscales (χ_D = 1)")
    dt_u = time_derivative(u, level.history, level.alpha0, level.dt)
    residual = momentum_residual(a, u, p, nu, dt_u, level.forcing)
    if chi.is_dynamic:
        residual = residual + u_sub_history / level.dt
    return tau_t * residual


def accelerated_subscale(
    u: Expr,
    p: Expr,
    nu: float,
    level: SolveLevel,
    tau_t: Expr,
    chi: Switches,
    u_sub_history: Expr | None,
) -> Expr:
    """ũ solving (I + τ_t ∇u) ũ = τ_t (R_u with a = u, + χ_D ũ_hist/δt) pointwise (spec step 3)."""
    right_side = velocity_subscale(u, p, u, nu, level, tau_t, chi, u_sub_history)
    identity = ufl.Identity(u.ufl_shape[0])
    return ufl.dot(ufl.inv(identity + tau_t * ufl.grad(u)), right_side)


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
    chi: Switches,
    u_sub_history: Expr | None,
) -> Form:
    """χ_A χ_D (D_t ũ, v) + (ũ, L*_a(v, q)) + (τ_c ∇·u, ∇·v) with ũ substituted (spec: Master)."""
    u_sub = velocity_subscale(u, p, a, nu, level, tau_t, chi, u_sub_history)
    form = (
        ufl.inner(u_sub, adjoint_operator(a, v, q, nu)) * dx + tau_c * ufl.div(u) * ufl.div(v) * dx
    )
    if chi.is_asgs and chi.is_dynamic:
        dt_u_sub = time_derivative(u_sub, u_sub_history, level.alpha0, level.dt)
        form += ufl.inner(dt_u_sub, v) * dx
    return form
