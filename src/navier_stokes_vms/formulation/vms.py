"""VMS terms of the master FE equation and the subscale expressions (spec: Master formulation)."""

from dataclasses import dataclass

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.config import StabilizationOptions
from navier_stokes_vms.formulation.convective_forms import convective_subscale_coupling
from navier_stokes_vms.formulation.galerkin import SolveLevel, time_derivative
from navier_stokes_vms.formulation.operators import momentum_residual, viscous_operator


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


def subscale_residual(
    u: Expr,
    p: Expr,
    a: Expr,
    nu: float,
    level: SolveLevel,
    chi: Switches,
    u_sub_history: Expr | None,
    convective_form: str = "skew_symmetric_1",
    viscous_form: str = "laplacian",
) -> Expr:
    """X = R + χ_D ũ_hist/δt, with R = R_u for ASGS and R̂_u = f − L_a(u, p) for OSS (χ_A = 0)."""
    if chi.is_dynamic != (u_sub_history is not None):
        raise ValueError("ũ_hist must be given exactly for dynamic subscales (χ_D = 1)")
    dt_u = (
        time_derivative(u, level.history, level.alpha0, level.dt)
        if chi.is_asgs
        else ufl.zero(u.ufl_shape)
    )
    residual = momentum_residual(a, u, p, nu, dt_u, level.forcing, convective_form, viscous_form)
    if chi.is_dynamic:
        residual = residual + u_sub_history / level.dt
    return residual


def velocity_subscale(
    u: Expr,
    p: Expr,
    a: Expr,
    nu: float,
    level: SolveLevel,
    tau_t: Expr,
    chi: Switches,
    u_sub_history: Expr | None,
    xi: Expr | None = None,
    convective_form: str = "skew_symmetric_1",
    viscous_form: str = "laplacian",
) -> Expr:
    """ũ = τ_t (X − ξ_h), with ξ_h = 0 for ASGS (spec: Master formulation)."""
    if chi.is_asgs != (xi is None):
        raise ValueError("ξ_h must be given exactly for orthogonal subscales (χ_A = 0)")
    residual = subscale_residual(
        u, p, a, nu, level, chi, u_sub_history, convective_form, viscous_form
    )
    if xi is not None:
        residual = residual - xi
    return tau_t * residual


def accelerated_subscale(
    u: Expr,
    p: Expr,
    nu: float,
    level: SolveLevel,
    tau_t: Expr,
    chi: Switches,
    u_sub_history: Expr | None,
    xi: Expr | None = None,
    convective_form: str = "skew_symmetric_1",
    viscous_form: str = "laplacian",
) -> Expr:
    """ũ solving (I + τ_t ∇u) ũ = τ_t (X with a = u, − ξ_h) pointwise (spec step 3)."""
    right_side = velocity_subscale(
        u, p, u, nu, level, tau_t, chi, u_sub_history, xi, convective_form, viscous_form
    )
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
    xi: Expr | None = None,
    zeta: Expr | None = None,
    convective_form: str = "skew_symmetric_1",
    viscous_form: str = "laplacian",
) -> Form:
    """χ_A χ_D (D_t ũ, v) + (ũ, L*_a(v, q)) + (τ_c (∇·u − ζ_h), ∇·v) with ũ substituted.

    The convective part of (ũ, L*_a(v, q)) is the selected coupling, −(T̃, ∇ˢv) for EMAC–VMS
    (spec: EMAC–VMS). ζ_h is None for ASGS and for OSS without the pressure subscale, where τ_c = 0.
    """
    if chi.is_asgs and zeta is not None:
        raise ValueError("ζ_h is only defined for orthogonal subscales (χ_A = 0)")
    u_sub = velocity_subscale(
        u, p, a, nu, level, tau_t, chi, u_sub_history, xi, convective_form, viscous_form
    )
    divergence = ufl.div(u) if zeta is None else ufl.div(u) - zeta
    form = convective_subscale_coupling(a, u, v, u_sub, dx, convective_form, chi.is_nonlinear)
    form += (
        -ufl.inner(u_sub, viscous_operator(v, nu, viscous_form) + ufl.grad(q)) * dx
        + tau_c * divergence * ufl.div(v) * dx
    )
    if chi.is_asgs and chi.is_dynamic:
        dt_u_sub = time_derivative(u_sub, u_sub_history, level.alpha0, level.dt)
        form += ufl.inner(dt_u_sub, v) * dx
    return form
