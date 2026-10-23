"""Elementwise viscous and momentum operators and residual (spec: Operators and Chapter 5)."""

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.formulation.convective_forms import residual_convection


def elementwise_laplacian(u: Expr) -> Expr:
    """Δ_h u inside cell integrals (spec: Operators and parameters)."""
    return ufl.div(ufl.grad(u))


def viscous_operator(u: Expr, nu: float, form: str) -> Expr:
    """L_ν u = νΔ_h u or 2ν∇·∇ˢu, including the divergence gradient for stress viscosity."""
    if form == "laplacian":
        return nu * elementwise_laplacian(u)
    if form == "stress":
        return 2 * nu * ufl.div(ufl.sym(ufl.grad(u)))
    raise ValueError(f"unknown viscous form: {form!r}")


def viscous_bilinear(u: Expr, v: Expr, nu: float, dx: Measure, form: str) -> Form:
    """Resolved viscous bilinear form (spec: Convective forms, viscous form)."""
    if form == "laplacian":
        return nu * ufl.inner(ufl.grad(u), ufl.grad(v)) * dx
    if form == "stress":
        return 2 * nu * ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v))) * dx
    raise ValueError(f"unknown viscous form: {form!r}")


def momentum_operator(
    a: Expr,
    u: Expr,
    p: Expr,
    nu: float,
    convective_form: str,
    viscous_form: str,
) -> Expr:
    """L_a(u, p) = N − L_ν u + ∇p (spec: Operators and Chapter 5 pressure convention)."""
    return (
        residual_convection(a, u, convective_form)
        - viscous_operator(u, nu, viscous_form)
        + ufl.grad(p)
    )


def momentum_residual(
    a: Expr,
    u: Expr,
    p: Expr,
    nu: float,
    dt_u: Expr,
    f: Expr,
    convective_form: str,
    viscous_form: str,
) -> Expr:
    """R_u = f − D_t u − L_a(u, p) (spec: Operators and parameters)."""
    return f - dt_u - momentum_operator(a, u, p, nu, convective_form, viscous_form)
