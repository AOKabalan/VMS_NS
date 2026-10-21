"""Elementwise momentum and adjoint operators (spec: Operators and parameters)."""

import ufl
from ufl.core.expr import Expr


def elementwise_laplacian(u: Expr) -> Expr:
    """Δ_h u inside cell integrals (spec: Operators and parameters)."""
    return ufl.div(ufl.grad(u))


def momentum_operator(a: Expr, u: Expr, p: Expr, nu: float) -> Expr:
    """L_a(u, p) = a·∇u − νΔ_h u + ∇p (spec: Operators and parameters)."""
    return ufl.dot(ufl.grad(u), a) - nu * elementwise_laplacian(u) + ufl.grad(p)


def adjoint_operator(a: Expr, v: Expr, q: Expr, nu: float) -> Expr:
    """L*_a(v, q) = −a·∇v − νΔ_h v − ∇q (spec: Operators and parameters)."""
    return -ufl.dot(ufl.grad(v), a) - nu * elementwise_laplacian(v) - ufl.grad(q)


def momentum_residual(a: Expr, u: Expr, p: Expr, nu: float, dt_u: Expr, f: Expr) -> Expr:
    """R_u = f − D_t u − L_a(u, p) (spec: Operators and parameters)."""
    return f - dt_u - momentum_operator(a, u, p, nu)
