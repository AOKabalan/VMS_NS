"""Pointwise spatial stabilization parameters (spec: Operators and parameters)."""

import ufl
from ufl.core.expr import Expr
from ufl.domain import extract_unique_domain


def stabilization_parameters(a: Expr, nu: float, k: int, c1: float, c2: float) -> tuple[Expr, Expr]:
    """τ_m = (c1 k⁴ ν/h² + c2 k |a|/h)⁻¹, τ_c = h²/(c1 k⁴ τ_m)."""
    h = ufl.MinCellEdgeLength(extract_unique_domain(a))
    speed = ufl.sqrt(ufl.dot(a, a))
    tau_m = 1 / (c1 * k**4 * nu / h**2 + c2 * k * speed / h)
    tau_c = h**2 / (c1 * k**4 * tau_m)
    return tau_m, tau_c


def dynamic_parameter(tau_m: Expr, alpha0: Expr, dt: Expr) -> Expr:
    """τ_t = (α₀/δt + 1/τ_m)⁻¹ (spec: Operators and parameters)."""
    return 1 / (alpha0 / dt + 1 / tau_m)
