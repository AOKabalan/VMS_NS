"""Constant-property Boussinesq buoyancy (thermal formulation B1)."""

import ufl
from ufl.core.expr import Expr


def buoyancy(temperature: Expr, beta: float, gravity: tuple[float, ...], reference: float) -> Expr:
    """Acceleration −β(T − T_ref)g with gravity pointing downwards (B1)."""
    return -beta * (temperature - reference) * ufl.as_vector(gravity)
