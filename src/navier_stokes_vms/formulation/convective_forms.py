"""Galerkin convective form: skew-symmetric type 1 with the outflow term (spec: Operators)."""

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure


def convection_form(
    a: Expr, u_h: Expr, u: Expr, v: Expr, dx: Measure, outflow_ds: Measure | None
) -> Form:
    """c(a; u, v) = ½(a·∇u, v) − ½(a·∇v, u) + ½⟨(u_h·n) u, v⟩ on Γ_N."""
    form = (ufl.inner(ufl.dot(ufl.grad(u), a), v) - ufl.inner(ufl.dot(ufl.grad(v), a), u)) / 2 * dx
    if outflow_ds is None:
        return form
    n = ufl.FacetNormal(outflow_ds.ufl_domain())
    return form + ufl.dot(u_h, n) * ufl.inner(u, v) / 2 * outflow_ds
