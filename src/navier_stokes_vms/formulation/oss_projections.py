"""τ-weighted L² projection equations of the OSS residuals (spec: OSS projections)."""

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure


def projection_equations(
    xi: Expr,
    zeta: Expr | None,
    w: Expr,
    r: Expr | None,
    residual: Expr,
    divergence: Expr,
    tau_t: Expr,
    tau_c: Expr,
    dx: Measure,
) -> Form:
    """(τ_t (ξ_h − X), w_h) + (τ_c (ζ_h − ∇·u_h), r_h), consistent mass with the common rule.

    X = R̂_u + χ_D ũ_hist/δt. The spaces of ξ_h and ζ_h carry no Dirichlet condition and no null
    space. ζ_h, r_h and the second equation are None when τ_c is off.
    """
    form = tau_t * ufl.inner(xi - residual, w) * dx
    if zeta is not None:
        form += tau_c * (zeta - divergence) * r * dx
    return form
