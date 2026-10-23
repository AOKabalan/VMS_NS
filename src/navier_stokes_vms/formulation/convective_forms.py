"""Chapter 5 convection, subscale coupling, residual convection and outflow terms."""

import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

PRESSURE_UNKNOWNS = {
    "nonconservative": "p",
    "conservative": "p",
    "skew_symmetric_1": "p",
    "skew_symmetric_2": "p",
    "emac": "P",
    "emac_vms": "P",
}


def divergence_pairing(u_h: Expr, u_sub: Expr, w: Expr, dx: Measure) -> Form:
    """(∇·(u_h + ũ), w) := (∇·u_h, w) − (ũ, ∇w) (spec: Convective forms)."""
    return (ufl.div(u_h) * w - ufl.inner(u_sub, ufl.grad(w))) * dx


def minimal_stress(u: Expr, u_sub: Expr, nonlinear: bool = True) -> Expr:
    """T̃ = ũ⊗u + u⊗ũ + χ_N ũ⊗ũ, with the same expression in every slot (spec: EMAC–VMS)."""
    stress = ufl.outer(u_sub, u) + ufl.outer(u, u_sub)
    return stress + ufl.outer(u_sub, u_sub) if nonlinear else stress


def resolved_convection(
    a: Expr, u_h: Expr, u: Expr, v: Expr, dx: Measure, form: str = "skew_symmetric_1"
) -> Form:
    """The selected resolved volume form; u_h is the resolved advector (spec: Chapter 5)."""
    if form == "skew_symmetric_1":
        return (
            (ufl.inner(ufl.dot(ufl.grad(u), a), v) - ufl.inner(ufl.dot(ufl.grad(v), a), u)) / 2 * dx
        )
    if form == "nonconservative":
        return ufl.inner(ufl.dot(ufl.grad(u), a), v) * dx
    if form == "conservative":
        return -ufl.inner(u, ufl.dot(ufl.grad(v), a)) * dx
    if form == "skew_symmetric_2":
        return ufl.inner(ufl.dot(ufl.grad(u), a), v) * dx + 0.5 * divergence_pairing(
            u_h, a - u_h, ufl.inner(u, v), dx
        )
    if form == "emac_vms":
        # Resolved EMAC uses u alone, including with nonlinear subscales.
        return (
            ufl.inner(2 * ufl.dot(ufl.sym(ufl.grad(u)), u), v) + ufl.div(u) * ufl.inner(u, v)
        ) * dx
    if form == "emac":
        return ufl.inner(2 * ufl.dot(ufl.sym(ufl.grad(u)), a), v) * dx + divergence_pairing(
            u_h, a - u_h, ufl.inner(u, v), dx
        )
    raise ValueError(f"unknown convective form: {form!r}")


def convective_subscale_coupling(
    a: Expr,
    u: Expr,
    v: Expr,
    u_sub: Expr,
    dx: Measure,
    form: str = "skew_symmetric_1",
    nonlinear: bool = True,
) -> Form:
    """Convective subscale feedback, replacing the standard coupling for minimal EMAC–VMS."""
    if form == "emac_vms":
        return -ufl.inner(minimal_stress(u, u_sub, nonlinear), ufl.sym(ufl.grad(v))) * dx
    if form not in PRESSURE_UNKNOWNS:
        raise ValueError(f"unknown convective form: {form!r}")
    return -ufl.inner(u_sub, ufl.dot(ufl.grad(v), a)) * dx


def residual_convection(a: Expr, u: Expr, form: str = "skew_symmetric_1") -> Expr:
    """N = (a·∇)u, plus the kinetic-pressure gradient for a P unknown (spec: Chapter 5)."""
    convection = ufl.dot(ufl.grad(u), a)
    if PRESSURE_UNKNOWNS[form] == "P":
        return convection + ufl.dot(u, ufl.grad(u))
    return convection


def outflow_term(
    u_h: Expr, u: Expr, v: Expr, outflow_ds: Measure | None, form: str = "skew_symmetric_1"
) -> Form:
    """Corrective facet term; EMAC retains its modified-pressure traction (spec: outflow)."""
    if form not in PRESSURE_UNKNOWNS:
        raise ValueError(f"unknown convective form: {form!r}")
    weight = {"conservative": 1, "skew_symmetric_1": 0.5}.get(form, 0)
    if outflow_ds is None or weight == 0:
        return Form([])
    n = ufl.FacetNormal(outflow_ds.ufl_domain())
    return weight * ufl.dot(u_h, n) * ufl.inner(u, v) * outflow_ds


def convection_form(
    a: Expr,
    u_h: Expr,
    u: Expr,
    v: Expr,
    dx: Measure,
    outflow_ds: Measure | None,
    form: str = "skew_symmetric_1",
) -> Form:
    """Resolved volume convection plus its outflow term (spec: Convective forms)."""
    return resolved_convection(a, u_h, u, v, dx, form) + outflow_term(u_h, u, v, outflow_ds, form)
