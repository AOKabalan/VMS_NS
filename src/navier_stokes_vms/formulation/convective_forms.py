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
    "emac_vms_full": "P",
}
EMAC_VMS_FORMS = ("emac_vms", "emac_vms_full")


def physical_pressure(pressure: Expr, u_h: Expr, form: str) -> Expr:
    """p_h = P_h + ½|u_h|² for a P unknown, else the unknown itself (spec: EMAC pressure)."""
    if PRESSURE_UNKNOWNS[form] == "P":
        return pressure + ufl.inner(u_h, u_h) / 2
    return pressure


def divergence_pairing(u_h: Expr, u_sub: Expr, w: Expr, dx: Measure) -> Form:
    """(∇·(u_h + ũ), w) := (∇·u_h, w) − (ũ, ∇w) (spec: Convective forms)."""
    return (ufl.div(u_h) * w - ufl.inner(u_sub, ufl.grad(w))) * dx


def minimal_stress(u: Expr, u_sub: Expr, nonlinear: bool) -> Expr:
    """T̃ = ũ⊗u + u⊗ũ + χ_N ũ⊗ũ, with the same expression in every slot (spec: EMAC–VMS)."""
    stress = ufl.outer(u_sub, u) + ufl.outer(u, u_sub)
    return stress + ufl.outer(u_sub, u_sub) if nonlinear else stress


def full_stress(u: Expr, u_sub: Expr, nonlinear: bool) -> Expr:
    """T̃_full = S(U) − S(u) = T̃ + (u·ũ + χ_N ½|ũ|²) I (chapter §5.8.5, Eq. 5.37)."""
    kinetic = ufl.inner(u, u_sub) + (ufl.inner(u_sub, u_sub) / 2 if nonlinear else 0)
    return minimal_stress(u, u_sub, nonlinear) + kinetic * ufl.Identity(u.ufl_shape[0])


def full_stress_matrix(u: Expr) -> Expr:
    """M_h = ∇ˢu + ½(∇·u) I, the subscale-dependent part of N_E (chapter Eq. 5.38)."""
    return ufl.sym(ufl.grad(u)) + ufl.div(u) * ufl.Identity(u.ufl_shape[0]) / 2


def resolved_convection(a: Expr, u_h: Expr, u: Expr, v: Expr, dx: Measure, form: str) -> Form:
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
    if form in EMAC_VMS_FORMS:
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
    form: str,
    nonlinear: bool,
) -> Form:
    """Convective subscale feedback, replacing the standard coupling for minimal EMAC–VMS."""
    if form == "emac_vms":
        return -ufl.inner(minimal_stress(u, u_sub, nonlinear), ufl.sym(ufl.grad(v))) * dx
    if form == "emac_vms_full":
        return -ufl.inner(full_stress(u, u_sub, nonlinear), ufl.sym(ufl.grad(v))) * dx
    if form not in PRESSURE_UNKNOWNS:
        raise ValueError(f"unknown convective form: {form!r}")
    return -ufl.inner(u_sub, ufl.dot(ufl.grad(v), a)) * dx


def residual_convection(a: Expr, u: Expr, form: str) -> Expr:
    """N = (a·∇)u, plus the kinetic-pressure gradient for a P unknown (spec: Chapter 5).

    For emac_vms_full, N_E = 2(∇ˢu)u + (∇·u)u + M_h ũ with ũ = a − u (chapter Eq. 5.38).
    """
    if form == "emac_vms_full":
        return (
            2 * ufl.dot(ufl.sym(ufl.grad(u)), u)
            + ufl.div(u) * u
            + ufl.dot(full_stress_matrix(u), a - u)
        )
    convection = ufl.dot(ufl.grad(u), a)
    if PRESSURE_UNKNOWNS[form] == "P":
        return convection + ufl.dot(u, ufl.grad(u))
    return convection


def outflow_term(u_h: Expr, u: Expr, v: Expr, outflow_ds: Measure | None, form: str) -> Form:
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


def convective_defect(
    a: Expr,
    u_h: Expr,
    w: Expr,
    u_sub: Expr | None,
    dx: Measure,
    outflow_ds: Measure | None,
    form: str,
    nonlinear: bool,
) -> Form:
    """D(w) = c(a; u_h, w) − (ũ, (a·∇)w), or G(w) for EMAC–VMS (spec: Conservation checks)."""
    defect = convection_form(a, u_h, u_h, w, dx, outflow_ds, form)
    if u_sub is None:
        return defect
    return defect + convective_subscale_coupling(a, u_h, w, u_sub, dx, form, nonlinear)


def energy_defect(
    a: Expr,
    u_h: Expr,
    u_sub: Expr | None,
    dx: Measure,
    outflow_ds: Measure | None,
    form: str,
    nonlinear: bool,
) -> Form:
    """𝒩 = D(u_h) + (ũ, N), the total convective energy defect (spec: Conservation checks)."""
    defect = convective_defect(a, u_h, u_h, u_sub, dx, outflow_ds, form, nonlinear)
    if u_sub is None:
        return defect
    return defect + ufl.inner(u_sub, residual_convection(a, u_h, form)) * dx


def predicted_energy_defect(
    u_h: Expr, u_sub: Expr, dx: Measure, form: str, nonlinear: bool
) -> Form:
    """Interior 𝒩 of chapter Tables 5.1–5.2, with D_χ(w) = (∇·u_h, w) − χ_N(ũ, ∇w)."""
    pairing_subscale = u_sub if nonlinear else ufl.zero(u_h.ufl_shape)
    pairing = divergence_pairing(u_h, pairing_subscale, ufl.inner(u_h, u_h), dx)
    if form == "nonconservative":
        return -0.5 * pairing
    if form == "conservative":
        return 0.5 * pairing
    if form == "emac":
        transfer = ufl.dot(ufl.grad(u_h), u_h) if nonlinear else ufl.dot(u_h, ufl.grad(u_h))
        return ufl.inner(u_sub, transfer) * dx
    if form not in PRESSURE_UNKNOWNS:
        raise ValueError(f"unknown convective form: {form!r}")
    return Form([])


def predicted_momentum_defect(
    u_h: Expr, u_sub: Expr | None, rigid: Expr, dx: Measure, form: str, nonlinear: bool
) -> Form:
    """Interior D_k or D_φ of chapter Tables 5.2–5.3 for the rigid test e_k or φ = e_k × x.

    With D_χ(w) = (∇·u_h, w) − χ_N(ũ, ∇w) and J = (ũ, (u_h·∇)r), which vanishes for e_k
    (spec: Conservation checks, predicted momentum defects).
    """
    if form not in PRESSURE_UNKNOWNS:
        raise ValueError(f"unknown convective form: {form!r}")
    zero = ufl.zero(u_h.ufl_shape)
    pairing_subscale = u_sub if nonlinear and u_sub is not None else zero
    pairing = divergence_pairing(u_h, pairing_subscale, ufl.inner(u_h, rigid), dx)
    transfer = Form([])
    if u_sub is not None and not nonlinear:
        transfer = -ufl.inner(u_sub, ufl.dot(ufl.grad(rigid), u_h)) * dx
    if form == "nonconservative":
        return -pairing + transfer
    if form in ("skew_symmetric_1", "skew_symmetric_2"):
        return -0.5 * pairing + transfer
    if form == "emac" and nonlinear and u_sub is not None:
        return ufl.inner(u_sub, ufl.dot(ufl.grad(u_h), rigid)) * dx
    if form in EMAC_VMS_FORMS:
        return Form([])
    return transfer
