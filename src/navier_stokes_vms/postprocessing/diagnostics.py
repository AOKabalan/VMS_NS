"""Registry of named per-step diagnostics; error norms use a rule of twice the common degree."""

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import firedrake as fd
import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.formulation.convective_forms import (
    convective_defect,
    energy_defect,
    predicted_energy_defect,
)
from navier_stokes_vms.formulation.operators import viscous_bilinear, viscous_operator
from navier_stokes_vms.problem.boundary_conditions import named_boundary_measure
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.solver.reusable_assembly import reusable_assembler

ENERGY_BUDGET_TERMS = (
    "rate",
    "viscous",
    "subscale",
    "viscous_cross",
    "pressure_subscale",
    "dissipation",
    "power",
    "defect_predicted",
    "defect_measured",
    "residual",
)


@dataclass(frozen=True)
class ConservationFields:
    """Solve-level and endpoint fields of the forms actually used (spec: Conservation checks).

    The solve level is t^(n+½) for Crank–Nicolson; a* = ū + χ_N ũ* and τ are those of the
    residual; the subscale is the substituted expression, None for Galerkin.
    """

    velocity: Expr
    advector: Expr
    subscale: Expr | None
    endpoint_velocity: fd.Function
    endpoint_subscale: fd.Function | None
    zeta: Expr | None
    tau_m: Expr | None
    tau_c: Expr | None
    forcing: Expr
    viscosity: float
    dt: float
    dx: Measure
    outflow_ds: Measure | None
    convective_form: str
    viscous_form: str
    is_nonlinear: bool
    is_asgs: bool


def _functional(form: Form) -> Callable[[], float]:
    """A scalar functional preprocessed once; one that vanishes identically evaluates to 0."""
    expanded = ufl.algorithms.expand_derivatives(form)
    if not any(integral.integrand() != 0 for integral in expanded.integrals()):
        return lambda: 0.0
    return reusable_assembler(form).assemble


class Conservation:
    """Convective defects and the Crank–Nicolson energy closure (spec: Conservation checks).

    `capture` runs after each converged solve, before the histories advance, so every term uses
    the solve-level state and the old endpoint; `energy_budget` then adds only the new endpoint
    energy.
    """

    def __init__(self, names: Sequence[str], fields: ConservationFields) -> None:
        mesh = fields.endpoint_velocity.function_space().mesh()
        dimension = mesh.geometric_dimension
        axes = "xyz"[:dimension]
        u, s, dx = fields.velocity, fields.subscale, fields.dx

        def defect(w: Expr) -> Callable[[], float]:
            return _functional(
                convective_defect(
                    fields.advector,
                    u,
                    w,
                    s,
                    dx,
                    fields.outflow_ds,
                    fields.convective_form,
                    fields.is_nonlinear,
                )
            )

        x = fd.SpatialCoordinate(mesh)
        # Rigid fields carry the mesh, so the coupling's ∇e_k is defined (and zero).
        rigid = [
            fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1)).assign(
                fd.Constant(tuple(float(i == k) for i in range(dimension)))
            )
            for k in range(dimension)
        ]
        rotations = (
            {"": ufl.as_vector((-x[1], x[0]))}
            if dimension == 2
            else {axis: ufl.cross(e, x) for axis, e in zip(axes, rigid, strict=True)}
        )
        self._defects: dict[str, dict[str, Callable[[], float]]] = {}
        if "linear_momentum_defect" in names:
            self._defects["linear_momentum_defect"] = {
                axis: defect(e) for axis, e in zip(axes, rigid, strict=True)
            }
        if "angular_momentum_defect" in names:
            self._defects["angular_momentum_defect"] = {
                axis: defect(phi) for axis, phi in rotations.items()
            }
        if "convective_energy_defect" in names:
            self._defects["convective_energy_defect"] = {
                "": _functional(
                    energy_defect(
                        fields.advector,
                        u,
                        s,
                        dx,
                        fields.outflow_ds,
                        fields.convective_form,
                        fields.is_nonlinear,
                    )
                )
            }
        self._budget: dict[str, Callable[[], float]] = {}
        self._energy: Callable[[], float] | None = None
        if "energy_budget" in names:
            nu, viscous = fields.viscosity, fields.viscous_form
            divergence = ufl.div(u) if fields.zeta is None else ufl.div(u) - fields.zeta
            self._budget = {
                name: _functional(form)
                for name, form in {
                    "viscous": viscous_bilinear(u, u, nu, dx, viscous),
                    # τ_m, not τ_t: the time part of τ_t belongs to the subscale rate.
                    "subscale": ufl.inner(s, s) / fields.tau_m * dx,
                    "viscous_cross": -2 * ufl.inner(s, viscous_operator(u, nu, viscous)) * dx,
                    "pressure_subscale": fields.tau_c * divergence * ufl.div(u) * dx,
                    "power": ufl.inner(fields.forcing, u + s) * dx,
                    "defect_predicted": predicted_energy_defect(
                        u, s, dx, fields.convective_form, fields.is_nonlinear
                    ),
                }.items()
            }
            u_end, s_end = fields.endpoint_velocity, fields.endpoint_subscale
            if fields.is_asgs:
                energy = ufl.inner(u_end + s_end, u_end + s_end) / 2 * dx
            else:
                energy = (ufl.inner(u_end, u_end) + ufl.inner(s_end, s_end)) / 2 * dx
            self._energy = _functional(energy)
        self._dt = fields.dt
        self._captured: dict[str, dict[str, float]] = {}
        self._terms: dict[str, float] | None = None
        self._old_energy = 0.0

    def capture(self, midpoint: bool) -> None:
        """Evaluate the solve-level terms and the old endpoint energy before the shift.

        The energy closure is only defined at midpoint steps; the BDF1 startup is excluded.
        """
        self._captured = {
            name: {axis: value() for axis, value in components.items()}
            for name, components in self._defects.items()
        }
        self._terms = None
        if self._energy is not None and midpoint:
            self._terms = {name: value() for name, value in self._budget.items()}
            self._old_energy = self._energy()

    def defect(self, name: str) -> float | dict[str, float]:
        """A captured defect: a scalar, or its components by axis."""
        components = self._captured[name]
        return components[""] if "" in components else components

    def energy_budget(self) -> dict[str, float | None]:
        """rate + 𝒟 + 𝒩_pred − power = r with 𝒩_meas = power − rate − 𝒟, after the shift."""
        if self._terms is None or self._energy is None:
            return dict.fromkeys(ENERGY_BUDGET_TERMS)
        terms = dict(self._terms)
        terms["rate"] = (self._energy() - self._old_energy) / self._dt
        terms["dissipation"] = sum(
            terms[name] for name in ("viscous", "subscale", "viscous_cross", "pressure_subscale")
        )
        terms["defect_measured"] = terms["power"] - terms["rate"] - terms["dissipation"]
        terms["residual"] = terms["defect_predicted"] - terms["defect_measured"]
        return {name: terms[name] for name in ENERGY_BUDGET_TERMS}


@dataclass(frozen=True)
class DiagnosticState:
    """What a diagnostic reads after a step: endpoint velocity and solve-level pressure.

    The pressure is the physical p_h, P_h + ½|u_h|² for EMAC variants (spec: EMAC pressure).
    """

    problem: Problem
    velocity: fd.Function
    pressure: Expr
    velocity_time: fd.Constant
    pressure_time: fd.Constant
    dx: Measure
    error_dx: Measure
    pressure_null_space: bool
    solve_velocity: fd.Function
    ds: Measure
    subscale_velocity: Expr | None
    subscale_time: fd.Constant | None
    conservation: Conservation | None = None


@dataclass(frozen=True)
class Diagnostic:
    """A functional, possibly with named components, and its exact-solution requirement."""

    evaluate: Callable[[DiagnosticState], float | dict[str, float | None]]
    needs_exact_solution: bool


def _kinetic_energy(state: DiagnosticState) -> float:
    return fd.assemble(ufl.inner(state.velocity, state.velocity) / 2 * state.dx)


def _subscale_kinetic_energy(state: DiagnosticState) -> float:
    return fd.assemble(ufl.inner(state.subscale_velocity, state.subscale_velocity) / 2 * state.dx)


def _components(expression: Expr, measure: Measure) -> dict[str, float]:
    return {
        axis: fd.assemble(expression[i] * measure)
        for i, axis in enumerate("xyz"[: expression.ufl_shape[0]])
    }


def _linear_momentum(state: DiagnosticState) -> dict[str, float]:
    return _components(state.velocity, state.dx)


def _angular_momentum(state: DiagnosticState) -> float | dict[str, float]:
    x = fd.SpatialCoordinate(state.problem.mesh)
    if state.problem.mesh.geometric_dimension == 2:
        return fd.assemble((x[0] * state.velocity[1] - x[1] * state.velocity[0]) * state.dx)
    return _components(ufl.cross(x, state.velocity), state.dx)


def _force(state: DiagnosticState) -> dict[str, float]:
    mesh = state.problem.mesh
    gradient = ufl.grad(state.solve_velocity)
    stress = -state.pressure * ufl.Identity(mesh.geometric_dimension)
    stress += state.problem.viscosity * (gradient + ufl.transpose(gradient))
    traction = -ufl.dot(stress, fd.FacetNormal(mesh))
    values = {}
    for name in state.problem.force_boundaries:
        ds = named_boundary_measure(state.ds, (name,), state.problem.boundary_tags)
        values.update(
            {f"{name}_{axis}": value for axis, value in _components(traction, ds).items()}
        )
    return values


def _divergence_l2(state: DiagnosticState) -> float:
    return math.sqrt(fd.assemble(ufl.div(state.velocity) ** 2 * state.dx))


def _velocity_error(state: DiagnosticState) -> Expr:
    x = fd.SpatialCoordinate(state.velocity.function_space().mesh())
    return state.velocity - state.problem.exact_velocity(x, state.velocity_time)


def _velocity_error_l2(state: DiagnosticState) -> float:
    error = _velocity_error(state)
    return math.sqrt(fd.assemble(ufl.inner(error, error) * state.error_dx))


def _velocity_gradient_error_l2(state: DiagnosticState) -> float:
    error = ufl.grad(_velocity_error(state))
    return math.sqrt(fd.assemble(ufl.inner(error, error) * state.error_dx))


def _pressure_error_l2(state: DiagnosticState) -> float:
    x = fd.SpatialCoordinate(state.problem.mesh)
    error = state.pressure - state.problem.exact_pressure(x, state.pressure_time)
    if state.pressure_null_space:
        area = fd.assemble(fd.Constant(1.0) * state.error_dx)
        error = error - fd.Constant(fd.assemble(error * state.error_dx) / area)
    return math.sqrt(fd.assemble(error**2 * state.error_dx))


def _conservation_defect(name: str) -> Callable[[DiagnosticState], float | dict[str, float]]:
    return lambda state: state.conservation.defect(name)


def _energy_budget(state: DiagnosticState) -> dict[str, float | None]:
    return state.conservation.energy_budget()


DIAGNOSTICS: dict[str, Diagnostic] = {
    "kinetic_energy": Diagnostic(_kinetic_energy, needs_exact_solution=False),
    "subscale_kinetic_energy": Diagnostic(_subscale_kinetic_energy, needs_exact_solution=False),
    "linear_momentum": Diagnostic(_linear_momentum, needs_exact_solution=False),
    "angular_momentum": Diagnostic(_angular_momentum, needs_exact_solution=False),
    "force": Diagnostic(_force, needs_exact_solution=False),
    "divergence_l2": Diagnostic(_divergence_l2, needs_exact_solution=False),
    "velocity_error_l2": Diagnostic(_velocity_error_l2, needs_exact_solution=True),
    "velocity_gradient_error_l2": Diagnostic(
        _velocity_gradient_error_l2, needs_exact_solution=True
    ),
    "pressure_error_l2": Diagnostic(_pressure_error_l2, needs_exact_solution=True),
    "linear_momentum_defect": Diagnostic(
        _conservation_defect("linear_momentum_defect"), needs_exact_solution=False
    ),
    "angular_momentum_defect": Diagnostic(
        _conservation_defect("angular_momentum_defect"), needs_exact_solution=False
    ),
    "convective_energy_defect": Diagnostic(
        _conservation_defect("convective_energy_defect"), needs_exact_solution=False
    ),
    "energy_budget": Diagnostic(_energy_budget, needs_exact_solution=False),
}

CONSERVATION_DIAGNOSTICS = (
    "linear_momentum_defect",
    "angular_momentum_defect",
    "convective_energy_defect",
    "energy_budget",
)


def check_diagnostics_supported(names: Sequence[str], problem: Problem) -> None:
    """Reject diagnostics whose required exact solution or force boundaries are absent."""
    if "force" in names and not problem.force_boundaries:
        raise ValueError("output.diagnostics: force requires the case to declare force_boundaries")
    needing = [name for name in names if DIAGNOSTICS[name].needs_exact_solution]
    if needing and problem.exact_velocity is None:
        raise ValueError(
            f"output.diagnostics: {', '.join(needing)} need an exact solution, "
            "which this case does not provide"
        )


def evaluate_diagnostics(names: Sequence[str], state: DiagnosticState) -> dict[str, float]:
    """The configured diagnostics of the current state, in configuration order."""
    values = {}
    for name in names:
        value = DIAGNOSTICS[name].evaluate(state)
        if isinstance(value, dict):
            values.update({f"{name}_{component}": entry for component, entry in value.items()})
        else:
            values[name] = value
    return values


def diagnostic_columns(names: Sequence[str], problem: Problem) -> tuple[str, ...]:
    """Scalar column names after expanding vectors and the case's measured boundaries."""
    axes = "xyz"[: problem.mesh.geometric_dimension]
    columns = []
    for name in names:
        if name in ("linear_momentum", "linear_momentum_defect") or (
            name in ("angular_momentum", "angular_momentum_defect") and len(axes) == 3
        ):
            columns.extend(f"{name}_{axis}" for axis in axes)
        elif name == "energy_budget":
            columns.extend(f"{name}_{term}" for term in ENERGY_BUDGET_TERMS)
        elif name == "force":
            columns.extend(
                f"force_{boundary}_{axis}" for boundary in problem.force_boundaries for axis in axes
            )
        else:
            columns.append(name)
    return tuple(columns)
