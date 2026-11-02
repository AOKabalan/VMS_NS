"""Solve-level thermal and buoyancy identities with actual discrete histories (B19–B24)."""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import firedrake as fd
import ufl
from ufl.core.expr import Expr
from ufl.measure import Measure

from navier_stokes_vms.config import Config, ThermalOptions
from navier_stokes_vms.formulation.convective_forms import convection_form
from navier_stokes_vms.formulation.galerkin import SolveLevel, time_derivative
from navier_stokes_vms.formulation.stabilization_parameters import transport_parameter
from navier_stokes_vms.formulation.thermal import (
    scalar_convection,
    scalar_galerkin_residual,
    scalar_vms_terms,
)
from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, Dirichlet
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.problem.scalar_problem import transported_scalars
from navier_stokes_vms.solver.reusable_assembly import reusable_assembler, scalar_functional

CONTENT_TERMS = (
    "content",
    "total_content",
    "subscale_moment",
    "rate",
    "source",
    "flux",
    "splitting_defect",
    "continuity_work",
    "iterate_work",
    "residual",
    "normalized_residual",
)
VARIANCE_TERMS = (
    "energy",
    "rate",
    "diffusion",
    "subscale",
    "diffusion_cross",
    "convection",
    "projection",
    "source",
    "flux",
    "reaction_work",
    "residual",
    "normalized_residual",
)
ENERGY_TERMS = (
    "kinetic",
    "potential",
    "kinetic_rate",
    "potential_rate",
    "viscous",
    "subscale",
    "viscous_cross",
    "pressure_subscale",
    "momentum_convection",
    "projection",
    "external_work",
    "thermal_source",
    "potential_diffusion",
    "potential_flux",
    "buoyancy_work",
    "potential_advection",
    "exchange_linear",
    "exchange_nonlinear",
    "exchange_convection",
    "exchange_reference",
    "exchange_iterate",
    "exchange_total",
    "residual",
    "normalized_residual",
)
THERMAL_BUDGETS = (
    "thermal_content_budget",
    "thermal_variance_budget",
    "boussinesq_energy_budget",
    "thermal_boundary_flux",
    "thermal_subscale_l2",
    "cavity_nusselt",
)
NORMALIZATION_FLOOR = 1e-12


@dataclass(frozen=True)
class ThermalBudgetFields:
    """Matched solve-level fields, parameters and history expressions used by the solver."""

    config: Config
    problem: Problem
    velocity: fd.Function
    advector: Expr
    velocity_subscale: Expr | None
    velocity_history: Expr | None
    scalar_fields: Mapping[str, fd.Function]
    scalar_subscales: Mapping[str, Expr]
    scalar_histories: Mapping[str, Expr]
    scalar_projections: Mapping[str, Expr]
    scalar_fluxes: Mapping[str, Mapping[str, Expr]]
    momentum_level: SolveLevel
    scalar_levels: Mapping[str, SolveLevel]
    xi: Expr | None
    zeta: Expr | None
    tau_u: Expr | None
    tau_c: Expr | None
    dx: Measure
    ds: Measure


def thermal_budget_columns(name: str, problem: Problem) -> tuple[str, ...]:
    """Components of the requested group, naming every transported scalar explicitly."""
    if name == "boussinesq_energy_budget":
        return ENERGY_TERMS
    if name == "cavity_nusselt":
        return ("volume", "hot_gradient", "cold_gradient", "hot_reaction", "cold_reaction")
    scalars = transported_scalars(problem)
    if name == "thermal_subscale_l2":
        return tuple(s.name for s in scalars)
    if name == "thermal_boundary_flux":
        return tuple(
            f"{s.name}_{boundary}_{kind}"
            for s in scalars
            for boundary in problem.boundary_tags
            for kind in ("gradient", "reaction")
        )
    terms = CONTENT_TERMS if name == "thermal_content_budget" else VARIANCE_TERMS
    return tuple(f"{s.name}_{term}" for s in scalars for term in terms)


def _boundary_masks(
    field: fd.Function, data: Mapping[str, object], tags: Mapping[str, int]
) -> dict[str, fd.Function]:
    space = field.function_space().collapse()
    masks = {}
    for name, bc in data.items():
        if isinstance(bc, Dirichlet):
            mask = fd.Function(space)
            fd.DirichletBC(space, 1.0, tags[name]).apply(mask)
            masks[name] = mask
    if masks:
        total = fd.Function(space).assign(sum(masks.values()))
        for mask in masks.values():
            mask.interpolate(ufl.conditional(total > 0, mask / total, 0.0))
    return masks


class ThermalBudgets:
    """Capture tested thermal identities before endpoint extrapolation or history shifting."""

    def __init__(self, names: Sequence[str], fields: ThermalBudgetFields) -> None:
        self._fields = fields
        names = tuple(
            dict.fromkeys(
                (*names, *(("thermal_boundary_flux",) if "cavity_nusselt" in names else ()))
            )
        )
        self._names = names
        self._forms = {}
        self._reactions = {}
        self._values = {}
        config, problem, dx, ds = fields.config, fields.problem, fields.dx, fields.ds
        thermal, stabilization = problem.thermal, config.stabilization
        vms = stabilization.method == "vms"
        A = int(vms and stabilization.subscale_space == "asgs")
        D = int(vms and stabilization.subscale_dynamics == "dynamic")
        N = int(vms and stabilization.scale_splitting == "nonlinear")
        form = (config.thermal or ThermalOptions()).convection
        lam = {"nonconservative": 1.0, "skew_symmetric": 0.5, "conservative": 0.0}[form]
        u, a = fields.velocity, fields.advector
        w = fields.velocity_subscale if vms else ufl.zero(u.ufl_shape)
        difference = a - u - N * w
        d_u = time_derivative(
            u, fields.momentum_level.history, fields.momentum_level.alpha0, fields.momentum_level.dt
        )
        d_w = (
            time_derivative(
                w, fields.velocity_history, fields.momentum_level.alpha0, fields.momentum_level.dt
            )
            if D
            else ufl.zero(u.ufl_shape)
        )
        g = ufl.as_vector(thermal.gravity)
        phi = ufl.dot(g, fd.SpatialCoordinate(problem.mesh))
        self._nusselt_scale = {}
        if "cavity_nusselt" in names:
            reference = problem.reference_values or {}
            if (
                not {"length_scale", "temperature_difference"} <= reference.keys()
                or not {"x_min", "x_max"} <= problem.boundary_tags.keys()
            ):
                raise ValueError(
                    "cavity_nusselt needs cavity length/temperature scales and x wall tags"
                )
            if any(
                not math.isfinite(reference[key]) or reference[key] <= 0
                for key in ("length_scale", "temperature_difference")
            ) or any(
                not isinstance(thermal.boundary_conditions[wall], Dirichlet)
                for wall in ("x_min", "x_max")
            ):
                raise ValueError(
                    "cavity_nusselt needs positive scales and prescribed x-wall temperatures"
                )
            factor = reference["length_scale"] / (
                thermal.diffusivity * reference["temperature_difference"]
            )
            for boundary in ("x_min", "x_max"):
                area = fd.assemble(fd.Constant(1.0) * ds(problem.boundary_tags[boundary]))
                if area <= 0:
                    raise ValueError("cavity_nusselt requires nonempty x walls")
                self._nusselt_scale[boundary] = factor / area
            volume = fd.assemble(fd.Constant(1.0) * dx)
            t = fields.scalar_fields["temperature"]
            self._add(
                "cavity_nusselt",
                {"volume": factor / volume * (u[0] * t - thermal.diffusivity * t.dx(0)) * dx},
            )
        if "boussinesq_energy_budget" in names:
            if any(isinstance(bc, Dirichlet) for bc in thermal.boundary_conditions.values()):
                raise ValueError(
                    "boussinesq_energy_budget needs an admissible potential test; "
                    "thermal Dirichlet walls are excluded"
                )
            space = fields.scalar_fields["temperature"].function_space()
            represented = fd.Function(space).interpolate(phi)
            error = fd.assemble((represented - phi) ** 2 * dx)
            if error > 1e-22 * max(1.0, scalar_functional(phi**2 * dx)()):
                raise ValueError(
                    "potential test g·x must be single-valued; "
                    "gravity-direction periodicity is inadmissible"
                )
        normal = fd.FacetNormal(problem.mesh)
        self._normal_flow = scalar_functional(ufl.dot(u, normal) ** 2 * ds)
        self._wall_speed = scalar_functional(ufl.inner(u, u) * ds)
        for scalar in transported_scalars(problem):
            name = scalar.name
            t, level = fields.scalar_fields[name], fields.scalar_levels[name]
            r = fields.scalar_subscales.get(name, ufl.as_ufl(0.0))
            d_t = time_derivative(t, level.history, level.alpha0, level.dt)
            d_r = (
                time_derivative(r, fields.scalar_histories[name], level.alpha0, level.dt)
                if D
                else ufl.as_ufl(0.0)
            )
            kappa = scalar.diffusivity
            masks = _boundary_masks(t, scalar.boundary_conditions, problem.boundary_tags)
            if "thermal_content_budget" in names and masks:
                raise ValueError(
                    "thermal_content_budget needs an admissible constant test; "
                    "scalar Dirichlet walls are excluded"
                )
            s = fd.TestFunction(t.function_space().collapse())
            residual = scalar_galerkin_residual(t, s, a, kappa, level, dx, form)
            if vms:
                residual += scalar_vms_terms(r, s, a, kappa, dx, d_r if A and D else None)
            flux_one, flux_t, flux_phi = 0 * dx, 0 * dx, 0 * dx
            for boundary, bc in scalar.boundary_conditions.items():
                measure = ds(problem.boundary_tags[boundary])
                if isinstance(bc, DiffusiveFlux) and bc.value is not None:
                    # The stepper supplies the same endpoint weights through the level data.
                    flux = fields.scalar_fluxes[name][boundary]
                    residual += flux * s * measure
                    flux_one += flux * measure
                    flux_t += flux * t * measure
                    flux_phi += flux * phi * measure
                if "thermal_boundary_flux" in names:
                    self._add(
                        "thermal_boundary_flux",
                        {
                            f"{name}_{boundary}_gradient": -kappa
                            * ufl.dot(ufl.grad(t), normal)
                            * measure,
                        },
                    )
            if masks and {"thermal_boundary_flux", "thermal_variance_budget"}.intersection(names):
                self._reactions[name] = (reusable_assembler(residual), masks, t)
            if "thermal_subscale_l2" in names:
                self._add("thermal_subscale_l2", {name: r * r * dx})
            if "thermal_content_budget" in names:
                j = ufl.dot(w, ufl.grad(t))
                self._add(
                    "thermal_content_budget",
                    {
                        f"{name}_{key}": value
                        for key, value in {
                            "content": t * dx,
                            "total_content": (t + A * D * r) * dx,
                            "subscale_moment": r * dx,
                            "rate": (d_t + A * D * d_r) * dx,
                            "source": level.forcing * dx,
                            "flux": flux_one,
                            "splitting_defect": lam * (1 - N) * j * dx,
                            "continuity_work": lam * (t * ufl.div(u) - j) * dx,
                            "iterate_work": -lam * ufl.dot(difference, ufl.grad(t)) * dx,
                        }.items()
                    },
                )
            if "thermal_variance_budget" in names:
                eta = ufl.as_ufl(0.0)
                if vms and not A:
                    eta = fields.scalar_projections[name]
                tau = (
                    transport_parameter(
                        a,
                        kappa,
                        config.discretization.velocity_degree,
                        stabilization.c1,
                        stabilization.c2,
                    )
                    if vms
                    else 1.0
                )
                self._add(
                    "thermal_variance_budget",
                    {
                        f"{name}_{key}": value
                        for key, value in {
                            "energy": (t * t + D * r * r + 2 * A * D * t * r) / 2 * dx,
                            "rate": (d_t * t + A * D * d_r * t + D * d_r * r + A * d_t * r) * dx,
                            "diffusion": kappa * ufl.inner(ufl.grad(t), ufl.grad(t)) * dx,
                            "subscale": r * r / tau * dx,
                            "diffusion_cross": -2 * kappa * r * ufl.div(ufl.grad(t)) * dx,
                            "convection": scalar_convection(a, t, t, dx, form),
                            "projection": eta * r * dx,
                            "source": level.forcing * (t + r) * dx,
                            "flux": flux_t,
                        }.items()
                    },
                )
            if name == "temperature" and "boussinesq_energy_budget" in names:
                beta, nu = thermal.expansion_coefficient, problem.viscosity
                t_ref = thermal.reference_temperature
                linear = beta * (N - 1) * t * ufl.dot(g, w)
                nonlinear = beta * N * r * ufl.dot(g, w)
                defect = (
                    beta * lam * (ufl.div(u) * t * phi - ufl.dot(a - u, phi * ufl.grad(t) + t * g))
                )
                reference = beta * t_ref * ufl.dot(g, u + w)
                iterate = beta * ufl.dot(g, difference) * (t + r)
                buoyancy = -beta * ((t - t_ref) * ufl.dot(g, u + w) + r * ufl.dot(g, u))
                advective = beta * (t + r) * ufl.dot(g, a) + defect
                divergence = ufl.div(u) - (fields.zeta if fields.zeta is not None else 0.0)
                self._add(
                    "boussinesq_energy_budget",
                    {
                        "kinetic": (
                            ufl.inner(u, u) + D * ufl.inner(w, w) + 2 * A * D * ufl.inner(u, w)
                        )
                        / 2
                        * dx,
                        "potential": beta * (t + A * D * r - t_ref) * phi * dx,
                        "kinetic_rate": (
                            ufl.inner(d_u, u)
                            + A * D * ufl.inner(d_w, u)
                            + D * ufl.inner(d_w, w)
                            + A * ufl.inner(d_u, w)
                        )
                        * dx,
                        "potential_rate": beta * (d_t + A * D * d_r) * phi * dx,
                        "viscous": nu * ufl.inner(ufl.grad(u), ufl.grad(u)) * dx,
                        "subscale": ufl.inner(w, w) / fields.tau_u * dx if vms else 0 * dx,
                        "viscous_cross": -2 * nu * ufl.inner(w, ufl.div(ufl.grad(u))) * dx,
                        "pressure_subscale": fields.tau_c * divergence * ufl.div(u) * dx
                        if vms
                        else 0 * dx,
                        "momentum_convection": convection_form(
                            a, u, u, u, dx, None, config.convection.form
                        ),
                        "projection": ufl.inner(fields.xi, w) * dx
                        if fields.xi is not None
                        else 0 * dx,
                        "external_work": ufl.inner(fields.momentum_level.forcing, u + w) * dx,
                        "thermal_source": beta * level.forcing * phi * dx,
                        "potential_diffusion": beta * kappa * ufl.dot(g, ufl.grad(t)) * dx,
                        "potential_flux": beta * flux_phi,
                        "buoyancy_work": buoyancy * dx,
                        "potential_advection": advective * dx,
                        "exchange_linear": linear * dx,
                        "exchange_nonlinear": nonlinear * dx,
                        "exchange_convection": defect * dx,
                        "exchange_reference": reference * dx,
                        "exchange_iterate": iterate * dx,
                        "exchange_total": (linear + nonlinear + defect + reference + iterate) * dx,
                    },
                )

    def _add(self, group: str, forms: Mapping[str, ufl.Form]) -> None:
        self._forms.setdefault(group, {}).update(
            {key: scalar_functional(value) for key, value in forms.items()}
        )

    def capture(self) -> None:
        """Assemble all rates, works and reactions using unshifted solve-level histories."""
        if {"thermal_content_budget", "boussinesq_energy_budget"}.intersection(
            self._names
        ) and self._normal_flow() > 1e-20:
            raise ValueError(
                "thermal content/potential budgets require impermeable velocity boundaries"
            )
        if "boussinesq_energy_budget" in self._names and self._wall_speed() > 1e-20:
            raise ValueError(
                "boussinesq_energy_budget requires homogeneous no-slip or periodic velocity"
            )
        self._values = {
            group: {key: evaluate() for key, evaluate in forms.items()}
            for group, forms in self._forms.items()
        }
        fluxes = self._values.get("thermal_boundary_flux", {})
        variance = self._values.get("thermal_variance_budget", {})
        for name, (assembler, masks, t) in self._reactions.items():
            residual = assembler.assemble()
            reaction_work = 0.0
            for boundary, mask in masks.items():
                with residual.dat.vec_ro as vector, mask.dat.vec_ro as lift:
                    fluxes[f"{name}_{boundary}_reaction"] = -vector.dot(lift)
                lift_t = fd.Function(mask.function_space()).interpolate(mask * t)
                with residual.dat.vec_ro as vector, lift_t.dat.vec_ro as lift:
                    reaction_work += vector.dot(lift)
            if variance:
                variance[f"{name}_reaction_work"] = reaction_work
        if "thermal_boundary_flux" in self._names:
            for key in thermal_budget_columns("thermal_boundary_flux", self._fields.problem):
                fluxes.setdefault(key, None)
            self._values["thermal_boundary_flux"] = fluxes
        if "cavity_nusselt" in self._names:
            nusselt = self._values["cavity_nusselt"]
            for side, boundary, sign in (("hot", "x_min", -1), ("cold", "x_max", 1)):
                for kind in ("gradient", "reaction"):
                    nusselt[f"{side}_{kind}"] = (
                        sign
                        * self._nusselt_scale[boundary]
                        * fluxes[f"temperature_{boundary}_{kind}"]
                    )
        for scalar in transported_scalars(self._fields.problem):
            name = scalar.name
            if variance:
                variance.setdefault(f"{name}_reaction_work", 0.0)
                signed = {
                    term: variance[f"{name}_{term}"]
                    * (-1 if term in ("source", "reaction_work") else 1)
                    for term in VARIANCE_TERMS
                    if term not in ("energy", "residual", "normalized_residual")
                }
                self._balance(variance, signed, f"{name}_")
            content = self._values.get("thermal_content_budget", {})
            if content:
                signed = {
                    term: content[f"{name}_{term}"]
                    * (
                        -1
                        if term in ("source", "splitting_defect", "continuity_work", "iterate_work")
                        else 1
                    )
                    for term in (
                        "rate",
                        "source",
                        "flux",
                        "splitting_defect",
                        "continuity_work",
                        "iterate_work",
                    )
                }
                self._balance(content, signed, f"{name}_")
        energy = self._values.get("boussinesq_energy_budget", {})
        if energy:
            signed = {
                term: energy[term]
                * (-1 if term in ("external_work", "thermal_source", "exchange_total") else 1)
                for term in (
                    "kinetic_rate",
                    "potential_rate",
                    "viscous",
                    "subscale",
                    "viscous_cross",
                    "pressure_subscale",
                    "momentum_convection",
                    "projection",
                    "external_work",
                    "thermal_source",
                    "potential_diffusion",
                    "potential_flux",
                    "exchange_total",
                )
            }
            self._balance(energy, signed)
        if "thermal_subscale_l2" in self._values:
            self._values["thermal_subscale_l2"] = {
                name: value**0.5 for name, value in self._values["thermal_subscale_l2"].items()
            }

    @staticmethod
    def _balance(values: dict[str, float], signed: Mapping[str, float], prefix: str = "") -> None:
        residual = sum(signed.values())
        values[f"{prefix}residual"] = residual
        values[f"{prefix}normalized_residual"] = residual / max(
            NORMALIZATION_FLOOR, sum(abs(v) for v in signed.values())
        )

    def values(self, name: str) -> dict[str, float | None]:
        """Previously captured values, safe to read after the histories have shifted."""
        return self._values[name]
