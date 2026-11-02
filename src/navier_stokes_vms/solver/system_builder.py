"""Composition of monolithic Galerkin and static VMS thermal systems from stateless kernels."""

from collections.abc import Mapping
from dataclasses import dataclass, replace

import firedrake as fd
import ufl
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.config import Config, ThermalOptions
from navier_stokes_vms.discretization.function_spaces import FieldLayout
from navier_stokes_vms.formulation.boussinesq import buoyancy
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.formulation.oss_projections import (
    projection_equations,
    scalar_projection_equation,
)
from navier_stokes_vms.formulation.stabilization_parameters import (
    stabilization_parameters,
    transport_parameter,
)
from navier_stokes_vms.formulation.thermal import (
    scalar_galerkin_residual,
    scalar_residual,
    scalar_subscale,
    scalar_vms_terms,
)
from navier_stokes_vms.formulation.vms import (
    subscale_residual,
    switches,
    velocity_subscale,
    vms_terms,
)
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.problem.scalar_problem import transported_scalars


@dataclass(frozen=True)
class ThermalSystem:
    """Coupled residual and eliminated subscales at a matched solve-level state."""

    residual: Form
    velocity_subscale: Expr | None
    scalar_subscales: Mapping[str, Expr]


def thermal_system(
    solution: fd.Function,
    layout: FieldLayout,
    problem: Problem,
    momentum_level: SolveLevel,
    scalar_levels: Mapping[str, SolveLevel],
    dx: Measure,
    config: Config,
    parameter_advector: Expr | None = None,
) -> ThermalSystem:
    """Compose Galerkin or static linear ASGS/OSS for temperature and passive scalars."""
    fields = layout.bind(ufl.split(solution))
    tests = layout.bind(fd.TestFunctions(solution.function_space()))
    u, p, v, q = (fields["velocity"], fields["pressure"], tests["velocity"], tests["pressure"])
    thermal = problem.thermal
    level = replace(
        momentum_level,
        forcing=momentum_level.forcing
        + buoyancy(
            fields["temperature"],
            thermal.expansion_coefficient,
            thermal.gravity,
            thermal.reference_temperature,
        ),
    )
    form = galerkin_residual(
        u, p, v, q, u, u, problem.viscosity, level, dx, None, config.convection.form
    )
    stabilization = config.stabilization
    chi = switches(stabilization) if stabilization.method == "vms" else None
    a_tau = u if parameter_advector is None else parameter_advector
    scalar_form = (config.thermal or ThermalOptions()).convection
    u_sub = None
    subscales = {}
    if chi is not None:
        k = config.discretization.velocity_degree
        tau_u, tau_c = stabilization_parameters(
            a_tau, problem.viscosity, k, stabilization.c1, stabilization.c2
        )
        if not stabilization.pressure_subscale:
            tau_c = ufl.zero()
        xi, zeta = fields.get("xi"), fields.get("zeta")
        form += vms_terms(
            u,
            p,
            v,
            q,
            u,
            problem.viscosity,
            level,
            tau_u,
            tau_c,
            dx,
            chi,
            None,
            xi,
            zeta,
            config.convection.form,
        )
        u_sub = velocity_subscale(
            u, p, u, problem.viscosity, level, tau_u, chi, None, xi, config.convection.form
        )
        if not chi.is_asgs:
            residual = subscale_residual(
                u, p, u, problem.viscosity, level, chi, None, config.convection.form
            )
            form += projection_equations(
                xi, zeta, tests["xi"], tests.get("zeta"), residual, ufl.div(u), tau_u, tau_c, dx
            )
    for scalar in transported_scalars(problem):
        field, test = fields[scalar.name], tests[scalar.name]
        scalar_level = scalar_levels[scalar.name]
        form += scalar_galerkin_residual(
            field, test, u, scalar.diffusivity, scalar_level, dx, scalar_form
        )
        if chi is not None:
            tau = transport_parameter(
                a_tau, scalar.diffusivity, k, stabilization.c1, stabilization.c2
            )
            residual = scalar_residual(field, u, scalar.diffusivity, scalar_level, chi.is_asgs)
            projection = fields.get(f"eta_{scalar.name}")
            subscale = scalar_subscale(residual, tau, projection)
            subscales[scalar.name] = subscale
            form += scalar_vms_terms(subscale, test, u, scalar.diffusivity, dx)
            if not chi.is_asgs:
                form += scalar_projection_equation(
                    projection, tests[f"eta_{scalar.name}"], residual, tau, dx
                )
    if chi is not None:
        form += (
            thermal.expansion_coefficient
            * subscales["temperature"]
            * ufl.dot(ufl.as_vector(thermal.gravity), v)
            * dx
        )
    return ThermalSystem(form, u_sub, subscales)
