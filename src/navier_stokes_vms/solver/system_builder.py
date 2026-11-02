"""Composition of the monolithic Galerkin thermal system from stateless kernels."""

from dataclasses import replace

import firedrake as fd
import ufl
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.discretization.function_spaces import FieldLayout
from navier_stokes_vms.formulation.boussinesq import buoyancy
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.formulation.thermal import scalar_galerkin_residual
from navier_stokes_vms.problem.thermal_problem import ThermalProblem


def galerkin_thermal_system(
    solution: fd.Function,
    layout: FieldLayout,
    viscosity: float,
    momentum_level: SolveLevel,
    scalar_level: SolveLevel,
    thermal: ThermalProblem,
    dx: Measure,
    momentum_convection: str,
    scalar_convection: str,
) -> Form:
    """Coupled momentum, continuity and temperature with buoyancy differentiated by Newton."""
    fields = layout.bind(ufl.split(solution))
    tests = layout.bind(fd.TestFunctions(solution.function_space()))
    u, p, temperature = (fields[name] for name in ("velocity", "pressure", "temperature"))
    v, q, test = (tests[name] for name in ("velocity", "pressure", "temperature"))
    level = replace(
        momentum_level,
        forcing=momentum_level.forcing
        + buoyancy(
            temperature,
            thermal.expansion_coefficient,
            thermal.gravity,
            thermal.reference_temperature,
        ),
    )
    return galerkin_residual(
        u, p, v, q, u, u, viscosity, level, dx, None, momentum_convection
    ) + scalar_galerkin_residual(
        temperature, test, u, thermal.diffusivity, scalar_level, dx, scalar_convection
    )
