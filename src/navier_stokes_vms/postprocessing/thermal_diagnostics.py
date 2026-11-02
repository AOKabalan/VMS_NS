"""Thermal Galerkin diagnostics with the shared discrete integration measure."""

from __future__ import annotations

from typing import TYPE_CHECKING

import firedrake as fd
import ufl

from navier_stokes_vms.formulation.boussinesq import buoyancy

if TYPE_CHECKING:
    from navier_stokes_vms.postprocessing.diagnostics import DiagnosticState


def temperature_l2(state: DiagnosticState) -> float:
    """L² norm of endpoint temperature."""
    return fd.assemble(state.temperature**2 * state.dx) ** 0.5


def thermal_content(state: DiagnosticState) -> float:
    """Integral of the resolved temperature."""
    return fd.assemble(state.temperature * state.dx)


def temperature_error_l2(state: DiagnosticState) -> float:
    """Endpoint temperature error against the declared exact solution."""
    x = fd.SpatialCoordinate(state.problem.mesh)
    exact = state.problem.thermal.exact_temperature(x, state.velocity_time)
    return fd.assemble((state.temperature - exact) ** 2 * state.error_dx) ** 0.5


def buoyancy_power(state: DiagnosticState) -> float:
    """Resolved BDF1 buoyancy work with the active thermal data."""
    thermal = state.problem.thermal
    force = buoyancy(
        state.temperature,
        thermal.expansion_coefficient,
        thermal.gravity,
        thermal.reference_temperature,
    )
    return fd.assemble(ufl.inner(force, state.solve_velocity) * state.dx)
