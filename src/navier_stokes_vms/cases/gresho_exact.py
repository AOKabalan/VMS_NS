"""Exact steady Euler Gresho vortex, separate from the legacy viscous/no-slip case."""

import ufl

from navier_stokes_vms.cases.gresho import LOWER, UPPER, frame, pressure, velocity
from navier_stokes_vms.config import Config
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.boundary_conditions import FreeSlip
from navier_stokes_vms.problem.problem_definition import Problem, at_initial_time

DEFAULT_CELLS = (48, 48)
DEFAULT_VISCOSITY = 0.0


def build(config: Config) -> Problem:
    """Unforced Euler on (−0.5, 0.5)² with only the normal velocity fixed on each wall.

    Keep ν exactly zero: positive viscosity would invalidate the exact steady reference.
    The prepared dynamic OSS configurations disable the optional pressure subscale to avoid
    zero-weight ζ equations in the resting exterior; no viscosity/speed floor is introduced.
    """
    unknown = set(config.case_parameters) - {"viscosity"}
    if unknown:
        raise ValueError(f"case_parameters: {', '.join(sorted(unknown))} not used by gresho_exact")
    viscosity = config.case_parameters.get("viscosity", DEFAULT_VISCOSITY)
    if isinstance(viscosity, bool) or not isinstance(viscosity, int | float) or viscosity != 0:
        raise ValueError("gresho_exact requires case_parameters.viscosity = 0 (exact Euler)")
    cells = mesh_loading.builtin_cells(config, DEFAULT_CELLS)
    mesh, tags = mesh_loading.rectangle(cells, LOWER, UPPER)
    return Problem(
        mesh=mesh,
        boundary_conditions={name: FreeSlip(0 if name.startswith("x_") else 1) for name in tags},
        boundary_tags=tags,
        viscosity=DEFAULT_VISCOSITY,
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
        exact_pressure=pressure,
        force_boundaries=tuple(tags),
        monitor_region=frame,
    )
