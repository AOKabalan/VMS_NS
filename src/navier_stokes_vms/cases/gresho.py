"""Gresho vortex: a steady piecewise Euler solution in a walled box (spec: Gresho vortex).

Viscosity is a case parameter (default 1e-6; viscous ASGS runs use 1e-3). At either positive
viscosity the Euler fields are reference fields, not exact solutions of the viscous problem.
"""

import math

import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.boundary_conditions import Dirichlet
from navier_stokes_vms.problem.problem_definition import Problem, at_initial_time, case_parameters

DEFAULT_CELLS = (32, 32)
LOWER = (-0.5, -0.5)
UPPER = (0.5, 0.5)
# The authorized fallback for ν = 0, where τ_c = 0 at rest leaves ζ rows of the frozen-τ OSS
# Jacobian identically zero, so the first Newton update does not exist (spec: Gresho vortex).
DEFAULT_PARAMETERS = {"viscosity": 1e-6}
CORE_RADIUS = 0.2
VORTEX_RADIUS = 0.4
FRAME = 0.45
RING_CONSTANT = -12.5 * VORTEX_RADIUS**2 + 20 * VORTEX_RADIUS - 4 * math.log(VORTEX_RADIUS)
CORE_CONSTANT = RING_CONSTANT - 20 * CORE_RADIUS + 4 * math.log(CORE_RADIUS)


def _radius(x: Expr) -> Expr:
    return ufl.sqrt(ufl.dot(x, x))


def _ring_radius(x: Expr) -> Expr:
    # Equals r wherever the ring and exterior branches apply, and keeps 1/r and ln r finite in
    # the core, where those branches are not selected.
    return ufl.max_value(_radius(x), CORE_RADIUS)


def velocity(x: Expr, _t: Expr) -> Expr:
    """u = (u_θ/r)(−y, x): u_θ/r = 5, 2/r − 5, 0; the core branch never divides by r."""
    r, ring = _radius(x), _ring_radius(x)
    angular_speed = ufl.conditional(
        r < CORE_RADIUS, 5.0, ufl.conditional(r < VORTEX_RADIUS, 2 / ring - 5, 0.0)
    )
    return angular_speed * ufl.as_vector((-x[1], x[0]))


def pressure(x: Expr, _t: Expr) -> Expr:
    """p = 12.5r² + C₁, 12.5r² − 20r + 4 ln r + C₂, 0, continuous at r = 0.2 and 0.4."""
    r, ring = _radius(x), _ring_radius(x)
    return ufl.conditional(
        r < CORE_RADIUS,
        12.5 * r**2 + CORE_CONSTANT,
        ufl.conditional(
            r < VORTEX_RADIUS, 12.5 * ring**2 - 20 * ring + 4 * ufl.ln(ring) + RING_CONSTANT, 0.0
        ),
    )


def frame(x: Expr) -> Expr:
    """Indicator of the near-wall frame max(|x|, |y|) ≥ 0.45."""
    return ufl.conditional(ufl.ge(ufl.max_value(abs(x[0]), abs(x[1])), FRAME), 1.0, 0.0)


def build(config: Config) -> Problem:
    """(−0.5, 0.5)², configurable ν, f = 0, four no-slip/force walls; Euler reference fields."""
    viscosity = case_parameters(config, DEFAULT_PARAMETERS)["viscosity"]
    cells = mesh_loading.builtin_cells(config, DEFAULT_CELLS)
    mesh, tags = mesh_loading.rectangle(cells, LOWER, UPPER)
    return Problem(
        mesh=mesh,
        boundary_conditions={name: Dirichlet(lambda _x, _t: ufl.zero(2)) for name in tags},
        boundary_tags=tags,
        viscosity=viscosity,
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
        exact_pressure=pressure,
        force_boundaries=tuple(tags),
        monitor_region=frame,
    )
