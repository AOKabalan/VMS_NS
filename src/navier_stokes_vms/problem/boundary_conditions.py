"""Boundary condition descriptions and their conversion to Firedrake by name (plan §8)."""

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from ufl.core.expr import Expr
from ufl.measure import Measure


@dataclass(frozen=True)
class Dirichlet:
    """Strong velocity data u_h = g(x, t) (spec: Boundary conditions and quadrature)."""

    value: Callable[[Expr, Expr], Expr]


@dataclass(frozen=True)
class Outflow:
    """Do-nothing outflow ν ∂n u − p n = 0 (spec: Boundary conditions and quadrature)."""


def dirichlet_bc(
    space: WithGeometry, name: str, value: Expr, boundary_tags: Mapping[str, int]
) -> fd.DirichletBC:
    """DirichletBC imposing value on the named boundary."""
    return fd.DirichletBC(space, value, boundary_tag(name, boundary_tags))


def named_boundary_measure(
    measure: Measure, names: Iterable[str], boundary_tags: Mapping[str, int]
) -> Measure:
    """The given ds restricted to the named boundaries, keeping its quadrature degree."""
    return measure(tuple(boundary_tag(name, boundary_tags) for name in names))


def boundary_tag(name: str, boundary_tags: Mapping[str, int]) -> int:
    """Mesh tag of a named boundary; an unknown name lists the available ones."""
    if name not in boundary_tags:
        available = ", ".join(sorted(boundary_tags)) or "none"
        raise ValueError(f"boundary {name!r} does not exist; available: {available}")
    return boundary_tags[name]


def check_dirichlet_flux(velocity: fd.Function, measure: Measure) -> None:
    """All-Dirichlet nodal data must satisfy the constant-test continuity equation."""
    normal = fd.FacetNormal(velocity.function_space().mesh())
    normal_velocity = fd.dot(velocity, normal)
    flux = fd.assemble(normal_velocity * measure)
    absolute_flux = fd.assemble(abs(normal_velocity) * measure)
    if abs(flux) > 1e-10 * absolute_flux:
        raise ValueError(
            "all-Dirichlet data must be discretely flux-free: interpolated boundary data "
            f"have net flux {flux:.16e}, exceeding 1e-10 times the absolute normal flux "
            f"{absolute_flux:.16e}; check the boundary values and their FE interpolation"
        )
