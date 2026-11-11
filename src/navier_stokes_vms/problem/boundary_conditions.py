"""Boundary condition descriptions and their conversion to Firedrake by name (plan §8)."""

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from firedrake.mesh import MeshGeometry
from ufl.core.expr import Expr
from ufl.measure import Measure


@dataclass(frozen=True)
class Dirichlet:
    """Strong field data prescribed as g(x, t) (spec: Boundary conditions and quadrature)."""

    value: Callable[[Expr, Expr], Expr]


@dataclass(frozen=True)
class DiffusiveFlux:
    """Prescribed outward scalar flux −κ∇T·n; None denotes insulation (B3)."""

    value: Callable[[Expr, Expr], Expr] | None = None


@dataclass(frozen=True)
class Outflow:
    """Do-nothing outflow ν ∂n u − p n = 0 (spec: Boundary conditions and quadrature)."""


@dataclass(frozen=True)
class FreeSlip:
    """Zero normal velocity and natural tangential traction on axis-aligned walls (G4)."""

    normal_component: int


def check_free_slip_boundaries(
    mesh: MeshGeometry,
    conditions: Mapping[str, Dirichlet | FreeSlip | Outflow],
    boundary_tags: Mapping[str, int],
) -> None:
    """Reject invalid components, empty tags and walls not normal to the requested axis."""
    dimension = mesh.geometric_dimension
    for name, condition in conditions.items():
        if not isinstance(condition, FreeSlip):
            continue
        component = condition.normal_component
        if (
            isinstance(component, bool)
            or not isinstance(component, int)
            or not 0 <= component < dimension
        ):
            raise ValueError(f"free-slip normal_component must be an integer in [0, {dimension})")
        measure = fd.ds(boundary_tag(name, boundary_tags), domain=mesh, degree=2)
        area = fd.assemble(fd.Constant(1.0) * measure)
        if area <= 0:
            raise ValueError(f"free-slip boundary {name!r} has no facets")
        normal = fd.FacetNormal(mesh)
        transverse = sum(normal[axis] ** 2 for axis in range(dimension) if axis != component)
        if fd.assemble(transverse * measure) > 1e-20 * area:
            raise ValueError(f"free-slip boundary {name!r} is not normal to axis {component}")


def free_slip_bc(
    space: WithGeometry, name: str, normal_component: int, boundary_tags: Mapping[str, int]
) -> fd.DirichletBC:
    """Build one component-only wall constraint using the shared geometry validation."""
    check_free_slip_boundaries(space.mesh(), {name: FreeSlip(normal_component)}, boundary_tags)
    return dirichlet_bc(space.sub(normal_component), name, 0.0, boundary_tags)


def velocity_bcs(
    space: WithGeometry,
    conditions: Mapping[str, Dirichlet | FreeSlip | Outflow],
    values: Mapping[str, Expr],
    boundary_tags: Mapping[str, int],
) -> list[fd.DirichletBC]:
    """Strong velocity constraints, preserving unconstrained free-slip tangent components."""
    bcs = []
    for name, condition in conditions.items():
        if isinstance(condition, Dirichlet):
            bcs.append(dirichlet_bc(space, name, values[name], boundary_tags))
        elif isinstance(condition, FreeSlip):
            bcs.append(
                dirichlet_bc(space.sub(condition.normal_component), name, 0.0, boundary_tags)
            )
    return bcs


def dirichlet_bc(
    space: WithGeometry, name: str, value: Expr | float, boundary_tags: Mapping[str, int]
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
    """Prescribed normal velocities must satisfy the constant-test continuity equation."""
    normal = fd.FacetNormal(velocity.function_space().mesh())
    normal_velocity = fd.dot(velocity, normal)
    flux = fd.assemble(normal_velocity * measure)
    absolute_flux = fd.assemble(abs(normal_velocity) * measure)
    if abs(flux) > 1e-10 * absolute_flux:
        raise ValueError(
            "prescribed velocity data must be discretely flux-free: interpolated boundary data "
            f"have net flux {flux:.16e}, exceeding 1e-10 times the absolute normal flux "
            f"{absolute_flux:.16e}; check the boundary values and their FE interpolation"
        )
