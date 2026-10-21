"""Boundary condition descriptions and their conversion to Firedrake by name (plan §8)."""

from collections.abc import Callable, Mapping
from dataclasses import dataclass

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from ufl.core.expr import Expr


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


def boundary_tag(name: str, boundary_tags: Mapping[str, int]) -> int:
    """Mesh tag of a named boundary; an unknown name lists the available ones."""
    if name not in boundary_tags:
        available = ", ".join(sorted(boundary_tags)) or "none"
        raise ValueError(f"boundary {name!r} does not exist; available: {available}")
    return boundary_tags[name]
