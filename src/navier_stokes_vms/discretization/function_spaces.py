"""Function spaces of the discretization and the constant-pressure null space."""

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from firedrake.mesh import MeshGeometry

from navier_stokes_vms.discretization.quadrature import quadrature_element


def mixed_space(mesh: MeshGeometry, velocity_degree: int, pressure_degree: int) -> WithGeometry:
    """W = V × Q with continuous Lagrange velocity and pressure of the configured degrees."""
    velocity = fd.VectorFunctionSpace(mesh, "CG", velocity_degree)
    pressure = fd.FunctionSpace(mesh, "CG", pressure_degree)
    return velocity * pressure


def oss_space(mesh: MeshGeometry, degree: int, pressure_projection: bool) -> WithGeometry:
    """W = V × Q × V̂ (× Q̂): u_h, p_h, ξ_h and, when τ_c is active, ζ_h (spec: OSS projections).

    V̂ and Q̂ have the elements of V and Q. They stay unconstrained because no Dirichlet condition
    or null space is ever applied to their components.
    """
    velocity = fd.VectorFunctionSpace(mesh, "CG", degree)
    pressure = fd.FunctionSpace(mesh, "CG", degree)
    spaces = [velocity, pressure, velocity] + ([pressure] if pressure_projection else [])
    return fd.MixedFunctionSpace(spaces)


def constant_pressure_null_space(space: WithGeometry) -> fd.MixedVectorSpaceBasis:
    """Null space of W spanned by the constant pressure, its component 1 (plan §8, §11)."""
    constant = fd.VectorSpaceBasis(constant=True, comm=space.comm)
    bases = [constant if index == 1 else space.sub(index) for index in range(len(space))]
    return fd.MixedVectorSpaceBasis(space, bases)


def subscale_space(mesh: MeshGeometry, degree: int) -> WithGeometry:
    """Vector space of values at the points of the dx rule of the given degree."""
    return fd.VectorFunctionSpace(mesh, quadrature_element(mesh.ufl_cell(), degree))
