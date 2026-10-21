"""Function spaces of the discretization and the constant-pressure null space."""

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from firedrake.mesh import MeshGeometry


def mixed_space(mesh: MeshGeometry, velocity_degree: int, pressure_degree: int) -> WithGeometry:
    """W = V × Q with continuous Lagrange velocity and pressure of the configured degrees."""
    velocity = fd.VectorFunctionSpace(mesh, "CG", velocity_degree)
    pressure = fd.FunctionSpace(mesh, "CG", pressure_degree)
    return velocity * pressure


def constant_pressure_null_space(space: WithGeometry) -> fd.MixedVectorSpaceBasis:
    """Null space of W spanned by the constant pressure (plan §8, §11)."""
    constant = fd.VectorSpaceBasis(constant=True, comm=space.comm)
    return fd.MixedVectorSpaceBasis(space, [space.sub(0), constant])
