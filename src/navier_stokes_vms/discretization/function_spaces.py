"""Function spaces of the discretization and the constant-pressure null space."""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TypeVar

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from firedrake.mesh import MeshGeometry

from navier_stokes_vms.config import StabilizationOptions
from navier_stokes_vms.discretization.quadrature import quadrature_element

_T = TypeVar("_T")


@dataclass(frozen=True)
class FieldLayout:
    """Ordered mixed fields with velocity first and the pressure gauge at component one."""

    names: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.names[:2] != ("velocity", "pressure"):
            raise ValueError("mixed fields must start with velocity and pressure")
        if len(set(self.names)) != len(self.names) or any(not name for name in self.names):
            raise ValueError("mixed field names must be nonempty and unique")

    def index(self, name: str) -> int:
        """Component index of the named field."""
        return self.names.index(name)

    def bind(self, components: Sequence[_T]) -> dict[str, _T]:
        """Name function, expression or test components, rejecting a mismatched space."""
        return dict(zip(self.names, components, strict=True))


def flow_layout(
    stabilization: StabilizationOptions, thermal: bool = False, passive_names: tuple[str, ...] = ()
) -> FieldLayout:
    """Named physical and projection fields of the configured method."""
    scalars = ("temperature", *passive_names) if thermal else ()
    names = ("velocity", "pressure", *scalars)
    if stabilization.method == "vms" and stabilization.subscale_space == "oss":
        names += ("xi", *(f"eta_{name}" for name in scalars))
        if stabilization.pressure_subscale:
            names += ("zeta",)
    return FieldLayout(names)


def mixed_space(
    mesh: MeshGeometry,
    velocity_degree: int,
    pressure_degree: int,
    temperature_degree: int | None = None,
    passive_scalars: int = 0,
) -> WithGeometry:
    """Physical mixed space with optional temperature and additional scalar components."""
    velocity = fd.VectorFunctionSpace(mesh, "CG", velocity_degree)
    pressure = fd.FunctionSpace(mesh, "CG", pressure_degree)
    spaces = [velocity, pressure]
    if temperature_degree is not None:
        spaces.extend(
            fd.FunctionSpace(mesh, "CG", temperature_degree) for _ in range(1 + passive_scalars)
        )
    return fd.MixedFunctionSpace(spaces)


def oss_space(
    mesh: MeshGeometry, degree: int, pressure_projection: bool, scalar_names: tuple[str, ...] = ()
) -> WithGeometry:
    """Physical fields and unconstrained OSS projections, including optional pressure projection."""
    velocity = fd.VectorFunctionSpace(mesh, "CG", degree)
    pressure = fd.FunctionSpace(mesh, "CG", degree)
    scalars = [fd.FunctionSpace(mesh, "CG", degree) for _ in scalar_names]
    spaces = [velocity, pressure, *scalars, velocity, *scalars]
    if pressure_projection:
        spaces.append(pressure)
    return fd.MixedFunctionSpace(spaces)


def constant_pressure_null_space(space: WithGeometry) -> fd.MixedVectorSpaceBasis:
    """Null space of W spanned by the constant pressure, its component 1 (plan §8, §11)."""
    constant = fd.VectorSpaceBasis(constant=True, comm=space.comm)
    bases = [constant if index == 1 else space.sub(index) for index in range(len(space))]
    return fd.MixedVectorSpaceBasis(space, bases)


def subscale_space(mesh: MeshGeometry, degree: int) -> WithGeometry:
    """Vector space of values at the points of the dx rule of the given degree."""
    return fd.VectorFunctionSpace(mesh, quadrature_element(mesh.ufl_cell(), degree))
