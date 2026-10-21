"""Built-in and gmsh meshes with their boundary name → tag tables (plan §8)."""

from collections.abc import Sequence
from pathlib import Path

import firedrake as fd
from firedrake.mesh import MeshGeometry

from navier_stokes_vms.config import Config

RECTANGLE_TAGS = {"x_min": 1, "x_max": 2, "y_min": 3, "y_max": 4}
BOX_TAGS = {**RECTANGLE_TAGS, "z_min": 5, "z_max": 6}


def rectangle(
    cells: Sequence[int], lower: Sequence[float], upper: Sequence[float]
) -> tuple[MeshGeometry, dict[str, int]]:
    """Triangle mesh of the rectangle [lower, upper]; boundaries x_min, x_max, y_min, y_max."""
    nx, ny = cells
    mesh = fd.RectangleMesh(nx, ny, upper[0] - lower[0], upper[1] - lower[1])
    _translate(mesh, lower)
    return mesh, dict(RECTANGLE_TAGS)


def periodic_rectangle(
    cells: Sequence[int], lower: Sequence[float], upper: Sequence[float]
) -> tuple[MeshGeometry, dict[str, int]]:
    """Triangle mesh of [lower, upper], periodic in both directions, so without boundaries."""
    nx, ny = cells
    mesh = fd.PeriodicRectangleMesh(nx, ny, upper[0] - lower[0], upper[1] - lower[1])
    _translate(mesh, lower)
    return mesh, {}


def box(
    cells: Sequence[int], lower: Sequence[float], upper: Sequence[float]
) -> tuple[MeshGeometry, dict[str, int]]:
    """Tetrahedral mesh of the box [lower, upper]; boundaries x_min, x_max, ..., z_max."""
    nx, ny, nz = cells
    lengths = [high - low for low, high in zip(lower, upper, strict=True)]
    mesh = fd.BoxMesh(nx, ny, nz, *lengths)
    _translate(mesh, lower)
    return mesh, dict(BOX_TAGS)


def gmsh_mesh(path: Path) -> tuple[MeshGeometry, dict[str, int]]:
    """Mesh from a .msh file; boundary names are its $PhysicalNames entries of facet dimension."""
    mesh = fd.Mesh(str(path))
    return mesh, _physical_names(path, dimension=mesh.topological_dimension - 1)


def builtin_cells(config: Config, default: tuple[int, ...]) -> tuple[int, ...]:
    """Cells per direction for a case that has only a built-in mesh: config override or default."""
    if config.mesh.source == "gmsh" or config.mesh.file is not None:
        raise ValueError(
            f"mesh.source and mesh.file: case {config.case!r} supports only builtin meshes"
        )
    cells = config.mesh.cells or default
    if len(cells) != len(default):
        raise ValueError(
            f"mesh.cells = {list(cells)}: case {config.case!r} needs {len(default)} entries"
        )
    return cells


def _translate(mesh: MeshGeometry, lower: Sequence[float]) -> None:
    mesh.coordinates.dat.data[:] += lower


def _physical_names(path: Path, dimension: int) -> dict[str, int]:
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if line.strip() == "$PhysicalNames":
                break
        else:
            return {}
        count = int(next(stream))
        entries = [next(stream).split(maxsplit=2) for _ in range(count)]
    return {
        name.strip().strip('"'): int(tag) for dim, tag, name in entries if int(dim) == dimension
    }
