"""The single quadrature rule: dx, ds, the matching quadrature element and its weight check."""

import firedrake as fd
import numpy as np
import ufl
from finat.element_factory import as_fiat_cell
from finat.quadrature import make_quadrature
from finat.ufl import FiniteElement
from firedrake.mesh import MeshGeometry
from ufl.measure import Measure

SCHEME = "default"


def quadrature_degree(velocity_degree: int) -> int:
    """q = 3k for velocity degree k (spec: Boundary conditions and quadrature)."""
    return 3 * velocity_degree


def volume_measure(mesh: MeshGeometry, degree: int) -> Measure:
    """dx with the rule of the given degree, after checking that its weights are positive."""
    check_positive_weights(mesh.ufl_cell(), degree)
    return fd.dx(domain=mesh, degree=degree, scheme=SCHEME)


def boundary_measure(mesh: MeshGeometry, degree: int) -> Measure:
    """ds with the facet rule of the given degree."""
    return fd.ds(domain=mesh, degree=degree, scheme=SCHEME)


def quadrature_element(cell: ufl.Cell, degree: int) -> FiniteElement:
    """Scalar element whose nodes are the points of the dx rule of the given degree."""
    return FiniteElement("Quadrature", cell, degree, quad_scheme=SCHEME)


def reference_rule(cell: ufl.Cell, degree: int) -> tuple[np.ndarray, np.ndarray]:
    """Points and weights on the reference cell of the rule that dx and the element use."""
    rule = make_quadrature(as_fiat_cell(cell), degree, scheme=SCHEME)
    return np.asarray(rule.point_set.points), np.asarray(rule.weights)


def check_positive_weights(cell: ufl.Cell, degree: int) -> None:
    """Raise unless every weight is positive, as the τ-weighted mass matrices require."""
    _, weights = reference_rule(cell, degree)
    if weights.min() <= 0:
        raise ValueError(
            f"quadrature rule of degree {degree} on {cell.cellname} has a non-positive weight "
            f"{weights.min():.3e}; the τ-weighted projections need positive weights"
        )
