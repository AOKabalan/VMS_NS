"""Skew symmetry of the convective form and its outflow energy flux (plan U3)."""

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.discretization.quadrature import (
    boundary_measure,
    quadrature_degree,
    volume_measure,
)
from navier_stokes_vms.formulation.convective_forms import convection_form
from navier_stokes_vms.problem.boundary_conditions import named_boundary_measure
from navier_stokes_vms.problem.mesh_loading import rectangle

RNG = np.random.default_rng(20261021)


def random_function(space: fd.functionspaceimpl.WithGeometry) -> fd.Function:
    function = fd.Function(space)
    function.dat.data[:] = RNG.standard_normal(function.dat.data.shape)
    return function


def quadratic_value(matrix: fd.Matrix, function: fd.Function) -> tuple[float, float]:
    with function.dat.vec_ro as vector:
        product = vector.duplicate()
        matrix.petscmat.mult(vector, product)
        return vector.dot(product), vector.norm() * product.norm()


@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("dimension", [2, 3])
def test_convective_form_is_skew_on_a_periodic_mesh(dimension: int, degree: int) -> None:
    mesh = fd.PeriodicUnitSquareMesh(4, 4) if dimension == 2 else fd.PeriodicUnitCubeMesh(2, 2, 2)
    space = fd.VectorFunctionSpace(mesh, "CG", degree)
    a = random_function(space)
    dx = volume_measure(mesh, quadrature_degree(degree))
    trial, test = fd.TrialFunction(space), fd.TestFunction(space)
    matrix = fd.assemble(convection_form(a, a, trial, test, dx, None))
    value, scale = quadratic_value(matrix, random_function(space))
    assert scale > 1.0
    assert abs(value) < 1e-12 * scale


@pytest.mark.parametrize("degree", [1, 2])
def test_outflow_term_gives_the_boundary_energy_flux(degree: int) -> None:
    mesh, tags = rectangle((4, 3), (0.0, 0.0), (1.0, 1.0))
    space = fd.VectorFunctionSpace(mesh, "CG", degree)
    q = quadrature_degree(degree)
    outflow = named_boundary_measure(boundary_measure(mesh, q), ["x_max", "y_min"], tags)
    assert outflow.metadata() == boundary_measure(mesh, q).metadata()
    u = random_function(space)
    trial, test = fd.TrialFunction(space), fd.TestFunction(space)
    matrix = fd.assemble(convection_form(u, u, trial, test, volume_measure(mesh, q), outflow))
    value, _ = quadratic_value(matrix, u)
    n = ufl.FacetNormal(mesh)
    flux = fd.assemble(ufl.dot(u, n) * ufl.inner(u, u) / 2 * outflow)
    assert abs(flux) > 1e-3
    assert value == pytest.approx(flux, rel=1e-12)
