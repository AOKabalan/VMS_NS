"""Eigenvalue checks of static ASGS and OSS dissipation, and the cell-diameter regression (U5)."""

import logging
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl
from ufl.form import Form

from navier_stokes_vms.discretization.function_spaces import mixed_space, oss_space
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.formulation.galerkin import SolveLevel
from navier_stokes_vms.formulation.oss_projections import projection_equations
from navier_stokes_vms.formulation.stabilization_parameters import stabilization_parameters
from navier_stokes_vms.formulation.vms import Switches, subscale_residual, vms_terms

LOGGER = logging.getLogger(__name__)

STATIC_LINEAR = Switches(is_asgs=True, is_dynamic=False, is_nonlinear=False)
STATIC_LINEAR_OSS = Switches(is_asgs=False, is_dynamic=False, is_nonlinear=False)


def free_matrix(
    form: Form, field: fd.Function, bc: fd.DirichletBC
) -> tuple[np.ndarray, np.ndarray]:
    marker = fd.Function(field.function_space()).assign(1.0)
    bc.zero(marker)
    with marker.dat.vec_ro as vector:
        free = np.flatnonzero(vector.getArray(readonly=True))
    matrix = fd.assemble(form).petscmat.convert("dense")
    restricted = matrix.getDenseArray()[np.ix_(free, free)].copy()
    return (restricted + restricted.T) / 2, free


@pytest.mark.parametrize("k", [1, 2])
@pytest.mark.parametrize("mesh_kind", ["left", "crossed", "box", "cylinder"])
@pytest.mark.parametrize("nu", [10.0, 1e-4], ids=["viscous", "convective"])
def test_static_asgs_energy_has_no_negative_eigenvalues(k: int, mesh_kind: str, nu: float) -> None:
    mesh = (
        fd.Mesh(str(Path(__file__).resolve().parents[2] / "meshes/cylinder_dfg_2d_coarse.msh"))
        if mesh_kind == "cylinder"
        else fd.UnitCubeMesh(2, 2, 2)
        if mesh_kind == "box"
        else fd.UnitSquareMesh(4, 4, diagonal=mesh_kind)
    )
    space = mixed_space(mesh, k, k)
    state = fd.Function(space)
    a = fd.Function(space.sub(0).collapse())
    rng = np.random.default_rng(51)
    a.dat.data[:] = rng.uniform(0.5, 1.5, a.dat.data.shape)
    u, p = fd.TrialFunctions(space)
    v, q = fd.TestFunctions(space)
    dx = volume_measure(mesh, 3 * k)
    level = SolveLevel(
        fd.Constant(0.0),
        fd.Constant(1.0),
        fd.Function(a.function_space()),
        ufl.zero(mesh.geometric_dimension),
    )
    tau_m, tau_c = stabilization_parameters(a, nu, k, 4.0, 2.0)
    energy = nu * ufl.inner(ufl.grad(u), ufl.grad(v)) * dx
    energy += vms_terms(u, p, v, q, a, nu, level, tau_m, tau_c, dx, STATIC_LINEAR, None)
    matrix, _ = free_matrix(fd.lhs(energy), state, fd.DirichletBC(space.sub(0), 0, "on_boundary"))
    eigenvalues = np.linalg.eigvalsh(matrix)
    smallest, largest = eigenvalues[0], eigenvalues[-1]
    LOGGER.info("U5 P%d/P%d %s ν=%g: λ_min=%g, λ_max=%g", k, k, mesh_kind, nu, smallest, largest)
    assert largest > 0
    assert smallest >= -1e-12 * largest


@pytest.mark.parametrize("k", [1, 2])
@pytest.mark.parametrize("mesh_kind", ["left", "crossed", "box", "cylinder"])
@pytest.mark.parametrize("nu", [10.0, 1e-4], ids=["viscous", "convective"])
def test_static_oss_energy_has_no_negative_eigenvalues(k: int, mesh_kind: str, nu: float) -> None:
    """E(u, p; v, q) with (ξ_h, ζ_h) eliminated: the Schur complement of the coupled form.

    The projections are linear in (u, p), so eliminating them from the coupled bilinear form gives
    the same matrix as assembling E column by column with each column's projections.
    """
    mesh = (
        fd.Mesh(str(Path(__file__).resolve().parents[2] / "meshes/cylinder_dfg_2d_coarse.msh"))
        if mesh_kind == "cylinder"
        else fd.UnitCubeMesh(2, 2, 2)
        if mesh_kind == "box"
        else fd.UnitSquareMesh(4, 4, diagonal=mesh_kind)
    )
    space = oss_space(mesh, k, pressure_projection=True)
    a = fd.Function(space.sub(0).collapse())
    rng = np.random.default_rng(51)
    a.dat.data[:] = rng.uniform(0.5, 1.5, a.dat.data.shape)
    u, p, xi, zeta = fd.TrialFunctions(space)
    v, q, w, r = fd.TestFunctions(space)
    dx = volume_measure(mesh, 3 * k)
    level = SolveLevel(
        fd.Constant(0.0),
        fd.Constant(1.0),
        fd.Function(a.function_space()),
        ufl.zero(mesh.geometric_dimension),
    )
    tau_m, tau_c = stabilization_parameters(a, nu, k, 4.0, 2.0)
    chi = STATIC_LINEAR_OSS
    residual = subscale_residual(u, p, a, nu, level, chi, None)
    form = nu * ufl.inner(ufl.grad(u), ufl.grad(v)) * dx
    form += vms_terms(u, p, v, q, a, nu, level, tau_m, tau_c, dx, chi, None, xi, zeta)
    form += projection_equations(xi, zeta, w, r, residual, ufl.div(u), tau_m, tau_c, dx)
    marker = fd.Function(space).assign(1.0)
    fd.DirichletBC(space.sub(0), 0, "on_boundary").zero(marker)
    physical = space.sub(0).dim() + space.sub(1).dim()
    with marker.dat.vec_ro as vector:
        free = np.flatnonzero(vector.getArray(readonly=True)[:physical])
    projected = np.arange(physical, space.dim())
    matrix = fd.assemble(fd.lhs(form)).petscmat.convert("dense").getDenseArray()
    coupled = matrix[np.ix_(free, projected)]
    mass = matrix[np.ix_(projected, projected)]
    schur = matrix[np.ix_(free, free)] - coupled @ np.linalg.solve(
        mass, matrix[np.ix_(projected, free)]
    )
    eigenvalues = np.linalg.eigvalsh((schur + schur.T) / 2)
    smallest, largest = eigenvalues[0], eigenvalues[-1]
    LOGGER.info(
        "U5 OSS P%d/P%d %s ν=%g: λ_min=%g, λ_max=%g", k, k, mesh_kind, nu, smallest, largest
    )
    assert largest > 0
    assert smallest >= -1e-12 * largest


def test_cell_diameter_detects_the_p2_counterexample() -> None:
    mesh = fd.UnitSquareMesh(4, 4, diagonal="left")
    space = fd.VectorFunctionSpace(mesh, "CG", 2)
    u, v = fd.TrialFunction(space), fd.TestFunction(space)
    dx = volume_measure(mesh, 6)
    h = ufl.CellDiameter(mesh)
    energy_at_rest = (
        ufl.inner(ufl.grad(u), ufl.grad(v))
        + ufl.div(u) * ufl.div(v)
        - h**2 / 64 * ufl.inner(ufl.div(ufl.grad(u)), ufl.div(ufl.grad(v)))
    ) * dx
    velocity = fd.Function(space)
    matrix, free = free_matrix(energy_at_rest, velocity, fd.DirichletBC(space, 0, "on_boundary"))
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    assert eigenvalues[0] < -1e-12 * eigenvalues[-1]
    with velocity.dat.vec as vector:
        vector.getArray()[free] = 1e-4 * eigenvectors[:, 0]
    pressure = fd.Function(fd.FunctionSpace(mesh, "CG", 2))
    level = SolveLevel(fd.Constant(0.0), fd.Constant(1.0), fd.Function(space), ufl.zero(2))
    tau_m, tau_c = stabilization_parameters(velocity, 1.0, 2, 4.0, 2.0)
    diameter_tau_m = 1 / (64 / h**2 + 4 * ufl.sqrt(ufl.dot(velocity, velocity)) / h)
    diameter_tau_c = h**2 / (64 * diameter_tau_m)
    viscous = ufl.inner(ufl.grad(velocity), ufl.grad(velocity)) * dx
    energies = []
    for momentum, continuity in ((diameter_tau_m, diameter_tau_c), (tau_m, tau_c)):
        stabilization = vms_terms(
            velocity,
            pressure,
            velocity,
            pressure,
            velocity,
            1.0,
            level,
            momentum,
            continuity,
            dx,
            STATIC_LINEAR,
            None,
        )
        energies.append(fd.assemble(viscous + stabilization))
    LOGGER.info("U5 P2 counterexample: diameter energy=%g, shortest-edge energy=%g", *energies)
    assert energies[0] < -1e-12 * fd.assemble(viscous)
    assert energies[1] > 0
