"""Weighted OSS projections: identity on the FE space and orthogonality of ũ (plan U6)."""

import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.formulation.galerkin import SolveLevel
from navier_stokes_vms.formulation.oss_projections import projection_equations
from navier_stokes_vms.formulation.stabilization_parameters import (
    dynamic_parameter,
    stabilization_parameters,
)
from navier_stokes_vms.formulation.vms import subscale_residual, switches
from navier_stokes_vms.solver.time_stepper import TimeStepper

LOGGER = logging.getLogger(__name__)
DIRECT = {"snes_type": "ksponly", "ksp_type": "preonly", "pc_type": "lu"}


@pytest.mark.parametrize("k", [1, 2])
@pytest.mark.parametrize("dimension", [2, 3])
def test_weighted_projection_returns_fe_functions_unchanged(k: int, dimension: int) -> None:
    mesh = fd.UnitSquareMesh(3, 3) if dimension == 2 else fd.UnitCubeMesh(2, 2, 2)
    velocity = fd.VectorFunctionSpace(mesh, "CG", k)
    pressure = fd.FunctionSpace(mesh, "CG", k)
    rng = np.random.default_rng(6)
    target_velocity, target_pressure, advector = (
        fd.Function(velocity),
        fd.Function(pressure),
        fd.Function(velocity),
    )
    # Random values everywhere, boundary nodes included.
    for field in (target_velocity, target_pressure):
        field.dat.data[:] = rng.standard_normal(field.dat.data.shape)
    advector.dat.data[:] = rng.uniform(-2.0, 2.0, advector.dat.data.shape)
    tau_m, tau_c = stabilization_parameters(advector, 0.01, k, 4.0, 2.0)
    tau_t = dynamic_parameter(tau_m, fd.Constant(1.5), fd.Constant(0.1))
    dx = volume_measure(mesh, 3 * k)
    projected = fd.Function(velocity * pressure)
    xi, zeta = ufl.split(projected)
    w, r = fd.TestFunctions(projected.function_space())
    equations = projection_equations(
        xi, zeta, w, r, target_velocity, target_pressure, tau_t, tau_c, dx
    )
    fd.solve(equations == 0, projected, solver_parameters=DIRECT)
    weights = [fd.assemble(weight * dx) for weight in (tau_t, tau_c)]
    assert max(weights) / min(weights) > 1  # the weights vary
    for value, target in zip(
        projected.subfunctions, (target_velocity, target_pressure), strict=True
    ):
        scale = np.abs(target.dat.data_ro).max()
        assert np.abs(value.dat.data_ro - target.dat.data_ro).max() <= 1e-12 * scale


def oss_config(k: int, dynamics: str, splitting: str) -> Config:
    return Config(
        case="kovasznay",
        run_name="u6",
        mesh=MeshOptions(cells=(12, 16)),
        discretization=DiscretizationOptions(velocity_degree=k, pressure_degree=k),
        time=TimeOptions(scheme="bdf2", dt=0.5, t_end=1.0),
        stabilization=StabilizationOptions(
            subscale_space="oss", subscale_dynamics=dynamics, scale_splitting=splitting
        ),
        solver=replace(SolverOptions(), inner_rtol=1e-12, inner_atol=1e-14),
        output=OutputOptions(vtk_every=0),
    )


@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
@pytest.mark.parametrize("k", [1, 2])
def test_converged_subscale_is_orthogonal_to_the_velocity_space(
    tmp_path: Path, k: int, dynamics: str, splitting: str
) -> None:
    """After step 2 of BDF2, whose dynamic ũ_hist = 2ũ¹ − ũ⁰/2 is non-zero."""
    config = oss_config(k, dynamics, splitting)
    problem = build_problem(config)
    stepper = TimeStepper(config, problem, tmp_path, LOGGER)
    stepper.run()
    chi = switches(config.stabilization)
    u, p, xi, _zeta = ufl.split(stepper._solution)
    x = fd.SpatialCoordinate(problem.mesh)
    dt = fd.Constant(config.time.dt)
    level = SolveLevel(
        stepper._alpha0, dt, stepper._velocity_history, problem.forcing(x, stepper._t_new)
    )
    a = u + stepper._iterate.value if chi.is_nonlinear else u
    tau_m, _ = stabilization_parameters(a, problem.viscosity, k, 4.0, 2.0)
    tau_t = dynamic_parameter(tau_m, stepper._alpha0, dt) if chi.is_dynamic else tau_m
    history = None
    if chi.is_dynamic:
        history = stepper._subscales.history
        assert fd.norm(history) > 1e-3 * fd.norm(stepper._subscales.latest)
    weighted = tau_t * subscale_residual(u, p, a, problem.viscosity, level, chi, history)
    w = fd.TestFunction(stepper.velocity.function_space())
    dx = stepper._dx
    scale = fd.assemble(ufl.inner(weighted, w) * dx).dat.norm
    subscales = {"expression": stepper._u_sub}
    if chi.is_dynamic:
        subscales["stored"] = stepper._subscales.latest
    for name, subscale in subscales.items():
        products = fd.assemble(ufl.inner(subscale, w) * dx).dat.norm
        LOGGER.info(
            "U6 P%d %s %s, %s ũ: |(ũ, w_h)| / |(τ_t X, w_h)| = %.3e",
            k,
            dynamics,
            splitting,
            name,
            products / scale,
        )
        assert products <= 1e-10 * scale, name
