"""The 3D exact formula is solenoidal and its symbolically derived forcing vanishes."""

import logging
import math
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import Config, MeshOptions, TimeOptions
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.problem.boundary_conditions import Dirichlet, Outflow
from navier_stokes_vms.solver.time_stepper import TimeStepper


@pytest.mark.parametrize("time", [0.0, 0.17, 1.0])
def test_ethier_steinman_has_zero_forcing(time: float) -> None:
    config = Config(case="ethier_steinman_3d", run_name="exact", mesh=MeshOptions(cells=(2, 2, 2)))
    problem = build_problem(config)
    x, t = fd.SpatialCoordinate(problem.mesh), fd.Constant(time)
    forcing = problem.forcing(x, t)
    velocity = problem.exact_velocity(x, t)
    dx = fd.dx(domain=problem.mesh, degree=6)
    assert fd.assemble(ufl.inner(forcing, forcing) * dx) ** 0.5 < 1e-13
    assert fd.assemble(ufl.div(velocity) ** 2 * dx) ** 0.5 < 1e-13
    assert len(problem.boundary_conditions) == len(problem.force_boundaries) == 6


@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("cells", [4, 8])
def test_ethier_dirichlet_interpolant_has_compatible_boundary_flux(degree: int, cells: int) -> None:
    config = Config(
        case="ethier_steinman_3d", run_name="flux", mesh=MeshOptions(cells=(cells,) * 3)
    )
    problem = build_problem(config)
    x = fd.SpatialCoordinate(problem.mesh)
    velocity = fd.Function(fd.VectorFunctionSpace(problem.mesh, "CG", degree)).interpolate(
        problem.exact_velocity(x, fd.Constant(0.01))
    )
    assert abs(fd.assemble(ufl.div(velocity) * fd.dx(degree=3 * degree))) < 1e-13


@pytest.mark.parametrize("d", [math.pi / 2, math.pi], ids=["incompatible", "compatible"])
@pytest.mark.parametrize("scheme", ["bdf2", "crank_nicolson"])
def test_time_stepper_checks_ethier_dirichlet_flux(tmp_path: Path, d: float, scheme: str) -> None:
    config = Config(
        case="ethier_steinman_3d",
        run_name="flux",
        mesh=MeshOptions(cells=(4, 4, 4)),
        time=TimeOptions(scheme=scheme, dt=0.01, t_end=0.03),
    )
    problem = build_problem(config)
    a = math.pi / 4

    def velocity(x: ufl.core.expr.Expr, t: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        return (
            -a
            * ufl.exp(-(d**2) * t)
            * ufl.as_vector(
                (
                    ufl.exp(a * x[0]) * ufl.sin(a * x[1] + d * x[2])
                    + ufl.exp(a * x[2]) * ufl.cos(a * x[0] + d * x[1]),
                    ufl.exp(a * x[1]) * ufl.sin(a * x[2] + d * x[0])
                    + ufl.exp(a * x[0]) * ufl.cos(a * x[1] + d * x[2]),
                    ufl.exp(a * x[2]) * ufl.sin(a * x[0] + d * x[1])
                    + ufl.exp(a * x[1]) * ufl.cos(a * x[2] + d * x[0]),
                )
            )
        )

    problem = replace(
        problem, boundary_conditions={name: Dirichlet(velocity) for name in problem.boundary_tags}
    )
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger(__name__))
    for step in (1, 2, 3):
        coefficients = step_coefficients(scheme, "lower_order", step)
        if d == math.pi / 2:
            with pytest.raises(ValueError, match=f"step {step}, solve time .*discretely flux-free"):
                stepper._set_level(step, coefficients)
        else:
            stepper._set_level(step, coefficients)


@pytest.mark.parametrize("outflow", [False, True])
def test_flux_guard_updates_boundary_data_at_each_solve_level(
    tmp_path: Path, outflow: bool
) -> None:
    config = Config(
        case="manufactured_temporal",
        run_name="flux",
        mesh=MeshOptions(cells=(2, 2)),
        time=TimeOptions(scheme="crank_nicolson", dt=0.1, t_end=0.2),
    )
    problem = build_problem(config)
    conditions = {name: Dirichlet(lambda x, t: t * (t - 0.1) * x) for name in problem.boundary_tags}
    if outflow:
        conditions["x_max"] = Outflow()
    problem = replace(problem, boundary_conditions=conditions)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger(__name__))
    stepper._set_level(1, step_coefficients("crank_nicolson", "lower_order", 1))
    coefficients = step_coefficients("crank_nicolson", "lower_order", 2)
    if outflow:
        stepper._set_level(2, coefficients)
    else:
        with pytest.raises(ValueError, match="step 2, solve time 0.15:.*discretely flux-free"):
            stepper._set_level(2, coefficients)
        x = fd.SpatialCoordinate(problem.mesh)
        error = stepper._boundary_velocity - 0.01 * x
        assert fd.assemble(ufl.inner(error, error) * stepper._diagnostic_state.ds) < 1e-28
