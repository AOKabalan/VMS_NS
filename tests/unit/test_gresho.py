"""Gresho vortex: reference fields, case data, and finite frozen-τ steps at rest (spec: Gresho).

The resting exterior has a = 0, where τ_m = ∞ at ν = 0 and |a| has no derivative (owner decision
K-D1). The case runs at the authorized fallback ν = 1e-6 because ν = 0 makes the frozen-τ OSS
Jacobian singular there.
"""

import csv
import logging
import math
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl
from firedrake.exceptions import ConvergenceError
from firedrake.petsc import PETSc

from navier_stokes_vms.cases import build_problem, gresho
from navier_stokes_vms.config import (
    Config,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
    ViscosityOptions,
)
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.problem.boundary_conditions import Dirichlet
from navier_stokes_vms.problem.problem_definition import (
    exact_solution_forcing,
    uses_pressure_null_space,
)
from navier_stokes_vms.solver.time_stepper import TimeStepper

LOGGER = logging.getLogger(__name__)

DT = 0.005
DIAGNOSTICS = (
    "kinetic_energy",
    "subscale_kinetic_energy",
    "linear_momentum",
    "angular_momentum",
    "energy_budget",
    "linear_momentum_defect",
    "angular_momentum_defect",
    "convective_energy_defect",
    "force",
    "velocity_error_l2",
    "pressure_error_l2",
    "divergence_l2",
)
STALLS = (
    "frozen-τ Newton stalls in the second outer iteration for nonlinear dynamic ASGS at "
    "ν ≤ 1e-4, including ν = 0 and the fallback; open owner question in the progress log"
)


def gresho_config(
    cells: int = 8,
    space: str = "asgs",
    splitting: str = "linear",
    newton_tau: str = "frozen",
    steps: int = 2,
) -> Config:
    """The run-matrix settings on a tiny mesh, with room for the linear frozen-τ contraction."""
    return Config(
        case="gresho",
        run_name="gresho",
        mesh=MeshOptions(cells=(cells, cells)),
        time=TimeOptions(scheme="crank_nicolson", dt=DT, t_end=steps * DT),
        stabilization=StabilizationOptions(subscale_space=space, scale_splitting=splitting),
        viscosity=ViscosityOptions(form="stress"),
        solver=SolverOptions(
            newton_tau=newton_tau,
            inner_rtol=1e-12,
            inner_max_iterations=100,
            outer_rtol=1e-11,
            outer_max_iterations=100,
            subscale_rtol=1e-13,
        ),
        output=OutputOptions(vtk_every=0, diagnostics=DIAGNOSTICS),
    )


def jacobian(stepper: TimeStepper) -> PETSc.Mat:
    """The Newton Jacobian at the current state, with the boundary rows."""
    problem = stepper._step_solver._inner._newton._problem
    return fd.assemble(problem.J, bcs=problem.bcs, mat_type="aij").petscmat


def jacobian_at_rest(stepper: TimeStepper) -> PETSc.Mat:
    """The Newton Jacobian at the interpolated initial state, first solve level."""
    stepper._initialise()
    stepper._set_level(1, step_coefficients("crank_nicolson", "lower_order", 1))
    return jacobian(stepper)


def entries(matrix: PETSc.Mat) -> np.ndarray:
    return matrix.getValuesCSR()[2]


def test_reference_fields_solve_the_steady_euler_equations() -> None:
    problem = build_problem(gresho_config(cells=16))
    x = fd.SpatialCoordinate(problem.mesh)
    dx = fd.dx(domain=problem.mesh, degree=8)
    t = fd.Constant(0.3)
    u = gresho.velocity(x, t)
    gradient = fd.assemble(ufl.inner(ufl.grad(u), ufl.grad(u)) * dx)
    assert fd.assemble(ufl.div(u) ** 2 * dx) <= 1e-26 * gradient
    euler = exact_solution_forcing(gresho.velocity, gresho.pressure, 0.0)(x, t)
    convection = ufl.dot(ufl.grad(u), u)
    scale = fd.assemble(ufl.inner(convection, convection) * dx)
    assert scale > 1.0
    assert fd.assemble(ufl.inner(euler, euler) * dx) <= 1e-26 * scale
    assert problem.forcing(x, t) == ufl.zero(2)
    later = gresho.velocity(x, fd.Constant(7.0)) - u
    assert fd.assemble(ufl.inner(later, later) * dx) == 0.0
    r = ufl.sqrt(ufl.dot(x, x))
    expected = ufl.conditional(r < 0.2, 10.0, ufl.conditional(r < 0.4, 2 / r - 10, 0.0))
    vorticity = u[1].dx(0) - u[0].dx(1)
    assert fd.assemble((vorticity - expected) ** 2 * dx) <= 1e-26 * gradient


@pytest.mark.parametrize("radius", [0.2, 0.4])
def test_velocity_and_pressure_are_continuous_at_the_interfaces(radius: float) -> None:
    mesh = fd.UnitSquareMesh(1, 1)
    x = fd.SpatialCoordinate(mesh)
    u, p = gresho.velocity(x, 0.0), gresho.pressure(x, 0.0)
    swirl = 5 * radius if radius == 0.2 else 0.0
    for angle in np.linspace(0, 2 * math.pi, 7):
        values = []
        for factor in (1 - 1e-12, 1 + 1e-12):
            point = (radius * factor * math.cos(angle), radius * factor * math.sin(angle))
            values.append((p(point), u[0](point), u[1](point)))
        (p_in, ux_in, uy_in), (p_out, ux_out, uy_out) = values
        assert p_out == pytest.approx(p_in, abs=1e-10)
        assert (ux_out, uy_out) == pytest.approx((ux_in, uy_in), abs=1e-10)
        assert (ux_in, uy_in) == pytest.approx(
            (-swirl * math.sin(angle), swirl * math.cos(angle)), abs=1e-10
        )
    ring_constant = -12.5 * 0.4**2 + 20 * 0.4 - 4 * math.log(0.4)
    assert p((0.0, 0.0)) == pytest.approx(ring_constant - 4 + 4 * math.log(0.2), rel=1e-15)
    assert p((0.3, 0.35)) == 0.0


def test_case_declares_walls_frame_and_fallback_viscosity() -> None:
    problem = build_problem(gresho_config(cells=20))
    names = ("x_min", "x_max", "y_min", "y_max")
    assert sorted(problem.boundary_conditions) == sorted(names)
    assert problem.force_boundaries == names
    assert all(isinstance(bc, Dirichlet) for bc in problem.boundary_conditions.values())
    x = fd.SpatialCoordinate(problem.mesh)
    assert all(bc.value(x, 0.0) == ufl.zero(2) for bc in problem.boundary_conditions.values())
    assert uses_pressure_null_space(problem)
    assert problem.viscosity == 1e-6
    # The 20×20 cell edges lie on the frame boundary, so the indicator integrates exactly.
    frame_area = fd.assemble(problem.monitor_region(x) * fd.dx(domain=problem.mesh, degree=3))
    assert frame_area == pytest.approx(1 - 0.9**2, rel=1e-13)
    with pytest.raises(ValueError, match="case_parameters: viscosity not used by case 'gresho'"):
        build_problem(replace(gresho_config(), case_parameters={"viscosity": 0.1}))


def test_zero_viscosity_makes_the_frozen_oss_jacobian_singular_at_rest(tmp_path: Path) -> None:
    """Why the case uses the fallback: τ_c = 0 at rest empties every ζ row supported there."""
    config = gresho_config(space="oss")
    counts = {}
    for viscosity in (0.0, gresho.VISCOSITY):
        problem = replace(build_problem(config), viscosity=viscosity)
        stepper = TimeStepper(config, problem, tmp_path / str(viscosity), LOGGER)
        matrix = jacobian_at_rest(stepper)
        assert np.all(np.isfinite(entries(matrix)))
        dense = matrix.convert("dense").getDenseArray()
        space = stepper._solution.function_space()
        first_zeta = sum(space.sub(i).dim() for i in range(3))
        counts[viscosity] = int(np.sum(np.all(dense[first_zeta:] == 0, axis=1)))
    assert counts[0.0] > 0
    assert counts[gresho.VISCOSITY] == 0


@pytest.mark.parametrize(
    "space,splitting",
    [
        ("asgs", "linear"),
        pytest.param(
            "asgs",
            "nonlinear",
            marks=pytest.mark.xfail(raises=ConvergenceError, strict=True, reason=STALLS),
        ),
        ("oss", "linear"),
        ("oss", "nonlinear"),
    ],
)
def test_frozen_tau_step_from_rest_is_finite_K_D1(
    tmp_path: Path, space: str, splitting: str
) -> None:
    """Residual, Jacobian, solve and every diagnostic stay finite through a Crank–Nicolson step."""
    config = gresho_config(space=space, splitting=splitting)
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    assert np.all(np.isfinite(entries(jacobian_at_rest(stepper))))
    assert math.isfinite(stepper._step_solver._inner.residual_norm())
    stepper.run()
    assert np.all(np.isfinite(entries(jacobian(stepper))))
    with (tmp_path / "time_steps.csv").open() as stream:
        assert [row["converged"] for row in csv.DictReader(stream)] == ["True", "True"]
    with (tmp_path / "diagnostics.csv").open() as stream:
        startup, midpoint = list(csv.DictReader(stream))
    assert startup["energy_budget_residual"] == ""
    for name, value in midpoint.items():
        assert value != "", name
        assert math.isfinite(float(value)), name


def test_exact_tau_fails_at_rest_K_D1(tmp_path: Path) -> None:
    config = gresho_config(newton_tau="exact")
    stepper = TimeStepper(config, build_problem(config), tmp_path / "assembly", LOGGER)
    assert not np.all(np.isfinite(entries(jacobian_at_rest(stepper))))
    (tmp_path / "run").mkdir()
    stepper = TimeStepper(config, build_problem(config), tmp_path / "run", LOGGER)
    with pytest.raises(ConvergenceError, match="solver.newton_tau = 'exact'"):
        stepper.run()
