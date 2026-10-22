"""VMS Jacobians (ASGS and coupled OSS), true residual logging, and the exact-τ singularity."""

import csv
import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl
from firedrake.exceptions import ConvergenceError

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
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.formulation.oss_projections import projection_equations
from navier_stokes_vms.formulation.stabilization_parameters import (
    dynamic_parameter,
    stabilization_parameters,
)
from navier_stokes_vms.formulation.vms import Switches, subscale_residual, vms_terms
from navier_stokes_vms.solver.time_stepper import TimeStepper

LOGGER = logging.getLogger(__name__)

STATIC_LINEAR = Switches(is_asgs=True, is_dynamic=False, is_nonlinear=False)


def asgs_config(
    degree: int,
    options: SolverOptions,
    dynamics: str = "static",
    splitting: str = "linear",
    space: str = "asgs",
) -> Config:
    return Config(
        case="manufactured_spatial",
        run_name="asgs",
        mesh=MeshOptions(cells=(4, 4)),
        discretization=DiscretizationOptions(velocity_degree=degree, pressure_degree=degree),
        time=TimeOptions(scheme="bdf1", dt=0.05, t_end=0.05),
        stabilization=StabilizationOptions(
            subscale_space=space, subscale_dynamics=dynamics, scale_splitting=splitting
        ),
        solver=options,
        output=OutputOptions(vtk_every=0),
    )


@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("newton_tau", ["frozen", "exact"])
def test_newton_jacobian_matches_finite_differences(
    tmp_path: Path, degree: int, newton_tau: str, dynamics: str, splitting: str, space: str
) -> None:
    """For OSS the unknown includes ξ_h, ζ_h and the residual their projection equations."""
    options = replace(SolverOptions(), newton_tau=newton_tau)
    config = asgs_config(degree, options, dynamics, splitting, space)
    problem = build_problem(config)
    stepper = TimeStepper(config, problem, tmp_path, LOGGER)
    stepper._initialise()
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    solution = stepper._solution
    mixed = solution.function_space()
    direction = fd.Function(mixed)
    rng = np.random.default_rng(56)
    for field in (*solution.subfunctions, *direction.subfunctions):
        field.dat.data[:] = rng.uniform(0.2, 1.0, field.dat.data.shape)
    base = solution.copy(deepcopy=True)
    frozen = base.subfunctions[0]
    u, p, *projections = ufl.split(solution)
    v, q, *projection_tests = fd.TestFunctions(mixed)
    level = SolveLevel(
        stepper._alpha0,
        fd.Constant(config.time.dt),
        stepper._velocity_history,
        problem.forcing(fd.SpatialCoordinate(problem.mesh), stepper._t_new),
    )
    a, a_tau = u, frozen if newton_tau == "frozen" else u
    if splitting == "nonlinear":
        u_sub_iterate = stepper._iterate.value
        u_sub_iterate.dat.data[:] = rng.uniform(-0.5, 0.5, u_sub_iterate.dat.data.shape)
        a, a_tau = a + u_sub_iterate, a_tau + u_sub_iterate
    tau_m, tau_c = stabilization_parameters(a_tau, problem.viscosity, degree, 4.0, 2.0)
    u_sub_history = None
    if dynamics == "dynamic":
        u_sub_history = stepper._subscales.history
        u_sub_history.dat.data[:] = rng.uniform(-1.0, 1.0, u_sub_history.dat.data.shape)
        tau_m = dynamic_parameter(tau_m, level.alpha0, level.dt)
    chi = Switches(space == "asgs", dynamics == "dynamic", splitting == "nonlinear")
    nu, dx = problem.viscosity, stepper._dx
    form = galerkin_residual(u, p, v, q, a, u, nu, level, dx, None)
    form += vms_terms(u, p, v, q, a, nu, level, tau_m, tau_c, dx, chi, u_sub_history, *projections)
    if space == "oss":
        xi, zeta = projections
        w, r = projection_tests
        residual = subscale_residual(u, p, a, nu, level, chi, u_sub_history)
        form += projection_equations(xi, zeta, w, r, residual, ufl.div(u), tau_m, tau_c, dx)
    jacobian = stepper._step_solver._inner._newton._problem.J
    action = fd.assemble(ufl.action(jacobian, direction))
    epsilon = 1e-6
    solution.assign(base + epsilon * direction)
    plus = fd.assemble(form)
    solution.assign(base - epsilon * direction)
    minus = fd.assemble(form)
    with plus.dat.vec_ro as vp, minus.dat.vec_ro as vm, action.dat.vec_ro as derivative:
        difference = vp.copy()
        difference.axpy(-1, vm)
        difference.scale(1 / (2 * epsilon))
        difference.axpy(-1, derivative)
        assert difference.norm() <= 1e-8 * derivative.norm()


@pytest.mark.parametrize(
    "linearization,newton_tau,warmup",
    [
        ("newton", "frozen", 0),
        ("newton", "exact", 0),
        ("picard", "frozen", 0),
        ("newton", "frozen", 2),
    ],
)
def test_each_iteration_logs_the_true_asgs_residual(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    linearization: str,
    newton_tau: str,
    warmup: int,
) -> None:
    options = replace(
        SolverOptions(),
        linearization=linearization,
        newton_tau=newton_tau,
        picard_iterations_before_newton=warmup,
    )
    config = asgs_config(1, options)
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    solver = stepper._step_solver._inner
    u, p = ufl.split(stepper._solution)
    v, q = fd.TestFunctions(stepper._solution.function_space())
    problem = stepper._problem
    level = SolveLevel(
        stepper._alpha0,
        fd.Constant(config.time.dt),
        stepper._velocity_history,
        problem.forcing(fd.SpatialCoordinate(problem.mesh), stepper._t_new),
    )
    tau_m, tau_c = stabilization_parameters(u, problem.viscosity, 1, 4.0, 2.0)
    reference = galerkin_residual(u, p, v, q, u, u, problem.viscosity, level, stepper._dx, None)
    reference += vms_terms(
        u, p, v, q, u, problem.viscosity, level, tau_m, tau_c, stepper._dx, STATIC_LINEAR, None
    )
    measured = []
    append = solver._append

    def checked_append(method: str, norm: float, update: float | None, linear: int | None) -> None:
        residual = fd.assemble(reference)
        for bc in solver._bcs:
            bc.zero(residual)
        with residual.dat.vec_ro as vector:
            measured.append((norm, vector.norm()))
        append(method, norm, update, linear)

    monkeypatch.setattr(solver, "_append", checked_append)
    stepper.run()
    for recorded, actual in measured:
        assert recorded == pytest.approx(actual, rel=1e-7, abs=1e-13)
    with (tmp_path / "inner_iterations.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert float(rows[-1]["residual_norm"]) == pytest.approx(
        solver.residual_norm(), rel=1e-7, abs=1e-13
    )
    assert float(rows[-1]["residual_norm"]) <= max(
        options.inner_atol, options.inner_rtol * float(rows[0]["residual_norm"])
    )


@pytest.mark.parametrize("newton_tau", ["frozen", "exact"])
def test_start_from_rest_handles_the_tau_derivative(tmp_path: Path, newton_tau: str) -> None:
    config = asgs_config(1, replace(SolverOptions(), newton_tau=newton_tau))
    problem = replace(build_problem(config), initial_velocity=lambda _x: ufl.zero(2))
    stepper = TimeStepper(config, problem, tmp_path, LOGGER)
    if newton_tau == "frozen":
        stepper.run()
        assert np.isfinite(stepper._step_solver._inner.residual_norm())
    else:
        with pytest.raises(ConvergenceError, match="solver.newton_tau = 'exact'"):
            stepper.run()
        assert (tmp_path / "fields" / "velocity.pvd").exists()
        with (tmp_path / "time_steps.csv").open() as stream:
            assert list(csv.DictReader(stream))[-1]["converged"] == "False"
