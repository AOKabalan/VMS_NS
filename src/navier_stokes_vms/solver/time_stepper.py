"""The time loop: steps of spec Time-step algorithm, logs and output."""

import logging
import time
from pathlib import Path

import firedrake as fd
import ufl
from firedrake.exceptions import ConvergenceError
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config, method_name
from navier_stokes_vms.discretization.function_spaces import (
    constant_pressure_null_space,
    mixed_space,
    subscale_space,
)
from navier_stokes_vms.discretization.quadrature import (
    boundary_measure,
    quadrature_degree,
    volume_measure,
)
from navier_stokes_vms.discretization.time_schemes import StepCoefficients, step_coefficients
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.formulation.stabilization_parameters import (
    dynamic_parameter,
    stabilization_parameters,
)
from navier_stokes_vms.formulation.vms import (
    Switches,
    accelerated_subscale,
    switches,
    velocity_subscale,
    vms_terms,
)
from navier_stokes_vms.postprocessing.diagnostics import DiagnosticState, evaluate_diagnostics
from navier_stokes_vms.postprocessing.output import CsvWriter, FieldWriter, open_diagnostics_csv
from navier_stokes_vms.postprocessing.run_logging import IterationLogs, open_iteration_logs
from navier_stokes_vms.problem.boundary_conditions import (
    Dirichlet,
    Outflow,
    dirichlet_bc,
    named_boundary_measure,
)
from navier_stokes_vms.problem.problem_definition import (
    Problem,
    SpaceTimeExpression,
    uses_pressure_null_space,
)
from navier_stokes_vms.solver.nonlinear_solver import InnerSolver, StepSolver, StepSolveResult
from navier_stokes_vms.solver.subscale_storage import SubscaleIterate, SubscaleStorage


class TimeStepper:
    """Integration of one problem: startup, histories, solves, logs and field output."""

    def __init__(
        self, config: Config, problem: Problem, folder: Path, logger: logging.Logger
    ) -> None:
        self._config = config
        self._problem = problem
        self._logger = logger
        mesh = problem.mesh
        degree = quadrature_degree(config.discretization.velocity_degree)
        self._dx = volume_measure(mesh, degree)
        space = mixed_space(
            mesh, config.discretization.velocity_degree, config.discretization.pressure_degree
        )
        velocity_space = space.sub(0).collapse()
        self._solution = fd.Function(space)
        self._velocity, self._pressure = self._solution.subfunctions
        self._pressure.rename("pressure")
        self._history = [
            fd.Function(velocity_space, name="velocity") for _ in range(self._history_levels())
        ]
        self._endpoint = fd.Function(velocity_space, name="velocity")
        self._velocity_history = fd.Function(velocity_space)
        self._frozen = fd.Function(velocity_space)
        self._alpha0 = fd.Constant(1.0)
        self._t_new = fd.Constant(0.0)
        self._t_old = fd.Constant(0.0)
        self._weight_new = fd.Constant(1.0)
        self._weight_old = fd.Constant(0.0)
        self._pressure_time = fd.Constant(0.0)
        x = fd.SpatialCoordinate(mesh)
        level = SolveLevel(
            self._alpha0,
            fd.Constant(config.time.dt),
            self._velocity_history,
            self._at_solve_level(problem.forcing, x),
        )
        conditions = problem.boundary_conditions
        outflow = [name for name, bc in conditions.items() if isinstance(bc, Outflow)]
        outflow_ds = None
        if outflow:
            outflow_ds = named_boundary_measure(
                boundary_measure(mesh, degree), outflow, problem.boundary_tags
            )
        u, p = ufl.split(self._solution)
        v, q = fd.TestFunctions(space)
        nu = problem.viscosity
        frozen = self._frozen
        jacobian = None
        self._subscales: SubscaleStorage | None = None
        self._u_sub: Expr | None = None
        self._iterate: SubscaleIterate | None = None
        stabilization = config.stabilization
        if stabilization.method == "none":
            newton = galerkin_residual(u, p, v, q, u, u, nu, level, self._dx, outflow_ds)
            picard = galerkin_residual(u, p, v, q, frozen, frozen, nu, level, self._dx, outflow_ds)
        else:
            chi = switches(stabilization)
            if chi.is_dynamic or chi.is_nonlinear:
                quadrature_space = subscale_space(mesh, degree)
            u_sub_history = None
            if chi.is_dynamic:
                self._subscales = SubscaleStorage(quadrature_space, len(self._history))
                u_sub_history = self._subscales.history
            a, a_hat = u, frozen
            if chi.is_nonlinear:
                u_sub_iterate = fd.Function(quadrature_space)
                a, a_hat = u + u_sub_iterate, frozen + u_sub_iterate
            tau_t, tau_c = self._parameters(a, nu, level, chi)
            tau_t_hat, tau_c_hat = self._parameters(a_hat, nu, level, chi)
            newton = galerkin_residual(u, p, v, q, a, u, nu, level, self._dx, outflow_ds)
            picard = galerkin_residual(u, p, v, q, a_hat, frozen, nu, level, self._dx, outflow_ds)
            if config.solver.linearization == "newton" and config.solver.newton_tau == "frozen":
                frozen_tau = newton + vms_terms(
                    u, p, v, q, a, nu, level, tau_t_hat, tau_c_hat, self._dx, chi, u_sub_history
                )
                # replace expands the derivative before restoring the current advector in τ.
                jacobian = ufl.replace(fd.derivative(frozen_tau, self._solution), {frozen: u})
            newton += vms_terms(
                u, p, v, q, a, nu, level, tau_t, tau_c, self._dx, chi, u_sub_history
            )
            picard += vms_terms(
                u, p, v, q, a_hat, nu, level, tau_t_hat, tau_c_hat, self._dx, chi, u_sub_history
            )
            u_sub = velocity_subscale(u, p, a, nu, level, tau_t, chi, u_sub_history)
            if chi.is_dynamic:
                self._u_sub = u_sub
            if chi.is_nonlinear:
                update = (
                    u_sub
                    if config.solver.subscale_update == "fixed_point"
                    else accelerated_subscale(u, p, nu, level, tau_t, chi, u_sub_history)
                )
                self._iterate = SubscaleIterate(
                    u_sub_iterate,
                    update,
                    config.solver.subscale_rtol,
                    config.solver.subscale_max_iterations,
                    is_static=not chi.is_dynamic,
                )
        bcs = [
            dirichlet_bc(
                space.sub(0), name, self._at_solve_level(bc.value, x), problem.boundary_tags
            )
            for name, bc in conditions.items()
            if isinstance(bc, Dirichlet)
        ]
        null_space = (
            constant_pressure_null_space(space) if uses_pressure_null_space(problem) else None
        )
        inner = InnerSolver(
            self._solution, newton, picard, frozen, bcs, null_space, config.solver, jacobian
        )
        self._step_solver = StepSolver(
            inner, self._solution, self._iterate, self._dx, config.solver
        )
        self._diagnostic_state = DiagnosticState(
            problem=problem,
            velocity=self._history[0],
            pressure=self._pressure,
            velocity_time=self._t_new,
            pressure_time=self._pressure_time,
            dx=self._dx,
            error_dx=fd.dx(domain=mesh, degree=2 * degree),
            pressure_null_space=null_space is not None,
        )
        self._fields = FieldWriter(folder)
        self._folder = folder
        self._logger.info(
            "%s P%d/P%d: %d unknowns, quadrature degree %d, pressure null space %s, "
            "%s Dirichlet and %s outflow boundaries",
            method_name(stabilization),
            config.discretization.velocity_degree,
            config.discretization.pressure_degree,
            space.dim(),
            degree,
            "on" if null_space is not None else "off",
            len(bcs),
            len(outflow),
        )

    @property
    def velocity(self) -> fd.Function:
        """Endpoint velocity u_h of the last completed step."""
        return self._history[0]

    @property
    def pressure(self) -> fd.Function:
        """Pressure of the last completed step, at its solve level."""
        return self._pressure

    def run(self) -> None:
        """Integrate from t = 0 to t_end, writing iteration logs, diagnostics and fields."""
        comm = self._problem.mesh.comm
        steps = round(self._config.time.t_end / self._config.time.dt)
        logs = open_iteration_logs(comm, self._folder)
        diagnostics = open_diagnostics_csv(comm, self._folder, self._config.output.diagnostics)
        try:
            self._initialise()
            for step in range(1, steps + 1):
                self._step(step, step == steps, logs, diagnostics)
        finally:
            logs.close()
            diagnostics.close()

    def _initialise(self) -> None:
        x = fd.SpatialCoordinate(self._problem.mesh)
        self._history[0].interpolate(self._problem.initial_velocity(x))
        timing = self._config.time
        if timing.startup == "exact_history":
            first = step_coefficients(timing.scheme, timing.startup, step=1)
            for level in range(1, len(first.history_weights)):
                earlier = fd.Constant(-level * timing.dt)
                self._history[level].interpolate(self._problem.exact_velocity(x, earlier))
        self._velocity.assign(self._history[0])
        self._pressure.assign(0.0)

    def _step(self, step: int, last: bool, logs: IterationLogs, diagnostics: CsvWriter) -> None:
        started = time.perf_counter()
        timing = self._config.time
        coefficients = step_coefficients(timing.scheme, timing.startup, step)
        self._set_level(step, coefficients)
        t = self._time(step)
        result = self._step_solver.solve()
        for outer in result.iterations:
            for row in outer.inner.iterations:
                logs.inner.write(
                    step=step,
                    time=t,
                    outer_iteration=outer.iteration,
                    inner_iteration=row.iteration,
                    method=row.method,
                    residual_norm=row.residual_norm,
                    relative_residual=row.relative_residual,
                    update_norm=row.update_norm,
                    linear_iterations=row.linear_iterations,
                )
            logs.outer.write(
                step=step,
                time=t,
                outer_iteration=outer.iteration,
                inner_iterations=len(outer.inner.iterations) - 1,
                velocity_change=outer.velocity_change,
                pressure_change=outer.pressure_change,
                subscale_change=outer.subscale_change,
                xi_change=None,
                zeta_change=None,
                fe_residual_norm=outer.fe_residual_norm,
                converged=outer.converged,
            )
        outer_iterations = len(result.iterations)
        inner_iterations = sum(len(outer.inner.iterations) - 1 for outer in result.iterations)
        local_iterations = sum(outer.local.iterations for outer in result.iterations if outer.local)
        if not result.converged:
            logs.time_steps.write(
                step=step,
                time=t,
                scheme_used=coefficients.scheme,
                dt=timing.dt,
                outer_iterations=outer_iterations,
                inner_iterations_total=inner_iterations,
                converged=False,
                wall_time=time.perf_counter() - started,
            )
            self._fail(step, result)
        self._advance(coefficients)
        values = evaluate_diagnostics(self._config.output.diagnostics, self._diagnostic_state)
        diagnostics.write(step=step, time=t, **values)
        logs.time_steps.write(
            step=step,
            time=t,
            scheme_used=coefficients.scheme,
            dt=timing.dt,
            outer_iterations=outer_iterations,
            inner_iterations_total=inner_iterations,
            converged=True,
            wall_time=time.perf_counter() - started,
        )
        summary = "".join(f", {name} {value:.6e}" for name, value in values.items())
        self._logger.info(
            "step %d, t = %.6g, %s: %d outer, %d inner, %d local iteration(s), residual %.3e%s",
            step,
            t,
            coefficients.scheme,
            outer_iterations,
            inner_iterations,
            local_iterations,
            result.iterations[-1].fe_residual_norm,
            summary,
        )
        vtk_every = self._config.output.vtk_every
        if vtk_every and (step % vtk_every == 0 or last):
            self._fields.write(self._history[0], self._pressure, t, float(self._pressure_time))

    def _set_level(self, step: int, coefficients: StepCoefficients) -> None:
        self._t_old.assign(self._time(step - 1))
        self._t_new.assign(self._time(step))
        self._pressure_time.assign(self._time(step - 1 + coefficients.solve_level))
        self._alpha0.assign(coefficients.alpha0)
        self._weight_new.assign(coefficients.data_weight_new)
        self._weight_old.assign(coefficients.data_weight_old)
        weighted = zip(coefficients.history_weights, self._history, strict=False)
        self._velocity_history.assign(sum(weight * velocity for weight, velocity in weighted))
        if self._subscales is not None:
            self._subscales.form_history(coefficients)
            if self._iterate is not None:
                self._iterate.value.assign(self._subscales.latest)

    def _advance(self, coefficients: StepCoefficients) -> None:
        if self._subscales is not None:
            self._subscales.store(self._u_sub, coefficients)
        self._endpoint.assign(
            coefficients.extrapolation_weight_solve * self._velocity
            + coefficients.extrapolation_weight_old * self._history[0]
        )
        for level in range(len(self._history) - 1, 0, -1):
            self._history[level].assign(self._history[level - 1])
        self._history[0].assign(self._endpoint)

    def _fail(self, step: int, result: StepSolveResult) -> None:
        solve_time = float(self._pressure_time)
        self._endpoint.assign(self._velocity)
        self._fields.write(self._endpoint, self._pressure, solve_time, solve_time)
        self._logger.error("Step %d: %s. Iteration history:", step, result.failure_reason)
        for outer in result.iterations:
            self._logger.error(
                "  outer iteration %d: velocity change %.3e, pressure change %.3e, subscale "
                "change %s, local sweeps %s, FE residual %.6e",
                outer.iteration,
                outer.velocity_change,
                outer.pressure_change,
                "-" if outer.subscale_change is None else f"{outer.subscale_change:.3e}",
                "-" if outer.local is None else outer.local.iterations,
                outer.fe_residual_norm,
            )
            for row in outer.inner.iterations:
                self._logger.error(
                    "    iteration %d (%s): residual %.6e, relative %.3e, update %s",
                    row.iteration,
                    row.method,
                    row.residual_norm,
                    row.relative_residual,
                    "-" if row.update_norm is None else f"{row.update_norm:.3e}",
                )
        raise ConvergenceError(
            f"step {step} did not converge: {result.failure_reason}; last state written to "
            f"{self._folder / 'fields'} at t = {solve_time:g}"
        )

    def _parameters(
        self, a: Expr, nu: float, level: SolveLevel, chi: Switches
    ) -> tuple[Expr, Expr]:
        stabilization = self._config.stabilization
        k = self._config.discretization.velocity_degree
        tau_m, tau_c = stabilization_parameters(a, nu, k, stabilization.c1, stabilization.c2)
        tau_t = dynamic_parameter(tau_m, level.alpha0, level.dt) if chi.is_dynamic else tau_m
        return tau_t, tau_c if stabilization.pressure_subscale else ufl.zero()

    def _time(self, steps: float) -> float:
        return float(f"{steps * self._config.time.dt:.15g}")

    def _history_levels(self) -> int:
        timing = self._config.time
        return max(
            len(step_coefficients(timing.scheme, timing.startup, step).history_weights)
            for step in (1, 2, 3)
        )

    def _at_solve_level(self, data: SpaceTimeExpression, x: Expr) -> Expr:
        return self._weight_new * data(x, self._t_new) + self._weight_old * data(x, self._t_old)
