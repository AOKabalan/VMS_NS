"""The time loop: steps of spec Time-step algorithm, logs and output."""

import logging
import time
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import ufl
from firedrake.exceptions import ConvergenceError
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config, ThermalOptions, method_name
from navier_stokes_vms.discretization.function_spaces import (
    constant_pressure_null_space,
    flow_layout,
    mixed_space,
    oss_space,
    scalar_subscale_space,
    subscale_space,
)
from navier_stokes_vms.discretization.quadrature import (
    boundary_measure,
    quadrature_degree,
    volume_measure,
)
from navier_stokes_vms.discretization.time_schemes import StepCoefficients, step_coefficients
from navier_stokes_vms.formulation.boussinesq import buoyancy
from navier_stokes_vms.formulation.convective_forms import PRESSURE_UNKNOWNS, physical_pressure
from navier_stokes_vms.formulation.galerkin import SolveLevel, galerkin_residual
from navier_stokes_vms.formulation.oss_projections import projection_equations
from navier_stokes_vms.formulation.stabilization_parameters import (
    dynamic_parameter,
    stabilization_parameters,
)
from navier_stokes_vms.formulation.vms import (
    Switches,
    accelerated_subscale,
    subscale_residual,
    switches,
    velocity_subscale,
    vms_terms,
)
from navier_stokes_vms.postprocessing.diagnostics import (
    CONSERVATION_DIAGNOSTICS,
    Conservation,
    ConservationFields,
    DiagnosticState,
    advective_cfl_sampler,
    diagnostic_columns,
    evaluate_diagnostics,
)
from navier_stokes_vms.postprocessing.horizontal_profiles import HorizontalProfiles
from navier_stokes_vms.postprocessing.output import (
    CsvWriter,
    FieldWriter,
    open_diagnostics_csv,
    write_initial_state,
    write_thermal_metadata,
)
from navier_stokes_vms.postprocessing.run_logging import IterationLogs, open_iteration_logs
from navier_stokes_vms.postprocessing.thermal_budgets import (
    THERMAL_BUDGETS,
    ThermalBudgetFields,
    ThermalBudgets,
)
from navier_stokes_vms.problem.boundary_conditions import (
    DiffusiveFlux,
    Dirichlet,
    FreeSlip,
    Outflow,
    check_dirichlet_flux,
    check_free_slip_boundaries,
    dirichlet_bc,
    named_boundary_measure,
    velocity_bcs,
)
from navier_stokes_vms.problem.problem_definition import (
    Problem,
    SpaceTimeExpression,
    uses_pressure_null_space,
)
from navier_stokes_vms.problem.scalar_problem import transported_scalars
from navier_stokes_vms.problem.thermal_problem import check_thermal_supported
from navier_stokes_vms.solver.field_history import FieldHistory
from navier_stokes_vms.solver.nonlinear_solver import InnerSolver, StepSolver, StepSolveResult
from navier_stokes_vms.solver.reusable_assembly import scalar_functional
from navier_stokes_vms.solver.subscale_storage import SubscaleIterate, SubscaleStorage
from navier_stokes_vms.solver.system_builder import ThermalSystem, thermal_system


class TimeStepper:
    """Integration of one problem: startup, histories, solves, logs and field output."""

    def __init__(
        self, config: Config, problem: Problem, folder: Path, logger: logging.Logger
    ) -> None:
        check_thermal_supported(config, problem)
        if config.output.horizontal_profiles_every and problem.horizontal_profile_geometry is None:
            raise ValueError(
                "output.horizontal_profiles_every requires horizontal profile geometry"
            )
        check_free_slip_boundaries(problem.mesh, problem.boundary_conditions, problem.boundary_tags)
        self._config = config
        self._problem = problem
        self._logger = logger
        mesh = problem.mesh
        degree = quadrature_degree(config.discretization.velocity_degree)
        self._dx = volume_measure(mesh, degree)
        self._ds = boundary_measure(mesh, degree)
        stabilization = config.stabilization
        k = config.discretization.velocity_degree
        thermal = problem.thermal
        scalar_data = transported_scalars(problem)
        self._scalar_data = {scalar.name: scalar for scalar in scalar_data}
        layout = flow_layout(
            stabilization,
            thermal=thermal is not None,
            passive_names=tuple(s.name for s in problem.passive_scalars),
        )
        self._layout = layout
        if stabilization.method == "vms" and stabilization.subscale_space == "oss":
            # ξ_h, ζ_h are unknowns of the inner solve; they start at zero and then carry over.
            space = oss_space(mesh, k, stabilization.pressure_subscale, tuple(self._scalar_data))
        else:
            space = mixed_space(
                mesh,
                k,
                config.discretization.pressure_degree,
                k if thermal is not None else None,
                len(problem.passive_scalars),
            )
        velocity_space = space.sub(layout.index("velocity")).collapse()
        self._solution = fd.Function(space)
        fields = layout.bind(self._solution.subfunctions)
        self._velocity, self._pressure = fields["velocity"], fields["pressure"]
        self._pressure.rename("pressure")
        self._velocity_states = FieldHistory(velocity_space, self._history_levels(), "velocity")
        self._temperature = fields.get("temperature")
        self._scalar_fields = {name: fields[name] for name in self._scalar_data}
        self._scalar_states = {}
        for name, field in self._scalar_fields.items():
            field.rename(name)
            self._scalar_states[name] = FieldHistory(
                space.sub(layout.index(name)).collapse(), self._history_levels(), name
            )
        self._scalar_subscales: dict[str, Expr] = {}
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
            self._velocity_states.history,
            self._at_solve_level(problem.forcing, x),
        )
        conditions = problem.boundary_conditions
        outflow = [name for name, bc in conditions.items() if isinstance(bc, Outflow)]
        outflow_ds = None
        if outflow:
            outflow_ds = named_boundary_measure(self._ds, outflow, problem.boundary_tags)
        unknowns = layout.bind(ufl.split(self._solution))
        tests = layout.bind(fd.TestFunctions(space))
        u, p = unknowns["velocity"], unknowns["pressure"]
        v, q = tests["velocity"], tests["pressure"]
        nu = problem.viscosity
        frozen = self._frozen
        jacobian = None
        self._subscales: SubscaleStorage | None = None
        self._scalar_subscale_storage: dict[str, SubscaleStorage] = {}
        self._u_sub: Expr | None = None
        self._iterate: SubscaleIterate | None = None
        forms = {"convective_form": config.convection.form, "viscous_form": config.viscosity.form}
        chi = switches(stabilization) if stabilization.method == "vms" else None
        u_sub_history = None
        a, a_hat = u, frozen
        tau_m = tau_t = tau_c = tau_t_hat = tau_c_hat = None
        xi, zeta = unknowns.get("xi"), unknowns.get("zeta")
        w, r = tests.get("xi"), tests.get("zeta")
        if chi is not None:
            if chi.is_dynamic or chi.is_nonlinear:
                quadrature_space = subscale_space(mesh, degree)
            if chi.is_dynamic:
                self._subscales = SubscaleStorage(
                    quadrature_space, len(self._velocity_states.levels)
                )
                u_sub_history = self._subscales.history
                if thermal is not None:
                    scalar_quadrature = scalar_subscale_space(mesh, degree)
                    self._scalar_subscale_storage = {
                        name: SubscaleStorage(scalar_quadrature, self._history_levels())
                        for name in self._scalar_data
                    }
            if chi.is_nonlinear:
                u_sub_iterate = fd.Function(quadrature_space)
                a, a_hat = u + u_sub_iterate, frozen + u_sub_iterate
            tau_m, tau_t, tau_c = self._parameters(a, nu, level, chi)
            _, tau_t_hat, tau_c_hat = self._parameters(a_hat, nu, level, chi)

        def coupled(a: Expr, u_h: Expr, tau_t: Expr | None, tau_c: Expr | None) -> ufl.Form:
            form = galerkin_residual(u, p, v, q, a, u_h, nu, level, self._dx, outflow_ds, **forms)
            if chi is None:
                return form
            form += vms_terms(
                u,
                p,
                v,
                q,
                a,
                nu,
                level,
                tau_t,
                tau_c,
                self._dx,
                chi,
                u_sub_history,
                xi,
                zeta,
                **forms,
            )
            if not chi.is_asgs:
                residual = subscale_residual(u, p, a, nu, level, chi, u_sub_history, **forms)
                form += projection_equations(
                    xi, zeta, w, r, residual, ufl.div(u), tau_t, tau_c, self._dx
                )
            return form

        if thermal is not None:
            scalar_levels = {
                name: SolveLevel(
                    self._alpha0,
                    level.dt,
                    self._scalar_states[name].history,
                    self._at_solve_level(data.source, x),
                )
                for name, data in self._scalar_data.items()
            }

            def thermal_at(advector: Expr, parameter_advector: Expr | None = None) -> ThermalSystem:
                return thermal_system(
                    self._solution,
                    layout,
                    problem,
                    level,
                    scalar_levels,
                    self._dx,
                    config,
                    parameter_advector,
                    velocity_history=u_sub_history,
                    scalar_histories={
                        name: storage.history
                        for name, storage in self._scalar_subscale_storage.items()
                    },
                    advector=advector,
                )

            system = thermal_at(a)
            flux = 0
            self._u_sub = system.velocity_subscale
            self._scalar_subscales = system.scalar_subscales
            for name, data in self._scalar_data.items():
                for boundary, bc in data.boundary_conditions.items():
                    if isinstance(bc, DiffusiveFlux) and bc.value is not None:
                        flux += (
                            self._at_solve_level(bc.value, x)
                            * tests[name]
                            * named_boundary_measure(self._ds, (boundary,), problem.boundary_tags)
                        )
            newton = system.residual + flux
            picard = thermal_at(a_hat).residual + flux
            if chi is not None and config.solver.newton_tau == "frozen":
                frozen_system = thermal_at(a, a_hat)
                jacobian = ufl.replace(
                    fd.derivative(frozen_system.residual, self._solution), {frozen: u}
                )
        else:
            newton = coupled(a, u, tau_t, tau_c)
            picard = coupled(a_hat, frozen, tau_t_hat, tau_c_hat)

        if chi is not None and thermal is None:
            if config.solver.linearization == "newton" and config.solver.newton_tau == "frozen":
                # τ on the frozen advector, then restored after differentiating: J has no τ′,
                # in the projection weights as elsewhere.
                frozen_tau = coupled(a, u, tau_t_hat, tau_c_hat)
                jacobian = ufl.replace(fd.derivative(frozen_tau, self._solution), {frozen: u})
            self._u_sub = velocity_subscale(
                u, p, a, nu, level, tau_t, chi, u_sub_history, xi, **forms
            )
        if chi is not None and chi.is_nonlinear:
            local_level = level
            if thermal is not None:
                local_level = replace(
                    level,
                    forcing=level.forcing
                    + buoyancy(
                        unknowns["temperature"],
                        thermal.expansion_coefficient,
                        thermal.gravity,
                        thermal.reference_temperature,
                    ),
                )
            update = (
                self._u_sub
                if config.solver.subscale_update == "fixed_point"
                else accelerated_subscale(
                    u, p, nu, local_level, tau_t, chi, u_sub_history, xi, **forms
                )
            )
            self._iterate = SubscaleIterate(
                u_sub_iterate,
                update,
                config.solver.subscale_rtol,
                config.solver.subscale_max_iterations,
                is_static=not chi.is_dynamic,
            )
        velocity_values = {
            name: self._at_solve_level(bc.value, x)
            for name, bc in conditions.items()
            if isinstance(bc, Dirichlet)
        }
        bcs = velocity_bcs(
            space.sub(layout.index("velocity")), conditions, velocity_values, problem.boundary_tags
        )
        slip_conditions = {name: bc for name, bc in conditions.items() if isinstance(bc, FreeSlip)}
        self._initial_slip_bcs = velocity_bcs(
            velocity_space, slip_conditions, {}, problem.boundary_tags
        )
        for name, data in self._scalar_data.items():
            bcs.extend(
                dirichlet_bc(
                    space.sub(layout.index(name)),
                    boundary,
                    self._at_solve_level(bc.value, x),
                    problem.boundary_tags,
                )
                for boundary, bc in data.boundary_conditions.items()
                if isinstance(bc, Dirichlet)
            )
        null_space = (
            constant_pressure_null_space(space) if uses_pressure_null_space(problem) else None
        )
        # With the null space the pressure is fixed by its zero spatial mean (spec: Boundary
        # conditions and quadrature), not by the solver's zero coefficient mean.
        self._zero_mean_pressure = null_space is not None
        self._boundary_velocity = None
        self._boundary_bcs = []
        if null_space is not None and any(
            isinstance(bc, Dirichlet | FreeSlip) for bc in conditions.values()
        ):
            self._boundary_velocity = fd.Function(velocity_space)
            self._boundary_bcs = velocity_bcs(
                velocity_space, conditions, velocity_values, problem.boundary_tags
            )
        self._area = fd.assemble(fd.Constant(1.0) * self._dx)
        # EMAC variants keep the unknown and its null space on P; every reported pressure is the
        # physical p_h at the solve-level velocity, with its own zero mean (spec: EMAC pressure).
        self._physical_pressure = physical_pressure(
            self._pressure, self._velocity, config.convection.form
        )
        self._pressure_offset = fd.Constant(0.0)
        self._reported_pressure = self._pressure
        reported_pressure: Expr = self._pressure
        if PRESSURE_UNKNOWNS[config.convection.form] == "P":
            # Degree max(k_p, 2k_u) contains P_h and ½|u_h|² on affine cells, so nothing is lost.
            output_degree = max(config.discretization.pressure_degree, 2 * k)
            self._reported_pressure = fd.Function(
                fd.FunctionSpace(mesh, "CG", output_degree), name="pressure"
            )
            reported_pressure = self._physical_pressure - self._pressure_offset
        thermal_options = config.thermal or ThermalOptions()
        residual_scales = None
        if thermal is not None:
            residual_scales = {}
            for name in layout.names:
                if name in ("velocity", "xi"):
                    scale = thermal_options.momentum_residual_scale
                elif name in ("pressure", "zeta"):
                    scale = thermal_options.continuity_residual_scale
                else:
                    scale = thermal_options.temperature_residual_scale
                residual_scales[name] = scale
        inner = InnerSolver(
            self._solution,
            newton,
            picard,
            frozen,
            bcs,
            null_space,
            config.solver,
            jacobian,
            residual_scales,
        )
        self._step_solver = StepSolver(
            inner,
            self._solution,
            self._iterate,
            self._dx,
            config.solver,
            layout,
            scalar_subscales=self._scalar_subscales if self._iterate is not None else None,
            scalar_space=scalar_subscale_space(mesh, degree)
            if thermal is not None and self._iterate is not None
            else None,
            velocity_subscale=self._u_sub
            if thermal is not None and self._iterate is not None
            else None,
        )
        self._conservation = None
        requested = [name for name in config.output.diagnostics if name in CONSERVATION_DIAGNOSTICS]
        if requested:
            self._conservation = Conservation(
                requested,
                ConservationFields(
                    velocity=u,
                    advector=a,
                    subscale=self._u_sub,
                    endpoint_velocity=self._velocity_states.levels[0],
                    endpoint_subscale=None if self._subscales is None else self._subscales.latest,
                    zeta=zeta,
                    tau_m=tau_m,
                    tau_c=tau_c,
                    forcing=level.forcing,
                    viscosity=nu,
                    dt=config.time.dt,
                    dx=self._dx,
                    outflow_ds=outflow_ds,
                    convective_form=config.convection.form,
                    viscous_form=config.viscosity.form,
                    is_nonlinear=chi is not None and chi.is_nonlinear,
                    is_asgs=chi is None or chi.is_asgs,
                ),
            )
        subscale_velocity = self._u_sub
        subscale_time = self._pressure_time if self._u_sub is not None else None
        if self._subscales is not None:
            subscale_velocity = self._subscales.latest
            subscale_time = self._t_new
        elif self._iterate is not None:
            subscale_velocity = self._iterate.value
        self._thermal_budgets = None
        thermal_requested = [name for name in config.output.diagnostics if name in THERMAL_BUDGETS]
        if thermal is not None and thermal_requested:
            self._thermal_budgets = ThermalBudgets(
                thermal_requested,
                ThermalBudgetFields(
                    config=config,
                    problem=problem,
                    velocity=self._velocity,
                    advector=a,
                    velocity_subscale=self._u_sub,
                    velocity_history=u_sub_history,
                    scalar_fields=self._scalar_fields,
                    scalar_subscales=self._scalar_subscales,
                    scalar_histories={
                        name: storage.history
                        for name, storage in self._scalar_subscale_storage.items()
                    },
                    scalar_projections={
                        name: unknowns[f"eta_{name}"]
                        for name in self._scalar_data
                        if f"eta_{name}" in unknowns
                    },
                    scalar_fluxes={
                        name: {
                            boundary: self._at_solve_level(bc.value, x)
                            for boundary, bc in data.boundary_conditions.items()
                            if isinstance(bc, DiffusiveFlux) and bc.value is not None
                        }
                        for name, data in self._scalar_data.items()
                    },
                    momentum_level=level,
                    scalar_levels=scalar_levels,
                    xi=xi,
                    zeta=zeta,
                    tau_u=tau_m,
                    tau_c=tau_c,
                    dx=self._dx,
                    ds=self._ds,
                ),
            )
        self._diagnostic_state = DiagnosticState(
            problem=problem,
            velocity=self._velocity_states.levels[0],
            pressure=reported_pressure,
            velocity_time=self._t_new,
            pressure_time=self._pressure_time,
            dx=self._dx,
            error_dx=fd.dx(domain=mesh, degree=2 * degree),
            pressure_null_space=null_space is not None,
            solve_velocity=self._velocity,
            ds=self._ds,
            subscale_velocity=subscale_velocity,
            subscale_time=subscale_time,
            conservation=self._conservation,
            temperature=self.temperature,
            solve_temperature=self._temperature,
            passive_scalars=self.passive_scalars,
            thermal_budgets=self._thermal_budgets,
            advective_cfl=advective_cfl_sampler(
                self.velocity, config.time.dt, config.discretization.velocity_degree, degree
            )
            if "advective_cfl" in config.output.diagnostics or problem.initial_pressure is not None
            else None,
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
        return self._velocity_states.levels[0]

    @property
    def pressure(self) -> fd.Function:
        """Reported pressure of the last completed step at its solve level, physical for EMAC."""
        return self._reported_pressure

    @property
    def temperature(self) -> fd.Function | None:
        """Endpoint temperature of the last successful step, when thermal physics is active."""
        return None if "temperature" not in self._scalar_states else self.scalar("temperature")

    def scalar(self, name: str) -> fd.Function:
        """Endpoint value of a named transported scalar after the last successful step."""
        return self._scalar_states[name].levels[0]

    @property
    def passive_scalars(self) -> dict[str, fd.Function]:
        """Named passive scalar endpoints for output and diagnostics."""
        return {scalar.name: self.scalar(scalar.name) for scalar in self._problem.passive_scalars}

    def run(self, after_step: Callable[[int, float], bool] | None = None) -> None:
        """Integrate from t = 0 to t_end, writing iteration logs, diagnostics and fields."""
        comm = self._problem.mesh.comm
        steps = round(self._config.time.t_end / self._config.time.dt)
        logs = open_iteration_logs(
            comm,
            self._folder,
            tuple(
                name
                for name in self._layout.names
                if name not in ("velocity", "pressure", "xi", "zeta")
            )
            if self._temperature is not None
            else (),
            self._layout.names if self._temperature is not None else (),
            tuple(self._scalar_data) if self._iterate is not None else (),
        )
        columns = diagnostic_columns(self._config.output.diagnostics, self._problem)
        diagnostics = open_diagnostics_csv(comm, self._folder, columns)
        profiles = None
        try:
            self._initialise()
            write_thermal_metadata(comm, self._folder, self._problem, self._config)
            write_initial_state(comm, self._folder, self._initial_state())
            profile_every = self._config.output.horizontal_profiles_every
            if profile_every:
                profiles = HorizontalProfiles(
                    self._folder,
                    self._problem.horizontal_profile_geometry,
                    self.velocity,
                    {name: self.scalar(name) for name in self._scalar_data},
                    self._config.output.horizontal_fourier_modes,
                )
                profiles.write(0, 0.0)
            for step in range(1, steps + 1):
                self._step(step, step == steps, logs, diagnostics)
                stopped = after_step is not None and after_step(step, self._time(step))
                if profiles is not None and (step % profile_every == 0 or step == steps or stopped):
                    profiles.write(step, self._time(step))
                if stopped:
                    if self._config.output.vtk_every and step != steps:
                        self._fields.write(
                            self.velocity,
                            self.pressure,
                            self._time(step),
                            float(self._pressure_time),
                            self.temperature,
                            self.passive_scalars,
                        )
                    break
        finally:
            if profiles is not None:
                profiles.close()
            logs.close()
            diagnostics.close()

    def _initialise(self) -> None:
        x = fd.SpatialCoordinate(self._problem.mesh)
        self._velocity_states.levels[0].interpolate(self._problem.initial_velocity(x))
        timing = self._config.time
        if timing.startup == "exact_history":
            first = step_coefficients(timing.scheme, timing.startup, step=1)
            for level in range(1, len(first.history_weights)):
                earlier = fd.Constant(-level * timing.dt)
                self._velocity_states.levels[level].interpolate(
                    self._problem.exact_velocity(x, earlier)
                )
        for state in self._velocity_states.levels:
            for bc in self._initial_slip_bcs:
                bc.apply(state)
        self._velocity.assign(self._velocity_states.levels[0])
        self._pressure.assign(0.0)
        for name, data in self._scalar_data.items():
            self._scalar_states[name].levels[0].interpolate(data.initial_value(x))
            self._scalar_fields[name].assign(self.scalar(name))
            if timing.startup == "exact_history":
                first = step_coefficients(timing.scheme, timing.startup, step=1)
                for history_level in range(1, len(first.history_weights)):
                    earlier = fd.Constant(-history_level * timing.dt)
                    self._scalar_states[name].levels[history_level].interpolate(
                        data.exact_solution(x, earlier)
                    )
        if self._problem.initial_pressure is not None:
            initial_pressure = self._problem.initial_pressure(x)
            if PRESSURE_UNKNOWNS[self._config.convection.form] == "P":
                initial_pressure -= ufl.inner(self._velocity, self._velocity) / 2
            self._pressure.interpolate(initial_pressure)
            self._normalize_pressure()

    def _initial_state(self) -> dict[str, float]:
        """|Ω|, d, and the energy and moments of u_h(0), which fix the conservation-run scales.

        Spec: Conservation-run normalization (owner decision K-D2).
        """
        names = ("kinetic_energy", "linear_momentum", "angular_momentum")
        values = {
            "time": 0.0,
            "domain_volume": self._area,
            "dimension": self._problem.mesh.geometric_dimension,
            **evaluate_diagnostics(names, self._diagnostic_state),
        }
        if self._problem.initial_pressure is not None:
            normal = fd.FacetNormal(self._problem.mesh)
            flux = ufl.dot(self.velocity, normal)
            values.update(
                evaluate_diagnostics(
                    ("divergence_l2", "vertical_kinetic_energy", "advective_cfl"),
                    self._diagnostic_state,
                )
            )
            values["pressure_mean"] = scalar_functional(self.pressure * self._dx)() / self._area
            values["boundary_flux"] = scalar_functional(flux * self._ds)()
            values["boundary_normal_l2"] = scalar_functional(flux**2 * self._ds)() ** 0.5
            if self.temperature is not None:
                thermal = self._problem.thermal
                mismatch = ufl.grad(self._physical_pressure) - buoyancy(
                    self._temperature,
                    thermal.expansion_coefficient,
                    thermal.gravity,
                    thermal.reference_temperature,
                )
                values["pressure_buoyancy_defect_l2"] = (
                    scalar_functional(ufl.inner(mismatch, mismatch) * self._dx)() ** 0.5
                )
                for name, data in self._scalar_data.items():
                    flux = data.diffusivity * ufl.dot(ufl.grad(self.scalar(name)), normal)
                    values[f"{name}_boundary_flux_l2"] = (
                        scalar_functional(flux**2 * self._ds)() ** 0.5
                    )
        return values

    def _step(self, step: int, last: bool, logs: IterationLogs, diagnostics: CsvWriter) -> None:
        started = time.perf_counter()
        timing = self._config.time
        coefficients = step_coefficients(timing.scheme, timing.startup, step)
        self._set_level(step, coefficients)
        t = self._time(step)
        result = self._step_solver.solve()
        self._normalize_pressure()
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
                    **(
                        {}
                        if row.field_residuals is None
                        else {
                            f"{name}_residual": value for name, value in row.field_residuals.items()
                        }
                    ),
                )
            logs.outer.write(
                step=step,
                time=t,
                outer_iteration=outer.iteration,
                inner_iterations=len(outer.inner.iterations) - 1,
                **(
                    {
                        f"{name}_change": change
                        for name, change in outer.field_changes.items()
                        if name not in ("velocity", "pressure", "xi", "zeta")
                    }
                    if self._temperature is not None
                    else {}
                ),
                velocity_change=outer.velocity_change,
                pressure_change=outer.pressure_change,
                subscale_change=outer.subscale_change,
                xi_change=outer.xi_change,
                zeta_change=outer.zeta_change,
                fe_residual_norm=outer.fe_residual_norm,
                converged=outer.converged,
                **(
                    {
                        "velocity_subscale_residual": outer.velocity_closure_residual,
                        **{
                            f"{name}_subscale_change": change
                            for name, change in outer.scalar_subscale_changes.items()
                        },
                        **{
                            f"{name}_residual": value
                            for name, value in outer.field_residuals.items()
                        },
                    }
                    if outer.scalar_subscale_changes is not None
                    else {}
                ),
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
        if self._conservation is not None:
            self._conservation.capture(coefficients.scheme == "crank_nicolson")
        if self._thermal_budgets is not None:
            self._thermal_budgets.capture()
        self._advance(coefficients)
        values = evaluate_diagnostics(self._config.output.diagnostics, self._diagnostic_state)
        subscale_time = self._diagnostic_state.subscale_time
        diagnostics.write(
            step=step,
            time=t,
            pressure_time=float(self._pressure_time),
            subscale_time=None if subscale_time is None else float(subscale_time),
            **values,
        )
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
        summary = "".join(
            f", {name} {value:.6e}" for name, value in values.items() if value is not None
        )
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
            self._fields.write(
                self.velocity,
                self.pressure,
                t,
                float(self._pressure_time),
                self.temperature,
                self.passive_scalars,
            )

    def _normalize_pressure(self) -> None:
        if self._zero_mean_pressure:
            mean = fd.assemble(self._pressure * self._dx) / self._area
            self._pressure.assign(self._pressure - mean)
        if self._reported_pressure is not self._pressure:
            if self._zero_mean_pressure:
                physical_mean = fd.assemble(self._physical_pressure * self._dx) / self._area
                self._pressure_offset.assign(physical_mean)
            self._reported_pressure.interpolate(self._physical_pressure - self._pressure_offset)

    def _set_level(self, step: int, coefficients: StepCoefficients) -> None:
        self._t_old.assign(self._time(step - 1))
        self._t_new.assign(self._time(step))
        self._pressure_time.assign(self._time(step - 1 + coefficients.solve_level))
        self._alpha0.assign(coefficients.alpha0)
        self._weight_new.assign(coefficients.data_weight_new)
        self._weight_old.assign(coefficients.data_weight_old)
        if self._boundary_velocity is not None:
            for bc in self._boundary_bcs:
                bc.apply(self._boundary_velocity)
            try:
                check_dirichlet_flux(self._boundary_velocity, self._ds)
            except ValueError as error:
                raise ValueError(
                    f"step {step}, solve time {float(self._pressure_time):g}: {error}"
                ) from error
        if self._problem.thermal is not None and self._boundary_velocity is not None:
            normal = fd.FacetNormal(self._problem.mesh)
            flux_square = fd.assemble(ufl.dot(self._boundary_velocity, normal) ** 2 * self._ds)
            if flux_square > 1e-20:
                raise ValueError("B1 thermal boundaries require zero normal velocity")
        self._velocity_states.form_history(coefficients)
        for state in self._scalar_states.values():
            state.form_history(coefficients)
        for storage in self._scalar_subscale_storage.values():
            storage.form_history(coefficients)
        if self._subscales is not None:
            self._subscales.form_history(coefficients)
            if self._iterate is not None:
                self._iterate.value.assign(self._subscales.latest)

    def _advance(self, coefficients: StepCoefficients) -> None:
        if self._subscales is not None:
            self._subscales.store(self._u_sub, coefficients)
        for name, storage in self._scalar_subscale_storage.items():
            storage.store(self._scalar_subscales[name], coefficients)
        self._velocity_states.advance(self._velocity, coefficients)
        for name, state in self._scalar_states.items():
            state.advance(self._scalar_fields[name], coefficients)

    def _fail(self, step: int, result: StepSolveResult) -> None:
        solve_time = float(self._pressure_time)
        self._velocity_states.endpoint.assign(self._velocity)
        self._fields.write(
            self._velocity_states.endpoint,
            self.pressure,
            solve_time,
            solve_time,
            self._temperature,
            {
                scalar.name: self._scalar_fields[scalar.name]
                for scalar in self._problem.passive_scalars
            },
        )
        self._logger.error("Step %d: %s. Iteration history:", step, result.failure_reason)
        for outer in result.iterations:
            self._logger.error(
                "  outer iteration %d: velocity change %.3e, pressure change %.3e, subscale "
                "change %s, local sweeps %s, ξ change %s, ζ change %s, FE residual %.6e",
                outer.iteration,
                outer.velocity_change,
                outer.pressure_change,
                _optional(outer.subscale_change),
                "-" if outer.local is None else outer.local.iterations,
                _optional(outer.xi_change),
                _optional(outer.zeta_change),
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
    ) -> tuple[Expr, Expr, Expr]:
        stabilization = self._config.stabilization
        k = self._config.discretization.velocity_degree
        tau_m, tau_c = stabilization_parameters(a, nu, k, stabilization.c1, stabilization.c2)
        tau_t = dynamic_parameter(tau_m, level.alpha0, level.dt) if chi.is_dynamic else tau_m
        return tau_m, tau_t, tau_c if stabilization.pressure_subscale else ufl.zero()

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


def _optional(change: float | None) -> str:
    return "-" if change is None else f"{change:.3e}"
