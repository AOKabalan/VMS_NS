"""Stable-wave reference signs, wall compatibility and discrete hydrostatic rest."""

import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.cases.internal_gravity_wave import WaveParameters, linear_fields
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.solver.time_stepper import TimeStepper


def configuration() -> Config:
    config = Config(
        case="internal_gravity_wave",
        run_name="gravity_wave_rest",
        case_parameters={"n_waves": 1, "displacement_amplitude": 0.0},
        mesh=MeshOptions(cells=(8, 16)),
        discretization=DiscretizationOptions(2, 2),
        stabilization=StabilizationOptions(subscale_dynamics="dynamic", scale_splitting="linear"),
        time=TimeOptions(scheme="bdf2", dt=0.05, t_end=0.1),
        output=OutputOptions(vtk_every=0, diagnostics=("kinetic_energy", "passive_scalar_content")),
    )
    return replace(
        config,
        solver=replace(config.solver, petsc={**config.solver.petsc, "mat_mumps_icntl_24": 1}),
    )


def test_linear_reference_solves_buoyancy_wave_equations_and_wall_data() -> None:
    config = configuration()
    problem = build_problem(config)
    assert problem.thermal.diffusivity == problem.viscosity == WaveParameters().diffusivity
    assert problem.passive_scalars[0].diffusivity == WaveParameters().moisture_diffusivity
    x = fd.SpatialCoordinate(problem.mesh)
    t = ufl.variable(fd.Constant(0.73))
    parameters = WaveParameters(n_waves=1)
    u, p, theta = linear_fields(x, t, parameters)
    dx = fd.dx(domain=problem.mesh, degree=12)
    momentum = ufl.diff(u, t) + ufl.grad(p) - parameters.diffusivity * ufl.div(ufl.grad(u))
    momentum -= ufl.as_vector((0.0, parameters.buoyancy_coefficient * theta))
    thermal = (
        ufl.diff(theta, t)
        + parameters.temperature_gradient * u[1]
        - parameters.diffusivity * ufl.div(ufl.grad(theta))
    )
    assert fd.assemble(ufl.inner(momentum, momentum) * dx) < 1e-24
    assert fd.assemble(thermal**2 * dx) < 1e-24
    assert fd.assemble(ufl.div(u) ** 2 * dx) < 1e-24
    walls = fd.ds(domain=problem.mesh, degree=10)
    for value in (
        u[1],
        u[0].dx(1),
        theta - parameters.temperature_gradient * x[1],
        problem.passive_scalars[0].initial_value(x).dx(1),
    ):
        assert fd.assemble(value**2 * walls) < 1e-24
    q = problem.passive_scalars[0].initial_value(x)
    assert (
        fd.assemble(ufl.max_value(q - ufl.exp(-parameters.saturation_height_rate * x[1]), 0.0) * dx)
        == 0.0
    )


def test_unperturbed_stratification_remains_at_rest(tmp_path: Path) -> None:
    config = configuration()
    problem = build_problem(config)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("gravity_wave_rest"))
    stepper.run()
    dx = fd.dx(domain=problem.mesh, degree=6)
    assert fd.assemble(ufl.inner(stepper.velocity, stepper.velocity) * dx) ** 0.5 < 1e-10
