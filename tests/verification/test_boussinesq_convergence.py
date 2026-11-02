"""B5 coupled spatial refinement and buoyancy-driven semidiscrete temporal refinement."""

import json
import logging
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl
from study import finest_rate, run_study
from ufl.core.expr import Expr

from navier_stokes_vms.config import (
    Config,
    MeshOptions,
    OutputOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.problem.thermal_problem import ThermalProblem
from navier_stokes_vms.solver.time_stepper import TimeStepper

pytestmark = pytest.mark.slow


@pytest.mark.parametrize("degree", [1, 2])
def test_coupled_manufactured_spatial_rates(tmp_path: Path, degree: int) -> None:
    base = f"""case = "manufactured_boussinesq"
[discretization]
velocity_degree = {degree}
pressure_degree = {degree}
[stabilization]
subscale_space = "oss"
subscale_dynamics = "static"
scale_splitting = "linear"
[time]
scheme = "bdf1"
dt = 0.05
t_end = 0.1
[output]
vtk_every = 0
diagnostics = [
  "velocity_error_l2", "velocity_gradient_error_l2", "pressure_error_l2",
  "temperature_error_l2", "temperature_gradient_error_l2"
]
"""
    table = run_study(tmp_path, base, f"b5_spatial_p{degree}", "--cells", "8x8", "16x16", "32x32")
    for name, order in {
        "velocity_error_l2": degree + 1,
        "velocity_gradient_error_l2": degree,
        "pressure_error_l2": degree,
        "temperature_error_l2": degree + 1,
        "temperature_gradient_error_l2": degree,
    }.items():
        assert finest_rate(table, name) >= order - 0.2, (name, table)


@pytest.mark.parametrize(
    "scheme,order", [("bdf1", 1), ("bdf2", 2), ("bdf3", 3), ("crank_nicolson", 2)]
)
def test_buoyant_discrete_shear_temporal_rates(tmp_path: Path, scheme: str, order: int) -> None:
    mesh = fd.PeriodicUnitSquareMesh(4, 4)
    S = fd.FunctionSpace(mesh, "CG", 1)
    x = fd.SpatialCoordinate(mesh)
    mode = fd.Function(S).interpolate(ufl.cos(2 * ufl.pi * x[0]))
    s = fd.TestFunction(S)
    mass = fd.assemble(mode * s * fd.dx)
    stiffness = fd.assemble(ufl.dot(ufl.grad(mode), ufl.grad(s)) * fd.dx)
    with mode.dat.vec_ro as v, mass.dat.vec_ro as m, stiffness.dat.vec_ro as k:
        eigenvalue = v.dot(k) / v.dot(m)
    nu, kappa, beta = 0.07, 0.1, 0.7

    def velocity(_x: Expr, t: Expr) -> Expr:
        amplitude = (
            beta
            * (ufl.exp(-kappa * eigenvalue * t) - ufl.exp(-nu * eigenvalue * t))
            / ((nu - kappa) * eigenvalue)
        )
        return ufl.as_vector((0.0 * mode, amplitude * mode))

    def temperature(_x: Expr, t: Expr) -> Expr:
        return ufl.exp(-kappa * eigenvalue * t) * mode

    def zero(_x: Expr, _t: Expr) -> Expr:
        return ufl.as_vector((0.0 * mode, 0.0 * mode))

    thermal = ThermalProblem(
        kappa,
        beta,
        (0.0, -1.0),
        0.0,
        lambda x, _t: 0.0 * x[0],
        lambda x: temperature(x, 0.0),
        {},
        temperature,
    )
    problem = Problem(
        mesh,
        {},
        {},
        nu,
        zero,
        lambda x: velocity(x, 0.0),
        exact_velocity=velocity,
        exact_pressure=lambda x, _t: 0.0 * x[0],
        thermal=thermal,
    )
    errors = {"velocity": [], "temperature": []}
    for dt in (0.02, 0.01, 0.005):
        folder = tmp_path / f"dt{dt}"
        folder.mkdir()
        config = Config(
            case="manufactured_boussinesq",
            run_name="b5_temporal",
            mesh=MeshOptions(cells=(4, 4)),
            time=TimeOptions(scheme=scheme, dt=dt, t_end=0.2, startup="exact_history"),
            stabilization=StabilizationOptions(scale_splitting="linear"),
            output=OutputOptions(vtk_every=0),
        )
        stepper = TimeStepper(config, problem, folder, logging.getLogger("b5_order"))
        stepper.run()
        eu = stepper.velocity - velocity(x, 0.2)
        et = stepper.temperature - temperature(x, 0.2)
        errors["velocity"].append(fd.assemble(ufl.inner(eu, eu) * stepper._dx) ** 0.5)
        errors["temperature"].append(fd.assemble(et * et * stepper._dx) ** 0.5)
    rates = {name: float(np.log2(values[-2] / values[-1])) for name, values in errors.items()}
    (tmp_path / "temporal_order.json").write_text(
        json.dumps({"scheme": scheme, "errors": errors, "rates": rates}, indent=2)
    )
    assert min(rates.values()) >= order - 0.2, rates
