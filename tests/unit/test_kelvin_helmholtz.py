"""Physical setup and component wall constraints of the KH experiment."""

import logging
import math
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    ConvectionOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    TimeOptions,
    ViscosityOptions,
)
from navier_stokes_vms.problem.boundary_conditions import free_slip_bc
from navier_stokes_vms.problem.problem_definition import uses_pressure_null_space
from navier_stokes_vms.solver.time_stepper import TimeStepper


def test_paper_initial_condition_and_channel_geometry() -> None:
    cfg = Config(case="kelvin_helmholtz", run_name="test", mesh=MeshOptions(cells=(12, 16)))
    p = build_problem(cfg)
    assert p.viscosity == pytest.approx(1 / 28000)
    assert uses_pressure_null_space(p)
    assert p.exact_velocity is None
    x = fd.SpatialCoordinate(p.mesh)
    ds, dx = fd.ds(domain=p.mesh, degree=4), fd.dx(domain=p.mesh, degree=12)
    assert fd.assemble(fd.Constant(1) * ds) == pytest.approx(2)
    assert fd.assemble(x[1] * fd.ds(3, domain=p.mesh)) == pytest.approx(0)
    assert fd.assemble(x[1] * fd.ds(4, domain=p.mesh)) == pytest.approx(1)
    s, delta = x[1] - 0.5, 1 / 28
    g = ufl.exp(-((s / delta) ** 2))
    by_hand = ufl.as_vector(
        (
            ufl.tanh(2 * s / delta)
            - 0.001
            * 2
            * s
            / delta**2
            * g
            * (ufl.cos(8 * math.pi * x[0]) + ufl.cos(20 * math.pi * x[0])),
            0.001
            * g
            * (
                8 * math.pi * ufl.sin(8 * math.pi * x[0])
                + 20 * math.pi * ufl.sin(20 * math.pi * x[0])
            ),
        )
    )
    initial = p.initial_velocity(x)
    assert fd.assemble(ufl.inner(initial - by_hand, initial - by_hand) * dx) < 1e-25
    assert fd.assemble(ufl.div(initial) ** 2 * dx) < 1e-23


def test_free_slip_changes_normal_component_only_in_mixed_space() -> None:
    p = build_problem(
        Config(case="kelvin_helmholtz", run_name="test", mesh=MeshOptions(cells=(4, 4)))
    )
    V = fd.VectorFunctionSpace(p.mesh, "CG", 2)
    Q = fd.FunctionSpace(p.mesh, "CG", 1)
    state = fd.Function(V * Q)
    state.subfunctions[0].assign(fd.Constant((2.0, 3.0)))
    for name in p.boundary_tags:
        free_slip_bc(state.function_space().sub(0), name, 1, p.boundary_tags).apply(state)
    u = state.subfunctions[0]
    assert fd.assemble(u[1] ** 2 * fd.ds(domain=p.mesh)) < 1e-24
    assert fd.assemble((u[0] - 2) ** 2 * fd.dx(domain=p.mesh)) < 1e-24
    with pytest.raises(ValueError, match="not normal to axis"):
        free_slip_bc(V, "y_min", 0, p.boundary_tags)


@pytest.mark.parametrize("form", ["skew_symmetric_1", "emac_vms"])
def test_free_slip_uniform_tangential_flow_is_preserved(form: str, tmp_path: Path) -> None:
    cfg = Config(
        case="kelvin_helmholtz",
        run_name="test",
        mesh=MeshOptions(cells=(4, 4)),
        time=TimeOptions(scheme="bdf2", dt=0.01, t_end=0.02),
        convection=ConvectionOptions(form=form),
        viscosity=ViscosityOptions(form="stress"),
        solver=SolverOptions(inner_max_iterations=200),
        output=OutputOptions(vtk_every=0, diagnostics=("kinetic_energy",)),
    )
    p = build_problem(cfg)
    exact = fd.Constant((1.0, 0.0))
    p = replace(p, initial_velocity=lambda _x: exact)
    stepper = TimeStepper(cfg, p, tmp_path, logging.getLogger("kh_test"))
    stepper.run()
    error = stepper.velocity - exact
    assert fd.assemble(ufl.inner(error, error) * fd.dx(domain=p.mesh)) < 1e-18
    assert fd.assemble(stepper.pressure**2 * fd.dx(domain=p.mesh)) < 1e-18
