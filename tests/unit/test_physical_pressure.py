"""Physical pressure p_h = P_h + ½|u_h|² reported for the EMAC variants (spec: EMAC pressure)."""

import csv
import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    ConvectionOptions,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    StabilizationOptions,
    TimeOptions,
    ViscosityOptions,
)
from navier_stokes_vms.solver.time_stepper import TimeStepper

LOGGER = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def config_for(
    case: str,
    mesh: MeshOptions,
    time: TimeOptions,
    form: str,
    galerkin: bool,
    diagnostics: tuple[str, ...],
) -> Config:
    k = 2 if galerkin else 1
    return Config(
        case=case,
        run_name="physical_pressure",
        mesh=mesh,
        discretization=DiscretizationOptions(velocity_degree=k, pressure_degree=1),
        time=time,
        stabilization=StabilizationOptions(method="none") if galerkin else StabilizationOptions(),
        convection=ConvectionOptions(form=form),
        viscosity=ViscosityOptions(form="stress"),
        output=OutputOptions(vtk_every=1, diagnostics=diagnostics),
    )


def norm(expression: Expr, mesh: fd.MeshGeometry) -> float:
    return np.sqrt(fd.assemble(expression**2 * fd.dx(domain=mesh, degree=12)))


def last_row(folder: Path) -> dict[str, str]:
    with (folder / "diagnostics.csv").open() as stream:
        return list(csv.DictReader(stream))[-1]


@pytest.mark.parametrize("scheme", ["bdf2", "crank_nicolson"])
@pytest.mark.parametrize(
    ("form", "galerkin"), [("emac", False), ("emac_vms", False), ("emac", True)]
)
def test_physical_pressure_uses_the_solve_level_velocity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scheme: str, form: str, galerkin: bool
) -> None:
    """With the null space, P_h and p_h each have zero mean; output, error and force use p_h."""
    config = config_for(
        "manufactured_temporal",
        MeshOptions(cells=(2, 2)),
        TimeOptions(scheme=scheme, dt=0.1, t_end=0.2),
        form,
        galerkin,
        ("pressure_error_l2", "force"),
    )
    problem = replace(build_problem(config), force_boundaries=("y_max",))
    stepper = TimeStepper(config, problem, tmp_path, LOGGER)
    written = []
    write = stepper._fields.write

    def record(velocity: fd.Function, pressure: fd.Function, *times: float) -> None:
        written.append((pressure, pressure.copy(deepcopy=True)))
        write(velocity, pressure, *times)

    monkeypatch.setattr(stepper._fields, "write", record)
    endpoints = []
    stepper.run(lambda _step, _time: endpoints.append(stepper.velocity.copy(deepcopy=True)))
    mesh = problem.mesh
    old, new = endpoints
    level = (old + new) / 2 if scheme == "crank_nicolson" else new
    pressure_time = 0.15 if scheme == "crank_nicolson" else 0.2
    assert norm(stepper._velocity - level, mesh) <= 1e-13 * norm(level, mesh)
    unknown = stepper._pressure
    area = fd.assemble(fd.Constant(1.0) * fd.dx(domain=mesh))
    assert abs(fd.assemble(unknown * fd.dx(domain=mesh, degree=12))) <= 1e-13 * area
    raw = unknown + ufl.inner(level, level) / 2
    expected = raw - fd.assemble(raw * fd.dx(domain=mesh, degree=12)) / area
    assert norm(expected - unknown, mesh) > 1e-6

    reported = stepper.pressure
    element = reported.function_space().ufl_element()
    assert element.family() == "Lagrange"
    assert element.degree() == 2 * config.discretization.velocity_degree
    assert norm(reported - expected, mesh) <= 1e-12 * norm(expected, mesh)
    assert len(written) == 2
    assert all(function is reported for function, _ in written)
    assert norm(written[-1][1] - expected, mesh) <= 1e-12 * norm(expected, mesh)
    if scheme == "crank_nicolson":
        averaged = unknown + (ufl.inner(old, old) + ufl.inner(new, new)) / 4
        averaged -= fd.assemble(averaged * fd.dx(domain=mesh, degree=12)) / area
        assert norm(averaged - expected, mesh) > 1e6 * 1e-12 * norm(expected, mesh)

    row = last_row(tmp_path)
    assert float(row["pressure_time"]) == pytest.approx(pressure_time)
    x = fd.SpatialCoordinate(mesh)
    error = expected - problem.exact_pressure(x, fd.Constant(pressure_time))
    error -= fd.assemble(error * fd.dx(domain=mesh, degree=12)) / area
    assert float(row["pressure_error_l2"]) == pytest.approx(norm(error, mesh), rel=1e-10)
    n = fd.FacetNormal(mesh)
    gradient = ufl.grad(level)
    stress = -expected * ufl.Identity(2) + problem.viscosity * (gradient + ufl.transpose(gradient))
    ds = fd.ds(domain=mesh, subdomain_id=problem.boundary_tags["y_max"], degree=12)
    for i, axis in enumerate("xy"):
        force = fd.assemble(-ufl.dot(stress, n)[i] * ds)
        assert float(row[f"force_y_max_{axis}"]) == pytest.approx(force, rel=1e-10, abs=1e-13)


def test_outflow_force_uses_the_unnormalized_physical_pressure(tmp_path: Path) -> None:
    """Without the null space only P_h's gauge is fixed by the outflow, so p_h is not shifted."""
    config = config_for(
        "cylinder_dfg_2d",
        MeshOptions(file=str(ROOT / "meshes/cylinder_dfg_2d_coarse.msh")),
        TimeOptions(scheme="crank_nicolson", dt=0.005, t_end=0.01),
        "emac_vms",
        False,
        ("force",),
    )
    problem = build_problem(config)
    stepper = TimeStepper(config, problem, tmp_path, LOGGER)
    stepper.run()
    assert float(stepper._pressure_offset) == 0.0
    u = stepper._velocity
    physical = stepper._pressure + ufl.inner(u, u) / 2
    assert norm(stepper.pressure - physical, problem.mesh) <= 1e-12 * norm(physical, problem.mesh)
    n = fd.FacetNormal(problem.mesh)
    gradient = ufl.grad(u)
    stress = -physical * ufl.Identity(2) + problem.viscosity * (gradient + ufl.transpose(gradient))
    ds = fd.ds(domain=problem.mesh, subdomain_id=problem.boundary_tags["cylinder"], degree=12)
    row = last_row(tmp_path)
    for i, axis in enumerate("xy"):
        force = fd.assemble(-ufl.dot(stress, n)[i] * ds)
        assert float(row[f"force_cylinder_{axis}"]) == pytest.approx(force, rel=1e-10)


@pytest.mark.parametrize(
    "form", ["nonconservative", "conservative", "skew_symmetric_1", "skew_symmetric_2"]
)
def test_the_other_forms_report_the_pressure_unknown(tmp_path: Path, form: str) -> None:
    config = config_for(
        "manufactured_temporal",
        MeshOptions(cells=(2, 2)),
        TimeOptions(scheme="crank_nicolson", dt=0.1, t_end=0.2),
        form,
        False,
        ("pressure_error_l2",),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    assert stepper.pressure is stepper._pressure
    assert stepper._diagnostic_state.pressure is stepper._pressure
    assert stepper.pressure.function_space().ufl_element().degree() == 1
