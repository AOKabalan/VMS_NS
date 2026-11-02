"""Solved Boussinesq budget identities across forms, scales and time histories."""

import csv
import logging
from dataclasses import replace
from pathlib import Path

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
    ThermalOptions,
    TimeOptions,
)
from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, Dirichlet
from navier_stokes_vms.problem.mesh_loading import box, periodic_rectangle
from navier_stokes_vms.problem.problem_definition import (
    Problem,
    at_initial_time,
    exact_solution_forcing,
)
from navier_stokes_vms.problem.scalar_problem import PassiveScalarProblem
from navier_stokes_vms.problem.thermal_problem import ThermalProblem
from navier_stokes_vms.solver.time_stepper import TimeStepper

BUDGETS = (
    "thermal_content_budget",
    "thermal_variance_budget",
    "boussinesq_energy_budget",
    "thermal_subscale_l2",
    "temperature_extrema",
)


def configuration(
    mode: str = "oss",
    dynamics: str = "dynamic",
    splitting: str = "nonlinear",
    form: str = "skew_symmetric",
    scheme: str = "bdf2",
    degree: int = 1,
) -> Config:
    return Config(
        case="manufactured_boussinesq",
        run_name="b5",
        mesh=MeshOptions(cells=(3, 3)),
        discretization=DiscretizationOptions(degree, degree),
        stabilization=StabilizationOptions(
            subscale_space=mode, subscale_dynamics=dynamics, scale_splitting=splitting
        ),
        time=TimeOptions(scheme=scheme, dt=0.02, t_end=0.06),
        thermal=ThermalOptions(convection=form),
        output=OutputOptions(vtk_every=0, diagnostics=BUDGETS),
    )


def insulated_problem(config: Config, *, reference: float = 0.0, passive: bool = True) -> Problem:
    original = build_problem(config)

    def velocity(x: Expr, t: Expr) -> Expr:
        psi = x[0] ** 2 * (1 - x[0]) ** 2 * x[1] ** 2 * (1 - x[1]) ** 2
        return (1 + t) * ufl.as_vector((psi.dx(1), -psi.dx(0)))

    def temperature(x: Expr, t: Expr) -> Expr:
        return 0.4 + (1 + t) * ufl.cos(ufl.pi * x[0]) * ufl.cos(ufl.pi * x[1])

    def pressure(x: Expr, t: Expr) -> Expr:
        return (1 + t) * (x[0] - 0.5)

    mechanical = exact_solution_forcing(velocity, pressure, original.viscosity)

    def forcing(x: Expr, t: Expr) -> Expr:
        return mechanical(x, t) + 0.7 * (temperature(x, t) - reference) * ufl.as_vector((0.0, -1.0))

    def source(x: Expr, t: Expr) -> Expr:
        time = ufl.variable(t)
        field = temperature(x, time)
        return (
            ufl.diff(field, time)
            + ufl.dot(velocity(x, time), ufl.grad(field))
            - 0.2 * ufl.div(ufl.grad(field))
        )

    thermal = ThermalProblem(
        0.2,
        0.7,
        (0.0, -1.0),
        reference,
        source,
        at_initial_time(temperature),
        {name: DiffusiveFlux() for name in original.boundary_tags},
        temperature,
    )
    tracer = PassiveScalarProblem(
        "q", 0.2, source, thermal.initial_temperature, thermal.boundary_conditions, temperature
    )
    return replace(
        original,
        forcing=forcing,
        initial_velocity=at_initial_time(velocity),
        exact_velocity=velocity,
        exact_pressure=pressure,
        thermal=thermal,
        passive_scalars=(tracer,) if passive else (),
    )


def read_budgets(config: Config, problem: Problem, folder: Path) -> list[dict[str, str]]:
    stepper = TimeStepper(config, problem, folder, logging.getLogger("b5"))
    stepper.run()
    with (folder / "diagnostics.csv").open() as stream:
        return list(csv.DictReader(stream))


@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize("dynamics", ["static", "dynamic"])
@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("form", ["conservative", "skew_symmetric", "nonconservative"])
@pytest.mark.parametrize("momentum", ["skew_symmetric_1", "nonconservative"])
def test_coupled_budget_variant_matrix(
    tmp_path: Path, mode: str, dynamics: str, splitting: str, form: str, momentum: str
) -> None:
    config = replace(
        configuration(mode, dynamics, splitting, form), convection=ConvectionOptions(momentum)
    )
    rows = read_budgets(config, insulated_problem(config), tmp_path)
    for row in rows:
        assert abs(float(row["boussinesq_energy_budget_residual"])) < 2e-10
        assert abs(float(row["boussinesq_energy_budget_normalized_residual"])) < 1e-8
        for name in ("temperature", "q"):
            for group in ("thermal_content_budget", "thermal_variance_budget"):
                assert abs(float(row[f"{group}_{name}_residual"])) < 2e-10
            assert abs(float(row[f"thermal_content_budget_{name}_continuity_work"])) < 2e-10
        exchange = float(row["boussinesq_energy_budget_exchange_total"])
        measured = float(row["boussinesq_energy_budget_buoyancy_work"]) + float(
            row["boussinesq_energy_budget_potential_advection"]
        )
        assert abs(exchange - measured) < 1e-12
        if form == "conservative" and mode == "oss" and splitting == "linear":
            assert abs(exchange) < 2e-10
        if form == "conservative" and mode == "asgs" and splitting == "linear":
            assert abs(float(row["boussinesq_energy_budget_exchange_linear"])) > 1e-8


@pytest.mark.parametrize("scheme", ["bdf1", "bdf2", "bdf3", "crank_nicolson"])
@pytest.mark.parametrize("mode", ["asgs", "oss"])
@pytest.mark.parametrize("degree", [1, 2])
def test_budgets_use_actual_history_and_higher_degree_diffusion(
    tmp_path: Path, scheme: str, mode: str, degree: int
) -> None:
    config = configuration(mode, scheme=scheme, degree=degree)
    rows = read_budgets(config, insulated_problem(config, reference=0.3), tmp_path)
    for row in rows:
        assert abs(float(row["boussinesq_energy_budget_residual"])) < 2e-10
        assert abs(float(row["thermal_variance_budget_temperature_residual"])) < 2e-10
    if degree == 2:
        assert abs(float(rows[-1]["thermal_variance_budget_temperature_diffusion_cross"])) > 1e-10


def test_nonhomogeneous_dirichlet_variance_reaction(tmp_path: Path) -> None:
    config = configuration()
    config = replace(
        config,
        output=replace(
            config.output, diagnostics=("thermal_variance_budget", "thermal_boundary_flux")
        ),
    )
    rows = read_budgets(config, build_problem(config), tmp_path)
    assert abs(float(rows[-1]["thermal_variance_budget_temperature_residual"])) < 1e-10
    assert abs(float(rows[-1]["thermal_variance_budget_temperature_reaction_work"])) > 1e-4


@pytest.mark.parametrize("group", ["thermal_content_budget", "boussinesq_energy_budget"])
def test_inadmissible_dirichlet_tests_are_rejected(tmp_path: Path, group: str) -> None:
    config = replace(configuration(), output=OutputOptions(diagnostics=(group,)))
    with pytest.raises(ValueError, match="admissible"):
        TimeStepper(config, build_problem(config), tmp_path, logging.getLogger("b5"))


def test_potential_rejects_gravity_direction_periodicity(tmp_path: Path) -> None:
    config = configuration()
    problem = insulated_problem(config, passive=False)
    mesh, tags = periodic_rectangle((4, 4), (0.0, 0.0), (1.0, 1.0))
    problem = replace(
        problem,
        mesh=mesh,
        boundary_tags=tags,
        boundary_conditions={},
        thermal=replace(problem.thermal, boundary_conditions={}),
    )
    with pytest.raises(ValueError, match="single-valued"):
        TimeStepper(config, problem, tmp_path, logging.getLogger("b5"))


def test_prescribed_flux_and_source_work_in_cn_budget(tmp_path: Path) -> None:
    config = configuration("asgs", scheme="crank_nicolson", form="conservative")
    original = insulated_problem(config, passive=False)

    def temperature(x: Expr, t: Expr) -> Expr:
        return (1 + t) * (1 + x[0] + x[1] + x[0] * x[1])

    def source(x: Expr, t: Expr) -> Expr:
        time = ufl.variable(t)
        field = temperature(x, time)
        return ufl.diff(field, time) + ufl.dot(original.exact_velocity(x, time), ufl.grad(field))

    normals = {"x_min": (-1.0, 0.0), "x_max": (1.0, 0.0), "y_min": (0.0, -1.0), "y_max": (0.0, 1.0)}

    def flux(normal: tuple[float, float]) -> DiffusiveFlux:
        return DiffusiveFlux(
            lambda x, t: -0.2 * ufl.dot(ufl.grad(temperature(x, t)), ufl.as_vector(normal))
        )

    mechanical = exact_solution_forcing(
        original.exact_velocity, original.exact_pressure, original.viscosity
    )

    def forcing(x: Expr, t: Expr) -> Expr:
        return mechanical(x, t) + 0.7 * temperature(x, t) * ufl.as_vector((0.0, -1.0))

    problem = replace(
        original,
        forcing=forcing,
        thermal=replace(
            original.thermal,
            source=source,
            initial_temperature=at_initial_time(temperature),
            exact_temperature=temperature,
            boundary_conditions={name: flux(n) for name, n in normals.items()},
        ),
    )
    rows = read_budgets(config, problem, tmp_path)
    for row in rows:
        assert abs(float(row["boussinesq_energy_budget_residual"])) < 1e-10
        assert abs(float(row["thermal_content_budget_temperature_residual"])) < 1e-10
        assert abs(float(row["thermal_variance_budget_temperature_residual"])) < 1e-10
    assert abs(float(rows[-1]["boussinesq_energy_budget_potential_flux"])) > 1e-3


def test_three_dimensional_coupled_budget_smoke(tmp_path: Path) -> None:
    config = configuration(scheme="crank_nicolson", splitting="linear")
    mesh, tags = box((2, 2, 2), (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))

    def velocity(x: Expr, t: Expr) -> Expr:
        psi = (
            x[0] ** 2 * (1 - x[0]) ** 2 * x[1] ** 2 * (1 - x[1]) ** 2 * x[2] ** 2 * (1 - x[2]) ** 2
        )
        return (1 + t) * ufl.as_vector((psi.dx(1), -psi.dx(0), 0.0 * x[2]))

    def pressure(x: Expr, t: Expr) -> Expr:
        return (1 + t) * (x[0] - 0.5)

    def temperature(x: Expr, t: Expr) -> Expr:
        return 0.4 + (1 + t) * ufl.cos(ufl.pi * x[0]) * ufl.cos(ufl.pi * x[1]) * ufl.cos(
            ufl.pi * x[2]
        )

    mechanical = exact_solution_forcing(velocity, pressure, 0.1)

    def forcing(x: Expr, t: Expr) -> Expr:
        return mechanical(x, t) + 0.7 * temperature(x, t) * ufl.as_vector((0.0, 0.0, -1.0))

    def source(x: Expr, t: Expr) -> Expr:
        time = ufl.variable(t)
        field = temperature(x, time)
        return (
            ufl.diff(field, time)
            + ufl.dot(velocity(x, time), ufl.grad(field))
            - 0.2 * ufl.div(ufl.grad(field))
        )

    problem = Problem(
        mesh,
        {name: Dirichlet(velocity) for name in tags},
        tags,
        0.1,
        forcing,
        at_initial_time(velocity),
        velocity,
        pressure,
        thermal=ThermalProblem(
            0.2,
            0.7,
            (0.0, 0.0, -1.0),
            0.0,
            source,
            at_initial_time(temperature),
            {name: DiffusiveFlux() for name in tags},
            temperature,
        ),
    )
    rows = read_budgets(config, problem, tmp_path)
    assert abs(float(rows[-1]["boussinesq_energy_budget_residual"])) < 1e-10
    assert abs(float(rows[-1]["thermal_variance_budget_temperature_residual"])) < 1e-10


@pytest.mark.parametrize("height", [1.0, 8.0])
def test_cavity_scaling_and_exact_conduction_nusselt(tmp_path: Path, height: float) -> None:
    config = replace(
        configuration(),
        case="differentially_heated_cavity",
        case_parameters={"rayleigh": 1000.0, "prandtl": 0.71, "aspect_ratio": height},
        output=OutputOptions(vtk_every=0, diagnostics=("cavity_nusselt", "temperature_extrema")),
    )
    problem = build_problem(config)
    assert problem.viscosity == pytest.approx((0.71 / 1000) ** 0.5)
    assert problem.thermal.diffusivity == pytest.approx(1 / (0.71 * 1000) ** 0.5)
    problem = replace(problem, thermal=replace(problem.thermal, expansion_coefficient=0.0))
    rows = read_budgets(config, problem, tmp_path)
    for row in rows:
        for name in ("volume", "hot_gradient", "cold_gradient", "hot_reaction", "cold_reaction"):
            assert float(row[f"cavity_nusselt_{name}"]) == pytest.approx(1.0, abs=1e-12)
        assert float(row["temperature_extrema_nodal_min"]) == pytest.approx(-0.5)
        assert float(row["temperature_extrema_nodal_max"]) == pytest.approx(0.5)


@pytest.mark.parametrize(
    "parameters",
    [
        {"rayleigh": -1},
        {"prandtl": False},
        {"aspect_ratio": 0},
        {"rayleigh": float("nan")},
        {"initial_profile": "invalid"},
        {"unknown": 1},
    ],
)
def test_invalid_cavity_data_are_rejected(parameters: dict[str, object]) -> None:
    config = replace(
        configuration(), case="differentially_heated_cavity", case_parameters=parameters
    )
    with pytest.raises(ValueError):
        build_problem(config)


def test_periodic_dynamic_oss_variance_identity(tmp_path: Path) -> None:
    config = configuration(scheme="crank_nicolson")
    mesh, tags = periodic_rectangle((4, 4), (0.0, 0.0), (1.0, 1.0))

    def zero(x: Expr, _t: Expr) -> Expr:
        return ufl.as_vector((0.0 * x[0], 0.0 * x[1]))

    thermal = ThermalProblem(
        0.2,
        0.7,
        (0.0, 0.0),
        0.0,
        lambda x, _t: ufl.sin(2 * ufl.pi * x[0]),
        lambda x: 0.4 + ufl.cos(2 * ufl.pi * x[0]),
        {},
    )
    problem = Problem(
        mesh, {}, tags, 0.1, zero, lambda x: zero(x, ufl.as_ufl(0.0)), thermal=thermal
    )
    rows = read_budgets(config, problem, tmp_path)
    assert float(rows[-1]["thermal_subscale_l2_temperature"]) > 1e-6
    assert abs(float(rows[-1]["thermal_variance_budget_temperature_residual"])) < 1e-10
    assert abs(float(rows[-1]["boussinesq_energy_budget_residual"])) < 1e-12
