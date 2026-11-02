"""G4 channel geometry, component constraints, coupled transport and natural traction."""

import csv
import logging
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import pytest
import ufl
from firedrake.exceptions import ConvergenceError
from ufl.core.expr import Expr

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
    ViscosityOptions,
)
from navier_stokes_vms.problem.boundary_conditions import (
    Dirichlet,
    FreeSlip,
    check_free_slip_boundaries,
    velocity_bcs,
)
from navier_stokes_vms.problem.mesh_loading import box, periodic_channel, rectangle
from navier_stokes_vms.problem.problem_definition import check_problem, uses_pressure_null_space
from navier_stokes_vms.solver.time_stepper import TimeStepper


def configuration(
    mode: str = "oss",
    dynamics: str = "dynamic",
    splitting: str = "nonlinear",
    thermal: bool = True,
    viscous: str = "laplacian",
    scheme: str = "bdf2",
) -> Config:
    options = SolverOptions()
    return Config(
        case="free_slip_channel",
        run_name="g4",
        case_parameters={"thermal": thermal},
        mesh=MeshOptions(cells=(4, 4)),
        discretization=DiscretizationOptions(2, 1)
        if mode == "none"
        else DiscretizationOptions(1, 1),
        stabilization=StabilizationOptions(method="none")
        if mode == "none"
        else StabilizationOptions(
            subscale_space=mode, subscale_dynamics=dynamics, scale_splitting=splitting
        ),
        viscosity=ViscosityOptions(viscous),
        time=TimeOptions(scheme=scheme, dt=0.05, t_end=0.15),
        solver=replace(options, petsc={**options.petsc, "mat_mumps_icntl_24": 1}),
        output=OutputOptions(vtk_every=0),
    )


def norm(expression: Expr, measure: ufl.Measure) -> float:
    return fd.assemble(ufl.inner(expression, expression) * measure) ** 0.5


def test_channel_translation_and_boundary_tags() -> None:
    mesh, tags = periodic_channel((4, 3), (-0.75, -1.5), (1.75, 2.5))
    assert tags == {"y_min": 3, "y_max": 4}
    x = fd.SpatialCoordinate(mesh)
    for tag, height in ((3, -1.5), (4, 2.5)):
        measure = fd.ds(tag, domain=mesh, degree=2)
        assert fd.assemble(fd.Constant(1.0) * measure) == pytest.approx(2.5)
        assert fd.assemble((x[1] - height) ** 2 * measure) < 1e-25
        normal = fd.FacetNormal(mesh)
        assert fd.assemble(normal[0] ** 2 * measure) == 0
    for tag in (1, 2):
        assert fd.assemble(fd.Constant(1.0) * fd.ds(tag, domain=mesh)) == 0


@pytest.mark.parametrize("component", [-1, 2, True, 1.0])
def test_invalid_normal_components_are_rejected(component: object) -> None:
    mesh, tags = periodic_channel((4, 3), (0.0, -1.0), (2.0, 1.0))
    with pytest.raises(ValueError, match="normal_component"):
        check_free_slip_boundaries(mesh, {"y_min": FreeSlip(component)}, tags)


def test_wrong_axis_and_tilted_wall_are_rejected() -> None:
    mesh, tags = rectangle((2, 2), (0.0, 0.0), (1.0, 1.0))
    with pytest.raises(ValueError, match="not normal"):
        check_free_slip_boundaries(mesh, {"y_min": FreeSlip(0)}, tags)
    mesh.coordinates.dat.data[:, 0] += 0.2 * mesh.coordinates.dat.data[:, 1]
    with pytest.raises(ValueError, match="not normal"):
        check_free_slip_boundaries(mesh, {"x_min": FreeSlip(0)}, tags)


def test_three_dimensional_axis_aligned_component_bc() -> None:
    mesh, tags = box((2, 2, 2), (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    conditions = {name: FreeSlip(2) for name in ("z_min", "z_max")}
    check_free_slip_boundaries(mesh, conditions, tags)
    space = fd.VectorFunctionSpace(mesh, "CG", 1)
    field = fd.Function(space).assign(fd.Constant((1.0, 2.0, 3.0)))
    for bc in velocity_bcs(space, conditions, {}, tags):
        bc.apply(field)
    assert norm(field[0] - 1, fd.dx(domain=mesh)) < 1e-14
    assert norm(field[1] - 2, fd.dx(domain=mesh)) < 1e-14
    assert fd.assemble(field[2] ** 2 * fd.ds((5, 6), domain=mesh)) == 0


def test_only_normal_velocity_rows_are_constrained(tmp_path: Path) -> None:
    config = configuration()
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger("g4"))
    stepper._solution.assign(7.0)
    inner = stepper._step_solver._inner
    assert len(inner._bcs) == 2
    for bc in inner._bcs:
        bc.apply(stepper._solution)
    fields = stepper._layout.bind(stepper._solution.subfunctions)
    assert norm(fields["velocity"][0] - 7, stepper._dx) < 1e-13
    assert fd.assemble(fields["velocity"][1] ** 2 * stepper._ds) == 0
    assert fd.assemble(fields["velocity"][1] ** 2 * stepper._dx) > 1
    for name in stepper._layout.names[1:]:
        field = fields[name]
        target = ufl.as_vector((7.0,) * field.ufl_shape[0]) if field.ufl_shape else 7.0
        assert norm(field - target, stepper._dx) < 1e-13
    residual = fd.Cofunction(stepper._solution.function_space().dual())
    with residual.dat.vec as vector:
        vector.set(1.0)
    for bc in inner._bcs:
        bc.zero(residual)
    u_residual = residual.subfunctions[0]
    assert (u_residual.dat.data_ro[:, 0] == 1.0).all()
    assert (u_residual.dat.data_ro[:, 1] == 0.0).any()
    assert (u_residual.dat.data_ro[:, 1] == 1.0).any()
    assert stepper._zero_mean_pressure
    assert stepper._layout.index("pressure") == 1
    jacobian = inner._newton._problem.J
    for name in ("pressure", "velocity", "xi"):
        direction = fd.Function(stepper._solution.function_space())
        target = direction.subfunctions[stepper._layout.index(name)]
        target.assign(1.0 if name == "pressure" else fd.Constant((1.0, 0.0)))
        action = fd.assemble(ufl.action(jacobian, direction))
        for bc in inner._bcs:
            bc.zero(action)
        with action.dat.vec_ro as vector:
            size = vector.norm()
        assert size < 1e-12 if name == "pressure" else size > 1e-5


VARIANTS = [("none", "static", "linear")] + [
    (mode, dynamics, splitting)
    for mode in ("asgs", "oss")
    for dynamics in ("static", "dynamic")
    for splitting in ("linear", "nonlinear")
]


@pytest.mark.parametrize("mode,dynamics,splitting", VARIANTS)
@pytest.mark.parametrize(
    "thermal,viscous", [(False, "laplacian"), (False, "stress"), (True, "laplacian")]
)
def test_uniform_channel_preserved_across_variants(
    tmp_path: Path, mode: str, dynamics: str, splitting: str, thermal: bool, viscous: str
) -> None:
    config = configuration(mode, dynamics, splitting, thermal, viscous)
    problem = build_problem(config)
    assert uses_pressure_null_space(problem)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("g4"))
    stepper.run()
    assert norm(stepper.velocity - ufl.as_vector((1.0, 0.0)), stepper._dx) < 1e-11
    assert norm(stepper.pressure, stepper._dx) < 1e-11
    assert (
        fd.assemble(ufl.dot(stepper.velocity, fd.FacetNormal(problem.mesh)) ** 2 * stepper._ds)
        < 1e-22
    )
    assert abs(fd.assemble(stepper.velocity[0] * stepper._dx) - 2.0) < 1e-11
    if thermal:
        assert norm(stepper.temperature - 1, stepper._dx) < 1e-11
        assert norm(stepper.scalar("q") - 2, stepper._dx) < 1e-11
    with (tmp_path / "time_steps.csv").open() as stream:
        assert all(row["converged"] == "True" for row in csv.DictReader(stream))


@pytest.mark.parametrize("scheme", ["bdf1", "bdf2", "bdf3", "crank_nicolson"])
@pytest.mark.parametrize("degree", [1, 2])
def test_thermal_channel_histories_and_unconstrained_scalar_fields(
    tmp_path: Path, scheme: str, degree: int
) -> None:
    config = replace(
        configuration(scheme=scheme),
        discretization=DiscretizationOptions(degree, degree),
        output=OutputOptions(
            vtk_every=1,
            diagnostics=(
                "thermal_content_budget",
                "thermal_variance_budget",
                "boussinesq_energy_budget",
            ),
        ),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, logging.getLogger("g4"))
    stepper.run()
    assert norm(stepper.velocity - ufl.as_vector((1.0, 0.0)), stepper._dx) < 1e-11
    assert norm(stepper.temperature - 1, stepper._dx) < 1e-11
    assert norm(stepper.scalar("q") - 2, stepper._dx) < 1e-11
    assert (tmp_path / "fields" / "q.pvd").exists()
    with (tmp_path / "diagnostics.csv").open() as stream:
        for row in csv.DictReader(stream):
            for key, value in row.items():
                if key.endswith("_residual") and not key.endswith("normalized_residual"):
                    assert abs(float(value)) < 1e-10


def test_shear_decay_has_natural_zero_tangential_traction(tmp_path: Path) -> None:
    config = replace(
        configuration(
            dynamics="static", splitting="linear", thermal=False, scheme="crank_nicolson"
        ),
        mesh=MeshOptions(cells=(8, 8)),
        time=TimeOptions(scheme="crank_nicolson", dt=0.02, t_end=0.1),
    )
    problem = build_problem(config)

    def velocity(x: Expr, t: Expr) -> Expr:
        return ufl.as_vector(
            (
                1
                + 0.2
                * ufl.cos(ufl.pi * (x[1] + 0.5))
                * ufl.exp(-problem.viscosity * ufl.pi**2 * t),
                0.0 * x[0],
            )
        )

    problem = replace(
        problem, initial_velocity=lambda x: velocity(x, ufl.as_ufl(0.0)), exact_velocity=velocity
    )
    check_problem(problem)
    steppers = []
    for viscous in ("laplacian", "stress"):
        folder = tmp_path / viscous
        folder.mkdir()
        stepper = TimeStepper(
            replace(config, viscosity=ViscosityOptions(viscous)),
            problem,
            folder,
            logging.getLogger("g4"),
        )
        stepper.run()
        exact = velocity(fd.SpatialCoordinate(problem.mesh), ufl.as_ufl(0.1))
        assert norm(stepper.velocity - exact, stepper._dx) < 0.005
        assert abs(fd.assemble(stepper.velocity[0] * stepper._dx) - 2.0) < 1e-10
        assert norm(stepper.velocity[1], stepper._dx) < 1e-10
        steppers.append(stepper)
    assert norm(steppers[0].velocity - steppers[1].velocity, steppers[0]._dx) < 1e-10


def test_free_slip_energy_budget_with_nonzero_tangential_wall_speed(tmp_path: Path) -> None:
    config = replace(
        configuration(dynamics="static", splitting="linear", scheme="crank_nicolson"),
        output=OutputOptions(vtk_every=0, diagnostics=("boussinesq_energy_budget",)),
    )
    problem = build_problem(config)

    def initial(x: Expr) -> Expr:
        return ufl.as_vector((1 + 0.2 * ufl.cos(ufl.pi * (x[1] + 0.5)), 0.0 * x[0]))

    problem = replace(problem, initial_velocity=initial, exact_velocity=None, exact_pressure=None)
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("g4"))
    stepper.run()
    assert fd.assemble(stepper.velocity[0] ** 2 * stepper._ds) > 1.0
    with (tmp_path / "diagnostics.csv").open() as stream:
        for row in csv.DictReader(stream):
            assert abs(float(row["boussinesq_energy_budget_residual"])) < 1e-10


def test_slip_initialization_clamps_only_the_normal_trace(tmp_path: Path) -> None:
    config = configuration(thermal=False)
    problem = build_problem(config)
    problem = replace(
        problem,
        initial_velocity=lambda x: ufl.as_vector((1 + 0.1 * x[1], 0.2 + 0.0 * x[0])),
        exact_velocity=None,
    )
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("g4"))
    stepper._initialise()
    x = fd.SpatialCoordinate(problem.mesh)
    assert norm(stepper.velocity[0] - (1 + 0.1 * x[1]), stepper._dx) < 1e-14
    assert fd.assemble(stepper.velocity[1] ** 2 * stepper._ds) == 0
    assert fd.assemble(stepper.velocity[1] ** 2 * stepper._dx) > 0


def test_mixed_dirichlet_slip_flux_incompatibility_is_rejected(tmp_path: Path) -> None:
    config = configuration(thermal=False)
    original = build_problem(config)
    mesh, tags = rectangle((3, 3), (0.0, 0.0), (1.0, 1.0))

    def inflow(x: Expr, _t: Expr) -> Expr:
        return ufl.as_vector((1.0 + 0.0 * x[0], 0.0 * x[1]))

    problem = replace(
        original,
        mesh=mesh,
        boundary_tags=tags,
        exact_velocity=None,
        exact_pressure=None,
        initial_velocity=lambda x: 0.0 * x,
        boundary_conditions={
            "x_min": Dirichlet(inflow),
            "x_max": FreeSlip(0),
            "y_min": FreeSlip(1),
            "y_max": FreeSlip(1),
        },
    )
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("g4"))
    with pytest.raises(ValueError, match="discretely flux-free"):
        stepper.run()


def test_failed_step_keeps_slip_tangent_and_writes_all_physical_fields(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = configuration()
    problem = build_problem(config)
    problem = replace(problem, forcing=lambda x, _t: ufl.as_vector((1.0 + 0.0 * x[0], 0.0 * x[1])))
    stepper = TimeStepper(config, problem, tmp_path, logging.getLogger("g4"))

    def fail() -> None:
        raise ConvergenceError("G4 synthetic factorization failure")

    monkeypatch.setattr(stepper._step_solver._inner._newton, "solve", fail)
    with pytest.raises(ConvergenceError, match="G4 synthetic"):
        stepper.run()
    assert norm(stepper._velocity[0] - 1, stepper._dx) < 1e-13
    assert fd.assemble(stepper._velocity[1] ** 2 * stepper._ds) == 0
    assert norm(stepper.velocity[0] - 1, stepper._dx) < 1e-13
    for name in ("velocity", "pressure", "temperature", "q"):
        assert (tmp_path / "fields" / f"{name}.pvd").exists()
