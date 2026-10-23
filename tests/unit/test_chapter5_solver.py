"""Chapter 5 forms in the solver: EMAC Newton Jacobians and one short run per form."""

import csv
import logging
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
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
    ViscosityOptions,
)
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.formulation.stabilization_parameters import stabilization_parameters
from navier_stokes_vms.solver.time_stepper import TimeStepper

LOGGER = logging.getLogger(__name__)

FORMS = (
    "nonconservative",
    "conservative",
    "skew_symmetric_1",
    "skew_symmetric_2",
    "emac",
    "emac_vms",
)

ROOT = Path(__file__).resolve().parents[2]


def chapter5_config(
    case: str,
    mesh: MeshOptions,
    degree: int,
    time: TimeOptions,
    form: str,
    viscous: str,
    space: str,
    splitting: str,
    solver: SolverOptions,
) -> Config:
    return Config(
        case=case,
        run_name="chapter5",
        mesh=mesh,
        discretization=DiscretizationOptions(velocity_degree=degree, pressure_degree=degree),
        time=time,
        stabilization=StabilizationOptions(subscale_space=space, scale_splitting=splitting),
        convection=ConvectionOptions(form=form),
        viscosity=ViscosityOptions(form=viscous),
        solver=solver,
        output=OutputOptions(vtk_every=0),
    )


def emac_reference(
    form: str,
    viscous: str,
    stepper: TimeStepper,
    a_tau: Expr,
    nonlinear: bool,
    asgs: bool,
) -> ufl.Form:
    """Dyn-ASGS/OSS residual written from spec: EMAC–VMS and Convective forms, not from vms.py.

    N and τ use the stored advector a* = u + χ_N ũ*; T̃ uses the substituted ũ in every slot.
    """
    problem = stepper._problem
    nu, dx = problem.viscosity, stepper._dx
    solution = stepper._solution
    u, pressure, *projections = ufl.split(solution)
    v, q, *projection_tests = fd.TestFunctions(solution.function_space())
    alpha0, dt = stepper._alpha0, fd.Constant(stepper._config.time.dt)
    f = problem.forcing(fd.SpatialCoordinate(problem.mesh), stepper._t_new)
    u_hist, s_hist = stepper._velocity_history, stepper._subscales.history
    s_star = stepper._iterate.value if nonlinear else ufl.zero(u.ufl_shape)
    a = u + s_star
    tau_m, tau_c = stabilization_parameters(a_tau, nu, 2, 4.0, 2.0)
    tau_t = 1 / (alpha0 / dt + 1 / tau_m)

    def strong_viscous(w: Expr) -> Expr:
        if viscous == "laplacian":
            return nu * ufl.div(ufl.grad(w))
        return 2 * nu * ufl.div(ufl.sym(ufl.grad(w)))

    if viscous == "laplacian":
        viscous_energy = nu * ufl.inner(ufl.grad(u), ufl.grad(v))
    else:
        viscous_energy = 2 * nu * ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v)))
    dt_u = (alpha0 * u - u_hist) / dt
    residual_convection = ufl.dot(ufl.grad(u), a) + ufl.dot(u, ufl.grad(u))
    x = f - residual_convection + strong_viscous(u) - ufl.grad(pressure) + s_hist / dt
    if asgs:
        x -= dt_u
    xi, zeta = (ufl.zero(u.ufl_shape), 0) if asgs else projections
    s = tau_t * (x - xi)
    if form == "emac":
        resolved = ufl.inner(2 * ufl.dot(ufl.sym(ufl.grad(u)), a), v)
        resolved += ufl.div(u) * ufl.inner(u, v) - ufl.inner(s_star, ufl.grad(ufl.inner(u, v)))
        coupling = -ufl.inner(s, ufl.dot(ufl.grad(v), a))
    else:
        resolved = ufl.inner(2 * ufl.dot(ufl.sym(ufl.grad(u)), u), v)
        resolved += ufl.div(u) * ufl.inner(u, v)
        stress = ufl.outer(s, u) + ufl.outer(u, s)
        if nonlinear:
            stress += ufl.outer(s, s)
        coupling = -ufl.inner(stress, ufl.sym(ufl.grad(v)))
    integrand = (
        ufl.inner(dt_u, v)
        + viscous_energy
        + resolved
        - pressure * ufl.div(v)
        + q * ufl.div(u)
        - ufl.inner(f, v)
        + coupling
        - ufl.inner(s, strong_viscous(v) + ufl.grad(q))
        + tau_c * (ufl.div(u) - zeta) * ufl.div(v)
    )
    if asgs:
        integrand += ufl.inner((alpha0 * s - s_hist) / dt, v)
    else:
        w, r = projection_tests
        integrand += tau_t * ufl.inner(xi - x, w) + tau_c * (zeta - ufl.div(u)) * r
    return integrand * dx


@pytest.mark.parametrize("newton_tau", ["frozen", "exact"])
@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize(
    ("form", "viscous"), [("emac", "laplacian"), ("emac", "stress"), ("emac_vms", "stress")]
)
def test_emac_newton_jacobian_matches_finite_differences(
    tmp_path: Path, form: str, viscous: str, space: str, splitting: str, newton_tau: str
) -> None:
    """P2/P2 dynamic subscales, random state; for OSS the rows include the projections."""
    config = chapter5_config(
        "manufactured_spatial",
        MeshOptions(cells=(3, 3)),
        2,
        TimeOptions(scheme="bdf1", dt=0.05, t_end=0.05),
        form,
        viscous,
        space,
        splitting,
        SolverOptions(newton_tau=newton_tau),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    stepper._initialise()
    stepper._set_level(1, step_coefficients("bdf1", "lower_order", 1))
    solution = stepper._solution
    direction = fd.Function(solution.function_space())
    rng = np.random.default_rng(57)
    for field in (*solution.subfunctions, *direction.subfunctions):
        field.dat.data[:] = rng.uniform(0.2, 1.0, field.dat.data.shape)
    stored = [stepper._velocity_history, stepper._subscales.history]
    if splitting == "nonlinear":
        stored.append(stepper._iterate.value)
    for field in stored:
        field.dat.data[:] = rng.uniform(-0.5, 0.5, field.dat.data.shape)
    base = solution.copy(deepcopy=True)
    u = ufl.split(solution)[0]
    nonlinear, asgs = splitting == "nonlinear", space == "asgs"
    a_tau = base.subfunctions[0] if newton_tau == "frozen" else u
    if nonlinear:
        a_tau = a_tau + stepper._iterate.value
    reference = emac_reference(form, viscous, stepper, a_tau, nonlinear, asgs)
    newton = stepper._step_solver._inner._newton._problem
    with fd.assemble(newton.F).dat.vec_ro as actual, fd.assemble(reference).dat.vec_ro as expected:
        assert actual.norm() > 0
        difference = actual.copy()
        difference.axpy(-1, expected)
        assert difference.norm() <= 1e-12 * expected.norm()
    action = fd.assemble(ufl.action(newton.J, direction))
    epsilon = 1e-6
    solution.assign(base + epsilon * direction)
    plus = fd.assemble(reference)
    solution.assign(base - epsilon * direction)
    minus = fd.assemble(reference)
    with plus.dat.vec_ro as vp, minus.dat.vec_ro as vm, action.dat.vec_ro as derivative:
        difference = vp.copy()
        difference.axpy(-1, vm)
        difference.scale(1 / (2 * epsilon))
        difference.axpy(-1, derivative)
        assert difference.norm() <= 1e-8 * derivative.norm()


def assert_every_step_converged(folder: Path, steps: int, stepper: TimeStepper) -> None:
    with (folder / "time_steps.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert [row["converged"] for row in rows] == ["True"] * steps
    for field in (stepper.velocity, stepper.pressure):
        assert np.all(np.isfinite(field.dat.data_ro))


@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize("form", FORMS)
def test_each_form_runs_and_converges_on_a_periodic_mesh(
    tmp_path: Path, form: str, space: str
) -> None:
    config = chapter5_config(
        "taylor_green_2d",
        MeshOptions(cells=(4, 4)),
        1,
        TimeOptions(scheme="crank_nicolson", dt=0.01, t_end=0.03),
        form,
        "stress",
        space,
        "nonlinear",
        SolverOptions(),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    stepper.run()
    assert_every_step_converged(tmp_path, 3, stepper)


@pytest.mark.parametrize("form", FORMS)
def test_each_form_runs_and_converges_with_an_outflow_boundary(tmp_path: Path, form: str) -> None:
    config = chapter5_config(
        "cylinder_dfg_2d",
        MeshOptions(file=str(ROOT / "meshes/cylinder_dfg_2d_coarse.msh")),
        1,
        TimeOptions(scheme="crank_nicolson", dt=0.005, t_end=0.01),
        form,
        "stress",
        "oss",
        "nonlinear",
        SolverOptions(),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    stepper.run()
    assert_every_step_converged(tmp_path, 2, stepper)


@pytest.mark.parametrize("form", FORMS[:4])
def test_picard_converges_for_the_first_four_forms(tmp_path: Path, form: str) -> None:
    """Picard freezes a in c and in the divergence pairing (spec: Convective forms)."""
    config = chapter5_config(
        "taylor_green_2d",
        MeshOptions(cells=(4, 4)),
        1,
        TimeOptions(scheme="crank_nicolson", dt=0.01, t_end=0.03),
        form,
        "laplacian",
        "oss",
        "nonlinear",
        SolverOptions(linearization="picard"),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    stepper.run()
    assert_every_step_converged(tmp_path, 3, stepper)
