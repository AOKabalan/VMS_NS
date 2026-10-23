"""Convective defects and the Crank–Nicolson energy closure (spec: Conservation checks, T5).

Every run has u_h = 0 on the whole boundary, so all convective fluxes vanish and the rigid tests
are admissible for the algebraic defects; no angular-momentum check runs on a periodic domain.
"""

import csv
import logging
from pathlib import Path
from typing import Any

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
DT = 0.05
K = 2
DIAGNOSTICS = (
    "energy_budget",
    "convective_energy_defect",
    "linear_momentum_defect",
    "angular_momentum_defect",
    "linear_momentum_defect_predicted",
    "angular_momentum_defect_predicted",
)


def budget_config(form: str, space: str, splitting: str, t_end: float) -> Config:
    return Config(
        case="manufactured_spatial",
        run_name="budget",
        mesh=MeshOptions(cells=(4, 4)),
        discretization=DiscretizationOptions(velocity_degree=K, pressure_degree=K),
        time=TimeOptions(scheme="crank_nicolson", dt=DT, t_end=t_end),
        stabilization=StabilizationOptions(subscale_space=space, scale_splitting=splitting),
        convection=ConvectionOptions(form=form),
        viscosity=ViscosityOptions(form="stress"),
        solver=SolverOptions(
            inner_rtol=1e-12,
            inner_atol=1e-14,
            outer_rtol=1e-11,
            outer_max_iterations=100,
            subscale_rtol=1e-13,
        ),
        output=OutputOptions(vtk_every=0, diagnostics=DIAGNOSTICS),
    )


def rows(folder: Path) -> list[dict[str, str]]:
    with (folder / "diagnostics.csv").open() as stream:
        return list(csv.DictReader(stream))


def energy(u: Expr, s: Expr, asgs: bool, dx: ufl.Measure) -> float:
    if asgs:
        return fd.assemble(ufl.inner(u + s, u + s) / 2 * dx)
    return fd.assemble((ufl.inner(u, u) + ufl.inner(s, s)) / 2 * dx)


def table_predictions(
    form: str, nonlinear: bool, u: Expr, s: Expr, dx: ufl.Measure
) -> dict[str, float]:
    """Chapter Tables 5.1–5.2 for zero boundary velocity, written out per entry."""

    def pairing(w: Expr) -> float:
        return fd.assemble((ufl.div(u) * w - (ufl.inner(s, ufl.grad(w)) if nonlinear else 0)) * dx)

    x = fd.SpatialCoordinate(u.function_space().mesh())
    phi = ufl.as_vector((-x[1], x[0]))
    factor = {"nonconservative": -1, "skew_symmetric_1": -0.5, "skew_symmetric_2": -0.5}
    predicted = {
        f"linear_momentum_defect_{axis}": factor.get(form, 0) * pairing(u[k])
        for k, axis in enumerate("xy")
    }
    predicted["angular_momentum_defect"] = factor.get(form, 0) * pairing(ufl.inner(u, phi))
    energy_factor = {"nonconservative": -0.5, "conservative": 0.5}.get(form, 0)
    predicted["energy_budget_defect_predicted"] = energy_factor * pairing(ufl.inner(u, u))
    if form == "emac":
        transfer = ufl.dot(ufl.grad(u), u) if nonlinear else ufl.dot(u, ufl.grad(u))
        predicted["energy_budget_defect_predicted"] = fd.assemble(ufl.inner(s, transfer) * dx)
        if nonlinear:
            for k, axis in enumerate("xy"):
                predicted[f"linear_momentum_defect_{axis}"] = fd.assemble(
                    ufl.inner(s, u.dx(k)) * dx
                )
            predicted["angular_momentum_defect"] = fd.assemble(
                ufl.inner(s, ufl.dot(ufl.grad(u), phi)) * dx
            )
    if not nonlinear and form != "emac_vms":
        predicted["angular_momentum_defect"] -= fd.assemble(
            ufl.inner(s, ufl.dot(ufl.grad(phi), u)) * dx
        )
    return predicted


def instrument(
    stepper: TimeStepper, monkeypatch: pytest.MonkeyPatch, midpoint: bool | None = None
) -> list[dict[str, Any]]:
    """Independent references from the state at each capture, before the histories advance.

    midpoint=True forces the closure at the BDF1 startup step as well.
    """
    problem, dx = stepper._problem, stepper._dx
    nu, asgs = problem.viscosity, stepper._config.stabilization.subscale_space == "asgs"
    nonlinear = stepper._iterate is not None
    subscales = stepper._subscales
    conservation = stepper._conservation
    capture = conservation.capture
    x = fd.SpatialCoordinate(problem.mesh)
    references: list[dict[str, Any]] = []

    def checked_capture(is_midpoint: bool) -> None:
        capture(is_midpoint if midpoint is None else midpoint)
        u = stepper._velocity.copy(deepcopy=True)
        s = fd.Function(subscales.latest.function_space()).interpolate(stepper._u_sub)
        a = u + stepper._iterate.value if nonlinear else u
        zeta = 0 if asgs else stepper._solution.subfunctions[3]
        h = ufl.MinCellEdgeLength(problem.mesh)
        tau_m = 1 / (4 * K**4 * nu / h**2 + 2 * K * ufl.sqrt(ufl.dot(a, a)) / h)
        tau_c = h**2 / (4 * K**4 * tau_m)
        tau_t = 1 / (2 / DT + 1 / tau_m)
        t_new, t_old = float(stepper._t_new), float(stepper._t_old)
        weight = 0.5 if is_midpoint else 1.0
        forcing = weight * problem.forcing(x, fd.Constant(t_new))
        forcing += (1 - weight) * problem.forcing(x, fd.Constant(t_old))
        reference = {
            "viscous": fd.assemble(
                2 * nu * ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(u))) * dx
            ),
            "subscale": fd.assemble(ufl.inner(s, s) / tau_m * dx),
            "subscale_with_tau_t": fd.assemble(ufl.inner(s, s) / tau_t * dx),
            "viscous_cross": fd.assemble(
                -2 * nu * ufl.inner(s, ufl.div(ufl.grad(u)) + ufl.grad(ufl.div(u))) * dx
            ),
            "pressure_subscale": fd.assemble(tau_c * (ufl.div(u) - zeta) * ufl.div(u) * dx),
            "power": fd.assemble(ufl.inner(forcing, u + s) * dx),
            "old_energy": energy(stepper.velocity, subscales.latest, asgs, dx),
            "old_velocity": stepper.velocity.copy(deepcopy=True),
            "old_subscale": subscales.latest.copy(deepcopy=True),
        }
        reference.update(table_predictions(stepper._config.convection.form, nonlinear, u, s, dx))
        references.append(reference)

    monkeypatch.setattr(conservation, "capture", checked_capture)
    return references


def close(actual: str, expected: float, scale: float) -> None:
    assert float(actual) == pytest.approx(expected, rel=1e-9, abs=1e-11 * scale)


@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize("form", FORMS)
def test_midpoint_closure_vanishes_from_step_2_T5(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, form: str, space: str, splitting: str
) -> None:
    config = budget_config(form, space, splitting, t_end=3 * DT)
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    references = instrument(stepper, monkeypatch)
    asgs = space == "asgs"
    new_energies = []
    stepper.run(
        lambda _step, _time: new_energies.append(
            energy(stepper.velocity, stepper._subscales.latest, asgs, stepper._dx)
        )
    )
    table = rows(tmp_path)
    assert len(table) == 3
    assert all(table[0][f"energy_budget_{term}"] == "" for term in ("rate", "residual", "power"))
    for row, reference, new_energy in zip(table[1:], references[1:], new_energies[1:], strict=True):
        scale = max(abs(float(row[f"energy_budget_{term}"])) for term in ("rate", "power"))
        rate = (new_energy - reference["old_energy"]) / DT
        close(row["energy_budget_rate"], rate, scale)
        for term in ("viscous", "subscale", "viscous_cross", "pressure_subscale", "power"):
            close(row[f"energy_budget_{term}"], reference[term], scale)
        # Each term is far above the closure tolerance: weighting by τ_t instead of τ_m, or
        # dropping the P2 viscous cross term or the pressure term, would show in r.
        significant = 1e-5 * scale
        assert abs(reference["subscale_with_tau_t"] - reference["subscale"]) > significant
        assert abs(reference["viscous_cross"]) > significant
        assert abs(reference["pressure_subscale"]) > significant
        dissipation = sum(
            reference[term]
            for term in ("viscous", "subscale", "viscous_cross", "pressure_subscale")
        )
        close(row["energy_budget_dissipation"], dissipation, scale)
        predicted = reference["energy_budget_defect_predicted"]
        close(row["energy_budget_defect_predicted"], predicted, scale)
        if form in ("nonconservative", "conservative", "emac"):
            assert abs(predicted) > significant
        else:
            assert float(row["energy_budget_defect_predicted"]) == 0.0
        measured = float(row["energy_budget_defect_measured"])
        assert measured == pytest.approx(predicted, rel=1e-9, abs=1e-10 * scale)
        assert measured == pytest.approx(
            float(row["convective_energy_defect"]), rel=1e-9, abs=1e-10 * scale
        )
        assert abs(float(row["energy_budget_residual"])) <= 1e-10 * scale
        for name in ("linear_momentum_defect_x", "linear_momentum_defect_y"):
            close(row[name], reference[name], scale)
            predicted_name = name.replace("_x", "_predicted_x").replace("_y", "_predicted_y")
            close(row[predicted_name], reference[name], scale)
        close(row["angular_momentum_defect"], reference["angular_momentum_defect"], scale)
        close(row["angular_momentum_defect_predicted"], reference["angular_momentum_defect"], scale)
        if form == "emac_vms":
            assert float(row["linear_momentum_defect_predicted_x"]) == 0.0
            assert float(row["angular_momentum_defect_predicted"]) == 0.0
        # J = (ũ, (u_h·∇)φ) vanishes for OSS: (u_h·∇)φ is a rotated u_h, in the FE space.
        nonzero = form in ("nonconservative", "skew_symmetric_1", "skew_symmetric_2")
        nonzero |= form == "emac" and (splitting == "nonlinear" or space == "asgs")
        nonzero |= form == "conservative" and splitting == "linear" and space == "asgs"
        if nonzero:
            assert abs(reference["angular_momentum_defect"]) > 1e-5 * scale
    if form == "emac_vms":
        final = table[-1]
        for name in ("linear_momentum_defect_x", "linear_momentum_defect_y"):
            assert abs(float(final[name])) <= 1e-12
        assert abs(float(final["angular_momentum_defect"])) <= 1e-12
        assert abs(float(final["convective_energy_defect"])) <= 1e-12


@pytest.mark.parametrize("space", ["asgs", "oss"])
def test_bdf1_startup_adds_the_increment_term(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, space: str
) -> None:
    """Forced at the startup step, the midpoint closure leaves −‖Δ(u + ũ)‖²/(2δt) for ASGS
    and −(‖Δu‖² + ‖Δũ‖²)/(2δt) for OSS, which is why that step is excluded."""
    config = budget_config("skew_symmetric_1", space, "nonlinear", t_end=DT)
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    references = instrument(stepper, monkeypatch, midpoint=True)
    stepper.run()
    (row,) = rows(tmp_path)
    (reference,) = references
    dx = stepper._dx
    du = stepper.velocity - reference["old_velocity"]
    ds = stepper._subscales.latest - reference["old_subscale"]
    assert np.max(np.abs(reference["old_subscale"].dat.data_ro)) == 0.0
    increment = energy(du, ds, space == "asgs", dx) * 2 / (2 * DT)
    scale = abs(float(row["energy_budget_power"]))
    assert increment > 1e-3 * scale
    assert float(row["energy_budget_residual"]) == pytest.approx(-increment, rel=1e-8)
