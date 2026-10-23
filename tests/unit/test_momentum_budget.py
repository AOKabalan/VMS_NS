"""Crank–Nicolson linear-momentum closure on a periodic domain (spec: Momentum budget).

e_k is an admissible test on the periodic Taylor–Green mesh, so r_M,k = rate + D_k − F_k
vanishes to solver tolerance for every form from step 2. A uniform background flow and a uniform
force make the rate, the defect and the force all nonzero.
"""

import csv
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
AREA = 4 * math.pi**2


def periodic_config(form: str, space: str, splitting: str) -> Config:
    return Config(
        case="taylor_green_2d",
        run_name="momentum",
        mesh=MeshOptions(cells=(8, 8)),
        time=TimeOptions(scheme="crank_nicolson", dt=DT, t_end=3 * DT),
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
        output=OutputOptions(
            vtk_every=0, diagnostics=("momentum_budget", "linear_momentum_defect")
        ),
    )


@pytest.mark.parametrize("splitting", ["linear", "nonlinear"])
@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize("form", FORMS)
def test_periodic_momentum_closure_vanishes_from_step_2(
    tmp_path: Path, form: str, space: str, splitting: str
) -> None:
    config = periodic_config(form, space, splitting)
    problem = build_problem(config)
    taylor_green = problem.initial_velocity
    problem = replace(
        problem,
        initial_velocity=lambda x: taylor_green(x) + ufl.as_vector((0.3, -0.2)),
        forcing=lambda _x, t: ufl.as_vector((1 + t, -0.5)),
    )
    stepper = TimeStepper(config, problem, tmp_path, LOGGER)
    asgs = space == "asgs"
    moments = []

    def record(_step: int, _time: float) -> bool:
        u, s, dx = stepper.velocity, stepper._subscales.latest, stepper._dx
        resolved = [fd.assemble(u[k] * dx) for k in range(2)]
        subscale = [fd.assemble(s[k] * dx) for k in range(2)]
        moments.append((resolved, subscale))
        return False

    stepper.run(record)
    with (tmp_path / "diagnostics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 3
    assert all(rows[0][f"momentum_budget_{term}_x"] == "" for term in ("rate", "residual"))
    for step in (1, 2):
        row = rows[step]
        (resolved_old, subscale_old), (resolved_new, subscale_new) = moments[step - 1 : step + 1]
        midpoint_time = (step + 0.5) * DT
        for k, axis in enumerate("xy"):
            resolved_rate = (resolved_new[k] - resolved_old[k]) / DT
            subscale_rate = (subscale_new[k] - subscale_old[k]) / DT
            rate = resolved_rate + subscale_rate if asgs else resolved_rate
            force = AREA * ((1 + midpoint_time) if axis == "x" else -0.5)
            defect = float(row[f"linear_momentum_defect_{axis}"])
            scale = max(abs(rate), abs(force), abs(defect))
            assert float(row[f"momentum_budget_rate_{axis}"]) == pytest.approx(
                rate, rel=1e-10, abs=1e-12 * scale
            )
            assert float(row[f"momentum_budget_forcing_{axis}"]) == pytest.approx(force, rel=1e-12)
            assert float(row[f"momentum_budget_defect_{axis}"]) == defect
            assert abs(float(row[f"momentum_budget_residual_{axis}"])) <= 1e-10 * scale
            if asgs:
                # The ASGS subscale moment is part of the discrete balance.
                assert abs(subscale_rate) > 1e-6 * scale
            else:
                # Matched OSS projection: e_k is in the unconstrained space, so ∫ũ = 0.
                assert abs(subscale_new[k]) <= 1e-10 * scale * DT
