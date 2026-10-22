"""Dynamic subscale integration against an exact subscale ODE solution (plan U8).

With the FE fields prescribed, u_h = U constant in time and space and p_h a fixed FE pressure, τ is
constant (|a| = |U| and h is the same in every cell), D_t u_h = 0 and the momentum residual is
R = f − ∇p_h. The forcing f = g(t) F(x) + ∇p_h makes the dynamic subscale equation
∂t ũ + ũ/τ = g(t) F(x), with ũ(0) = 0, whose solution is ũ = y(t) F(x), y′ + y/τ = g, y(0) = 0.
The test marches the time stepper's own subscale path (step coefficients, history, the subscale
expression with τ_t, storage and Crank–Nicolson extrapolation) without solving for u_h, p_h.
"""

import logging
import math
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
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.solver.time_stepper import TimeStepper

LOGGER = logging.getLogger(__name__)

T_END = 1.0
STEPS = (10, 20, 40, 80)
SPEED = 0.2
CELLS = 2
VISCOSITY = 0.01  # manufactured_spatial's default
# τ_m = (4ν/h² + 2|U|/h)⁻¹ for P1 with h = 1/CELLS, the shortest edge of every cell.
TAU = 1 / (4 * VISCOSITY * CELLS**2 + 2 * SPEED * CELLS)


def generic(t: Expr | float) -> tuple[Expr, float]:
    """g = cos 2t, so y(0) = 0 but y′(0) = 1."""
    a = TAU / (1 + 4 * TAU**2)
    y = a * (math.cos(2 * T_END) + 2 * TAU * math.sin(2 * T_END) - math.exp(-T_END / TAU))
    return ufl.cos(2 * t), y


def smooth_start(t: Expr | float) -> tuple[Expr, float]:
    """y = sin³ t, which vanishes with its first two derivatives at t = 0."""
    g = 3 * ufl.sin(t) ** 2 * ufl.cos(t) + ufl.sin(t) ** 3 / TAU
    return g, math.sin(T_END) ** 3


def shape(x: Expr) -> Expr:
    return ufl.as_vector((ufl.sin(math.pi * x[0]) + x[1], ufl.cos(x[0] * x[1])))


def final_subscale_error(folder: Path, scheme: str, startup: str, data: str, steps: int) -> float:
    data_function = {"generic": generic, "smooth_start": smooth_start}[data]
    config = Config(
        case="manufactured_spatial",
        run_name="u8",
        mesh=MeshOptions(cells=(CELLS, CELLS)),
        discretization=DiscretizationOptions(velocity_degree=1, pressure_degree=1),
        time=TimeOptions(scheme=scheme, dt=T_END / steps, t_end=T_END, startup=startup),
        stabilization=StabilizationOptions(
            subscale_space="asgs", subscale_dynamics="dynamic", scale_splitting="linear"
        ),
        output=OutputOptions(vtk_every=0),
    )
    base = build_problem(config)
    assert base.viscosity == VISCOSITY
    x = fd.SpatialCoordinate(base.mesh)
    pressure = fd.Function(fd.FunctionSpace(base.mesh, "CG", 1)).interpolate(x[0] * x[1])
    velocity = ufl.as_vector((SPEED, 0.0))
    problem = replace(
        base,
        initial_velocity=lambda _x: velocity,
        exact_velocity=lambda _x, _t: velocity,
        forcing=lambda y, t: data_function(t)[0] * shape(y) + ufl.grad(pressure),
    )
    stepper = TimeStepper(config, problem, folder, LOGGER)
    stepper._initialise()
    for step in range(1, steps + 1):
        coefficients = step_coefficients(scheme, startup, step)
        stepper._set_level(step, coefficients)
        stepper._velocity.interpolate(velocity)
        stepper._pressure.assign(pressure)
        stepper._advance(coefficients)
    computed = stepper._subscales.latest
    exact = fd.Function(computed.function_space()).interpolate(data_function(T_END)[1] * shape(x))
    scale = np.abs(exact.dat.data_ro).max()
    return np.abs(computed.dat.data_ro - exact.dat.data_ro).max() / scale


@pytest.mark.parametrize(
    ("scheme", "startup", "data", "order"),
    [
        ("bdf1", "lower_order", "generic", 1),
        ("bdf2", "lower_order", "generic", 2),
        ("bdf3", "lower_order", "generic", 2),
        # ũ is zero at negative times, while the smooth extension of y has y(−jδt) ≈ −j y′(0) δt.
        # BDF3's first step therefore sees a history wrong by O(δt), which gives ũ¹ an O(δt)
        # error that the stable recursion carries to T: first order, unless y′(0) = y″(0) = 0.
        ("bdf3", "exact_history", "generic", 1),
        ("bdf3", "exact_history", "smooth_start", 3),
        ("crank_nicolson", "lower_order", "generic", 2),
    ],
)
def test_dynamic_subscale_reaches_the_order_of_its_scheme(
    tmp_path: Path, scheme: str, startup: str, data: str, order: int
) -> None:
    errors = [
        final_subscale_error(tmp_path / str(steps), scheme, startup, data, steps) for steps in STEPS
    ]
    rates = [math.log2(coarse / fine) for coarse, fine in zip(errors, errors[1:], strict=False)]
    LOGGER.info("U8 %s %s %s: errors %s, rates %s", scheme, startup, data, errors, rates)
    assert abs(rates[-1] - order) <= 0.2
