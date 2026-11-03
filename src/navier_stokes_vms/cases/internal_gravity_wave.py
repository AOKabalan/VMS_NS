"""Stable Boussinesq standing wave with a passive moist layer and no imposed shear."""

import math
from dataclasses import asdict, dataclass, replace

import firedrake as fd
import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, Dirichlet, FreeSlip
from navier_stokes_vms.problem.mesh_loading import builtin_cells, periodic_channel
from navier_stokes_vms.problem.problem_definition import HorizontalProfileGeometry, Problem
from navier_stokes_vms.problem.scalar_problem import PassiveScalarProblem
from navier_stokes_vms.problem.thermal_problem import ThermalProblem


@dataclass(frozen=True)
class WaveParameters:
    """Dimensionless stable gradient, wave kick, common diffusion and moisture profile."""

    n_waves: int = 4
    wave_number: float = 1.0
    height: float = 2 * math.pi
    temperature_gradient: float = 1.0
    buoyancy_coefficient: float = 1.0
    displacement_amplitude: float = 0.2
    diffusivity: float = 0.01
    moisture_diffusivity: float = 0.0002
    moisture: bool = True
    rh_peak: float = 0.98
    moisture_width: float = 1.2
    saturation_height_rate: float = 0.3


def checked_parameters(config: Config) -> WaveParameters:
    """Check every physical input before constructing the wave or mesh."""
    defaults = WaveParameters()
    if set(config.case_parameters) - asdict(defaults).keys():
        raise ValueError("unknown internal_gravity_wave case_parameters")
    values = replace(defaults, **config.case_parameters)
    if (
        isinstance(values.n_waves, bool)
        or not isinstance(values.n_waves, int)
        or values.n_waves < 1
    ):
        raise ValueError("case_parameters.n_waves must be a positive integer")
    if not isinstance(values.moisture, bool):
        raise ValueError("case_parameters.moisture must be boolean")
    for name, value in asdict(values).items():
        if name in ("n_waves", "moisture"):
            continue
        if (
            isinstance(value, bool)
            or not isinstance(value, int | float)
            or not math.isfinite(value)
        ):
            raise ValueError(f"case_parameters.{name} must be finite")
        if value < 0 or (value == 0 and name != "displacement_amplitude"):
            raise ValueError(f"case_parameters.{name} has an invalid range")
    if values.rh_peak > 1:
        raise ValueError("case_parameters.rh_peak must be at most one")
    if values.displacement_amplitude * math.pi / values.height >= 0.5:
        raise ValueError("wave kick must satisfy amplitude*pi/height < 0.5")
    return values


def frequency(parameters: WaveParameters) -> float:
    """Nonrotating linear wave frequency ω²=N²k²/(k²+m²), m=π/H."""
    k, m = parameters.wave_number, math.pi / parameters.height
    return (
        math.sqrt(parameters.buoyancy_coefficient * parameters.temperature_gradient)
        * k
        / math.hypot(k, m)
    )


def moisture_upper_bound(parameters: WaveParameters) -> float:
    """Continuum maximum of the log-concave Gaussian/cosine initial tracer profile."""
    low, high = -parameters.height / 2, parameters.height / 2
    m = math.pi / parameters.height
    for _ in range(64):
        z = (low + high) / 2
        slope = (
            -2 * z / parameters.moisture_width**2
            - parameters.saturation_height_rate
            - 2 * m * math.tan(m * z)
        )
        if slope > 0:
            low = z
        else:
            high = z
    z = (low + high) / 2
    return (
        parameters.rh_peak
        * math.exp(-((z / parameters.moisture_width) ** 2) - parameters.saturation_height_rate * z)
        * math.cos(m * z) ** 2
    )


def linear_fields(x: Expr, t: Expr, parameters: WaveParameters) -> tuple[Expr, Expr, Expr]:
    """Damped linear reference (u,p,T); not an exact solution of the full nonlinear problem."""
    k, m = parameters.wave_number, math.pi / parameters.height
    omega = frequency(parameters)
    amplitude = parameters.displacement_amplitude
    decay = ufl.exp(-parameters.diffusivity * (k * k + m * m) * t)
    phase = ufl.cos(k * x[0]) * ufl.cos(m * x[1])
    velocity = (
        amplitude
        * omega
        * decay
        * ufl.cos(omega * t)
        * ufl.as_vector(
            (
                m / k * ufl.sin(k * x[0]) * ufl.sin(m * x[1]),
                phase,
            )
        )
    )
    temperature = parameters.temperature_gradient * (
        x[1] - amplitude * decay * ufl.sin(omega * t) * phase
    )
    pressure = parameters.buoyancy_coefficient * parameters.temperature_gradient * x[1] ** 2 / 2
    pressure -= (
        amplitude
        * omega**2
        * m
        / k**2
        * decay
        * ufl.sin(omega * t)
        * ufl.cos(k * x[0])
        * ufl.sin(m * x[1])
    )
    return velocity, pressure, temperature


def build(config: Config) -> Problem:
    """Build a free-slip periodic channel, fixed-temperature walls and source-free moist tracer."""
    parameters = checked_parameters(config)
    length = parameters.n_waves * 2 * math.pi / parameters.wave_number
    cells = builtin_cells(config, (24 * parameters.n_waves, 32))
    mesh, tags = periodic_channel(
        cells, (0.0, -parameters.height / 2), (length, parameters.height / 2)
    )
    x = fd.SpatialCoordinate(mesh)
    velocity, pressure, temperature = linear_fields(x, ufl.as_ufl(0.0), parameters)
    m = math.pi / parameters.height
    moisture = (
        parameters.rh_peak
        * ufl.exp(
            -((x[1] / parameters.moisture_width) ** 2) - parameters.saturation_height_rate * x[1]
        )
        * ufl.cos(m * x[1]) ** 2
    )

    def scalar_zero(position: Expr, _time: Expr) -> Expr:
        return 0.0 * position[0]

    def background(position: Expr, _time: Expr) -> Expr:
        return parameters.temperature_gradient * position[1]

    def forcing(position: Expr, _time: Expr) -> Expr:
        return ufl.as_vector((0.0 * position[0], 0.0 * position[1]))

    omega = frequency(parameters)
    metadata = {
        "parameters": {**asdict(parameters), "saturation_temperature_rate": 0.0},
        "domain_length": length,
        "mesh_cells": list(cells),
        "buoyancy_frequency_squared": parameters.buoyancy_coefficient
        * parameters.temperature_gradient,
        "linear_wave_frequency": omega,
        "linear_wave_period": 2 * math.pi / omega,
        "linear_decay_rate": parameters.diffusivity * (parameters.wave_number**2 + m * m),
        "initial_moisture_upper_bound": moisture_upper_bound(parameters),
        "background": "T=gradient*z; p=beta*gradient*z^2/2; no mean flow or shear",
        "thermal_walls": "fixed background temperature; maintains the stable linear gradient",
        "moisture_profile": "RH_peak*exp(-(z/width)^2-gamma*z)*cos(pi*z/H)^2",
        "moisture_wall_adjustment": "cosine-squared envelope gives zero value and normal gradient",
        "initial_phase": "zero displacement, maximum wave velocity; subscales zero",
        "visualization": {
            "title": "Stable internal gravity wave · four-wave moisture pilot"
            if parameters.n_waves == 4
            else "Stable internal gravity wave",
            "temperature_levels": [-1.0, -0.5, 0.0, 0.5, 1.0],
            "vertical_limits": [-2.0, 2.0],
            "vorticity_limit": 0.3,
            "motion_title": "Oscillating wave motion: vorticity",
        },
    }
    return Problem(
        mesh=mesh,
        boundary_conditions={name: FreeSlip(1) for name in tags},
        boundary_tags=tags,
        viscosity=parameters.diffusivity,
        forcing=forcing,
        initial_velocity=lambda _x: velocity,
        initial_pressure=lambda _x: pressure,
        thermal=ThermalProblem(
            parameters.diffusivity,
            parameters.buoyancy_coefficient,
            (0.0, -1.0),
            0.0,
            scalar_zero,
            lambda _x: temperature,
            {name: Dirichlet(background) for name in tags},
            scaling="nondimensional",
        ),
        passive_scalars=(
            PassiveScalarProblem(
                "q",
                parameters.moisture_diffusivity,
                scalar_zero,
                lambda _x: moisture,
                {name: DiffusiveFlux() for name in tags},
            ),
        )
        if parameters.moisture
        else (),
        case_metadata=metadata,
        reference_values={
            "length_scale": 1 / parameters.wave_number,
            "velocity_scale": math.sqrt(
                parameters.buoyancy_coefficient * parameters.temperature_gradient
            )
            / parameters.wave_number,
            "temperature_scale": parameters.temperature_gradient / parameters.wave_number,
            "wave_frequency": omega,
            "wave_period": 2 * math.pi / omega,
        },
        horizontal_profile_geometry=HorizontalProfileGeometry(
            length,
            cells[0],
            tuple(
                -parameters.height / 2 + j * parameters.height / cells[1]
                for j in range(cells[1] + 1)
            ),
        ),
    )
