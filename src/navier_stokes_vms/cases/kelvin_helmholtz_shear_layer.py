"""Stratified tanh shear layer with compatible walls and a passive moisture proxy (B7.1)."""

import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace

import firedrake as fd
import numpy as np
import ufl
from ufl.core.expr import Expr

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem.boundary_conditions import DiffusiveFlux, FreeSlip
from navier_stokes_vms.problem.mesh_loading import builtin_cells, periodic_channel
from navier_stokes_vms.problem.problem_definition import HorizontalProfileGeometry, Problem
from navier_stokes_vms.problem.scalar_problem import PassiveScalarProblem
from navier_stokes_vms.problem.thermal_problem import ThermalProblem


@dataclass(frozen=True)
class ShearLayerParameters:
    """Dimensionless physical, perturbation and wall-profile parameters."""

    reynolds: float = 1000.0
    prandtl: float = 1.0
    schmidt: float = 1.0
    richardson: float = 0.1
    n_waves: int = 1
    wave_number: float = 0.4446
    height: float = 20.0
    epsilon: float = 1e-4
    epsilon_noise: float = 0.0
    random_seed: int = 1729
    noise_modes: tuple[int, ...] = (1, 2, 3)
    wall_buffer: float = 1.0
    wall_transition_width: float = 2.0
    moisture: bool = True
    rh_peak: float = 0.97
    moisture_width: float = 1.5
    saturation_height_rate: float = 0.03
    saturation_temperature_rate: float = 0.0
    clustering: float = 0.0


def checked_parameters(values: Mapping[str, object]) -> ShearLayerParameters:
    """Validate each parameter before constructing a mesh or evaluating profile formulas."""
    defaults = ShearLayerParameters()
    if set(values) - asdict(defaults).keys():
        raise ValueError("unknown kelvin_helmholtz_shear_layer case_parameters")
    parameters = replace(defaults, **values)
    for name in (
        "reynolds",
        "prandtl",
        "schmidt",
        "richardson",
        "wave_number",
        "height",
        "epsilon",
        "epsilon_noise",
        "wall_buffer",
        "wall_transition_width",
        "rh_peak",
        "moisture_width",
        "saturation_height_rate",
        "saturation_temperature_rate",
        "clustering",
    ):
        value = getattr(parameters, name)
        if (
            isinstance(value, bool)
            or not isinstance(value, int | float)
            or not math.isfinite(value)
        ):
            raise ValueError(f"case_parameters.{name} must be finite")
        if name in (
            "richardson",
            "epsilon",
            "epsilon_noise",
            "saturation_height_rate",
            "clustering",
        ):
            valid = value >= 0
        elif name == "saturation_temperature_rate":
            valid = True
        else:
            valid = value > 0
        if not valid:
            raise ValueError(f"case_parameters.{name} has an invalid range")
    for name, minimum in (("n_waves", 1), ("random_seed", 0)):
        value = getattr(parameters, name)
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"case_parameters.{name} must be an integer >= {minimum}")
    if not isinstance(parameters.moisture, bool):
        raise ValueError("case_parameters.moisture must be boolean")
    if (
        not isinstance(parameters.noise_modes, list | tuple)
        or any(
            isinstance(mode, bool) or not isinstance(mode, int) or mode <= 0
            for mode in parameters.noise_modes
        )
        or len(set(parameters.noise_modes)) != len(parameters.noise_modes)
    ):
        raise ValueError("noise_modes must be distinct positive integers")
    if parameters.epsilon_noise > 0 and not parameters.noise_modes:
        raise ValueError("nonzero epsilon_noise requires noise_modes")
    if parameters.rh_peak > 1:
        raise ValueError("rh_peak must be at most one")
    if parameters.wall_transition_width > 2.5:
        raise ValueError("wall_transition_width must be <= 2.5 for the monotone degree-five taper")
    if parameters.height / 2 <= parameters.wall_buffer + parameters.wall_transition_width:
        raise ValueError(
            "height must leave a positive unmodified core between the wall transitions"
        )
    return replace(parameters, noise_modes=tuple(parameters.noise_modes))


def vertical_mesh_vertices(parameters: ShearLayerParameters, cells: int) -> np.ndarray:
    """Uniform or sinh-clustered rows, checked for collapse and resolved wall plateaus."""
    half_height = parameters.height / 2
    vertices = np.linspace(-half_height, half_height, cells + 1)
    strength = parameters.clustering
    if strength > 0:
        radius = np.abs(vertices) / half_height
        vertices = (
            half_height
            * np.sign(vertices)
            * np.exp(strength * (radius - 1))
            * (-np.expm1(-2 * strength * radius) / -math.expm1(-2 * strength))
        )
        spacing = np.diff(vertices)
        if not np.isfinite(vertices).all() or not (spacing > 0).all():
            raise ValueError("clustering collapses mesh rows; reduce clustering or change nz")
        if max(spacing[0], spacing[-1]) > parameters.wall_buffer * (1 + 1e-12):
            raise ValueError("clustering needs a full wall cell in each plateau; increase nz")
    return vertices


def mesh_quality(mesh: fd.MeshGeometry, reference_jacobian: float) -> dict[str, float]:
    """MPI-global cell measures, vertical-map Jacobian and longest-edge/altitude ratio."""
    space = fd.FunctionSpace(mesh, "DG", 0)
    quality = {}
    for name, expression in {
        "absolute_jacobian": abs(ufl.JacobianDeterminant(mesh)),
        "map_jacobian": abs(ufl.JacobianDeterminant(mesh)) / reference_jacobian,
        "edge": fd.MinCellEdgeLength(mesh),
        "aspect": fd.MaxCellEdgeLength(mesh) ** 2 / (2 * fd.CellVolume(mesh)),
    }.items():
        field = fd.Function(space).interpolate(expression)
        with field.dat.vec_ro as vector:
            quality[f"minimum_{name}"] = float(vector.min()[1])
            quality[f"maximum_{name}"] = float(vector.max()[1])
    if not all(math.isfinite(value) and value > 0 for value in quality.values()):
        raise ValueError(
            "channel mesh must have finite positive Jacobians, edges and aspect ratios"
        )
    return quality


def wall_profiles(z: Expr, parameters: ShearLayerParameters) -> tuple[Expr, Expr, Expr]:
    """Return base velocity/temperature, its even primitive and a C² compact envelope."""
    core = parameters.height / 2 - parameters.wall_buffer - parameters.wall_transition_width
    flat = parameters.height / 2 - parameters.wall_buffer
    width = parameters.wall_transition_width
    value = math.tanh(core)
    tail = math.exp(-2 * core)
    gradient = 4 * tail / (1 + tail) ** 2
    bernstein = gradient * (1 - 2 * width * value / 5)
    potential = core + math.log1p(tail) - math.log(2)
    radius = abs(z)
    t = (radius - core) / width
    r = 1 - t
    integral = gradient * (1 - r**6) / 6 + bernstein * (1 / 6 - r**5 + 5 * r**6 / 6)
    twice_integrated = gradient / 6 * (t - (1 - r**7) / 7) + bernstein * (
        t / 6 - (1 - r**6) / 6 + 5 * (1 - r**7) / 42
    )
    plateau = value + width * (gradient + bernstein) / 6
    potential_flat = potential + width * value + width**2 * (gradient / 7 + 5 * bernstein / 42)
    sign = ufl.conditional(z >= 0, 1.0, -1.0)
    base = ufl.conditional(
        radius <= core,
        ufl.tanh(z),
        sign * ufl.conditional(radius < flat, value + width * integral, plateau),
    )
    core_potential = ufl.conditional(
        z >= 0,
        z + ufl.ln(1 + ufl.exp(-2 * z)) - math.log(2),
        -z + ufl.ln(1 + ufl.exp(2 * z)) - math.log(2),
    )
    primitive = ufl.conditional(
        radius <= core,
        core_potential,
        ufl.conditional(
            radius < flat,
            potential + width * value * t + width**2 * twice_integrated,
            potential_flat + plateau * (radius - flat),
        ),
    )
    envelope = ufl.conditional(
        radius <= core, 1.0, ufl.conditional(radius < flat, r**3 * (1 + 3 * t + 6 * t * t), 0.0)
    )
    return base, primitive, envelope


def build(config: Config) -> Problem:
    """Build the shear-layer case on checked affine rows with zero force and scalar sources."""
    parameters = checked_parameters(config.case_parameters)
    length = parameters.n_waves * 2 * math.pi / parameters.wave_number
    cells = builtin_cells(config, (128 * parameters.n_waves, 192))
    vertices = vertical_mesh_vertices(parameters, cells[1])
    mesh, tags = periodic_channel(
        cells, (0.0, -parameters.height / 2), (length, parameters.height / 2)
    )
    if parameters.clustering > 0:
        coordinates = mesh.coordinates.dat.data
        rows = np.rint((coordinates[:, 1] / parameters.height + 0.5) * cells[1]).astype(int)
        coordinates[:, 1] = vertices[rows]
    quality = mesh_quality(mesh, length / cells[0] * parameters.height / cells[1])
    x = ufl.SpatialCoordinate(mesh)
    base, primitive, envelope = wall_profiles(x[1], parameters)
    phases = (
        np.random.Generator(np.random.PCG64(parameters.random_seed))
        .uniform(0.0, 2 * math.pi, len(parameters.noise_modes))
        .tolist()
    )
    normalization = 1 / math.sqrt(len(phases)) if phases else 0.0
    noise = sum(
        normalization * ufl.cos(2 * math.pi * mode * x[0] / length + phase)
        for mode, phase in zip(parameters.noise_modes, phases, strict=True)
    )
    perturbation = ufl.zero((2,))
    if parameters.epsilon > 0:
        stream = (
            parameters.epsilon
            * ufl.exp(-(x[1] ** 2))
            * envelope
            * (ufl.cos(parameters.wave_number * x[0]) + parameters.epsilon_noise * noise)
        )
        perturbation = ufl.as_vector((stream.dx(1), -stream.dx(0)))
    initial = ufl.as_vector((base, 0.0)) + perturbation
    moisture = (
        parameters.rh_peak
        * ufl.exp(
            -((x[1] / parameters.moisture_width) ** 2)
            + parameters.saturation_temperature_rate * base
            - parameters.saturation_height_rate * x[1]
        )
        * envelope
    )
    boundaries = {name: DiffusiveFlux() for name in tags}

    def zero(position: Expr, _time: Expr) -> Expr:
        return 0.0 * position[0]

    def forcing(position: Expr, _time: Expr) -> Expr:
        return ufl.as_vector((0.0 * position[0], 0.0 * position[1]))

    metadata = {
        "parameters": asdict(parameters),
        "length_scale": "delta",
        "velocity_scale": "U0",
        "temperature": "theta=(T-T_ref)/DeltaT",
        "domain_length": length,
        "mesh_cells": list(cells),
        "minimum_edge": quality["minimum_edge"],
        "mesh_quality": quality,
        "vertical_vertices": vertices.tolist(),
        "wall_core": parameters.height / 2
        - parameters.wall_buffer
        - parameters.wall_transition_width,
        "wall_flat": parameters.height / 2 - parameters.wall_buffer,
        "profile": "tanh core; monotone degree-five gradient taper; constant wall plateau",
        "hydrostatic_pressure": "Ri times analytic even primitive of the adjusted temperature",
        "streamfunction": "epsilon*exp(-z^2)*envelope*(cos(k*x)+epsilon_noise*noise)",
        "noise_generator": "PCG64",
        "noise_phases": phases,
        "noise_amplitude_per_mode": normalization,
        "noise_normalization": "same spatial RMS as cos(k*x) for distinct integer modes",
        "moisture_wall_adjustment": "multiply core RH/saturation profile by the compact envelope",
        "subscale_initialization": "zero",
        "mesh_clustering": "uniform"
        if parameters.clustering == 0
        else "sinh vertex map; affine cells",
    }
    passive = (
        (
            PassiveScalarProblem(
                "q",
                1 / (parameters.reynolds * parameters.schmidt),
                zero,
                lambda _x: moisture,
                boundaries,
            ),
        )
        if parameters.moisture
        else ()
    )
    return Problem(
        mesh=mesh,
        boundary_tags=tags,
        boundary_conditions={name: FreeSlip(1) for name in tags},
        viscosity=1 / parameters.reynolds,
        forcing=forcing,
        initial_velocity=lambda _x: initial,
        thermal=ThermalProblem(
            1 / (parameters.reynolds * parameters.prandtl),
            parameters.richardson,
            (0.0, -1.0),
            0.0,
            zero,
            lambda _x: base,
            boundaries,
            scaling="nondimensional",
        ),
        passive_scalars=passive,
        initial_pressure=lambda _x: parameters.richardson * primitive,
        reference_values={
            "reynolds": parameters.reynolds,
            "prandtl": parameters.prandtl,
            "schmidt": parameters.schmidt,
            "richardson": parameters.richardson,
            "length_scale": 1.0,
            "velocity_scale": 1.0,
            "temperature_scale": 1.0,
        },
        case_metadata=metadata,
        horizontal_profile_geometry=HorizontalProfileGeometry(length, cells[0], tuple(vertices)),
    )
