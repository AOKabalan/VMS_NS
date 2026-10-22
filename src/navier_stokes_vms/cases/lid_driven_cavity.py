"""Unit cavity: Ghia, Ghia & Shin (1982), JCP 48:387–411, Tables I–II."""

from pathlib import Path

import firedrake as fd
import numpy as np
import ufl

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.boundary_conditions import Dirichlet
from navier_stokes_vms.problem.problem_definition import Problem, case_parameters

REFERENCES = Path(__file__).resolve().parents[3] / "data/references/lid_driven_cavity"
STEADY_TOLERANCE = 1e-8


def reference_profiles(reynolds: float) -> dict[str, np.ndarray]:
    """The paper's 17 rounded stations and velocities, without interpolation or smoothing."""
    if reynolds not in (100, 1000):
        raise ValueError("lid_driven_cavity: reference Reynolds numbers are 100 and 1000")
    column = 1 if reynolds == 100 else 2
    return {
        name: np.loadtxt(REFERENCES / f"{name}.csv", delimiter=",", skiprows=1)[:, [0, column]]
        for name in ("u", "v")
    }


def build(config: Config) -> Problem:
    """Unit-speed lid, with both top corner nodes assigned to the stationary walls."""
    reynolds = case_parameters(config, {"reynolds": 100.0})["reynolds"]
    profiles = reference_profiles(reynolds)
    mesh, tags = mesh_loading.rectangle(
        mesh_loading.builtin_cells(config, (64, 64)), (0.0, 0.0), (1.0, 1.0)
    )
    conditions = {name: Dirichlet(lambda _x, _t: ufl.zero(2)) for name in tags}
    # Matching corner values on both incident boundaries avoids dependence on BC ordering.
    conditions["y_max"] = Dirichlet(
        lambda x, _t: ufl.as_vector((ufl.conditional(ufl.And(x[0] > 0, x[0] < 1), 1.0, 0.0), 0.0))
    )
    return Problem(
        mesh=mesh,
        boundary_conditions=conditions,
        boundary_tags=tags,
        viscosity=1 / reynolds,
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=lambda _x: ufl.zero(2),
        reference_values={
            f"{name}_{coordinate:g}": value
            for name, profile in profiles.items()
            for coordinate, value in profile
        },
    )


def sample_profiles(velocity: fd.Function, reynolds: float) -> dict[str, np.ndarray]:
    """Evaluate the two centrelines collectively; every rank receives all values.

    PointEvaluator gives each point a single owning cell. Function.at evaluates a point on a
    partition boundary on every rank that touches it and rejects round-off differences.
    """
    mesh = velocity.function_space().mesh()
    return {
        name: fd.PointEvaluator(
            mesh, [(0.5, s) if name == "u" else (s, 0.5) for s in profile[:, 0]]
        ).evaluate(velocity)[:, component]
        for component, (name, profile) in enumerate(reference_profiles(reynolds).items())
    }


def compare(samples: dict[str, np.ndarray], reynolds: float) -> dict[str, float]:
    """Maximum absolute deviations from Ghia over all 17 stations on each centreline."""
    return {
        f"{name}_max_absolute_deviation": float(np.max(np.abs(samples[name] - profile[:, 1])))
        for name, profile in reference_profiles(reynolds).items()
    }
