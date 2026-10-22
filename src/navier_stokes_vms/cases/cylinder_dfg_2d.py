"""DFG 2D-2: Schäfer & Turek (1996), Table 4; FEATFLOW's Re=100 periodic benchmark."""

import json
from pathlib import Path

import firedrake as fd
import numpy as np
import ufl

from navier_stokes_vms.config import Config
from navier_stokes_vms.problem import mesh_loading
from navier_stokes_vms.problem.boundary_conditions import Dirichlet, Outflow
from navier_stokes_vms.problem.problem_definition import Problem, case_parameters

ROOT = Path(__file__).resolve().parents[3]
REFERENCES = ROOT / "data/references/cylinder_dfg_2d"
DIAMETER = 0.1
MEAN_SPEED = 1.0


def build(config: Config) -> Problem:
    """No-slip cylinder and walls, parabolic inflow, gradient-form do-nothing outlet."""
    case_parameters(config, {})
    if config.mesh.source == "builtin" or config.mesh.cells is not None:
        raise ValueError("cylinder_dfg_2d requires a gmsh mesh, not mesh.cells")
    path = Path(config.mesh.file) if config.mesh.file else ROOT / "meshes/cylinder_dfg_2d.msh"
    mesh, tags = mesh_loading.gmsh_mesh(path)
    if set(tags) != {"inlet", "outlet", "walls", "cylinder"}:
        raise ValueError("cylinder mesh requires exactly inlet, outlet, walls, cylinder")
    # For an incompressible no-slip wall, (∇u)^T n=0: FEATFLOW's gradient stress
    # and our symmetric stress give the same continuum force (discrete errors can differ).
    return Problem(
        mesh=mesh,
        boundary_conditions={
            "inlet": Dirichlet(
                lambda x, _t: ufl.as_vector((6 * x[1] * (0.41 - x[1]) / 0.41**2, 0))
            ),
            "outlet": Outflow(),
            "walls": Dirichlet(lambda _x, _t: ufl.zero(2)),
            "cylinder": Dirichlet(lambda _x, _t: ufl.zero(2)),
        },
        boundary_tags=tags,
        viscosity=1e-3,
        forcing=lambda _x, _t: ufl.zero(2),
        initial_velocity=lambda _x: ufl.zero(2),
        reference_values=json.loads((REFERENCES / "featflow.json").read_text()),
        force_boundaries=("cylinder",),
    )


def lift_minima(times: np.ndarray, lift: np.ndarray) -> np.ndarray:
    """Interior minima, refined by the local quadratic on three samples."""
    indices = np.flatnonzero((lift[1:-1] < lift[:-2]) & (lift[1:-1] <= lift[2:])) + 1
    minima = []
    for i in indices:
        offsets = times[i - 1 : i + 2] - times[i]
        a, b, _ = np.polyfit(offsets, lift[i - 1 : i + 2], 2)
        vertex = -b / (2 * a)
        minima.append(times[i] + np.clip(vertex, offsets[0], offsets[-1]))
    return np.asarray(minima)


def last_cycle(times: np.ndarray, drag: np.ndarray, lift: np.ndarray) -> np.ndarray:
    """Last complete lift-minimum-to-minimum cycle, with elapsed physical time in column 0."""
    minima = lift_minima(times, lift)
    if len(minima) < 2:
        raise ValueError("insufficient lift minima: no complete shedding cycle")
    start, end = minima[-2:]
    t = np.r_[start, times[(times > start) & (times < end)], end]
    return np.column_stack((t - start, np.interp(t, times, drag), np.interp(t, times, lift)))


def compare(times: np.ndarray, drag: np.ndarray, lift: np.ndarray) -> dict[str, float]:
    """FEATFLOW cycle extrema, midrange mean, peak-to-peak amplitude, St and periodicity."""
    cycle = last_cycle(times, drag, lift)
    values = {"st": DIAMETER / (MEAN_SPEED * cycle[-1, 0])}
    minima = lift_minima(times, lift)
    values.update(cycle_start=float(minima[-2]), cycle_end=float(minima[-1]))
    for column, name in ((1, "cd"), (2, "cl")):
        low, high = float(np.min(cycle[:, column])), float(np.max(cycle[:, column]))
        values.update(
            {
                f"{name}_min": low,
                f"{name}_max": high,
                f"{name}_mean": (low + high) / 2,
                f"{name}_amplitude": high - low,
            }
        )
    if len(minima) >= 4:
        peaks = [
            float(np.max(lift[(times >= start) & (times <= end)]))
            for start, end in zip(minima[-4:-1], minima[-3:], strict=True)
        ]
        values["cl_max_change_last_three_cycles"] = max(
            abs(b - a) for a, b in zip(peaks[:-1], peaks[1:], strict=True)
        )
    return values


def reference_curve() -> np.ndarray:
    """Newest high-accuracy level-6 data, columns elapsed time, C_D, C_L."""
    data = np.loadtxt(REFERENCES / "featflow_curve.csv", delimiter=",", skiprows=1)
    return last_cycle(data[:, 0], data[:, 1], data[:, 2])


def compare_curves(cycle: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    """Max errors on the common elapsed-time interval, aligned at C_L minima, no time scaling."""
    end = min(cycle[-1, 0], reference[-1, 0])
    times = np.unique(
        np.r_[cycle[cycle[:, 0] <= end, 0], reference[reference[:, 0] <= end, 0], end]
    )
    return {
        f"{name}_curve_max_absolute_deviation": float(
            np.max(
                np.abs(
                    np.interp(times, cycle[:, 0], cycle[:, column])
                    - np.interp(times, reference[:, 0], reference[:, column])
                )
            )
        )
        for column, name in ((1, "cd"), (2, "cl"))
    }


def outlet_evaluator(problem: Problem) -> fd.PointEvaluator:
    """Collective evaluation at the outlet nodes and edge midpoints, in boundary order."""
    return fd.PointEvaluator(problem.mesh, outlet_points(problem))


def outlet_points(problem: Problem) -> np.ndarray:
    """Global endpoints and midpoints of the straight outlet edges, in boundary order."""
    space = fd.FunctionSpace(problem.mesh, "CG", 1)
    y = fd.Function(space).interpolate(fd.SpatialCoordinate(problem.mesh)[1])
    nodes = fd.DirichletBC(space, 0, problem.boundary_tags["outlet"]).nodes
    coordinates = problem.mesh.comm.allgather(y.dat.data_ro_with_halos[nodes].tolist())
    stations = np.unique(np.concatenate(coordinates))
    samples = np.sort(np.r_[stations, (stations[:-1] + stations[1:]) / 2])
    return np.column_stack((np.full(len(samples), 2.2), samples))


def outlet_minimum(velocity: fd.Function, evaluator: fd.PointEvaluator) -> float:
    """Global minimum of the P1/P2 normal trace, including each edge's interior extremum."""
    values = evaluator.evaluate(velocity)[:, 0]
    left, middle, right = values[:-2:2], values[1::2], values[2::2]
    a = 2 * (left + right - 2 * middle)
    b = right - left - a
    candidates = list(values)
    for quadratic, linear, constant in zip(a, b, left, strict=True):
        if quadratic > 0:
            vertex = -linear / (2 * quadratic)
            if 0 < vertex < 1:
                candidates.append(constant - linear**2 / (4 * quadratic))
    return float(min(candidates))
