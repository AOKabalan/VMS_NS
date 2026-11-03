"""Collective horizontal means of CG fields in a structured affine triangle channel."""

import json
import math
from collections.abc import Mapping
from pathlib import Path

import firedrake as fd
import numpy as np

from navier_stokes_vms.postprocessing.output import CsvWriter
from navier_stokes_vms.problem.problem_definition import HorizontalProfileGeometry


def profile_quadrature(
    geometry: HorizontalProfileGeometry, degree: int = 2
) -> tuple[np.ndarray, np.ndarray]:
    """Gauss integration per horizontal half-cell, on every row and midrow."""
    vertices = np.asarray(geometry.vertical_vertices)
    if (
        not math.isfinite(geometry.length)
        or geometry.length <= 0
        or isinstance(geometry.horizontal_cells, bool)
        or not isinstance(geometry.horizontal_cells, int)
        or geometry.horizontal_cells <= 0
        or vertices.ndim != 1
        or len(vertices) < 2
        or not np.isfinite(vertices).all()
        or not (np.diff(vertices) > 0).all()
    ):
        raise ValueError(
            "horizontal profiles require positive length/cells and ordered finite rows"
        )
    heights = np.sort(np.r_[vertices, (vertices[:-1] + vertices[1:]) / 2])
    half_cells = 2 * geometry.horizontal_cells
    nodes, _weights = np.polynomial.legendre.leggauss(max(2, (degree + 2) // 2))
    offsets = (nodes + 1) / 2
    horizontal = ((np.arange(half_cells)[:, None] + offsets) * geometry.length / half_cells).ravel()
    points = np.column_stack(
        (np.tile(horizontal, len(heights)), np.repeat(heights, len(horizontal)))
    )
    return heights, points


class HorizontalProfiles:
    """Reusable collective point evaluator and named endpoint profile CSV, including t=0."""

    def __init__(
        self,
        folder: Path,
        geometry: HorizontalProfileGeometry,
        velocity: fd.Function,
        scalars: Mapping[str, fd.Function],
        modes: tuple[int, ...] = (),
    ) -> None:
        mesh = velocity.function_space().mesh()
        if mesh.geometric_dimension != 2:
            raise ValueError("horizontal profiles require a two-dimensional channel")
        if {"step", "time", "z", "velocity_x", "velocity_y"}.intersection(scalars):
            raise ValueError("horizontal profile scalar names conflict with reserved CSV columns")
        fields = (velocity, *scalars.values())
        degree = max(field.ufl_element().degree() for field in fields)
        # Fourier weights are not polynomials; overintegrate their product with FE traces.
        integration_degree = max(degree, 14) if modes else degree
        self._heights, points = profile_quadrature(geometry, integration_degree)
        _nodes, weights = np.polynomial.legendre.leggauss(max(2, (integration_degree + 2) // 2))
        self._weights = np.tile(weights / 2, 2 * geometry.horizontal_cells) / (
            2 * geometry.horizontal_cells
        )
        self._evaluator = fd.PointEvaluator(mesh, points)
        self._velocity = velocity
        self._scalars = dict(scalars)
        self._names = ("velocity_x", "velocity_y", *scalars)
        self._fourier_weights = {
            f"mode_{mode}_{part}": self._weights
            * function(2 * np.pi * mode * points[: len(self._weights), 0] / geometry.length)
            for mode in modes
            for part, function in (("cos", np.cos), ("sin", np.sin))
        }
        names = (
            *self._names,
            *(f"{name}_{suffix}" for suffix in self._fourier_weights for name in self._names),
        )
        if len(names) != len(set(names)):
            raise ValueError("horizontal profile scalar names conflict with Fourier columns")
        self._writer = CsvWriter(
            mesh.comm, folder / "horizontal_profiles.csv", ("step", "time", "z", *names)
        )
        if mesh.comm.rank == 0:
            (folder / "horizontal_profiles.json").write_text(
                json.dumps(
                    {
                        "heights": self._heights.tolist(),
                        "horizontal_cells": geometry.horizontal_cells,
                        "length": geometry.length,
                        "quadrature": "Gauss points per horizontal half-cell; divide by length",
                        "gauss_points_per_half_cell": len(weights),
                        "exactness": "polynomial traces at structured mesh rows and midrows",
                        "time_level": "endpoint; includes initialization",
                        "fields": self._names,
                        **(
                            {
                                "fourier_modes": modes,
                                "fourier_normalization": (
                                    "C=(1/L) integral f cos; S=(1/L) integral f sin; "
                                    "f_mode=2(C cos+S sin)"
                                ),
                                "fourier_quadrature": (
                                    "8 or more Gauss points per half-cell; nonpolynomial weights"
                                ),
                            }
                            if modes
                            else {}
                        ),
                    },
                    indent=2,
                )
                + "\n"
            )

    def write(self, step: int, time: float) -> None:
        """Integrate line traces collectively, then write one row per height on rank zero."""
        vector = self._evaluator.evaluate(self._velocity).reshape(len(self._heights), -1, 2)
        samples = {
            "velocity_x": vector[:, :, 0],
            "velocity_y": vector[:, :, 1],
        }
        for name, field in self._scalars.items():
            samples[name] = self._evaluator.evaluate(field).reshape(len(self._heights), -1)
        means = {name: values @ self._weights for name, values in samples.items()}
        means.update(
            {
                f"{name}_{suffix}": values @ weights
                for suffix, weights in self._fourier_weights.items()
                for name, values in samples.items()
            }
        )
        if not all(np.isfinite(values).all() for values in means.values()):
            raise ValueError("horizontal profile sampling returned nonfinite values")
        for index, height in enumerate(self._heights):
            self._writer.write(
                step=step,
                time=time,
                z=height,
                **{name: values[index] for name, values in means.items()},
            )

    def close(self) -> None:
        """Close the profile CSV."""
        self._writer.close()
