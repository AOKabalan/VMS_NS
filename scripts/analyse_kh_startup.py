"""Compare completed zero-seed/seeded KH startup pairs; never run or fit a simulation."""

import argparse
import csv
import json
import tomllib
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a completed CSV."""
    with path.open() as stream:
        return list(csv.DictReader(stream))


def read_run(folder: Path) -> tuple[str, dict[str, object]]:
    """Validate one completed run and extract raw energy and sampled mean-profile drift."""
    info = json.loads((folder / "run_info.json").read_text())
    if info["exit_status"] != 0:
        raise ValueError(f"{folder}: incomplete or failed run")
    config = tomllib.loads((folder / "config_resolved.toml").read_text())
    if config["case"] != "kelvin_helmholtz_shear_layer":
        raise ValueError(f"{folder}: expected kelvin_helmholtz_shear_layer")
    metadata = json.loads((folder / "thermal_metadata.json").read_text())["case_metadata"]
    initial = json.loads((folder / "initial_state.json").read_text())
    diagnostics = read_rows(folder / "diagnostics.csv")
    profiles = read_rows(folder / "horizontal_profiles.csv")
    expected_steps = round(config["time"]["t_end"] / config["time"]["dt"])
    if (
        len(diagnostics) != expected_steps
        or int(diagnostics[-1]["step"]) != expected_steps
        or not np.isclose(
            float(diagnostics[-1]["time"]), config["time"]["t_end"], atol=1e-13, rtol=0
        )
    ):
        raise ValueError(f"{folder}: expected all configured timesteps")
    times = np.array([float(row["time"]) for row in diagnostics])
    if [int(row["step"]) for row in diagnostics] != list(
        range(1, expected_steps + 1)
    ) or not np.allclose(
        times, np.arange(1, expected_steps + 1) * config["time"]["dt"], atol=1e-13, rtol=0
    ):
        raise ValueError(f"{folder}: inconsistent diagnostic steps/times")
    energy = np.array([float(row["vertical_kinetic_energy"]) for row in diagnostics])
    heights = np.array([float(row["z"]) for row in profiles if int(row["step"]) == 0])
    if len(heights) < 2 or not np.all(np.diff(heights) > 0):
        raise ValueError(f"{folder}: missing or unordered initial profiles")
    profile_steps = sorted({int(row["step"]) for row in profiles})
    if profile_steps[0] != 0 or profile_steps[-1] != expected_steps:
        raise ValueError(f"{folder}: profiles must include initial and final states")
    profile_fields = tuple(name for name in profiles[0] if name not in ("step", "time", "z"))
    profile_values = {}
    profile_times = []
    for step in profile_steps:
        rows = [row for row in profiles if int(row["step"]) == step]
        if len(rows) != len(heights) or not np.array_equal(
            [float(row["z"]) for row in rows], heights
        ):
            raise ValueError(f"{folder}: inconsistent profile heights")
        time = float(rows[0]["time"])
        if any(float(row["time"]) != time for row in rows) or not np.isclose(
            time, step * config["time"]["dt"], atol=1e-13, rtol=0
        ):
            raise ValueError(f"{folder}: inconsistent profile times")
        profile_times.append(time)
        profile_values[step] = {
            name: np.array([float(row[name]) for row in rows]) for name in profile_fields
        }
    if not (
        np.isfinite(energy).all()
        and (energy >= 0).all()
        and all(
            np.isfinite(values).all()
            for fields in profile_values.values()
            for values in fields.values()
        )
    ):
        raise ValueError(f"{folder}: nonfinite profile/energy data")
    drifts = {
        name: [
            float(np.max(abs(profile_values[step][name] - profile_values[0][name])))
            for step in profile_steps
        ]
        for name in profile_fields
    }
    normalized_config = dict(config)
    normalized_config.pop("run_name")
    normalized_config.pop("output")
    parameters = dict(metadata["parameters"])
    epsilon = parameters.pop("epsilon")
    normalized_config["case_parameters"] = parameters
    key = json.dumps(normalized_config, sort_keys=True)
    return key, {
        "folder": str(folder),
        "run_name": config["run_name"],
        "cells": metadata["mesh_cells"],
        "degree": config["discretization"]["velocity_degree"],
        "dt": config["time"]["dt"],
        "t_end": config["time"]["t_end"],
        "epsilon": epsilon,
        "parameters": metadata["parameters"],
        "stabilization": config["stabilization"],
        "scheme": config["time"]["scheme"],
        "mesh_quality": metadata["mesh_quality"],
        "initial": initial,
        "times": times.tolist(),
        "vertical_energy": energy.tolist(),
        "maximum_cfl": max(float(row["advective_cfl"]) for row in diagnostics),
        "profile_times": profile_times,
        "profile_maximum_absolute_drift": drifts,
        "profile_heights": heights.tolist(),
        "final_minus_initial_profiles": {
            name: (profile_values[expected_steps][name] - profile_values[0][name]).tolist()
            for name in profile_fields
        },
        "maximum_absolute_energy_budget_residual": max(
            abs(float(row["boussinesq_energy_budget_residual"])) for row in diagnostics
        ),
        "wall_seconds": info["wall_time_seconds"],
        "provenance": {
            name: info.get(name)
            for name in (
                "git_commit",
                "git_dirty",
                "firedrake_version",
                "petsc_version",
                "python_version",
            )
        },
    }


def compare(folders: list[Path]) -> dict[str, object]:
    """Match numerically identical pairs differing only in epsilon; retain both raw energies."""
    groups = defaultdict(list)
    for folder in folders:
        key, run = read_run(folder)
        groups[key].append(run)
    pairs = []
    for runs in groups.values():
        controls = [run for run in runs if run["epsilon"] == 0]
        seeds = [run for run in runs if run["epsilon"] > 0]
        if len(controls) != 1 or len(seeds) != 1:
            raise ValueError("each setup requires exactly one zero-epsilon and one seeded run")
        control, seeded = controls[0], seeds[0]
        if control["times"] != seeded["times"]:
            raise ValueError("paired diagnostic times differ")
        seed_initial = seeded["initial"]["vertical_kinetic_energy"]
        if seed_initial <= 0:
            raise ValueError("seeded run needs positive initial vertical energy")
        if seeded["vertical_energy"][-1] <= 0:
            raise ValueError("seeded run needs positive final vertical energy for the ratio")
        pairs.append(
            {
                "control": control,
                "seeded": seeded,
                "control_final_over_seeded_final": control["vertical_energy"][-1]
                / seeded["vertical_energy"][-1],
                "maximum_control_over_initial_seed_energy": max(control["vertical_energy"])
                / seed_initial,
            }
        )
    return {
        "phase": "B7.2",
        "pairs": pairs,
        "interpretation": {
            "energies": "raw vertical energies; their difference is not perturbation energy",
            "profile_drift": "max at sampled heights of each mean minus its own initial mean",
            "limits": "startup/mesh checks only; no growth fit or resolution certificate",
        },
    }


def plot_comparison(report: dict[str, object], path: Path) -> None:
    """Export the startup energies and mean-profile drift without fitting growth rates."""
    import matplotlib.pyplot as plt

    pairs = report["pairs"]
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for pair in pairs:
        control, seeded = pair["control"], pair["seeded"]
        label = (
            f"P{control['degree']} {control['cells'][0]}×{control['cells'][1]}, "
            f"λ={control['parameters']['clustering']:g}, Δt={control['dt']:g}"
        )
        (line,) = axes[0, 0].semilogy(control["times"], control["vertical_energy"], label=label)
        axes[0, 0].semilogy(
            seeded["times"], seeded["vertical_energy"], "--", color=line.get_color()
        )
        heights = control["profile_heights"]
        axes[0, 1].plot(
            control["final_minus_initial_profiles"]["velocity_x"], heights, color=line.get_color()
        )
        axes[1, 0].plot(
            control["final_minus_initial_profiles"]["temperature"], heights, color=line.get_color()
        )
        axes[1, 1].plot(
            control["times"],
            np.array(control["vertical_energy"]) / np.array(seeded["vertical_energy"]),
            color=line.get_color(),
        )
    axes[0, 0].set(
        xlabel="Time",
        ylabel="Resolved vertical kinetic energy",
        title="Control (solid) and seeded (dashed)",
    )
    axes[0, 0].legend(fontsize=7)
    axes[0, 1].set(
        xlabel="Change in horizontal-mean velocity",
        ylabel="z",
        title="Zero-seed control: final minus initial",
    )
    axes[1, 0].set(
        xlabel="Change in horizontal-mean temperature",
        ylabel="z",
        title="Zero-seed control: final minus initial",
    )
    axes[1, 1].set(
        xlabel="Time",
        ylabel="Control energy / seeded energy",
        title="Contamination diagnostic; not a growth rate",
    )
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    """Analyse existing outputs only; all integrations use scripts/run_simulation.py."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figure", type=Path)
    arguments = parser.parse_args(argv)
    report = compare(arguments.folders)
    if arguments.figure is not None:
        plot_comparison(report, arguments.figure)
    arguments.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
