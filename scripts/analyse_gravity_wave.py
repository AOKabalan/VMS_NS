"""Measure stable wave period, isentrope displacement and overturning from original FE fields."""

import argparse
import csv
import json
import logging
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from visualize_clouds import read_frame, read_series

from navier_stokes_vms.postprocessing.cloud_visualization import (
    triangle_basis,
    triangle_gradients,
    triangle_rule,
)


def analyse(folder: Path, output: Path) -> dict:
    """Compare the standing mode with linear theory and check vertex gradients for overturning."""
    info = json.loads((folder / "run_info.json").read_text())
    if info["exit_status"] != 0:
        raise ValueError("wave analysis requires a completed run")
    thermal = json.loads((folder / "thermal_metadata.json").read_text())
    if thermal["viscosity"] != thermal["diffusivity"]:
        raise ValueError("the damped linear reference assumes equal momentum/thermal diffusion")
    metadata = thermal["case_metadata"]
    parameters = metadata["parameters"]
    fields = ["velocity", "temperature"] + (["q"] if parameters["moisture"] else [])
    series = {name: read_series(folder / f"fields/{name}.pvd") for name in fields}
    times = list(series["velocity"])
    if any(list(s) != times for s in series.values()):
        raise ValueError("wave field times differ")
    bary, weights = triangle_rule(10)
    k, m = parameters["wave_number"], np.pi / parameters["height"]
    gradient = parameters["temperature_gradient"]
    area = metadata["domain_length"] * parameters["height"]
    records = []
    for time in times:
        _, points, cells, arrays = read_frame({n: s[time] for n, s in series.items()}, area)
        corners = points[cells[:, :3], :2]
        edges = corners[:, 1:] - corners[:, :1]
        cell_area = abs(edges[:, 0, 0] * edges[:, 1, 1] - edges[:, 1, 0] * edges[:, 0, 1]) / 2
        basis, _ = triangle_basis(bary, cells.shape[1])
        x = corners[:, :, 0] @ bary
        z = corners[:, :, 1] @ bary
        mode = np.cos(k * x) * np.cos(m * z)
        w = arrays["velocity"][cells, 1] @ basis
        theta = arrays["temperature"][cells] @ basis
        displacement = -(theta - gradient * z) / gradient
        normalization = np.dot(cell_area, mode**2 @ weights)
        derivatives = triangle_gradients(points, cells, arrays["temperature"][:, None], np.eye(3))
        record = {
            "time": time,
            "vertical_velocity_mode": float(
                np.dot(cell_area, (w * mode) @ weights) / normalization
            ),
            "displacement_mode": float(
                np.dot(cell_area, (displacement * mode) @ weights) / normalization
            ),
            "minimum_temperature_gradient": float(np.min(derivatives[:, :, 0, 1])),
        }
        if "q" in arrays:
            q = arrays["q"][cells] @ basis
            cloud = q - np.exp(-parameters["saturation_height_rate"] * z) > 1e-6
            cloud_area = np.dot(cell_area, cloud @ weights)
            lifted_area = np.dot(cell_area, (cloud & (displacement > 0)) @ weights)
            record["cloud_fraction"] = float(cloud_area / area)
            record["cloud_fraction_in_lifted_air"] = (
                float(lifted_area / cloud_area) if cloud_area > 0 else None
            )
        records.append(record)
    t = np.array(times)
    w_mode = np.array([r["vertical_velocity_mode"] for r in records])
    sign_changes = np.flatnonzero(w_mode[:-1] * w_mode[1:] < 0)
    crossings = [
        float(t[j] - w_mode[j] * (t[j + 1] - t[j]) / (w_mode[j + 1] - w_mode[j]))
        for j in sign_changes
    ]
    periods = np.array(crossings[2:]) - np.array(crossings[:-2])
    if len(periods) < 2:
        raise ValueError("wave period measurement needs at least four signed zero crossings")
    observed = float(np.mean(periods))
    expected = metadata["linear_wave_period"]
    summary = {
        "run": str(folder),
        "run_info": info,
        "linear_reference_period": expected,
        "measured_period": observed,
        "period_relative_error": abs(observed / expected - 1),
        "period_estimates": periods.tolist(),
        "velocity_zero_crossings": crossings,
        "minimum_temperature_gradient": min(r["minimum_temperature_gradient"] for r in records),
        "maximum_displacement_mode": max(abs(r["displacement_mode"]) for r in records),
        "reference_scope": "linear damped mode; full solver retains nonlinear advection",
        "overturning_check": "minimum P1/P2 dT/dz at vertices of every affine parent triangle",
    }
    amplitude = parameters["displacement_amplitude"]
    omega = metadata["linear_wave_frequency"]
    decay = np.exp(-metadata["linear_decay_rate"] * t)
    velocity_reference = amplitude * omega * decay * np.cos(omega * t)
    displacement_reference = amplitude * decay * np.sin(omega * t)
    velocity_error = float(
        np.sqrt(np.mean((w_mode - velocity_reference) ** 2)) / (amplitude * omega)
    )
    displacement_error = float(
        np.sqrt(
            np.mean(
                (np.array([r["displacement_mode"] for r in records]) - displacement_reference) ** 2
            )
        )
        / amplitude
    )
    summary["velocity_mode_rms_error_over_initial_amplitude"] = velocity_error
    summary["displacement_mode_rms_error_over_initial_amplitude"] = displacement_error
    summary["acceptance"] = {
        "period_within_two_percent": summary["period_relative_error"] < 0.02,
        "modal_errors_within_five_percent": max(velocity_error, displacement_error) < 0.05,
        "no_overturning": summary["minimum_temperature_gradient"] > 0,
    }
    if parameters["moisture"]:
        lifted = [
            r["cloud_fraction_in_lifted_air"]
            for r in records
            if r["cloud_fraction_in_lifted_air"] is not None
        ]
        summary["minimum_cloud_area_fraction_in_lifted_air"] = min(lifted) if lifted else None
        summary["cloud_fraction_quadrature_order"] = 10
    output.mkdir(parents=True, exist_ok=True)
    (output / "wave_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (output / "wave_modes.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.8), sharex=True)
    axes[0].plot(t, w_mode, label="Measured velocity mode")
    axes[0].plot(t, amplitude * omega * decay * np.cos(omega * t), "--", label="Linear reference")
    axes[0].set_ylabel("Vertical velocity amplitude")
    axes[1].plot(t, [r["displacement_mode"] for r in records], label="Inferred displacement mode")
    axes[1].plot(t, amplitude * decay * np.sin(omega * t), "--", label="Linear reference")
    axes[1].set_ylabel("Displacement amplitude")
    axes[1].set_xlabel("Dimensionless time")
    for ax in axes:
        ax.legend()
        ax.grid(alpha=0.25)
    fig.suptitle(f"Internal gravity wave · measured period {observed:.4f} / linear {expected:.4f}")
    fig.tight_layout()
    fig.savefig(output / "wave_modes.png", dpi=180)
    fig.savefig(output / "wave_modes.pdf")
    plt.close(fig)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    summary = analyse(args.run, args.output or args.run / "wave_analysis")
    logging.info(
        "Measured period %.6g; relative error %.3g; minimum dT/dz %.6g",
        summary["measured_period"],
        summary["period_relative_error"],
        summary["minimum_temperature_gradient"],
    )


if __name__ == "__main__":
    main()
