"""Measure a KH Fourier mode in completed matched runs, with an explicit fit window."""

import argparse
import json
from pathlib import Path

import numpy as np
from analyse_kh_startup import compare as compare_startup
from analyse_kh_startup import read_rows


def vertical_integral(values: np.ndarray, heights: np.ndarray) -> np.ndarray:
    """Composite Simpson integration over each actual mesh row and its midpoint."""
    if (
        len(heights) < 3
        or len(heights) % 2 != 1
        or not np.all(np.diff(heights) > 0)
        or not np.allclose(heights[1::2], (heights[:-2:2] + heights[2::2]) / 2, atol=1e-12, rtol=0)
    ):
        raise ValueError("Fourier energy requires ordered mesh rows and midrows")
    return np.sum(
        np.diff(heights[::2])
        * (values[..., :-2:2] + 4 * values[..., 1::2] + values[..., 2::2])
        / 6,
        axis=-1,
    )


def mode_profiles(folder: Path, mode: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Read real sine/cosine coefficients of vertical velocity, including initialization."""
    metadata = json.loads((folder / "horizontal_profiles.json").read_text())
    if mode not in metadata.get("fourier_modes", ()):
        raise ValueError(f"{folder}: mode {mode} was not recorded")
    rows = read_rows(folder / "horizontal_profiles.csv")
    steps = sorted({int(row["step"]) for row in rows})
    heights = np.asarray(metadata["heights"])
    coefficients, times = [], []
    for step in steps:
        records = [row for row in rows if int(row["step"]) == step]
        if not np.array_equal([float(row["z"]) for row in records], heights):
            raise ValueError(f"{folder}: inconsistent Fourier profile heights")
        coefficients.append(
            [
                [float(row[f"velocity_y_mode_{mode}_{part}"]) for row in records]
                for part in ("cos", "sin")
            ]
        )
        times.append(float(records[0]["time"]))
    coefficients = np.asarray(coefficients)
    if not np.isfinite(coefficients).all():
        raise ValueError(f"{folder}: nonfinite Fourier coefficients")
    return np.asarray(times), heights, coefficients, metadata["length"]


def fit_energy(times: np.ndarray, energy: np.ndarray, window: tuple[float, float]) -> dict:
    """Fit log energy = intercept + 2 sigma t; retain curvature and OLS residual uncertainty."""
    start, end = window
    if not np.isfinite(window).all() or start < 0 or end <= start:
        raise ValueError("growth window must have finite 0 <= start < end")
    if start < times[0] - 1e-12 or end > times[-1] + 1e-12:
        raise ValueError("growth window lies outside recorded times")
    mask = (times >= start - 1e-12) & (times <= end + 1e-12)
    t, values = times[mask], energy[mask]
    if len(t) < 6 or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("growth fit requires at least six positive finite energy samples")
    centered = t - np.mean(t)
    logarithm = np.log(values)
    slope = np.dot(centered, logarithm) / np.dot(centered, centered)
    residual = logarithm - np.mean(logarithm) - slope * centered
    half = len(t) // 2
    half_slopes = [
        float(np.polyfit(t[s], logarithm[s], 1)[0] / 2)
        for s in (slice(None, half), slice(half, None))
    ]
    return {
        "window": [float(t[0]), float(t[-1])],
        "samples": len(t),
        "sigma": float(slope / 2),
        "rms_log_energy_residual": float(np.sqrt(np.mean(residual**2))),
        "maximum_log_energy_residual": float(np.max(abs(residual))),
        "ols_sigma_standard_error": float(
            np.sqrt(np.dot(residual, residual) / (len(t) - 2) / np.dot(centered, centered)) / 2
        ),
        "half_window_sigmas": half_slopes,
        "half_window_sigma_difference": abs(half_slopes[1] - half_slopes[0]),
    }


def compare(folders: list[Path], mode: int, window: tuple[float, float] | None) -> dict:
    """Retain raw and mode energies; subtract matched coefficients before squaring."""
    startup = compare_startup(folders)
    pairs = []
    for pair in startup["pairs"]:
        control, seeded = pair["control"], pair["seeded"]
        times, heights, c, length = mode_profiles(Path(control["folder"]), mode)
        t_seed, z_seed, s, l_seed = mode_profiles(Path(seeded["folder"]), mode)
        if not (
            np.array_equal(times, t_seed) and np.array_equal(heights, z_seed) and length == l_seed
        ):
            raise ValueError("matched Fourier profiles require identical times and geometry")
        energies = {
            name: length * vertical_integral(np.sum(value**2, axis=1), heights)
            for name, value in (("control_mode", c), ("seeded_mode", s), ("difference_mode", s - c))
        }
        if np.any(energies["difference_mode"] <= 0):
            raise ValueError("seeded-minus-control Fourier energy must be positive")
        raw = np.r_[seeded["initial"]["vertical_kinetic_energy"], seeded["vertical_energy"]]
        raw_times = np.r_[0, seeded["times"]]
        difference = s - c
        entry = {
            "control": control,
            "seeded": seeded,
            "times": times.tolist(),
            "energies": {name: values.tolist() for name, values in energies.items()},
            "maximum_control_over_difference_mode": float(
                np.max(energies["control_mode"] / energies["difference_mode"])
            ),
            "maximum_difference_mode_velocity": float(
                np.max(2 * np.sqrt(np.sum(difference**2, axis=1)))
            ),
            "fits": {},
        }
        profile_metadata = json.loads(
            (Path(seeded["folder"]) / "horizontal_profiles.json").read_text()
        )
        if 2 * mode in profile_metadata["fourier_modes"]:
            _, _, harmonic_c, _ = mode_profiles(Path(control["folder"]), 2 * mode)
            _, _, harmonic_s, _ = mode_profiles(Path(seeded["folder"]), 2 * mode)
            harmonic_energy = length * vertical_integral(
                np.sum((harmonic_s - harmonic_c) ** 2, axis=1), heights
            )
            entry["second_harmonic_over_difference_mode"] = (
                harmonic_energy / energies["difference_mode"]
            ).tolist()
        if window is not None:
            entry["fits"] = {
                name: fit_energy(times, values, window)
                for name, values in energies.items()
                if name != "control_mode"
            }
            entry["fits"]["raw_vertical"] = fit_energy(raw_times, raw, window)
            mask = (times >= window[0] - 1e-12) & (times <= window[1] + 1e-12)
            interpolated_raw = np.interp(times[mask], raw_times, raw)
            entry["window_mode_fraction_of_raw_energy"] = (
                energies["seeded_mode"][mask] / interpolated_raw
            ).tolist()
            v = difference[mask, 0] - 1j * difference[mask, 1]
            overlap = vertical_integral(np.conj(v[0]) * v[-1], heights)
            norm_product = np.sqrt(
                vertical_integral(abs(v[0]) ** 2, heights)
                * vertical_integral(abs(v[-1]) ** 2, heights)
            )
            entry["window_endpoint_mode_shape_coherence"] = float(abs(overlap) / norm_product)
            entry["window_maximum_control_over_difference_mode"] = float(
                np.max(energies["control_mode"][mask] / energies["difference_mode"][mask])
            )
        pairs.append(entry)
    return {
        "phase": "B7.3",
        "mode": mode,
        "pairs": pairs,
        "definitions": {
            "coefficients": "C=<w cos>, S=<w sin>; horizontal average over L",
            "energy": "L integral_z (C^2+S^2); vertical composite Simpson on rows/midrows",
            "difference": "subtract matched real coefficients before computing mode energy",
            "fit": "sigma=0.5 slope(log energy); window is supplied explicitly",
            "uncertainty": "OLS error excludes discretization error and base drift",
            "acceptance": "fits alone do not certify an eigenmode, resolution or linear regime",
        },
    }


def plot_comparison(report: dict, path: Path) -> None:
    """Plot raw/mode energies, local logarithmic slopes, contamination and base drift."""
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for pair in report["pairs"]:
        seeded = pair["seeded"]
        label = (
            f"Ri={seeded['parameters']['richardson']:g}, "
            f"Re={seeded['parameters']['reynolds']:g}, "
            f"{seeded['cells'][0]}×{seeded['cells'][1]}, dt={seeded['dt']:g}"
        )
        times = np.asarray(pair["times"])
        energies = pair["energies"]
        difference = np.asarray(energies["difference_mode"])
        (line,) = axes[0, 0].semilogy(times, difference, label=label)
        color = line.get_color()
        axes[0, 0].semilogy(seeded["times"], seeded["vertical_energy"], "--", color=color)
        axes[0, 1].plot(times, 0.5 * np.gradient(np.log(difference), times), color=color)
        axes[1, 0].plot(times, np.asarray(energies["control_mode"]) / difference, color=color)
        axes[1, 1].plot(
            seeded["profile_times"],
            seeded["profile_maximum_absolute_drift"]["velocity_x"],
            color=color,
        )
        fit = pair["fits"].get("difference_mode")
        if fit is not None:
            for axis in axes[0]:
                axis.axvspan(*fit["window"], color=color, alpha=0.08)
    axes[0, 0].set(
        xlabel="Time",
        ylabel="Vertical kinetic energy",
        title="Difference mode (solid); total seeded (dashed)",
    )
    axes[0, 0].legend(fontsize=7)
    axes[0, 1].set(
        xlabel="Time",
        ylabel="Half local slope of log mode energy",
        title="Local slope diagnostic; shaded explicit fit window",
    )
    axes[1, 0].set(
        xlabel="Time",
        ylabel="Control mode / difference mode energy",
        title="Contamination at the imposed wavelength",
    )
    axes[1, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[1, 1].set(
        xlabel="Time", ylabel="Maximum sampled mean-velocity change", title="Base-profile evolution"
    )
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    """Analyse completed runs; all simulation integration uses run_simulation.py."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", nargs="+", type=Path)
    parser.add_argument("--mode", type=int, default=1)
    parser.add_argument("--window", nargs=2, type=float, metavar=("START", "END"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figure", type=Path)
    arguments = parser.parse_args(argv)
    if arguments.mode < 1:
        parser.error("--mode must be positive")
    report = compare(arguments.folders, arguments.mode, arguments.window)
    if arguments.figure is not None:
        plot_comparison(report, arguments.figure)
    arguments.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
