"""Plot recorded Boussinesq refinement errors without rerunning the solver."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = ROOT / "docs/verification/boussinesq_mms_recheck_20261103.json"
DEFAULT_OUTPUT = ROOT / "results/bousinesq/verification/manufactured_solution/2026-11-03"
COLORS = {"velocity": "#1766a1", "pressure": "#7351a6", "temperature": "#cf6216"}
SCHEMES = (
    ("bdf1", "BDF1", 1),
    ("bdf2", "BDF2", 2),
    ("bdf3", "BDF3", 3),
    ("crank_nicolson", "Crank–Nicolson", 2),
)


def refinement_axis(ax: Axes, steps: np.ndarray, label: str) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(max(steps) * 1.15, min(steps) / 1.15)
    ax.xaxis.set_major_locator(FixedLocator(steps))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(label)
    ax.grid(which="major", alpha=0.23)
    ax.spines[["top", "right"]].set_visible(False)


def save_figure(fig: Figure, output: Path, name: str) -> None:
    fig.savefig(output / f"{name}.png", dpi=180, facecolor="white")
    fig.savefig(output / f"{name}.pdf", facecolor="white")
    plt.close(fig)


def plot_spatial(summary: dict, output: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
    for row, (degree, result) in enumerate(summary["spatial"].items()):
        k = int(degree[1:])
        h = np.array([1 / float(record["refinement"]) for record in result["table"]])
        for ax, field in zip(axes[row], ("velocity", "pressure", "temperature"), strict=True):
            name = f"{field}_error_l2"
            errors = np.array([float(record[name]) for record in result["table"]])
            order = k if field == "pressure" else k + 1
            rate = result["finest_rates"][name]
            ax.plot(h, errors, "o-", color=COLORS[field], label=f"Measured; rate {rate:.3f}")
            ax.plot(
                h,
                1.8 * errors[0] * (h / h[0]) ** order,
                "--",
                color="0.5",
                label=rf"Target slope $h^{order}$",
            )
            refinement_axis(ax, h, r"Mesh parameter $h=1/N$")
            ax.set_title(f"{degree} · {field}")
            ax.set_ylabel(r"$L^2$ error")
            ax.legend(fontsize=8, loc="best")
    fig.suptitle("Boussinesq manufactured solution: spatial convergence", fontsize=15)
    fig.text(
        0.5,
        0.932,
        "Static linear OSS · 8², 16², 32² meshes · BDF1 · t = 0.1",
        ha="center",
        fontsize=10,
    )
    fig.text(
        0.5,
        0.022,
        "Rates use the finest mesh pair. Exact T is bilinear and belongs to P2; "
        "its observed rate ≈ 4 is case specific.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.047, 1, 0.915), h_pad=2)
    save_figure(fig, output, "spatial_convergence")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
    for ax, field in zip(axes, ("velocity", "temperature"), strict=True):
        for degree, marker, color in (("P1", "o", "#1766a1"), ("P2", "s", "#cf6216")):
            result = summary["spatial"][degree]
            h = np.array([1 / float(record["refinement"]) for record in result["table"]])
            name = f"{field}_gradient_error_l2"
            errors = np.array([float(record[name]) for record in result["table"]])
            rate = result["finest_rates"][name]
            ax.plot(h, errors, marker=marker, color=color, label=f"{degree}; rate {rate:.3f}")
            k = int(degree[1:])
            ax.plot(
                h,
                1.8 * errors[0] * (h / h[0]) ** k,
                "--",
                color=color,
                alpha=0.6,
                label=rf"Target $h^{k}$",
            )
        refinement_axis(ax, h, r"Mesh parameter $h=1/N$")
        ax.set_title(field.capitalize())
        ax.set_ylabel(r"$H^1$ seminorm error")
        ax.legend(fontsize=9)
    fig.suptitle("Boussinesq manufactured solution: gradient convergence", fontsize=14)
    fig.text(
        0.5,
        0.025,
        "Finest-pair rates · bilinear exact T is representable in P2",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.95))
    save_figure(fig, output, "spatial_gradient_convergence")


def plot_temporal(summary: dict, output: Path) -> None:
    dt = np.array(summary["temporal_steps"])
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    for ax, (scheme, title, order) in zip(axes.flat, SCHEMES, strict=True):
        result = summary["temporal"][scheme]
        for field, marker in (("velocity", "o"), ("temperature", "s")):
            errors = np.array(result["errors"][field])
            rate = result["rates"][field]
            ax.plot(
                dt,
                errors,
                marker=marker,
                color=COLORS[field],
                label=f"{field.capitalize()}; rate {rate:.3f}",
            )
        anchor = 1.8 * result["errors"]["temperature"][0]
        ax.plot(
            dt,
            anchor * (dt / dt[0]) ** order,
            "--",
            color="0.5",
            label=rf"Expected slope $\Delta t^{order}$",
        )
        refinement_axis(ax, dt, r"Time step $\Delta t$")
        ax.set_title(title)
        ax.set_ylabel(r"$L^2$ error at $t=0.2$")
        ax.legend(fontsize=8)
    fig.suptitle("Buoyancy-driven shear: temporal convergence", fontsize=15)
    fig.text(
        0.5,
        0.932,
        "Analytic semidiscrete eigenmode · dynamic linear OSS · no spatial-error floor",
        ha="center",
        fontsize=10,
    )
    fig.text(
        0.5,
        0.022,
        "Finest-pair rates · BDF uses exact histories; CN uses a BDF1 start · "
        "compatible OSS subscales vanish",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.047, 1, 0.915), h_pad=2)
    save_figure(fig, output, "temporal_convergence")


def write_tables(summary: dict, output: Path) -> None:
    (output / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    for degree, result in summary["spatial"].items():
        table = result["table"]
        with (output / f"spatial_{degree.lower()}.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(table[0]))
            writer.writeheader()
            writer.writerows(table)
    with (output / "temporal.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("scheme", "dt", "velocity_error_l2", "temperature_error_l2"))
        for scheme, _, _ in SCHEMES:
            errors = summary["temporal"][scheme]["errors"]
            writer.writerows(
                (scheme, dt, u, t)
                for dt, u, t in zip(
                    summary["temporal_steps"],
                    errors["velocity"],
                    errors["temperature"],
                    strict=True,
                )
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"font.size": 10, "axes.titleweight": "bold", "pdf.fonttype": 42}):
        plot_spatial(summary, args.output)
        plot_temporal(summary, args.output)
    write_tables(summary, args.output)


if __name__ == "__main__":
    main()
