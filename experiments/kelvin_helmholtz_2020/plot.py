"""Plot absolute vorticity at Figure 4 times and compare saved KH time histories."""

import argparse
import csv
import json
import math
import tomllib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
TIMES = (1.0, 4.5, 5.0, 5.5, 6.0, 6.5, 10.0)


def read_csv(path: Path) -> list[dict[str, float]]:
    """Read numeric diagnostics, preserving missing values as NaN."""
    with path.open() as stream:
        return [
            {k: float(v) if v else math.nan for k, v in r.items()} for r in csv.DictReader(stream)
        ]


def load_run(folder: Path, form: str) -> dict:
    """Check the selected formulation and collect only existing snapshots."""
    cfg = tomllib.loads((folder / "config_resolved.toml").read_text())
    if cfg["case"] != "kelvin_helmholtz" or cfg["convection"]["form"] != form:
        raise ValueError(f"{folder}: expected Kelvin–Helmholtz with {form}")
    geometry = np.load(folder / "mesh.npz")["vertices"]
    initial = json.loads((folder / "initial_state.json").read_text())
    enstrophy = read_csv(folder / "enstrophy.csv")
    a, b = geometry[:, 1] - geometry[:, 0], geometry[:, 2] - geometry[:, 0]
    areas = np.abs(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]) / 2
    if not np.isclose(areas.sum(), initial["domain_volume"], rtol=1e-12, atol=1e-12):
        raise ValueError(f"{folder}: snapshot triangles do not cover the domain")
    flat = geometry.reshape(-1, 2)
    tri = mtri.Triangulation(flat[:, 0], flat[:, 1], np.arange(len(flat)).reshape(-1, 3))
    snapshots = {}
    for path in sorted(folder.glob("snapshot_*.npz")):
        with np.load(path) as data:
            if data["vorticity"].shape != geometry.shape[:2]:
                raise ValueError(f"snapshot/mesh shape mismatch: {path}")
            values = data["vorticity"]
            # Exact triangle integral for piecewise linear vorticity. This
            # also detects missing/duplicated MPI cells and non-vertex exports.
            reconstructed = float(
                np.sum(areas * (values.sum(axis=1) ** 2 + (values**2).sum(axis=1)) / 24)
            )
            t = float(data["time"])
            measured = next((r["enstrophy"] for r in enstrophy if abs(r["time"] - t) < 1e-9), None)
            if measured is None or not np.isclose(reconstructed, measured, rtol=1e-9, atol=1e-12):
                raise ValueError(f"{path}: snapshot does not reproduce the enstrophy diagnostic")
            snapshots[float(data["time"])] = data["vorticity"].copy().ravel()
    return dict(
        folder=folder,
        config=cfg,
        triangulation=tri,
        snapshots=snapshots,
        initial=initial,
        info=json.loads((folder / "run_info.json").read_text()),
        diagnostics=read_csv(folder / "diagnostics.csv"),
        enstrophy=enstrophy,
    )


def main() -> None:
    """No interpolation in time; an unavailable requested frame is explicitly marked."""
    parser = argparse.ArgumentParser(description=__doc__)
    base = ROOT / "results/kelvin_helmholtz"
    parser.add_argument("--emac", type=Path, default=base / "kh2020_emac_vms_p1_n196")
    parser.add_argument("--skew", type=Path, default=base / "kh2020_skew1_p1_n196")
    parser.add_argument("--output", type=Path, default=base / "comparison_p1_n196")
    parser.add_argument("--times", nargs="+", type=float, default=list(TIMES))
    parser.add_argument("--vmax", type=float, default=80.0)
    args = parser.parse_args()
    if not math.isfinite(args.vmax) or args.vmax <= 0:
        parser.error("--vmax must be finite and positive")
    if any(not math.isfinite(t) or t < 0 for t in args.times):
        parser.error("--times must be finite and nonnegative")
    runs = [load_run(args.emac, "emac_vms"), load_run(args.skew, "skew_symmetric_1")]
    configs = []
    for r in runs:
        c = dict(r["config"])
        c.pop("run_name")
        c.pop("convection")
        configs.append(c)
    if configs[0] != configs[1]:
        parser.error("the two runs differ in settings beyond convection form and run name")
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42})
    labels, colors = ("EMAC–VMS", "SKEW1"), ("#bf4b0b", "#007e70")
    norm = Normalize(0, args.vmax)
    fig, axes = plt.subplots(
        len(args.times),
        2,
        figsize=(7.5, 3.2 * len(args.times)),
        squeeze=False,
        constrained_layout=True,
    )
    missing = []
    for row, t in enumerate(args.times):
        for col, run in enumerate(runs):
            ax = axes[row, col]
            found = next((k for k in run["snapshots"] if abs(k - t) < 1e-9), None)
            if found is None:
                ax.text(0.5, 0.5, "Snapshot not available", ha="center", va="center")
                missing.append(dict(form=labels[col], time=t))
            else:
                ax.tripcolor(
                    run["triangulation"],
                    np.abs(run["snapshots"][found]),
                    shading="gouraud",
                    cmap="magma",
                    norm=norm,
                    rasterized=True,
                )
            ax.set(xlim=(0, 1), ylim=(0, 1), aspect="equal", xticks=[], yticks=[])
            if row == 0:
                ax.set_title(labels[col])
            if col == 0:
                ax.set_ylabel(f"t = {t:g}")
    partial = missing or any(r["info"]["exit_status"] != 0 for r in runs)
    suffix = " · PARTIAL OUTPUT" if partial else ""
    fig.suptitle("Kelvin–Helmholtz · absolute resolved vorticity |ωₕ|" + suffix)
    fig.colorbar(
        ScalarMappable(norm=norm, cmap="magma"),
        ax=axes.ravel().tolist(),
        shrink=0.65,
        label=f"|ωₕ| (values above {args.vmax:g} saturated)",
        extend="max",
    )
    for ext in ("png", "pdf"):
        fig.savefig(args.output / f"figure4_vorticity.{ext}", dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    summary = []
    for label, color, run in zip(labels, colors, runs, strict=True):
        d, z, initial = run["diagnostics"], run["enstrophy"], run["initial"]
        td, tz = [0.0] + [r["time"] for r in d], [r["time"] for r in z]
        energy = [1.0] + [r["kinetic_energy"] / initial["kinetic_energy"] for r in d]
        enstrophy = [r["enstrophy"] / z[0]["enstrophy"] for r in z]
        momentum = [math.hypot(initial["linear_momentum_x"], initial["linear_momentum_y"])] + [
            math.hypot(r["linear_momentum_x"], r["linear_momentum_y"]) for r in d
        ]
        angular = [initial["angular_momentum"]] + [r["angular_momentum"] for r in d]
        axes[0, 0].plot(td, energy, color=color, label=label)
        axes[0, 1].plot(tz, enstrophy, color=color, label=label)
        axes[0, 2].plot(tz, [r["vorticity_max_abs"] for r in z], color=color, label=label)
        axes[1, 0].plot(td, momentum, color=color, label=label)
        axes[1, 1].plot(td, angular, color=color, label=label)
        axes[1, 2].plot(td[1:], [r["divergence_l2"] for r in d], color=color, label=label)
        window = [r for r in z if 5 <= r["time"] <= 7]
        peak = max(window, key=lambda r: r["enstrophy"]) if window else None
        summary.append(
            dict(
                form=label,
                completed_time=td[-1],
                exit_status=run["info"]["exit_status"],
                final_energy_fraction=energy[-1],
                final_enstrophy_fraction=enstrophy[-1],
                peak_enstrophy_t5_to_t7=None if peak is None else peak["enstrophy"],
                peak_time_t5_to_t7=None if peak is None else peak["time"],
                max_normal_wall_velocity=max(r["normal_velocity_wall_l2"] for r in z),
            )
        )
    titles = (
        "Resolved energy E/E(0)",
        "Resolved enstrophy Z/Z(0)",
        "Maximum |ωₕ|",
        "Total momentum magnitude",
        "Angular momentum about origin",
        "L2 divergence",
    )
    for ax, title in zip(axes.ravel(), titles, strict=True):
        ax.set(title=title, xlabel="Time")
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle("Matched VMS comparison · raw histories, no smoothing" + suffix)
    for ext in ("png", "pdf"):
        fig.savefig(args.output / f"histories.{ext}", dpi=170)
    plt.close(fig)
    (args.output / "summary.json").write_text(
        json.dumps(dict(runs=summary, missing_snapshots=missing, partial=bool(partial)), indent=2)
        + "\n"
    )


if __name__ == "__main__":
    main()
