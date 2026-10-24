"""Figure 10: rigid-rotation and shear robustness (E6, E7): python robustness_figure.py."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).parent
RESULTS = HERE / "results"
SPEEDS = (1, 3, 10, 30, 100)
FORMS = ("nonconservative", "conservative", "skew_symmetric_1", "emac", "emac_vms")
LABEL = {
    "nonconservative": "CONV",
    "conservative": "CONS",
    "skew_symmetric_1": "SKEW",
    "emac": "EMAC",
    "emac_vms": "EMAC–VMS",
}
COLOR = {
    "nonconservative": "#1f6fb2",
    "conservative": "#e0a020",
    "skew_symmetric_1": "#169c78",
    "emac": "#a0609a",
    "emac_vms": "#d4570f",
}
# Relative error: ‖u_h − u‖ / ‖u‖ on the unit box, ‖C(−y, x)‖ = C/√6 and ‖C(y, 0)‖ = C/√12.
NORM = {"E6": 1 / math.sqrt(6), "E7": 1 / math.sqrt(12)}
FLOOR = 1e-17
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})


def relative_error(experiment: str, form: str, speed: int, cells: int) -> float | None:
    """Final relative error, or None when the run failed."""
    folder = RESULTS / f"{experiment}_{form}_C{speed}_{cells}"
    status = (folder / "status.txt").read_text().strip() if (folder / "status.txt").exists() else ""
    with (folder / "diagnostics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if status != "0" or not rows:
        return None
    return float(rows[-1]["velocity_error_l2"]) / (speed * NORM[experiment])


def flow_inset(axis: plt.Axes, experiment: str) -> None:
    inset = axis.inset_axes((0.04, 0.33, 0.2, 0.27))
    grid = np.linspace(-0.4, 0.4, 7)
    x, y = np.meshgrid(grid, grid)
    u, v = (-y, x) if experiment == "E6" else (y, 0 * x)
    inset.quiver(x, y, u, v, color="#333", scale=4, width=0.012)
    inset.set_xlim(-0.5, 0.5)
    inset.set_ylim(-0.5, 0.5)
    inset.set_xticks([])
    inset.set_yticks([])
    inset.set_aspect("equal")
    for spine in inset.spines.values():
        spine.set_visible(True)
        spine.set_color("#999")


def main() -> None:
    figure, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    panels = (
        (
            "E6",
            "Rigid rotation  u = C(−y, x),  p = C²r²/2",
            "classical forms: P1 pressure must balance the centripetal gradient\n"
            "EMAC forms: ∇ˢu = 0 and P = p − ½|u|² = 0  →  exact",
        ),
        (
            "E7",
            "Parallel shear  u = C(y, 0),  p = const",
            "classical forms: (u·∇)u = 0  →  exact\n"
            "EMAC forms: P1 pressure must balance ∇½|u|² = C²y e_y",
        ),
    )
    for axis, (experiment, title, note) in zip(axes, panels, strict=True):
        axis.axhspan(FLOOR, 1e-13, color="#e6e6e6", zorder=0)
        axis.text(
            0.98,
            0.17,
            "round-off level = exact",
            transform=axis.transAxes,
            ha="right",
            fontsize=9,
            color="#555",
        )
        for form in FORMS:
            for cells, style, marker in ((16, "-", "o"), (32, "--", "s")):
                values = [relative_error(experiment, form, c, cells) for c in SPEEDS]
                ok = [
                    (c, max(e, FLOOR)) for c, e in zip(SPEEDS, values, strict=True) if e is not None
                ]
                width = 2.8 if form == "emac_vms" else 1.6
                axis.plot(
                    [c for c, _ in ok],
                    [e for _, e in ok],
                    style,
                    marker=marker,
                    ms=5,
                    color=COLOR[form],
                    lw=width,
                    mfc=COLOR[form] if cells == 16 else "white",
                    label=f"{LABEL[form]}, {cells}²",
                )
                failed = [c for c, e in zip(SPEEDS, values, strict=True) if e is None]
                for c in failed:
                    axis.plot(c, 2e-1, "x", color=COLOR[form], ms=11, mew=2.5)
                    axis.annotate(
                        "no convergence",
                        (c, 2e-1),
                        textcoords="offset points",
                        xytext=(-90, 6),
                        fontsize=8,
                        color=COLOR[form],
                    )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_ylim(FLOOR / 3, 1)
        axis.set_xlabel("rotation or shear rate C")
        axis.set_title(title, loc="left", fontweight="bold", fontsize=11)
        axis.text(
            0.98,
            0.47,
            note,
            transform=axis.transAxes,
            ha="right",
            va="center",
            fontsize=9,
            color="#333",
            bbox={"boxstyle": "round", "fc": "white", "ec": "#bbb"},
        )
        axis.grid(alpha=0.3, which="major")
        flow_inset(axis, experiment)
    axes[0].set_ylabel("relative L² velocity error after 10 steps")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=5, fontsize=9, frameon=False)
    figure.suptitle(
        "Robustness: each form is exact for a different flow type "
        "(P1/P1 dynamic OSS, ν = 1e-6, exact Dirichlet data)",
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0.1, 1, 0.97))
    out = HERE / "presentation" / "10_rotation_shear_robustness.png"
    figure.savefig(out, dpi=160)
    plt.close(figure)
    print(out)  # noqa: T201 — command-line summary


if __name__ == "__main__":
    main()
