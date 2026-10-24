"""Figures and summary tables for E1, E2 and E3: python report.py."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
FORMS = (
    "nonconservative",
    "conservative",
    "skew_symmetric_1",
    "skew_symmetric_2",
    "emac",
    "emac_vms",
)
LABEL = {
    "nonconservative": "Nonconservative",
    "conservative": "Conservative",
    "skew_symmetric_1": "Skew-symmetric 1",
    "skew_symmetric_2": "Skew-symmetric 2",
    "emac": "EMAC",
    "emac_vms": "EMAC–VMS",
}
STYLE = {
    "nonconservative": dict(color="#1f6fb2", lw=1.6),
    "conservative": dict(color="#e0a020", lw=1.6, ls="--"),
    "skew_symmetric_1": dict(color="#169c78", lw=1.6),
    "skew_symmetric_2": dict(color="#5bb8e0", lw=1.3, ls="--"),
    "emac": dict(color="#a0609a", lw=1.6, ls="-."),
    "emac_vms": dict(color="#d4570f", lw=2.6),
}
VARIANT_TITLE = {
    "oss": "dynamic OSS, default τ (c2 = 2, τ_c on)",
    "oss_weak": "dynamic OSS, weak τ (c2 = 32, τ_c off)",
    "asgs": "dynamic ASGS, default τ",
}


def table(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return {}
    data = {}
    for key in rows[0]:
        values = []
        for row in rows:
            try:
                values.append(float(row[key]))
            except (TypeError, ValueError):
                values.append(math.nan)
        data[key] = np.array(values)
    return data


def cumulative(values: np.ndarray, dt: float) -> np.ndarray:
    return np.cumsum(np.nan_to_num(values)) * dt


def status(run: str) -> str:
    path = RESULTS / run / "status.txt"
    if not path.exists():
        return "missing"
    extra = table(RESULTS / run / "extra.csv")
    reached = extra["time"][-1] if extra else 0.0
    return "ok" if path.read_text().strip() == "0" else f"failed at t = {reached:.2f}"


def panel_figure(
    experiment: str,
    variant: str,
    panels: list[tuple[str, str, object]],
    note: str,
    dt: float,
    suffix: str = "",
) -> Path | None:
    runs = {form: f"{experiment}_{form}_{variant}{suffix}" for form in FORMS}
    if not any((RESULTS / run).exists() for run in runs.values()):
        return None
    figure, axes = plt.subplots(2, (len(panels) + 1) // 2, figsize=(14, 8.2))
    axes = axes.ravel()
    for axis, (title, ylabel, quantity) in zip(axes, panels, strict=False):
        for form, run in runs.items():
            extra = table(RESULTS / run / "extra.csv")
            diagnostics = table(RESULTS / run / "diagnostics.csv")
            if not extra:
                continue
            values = quantity(extra, diagnostics, dt)
            if values is None:
                continue
            time = extra["time"][: len(values)]
            axis.plot(time, values[: len(time)], label=LABEL[form], **STYLE[form])
        axis.set_title(title, fontsize=11, loc="left", fontweight="bold")
        axis.set_ylabel(ylabel)
        axis.set_xlabel("time t")
        axis.grid(alpha=0.3)
    for axis in axes[len(panels) :]:
        axis.axis("off")
    axes[0].legend(fontsize=8, ncol=2)
    figure.suptitle(
        f"{experiment} · {VARIANT_TITLE.get(variant, variant)}{suffix}",
        x=0.01,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    figure.text(0.01, 0.005, note, fontsize=8, color="#444")
    figure.tight_layout(rect=(0, 0.03, 1, 0.96))
    FIGURES.mkdir(exist_ok=True)
    path = FIGURES / f"{experiment}_{variant}{suffix}.png"
    figure.savefig(path, dpi=130)
    plt.close(figure)
    return path


def e1() -> list[str]:
    def energy(e: dict, _d: dict, _dt: float) -> np.ndarray:
        return e["energy"] / e["energy"][0]

    def defect(_e: dict, d: dict, dt: float) -> np.ndarray | None:
        if "convective_energy_defect" not in d:
            return None
        return cumulative(d["convective_energy_defect"], dt) / 0.70703

    panels = [
        ("Resolved kinetic energy", "E(t)/E(0)", energy),
        ("L² velocity error vs exact", "‖u_h − u‖", lambda e, _d, _t: e["error_l2"]),
        ("Vortex-centre position error", "|x_c − U₀t|", lambda e, _d, _t: e["centre_error"]),
        ("Cumulative convective energy defect", "∫𝒩 dt / E(0)", defect),
        ("Resolved momentum drift (x)", "M_x(t) − M_x(0)", lambda e, _d, _t: e["momentum_x"] - 1.0),
        ("Enstrophy", "½‖ω_h‖²", lambda e, _d, _t: e["enstrophy"]),
    ]
    note = (
        "Translating Gresho vortex, periodic unit box, U₀ = (1, ½), ν = 1e-6, P1/P1 32×32, "
        "CN Δt = 0.01. Exact solution u = U₀ + u_G(x − U₀t). Resolved fields."
    )
    lines = [
        "## E1 — translating Gresho vortex (energy and momentum invariants)",
        "",
        "| Variant | Form | Status | E(T)/E(0) | ‖u_h − u‖(T) | centre error (T) | "
        "∫𝒩 dt / E(0) | max |M_x − 1| |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for variant in ("oss", "oss_weak", "asgs"):
        panel_figure("E1", variant, panels, note, 0.01)
        for form in FORMS:
            run = f"E1_{form}_{variant}"
            e, d = table(RESULTS / run / "extra.csv"), table(RESULTS / run / "diagnostics.csv")
            if not e:
                lines.append(f"| {variant} | {LABEL[form]} | {status(run)} | | | | | |")
                continue
            n = cumulative(d.get("convective_energy_defect", np.zeros(1)), 0.01)[-1] / 0.70703
            lines.append(
                f"| {variant} | {LABEL[form]} | {status(run)} | "
                f"{e['energy'][-1] / e['energy'][0]:.4f} | "
                f"{e['error_l2'][-1]:.4f} | {e['centre_error'][-1]:.4f} | {n:+.2e} | "
                f"{np.nanmax(np.abs(e['momentum_x'] - 1.0)):.1e} |"
            )
    return lines


def e3() -> list[str]:
    def angular(_e: dict, d: dict, _dt: float) -> np.ndarray:
        return d["angular_momentum"] / d["angular_momentum"][0] - 1

    def torque(_e: dict, d: dict, dt: float) -> np.ndarray:
        return cumulative(d["angular_momentum_defect"], dt) / abs(d["angular_momentum"][0])

    panels = [
        ("Resolved angular momentum drift", "A(t)/A(0) − 1", angular),
        ("Cumulative spurious convective torque", "∫D_φ dt / A(0)", torque),
        ("Resolved kinetic energy", "E(t)/E(0)", lambda e, _d, _t: e["energy"] / e["energy"][0]),
        ("L² velocity error vs exact", "‖u_h − u‖", lambda _e, d, _t: d["velocity_error_l2"]),
        (
            "Velocity in the wall frame (|x| or |y| ≥ 0.9)",
            "‖u_h‖_frame",
            lambda _e, d, _t: d["region_velocity_l2"],
        ),
        (
            "Cumulative convective energy defect",
            "∫𝒩 dt / E(0)",
            lambda e, d, dt: cumulative(d["convective_energy_defect"], dt) / e["energy"][0],
        ),
    ]
    note = (
        "Off-centre Gresho vortex (centre (0.1, 0.05)) in a no-slip box (−1, 1)², ν = 1e-6, "
        "P1/P1 48×48, CN Δt = 0.01. Exact steady Euler field; A about the origin."
    )
    lines = [
        "",
        "## E3 — Gresho in a large box (angular-momentum invariant)",
        "",
        "| Variant | Form | Status | A(T)/A(0) − 1 | ∫D_φ dt / A(0) | E(T)/E(0) | "
        "‖u_h − u‖(T) | "
        "frame ‖u_h‖ (T) |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for variant in ("oss", "oss_weak"):
        panel_figure("E3", variant, panels, note, 0.01)
        for form in FORMS:
            run = f"E3_{form}_{variant}"
            e, d = table(RESULTS / run / "extra.csv"), table(RESULTS / run / "diagnostics.csv")
            if not d:
                lines.append(f"| {variant} | {LABEL[form]} | {status(run)} | | | | | |")
                continue
            a0 = d["angular_momentum"][0]
            lines.append(
                f"| {variant} | {LABEL[form]} | {status(run)} | "
                f"{d['angular_momentum'][-1] / a0 - 1:+.2e} | "
                f"{cumulative(d['angular_momentum_defect'], 0.01)[-1] / abs(a0):+.2e} | "
                f"{e['energy'][-1] / e['energy'][0]:.4f} | {d['velocity_error_l2'][-1]:.4f} | "
                f"{d['region_velocity_l2'][-1]:.2e} |"
            )
    return lines


def e2() -> list[str]:
    lines = [
        "",
        "## E2 — thin double shear layer (under-resolved robustness)",
        "",
        "| Run | Status | E(T)/E(0) | Z(T)/Z(0) | ∫𝒩 dt / E(0) |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    runs = sorted(path.name for path in RESULTS.glob("E2_*"))
    for run in runs:
        e, d = table(RESULTS / run / "extra.csv"), table(RESULTS / run / "diagnostics.csv")
        if not e:
            lines.append(f"| {run} | {status(run)} | | | |")
            continue
        n = cumulative(d.get("convective_energy_defect", np.zeros(1)), 0.005)[-1] / e["energy"][0]
        lines.append(
            f"| {run} | {status(run)} | {e['energy'][-1] / e['energy'][0]:.4f} | "
            f"{e['enstrophy'][-1] / e['enstrophy'][0]:.4f} | {n:+.2e} |"
        )
    for cells in sorted({int(run.rsplit("_", 1)[1]) for run in runs}):
        group = [run for run in runs if run.endswith(f"_{cells}")]
        if not group:
            continue
        figure, axes = plt.subplots(1, 2, figsize=(13, 4.6))
        for run in group:
            form = next(
                f for f in sorted(FORMS, key=len, reverse=True) if run.startswith(f"E2_{f}_")
            )
            e = table(RESULTS / run / "extra.csv")
            if not e:
                continue
            axes[0].plot(e["time"], e["energy"] / e["energy"][0], label=LABEL[form], **STYLE[form])
            axes[1].plot(e["time"], e["enstrophy"] / e["enstrophy"][0], **STYLE[form])
        axes[0].set_title("Resolved kinetic energy E/E(0)", loc="left", fontweight="bold")
        axes[1].set_title("Resolved enstrophy Z/Z(0)", loc="left", fontweight="bold")
        for axis in axes:
            axis.grid(alpha=0.3)
            axis.set_xlabel("time t")
        axes[0].legend(fontsize=8)
        figure.suptitle(
            f"E2 · double shear layer · {cells}×{cells}", x=0.01, ha="left", fontweight="bold"
        )
        figure.tight_layout()
        figure.savefig(FIGURES / f"E2_{cells}_history.png", dpi=130)
        plt.close(figure)
        snapshot_figure(group, cells)
    return lines


def snapshot_figure(group: list[str], cells: int) -> None:
    keys = ("t1.00", "t1.20", "t1.50")
    figure, axes = plt.subplots(
        len(keys), len(group), figsize=(2.6 * len(group), 2.7 * len(keys)), squeeze=False
    )
    for column, run in enumerate(group):
        path = RESULTS / run / "vorticity.npz"
        if not path.exists():
            continue
        data = np.load(path)
        centres = data["centres"]
        for row, key in enumerate(keys):
            axis = axes[row, column]
            axis.set_xticks([])
            axis.set_yticks([])
            if key in data:
                axis.tripcolor(
                    centres[:, 0],
                    centres[:, 1],
                    data[key],
                    cmap="RdBu_r",
                    vmin=-40,
                    vmax=40,
                    shading="gouraud",
                )
            axis.set_aspect("equal")
            if row == 0:
                axis.set_title(run.removeprefix("E2_").rsplit("_", 2)[0], fontsize=8)
            if column == 0:
                axis.set_ylabel(key.replace("t", "t = "), fontsize=9)
    figure.suptitle(f"E2 · vorticity ω_h, {cells}×{cells}, colour scale ±40", fontweight="bold")
    figure.tight_layout()
    figure.savefig(FIGURES / f"E2_{cells}_vorticity.png", dpi=110)
    plt.close(figure)


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    lines = ["# EMAC–VMS advantage experiments", ""] + e1() + e3() + e2()
    (HERE / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))  # noqa: T201 — command-line summary


if __name__ == "__main__":
    main()
