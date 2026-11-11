"""Matched-dissipation comparison: error against energy loss as τ_m varies (python matched.py)."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).parent
RESULTS = HERE / "results"
OUT = HERE / "presentation"
SCALES = {
    "emac_vms": (("oss", 1.0), ("oss_s0.5", 0.5), ("oss_s0.25", 0.25)),
    "skew_symmetric_1": (
        ("oss_s0.25", 0.25),
        ("oss_s0.5", 0.5),
        ("oss", 1.0),
        ("oss_s2", 2.0),
        ("oss_s4", 4.0),
    ),
}
LABEL = {"emac_vms": "EMAC–VMS", "skew_symmetric_1": "SKEW"}
COLOR = {"emac_vms": "#d4570f", "skew_symmetric_1": "#169c78"}
# (title, run pattern with {form} and {variant}, table, energy column, error column or None)
CASES = (
    (
        "Translating smooth vortex (E4, 32², T = 4)",
        "E4_{form}_{variant}",
        "extra",
        "energy",
        "error_l2",
    ),
    (
        "Smooth vortex in a box (E5, 48², T = 3)",
        "E5_{form}_{variant}",
        "diagnostics",
        "kinetic_energy",
        "velocity_error_l2",
    ),
    (
        "Gresho vortex in a box (E3, 96², T = 3)",
        "E3_{form}_{variant}_96",
        "diagnostics",
        "kinetic_energy",
        "velocity_error_l2",
    ),
)


def rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream))


def number(text: str | None) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def point(
    pattern: str, form: str, variant: str, table: str, energy: str, error: str
) -> tuple[float, float] | None:
    data = rows(RESULTS / pattern.format(form=form, variant=variant) / f"{table}.csv")
    if not data:
        return None
    e = [number(r[energy]) for r in data]
    return 1 - e[-1] / e[0], number(data[-1][error])


def shear_point(form: str, variant: str, time: float = 1.2) -> tuple[float, float] | None:
    data = rows(RESULTS / f"E2_{form}_{variant}_64" / "extra.csv")
    if not data:
        return None
    e0, z0 = number(data[0]["energy"]), number(data[0]["enstrophy"])
    row = min(data, key=lambda r: abs(number(r["time"]) - time))
    return 1 - number(row["energy"]) / e0, number(row["enstrophy"]) / z0


def main() -> None:
    OUT.mkdir(exist_ok=True)
    lines = [
        "# Matched-dissipation comparison",
        "",
        "τ_m multiplied by s (c1 and c2 divided by s), τ_c unchanged; s = 1 is the default.",
        "",
    ]
    figure, axes = plt.subplots(1, len(CASES) + 1, figsize=(18, 4.4))
    for axis, (title, pattern, table, energy, error) in zip(axes, CASES, strict=False):
        lines += [
            f"## {title}",
            "",
            "| Form | s | energy lost | L² error |",
            "| --- | ---: | ---: | ---: |",
        ]
        for form, variants in SCALES.items():
            points = []
            for variant, scale in variants:
                found = point(pattern, form, variant, table, energy, error)
                if found is None:
                    lines.append(f"| {LABEL[form]} | {scale:g} | missing | |")
                    continue
                points.append((scale, *found))
                lines.append(f"| {LABEL[form]} | {scale:g} | {found[0]:.4f} | {found[1]:.4f} |")
            if points:
                points.sort(key=lambda p: p[1])
                axis.plot(
                    [p[1] for p in points],
                    [p[2] for p in points],
                    "o-",
                    color=COLOR[form],
                    lw=2,
                    label=LABEL[form],
                )
                for scale, loss, err in points:
                    axis.annotate(
                        f"s={scale:g}",
                        (loss, err),
                        textcoords="offset points",
                        xytext=(4, 4),
                        fontsize=8,
                        color=COLOR[form],
                    )
        axis.set_title(title, fontsize=10, loc="left", fontweight="bold")
        axis.set_xlabel("energy lost, 1 − E(T)/E(0)")
        axis.set_ylabel("L² velocity error at T")
        axis.grid(alpha=0.3)
        lines.append("")
    axis = axes[-1]
    lines += [
        "## Thin shear layer (E2, 64², t = 1.2)",
        "",
        "| Form | s | energy lost | Z/Z(0) |",
        "| --- | ---: | ---: | ---: |",
    ]
    for form, variants in SCALES.items():
        points = []
        for variant, scale in variants:
            found = shear_point(form, variant)
            if found is None:
                lines.append(f"| {LABEL[form]} | {scale:g} | missing | |")
                continue
            points.append((scale, *found))
            lines.append(f"| {LABEL[form]} | {scale:g} | {found[0]:.4f} | {found[1]:.4f} |")
        if points:
            points.sort(key=lambda p: p[1])
            axis.plot(
                [p[1] for p in points],
                [p[2] for p in points],
                "o-",
                color=COLOR[form],
                lw=2,
                label=LABEL[form],
            )
            for scale, loss, z in points:
                axis.annotate(
                    f"s={scale:g}",
                    (loss, z),
                    textcoords="offset points",
                    xytext=(4, 4),
                    fontsize=8,
                    color=COLOR[form],
                )
    references = []
    for form in SCALES:
        data = rows(RESULTS / f"E2_{form}_oss_128" / "extra.csv")
        if data:
            row = min(data, key=lambda r: abs(number(r["time"]) - 1.2))
            references.append(
                (
                    1 - number(row["energy"]) / number(data[0]["energy"]),
                    number(row["enstrophy"]) / number(data[0]["enstrophy"]),
                )
            )
    for loss, z in references:
        axis.plot(loss, z, "k*", ms=12)
    if references:
        axis.plot([], [], "k*", ms=10, label="128² runs")
    axis.set_title(
        "Thin shear layer (E2, 64², t = 1.2)", fontsize=10, loc="left", fontweight="bold"
    )
    axis.set_xlabel("energy lost, 1 − E/E(0)")
    axis.set_ylabel("enstrophy Z/Z(0)")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8)
    axes[0].legend(fontsize=9)
    figure.suptitle(
        "Matched dissipation: τ_m × s with τ_c fixed. Below and left is better; "
        "compare the curves at equal energy loss",
        fontweight="bold",
    )
    figure.tight_layout()
    figure.savefig(OUT / "8_matched_dissipation.png", dpi=150)
    plt.close(figure)
    shear_snapshots()
    (HERE / "matched.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))  # noqa: T201 — command-line summary


def shear_snapshots() -> None:
    """Vorticity at t = 1.0 and 1.2 for every τ_m scaling, next to the 128² runs."""
    runs = [
        ("E2_skew_symmetric_1_oss_128", "SKEW 128²"),
        ("E2_emac_vms_oss_128", "EMAC–VMS 128²"),
        ("E2_skew_symmetric_1_oss_s0.25_64", "SKEW s=0.25"),
        ("E2_skew_symmetric_1_oss_s0.5_64", "SKEW s=0.5"),
        ("E2_skew_symmetric_1_oss_64", "SKEW s=1"),
        ("E2_skew_symmetric_1_oss_s2_64", "SKEW s=2"),
        ("E2_skew_symmetric_1_oss_s4_64", "SKEW s=4"),
        ("E2_emac_vms_oss_64", "EMAC–VMS s=1"),
        ("E2_emac_vms_oss_s0.5_64", "EMAC–VMS s=0.5"),
        ("E2_emac_vms_oss_s0.25_64", "EMAC–VMS s=0.25"),
    ]
    keys = ("t1.00", "t1.20")
    figure, axes = plt.subplots(
        len(keys), len(runs), figsize=(2.5 * len(runs), 2.8 * len(keys)), squeeze=False
    )
    for column, (run, label) in enumerate(runs):
        path = RESULTS / run / "vorticity.npz"
        data = np.load(path) if path.exists() else None
        for row, key in enumerate(keys):
            axis = axes[row, column]
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_aspect("equal")
            if data is not None and key in data:
                centres = data["centres"]
                axis.tripcolor(
                    centres[:, 0],
                    centres[:, 1],
                    data[key],
                    cmap="RdBu_r",
                    vmin=-40,
                    vmax=40,
                    shading="gouraud",
                )
            if row == 0:
                axis.set_title(label, fontsize=9)
            if column == 0:
                axis.set_ylabel(key.replace("t", "t = "))
    figure.suptitle(
        "Thin shear layer, 64², vorticity ±40: does matched dissipation change the braid topology?",
        fontweight="bold",
    )
    figure.tight_layout()
    figure.savefig(OUT / "9_shear_layer_matched.png", dpi=110)
    plt.close(figure)


if __name__ == "__main__":
    main()
