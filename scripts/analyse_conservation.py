"""Tables and plots of the Phase 12b conservation runs: python scripts/analyse_conservation.py.

Every residual, defect and drift is reported in absolute terms and divided by the fixed scales of
its run (spec: Conservation-run normalization, owner decision K-D2). The script only reads the
results tree; a run without its files is reported as missing or incomplete.
"""

import argparse
import csv
import json
import logging
import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from matplotlib.figure import Figure

from navier_stokes_vms.config import Config, load_config

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ("gresho", "gresho_long", "gresho_refined", "taylor_green")
CLOSURE_LIMIT = 1e-9
EMAC_VMS_LIMIT = 1e-12
DEFECT_LIMIT = 1e-9
Row = Mapping[str, str]


@dataclass(frozen=True)
class Scales:
    """U_c = √(2E(0)/|Ω|), L_c = |Ω|^(1/d) and the scales built from them (owner decision K-D2)."""

    energy: float
    momentum: float
    angular_momentum: float
    force: float
    torque: float
    power: float


def run_scales(initial: Mapping[str, float]) -> Scales | None:
    """The scales of one run; None when E(0) = 0, so that nothing is divided by zero."""
    volume, energy = initial["domain_volume"], initial["kinetic_energy"]
    if energy <= 0:
        return None
    velocity = math.sqrt(2 * energy / volume)
    length = volume ** (1 / initial["dimension"])
    force = volume * velocity**2 / length
    return Scales(
        energy=energy,
        momentum=volume * velocity,
        angular_momentum=volume * velocity * length,
        force=force,
        torque=force * length,
        power=volume * velocity**3 / length,
    )


def family(run_name: str) -> str:
    """gresho_long and gresho_refined before gresho; taylor_green for the periodic runs."""
    for prefix in ("gresho_long", "gresho_refined", "gresho", "taylor_green"):
        if run_name.startswith(f"{prefix}_"):
            return prefix
    raise ValueError(f"run {run_name!r} belongs to no conservation family")


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def value(row: Row, column: str) -> float | None:
    """A CSV entry as a number; blank or absent entries are None."""
    text = row.get(column, "")
    return float(text) if text else None


def run_status(folder: Path) -> str:
    """missing, incomplete (no recorded exit status), completed or failed."""
    if not folder.is_dir():
        return "missing"
    info_path = folder / "run_info.json"
    info = json.loads(info_path.read_text(encoding="utf-8")) if info_path.exists() else {}
    exit_status = info.get("exit_status")
    if exit_status is None:
        return "incomplete"
    return "completed" if exit_status == 0 else "failed"


def _maximum(values: Sequence[float]) -> float | None:
    return max(values) if values else None


def _scaled(number: float | None, scale: float | None) -> float | None:
    return None if number is None or scale is None else number / scale


def _defect_columns(axes: str, angular: bool) -> list[tuple[str, list[tuple[str, str]], str]]:
    """(label, [(measured, predicted)] column pairs, scale name) for each defect family."""
    columns = [
        (
            "linear_defect",
            [
                (f"linear_momentum_defect_{a}", f"linear_momentum_defect_predicted_{a}")
                for a in axes
            ],
            "force",
        )
    ]
    if angular:
        columns.append(
            (
                "angular_defect",
                [("angular_momentum_defect", "angular_momentum_defect_predicted")],
                "torque",
            )
        )
    columns.append(
        (
            "energy_defect",
            [("convective_energy_defect", "energy_budget_defect_predicted")],
            "power",
        )
    )
    return columns


def summarise(config: Config, folder: Path) -> dict[str, object]:
    """One table row: status, iterations, closures, defects, drifts and final monitors."""
    stabilization = config.stabilization
    summary: dict[str, object] = {
        "run": config.run_name,
        "form": config.convection.form,
        "space": stabilization.subscale_space,
        "splitting": stabilization.scale_splitting,
    }
    status = run_status(folder)
    steps = read_rows(folder / "time_steps.csv")
    diagnostics = read_rows(folder / "diagnostics.csv")
    converged = [row for row in steps if row["converged"] == "True"]
    summary.update(
        status=status,
        steps_expected=round(config.time.t_end / config.time.dt),
        steps_converged=len(converged),
        failed_step=next((int(row["step"]) for row in steps if row["converged"] != "True"), None),
        outer_iterations_max=_maximum([int(row["outer_iterations"]) for row in steps]),
        inner_iterations_max=_maximum([int(row["inner_iterations_total"]) for row in steps]),
        inner_iterations_total=sum(int(row["inner_iterations_total"]) for row in steps),
    )
    initial_path = folder / "initial_state.json"
    initial = (
        json.loads(initial_path.read_text(encoding="utf-8")) if initial_path.exists() else None
    )
    scales = None if initial is None else run_scales(initial)
    axes = "xyz"[: int(initial["dimension"])] if initial else "xy"
    periodic = config.case == "taylor_green_2d"

    energy_residuals = [
        abs(number)
        for row in diagnostics
        if (number := value(row, "energy_budget_residual")) is not None
    ]
    momentum_residuals = [
        abs(number)
        for row in diagnostics
        for axis in axes
        if (number := value(row, f"momentum_budget_residual_{axis}")) is not None
    ]
    summary["energy_closure"] = _maximum(energy_residuals)
    summary["energy_closure_scaled"] = _scaled(summary["energy_closure"], scales and scales.power)
    summary["momentum_closure"] = _maximum(momentum_residuals)
    summary["momentum_closure_scaled"] = _scaled(
        summary["momentum_closure"], scales and scales.force
    )

    defects_within = []
    for label, pairs, scale_name in _defect_columns(axes, angular=not periodic):
        scale = None if scales is None else getattr(scales, scale_name)
        measured, predicted, differences = [], [], []
        for row in diagnostics:
            for measured_column, predicted_column in pairs:
                m, p = value(row, measured_column), value(row, predicted_column)
                if m is None or p is None:
                    continue
                measured.append(abs(m))
                predicted.append(abs(p))
                differences.append(abs(m - p))
        for kind, numbers in (
            ("measured", measured),
            ("predicted", predicted),
            ("difference", differences),
        ):
            summary[f"{label}_{kind}"] = _maximum(numbers)
            summary[f"{label}_{kind}_scaled"] = _scaled(_maximum(numbers), scale)
        if config.convection.form == "emac_vms":
            check, limit = summary[f"{label}_measured_scaled"], EMAC_VMS_LIMIT
        else:
            check, limit = summary[f"{label}_difference_scaled"], DEFECT_LIMIT
        defects_within.append(None if check is None else check <= limit)

    summary["energy_closure_within"] = _within(summary["energy_closure_scaled"], CLOSURE_LIMIT)
    summary["momentum_closure_within"] = (
        _within(summary["momentum_closure_scaled"], CLOSURE_LIMIT) if periodic else None
    )
    summary["defects_within"] = (
        None if None in defects_within or not defects_within else all(defects_within)
    )
    summary.update(_drifts(diagnostics, initial, scales, axes, angular=not periodic))
    final = diagnostics[-1] if diagnostics else {}
    summary["velocity_error_final"] = value(final, "velocity_error_l2")
    if not periodic:
        summary["region_velocity_final"] = value(final, "region_velocity_l2")
        summary["region_subscale_final"] = value(final, "region_subscale_l2")
        forces = [
            abs(number)
            for column in final
            if column.startswith("force_") and (number := value(final, column)) is not None
        ]
        summary["wall_force_final_max"] = _maximum(forces)
        summary["wall_force_final_max_scaled"] = _scaled(
            summary["wall_force_final_max"], scales and scales.force
        )
        for column in sorted(column for column in final if column.startswith("force_")):
            summary[f"{column}_final"] = value(final, column)
    return summary


def _within(number: object, limit: float) -> bool | None:
    return None if number is None else number <= limit


def _drifts(
    diagnostics: Sequence[Row],
    initial: Mapping[str, float] | None,
    scales: Scales | None,
    axes: str,
    angular: bool,
) -> dict[str, float | None]:
    """Final and largest drifts of M, A and E from the recorded initial state."""
    if initial is None or not diagnostics:
        return {}
    series: dict[str, tuple[list[float], float | None]] = {}
    momentum = []
    for row in diagnostics:
        components = [value(row, f"linear_momentum_{axis}") for axis in axes]
        if None not in components:
            momentum.append(
                math.sqrt(
                    sum(
                        (component - initial[f"linear_momentum_{axis}"]) ** 2
                        for component, axis in zip(components, axes, strict=True)
                    )
                )
            )
    series["momentum_drift"] = (momentum, scales and scales.momentum)
    if angular:
        series["angular_momentum_drift"] = (
            [
                number - initial["angular_momentum"]
                for row in diagnostics
                if (number := value(row, "angular_momentum")) is not None
            ],
            scales and scales.angular_momentum,
        )
    series["energy_drift"] = (
        [
            number - initial["kinetic_energy"]
            for row in diagnostics
            if (number := value(row, "kinetic_energy")) is not None
        ],
        scales and scales.energy,
    )
    drifts: dict[str, float | None] = {}
    for name, (numbers, scale) in series.items():
        final = numbers[-1] if numbers else None
        largest = _maximum([abs(number) for number in numbers])
        drifts.update(
            {
                f"{name}_final": final,
                f"{name}_final_scaled": _scaled(final, scale),
                f"{name}_max": largest,
                f"{name}_max_scaled": _scaled(largest, scale),
            }
        )
    return drifts


MARKDOWN_COLUMNS = (
    ("run", "run"),
    ("status", "status"),
    ("steps_converged", "steps"),
    ("outer_iterations_max", "max outer"),
    ("inner_iterations_max", "max inner"),
    ("energy_closure_scaled", "max |r|/P_c"),
    ("momentum_closure_scaled", "max |r_M|/F_c"),
    ("linear_defect_difference_scaled", "max |D_k − pred|/F_c"),
    ("linear_defect_predicted_scaled", "max |pred D_k|/F_c"),
    ("angular_defect_difference_scaled", "max |D_φ − pred|/(F_c L_c)"),
    ("angular_defect_predicted_scaled", "max |pred D_φ|/(F_c L_c)"),
    ("energy_defect_difference_scaled", "max |𝒩 − pred|/P_c"),
    ("energy_defect_predicted_scaled", "max |pred 𝒩|/P_c"),
    ("momentum_drift_final_scaled", "|M(T) − M(0)|/M_c"),
    ("angular_momentum_drift_final_scaled", "(A(T) − A(0))/A_c"),
    ("energy_drift_final_scaled", "(E(T) − E(0))/E(0)"),
    ("region_velocity_final", "frame ‖u_h‖"),
    ("region_subscale_final", "frame ‖ũ‖"),
    ("wall_force_final_max_scaled", "max wall force/F_c"),
    ("velocity_error_final", "‖u_h − u‖"),
    ("energy_closure_within", "r ok"),
    ("momentum_closure_within", "r_M ok"),
    ("defects_within", "defects ok"),
)


def _text(entry: object) -> str:
    if entry is None:
        return "—"
    if isinstance(entry, bool):
        return "yes" if entry else "no"
    if isinstance(entry, float):
        return f"{entry:.3e}"
    return str(entry)


def write_tables(name: str, rows: Sequence[Mapping[str, object]], folder: Path) -> str:
    """<family>.csv with every column and <family>.md with the main ones; returns the Markdown."""
    columns = list(dict.fromkeys(column for row in rows for column in row))
    with (folder / f"{name}.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {column: "" if row.get(column) is None else row.get(column) for column in columns}
            )
    shown = [(key, title) for key, title in MARKDOWN_COLUMNS if key in columns]
    lines = [
        f"## {name}",
        "",
        "| " + " | ".join(title for _, title in shown) + " |",
        "| " + " | ".join("---" for _ in shown) + " |",
    ]
    lines += ["| " + " | ".join(_text(row.get(key)) for key, _ in shown) + " |" for row in rows]
    text = "\n".join(lines) + "\n"
    (folder / f"{name}.md").write_text(text, encoding="utf-8")
    return text


def failure_report(
    rows: Sequence[Mapping[str, object]], results: Path, configs: Mapping[str, Config]
) -> str:
    """Every run that did not complete, with its iteration counts and the logged reason."""
    lines = ["## Runs that did not complete", ""]
    for row in rows:
        if row["status"] == "completed":
            continue
        config = configs[str(row["run"])]
        log = results / config.case / config.run_name / "run.log"
        reasons = [
            line.split("Run stopped: ", 1)[1]
            for line in (log.read_text(encoding="utf-8").splitlines() if log.exists() else [])
            if "Run stopped: " in line
        ]
        lines.append(
            f"- {row['run']}: {row['status']}, {row['steps_converged']} of {row['steps_expected']}"
            f" steps converged, failed step {_text(row['failed_step'])}, max outer"
            f" {_text(row['outer_iterations_max'])}, max inner {_text(row['inner_iterations_max'])}"
            + (f"; {reasons[-1]}" if reasons else "")
        )
    if len(lines) == 2:
        lines.append("None.")
    return "\n".join(lines) + "\n"


def plot_long_gresho(configs: Sequence[Config], results: Path, path: Path) -> bool:
    """Energy, both momentum drifts and the velocity error of the long Gresho runs, per form."""
    figure = Figure(figsize=(11, 8), layout="constrained")
    axes = figure.subplots(2, 2, sharex=True)
    panels = (
        ("E(t)/E(0)", axes[0, 0]),
        ("|M(t) − M(0)|/M_c", axes[0, 1]),
        ("(A(t) − A(0))/A_c", axes[1, 0]),
        ("‖u_h − u‖ (L²)", axes[1, 1]),
    )
    plotted = False
    for config in configs:
        folder = results / config.case / config.run_name
        rows = read_rows(folder / "diagnostics.csv")
        initial_path = folder / "initial_state.json"
        if not rows or not initial_path.exists():
            continue
        initial = json.loads(initial_path.read_text(encoding="utf-8"))
        scales = run_scales(initial)
        if scales is None:
            continue
        times = [float(row["time"]) for row in rows]
        momentum = [
            math.hypot(
                float(row["linear_momentum_x"]) - initial["linear_momentum_x"],
                float(row["linear_momentum_y"]) - initial["linear_momentum_y"],
            )
            / scales.momentum
            for row in rows
        ]
        # E, M and A start from the recorded initial state; the error is known from step 1.
        series = (
            ([0.0, *times], [1.0] + [float(row["kinetic_energy"]) / scales.energy for row in rows]),
            ([0.0, *times], [0.0, *momentum]),
            (
                [0.0, *times],
                [0.0]
                + [
                    (float(row["angular_momentum"]) - initial["angular_momentum"])
                    / scales.angular_momentum
                    for row in rows
                ],
            ),
            (times, [float(row["velocity_error_l2"]) for row in rows]),
        )
        for (title, axis), (abscissae, values) in zip(panels, series, strict=True):
            axis.plot(abscissae, values, label=config.convection.form)
            axis.set_title(title)
        plotted = True
    if not plotted:
        return False
    for _, axis in panels:
        axis.set_xlabel("t")
        axis.grid(True, alpha=0.3)
    axes[0, 0].legend()
    figure.suptitle("Long Gresho runs: Dyn-Nl-OSS, P1/P1, 32×32, δt = 0.01")
    figure.savefig(path, dpi=150)
    return True


def analyse(
    configs_folder: Path, results: Path, output: Path
) -> dict[str, list[dict[str, object]]]:
    """Write every family table, the failure report, the long-Gresho plot and summary.md."""
    configs = {path.stem: load_config(path) for path in sorted(configs_folder.glob("*.toml"))}
    output.mkdir(parents=True, exist_ok=True)
    tables: dict[str, list[dict[str, object]]] = {name: [] for name in FAMILIES}
    for name, config in configs.items():
        folder = results / config.case / config.run_name
        tables[family(name)].append(summarise(config, folder))
    sections = [
        "# Conservation runs (Phase 12b)",
        "",
        "Scales per run (owner decision K-D2): U_c = √(2E(0)/|Ω|), L_c = |Ω|^(1/d), "
        "M_c = |Ω|U_c, A_c = M_c L_c, F_c = |Ω|U_c²/L_c, torque F_c L_c, "
        "P_c = |Ω|U_c³/L_c, from initial_state.json. Closures use steps 2 onward. Expectations: "
        f"|r| ≤ {CLOSURE_LIMIT:g} P_c; |r_M| ≤ {CLOSURE_LIMIT:g} F_c on Taylor–Green; EMAC–VMS "
        f"measured defects ≤ {EMAC_VMS_LIMIT:g} of their scales; other forms "
        f"|measured − predicted| ≤ {DEFECT_LIMIT:g} of their scales. On Gresho, r_M is the "
        "wall reaction and is only reported. Drifts use the resolved diagnostics.",
        "",
    ]
    for name, rows in tables.items():
        if rows:
            sections.append(write_tables(name, rows, output))
    all_rows = [row for rows in tables.values() for row in rows]
    sections.append(failure_report(all_rows, results, configs))
    long_configs = [configs[str(row["run"])] for row in tables["gresho_long"]]
    if plot_long_gresho(long_configs, results, output / "gresho_long.png"):
        sections.append("![Long Gresho runs](gresho_long.png)\n")
    (output / "summary.md").write_text("\n".join(sections), encoding="utf-8")
    return tables


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--configs", type=Path, default=ROOT / "configs" / "chapter5" / "conservation"
    )
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--output", type=Path, default=Path("results") / "conservation_analysis")
    arguments = parser.parse_args(argv)
    # The format of run logs, without importing Firedrake for a file-only task.
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(Path(__file__).stem)
    tables = analyse(arguments.configs, arguments.results, arguments.output)
    statuses = [row["status"] for rows in tables.values() for row in rows]
    logger.info(
        "Analysed %d runs (%s); tables in %s",
        len(statuses),
        ", ".join(f"{statuses.count(s)} {s}" for s in sorted(set(statuses))),
        arguments.output,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
