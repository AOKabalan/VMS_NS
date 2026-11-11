"""Summarise the blow-up screening: one row per run, plus per-setting Markdown tables."""

import csv
import math
from pathlib import Path

HERE = Path(__file__).parent
RESULTS = HERE / "results" / "gresho"
T_END = 3.0
BLOW_UP = 1.5  # E(t)/E(0) above this counts as blow-up
FORMS = (
    "nonconservative",
    "conservative",
    "skew_symmetric_1",
    "skew_symmetric_2",
    "emac",
    "emac_vms",
)
SHORT = {
    "nonconservative": "CONV",
    "conservative": "CONS",
    "skew_symmetric_1": "SKEW1",
    "skew_symmetric_2": "SKEW2",
    "emac": "EMAC",
    "emac_vms": "EMAC–VMS",
}


def rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream))


def number(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def summarise(name: str, status: int | None, wall: int | None) -> dict[str, object]:
    folder = RESULTS / name
    diagnostics = rows(folder / "diagnostics.csv")
    steps = rows(folder / "time_steps.csv")
    times = [number(row["time"]) for row in diagnostics]
    energy = [number(row["kinetic_energy"]) for row in diagnostics]
    error = [number(row["velocity_error_l2"]) for row in diagnostics]
    defect = [number(row["convective_energy_defect"]) for row in diagnostics]
    energy0 = energy[0] if energy else math.nan
    ratio = [value / energy0 for value in energy]
    finite = all(math.isfinite(value) for value in ratio)
    reached = times[-1] if times else 0.0
    blow_up_time = next(
        (t for t, value in zip(times, ratio, strict=True) if not value <= BLOW_UP), None
    )
    completed = status == 0 and reached >= T_END - 1e-9
    if blow_up_time is not None:
        outcome = f"blow-up t={blow_up_time:.2f}"
    elif completed:
        outcome = "ok"
    elif status == 124:
        outcome = f"timeout t={reached:.2f}"
    elif status is None:
        outcome = "not run"
    else:
        outcome = f"failed t={reached:.2f}"
    dt = 0.01
    cumulative_defect = sum(value for value in defect if math.isfinite(value)) * dt
    outer = [int(row["outer_iterations"]) for row in steps if row["outer_iterations"]]
    parts = name.split("_")
    return {
        "run": name,
        "outcome": outcome,
        "exit": status,
        "wall_s": wall,
        "t_reached": reached,
        "energy_ratio_final": ratio[-1] if ratio else math.nan,
        "energy_ratio_max": max(ratio) if ratio and finite else math.nan,
        "error_final": error[-1] if error else math.nan,
        "cumulative_defect_over_E0": cumulative_defect / energy0 if energy else math.nan,
        "divergence_max": max(
            (number(row["divergence_l2"]) for row in diagnostics), default=math.nan
        ),
        "outer_max": max(outer, default=0),
        "family": parts[0],
    }


def setting_of(name: str) -> tuple[str, str]:
    """(form, setting suffix) for a VMS run name vms_<form>_<cells>_tauc<on|off>_c2_<c2>."""
    body = name.removeprefix("vms_")
    for form in sorted(FORMS, key=len, reverse=True):
        if body.startswith(form + "_"):
            return form, body[len(form) + 1 :]
    raise ValueError(name)


def cell(row: dict[str, object]) -> str:
    if row["outcome"] == "ok":
        return f"{row['error_final']:.3f} (E {100 * (1 - row['energy_ratio_final']):.1f}%)"
    return f"**{row['outcome']}**"


def main() -> None:
    log = {}
    for line in (HERE / "runs.log").read_text().splitlines():
        name, status, wall = line.split()
        log[name] = (int(status), int(wall))
    names = sorted(path.stem for path in (HERE / "configs").glob("*.toml"))
    table = [summarise(name, *log.get(name, (None, None))) for name in names]
    with (HERE / "screening.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    by_name = {row["run"]: row for row in table}
    lines = [
        "# Blow-up screening: Gresho, ν = 1e-6, Δt = 0.01, Crank–Nicolson, T = 3",
        "",
        "Cells: final L² velocity error (resolved energy lost). Blow-up: E/E(0) > 1.5.",
        "",
        "## Galerkin Taylor–Hood P2/P1, 48×48 (Li & Rui 2023 §4.2 resolution)",
        "",
        "| Form | Outcome | Error at end | E/E(0) max | ∫𝒩 dt / E(0) |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for form in FORMS[:-1]:
        row = by_name.get(f"galerkin_th_{form}_48")
        if row:
            lines.append(
                f"| {SHORT[form]} | {row['outcome']} | {row['error_final']:.3g} | "
                f"{row['energy_ratio_max']:.3g} | {row['cumulative_defect_over_E0']:.3g} |"
            )
    settings = sorted(
        {setting_of(name)[1] for name in names if name.startswith("vms_")},
        key=lambda s: (int(s.split("_")[0]), s.split("_")[1], int(s.split("_")[-1])),
    )
    lines += [
        "",
        "## P1/P1 dynamic OSS, nonlinear splitting",
        "",
        "| Mesh, τ_c, c2 | " + " | ".join(SHORT[form] for form in FORMS) + " |",
        "| --- |" + " --- |" * len(FORMS),
    ]
    for setting in settings:
        cells, tau, _, c2 = setting.split("_")
        label = f"{cells}², {tau.removeprefix('tauc')}, {c2}"
        entries = [
            cell(by_name[f"vms_{form}_{setting}"]) if f"vms_{form}_{setting}" in by_name else ""
            for form in FORMS
        ]
        lines.append(f"| {label} | " + " | ".join(entries) + " |")
    (HERE / "screening.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))  # noqa: T201 — command-line summary


if __name__ == "__main__":
    main()
