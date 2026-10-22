"""Convergence study: run_convergence_study.py BASE.toml --cells 8x8 16x16 | --steps 20 40 80."""

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

from mpi4py import MPI

from navier_stokes_vms.config import Config, config_to_toml, load_config, validate_config
from navier_stokes_vms.postprocessing.run_logging import configure_logging

RESULTS_ROOT = Path("results")
RUN_SIMULATION = Path(__file__).resolve().with_name("run_simulation.py")


def main(argv: list[str] | None = None) -> int:
    """Run the base configuration on every level and tabulate final diagnostics and rates."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="base TOML configuration")
    levels = parser.add_mutually_exclusive_group(required=True)
    levels.add_argument("--cells", nargs="+", help="mesh levels such as 8x8 16x16 32x32")
    levels.add_argument("--steps", nargs="+", type=int, help="numbers of time steps to t_end")
    arguments = parser.parse_args(argv)
    logger = configure_logging(MPI.COMM_WORLD, "INFO")
    try:
        base = load_config(arguments.config)
        studies = _level_configs(base, arguments.cells, arguments.steps)
    except ValueError as error:
        logger.error("Invalid study %s: %s", arguments.config, error)
        return 1
    folder = RESULTS_ROOT / base.case / f"{base.run_name}_study"
    try:
        folder.mkdir(parents=True)
    except FileExistsError:
        logger.error("study folder %s already exists; remove it or set another run_name", folder)
        return 1
    finals = []
    for label, refinement, config in studies:
        path = folder / f"{label}.toml"
        path.write_text(config_to_toml(config), encoding="utf-8")
        logger.info("Level %s: running %s", label, path)
        if subprocess.run([sys.executable, str(RUN_SIMULATION), str(path)], check=False).returncode:
            logger.error(
                "Level %s failed; see %s", label, RESULTS_ROOT / config.case / config.run_name
            )
            return 1
        diagnostics = RESULTS_ROOT / config.case / config.run_name / "diagnostics.csv"
        with diagnostics.open(encoding="utf-8") as stream:
            last = list(csv.DictReader(stream))[-1]
        names = tuple(
            name for name in last if name not in ("step", "time", "pressure_time", "subscale_time")
        )
        finals.append((label, refinement, {name: float(last[name]) for name in names}))
    _write_table(folder / "convergence.csv", names, finals)
    for label, refinement, values, rates in _with_rates(names, finals):
        cells = ", ".join(
            f"{name} {values[name]:.4e} ({'-' if rates[name] is None else f'{rates[name]:.2f}'})"
            for name in names
        )
        logger.info("%s (N = %d): %s", label, refinement, cells)
    logger.info("Table written to %s", folder / "convergence.csv")
    return 0


def _level_configs(
    base: Config, cells: list[str] | None, steps: list[int] | None
) -> list[tuple[str, int, Config]]:
    levels = []
    if cells is not None:
        for text in cells:
            counts = tuple(int(part) for part in text.split("x"))
            label = f"cells{text}"
            mesh = replace(base.mesh, cells=counts)
            levels.append((label, counts[0], replace(base, mesh=mesh)))
    else:
        for count in steps:
            label = f"steps{count}"
            timing = replace(base.time, dt=base.time.t_end / count)
            levels.append((label, count, replace(base, time=timing)))
    named = [
        (label, refinement, replace(config, run_name=f"{base.run_name}_{label}"))
        for label, refinement, config in levels
    ]
    for _, _, config in named:
        validate_config(config)
    return named


def _with_rates(
    names: tuple[str, ...], finals: list[tuple[str, int, dict[str, float]]]
) -> list[tuple[str, int, dict[str, float], dict[str, float | None]]]:
    rows = []
    for index, (label, refinement, values) in enumerate(finals):
        rates: dict[str, float | None] = dict.fromkeys(names)
        if index > 0:
            _, previous_refinement, previous = finals[index - 1]
            ratio = math.log(refinement / previous_refinement)
            for name in names:
                if previous[name] > 0 and values[name] > 0:
                    rates[name] = math.log(previous[name] / values[name]) / ratio
        rows.append((label, refinement, values, rates))
    return rows


def _write_table(
    path: Path, names: tuple[str, ...], finals: list[tuple[str, int, dict[str, float]]]
) -> None:
    columns = ["level", "refinement"] + [
        f"{name}{suffix}" for name in names for suffix in ("", "_rate")
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        for label, refinement, values, rates in _with_rates(names, finals):
            row = [label, refinement]
            for name in names:
                row += [values[name], "" if rates[name] is None else rates[name]]
            writer.writerow(row)


if __name__ == "__main__":
    sys.exit(main())
