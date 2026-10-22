"""Run a validation case and write comparisons, or analyse its existing sampled outputs."""

import argparse
import csv
import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
from mpi4py import MPI

from navier_stokes_vms.cases import (
    build_problem,
    cylinder_dfg_2d,
    lattice_vortex,
    lid_driven_cavity,
)
from navier_stokes_vms.config import Config, check_diagnostics_known, load_config
from navier_stokes_vms.discretization.quadrature import volume_measure
from navier_stokes_vms.postprocessing.diagnostics import DIAGNOSTICS, check_diagnostics_supported
from navier_stokes_vms.postprocessing.output import (
    CsvWriter,
    collect_run_info,
    create_results_folder,
    finish_run_info,
    write_resolved_config,
    write_run_info,
)
from navier_stokes_vms.postprocessing.run_logging import configure_logging, log_setup_summary
from navier_stokes_vms.solver.time_stepper import TimeStepper


def run(config: Config, root: Path) -> int:
    """Integrate, retaining the small samples needed for analysis without reloading FE fields."""
    comm = MPI.COMM_WORLD
    clock = time.perf_counter()
    logger = configure_logging(comm, config.output.console_level)
    check_diagnostics_known(config, DIAGNOSTICS)
    problem = build_problem(config)
    check_diagnostics_supported(config.output.diagnostics, problem)
    folder = create_results_folder(comm, root, config)
    logger = configure_logging(comm, config.output.console_level, folder / "run.log")
    info = collect_run_info(comm)
    write_resolved_config(comm, folder, config)
    write_run_info(comm, folder, info)
    log_setup_summary(logger, config, folder, info)
    monitor = None
    status = 1
    try:
        stepper = TimeStepper(config, problem, folder, logger)
        after_step = None
        if config.case == "lid_driven_cavity":
            previous = fd.Function(stepper.velocity.function_space())
            dx = volume_measure(problem.mesh, 3 * config.discretization.velocity_degree)
            monitor = CsvWriter(comm, folder / "steadiness.csv", ("time", "relative_change"))
            steady_steps = 0

            def after_step(_step: int, time: float) -> bool:
                nonlocal steady_steps
                difference = stepper.velocity - previous
                change = math.sqrt(fd.assemble(fd.inner(difference, difference) * dx))
                norm = math.sqrt(fd.assemble(fd.inner(stepper.velocity, stepper.velocity) * dx))
                relative = change / max(norm, 1e-30)
                monitor.write(time=time, relative_change=relative)
                previous.assign(stepper.velocity)
                steady_steps = (
                    steady_steps + 1 if relative < lid_driven_cavity.STEADY_TOLERANCE else 0
                )
                return steady_steps >= 3

        elif config.case == "cylinder_dfg_2d":
            evaluator = cylinder_dfg_2d.outlet_evaluator(problem)
            monitor = CsvWriter(comm, folder / "outlet.csv", ("time", "minimum_normal_velocity"))

            def after_step(_step: int, time: float) -> bool:
                monitor.write(
                    time=time,
                    minimum_normal_velocity=cylinder_dfg_2d.outlet_minimum(
                        stepper.velocity, evaluator
                    ),
                )
                return False

        stepper.run(after_step)
        if config.case == "lid_driven_cavity":
            reynolds = config.case_parameters.get("reynolds", 100)
            samples = lid_driven_cavity.sample_profiles(stepper.velocity, reynolds)
            writer = CsvWriter(
                comm, folder / "centreline_samples.csv", ("component", "coordinate", "value")
            )
            try:
                for name, profile in lid_driven_cavity.reference_profiles(reynolds).items():
                    for coordinate, value in zip(profile[:, 0], samples[name], strict=True):
                        writer.write(component=name, coordinate=coordinate, value=value)
            finally:
                writer.close()
        status = 0
    except Exception:
        logger.exception("Validation simulation failed")
    finally:
        if monitor is not None:
            monitor.close()
        write_run_info(comm, folder, finish_run_info(info, status, time.perf_counter() - clock))
    return status


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a completed CSV, preserving strings for component names and optional columns."""
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    """Write an analysis table on the calling rank (analysis runs on rank zero only)."""
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyse(config: Config, folder: Path) -> None:
    """Case-specific comparisons; validation.csv is the completion marker, written last."""
    info = json.loads((folder / "run_info.json").read_text())
    if info["exit_status"] != 0:
        raise ValueError("Cannot analyse an incomplete or failed simulation")
    config = load_config(folder / "config_resolved.toml")
    data = read_rows(folder / "diagnostics.csv")
    times = np.asarray([float(row["time"]) for row in data])
    references = {}
    notes = [
        f"# {config.case}: {config.run_name}",
        f"Completed {len(times)} steps to t={times[-1]:g}.",
    ]
    if config.case == "lid_driven_cavity":
        reynolds = config.case_parameters.get("reynolds", 100)
        raw = read_rows(folder / "centreline_samples.csv")
        samples = {
            name: np.asarray([float(row["value"]) for row in raw if row["component"] == name])
            for name in ("u", "v")
        }
        metrics = lid_driven_cavity.compare(samples, reynolds)
        changes = [float(row["relative_change"]) for row in read_rows(folder / "steadiness.csv")]
        metrics["last_time_level_relative_change"] = changes[-1]
        steady = len(changes) >= 3 and max(changes[-3:]) < lid_driven_cavity.STEADY_TOLERANCE
        metrics["steady"] = float(steady)
        rows = []
        for name, profile in lid_driven_cavity.reference_profiles(reynolds).items():
            for (coordinate, reference), value in zip(profile, samples[name], strict=True):
                rows.append(
                    dict(
                        component=name,
                        coordinate=coordinate,
                        value=value,
                        reference=reference,
                        absolute_deviation=abs(value - reference),
                    )
                )
        write_table(folder / "profile_comparison.csv", rows)
        notes += [
            f"Steady between time levels: {steady}; criterion is relative L² change "
            "<1e-8 for three consecutive steps.",
            "Both top corner nodes belong to the stationary walls. Ghia’s finite-difference "
            "accuracy limits the comparison near the walls at Re=1000.",
            "Reference: Ghia, Ghia & Shin (1982), Tables I–II. See profile_comparison.csv "
            "for all 17 stations per profile.",
        ]
    elif config.case == "lattice_vortex":
        viscosity = config.case_parameters.get("viscosity", 7e-5)
        energy = np.asarray([float(row["kinetic_energy"]) for row in data])
        exact = 0.25 * np.exp(-16 * viscosity * math.pi**2 * times)
        metrics = lattice_vortex.compare(times, energy, viscosity)
        write_table(
            folder / "energy_comparison.csv",
            [
                dict(time=t, energy=e, reference=r, relative_error=abs(e - r) / r)
                for t, e, r in zip(times, energy, exact, strict=True)
            ],
        )
        notes += [
            "Reference: thesis §6.2, equations (6.3)–(6.5); E(t)=¼ exp(−16νπ²t). "
            "See energy_comparison.csv for the full history."
        ]
    else:
        force_times = np.asarray([float(row["pressure_time"]) for row in data])
        drag = np.asarray([float(row["force_cylinder_x"]) for row in data]) * 20
        lift = np.asarray([float(row["force_cylinder_y"]) for row in data]) * 20
        minima = cylinder_dfg_2d.lift_minima(force_times, lift)
        metrics = {
            "complete_cycles": float(max(0, len(minima) - 1)),
            "outlet_minimum_normal_velocity": min(
                float(row["minimum_normal_velocity"]) for row in read_rows(folder / "outlet.csv")
            ),
        }
        references = json.loads((cylinder_dfg_2d.REFERENCES / "featflow.json").read_text())
        if len(minima) >= 2:
            metrics.update(cylinder_dfg_2d.compare(force_times, drag, lift))
            cycle = cylinder_dfg_2d.last_cycle(force_times, drag, lift)
            reference = cylinder_dfg_2d.reference_curve()
            metrics.update(cylinder_dfg_2d.compare_curves(cycle, reference))
            end = min(cycle[-1, 0], reference[-1, 0])
            aligned = np.unique(
                np.r_[cycle[cycle[:, 0] <= end, 0], reference[reference[:, 0] <= end, 0], end]
            )
            write_table(
                folder / "cycle_comparison.csv",
                [
                    dict(
                        elapsed_time=t,
                        cd=np.interp(t, cycle[:, 0], cycle[:, 1]),
                        cl=np.interp(t, cycle[:, 0], cycle[:, 2]),
                        reference_cd=np.interp(t, reference[:, 0], reference[:, 1]),
                        reference_cl=np.interp(t, reference[:, 0], reference[:, 2]),
                    )
                    for t in aligned
                ],
            )
        else:
            notes.append(
                "Insufficient lift minima for a complete shedding cycle: cycle statistics "
                "are unavailable (expected for smoke/probe runs)."
            )
        if len(minima) < 4:
            notes.append("Fewer than three complete cycles: periodicity change is unavailable.")
        intervals = json.loads((cylinder_dfg_2d.REFERENCES / "dfg_1996.json").read_text())
        notes += [
            "Primary scalar reference: FEATFLOW finest listed level. Newer 20141124 level-6 "
            "curves are compared over the common elapsed-time interval, aligned at a lift "
            "minimum without stretching time.",
            "FEATFLOW mean=(max+min)/2 and amplitude=max−min. St=D/(Ū×cycle length), "
            "D=0.1, Ū=1. Minima use local quadratic interpolation; "
            "cycle extrema use sampled values.",
            "The outlet minimum includes interior extrema of every P1/P2 boundary edge, "
            "evaluated collectively. Negative values indicate backflow.",
            f"1996 contextual intervals: {intervals}",
            "Reference: https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark2_re100.html",
        ]
    rows = []
    for name, value in metrics.items():
        ref = references.get(name)
        rows.append(
            dict(
                metric=name,
                value=value,
                reference=ref,
                absolute_deviation=None if ref is None else abs(value - ref),
                relative_deviation=None if ref is None or ref == 0 else abs(value - ref) / abs(ref),
            )
        )
    table = ["| Metric | Value | Reference |", "| --- | ---: | ---: |"]
    table += [
        f"| {row['metric']} | {row['value']:.8g} | "
        f"{row['reference'] if row['reference'] is not None else '—'} |"
        for row in rows
    ]
    notes.append("\n".join(table))
    (folder / "validation.md").write_text("\n\n".join(notes) + "\n")
    temporary = folder / "validation.csv.tmp"
    write_table(temporary, rows)
    temporary.replace(folder / "validation.csv")


def main(argv: list[str] | None = None) -> int:
    """Run/observe or reanalyse; smoke and step limits preserve mesh, dt and solver settings."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--analyse-only", action="store_true")
    mode.add_argument("--run-only", action="store_true")
    parser.add_argument("--results-root", type=Path, default=Path("results"))
    limit = parser.add_mutually_exclusive_group()
    limit.add_argument("--smoke", action="store_true")
    limit.add_argument("--steps", type=int, help="Short runtime probe with at most this many steps")
    args = parser.parse_args(argv)
    comm = MPI.COMM_WORLD
    logger = configure_logging(comm, "INFO")
    try:
        config = load_config(args.config)
        if config.case not in ("lid_driven_cavity", "lattice_vortex", "cylinder_dfg_2d"):
            raise ValueError("run_validation supports only the three Phase 10 validation cases")
        if args.steps is not None and args.steps <= 0:
            raise ValueError("--steps must be positive")
        if args.smoke or args.steps:
            steps = 3 if args.smoke else args.steps
            config = replace(
                config,
                time=replace(config.time, t_end=min(config.time.t_end, steps * config.time.dt)),
            )
        if not args.analyse_only and run(config, args.results_root):
            return 1
        status = 0
        if not args.run_only and comm.rank == 0:
            try:
                analyse(config, args.results_root / config.case / config.run_name)
            except Exception:
                logger.exception("Validation analysis failed")
                status = 1
        return comm.bcast(status, root=0)
    except (OSError, ValueError) as error:
        logger.error("%s", error)
        return 1


if __name__ == "__main__":
    sys.exit(main())
