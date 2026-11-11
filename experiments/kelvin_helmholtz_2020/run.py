"""Run one KH configuration with enstrophy and Figure 4 vorticity snapshots (MPI-safe)."""

import argparse
import json
import logging
import math
import resource
import sys
import time
from dataclasses import replace
from pathlib import Path

import firedrake as fd
import numpy as np
from firedrake.petsc import PETSc
from kh_profiling import PhaseProfiler
from mpi4py import MPI

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import Config, load_config, validate_config
from navier_stokes_vms.postprocessing.diagnostics import check_diagnostics_supported
from navier_stokes_vms.postprocessing.output import (
    CsvWriter,
    collect_run_info,
    create_results_folder,
    finish_run_info,
    write_resolved_config,
    write_run_info,
)
from navier_stokes_vms.postprocessing.run_logging import configure_logging, log_setup_summary
from navier_stokes_vms.problem.problem_definition import Problem
from navier_stokes_vms.solver.time_stepper import TimeStepper

FIGURE_TIMES = (1.0, 4.5, 5.0, 5.5, 6.0, 6.5, 10.0)


class SnapshotStepper(TimeStepper):
    """Keep DG1 vorticity on each triangle, without smoothing across element boundaries."""

    def __init__(
        self,
        cfg: Config,
        problem: Problem,
        folder: Path,
        logger: logging.Logger,
        snapshot_times: list[float],
    ) -> None:
        super().__init__(cfg, problem, folder, logger)
        self.comm = problem.mesh.comm
        self.snapshot_steps = {round(t / cfg.time.dt) for t in snapshot_times}
        self.snapshot_steps.add(round(cfg.time.t_end / cfg.time.dt))
        mesh = problem.mesh
        self.measure = fd.dx(domain=mesh, degree=max(6, 3 * cfg.discretization.velocity_degree))
        self.wall_measure = fd.ds((3, 4), domain=mesh, degree=6)
        self.curl = self.velocity[1].dx(0) - self.velocity[0].dx(1)
        # The default DG variant has interior nodes. Equispaced DG1 puts its
        # nodes at the triangle vertices, as required by the saved connectivity.
        self.omega = fd.Function(fd.FunctionSpace(mesh, "DG", 1, variant="equispaced"))
        coordinates = fd.Function(
            fd.VectorFunctionSpace(mesh, "DG", 1, variant="equispaced")
        ).interpolate(fd.SpatialCoordinate(mesh))
        self.cell_nodes = self.omega.function_space().cell_node_map().values[: mesh.cell_set.size]
        vertices = coordinates.dat.data_ro_with_halos[self.cell_nodes].copy()
        parts = self.comm.gather(vertices, root=0)
        if self.comm.rank == 0:
            np.savez_compressed(folder / "mesh.npz", vertices=np.concatenate(parts))
        self.enstrophy = CsvWriter(
            self.comm,
            folder / "enstrophy.csv",
            ("step", "time", "enstrophy", "vorticity_max_abs", "normal_velocity_wall_l2"),
        )
        self.started = time.perf_counter()

    def _initialise(self) -> None:
        super()._initialise()
        self.record(0, 0.0)

    def record(self, step: int, t: float) -> bool:
        """Write durable diagnostics each step and snapshots at the selected endpoints."""
        self.omega.interpolate(self.curl)
        max_local = float(np.max(np.abs(self.omega.dat.data_ro), initial=0.0))
        self.enstrophy.write(
            step=step,
            time=t,
            enstrophy=float(fd.assemble(self.curl**2 / 2 * self.measure)),
            vorticity_max_abs=self.comm.allreduce(max_local, op=MPI.MAX),
            normal_velocity_wall_l2=math.sqrt(
                abs(float(fd.assemble(self.velocity[1] ** 2 * self.wall_measure)))
            ),
        )
        if step == 0 or step in self.snapshot_steps:
            local = self.omega.dat.data_ro_with_halos[self.cell_nodes].copy()
            parts = self.comm.gather(local, root=0)
            if self.comm.rank == 0:
                target = self._folder / f"snapshot_{step:08d}.npz"
                temporary = self._folder / f".snapshot_{step:08d}.tmp.npz"
                np.savez_compressed(
                    temporary,
                    time=t,
                    step=step,
                    vorticity=np.concatenate(parts),
                )
                temporary.replace(target)
        if step and (
            step == 1
            or step % 100 == 0
            or step == round(self._config.time.t_end / self._config.time.dt)
        ):
            if self.comm.rank == 0:
                elapsed = time.perf_counter() - self.started
                sys.stdout.write(
                    f"{self._config.run_name}: step {step}, t={t:g}, "
                    f"elapsed {elapsed / 60:.1f} min\n"
                )
                sys.stdout.flush()
        return False


def main() -> int:
    """Use the production solver and preserve normal run provenance and failure status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output-root", type=Path, default=Path("results"))
    parser.add_argument(
        "--max-steps", type=int, help="timing/smoke prefix, saved under a distinct name"
    )
    parser.add_argument("--snapshot-times", nargs="+", type=float, default=list(FIGURE_TIMES))
    parser.add_argument("--petsc-profile", action="store_true", help="save PETSc event timings")
    parser.add_argument(
        "--profile-phases", action="store_true", help="split startup and measured solver phases"
    )
    parser.add_argument(
        "--save-final-state",
        action="store_true",
        help="save owned FE arrays per rank for comparison",
    )
    args = parser.parse_args()
    cfg = load_config(args.config)
    if cfg.case != "kelvin_helmholtz":
        parser.error("this driver requires case = 'kelvin_helmholtz'")
    if cfg.discretization.velocity_degree not in (1, 2):
        parser.error("DG1 snapshots support velocity degree 1 or 2")
    if args.save_final_state and (
        cfg.discretization.velocity_degree != cfg.discretization.pressure_degree
    ):
        parser.error("coordinate-based final-state exports require equal-order velocity/pressure")
    if args.max_steps is not None:
        if not 0 < args.max_steps <= round(cfg.time.t_end / cfg.time.dt):
            parser.error("--max-steps must be positive and no longer than the configured run")
        cfg = replace(
            cfg,
            run_name=cfg.run_name + f"_prefix{args.max_steps}",
            time=replace(cfg.time, t_end=args.max_steps * cfg.time.dt),
        )
    if any(not math.isfinite(t) or t < 0 for t in args.snapshot_times):
        parser.error("snapshot times must be finite and nonnegative")
    times = [t for t in args.snapshot_times if t <= cfg.time.t_end]
    if any(not math.isclose(t / cfg.time.dt, round(t / cfg.time.dt), abs_tol=1e-8) for t in times):
        parser.error("snapshot times must coincide with time-step endpoints")
    validate_config(cfg)
    if args.petsc_profile or args.profile_phases:
        PETSc.Log.begin()
    comm = MPI.COMM_WORLD
    logger = configure_logging(comm, "INFO")
    folder = create_results_folder(comm, args.output_root, cfg)
    info = collect_run_info(comm)
    write_resolved_config(comm, folder, cfg)
    write_run_info(comm, folder, info)
    logger = configure_logging(comm, cfg.output.console_level, folder / "run.log")
    if comm.rank == 0:
        sys.stdout.write(
            f"{cfg.run_name}: {round(cfg.time.t_end / cfg.time.dt)} steps, "
            f"T={cfg.time.t_end:g}, {comm.size} MPI rank(s); output {folder}\n"
        )
        sys.stdout.flush()
    started = time.perf_counter()
    status, stepper = 1, None
    profiler = None
    try:
        problem = build_problem(cfg)
        check_diagnostics_supported(cfg.output.diagnostics, problem)
        log_setup_summary(logger, cfg, folder, info)
        stepper = SnapshotStepper(cfg, problem, folder, logger, times)
        if args.profile_phases:
            profiler = PhaseProfiler()
            profiler.install(stepper)
        if comm.rank == 0:
            metadata = dict(
                source="Olshanskii & Rebholz 2020, section 4.3",
                source_url="https://arxiv.org/abs/2002.01416",
                viscosity=problem.viscosity,
                figure_times=list(FIGURE_TIMES),
                requested_snapshot_times=times,
                velocity_dofs=stepper.velocity.function_space().dim(),
                mixed_dofs=stepper._solution.function_space().dim(),
                vorticity_representation="discontinuous linear, triangle-local nodes",
                prefix=args.max_steps is not None,
            )
            (folder / "experiment.json").write_text(json.dumps(metadata, indent=2) + "\n")
        stepper.run(stepper.record)
        if args.save_final_state:
            fields = {
                f"field_{i}": field.dat.data_ro.copy()
                for i, field in enumerate(stepper._solution.subfunctions)
            }
            # Rank partitions may differ between launches. Coordinates allow
            # comparison by physical location rather than local array order.
            x = fd.SpatialCoordinate(problem.mesh)
            coordinates = fd.Function(stepper.velocity.function_space()).interpolate(x)
            fields["fe_coordinates"] = coordinates.dat.data_ro.copy()
            if stepper._iterate is not None:
                fields["subscale"] = stepper._iterate.value.dat.data_ro.copy()
                subspace = stepper._iterate.value.function_space()
                subcoordinates = fd.Function(subspace).interpolate(x)
                centres = fd.Function(fd.VectorFunctionSpace(problem.mesh, "DG", 0)).interpolate(x)
                subcentres = fd.Function(subspace).interpolate(centres)
                fields["subscale_coordinates"] = subcoordinates.dat.data_ro.copy()
                fields["subscale_cell_centres"] = subcentres.dat.data_ro.copy()
            np.savez_compressed(folder / f"final_state_rank{comm.rank:04d}.npz", **fields)
        status = 0
    except Exception:
        logger.exception("Kelvin–Helmholtz run stopped; completed snapshots and rows are retained")
    finally:
        if profiler is not None:
            profiler.close(comm, folder / "phase_timings.json")
        if stepper is not None:
            stepper.enstrophy.close()
        write_run_info(comm, folder, finish_run_info(info, status, time.perf_counter() - started))
        if args.petsc_profile or args.profile_phases:
            peaks = comm.gather(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, root=0)
            if comm.rank == 0:
                (folder / "memory_peaks.json").write_text(
                    json.dumps(
                        {
                            "per_rank_max_rss_kib": peaks,
                            "note": "Linux high-water marks; rank peaks need not be simultaneous",
                        },
                        indent=2,
                    )
                    + "\n"
                )
            if stepper is not None and stepper._step_solver._inner._newton is not None:
                viewer = PETSc.Viewer().createASCII(str(folder / "solver_view.txt"), comm=comm)
                try:
                    stepper._step_solver._inner._newton.snes.view(viewer)
                finally:
                    viewer.destroy()
            viewer = PETSc.Viewer().createASCII(str(folder / "petsc_profile.txt"), comm=comm)
            try:
                PETSc.Log.view(viewer)
            finally:
                viewer.destroy()
    return status


if __name__ == "__main__":
    sys.exit(main())
