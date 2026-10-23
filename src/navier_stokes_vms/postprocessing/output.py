"""Results folder, resolved configuration, run information, VTK fields and CSV files (plan §9)."""

import csv
import json
import platform
import socket
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from importlib import metadata
from pathlib import Path

import firedrake as fd
from mpi4py import MPI
from petsc4py import PETSc

from navier_stokes_vms.config import Config, config_to_toml


class CsvWriter:
    """One CSV file with fixed columns, written and flushed row by row on rank 0 only."""

    def __init__(self, comm: MPI.Comm, path: Path, columns: Sequence[str]) -> None:
        self._columns = tuple(columns)
        self._stream = path.open("w", newline="", encoding="utf-8") if comm.rank == 0 else None
        self._write(self._columns)

    def write(self, **values: object) -> None:
        """Append one row; every column must be given, None is written as an empty field."""
        if set(values) != set(self._columns):
            raise ValueError(f"CSV row has keys {sorted(values)}, expected {list(self._columns)}")
        self._write(["" if values[name] is None else values[name] for name in self._columns])

    def close(self) -> None:
        """Close the file."""
        if self._stream is not None:
            self._stream.close()

    def _write(self, row: Sequence[object]) -> None:
        if self._stream is not None:
            csv.writer(self._stream).writerow(row)
            self._stream.flush()


class FieldWriter:
    """VTK series fields/velocity.pvd and fields/pressure.pvd, each with its own time."""

    def __init__(self, folder: Path) -> None:
        self._folder = folder / "fields"
        self._files: tuple[fd.VTKFile, fd.VTKFile] | None = None

    def write(
        self,
        velocity: fd.Function,
        pressure: fd.Function,
        velocity_time: float,
        pressure_time: float,
    ) -> None:
        """Append the velocity at velocity_time and the pressure at pressure_time."""
        if self._files is None:
            self._files = (
                fd.VTKFile(str(self._folder / "velocity.pvd")),
                fd.VTKFile(str(self._folder / "pressure.pvd")),
            )
        self._files[0].write(velocity, time=velocity_time)
        self._files[1].write(pressure, time=pressure_time)


def open_diagnostics_csv(comm: MPI.Comm, folder: Path, names: Sequence[str]) -> CsvWriter:
    """diagnostics.csv with field times and the expanded diagnostic columns."""
    return CsvWriter(
        comm, folder / "diagnostics.csv", ("step", "time", "pressure_time", "subscale_time", *names)
    )


@dataclass(frozen=True)
class RunInfo:
    """Contents of run_info.json; the end fields stay None until the run finishes."""

    git_commit: str | None
    git_dirty: bool | None
    firedrake_version: str
    petsc_version: str
    python_version: str
    mpi_size: int
    host: str
    start_time: str
    end_time: str | None = None
    wall_time_seconds: float | None = None
    exit_status: int | None = None


def create_results_folder(comm: MPI.Comm, root: Path, config: Config) -> Path:
    """Create root/<case>/<run_name>/ on rank 0, refusing to reuse an existing folder."""
    folder = root / config.case / config.run_name
    error = None
    if comm.rank == 0:
        try:
            folder.mkdir(parents=True)
        except FileExistsError:
            error = f"results folder {folder} already exists; remove it or set another run_name"
    error = comm.bcast(error, root=0)
    if error is not None:
        raise FileExistsError(error)
    return folder


def write_resolved_config(comm: MPI.Comm, folder: Path, config: Config) -> None:
    """Write config_resolved.toml on rank 0."""
    if comm.rank == 0:
        (folder / "config_resolved.toml").write_text(config_to_toml(config), encoding="utf-8")


def collect_run_info(comm: MPI.Comm) -> RunInfo:
    """Code and software versions, MPI size, host and start time, gathered on rank 0."""
    info = None
    if comm.rank == 0:
        commit, dirty = _git_state()
        info = RunInfo(
            git_commit=commit,
            git_dirty=dirty,
            firedrake_version=metadata.version("firedrake"),
            petsc_version=".".join(str(part) for part in PETSc.Sys.getVersion()),
            python_version=platform.python_version(),
            mpi_size=comm.size,
            host=socket.gethostname(),
            start_time=_now(),
        )
    return comm.bcast(info, root=0)


def finish_run_info(info: RunInfo, exit_status: int, wall_time_seconds: float) -> RunInfo:
    """Add the end time, wall time and process exit status."""
    return replace(
        info, end_time=_now(), wall_time_seconds=wall_time_seconds, exit_status=exit_status
    )


def write_run_info(comm: MPI.Comm, folder: Path, info: RunInfo) -> None:
    """Write run_info.json on rank 0."""
    if comm.rank == 0:
        text = json.dumps(asdict(info), indent=2)
        (folder / "run_info.json").write_text(text + "\n", encoding="utf-8")


def write_initial_state(comm: MPI.Comm, folder: Path, values: Mapping[str, float]) -> None:
    """Write initial_state.json on rank 0: domain size and moments of the initial velocity."""
    if comm.rank == 0:
        text = json.dumps(dict(values), indent=2)
        (folder / "initial_state.json").write_text(text + "\n", encoding="utf-8")


def _git_state() -> tuple[str | None, bool | None]:
    source = Path(__file__).resolve().parent
    try:
        commit = _git(source, "rev-parse", "HEAD").strip()
        status = _git(source, "status", "--porcelain")
    except (OSError, subprocess.CalledProcessError):
        return None, None
    return commit, bool(status.strip())


def _git(cwd: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", *arguments], cwd=cwd, capture_output=True, text=True, check=True
    )
    return result.stdout


def _now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")
