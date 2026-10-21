"""Results folder, resolved configuration and run information (plan §9)."""

import json
import platform
import socket
import subprocess
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from importlib import metadata
from pathlib import Path

from mpi4py import MPI
from petsc4py import PETSc

from navier_stokes_vms.config import Config, config_to_toml


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
