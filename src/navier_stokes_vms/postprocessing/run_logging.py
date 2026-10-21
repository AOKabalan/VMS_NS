"""Console and run.log logging on MPI rank 0 (plan §9)."""

import logging
from pathlib import Path

from mpi4py import MPI

from navier_stokes_vms.config import Config, config_as_dict
from navier_stokes_vms.postprocessing.output import RunInfo

LOGGER_NAME = "navier_stokes_vms"

_FORMAT = "%(asctime)s %(levelname)-7s %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def configure_logging(
    comm: MPI.Comm, console_level: str, log_file: Path | None = None
) -> logging.Logger:
    """Route the package log to the console and, if given, to run.log; other ranks stay silent."""
    logger = logging.getLogger(LOGGER_NAME)
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    if comm.rank != 0:
        logger.addHandler(logging.NullHandler())
        return logger
    formatter = logging.Formatter(_FORMAT, _DATE_FORMAT)
    console = logging.StreamHandler()
    console.setLevel(console_level)
    console.setFormatter(formatter)
    logger.addHandler(console)
    if log_file is not None:
        run_log = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        run_log.setLevel(logging.INFO)
        run_log.setFormatter(formatter)
        logger.addHandler(run_log)
    return logger


def log_setup_summary(logger: logging.Logger, config: Config, folder: Path, info: RunInfo) -> None:
    """Log the run identity, every configuration section and the software environment."""
    logger.info("Case %s, run %s, results in %s", config.case, config.run_name, folder)
    for section, values in config_as_dict(config).items():
        if isinstance(values, dict):
            settings = ", ".join(f"{key} = {value}" for key, value in values.items())
            logger.info("[%s] %s", section, settings or "(none)")
    state = "unknown" if info.git_dirty is None else ("dirty" if info.git_dirty else "clean")
    logger.info(
        "Git %s (%s), Firedrake %s, PETSc %s, Python %s, %d MPI rank(s) on %s",
        info.git_commit or "unavailable",
        state,
        info.firedrake_version,
        info.petsc_version,
        info.python_version,
        info.mpi_size,
        info.host,
    )
