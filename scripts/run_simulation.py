"""Run one simulation: python scripts/run_simulation.py configs/<file>.toml."""

import argparse
import sys
import time
from pathlib import Path

from mpi4py import MPI

from navier_stokes_vms.config import load_config
from navier_stokes_vms.postprocessing.output import (
    collect_run_info,
    create_results_folder,
    finish_run_info,
    write_resolved_config,
    write_run_info,
)
from navier_stokes_vms.postprocessing.run_logging import configure_logging, log_setup_summary

RESULTS_ROOT = Path("results")


def main(argv: list[str] | None = None) -> int:
    """Validate the configuration, prepare results/<case>/<run_name>/ and log the setup."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="TOML configuration file")
    arguments = parser.parse_args(argv)
    clock = time.perf_counter()
    comm = MPI.COMM_WORLD
    logger = configure_logging(comm, "INFO")
    info = collect_run_info(comm)
    try:
        config = load_config(arguments.config)
        folder = create_results_folder(comm, RESULTS_ROOT, config)
    except ValueError as error:
        logger.error("Invalid configuration %s: %s", arguments.config, error)
        return 1
    except OSError as error:
        logger.error("%s", error)
        return 1
    logger = configure_logging(comm, config.output.console_level, folder / "run.log")
    write_resolved_config(comm, folder, config)
    write_run_info(comm, folder, info)
    exit_status = 1
    try:
        log_setup_summary(logger, config, folder, info)
        logger.info("No solver is implemented yet; stopping after setup.")
        exit_status = 0
    except Exception:
        logger.exception("Run failed")
    finally:
        wall_time = time.perf_counter() - clock
        write_run_info(comm, folder, finish_run_info(info, exit_status, wall_time))
    return exit_status


if __name__ == "__main__":
    sys.exit(main())
