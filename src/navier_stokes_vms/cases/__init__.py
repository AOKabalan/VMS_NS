"""Registry of cases, case name → build function, and the checked problem build (plan §8)."""

from collections.abc import Callable

from navier_stokes_vms.cases import (
    ethier_steinman_3d,
    kovasznay,
    manufactured_spatial,
    manufactured_temporal,
    taylor_green_2d,
)
from navier_stokes_vms.config import Config, check_startup_supported
from navier_stokes_vms.problem.problem_definition import Problem, check_problem

CASES: dict[str, Callable[[Config], Problem]] = {
    "ethier_steinman_3d": ethier_steinman_3d.build,
    "kovasznay": kovasznay.build,
    "manufactured_spatial": manufactured_spatial.build,
    "manufactured_temporal": manufactured_temporal.build,
    "taylor_green_2d": taylor_green_2d.build,
}


def build_problem(config: Config) -> Problem:
    """Build the configured case, check it, and apply plan §7 rule 3, which needs the case."""
    if config.case not in CASES:
        raise ValueError(f"case = {config.case!r}: unknown; available: {', '.join(sorted(CASES))}")
    problem = CASES[config.case](config)
    check_problem(problem)
    check_startup_supported(config, has_exact_solution=problem.exact_velocity is not None)
    return problem
