"""Write the Phase 12b conservation matrix: python scripts/generate_conservation_configs.py."""

import argparse
import logging
import sys
from pathlib import Path
from types import MappingProxyType

from navier_stokes_vms.config import (
    Config,
    ConvectionOptions,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
    ViscosityOptions,
    config_to_toml,
    validate_config,
)

FOLDER = Path(__file__).resolve().parents[1] / "configs" / "chapter5" / "conservation"
FORMS = (
    "nonconservative",
    "conservative",
    "skew_symmetric_1",
    "skew_symmetric_2",
    "emac",
    "emac_vms",
)
LONG_FORMS = ("nonconservative", "conservative", "skew_symmetric_1", "emac", "emac_vms")
SPACES = ("asgs", "oss")
SPLITTINGS = ("linear", "nonlinear")
GRESHO_DIAGNOSTICS = (
    "kinetic_energy",
    "subscale_kinetic_energy",
    "linear_momentum",
    "angular_momentum",
    "energy_budget",
    "momentum_budget",
    "linear_momentum_defect",
    "linear_momentum_defect_predicted",
    "angular_momentum_defect",
    "angular_momentum_defect_predicted",
    "convective_energy_defect",
    "force",
    "region_velocity_l2",
    "region_subscale_l2",
    "velocity_error_l2",
    "pressure_error_l2",
    "divergence_l2",
)
# No angular quantity on the periodic domain: e_k × x is not periodic.
TAYLOR_GREEN_DIAGNOSTICS = (
    "kinetic_energy",
    "subscale_kinetic_energy",
    "linear_momentum",
    "energy_budget",
    "momentum_budget",
    "linear_momentum_defect",
    "linear_momentum_defect_predicted",
    "convective_energy_defect",
    "velocity_error_l2",
)
SOLVER = SolverOptions(
    linearization="newton",
    newton_tau="frozen",
    picard_iterations_before_newton=0,
    inner_rtol=1e-12,
    inner_max_iterations=200,
    outer_rtol=1e-11,
    outer_max_iterations=100,
    subscale_rtol=1e-13,
)


def run_config(
    run_name: str,
    case: str,
    cells: int,
    dt: float,
    t_end: float,
    form: str,
    space: str,
    splitting: str,
    diagnostics: tuple[str, ...],
    *,
    viscosity: float = 1e-6,
) -> Config:
    """One run with the common P1/P1 stress, dynamic-subscale, Crank–Nicolson settings."""
    return Config(
        case=case,
        run_name=run_name,
        case_parameters=MappingProxyType({"viscosity": viscosity}),
        mesh=MeshOptions(cells=(cells, cells)),
        discretization=DiscretizationOptions(velocity_degree=1, pressure_degree=1),
        time=TimeOptions(scheme="crank_nicolson", dt=dt, t_end=t_end, startup="lower_order"),
        stabilization=StabilizationOptions(
            method="vms",
            subscale_space=space,
            subscale_dynamics="dynamic",
            scale_splitting=splitting,
            pressure_subscale=True,
        ),
        convection=ConvectionOptions(form=form),
        viscosity=ViscosityOptions(form="stress"),
        solver=SOLVER,
        output=OutputOptions(diagnostics=diagnostics),
    )


def conservation_configs() -> list[Config]:
    """The 55 runs of plan Phase 12b: 12 OSS + 12 viscous ASGS + 5 long + 1 + 1 refined + 24."""
    matrix = [
        (form, space, splitting) for form in FORMS for space in SPACES for splitting in SPLITTINGS
    ]
    configs = [
        run_config(
            f"{prefix}_{form}_{space}_{splitting}",
            "gresho",
            32,
            0.005,
            0.5,
            form,
            space,
            splitting,
            GRESHO_DIAGNOSTICS,
            viscosity=viscosity,
        )
        for prefix, space, viscosity in (("gresho", "oss", 1e-6), ("gresho_visc", "asgs", 1e-3))
        for form in FORMS
        for splitting in SPLITTINGS
    ]
    configs += [
        run_config(
            f"gresho_long_{form}",
            "gresho",
            32,
            0.01,
            10.0,
            form,
            "oss",
            "nonlinear",
            GRESHO_DIAGNOSTICS,
        )
        for form in LONG_FORMS
    ]
    configs += [
        run_config(
            f"{prefix}_refined_emac_vms_{space}",
            "gresho",
            64,
            0.0025,
            0.5,
            "emac_vms",
            space,
            "nonlinear",
            GRESHO_DIAGNOSTICS,
            viscosity=viscosity,
        )
        for prefix, space, viscosity in (("gresho", "oss", 1e-6), ("gresho_visc", "asgs", 1e-3))
    ]
    configs += [
        run_config(
            f"taylor_green_{form}_{space}_{splitting}",
            "taylor_green_2d",
            32,
            0.005,
            0.5,
            form,
            space,
            splitting,
            TAYLOR_GREEN_DIAGNOSTICS,
            viscosity=0.01,
        )
        for form, space, splitting in matrix
    ]
    for config in configs:
        validate_config(config)
    return configs


def write_configs(folder: Path) -> list[Path]:
    """Write the matrix and remove the 13 retired near-inviscid ASGS configurations."""
    folder.mkdir(parents=True, exist_ok=True)
    paths = []
    for config in conservation_configs():
        path = folder / f"{config.run_name}.toml"
        path.write_text(config_to_toml(config), encoding="utf-8")
        paths.append(path)
    retired = [f"gresho_{form}_asgs_{splitting}" for form in FORMS for splitting in SPLITTINGS]
    retired.append("gresho_refined_emac_vms_asgs")
    for name in retired:
        (folder / f"{name}.toml").unlink(missing_ok=True)
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type=Path, default=FOLDER, help="output directory")
    arguments = parser.parse_args(argv)
    # The format of run logs, without importing Firedrake for a file-only task.
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(Path(__file__).stem)
    paths = write_configs(arguments.folder)
    logger.info("Wrote %d configurations to %s", len(paths), arguments.folder)
    return 0


if __name__ == "__main__":
    sys.exit(main())
