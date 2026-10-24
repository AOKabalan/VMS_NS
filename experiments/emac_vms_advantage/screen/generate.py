"""Blow-up screening: long-Gresho settings shortened to T = 3, writing configs to ./configs."""

from dataclasses import replace
from pathlib import Path
from types import MappingProxyType

from navier_stokes_vms.cases import build_problem
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
from navier_stokes_vms.postprocessing.diagnostics import check_diagnostics_supported

FOLDER = Path(__file__).parent / "configs"
FORMS = (
    "nonconservative",
    "conservative",
    "skew_symmetric_1",
    "skew_symmetric_2",
    "emac",
    "emac_vms",
)
SOLVER = SolverOptions(
    linearization="newton",
    newton_tau="frozen",
    inner_rtol=1e-10,
    inner_max_iterations=200,
    outer_rtol=1e-9,
    outer_max_iterations=100,
    subscale_rtol=1e-11,
)
VMS_DIAGNOSTICS = (
    "kinetic_energy",
    "subscale_kinetic_energy",
    "linear_momentum",
    "angular_momentum",
    "energy_budget",
    "convective_energy_defect",
    "velocity_error_l2",
    "divergence_l2",
)
GALERKIN_DIAGNOSTICS = (
    "kinetic_energy",
    "linear_momentum",
    "angular_momentum",
    "convective_energy_defect",
    "velocity_error_l2",
    "divergence_l2",
)


def config(
    name: str,
    form: str,
    cells: int,
    stabilization: StabilizationOptions,
    degrees: tuple[int, int],
    diagnostics: tuple[str, ...],
) -> Config:
    return Config(
        case="gresho",
        run_name=name,
        case_parameters=MappingProxyType({"viscosity": 1e-6}),
        mesh=MeshOptions(cells=(cells, cells)),
        discretization=DiscretizationOptions(*degrees),
        time=TimeOptions(scheme="crank_nicolson", dt=0.01, t_end=3.0, startup="lower_order"),
        stabilization=stabilization,
        convection=ConvectionOptions(form=form),
        viscosity=ViscosityOptions(form="stress"),
        solver=SOLVER,
        output=OutputOptions(vtk_every=0, diagnostics=diagnostics),
    )


def screening() -> list[Config]:
    configs = []
    for form in FORMS:
        for cells in (16, 32):
            for tau_c in (True, False):
                for c2 in (2.0, 32.0, 512.0):
                    stabilization = StabilizationOptions(
                        method="vms",
                        subscale_space="oss",
                        subscale_dynamics="dynamic",
                        scale_splitting="nonlinear",
                        c2=c2,
                        pressure_subscale=tau_c,
                    )
                    name = f"vms_{form}_{cells}_tauc{'on' if tau_c else 'off'}_c2_{int(c2)}"
                    configs.append(
                        config(name, form, cells, stabilization, (1, 1), VMS_DIAGNOSTICS)
                    )
    # Anchor: Galerkin Taylor–Hood at the Li & Rui (2023) §4.2 resolution and step.
    for form in FORMS[:-1]:
        configs.append(
            config(
                f"galerkin_th_{form}_48",
                form,
                48,
                StabilizationOptions(method="none"),
                (2, 1),
                GALERKIN_DIAGNOSTICS,
            )
        )
    for item in configs:
        validate_config(item)
        problem = build_problem(replace(item, mesh=MeshOptions(cells=(2, 2))))
        check_diagnostics_supported(item.output.diagnostics, problem)
    return configs


if __name__ == "__main__":
    FOLDER.mkdir(exist_ok=True)
    items = screening()
    for item in items:
        (FOLDER / f"{item.run_name}.toml").write_text(config_to_toml(item), encoding="utf-8")
    print(len(items), "configs")  # noqa: T201 — command-line summary
