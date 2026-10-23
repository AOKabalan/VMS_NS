"""The committed Phase 12b conservation matrix (plan 12b-1 step 3)."""

import importlib.util
from collections import Counter
from dataclasses import replace
from pathlib import Path

import pytest

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import MeshOptions, check_diagnostics_known, load_config
from navier_stokes_vms.postprocessing.diagnostics import DIAGNOSTICS, check_diagnostics_supported

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "generate_conservation_configs.py"
_SPEC = importlib.util.spec_from_file_location("generate_conservation_configs", SCRIPT)
GENERATOR = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(GENERATOR)
FOLDER = ROOT / "configs" / "chapter5" / "conservation"
FORMS = (
    "nonconservative",
    "conservative",
    "skew_symmetric_1",
    "skew_symmetric_2",
    "emac",
    "emac_vms",
)


def expected_names() -> set[str]:
    combinations = [
        f"{f}_{s}_{p}" for f in FORMS for s in ("asgs", "oss") for p in ("linear", "nonlinear")
    ]
    long_forms = ("nonconservative", "conservative", "skew_symmetric_1", "emac", "emac_vms")
    return {
        *(f"gresho_{c}" for c in combinations),
        *(f"gresho_long_{form}" for form in long_forms),
        *(f"gresho_refined_emac_vms_{space}" for space in ("asgs", "oss")),
        *(f"taylor_green_{c}" for c in combinations),
    }


def test_committed_files_are_exactly_the_generated_matrix(tmp_path: Path) -> None:
    paths = GENERATOR.write_configs(tmp_path)
    committed = sorted(FOLDER.glob("*"))
    assert [path.name for path in committed] == sorted(path.name for path in paths)
    assert {path.stem for path in committed} == expected_names()
    assert len(committed) == 55
    for path in committed:
        assert path.read_bytes() == (tmp_path / path.name).read_bytes(), path.name


def test_every_configuration_is_valid_with_a_unique_run_name() -> None:
    configs = [load_config(path) for path in sorted(FOLDER.glob("*.toml"))]
    names = [config.run_name for config in configs]
    assert len(set(names)) == len(names) == 55
    steps = Counter()
    for path, config in zip(sorted(FOLDER.glob("*.toml")), configs, strict=True):
        assert config.run_name == path.stem
        family = next(
            prefix
            for prefix in ("gresho_long", "gresho_refined", "gresho", "taylor_green")
            if config.run_name.startswith(f"{prefix}_")
        )
        steps[family] += round(config.time.t_end / config.time.dt)
        assert config.discretization.velocity_degree == config.discretization.pressure_degree == 1
        assert config.viscosity.form == "stress"
        stabilization = config.stabilization
        assert stabilization.method == "vms" and stabilization.pressure_subscale
        assert stabilization.subscale_dynamics == "dynamic"
        assert (config.time.scheme, config.time.startup) == ("crank_nicolson", "lower_order")
        solver = config.solver
        assert (solver.linearization, solver.newton_tau) == ("newton", "frozen")
        assert solver.picard_iterations_before_newton == 0
        assert (solver.outer_rtol, solver.outer_max_iterations) == (1e-11, 100)
        assert (solver.inner_rtol, solver.subscale_rtol) == (1e-12, 1e-13)
        check_diagnostics_known(config, DIAGNOSTICS)
    assert steps == {
        "gresho": 2400,
        "gresho_long": 5000,
        "gresho_refined": 400,
        "taylor_green": 2400,
    }


@pytest.mark.parametrize(
    "name",
    [
        "gresho_emac_vms_oss_nonlinear",
        "gresho_long_emac",
        "gresho_refined_emac_vms_asgs",
        "taylor_green_conservative_asgs_linear",
    ],
)
def test_family_settings_and_case_support(name: str) -> None:
    config = load_config(FOLDER / f"{name}.toml")
    family = {
        "gresho_emac_vms_oss_nonlinear": ("gresho", 32, 0.005, 0.5),
        "gresho_long_emac": ("gresho", 32, 0.01, 10.0),
        "gresho_refined_emac_vms_asgs": ("gresho", 64, 0.0025, 0.5),
        "taylor_green_conservative_asgs_linear": ("taylor_green_2d", 32, 0.005, 0.5),
    }[name]
    assert (config.case, config.mesh.cells[0], config.time.dt, config.time.t_end) == family
    assert config.mesh.cells == (family[1], family[1])
    if name == "gresho_long_emac":
        assert (config.stabilization.subscale_space, config.stabilization.scale_splitting) == (
            "oss",
            "nonlinear",
        )
    problem = build_problem(replace(config, mesh=MeshOptions(cells=(2, 2))))
    check_diagnostics_supported(config.output.diagnostics, problem)
    angular = [name for name in config.output.diagnostics if name.startswith("angular")]
    if config.case == "taylor_green_2d":
        assert dict(config.case_parameters) == {"viscosity": 0.01}
        assert angular == []
    else:
        assert problem.monitor_region is not None and len(problem.force_boundaries) == 4
        assert len(angular) == 3
