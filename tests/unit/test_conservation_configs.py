"""The committed Phase 12b conservation matrix (plan 12b-1 step 3)."""

import importlib.util
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
        *(f"gresho_{f}_oss_{p}" for f in FORMS for p in ("linear", "nonlinear")),
        *(f"gresho_visc_{f}_asgs_{p}" for f in FORMS for p in ("linear", "nonlinear")),
        *(f"gresho_long_{form}" for form in long_forms),
        "gresho_refined_emac_vms_oss",
        "gresho_visc_refined_emac_vms_asgs",
        *(f"taylor_green_{c}" for c in combinations),
    }


def test_committed_files_are_exactly_the_generated_matrix(tmp_path: Path) -> None:
    # Regeneration must remove the old ASGS runs so the queue cannot pick them up.
    retired = [f"gresho_{f}_asgs_{p}" for f in FORMS for p in ("linear", "nonlinear")]
    for name in [*retired, "gresho_refined_emac_vms_asgs"]:
        (tmp_path / f"{name}.toml").write_text("retired")
    paths = GENERATOR.write_configs(tmp_path)
    assert {path.stem for path in tmp_path.glob("*.toml")} == expected_names()
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
    for path, config in zip(sorted(FOLDER.glob("*.toml")), configs, strict=True):
        assert config.run_name == path.stem
        check_diagnostics_known(config, DIAGNOSTICS)


@pytest.mark.parametrize(
    "name",
    [
        "gresho_emac_vms_oss_nonlinear",
        "gresho_visc_emac_vms_asgs_nonlinear",
        "taylor_green_conservative_asgs_linear",
    ],
)
def test_conservation_case_diagnostics_are_supported(name: str) -> None:
    config = load_config(FOLDER / f"{name}.toml")
    problem = build_problem(replace(config, mesh=MeshOptions(cells=(2, 2))))
    check_diagnostics_supported(config.output.diagnostics, problem)
    angular = [name for name in config.output.diagnostics if name.startswith("angular")]
    if config.case == "taylor_green_2d":
        assert dict(config.case_parameters) == {"viscosity": 0.01}
        assert angular == []
    else:
        assert problem.viscosity == (1e-3 if name.startswith("gresho_visc_") else 1e-6)
        assert problem.monitor_region is not None and len(problem.force_boundaries) == 4
        assert len(angular) == 3
