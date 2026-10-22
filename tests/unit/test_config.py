"""Configuration loading, validation, resolved dump and the run script (plan §7, §9)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from navier_stokes_vms.config import (
    Config,
    MeshOptions,
    check_diagnostics_known,
    check_startup_supported,
    config_to_toml,
    load_config,
    method_name,
)

RUN_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "run_simulation.py"

GALERKIN = """
case = "manufactured_spatial"

[discretization]
velocity_degree = 2
pressure_degree = 1

[stabilization]
method = "none"
"""

FULL = """
case = "cylinder_dfg_2d"
run_name = "cn_picard"

[case_parameters]
reynolds = 100
inflow = { profile = "parabolic", peak = 1.5 }

[mesh]
source = "gmsh"
file = "meshes/cylinder_dfg_2d.msh"

[discretization]
velocity_degree = 3
pressure_degree = 2

[time]
scheme = "crank_nicolson"
dt = 1
t_end = 8.0
startup = "exact_history"

[stabilization]
method = "none"

[convection]
form = "skew_symmetric_1"

[solver]
linearization = "picard"
newton_tau = "exact"
inner_rtol = 1e-9
inner_atol = 0
inner_max_iterations = 40
outer_rtol = 1e-7
outer_max_iterations = 20
subscale_update = "fixed_point"
subscale_rtol = 1e-9
subscale_max_iterations = 60
petsc = { ksp_type = "gmres", ksp_rtol = 1e-8, ksp_monitor = true }

[output]
vtk_every = 0
diagnostics = ["kinetic_energy"]
console_level = "DEBUG"
"""

SHORT_RUN = (
    GALERKIN
    + "\n[mesh]\ncells = [2, 2]\n\n[time]\ndt = 0.1\nt_end = 0.1\n\n[output]\nvtk_every = 0\n"
)

STATIC_ASGS = """case = "kovasznay"
[stabilization]
method = "vms"
subscale_space = "asgs"
subscale_dynamics = "static"
scale_splitting = "linear"
"""


def write_config(folder: Path, text: str, name: str = "run.toml") -> Path:
    path = folder / name
    path.write_text(text, encoding="utf-8")
    return path


def load_text(folder: Path, text: str) -> Config:
    return load_config(write_config(folder, text))


def run_script(config_path: Path, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(RUN_SCRIPT), str(config_path)],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def test_minimal_config_takes_defaults(tmp_path: Path) -> None:
    config = load_text(tmp_path, GALERKIN)
    assert config.case == "manufactured_spatial"
    assert config.run_name == "run"
    assert dict(config.case_parameters) == {}
    assert config.mesh == MeshOptions(source=None, cells=None, file=None)
    assert config.time.scheme == "bdf2"
    assert config.time.startup == "lower_order"
    assert config.solver.newton_tau == "frozen"
    assert dict(config.solver.petsc) == {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    assert config.output.diagnostics == ("kinetic_energy", "divergence_l2")


def test_full_config_loads_every_field(tmp_path: Path) -> None:
    config = load_text(tmp_path, FULL)
    assert config.run_name == "cn_picard"
    assert config.case_parameters["inflow"] == {"profile": "parabolic", "peak": 1.5}
    assert config.mesh.file == "meshes/cylinder_dfg_2d.msh"
    assert (config.discretization.velocity_degree, config.discretization.pressure_degree) == (3, 2)
    assert config.time.dt == 1.0
    assert isinstance(config.time.dt, float)
    assert config.time.startup == "exact_history"
    assert config.solver.linearization == "picard"
    assert config.solver.inner_atol == 0.0
    assert dict(config.solver.petsc) == {"ksp_type": "gmres", "ksp_rtol": 1e-8, "ksp_monitor": True}
    assert config.output.diagnostics == ("kinetic_energy",)
    assert config.output.console_level == "DEBUG"


@pytest.mark.parametrize(
    ("table", "expected"),
    [
        ("source = 'gmsh'", MeshOptions(source="gmsh")),
        ("file = 'a.msh'", MeshOptions(file="a.msh")),
        ("cells = [8, 8, 8]", MeshOptions(cells=(8, 8, 8))),
        ("source = 'builtin'\ncells = [4, 6]", MeshOptions(source="builtin", cells=(4, 6))),
    ],
    ids=["source-only", "file-only", "cells-only", "builtin-cells"],
)
def test_mesh_fields_override_the_case_default_independently(
    tmp_path: Path, table: str, expected: MeshOptions
) -> None:
    assert load_text(tmp_path, GALERKIN + f"\n[mesh]\n{table}\n").mesh == expected


@pytest.mark.parametrize("text", [GALERKIN, FULL, STATIC_ASGS], ids=["minimal", "full", "asgs"])
def test_resolved_dump_reloads_to_the_same_config(tmp_path: Path, text: str) -> None:
    config = load_text(tmp_path, text)
    dumped = write_config(tmp_path, config_to_toml(config), name="other_name.toml")
    assert load_config(dumped) == config


def test_rule1_galerkin_needs_taylor_hood(tmp_path: Path) -> None:
    text = GALERKIN.replace("velocity_degree = 2", "velocity_degree = 1")
    with pytest.raises(ValueError, match="discretization.velocity_degree = 1, .*pressure_degree"):
        load_text(tmp_path, text)


def test_rule1_galerkin_rejects_stabilization_fields(tmp_path: Path) -> None:
    text = GALERKIN + 'c1 = 12.0\nsubscale_space = "asgs"\n'
    with pytest.raises(ValueError, match="stabilization.subscale_space, stabilization.c1"):
        load_text(tmp_path, text)


def test_rule2_vms_needs_equal_order(tmp_path: Path) -> None:
    text = GALERKIN.replace('method = "none"', 'method = "vms"')
    with pytest.raises(ValueError, match="method = 'vms' needs equal orders"):
        load_text(tmp_path, text)


def test_rule3_exact_history_needs_an_exact_solution(tmp_path: Path) -> None:
    exact_history = load_text(tmp_path, GALERKIN + '\n[time]\nstartup = "exact_history"\n')
    with pytest.raises(ValueError, match="time.startup = 'exact_history'"):
        check_startup_supported(exact_history, has_exact_solution=False)
    check_startup_supported(exact_history, has_exact_solution=True)
    check_startup_supported(load_text(tmp_path, GALERKIN), has_exact_solution=False)


def test_rule4_picard_warmup_needs_newton(tmp_path: Path) -> None:
    text = GALERKIN + '\n[solver]\nlinearization = "picard"\npicard_iterations_before_newton = 2\n'
    with pytest.raises(ValueError, match="solver.picard_iterations_before_newton = 2 requires"):
        load_text(tmp_path, text)


def test_rule6_diagnostics_must_be_registered(tmp_path: Path) -> None:
    config = load_text(tmp_path, GALERKIN)
    with pytest.raises(ValueError, match="output.diagnostics: unknown divergence_l2"):
        check_diagnostics_known(config, available={"kinetic_energy"})
    check_diagnostics_known(config, available={"kinetic_energy", "divergence_l2"})


@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("pressure_subscale", ["true", "false"])
@pytest.mark.parametrize("space", ["asgs", "oss"])
@pytest.mark.parametrize(("dynamics", "name"), [("static", "Sta"), ("dynamic", "Dyn")])
@pytest.mark.parametrize(("splitting", "label"), [("linear", "Lin"), ("nonlinear", "Nl")])
def test_rule5_accepts_all_eight_cases(
    tmp_path: Path,
    degree: int,
    pressure_subscale: str,
    space: str,
    dynamics: str,
    name: str,
    splitting: str,
    label: str,
) -> None:
    text = (
        STATIC_ASGS.replace('"asgs"', f'"{space}"')
        .replace('"static"', f'"{dynamics}"')
        .replace('"linear"', f'"{splitting}"')
        + f"pressure_subscale = {pressure_subscale}\n[discretization]\n"
        f"velocity_degree = {degree}\npressure_degree = {degree}\n"
    )
    config = load_text(tmp_path, text)
    assert method_name(config.stabilization) == f"{name}-{label}-{space.upper()}"
    assert config.discretization.velocity_degree == degree


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("[stabilization]\nmethod = 'none'\n", "case: required"),
        (GALERKIN + "\n[physics]\nnu = 1.0\n", "unknown top-level keys: physics"),
        (GALERKIN + "\n[time]\ndelta_t = 0.1\n", r"unknown keys in \[time\]: delta_t"),
        (GALERKIN + "\n[time]\ndt = 'small'\n", "time.dt = 'small': expected float"),
        (GALERKIN + "\n[output]\nvtk_every = true\n", "output.vtk_every = True: expected int"),
        (GALERKIN + "\n[time]\nscheme = 'rk4'\n", "time.scheme = 'rk4': expected one of"),
        (GALERKIN + "\n[time]\ndt = 0.0\n", "time.dt = 0.0: must be > 0"),
        (GALERKIN + "\n[time]\nt_end = inf\n", "time.t_end = inf: must be > 0"),
        (GALERKIN + "\n[output]\nvtk_every = -1\n", "output.vtk_every = -1: must be >= 0"),
        (GALERKIN + "\n[time]\ndt = 0.3\nt_end = 1.0\n", "whole number of uniform steps"),
        (GALERKIN + "\n[mesh]\nsource = 'builtin'\nfile = 'a.msh'\n", "mesh.file: only used when"),
        (
            GALERKIN + "\n[mesh]\nsource = 'gmsh'\ncells = [8, 8]\n",
            "mesh.cells: only used when mesh.source = 'builtin'",
        ),
        (GALERKIN + "\n[mesh]\ncells = [8, 8]\nfile = 'a.msh'\n", "mesh.cells and mesh.file"),
        (GALERKIN + "\n[mesh]\ncells = [8]\n", r"mesh.cells = \[8\]: needs 2 or 3 entries"),
        (GALERKIN + "\n[mesh]\ncells = [0, 4]\n", "each >= 1"),
        (GALERKIN + "\n[mesh]\nlengths = [1.0, 1.0]\n", r"unknown keys in \[mesh\]: lengths"),
        (GALERKIN + "\n[solver]\npetsc = { a = [1, 2] }\n", "solver.petsc: a must be"),
        (GALERKIN.replace("manufactured_spatial", "../elsewhere"), "case = '../elsewhere'"),
        ('run_name = ".hidden"\n' + GALERKIN, "run_name = '.hidden'"),
        (GALERKIN + "\n[time\n", "Expected ']'"),
    ],
    ids=[
        "missing-case",
        "unknown-section",
        "unknown-key",
        "wrong-type",
        "bool-for-int",
        "unknown-choice",
        "zero-dt",
        "infinite-t_end",
        "negative-vtk_every",
        "t_end-not-whole-steps",
        "builtin-with-file",
        "gmsh-with-cells",
        "cells-and-file",
        "cells-entries",
        "zero-cells",
        "geometry-not-in-config",
        "petsc-array",
        "unsafe-case",
        "unsafe-run_name",
        "toml-syntax",
    ],
)
def test_invalid_input_is_rejected_naming_the_field(
    tmp_path: Path, text: str, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        load_text(tmp_path, text)


def test_run_script_writes_resolved_config_run_info_and_log(tmp_path: Path) -> None:
    config_path = write_config(tmp_path, SHORT_RUN, name="galerkin_p2p1.toml")
    result = run_script(config_path, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    folder = tmp_path / "results" / "manufactured_spatial" / "galerkin_p2p1"
    assert sorted(path.name for path in folder.iterdir()) == [
        "config_resolved.toml",
        "diagnostics.csv",
        "inner_iterations.csv",
        "outer_iterations.csv",
        "run.log",
        "run_info.json",
        "time_steps.csv",
    ]
    assert load_config(folder / "config_resolved.toml") == load_config(config_path)
    info = json.loads((folder / "run_info.json").read_text(encoding="utf-8"))
    assert info["exit_status"] == 0
    assert info["mpi_size"] == 1
    assert info["wall_time_seconds"] > 0
    assert info["start_time"] <= info["end_time"]
    expected = {"git_commit", "git_dirty", "firedrake_version", "petsc_version", "python_version"}
    assert expected | {"host"} <= info.keys()
    log = (folder / "run.log").read_text(encoding="utf-8")
    assert "Case manufactured_spatial, run galerkin_p2p1" in log
    assert "[stabilization] method = none" in log
    assert "step 1, t = 0.1, bdf1" in log


def test_run_script_refuses_an_existing_results_folder(tmp_path: Path) -> None:
    config_path = write_config(tmp_path, SHORT_RUN)
    assert run_script(config_path, cwd=tmp_path).returncode == 0
    second = run_script(config_path, cwd=tmp_path)
    assert second.returncode == 1
    assert "already exists" in second.stderr


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (
            SHORT_RUN.replace("vtk_every = 0", "vtk_every = 0\ndiagnostics = ['vorticity']"),
            "output.diagnostics: unknown vorticity; available: angular_momentum",
        ),
        (
            SHORT_RUN.replace(
                "vtk_every = 0", "vtk_every = 0\ndiagnostics = ['subscale_kinetic_energy']"
            ),
            "subscale_kinetic_energy requires stabilization.method = 'vms'",
        ),
        (
            SHORT_RUN.replace("vtk_every = 0", "vtk_every = 0\ndiagnostics = ['force']"),
            "force requires the case to declare force_boundaries",
        ),
        (SHORT_RUN.replace("manufactured_spatial", "cavity"), "case = 'cavity': unknown"),
    ],
    ids=["rule-6", "subscale-without-vms", "force-without-boundary", "unknown-case"],
)
def test_run_script_rejects_an_invalid_setup_before_writing(
    tmp_path: Path, text: str, message: str
) -> None:
    result = run_script(write_config(tmp_path, text), cwd=tmp_path)
    assert result.returncode == 1
    assert message in result.stderr
    assert not (tmp_path / "results").exists()
