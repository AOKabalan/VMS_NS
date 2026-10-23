"""The resumable run queue, with stub runners in a temporary copy of the repository layout.

A stub mpiexec runs its command once, so no simulation is launched; the stubs record each call.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

QUEUE = Path(__file__).resolve().parents[2] / "scripts" / "run_validation_queue.sh"
SIMULATION_STUB = """import json, sys, tomllib
from pathlib import Path
config = tomllib.loads(Path(sys.argv[1]).read_text())
folder = Path("results") / config["case"] / config["run_name"]
folder.mkdir(parents=True)
status = 1 if "fail" in config["run_name"] else 0
(folder / "run_info.json").write_text(json.dumps({"exit_status": status}, indent=2))
with open("calls.txt", "a") as stream:
    stream.write(f"simulation {sys.argv[1]}\\n")
sys.exit(status)
"""
VALIDATION_STUB = """import sys, tomllib
from pathlib import Path
config = tomllib.loads(Path(sys.argv[1]).read_text())
root = Path(sys.argv[sys.argv.index("--results-root") + 1])
folder = root / config["case"] / config["run_name"]
if "--run-only" in sys.argv:
    folder.mkdir(parents=True)
else:
    (folder / "validation.csv").write_text("ok\\n")
with open("calls.txt", "a") as stream:
    stream.write(" ".join(["validation", *sys.argv[1:]]) + "\\n")
"""


def config_text(case: str, run_name: str, cells: int, steps: int) -> str:
    return (
        f'case = "{case}"\nrun_name = "{run_name}"\n[mesh]\ncells = [{cells}, {cells}]\n'
        "[discretization]\nvelocity_degree = 1\npressure_degree = 1\n"
        '[stabilization]\nmethod = "vms"\nsubscale_space = "oss"\nscale_splitting = "linear"\n'
        f"[time]\ndt = 0.1\nt_end = {0.1 * steps}\n"
    )


@pytest.fixture
def layout(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    """A repository layout with the real queue, stub runners and a stub mpiexec on PATH."""
    (tmp_path / "scripts").mkdir()
    shutil.copy(QUEUE, tmp_path / "scripts" / QUEUE.name)
    (tmp_path / "scripts" / "run_simulation.py").write_text(SIMULATION_STUB)
    (tmp_path / "scripts" / "run_validation.py").write_text(VALIDATION_STUB)
    bin_folder = tmp_path / "bin"
    bin_folder.mkdir()
    mpiexec = bin_folder / "mpiexec"
    mpiexec.write_text('#!/bin/sh\nshift 2\nexec "$@"\n')
    mpiexec.chmod(0o755)
    environment = {**os.environ, "PATH": f"{bin_folder}{os.pathsep}{os.environ['PATH']}"}
    return tmp_path, environment


def run_queue(folder: Path, environment: dict[str, str], *arguments: str) -> int:
    script = folder / "scripts" / QUEUE.name
    command = ["bash", str(script), *arguments]
    return subprocess.run(command, cwd=folder, env=environment, check=False).returncode


def calls(folder: Path) -> list[str]:
    path = folder / "calls.txt"
    return path.read_text().splitlines() if path.exists() else []


def test_simulation_queue_resumes_and_keeps_recorded_failures(
    layout: tuple[Path, dict[str, str]],
) -> None:
    folder, environment = layout
    configs = folder / "configs" / "chapter5" / "conservation"
    configs.mkdir(parents=True)
    runs = {
        "taylor_green_e": ("taylor_green_2d", 4, 10),
        "gresho_b_fail": ("gresho", 8, 10),
        "gresho_a": ("gresho", 4, 10),
        "gresho_c_done": ("gresho", 4, 20),
        "gresho_d_interrupted": ("gresho", 8, 20),
    }
    for name, (case, cells, steps) in runs.items():
        (configs / f"{name}.toml").write_text(config_text(case, name, cells, steps))
    results = folder / "results"
    (results / "gresho" / "gresho_c_done").mkdir(parents=True)
    (results / "gresho" / "gresho_c_done" / "run_info.json").write_text('{"exit_status": 0}')
    (results / "gresho" / "gresho_d_interrupted").mkdir()
    (results / "gresho" / "gresho_d_interrupted" / "run_info.json").write_text(
        '{"exit_status": null}'
    )
    (results / "VALIDATION_DONE").write_text("validation summary\n")
    arguments = ("--configs", "configs/chapter5/conservation/", "--runner", "simulation")

    assert run_queue(folder, environment, *arguments, "--dry-run") == 0
    assert calls(folder) == []
    assert not (results / "conservation_queue.log").exists()

    assert run_queue(folder, environment, *arguments) == 1
    prefix = "simulation configs/chapter5/conservation/"
    assert calls(folder) == [
        f"{prefix}gresho_a.toml",
        f"{prefix}gresho_b_fail.toml",
        f"{prefix}gresho_d_interrupted.toml",
        f"{prefix}taylor_green_e.toml",
    ]
    moved = list((results / "gresho").glob("gresho_d_interrupted.incomplete-*"))
    assert len(moved) == 1
    log = [
        line.split("\t") for line in (results / "conservation_queue.log").read_text().splitlines()
    ]
    assert [(Path(row[0]).stem, row[1], row[3]) for row in log] == [
        ("gresho_a", "0", "ran"),
        ("gresho_c_done", "0", "skipped"),
        ("gresho_b_fail", "1", "ran"),
        ("gresho_d_interrupted", "0", "ran"),
        ("taylor_green_e", "0", "ran"),
    ]
    assert (results / "CONSERVATION_DONE").read_text() == (
        "Conservation queue: 3 completed, 1 skipped, 1 failed\n"
    )
    assert (results / "VALIDATION_DONE").read_text() == "validation summary\n"

    # A finished run that failed to converge is a recorded result, not repeated.
    assert run_queue(folder, environment, *arguments) == 0
    assert len(calls(folder)) == 4
    assert (results / "CONSERVATION_DONE").read_text() == (
        "Conservation queue: 0 completed, 5 skipped, 0 failed\n"
    )


@pytest.mark.parametrize("smoke", [False, True])
def test_default_queue_runs_and_analyses_validation_configs(
    layout: tuple[Path, dict[str, str]], smoke: bool
) -> None:
    folder, environment = layout
    configs = folder / "configs" / "validation"
    configs.mkdir(parents=True)
    (configs / "cavity_a.toml").write_text(config_text("lid_driven_cavity", "cavity_a", 4, 5))
    (configs / "cavity_b.toml").write_text(config_text("lid_driven_cavity", "cavity_b", 8, 5))
    root = folder / ("results/validation_smoke" if smoke else "results")
    (root / "lid_driven_cavity" / "cavity_b").mkdir(parents=True)
    (root / "lid_driven_cavity" / "cavity_b" / "validation.csv").write_text("ok\n")
    assert run_queue(folder, environment, *(["--smoke"] if smoke else [])) == 0
    root_name = "results/validation_smoke" if smoke else "results"
    smoke_flag = " --smoke" if smoke else ""
    assert calls(folder) == [
        f"validation configs/validation/cavity_a.toml --run-only --results-root {root_name}"
        + smoke_flag,
        f"validation configs/validation/cavity_a.toml --analyse-only --results-root {root_name}",
    ]
    assert (root / "VALIDATION_DONE").read_text() == (
        "Validation queue: 1 completed, 1 skipped, 0 failed\n"
    )
    assert len((root / "validation_queue.log").read_text().splitlines()) == 2


@pytest.mark.parametrize(
    "arguments",
    [
        ("--runner", "other"),
        ("--configs", "configs/missing"),
        ("--smoke", "--runner", "simulation"),
        ("--configs",),
        ("--unknown",),
    ],
)
def test_invalid_usage_is_rejected(
    layout: tuple[Path, dict[str, str]], arguments: tuple[str, ...]
) -> None:
    folder, environment = layout
    (folder / "configs" / "validation").mkdir(parents=True)
    assert run_queue(folder, environment, *arguments) == 2
    assert not (folder / "results").exists()
