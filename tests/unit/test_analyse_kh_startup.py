"""Startup analysis rejects incomplete/mismatched runs and preserves raw energy comparisons."""

import csv
import importlib.util
import json
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location(
    "analyse_kh_startup", Path(__file__).resolve().parents[2] / "scripts/analyse_kh_startup.py"
)
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def write_table(path: Path, records: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def fixture_run(folder: Path, epsilon: float) -> Path:
    folder.mkdir()
    (folder / "run_info.json").write_text(json.dumps({"exit_status": 0, "wall_time_seconds": 1.0}))
    (folder / "config_resolved.toml").write_text(f'''case = "kelvin_helmholtz_shear_layer"
run_name = "{folder.name}"
[time]
dt = 0.1
t_end = 0.2
scheme = "bdf2"
[discretization]
velocity_degree = 1
[stabilization]
subscale_dynamics = "dynamic"
[output]
horizontal_profiles_every = 1
''')
    metadata = {
        "mesh_cells": [4, 4],
        "parameters": {"epsilon": epsilon, "clustering": 0},
        "mesh_quality": {"minimum_edge": 0.25},
    }
    (folder / "thermal_metadata.json").write_text(json.dumps({"case_metadata": metadata}))
    (folder / "initial_state.json").write_text(json.dumps({"vertical_kinetic_energy": epsilon**2}))
    write_table(
        folder / "diagnostics.csv",
        [
            {
                "step": step,
                "time": step * 0.1,
                "vertical_kinetic_energy": 0.01 + epsilon**2,
                "advective_cfl": 0.1,
                "boussinesq_energy_budget_residual": 1e-13,
            }
            for step in (1, 2)
        ],
    )
    write_table(
        folder / "horizontal_profiles.csv",
        [
            {
                "step": step,
                "time": step * 0.1,
                "z": z,
                "velocity_x": z + 0.1 * step,
                "velocity_y": 0,
                "temperature": z - 0.2 * step,
            }
            for step in (0, 1, 2)
            for z in (-1, 0, 1)
        ],
    )
    return folder


def test_raw_energies_and_profile_drift_are_retained(tmp_path: Path) -> None:
    control = fixture_run(tmp_path / "control", 0)
    seeded = fixture_run(tmp_path / "seeded", 0.1)
    report = ANALYSIS.compare([seeded, control])
    pair = report["pairs"][0]
    assert pair["control_final_over_seeded_final"] == pytest.approx(0.5)
    assert pair["maximum_control_over_initial_seed_energy"] == pytest.approx(1)
    assert pair["control"]["profile_maximum_absolute_drift"]["velocity_x"] == pytest.approx(
        [0, 0.1, 0.2]
    )
    assert pair["seeded"]["profile_maximum_absolute_drift"]["temperature"] == pytest.approx(
        [0, 0.2, 0.4]
    )
    assert pair["control"]["vertical_energy"] == [0.01, 0.01]
    output = tmp_path / "report.json"
    assert ANALYSIS.main([str(control), str(seeded), "--output", str(output)]) == 0
    assert json.loads(output.read_text())["phase"] == "B7.2"


@pytest.mark.parametrize(
    "damage",
    [
        "failed",
        "missing_initial",
        "missing_final",
        "nan_profile",
        "negative_energy",
        "missing_step",
        "different_rows",
        "different_time",
    ],
)
def test_damaged_runs_are_rejected(tmp_path: Path, damage: str) -> None:
    folder = fixture_run(tmp_path / "control", 0)
    if damage == "failed":
        (folder / "run_info.json").write_text(json.dumps({"exit_status": 1}))
    elif damage in ("negative_energy", "missing_step"):
        rows = ANALYSIS.read_rows(folder / "diagnostics.csv")
        if damage == "missing_step":
            rows.pop()
        else:
            rows[0]["vertical_kinetic_energy"] = -1
        write_table(folder / "diagnostics.csv", rows)
    else:
        rows = ANALYSIS.read_rows(folder / "horizontal_profiles.csv")
        if damage == "missing_initial":
            rows = rows[3:]
        elif damage == "missing_final":
            rows = rows[:-3]
        elif damage == "nan_profile":
            rows[-1]["temperature"] = float("nan")
        elif damage == "different_rows":
            rows[-1]["z"] = 2
        else:
            rows[-1]["time"] = 0.3
        write_table(folder / "horizontal_profiles.csv", rows)
    with pytest.raises(ValueError):
        ANALYSIS.read_run(folder)


@pytest.mark.parametrize("mismatch", ["missing_pair", "duplicate_control", "physics"])
def test_unmatched_pairs_are_rejected(tmp_path: Path, mismatch: str) -> None:
    control = fixture_run(tmp_path / "control", 0)
    seeded = fixture_run(tmp_path / "seeded", 0.1)
    paths = [control, seeded]
    if mismatch == "missing_pair":
        paths.pop()
    elif mismatch == "duplicate_control":
        paths.append(fixture_run(tmp_path / "extra", 0))
    else:
        path = seeded / "thermal_metadata.json"
        data = json.loads(path.read_text())
        data["case_metadata"]["parameters"]["clustering"] = 1
        path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="exactly one"):
        ANALYSIS.compare(paths)


@pytest.fixture
def growth(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / "scripts"))
    return importlib.import_module("analyse_kh_growth")


def test_fourier_growth_uses_difference_coefficients_before_squaring(
    tmp_path: Path, growth: ModuleType
) -> None:
    times = np.arange(7) * 0.1
    rate = 0.17
    folders = []
    for name, epsilon, amplitude in (("control", 0, 0.3), ("seeded", 0.1, 0.4)):
        folder = fixture_run(tmp_path / name, epsilon)
        folders.append(folder)
        config = folder / "config_resolved.toml"
        config.write_text(config.read_text().replace("t_end = 0.2", "t_end = 0.6"))
        (folder / "horizontal_profiles.json").write_text(
            json.dumps(
                {
                    "fourier_modes": [1],
                    "heights": [-1, 0, 1],
                    "length": 2 * np.pi,
                }
            )
        )
        write_table(
            folder / "horizontal_profiles.csv",
            [
                {
                    "step": step,
                    "time": time,
                    "z": z,
                    "velocity_x": z,
                    "velocity_y": 0,
                    "temperature": z,
                    "velocity_y_mode_1_cos": amplitude * np.exp(rate * time) * (z + 2),
                    "velocity_y_mode_1_sin": 0,
                }
                for step, time in enumerate(times)
                for z in (-1, 0, 1)
            ],
        )
        write_table(
            folder / "diagnostics.csv",
            [
                {
                    "step": step,
                    "time": time,
                    "vertical_kinetic_energy": 20 + amplitude**2 * np.exp(2 * rate * time),
                    "advective_cfl": 0.1,
                    "boussinesq_energy_budget_residual": 1e-13,
                }
                for step, time in enumerate(times[1:], 1)
            ],
        )
    report = growth.compare(folders, 1, (0, 0.6))
    pair = report["pairs"][0]
    expected = 2 * np.pi * (26 / 3) * 0.1**2 * np.exp(2 * rate * times)
    np.testing.assert_allclose(pair["energies"]["difference_mode"], expected, rtol=1e-13)
    assert pair["fits"]["difference_mode"]["sigma"] == pytest.approx(rate, abs=1e-14)
    assert pair["fits"]["difference_mode"]["rms_log_energy_residual"] < 1e-14
    assert pair["maximum_control_over_difference_mode"] == pytest.approx(9)
    assert pair["window_endpoint_mode_shape_coherence"] == pytest.approx(1)
    output = tmp_path / "growth.json"
    assert (
        growth.main([*(str(p) for p in folders), "--window", "0", "0.6", "--output", str(output)])
        == 0
    )
    assert json.loads(output.read_text())["phase"] == "B7.3"
    with pytest.raises(ValueError, match="not recorded"):
        growth.compare(folders, 2, None)


def test_growth_fits_reject_bad_windows_and_report_curvature(growth: ModuleType) -> None:
    times = np.linspace(0, 1, 11)
    exponential = np.exp(0.4 * times)
    curved = np.exp(0.4 * times + times**2)
    assert growth.fit_energy(times, curved, (0, 1))["half_window_sigma_difference"] > 0.4
    for window in ((0.2, 0.2), (0, 2), (0, float("nan")), (0, 0.1)):
        with pytest.raises(ValueError):
            growth.fit_energy(times, exponential, window)
    with pytest.raises(ValueError, match="positive finite"):
        growth.fit_energy(times, np.zeros_like(times), (0, 1))
