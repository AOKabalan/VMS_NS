"""The conservation analysis on a synthetic results tree (plan 12b-1 step 5, owner decision K-D2).

The Gresho runs use |Ω| = 1, d = 2 and E(0) = 0.08, so U_c = 0.4, L_c = 1, M_c = A_c = 0.4,
F_c = 0.16 (also the torque scale) and P_c = 0.064.
"""

import csv
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "analyse_conservation.py"
_SPEC = importlib.util.spec_from_file_location("analyse_conservation", SCRIPT)
ANALYSIS = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ANALYSIS)
CONFIGS = ROOT / "configs" / "chapter5" / "conservation"
GRESHO_INITIAL = {
    "time": 0.0,
    "domain_volume": 1.0,
    "dimension": 2,
    "kinetic_energy": 0.08,
    "linear_momentum_x": 0.0,
    "linear_momentum_y": 0.0,
    "angular_momentum": 0.05,
}


def write_run(
    results: Path,
    case: str,
    name: str,
    exit_status: int | None,
    steps: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
    initial: dict[str, float] | None,
    log: str = "",
) -> None:
    folder = results / case / name
    folder.mkdir(parents=True)
    (folder / "run_info.json").write_text(json.dumps({"exit_status": exit_status}, indent=2))
    for file_name, rows in (("time_steps.csv", steps), ("diagnostics.csv", diagnostics)):
        columns = list(dict.fromkeys(column for row in rows for column in row))
        with (folder / file_name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, columns)
            writer.writeheader()
            writer.writerows(rows)
    if initial is not None:
        (folder / "initial_state.json").write_text(json.dumps(initial))
    (folder / "run.log").write_text(log)


def steps(count: int, converged: bool = True) -> list[dict[str, object]]:
    return [
        {
            "step": step,
            "converged": converged,
            "outer_iterations": 5 + step,
            "inner_iterations_total": 10 * step,
        }
        for step in range(1, count + 1)
    ]


def gresho_rows(measured: float, predicted: float) -> list[dict[str, object]]:
    """Two steps: blank budgets at the BDF1 startup, then the midpoint closures."""
    common = {
        "linear_momentum_defect_x": measured,
        "linear_momentum_defect_y": -2 * measured,
        "linear_momentum_defect_predicted_x": predicted,
        "linear_momentum_defect_predicted_y": -2 * predicted,
        "angular_momentum_defect": measured,
        "angular_momentum_defect_predicted": predicted,
        "convective_energy_defect": measured / 2.5,
    }
    return [
        {
            "step": 1,
            "time": 0.005,
            "kinetic_energy": 0.0799,
            "linear_momentum_x": 0.001,
            "linear_momentum_y": 0.0,
            "angular_momentum": 0.0495,
            "energy_budget_residual": "",
            "energy_budget_defect_predicted": "",
            "momentum_budget_residual_x": "",
            "momentum_budget_residual_y": "",
            **common,
            "region_velocity_l2": 1e-4,
            "region_subscale_l2": 2e-5,
            "velocity_error_l2": 0.01,
            "force_x_min_x": 0.0,
            "force_x_min_y": 0.0,
        },
        {
            "step": 2,
            "time": 0.01,
            "kinetic_energy": 0.0796,
            "linear_momentum_x": 0.004,
            "linear_momentum_y": -0.003,
            "angular_momentum": 0.049,
            "energy_budget_residual": -6.4e-12,
            "energy_budget_defect_predicted": predicted / 2.5,
            "momentum_budget_residual_x": 1.6e-3,
            "momentum_budget_residual_y": -3.2e-3,
            **common,
            "region_velocity_l2": 3e-4,
            "region_subscale_l2": 5e-5,
            "velocity_error_l2": 0.02,
            "force_x_min_x": 0.0016,
            "force_x_min_y": -0.0032,
        },
    ]


def long_rows() -> list[dict[str, object]]:
    return [
        {
            "step": step,
            "time": 0.01 * step,
            "kinetic_energy": 0.08 - 1e-4 * step,
            "linear_momentum_x": 1e-5 * step,
            "linear_momentum_y": 0.0,
            "angular_momentum": 0.05 - 1e-4 * step,
            "velocity_error_l2": 1e-3 * step,
        }
        for step in range(1, 4)
    ]


def test_scales_follow_the_owner_definition() -> None:
    scales = ANALYSIS.run_scales(GRESHO_INITIAL)
    assert scales.momentum == pytest.approx(0.4)
    assert scales.angular_momentum == pytest.approx(0.4)
    assert scales.force == pytest.approx(0.16)
    assert scales.torque == pytest.approx(0.16)
    assert scales.power == pytest.approx(0.064)
    taylor_green = {**GRESHO_INITIAL, "domain_volume": 4 * math.pi**2, "kinetic_energy": math.pi**2}
    scales = ANALYSIS.run_scales(taylor_green)
    assert scales.force == pytest.approx(math.pi)
    assert scales.power == pytest.approx(2 * math.pi * 2**-1.5)
    assert ANALYSIS.run_scales({**GRESHO_INITIAL, "kinetic_energy": 0.0}) is None


def test_tables_failures_and_plot_on_a_synthetic_tree(tmp_path: Path) -> None:
    results, output = tmp_path / "results", tmp_path / "analysis"
    write_run(
        results,
        "gresho",
        "gresho_emac_vms_oss_nonlinear",
        0,
        steps(2),
        gresho_rows(3.2e-14, 0.0),
        GRESHO_INITIAL,
    )
    write_run(
        results,
        "gresho",
        "gresho_nonconservative_oss_linear",
        0,
        steps(2),
        gresho_rows(1e-3, 1e-3 + 1.6e-13),
        GRESHO_INITIAL,
    )
    write_run(
        results,
        "gresho",
        "gresho_conservative_oss_linear",
        0,
        steps(2),
        gresho_rows(1e-3, 1e-3 + 1.6e-6),
        GRESHO_INITIAL,
    )
    write_run(
        results,
        "gresho",
        "gresho_visc_skew_symmetric_1_asgs_nonlinear",
        1,
        steps(1, converged=False),
        [],
        GRESHO_INITIAL,
        log="2026-10-23 12:00:00 ERROR   Run stopped: step 1 did not converge: stalled\n",
    )
    write_run(results, "gresho", "gresho_emac_oss_linear", None, [], [], None)
    for form in ("emac", "emac_vms"):
        write_run(
            results, "gresho", f"gresho_long_{form}", 0, steps(3), long_rows(), GRESHO_INITIAL
        )
    taylor_green_initial = {
        **GRESHO_INITIAL,
        "domain_volume": 4 * math.pi**2,
        "kinetic_energy": math.pi**2,
        "angular_momentum": 0.0,
    }
    write_run(
        results,
        "taylor_green_2d",
        "taylor_green_conservative_asgs_linear",
        0,
        steps(2),
        [
            {
                "step": 1,
                "time": 0.005,
                "energy_budget_residual": "",
                "momentum_budget_residual_x": "",
                "momentum_budget_residual_y": "",
                "kinetic_energy": math.pi**2,
            },
            {
                "step": 2,
                "time": 0.01,
                "energy_budget_residual": 1e-12,
                "momentum_budget_residual_x": 2e-9 * math.pi,
                "momentum_budget_residual_y": 0.0,
                "kinetic_energy": math.pi**2 * 0.99,
            },
        ],
        taylor_green_initial,
    )

    tables = ANALYSIS.analyse(CONFIGS, results, output)

    statuses = [row["status"] for rows in tables.values() for row in rows]
    assert len(statuses) == 55
    assert (statuses.count("completed"), statuses.count("failed")) == (6, 1)
    assert (statuses.count("incomplete"), statuses.count("missing")) == (1, 47)
    gresho = {row["run"]: row for row in tables["gresho"]}
    emac_vms = gresho["gresho_emac_vms_oss_nonlinear"]
    assert emac_vms["steps_converged"] == 2 and emac_vms["outer_iterations_max"] == 7
    assert emac_vms["energy_closure"] == pytest.approx(6.4e-12)
    assert emac_vms["energy_closure_scaled"] == pytest.approx(1e-10)
    assert emac_vms["momentum_closure_scaled"] == pytest.approx(0.02)
    assert emac_vms["momentum_closure_within"] is None  # the Gresho reaction is only reported
    assert emac_vms["linear_defect_measured_scaled"] == pytest.approx(6.4e-14 / 0.16)
    assert emac_vms["angular_defect_measured_scaled"] == pytest.approx(3.2e-14 / 0.16)
    assert emac_vms["energy_defect_measured_scaled"] == pytest.approx(1.28e-14 / 0.064)
    assert emac_vms["energy_closure_within"] is True and emac_vms["defects_within"] is True
    assert emac_vms["momentum_drift_final_scaled"] == pytest.approx(0.005 / 0.4)
    assert emac_vms["angular_momentum_drift_final_scaled"] == pytest.approx(-0.001 / 0.4)
    assert emac_vms["energy_drift_final_scaled"] == pytest.approx(-0.0004 / 0.08)
    assert emac_vms["energy_drift_max"] == pytest.approx(0.0004)
    assert emac_vms["wall_force_final_max_scaled"] == pytest.approx(0.02)
    assert emac_vms["region_velocity_final"] == 3e-4
    nonconservative = gresho["gresho_nonconservative_oss_linear"]
    assert nonconservative["linear_defect_difference_scaled"] == pytest.approx(3.2e-13 / 0.16)
    assert nonconservative["linear_defect_predicted_scaled"] == pytest.approx(2e-3 / 0.16, rel=1e-6)
    assert nonconservative["defects_within"] is True
    assert gresho["gresho_conservative_oss_linear"]["defects_within"] is False
    (failed,) = [row for row in tables["gresho_visc"] if row["status"] == "failed"]
    assert failed["run"] == "gresho_visc_skew_symmetric_1_asgs_nonlinear"
    assert (failed["status"], failed["failed_step"], failed["inner_iterations_max"]) == (
        "failed",
        1,
        10,
    )
    (taylor_green,) = [row for row in tables["taylor_green"] if row["status"] == "completed"]
    assert taylor_green["energy_closure_scaled"] == pytest.approx(1e-12 / (2 * math.pi * 2**-1.5))
    assert taylor_green["momentum_closure_scaled"] == pytest.approx(2e-9)
    assert taylor_green["momentum_closure_within"] is False
    assert "angular_momentum_drift_final" not in taylor_green

    with (output / "gresho.csv").open() as stream:
        written = {row["run"]: row for row in csv.DictReader(stream)}
    assert float(
        written["gresho_emac_vms_oss_nonlinear"]["energy_closure_scaled"]
    ) == pytest.approx(1e-10)
    assert written["gresho_emac_oss_linear"]["status"] == "incomplete"
    markdown = (output / "gresho.md").read_text()
    assert "| gresho_conservative_oss_linear | completed | 2 |" in markdown
    summary = (output / "summary.md").read_text()
    assert (
        "gresho_visc_skew_symmetric_1_asgs_nonlinear: failed, 0 of 100 steps converged" in summary
    )
    assert "step 1 did not converge: stalled" in summary
    assert "![Long Gresho runs](gresho_long.png)" in summary
    assert (output / "gresho_long.png").stat().st_size > 10_000
    assert [row["run"] for row in tables["gresho_visc_refined"]] == [
        "gresho_visc_refined_emac_vms_asgs"
    ]
    assert [row["run"] for row in tables["gresho_refined"]] == ["gresho_refined_emac_vms_oss"]
    assert {name: len(rows) for name, rows in tables.items()} == {
        "gresho": 12,
        "gresho_visc": 12,
        "gresho_long": 5,
        "gresho_refined": 1,
        "gresho_visc_refined": 1,
        "taylor_green": 24,
    }
    for name in ANALYSIS.FAMILIES:
        assert (output / f"{name}.csv").exists() and (output / f"{name}.md").exists()
