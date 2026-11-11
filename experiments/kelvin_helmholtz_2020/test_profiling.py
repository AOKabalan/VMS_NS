"""Check nested accounting without running a fluid simulation."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from kh_profiling import PhaseProfiler


def test_nested_times_are_not_double_counted() -> None:
    ticks = iter((0.0, 2.0, 9.0, 12.0))
    profiler = PhaseProfiler(lambda: next(ticks))
    profiler.step = 3
    child = profiler.timed("child", lambda value: value + 1)
    parent = profiler.timed("parent", child)
    assert parent(value=41) == 42
    rows = profiler.rows[3]
    assert rows["parent"]["inclusive_seconds"] == 12
    assert rows["parent"]["exclusive_seconds"] == 5
    assert rows["child"]["exclusive_seconds"] == 7
    assert sum(row["exclusive_seconds"] for row in rows.values()) == 12
    assert not profiler._stack


def test_exception_unwinds_all_timers() -> None:
    ticks = iter((0.0, 2.0, 9.0, 12.0))
    profiler = PhaseProfiler(lambda: next(ticks))

    def fail() -> None:
        raise ValueError("intentional")

    call = profiler.timed("parent", profiler.timed("child", fail))
    with pytest.raises(ValueError, match="intentional"):
        call()
    assert not profiler._stack
    assert profiler.rows[0]["parent"]["exclusive_seconds"] == 5


def test_restoration_and_serial_export(tmp_path: Path) -> None:
    profiler = PhaseProfiler()
    target = SimpleNamespace(value="original")
    profiler._replace(target, "value", "first")
    profiler._replace(target, "value", "second")
    comm = SimpleNamespace(rank=0, gather=lambda value, root: [value] if root == 0 else None)
    path = tmp_path / "profile.json"
    profiler.close(comm, path)
    assert target.value == "original"
    assert not profiler._restore
    assert json.loads(path.read_text())["ranks"] == [{}]


def test_profiled_solver_preserves_uniform_flow_and_restores_checks(tmp_path: Path) -> None:
    """Exercise the profiler against the current named-field solver interface."""
    import logging
    from dataclasses import replace

    import firedrake as fd
    from firedrake.petsc import PETSc
    from run import SnapshotStepper

    from navier_stokes_vms.cases import build_problem
    from navier_stokes_vms.config import Config, MeshOptions, OutputOptions, TimeOptions

    config = Config(
        case="kelvin_helmholtz",
        run_name="profile_smoke",
        mesh=MeshOptions(cells=(4, 4)),
        time=TimeOptions(dt=0.01, t_end=0.02),
        output=OutputOptions(vtk_every=0),
    )
    problem = build_problem(config)
    problem = replace(problem, initial_velocity=lambda _x: fd.Constant((1.0, 0.0)))
    stepper = SnapshotStepper(config, problem, tmp_path, logging.getLogger(__name__), [])
    original_checks = stepper._step_solver._changes
    PETSc.Log.begin()
    profiler = PhaseProfiler()
    profiler.install(stepper)
    path = tmp_path / "phase_profile.json"
    try:
        stepper.run(stepper.record)
    finally:
        profiler.close(problem.mesh.comm, path)
        stepper.enstrophy.close()
    assert stepper._step_solver._changes is original_checks
    error = stepper.velocity - fd.Constant((1.0, 0.0))
    assert fd.assemble(fd.inner(error, error) * fd.dx) < 1e-18
    rows = json.loads(path.read_text())["ranks"][0]
    assert set(rows) == {"0", "1", "2"}
    assert rows["2"]["change_norms"]["calls"] > 0
