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
