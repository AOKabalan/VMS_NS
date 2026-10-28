"""Optional nested wall timings and PETSc stages for the KH experiment only."""

import json
import time
from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from mpi4py import MPI

    from navier_stokes_vms.solver.time_stepper import TimeStepper


class PhaseProfiler:
    """Measure inclusive/exclusive local wall time without barriers inside the solve.

    Timers nest, so their exclusive times can be added within a rank. MPI rank
    maxima of separate phases must not be summed as if they were one timeline.
    Instrumentation is removed when writing the result, including on failure.
    """

    def __init__(self, clock: Callable[[], float] = time.perf_counter) -> None:
        self.step = 0
        self._clock = clock
        self._stack: list[list[float]] = []
        self.rows: dict[int, dict[str, dict[str, float]]] = {}
        self._restore: list[tuple[object, str, object]] = []

    def timed(self, name: str, function: Callable[..., object]) -> Callable[..., object]:
        """Wrap one callable, tracking child time and retaining its return/exception."""

        @wraps(function)
        def wrapped(*args: object, **kwargs: object) -> object:
            step = self.step
            frame = [self._clock(), 0.0]
            self._stack.append(frame)
            try:
                return function(*args, **kwargs)
            finally:
                elapsed = self._clock() - frame[0]
                self._stack.pop()
                if self._stack:
                    self._stack[-1][1] += elapsed
                row = self.rows.setdefault(step, {}).setdefault(
                    name, {"calls": 0, "inclusive_seconds": 0.0, "exclusive_seconds": 0.0}
                )
                row["calls"] += 1
                row["inclusive_seconds"] += elapsed
                row["exclusive_seconds"] += elapsed - frame[1]

        return wrapped

    def _replace(self, target: object, attribute: str, value: object) -> None:
        self._restore.append((target, attribute, getattr(target, attribute)))
        setattr(target, attribute, value)

    def install(self, stepper: "TimeStepper") -> None:
        """Instrument an existing stepper; PETSc logging must already be enabled."""
        from firedrake.petsc import PETSc

        import navier_stokes_vms.solver.time_stepper as module

        stages = (PETSc.Log.Stage("KH_initial_steps"), PETSc.Log.Stage("KH_measured_steps"))

        def staged(function: Callable[..., object]) -> Callable[..., object]:
            @wraps(function)
            def wrapped(step: int, *args: object, **kwargs: object) -> object:
                self.step = step
                stage = stages[int(step > 2)]
                stage.push()
                try:
                    return function(step, *args, **kwargs)
                finally:
                    stage.pop()

            return wrapped

        outer = stepper._step_solver
        inner = outer._inner
        for target, attribute, label in (
            (stepper, "_step", "time_step"),
            (stepper, "record", "enstrophy_and_snapshots"),
            (stepper, "_set_level", "time_level"),
            (stepper, "_normalize_pressure", "pressure_normalization"),
            (stepper, "_advance", "history_update"),
            (outer, "solve", "nonlinear_solve"),
            (outer, "_row", "outer_checks"),
            (inner, "solve", "inner_solve"),
            (inner, "_newton_iterations", "snes_solve"),
            (inner, "residual_norm", "explicit_residual_norm"),
            (module, "evaluate_diagnostics", "diagnostics"),
        ):
            self._replace(target, attribute, self.timed(label, getattr(target, attribute)))
        if outer._subscale is not None:
            self._replace(
                outer._subscale, "update", self.timed("local_subscale", outer._subscale.update)
            )
            self._replace(
                outer, "_subscale_change", self.timed("change_norms", outer._subscale_change)
            )
        self._replace(
            outer, "_changes", [self.timed("change_norms", change) for change in outer._changes]
        )
        self._replace(stepper, "_step", staged(stepper._step))
        self._replace(stepper, "record", staged(stepper.record))

    def close(self, comm: "MPI.Comm", path: Path) -> None:
        """Restore methods and gather per-rank data once, after the simulation."""
        for target, attribute, original in reversed(self._restore):
            setattr(target, attribute, original)
        self._restore.clear()
        ranks = comm.gather(self.rows, root=0)
        if comm.rank == 0:
            path.write_text(
                json.dumps(
                    {
                        "measured_steps": "steps >= 3; initial diagnostics are step 0",
                        "timing": "wall seconds; inclusive events nest; exclusive add per rank",
                        "ranks": ranks,
                    },
                    indent=2,
                )
                + "\n"
            )
