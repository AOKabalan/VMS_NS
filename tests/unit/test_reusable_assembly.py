"""Assemblers preprocessed once give exactly what fd.assemble and interpolate give.

`reusable_assembler` relies on Firedrake internals (`BaseFormAssembler.preprocess_base_form` and
the `is_base_form_preprocessed` flag of `get_assembler`), so this test guards a Firedrake upgrade.
"""

import logging
import math
from pathlib import Path

import firedrake as fd
import numpy as np
import pytest
import ufl

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import (
    Config,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
)
from navier_stokes_vms.discretization.time_schemes import step_coefficients
from navier_stokes_vms.solver.nonlinear_solver import RelativeChange
from navier_stokes_vms.solver.reusable_assembly import reusable_assembler
from navier_stokes_vms.solver.time_stepper import TimeStepper

LOGGER = logging.getLogger(__name__)


def values(tensor: fd.Function | fd.Cofunction) -> np.ndarray:
    with tensor.dat.vec_ro as vector:
        return vector.getArray(readonly=True).copy()


@pytest.mark.parametrize("space", ["asgs", "oss"])
def test_reusable_assembly_matches_fd_assemble_and_interpolate(tmp_path: Path, space: str) -> None:
    config = Config(
        case="kovasznay",
        run_name="cache",
        mesh=MeshOptions(cells=(3, 4)),
        discretization=DiscretizationOptions(velocity_degree=2, pressure_degree=2),
        time=TimeOptions(scheme="bdf2", dt=0.5, t_end=1.0),
        stabilization=StabilizationOptions(
            subscale_space=space, subscale_dynamics="dynamic", scale_splitting="nonlinear"
        ),
        solver=SolverOptions(subscale_update="fixed_point", subscale_max_iterations=1),
        output=OutputOptions(vtk_every=0),
    )
    stepper = TimeStepper(config, build_problem(config), tmp_path, LOGGER)
    stepper._initialise()
    stepper._set_level(2, step_coefficients("bdf2", "lower_order", 2))
    solution, iterate, storage = stepper._solution, stepper._iterate, stepper._subscales
    previous = fd.Function(solution.function_space())
    residual_form = stepper._step_solver._inner._newton_residual
    residual = reusable_assembler(residual_form)
    dx = stepper._dx
    pairs = list(zip(solution.subfunctions, previous.subfunctions, strict=True))
    changes = [RelativeChange(new, old, dx) for new, old in pairs]
    rng = np.random.default_rng(11)
    # Redrawing the inputs between checks shows that the cached assemblers see new values.
    for _ in range(2):
        fields = (*solution.subfunctions, *previous.subfunctions, iterate.value, storage.history)
        for field in fields:
            field.dat.data[:] = rng.uniform(-1.0, 1.0, field.dat.data.shape)
        cached = residual.assemble(tensor=fd.Cofunction(solution.function_space().dual()))
        assert np.array_equal(values(cached), values(fd.assemble(residual_form)))
        for change, (new, old) in zip(changes, pairs, strict=True):
            direct = math.sqrt(fd.assemble(ufl.inner(new - old, new - old) * dx)) / math.sqrt(
                fd.assemble(ufl.inner(new, new) * dx)
            )
            assert change() == direct
        # Fixed point: the update is the subscale expression with the current ũ*.
        expected = fd.Function(iterate.value.function_space()).interpolate(stepper._u_sub)
        iterate.update()
        assert np.array_equal(values(iterate.value), values(expected))
        expected = fd.Function(iterate.value.function_space()).interpolate(stepper._u_sub)
        storage.store(stepper._u_sub, step_coefficients("bdf2", "lower_order", 2))
        assert np.array_equal(values(storage.latest), values(expected))
    assert len(changes) == (4 if space == "oss" else 2)
