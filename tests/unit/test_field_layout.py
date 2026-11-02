"""Mixed field identities and convergence reporting independent of component count."""

from unittest.mock import Mock

import firedrake as fd
import pytest
import ufl

from navier_stokes_vms.config import SolverOptions, StabilizationOptions
from navier_stokes_vms.discretization.function_spaces import (
    FieldLayout,
    flow_layout,
    mixed_space,
    oss_space,
)
from navier_stokes_vms.solver.nonlinear_solver import (
    InnerIteration,
    InnerSolver,
    InnerSolveResult,
    StepSolver,
)
from navier_stokes_vms.solver.subscale_storage import LocalUpdateResult, SubscaleIterate


@pytest.mark.parametrize(
    ("method", "subscale_space", "pressure_projection", "names"),
    [
        ("none", "oss", True, ("velocity", "pressure")),
        ("vms", "asgs", True, ("velocity", "pressure")),
        ("vms", "oss", False, ("velocity", "pressure", "xi")),
        ("vms", "oss", True, ("velocity", "pressure", "xi", "zeta")),
    ],
)
def test_flow_layout_matches_space(
    method: str, subscale_space: str, pressure_projection: bool, names: tuple[str, ...]
) -> None:
    options = StabilizationOptions(
        method=method, subscale_space=subscale_space, pressure_subscale=pressure_projection
    )
    layout = flow_layout(options)
    mesh = fd.UnitSquareMesh(2, 2)
    space = (
        oss_space(mesh, 1, pressure_projection)
        if method == "vms" and subscale_space == "oss"
        else mixed_space(mesh, 2 if method == "none" else 1, 1)
    )
    solution = fd.Function(space)
    assert layout.names == names
    assert layout.index("pressure") == 1
    for components in (solution.subfunctions, ufl.split(solution), fd.TestFunctions(space)):
        fields = layout.bind(components)
        assert tuple(fields) == names
        assert all(
            fields[name] is component for name, component in zip(names, components, strict=True)
        )
    with pytest.raises(ValueError):
        layout.bind(solution.subfunctions[:-1])


@pytest.mark.parametrize(
    "names",
    [
        ("pressure", "velocity"),
        ("velocity",),
        ("velocity", "pressure", "pressure"),
        ("velocity", "pressure", ""),
    ],
)
def test_invalid_layout(names: tuple[str, ...]) -> None:
    with pytest.raises(ValueError):
        FieldLayout(names)


def test_step_solver_names_changes_with_an_intervening_field() -> None:
    mesh = fd.UnitSquareMesh(2, 2)
    velocity = fd.VectorFunctionSpace(mesh, "CG", 1)
    scalar = fd.FunctionSpace(mesh, "CG", 1)
    space = fd.MixedFunctionSpace([velocity, scalar, scalar, velocity])
    layout = FieldLayout(("velocity", "pressure", "marker", "xi"))
    solution = fd.Function(space).assign(1.0)
    inner = Mock(spec=InnerSolver)
    result = InnerSolveResult(True, (InnerIteration(0, "newton", 0.0, 0.0, None, None),))

    def update() -> InnerSolveResult:
        for component, value in zip(solution.subfunctions, (2.0, 3.0, 5.0, 7.0), strict=True):
            component.assign(value)
        return result

    inner.solve.side_effect = update
    solver = StepSolver(inner, solution, None, fd.dx(domain=mesh), SolverOptions(), layout)
    row = solver.solve().iterations[0]
    assert tuple(row.field_changes) == layout.names
    assert row.velocity_change == pytest.approx(1 / 2)
    assert row.pressure_change == pytest.approx(2 / 3)
    assert row.field_changes["marker"] == pytest.approx(4 / 5)
    assert row.xi_change == pytest.approx(6 / 7)
    assert row.zeta_change is None
    with pytest.raises(TypeError):
        row.field_changes["xi"] = 0.0


def test_additional_field_participates_in_outer_convergence() -> None:
    mesh = fd.UnitSquareMesh(2, 2)
    velocity = fd.VectorFunctionSpace(mesh, "CG", 1)
    scalar = fd.FunctionSpace(mesh, "CG", 1)
    space = fd.MixedFunctionSpace([velocity, scalar, scalar])
    layout = FieldLayout(("velocity", "pressure", "marker"))
    solution = fd.Function(space).assign(1.0)
    inner = Mock(spec=InnerSolver)
    result = InnerSolveResult(True, (InnerIteration(0, "newton", 1.0, 1.0, None, None),))

    def update() -> InnerSolveResult:
        solution.subfunctions[2].assign(2.0)
        return result

    inner.solve.side_effect = update
    inner.residual_norm.return_value = 0.0
    subscale = Mock(spec=SubscaleIterate)
    subscale.value = fd.Function(velocity)
    subscale.update.return_value = LocalUpdateResult(True, 1)
    solver = StepSolver(inner, solution, subscale, fd.dx(domain=mesh), SolverOptions(), layout)
    solved = solver.solve()
    assert solved.converged
    assert len(solved.iterations) == 2
    assert not solved.iterations[0].converged
    assert solved.iterations[0].velocity_change == 0.0
    assert solved.iterations[0].pressure_change == 0.0
    assert solved.iterations[0].field_changes["marker"] == pytest.approx(0.5)
    assert solved.iterations[1].field_changes["marker"] == 0.0
