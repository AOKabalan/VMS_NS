"""Inner Newton or Picard solve with its convergence test and iteration records."""

import math
from collections.abc import Sequence
from dataclasses import dataclass

import firedrake as fd
import ufl
from firedrake.exceptions import ConvergenceError
from firedrake.petsc import PETSc
from ufl.core.expr import Expr
from ufl.form import Form
from ufl.measure import Measure

from navier_stokes_vms.config import SolverOptions

_REASON = PETSc.SNES.ConvergedReason


@dataclass(frozen=True)
class InnerIteration:
    """One row of inner_iterations.csv: the state after iteration `iteration` (0 is the guess)."""

    iteration: int
    method: str
    residual_norm: float
    relative_residual: float
    update_norm: float | None
    linear_iterations: int | None


@dataclass(frozen=True)
class InnerSolveResult:
    """Whether the inner solve met its tolerance, and every iteration it took."""

    converged: bool
    iterations: tuple[InnerIteration, ...]
    failure_reason: str | None = None


class InnerSolver:
    """Newton with full steps and Picard on the frozen-advector residual, for one unknown."""

    def __init__(
        self,
        solution: fd.Function,
        newton_residual: Form,
        picard_residual: Form,
        frozen_advector: fd.Function,
        bcs: Sequence[fd.DirichletBC],
        null_space: fd.MixedVectorSpaceBasis | None,
        options: SolverOptions,
        newton_jacobian: Form | None = None,
    ) -> None:
        self._solution = solution
        self._newton_residual = newton_residual
        self._frozen_advector = frozen_advector
        self._bcs = list(bcs)
        self._options = options
        self._previous = fd.Function(solution.function_space())
        self._rows: list[InnerIteration] = []
        self._initial_norm = 0.0
        self._tolerance = 0.0
        self._newton_max_iterations = 0
        self._failure_reason: str | None = None
        self._picard_iterations = (
            options.inner_max_iterations
            if options.linearization == "picard"
            else options.picard_iterations_before_newton
        )
        self._picard = None
        if self._picard_iterations > 0:
            self._picard = self._variational_solver(picard_residual, null_space, "ksponly")
        self._newton = None
        if options.linearization == "newton":
            self._newton = self._variational_solver(
                newton_residual, null_space, "newtonls", newton_jacobian
            )
            self._newton.snes.setMonitor(self._record_newton_iteration)
            self._newton.snes.setConvergenceTest(self._newton_convergence)

    def solve(self) -> InnerSolveResult:
        """Iterate from the current solution until the residual meets the tolerance or the limit."""
        for bc in self._bcs:
            bc.apply(self._solution)
        self._rows = []
        self._failure_reason = None
        self._initial_norm = self.residual_norm()
        self._tolerance = max(
            self._options.inner_atol, self._options.inner_rtol * self._initial_norm
        )
        first_method = "picard" if self._picard is not None else "newton"
        self._append(first_method, self._initial_norm, None, None)
        if not math.isfinite(self._initial_norm):
            return InnerSolveResult(False, tuple(self._rows), self._nonfinite_message("residual"))
        converged = self._initial_norm <= self._tolerance
        if self._picard is not None:
            while not converged and self._iteration() < self._picard_iterations:
                converged = self._picard_iteration(self._picard)
        remaining = self._options.inner_max_iterations - self._iteration()
        if not converged and self._newton is not None and remaining > 0:
            converged = self._newton_iterations(self._newton, remaining)
        return InnerSolveResult(converged, tuple(self._rows), self._failure_reason)

    def residual_norm(self) -> float:
        """ℓ² norm of the full nonlinear residual vector with the Dirichlet rows excluded."""
        residual = fd.assemble(self._newton_residual)
        for bc in self._bcs:
            bc.zero(residual)
        with residual.dat.vec_ro as vector:
            return vector.norm()

    def _picard_iteration(self, solver: fd.NonlinearVariationalSolver) -> bool:
        self._frozen_advector.assign(self._solution.subfunctions[0])
        self._previous.assign(self._solution)
        try:
            solver.solve()
        except ConvergenceError as error:
            self._failure_reason = str(error)
            return False
        residual = self.residual_norm()
        with self._solution.dat.vec_ro as new, self._previous.dat.vec as update:
            update.aypx(-1.0, new)
            update_norm = update.norm()
        linear = solver.snes.getKSP().getIterationNumber()
        self._append("picard", residual, update_norm, linear)
        return residual <= self._tolerance

    def _newton_iterations(self, solver: fd.NonlinearVariationalSolver, limit: int) -> bool:
        self._newton_max_iterations = limit
        try:
            solver.solve()
        except PETSc.Error:
            # petsc4py wraps exceptions raised by the Jacobian callback.
            if self._failure_reason is None:
                raise
            return False
        except ConvergenceError as error:
            if self._failure_reason is None:
                self._failure_reason = str(error)
            return False
        return solver.snes.getConvergedReason() > 0

    def _record_newton_iteration(self, snes: PETSc.SNES, iteration: int, norm: float) -> None:
        if iteration == 0:
            return
        update = snes.getSolutionUpdate().norm()
        linear = snes.getKSP().getIterationNumber()
        self._append("newton", norm, update, linear)

    def _newton_convergence(
        self, _snes: PETSc.SNES, iteration: int, norms: tuple[float, float, float]
    ) -> int:
        norm = norms[2]
        if not math.isfinite(norm):
            self._failure_reason = self._nonfinite_message("residual")
            return _REASON.DIVERGED_FNORM_NAN
        if norm <= self._tolerance:
            return _REASON.CONVERGED_FNORM_ABS
        if iteration >= self._newton_max_iterations:
            return _REASON.DIVERGED_MAX_IT
        return _REASON.CONVERGED_ITERATING

    def _variational_solver(
        self,
        residual: Form,
        null_space: fd.MixedVectorSpaceBasis | None,
        snes_type: str,
        jacobian: Form | None = None,
    ) -> fd.NonlinearVariationalSolver:
        parameters = {
            **self._options.petsc,
            "snes_type": snes_type,
            "snes_linesearch_type": "basic",
        }
        problem = fd.NonlinearVariationalProblem(residual, self._solution, self._bcs, J=jacobian)
        return fd.NonlinearVariationalSolver(
            problem,
            nullspace=null_space,
            transpose_nullspace=null_space,
            solver_parameters=parameters,
            post_jacobian_callback=(
                self._check_jacobian
                if snes_type == "newtonls" and self._options.newton_tau == "exact"
                else None
            ),
        )

    def _check_jacobian(self, _state: PETSc.Vec, matrix: PETSc.Mat) -> None:
        if not math.isfinite(matrix.norm(PETSc.NormType.FROBENIUS)):
            self._failure_reason = self._nonfinite_message("Jacobian")
            raise ConvergenceError(self._failure_reason)

    def _nonfinite_message(self, quantity: str) -> str:
        message = f"non-finite {quantity} (NaN or infinity)"
        if self._options.newton_tau == "exact":
            message += (
                "; solver.newton_tau = 'exact' differentiates |a|, which is undefined at a = 0"
            )
        return message

    def _append(self, method: str, norm: float, update: float | None, linear: int | None) -> None:
        relative = norm / self._initial_norm if self._initial_norm > 0 else 0.0
        self._rows.append(InnerIteration(len(self._rows), method, norm, relative, update, linear))

    def _iteration(self) -> int:
        return len(self._rows) - 1


def relative_change(new: Expr, old: Expr, dx: Measure) -> float:
    """‖new − old‖ / ‖new‖ in L² with the common rule; the absolute change when new is zero."""
    change = math.sqrt(fd.assemble(ufl.inner(new - old, new - old) * dx))
    size = math.sqrt(fd.assemble(ufl.inner(new, new) * dx))
    return change / size if size > 0 else change
