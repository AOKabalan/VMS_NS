"""Stored subscales at the quadrature points: histories and the iterate ũ* (spec: steps 3, 6)."""

from dataclasses import dataclass

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from firedrake.petsc import PETSc
from ufl.core.expr import Expr

from navier_stokes_vms.discretization.time_schemes import StepCoefficients
from navier_stokes_vms.solver.reusable_assembly import reusable_assembler


class SubscaleStorage:
    """ũⁿ, ũⁿ⁻¹, … starting from zero, and their combination ũ_hist for the current step."""

    def __init__(self, space: WithGeometry, levels: int) -> None:
        self._space = space
        self._levels = [fd.Function(space) for _ in range(levels)]
        self._solve_level = fd.Function(space)
        self._sampled: Expr | None = None
        self._endpoint = fd.Function(space)
        self.history = fd.Function(space)

    @property
    def latest(self) -> fd.Function:
        """ũⁿ, the subscale at the start of the current step."""
        return self._levels[0]

    def form_history(self, coefficients: StepCoefficients) -> None:
        """Set ũ_hist = Σ_j history_weights[j] ũ^(n−j) for this step."""
        weighted = zip(coefficients.history_weights, self._levels, strict=False)
        self.history.assign(sum(weight * level for weight, level in weighted))

    def store(self, u_sub: Expr, coefficients: StepCoefficients) -> None:
        """Sample the converged subscale, extrapolate it to t^(n+1) and shift the histories."""
        if u_sub is not self._sampled:
            # The time stepper passes the same expression every step: preprocess it once.
            self._sample = reusable_assembler(fd.interpolate(u_sub, self._space))
            self._sampled = u_sub
        self._sample.assemble(tensor=self._solve_level)
        self._endpoint.assign(
            coefficients.extrapolation_weight_solve * self._solve_level
            + coefficients.extrapolation_weight_old * self._levels[0]
        )
        for level in range(len(self._levels) - 1, 0, -1):
            self._levels[level].assign(self._levels[level - 1])
        self._levels[0].assign(self._endpoint)


@dataclass(frozen=True)
class LocalUpdateResult:
    """Whether the pointwise subscale iteration met its tolerance, and how many sweeps it took."""

    converged: bool
    iterations: int


class SubscaleIterate:
    """The stored iterate ũ* inside a = u_h + ũ* and its pointwise update (spec: step 3)."""

    def __init__(
        self,
        value: fd.Function,
        update: Expr,
        rtol: float,
        max_iterations: int,
        is_static: bool = False,
    ) -> None:
        self.value = value
        self.is_static = is_static
        self._rtol = rtol
        self._max_iterations = max_iterations
        self._candidate = fd.Function(value.function_space())
        # The update is preprocessed once and evaluated every sweep.
        self._update = reusable_assembler(fd.interpolate(update, value.function_space()))

    def update(self) -> LocalUpdateResult:
        """Repeat ũ* ← update(ũ*) until max|change| ≤ rtol · max|ũ*| over all points."""
        for iteration in range(1, self._max_iterations + 1):
            self._update.assemble(tensor=self._candidate)
            with self._candidate.dat.vec_ro as new, self.value.dat.vec as old:
                size = new.norm(PETSc.NormType.NORM_INFINITY)
                old.aypx(-1.0, new)
                change = old.norm(PETSc.NormType.NORM_INFINITY)
            self.value.assign(self._candidate)
            if change <= self._rtol * size:
                return LocalUpdateResult(True, iteration)
        return LocalUpdateResult(False, self._max_iterations)
