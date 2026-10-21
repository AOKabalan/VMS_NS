"""Stored subscale histories at the quadrature points (spec: Time-step algorithm step 6)."""

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from ufl.core.expr import Expr

from navier_stokes_vms.discretization.time_schemes import StepCoefficients


class SubscaleStorage:
    """ũⁿ, ũⁿ⁻¹, … starting from zero, and their combination ũ_hist for the current step."""

    def __init__(self, space: WithGeometry, levels: int) -> None:
        self._levels = [fd.Function(space) for _ in range(levels)]
        self._solve_level = fd.Function(space)
        self._endpoint = fd.Function(space)
        self.history = fd.Function(space)

    def form_history(self, coefficients: StepCoefficients) -> None:
        """Set ũ_hist = Σ_j history_weights[j] ũ^(n−j) for this step."""
        weighted = zip(coefficients.history_weights, self._levels, strict=False)
        self.history.assign(sum(weight * level for weight, level in weighted))

    def store(self, u_sub: Expr, coefficients: StepCoefficients) -> None:
        """Sample the converged subscale, extrapolate it to t^(n+1) and shift the histories."""
        self._solve_level.interpolate(u_sub)
        self._endpoint.assign(
            coefficients.extrapolation_weight_solve * self._solve_level
            + coefficients.extrapolation_weight_old * self._levels[0]
        )
        for level in range(len(self._levels) - 1, 0, -1):
            self._levels[level].assign(self._levels[level - 1])
        self._levels[0].assign(self._endpoint)
