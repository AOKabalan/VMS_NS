"""Resolved field histories and endpoint recovery shared by velocity and scalars."""

import firedrake as fd
from firedrake.functionspaceimpl import WithGeometry
from ufl.core.expr import Expr

from navier_stokes_vms.discretization.time_schemes import StepCoefficients


class FieldHistory:
    """Persistent endpoint states, solve-level history and endpoint workspace of one field."""

    def __init__(self, space: WithGeometry, levels: int, name: str) -> None:
        self.levels = [fd.Function(space, name=name) for _ in range(levels)]
        self.history = fd.Function(space)
        self.endpoint = fd.Function(space, name=name)

    def form_history(self, coefficients: StepCoefficients) -> None:
        """Combine previous endpoints with the integrator's actual history weights."""
        weighted = zip(coefficients.history_weights, self.levels, strict=False)
        self.history.assign(sum(weight * value for weight, value in weighted))

    def advance(self, value: Expr, coefficients: StepCoefficients) -> None:
        """Recover the endpoint and shift histories only after a successful coupled solve."""
        self.endpoint.assign(
            coefficients.extrapolation_weight_solve * value
            + coefficients.extrapolation_weight_old * self.levels[0]
        )
        for level in range(len(self.levels) - 1, 0, -1):
            self.levels[level].assign(self.levels[level - 1])
        self.levels[0].assign(self.endpoint)
