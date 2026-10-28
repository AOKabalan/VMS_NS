"""Time-scheme coefficients and startup sequences, without Firedrake (spec: Time schemes)."""

from dataclasses import dataclass

from navier_stokes_vms.config import Startup, TimeScheme

_BDF = {
    "bdf1": (1.0, (1.0,)),
    "bdf2": (3 / 2, (2.0, -1 / 2)),
    "bdf3": (11 / 6, (3.0, -3 / 2, 1 / 3)),
}
_BDF_ORDER = {"bdf1": 1, "bdf2": 2, "bdf3": 3}


@dataclass(frozen=True)
class StepCoefficients:
    """One step: D_t g = (α₀ g − Σ_j history_weights[j] g^(n−j))/δt at t^n + solve_level·δt."""

    scheme: TimeScheme
    alpha0: float
    history_weights: tuple[float, ...]
    solve_level: float
    data_weight_new: float
    data_weight_old: float
    extrapolation_weight_solve: float
    extrapolation_weight_old: float


def scheme_at_step(scheme: TimeScheme, startup: Startup, step: int) -> TimeScheme:
    """Scheme used at step n+1 = step (1-based), after the startup of spec: Time schemes."""
    if step < 1:
        raise ValueError(f"step = {step}: steps are numbered from 1")
    if startup == "none":
        if scheme != "crank_nicolson":
            raise ValueError("time.startup = 'none' requires time.scheme = 'crank_nicolson'")
        return scheme
    if scheme == "crank_nicolson":
        return "bdf1" if step == 1 else "crank_nicolson"
    if startup == "exact_history":
        return scheme
    order = min(_BDF_ORDER[scheme], step)
    return next(name for name, value in _BDF_ORDER.items() if value == order)


def step_coefficients(scheme: TimeScheme, startup: Startup, step: int) -> StepCoefficients:
    """α₀, history weights, solve level, data weights and extrapolation for one step."""
    used = scheme_at_step(scheme, startup, step)
    if used == "crank_nicolson":
        return StepCoefficients(
            scheme=used,
            alpha0=2.0,
            history_weights=(2.0,),
            solve_level=1 / 2,
            data_weight_new=1 / 2,
            data_weight_old=1 / 2,
            extrapolation_weight_solve=2.0,
            extrapolation_weight_old=-1.0,
        )
    alpha0, history_weights = _BDF[used]
    return StepCoefficients(
        scheme=used,
        alpha0=alpha0,
        history_weights=history_weights,
        solve_level=1.0,
        data_weight_new=1.0,
        data_weight_old=0.0,
        extrapolation_weight_solve=1.0,
        extrapolation_weight_old=0.0,
    )
