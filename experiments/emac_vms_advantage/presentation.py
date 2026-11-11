"""Presentation figures from the EMAC–VMS experiments: python presentation.py."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).parent
RESULTS = HERE / "results"
SCREEN = HERE / "screen" / "results" / "gresho"
DTSWEEP = HERE / "dtsweep" / "results" / "gresho"
OUT = HERE / "presentation"
FORMS = ("nonconservative", "conservative", "skew_symmetric_1", "emac", "emac_vms")
SHORT = {
    "nonconservative": "CONV",
    "conservative": "CONS",
    "skew_symmetric_1": "SKEW",
    "skew_symmetric_2": "SKEW2",
    "emac": "EMAC",
    "emac_vms": "EMAC–VMS",
    "emac_vms_full": "EMAC–VMS full",
}
COLOR = {
    "nonconservative": "#1f6fb2",
    "conservative": "#e0a020",
    "skew_symmetric_1": "#169c78",
    "emac": "#a0609a",
    "emac_vms": "#d4570f",
    "emac_vms_full": "#7a3a10",
}
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})


def rows(path: Path) -> list[dict[str, str]]:
    with path.open() as stream:
        return list(csv.DictReader(stream))


def column(data: list[dict[str, str]], key: str) -> np.ndarray:
    return np.array([float(r[key]) if r[key] not in ("", None) else math.nan for r in data])


def save(figure: plt.Figure, name: str) -> None:
    figure.tight_layout()
    figure.savefig(OUT / f"{name}.png", dpi=160)
    plt.close(figure)


def defect_matrix() -> None:
    """Three structural defects per form: only EMAC–VMS is at round-off in all three."""
    energy, momentum, torque = [], [], []
    for form in FORMS:
        e3 = rows(RESULTS / f"E3_{form}_oss" / "diagnostics.csv")
        e1 = rows(RESULTS / f"E1_{form}_oss" / "extra.csv")
        e0 = column(e3, "kinetic_energy")[0]
        a0 = abs(column(e3, "angular_momentum")[0])
        energy.append(abs(np.nansum(column(e3, "convective_energy_defect")) * 0.01) / e0)
        torque.append(abs(np.nansum(column(e3, "angular_momentum_defect")) * 0.01) / a0)
        momentum.append(np.nanmax(np.abs(column(e1, "momentum_x") - 1.0)))
    floor = 1e-17
    figure, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for axis, values, title in zip(
        axes,
        (energy, momentum, torque),
        (
            "Spurious convective energy\n|∫N dt| / E(0)  (box vortex, T = 10)",
            "Linear-momentum drift\nmax |M(t) − M(0)| / M(0)  (translating vortex)",
            "Spurious convective torque\n|∫D_φ dt| / A(0)  (box vortex, T = 10)",
        ),
        strict=True,
    ):
        bars = np.maximum(values, floor)
        axis.bar([SHORT[f] for f in FORMS], bars, color=[COLOR[f] for f in FORMS])
        axis.set_yscale("log")
        axis.set_ylim(floor / 2, 1)
        axis.axhspan(floor / 2, 1e-12, color="#ddd", alpha=0.6, zorder=0)
        axis.set_title(title, fontsize=10)
        axis.tick_params(axis="x", rotation=30)
    axes[0].set_ylabel("relative defect (log scale)")
    axes[0].text(-0.4, 3e-15, "round-off", fontsize=9, color="#555")
    figure.suptitle(
        "Only EMAC–VMS has all three convective defects at round-off (P1/P1 dynamic OSS)",
        fontweight="bold",
    )
    save(figure, "1_defect_matrix")


def emac_vs_emac_vms() -> None:
    """Standard EMAC (fixed coupling) against EMAC–VMS on the translating vortex."""
    figure, axes = plt.subplots(1, 2, figsize=(11, 4))
    for form in ("emac", "emac_vms", "skew_symmetric_1"):
        e = rows(RESULTS / f"E1_{form}_oss" / "extra.csv")
        d = rows(RESULTS / f"E1_{form}_oss" / "diagnostics.csv")
        t = column(e, "time")
        lw = 2.6 if form == "emac_vms" else 1.6
        axes[0].plot(t, column(e, "momentum_x") - 1.0, color=COLOR[form], lw=lw, label=SHORT[form])
        n = np.cumsum(np.nan_to_num(column(d, "convective_energy_defect"))) * 0.01
        axes[1].plot(column(d, "time"), n / column(e, "energy")[0], color=COLOR[form], lw=lw)
    axes[0].set_title("Momentum drift M_x(t) − M_x(0)", loc="left", fontweight="bold")
    axes[1].set_title("Energy injected by convection ∫N dt / E(0)", loc="left", fontweight="bold")
    for axis in axes:
        axis.set_xlabel("time t")
        axis.grid(alpha=0.3)
    axes[0].legend()
    figure.suptitle(
        "Translating Gresho vortex (periodic, U₀ = (1, ½)): the VMS coupling, "
        "not EMAC itself, decides conservation",
        fontweight="bold",
        fontsize=11,
    )
    save(figure, "2_emac_vs_emac_vms")


def galerkin_anchor() -> None:
    """Without stabilization the papers' blow-up is reproduced."""
    figure, axis = plt.subplots(figsize=(6.5, 4))
    for form in FORMS[:-1]:
        d = rows(SCREEN / f"galerkin_th_{form}_48" / "diagnostics.csv")
        e = column(d, "kinetic_energy")
        axis.plot(column(d, "time"), e / e[0], color=COLOR[form], lw=2, label=SHORT[form])
    axis.set_ylim(0.9, 2.0)
    axis.set_xlabel("time t")
    axis.set_ylabel("E(t)/E(0)")
    axis.legend()
    axis.grid(alpha=0.3)
    axis.set_title(
        "Galerkin Taylor–Hood, Gresho 48², ν = 1e-6 (Li & Rui 2023 setup)\n"
        "CONV, CONS blow up; SKEW fills with grid-scale noise (‖∇·u_h‖ → 90)\n"
        "that even ν = 1e-6 dissipates; EMAC stays clean. With VMS no form blows up.",
        fontsize=9,
    )
    save(figure, "3_galerkin_blowup_anchor")


def smooth_vs_gresho() -> None:
    """Accuracy at t = 3 in the no-slip box: Gresho (coarse, fine) against a smooth vortex."""
    cases = (
        ("Gresho vortex, 48²", "E3_{}_oss", 299),
        ("Gresho vortex, 96²", "E3_{}_oss_96", -1),
        ("Smooth vortex, 48²", "E5_{}_oss", -1),
    )
    figure, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for axis, (title, pattern, index) in zip(axes, cases, strict=True):
        errors = [
            column(rows(RESULTS / pattern.format(form) / "diagnostics.csv"), "velocity_error_l2")[
                index
            ]
            for form in FORMS
        ]
        axis.bar([SHORT[f] for f in FORMS], errors, color=[COLOR[f] for f in FORMS])
        axis.set_title(title, loc="left", fontweight="bold")
        axis.tick_params(axis="x", rotation=30)
    axes[0].set_ylabel("L² velocity error at t = 3")
    figure.suptitle(
        "Vortex in a no-slip box, P1/P1 OSS: EMAC forms are even on coarse Gresho, "
        "worse on refined Gresho (kinked |u|²), best on a smooth vortex",
        fontweight="bold",
        fontsize=11,
    )
    save(figure, "4_smooth_vs_gresho_accuracy")


def dissipation_budget() -> None:
    """The extra EMAC–VMS energy loss is velocity-subscale dissipation and independent of δt."""
    figure, axis = plt.subplots(figsize=(8, 4))
    dts = ("0.02", "0.01", "0.005", "0.0025")
    width = 0.35
    for k, form in enumerate(("skew_symmetric_1", "emac_vms")):
        sub, pres = [], []
        for dt in dts:
            d = rows(DTSWEEP / f"dt_{form}_{dt}" / "diagnostics.csv")
            step = float(d[1]["time"]) - float(d[0]["time"])
            e0 = float(d[0]["kinetic_energy"])
            sub.append(np.nansum(column(d, "energy_budget_subscale")) * step / e0)
            pres.append(np.nansum(column(d, "energy_budget_pressure_subscale")) * step / e0)
        x = np.arange(len(dts)) + (k - 0.5) * width
        axis.bar(x, sub, width, color=COLOR[form], label=f"{SHORT[form]}: velocity subscale")
        axis.bar(
            x,
            pres,
            width,
            bottom=sub,
            color=COLOR[form],
            alpha=0.35,
            label=f"{SHORT[form]}: pressure subscale",
        )
    axis.set_xticks(np.arange(len(dts)), [f"δt = {dt}" for dt in dts])
    axis.set_ylabel("energy dissipated by t = 2, / E(0)")
    axis.set_ylim(0, 0.105)
    axis.legend(fontsize=8, ncol=2)
    axis.set_title(
        "Gresho 32², CN: extra EMAC–VMS loss is all velocity-subscale "
        "dissipation, independent of δt",
        fontsize=10,
    )
    save(figure, "5_dissipation_budget_dt_sweep")


def full_variant() -> None:
    """Minimal against full-stress EMAC–VMS (chapter §5.8.5)."""
    cases = [
        ("Translating Gresho\n(32², T = 4)", "E1_{}_oss", "extra", "energy", "error_l2"),
        (
            "Box Gresho\n(96², T = 3)",
            "E3_{}_oss_96",
            "diagnostics",
            "kinetic_energy",
            "velocity_error_l2",
        ),
        (
            "Box Gresho\n(48², T = 10)",
            "E3_{}_oss",
            "diagnostics",
            "kinetic_energy",
            "velocity_error_l2",
        ),
    ]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(len(cases))
    for k, form in enumerate(("emac_vms", "emac_vms_full")):
        loss, error = [], []
        for _, pattern, table, energy_key, error_key in cases:
            d = rows(RESULTS / pattern.format(form) / f"{table}.csv")
            e = column(d, energy_key)
            loss.append(1 - e[-1] / e[0])
            error.append(column(d, error_key)[-1])
        offset = (k - 0.5) * 0.35
        axes[0].bar(x + offset, loss, 0.35, color=COLOR[form], label=SHORT[form])
        axes[1].bar(x + offset, error, 0.35, color=COLOR[form], label=SHORT[form])
    for axis, title in zip(
        axes, ("energy lost, 1 − E(T)/E(0)", "L² velocity error at T"), strict=True
    ):
        axis.set_xticks(x, [c[0] for c in cases], fontsize=9)
        axis.set_title(title, loc="left", fontweight="bold")
    axes[0].legend()
    figure.suptitle(
        "Full-stress variant: same exact conservation, more dissipation, no accuracy gain",
        fontweight="bold",
        fontsize=11,
    )
    save(figure, "6_full_stress_variant")


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for make in (
        defect_matrix,
        emac_vs_emac_vms,
        galerkin_anchor,
        smooth_vs_gresho,
        dissipation_budget,
        full_variant,
    ):
        make()
    print(sorted(path.name for path in OUT.glob("*.png")))  # noqa: T201 — command-line summary


if __name__ == "__main__":
    main()
