"""Run configuration: option dataclasses, TOML loading, validation and resolved dump (plan §7)."""

import math
import re
import tomllib
from collections.abc import Collection, Mapping
from dataclasses import dataclass, field, fields
from pathlib import Path
from types import MappingProxyType, UnionType
from typing import Any, Literal, get_args, get_origin, get_type_hints

import tomli_w

MeshSource = Literal["builtin", "gmsh"]
TimeScheme = Literal["bdf1", "bdf2", "bdf3", "crank_nicolson"]
Startup = Literal["lower_order", "exact_history"]
Method = Literal["none", "vms"]
SubscaleSpace = Literal["asgs", "oss"]
SubscaleDynamics = Literal["static", "dynamic"]
ScaleSplitting = Literal["linear", "nonlinear"]
ConvectiveForm = Literal[
    "nonconservative", "conservative", "skew_symmetric_1", "skew_symmetric_2", "emac", "emac_vms"
]
ViscousForm = Literal["laplacian", "stress"]
Linearization = Literal["newton", "picard"]
NewtonTau = Literal["frozen", "exact"]
SubscaleUpdate = Literal["fixed_point", "accelerated"]
LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]
PetscValue = str | int | float | bool


@dataclass(frozen=True)
class MeshOptions:
    """Overrides of the case's default mesh; None keeps the case's choice (plan §8)."""

    source: MeshSource | None = None
    cells: tuple[int, ...] | None = None
    file: str | None = None


@dataclass(frozen=True)
class DiscretizationOptions:
    """Polynomial degrees of the velocity and pressure spaces."""

    velocity_degree: int = 1
    pressure_degree: int = 1


@dataclass(frozen=True)
class TimeOptions:
    """Time scheme, uniform step, final time and startup mode (spec: Time schemes)."""

    scheme: TimeScheme = "bdf2"
    dt: float = 0.01
    t_end: float = 8.0
    startup: Startup = "lower_order"


@dataclass(frozen=True)
class StabilizationOptions:
    """Galerkin or VMS, the three VMS switches and the τ constants (spec: The eight cases)."""

    method: Method = "vms"
    subscale_space: SubscaleSpace = "oss"
    subscale_dynamics: SubscaleDynamics = "dynamic"
    scale_splitting: ScaleSplitting = "nonlinear"
    c1: float = 4.0
    c2: float = 2.0
    pressure_subscale: bool = True


@dataclass(frozen=True)
class ConvectionOptions:
    """Convective form of the resolved equation and subscale coupling (spec: Convective forms)."""

    form: ConvectiveForm = "skew_symmetric_1"


@dataclass(frozen=True)
class ViscosityOptions:
    """Viscous bilinear form and its elementwise operator L_ν (spec: Convective forms)."""

    form: ViscousForm = "laplacian"


@dataclass(frozen=True)
class SolverOptions:
    """Inner linearization, outer fixed point, subscale update and PETSc options."""

    linearization: Linearization = "newton"
    newton_tau: NewtonTau = "frozen"
    picard_iterations_before_newton: int = 0
    inner_rtol: float = 1e-10
    inner_atol: float = 1e-12
    inner_max_iterations: int = 30
    outer_rtol: float = 1e-8
    outer_max_iterations: int = 30
    subscale_update: SubscaleUpdate = "accelerated"
    subscale_rtol: float = 1e-10
    subscale_max_iterations: int = 50
    petsc: Mapping[str, PetscValue] = field(
        default_factory=lambda: MappingProxyType(
            {"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}
        )
    )


@dataclass(frozen=True)
class OutputOptions:
    """Field output frequency, per-step diagnostics and console log level."""

    vtk_every: int = 10
    diagnostics: tuple[str, ...] = ("kinetic_energy", "divergence_l2")
    console_level: LogLevel = "INFO"


@dataclass(frozen=True)
class Config:
    """A fully resolved run configuration."""

    case: str
    run_name: str
    case_parameters: Mapping[str, Any] = field(default_factory=lambda: MappingProxyType({}))
    mesh: MeshOptions = field(default_factory=MeshOptions)
    discretization: DiscretizationOptions = field(default_factory=DiscretizationOptions)
    time: TimeOptions = field(default_factory=TimeOptions)
    stabilization: StabilizationOptions = field(default_factory=StabilizationOptions)
    convection: ConvectionOptions = field(default_factory=ConvectionOptions)
    viscosity: ViscosityOptions = field(default_factory=ViscosityOptions)
    solver: SolverOptions = field(default_factory=SolverOptions)
    output: OutputOptions = field(default_factory=OutputOptions)


_SECTIONS: dict[str, type] = {
    "mesh": MeshOptions,
    "discretization": DiscretizationOptions,
    "time": TimeOptions,
    "stabilization": StabilizationOptions,
    "convection": ConvectionOptions,
    "viscosity": ViscosityOptions,
    "solver": SolverOptions,
    "output": OutputOptions,
}

_NAME = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]*")


def load_config(path: Path) -> Config:
    """Read a TOML file, fill defaults and validate; run_name defaults to the file stem."""
    with path.open("rb") as stream:
        data = tomllib.load(stream)
    config = _config_from_dict(data, default_run_name=path.stem)
    validate_config(config)
    return config


def _config_from_dict(data: Mapping[str, object], default_run_name: str) -> Config:
    """Build a Config from parsed TOML, checking keys and value types."""
    unknown = sorted(set(data) - {"case", "run_name", "case_parameters", *_SECTIONS})
    if unknown:
        raise ValueError(f"unknown top-level keys: {', '.join(unknown)}")
    if "case" not in data:
        raise ValueError("case: required field is missing")
    case_parameters = data.get("case_parameters", {})
    if not isinstance(case_parameters, dict):
        raise ValueError("case_parameters: expected a table")
    sections = {
        name: _parse_section(name, options, data.get(name, {}))
        for name, options in _SECTIONS.items()
    }
    return Config(
        case=_parse_value("case", data["case"], str),
        run_name=_parse_value("run_name", data.get("run_name", default_run_name), str),
        case_parameters=MappingProxyType(dict(case_parameters)),
        **sections,
    )


def validate_config(config: Config) -> None:
    """Raise ValueError for an illegal option or combination (plan §7, rules 1, 2, 4).

    Rule 5 has nothing left to reject: all eight VMS cases are implemented. The Chapter 5 forms
    add their own restrictions (spec: Convective forms and EMAC–VMS).
    """
    for name, value in (("case", config.case), ("run_name", config.run_name)):
        if not _NAME.fullmatch(value):
            raise ValueError(
                f"{name} = {value!r}: use letters, digits, '_', '-' and '.', "
                "starting with a letter, digit or '_'"
            )
    _check_bounds(config)
    steps = config.time.t_end / config.time.dt
    if abs(steps - round(steps)) > 1e-9 * steps:
        raise ValueError(
            f"time.t_end = {config.time.t_end}, time.dt = {config.time.dt}: "
            "t_end must be a whole number of uniform steps"
        )
    _check_mesh(config.mesh)
    _check_stabilization(config)
    _check_forms(config)
    for name in ("subscale_kinetic_energy", "region_subscale_l2"):
        if name in config.output.diagnostics and config.stabilization.method == "none":
            raise ValueError(f"output.diagnostics: {name} requires stabilization.method = 'vms'")
    for name in ("energy_budget", "momentum_budget"):
        if name in config.output.diagnostics and (
            config.time.scheme != "crank_nicolson"
            or config.stabilization.method != "vms"
            or config.stabilization.subscale_dynamics != "dynamic"
        ):
            # The per-step closures are midpoint identities of both scales (spec: Conservation
            # checks).
            raise ValueError(
                f"output.diagnostics: {name} with time.scheme = {config.time.scheme!r}, "
                f"stabilization.method = {config.stabilization.method!r}, "
                f"stabilization.subscale_dynamics = {config.stabilization.subscale_dynamics!r}: "
                "needs 'crank_nicolson' with dynamic VMS subscales"
            )
    solver = config.solver
    if solver.picard_iterations_before_newton > 0 and solver.linearization != "newton":
        raise ValueError(
            f"solver.picard_iterations_before_newton = {solver.picard_iterations_before_newton} "
            f"requires solver.linearization = 'newton', not {solver.linearization!r}"
        )


def check_startup_supported(config: Config, has_exact_solution: bool) -> None:
    """Rule 3 of plan §7: exact-history startup needs a case with an exact solution."""
    if config.time.startup == "exact_history" and not has_exact_solution:
        raise ValueError(
            "time.startup = 'exact_history' requires an exact solution, "
            f"which case {config.case!r} does not provide"
        )


def check_diagnostics_known(config: Config, available: Collection[str]) -> None:
    """Rule 6 of plan §7: every requested diagnostic must exist in the registry."""
    unknown = [name for name in config.output.diagnostics if name not in available]
    if unknown:
        raise ValueError(
            f"output.diagnostics: unknown {', '.join(unknown)}; "
            f"available: {', '.join(sorted(available))}"
        )


def config_as_dict(config: Config) -> dict[str, object]:
    """Every field as plain TOML-ready data; unset optional fields are left out."""
    data: dict[str, object] = {
        "case": config.case,
        "run_name": config.run_name,
        "case_parameters": dict(config.case_parameters),
    }
    for name in _SECTIONS:
        section = getattr(config, name)
        data[name] = {
            item.name: _plain(getattr(section, item.name))
            for item in fields(section)
            if getattr(section, item.name) is not None
        }
    return data


def config_to_toml(config: Config) -> str:
    """The fully resolved configuration as TOML text."""
    return tomli_w.dumps(config_as_dict(config))


def method_name(stabilization: StabilizationOptions) -> str:
    """Galerkin, or the spec name of the VMS case such as Dyn-Nl-OSS."""
    if stabilization.method == "none":
        return "Galerkin"
    return _vms_case_name(
        stabilization.subscale_space, stabilization.subscale_dynamics, stabilization.scale_splitting
    )


def _vms_case_name(space: str, dynamics: str, splitting: str) -> str:
    return "-".join(
        (
            {"static": "Sta", "dynamic": "Dyn"}[dynamics],
            {"linear": "Lin", "nonlinear": "Nl"}[splitting],
            space.upper(),
        )
    )


def _parse_section(name: str, options: type, table: object) -> object:
    if not isinstance(table, dict):
        raise ValueError(f"{name}: expected a table")
    hints = get_type_hints(options)
    unknown = sorted(set(table) - set(hints))
    if unknown:
        raise ValueError(f"unknown keys in [{name}]: {', '.join(unknown)}")
    return options(
        **{key: _parse_value(f"{name}.{key}", value, hints[key]) for key, value in table.items()}
    )


def _parse_value(name: str, value: object, annotation: object) -> object:
    origin = get_origin(annotation)
    if origin is Literal:
        allowed = get_args(annotation)
        if value not in allowed:
            choices = ", ".join(repr(choice) for choice in allowed)
            raise ValueError(f"{name} = {value!r}: expected one of {choices}")
        return value
    if origin is UnionType:
        return _parse_value(name, value, get_args(annotation)[0])
    if origin is tuple:
        if not isinstance(value, list):
            raise ValueError(f"{name} = {value!r}: expected an array")
        item = get_args(annotation)[0]
        return tuple(_parse_value(f"{name}[{i}]", entry, item) for i, entry in enumerate(value))
    if origin is Mapping:
        if not isinstance(value, dict):
            raise ValueError(f"{name} = {value!r}: expected a table")
        scalar = get_args(get_args(annotation)[1])
        invalid = sorted(key for key, entry in value.items() if not isinstance(entry, scalar))
        if invalid:
            raise ValueError(f"{name}: {', '.join(invalid)} must be strings, numbers or booleans")
        return MappingProxyType(dict(value))
    if annotation is float and isinstance(value, int | float) and not isinstance(value, bool):
        return float(value)
    if isinstance(value, annotation) and not (annotation is int and isinstance(value, bool)):
        return value
    raise ValueError(f"{name} = {value!r}: expected {annotation.__name__}")


def _plain(value: object) -> object:
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, Mapping):
        return dict(value)
    return value


def _check_bounds(config: Config) -> None:
    degrees, time, stabilization = config.discretization, config.time, config.stabilization
    solver, output = config.solver, config.output
    limits = (
        ("discretization.velocity_degree", degrees.velocity_degree, ">=", 1),
        ("discretization.pressure_degree", degrees.pressure_degree, ">=", 1),
        ("time.dt", time.dt, ">", 0),
        ("time.t_end", time.t_end, ">", 0),
        ("stabilization.c1", stabilization.c1, ">", 0),
        ("stabilization.c2", stabilization.c2, ">", 0),
        ("solver.picard_iterations_before_newton", solver.picard_iterations_before_newton, ">=", 0),
        ("solver.inner_rtol", solver.inner_rtol, ">", 0),
        ("solver.inner_atol", solver.inner_atol, ">=", 0),
        ("solver.inner_max_iterations", solver.inner_max_iterations, ">=", 1),
        ("solver.outer_rtol", solver.outer_rtol, ">", 0),
        ("solver.outer_max_iterations", solver.outer_max_iterations, ">=", 1),
        ("solver.subscale_rtol", solver.subscale_rtol, ">", 0),
        ("solver.subscale_max_iterations", solver.subscale_max_iterations, ">=", 1),
        ("output.vtk_every", output.vtk_every, ">=", 0),
    )
    for name, value, relation, bound in limits:
        within = value > bound if relation == ">" else value >= bound
        if not (math.isfinite(value) and within):
            raise ValueError(f"{name} = {value}: must be {relation} {bound}")


def _check_mesh(mesh: MeshOptions) -> None:
    if mesh.cells is not None and mesh.file is not None:
        raise ValueError("mesh.cells and mesh.file: cells is for builtin meshes, file for gmsh")
    if mesh.source == "gmsh" and mesh.cells is not None:
        raise ValueError("mesh.cells: only used when mesh.source = 'builtin'")
    if mesh.source == "builtin" and mesh.file is not None:
        raise ValueError("mesh.file: only used when mesh.source = 'gmsh'")
    if mesh.cells is not None and (len(mesh.cells) not in (2, 3) or min(mesh.cells) < 1):
        raise ValueError(f"mesh.cells = {list(mesh.cells)}: needs 2 or 3 entries, each >= 1")


def _check_stabilization(config: Config) -> None:
    degrees, stabilization = config.discretization, config.stabilization
    velocity, pressure = degrees.velocity_degree, degrees.pressure_degree
    if stabilization.method == "none":
        if velocity != pressure + 1:
            raise ValueError(
                f"discretization.velocity_degree = {velocity}, "
                f"discretization.pressure_degree = {pressure}: stabilization.method = 'none' "
                "needs a Taylor–Hood pair, velocity_degree = pressure_degree + 1"
            )
        changed = [
            f"stabilization.{item.name}"
            for item in fields(stabilization)
            if item.name != "method"
            and getattr(stabilization, item.name) != getattr(StabilizationOptions(), item.name)
        ]
        if changed:
            raise ValueError(
                f"{', '.join(changed)}: not allowed with stabilization.method = 'none'"
            )
        return
    if velocity != pressure:
        raise ValueError(
            f"discretization.velocity_degree = {velocity}, "
            f"discretization.pressure_degree = {pressure}: stabilization.method = 'vms' "
            "needs equal orders, velocity_degree = pressure_degree"
        )


def _check_forms(config: Config) -> None:
    form, solver = config.convection.form, config.solver
    if form == "emac_vms":
        for name, value, required in (
            ("stabilization.method", config.stabilization.method, "vms"),
            ("viscosity.form", config.viscosity.form, "stress"),
        ):
            if value != required:
                raise ValueError(
                    f"convection.form = 'emac_vms', {name} = {value!r}: "
                    f"EMAC–VMS needs {name} = {required!r}"
                )
    if form in ("emac", "emac_vms"):
        # Owner decision: no Picard linearization of the kinetic pressure or substituted stress.
        if solver.linearization != "newton":
            raise ValueError(
                f"convection.form = {form!r}, solver.linearization = {solver.linearization!r}: "
                "EMAC forms need solver.linearization = 'newton'"
            )
        if solver.picard_iterations_before_newton > 0:
            raise ValueError(
                f"convection.form = {form!r}, solver.picard_iterations_before_newton = "
                f"{solver.picard_iterations_before_newton}: EMAC forms need no Picard warm-up"
            )
