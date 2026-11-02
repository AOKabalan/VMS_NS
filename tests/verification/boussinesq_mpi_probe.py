"""Serial/MPI comparison probe for global thermal diagnostics and OSS orthogonality."""

import json
import logging
import sys
from pathlib import Path

import firedrake as fd
import ufl
from mpi4py import MPI

from navier_stokes_vms.cases import build_problem
from navier_stokes_vms.config import Config, MeshOptions, OutputOptions, TimeOptions
from navier_stokes_vms.postprocessing.diagnostics import evaluate_diagnostics
from navier_stokes_vms.solver.time_stepper import TimeStepper


def probe(folder: Path) -> None:
    """Write globally reduced solve diagnostics and orthogonality norms on rank zero."""
    comm = MPI.COMM_WORLD
    if comm.rank == 0:
        folder.mkdir(parents=True, exist_ok=True)
    comm.barrier()
    config = Config(
        case="differentially_heated_cavity",
        run_name="mpi_probe",
        mesh=MeshOptions(cells=(4, 4)),
        time=TimeOptions(dt=0.05, t_end=0.1),
        output=OutputOptions(
            vtk_every=0,
            diagnostics=(
                "temperature_l2",
                "thermal_content",
                "thermal_variance_budget",
                "thermal_subscale_l2",
                "cavity_nusselt",
            ),
        ),
    )
    stepper = TimeStepper(config, build_problem(config), folder, logging.getLogger("mpi_probe"))
    stepper.run()
    values = evaluate_diagnostics(config.output.diagnostics, stepper._diagnostic_state)
    for name, field in (("velocity", stepper._velocity), ("temperature", stepper._temperature)):
        test = fd.TestFunction(field.function_space().collapse())
        sub = stepper._u_sub if name == "velocity" else stepper._scalar_subscales[name]
        residual = fd.assemble(ufl.inner(sub, test) * stepper._dx)
        with residual.dat.vec_ro as vector:
            values[f"{name}_orthogonality"] = vector.norm()
    if comm.rank == 0:
        (folder / "global_diagnostics.json").write_text(json.dumps(values, indent=2))


if __name__ == "__main__":
    probe(Path(sys.argv[1]))
