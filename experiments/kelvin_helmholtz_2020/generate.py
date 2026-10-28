"""Generate the matched SKEW1 and EMAC–VMS Kelvin–Helmholtz configurations."""

import argparse
from pathlib import Path

from navier_stokes_vms.config import (
    Config,
    ConvectionOptions,
    DiscretizationOptions,
    MeshOptions,
    OutputOptions,
    SolverOptions,
    StabilizationOptions,
    TimeOptions,
    ViscosityOptions,
    config_to_toml,
    validate_config,
)

HERE = Path(__file__).resolve().parent


def main() -> None:
    """Write two configurations; use another output directory for parameter sweeps."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=196)
    parser.add_argument("--degree", type=int, choices=(1, 2), default=1)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--end", type=float, default=14.0)
    parser.add_argument("--tau-scale", type=float, default=1.0)
    parser.add_argument(
        "--linear-solver",
        choices=("direct", "reuse_lu", "reuse_lu_tuned", "block_oss"),
        default="direct",
    )
    parser.add_argument("--output", type=Path, default=HERE / "configs")
    args = parser.parse_args()
    if args.tau_scale <= 0:
        parser.error("--tau-scale must be positive")
    for short, form in (("skew1", "skew_symmetric_1"), ("emac_vms", "emac_vms")):
        name = f"kh2020_{short}_p{args.degree}_n{args.cells}"
        if args.tau_scale != 1:
            name += "_tau" + f"{args.tau_scale:g}".replace(".", "p")
        petsc = dict(SolverOptions().petsc)
        if args.linear_solver != "direct":
            name += "_" + args.linear_solver
            # Reassemble the current Jacobian; reuse only its LU preconditioner.
            # Persist the rebuild counter across the outer subscale solves.
            petsc.update(
                ksp_type="gmres",
                ksp_pc_side="right",
                ksp_norm_type="unpreconditioned",
                ksp_rtol=1e-10,
                ksp_atol=1e-14,
                ksp_max_it=100,
                ksp_gmres_restart=50,
                snes_lag_jacobian=1,
                snes_lag_preconditioner=5,
                snes_lag_preconditioner_persists=True,
            )
            if args.linear_solver in ("reuse_lu_tuned", "block_oss"):
                petsc.update(
                    snes_lag_preconditioner=20,
                    snes_ksp_ew=True,
                    snes_ksp_ew_version=2,
                    snes_ksp_ew_rtol0=0.1,
                    snes_ksp_ew_rtolmax=0.5,
                )
            if args.linear_solver == "block_oss":
                # Keep the full consistent-mass Jacobian as the operator.
                # Approximate only its inverse: mass blocks first, then (u,p).
                petsc.pop("pc_factor_mat_solver_type")
                petsc.update(
                    mat_type="aij",
                    ksp_type="fgmres",
                    ksp_max_it=200,
                    pc_type="fieldsplit",
                    pc_fieldsplit_type="schur",
                    pc_fieldsplit_schur_fact_type="full",
                    pc_fieldsplit_schur_precondition="a11",
                    pc_fieldsplit_0_fields="2,3",
                    pc_fieldsplit_1_fields="0,1",
                    fieldsplit_0_ksp_type="cg",
                    fieldsplit_0_ksp_rtol=1e-2,
                    fieldsplit_0_ksp_atol=1e-14,
                    fieldsplit_0_ksp_max_it=50,
                    fieldsplit_0_pc_type="jacobi",
                    fieldsplit_1_ksp_type="preonly",
                    fieldsplit_1_pc_type="lu",
                    fieldsplit_1_pc_factor_mat_solver_type="mumps",
                )
        cfg = Config(
            case="kelvin_helmholtz",
            run_name=name,
            case_parameters={
                "reynolds": 1000.0,
                "delta_0": 1 / 28,
                "u_infinity": 1.0,
                "noise": 0.001,
            },
            mesh=MeshOptions(cells=(args.cells, args.cells)),
            discretization=DiscretizationOptions(args.degree, args.degree),
            time=TimeOptions(scheme="bdf2", dt=args.dt, t_end=args.end, startup="lower_order"),
            stabilization=StabilizationOptions(c1=4 / args.tau_scale, c2=2 / args.tau_scale),
            convection=ConvectionOptions(form=form),
            viscosity=ViscosityOptions(form="stress"),
            solver=SolverOptions(
                inner_max_iterations=200,
                outer_rtol=1e-9,
                outer_max_iterations=100,
                subscale_rtol=1e-11,
                subscale_max_iterations=500,
                petsc=petsc,
            ),
            output=OutputOptions(
                vtk_every=0,
                console_level="WARNING",
                diagnostics=(
                    "kinetic_energy",
                    "subscale_kinetic_energy",
                    "linear_momentum",
                    "angular_momentum",
                    "divergence_l2",
                    "convective_energy_defect",
                ),
            ),
        )
        validate_config(cfg)
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / f"{name}.toml").write_text(config_to_toml(cfg))


if __name__ == "__main__":
    main()
