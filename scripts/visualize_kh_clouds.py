"""Gather saved P1 KH fields into cloud-proxy VTK, figures, an animation and diagnostic CSV."""

import argparse
import csv
import json
import logging
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from PIL import Image
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
from vtkmodules.vtkIOXML import (
    vtkXMLPUnstructuredGridReader,
    vtkXMLUnstructuredGridReader,
    vtkXMLUnstructuredGridWriter,
)

from navier_stokes_vms.postprocessing.cloud_visualization import cloud_fields, cloud_statistics

LOGGER = logging.getLogger(__name__)


def read_series(path: Path) -> dict[float, Path]:
    """Read the actual output times and piece filenames, including MPI PVTU collections."""
    series = {}
    for entry in ET.parse(path).findall(".//DataSet"):
        time = float(entry.attrib["timestep"])
        if not np.isfinite(time) or time in series:
            raise ValueError(f"invalid or duplicate field time in {path}")
        series[time] = path.parent / entry.attrib["file"]
    if not series or list(series) != sorted(series):
        raise ValueError(f"empty or unordered field series: {path}")
    return series


def read_frame(paths: dict[str, Path], domain_area: float) -> tuple:
    """Read the whole gathered mesh and reject mismatched, non-P1 or duplicated cells."""
    fields = {}
    grid = None
    for name, path in paths.items():
        reader = (
            vtkXMLPUnstructuredGridReader()
            if path.suffix == ".pvtu"
            else vtkXMLUnstructuredGridReader()
        )
        reader.SetFileName(str(path))
        reader.Update()
        current = reader.GetOutput()
        points = vtk_to_numpy(current.GetPoints().GetData())
        cells = current.GetCells()
        offsets = vtk_to_numpy(cells.GetOffsetsArray())
        if not np.all(np.diff(offsets) == 3):
            raise ValueError("cloud visualization currently supports P1 triangles only")
        triangles = vtk_to_numpy(cells.GetConnectivityArray()).reshape(-1, 3)
        if grid is None:
            grid = vtkUnstructuredGrid()
            grid.DeepCopy(current)
            coordinates, topology = points.copy(), triangles.copy()
        elif not np.array_equal(coordinates, points) or not np.array_equal(topology, triangles):
            raise ValueError("saved velocity, temperature and moisture meshes differ")
        array = current.GetPointData().GetArray(name)
        if array is None:
            raise ValueError(f"missing point field {name} in {path}")
        fields[name] = vtk_to_numpy(array).copy()
        if not np.isfinite(fields[name]).all():
            raise ValueError(f"nonfinite field {name} in {path}")
    edges = coordinates[topology[:, 1:]] - coordinates[topology[:, :1]]
    area = np.sum(abs(edges[:, 0, 0] * edges[:, 1, 1] - edges[:, 1, 0] * edges[:, 0, 1])) / 2
    if not np.isclose(area, domain_area, rtol=1e-10, atol=1e-10):
        raise ValueError("gathered VTK cell areas disagree with the full domain area")
    return grid, coordinates, topology, fields


def write_cloud_vtu(
    path: Path,
    grid: vtkUnstructuredGrid,
    fields: dict[str, np.ndarray],
    cloud: np.ndarray,
    vorticity: np.ndarray,
) -> None:
    """Write combined unmodified transport fields, nodal q_c and cell vorticity."""
    for name, values in {**fields, "q_c": cloud}.items():
        array = numpy_to_vtk(np.ascontiguousarray(values), deep=True)
        array.SetName(name)
        grid.GetPointData().AddArray(array)
    array = numpy_to_vtk(np.ascontiguousarray(vorticity), deep=True)
    array.SetName("vorticity")
    grid.GetCellData().AddArray(array)
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    if writer.Write() != 1:
        raise OSError(f"could not write {path}")


def plot_frame(
    path: Path,
    time: float,
    points: np.ndarray,
    triangles: np.ndarray,
    theta: np.ndarray,
    cloud: np.ndarray,
    vorticity: np.ndarray,
    cloud_limit: float,
) -> None:
    """Plot the central cross-section with common color scales across animation frames."""
    triangulation = Triangulation(points[:, 0], points[:, 1], triangles)
    fig, axes = plt.subplots(2, 1, figsize=(12, 5.8), sharex=True, sharey=True)
    image = axes[0].tripcolor(
        triangulation, cloud, shading="gouraud", cmap="Blues", vmin=0, vmax=cloud_limit
    )
    axes[0].tricontour(
        triangulation, theta, levels=[-0.75, -0.5, 0, 0.5, 0.75], colors="#4b4b4b", linewidths=0.6
    )
    fig.colorbar(image, ax=axes[0], label=r"Cloud proxy $q_c$", pad=0.02)
    axes[0].set_title(r"Cloud proxy over temperature contours ($\theta=-0.75,-0.5,0,0.5,0.75$)")
    image = axes[1].tripcolor(
        triangulation, facecolors=vorticity, cmap="RdBu_r", vmin=-1.5, vmax=1.5
    )
    fig.colorbar(image, ax=axes[1], label="Vorticity", pad=0.02, extend="both")
    axes[1].set_title("Velocity roll-up")
    for ax in axes:
        ax.set(xlim=(np.min(points[:, 0]), np.max(points[:, 0])), ylim=(-5, 5), ylabel="z")
        ax.set_aspect("equal")
    axes[1].set_xlabel("x")
    fig.suptitle(f"Stratified Kelvin–Helmholtz · four-wave visualization pilot · t = {time:g}")
    fig.text(
        0.5,
        0.015,
        "Passive moisture / diagnostic saturation proxy · no latent heat · "
        "central |z| ≤ 5 of the height-20 domain",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    fig.savefig(path, dpi=140, facecolor="white")
    plt.close(fig)


def make_animation(folder: Path, files: list[Path]) -> None:
    """Produce a browser-viewable GIF and, when ffmpeg is present, an MP4."""
    images = []
    for path in files:
        with Image.open(path) as image:
            images.append(image.convert("RGB").quantize(colors=128))
    images[0].save(
        folder / "clouds.gif",
        save_all=True,
        append_images=images[1:],
        duration=100,
        loop=0,
        optimize=False,
    )
    for image in images:
        image.close()
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        subprocess.run(
            [
                ffmpeg,
                "-y",
                "-loglevel",
                "error",
                "-framerate",
                "10",
                "-i",
                str(folder / "frames/frame_%04d.png"),
                "-vf",
                "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                "-c:v",
                "libx264",
                "-pix_fmt",
                "yuv420p",
                str(folder / "clouds.mp4"),
            ],
            check=True,
            capture_output=True,
            text=True,
        )


def plot_history(folder: Path, rows: list[dict[str, float]], initial_bound: float | None) -> None:
    """Plot cloud coverage, vertical motion and moisture extrema without clipping q."""
    time = [r["time"] for r in rows]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    axes[0].plot(time, [100 * r["cloud_fraction"] for r in rows])
    axes[0].set_ylabel("Cloud fraction of whole domain (%)")
    axes[1].semilogy(time, [r["vertical_kinetic_energy"] for r in rows])
    axes[1].set_ylabel("Vertical kinetic energy")
    for name in ("moisture_nodal_min", "moisture_nodal_max"):
        axes[2].plot(time, [r[name] for r in rows], label=name.removeprefix("moisture_nodal_"))
    axes[2].set_ylabel("Unmodified moisture extrema")
    if initial_bound is not None:
        axes[2].axhline(initial_bound, color="0.4", linestyle="--", label="Initial upper bound")
    axes[2].legend()
    for ax in axes:
        ax.set_xlabel("Time")
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(folder / "cloud_history.png", dpi=180)
    fig.savefig(folder / "cloud_history.pdf")
    plt.close(fig)


def compare_solver_diagnostics(folder: Path, rows: list[dict[str, float]]) -> dict:
    """Check saved-field integrals against assembled diagnostics at the same endpoint times."""
    with (folder / "diagnostics.csv").open() as stream:
        records = list(csv.DictReader(stream))
    by_time = {float(record["time"]): record for record in records}
    errors = {}
    for plotted, assembled in (
        ("kinetic_energy", "kinetic_energy"),
        ("vertical_kinetic_energy", "vertical_kinetic_energy"),
        ("temperature_content", "thermal_content"),
        ("moisture_content", "passive_scalar_content_q"),
    ):
        if assembled not in records[0]:
            continue
        maximum = 0.0
        for row in rows:
            reference = float(by_time[row["time"]][assembled])
            error = abs(reference - row[plotted])
            if error > 1e-10 * max(1.0, abs(reference)):
                raise ValueError(f"saved-field {plotted} disagrees with assembled FE diagnostic")
            maximum = max(maximum, error)
        errors[plotted] = maximum
    result = {"saved_integrals_max_absolute_difference": errors, "solver_steps": len(records)}
    for name in ("thermal_content_budget_q_residual", "boussinesq_energy_budget_residual"):
        if name in records[0]:
            result[f"maximum_absolute_{name}"] = max(abs(float(r[name])) for r in records)
    if "passive_scalar_content_q" in records[0]:
        initial = float(records[0]["passive_scalar_content_q"])
        result["all_steps_moisture_relative_drift"] = max(
            abs(float(r["passive_scalar_content_q"]) / initial - 1) for r in records
        )
    return result


def visualize(folder: Path, output: Path, allow_incomplete: bool = False) -> dict:
    """Post-process one saved trajectory; permit an explicitly labelled in-progress preview."""
    info = json.loads((folder / "run_info.json").read_text())
    if info["exit_status"] != 0 and not allow_incomplete:
        raise ValueError("run is not complete; --allow-incomplete is only for a labelled preview")
    metadata = json.loads((folder / "thermal_metadata.json").read_text())
    case = metadata["case_metadata"]
    parameters = case["parameters"]
    if not parameters["moisture"]:
        raise ValueError("cloud visualization requires a transported q field")
    if metadata["temperature_degree"] != 1:
        raise ValueError("offline cloud diagnostics currently require P1 output")
    series = {
        name: read_series(folder / f"fields/{name}.pvd")
        for name in ("velocity", "temperature", "q")
    }
    times = list(series["velocity"])
    if any(list(values) != times for values in series.values()):
        raise ValueError("temperature, moisture and velocity output times differ")
    output.mkdir(parents=True, exist_ok=True)
    (output / "frames").mkdir(exist_ok=True)
    (output / "fields").mkdir(exist_ok=True)
    alpha, gamma = parameters["saturation_temperature_rate"], parameters["saturation_height_rate"]
    area = case["domain_length"] * parameters["height"]
    rows = []
    collection = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
    datasets = ET.SubElement(collection, "Collection")
    for index, time in enumerate(times):
        grid, points, triangles, fields = read_frame(
            {name: values[time] for name, values in series.items()}, area
        )
        arguments = (
            points,
            triangles,
            fields["temperature"],
            fields["q"],
            fields["velocity"],
            alpha,
            gamma,
        )
        cloud, vorticity = cloud_fields(*arguments)
        values = cloud_statistics(*arguments, parameters["richardson"], quadrature_order=12)
        coarse = cloud_statistics(*arguments, parameters["richardson"], quadrature_order=8)
        rows.append(
            {
                "time": time,
                **values,
                "cloud_fraction_quadrature_change": abs(
                    values["cloud_fraction"] - coarse["cloud_fraction"]
                ),
            }
        )
        filename = f"fields/cloud_{index:04d}.vtu"
        write_cloud_vtu(output / filename, grid, fields, cloud, vorticity)
        ET.SubElement(datasets, "DataSet", timestep=str(time), file=filename)
    ET.ElementTree(collection).write(output / "clouds.pvd", encoding="utf-8", xml_declaration=True)
    with (output / "cloud_diagnostics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    limit = max(max(r["cloud_proxy_nodal_max"] for r in rows), 1e-5)
    images = []
    for index, time in enumerate(times):
        _, points, triangles, fields = read_frame(
            {name: values[time] for name, values in series.items()}, area
        )
        cloud, vorticity = cloud_fields(
            points, triangles, fields["temperature"], fields["q"], fields["velocity"], alpha, gamma
        )
        path = output / f"frames/frame_{index:04d}.png"
        plot_frame(path, time, points, triangles, fields["temperature"], cloud, vorticity, limit)
        images.append(path)
    make_animation(output, images)
    initial_bound = (
        parameters["rh_peak"] * np.exp(gamma**2 * parameters["moisture_width"] ** 2 / 4)
        if alpha == 0
        else None
    )
    plot_history(output, rows, initial_bound)
    maximum = max(rows, key=lambda row: row["cloud_proxy_nodal_max"])
    shutil.copyfile(images[times.index(maximum["time"])], output / "maximum_cloud.png")
    maximum_point = None
    if maximum["cloud_proxy_nodal_max"] > 0:
        _, points, triangles, fields = read_frame(
            {name: values[maximum["time"]] for name, values in series.items()}, area
        )
        cloud, _ = cloud_fields(
            points, triangles, fields["temperature"], fields["q"], fields["velocity"], alpha, gamma
        )
        index = int(np.argmax(cloud))
        maximum_point = {
            "position": points[index, :2].tolist(),
            "moisture": float(fields["q"][index]),
            "saturation": float(
                np.exp(alpha * fields["temperature"][index] - gamma * points[index, 1])
            ),
        }
    summary = {
        "run": str(folder),
        "run_info": info,
        "parameters": parameters,
        "status": "complete" if info["exit_status"] == 0 else "in_progress_preview",
        "frames": len(times),
        "first_time": times[0],
        "last_time": times[-1],
        "first_cloud_time": next((r["time"] for r in rows if r["cloud_fraction"] > 0), None),
        "maximum_cloud": maximum,
        "maximum_cloud_point": maximum_point,
        "maximum_cloud_fraction": max(r["cloud_fraction"] for r in rows),
        "maximum_fraction_quadrature_change": max(
            r["cloud_fraction_quadrature_change"] for r in rows
        ),
        "moisture_content_relative_drift": max(
            abs(r["moisture_content"] / rows[0]["moisture_content"] - 1) for r in rows
        ),
        "moisture_nodal_min": min(r["moisture_nodal_min"] for r in rows),
        "moisture_nodal_max": max(r["moisture_nodal_max"] for r in rows),
        "initial_continuum_moisture_upper_bound": initial_bound,
        "maximum_moisture_bound_excess": (
            max(0.0, max(r["moisture_nodal_max"] for r in rows) - initial_bound)
            if initial_bound is not None
            else None
        ),
        "cloud_threshold": 1e-6,
        "solver_comparison": compare_solver_diagnostics(folder, rows),
        "quadrature_orders_per_direction": [8, 12],
        "display": "nodal positive part; global common cloud scale; central |z|<=5",
        "fraction": "raw P1 q-exp(alpha*T-gamma*z) at triangle quadrature points; whole domain",
        "limits": (
            "Coarse visualization pilot; no latent heat, condensate removal "
            "or production convergence claim."
        ),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    output = args.output or args.run / "cloud_visualization"
    summary = visualize(args.run, output, args.allow_incomplete)
    LOGGER.info(
        "Wrote %d frames through t=%g to %s", summary["frames"], summary["last_time"], output
    )


if __name__ == "__main__":
    main()
