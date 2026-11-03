"""Cloud diagnostics distinguish raw saturation from interpolated display fields."""

import numpy as np

from navier_stokes_vms.postprocessing.cloud_visualization import cloud_fields, cloud_statistics


def test_raw_cloud_can_exist_between_cloud_free_vertices() -> None:
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    triangles = np.array([[0, 1, 2], [1, 2, 3]])
    theta = np.zeros(4)
    q = np.exp(-points[:, 1])
    original = q.copy()
    velocity = np.column_stack((-points[:, 1], points[:, 0]))
    cloud, vorticity = cloud_fields(points, triangles, theta, q, velocity, 0.0, 1.0)
    assert np.array_equal(cloud, np.zeros(4))
    assert np.allclose(vorticity, 2.0)
    values = cloud_statistics(points, triangles, theta, q, velocity, 0.0, 1.0, 0.1)
    assert values["cloud_fraction"] > 0.999
    assert values["cloud_proxy_quadrature_max"] > 0.07
    assert np.isclose(values["moisture_content"], (1 + np.exp(-1)) / 2)
    assert np.isclose(values["kinetic_energy"], 1 / 3)
    assert np.isclose(values["vertical_kinetic_energy"], 1 / 6)
    assert np.array_equal(q, original)


def test_cloud_fraction_against_analytic_height_threshold() -> None:
    rows = np.linspace(0, 1, 21)
    points = np.array([(x, z) for z in rows for x in (0.0, 1.0)])
    triangles = np.array(
        [
            tri
            for j in range(20)
            for tri in ((2 * j, 2 * j + 1, 2 * j + 2), (2 * j + 1, 2 * j + 3, 2 * j + 2))
        ]
    )
    threshold = 1e-6
    values = cloud_statistics(
        points,
        triangles,
        points[:, 1],
        np.full(42, 0.9),
        np.zeros((42, 2)),
        0.0,
        0.2,
        0.1,
        12,
        threshold,
    )
    exact_fraction = 1 + np.log(0.9 - threshold) / 0.2
    assert abs(values["cloud_fraction"] - exact_fraction) < 0.002
    assert np.isclose(values["potential_energy"], -0.1 / 3)


def test_quadratic_vtk_fields_keep_original_derivatives_and_integrals() -> None:
    points = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5],
            [0.5, 1.0],
            [1.0, 0.5],
        ]
    )
    triangles = np.array([[0, 1, 2, 4, 5, 6], [1, 2, 3, 5, 7, 8]])
    theta = np.sum(points**2, axis=1)
    velocity = np.column_stack((-(points[:, 1] ** 2), points[:, 0] ** 2))
    q = np.ones(len(points))
    _, vorticity = cloud_fields(points, triangles, theta, q, velocity, 0.0, 0.2)
    assert np.allclose(vorticity, [4 / 3, 8 / 3])
    values = cloud_statistics(points, triangles, theta, q, velocity, 0.0, 0.2, 0.1)
    assert np.isclose(values["temperature_content"], 2 / 3)
    assert np.isclose(values["kinetic_energy"], 1 / 5)
    assert np.isclose(values["vertical_kinetic_energy"], 1 / 10)
    assert np.isclose(values["potential_energy"], -0.1 * 5 / 12)
