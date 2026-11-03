"""B25 cloud proxy and whole-mesh P1 triangle diagnostics from gathered VTK fields."""

import numpy as np
from numpy.polynomial.legendre import leggauss


def cloud_fields(
    points: np.ndarray,
    triangles: np.ndarray,
    temperature: np.ndarray,
    moisture: np.ndarray,
    velocity: np.ndarray,
    alpha: float,
    gamma: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Nodal positive part q_c and cell vorticity ∂x w − ∂z u, without modifying q."""
    saturation = np.exp(alpha * temperature - gamma * points[:, 1])
    cloud = np.maximum(moisture - saturation, 0)
    edges = points[triangles[:, 1:]] - points[triangles[:, :1]]
    differences = velocity[triangles[:, 1:]] - velocity[triangles[:, :1]]
    determinant = edges[:, 0, 0] * edges[:, 1, 1] - edges[:, 1, 0] * edges[:, 0, 1]
    vertical_x = (
        differences[:, 0, 1] * edges[:, 1, 1] - differences[:, 1, 1] * edges[:, 0, 1]
    ) / determinant
    horizontal_z = (
        edges[:, 0, 0] * differences[:, 1, 0] - edges[:, 1, 0] * differences[:, 0, 0]
    ) / determinant
    return cloud, vertical_x - horizontal_z


def cloud_statistics(
    points: np.ndarray,
    triangles: np.ndarray,
    temperature: np.ndarray,
    moisture: np.ndarray,
    velocity: np.ndarray,
    alpha: float,
    gamma: float,
    richardson: float,
    quadrature_order: int = 8,
    threshold: float = 1e-6,
) -> dict[str, float]:
    """Integrate raw B25 at Duffy–Gauss points on every owned triangle of the gathered mesh."""
    coordinates = points[triangles, :2]
    edges = coordinates[:, 1:] - coordinates[:, :1]
    areas = abs(edges[:, 0, 0] * edges[:, 1, 1] - edges[:, 1, 0] * edges[:, 0, 1]) / 2
    if not np.isfinite(areas).all() or np.any(areas <= 0):
        raise ValueError("cloud visualization requires nondegenerate planar triangles")
    volume = np.sum(areas)
    nodes, weights = leggauss(quadrature_order)
    r, s = np.meshgrid((nodes + 1) / 2, (nodes + 1) / 2, indexing="ij")
    barycentric = np.array([1 - r.ravel(), (r * (1 - s)).ravel(), (r * s).ravel()])
    normalized_weights = (np.outer(weights, weights) * r / 2).ravel()
    temperatures = temperature[triangles] @ barycentric
    tracers = moisture[triangles] @ barycentric
    heights = coordinates[:, :, 1] @ barycentric
    raw_cloud = np.maximum(tracers - np.exp(alpha * temperatures - gamma * heights), 0)
    cloud, _ = cloud_fields(points, triangles, temperature, moisture, velocity, alpha, gamma)
    nodal_velocity = velocity[triangles, :2]
    square_integrals = (
        areas[:, None]
        / 12
        * (np.sum(nodal_velocity, axis=1) ** 2 + np.sum(nodal_velocity**2, axis=1))
    )
    theta = temperature[triangles]
    z = coordinates[:, :, 1]
    potential = (
        -richardson
        * areas
        / 12
        * (np.sum(theta, axis=1) * np.sum(z, axis=1) + np.sum(theta * z, axis=1))
    )
    return {
        "kinetic_energy": float(np.sum(square_integrals) / 2),
        "vertical_kinetic_energy": float(np.sum(square_integrals[:, 1]) / 2),
        "potential_energy": float(np.sum(potential)),
        "temperature_content": float(np.dot(areas, np.mean(theta, axis=1))),
        "moisture_content": float(np.dot(areas, np.mean(moisture[triangles], axis=1))),
        "moisture_nodal_min": float(np.min(moisture)),
        "moisture_nodal_max": float(np.max(moisture)),
        "cloud_proxy_nodal_max": float(np.max(cloud)),
        "cloud_proxy_quadrature_max": float(np.max(raw_cloud)),
        "cloud_fraction": float(
            np.dot(areas, (raw_cloud > threshold) @ normalized_weights) / volume
        ),
    }
