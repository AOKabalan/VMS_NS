"""B25 cloud proxy and whole-mesh diagnostics of affine P1/P2 triangle VTK fields."""

import numpy as np
from numpy.polynomial.legendre import leggauss


def triangle_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Barycentric Duffy–Gauss points and weights normalized to unit triangle area."""
    nodes, weights = leggauss(order)
    r, s = np.meshgrid((nodes + 1) / 2, (nodes + 1) / 2, indexing="ij")
    barycentric = np.array([1 - r.ravel(), (r * (1 - s)).ravel(), (r * s).ravel()])
    return barycentric, (np.outer(weights, weights) * r / 2).ravel()


def triangle_basis(barycentric: np.ndarray, nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Affine triangle Lagrange basis and reference gradients in VTK vertex/edge order."""
    reference = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
    if nodes == 3:
        return barycentric, np.broadcast_to(reference[:, None, :], (3, barycentric.shape[1], 2))
    if nodes != 6:
        raise ValueError("cloud diagnostics require P1/P2 triangles with 3 or 6 nodes")
    basis = [barycentric[i] * (2 * barycentric[i] - 1) for i in range(3)]
    gradient = [(4 * barycentric[i, :, None] - 1) * reference[i] for i in range(3)]
    for i, j in ((0, 1), (1, 2), (2, 0)):
        basis.append(4 * barycentric[i] * barycentric[j])
        gradient.append(
            4 * (barycentric[i, :, None] * reference[j] + barycentric[j, :, None] * reference[i])
        )
    return np.array(basis), np.array(gradient)


def triangle_gradients(
    points: np.ndarray, triangles: np.ndarray, values: np.ndarray, barycentric: np.ndarray
) -> np.ndarray:
    """Physical gradients [cell, point, field, direction] of one or more nodal fields."""
    _, derivative = triangle_basis(barycentric, triangles.shape[1])
    reference = np.einsum("niv,iqr->nqvr", values[triangles], derivative)
    edges = points[triangles[:, 1:3], :2] - points[triangles[:, :1], :2]
    determinant = edges[:, 0, 0] * edges[:, 1, 1] - edges[:, 1, 0] * edges[:, 0, 1]
    scale = determinant[:, None, None]
    dx = (
        reference[..., 0] * edges[:, 1, 1, None, None]
        - reference[..., 1] * edges[:, 0, 1, None, None]
    ) / scale
    dz = (
        reference[..., 1] * edges[:, 0, 0, None, None]
        - reference[..., 0] * edges[:, 1, 0, None, None]
    ) / scale
    return np.stack((dx, dz), axis=-1)


def cloud_fields(
    points: np.ndarray,
    triangles: np.ndarray,
    temperature: np.ndarray,
    moisture: np.ndarray,
    velocity: np.ndarray,
    alpha: float,
    gamma: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Nodal positive part q_c and exact cell-mean vorticity, without modifying q."""
    saturation = np.exp(alpha * temperature - gamma * points[:, 1])
    cloud = np.maximum(moisture - saturation, 0)
    gradient = triangle_gradients(points, triangles, velocity[:, :2], np.full((3, 1), 1 / 3))
    return cloud, gradient[:, 0, 1, 0] - gradient[:, 0, 0, 1]


def cloud_statistics(
    points: np.ndarray,
    triangles: np.ndarray,
    temperature: np.ndarray,
    moisture: np.ndarray,
    velocity: np.ndarray,
    alpha: float,
    gamma: float,
    buoyancy_coefficient: float,
    quadrature_order: int = 8,
    threshold: float = 1e-6,
) -> dict[str, float]:
    """Integrate original P1/P2 FE fields and raw B25 on every gathered affine triangle."""
    coordinates = points[triangles[:, :3], :2]
    edges = coordinates[:, 1:] - coordinates[:, :1]
    areas = abs(edges[:, 0, 0] * edges[:, 1, 1] - edges[:, 1, 0] * edges[:, 0, 1]) / 2
    if not np.isfinite(areas).all() or np.any(areas <= 0):
        raise ValueError("cloud visualization requires nondegenerate affine triangles")
    barycentric, weights = triangle_rule(quadrature_order)
    basis, _ = triangle_basis(barycentric, triangles.shape[1])
    theta = temperature[triangles] @ basis
    tracers = moisture[triangles] @ basis
    heights = coordinates[:, :, 1] @ barycentric
    velocities = np.einsum("niv,iq->nqv", velocity[triangles, :2], basis)
    raw_cloud = np.maximum(tracers - np.exp(alpha * theta - gamma * heights), 0)
    cloud, _ = cloud_fields(points, triangles, temperature, moisture, velocity, alpha, gamma)
    volume = np.sum(areas)

    def integral(values: np.ndarray) -> float:
        return float(np.dot(areas, values @ weights))

    return {
        "kinetic_energy": integral(np.sum(velocities**2, axis=-1) / 2),
        "vertical_kinetic_energy": integral(velocities[..., 1] ** 2 / 2),
        "potential_energy": integral(-buoyancy_coefficient * theta * heights),
        "temperature_content": integral(theta),
        "moisture_content": integral(tracers),
        "moisture_nodal_min": float(np.min(moisture)),
        "moisture_nodal_max": float(np.max(moisture)),
        "moisture_quadrature_min": float(np.min(tracers)),
        "moisture_quadrature_max": float(np.max(tracers)),
        "cloud_proxy_nodal_max": float(np.max(cloud)),
        "cloud_proxy_quadrature_max": float(np.max(raw_cloud)),
        "cloud_fraction": integral(raw_cloud > threshold) / volume,
    }
