"""Tools for surface integrals and hydrostatics."""
# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

from abc import ABC
import numpy as np

# NEWER IMPLEMENTATIONS FROM CAPYTAINE

class SurfaceIntegralsMixin(ABC):
    # These methods need to be defined for both Mesh and CollectionOfMeshes with the exact same definitions.
    # To avoid redunduncy, they are defined here in a mixin inherited by both Mesh and CollectionOfMeshes.

    def surface_integral(self, data, **kwargs):
        """Returns integral of given data along wet surface area."""
        return np.sum(data * self.faces_areas, **kwargs)

    def waterplane_integral(self, data, **kwargs):
        """Returns integral of given data along water plane area."""
        return self.surface_integral(self.faces_normals[:,2] * data, **kwargs)

    @property
    def wet_surface_area(self):
        """Returns wet surface area."""
        return self.surface_integral(1)

    @property
    def volumes(self):
        """Returns volumes using x, y, z components of the mesh."""
        norm_coord = self.faces_normals * self.faces_centers
        return self.surface_integral(norm_coord.T, axis=1)

    @property
    def volume(self):
        """Returns volume of the mesh."""
        return np.mean(self.volumes)

    def disp_mass(self, *, rho=1000):
        return rho * self.volume

    @property
    def center_of_buoyancy(self):
        """Returns center of buoyancy of the mesh."""
        coords_sq_norm = self.faces_normals * self.faces_centers**2
        return self.surface_integral(coords_sq_norm.T, axis=1) / (2*self.volume)

    @property
    def waterplane_area(self):
        """Returns water plane area of the mesh."""
        waterplane_area = -self.waterplane_integral(1)
        return waterplane_area

    @property
    def waterplane_center(self):
        """Returns water plane center of the mesh.

        Note: Returns None if the mesh is full submerged.
        """
        waterplane_area = self.waterplane_area
        if abs(waterplane_area) < 1e-10:
            return None
        else:
            waterplane_center = -self.waterplane_integral(
                self.faces_centers.T, axis=1) / waterplane_area
            return waterplane_center[:-1]


# FORMER IMPLEMENTATIONS
# Based on meshmagick <https://github.com/LHEEA/meshmagick> by François Rongère.
# (to be deleted in 2.0?)

def compute_faces_integrals(mesh):
    """Compute intergrals on the faces."""

    surface_integrals = np.zeros((15, mesh.nb_faces), dtype=float)

    # First triangles
    if mesh.nb_triangles > 0:
        triangles_ids = mesh.triangles_ids
        triangles_vertices = mesh.vertices[mesh._faces[triangles_ids][:, :3]]
        surface_integrals[:, triangles_ids] = _compute_triangles_integrals(triangles_vertices)

    # Now quadrangles by splitting them up
    if mesh.nb_quadrangles > 0:
        quadrangles_ids = mesh.quadrangles_ids
        quadrangles = mesh.faces[quadrangles_ids]

        # First pass
        surface_integrals[:, quadrangles_ids] = \
            _compute_triangles_integrals(mesh._vertices[quadrangles[:, (0, 1, 2)]])

        # Second pass
        surface_integrals[:, quadrangles_ids] += \
            _compute_triangles_integrals(mesh._vertices[quadrangles[:, (0, 2, 3)]])

    return surface_integrals


def _compute_triangles_integrals(triangles_vertices):
    """Performs the computation of the various interesting surface integrals.

    Notes
    -----
    triangles vertices must describe by increasing dimension from the general to the particular:
    dimension 0 : information about each facet -- triangles_vertices[0] -> facet 0)
    dimension 1 : information about each facet vertex -- triangles_vertices[0, 1] -> vertex 1 of facet 0
    dimension 2 : information on each vertex coordinate -- triangles_vertices[0, 1, 2] -> z coordinate of vertex 1 of facet 0

    Todo
    ----
    Explicit the integrals
    """

    nb_triangles = triangles_vertices.shape[0]
    s_int = np.zeros((15, nb_triangles), dtype=float)

    point_0, point_1, point_2 = list(map(_3DPointsArray, np.rollaxis(triangles_vertices, 1, 0)))

    t0 = point_0 + point_1                    # t0 = p0 + p1
    f1 = t0 + point_2                         # f1 = p0 + p1 + p2
    t1 = point_0 * point_0                    # t1 = p0**2
    t2 = t1 + point_1*t0                      # t2 = p0**2 + p1**2 + p1p0
    f2 = t2 + point_2*f1                      # f2 = p0**2 + p1**2 + p2**2 + p0*p1 + p0*p2 + p1*p2
    f3 = point_0*t1 + point_1*t2 + point_2*f2 # f3 = p0**3 + p1**3 + p2**3 + p1*p0**2 + ...
    g0 = f2 + point_0 * (f1 + point_0)        # g0 = 3p0**2 + p1**2  + p2**2  + 2p0*p1 + 2p0*p2 + p1*p2
    g1 = f2 + point_1 * (f1 + point_1)        # g1 = p0**2  + 3p1**2 + p2**2  + 2p0*p1 + p0*p2  + 2p1*p2
    g2 = f2 + point_2 * (f1 + point_2)        # g2 = p0**2  + p1**2  + 3p2**2 + p0*p1  + 2p0*p2 + 2p1*p2

    e1 = point_1 - point_0
    e2 = point_2 - point_0

    delta = np.linalg.norm(np.cross(e1, e2), axis=1)

    s_int[0:3] = np.einsum('i, ij -> ji', delta, f1) / 6.

    s_int[3] = delta * (6.*point_0.y*point_0.z + 3*(point_1.y*point_1.z + point_2.y*point_2.z) - point_0.y*f1[:, 2] - point_0.z*f1[:, 1]) / 12.
    s_int[4] = delta * (6.*point_0.x*point_0.z + 3*(point_1.x*point_1.z + point_2.x*point_2.z) - point_0.x*f1[:, 2] - point_0.z*f1[:, 0]) / 12.
    s_int[5] = delta * (6.*point_0.x*point_0.y + 3*(point_1.x*point_1.y + point_2.x*point_2.y) - point_0.x*f1[:, 1] - point_0.y*f1[:, 0]) / 12.

    s_int[6:9] = np.einsum('i, ij -> ji', delta, f2) / 12.
    s_int[9:12] = np.einsum('i, ij -> ji', delta, f3) / 20.

    s_int[12] = delta * (point_0.y*g0[:, 0] + point_1.y*g1[:, 0] + point_2.y*g2[:, 0]) / 60.
    s_int[13] = delta * (point_0.z*g0[:, 1] + point_1.z*g1[:, 1] + point_2.z*g2[:, 1]) / 60.
    s_int[14] = delta * (point_0.x*g0[:, 2] + point_1.x*g1[:, 2] + point_2.x*g2[:, 2]) / 60.

    return s_int


class _3DPointsArray(np.ndarray):
    def __new__(cls, points):
        obj = np.asarray(points).view(cls)
        cls.x = property(fget=lambda cls: cls[:, 0])
        cls.y = property(fget=lambda cls: cls[:, 1])
        cls.z = property(fget=lambda cls: cls[:, 2])
        return obj
