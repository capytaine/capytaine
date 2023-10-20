import numpy as np
from capytaine.bodies import FloatingBody
from capytaine.post_pro.free_surfaces import FreeSurface
from capytaine.meshes import Mesh, CollectionOfMeshes


def _normalize_points(points, keep_mesh=False):
    if isinstance(points, (FloatingBody, FreeSurface)):
        if keep_mesh:
            return points.mesh, (points.mesh.nb_faces,)
        else:
            return points.mesh.faces_centers, (points.mesh.nb_faces,)

    if isinstance(points, (Mesh, CollectionOfMeshes)):
        if keep_mesh:
            return points, (points.nb_faces,)
        else:
            return points.faces_centers, (points.nb_faces,)

    points = np.asarray(points)

    if points.ndim == 1:  # A single point has been provided
        output_shape = (1,)
        points = points.reshape((1, points.shape[0]))

    elif points.ndim == 2:
        output_shape = (points.shape[0],)

    elif points.ndim > 2:
        # `points` is expected to be the results of a meshgrid. Points has shape (d, nx, ny, ...)
        output_shape = points.shape[1:]
        points = points.reshape(points.shape[0], -1).transpose()
        # points is now a (nx*ny*... , d) array

    else:
        raise ValueError(f"Expected a list of points or a mesh, but got instead: {points}")

    return points, output_shape

def _normalize_free_surface_points(points, keep_mesh=False):
    if keep_mesh and isinstance(points, (FloatingBody, FreeSurface)):
        return points.mesh, (points.mesh.nb_faces,)

    if keep_mesh and isinstance(points, (Mesh, CollectionOfMeshes)):
        return points, (points.nb_faces,)

    points, output_shape = _normalize_points(points, keep_mesh)

    if points.ndim == 2 and points.shape[1] == 2:  # Only x and y have been provided
        points = np.concatenate([points, np.zeros((points.shape[0], 1))], axis=1)

    return points, output_shape
