import logging

import numpy as np

from capytaine.tools.optional_imports import silently_import_optional_dependency

LOG = logging.getLogger(__name__)

# The builtin methods are stored as a list of 2D-points in [-1, 1]² and a list
# of corresponding weights. The 2D points will be remapped to the actual shape
# of the faces. They are only defined for quadrilaterals. They also work for
# triangles (although they might be subobtimal).

builtin_methods = {
        "First order": (np.array([(0.0, 0.0)]), np.array([1.0])),
        "Gauss-Legendre 2": (
            np.array([(+1/np.sqrt(3), +1/np.sqrt(3)),
             (+1/np.sqrt(3), -1/np.sqrt(3)),
             (-1/np.sqrt(3), +1/np.sqrt(3)),
             (-1/np.sqrt(3), -1/np.sqrt(3))]),
            np.array([1/4, 1/4, 1/4, 1/4])
            )
        }


def compute_quadrature_on_faces(faces, method):
    """
    Compute the quadrature points and weight for numerical integration over the faces.

    Parameters
    ----------
    faces: array of shape (nb_faces, 4, 3)
        The 3D-coordinates of each of the 4 corners of each face.
    method: string or quadpy object
        The method used to compute the quadrature scheme

    Returns
    -------
    points: array of shape (nb_faces, nb_quad_points, 3)
        The 3D-coordinates of each of the quadrature points (their number depends on the method) of each face.
    weights: array of shape (nb_faces, nb_quad_points)
        Weights associated to each of the quadrature points.
    """

    if method in builtin_methods:
        LOG.debug("Quadrature method found in builtin methods.")
        local_points, local_weights = builtin_methods[method]

    elif ((quadpy := silently_import_optional_dependency("quadpy")) is not None
                and isinstance(method, quadpy.c2._helpers.C2Scheme)):
        LOG.debug("Quadrature method is a Quadpy scheme: %s", method.name)
        local_points = method.points.T
        local_weights = method.weights

    else:
        raise ValueError(f"Unrecognized quadrature scheme: {method}.\n"
                         f"Consider using one of the following: {set(builtin_methods.keys())}")

    nb_faces = faces.shape[0]
    nb_quad_points = len(local_weights)
    points = np.empty((nb_faces, nb_quad_points, 3))
    weights = np.empty((nb_faces, nb_quad_points))

    for i_face in range(nb_faces):
        for k_quad in range(nb_quad_points):
            xk, yk = local_points[k_quad, :]
            points[i_face, k_quad, :] = (
                      (1+xk)*(1+yk) * faces[i_face, 0, :]
                    + (1+xk)*(1-yk) * faces[i_face, 1, :]
                    + (1-xk)*(1-yk) * faces[i_face, 2, :]
                    + (1-xk)*(1+yk) * faces[i_face, 3, :]
                    )/4
            dxidx = ((1+yk)*faces[i_face, 0, :] + (1-yk)*faces[i_face, 1, :]
                     - (1-yk)*faces[i_face, 2, :] - (1+yk)*faces[i_face, 3, :])/4
            dxidy = ((1+xk)*faces[i_face, 0, :] - (1+xk)*faces[i_face, 1, :]
                     - (1-xk)*faces[i_face, 2, :] + (1-xk)*faces[i_face, 3, :])/4
            detJ = np.linalg.norm(np.cross(dxidx, dxidy))
            weights[i_face, k_quad] = local_weights[k_quad] * 4 * detJ

    return points, weights
