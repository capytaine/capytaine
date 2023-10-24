import logging

import numpy as np

from capytaine.tools.optional_imports import silently_import_optional_dependency

LOG = logging.getLogger(__name__)

# Quadrature methods are only defined for quadrilaterals. They also work for
# triangles (although they might be subobtimal).

builtin_methods = {
        "Gauss-Legendre 2": (
            np.array([(+1/np.sqrt(3), +1/np.sqrt(3)),
             (+1/np.sqrt(3), -1/np.sqrt(3)),
             (-1/np.sqrt(3), +1/np.sqrt(3)),
             (-1/np.sqrt(3), -1/np.sqrt(3))]),
            np.array([1/4, 1/4, 1/4, 1/4])
            )
        }

DEFAULT_QUADRATURE = "Gauss-Legendre 2"


def compute_quadrature(mesh, method):
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

    nb_points = len(local_weights)
    points = np.empty((mesh.nb_faces, nb_points, 3))
    weights = np.empty((mesh.nb_faces, nb_points))

    for i_face in range(mesh.nb_faces):
        corners = mesh.vertices[mesh.faces[i_face, :]]

        for k_quad in range(nb_points):
            xk, yk = local_points[k_quad, :]
            points[i_face, k_quad, :] = (
                      (1 + xk)*(1 + yk) * corners[0, :]
                    + (1 + xk)*(1 - yk) * corners[1, :]
                    + (1 - xk)*(1 - yk) * corners[2, :]
                    + (1 - xk)*(1 + yk) * corners[3, :]
                    )/4
            dxidx = ((1+yk)*corners[0, :] + (1-yk)*corners[1, :] - (1-yk)*corners[2, :] - (1+yk)*corners[3, :])/4
            dxidy = ((1+xk)*corners[0, :] - (1+xk)*corners[1, :] - (1-xk)*corners[2, :] + (1-xk)*corners[3, :])/4
            detJ = np.linalg.norm(np.cross(dxidx, dxidy))
            weights[i_face, k_quad] = local_weights[k_quad] * 4 * detJ

    return points, weights
