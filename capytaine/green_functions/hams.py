from importlib import import_module
import numpy as np

from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.green_functions.abstract_green_function import AbstractGreenFunction


class LiangWuNoblesseGF(AbstractGreenFunction):
    fortran_core = import_module("capytaine.green_functions.libs.Delhommeau_float64")
    tabulation_grid_shape_index = fortran_core.constants.liang_wu_noblesse
    gf_singularities_index = fortran_core.constants.low_freq
    exportable_settings = {'green_function': "LiangWuNoblesseGF"}

    # Dummy arrays that won't actually be used by the fortran code.
    a_exp, lamda_exp = np.empty(1), np.empty(1)
    tabulation_nb_integration_points = 1
    tabulated_r_range = np.empty(1)
    tabulated_z_range = np.empty(1)
    tabulated_integrals = np.empty(1)

    def evaluate(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=True, early_dot_product=True):
        if free_surface == np.inf: # No free surface, only a single Rankine source term
            coeffs = np.array((1.0, 0.0, 0.0))
        elif free_surface == 0.0 and water_depth < np.inf:
            raise NotImplementedError()
        else:  # Infinite water depth
            if wavenumber == 0.0:
                coeffs = np.array((1.0, 1.0, 0.0))
            elif wavenumber == np.inf:
                coeffs = np.array((1.0, -1.0, 0.0))
            else:
                coeffs = np.array((1.0, 1.0, 1.0))

        if isinstance(mesh1, Mesh) or isinstance(mesh1, CollectionOfMeshes):
            collocation_points = mesh1.faces_centers
            nb_collocation_points = mesh1.nb_faces
            if not adjoint_double_layer: # Computing the D matrix
                early_dot_product_normals = mesh2.faces_normals
            else: # Computing the K matrix
                early_dot_product_normals = mesh1.faces_normals

        elif isinstance(mesh1, np.ndarray) and mesh1.ndim == 2 and mesh1.shape[1] == 3:
            # This is used when computing potential or velocity at given points in postprocessing
            collocation_points = mesh1
            nb_collocation_points = mesh1.shape[0]
            if not adjoint_double_layer: # Computing the D matrix
                early_dot_product_normals = mesh2.faces_normals
            else: # Computing the K matrix
                early_dot_product_normals = np.zeros((nb_collocation_points, 3))
                # Dummy argument since this method is meant to be used either
                # - to compute potential, then only S is needed and early_dot_product_normals is irrelevant,
                # - to compute velocity, then the adjoint full gradient is needed and early_dot_product is False and this value is unused.
                # TODO: add an only_S argument and return an error here if (early_dot_product and not only_S)

        else:
            raise ValueError(f"Unrecognized first input for {self.__class__.__name__}.evaluate:\n{mesh1}")

        dtype = "complex128"
        S = np.empty((nb_collocation_points, mesh2.nb_faces), order="F", dtype=dtype)
        K = np.empty((nb_collocation_points, mesh2.nb_faces, 1 if early_dot_product else 3), order="F", dtype=dtype)

        self.fortran_core.matrices.build_matrices(
            collocation_points,  early_dot_product_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
            *mesh2.quadrature_points,
            wavenumber, np.inf,
            coeffs,
            self.tabulation_nb_integration_points, self.tabulation_grid_shape_index,
            self.tabulated_r_range, self.tabulated_z_range, self.tabulated_integrals,
            self.lamda_exp, self.a_exp,
            mesh1 is mesh2, self.gf_singularities_index, adjoint_double_layer,
            S, K
        )

        if np.any(np.isnan(S)) or np.any(np.isnan(K)):
            raise RuntimeError("Green function returned a NaN in the interaction matrix.\n"
                    "It could be due to overlapping panels.")

        if early_dot_product: K = K.reshape((nb_collocation_points, mesh2.nb_faces))

        return S, K

