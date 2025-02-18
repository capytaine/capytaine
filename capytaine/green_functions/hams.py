from importlib import import_module
import numpy as np

from capytaine.green_functions.abstract_green_function import AbstractGreenFunction


class LiangWuNoblesseGF(AbstractGreenFunction):
    """Wrapper for the Fortran implementation of the infinite depth Green function of [Liang, Wu, Noblesse, 2018].

    Uses the same implementation as Delhommeau() for the Rankine and reflected Rankine terms.

    """
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

    def __str__(self):
        return "LiangWuNoblesseGF()"

    def __repr__(self):
        return "LiangWuNoblesseGF()"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__repr__())

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

        collocation_points, early_dot_product_normals = self._get_colocation_points_and_normals(mesh1, mesh2, adjoint_double_layer)

        S, K = self._init_matrices((collocation_points.shape[0], mesh2.nb_faces), "complex128", early_dot_product=early_dot_product)

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

        if early_dot_product: K = K.reshape((collocation_points.shape[0], mesh2.nb_faces))

        return S, K

