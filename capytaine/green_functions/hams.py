from importlib import import_module
from scipy.optimize import brentq
import numpy as np

from capytaine.green_functions.abstract_green_function import AbstractGreenFunction, GreenFunctionEvaluationError


class LiangWuNoblesseGF(AbstractGreenFunction):
    """Wrapper for the Fortran implementation of the infinite depth Green function of [Liang, Wu, Noblesse, 2018].

    Uses the same implementation as Delhommeau() for the Rankine and reflected Rankine terms.

    """
    floating_point_precision = "float64"

    fortran_core = import_module("capytaine.green_functions.libs.Delhommeau_float64")
    tabulation_grid_shape_index = fortran_core.constants.liang_wu_noblesse
    exportable_settings = {'green_function': "LiangWuNoblesseGF"}

    # Dummy arrays that won't actually be used by the fortran code.
    prony_decomposition = np.zeros((1, 1))
    dispersion_relation_roots = np.empty(1)
    finite_depth_method_index = -9999
    tabulation_nb_integration_points = 1
    tabulated_r_range = np.empty(1)
    tabulated_z_range = np.empty(1)
    tabulated_integrals = np.empty(1)
    dummy_param = -999

    def __str__(self):
        return "LiangWuNoblesseGF()"

    def __repr__(self):
        return "LiangWuNoblesseGF()"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__repr__())

    def evaluate(self,
                 mesh1, mesh2,
                 free_surface=0.0, water_depth=np.inf, wavenumber=1.0,
                 adjoint_double_layer=True, early_dot_product=True
                 ):

        if free_surface == np.inf or water_depth < np.inf:
            raise NotImplementedError("LiangWuNoblesseGF() is only implemented for infinite depth with a free surface")

        if wavenumber == np.inf:
            gf_singularities_index = self.fortran_core.constants.high_freq
        else:
            gf_singularities_index = self.fortran_core.constants.low_freq

        collocation_points, early_dot_product_normals = \
                self._get_colocation_points_and_normals(mesh1, mesh2, adjoint_double_layer)

        S, K = self._init_matrices(
            (collocation_points.shape[0], mesh2.nb_faces), early_dot_product=early_dot_product
        )

        self.fortran_core.matrices.build_matrices(
            collocation_points,  early_dot_product_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
            *mesh2.quadrature_points,
            wavenumber, np.inf,
            self.tabulation_nb_integration_points, self.tabulation_grid_shape_index,
            self.tabulated_r_range, self.tabulated_z_range, self.tabulated_integrals,
            self.dummy_param, self.prony_decomposition, self.dispersion_relation_roots,
            gf_singularities_index, adjoint_double_layer,
            S, K
        )

        if mesh1 is mesh2:
            self.fortran_core.matrices.add_diagonal_term(
                    mesh2.faces_centers, early_dot_product_normals, free_surface, K,
                    )

        if np.any(np.isnan(S)) or np.any(np.isnan(K)):
            raise GreenFunctionEvaluationError(
                    "Green function returned a NaN in the interaction matrix.\n"
                    "It could be due to overlapping panels.")

        if early_dot_product:
            K = K.reshape((collocation_points.shape[0], mesh2.nb_faces))

        return S, K


class FinGreen3D(AbstractGreenFunction):
    """Wrapper for the Fortran implementation of the finite depth Green function of [Liu et al.].

    Uses the same implementation as Delhommeau() for the Rankine and reflected Rankine terms.

    """
    floating_point_precision = "float64"

    fortran_core = import_module("capytaine.green_functions.libs.Delhommeau_float64")
    finite_depth_method_index = fortran_core.constants.fingreen3d_method
    gf_singularities_index = fortran_core.constants.low_freq

    # Dummy arrays that won't actually be used by the fortran code.
    prony_decomposition = np.zeros((1, 1))
    tabulation_nb_integration_points = 1
    tabulated_r_range = np.empty(1)
    tabulated_z_range = np.empty(1)
    tabulated_integrals = np.empty(1)
    dummy_param = -999

    def __init__(self, *, nb_dispersion_roots=200):
        self.nb_dispersion_roots = nb_dispersion_roots
        self.exportable_settings = {
            'green_function': "FinGreen3D",
            'nb_dispersion_roots': nb_dispersion_roots
        }

    def __str__(self):
        return f"FinGreen3D(nb_dispersion_roots={self.nb_dispersion_roots})"

    def __repr__(self):
        return f"FinGreen3D(nb_dispersion_roots={self.nb_dispersion_roots})"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__repr__())

    def compute_dispersion_relation_roots(self, nk, wavenumber, depth):
        omega2_h_over_g = wavenumber*np.tanh(wavenumber*depth)*depth
        def root(i_root):
            return brentq(lambda y: omega2_h_over_g + y*np.tan(y), (2*i_root+1)*np.pi/2 + 1e-10, (2*i_root+2)*np.pi/2 - 1e-10)/depth
        return np.array([wavenumber] + [root(i_root) for i_root in range(nk-1)])

    def evaluate(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=True, early_dot_product=True):

        if free_surface == np.inf or water_depth == np.inf:
            raise NotImplementedError("FinGreen3D is only implemented for finite depth with a free surface.")
        if wavenumber == 0.0 or wavenumber == np.inf:
            raise NotImplementedError("FinGreen3D is only implemented for non-zero and non-infinite frequencies")

        dispersion_relation_roots = self.compute_dispersion_relation_roots(
            self.nb_dispersion_roots,
            wavenumber,
            water_depth
        )

        collocation_points, early_dot_product_normals = \
                self._get_colocation_points_and_normals(mesh1, mesh2, adjoint_double_layer)

        S, K = self._init_matrices(
            (collocation_points.shape[0], mesh2.nb_faces), early_dot_product=early_dot_product
        )

        self.fortran_core.matrices.build_matrices(
            collocation_points,  early_dot_product_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
            *mesh2.quadrature_points,
            wavenumber, water_depth,
            self.tabulation_nb_integration_points, self.dummy_param,
            self.tabulated_r_range, self.tabulated_z_range, self.tabulated_integrals,
            self.finite_depth_method_index, self.prony_decomposition, dispersion_relation_roots,
            self.gf_singularities_index, adjoint_double_layer,
            S, K
        )

        if mesh1 is mesh2:
            self.fortran_core.matrices.add_diagonal_term(
                    mesh2.faces_centers, early_dot_product_normals, free_surface, K,
                    )

        if np.any(np.isnan(S)) or np.any(np.isnan(K)):
            raise GreenFunctionEvaluationError(
                    "Green function returned a NaN in the interaction matrix.\n"
                    "It could be due to overlapping panels.")

        if early_dot_product:
            K = K.reshape((collocation_points.shape[0], mesh2.nb_faces))

        return S, K


class HAMS_GF(AbstractGreenFunction):
    floating_point_precision = "float64"

    exportable_settings = {'green_function': "HAMS_GF"}

    def __init__(self):
        self.infinite_depth_gf = LiangWuNoblesseGF()
        self.finite_depth_gf = FinGreen3D(nb_dispersion_roots=200)

    def __str__(self):
        return "HAMS_GF()"

    def __repr__(self):
        return "HAMS_GF()"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__repr__())

    def evaluate(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=True, early_dot_product=True):
        if water_depth == np.inf:
            return self.infinite_depth_gf.evaluate(mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer, early_dot_product)
        else:
            return self.finite_depth_gf.evaluate(mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer, early_dot_product)
