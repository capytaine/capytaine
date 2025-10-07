"""Variants of Delhommeau's method for the computation of the Green function."""
# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import os
import logging
from functools import lru_cache
from importlib import import_module

import numpy as np

from capytaine.tools.prony_decomposition import find_best_exponential_decomposition, PronyDecompositionFailure
from capytaine.tools.cache_on_disk import cache_directory

from capytaine.green_functions.abstract_green_function import AbstractGreenFunction, GreenFunctionEvaluationError

LOG = logging.getLogger(__name__)

_default_parameters = dict(
    tabulation_nr=676,
    tabulation_rmax=100.0,
    tabulation_nz=372,
    tabulation_zmin=-251.0,
    tabulation_nb_integration_points=1001,
    tabulation_grid_shape="scaled_nemoh3",
    finite_depth_method="newer",
    finite_depth_prony_decomposition_method="python",
    floating_point_precision="float64",
    gf_singularities="low_freq",
)


class Delhommeau(AbstractGreenFunction):
    """The Green function as implemented in Aquadyn and Nemoh.

    Parameters
    ----------
    tabulation_nr: int, optional
        Number of tabulation points for horizontal coordinate.
        If 0 is given, no tabulation is used at all.
        Default: 676
    tabulation_rmax: float, optional
        Maximum value of r range for the tabulation. (Minimum is zero.)
        Only used with the :code:`"scaled_nemoh3"` method.
        Default: 100.0
    tabulation_nz: int, optional
        Number of tabulation points for vertical coordinate.
        If 0 is given, no tabulation is used at all.
        Default: 372
    tabulation_zmin: float, optional
        Minimum value of z range for the tabulation. (Maximum is zero.)
        Only used with the :code:`"scaled_nemoh3"` method.
        Default: -251.0
    tabulation_nb_integration_points: int, optional
        Number of points for the numerical integration w.r.t. :math:`theta` of
        Delhommeau's integrals
        Default: 1000
    tabulation_grid_shape: string, optional
        Either :code:`"legacy"` or :code:`"scaled_nemoh3"`, which are the two
        methods currently implemented.
        Default: :code:`"scaled_nemoh3"`
    tabulation_cache_dir: str or None, optional
        Directory in which to save the tabulation file(s).
        If None, the tabulation is not saved on disk.
        Default: calls capytaine.tools.cache_on_disk.cache_directory(), which
        returns the value of the environment variable CAPYTAINE_CACHE_DIR if
        set, or else the default cache directory on your system.
    finite_depth_method: string, optional
        The method used to compute the finite depth Green function.
    finite_depth_prony_decomposition_method: string, optional
        The implementation of the Prony decomposition used to compute the
        finite water_depth Green function. Accepted values: :code:`'fortran'`
        for Nemoh's implementation (by default), :code:`'python'` for an
        experimental Python implementation.
        See :func:`find_best_exponential_decomposition`.
    floating_point_precision: string, optional
        Either :code:`'float32'` for single precision computations or
        :code:`'float64'` for double precision computations.
        Default: :code:`'float64'`.
    gf_singularities: string, optional
        Chose of the variant among the ways singularities can be extracted from
        the Green function. Currently only affects the infinite depth Green
        function.
        Default: "low_freq".

    Attributes
    ----------
    fortran_core:
        Compiled Fortran module with functions used to compute the Green
        function.
    tabulation_grid_shape_index: int
    gf_singularities_index: int
    finite_depth_method_index: int
        Integers passed to Fortran code to describe which method is used.
    tabulated_r_range: numpy.array of shape (tabulation_nr,) and type floating_point_precision
    tabulated_z_range: numpy.array of shape (tabulation_nz,) and type floating_point_precision
        Coordinates of the tabulation points.
    tabulated_integrals: numpy.array of shape (tabulation_nr, tabulation_nz, nb_tabulated_values) and type floating_point_precision
        Tabulated Delhommeau integrals.
    """

    dispersion_relation_roots = np.empty(1)  # dummy array


    def __init__(self, *,
                 tabulation_nr=_default_parameters["tabulation_nr"],
                 tabulation_rmax=_default_parameters["tabulation_rmax"],
                 tabulation_nz=_default_parameters["tabulation_nz"],
                 tabulation_zmin=_default_parameters["tabulation_zmin"],
                 tabulation_nb_integration_points=_default_parameters["tabulation_nb_integration_points"],
                 tabulation_grid_shape=_default_parameters["tabulation_grid_shape"],
                 tabulation_cache_dir=cache_directory(),
                 finite_depth_method=_default_parameters["finite_depth_method"],
                 finite_depth_prony_decomposition_method=_default_parameters["finite_depth_prony_decomposition_method"],
                 floating_point_precision=_default_parameters["floating_point_precision"],
                 gf_singularities=_default_parameters["gf_singularities"],
                 ):

        self.fortran_core = import_module(f"capytaine.green_functions.libs.Delhommeau_{floating_point_precision}")

        self.tabulation_grid_shape = tabulation_grid_shape
        fortran_enum = {
                'legacy': self.fortran_core.constants.legacy_grid,
                'scaled_nemoh3': self.fortran_core.constants.scaled_nemoh3_grid,
                              }
        self.tabulation_grid_shape_index = fortran_enum[tabulation_grid_shape]

        self.gf_singularities = gf_singularities
        self.gf_singularities_fortran_enum = {
                'high_freq': self.fortran_core.constants.high_freq,
                'low_freq': self.fortran_core.constants.low_freq,
                'low_freq_with_rankine_part': self.fortran_core.constants.low_freq_with_rankine_part,
                }

        self.finite_depth_method = finite_depth_method
        fortran_enum = {
                'legacy': self.fortran_core.constants.legacy_finite_depth,
                'newer': self.fortran_core.constants.newer_finite_depth,
                              }
        self.finite_depth_method_index = fortran_enum[finite_depth_method]

        self.floating_point_precision = floating_point_precision
        self.tabulation_nb_integration_points = tabulation_nb_integration_points

        self.tabulation_cache_dir = tabulation_cache_dir
        if tabulation_cache_dir is None:
            self._create_tabulation(tabulation_nr, tabulation_rmax,
                                    tabulation_nz, tabulation_zmin,
                                    tabulation_nb_integration_points)
        else:
            self._create_or_load_tabulation(tabulation_nr, tabulation_rmax,
                                            tabulation_nz, tabulation_zmin,
                                            tabulation_nb_integration_points,
                                            tabulation_cache_dir)

        self.finite_depth_prony_decomposition_method = finite_depth_prony_decomposition_method

        self.exportable_settings = {
            'green_function': self.__class__.__name__,
            'tabulation_nr': tabulation_nr,
            'tabulation_rmax': tabulation_rmax,
            'tabulation_nz': tabulation_nz,
            'tabulation_zmin': tabulation_zmin,
            'tabulation_nb_integration_points': tabulation_nb_integration_points,
            'tabulation_grid_shape': tabulation_grid_shape,
            'finite_depth_method': finite_depth_method,
            'finite_depth_prony_decomposition_method': finite_depth_prony_decomposition_method,
            'floating_point_precision': floating_point_precision,
            'gf_singularities': gf_singularities,
        }

        self._hash = hash(self.exportable_settings.values())

    def __hash__(self):
        return self._hash

    def __str__(self):
        # Print only the non-default values.
        to_be_printed = []
        for name, value in self.exportable_settings.items():
            if name in _default_parameters and value != _default_parameters[name]:
                to_be_printed.append(f"{name}={repr(value)}")
        return f"{self.__class__.__name__}({', '.join(to_be_printed)})"

    def __repr__(self):
        # Same as __str__ except all values are printed even when they are the
        # default value.
        to_be_printed = []
        for name, value in self.exportable_settings.items():
            if name in _default_parameters:
                to_be_printed.append(f"{name}={repr(value)}")
        return f"{self.__class__.__name__}({', '.join(to_be_printed)})"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__repr__())

    def _create_or_load_tabulation(self, tabulation_nr, tabulation_rmax,
                                   tabulation_nz, tabulation_zmin,
                                   tabulation_nb_integration_points,
                                   tabulation_cache_dir):
        """This method either:
            - loads an existing tabulation saved on disk
            - generates a new tabulation with the data provided as argument and save it on disk.
        """

        # Normalize inputs
        tabulation_rmax = float(tabulation_rmax)
        tabulation_zmin = float(tabulation_zmin)

        filename = "tabulation_{}_{}_{}_{}_{}_{}_{}.npz".format(
            self.floating_point_precision, self.tabulation_grid_shape,
            tabulation_nr, tabulation_rmax, tabulation_nz, tabulation_zmin,
            tabulation_nb_integration_points
        )
        filepath = os.path.join(tabulation_cache_dir, filename)

        if os.path.exists(filepath):
            try:
                LOG.info("Loading tabulation from %s", filepath)
                loaded_arrays = np.load(filepath)
                self.tabulated_r_range = loaded_arrays["r_range"]
                self.tabulated_z_range = loaded_arrays["z_range"]
                self.tabulated_integrals = loaded_arrays["values"]
                return filename
            except (EOFError, FileNotFoundError, KeyError, ValueError):
                LOG.warning("Error loading tabulation from %s", filepath)

        self._create_tabulation(tabulation_nr, tabulation_rmax,
                                tabulation_nz, tabulation_zmin,
                                tabulation_nb_integration_points)
        LOG.debug("Saving tabulation in %s", filepath)
        np.savez_compressed(
            filepath, r_range=self.tabulated_r_range, z_range=self.tabulated_z_range,
            values=self.tabulated_integrals
        )
        return filename

    def _create_tabulation(self, tabulation_nr, tabulation_rmax,
                                   tabulation_nz, tabulation_zmin,
                                   tabulation_nb_integration_points):
        LOG.warning("Precomputing tabulation, it may take a few seconds.")
        self.tabulated_r_range = self.fortran_core.delhommeau_integrals.default_r_spacing(
                tabulation_nr, tabulation_rmax, self.tabulation_grid_shape_index
                )
        self.tabulated_z_range = self.fortran_core.delhommeau_integrals.default_z_spacing(
                tabulation_nz, tabulation_zmin, self.tabulation_grid_shape_index
                )
        self.tabulated_integrals = self.fortran_core.delhommeau_integrals.construct_tabulation(
                self.tabulated_r_range, self.tabulated_z_range, tabulation_nb_integration_points,
                )

    @property
    def all_tabulation_parameters(self):
        """An alias meant to pass to the Fortran functions all the parameters controlling the tabulation in a single item."""
        return (self.tabulation_nb_integration_points, self.tabulation_grid_shape_index,
                self.tabulated_r_range, self.tabulated_z_range, self.tabulated_integrals)

    @lru_cache(maxsize=128)
    def find_best_exponential_decomposition(self, dimensionless_wavenumber, *, method=None):
        """Compute the decomposition of a part of the finite water_depth Green function as a sum of exponential functions.

        Two implementations are available: the legacy Fortran implementation from Nemoh and a newer one written in Python.
        For some still unexplained reasons, the two implementations do not always give the exact same result.
        Until the problem is better understood, the Fortran implementation is the default one, to ensure consistency with Nemoh.
        The Fortran version is also significantly faster...

        Results are cached.

        Parameters
        ----------
        dimensionless_wavenumber: float
            dimensionless wavenumber: :math:`kh`
        method: str, optional
            "python" or "fortran". If not provided, uses self.finite_depth_prony_decomposition_method.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            the amplitude and growth rates of the exponentials
        """
        kh = dimensionless_wavenumber

        if method is None:
            method = self.finite_depth_prony_decomposition_method

        LOG.debug("\tCompute Prony decomposition in finite water_depth Green function "
                  "for dimensionless_wavenumber=%.2e", dimensionless_wavenumber)

        if method.lower() == 'python':
            if kh <= 0.1:
                raise NotImplementedError(
                    f"{self} cannot evaluate finite depth Green function "
                    f"for kh<0.1 (kh={kh})"
                    )
            elif kh < 1e5:
                # The function that will be approximated.
                sing_coef = (1 + np.tanh(kh))**2/(1 - np.tanh(kh)**2 + np.tanh(kh)/kh)
                def ref_function(x):
                    """The function that should be approximated by a sum of exponentials."""
                    return ((x + kh*np.tanh(kh)) * np.exp(x))/(x*np.sinh(x) - kh*np.tanh(kh)*np.cosh(x)) - sing_coef/(x - kh) - 2
            else:
                # Asymptotic approximation of the function for large kh, including infinite frequency
                def ref_function(x):
                    return -2/(1 + np.exp(-2*x)) + 2

            try:
                a, lamda = find_best_exponential_decomposition(ref_function, x_min=-0.1, x_max=20.0, n_exp_range=range(4, 31, 2), tol=1e-4)
                return np.stack([lamda, a])
            except PronyDecompositionFailure as e:
                raise GreenFunctionEvaluationError(
                    f"{self} cannot evaluate finite depth Green function "
                    f"for kh={dimensionless_wavenumber}"
                ) from e

        elif method.lower() == 'fortran':
            if kh > 1e5:
                raise NotImplementedError("Fortran implementation of the Prony decomposition does not support infinite frequency")
            omega2_h_over_g = kh*np.tanh(kh)
            nexp, pr_d = self.fortran_core.old_prony_decomposition.lisc(omega2_h_over_g, kh)
            return pr_d[0:2, :nexp]

        else:
            raise ValueError(f"Unrecognized name for the Prony decomposition method: {repr(method)}. Expected 'python' or 'fortran'.")

    def evaluate_rankine_only(self,
                              mesh1, mesh2,
                              adjoint_double_layer=True, early_dot_product=True
                              ):
        r"""Construct the matrices between mesh1 (that can also be a list of points)
        and mesh2 for a Rankine kernel.

        Parameters
        ----------
        mesh1: Mesh or CollectionOfMeshes or list of points
            mesh of the receiving body (where the potential is measured)
            if only S is wanted or early_dot_product is False, then only a list
            of points as an array of shape (n, 3) can be passed.
        mesh2: Mesh or CollectionOfMeshes
            mesh of the source body (over which the source distribution is integrated)
        adjoint_double_layer: bool, optional
            compute double layer for direct method (F) or adjoint double layer
            for indirect method (T) matrices (default: True)
        early_dot_product: boolean, optional
            if False, return K as a (n, m, 3) array storing ∫∇G
            if True, return K as a (n, m) array storing ∫∇G·n

        Returns
        -------
        tuple of real-valued numpy arrays
            the matrices :math:`S` and :math:`K`
        """
        collocation_points, early_dot_product_normals = \
                self._get_colocation_points_and_normals(mesh1, mesh2, adjoint_double_layer)

        S, K = self._init_matrices(
            (collocation_points.shape[0], mesh2.nb_faces), early_dot_product
        )

        self.fortran_core.matrices.add_rankine_term_only(
                collocation_points,  early_dot_product_normals,
                mesh2.vertices,      mesh2.faces + 1,
                mesh2.faces_centers, mesh2.faces_normals,
                mesh2.faces_areas,   mesh2.faces_radiuses,
                *mesh2.quadrature_points,
                adjoint_double_layer,
                S, K)

        if mesh1 is mesh2:
            self.fortran_core.matrices.add_diagonal_term(
                    mesh2.faces_centers, early_dot_product_normals, np.inf, K,
                    )

        S, K = np.real(S), np.real(K)

        if np.any(np.isnan(S)) or np.any(np.isnan(K)):
            raise GreenFunctionEvaluationError(
                    "Green function returned a NaN in the interaction matrix.\n"
                    "It could be due to overlapping panels.")

        if early_dot_product:
            K = K.reshape((collocation_points.shape[0], mesh2.nb_faces))

        return S, K


    def evaluate(self,
                 mesh1, mesh2,
                 free_surface=0.0, water_depth=np.inf, wavenumber=1.0,
                 adjoint_double_layer=True, early_dot_product=True
                 ):
        r"""The main method of the class, called by the engine to assemble the influence matrices.

        Parameters
        ----------
        mesh1: MeshLike or list of points
            mesh of the receiving body (where the potential is measured)
            if only S is wanted or early_dot_product is False, then only a list of points as an array of shape (n, 3) can be passed.
        mesh2: MeshLike
            mesh of the source body (over which the source distribution is integrated)
        free_surface: float, optional
            position of the free surface (default: :math:`z = 0`)
        water_depth: float, optional
            constant depth of water (default: :math:`+\infty`)
        wavenumber: float, optional
            wavenumber (default: 1.0)
        adjoint_double_layer: bool, optional
            compute double layer for direct method (F) or adjoint double layer for indirect method (T) matrices (default: True)
        early_dot_product: boolean, optional
            if False, return K as a (n, m, 3) array storing ∫∇G
            if True, return K as a (n, m) array storing ∫∇G·n

        Returns
        -------
        tuple of numpy arrays
            the matrices :math:`S` and :math:`K`
            the dtype of the matrix can be real or complex and depends on self.floating_point_precision
        """

        if free_surface == np.inf:  # No free surface, only a single Rankine source term
            if water_depth != np.inf:
                raise ValueError("When setting free_surface=inf, "
                                 "the water depth should also be infinite "
                                 f"(got water_depth={water_depth})")

            return self.evaluate_rankine_only(
                mesh1, mesh2,
                adjoint_double_layer=adjoint_double_layer,
                early_dot_product=early_dot_product,
            )

        # Main case:
        collocation_points, early_dot_product_normals = \
                self._get_colocation_points_and_normals(mesh1, mesh2, adjoint_double_layer)

        S, K = self._init_matrices(
            (collocation_points.shape[0], mesh2.nb_faces), early_dot_product
        )

        wavenumber = float(wavenumber)

        # Overrides gf_singularities setting in some specific cases, else use the class one.
        if water_depth < np.inf and self.finite_depth_method == 'legacy' and not self.gf_singularities == 'low_freq':
            gf_singularities = "low_freq"  # Reproduce legacy method behavior
            LOG.debug(
                f"Overriding gf_singularities='{self.gf_singularities}' because of finite_depth_method=='legacy'"
            )
        elif wavenumber == 0.0 and not self.gf_singularities == 'low_freq':
            gf_singularities = "low_freq"
            LOG.debug(
                f"Overriding gf_singularities='{self.gf_singularities}' because of wavenumber==0.0"
            )
        elif wavenumber == np.inf and not self.gf_singularities == 'high_freq':
            gf_singularities = "high_freq"
            LOG.debug(
                f"Overriding gf_singularities='{self.gf_singularities}' because of wavenumber==np.inf"
            )
        elif np.any(abs(mesh2.faces_centers[:, 2]) < 1e-6) and not self.gf_singularities == 'low_freq':
            gf_singularities = "low_freq"
            LOG.warning(
                f"Overriding gf_singularities='{self.gf_singularities}' because of free surface panels, "
                "which are currently only supported by gf_singularities='low_freq'"
            )
        else:
            gf_singularities = self.gf_singularities
        gf_singularities_index = self.gf_singularities_fortran_enum[gf_singularities]

        if water_depth == np.inf:
            prony_decomposition = np.zeros((1, 1))  # Dummy array that won't actually be used by the fortran code.
        else:
            prony_decomposition = self.find_best_exponential_decomposition(wavenumber*water_depth)

        # Main call to Fortran code
        self.fortran_core.matrices.build_matrices(
            collocation_points,  early_dot_product_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
            *mesh2.quadrature_points,
            wavenumber, water_depth,
            *self.all_tabulation_parameters,
            self.finite_depth_method_index, prony_decomposition, self.dispersion_relation_roots,
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

################################

class XieDelhommeau(Delhommeau):
    """Legacy way to call the gf_singularities="low_freq" variant."""

    def __init__(self, **kwargs):
        super().__init__(gf_singularities="low_freq", **kwargs)
