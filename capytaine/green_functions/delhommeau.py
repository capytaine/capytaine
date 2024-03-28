"""Variants of Delhommeau's method for the computation of the Green function."""
# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import os
import logging
from functools import lru_cache
from importlib import import_module

import numpy as np

from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.tools.prony_decomposition import exponential_decomposition, error_exponential_decomposition
from capytaine.tools.cache_on_disk import cache_directory

from capytaine.green_functions.abstract_green_function import AbstractGreenFunction

LOG = logging.getLogger(__name__)

_default_parameters = dict(
    tabulation_nr=676,
    tabulation_rmax=100.0,
    tabulation_nz=372,
    tabulation_zmin=-251.0,
    tabulation_nb_integration_points=1000,
    tabulation_method="scaled_nemoh3",
    finite_depth_prony_decomposition_method="fortran",
    floating_point_precision="float64"
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
    tabulation_method: string, optional
        Either :code:`"legacy"` or :code:`"scaled_nemoh3"`, which are the two
        methods currently implemented.
        Default: :code:`"scaled_nemoh3"`
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

    Attributes
    ----------
    fortran_core:
        Compiled Fortran module with functions used to compute the Green
        function.
    tabulation_method_index: int
        Integer passed to Fortran code to describe which method is used.
    tabulated_r_range: numpy.array of shape (tabulation_nr,) and type floating_point_precision
    tabulated_z_range: numpy.array of shape (tabulation_nz,) and type floating_point_precision
        Coordinates of the tabulation points.
    tabulated_integrals: numpy.array of shape (tabulation_nr, tabulation_nz, 2, 2) and type floating_point_precision
        Tabulated Delhommeau integrals.
    """

    legacy_wave_part = True


    def __init__(self, *,
                 tabulation_nr=_default_parameters["tabulation_nr"],
                 tabulation_rmax=_default_parameters["tabulation_rmax"],
                 tabulation_nz=_default_parameters["tabulation_nz"],
                 tabulation_zmin=_default_parameters["tabulation_zmin"],
                 tabulation_nb_integration_points=_default_parameters["tabulation_nb_integration_points"],
                 tabulation_method=_default_parameters["tabulation_method"],
                 finite_depth_prony_decomposition_method=_default_parameters["finite_depth_prony_decomposition_method"],
                 floating_point_precision=_default_parameters["floating_point_precision"],
                 ):


        self.fortran_core = import_module(f"capytaine.green_functions.libs.Delhommeau_{floating_point_precision}")

        self.tabulation_method = tabulation_method
        fortran_indices_for_methods = {
                'legacy': self.fortran_core.delhommeau_integrals.legacy_method,
                'scaled_nemoh3': self.fortran_core.delhommeau_integrals.scaled_nemoh3_method,
                              }
        self.tabulation_method_index = fortran_indices_for_methods[tabulation_method]

        self.floating_point_precision = floating_point_precision

        self._create_or_load_tabulation(tabulation_nr, tabulation_rmax, tabulation_nz, tabulation_zmin, tabulation_nb_integration_points)

        self.finite_depth_prony_decomposition_method = finite_depth_prony_decomposition_method

        self.exportable_settings = {
            'green_function': self.__class__.__name__,
            'tabulation_nr': tabulation_nr,
            'tabulation_rmax': tabulation_rmax,
            'tabulation_nz': tabulation_nz,
            'tabulation_zmin': tabulation_zmin,
            'tabulation_nb_integration_points': tabulation_nb_integration_points,
            'tabulation_method': tabulation_method,
            'finite_depth_prony_decomposition_method': finite_depth_prony_decomposition_method,
            'floating_point_precision': floating_point_precision,
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
                                   tabulation_nb_integration_points):
        """This method either:
            - loads an existing tabulation saved on disk
            - generates a new tabulation with the data provided as argument and save it on disk.
        """

        # Normalize inputs
        tabulation_rmax = float(tabulation_rmax)
        tabulation_zmin = float(tabulation_zmin)

        filename = "tabulation_{}_{}_{}_{}_{}_{}_{}_{}.npz".format(
            self.__class__.__name__, self.floating_point_precision,
            self.tabulation_method,
            tabulation_nr, tabulation_rmax, tabulation_nz, tabulation_zmin,
            tabulation_nb_integration_points
        )
        filepath = os.path.join(cache_directory(), filename)

        if os.path.exists(filepath):
            LOG.info("Loading tabulation from %s", filepath)
            loaded_arrays = np.load(filepath)
            self.tabulated_r_range = loaded_arrays["r_range"]
            self.tabulated_z_range = loaded_arrays["z_range"]
            self.tabulated_integrals = loaded_arrays["values"]

        else:
            LOG.warning("Precomputing tabulation, it may take a few seconds.")
            self.tabulated_r_range = self.fortran_core.delhommeau_integrals.default_r_spacing(
                    tabulation_nr, tabulation_rmax, self.tabulation_method_index
                    )
            self.tabulated_z_range = self.fortran_core.delhommeau_integrals.default_z_spacing(
                    tabulation_nz, tabulation_zmin, self.tabulation_method_index
                    )
            self.tabulated_integrals = self.fortran_core.delhommeau_integrals.construct_tabulation(
                    self.tabulated_r_range, self.tabulated_z_range, tabulation_nb_integration_points,
                    self.legacy_wave_part,
                    )
            LOG.debug("Saving tabulation in %s", filepath)
            np.savez_compressed(
                filepath, r_range=self.tabulated_r_range, z_range=self.tabulated_z_range,
                values=self.tabulated_integrals
            )

    @lru_cache(maxsize=128)
    def find_best_exponential_decomposition(self, dimensionless_omega, dimensionless_wavenumber):
        """Compute the decomposition of a part of the finite water_depth Green function as a sum of exponential functions.

        Two implementations are available: the legacy Fortran implementation from Nemoh and a newer one written in Python.
        For some still unexplained reasons, the two implementations do not always give the exact same result.
        Until the problem is better understood, the Fortran implementation is the default one, to ensure consistency with Nemoh.
        The Fortran version is also significantly faster...

        Results are cached.

        Parameters
        ----------
        dimensionless_omega: float
            dimensionless angular frequency: :math:`kh \\tanh (kh) = \\omega^2 h/g`
        dimensionless_wavenumber: float
            dimensionless wavenumber: :math:`kh`
        method: string, optional
            the implementation that should be used to compute the Prony decomposition

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            the amplitude and growth rates of the exponentials
        """

        LOG.debug(f"\tCompute Prony decomposition in finite water_depth Green function "
                  f"for dimless_omega=%.2e and dimless_wavenumber=%.2e",
                  dimensionless_omega, dimensionless_wavenumber)

        if self.finite_depth_prony_decomposition_method.lower() == 'python':
            # The function that will be approximated.
            @np.vectorize
            def f(x):
                return self.fortran_core.initialize_green_wave.ff(x, dimensionless_omega, dimensionless_wavenumber)

            # Try different increasing number of exponentials
            for n_exp in range(4, 31, 2):

                # The coefficients are computed on a resolution of 4*n_exp+1 ...
                X = np.linspace(-0.1, 20.0, 4*n_exp+1)
                a, lamda = exponential_decomposition(X, f(X), n_exp)

                # ... and they are evaluated on a finer discretization.
                X = np.linspace(-0.1, 20.0, 8*n_exp+1)
                if error_exponential_decomposition(X, f(X), a, lamda) < 1e-4:
                    break

            else:
                LOG.warning("No suitable exponential decomposition has been found"
                            "for dimless_omega=%.2e and dimless_wavenumber=%.2e",
                            dimensionless_omega, dimensionless_wavenumber)

        elif self.finite_depth_prony_decomposition_method.lower() == 'fortran':
            lamda, a, nexp = self.fortran_core.old_prony_decomposition.lisc(dimensionless_omega, dimensionless_wavenumber)
            lamda = lamda[:nexp]
            a = a[:nexp]

        else:
            raise ValueError("Unrecognized method name for the Prony decomposition.")

        # Add one more exponential function (actually a constant).
        # It is not clear where it comes from exactly in the theory...
        a = np.concatenate([a, np.array([2])])
        lamda = np.concatenate([lamda, np.array([0.0])])

        return a, lamda

    def evaluate(self, mesh1, mesh2, free_surface=0.0, water_depth=np.inf, wavenumber=1.0, adjoint_double_layer=True, early_dot_product=True):
        r"""The main method of the class, called by the engine to assemble the influence matrices.

        Parameters
        ----------
        mesh1: Mesh or CollectionOfMeshes or list of points
            mesh of the receiving body (where the potential is measured)
            if only S is wanted or early_dot_product is False, then only a list of points as an array of shape (n, 3) can be passed.
        mesh2: Mesh or CollectionOfMeshes
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
        """

        wavenumber = float(wavenumber)

        if free_surface == np.inf: # No free surface, only a single Rankine source term

            a_exp, lamda_exp = np.empty(1), np.empty(1)  # Dummy arrays that won't actually be used by the fortran code.

            coeffs = np.array((1.0, 0.0, 0.0))

        elif water_depth == np.inf:

            a_exp, lamda_exp = np.empty(1), np.empty(1)  # Idem

            if wavenumber == 0.0:
                coeffs = np.array((1.0, 1.0, 0.0))
            elif wavenumber == np.inf:
                coeffs = np.array((1.0, -1.0, 0.0))
            else:
                if self.legacy_wave_part:
                    coeffs = np.array((1.0, -1.0, 1.0))
                else:
                    coeffs = np.array((1.0, 1.0, 1.0))

        else:  # Finite water_depth
            if wavenumber == 0.0 or wavenumber == np.inf:
                raise NotImplementedError("Zero or infinite frequencies not implemented for finite depth.")
            else:
                a_exp, lamda_exp = self.find_best_exponential_decomposition(
                    wavenumber*water_depth*np.tanh(wavenumber*water_depth),
                    wavenumber*water_depth,
                )
                coeffs = np.array((1.0, 1.0, 1.0))

        if isinstance(mesh1, Mesh) or isinstance(mesh1, CollectionOfMeshes):
            collocation_points = mesh1.faces_centers
            nb_collocation_points = mesh1.nb_faces
            if ( adjoint_double_layer == False ):
                early_dot_product_normals = np.zeros((nb_collocation_points, 3))  # Should not be used
            else:
                early_dot_product_normals = mesh1.faces_normals
        elif isinstance(mesh1, np.ndarray) and mesh1.ndim ==2 and mesh1.shape[1] == 3:
            # This is used when computing potential or velocity at given points in postprocessing
            collocation_points = mesh1
            nb_collocation_points = mesh1.shape[0]
            early_dot_product_normals = np.zeros((nb_collocation_points, 3))  # Should not be used
            if ( adjoint_double_layer == False ):
                raise NotImplementedError("Using a list of points as collocation points is not supported in computing adjoint double layer matrices.")
        else:
            raise ValueError(f"Unrecognized input for {self.__class__.__name__}.evaluate")

        if self.floating_point_precision == "float32":
            dtype = "complex64"
        elif self.floating_point_precision == "float64":
            dtype = "complex128"
        else:
            raise NotImplementedError

        S = np.empty((nb_collocation_points, mesh2.nb_faces), order="F", dtype=dtype)
        K = np.empty((nb_collocation_points, mesh2.nb_faces, 1 if early_dot_product else 3), order="F", dtype=dtype)

        # Main call to Fortran code
        self.fortran_core.matrices.build_matrices(
            collocation_points,  early_dot_product_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
            *mesh2.quadrature_points,
            wavenumber, water_depth,
            coeffs,
            self.tabulation_method_index, self.tabulated_r_range, self.tabulated_z_range, self.tabulated_integrals,
            lamda_exp, a_exp,
            mesh1 is mesh2, self.legacy_wave_part, adjoint_double_layer,
            S, K
        )

        if np.any(np.isnan(S)) or np.any(np.isnan(K)):
            raise RuntimeError("Green function returned a NaN in the interaction matrix.\n"
                    "It could be due to overlapping panels.")

        if early_dot_product: K = K.reshape((nb_collocation_points, mesh2.nb_faces))

        return S, K

################################

class XieDelhommeau(Delhommeau):
    """Variant of Nemoh's Green function, more accurate near the free surface.

    Same arguments and methods as :class:`Delhommeau`.
    """

    legacy_wave_part = False
