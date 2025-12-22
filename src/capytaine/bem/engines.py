"""Definition of the methods to build influence matrices, using possibly some sparse structures."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
from abc import ABC, abstractmethod
from typing import Tuple, Union, Optional, Callable

import numpy as np
import scipy.sparse.linalg as ssl

from capytaine.meshes.symmetric import ReflectionSymmetricMesh as OldReflectionSymmetricMesh
from capytaine.new_meshes.symmetric_meshes import ReflectionSymmetricMesh

from capytaine.green_functions.abstract_green_function import AbstractGreenFunction
from capytaine.green_functions.delhommeau import Delhommeau

from capytaine.tools.block_circulant_matrices import (
        BlockCirculantMatrix, lu_decompose, has_been_lu_decomposed,
        MatrixLike, LUDecomposedMatrixLike
        )

LOG = logging.getLogger(__name__)


####################
#  ABSTRACT CLASS  #
####################

class MatrixEngine(ABC):
    """Abstract method to build a matrix."""

    @abstractmethod
    def build_matrices(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer):
        pass

    @abstractmethod
    def build_S_matrix(self, mesh1, mesh2, free_surface, water_depth, wavenumber):
        pass

    @abstractmethod
    def build_fullK_matrix(self, mesh1, mesh2, free_surface, water_depth, wavenumber):
        pass


##################
#  BASIC ENGINE  #
##################

class Counter:
    def __init__(self):
        self.nb_iter = 0

    def __call__(self, *args, **kwargs):
        self.nb_iter += 1


def solve_gmres(A, b):
    LOG.debug(f"Solve with GMRES for {A}.")

    if LOG.isEnabledFor(logging.INFO):
        counter = Counter()
        x, info = ssl.gmres(A, b, atol=1e-6, callback=counter)
        LOG.info(f"End of GMRES after {counter.nb_iter} iterations.")

    else:
        x, info = ssl.gmres(A, b, atol=1e-6)

    if info > 0:
        raise RuntimeError(f"No convergence of the GMRES after {info} iterations.\n"
                            "This can be due to overlapping panels or irregular frequencies.")

    return x


LUDecomposedMatrixOrNot = Union[MatrixLike, LUDecomposedMatrixLike]


class BasicMatrixEngine(MatrixEngine):
    """
    Default matrix engine.

    Features:
        - Caching of the last computed matrices.
        - Supports plane symmetries and nested plane symmetries.
        - Linear solver can be customized. Default is `lu_decomposition` with caching of the LU decomposition.

    Parameters
    ----------
    green_function: AbstractGreenFunction
        the low level implementation used to compute the coefficients of the matrices.
    linear_solver: str or function, optional
        Setting of the numerical solver for linear problems Ax = b.
        It can be set with the name of a preexisting solver
        (available: "lu_decomposition", "lu_decompositon_with_overwrite" and "gmres", the former is the default choice)
        or by passing directly a solver function.
    """

    green_function: AbstractGreenFunction
    _linear_solver: Union[str, Callable]
    last_computed_matrices: Optional[Tuple[MatrixLike, LUDecomposedMatrixOrNot]]

    def __init__(self, *, green_function=None, linear_solver='lu_decomposition'):

        self.green_function = Delhommeau() if green_function is None else green_function

        self._linear_solver = linear_solver

        self.last_computed_inputs = None
        self.last_computed_matrices = None

        self.exportable_settings = {
            'engine': 'BasicMatrixEngine',
            'linear_solver': str(linear_solver),
            **self.green_function.exportable_settings,
        }

    def __str__(self):
        params= [f"green_function={self.green_function}", f"linear_solver={repr(self._linear_solver)}"]
        return f"BasicMatrixEngine({', '.join(params)})"

    def __repr__(self):
        return self.__str__()

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def build_S_matrix(
            self, mesh1, mesh2, free_surface, water_depth, wavenumber
            ) -> np.ndarray:
        """Similar to :code:`build_matrices`, but returning only :math:`S`"""
        # Calls directly evaluate instead of build_matrices because the caching
        # mechanism of build_matrices is not compatible with giving mesh1 as a
        # list of points, but we need that for post-processing
        S, _ = self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber
            )
        return S

    def build_fullK_matrix(
            self, mesh1, mesh2, free_surface, water_depth, wavenumber
            ) -> np.ndarray:
        """Similar to :code:`build_matrices`, but returning only full :math:`K`
        (that is the three components of the gradient, not just the normal one)"""
        # TODO: could use symmetries. In particular for forward, we compute the
        # full velocity on the same mesh so symmetries could be used.
        _, fullK = self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber,
                adjoint_double_layer=True, early_dot_product=False,
            )
        return fullK

    def _build_matrices_with_symmetries(
            self, mesh1, mesh2, *, free_surface, water_depth, wavenumber, adjoint_double_layer
            ) -> Tuple[MatrixLike, MatrixLike]:
        if (isinstance(mesh1, OldReflectionSymmetricMesh)
                and isinstance(mesh2, OldReflectionSymmetricMesh)
                and mesh1.plane == mesh2.plane):

            S_a, K_a = self._build_matrices_with_symmetries(
                mesh1[0], mesh2[0],
                free_surface=free_surface, water_depth=water_depth,
                wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
                )
            S_b, K_b = self._build_matrices_with_symmetries(
                mesh1[0], mesh2[1],
                free_surface=free_surface, water_depth=water_depth,
                wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
                )

            return BlockCirculantMatrix([S_a, S_b]), BlockCirculantMatrix([K_a, K_b])

        elif (isinstance(mesh1, ReflectionSymmetricMesh)
                and isinstance(mesh2, ReflectionSymmetricMesh)
                and mesh1.plane == mesh2.plane):

            S_a, K_a = self._build_matrices_with_symmetries(
                mesh1.half, mesh2.half,
                free_surface=free_surface, water_depth=water_depth,
                wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
                )
            S_b, K_b = self._build_matrices_with_symmetries(
                mesh1.half, mesh2.other_half,
                free_surface=free_surface, water_depth=water_depth,
                wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
                )

            return BlockCirculantMatrix([S_a, S_b]), BlockCirculantMatrix([K_a, K_b])

        else:
            return self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber,
                adjoint_double_layer=adjoint_double_layer, early_dot_product=True,
            )

    def _build_and_cache_matrices_with_symmetries(
            self, mesh1, mesh2, *, free_surface, water_depth, wavenumber, adjoint_double_layer
            ) -> Tuple[MatrixLike, LUDecomposedMatrixOrNot]:
        if (mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer) == self.last_computed_inputs:
            LOG.debug("%s: reading cache.", self.__class__.__name__)
            return self.last_computed_matrices
        else:
            LOG.debug("%s: computing new matrices.", self.__class__.__name__)
            self.last_computed_matrices = None  # Unlink former cached values, so the memory can be freed to compute new matrices.
            S, K = self._build_matrices_with_symmetries(
                    mesh1, mesh2,
                    free_surface=free_surface, water_depth=water_depth,
                    wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
                    )
            self.last_computed_inputs = (mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer)
            self.last_computed_matrices = (S, K)
            return self.last_computed_matrices

    # Main interface for compliance with AbstractGreenFunction interface
    def build_matrices(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=True):
        r"""Build the influence matrices between mesh1 and mesh2.

        Parameters
        ----------
        mesh1: MeshLike or list of points
            mesh of the receiving body (where the potential is measured)
        mesh2: MeshLike
            mesh of the source body (over which the source distribution is integrated)
        free_surface: float
            position of the free surface (default: :math:`z = 0`)
        water_depth: float
            position of the sea bottom (default: :math:`z = -\infty`)
        wavenumber: float
            wavenumber (default: 1.0)
        adjoint_double_layer: bool, optional
            compute double layer for direct method (F) or adjoint double layer for indirect method (T) matrices (default: True)

        Returns
        -------
        tuple of matrix-like (Numpy arrays or BlockCirculantMatrix)
            the matrices :math:`S` and :math:`K`
        """
        return self._build_and_cache_matrices_with_symmetries(
            mesh1, mesh2,
            free_surface=free_surface, water_depth=water_depth,
            wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
        )

    def linear_solver(self, A: LUDecomposedMatrixOrNot, b: np.ndarray) -> np.ndarray:
        """Solve a linear system with left-hand side A and right-hand-side b

        Parameters
        ----------
        A: matrix-like
            Expected to be the second output of `build_matrices`
        b: np.ndarray
            Vector of the correct length

        Returns
        -------
        x: np.ndarray
            Vector such that A@x = b
        """
        if not isinstance(self._linear_solver, str):
            # If not a string, it is expected to be a custom function that can
            # be called to solve the system
            x = self._linear_solver(A, b)

            if not x.shape == b.shape:
                raise ValueError(f"Error in linear solver of {self}: the shape of the output ({x.shape}) "
                                 f"does not match the expected shape ({b.shape})")

            return x

        elif self._linear_solver in ("lu_decomposition", "lu_decomposition_with_overwrite") :
            overwrite_a = (self._linear_solver == "lu_decompositon_with_overwrite")
            if not has_been_lu_decomposed(A):
                luA = lu_decompose(A, overwrite_a=overwrite_a)
                if A is self.last_computed_matrices[1]:
                    # In normal operation of Capytaine, `A` is always the $D$
                    # or $K$ matrix stored in the cache of the solver.
                    # Here we replace the matrix by its LU decomposition in the
                    # cache to avoid doing the decomposition again.
                    self.last_computed_matrices = (self.last_computed_matrices[0], luA)
            else:
                luA: LUDecomposedMatrixLike = A
            return luA.solve(b)

        elif self._linear_solver == "gmres":
            return solve_gmres(A, b)

        else:
            raise NotImplementedError(
                f"Unknown `linear_solver` in BasicMatrixEngine: {self._linear_solver}"
            )

    def compute_ram_estimation(self, problem):
        if isinstance(problem.body.mesh, ReflectionSymmetricMesh):
            if isinstance(problem.body.mesh.half, ReflectionSymmetricMesh):
                # Should not go deeper than that, there is currently only two
                # symmetries available
                symmetry_factor = 0.25
            else:
                symmetry_factor = 0.5
        else:
            symmetry_factor = 1.0


        nb_faces = problem.body.mesh.nb_faces
        nb_matrices = 2
        nb_bytes = 16

        if self._linear_solver == "lu_decomposition":
            nb_matrices += 1

        if self.green_function.floating_point_precision == "float32":
            nb_bytes = 8

        peak_memory = symmetry_factor * nb_faces**2 * nb_matrices * nb_bytes/1e9
        return peak_memory
